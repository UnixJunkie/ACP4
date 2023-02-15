(* Copyright (C) 2022, Francois Berenger
   Tsuda laboratory, The University of Tokyo,
   5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.

   ACP4: Encode and compare molecules in the pharmacophore space. *)

open Printf

module A = BatArray
module BA = Bigarray
module BA1 = BA.Array1
module CLI = Minicli.CLI
module Fn = Filename
module Ht = BatHashtbl
module IntMap = BatMap.Int
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module S = BatString

module Bst_point = struct
  type t = string * Common.sparse_float_fp
  let dist = Common.tani_dist
end

module BST = Bst.Bisec_tree.Make(Bst_point)

(* FBR: I should check query performance Vs One_band OR Two_bands *)

(* FBR: this can get silently OOM killed if nprocs > 1 *)
let index_fps nprocs db_fn chunk_size fps =
  let feedback =
    if nprocs = 1 then (fun curr tot -> printf "BST: %d/%d\r%!" curr tot)
    else (fun _ _ -> ()) in
  let size_prfx = chunk_size / 1000 in
  let n = A.length fps in
  (if n >= 1_000_000
   then Log.warn
   else Log.info) "indexing %d molecules" n;
  let index_chunk (curr_chunk_size, j, bst_j) =
    let bst_fn = sprintf "%s.bst%dk.%d" db_fn size_prfx bst_j in
    Log.info "size: %d creating %s" curr_chunk_size bst_fn;
    (* chunk feedback OFF in parallel case *)
    let bst =
      let to_index = A.sub fps j curr_chunk_size in
      (* Warning: Two_bands is too slow w/ 100k chunks *)
      BST.create
        ~progress_callback:feedback 1 Bst.Bisec_tree.One_band to_index in
    LO.with_out_file bst_fn (fun out ->
        Marshal.(to_channel out bst [No_sharing])
      ) in
  Parany.run nprocs
    ~demux:(
      let i = ref 0 in
      let bst_i = ref 0 in
      fun () ->
        if !i < n then
          let size = min chunk_size (n - !i) in
          let res = (size, !i, !bst_i) in
          i := !i + size;
          incr bst_i;
          res
        else
          let () = assert(n = !i) in
          raise Parany.End_of_input)
    ~work:index_chunk
    ~mux:(fun () -> ())

module Ligand_defaults = struct
  let radial_cutoff = 5.0
  let dx = 0.5
  let bst_chunk_size = 100_000
end

(* protein Binding-Sites *)
module BS_defaults = struct
  let radial_cutoff = 40.0
  let dx = 0.9
end

type mode =
  (* default: ph4 to liblinear classification format *)
  | Encode of string
  (* (query_fn, db_fn) -q and -db CLI options *)
  | Single_query of string * string
  (* (query_fns, db_fn) --queries and -db CLI options *)
  | Multiple_queries of string * string

let decode_mode maybe_input_fn maybe_db_fn maybe_query_fn maybe_queries_fn =
  match (maybe_input_fn, maybe_db_fn, maybe_query_fn, maybe_queries_fn) with
  | (Some input_fn, None, None, None) -> Encode input_fn
  | (None, Some db_fn, Some query_fn, None) -> Single_query (query_fn, db_fn)
  | (None, Some db_fn, None, Some queries_fn) ->
    Multiple_queries (queries_fn, db_fn)
  | _ -> failwith "Acp4.decode_mode: unsupported mode: \
                   use EITHER -i OR (-q and -db) OR (--queries and -db)"

let multiple_queries_asked = function
  | Multiple_queries _ -> true
  | _ -> false

let parse_score_line line =
  try Scanf.sscanf line "%f %s" (fun score name -> (score, name))
  with exn ->
    let () = Log.fatal "Acp4.parse_score_line: cannot parse: '%s'" line in
    raise exn

let multi_conf_post_proc_scores_fn scores_fn =
  Log.info "resolving multi-conformer scores file";
  let num_scores = LO.count scores_fn in
  let ht = Ht.create num_scores in
  LO.iter scores_fn (fun line ->
      let score, name_confid = parse_score_line line in
      let last_underscore_i = S.rfind name_confid "_" in
      let name = S.sub name_confid 0 last_underscore_i in
      try
        let prev_score = Ht.find ht name in
        if score > prev_score then
          Ht.replace ht name score
      with Not_found ->
        Ht.add ht name score
    );
  Log.warn "decr. sort of scores file";
  let key_values_a = A.of_list (Ht.to_list ht) in
  A.sort (fun (_n1, s1) (_n2, s2) -> BatFloat.compare s2 s1) key_values_a;
  LO.with_out_file scores_fn (fun out ->
      A.iter (fun (name, score) ->
          fprintf out "%f %s\n" score name
        ) key_values_a
    )

(* query using Tanimoto distance BSTs stored on disk *)
let distance_query nprocs output_fn bst_fns dist query =
  LO.with_out_file output_fn (fun out ->
      Parany.run nprocs
        ~demux:(
          let todo = ref bst_fns in
          fun () -> match !todo with
            | [] -> raise Parany.End_of_input
            | x :: xs -> (todo := xs; x)
        )
        ~work:(fun bst_fn ->
            Log.info "chunk: %s" bst_fn;
            LO.with_in_file bst_fn (fun input ->
                let (bst: BST.t) = Marshal.from_channel input in
                BST.neighbors query dist bst
              )
          )
        ~mux:(L.iter (fun (name, _fp) ->
            fprintf out "%s\n" name)
          )
    )

let main () =
  let start = Unix.gettimeofday () in
  Log.(set_prefix_builder short_prefix_builder);
  Log.color_on ();
  Log.(set_log_level INFO);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              [-i <filename.ph4>]: input file to encode\n  \
              [-q <filename.ph4>]: query molecule\n  \
              [--queries <filename.ph4>]: query molecules\n  \
              [-db <filename.ph4>]: database to screen\n  \
              -o <filename.{csv|scores}>: output file \
              (encoded db or query scores)\n  \
              [-np <int>]: nprocs (default=1)\n  \
              [-s <int>]: BST chunk size (default=%d)\n  \
              [-c <float>]: cutoff distance (default=%.2f)\n  \
              [-csize <int>]: chunk size (for better parallelization)\n  \
              [-f]: force overwriting binary cache files, if any\n  \
              use -f if you are changing -dx or -c compared to previous \
              queries\n  \
              [-dx <float>]: radial discretization step (default=%g)\n  \
              [--no-plot]: turn OFF gnuplot\n  \
              [--no-tap]: \
              neither read nor write from/to binary cache files\n  \
              [--confs]: \
              if the DB has several _consecutive_ conformers\n  \
              [--quick]: \
              quick; exit right after scoring; no perf. metrics\n  \
              [--index]: \
              create a BST index for the DB being screened\n  \
              [--BSTs <filename.txt>]: \
              file containing a list of serialized BST filenames\n  \
              [-td <float>]: \
              maximum Tanimoto _distance_ to (single) query\n  \
              [--BS]: \
              load optimal defaults for binding-sites (ignores -c and -dx)\n  \
              [-v]: verbose/debug mode\n"
       Sys.argv.(0)
       Ligand_defaults.bst_chunk_size
       Ligand_defaults.radial_cutoff
       Ligand_defaults.dx;
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  if verbose then Log.(set_log_level DEBUG);
  let maybe_input_fn = CLI.get_string_opt ["-i"] args in
  let maybe_db_fn = CLI.get_string_opt ["-db"] args in
  let maybe_query_fn = CLI.get_string_opt ["-q"] args in
  let maybe_queries_fn = CLI.get_string_opt ["--queries"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let force = CLI.get_set_bool ["-f"] args in
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  let csize = CLI.get_int_def ["-csize"] args 1 in
  let bst_chunk_size =
    CLI.get_int_def ["-s"] args Ligand_defaults.bst_chunk_size in
  let binding_site_mode = CLI.get_set_bool ["--BS"] args in
  let cutoff =
    CLI.get_float_def ["-c"] args
      (if binding_site_mode
       then BS_defaults.radial_cutoff
       else Ligand_defaults.radial_cutoff) in
  let dx =
    CLI.get_float_def ["-dx"] args
      (if binding_site_mode
       then BS_defaults.dx
       else Ligand_defaults.dx) in
  let no_tap = CLI.get_set_bool ["--no-tap"] args in
  let multi_conf = CLI.get_set_bool ["--confs"] args in
  let quick = CLI.get_set_bool ["--quick"] args in
  let index = CLI.get_set_bool ["--index"] args in
  let maybe_bst_fn = CLI.get_string_opt ["--BSTs"] args in
  let maybe_dist = match CLI.get_float_opt ["-td"] args with
    | None -> None
    | Some td ->
      assert(td >= 0.0 && td < 1.0);
      Some td in
  let maybe_roc_curve_fn =
    if CLI.get_set_bool ["--no-plot"] args then None
    else Some (Fn.temp_file ~temp_dir:"/tmp" "acp4_" ".roc") in
  CLI.finalize (); (* -------------------------------------------------- *)
  (if multi_conf then
     if output_fn = "/dev/null" then
       (Log.fatal "--confs needs an actual output file; not /dev/null";
        exit 1)
  );
  let selected_mode =
    decode_mode maybe_input_fn maybe_db_fn maybe_query_fn maybe_queries_fn in
  let mol_count = ref 0 in
  let nb_dx = 1 + BatFloat.round_to_int (cutoff /. dx) in
  (match selected_mode with
   | Encode input_fn ->
     LO.with_infile_outfile input_fn output_fn (fun input output ->
         Parany.run ~preserve:true nprocs ~csize:csize
           ~demux:(fun () ->
               try Ph4.preread_one_ph4_encoded_molecule input
               with End_of_file -> raise Parany.End_of_input)
           ~work:(fun name_lines ->
               let encoded =
                 Ph4.encode verbose cutoff dx
                   (Ph4.parse_one_ph4_encoded_molecule name_lines) in
               Common.sprintf_code nb_dx encoded)
           ~mux:(fun line ->
               output_string output line;
               incr mol_count)
       )
   | Single_query (queries_fn, db_fn)
   | Multiple_queries (queries_fn, db_fn) ->
     let maybe_tap_out, maybe_tap_in =
       if no_tap then
         (None, None)
       else
         let bin_tap_fn = db_fn ^ ".bin" in
         if not force && Sys.file_exists bin_tap_fn then
           let () = Log.warn "%s exists" bin_tap_fn in
           (None, Some (open_in_bin bin_tap_fn))
         else
           (Some (open_out_bin bin_tap_fn), None) in
     let queries =
       assert(S.ends_with queries_fn ".ph4");
       A.of_list (Common.parse_all verbose cutoff dx nb_dx queries_fn) in
     let num_queries = A.length queries in
     Log.info "%d queries in %s" num_queries queries_fn;
     assert(num_queries >= 1);
     (if multiple_queries_asked selected_mode && num_queries = 1 then
        Log.warn "single query in %s" queries_fn
     );
     begin match (maybe_bst_fn, maybe_dist) with
       | (Some fn, Some dist) ->
         begin
           if num_queries > 1 then
             Log.warn "distance query only uses 1st query"
           ;
           let bst_fns = LO.lines_of_file fn in
           distance_query nprocs output_fn bst_fns dist queries.(0)
         end
       | _ ->
         begin
           (* assert(A.for_all (fun x -> is_active (fst x)) queries);
              (\* active labels *\) *)
           assert(S.ends_with db_fn ".ph4");
           LO.with_infile_outfile db_fn output_fn (fun input scores_out ->
               let reader = match maybe_tap_in with
                 | None -> (fun () -> Common.parse_one verbose cutoff dx nb_dx input)
                 | Some bin_tap_in -> (fun () -> Common.read_one_bin bin_tap_in) in
               let writer = match maybe_tap_out with
                 | None -> (fun _x -> ())
                 | Some bin_out -> Common.write_one_bin bin_out in
               let name_scores = ref [] in
               let to_index = ref [] in
               try
                 while true do
                   let (cand_name, _cand_fp) as cand = reader () in
                   (if index then
                      to_index := cand :: !to_index
                   );
                   writer cand;
                   let scores = A.map (Common.tanimoto cand) queries in
                   let score = A.max scores in (* MAX consensus *)
                   name_scores := (Common.is_active cand_name, score) :: !name_scores;
                   fprintf scores_out "%f %s\n" score cand_name;
                   incr mol_count
                 done
               with End_of_file ->
                 begin
                   BatOption.may close_out maybe_tap_out;
                   BatOption.may close_in maybe_tap_in;
                   (* put them back in the right order *)
                   name_scores := L.rev !name_scores;
                   let score_labels_a = A.of_list !name_scores in
                   name_scores := [];
                   let num_entries = A.length score_labels_a in
                   Log.info "%d molecules in %s" num_entries db_fn;
                   (if !to_index <> [] then
                      let fps = A.of_list !to_index in
                      to_index := [];
                      index_fps nprocs db_fn bst_chunk_size fps
                   );
                   if not quick && not multi_conf then
                     Common.performance_metrics
                       maybe_roc_curve_fn db_fn score_labels_a cutoff dx
                 end
             );
           if not quick && multi_conf then
             begin
               multi_conf_post_proc_scores_fn output_fn;
               let score_labels_a =
                 let score_labels = LO.map output_fn parse_score_line in
                 A.map (fun (score, name) -> (Common.is_active name, score))
                   (A.of_list score_labels) in
               Common.performance_metrics
                 maybe_roc_curve_fn db_fn score_labels_a cutoff dx
             end
         end
     end
  );
  let max_feat = 1 + (Ph4.nb_channels * nb_dx) in
  let now = Unix.gettimeofday () in
  let dt = now -. start in
  let freq = (float !mol_count) /. dt in
  Log.info "%d molecules @ %.2f Hz max_feat=%d" !mol_count freq max_feat

let () = main ()
