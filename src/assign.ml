(* Find a good assignment between two sets of ph4 features *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module V3 = Vector3

(* small constant to avoid division by 0.0 *)
let epsilon = 0.001

let all_distinct_pairs l =
  let rec loop acc = function
    | [] -> acc
    | x :: ys ->
      let acc' = L.rev_append (L.rev_map (fun y -> (x, y)) ys) acc in
      loop acc' ys in
  loop [] l

(* Heuristic: a good assignment means a maximum of points have been assigned
   a correspondant in the other molecule
   AND the sum of internal distance discrepancies is minimal *)
let score_assignment m1 m2 pairs' =
  (* y = -1 <=> x was not assigned *)
  let pairs = L.filter (fun (_x, y) -> y <> -1) pairs' in
  (* maximize number of assignments *)
  let n = L.length pairs in
  let assign1, assign2 = L.split pairs in
  let error = ref epsilon in
  let coords1 = Ph4.get_coords m1 in
  let coords2 = Ph4.get_coords m2 in
  (* compute all internal dists in m1 and m2 *)
  let pairs1 = all_distinct_pairs assign1 in
  let pairs2 = all_distinct_pairs assign2 in
  (* FBR: we need to use a distance cache in there
   *      it can be initialized the first time molecules are read in *)
  let dists1 =
    L.rev_map (fun (i, j) -> V3.dist coords1.(i) coords1.(j)) pairs1 in
  let dists2 =
    L.rev_map (fun (i, j) -> V3.dist coords2.(i) coords2.(j)) pairs2 in
  (* minimize "disagreement over distances" *)
  L.iter2 (fun d1 d2 ->
      error := !error +. abs_float (d1 -. d2)
    ) dists1 dists2;
  (n, !error)

(* which feature indexes in [m1] are compatible with
   which feature indexes in [m2] *)
let compatible_feats m1 m2 =
  let all_feats = Ph4.[ARO; HYD; HBA; HBD; POS; NEG] in
  L.fold_left (fun acc feat ->
      let xs = ref [] in
      let m1_feats = Ph4.get_features m1 in
      A.iteri (fun i x ->
          if feat = x then
            xs := i :: !xs
        ) m1_feats;
      let ys = ref [] in
      let m2_feats = Ph4.get_features m2 in
      A.iteri (fun i x ->
          if feat = x then
            ys := i :: !ys
        ) m2_feats;
      let res = (!xs, !ys) in
      if res <> ([], []) then
        res :: acc
      else
        acc
    ) [] all_feats

(* if [i] is assigned to [j] (i, j), then no one else can be assigned to [j] *)
let assign_only_once l =
  let rec loop acc = function
    | [] -> acc
    | ((i: int), (j: int)) :: rest ->
      let acc' =
        let l' = (i, j) :: (L.filter (fun (_i', j') -> j <> j') rest) in
        L.rev_append l' acc in
      loop acc' rest in
  loop [] l

let all_possible_assignments compats =
  let intermediate =
    L.fold_left (fun acc (feats1, feats2) ->
        let prod = L.cartesian_product feats1 feats2 in
        (assign_only_once prod) :: acc
      ) [] compats in
  failwith "not implemented yet"

let string_of_list to_str l =
  let strings = L.map to_str l in
  let body = String.concat "; " strings in
  String.concat "" ["["; body; "]"]

let string_of_int_pair (x, y) =
  sprintf "(%d, %d)" x y

let main () =
  Log.(set_prefix_builder short_prefix_builder);
  Log.color_on ();
  Log.(set_log_level INFO);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              -r <ref.ph4>: reference set of ph4 features\n  \
              -c <cand.ph4>: candidate set of ph4 features\n  \
              [-v]: verbose/debug mode\n"
       Sys.argv.(0);
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  if verbose then Log.(set_log_level DEBUG);
  let ref_fn = CLI.get_string ["-r"] args in
  let cand_fn = CLI.get_string ["-c"] args in
  CLI.finalize (); (* ----------------------------------------------------- *)
  let reference = LO.with_in_file ref_fn  Ph4.read_one_ph4_encoded_molecule in
  let candidate = LO.with_in_file cand_fn Ph4.read_one_ph4_encoded_molecule in
  Log.info "|ref|=%d" (Ph4.num_features reference);
  Log.info "|cand|=%d" (Ph4.num_features candidate);
  (* test assignments *)
  let compats = compatible_feats reference candidate in
  L.iter (fun (xs, ys) ->
      Log.info "%s - %s"
        (string_of_list string_of_int xs)
        (string_of_list string_of_int ys)
    ) compats;
  let all_assignments = all_possible_assignments compats in
  Log.warn "possible assignments: %d" (L.length all_assignments);
  L.iter (fun x ->
      let n, error = score_assignment reference candidate x in
      let score = (float n) /. error in
      Log.info "score: %f N: %d %s"
        score n (string_of_list string_of_int_pair x)
    ) all_assignments

let () = main ()
