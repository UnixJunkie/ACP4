(* Copyright (C) 2023, Francois Berenger
   Tsuda laboratory, The University of Tokyo,
   5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.

   Minimum Spanning Tree (MST) over a dataset's Gram matrix.
   Nodes are connected all to all (undirected graph).
   The edge weight between two nodes is the Tanimoto distance between them.

   Probst, D., & Reymond, J. L. (2020).
   Visualization of very large high-dimensional data sets as minimum spanning
   trees. Journal of Cheminformatics, 12(1), 1-13. *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module Ht = Hashtbl
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log

(* Undirected graph with integer labels on vertices and float labels (weights)
   on edges *)
module Node = struct
  type t = int
  let compare = BatInt.compare
  let hash = Ht.hash
  let equal = BatInt.equal
  let default = 0
end

module Edge = struct
  type t = float
  let compare = BatFloat.compare
  let hash = Ht.hash
  let equal = BatFloat.equal
  let default = 0.0
end

(* JCF said we should use Graph.Imperative.Matrix.Graph for performance reasons
   because our graphs are dense
   JS said we should use Prim's algorithm instead of Kruskal; also for performance
   reasons since we have a dense graph *)
module G = Graph.Imperative.Graph.ConcreteLabeled(Node)(Edge)

module W = struct
  type label = G.E.label
  type edge = G.E.t
  type t = float (* mandatory *)
  let weight x = G.E.label x
  let zero = 0.0
  let add = (+.)
  let compare = BatFloat.compare (* mandatory *)
end

module Kruskal = Graph.Kruskal.Make(G)(W)

(* write graph to file in graphviz dot format *)
let graph_to_dot fn g =
  LO.with_out_file fn (fun out ->
      fprintf out "graph all_to_all {\n";
      G.iter_edges_e (fun e ->
          fprintf out "%d -- %d [label=\"%.2f\"]\n"
            (G.E.src e) (G.E.dst e) (G.E.label e)
        ) g;
      fprintf out "}\n";
    )

(* write the MST edges to file in dot format *)
let mst_edges_to_dot fn names edges =
  LO.with_out_file fn (fun out ->
      fprintf out "graph min_span_tree {\n";
      let nb_nodes = A.length names in
      for i = 0 to nb_nodes - 1 do
        (* let ic50 = pIC50s.(i) in *)
        let name = names.(i) in
        (* node color *)
        (* FBR: maybe later; use one color per EC class? *)
        (* let red, green, blue = rgb_triplet min_pIC50 delta_pIC50 ic50 in *)
        let red, green, blue = 255, 255, 255 in
        assert(0 <= red && red <= 255 &&
               0 <= green && green <= 255 &&
               0 <= blue && blue <= 255);
        fprintf out "\"%d\" [label=\"\" style=\"filled\" \
                     color=\"#%02x%02x%02x\" \
                     image=\"pix/%s_lig.png\"]\n"
          i red green blue name
      done;
      L.iter (fun e ->
          fprintf out "%d -- %d [label=\"%.2f\"]\n"
            (G.E.src e) (G.E.dst e) (G.E.label e)
        ) edges;
      fprintf out "}\n"
    )

let minimum_spanning_tree g =
  Kruskal.spanningtree g

(* mean and stddev from n samples Tanimoto *)
let tanimoto_mean_std n a =
  let rng = Random.State.make_self_init () in
  let len = A.length a in
  let arr =
    A.init n (fun _i ->
        let i = Random.State.int rng len in
        let j = Random.State.int rng len in
        Common.tanimoto' a.(i) a.(j)
      ) in
  Common.(average arr, stddev arr)

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s\n  \
              -i <filename.ph4>: input file\n  \
              -o <filename.dot>: output MST to dot file\n  \
              [-go <filename.dot>]: output fully connected graph to dot file\n  \
              [-np <int>]: maximum number of CPU cores (default=1)\n  \
              [-cs <int>]: parallel job chunk size (default=1)\n  \
              [--BS]: optimal defaults for binding-sites (ignores -c and -dx)\n  \
              [-c <float>]: cutoff distance (default=%.2f)\n  \
              [-dx <float>]: radial discretization step (default=%g)\n  \
              [-v]: verbose/debug mode\n"
       Sys.argv.(0)
       Common.Ligand_defaults.radial_cutoff
       Common.Ligand_defaults.dx;
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let maybe_full_graph_fn = CLI.get_string_opt ["-go"] args in
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  let csize = CLI.get_int_def ["-cs"] args 1 in
  let verbose = CLI.get_set_bool ["-v"] args in
  let binding_site_mode = CLI.get_set_bool ["--BS"] args in
  let cutoff =
    CLI.get_float_def ["-c"] args
      (if binding_site_mode
       then Common.BS_defaults.radial_cutoff
       else Common.Ligand_defaults.radial_cutoff) in
  let dx =
    CLI.get_float_def ["-dx"] args
      (if binding_site_mode
       then Common.BS_defaults.dx
       else Common.Ligand_defaults.dx) in
  CLI.finalize (); (* ------------------------------------------------------ *)
  let nb_dx = 1 + BatFloat.round_to_int (cutoff /. dx) in
  Log.info "reading molecules...";
  let names, all_mols =
    A.split (A.of_list (Common.parse_all verbose cutoff dx nb_dx input_fn)) in
  let nb_mols = A.length all_mols in
  Log.info "read %d" nb_mols;
  let g = G.create ~size:nb_mols () in
  (* add all nodes to graph *)
  for i = 0 to nb_mols - 1 do
    G.add_vertex g i
  done;
  (* compute Gram matrix in // *)
  let matrix = A.make_matrix nb_mols nb_mols 0.0 in
  Log.info "Gram matrix initialization...";
  Molenc.Gram.initialize_matrix Common.tani_dist' nprocs csize all_mols matrix;
  Molenc.Gram.print_corners matrix;
  Log.info "Adding edges to graph...";
  (* add all edges to graph *)
  let disconnected = ref 0 in
  for i = 0 to nb_mols - 1 do
    (* WARNING: we don't initialize the diagonal
       (it is supposed to be all 0s) *)
    for j = i + 1 to nb_mols - 1 do
      let w = matrix.(i).(j) in
      if w = 1.0 then
        incr disconnected
      else
        let edge = G.E.create i w j in
        G.add_edge_e g edge
    done;
    printf "done: %d/%d\r%!" (i + 1) nb_mols;
  done;
  printf "\n%!";
  (if !disconnected > 0 then
     Log.info "disconnected molecules: %d" !disconnected);
  BatOption.may (fun fn -> graph_to_dot fn g) maybe_full_graph_fn;
  (* MST *)
  Log.info "MST...";
  let mst = minimum_spanning_tree g in
  (* dump to file *)
  Log.info "writing file %s ..." output_fn;
  mst_edges_to_dot output_fn names mst

let () = main ()
