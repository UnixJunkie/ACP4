
module CLI = Minicli.CLI
module Log = Dolog.Log
module LO = Line_oriented

open Printf

(* FBR: TODO: turn on clustering of features (Ph4.cluster_features)
   via a CLI option *)

let read_one_ph4_molecule fn =
  LO.with_in_file fn Ph4.read_one_ph4_encoded_molecule

let main () =
  Log.(set_prefix_builder short_prefix_builder);
  Log.color_on ();
  Log.(set_log_level INFO);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              [-ref <filename.ph4>]: reference structure\n  \
              [-cand <filename.ph4>]: candidate structure\n  \
              [-v]: verbose/debug mode\n"
       Sys.argv.(0);
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  if verbose then Log.(set_log_level DEBUG);
  (* given two .ph4 files, assume: *)
  (* - LHS = A: the reference structure *)
  (* - RHS = B: the candidate structure *)
  let ref_fn = CLI.get_string ["-ref"] args in
  let cand_fn = CLI.get_string ["-cand"] args in
  CLI.finalize (); (* -------------------------------------------------- *)
  let _ref_mol = read_one_ph4_molecule ref_fn in
  let _cand_mol = read_one_ph4_molecule cand_fn in
  (* Each ph4 type in A is BST indexed first. *)

  (* - find the best rotation and translation of B that minimizes error(A, B). *)
  (* To evaluate one rotation and one translation, the error function is: *)
  (* - for each ph4 point in B, add the distance to its nearest neighbor in A *)
  (*   which has the same type into the error sum *)
  (* At the end, output the (rot, trans) w/ the smallest error. *)
  ()
