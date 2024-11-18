(* Find a good assignment between two sets of ph4 features *)

open Printf

module A = BatArray
module BA = Bigarray
module BA1 = BA.Array1
module CLI = Minicli.CLI
module Fn = Filename
module IntMap = BatMap.Int
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module V3 = Vector3

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
  let reference = LO.with_in_file ref_fn Ph4.read_one_ph4_encoded_molecule in
  let candidate = LO.with_in_file cand_fn Ph4.read_one_ph4_encoded_molecule in
  Log.info "|ref|=%d" (Ph4.num_features reference);
  Log.info "|cand|=%d" (Ph4.num_features candidate)

let () = main ()
