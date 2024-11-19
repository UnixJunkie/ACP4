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

(* small constant to avoid division by 0.0 *)
let epsilon = 0.001

let all_distinct_pairs l =
  let rec loop acc = function
    | [] -> acc
    | x :: ys ->
      let acc' = L.rev_append (L.rev_map (fun y -> (x, y)) ys) acc in
      loop acc' ys in
  loop [] l

(* intuitively, a good assignment means a maximum of points have been assigned
   a correspondant in the other molecule AND the sum of distance errors is
   minimal *)
let score_assignment m1 m2 pairs =
  let n = L.length pairs in
  let assign1, assign2 = L.split pairs in
  let error = ref epsilon in
  let coords1 = Ph4.get_coords m1 in
  let coords2 = Ph4.get_coords m2 in
  (* compute all internal dists in m1 *)
  
  (* compute all internal dists in m2 *)
  (* measure the "disagreement error" *)
  (float n) /. !error

(* for each feature in [m1], find all possible assignments from [m2];
   including no assignment (-1);
   the smaller molecule is supposed to be [m1] *)
let possible_assignments m1 m2 =
  let n1 = Ph4.num_features m1 in
  let n2 = Ph4.num_features m2 in
  assert(n1 <= n2);
  let res = ref [] in
  for i = 0 to n1 - 1 do
    let feat1 = m1.features.(i) in
    (* no match is always a possibility *)
    let i_matches = ref [-1] in
    for j = 0 to n2 - 1 do
      let feat2 = m2.features.(j) in
      if feat1 = feat2 then
        i_matches := j :: !i_matches
    done;
    res := (i, !i_matches) :: !res
  done;
  !res

(* let rec enumerate_assignments acc m1 m2 m2_unassigned possibilities = *)
(*   match possibilities with *)
(*   | [] -> acc *)
(*   | x :: xs -> *)
(*     begin *)
(*       match x with *)
(*       | (i, js) -> *)
(*         L.map (fun j -> *)
(*             enumerate_assignments () *)

(*     end *)

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
