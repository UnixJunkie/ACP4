
module A = BatArray
module CLI = Minicli.CLI
module Log = Dolog.Log
module LO = Line_oriented
module V3 = Mmo.V3

open Printf

(* FBR: TODO: turn on clustering of features (Ph4.cluster_features)
   via a CLI option *)

let read_one_ph4_molecule fn =
  LO.with_in_file fn Ph4.read_one_ph4_encoded_molecule

module P = struct
  type t = V3.t
  let dist = V3.dist
end

module BST = Bst.Bisec_tree.Make(P)

(* score already rotated and translated molecule [m] *)
let score_molecule_pose bsts m =
  let error = ref 0.0 in
  A.iteri (fun i feat ->
      let bst = bsts.(i) in
      let cand_coords = Ph4.get_coords_with_feat m feat in
      A.iter (fun xyz ->
          let _nn, dist = BST.nearest_neighbor xyz bst in
          error := !error +. dist
        ) cand_coords
    ) Ph4.all_features;
  !error

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
  let ref_mol = read_one_ph4_molecule ref_fn in
  let _cand_mol = read_one_ph4_molecule cand_fn in
  (* Each ph4 type in A is BST indexed first *)
  let aro_coords = Ph4.get_coords_with_feat ref_mol ARO in
  let hyd_coords = Ph4.get_coords_with_feat ref_mol HYD in
  let hba_coords = Ph4.get_coords_with_feat ref_mol HBA in
  let hbd_coords = Ph4.get_coords_with_feat ref_mol HBD in
  let pos_coords = Ph4.get_coords_with_feat ref_mol POS in
  let neg_coords = Ph4.get_coords_with_feat ref_mol NEG in
  let _aro_bst = BST.(create 1 Two_bands) aro_coords in
  let _hyd_bst = BST.(create 1 Two_bands) hyd_coords in
  let _hba_bst = BST.(create 1 Two_bands) hba_coords in
  let _hbd_bst = BST.(create 1 Two_bands) hbd_coords in
  let _pos_bst = BST.(create 1 Two_bands) pos_coords in
  let _neg_bst = BST.(create 1 Two_bands) neg_coords in
  (* - find the best rotation and translation of B that minimizes error(A, B). *)
  (* To evaluate one rotation and one translation, the error function is: *)
  (* - for each ph4 point in B, add the distance to its nearest neighbor in A *)
  (*   which has the same type into the error sum *)
  (* At the end, output the (rot, trans) w/ the smallest error. *)
  ()
