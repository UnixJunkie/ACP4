(* Copyright (C) 2022, Francois Berenger
   Tsuda laboratory, The University of Tokyo,
   5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.

   Extract ligand-defined binding site. *)

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

module Point_3D = struct
  type t = V3.t
  let dist = V3.dist
end

module BST = Bst.Bisec_tree.Make(Point_3D)

(* carve ligand-defined binding-site out of the protein receptor *)
let carve ligand_bst cutoff protein =
  let feats = ref [] in
  let coords = ref [] in
  let open Ph4 in
  A.iter2 (fun feat xyz ->
      let _nearest, dist = BST.nearest_neighbor xyz ligand_bst in
      if dist <= cutoff then
        (feats := feat :: !feats;
         coords := xyz :: !coords)
    ) protein.features protein.coords;
  { protein with features = A.of_list (L.rev !feats);
                 coords = A.of_list (L.rev !coords) }

let main () =
  Log.(set_prefix_builder short_prefix_builder);
  Log.color_on ();
  Log.(set_log_level INFO);
  let default_cutoff = 5.0 in (* (Angstrom) around ligand's heavy atoms *)
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              -l <ligand.sdf>: binding-site ligand input file\n  \
              -p <protein.ph4>: receptor protein input file\n  \
              [-d <float>]: distance cutoff (default=%.2f)\n  \
              -o <output.ph4>: ligand-defined binding site output file\n  \
              [-v]: verbose/debug mode\n"
       Sys.argv.(0) default_cutoff;
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  if verbose then Log.(set_log_level DEBUG);
  let ligand_fn = CLI.get_string ["-l"] args in
  let protein_fn = CLI.get_string ["-p"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let cutoff = CLI.get_float_def ["-d"] args default_cutoff in
  CLI.finalize (); (* ----------------------------------------------------- *)
  LO.with_out_file output_fn (fun output ->
      (* read 3D ligand from SDF *)
      let ligand = LO.with_in_file ligand_fn Sdf_3D.read_one_molecule in
      (* keep only heavy atoms *)
      let heavies_3D = Sdf_3D.heavy_atom_coords ligand in
      (* index ligand using a BST *)
      let bst = BST.(create 1 Two_bands heavies_3D) in
      (* read protein from a .ph4 file *)
      let protein =
        LO.with_in_file protein_fn Ph4.read_one_ph4_encoded_molecule in
      (* prune points which are further than cutoff distance
         from any HA ligand atom *)
      let binding_site = carve bst cutoff protein in
      (* output pruned result to file *)
      Ph4.write_out output binding_site
    )

let () = main ()
