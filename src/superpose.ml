
module A = BatArray
module CLI = Minicli.CLI
module Log = Dolog.Log
module LO = Line_oriented
module V3 = Mmo.V3
module Rot = Mmo.Rot

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
(* To evaluate one pose (rotation-translation), the error function is: *)
(* - for each ph4 point in B, add the distance to its nearest neighbor in A *)
(*   which has the same type into the total error *)
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

let geometric_center vects =
  V3.div
    (A.fold_left V3.add V3.origin vects)
    (float (A.length vects))

let translate_by vects delta =
  let dx, dy, dz = V3.to_triplet delta in
  A.map (fun V3.{x; y; z} ->
      V3.make (x +. dx) (y +. dy) (z +. dz)
    ) vects

let translate_to vects p =
  translate_by vects (V3.diff p (geometric_center vects))

let center_molecule m =
  Ph4.{ m with coords = translate_to m.coords V3.origin }

(* CONSTRAINT: vects must already be centered at origin *)
let centered_rotate vects r =
  A.map (Rot.rotate r) vects

(* rotate then translate already centered molecule [centered_m] *)
let place_molecule centered_m params =
  let rotated =
    centered_rotate
      Ph4.(centered_m.coords)
      (Rot.r_xyz params.(0) params.(1) params.(2)) in
  Ph4.{ centered_m with
        coords =
          translate_to rotated (V3.make params.(3) params.(4) params.(5)) }

(* we don't have a gradient --> _gradient *)
let nlopt_eval_solution verbose bsts centered_m params _gradient =
  let posed = place_molecule centered_m params in
  let error = score_molecule_pose bsts posed in
  (if verbose then
     Log.info "err: %f" error
  );
  error

let two_pi = Mmo.Math.two_pi

let nlopt_optimize verbose bsts centered_m max_evals =
  let ndims = 6 in (* DOFs: 3 rotations + 3 translations *)
  (* local optimizer that will be passed to the global one *)
  let local = Nlopt.(create sbplx ndims) in (* local optimizer: gradient-free *)
  Nlopt.set_min_objective local (nlopt_eval_solution verbose bsts centered_m);
  (* I don't set parameter bounds on the local optimizer, I guess
   * the global optimizer handles this *)
  (* hard stop conditions *)
  Nlopt.set_stopval local 0.0; (* absolute minimum total error *)
  (* smart stop conditions *)
  Nlopt.set_ftol_abs local 0.0001; (* might need to be tweaked *)
  let global = Nlopt.(create auglag ndims) in (* global optimizer *)
  Nlopt.set_local_optimizer global local;
  Nlopt.set_min_objective global (nlopt_eval_solution verbose bsts centered_m);
  (* bounds for parameters *)
  (* FBR: found tx_min, ty_min, tz_min using bounding boxes? *)
  (*                              rz      ry      rz       tx       ty       tz *)
  Nlopt.set_lower_bounds global [|0.0;    0.0;    0.0;    -1000.0; -1000.0; -1000.0|];
  Nlopt.set_upper_bounds global [|two_pi; two_pi; two_pi;  1000.0;  1000.0;  1000.0|];
  (* hard stop conditions *)
  Nlopt.set_stopval global 0.0;
  (* max number of iterations *)
  Nlopt.set_maxeval global max_evals;
  (* starting solution *)
  let start_sol = Array.make ndims 0.0 in (* no rot., no trans. *)
  let stop_cond, params, error = Nlopt.optimize global start_sol in
  Log.info "NLopt optimize global: %s" (Nlopt.string_of_result stop_cond);
  Log.info "final error: %f" error;
  params

let random_start m =
  Random.self_init ();
  place_molecule m [|Random.float two_pi;
                     Random.float two_pi;
                     Random.float two_pi;
                     -1000.0 +. Random.float 2000.;
                     -1000.0 +. Random.float 2000.;
                     -1000.0 +. Random.float 2000.|]

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
              [-o <filename.ph4>]: optimally superposed cand. structure\n  \
              [--rand-start]: randomize starting position of cand. before\n  \
              starting optimization (for tests).\n  \
              [-v]: verbose/debug mode\n"
       Sys.argv.(0);
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  if verbose then Log.(set_log_level DEBUG);
  let rand_start = CLI.get_set_bool ["--rand-start"] args in
  (* given two .ph4 files, assume: *)
  (* - LHS = A: the reference structure *)
  (* - RHS = B: the candidate structure *)
  let ref_fn = CLI.get_string ["-ref"] args in
  let cand_fn = CLI.get_string ["-cand"] args in
  let output_fn = CLI.get_string ["-o"] args in
  CLI.finalize (); (* -------------------------------------------------- *)
  let ref_mol = read_one_ph4_molecule ref_fn in
  let cand_mol =
    center_molecule
      (if rand_start then
         random_start (read_one_ph4_molecule cand_fn)         
       else
         read_one_ph4_molecule cand_fn) in
  (* Each ph4 type in A is BST indexed first *)
  let aro_coords = Ph4.get_coords_with_feat ref_mol ARO in
  let hyd_coords = Ph4.get_coords_with_feat ref_mol HYD in
  let hba_coords = Ph4.get_coords_with_feat ref_mol HBA in
  let hbd_coords = Ph4.get_coords_with_feat ref_mol HBD in
  let pos_coords = Ph4.get_coords_with_feat ref_mol POS in
  let neg_coords = Ph4.get_coords_with_feat ref_mol NEG in
  let aro_bst = BST.(create 1 Two_bands) aro_coords in
  let hyd_bst = BST.(create 1 Two_bands) hyd_coords in
  let hba_bst = BST.(create 1 Two_bands) hba_coords in
  let hbd_bst = BST.(create 1 Two_bands) hbd_coords in
  let pos_bst = BST.(create 1 Two_bands) pos_coords in
  let neg_bst = BST.(create 1 Two_bands) neg_coords in
  let bsts = [|aro_bst; hyd_bst; hba_bst; hbd_bst; pos_bst; neg_bst|] in
  (* - find the best rotation and translation of B that minimizes error(A, B). *)
  (* At the end, output the (rot, trans) w/ the smallest error. *)
  let best_params = nlopt_optimize verbose bsts cand_mol 1000 in
  Log.info "best_params: %f %f %f %f %f %f"
    best_params.(0)
    best_params.(1)
    best_params.(2)
    best_params.(3)
    best_params.(4)
    best_params.(5);
  let cand_opt = place_molecule cand_mol best_params in
  Ph4.to_file output_fn cand_opt

let () = main ()
