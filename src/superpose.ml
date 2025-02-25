
module A = BatArray
module CLI = Minicli.CLI
module Log = Dolog.Log
module LO = Line_oriented
module V3 = Mmo.V3
module Rot = Mmo.Rot
module SO3 = Mmo.SO3

open Printf

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
          try
            let _nn, dist = BST.nearest_neighbor xyz bst in
            error := !error +. dist
          with Not_found ->
            () (* this feature is unmatched *)
        ) cand_coords
    ) Ph4.all_features;
  !error

let geometric_center vects =
  V3.div
    (A.fold_left V3.add V3.origin vects)
    (float (A.length vects))

let center_of m =
  geometric_center Ph4.(m.coords)

let translate_by vects delta =
  let dx, dy, dz = V3.to_triplet delta in
  A.map (fun V3.{x; y; z} ->
      V3.make (x +. dx) (y +. dy) (z +. dz)
    ) vects

let translate_to vects p =
  translate_by vects (V3.diff p (geometric_center vects))

let center_molecule m =
  Ph4.{ m with coords = translate_to m.coords V3.origin }

let translate_molecule m p =
  Ph4.{ m with coords = translate_to m.coords p }

(* CONSTRAINT: vects must already be centered at origin *)
let centered_rotate vects r =
  A.map (Rot.rotate r) vects

let rotate_centered_molecule m r =
  Ph4.{ m with coords = centered_rotate m.coords r }

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
let nlopt_eval_solution bsts centered_m params _gradient =
  let posed = place_molecule centered_m params in
  let error = score_molecule_pose bsts posed in
  Log.debug "err: %f" error;
  error

let pi = Mmo.Math.pi
let two_pi = Mmo.Math.two_pi

let nlopt_optimize init_params bsts centered_m max_evals =
  let ndims = 6 in (* DOFs: 3 rotations + 3 translations *)
  (* local optimizer that will be passed to the global one *)
  let local = Nlopt.(create sbplx ndims) in (* local optimizer: gradient-free *)
  Nlopt.set_min_objective local (nlopt_eval_solution bsts centered_m);
  (* I don't set parameter bounds on the local optimizer, I guess
   * the global optimizer handles this *)
  (* hard stop conditions *)
  Nlopt.set_stopval local 0.0; (* absolute minimum total error *)
  (* (\* smart stop conditions *\) *)
  (* Nlopt.set_ftol_abs local 0.0001; (\* might need to be tweaked *\) *)
  let global = Nlopt.(create auglag ndims) in (* global optimizer *)
  Nlopt.set_local_optimizer global local;
  Nlopt.set_min_objective global (nlopt_eval_solution bsts centered_m);
  (* bounds for parameters *)
  (* FBR: found tx_min, ty_min, tz_min using bounding boxes? *)
  (*                              rz    ry    rz    tx       ty       tz *)
  Nlopt.set_lower_bounds global [|-.pi; -.pi; -.pi; -1000.0; -1000.0; -1000.0|];
  Nlopt.set_upper_bounds global [|  pi;   pi;   pi;  1000.0;  1000.0;  1000.0|];
  (* hard stop conditions *)
  Nlopt.set_stopval global 0.0;
  (* max number of iterations *)
  Nlopt.set_maxeval global max_evals;
  let stop_cond, params, error = Nlopt.optimize global init_params in
  Log.info "NLopt optimize global: %s" (Nlopt.string_of_result stop_cond);
  Log.info "error: %f" error;
  (params, error)

let random_start m =
  Random.self_init ();
  place_molecule m [|-.pi +. Random.float two_pi;
                     -.pi +. Random.float two_pi;
                     -.pi +. Random.float two_pi;
                     -1000.0 +. Random.float 2000.;
                     -1000.0 +. Random.float 2000.;
                     -1000.0 +. Random.float 2000.|]

let bsts_for_molecule m =
  [|BST.(create 1 Two_bands) (Ph4.get_coords_with_feat m ARO);
    BST.(create 1 Two_bands) (Ph4.get_coords_with_feat m HYD);
    BST.(create 1 Two_bands) (Ph4.get_coords_with_feat m HBA);
    BST.(create 1 Two_bands) (Ph4.get_coords_with_feat m HBD);
    BST.(create 1 Two_bands) (Ph4.get_coords_with_feat m POS);
    BST.(create 1 Two_bands) (Ph4.get_coords_with_feat m NEG)|]

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
              [-oph4 <filename.ph4>]: optimally superposed cand. structure\n  \
              [-o <filename.txt>]: 6 optimal parameters\n  \
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
  let txt_out_fn = CLI.get_string ["-o"] args in
  let ph4_out_fn = CLI.get_string ["-oph4"] args in
  CLI.finalize (); (* -------------------------------------------------- *)
  let ref_mol = read_one_ph4_molecule ref_fn in
  let cand_mol =
    center_molecule
      (let mol = read_one_ph4_molecule cand_fn in
       if rand_start then
         random_start mol
       else
         mol) in
  (* ROTATIONAL SEARCH ONLY ------------------------------------------------ *)
  Log.info "discrete rotational search...";
  (* center both molecules, then brute force search over a discretized
     space of all possible rotations *)
  let centered_ref_mol = center_molecule ref_mol in
  let centered_ref_mol_bsts = bsts_for_molecule centered_ref_mol in
  (* at the end, translate the best rotated candidate there *)
  let ref_mol_init_center = center_of ref_mol in
  (* brute force search of the best rotation *)
  let rotations = SO3.rotations 100_000 in
  let best_score = ref infinity in
  let best_rot = ref (Rot.id ()) in
  A.iter (fun rot ->
      let rotated = rotate_centered_molecule cand_mol rot in
      let score = score_molecule_pose centered_ref_mol_bsts rotated in
      if score < !best_score then
        (Log.debug "score: %f" score;
         best_score := score;
         best_rot := rot)
    ) rotations;
  Log.info "error: %f" !best_score;
  (* rotate centered cand. *)
  let rotated = rotate_centered_molecule cand_mol !best_rot in
  (* translate it back *)
  let translated = translate_molecule rotated ref_mol_init_center in
  Log.info "writing: %s" ph4_out_fn;
  Ph4.to_file ph4_out_fn translated;
  (* try global optim. of this initial solution
     if it improves the current solution, we will report it *)
  (* GLOBAL SEARCH --------------------------------------------------------------- *)
  Log.info "global optimization...";
  (* Each ph4 type in A is BST indexed first *)
  let ref_mol_bsts = bsts_for_molecule ref_mol in
  (* - find the best rotation and translation of B that minimizes error(A, B). *)
  (* At the end, output the (rot, trans) w/ the smallest error. *)
  let rx, ry, rz = Rot.decompose !best_rot in
  let cx, cy, cz = V3.to_triplet ref_mol_init_center in
  (* start solution: best rotation found so far; translated at ref_mol's center *)
  let init_params = [|rx; ry; rz; cx; cy; cz|] in
  Log.info "writing: %s" txt_out_fn;
  LO.with_out_file txt_out_fn (fun out ->
      fprintf out "%f\n%f\n%f\n%f\n%f\n%f\n"
        init_params.(0)
        init_params.(1)
        init_params.(2)
        init_params.(3)
        init_params.(4)
        init_params.(5)
    );
  let best_params, error = nlopt_optimize init_params ref_mol_bsts cand_mol 100_000 in
  if error < !best_score then
    (Log.info "score improved: %f" error;
     LO.with_out_file txt_out_fn (fun out ->
         fprintf out "%f\n%f\n%f\n%f\n%f\n%f\n"
           best_params.(0)
           best_params.(1)
           best_params.(2)
           best_params.(3)
           best_params.(4)
           best_params.(5)
       );
     let cand_opt = place_molecule cand_mol best_params in
     Ph4.to_file ph4_out_fn cand_opt
    )

let () = main ()
