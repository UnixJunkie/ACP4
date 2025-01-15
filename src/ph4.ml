(* Copyright (C) 2022, Francois Berenger
   Tsuda laboratory, The University of Tokyo,
   5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.

   Read ph4 features and their 3D coordinates
   from a .sdf file with 3D conformers. *)

module A = BatArray
module IntMap = BatMap.Int
module IntSet = BatSet.Int
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module V3 = Vector3

open Printf

type feature = ARO
             | HYD
             | HBA
             | HBD
             | POS
             | NEG

let all_features = [|ARO; HYD; HBA; HBD; POS; NEG|]

let nb_features = A.length all_features
let nb_channels = 21 (* 6 + 6*5/2 *)

type molecule_ph4 = { name: string;
                      features: feature array;
                      coords: Vector3.t array }

let get_features m =
  m.features

let get_coords m =
  m.coords

(* all coords in [m] w/ the given ph4 feature type *)
let get_coords_with_feat m f =
  A.filteri (fun i _xyz ->
      f = m.features.(i)
    ) m.coords

(* a sparse float array *)
type channels = (float IntMap.t) array

let num_features m =
  A.length m.features

let feature_of_string = function
  | "ARO" -> ARO
  | "HYD" -> HYD
  | "HBA" -> HBA
  | "HBD" -> HBD
  | "POS" -> POS
  | "NEG" -> NEG
  | unk ->
    let () = Log.fatal "Sdf_3D.feature_of_string: %s" unk in
    exit 1

let string_of_feature = function
  | ARO -> "ARO"
  | HYD -> "HYD"
  | HBA -> "HBA"
  | HBD -> "HBD"
  | POS -> "POS"
  | NEG -> "NEG"

(* write molecule to .ph4 file *)
let write_out out mol =
  fprintf out "%d:%s\n" (A.length mol.features) mol.name;
  A.iter2 (fun feat xyz ->
      let x, y, z = Vector3.to_triplet xyz in
      fprintf out "%s %f %f %f\n" (string_of_feature feat) x y z
    ) mol.features mol.coords

let sprintf_feat i mol =
  let feat = mol.features.(i) in
  let x, y, z = Vector3.to_triplet mol.coords.(i) in
  sprintf "%s %f %f %f" (string_of_feature feat) x y z

let average_vects l =
  let n = float (L.length l) in
  let sum = L.fold_left (fun acc v -> V3.add acc v) (V3.make 0. 0. 0.) l in
  Vector3.div sum n

(* if we want to simplify a complex query: average features which are within
   2A of each others *)
let cluster_features mol =
  Log.info "feats before cluster: %d" (A.length mol.features);
  A.iteri (fun i _feat ->
      Log.info "%s" (sprintf_feat i mol)
    ) mol.features;
  let threshold = 2.0 in (* Angstrom *)
  let n = A.length mol.features in
  let idx2cluster = A.init n (fun i -> i) in
  (* clustering / grouping *)
  for i = 0 to n - 2 do
    let feat_i = mol.features.(i) in
    let pos_i = mol.coords.(i) in
    let clust_i = idx2cluster.(i) in
    for j = i + 1 to n - 1 do
      let feat_j = mol.features.(j) in
      if feat_i = feat_j then
        let pos_j = mol.coords.(j) in
        if V3.dist pos_i pos_j <= threshold then
          (* same cluster *)
          idx2cluster.(j) <- clust_i
    done
  done;
  (* they might be non consecutive integers *)
  let cluster_ids = IntSet.of_array idx2cluster in
  let m = IntSet.cardinal cluster_ids in
  (* Log.info "%d clusters" m; *)
  let features = A.make m HYD in
  let coords = A.make m (V3.make 0. 0. 0.) in
  let cid_i = ref 0 in
  IntSet.iter (fun cid ->
      (* Log.info "cid: %d" cid; *)
      features.(!cid_i) <- mol.features.(cid);
      (* average members of same cluster *)
      let res = ref [] in
      A.iteri (fun i cid' ->
          if cid = cid' then
            res := mol.coords.(i) :: !res
        ) idx2cluster;
      coords.(!cid_i) <- average_vects !res;
      incr cid_i
    ) cluster_ids;
  Log.info "feats after cluster: %d" (A.length features);
  let res = { name = mol.name; features; coords } in
  A.iteri (fun i _feat ->
      Log.info "%s" (sprintf_feat i res)
    ) features;
  res

(* canonicalization of the feature pair happens here *)
let channel_of_features = function
  | (ARO, ARO) -> 0
  | (HYD, HYD) -> 1
  | (HBA, HBA) -> 2
  | (HBD, HBD) -> 3
  | (POS, POS) -> 4
  | (NEG, NEG) -> 5
  | (ARO, HYD) | (HYD, ARO) ->  6
  | (ARO, HBA) | (HBA, ARO) ->  7
  | (ARO, HBD) | (HBD, ARO) ->  8
  | (ARO, POS) | (POS, ARO) ->  9
  | (ARO, NEG) | (NEG, ARO) -> 10
  | (HYD, HBA) | (HBA, HYD) -> 11
  | (HYD, HBD) | (HBD, HYD) -> 12
  | (HYD, POS) | (POS, HYD) -> 13
  | (HYD, NEG) | (NEG, HYD) -> 14
  | (HBA, HBD) | (HBD, HBA) -> 15
  | (HBA, POS) | (POS, HBA) -> 16
  | (HBA, NEG) | (NEG, HBA) -> 17
  | (HBD, POS) | (POS, HBD) -> 18
  | (HBD, NEG) | (NEG, HBD) -> 19
  | (POS, NEG) | (NEG, POS) -> 20

let channels = [|(ARO, ARO);
                 (HYD, HYD);
                 (HBA, HBA);
                 (HBD, HBD);
                 (POS, POS);
                 (NEG, NEG);
                 (ARO, HYD);
                 (ARO, HBA);
                 (ARO, HBD);
                 (ARO, POS);
                 (ARO, NEG);
                 (HYD, HBA);
                 (HYD, HBD);
                 (HYD, POS);
                 (HYD, NEG);
                 (HBA, HBD);
                 (HBA, POS);
                 (HBA, NEG);
                 (HBD, POS);
                 (HBD, NEG);
                 (POS, NEG)|]

let channel_strings = [|"ARO-ARO";
                        "HYD-HYD";
                        "HBA-HBA";
                        "HBD-HBD";
                        "POS-POS";
                        "NEG-NEG";
                        "ARO-HYD";
                        "ARO-HBA";
                        "ARO-HBD";
                        "ARO-POS";
                        "ARO-NEG";
                        "HYD-HBA";
                        "HYD-HBD";
                        "HYD-POS";
                        "HYD-NEG";
                        "HBA-HBD";
                        "HBA-POS";
                        "HBA-NEG";
                        "HBD-POS";
                        "HBD-NEG";
                        "POS-NEG"|]

let features_of_channel c =
  channels.(c)

let string_of_channel c =
  channel_strings.(c)

let parse_ph4_line to_parse =
  (* "^HYD 1.0 2.0 3.0$" *)
  try
    Scanf.sscanf to_parse "%s@ %f %f %f"
      (fun ph4 x y z ->
         (feature_of_string ph4, Vector3.make x y z))
  with exn ->
    let () = Log.fatal "Sdf_3D.parse_ph4_line: cannot parse: '%s'"
        to_parse in
    raise exn

let read_and_parse_ph4_line input =
  parse_ph4_line (input_line input)

let get_name_num_features line =
  Scanf.sscanf line "%d:%s" (fun num_feats name ->
      (name, num_feats)
    )

let read_name_num_ph4_features input =
  get_name_num_features (input_line input)

(* reader for the output of molenc_ph4.py *)
let read_one_ph4_encoded_molecule input =
  let name, num_features = read_name_num_ph4_features input in
  (if num_features <= 1 then
     Log.warn "mol %s: %d ph4-feats" name num_features
  );
  let features = Array.make num_features ARO in
  let coords =
    Array.init num_features
      (fun i ->
         let (feat, xyz) = read_and_parse_ph4_line input in
         features.(i) <- feat;
         xyz
      ) in
  { name; features; coords }

(* reader for the parallel encoder (does less work) *)
let preread_one_ph4_encoded_molecule input: string * (string array) =
  let name, num_features = read_name_num_ph4_features input in
  (if num_features <= 1 then
     Log.warn "mol %s: %d ph4-feats" name num_features
  );
  let lines = Array.init num_features (fun _i -> input_line input) in
  (name, lines)

(* parser for the parallel encoder (in worker's loop) *)
let parse_one_ph4_encoded_molecule (name, lines): molecule_ph4 =
  let num_features = A.length lines in
  let features = Array.make num_features ARO in
  let coords =
    Array.init num_features
      (fun i ->
         let (feat, xyz) = parse_ph4_line lines.(i) in
         features.(i) <- feat;
         xyz
      ) in
  { name; features; coords }

(* [cutoff]: longest interatomic distance allowed (A)
   [dx]: axis discretization step
   [mol]: molecule to encode *)
let encode verbose (cutoff: float) (dx: float) (mol: molecule_ph4) =
  let res = A.make nb_channels IntMap.empty in
  let n = A.length mol.features in
  for i = 0 to n - 2 do
    let feat_i = A.unsafe_get mol.features i in
    let xyz_i = A.unsafe_get mol.coords i in
    (* WARNING: this way, molecules w/ a single ph4 feature are all 0s *)
    for j = i + 1 to n - 1 do
      let feat_j = A.unsafe_get mol.features j in
      let xyz_j = A.unsafe_get mol.coords j in
      let feat_pair = (feat_i, feat_j) in
      let chan = channel_of_features feat_pair in
      let dist = V3.dist xyz_i xyz_j in
      if verbose && dist > cutoff then
        Log.warn "%g > %g" dist cutoff
      else (* dist <= cutoff *)
        let bin_before = int_of_float (dist /. dx) in
        let bin_after = bin_before + 1 in
        let before = dx *. (float bin_before) in
        let after = before +. dx in
        (* linear binning *)
        let w_l = (after -. dist) /. dx in
        let w_r = 1.0 -. w_l in
        (if verbose then
           Log.debug "chan: %02d:%s left: %g right: %g x_l: %d x_r: %d \
                      w_l: %g w_r: %g dist: %g"
             chan (string_of_channel chan) before after bin_before bin_after
             w_l w_r dist
        );
        (* radial.(bin_before).(chan) <- radial.(bin_before).(chan) +. w_l  ;
           radial.(bin_after).(chan)  <- radial.(bin_after).(chan)  +. w_r *)
        let res_chan = A.unsafe_get res chan in
        let prev_before = IntMap.find_default 0.0 bin_before res_chan in
        let prev_after  = IntMap.find_default 0.0 bin_after  res_chan in
        (* update map *)
        res.(chan) <-
          IntMap.add bin_before (prev_before +. w_l)
            (IntMap.add bin_after (prev_after +. w_r) res_chan)
    done
  done;
  (mol.name, res)

(* read then encode all molecules from file *)
let read_then_encode_ph4_molecules input_fn verbose cutoff dx cluster =
  LO.with_in_file input_fn (fun input ->
      let res, exn =
        L.unfold_exn (fun () ->
            let cand_ph4_mol =
              let cand = read_one_ph4_encoded_molecule input in
              if cluster then
                cluster_features cand
              else
                cand in
            encode verbose cutoff dx cand_ph4_mol) in
      if exn = End_of_file then ()
      else raise exn;
      res
    )
