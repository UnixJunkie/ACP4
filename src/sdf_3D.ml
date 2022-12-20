(* Copyright (C) 2022, Francois Berenger
   Tsuda laboratory, The University of Tokyo,
   5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.

   Read molecule name, element symbols and 3D coordinates
   from a .sdf file holding 3D conformers. *)

module A = BatArray
module IntMap = BatMap.Int
module L = BatList
module Log = Dolog.Log
module S = BatString
module V3 = Vector3

type atoms_3D = { name: string;
                  elements: int array;
                  coords: Vector3.t array;
                  (* just which atom is connected to which other;
                     no bond order info *)
                  bonds: int list array }

let heavy_atom_coords mol =
  let res = ref [] in
  A.iter2 (fun anum xyz ->
      if anum > 1 then (* heavy atom? *)
        res := xyz :: !res
    ) mol.elements mol.coords;
  A.of_list !res

let anum_of_symbol = function
  | "C" -> 6
  | "H" -> 1
  | "N" -> 7
  | "O" -> 8
  | "P" -> 15
  | "S" -> 16
  | "F" -> 9
  | "Cl" -> 17
  | "Br" -> 35
  | "I" -> 53
  | unk ->
    let () = Log.warn "unsuported elt: %s" unk in
    -1

let read_name input =
  input_line input

let skip_header_lines input =
  let (_: string) = input_line input in
  let (_: string) = input_line input in
  ()

let read_atom_bonds_header input =
  let to_parse = input_line input in
  (* first two integers are 3-char fixed width each *)
  let num_atoms = int_of_string (S.strip (S.sub to_parse 0 3)) in
  let num_bonds = int_of_string (S.strip (S.sub to_parse 3 3)) in
  assert(S.ends_with to_parse "V2000");
  (num_atoms, num_bonds)

let read_bond_line input =
  let to_parse = input_line input in
  (* first two integers are 3-char fixed width each *)
  let src = int_of_string (S.strip (S.sub to_parse 0 3)) in
  let dst = int_of_string (S.strip (S.sub to_parse 3 3)) in
  (* indexes start at 1 in the SDF but atom indexes in the atoms array
   * start at 0 *)
  (src - 1, dst - 1)

let parse_atom_line input =
  let to_parse = input_line input in
  (* "^    5.0751   -3.8284   -4.0739 Br  0  0  0  ...$" *)
  try
    Scanf.sscanf to_parse
      " %f %f %f %s@ " (* ignore all the rest of the line *)
      (fun x y z elt_symbol ->
         (anum_of_symbol elt_symbol, Vector3.make x y z))
  with exn ->
    let () = Log.fatal "Sdf_3D.parse_atom_line: cannot parse: '%s'"
        to_parse in
    raise exn

exception Four_dollars

let read_one_molecule input =
  let name = read_name input in
  (skip_header_lines input);
  let num_atoms, num_bonds = read_atom_bonds_header input in
  let elements = Array.make num_atoms 0 in
  let coords =
    Array.init num_atoms
      (fun i ->
         let (anum, xyz) = parse_atom_line input in
         elements.(i) <- anum;
         xyz
      ) in
  (* read all bonds *)
  let bonds = A.create num_atoms [] in
  for _i = 1 to num_bonds do
    let src, dst = read_bond_line input in
    bonds.(src) <- dst :: bonds.(src)
  done;
  (* put them back in the same order than what was read *)
  for i = 0 to num_atoms - 1 do
    bonds.(i) <- L.rev bonds.(i)
  done;
  (try
     (* look for end of this molecule's record *)
     while true do
       if input_line input = "$$$$" then
         raise Four_dollars
     done;
     assert(false)
   with Four_dollars -> ()
  );
  { name; elements; coords; bonds }
