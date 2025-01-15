
let main () =
  (*
  given two .ph4 files, assume:
  - LHS = A: the reference structure
  - RHS = B: the candidate structure
  - find the best rotation and translation of B that minimizes error(A, B).
  Each ph4 type in A is BST indexed first.
  To evaluate one rotation and one translation, the error function is:
  - for each ph4 point in B, add the distance to its nearest neighbor in A
    which has the same type into the error sum
  At the end, output the (rot, trans) w/ the smallest error.
  *)
  failwith "not implemented yet"
