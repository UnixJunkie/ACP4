# ACP4: AutoCorrelation of Pharmacophore Features

```
acp4 program usage:
  ./acp4
  [-i <filename.ph4>]: input file to encode
  [-q <filename.ph4>]: query molecule
  [--queries <filename.ph4>]: query molecules
  [-db <filename.ph4>]: database to screen
  -o <filename.{csv|scores}>: output file (encoded db or query scores)
  [-np <int>]: nprocs (default=1)
  [-s <int>]: BST chunk size (default=100000)
  [-c <float>]: cutoff distance (default=5.00)
  [-f]: force overwriting binary cache files, if any
  use -f if you are changing -dx or -c compared to previous queries
  [-dx <float>]: radial discretization step (default=0.5)
  [--no-plot]: turn OFF gnuplot
  [--no-tap]: neither read nor write from/to binary cache files
  [--confs]: if the DB has several _consecutive_ conformers
  [--quick]: quick; exit right after scoring; no perf. metrics
  [--index]: create a BST index for the DB being screened
  [--BSTs <filename.txt>]: file containing a list of serialized BST filenames
  [-td <float>]: maximum Tanimoto _distance_ to (single) query
  [--BS]: load optimal defaults for binding-sites (ignores -c and -dx)
  [-v]: verbose/debug mode
```

```
acp4_scissors program usage:
  ./scissors
  -l <ligand.sdf>: binding-site ligand input file
  -p <protein.ph4>: receptor protein input file
  [-d <float>]: distance cutoff (default=5.00)
  -o <output.ph4>: ligand-defined binding site output file
  [-v]: verbose/debug mode
```
