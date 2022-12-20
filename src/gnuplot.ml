(* Copyright (C) 2022, Francois Berenger
   Tsuda laboratory, The University of Tokyo,
   5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

open Printf

module Fn = Filename
module LO = Line_oriented
module Log = Dolog.Log
module Stats = Cpm.RegrStats

let protect_underscores title =
  BatString.nreplace ~str:title ~sub:"_" ~by:"\\\\_"

let run_command ?(debug = false) (cmd: string): unit =
  if debug then Log.info "run_command: %s" cmd;
  match Unix.system cmd with
  | Unix.WSIGNALED _ -> (Log.fatal "run_command: signaled: %s" cmd; exit 1)
  | Unix.WSTOPPED _ -> (Log.fatal "run_command: stopped: %s" cmd; exit 1)
  | Unix.WEXITED i when i <> 0 ->
    (Log.fatal "run_command: exit %d: %s" i cmd; exit 1)
  | Unix.WEXITED _ (* i = 0 then *) -> ()

(* comes from RanKers Gnuplot module *)
let roc_curve title roc_curve_fn nb_actives nb_decoys =
  (* Utls.run_command
   *   (sprintf "cat %s | time croc-curve 2>/dev/null > %s"
   *      score_labels_fn roc_curve_fn); *)
  let gnuplot_script_fn = Fn.temp_file ~temp_dir:"/tmp" "linwrap_" ".gpl" in
  LO.with_out_file gnuplot_script_fn (fun out ->
      fprintf out
        "set title \"|A|:|D|=%d:%d %s\"\n\
         set xtics out nomirror\n\
         set ytics out nomirror\n\
         set size square\n\
         set xrange [0:1]\n\
         set yrange [0:1]\n\
         set xlabel 'FPR'\n\
         set ylabel 'TPR'\n\
         set key outside right\n\
         f(x) = x\n\
         plot '%s' u 1:2 w lines t 'ROC' , \
         f(x) lc rgb 'black' not\n"
        nb_actives nb_decoys (protect_underscores title)
        roc_curve_fn
    );
  let gnuplot_log = Fn.temp_file ~temp_dir:"/tmp" "gnuplot_" ".log" in
  run_command (sprintf "(gnuplot -persist %s 2>&1) > %s"
                 gnuplot_script_fn gnuplot_log)
