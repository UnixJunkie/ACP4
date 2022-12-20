(* Copyright (C) 2022, Francois Berenger
   Tsuda laboratory, The University of Tokyo,
   5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

open Printf

module A = BatArray
module LO = Line_oriented
module Log = Dolog.Log

module SL = struct
  type t = bool * float (* (label, score) *)
  let create label score =
    (label, score)
  let get_score (_, s) = s
  let get_label (l, _) = l
  (* to do a decreasing sort of a score labels list *)
  let high_score_first_cmp (_, s1) (_, s2) =
    BatFloat.compare s2 s1
  let to_string (label, score) =
    sprintf "%d %f" label score
end

module ROC = Cpm.MakeROC.Make(SL)

let performance_metrics maybe_out_fn db_fn score_labels d_cut dx =
  ROC.rank_order_by_score_a score_labels;
  let auc = ROC.fast_auc_a score_labels in
  let ef_05 = ROC.fast_enrichment_factor 0.05 score_labels in
  let pm_05 = ROC.fast_power_metric_a 0.05 score_labels in
  let log_title =
    sprintf "AUC=%.3f EF5=%.2f PM5=%.2f DB=%s d_c=%g dx=%g"
      auc ef_05 pm_05 db_fn d_cut dx in
  let plot_title =
    sprintf "AUC=%.3f EF5=%.2f PM5=%.2f\\nDB=%s" auc ef_05 pm_05 db_fn in
  Log.info "%s" log_title;
  BatOption.may (fun out_fn ->
      Log.info "writing ROC to %s" out_fn;
      LO.with_out_file out_fn (fun out ->
          let roc = ROC.fast_roc_curve_a score_labels in
          A.iter (fun (fpr, tpr) ->
              fprintf out "%g %g\n" fpr tpr
            ) roc
        );
      let nb_actives = A.count_matching fst score_labels in
      let nb_decoys = (A.length score_labels) - nb_actives in
      Gnuplot.roc_curve plot_title out_fn nb_actives nb_decoys
    ) maybe_out_fn

let int_of_bool = function
  | true -> 1
  | false -> 0
