# ---- Packages ----
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(rlang)

# ---- EXPECTED INPUT ----
# A data frame at the FOV level with (at minimum) these columns:
#   molecule_id : chr/int (10,000 unique values)
#   well_id     : chr/int (3 per molecule)
#   fov_id      : chr/int (4 per well)
#   measurement : numeric  (assay readout per FOV)
#   label       : int/logical/factor (ground truth at MOLECULE level: 1=active, 0=inactive)
#
# Notes:
# - label must be constant within a molecule (copied to all rows for convenience).
# - If your labels are only at molecule level, just join/merge them onto this FOV table.

# ---- TUNABLE AGGREGATION STRATEGY ----
# You can customize how you aggregate FOV -> Well -> Molecule.
# By default: mean at each level (robust alternatives provided).
agg_fov_to_well <- function(.x) mean(.x, na.rm = TRUE)
agg_well_to_molecule <- function(.x) mean(.x, na.rm = TRUE)

# Optional robust alternatives you might try:
# agg_fov_to_well <- function(.x) median(.x, na.rm = TRUE)
# agg_well_to_molecule <- function(.x) median(.x, na.rm = TRUE)
# Or “signal enhancing”: max/trimmed mean etc.

# ---- PR CURVE CORE ----
# Compute precision-recall pairs over a grid of thresholds.
# scores: numeric vector of prediction scores (higher = more likely positive)
# labels: logical or 0/1 vector (1/TRUE = positive class)
# n_thresholds: number of thresholds to evaluate (use length(unique(scores)) for exact)
compute_pr_curve <- function(scores, labels, n_thresholds = 200) {
  labels <- as.integer(labels) # 1 for positive, 0 for negative
  stopifnot(length(scores) == length(labels))
  
  # Threshold grid (descending so recall moves from 0 -> 1)
  thr <- quantile(scores, probs = seq(0, 1, length.out = n_thresholds), na.rm = TRUE) %>%
    unique() %>%
    sort(decreasing = TRUE)
  
  # Vectorized computation over thresholds
  tp <- fp <- fn <- numeric(length(thr))
  
  for (i in seq_along(thr)) {
    pred_pos <- scores >= thr[i]
    tp[i] <- sum(pred_pos & labels == 1, na.rm = TRUE)
    fp[i] <- sum(pred_pos & labels == 0, na.rm = TRUE)
    fn[i] <- sum(!pred_pos & labels == 1, na.rm = TRUE)
  }
  
  precision <- ifelse(tp + fp > 0, tp / (tp + fp), 1) # define P=1 when no predicted positives
  recall    <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
  
  tibble(
    threshold = thr,
    precision = precision,
    recall    = recall,
    tp = tp, fp = fp, fn = fn
  ) %>%
    arrange(recall, precision)
}

# ---- (Optional) PR-AUC (Average Precision-ish via trapezoid) ----
pr_auc_trapezoid <- function(pr_tbl) {
  pr_tbl <- pr_tbl %>% arrange(recall)
  # Trapezoidal area under P(R)
  r <- pr_tbl$recall
  p <- pr_tbl$precision
  sum(diff(r) * (head(p, -1) + tail(p, -1)) / 2, na.rm = TRUE)
}

# ---- PIPELINE: FROM FOV -> MOLECULE SCORE -> PR CURVE ----
# 'fov_df' is your input data frame (FOV-level)
build_pr_from_fov <- function(fov_df,
                              agg1 = agg_fov_to_well,
                              agg2 = agg_well_to_molecule,
                              n_thresholds = 200) {
  
  # Sanity checks
  required_cols <- c("molecule_id", "well_id", "fov_id", "measurement", "label")
  missing_cols <- setdiff(required_cols, names(fov_df))
  if (length(missing_cols) > 0) {
    abort(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Ensure label is at molecule level (warn if not constant)
  chk <- fov_df %>%
    group_by(molecule_id) %>%
    summarize(label_range = n_distinct(label), .groups = "drop")
  if (any(chk$label_range != 1)) {
    warning("Some molecules have multiple label values within them. Using first encountered per molecule.")
  }
  
  # 1) Aggregate FOV -> WELL
  well_level <- fov_df %>%
    group_by(molecule_id, well_id) %>%
    summarize(well_score = agg1(measurement),
              label = first(label), .groups = "drop")
  
  # 2) Aggregate WELL -> MOLECULE
  mol_level <- well_level %>%
    group_by(molecule_id) %>%
    summarize(score = agg2(well_score),
              label = first(label), .groups = "drop")
  
  # Convert label to 0/1
  lab <- mol_level$label
  if (is.logical(lab))       lab <- as.integer(lab)
  if (is.factor(lab))        lab <- as.integer(lab) - 1L
  if (!all(lab %in% c(0, 1))) abort("Labels must be logical or encoded as 0/1 at the molecule level.")
  
  pr_tbl <- compute_pr_curve(scores = mol_level$score, labels = lab, n_thresholds = n_thresholds)
  attr(pr_tbl, "pr_auc_trap") <- pr_auc_trapezoid(pr_tbl)
  
  list(
    pr = pr_tbl,
    mol_scores = mol_level
  )
}

# ---- PLOTTING ----
plot_pr_curve <- function(pr_tbl, title = "Precision–Recall Curve (Molecule-Level)") {
  auc_est <- attr(pr_tbl, "pr_auc_trap")
  ggplot(pr_tbl, aes(x = recall, y = precision)) +
    geom_path(linewidth = 1) +
    geom_point(size = 0.8, alpha = 0.4) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(
      title = title,
      subtitle = if (!is.null(auc_est)) sprintf("Trapezoidal PR-AUC ≈ %.3f", auc_est) else NULL,
      x = "Recall (TP / (TP + FN))",
      y = "Precision (TP / (TP + FP))"
    ) +
    theme_minimal(base_size = 12)
}

# ---- EXAMPLE USAGE (Replace this with your real data) ----
# Suppose you already have a data frame `fov_df` with the required columns.
# Below is a minimal schematic generator you can REMOVE when using real data.

# set.seed(1)
# n_mol  <- 10000
# n_well <- 3
# n_fov  <- 4
# mol_df <- tibble(
#   molecule_id = rep(sprintf("mol_%05d", 1:n_mol), each = n_well * n_fov),
#   well_id     = rep(rep(paste0("w", 1:n_well), each = n_fov), times = n_mol),
#   fov_id      = rep(paste0("f", 1:n_fov), times = n_mol * n_well)
# ) %>%
#   group_by(molecule_id) %>%
#   mutate(label = rbinom(1, 1, 0.2)) %>%  # 20% actives, constant per molecule
#   ungroup() %>%
#   mutate(
#     # Simulate FOV-level measurements: actives shifted higher
#     measurement = rnorm(n(), mean = 0 + 2 * label, sd = 1)
#   )
#
# out <- build_pr_from_fov(mol_df)
# plot_pr_curve(out$pr)

# ---- (Optional) Compare aggregation strategies in one plot ----
# If you want to see how different aggregation choices change the PR curve:
# strategies <- tribble(
#   ~name,       ~agg1,                         ~agg2,
#   "mean-mean", agg_fov_to_well,               agg_well_to_molecule,
#   "median-med", function(x) median(x, TRUE),  function(x) median(x, TRUE),
#   "max-mean",  function(x) max(x, na.rm=TRUE),agg_well_to_molecule
# )
#
# pr_list <- pmap(strategies, function(name, agg1, agg2) {
#   res <- build_pr_from_fov(fov_df = mol_df, agg1 = agg1, agg2 = agg2)
#   res$pr %>% mutate(strategy = name, auc = attr(res$pr, "pr_auc_trap"))
# }) %>% bind_rows()
#
# ggplot(pr_list, aes(recall, precision, color = strategy)) +
#   geom_path() +
#   coord_equal(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
#   labs(title = "PR Curves by Aggregation Strategy",
#        x = "Recall", y = "Precision", color = "Strategy") +
#   theme_minimal(base_size = 12)
