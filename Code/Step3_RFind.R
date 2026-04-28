# ============================================================
# Step3_RFind.R : RFscore-G 計算 + 診断プロット (ROSMAP)
#
# 入力 (Step1_SNPs/, Step2_Analysis/):
#   ROSMAP_metadata.csv
#   SNPs_weights_<PGS_ID>.csv  (Step2 出力)
#   SNPs_flipped_<PGS_ID>.csv  (Step2 出力, allele 整合済 dosage)
# 出力:
#   RFinD_G_scores_<PGS_ID>.csv         (individualID, RFscore, RFscore_z + clinical)
#   Step3_RFscore_plots_<DATE>_<var>.png × 7
#
# アルゴリズム:
#   g1 : PGS β を方向別 |β| 順で並べた順位
#   g2 : 個人 dosage 偏差を方向別 |fc| 順で並べた順位
#   RFscore = RF(g1↑,g2↑) + RF(g1↓,g2↓) − RF(g1↑,g2↓) − RF(g1↓,g2↑)
#       RF = 上位 N×M グリッドでの最大 -log10(hypergeom p)
# ============================================================

rm(list = ls())
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(dplyr) })
# parallel は :: prefix で利用 (library 不要)

pgs_id <- rstudioapi::showPrompt("PGS ID", "PGS ID:", "PGS004228")
if (is.null(pgs_id) || trimws(pgs_id) == "") stop("PGS ID 未入力")
pgs_id <- trimws(pgs_id)
cat("PGS ID:", pgs_id, "\n")

# ---- パス & パラメータ ----
# 実行前に setwd("/path/to/RFind-G") (project root) してください
base_dir  <- normalizePath(getwd(), mustWork = TRUE)
step1_dir <- file.path(base_dir, "Step1_SNPs")
step2_dir <- file.path(base_dir, "Step2_Analysis")
workdir   <- file.path(base_dir, "Step3_RFind")
if (!dir.exists(workdir)) dir.create(workdir, recursive = TRUE)
setwd(workdir)
GRID_STEP <- 50L
TOP_K     <- 10000L
n_cores   <- max(1L, parallel::detectCores() - 1L)

# ============================================================
# (1) 重み (g1) 構築
# ============================================================
w <- fread(file.path(step2_dir, paste0("SNPs_weights_", pgs_id, ".csv")))
w <- w[present_in_Dosage == TRUE & mismatch == FALSE, .(SNP = id, beta = beta)]
P <- nrow(w)
cat("Effective SNPs P =", P, "\n")

grid_seq <- function(n, step) {
  if (n < step) return(integer(0))
  c(seq.int(step, n, step), if (n %% step) n else integer(0))
}
build_g1 <- function(idx, abs_beta, P, step) {
  if (length(idx) > 1L) idx <- idx[order(-abs_beta[idx])]
  rank_map <- integer(P); rank_map[idx] <- seq_along(idx)
  list(n_idx = length(idx), N_seq = grid_seq(length(idx), step),
       rank_map = rank_map)
}
g1_up   <- build_g1(which(w$beta > 0), abs(w$beta), P, GRID_STEP)
g1_down <- build_g1(which(w$beta < 0), abs(w$beta), P, GRID_STEP)

# ============================================================
# (2) SNPs_flipped を一括ロード → P × N_samp 行列
# ============================================================
flip <- fread(file.path(step2_dir, paste0("SNPs_flipped_", pgs_id, ".csv")),
              na.strings = c("NA","NaN",""))
sample_names <- setdiff(names(flip), c("chr","id","allele1","allele2"))
N <- length(sample_names)
cat("Samples:", N, "\n")

pos  <- match(flip$id, w$SNP)
flip <- flip[!is.na(pos)]; pos <- pos[!is.na(pos)]
M    <- matrix(NA_real_, nrow = P, ncol = N,
               dimnames = list(NULL, sample_names))
M[pos, ] <- as.matrix(flip[, ..sample_names])
mode(M)  <- "numeric"
rm(flip); invisible(gc(verbose = FALSE))

# metadata は Control 同定 (Step3 中盤) と最終 merge (Step3 末尾) で利用するので一度だけ読む
# ROSMAP では dcfdx_lv == 1 (NCI) を Control とする
meta <- fread(file.path(step1_dir, "ROSMAP_metadata.csv"), na.strings = c("NA",""))
meta[, individualID := trimws(as.character(individualID))]

ctrl_ids  <- meta$individualID[!is.na(meta$dcfdx_lv) & meta$dcfdx_lv == 1L]
ctrl_cols <- which(sample_names %in% ctrl_ids)
if (length(ctrl_cols) < 10L) {
  warning("Control < 10. Falling back to all-sample mean.")
  ref_mean <- rowMeans(M, na.rm = TRUE)
} else {
  cat("Reference = mean over", length(ctrl_cols), "Control samples (dcfdx_lv == 1)\n")
  ref_mean <- rowMeans(M[, ctrl_cols, drop = FALSE], na.rm = TRUE)
}
FC <- M - ref_mean
rm(M); invisible(gc(verbose = FALSE))

# ============================================================
# (3) Running Fisher + サンプルごとスコアリング（並列）
# ============================================================
running_fisher <- function(g1, g2_idx, P, step) {
  if (g1$n_idx == 0L || length(g2_idx) == 0L) return(0)
  M_seq <- grid_seq(length(g2_idx), step)
  if (!length(M_seq)) return(0)
  r <- g1$rank_map[g2_idx]; hit <- (r > 0L)
  best <- 0
  for (Nn in g1$N_seq) {
    overlaps <- cumsum(hit & (r <= Nn))[M_seq]
    pv <- phyper(overlaps - 1L, Nn, P - Nn, M_seq, lower.tail = FALSE)
    lb <- max(-log10(pv))
    if (lb > best) { best <- lb; if (is.infinite(best)) return(best) }
  }
  best
}
score_one <- function(j) {
  fc <- FC[, j]; ok <- !is.na(fc)
  i_up <- which(ok & fc > 0); if (length(i_up) > 1L) i_up <- i_up[order(-fc[i_up])]
  i_dn <- which(ok & fc < 0); if (length(i_dn) > 1L) i_dn <- i_dn[order( fc[i_dn])]
  if (length(i_up) > TOP_K) i_up <- i_up[seq_len(TOP_K)]
  if (length(i_dn) > TOP_K) i_dn <- i_dn[seq_len(TOP_K)]
  running_fisher(g1_up,   i_up, P, GRID_STEP) +
  running_fisher(g1_down, i_dn, P, GRID_STEP) -
  running_fisher(g1_up,   i_dn, P, GRID_STEP) -
  running_fisher(g1_down, i_up, P, GRID_STEP)
}

cat("Scoring", N, "samples on", n_cores, "cores ...\n")
t0 <- Sys.time()
rf <- (
  if (n_cores > 1L)
    unlist(parallel::mclapply(seq_len(N), score_one,
                              mc.cores = n_cores, mc.preschedule = TRUE))
  else
    vapply(seq_len(N), score_one, numeric(1))
)
cat("Scoring time:",
    round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "sec\n")

# ============================================================
# (4) 出力 RFinD_G_scores
# ============================================================
sd_rf <- sd(rf)
rf_z  <- (
  if (is.na(sd_rf) || sd_rf == 0) rep(0, length(rf))
  else as.numeric(scale(rf))
)

out <- merge(data.table(individualID = sample_names, RFscore = rf, RFscore_z = rf_z),
             meta, by = "individualID", sort = FALSE)
fwrite(out, paste0("RFinD_G_scores_", pgs_id, ".csv"), na = "NA")
cat("Saved RFinD_G_scores. RFscore raw range = [",
    round(min(rf), 2), ",", round(max(rf), 2), "]\n\n")

# ============================================================
# (5) 診断プロット 7 種
# ============================================================
d <- as.data.frame(out, stringsAsFactors = FALSE) %>%
  mutate(across(c(msex, braaksc, ceradsc, cogdx, dcfdx_lv),
                ~ trimws(as.character(.x)))) %>%
  mutate(
    msex     = factor(msex,     levels = c("0","1"),       labels = c("Female","Male")),
    braaksc  = factor(braaksc,  levels = as.character(0:6)),
    ceradsc  = factor(ceradsc,  levels = as.character(1:4),
                      labels = c("1.Definite","2.Probable","3.Possible","4.No AD")),
    cogdx    = factor(cogdx,    levels = as.character(1:6),
                      labels = c("NCI","MCI","MCI+","AD","AD+","Dementia")),
    dcfdx_lv = factor(dcfdx_lv, levels = as.character(1:6),
                      labels = c("NCI","MCI","MCI+","AD","AD+","Dementia")),
    cts_mmse30_lv = as.numeric(cts_mmse30_lv),
    mmse_cat = cut(cts_mmse30_lv, breaks = c(-Inf, 9, 18, 24, Inf),
                   labels = c("0-9","10-18","19-24",">=25"), ordered_result = TRUE)
  )

theme_p <- theme_bw(base_size = 14) + theme(
  legend.position  = "none",
  axis.title       = element_text(size = 24, color = "black"),
  axis.text        = element_text(size = 18, color = "black"),
  panel.grid.major = element_line(color = "grey85"),
  panel.grid.minor = element_blank(),
  aspect.ratio     = 1
)
plot_box <- function(xv, xlab) {
  d2 <- d[!is.na(d[[xv]]) & !is.na(d$RFscore_z), , drop = FALSE]
  ggplot(d2, aes(.data[[xv]], RFscore_z, fill = .data[[xv]])) +
    geom_boxplot(width = 0.6, alpha = 0.70, outlier.size = 1.5) +
    geom_point(position = position_jitter(width = 0.12, height = 0),
               size = 1.6, alpha = 0.35, color = "black") +
    labs(x = xlab, y = "RFscore z-score") + theme_p
}
plot_lm <- function(xv, xlab) {
  d2 <- d[!is.na(d[[xv]]) & !is.na(d$RFscore_z), , drop = FALSE]
  ggplot(d2, aes(.data[[xv]], RFscore_z)) +
    geom_point(alpha = 0.35, size = 1.6) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
    labs(x = xlab, y = "RFscore z-score") + theme_p
}
plots <- list(
  sex       = plot_box("msex",          "Sex"),
  braak     = plot_box("braaksc",       "Braak stage"),
  cerad     = plot_box("ceradsc",       "CERAD score"),
  cogdx     = plot_box("cogdx",         "Final consensus cognitive diagnosis"),
  dcfdx_lv  = plot_box("dcfdx_lv",      "Last valid clinical cognitive diagnosis"),
  mmse_cont = plot_lm ("cts_mmse30_lv", "MMSE score"),
  mmse_cat  = plot_box("mmse_cat",      "MMSE category")
)

prefix  <- paste0("Step3_RFscore_plots_", format(Sys.Date(), "%Y%m%d"))
for (nm in names(plots)) {
  ggsave(file.path(workdir, paste0(prefix, "_", nm, ".png")),
         plots[[nm]], width = 6, height = 6, dpi = 300)
}
cat("Saved 7 plots in", workdir, "with prefix:", prefix, "\n")
