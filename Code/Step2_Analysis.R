# ============================================================
# Step2_Analysis.R : PRS 計算 + 診断プロット (ROSMAP)
#
# 入力 (Step1_SNPs/):
#   ROSMAP_metadata.csv         (specimenID, individualID + clinical)
#   <PGS_ID>.txt                (PGS Catalog 形式: rsID, effect_allele, effect_weight)
#   effectSNPs_<PGS_ID>.csv     (chr,id,allele1,allele2 + 個人 dosage 列)
#
# 出力:
#   SNPs_flipped_<PGS_ID>.csv   (allele 整合済 dosage)
#   SNPs_weights_<PGS_ID>.csv   (PGS 監査表: present/flipped/mismatch/used)
#   PRS_<PGS_ID>.csv            (individualID, PRS_raw, PRS_z + clinical)
#   Step2_PRSz_plots_<DATE>_<var>.png × 7
# ============================================================

rm(list = ls())
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(dplyr) })

# ---- PGS ID 入力 ----
pgs_id <- rstudioapi::showPrompt("PGS ID", "使用する PGS ID:", "PGS004228")
if (is.null(pgs_id) || trimws(pgs_id) == "") stop("PGS ID 未入力")
pgs_id <- trimws(pgs_id)
cat("PGS ID:", pgs_id, "\n")

# ---- パス ----
base_dir  <- "~/Miyakawa Lab Dropbox/Murano Tomoyuki/RFind-G_260428"
raw_dir   <- file.path(base_dir, "RawData")        # PGS txt 等の原データ
step1_dir <- file.path(base_dir, "Step1_SNPs")     # Step1 出力 (metadata, effectSNPs)
workdir   <- file.path(base_dir, "Step2_Analysis") # Step2 出力
if (!dir.exists(workdir)) dir.create(workdir, recursive = TRUE)
setwd(workdir)

# ============================================================
# (1) PGS と metadata 読み込み
# ============================================================
meta <- fread(file.path(step1_dir, "ROSMAP_metadata.csv"), na.strings = c("NA",""))
meta[, individualID := trimws(as.character(individualID))]
if (anyDuplicated(meta$individualID)) stop("metadata: individualID 重複")

w <- fread(file.path(raw_dir, paste0(pgs_id, ".txt")),
           sep = "\t", na.strings = c("NA",""))
miss <- setdiff(c("rsID","effect_allele","effect_weight"), names(w))
if (length(miss)) stop("PGS missing columns: ", paste(miss, collapse=", "))
w[, `:=`(rsID = trimws(as.character(rsID)),
         effect_allele = toupper(trimws(as.character(effect_allele))),
         effect_weight = as.numeric(effect_weight))]
w <- w[!is.na(rsID) & rsID != ""][!duplicated(rsID)]
chr_col <- intersect(c("chr_name","chr"), names(w))[1]
w[, chr_int := if (is.na(chr_col)) NA_integer_
               else suppressWarnings(as.integer(
                 gsub("^chr","", as.character(get(chr_col)), ignore.case=TRUE)))]
cat("PGS rsIDs:", nrow(w), "\n")

# ============================================================
# (2) effectSNPs を一括読み込み + allele 整合 + flip + PRS
# ============================================================
sn <- fread(file.path(step1_dir, paste0("effectSNPs_", pgs_id, ".csv")),
            na.strings = c("NA","NaN",""))
req <- c("chr","id","allele1","allele2")
if (!all(req %in% names(sn))) stop("effectSNPs missing leading 4 columns")
sn[, `:=`(allele1 = toupper(trimws(allele1)), allele2 = toupper(trimws(allele2)))]

sample_cols <- setdiff(names(sn), req)
if (!setequal(sample_cols, meta$individualID))
  stop("sample 列 vs metadata individualID 不整合")
setcolorder(sn, c(req, meta$individualID))
sample_cols <- meta$individualID
cat("rows =", nrow(sn), " samples =", length(sample_cols), "\n")

# allele 比較（vectorized）
m_idx     <- match(sn$id, w$rsID)
in_pgs    <- !is.na(m_idx)
sn_in     <- sn[in_pgs]
m_in      <- m_idx[in_pgs]
ea        <- w$effect_allele[m_in]
beta      <- w$effect_weight[m_in]
match_a1  <- (ea == sn_in$allele1)
match_a2  <- (ea == sn_in$allele2)
mismatch  <- !match_a1 & !match_a2
used_rows <- match_a1 | match_a2

cat("present:", length(m_in), " | matches a1:", sum(match_a1),
    " | flipped:", sum(match_a2), " | mismatch:", sum(mismatch), "\n")

# dosage 行列化 → flip → PRS
dose <- as.matrix(sn_in[, ..sample_cols]); mode(dose) <- "numeric"
dose[match_a2, ] <- 2 - dose[match_a2, ]
dose_used  <- dose[used_rows, , drop = FALSE]
dose_used[is.na(dose_used)] <- 0
prs_raw    <- as.numeric(crossprod(dose_used, beta[used_rows]))
names(prs_raw) <- sample_cols
prs_sd     <- sd(prs_raw)
prs_z      <- (
  if (is.na(prs_sd) || prs_sd == 0) rep(0, length(prs_raw))
  else as.numeric(scale(prs_raw))
)
cat("PRS_raw mean =", round(mean(prs_raw), 4),
    " sd =", round(prs_sd, 4), "\n")

# ============================================================
# (3) 出力 3 ファイル
# ============================================================
# SNPs_flipped: allele 整合済 dosage
flipped <- copy(sn_in)
flipped[, (sample_cols) := as.data.table(dose)]
flipped[match_a2, c("allele1","allele2") := list(
  sn_in$allele2[match_a2], sn_in$allele1[match_a2]
)]
fwrite(flipped, paste0("SNPs_flipped_", pgs_id, ".csv"), na = "NA")

# SNPs_weights: 監査用フラグ
n_w  <- nrow(w)
flag <- function(idx) { x <- logical(n_w); x[idx] <- TRUE; x }
fwrite(data.table(
  chr = w$chr_int, id = w$rsID, A1 = w$effect_allele, beta = w$effect_weight,
  present_in_Dosage = flag(m_in),
  flipped           = flag(m_in[match_a2]),
  mismatch          = flag(m_in[mismatch]),
  used_in_PRS       = flag(m_in[used_rows])
), paste0("SNPs_weights_", pgs_id, ".csv"), na = "NA")

# PRS.csv: PRS + metadata
prs_dt <- data.table(individualID = sample_cols,
                     PRS_raw = prs_raw, PRS_z = prs_z)
prs <- merge(prs_dt, meta, by = "individualID", sort = FALSE)
fwrite(prs, paste0("PRS_", pgs_id, ".csv"), na = "NA")
cat("Saved PRS, SNPs_flipped, SNPs_weights\n\n")

# ============================================================
# (4) 診断プロット 7 種
# ============================================================
d <- prs %>%
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
  d2 <- d[!is.na(d[[xv]]) & !is.na(d$PRS_z), , drop = FALSE]
  ggplot(d2, aes(.data[[xv]], PRS_z, fill = .data[[xv]])) +
    geom_boxplot(width = 0.6, alpha = 0.70, outlier.size = 1.5) +
    geom_point(position = position_jitter(width = 0.12, height = 0),
               size = 1.6, alpha = 0.35, color = "black") +
    labs(x = xlab, y = "PRS z-score") + theme_p
}
plot_lm <- function(xv, xlab) {
  d2 <- d[!is.na(d[[xv]]) & !is.na(d$PRS_z), , drop = FALSE]
  ggplot(d2, aes(.data[[xv]], PRS_z)) +
    geom_point(alpha = 0.35, size = 1.6) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
    labs(x = xlab, y = "PRS z-score") + theme_p
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

prefix  <- paste0("Step2_PRSz_plots_", format(Sys.Date(), "%Y%m%d"))
for (nm in names(plots)) {
  ggsave(file.path(workdir, paste0(prefix, "_", nm, ".png")),
         plots[[nm]], width = 6, height = 6, dpi = 300)
}
cat("Saved 7 plots in", workdir, "with prefix:", prefix, "\n")
