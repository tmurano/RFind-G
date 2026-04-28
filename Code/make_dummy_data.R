# ============================================================
# make_dummy_data.R : RawData/ に format demo 用 dummy data を生成
#
# 目的: 共同研究者が clone してすぐ Step1-4 の format / 動作確認できる
#       small-scale ダミーデータを作成。実際のヒト由来情報は一切含まない。
#
# 出力 (RawData/):
#   AMP-AD_..._chr{1..22}.dosage.gz  : PGS rsID × N_DUMMY 個人 (random dosage)
#   AMP-AD_..._Imputed.fam           : N_DUMMY dummy individuals
#   ROSMAP_biospecimen_metadata.csv  : N_DUMMY specimen-individual mapping
#   ROSMAP_clinical.csv              : N_DUMMY random clinical records
#   PGS004228.txt                    : (real PGS Catalog file, 公開データ)
#
# 注意:
#   - dosage 値は Uniform(0,2) random で実際の SNP 多型を反映しない
#   - clinical 値も random、実在の患者と無関係
#   - 結果として RFscore / PRS / AUC は random ≈ chance level となる
#     (Step4 出力は意味のない数値、format demo のみ)
# ============================================================

suppressPackageStartupMessages({ library(data.table) })

set.seed(42)

# ---- パス ----
base_dir <- normalizePath(getwd(), mustWork = TRUE)
src_pgs  <- file.path(base_dir, "RawData", "PGS004228.txt")  # 既存ならそれ使う
out_dir  <- file.path(base_dir, "RawData")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# PGS file が無い場合は退避先からコピー (一度きり)
if (!file.exists(src_pgs)) {
  fallback <- file.path(dirname(base_dir), "RFind-G_real_RawData", "PGS004228.txt")
  if (file.exists(fallback)) {
    file.copy(fallback, src_pgs)
    cat("Copied PGS file from", fallback, "->", src_pgs, "\n")
  } else {
    stop("PGS004228.txt が見つかりません。実在する PGS Catalog file を ",
         out_dir, "/ に置いてから再実行してください。")
  }
}

# ---- パラメータ ----
N_DUMMY <- 300L    # dummy individual 数 (~10 MB の demo data になる)
DUMMY_PREFIX <- "dummy_"

# ---- PGS file から rsID + chr を取得 ----
pgs <- fread(src_pgs, sep = "\t", na.strings = c("NA",""), skip = "rsID")
pgs[, chr_int := suppressWarnings(as.integer(
  gsub("^chr", "", as.character(chr_name), ignore.case = TRUE)))]
pgs <- pgs[!is.na(rsID) & substr(rsID, 1, 2) == "rs" & !is.na(chr_int) & chr_int %in% 1:22]
cat("PGS rsIDs (chr1-22):", nrow(pgs), "\n")

# ---- dummy individual ID 群 ----
ids <- sprintf("%s%03d", DUMMY_PREFIX, seq_len(N_DUMMY))

# ---- (1) dosage.gz × 22 chr ----
alleles <- c("A","T","C","G")
dos_tpl <- "AMP-AD_ROSMAP_Rush-Broad_AffymetrixGenechip6_Imputed_chr%d.dosage.gz"

for (chr in 1:22) {
  rs_chr <- pgs$rsID[pgs$chr_int == chr]
  n_snp  <- length(rs_chr)
  if (n_snp == 0) {
    # 空 chr 用に 1 行だけ dummy (script が空 chr を想定してない場合の保険)
    rs_chr <- paste0("rs_dummy_chr", chr, "_001")
    n_snp <- 1L
  }
  a1 <- sample(alleles, n_snp, replace = TRUE)
  a2 <- vapply(seq_len(n_snp), function(i) sample(setdiff(alleles, a1[i]), 1L),
               character(1))
  # dosage matrix: random Uniform(0,2)
  dose <- matrix(sprintf("%.4f", runif(n_snp * N_DUMMY, min = 0, max = 2)),
                 nrow = n_snp, ncol = N_DUMMY)
  # space-delimited rows: rsID a1 a2 dose1 dose2 ...
  lines <- paste(rs_chr, a1, a2,
                 apply(dose, 1L, paste, collapse = " "))
  out_path <- file.path(out_dir, sprintf(dos_tpl, chr))
  con <- gzfile(out_path, "w")
  writeLines(lines, con)
  close(con)
  cat(sprintf("  chr %2d : %5d SNPs -> %s\n", chr, n_snp, basename(out_path)))
}

# ---- (2) fam (PLINK 6 列) ----
fam_path <- file.path(out_dir, "AMP-AD_ROSMAP_Rush-Broad_AffymetrixGenechip6_Imputed.fam")
fam <- data.table(
  FID   = ids,
  IID   = ids,
  PID   = 0L,
  MID   = 0L,
  SEX   = sample(1:2, N_DUMMY, replace = TRUE),  # PLINK: 1=M, 2=F
  PHENO = -9L
)
fwrite(fam, fam_path, sep = " ", col.names = FALSE)
cat("Saved", basename(fam_path), "(", nrow(fam), "rows )\n")

# ---- (3) biospecimen metadata (specimenID == individualID for simplicity) ----
bios_path <- file.path(out_dir, "ROSMAP_biospecimen_metadata.csv")
bios <- data.table(
  specimenID   = ids,
  individualID = ids,
  organ        = "brain",
  tissue       = "DLPFC"
)
fwrite(bios, bios_path)
cat("Saved", basename(bios_path), "(", nrow(bios), "rows )\n")

# ---- (4) clinical (random, ROSMAP 風 column 構成) ----
clin_path <- file.path(out_dir, "ROSMAP_clinical.csv")
# dcfdx_lv: 1=NCI 50%, 4=AD 35%, その他 15% (合計 N_DUMMY)
n_nci   <- round(N_DUMMY * 0.50)
n_ad    <- round(N_DUMMY * 0.35)
n_other <- N_DUMMY - n_nci - n_ad
dcfdx_pool <- c(rep(1L, n_nci),
                rep(4L, n_ad),
                rep(c(2L,3L,5L,6L), length.out = n_other))
dcfdx_pool <- sample(dcfdx_pool)
stopifnot(length(dcfdx_pool) == N_DUMMY, !anyNA(dcfdx_pool))
clin <- data.table(
  individualID  = ids,
  msex          = sample(0:1, N_DUMMY, replace = TRUE),
  age_death     = sprintf("%.1f", runif(N_DUMMY, 70, 95)),
  braaksc       = sample(0:6, N_DUMMY, replace = TRUE),
  ceradsc       = sample(1:4, N_DUMMY, replace = TRUE),
  cogdx         = sample(1:6, N_DUMMY, replace = TRUE),
  dcfdx_lv      = dcfdx_pool,
  cts_mmse30_lv = round(runif(N_DUMMY, 0, 30))
)
fwrite(clin, clin_path)
cat("Saved", basename(clin_path), "(", nrow(clin), "rows )\n")
cat("  dcfdx_lv == 1 (NCI):", sum(clin$dcfdx_lv == 1), "\n")
cat("  dcfdx_lv == 4 (AD): ", sum(clin$dcfdx_lv == 4), "\n")

cat("\nDummy data generated in:", out_dir, "\n")
cat("Total size:", format(sum(file.info(list.files(out_dir, full.names = TRUE))$size) / 1024^2,
                          digits = 3), "MB\n")
