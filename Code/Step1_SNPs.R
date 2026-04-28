# ============================================================
# Step1_SNPs.R : ROSMAP dosage.gz から PGS 対象 SNP を抽出 + metadata 構築
#
# 入力 (Step1_SNPs/):
#   AMP-AD_*.fam, AMP-AD_chr*.dosage.gz × 22
#   ROSMAP_biospecimen_metadata.csv, ROSMAP_clinical.csv 
#   <PGS_ID>.txt (PGS Catalog 形式)
# 出力:
#   ROSMAP_metadata.csv          (specimenID + individualID + clinical)
#   effectSNPs_<PGS_ID>.csv      (chr,id,allele1,allele2 + 個人 dosage)
#
# 高速化: pigz 並列解凍 + awk フィルタ + 22 chr を mclapply 並列。
#         pigz が無ければ gzip にフォールバック。
# ============================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(parallel)
})

pgs_id <- rstudioapi::showPrompt("PGS ID", "PGS ID:", "PGS004228")
if (is.null(pgs_id) || trimws(pgs_id) == "") stop("PGS ID 未入力")
pgs_id <- trimws(pgs_id)
cat("PGS ID:", pgs_id, "\n")

# ---- パス & パラメータ ----
base_dir  <- path.expand("~/Miyakawa Lab Dropbox/Murano Tomoyuki/RFind-G_260428")
input_dir <- file.path(base_dir, "RawData")     # 読み込み (dosage.gz, fam, csv, txt)
workdir   <- file.path(base_dir, "Step1_SNPs")  # 書き出し (metadata, effectSNPs, tmp_)
if (!dir.exists(workdir)) dir.create(workdir, recursive = TRUE)
setwd(workdir)
dos_tpl  <- "AMP-AD_ROSMAP_Rush-Broad_AffymetrixGenechip6_Imputed_chr%d.dosage.gz"
chrs     <- 1:22
n_par    <- 6L  # 並列数（Mac 8 cores → gzip 律速で ~6 が最適）

GZ_DEC <- if (nzchar(Sys.which("pigz"))) "pigz -dc" else "gzip -dc"
cat("Decompressor:", GZ_DEC,
    if (GZ_DEC == "gzip -dc") "  (高速化なら brew install pigz)" else "", "\n")

# tmp_dir: /tmp 直下に固定 (Dropbox 配下を回避、R tempdir() の不安定挙動も回避)
# 完了後に自動削除はしない (debug 用に残す)。再実行時は冒頭で unlink。
tmp_dir <- file.path("/tmp", paste0("RFind-G_step1_", pgs_id))
if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE, force = TRUE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = TRUE, mode = "0755")
if (!dir.exists(tmp_dir)) stop("Failed to create tmp_dir: ", tmp_dir)
cat("tmp_dir:", tmp_dir, "\n")

# ============================================================
# (A) FAM / biospecimen / clinical → ROSMAP_metadata.csv
# ============================================================
clean_id <- function(x) { x <- trimws(as.character(x)); x[x == ""] <- NA; x }

fam  <- fread(file.path(input_dir, "AMP-AD_ROSMAP_Rush-Broad_AffymetrixGenechip6_Imputed.fam"),
              header = FALSE)
setnames(fam, c("FID","IID","PID","MID","SEX","PHENO"))
bios <- fread(file.path(input_dir, "ROSMAP_biospecimen_metadata.csv"))
clin <- fread(file.path(input_dir, "ROSMAP_clinical.csv"))

fam$IID           <- clean_id(fam$IID)
bios$specimenID   <- clean_id(bios$specimenID)
bios$individualID <- clean_id(bios$individualID)
clin$individualID <- clean_id(clin$individualID)

iid_all <- fam$IID                  # 元 FAM 順（dosage 列順と一致）
n_fam   <- length(iid_all)

id_map <- bios %>%
  filter(!is.na(specimenID), !is.na(individualID)) %>%
  distinct(specimenID, individualID)
if (nrow(id_map %>% count(specimenID) %>% filter(n > 1)) > 0)
  stop("FATAL: specimenID -> 複数 individualID")

mapped     <- setNames(id_map$individualID, id_map$specimenID)[iid_all]
keep_pos   <- which(!is.na(mapped))
ind_strict <- mapped[keep_pos]
if (anyDuplicated(ind_strict)) stop("FATAL: individualID 重複")

cat("FAM samples:", n_fam, " | kept (mapped):", length(ind_strict), "\n")

# clinical を結合
clin_df <- as.data.frame(clin, stringsAsFactors = FALSE)
if ("age_death" %in% colnames(clin_df))
  clin_df$age_death <- as.numeric(gsub("\\+", "", clin_df$age_death))
if (anyDuplicated(clin_df$individualID))
  stop("FATAL: clinical の individualID 重複")

metadata <- data.frame(
  specimenID = iid_all[keep_pos], individualID = ind_strict,
  stringsAsFactors = FALSE
) %>% left_join(clin_df, by = "individualID")
write.csv(metadata, "ROSMAP_metadata.csv", row.names = FALSE)
cat("Saved ROSMAP_metadata.csv (rows =", nrow(metadata), ")\n")

# ============================================================
# (B) PGS rsID リスト + awk 用ヘルパファイル
# ============================================================
score_df <- read.table(file.path(input_dir, paste0(pgs_id, ".txt")),
                       header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, comment.char = "#",
                       quote = "", fill = TRUE)
rs_set <- unique(trimws(as.character(score_df$rsID)))
rs_set <- rs_set[!is.na(rs_set) & rs_set != "" & substr(rs_set, 1, 2) == "rs"]
cat("PGS rsID count:", length(rs_set), "\n")

writeLines(rs_set,                       file.path(tmp_dir, "rsids.txt"))
writeLines(as.character(keep_pos),       file.path(tmp_dir, "strict_idx.txt"))
writeLines(paste(c("chr","id","allele1","allele2", ind_strict), collapse = ","),
           file.path(tmp_dir, "header.csv"))

# awk: rsID ハッシュで filter + strict_idx で列選択 → CSV 出力
writeLines('
BEGIN {
  FS = " "
  while ((getline l < rsfile)  > 0) keep[l] = 1; close(rsfile)
  ni = 0
  while ((getline l < idxfile) > 0) { ni++; idx[ni] = l + 3 }
  close(idxfile)
}
{
  if ($1 in keep) {
    out = chr "," $1 "," $2 "," $3
    for (i = 1; i <= ni; i++) out = out "," $idx[i]
    print out
  }
}', file.path(tmp_dir, "filter.awk"))

# ============================================================
# (C) chr ごと並列実行（gzip|awk）+ 進捗バー
# ============================================================
run_chr <- function(chr) {
  dos    <- file.path(input_dir, sprintf(dos_tpl, chr))
  outcsv <- file.path(tmp_dir, sprintf("chr_%02d.csv", chr))
  cmd <- sprintf("%s %s | awk -v rsfile=%s -v idxfile=%s -v chr=%d -f %s > %s",
                 GZ_DEC, shQuote(dos),
                 shQuote(file.path(tmp_dir, "rsids.txt")),
                 shQuote(file.path(tmp_dir, "strict_idx.txt")),
                 chr, shQuote(file.path(tmp_dir, "filter.awk")),
                 shQuote(outcsv))
  list(chr = chr, ret = system(cmd))
}

cat("Running awk on", length(chrs), "chrs (P =", n_par, ") ...\n")
t0 <- Sys.time()
jobs <- vector("list", length(chrs))
done <- rep(FALSE, length(chrs))
results <- vector("list", length(chrs))
nxt <- 1L; running <- 0L
n_expected <- length(rs_set); BAR <- 30L

while (any(!done)) {
  while (nxt <= length(chrs) && running < n_par) {
    jobs[[nxt]] <- mcparallel(run_chr(chrs[nxt]),
                              name = sprintf("chr%02d", chrs[nxt]))
    running <- running + 1L; nxt <- nxt + 1L
  }
  for (i in seq_along(chrs)) {
    if (done[i] || is.null(jobs[[i]])) next
    r <- mccollect(jobs[[i]], wait = FALSE, timeout = 0)
    if (!is.null(r)) {
      results[[i]] <- r[[1]]; done[i] <- TRUE; running <- running - 1L
    }
  }
  total <- sum(vapply(seq_along(chrs), function(i) {
    f <- file.path(tmp_dir, sprintf("chr_%02d.csv", chrs[i]))
    if (file.exists(f)) length(readLines(f, warn = FALSE)) else 0L
  }, integer(1)))
  pct      <- min(1, total / max(n_expected, 1))
  filled   <- as.integer(round(pct * BAR))
  bar      <- paste0(strrep("=", filled), strrep(".", BAR - filled))
  elapsed  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  rate     <- if (elapsed > 5) total / elapsed else 0
  eta_str  <- if (rate > 0 && total < n_expected)
                sprintf("%.0f min", (n_expected - total) / rate / 60)
              else "--"
  cat(sprintf("\r  [%4.0fs] [%s] %3.0f%%  chr %2d/%d  SNPs %5d/~%d  ETA %s",
              elapsed, bar, pct * 100, sum(done), length(chrs),
              total, n_expected, eta_str))
  flush.console()
  Sys.sleep(2)
}
cat("\n")

bad <- which(vapply(results, function(r) is.null(r) || r$ret != 0L, logical(1)))
if (length(bad)) stop("FAILED chrs: ", paste(chrs[bad], collapse = ", "))
cat("All chrs done in",
    round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "sec\n")

# ============================================================
# (D) 連結 → effectSNPs_<PGS_ID>.csv
# ============================================================
total_rows <- 0L
for (chr in chrs) {
  f <- file.path(tmp_dir, sprintf("chr_%02d.csv", chr))
  n <- as.integer(system2("wc", c("-l", shQuote(f)), stdout = TRUE) |>
                    sub(pattern = "^\\s*([0-9]+).*$", replacement = "\\1"))
  cat(sprintf("CHR %2d : rows = %d\n", chr, n))
  total_rows <- total_rows + n
}
cat("Total rows written:", total_rows, "\n")

snps_out <- paste0("effectSNPs_", pgs_id, ".csv")
system(sprintf(
  "cat %s %s > %s",
  shQuote(file.path(tmp_dir, "header.csv")),
  paste(shQuote(file.path(tmp_dir, sprintf("chr_%02d.csv", chrs))), collapse = " "),
  shQuote(snps_out)
))
cat("Saved:", snps_out, "\n")
