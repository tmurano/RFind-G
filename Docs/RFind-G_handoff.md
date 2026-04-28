# RFind-G 共同解析ご依頼資料

**作成**: 村野 (宮川研究室) / 2026-04-28
**宛先**: Genequest 技術担当者様

---

## 1. 依頼概要

宮川研では PRS の代替個人リスクスコア **RFscore-G** を開発し、ROSMAP cohort (Alzheimer's Disease, N=1,637) で検証しました。同手法を **貴社が保有する DM (糖尿病) / 高血圧の genotype + 表現型データ** に適用し、PRS 単独との比較評価をご依頼いたします。**コードのチューニング・最適化のご依頼ではありません**。本資料の pipeline をそのまま (or 必要最小の入力アダプタのみ追加で) 流していただく想定です。

### ご提供いただきたい結果

| # | 内容 (DM, HTN 各 cohort で 1 set) |
|---|---|
| 1 | `RFinD_G_scores_<PGS>.csv` (個人 RFscore-G + clinical) |
| 2 | `PRS_<PGS>.csv` (個人 PRS + clinical) |
| 3 | AUC ± SD (RFscore-G / PRS / Combined), Spearman ρ (vs PRS), Cohen's d, Wilcoxon p |
| 4 | 主要 3 図 (`G_violin.png`, `G_scatter_PRS.png`, `G_ROC.png`) |
| 5 | 使用 PGS Catalog ID + case/control 定義 (簡潔な README で OK) |

ROSMAP 参考結果 (PGS004228, AD vs NCI, N=1,170): AUC RFscore-G **0.730 ± 0.041** / PRS **0.755 ± 0.054** / Spearman ρ **0.922** / R² (RFscore_z ~ PRS_z) **0.859**。

---

## 2. RFscore-G アルゴリズム

PGS の β を **線形和 (PRS)** ではなく **方向別ランクの濃縮検定 (Running Fisher)** で集計するアプローチ。スケール不変・閾値不要・cohort 構成頑健、という性質を持つ。

| | PRS | RFscore-G |
|---|---|---|
| 集計 | β · dosage の **線形和** | β を順位化、個人 dosage 偏差分布との **濃縮度** |
| 出力 | 連続値 (β スケール依存) | -log10(p) ベースの 4 方向符号付き和 |
| スケール依存 | あり | **β の順位のみ**使用 → 不変 |
| 閾値選択 | clumping/thresholding 必要 | グリッド走査で自動最大化 |

### 2.1 g1 シグネチャ (PGS 由来、サンプル非依存)

PGS 重み β を **方向別・絶対値降順** で並べる:
```
g1_up   = { i : β_i > 0 } を |β| 降順で並べた順位列
g1_down = { i : β_i < 0 } を |β| 降順で並べた順位列
```

### 2.2 g2 シグネチャ (個人 j ごと)

各 SNP について **Control 平均 dosage** からの偏差を計算:
```
fc_{ij} = dose_{ij} - mean_Control(dose_i)
g2_up_j = { i : fc_{ij} > 0 } を fc 降順
g2_dn_j = { i : fc_{ij} < 0 } を |fc| 降順
```
**Reference を Control 平均にした理由**: cohort 全平均だと case 比率に応じて baseline が歪む。Control 基準なら「case が control から逸脱している方向 / 強さ」を直接測れる。Control N < 10 のときのみ自動で全 sample 平均にフォールバック。

### 2.3 Running Fisher 検定 RF(g1, g2)

`g1` 上位 N と `g2` 上位 M の重複を hypergeometric 検定で評価し、N×M グリッドで **最大 -log10(p)** を採用 (GSEA 様の non-parametric enrichment):
```
RF(g1, g2) = max_{N ∈ N_seq, M ∈ M_seq}
              -log10( phyper(overlap_NM - 1, N, P-N, M, lower.tail=FALSE) )
```
- `N_seq, M_seq = (50, 100, ..., n)` (`grid_step=50`)
- `P` = effective SNP 数 (PGS ∩ dosage、mismatch 除外後)
- `overlap_NM` = g1 top-N と g2 top-M の rsID 重複数
- グリッド走査により「best cutoff」を後決め → 任意閾値に依らない

### 2.4 最終 RFscore (4 方向の符号付き和)

```
RFscore_j = + RF(g1_up,   g2_up_j)    PGS↑ × 個人↑   ┐ concordant
            + RF(g1_down, g2_dn_j)    PGS↓ × 個人↓   ┘
            − RF(g1_up,   g2_dn_j)    逆方向         ┐ discordant
            − RF(g1_down, g2_up_j)    逆方向         ┘
```
**解釈**: PGS が示すリスク方向に個人偏差が強く濃縮しているほど高スコア。逆方向の濃縮は protective として減点。直感的には「PGS の SNP セットに対する個人 GSEA-like signed enrichment」。

### 2.5 パラメータ (既定で OK)

| 名前 | 既定 | 意味 |
|---|---|---|
| `GRID_STEP` | 50 | Running Fisher の N, M 走査刻み (小さいほど精細、計算コスト増) |
| `TOP_K` | 10,000 | 個人 g2 SNP 上限 (P が小さい場合は無効、underflow 抑制用) |

### 2.6 補足: PRS との関係 (重要)

ROSMAP では `lm(RFscore_z ~ PRS_z)` の **R² = 0.859** で、両者は強相関 (Spearman ρ = 0.922)。RF 残差の AUC = 0.533 ≈ chance。すなわち:

> **RFscore-G は PRS と独立した予測情報を持つわけではなく、同じ情報を別パラダイムで再表現したもの**

しかし以下の固有性質は PRS と異なる:
- **β スケール不変性**: PGS Catalog 重みの絶対スケール変動 (異なる study 間の β 単位差) に robust
- **閾値不要**: clumping/thresholding (C+T) の任意性に依らない
- **コホート構成頑健**: case 比率による baseline 変動を Control 基準で吸収

→ **DM / HTN でも同様の傾向 (高 ρ + AUC 同等) が再現するか** が本依頼の主要検証ポイント。

---

## 3. 必要な入力データ

### 3.1 SNP genotype

**貴社データ形式 (VCF / PLINK bed / dosage / BGEN 等) をまずご教示ください**。本 pipeline は ROSMAP の **PLINK dosage (.dosage.gz, 22 chr 分割)** を canonical input として実装していますが、貴社形式に応じて Step1 のローダ部分を調整します (本体ロジックは形式非依存)。

ROSMAP 形式 (参考):
- `*.dosage.gz` (chr 別、space-delimited、header なし):
  ```
  rsID a1 a2 dose_s1 dose_s2 ... dose_sN
  rs55998931 C T 1.93 1.94 ... 1.93
  ```
  - dose は {0, 1, 2} の連続値 (imputed dosage)
  - `a1` = effect-allele (PGS 側 effect_allele との一致 or flip 判定に使用)
- `*.fam` (PLINK 標準 6 列): `FID IID PID MID SEX PHENO`
  - `IID` 順序が dosage の sample 列順と **完全一致** すること

VCF / BGEN の場合、`PLINK2 --recode A-transpose` 等で dosage 形式へ変換してから入力するのが最短です。

### 3.2 PGS Catalog ファイル

DM / 高血圧の PGS は **貴社にて PGS Catalog から選定** をお願いします。形式は PGS Catalog 標準 (tab-delimited, `#` で始まる metadata header):
```
#PGS_ID=PGS00xxxx
#trait_reported=Type 2 diabetes
...
rsID    chr_name    effect_allele    effect_weight    other_allele
rs7903146    10    T    0.342    C
...
```
- 必須列: `rsID`, `effect_allele`, `effect_weight`
- chr 列名は `chr_name` または `chr` を許容

### 3.3 表現型 (case/control)

CSV、必須列:
- `individualID`: dosage の sample IID と一致する識別子
- **case/control 列**: 二値。**DM / 高血圧の case/control 判定は貴社の臨床定義 (ICD コード / 服薬歴 / 検査値 等) に従って決定** してください。スクリプト内の参照箇所 (`Step3_RFind.R` 80-91 行目: ROSMAP では `meta$dcfdx_lv == 1` で control を抽出) を、貴社の判定列名 + control 値に置換いただきます。

その他、サブグループ解析 (年齢 / 性別 / 既往歴 等) で使いたい列があれば任意で含めてください (Step2/3 の plot 軸変数として利用可)。

### 3.4 補助データ (任意)

dosage の sample ID が **specimen 単位** (1 個人 = 複数 sample) の場合、`(specimenID, individualID)` の対応 CSV を別途ご用意ください (ROSMAP では `ROSMAP_biospecimen_metadata.csv`)。1 個人 = 1 sample なら不要。

---

## 4. パイプライン構成 + 実行手順

### 4.1 ディレクトリ構成

```
RFind-G_260428/
├── Code/
│   ├── Step1_SNPs.R         入力読み込み + PGS rsID 抽出 (重い: ~80 分 @ 8 cores)
│   ├── Step2_Analysis.R     allele flip + PRS 計算 + 7 plot (~10 秒)
│   ├── Step3_RFind.R        RFscore-G 計算 + 7 plot (~5 秒, 並列)
│   └── Step4_figures.R      主要 3 図 (~3 秒)
├── RawData/                 入力 (genotype, fam, csv, PGS txt)
├── Step1_SNPs/              出力: ROSMAP_metadata.csv, effectSNPs_<PGS>.csv
├── Step2_Analysis/          出力: PRS_<PGS>.csv, SNPs_flipped, SNPs_weights, plots
├── Step3_RFind/             出力: RFinD_G_scores_<PGS>.csv, plots
└── Step4_figures/           出力: G_violin.png, G_scatter_PRS.png, G_ROC.png
```

### 4.2 各 Step の入出力

| Step | 入力 | 主な出力 |
|---|---|---|
| 1 | dosage.gz × 22, fam, biospecimen.csv, clinical.csv, `<PGS>.txt` | `effectSNPs_<PGS>.csv` (chr,id,a1,a2 + 個人 dosage 列), `metadata.csv` |
| 2 | Step1 出力 + `<PGS>.txt` | `PRS_<PGS>.csv` (individualID, PRS_raw, PRS_z + clinical), `SNPs_flipped`, `SNPs_weights` (PGS 監査表), 7 PNG |
| 3 | Step1/2 出力 | `RFinD_G_scores_<PGS>.csv` (individualID, RFscore, RFscore_z + clinical), 7 PNG |
| 4 | Step2/3 出力 | `G_violin.png`, `G_scatter_PRS.png`, `G_ROC.png` |

### 4.3 実行

```r
# project root に setwd するだけで全 script 動きます (base_dir は getwd() ベース)
setwd("/path/to/RFind-G_workspace")
source("Code/Step1_SNPs.R")        # PGS ID プロンプト
source("Code/Step2_Analysis.R")
source("Code/Step3_RFind.R")
source("Code/Step4_figures.R")
```

**環境**: R ≥ 4.2 / `data.table, dplyr, ggplot2, tidyverse, pROC, ggforce, patchwork, glue, parallel, rstudioapi` / 推奨 16 GB RAM / `pigz` (任意, `brew install pigz` で Step1 高速化)。

**Step4 の 5-fold fold は ROSMAP のハードコード** (Python sklearn 一致確認用) のため、**貴社 cohort では再生成必須**です。例:
```r
set.seed(42)
fold_vec <- caret::createFolds(label, k = 5, list = FALSE)
```
で `Step4_figures.R` の `folds_G_full` 部分を置き換えてください。

---

## 5. 補足

**Q. ちょうど良い PGS が複数ある。どれを選ぶ?**
→ (1) 貴社 cohort と **祖先構成が近い** publication ベース、(2) **N_SNP ~10⁴ 程度** (極端に小さい / 大きいものは除外)、(3) **multi-ancestry validation 報告あり** のものを推奨。最終判断は貴社にお任せします。複数 PGS で比較いただいても OK。

**Q. 数値的問題 (Inf / underflow) は?**
→ ROSMAP (P=8,193, N=1,637) では発生せず。万一 hypergeom p が 0 にアンダーフロー (`-log10(0) = Inf`) したら `TOP_K` を下げて (例 10,000 → 5,000) g2 サイズを制限してください。

**Q. Step3 の Reference (Control) の定義は?**
→ 貴社 cohort の **case/control 列で control に分類された個人の dosage 平均** を使用。Control N < 10 で自動的に全 sample 平均にフォールバックします (`Step3_RFind.R` 85-91 行目)。

**Q. Step4 図の軸変数は何にすべき?**
→ ROSMAP では `dcfdx_lv == 1 (NCI) vs == 4 (AD)` を case/control としています。DM / HTN では貴社の臨床判定二値変数 (例 `t2d_status`, `htn_status`) で `g_an <- g_raw %>% filter(...)` の filter 条件を書き換えてください。

**Q. Step1 の SNP マッチ率が低い (< 50%) 場合は?**
→ PGS の rsID が dbSNP build 違いで genotype 側と不一致の可能性。`Step1_SNPs.R` で extract された SNP 数と PGS rsID 総数の比 (Step1 完了時に表示) を確認し、明らかに低い場合は別 PGS の選定 or rsID liftover が必要です。

**Q. 結果ファイル受け渡し**
→ CSV (Step1-3 全部 or 最低限 PRS と RFinD_G_scores)、PNG (Step4 の 3 図 + 任意で Step2/3 の診断 plot)、統計サマリ (text or md) を zip 化して共有 drive / メール添付で。**生 genotype の送付は不要**です。
