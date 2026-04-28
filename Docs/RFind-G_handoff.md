# RFind-G 共同解析ご依頼資料

**作成**: 村野 / 2026-04-28

---

## 1. ROSMAP 結果概要

PRS の代替個人リスクスコア **RFscore-G** を開発し、ROSMAP cohort (Alzheimer's Disease, N=1,637 → AD vs NCI 評価対象 N=1,170, PGS004228) で検証。

### 主要結果 (5-fold StratifiedKFold CV, glm logistic regression)

| 指標 | RFscore-G | PRS | Combined (PRS + RFscore-G) |
|---|---|---|---|
| **AUC** (mean ± SD) | **0.728 ± 0.040** | **0.755 ± 0.054** | **0.758 ± 0.057** |
| Cohen's d (AD vs NCI) | 0.86 | 0.99 | — |
| Wilcoxon p | < 0.001 | < 0.001 | — |
| Spearman ρ (vs PRS) | **0.925** | — | — |
| R² (`lm(RFscore_z ~ PRS_z)`) | **0.859** | — | — |

### 図 (Step4 出力)

[Step4_figures/G_violin.png](../Step4_figures/G_violin.png) / [G_scatter_PRS.png](../Step4_figures/G_scatter_PRS.png) / [G_ROC.png](../Step4_figures/G_ROC.png)

### 解釈サマリ

- **PRS と RFscore-G は強相関** (ρ = 0.925, R² = 0.859) → 両者は同じ情報を別パラダイムで再表現
- AUC は **PRS がやや勝るが (0.755 vs 0.728)、Combined で改善は小** (0.758) → 独立予測情報は限定的
- **β スケール不変・閾値不要・Control 基準** という性質で PRS とは異なる方法論的位置づけ
- 別 cohort / 別表現型でも同様の傾向 (高 ρ + AUC 同等) が再現するか、別データで要検証

---

## 2. RFscore-G アルゴリズム

PGS の β を **線形和 (PRS)** ではなく **方向別ランクの濃縮検定 (Running Fisher)** で集計するアプローチ。

| | PRS | RFscore-G |
|---|---|---|
| 集計 | β · dosage の **線形和** | β を順位化、個人 dosage 偏差分布との **濃縮度** |
| 出力 | 連続値 (β スケール依存) | -log10(p) ベースの 4 方向符号付き和 |
| スケール依存 | あり | **β の順位のみ**使用 → 不変 |
| 閾値選択 | clumping/thresholding 必要 | グリッド走査で自動最大化 |

### 2.1 Effective SNP set (universe P) の決定

g1 と g2 の hypergeometric 検定の **母集団 (universe)** となる SNP セット P は、以下のフィルタで決定する:

```
PGS Catalog rsID (8,863 for PGS004228)
        ∩
ROSMAP dosage に存在する rsID (Step1 awk で抽出 → 8,193 個、92.4% マッチ)
        −
allele mismatch (effect_allele が dosage の a1/a2 どちらとも一致しない SNP, ROSMAP では 0 個)
        ↓
   P = 8,193  ← effective SNP universe
```

**P の意義**:
- g1, g2 の **両方が定義される空間**。g1 (PGS β 由来) も g2 (個人 dosage 偏差由来) も全て 1:P の index で表現される
- Hypergeometric の母集団 size として `phyper(overlap-1, N, P-N, M, ...)` の `P` に直接入る
- マッチ率が低い (< 50% 程度) と P が小さくなり、検定 power が低下 → 別 PGS の選定 or rsID liftover (dbSNP build 違い解消) を検討

**マッチ率の典型値**:
- PGS が dataset と同じ genotyping array / imputation panel ベースなら 80-95% (PGS004228 × ROSMAP は 92.4%)
- dbSNP build 違い (例 GRCh37 vs GRCh38) の場合 40-60% に低下することあり
- PGS 公開時の rsID と dataset の rsID 命名規則が異なる場合も低下

### 2.2 g1 シグネチャ (PGS 由来、サンプル非依存)

PGS 重み β を **方向別・絶対値降順** で並べる:
```
g1_up   = { i : β_i > 0 } を |β| 降順で並べた順位列
g1_down = { i : β_i < 0 } を |β| 降順で並べた順位列
```

### 2.3 g2 シグネチャ (個人 j ごと)

各 SNP について **Control 平均 dosage** からの偏差を計算:
```
fc_{ij} = dose_{ij} - mean_Control(dose_i)
g2_up_j = { i : fc_{ij} > 0 } を fc 降順
g2_dn_j = { i : fc_{ij} < 0 } を |fc| 降順
```
**Reference を Control 平均にした理由**: cohort 全平均だと case 比率に応じて baseline が歪む。Control 基準なら「case が control から逸脱している方向 / 強さ」を直接測れる。Control N < 10 のときのみ自動で全 sample 平均にフォールバック。

### 2.4 Running Fisher 検定 RF(g1, g2)

`g1` 上位 N と `g2` 上位 M の重複を hypergeometric 検定で評価し、N×M グリッドで **最大 -log10(p)** を採用 (GSEA 様の non-parametric enrichment):
```
RF(g1, g2) = max_{N ∈ N_seq, M ∈ M_seq}
              -log10( phyper(overlap_NM - 1, N, P-N, M, lower.tail=FALSE) )
```
- `N_seq, M_seq = (50, 100, ..., n)` (`grid_step=50`)
- `P` = effective SNP 数 (PGS ∩ dosage、mismatch 除外後)
- `overlap_NM` = g1 top-N と g2 top-M の rsID 重複数
- グリッド走査により「best cutoff」を後決め → 任意閾値に依らない

### 2.5 最終 RFscore (4 方向の符号付き和)

```
RFscore_j = + RF(g1_up,   g2_up_j)    PGS↑ × 個人↑   ┐ concordant
            + RF(g1_down, g2_dn_j)    PGS↓ × 個人↓   ┘
            − RF(g1_up,   g2_dn_j)    逆方向         ┐ discordant
            − RF(g1_down, g2_up_j)    逆方向         ┘
```
**解釈**: PGS が示すリスク方向に個人偏差が強く濃縮しているほど高スコア。逆方向の濃縮は protective として減点。直感的には「PGS の SNP セットに対する個人 GSEA-like signed enrichment」。

### 2.6 パラメータ (既定)

| 名前 | 既定 | 意味 |
|---|---|---|
| `GRID_STEP` | 50 | Running Fisher の N, M 走査刻み (小さいほど精細、計算コスト増) |

### 2.7 補足: PRS との関係 (重要)

ROSMAP では `lm(RFscore_z ~ PRS_z)` の **R² = 0.859** で、両者は強相関 (Spearman ρ = 0.922)。すなわち:

> **RFscore-G は PRS と独立した予測情報を持つわけではなく、同じ情報を別パラダイムで再表現したもの**

しかし以下の固有性質は PRS と異なる:
- **β スケール不変性**: PGS Catalog 重みの絶対スケール変動 (異なる study 間の β 単位差) に robust
- **閾値不要**: clumping/thresholding (C+T) の任意性に依らない
- **コホート構成頑健**: case 比率による baseline 変動を Control 基準で吸収

→ **別 cohort / 別表現型でも同様の傾向 (高 ρ + AUC 同等) が再現するか** が次の検証ポイント。

---

## 3. 入力データ (ROSMAP の場合)

ROSMAP cohort で実際に使用したデータの仕様。本 pipeline はこの形式を canonical input として実装している。

### 3.1 SNP genotype

[AD Knowledge Portal (synapse.org)](https://www.synapse.org/) から取得した imputed PLINK dosage:

- `AMP-AD_ROSMAP_Rush-Broad_AffymetrixGenechip6_Imputed_chr{1..22}.dosage.gz` (chr 別、22 ファイル、計 ~36 GB)
  - space-delimited、header なし
  - 列: `rsID a1 a2 dose_s1 dose_s2 ... dose_sN`
  - 例: `rs55998931 C T 1.93 1.94 ... 1.93`
  - dose は {0, 1, 2} の連続値 (Affymetrix Genechip 6 → imputation)
  - `a1` = effect-allele 候補 (PGS 側 effect_allele との一致 or flip 判定に使用)
- `AMP-AD_ROSMAP_Rush-Broad_AffymetrixGenechip6_Imputed.fam` (PLINK 標準 6 列): `FID IID PID MID SEX PHENO`
  - `IID` 順序が dosage の sample 列順と完全一致 (1,708 sample)

### 3.2 PGS Catalog ファイル

[PGS Catalog](https://www.pgscatalog.org/) から取得した **PGS004228** (Alzheimer's Disease):

```
#PGS_ID=PGS004228
#trait_reported=Alzheimer's disease
...
rsID    chr_name    effect_allele    effect_weight    other_allele
rs...   ...         ...              ...              ...
```
- tab-delimited, `#` で始まる metadata header
- 必須列: `rsID`, `effect_allele`, `effect_weight`, `chr_name` (or `chr`)
- 8,863 rsID (chr1-22)

### 3.3 表現型 (case/control)

`ROSMAP_clinical.csv` (AD Knowledge Portal):
- `individualID`: dosage の sample IID と biospecimen 経由で対応
- 主要列: `dcfdx_lv` (last valid clinical cognitive diagnosis), `cogdx`, `braaksc`, `ceradsc`, `cts_mmse30_lv`, `msex`, `age_death`, etc.
- **本解析の case/control 定義**: `dcfdx_lv == 1` (NCI, control) vs `dcfdx_lv == 4` (AD, case)
  - N_NCI = 665, N_AD = 505 (合計 N = 1,170 が Step4 評価対象)

### 3.4 補助データ

`ROSMAP_biospecimen_metadata.csv`: specimen-individual 対応表
- ROSMAP では 1 individual に複数 specimen があり得るため必要
- Step1 で fam の `IID` (specimen 単位) → `individualID` 変換に利用
- 1,637 sample が individual 単位に正しく mapping できた (元 1,708 から 71 失敗 = 5.4 cohort 構成上 specimen の重複等)

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

### 4.2 各 Step の入出力ファイル

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

**Step4 の分類評価フロー**:
- `dcfdx_lv ∈ {1, 4}` で AD vs NCI に絞る (N=1,170)
- 5-fold StratifiedKFold で個人を分割 (`y` 比率を保ったまま)
- 各 fold で **training set の z-score で標準化** → **`glm(family = binomial())` (logistic regression)** を fit → test set で予測確率を出す → `pROC::auc` で AUC 計算
- 5 fold の AUC を mean ± SD で報告
- 単 feature (PRS_z 単独 / RFscore_z 単独) の場合 logistic は単調変換なので AUC = 単変量 AUC、**Combined (PRS_z + RFscore_z) のみ logistic regression が 2 feature の最適線形結合を学習**
- ROC は 5 fold pooled prediction probability から作成

**Fold の再現性**: 現在の `folds_G_full` は Python sklearn `StratifiedKFold(seed=42)` と完全一致を確認するためにハードコードされた個人 ID → fold 番号のマッピング。**ROSMAP 個人 ID と一致しない別 cohort では再生成が必要**。例:
```r
set.seed(42)
fold_vec <- caret::createFolds(label, k = 5, list = FALSE)
```
で `Step4_figures.R` の `folds_G_full` 部分を置き換えてください (現在は dummy data でも動くよう auto-regen 機能を実装済)。

