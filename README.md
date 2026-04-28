# RFind-G

**Rank-based Fisher enrichment score for individual polygenic risk (RFscore-G)** — a PRS-alternative individual risk score using Running Fisher enrichment of dosage deviation against PGS β ranks.

Validated on ROSMAP (Alzheimer's Disease, N = 1,637, PGS004228):

| | RFscore-G | PRS | Combined |
|---|---|---|---|
| AUC (5-fold CV mean ± SD) | 0.728 ± 0.040 | 0.755 ± 0.054 | 0.758 ± 0.057 |
| Spearman ρ (vs PRS) | 0.925 | — | — |

## Algorithm at a glance

For each individual *j*:
```
RFscore_j = + RF(g1_up,   g2_up_j)    PGS↑ × 個人↑   ┐ concordant
            + RF(g1_down, g2_dn_j)    PGS↓ × 個人↓   ┘
            − RF(g1_up,   g2_dn_j)    逆方向         ┐ discordant
            − RF(g1_down, g2_up_j)    逆方向         ┘
```
where `RF` is the maximum −log10(hypergeometric p) over an N×M grid of top-N PGS ranks vs top-M individual dosage deviation ranks.

Full algorithm spec: [Docs/RFind-G_handoff.md](Docs/RFind-G_handoff.md) (or [Docs/RFind-G_handoff.html](Docs/RFind-G_handoff.html))

## Pipeline

```
RFind-G/
├── Code/
│   ├── Step1_SNPs.R         dosage.gz × 22 → effectSNPs抽出 (~80 min @ 8 cores + pigz)
│   ├── Step2_Analysis.R     allele flip + PRS 計算 + 7 plot (~10 sec)
│   ├── Step3_RFind.R        RFscore-G 計算 + 7 plot (~5 sec, 並列)
│   ├── Step4_figures.R      主要 3 図 (violin / scatter / ROC) (~3 sec)
│   └── make_dummy_data.R    RawData/ の demo 用 dummy data 生成
├── Docs/                    handoff document (md / html / pptx)
├── RawData/                 dummy demo data (300 samples, ~7 MB) — 実運用時は自前データに置換
├── Step1_SNPs/              outputs (gitignored)
├── Step2_Analysis/          outputs (gitignored)
├── Step3_RFind/             outputs (gitignored)
└── Step4_figures/           主要図 (commit 済 example)
```

## Quick start

```r
# 1. project root に setwd (script は getwd() を base_dir として使う)
setwd("/path/to/RFind-G")

# 2. dummy demo data (RawData/) で動作確認
source("Code/Step1_SNPs.R")        # PGS ID プロンプト (デフォルト PGS004228)
source("Code/Step2_Analysis.R")
source("Code/Step3_RFind.R")
# Step4 は ROSMAP 個人 ID の fold ハードコードのため、自前 cohort では fold 再生成必須
# (詳細 Docs/RFind-G_handoff.md §4.3)
```

> **DEMO 注意**: `RawData/` の dummy は random dosage / random clinical なので、
> RFscore / PRS / AUC は意味のない値 (chance level) になります。format / 動作確認のみ。

### 環境

- R ≥ 4.2
- `data.table`, `dplyr`, `ggplot2`, `tidyverse`, `pROC`, `ggforce`, `patchwork`, `glue`, `parallel`, `rstudioapi`
- 推奨 16 GB RAM
- `pigz` (任意、`brew install pigz` で Step1 高速化)

### 入力データ仕様

- **Genotype**: PLINK dosage 形式 (`.dosage.gz`, 22 chr 分割) + `.fam` (詳細は [Docs/RFind-G_handoff.md §3](Docs/RFind-G_handoff.md))
- **PGS**: PGS Catalog 標準 (`rsID`, `effect_allele`, `effect_weight` 必須)
- **表現型**: `individualID` + case/control 列の CSV

> ROSMAP データは [AD Knowledge Portal (synapse.org)](https://www.synapse.org/) から取得できますが、controlled access のため別途承認が必要です。

## 適用例 (依頼ベース)

本 pipeline は元々 ROSMAP (Alzheimer's Disease) で検証されましたが、原理上 **任意の二値表現型 + PGS** に適用可能です。本 repo は DM / 高血圧 / その他疾患への横展開を前提に整備されています。詳細は [Docs/RFind-G_handoff.md](Docs/RFind-G_handoff.md) を参照。

## Citation

未発表 (2026-04 時点)。本 pipeline を利用された場合はその旨ご連絡いただけると助かります。

## Contact

- 村野 智之 (murano.mg@gmail.com)
- 宮川研究室 / 藤田医科大学 総合医科学研究所

## License

MIT (see [LICENSE](LICENSE))
