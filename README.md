[🇬🇧 English](README.md) | [🇨🇳 中文](README-Chinese.md)

# jwtools (Jian Wu's tools) <img src="man/figures/logo_hex.png" align="right" height="150" alt="jwtools logo" />

[![Website](https://img.shields.io/badge/website-jianwu--atlas.github.io-blue.svg)](https://jianwu-atlas.github.io/index.html)
![R-CMD-check](https://img.shields.io/badge/R-package-blue.svg)
![Version](https://img.shields.io/badge/version-0.2.0-brightgreen.svg)
![License](https://img.shields.io/badge/license-MIT-lightgrey.svg)
[![Views](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/jianwu-atlas/jwtools/main/traffic_report/badge_views.json)](https://github.com/jianwu-atlas/jwtools/graphs/traffic)
[![Clones](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/jianwu-atlas/jwtools/main/traffic_report/badge_clones.json)](https://github.com/jianwu-atlas/jwtools/graphs/traffic)

`Jian Wu`'s Personal collection of reusable R utility functions for computational biology workflows.

## Motivation & Design Philosophy

`jwtools` is designed to **streamline downstream analysis and visualization** of
**bulk RNA-seq, bulk ATAC-seq, ChIP-seq, spatial 10x (Visium/Xenium), and
single-cell/single-nucleus sequencing data**. Rather than integrating every
possible bioinformatics method at once, the package grows organically: every
time a useful helper function is written for a real project.
it gets added here as a documented, tested, reusable function — instead of
being copy-pasted between scripts and slowly diverging into inconsistent
versions.

The long-term goal is to **integrate commonly used bioinformatics routines
into modular and reusable scripts**, so that:

- the same QC / normalization / visualization logic is applied consistently
  across projects and assay types;
- large single-cell/spatial objects (e.g. `Seurat` objects with 200k+ cells)
  can be saved, loaded, and manipulated in a **memory-efficient** way;
- new team members (or future-me, six months later) can understand *why*
  a function exists, not just *how* to call it.

> **Current status (v0.2.0):** in addition to the original workspace I/O
> utilities, this release adds a donor-level **cell-type proportion analysis
> toolkit** (`ct_proportion_analysis()`), a **dataset-origin naming
> standardizer** (`rename_dataset_origin()`), and a **dual PDF+PNG
> figure-saving helper** (`save_dual()`) — see below.
> Additional modules for bulk RNA-seq, ATAC-seq/ChIP-seq, spatial
> transcriptomics, and further single-cell analysis are planned — see
> [Roadmap](#roadmap).

---

## Table of Contents

- [jwtools (Jian Wu's tools) ](#jwtools-jian-wus-tools-)
  - [Motivation \& Design Philosophy](#motivation--design-philosophy)
  - [Table of Contents](#table-of-contents)
  - [Currently Implemented Functions](#currently-implemented-functions)
    - [Why `.qs2` instead of base R `.RData`?](#why-qs2-instead-of-base-r-rdata)
  - [Installation](#installation)
    - [Option 1 — Local installation (works immediately, no GitHub required)](#option-1--local-installation-works-immediately-no-github-required)
    - [Option 2 — Install from GitHub (recommended — keeps server \& laptop in sync)](#option-2--install-from-github-recommended--keeps-server--laptop-in-sync)
  - [Usage](#usage)
  - [Roadmap](#roadmap)
  - [Adding New Functions (Developer Guide)](#adding-new-functions-developer-guide)
  - [Dependencies](#dependencies)
  - [Python 工具](#python-工具)
  - [Repo Traffic](#repo-traffic)
  - [Author \& Contact](#author--contact)
  - [License](#license)

---

## Currently Implemented Functions

| Function | Purpose | Equivalent to your old code |
|---|---|---|
| `qs_save_workspace()` | Serialize **all objects in the current environment** into a single `.qs2` file | `save(list = ls(all = TRUE), file = "xxx.RData")` |
| `qs_load_workspace()` | Restore **all objects** from a `.qs2` file back into the environment | `load("xxx.RData")` |
| `ct_proportion_analysis()` | Compute **donor-level cell-type proportions**, run Kruskal–Wallis (overall) + pairwise Wilcoxon (BH-adjusted) tests across a grouping variable (e.g. age group), and generate a faceted boxplot+jitter figure with optional dataset/batch-origin point overlay for batch-effect diagnostics | Manual `rstatix`/`ggpubr` scripts copy-pasted per project |
| `rename_dataset_origin()` | Standardize the `dataset_origin` column in Seurat object metadata to a consistent `"FirstAuthor_Year"` citation format | Manual `case_when()`/`recode()` lookup tables rewritten each time |
| `save_dual()` | Save a `ggplot2` plot **simultaneously as `.pdf` and `.png`**, with a single function call and consistent width/height/DPI settings | Two separate `ggsave()` calls with duplicated arguments |

### Why `.qs2` instead of base R `.RData`?

The [`qs2`](https://cran.r-project.org/package=qs2) format uses a
multi-threaded, high-compression serialization backend that is
**substantially faster and produces smaller files** than base R's
`save()`/`load()` — a difference that becomes very noticeable when the
environment contains large single-cell objects (e.g. a `Seurat` object with
a `scale.data` slot for tens of thousands of cells × hundreds of features).
This directly matters for reproducible, HPC/SLURM-based single-cell
workflows where I/O time and disk quota are real constraints.

---

## Installation

### Option 1 — Local installation (works immediately, no GitHub required)

Copy the `jwtools/` folder to any location on your server or laptop, then:

```r
install.packages("devtools")   # skip if already installed
devtools::install_local("path/to/jwtools")
```

### Option 2 — Install from GitHub (recommended — keeps server & laptop in sync)

1. Create a GitHub repository (e.g. `jwtools`).
2. Push the entire folder content, keeping the package structure intact at
   the repository root (`DESCRIPTION`, `NAMESPACE`, `R/`, `man/`, etc.).
3. On any machine:

```r
install.packages("remotes")   # skip if already installed
remotes::install_github("jianwu-atlas/jwtools")
```

Whenever a function is modified and pushed to GitHub, simply re-run
`remotes::install_github("jianwu-atlas/jwtools")` on any other machine to pull
the latest version — no manual file copying needed.

---

## Usage

```r
library(jwtools)

# =======================================================================
# Workspace I/O
# =======================================================================

# 1. Save the ENTIRE current workspace
#    (equivalent to base R: save(list = ls(all = TRUE), file = ...))
#    Recommended at natural checkpoints in a long analysis pipeline
#    (e.g. right after cell-type annotation, before starting DE testing).
qs_save_workspace("core_WT_YMO_workspace.qs2", nthreads = 14)

# 2. Exclude specific large objects that are already saved separately
#    (avoids duplicating disk space for objects like a full Seurat object
#    that has its own dedicated .qs2/.rds checkpoint).
qs_save_workspace(
  "core_WT_YMO_workspace.qs2",
  nthreads = 14,
  exclude = c("immune.combined")
)

# 3. Exclude a batch of objects by regex pattern
#    (e.g. transient ggplot objects "p_xxx", scratch variables "tmp_xxx",
#    or intermediate figure objects "fig1", "fig2", ...).
qs_save_workspace(
  "core_WT_YMO_workspace.qs2",
  nthreads = 14,
  exclude_pattern = "^p_|^tmp_|^fig[0-9]"
)

# 4. Restore the workspace on the same or a different machine
#    (e.g. moving from an HPC/SLURM node to a local Mac for plotting).
qs_load_workspace("core_WT_YMO_workspace.qs2", nthreads = 14)


# =======================================================================
# Cell-type proportion analysis (donor-level, with pairwise stats
# and dataset/batch-origin overlay)
# =======================================================================

res <- ct_proportion_analysis(
  seurat_obj      = immune.combined,
  celltype_col    = "Manuscript_Identity",
  sample_col      = "sample",
  group_col       = "AgeGroup_short",
  group_levels    = c("Y", "M", "O"),
  dataset_col     = "dataset_origin",
  celltype_order  = c("EC", "Epi", "B", "DC", "Mac", "Mono", "NK", "T", "Fib", "SMC"),
  celltype_colors = cols,               # reuse the same color vector as your main UMAP
  output_prefix   = "28_ExtFig1F"
)

# Access individual results:
res$plot            # main boxplot+jitter ggplot object (also auto-saved as PDF+PNG)
res$prop_df          # donor-level proportion table (also auto-saved as .csv)
res$kw_test          # Kruskal-Wallis overall test per cell type
res$pairwise_test    # pairwise Wilcoxon test (BH-adjusted) per cell type
res$stat_test         # full stats table incl. bracket coordinates for plotting


# =======================================================================
# Standardizing dataset_origin metadata
# =======================================================================

# Convert raw/inconsistent dataset labels into "FirstAuthor_Year" format
seurat_obj$dataset_origin <- rename_dataset_origin(seurat_obj$dataset_origin)


# =======================================================================
# Saving figures in both PDF and PNG
# =======================================================================

p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()

save_dual(
  filename_prefix = "figures/scatter_wt_mpg",   # no extension needed
  plot            = p,
  width           = 6,
  height          = 5,
  dpi             = 300
)
# → writes figures/scatter_wt_mpg.pdf AND figures/scatter_wt_mpg.png
```

Full function documentation is available via:

```r
?qs_save_workspace
?qs_load_workspace
?ct_proportion_analysis
?rename_dataset_origin
?save_dual
```

---

## Roadmap

The following modules are **planned but not yet implemented**. They reflect
the intended long-term scope of `jwtools` across the assay types used in
the author's ongoing atlas projects. Contributions/additions will be listed
here first, then moved to "Currently Implemented Functions" once merged.

| Planned module | Intended scope |
|---|---|
| **Bulk RNA-seq utilities** | `DESeq2` design/contrast wrappers, GO/KEGG batch enrichment, family-level count aggregation, publication-ready volcano plots |
| **Bulk ATAC-seq / ChIP-seq utilities** | Peak annotation (gene), differential accessibility wrappers, TSS/metagene signal profile plotting, KRAB-ZNF binding site analysis |
| **Spatial transcriptomics (10x Visium/Xenium)** | Standardized loading & QC, spatial feature/module-score overlay plotting, cell-type deconvolution visualization |
| **Single-cell / single-nucleus utilities** | QC filtering helpers, Harmony/Seurat v5 integration wrappers, marker-based lineage annotation, signed-log `AddModuleScore()` wrapper, pseudobulk DE aggregation |
| **Shared visualization utilities** | Themed `ggplot2`/`ComplexHeatmap`/`pheatmap` wrappers with consistent color palettes across projects (`save_dual()` for PDF+PNG export already implemented in v0.2.0) |

> Each module will follow the same principle: a function is added here
> **only after being used and validated in a real analysis**, ensuring every
> exported function in `jwtools` is battle-tested rather than speculative.

---

## Adding New Functions (Developer Guide)

To keep the package structure and documentation self-consistent:

1. Create a new `.R` file under `R/` (e.g. `R/slim_and_save.R`) and write
   the function there. Add `roxygen2` comments directly above the function
   (lines starting with `#'`), including `@param`, `@return`, and `@export`
   tags — and where relevant, a short note on the **biological or technical
   rationale** (e.g. *why* a signed-log transform is used instead of a plain
   log transform for module scores with negative values).

2. Regenerate documentation and test locally:

```r
devtools::document()   # auto-generates/updates NAMESPACE and man/*.Rd from roxygen comments
devtools::load_all()   # reloads the package in the current R session for testing
devtools::check()      # optional: runs R CMD check to catch structural issues
```

3. Commit and push to GitHub. On any other machine, re-run
   `remotes::install_github("jianwu-atlas/jwtools")` to pull the update.

This workflow means `NAMESPACE` and help files never need to be hand-edited —
`devtools::document()` generates them automatically from the roxygen comments.

---

## Dependencies

Core dependencies (required, from `DESCRIPTION`):

- [`qs2`](https://cran.r-project.org/package=qs2) — fast, multi-threaded
  serialization backend used by `qs_save_workspace()` / `qs_load_workspace()`
- `dplyr`, `tibble`, `tidyr` — data wrangling for `ct_proportion_analysis()`
  and `rename_dataset_origin()`
- `ggplot2`, `ggpubr`, `ggh4x`, `ggrepel`, `RColorBrewer`, `scales`, `grid` —
  plotting backbone for `ct_proportion_analysis()` and `save_dual()`
- `rstatix` — Kruskal–Wallis / pairwise Wilcoxon statistical testing in
  `ct_proportion_analysis()`
- `stats`, `utils` — base R statistical/utility functions

Optional (for running package tests):

- `testthat` (>= 3.0.0)

```r
install.packages(c(
  "qs2", "dplyr", "ggh4x", "ggplot2", "ggpubr", "ggrepel",
  "RColorBrewer", "rstatix", "scales", "tibble", "tidyr"
))
```

Future modules (see [Roadmap](#roadmap)) will introduce additional
dependencies as they are implemented, e.g. `Seurat` (>= 5.0.0), `DESeq2`,
`ChIPseeker`, `Squidpy`/`Giotto` interfaces, and `ComplexHeatmap`. These
will be documented per-function via roxygen2 `@importFrom` tags and
reflected in `DESCRIPTION` at that time.

---

## Python 工具

`python/` 目录提供配套的 Python 可视化模块，详见
[python/README.md](python/README.md)。

---

## Repo Traffic

Views and clones are tracked automatically via a scheduled GitHub Action
(see `.github/workflows/traffic-tracker.yml`) that pulls the GitHub Traffic
API every 10 days and accumulates the history in `traffic_report/`,
overcoming the API's native 14-day retention limit.

![Monthly Traffic Trend](traffic_report/monthly_trend.png)

---

## Author & Contact

**Jian Wu (Jack Ng)**
Cardio-Pulmonary Institute (CPI) & Max Planck Institute for Heart and Lung Research & DZL DataLung School

- ✉️ Email: [Jian.Wu@mpi-bn.mpg.de](mailto:Jian.Wu@mpi-bn.mpg.de)
- 🌐 Homepage: [jianwu-atlas.github.io](https://jianwu-atlas.github.io/index.html)
- GitHub: [github.com/jianwu-atlas](https://github.com/jianwu-atlas)
- Google Scholar: [scholar.google.com/citations?user=-pYIKQkAAAAJ](https://scholar.google.com/citations?user=-pYIKQkAAAAJ&hl)
- ORCID: [0000-0003-4720-2374](https://orcid.org/0000-0003-4720-2374)
- LinkedIn: [linkedin.com/in/jackng833](https://www.linkedin.com/in/jackng833/)
- X (Twitter): [@jackng831](https://x.com/jackng831)

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE)
file for details.
