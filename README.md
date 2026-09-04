# paocseq 0.1.0: A python package for mixed cell type annotation following the activation of cells & sequencing.

**paocseq** is a suite of statistical tools for the analysis of multimodal immunosequencing data. 

The purpose of this statistical tool is to quantify single cell data and provide **1)** A sample list of mixed cell types with accompanying data sheets to collate detectable clonotypes that are present in the blood or tissue and **2)** To detect, score and rank cells with optimal or desirable characteristics. Characteristics are defined by the gene expression of cells that are considered of interest because of their response to perturbation. 

# Research applications

The software goals are to combine probability and control theory with existing single cell analysis software for the classification of cells. Ultimately the aim is to provide a flexible framework that can be used to analyze a broad range of functional experiments that produce multi-modal data consisting of multiple different cell types. For a list of applications, see the section **Getting started** for the vignettes or the [documentation](https://github.com/MaeWoods/aocseq/raw/main/aocseq.pdf) listing the functions. 

# Rationale for initial package development

The sofware was initially developed to identify receptor sequences from single cell RNA sequencing coupled with amplification of the receptor sequence. Version 0.1.0 can be used with single cell sequencing data of cells that have been activated through the binding of their T cell receptors (TCRs) to cognate antigens. The tool provides a method to identify full length TCRs with specificity directed toward a target antigen. The long term aim of the package is to develop functionality for the analysis of a broad range of sequencing techniques applicable to polyclonal cells.

**Introduction to TCR sequencing**

T cells are part of the immune system and these cells recognize targets on their cell surface or other cell surfaces via expression of the TCR. The TCR sequence is important to study because T cell specificity depends on a variable region that differs between different people and single cells, equipping T cells with capacity to mount a response against a broad range of targets. T cells that are specific to a particular target can be grouped into clonotypes that share a common CDR3 beta chain and this way, used to estimate the frequency of target specific T cells in the blood. Harnessing this heterogeneity in sequence between T cells for a quantitative analysis of adaptive immunology has broad applicability in immunology and immunotherapy because the clearing of infection and cancer depends on availability of immune cells (including T cells) with capacity to mount a response. 
Immunosequencing is a PCR-based based method that exploits the capacity of high-throughput sequencing technology to characterize tens of thousands of TCR CDR3 chains simultaneously and **aocseq** has been developed to analyse and annotate this data.

Specifically, this is a software package of statistical tools that can be used to trace, analyse, annotate and query clonotypes subject to amplification or reduction in frequency following antigen stimulation or between experimental conditions. **aocseq** has initially been applied to Virus specific T cells (VSTs) because these clinical blood products contain non specific bystander T cells in addition to potent virus specific clonotypes. However, the tool can be adapted to model other barcoded time series frequency data and the accompanying vignette demonstrates how the tool could be used to track clonotypes *in vivo*, using Adaptive's ImmuneAccess database. 

# initial porting using Claude Sonnet 5

The initial Python port of the R package [`aocseq`](.) was generated using Claude desktop.

Where `aocseq` stored data in Seurat objects, `paocseq` uses
[`AnnData`](https://anndata.readthedocs.io/) (via
[`scanpy`](https://scanpy.readthedocs.io/) for I/O), the standard Python
single-cell data container -- the direct analogue of the Seurat object.
Numerical routines that were written by hand in base R or Rcpp/C++ in
`aocseq` (the EM demultiplexing algorithm, the Mahalanobis-distance solver,
the isolation-forest outlier-scoring solver) are re-implemented as
vectorized `numpy` array operations rather than swapped out for third-party
statistics packages, to keep the same algorithmic behaviour and keep the
dependency footprint minimal: `numpy`, `scipy`, `pandas`, `anndata`,
`scanpy`, `matplotlib`.

## Install

```bash
cd paocseq
uv venv --python 3.13
uv sync
# or, without uv:
pip install -e .
```

## Quick start

```python
import paocseq as pq

# 1. Import 10x feature-barcode matrices + VDJ contig annotations
adata_list = pq.combine_data(
    gex_path=["SeqData/gex/Control", "SeqData/gex/CMV"],
    vdj_path=[
        "SeqData/vdj/Control/all_contig_annotations.csv",
        "SeqData/vdj/CMV/all_contig_annotations.csv",
    ],
    marker_gene=["IFNG", "TNF"],
)

# 2. Threshold marker genes to "high"/"unassigned" per cell
for adata in adata_list:
    pq.threshold_marker_genes(adata, ["IFNG", "TNF"], threshold_cutoff=0.975)

# 3. Summarize clonotype frequency + marker-gene positivity across samples
summary = pq.annotate_cell_types(adata_list, goi="IFNG", path="SummaryTable.csv")

# 4. (Optional) CITE-seq / hashtag demultiplexing
demuxed = pq.combine_data(
    gex_path=["SeqData/gex/CD5"],
    marker_gene=["IFNG", "CRTAM", "GZMB"],
    demultiplex=True,
    demultiplex_index=[8, 9, 10],
    nameshashtags=["donor_1", "donor_2", "donor_3"],
    n_ht_per_sample=3,
)

# 5. Build a reference signature and classify clonotypes / cells against it
reference = pq.make_reference(adata_list[0], gene_list=["IFNG", "TNF"])
scored = pq.classify_cells(adata_list[1], reference, gene_list=["IFNG", "TNF"], distance=0)
```

See `examples/getting_started.py` for a runnable end-to-end example against
synthetic data (no real 10x data required), and the module docstrings in
`src/paocseq/*.py` for the mapping from each Python function back to its
original R source file/function name.

## Module map (R -> Python)

| R file (`aocseq`)        | R function(s)                                   | Python module            |
|---------------------------|--------------------------------------------------|---------------------------|
| `Import.R`                | `CombineData`                                    | `paocseq.io.combine_data` |
| `Import.R`                | `EStep`, `MStep`, `GMMDemux`, `PhenotypeMarkers` | `paocseq.demux`           |
| `Annotation.R`            | `AnnotateCellTypes` (aka `AnnotateClonotypes`)   | `paocseq.annotate.annotate_cell_types` |
| `Classification.R`        | `ClassifyCellTypes`, `ConstructTree`/`IsoForest` | `paocseq.classify`, `paocseq._solvers` |
| `Distances.R`             | `ClassifyCells`, `AddDistances`, `PercentOutlier`| `paocseq.distances`       |
| `Expression.R`            | `GetSpecificCells`, `GetGeneSignature`           | `paocseq.expression`      |
| `Genomics.R`              | `RankTCRs`                                       | `paocseq.genomics`        |
| `Plotting.R`              | `QCPlot`, `UMAPReduce`, `SegmentPlot`, `SaveHeatmap` | `paocseq.plotting`    |
| `References.R`            | `MakeReference`                                  | `paocseq.references`      |
| `solvers.cpp`             | `GetMahalanobis`, `isoForest`/`tree`             | `paocseq._solvers`        |

### What was intentionally *not* ported

* **`DEsingle`** (`Classification.R`) is a vendored copy of the Bioconductor
  zero-inflated-negative-binomial DE test. Rather than reimplement that
  statistical model from scratch, `paocseq.classify.classify_cell_types`
  offers `method="wilcoxon"`, using `scanpy.tl.rank_genes_groups` for the
  same "is this clonotype specific for the gene of interest" question.
* **`GetTCRs`** and **`MakeReceptor`** (`Genomics.R`) parse Cell Ranger
  `all_contig_annotations.json` files and make live network calls to
  UniProt to stitch full-length receptor sequences. These are one-off,
  network-dependent utilities orthogonal to the matrix/array-based analysis
  this package focuses on (10x feature-barcode matrices +
  `all_contig_annotations.csv`), and were left out of this first pass.
* `SegmentPlot`'s circos/chord diagram (via the R `circlize` package) is
  reproduced as a paired bar chart in `paocseq.plotting.segment_plot`
  instead, to avoid adding a circos-plotting dependency; and
  `SaveHeatmap`'s dataset-specific hard-coded column indices were
  generalized into `paocseq.plotting.expression_heatmap`.

This is a first, testable pass at the port -- please try it against your
real 10x + VDJ data and let me know what breaks or what's missing so the
next iteration can fix it.
