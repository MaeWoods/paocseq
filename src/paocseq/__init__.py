"""paocseq: a Python port of the R package ``aocseq``.

aocseq traces functional-marker-gene co-expression across heterogeneous
single cells (typically T cells), by combining 10x Genomics gene-expression
feature-barcode matrices with VDJ contig annotations, optionally
demultiplexing hashtag/CITE-seq data, and classifying clonotypes/cell types
by marker-gene expression relative to a reference.

Where the R package stored data in Seurat objects, paocseq uses
:class:`anndata.AnnData` (the Python single-cell ecosystem's equivalent
container), read via :mod:`scanpy`. The numerical routines that were
hand-written in base R or Rcpp/C++ in aocseq (the EM demultiplexing
algorithm, the Mahalanobis-distance solver, the isolation-forest solver) are
re-implemented directly as vectorized numpy array operations in
:mod:`paocseq.demux` / :mod:`paocseq._solvers`, matching the original
algorithms rather than substituting off-the-shelf equivalents, per the
package's original design.

Typical pipeline::

    import paocseq as pq

    adata_list = pq.combine_data(
        gex_path=["SeqData/gex/Control", "SeqData/gex/CMV"],
        vdj_path=["SeqData/vdj/Control/all_contig_annotations.csv",
                  "SeqData/vdj/CMV/all_contig_annotations.csv"],
        marker_gene=["IFNG", "TNF"],
    )

    for adata in adata_list:
        pq.threshold_marker_genes(adata, ["IFNG", "TNF"], threshold_cutoff=0.975)

    summary = pq.annotate_cell_types(adata_list, goi="IFNG")
"""

from .io import combine_data, read_10x_sample
from .demux import gmm_demux, phenotype_markers, e_step, m_step, fit_two_component_em
from .classify import threshold_marker_genes, classify_cell_types
from .annotate import annotate_cell_types
from .distances import classify_cells, add_distances, percent_outlier, pf_log1p_pf
from .expression import get_specific_cells, get_gene_signature
from .genomics import rank_tcrs
from .plotting import qc_plot, umap_reduce, segment_plot, expression_heatmap
from .references import make_reference
from . import _solvers as solvers

__version__ = "0.1.0"

__all__ = [
    "combine_data",
    "read_10x_sample",
    "gmm_demux",
    "phenotype_markers",
    "e_step",
    "m_step",
    "fit_two_component_em",
    "threshold_marker_genes",
    "classify_cell_types",
    "annotate_cell_types",
    "classify_cells",
    "add_distances",
    "percent_outlier",
    "pf_log1p_pf",
    "get_specific_cells",
    "get_gene_signature",
    "rank_tcrs",
    "qc_plot",
    "umap_reduce",
    "segment_plot",
    "expression_heatmap",
    "make_reference",
    "solvers",
]
