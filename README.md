# TNBC target discovery (work in progress)

Differential expression of basal-like breast tumors (n=190) vs normal
breast tissue (n=113) from TCGA-BRCA, to prioritise candidate targets
for RNA-based therapeutics.

## Method
TCGAbiolinks for data, DESeq2 for differential expression.

## Limitations
- Basal-like (PAM50) is a proxy for triple-negative, not identical.
- Bulk RNA-seq mixes tumor and immune cells.
- Computational prioritisation only; no experimental validation.

## Next steps
Check expression in normal organs (GEPIA2), survival, DepMap
dependency, then siRNA design.
