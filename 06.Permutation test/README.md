# Statistical Permutation Test for Gene Overlaps

We performed permutation tests to determine if the overlaps between cancer and hemochorial up-regulated genes for each cell type were significantly different from those expected at random. This analysis uses random sampling to generate null distributions for comparison with observed overlap counts.

## Method Overview
To assess statistical significance, we used the `sample` function in R to generate simulated gene sets in the human genome by randomly selecting genes equal in number to the observed up-regulated genes in hemochorial cell types. This process was repeated 10,000 times. We then compared the number of overlaps between the observed hemochorial up-regulated and cancer up-regulated genes with the overlaps between the simulated hemochorial up-regulated and cancer up-regulated genes. The statistical significance (P-values) represents the probability that a higher number of overlaps would be observed by chance.

## Script: `permutation_test.r`

### Usage
```bash
# Run the permutation test for a specific cell type against a directory of cancer gene lists
Rscript permutation_test.r <hemo_genes.txt> <bg_genes.txt> <cancer_dir> <out_summary.csv> <out_intersections.csv> [iterations] [seed]
```

### Example
```bash
Rscript permutation_test.r \
  data/EVT_upregulated_genes.txt \
  data/EVT_all_genes.txt \
  data/cancer_upragulated_genes/ \
  output/EVT_summary.csv \
  output/EVT_cancer_overlaps.csv
```

### Parameters
- **`<hemo_genes.txt>`**: A single-column text file containing the up-regulated genes in the hemochorial cell type (e.g., T cells or EVT).
- **`<bg_genes.txt>`**: A single-column text file containing all background genes expressed in hemochorial species.
- **`<cancer_dir>`**: A directory containing cancer up-regulated gene lists. Each file should represent a specific cancer type and contain a single-column list of its up-regulated genes.
- **`<out_summary.csv>`**: The output file path for the statistical summary and P-values.
- **`<out_intersections.csv>`**: The output file path for the extracted overlapping genes matrix.
- **`[iterations]`**: (Optional) The number of random samplings to perform. Default is 10,000.
- **`[seed]`**: (Optional) Random seed for reproducibility. Default is 42.

## Output Files
This script generates 2 files:

1. **`<out_summary.csv>`**: Contains the following information for each cancer type:
   - Overlap counts between observed hemochorial genes and cancer genes
   - Extreme sweep counts from permutations
   - P-values (raw statistical significance)
   - Adjusted P-values (FDR/Benjamini-Hochberg correction)

2. **`<out_intersections.csv>`**: A matrix table listing the exact overlapping genes identified between the hemochorial cells and each cancer type.
