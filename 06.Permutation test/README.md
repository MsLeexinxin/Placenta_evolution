Statistical Permutation Test for Gene Overlaps

We performed permutation test to check if the overlaps between cancer and hemochorial up-regulated genes for each cell type were significantly different from those expected at random. To this end, we used sample function in R to generate simulated gene sets in the human genome by randomly selecting genes of equal number to the observed up-regulated genes in hemochorial, and we replicated this process 10,000 times. We compared the number of overlaps between the observed hemochorial up-regulated and cancer up-regulated genes with the overlaps between the simulated hemochorial up-regulated and cancer up-regulated genes, and calculated the statistical significance of P values (i.e., the probability that a higher number of overlaps would be observed by chance).

Script: permutation_test.r

Usage:
    # Run the permutation test for a specific cell type against a directory of cancer gene lists
    Rscript permutation_test.r <hemo_genes.txt> <bg_genes.txt> <cancer_dir> <out_summary.csv> <out_intersections.csv> [iterations] [seed]
Example:
    Rscript permutation_test.r \
    data/EVT_upregulated_genes.txt \
    data/EVT_all_genes.txt \
    data/cancer_upragulated_genes/ \
    output/EVT_summary.csv \
    output/EVT_cancer_overlaps.csv
    
Parameters:
<hemo_genes.txt>: A single-column text file containing the up-regulated genes in the hemochorial cell type (e.g., T cells or EVT).
<bg_genes.txt>: A single-column text file containing all background genes expessed in hemochorial species.
<cancer_dir>: A directory containing cancer up-regulated gene lists. Each file should represent a specific cancer type and contain a single-column list of its up-regulated genes.
<out_summary.csv>: The output file path for the statistical summary and P-values.
<out_intersections.csv>: The output file path for the extracted overlapping genes matrix.
[iterations]: (Optional) The number of random samplings to perform. Default is 10000.
[seed]: (Optional) Random seed for reproducibility. Default is 42.

Output Files:
This script generates 2 files:
<out_summary.csv>: Contains the overlap counts, extreme sweep counts, P-values, and adjusted P-values (FDR/BH) for each cancer type.
<out_intersections.csv>: A matrix table listing the exact overlapping genes identified between the hemochorial cells and each cancer type.

