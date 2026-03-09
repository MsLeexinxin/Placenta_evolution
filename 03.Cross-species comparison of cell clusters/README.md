Cross-Species Single-Cell Integration and Homology Analysis Pipeline Using SAMap



This README provides detailed instructions for performing cross-species single-cell integration and homology analysis using SAMap. The pipeline involves pairwise BLAST alignments, SAMap integration, mapping score extraction, and visualization of cell type correspondence across multiple species.



Step 1: Prepare Pairwise BLAST Results



First, create a directory structure for storing BLAST alignment results between all species pairs.



Directory Structure



Create a folder named 01.blastout. Inside this folder, create subdirectories for each species pair, where the directory name consists of two-letter species codes. For example, HOMU represents the human (HO) and cynomolgus monkeys (MA) pair.



01.blastout/

├── HOMA/

├── MAMU/

├── TUDO/

└── ... (other species pairs)





Species Codes



The following two-letter species codes are used in this analysis:

• HO: Human (Homo sapiens)



• MA: Cynomolgus monkeys (Macaca fascicularis)



• MU: Mouse (Mus musculus)



• TU: Tree shrew (Tupaia chinensis)



• DO: Dog (Canis lupus familiaris)



• OV: Sheep (Ovis aries)



• GO: Goat (Capra hircus)



• SU: Pig (Sus scrofa)



BLAST Alignment Commands



For each species pair directory, you need to generate two reciprocal BLAST alignment files. Each file should be named in the format {species1}\_to\_{species2}.txt.



Example: Human-Macaca (HOMA) analysis

\# Set up sequence directory and output directory

seqdir="path/to/pep/sequences"

outdir="./01.blastout/HOMA"



\# Create output directory

mkdir -p $outdir

mkdir -p logs



\# Human (query) vs Macaca(database)

blastp -query $seqdir/HO.pep.fa \\

&nbsp;      -db $seqdir/Ma.pep.fa \\

&nbsp;      -outfmt 6 \\

&nbsp;      -out $outdir/HO\_to\_MA.txt \\

&nbsp;      -num\_threads $cpus \\

&nbsp;      -max\_hsps 1 \\

&nbsp;      -evalue 1e-6 >> logs/blast.HOtoMA.log 2>\&1



\# Macaca (query) vs Human (database)

blastp -query $seqdir/MA.pep.fa \\

&nbsp;      -db $seqdir/HO.pep.fa \\

&nbsp;      -outfmt 6 \\

&nbsp;      -out $outdir/MA\_to\_HO.txt \\

&nbsp;      -num\_threads $cpus \\

&nbsp;      -max\_hsps 1 \\

&nbsp;      -evalue 1e-6 >> logs/blast.MAtoHO.log 2>\&1





Important: The query species comes first in the filename (e.g., HO\_to\_MA.txt indicates human proteins were used as queries against the macaca protein database).



Repeat this process for all species pairs of interest. Each pair directory should contain two reciprocal BLAST files.



Step 2: Run SAMap Integration



With the BLAST results prepared, you can now perform SAMap integration for each species pair. SAMap is a method for aligning single-cell RNA-seq data across species.



Prepare Input Files



Ensure you have the following for each species:

• H5AD files: Single-cell expression data in h5ad format (e.g., HO.h5ad, MU.h5ad)



• Cluster information: Ensure the h5ad files contain cluster labels in the seurat\_clusters column (or other specified annotation column)



SAMap Command



Run SAMap for each species pair using the following command structure:

python /share/home/zhanglab/user/lixin/placenta\_analysis/snRNA-analysis/Samap/run\_samap.py \\

&nbsp; -a HO \\              # First species code

&nbsp; -A MA \\              # Second species code

&nbsp; -k seurat\_clusters \\ # Cluster column for first species

&nbsp; -K seurat\_clusters \\ # Cluster column for second species

&nbsp; -f ./h5ad\_file/HO.h5ad \\     # h5ad file for first species

&nbsp; -F ./h5ad\_file/MU.h5ad \\     # h5ad file for second species

&nbsp; -b ./01.blastout \\           # BLAST results directory

&nbsp; -o HOMA                     # Output directory prefix



Parameter Explanation:

• -a / -A: Species codes (must match the BLAST directory names)

• -k / -K: Annotation column names in the h5ad files containing cell type/cluster labels

• -f / -F: Paths to h5ad files for the two species

• -b: Path to the parent directory containing BLAST results

• -o: Output prefix for SAMap results



Repeat this step for all species pairs to generate SAMap integration results.



Step 3: Generate SAMap Mapping Scores

After running SAMap for each species pair, you need to extract the mapping scores from the results.

For each species pair, run the script to generate the mapping score results: 

python twospecies.readpkl\_plot.py

\# This generates: HOMA.results.txt, MAMU.results.txt, TUDO.results.txt, DOOV.results.txt, OVGO.results.txt, GOSU.results.txt



Each text file ({species\_pair}.results.txt) containing the SAMap mapping scores between cell types/clusters of the two species.



Step 4: Create Combined Matrix and Generate Sankey Diagram

After obtaining mapping scores for all species pairs, combine them into a single matrix and visualize the relationships using a Sankey diagram.



Step 4.1: Create Combined Matrix

Combine all pairwise mapping score files into a single matrix:



python step1.creat.merge.HO\_MA\_MU\_TU\_DO\_OV\_GO\_SU.py



Input files required: HOMA.results.txt, MAMU.results.txt, TUDO.results.txt, DOOV.results.txt, OVGO.results.txt, GOSU.results.txt

Output: This script generates a combined matrix file that contains mapping scores across all eight species.



Step 4.2: Generate Sankey Diagram Visualization



Create a Sankey diagram to visualize the cell type correspondence relationships across all species:

python step2.sanky.HO\_MA\_MU\_TU\_DO\_OV\_GO\_SU.py



Output: The script generates a Sankey diagram visualization (typically as a PDF or PNG file) showing the flow and strength of cell type correspondences across the eight species.



Workflow Summary



1\. Prepare BLAST results for all species pairs in 01.blastout/

2\. Run SAMap integration for each species pair using the BLAST results

3\. Extract mapping scores from SAMap results for each pair

4\. Combine all mapping scores into a single matrix

5\. Visualize cross-species relationships using a Sankey diagram



Note: Ensure all required Python scripts (HOMU.readpkl\_plot.py, step1.creat.merge.py, step2.sanky.py) are in your working directory or accessible in your PATH before running the pipeline.

