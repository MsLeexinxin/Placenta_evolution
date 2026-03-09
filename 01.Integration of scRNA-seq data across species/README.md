Prepare the Longest Transcripts for Protein Sequence Analysis Across Nine Species

This README outlines the data preparation pipeline to generate gene-level embeddings for cross-species analysis. The process involves extracting the longest protein-coding transcripts from the NCBI genome annotations of nine species, mapping protein IDs to gene symbols, generating protein embeddings using the ESM-2 model, and finally converting them to gene-level embeddings.

Step 1: Prepare Representative Protein Sequences

Download the genome annotation feature tables for the following nine species and extract the longest protein-coding transcript for each gene.

List of Species (NCBI Assembly Accession):
•   PetGlider (Petaurus breviceps): PetGlider_PUasm1.0 (GCF_028583685.1)

•   Human (Homo sapiens): GRCh38.p14 (GCF_000001405.40)

•   Cynomolgus monkeys (Macaca fascicularis): MFA1912RKSv2 (GCF_012559485.2)

•   Mouse (Mus musculus): GRCm39 (GCF_000001635.27)

•   Chinese Tree Shrew (Tupaia chinensis): TupChi_1.0 (GCF_000334495.1)

•   Dog (Canis lupus familiaris): ROS_Cfam_1.0 (GCF_014441545.1)

•   Rambouillet Sheep (Ovis aries): ARS_UI_Ramb_v3.0 (GCF_016772045.2)

•   Cattle (Bos taurus): ARS1.2 (GCF_001704415.2)

•   Pig (Sus scrofa): Sscrofa11.1 (GCF_000003025.6)

Example: Human Genome (GRCh38.p14)

1.  Download and Filter Feature Table:
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_feature_table.txt.gz
    zcat GCF_000001405.40_GRCh38.p14_feature_table.txt.gz | awk '{if($1=="mRNA")print}' | cut -f 13,15 | sort -n | uniq > pro2sym.txt
    
    This command downloads the human genome feature table, filters for mRNA entries, extracts the protein ID (column 13) and Gene ID (column 15), sorts them, removes duplicates, and saves the mapping to pro2sym.txt.

2.  Map Protein IDs to Gene IDs in the FASTA file:
    Use the provided script com.aa.sedname.pl to transfer the identifiers in your representative protein FASTA file (human_hg38.representative.pep.fa) using the mapping file (pro2sym.txt).
    perl com.aa.sedname.pl human_hg38.representative.pep.fa pro2sym.txt > human_hg38.pep.all.fa
    

3.  Filter the FASTA File:
    Run the clean_fasta.py script to standardize/clean the FASTA sequences.
    python clean_fasta.py --data_path human_hg38.pep.all.fa --save_path human_hg38.pep.all_clean.fa
    

Step 2: Generate Protein Embeddings with ESM-2

Use the Evolutionary Scale Modeling (ESM) framework to compute embeddings for the cleaned protein sequences.

1.  Model Preparation: Download the pre-trained ESM-2 model checkpoint esm2_t48_15B_UR50D.pt from the official https://github.com/facebookresearch/esm.

2.  Extract Embeddings: Run the extract.py script. This will generate per-residue or per-protein embeddings (here, the mean pooling of residue embeddings is used).
    python3 -u extract.py esm2_t48_15B_UR50D.pt human_hg38.pep.all_clean.fa hg38_pep.all_clean.fa_esm1b --include mean --truncate
    

Step 3: Convert Protein Embeddings to Gene-Level Embeddings

Since multiple protein isoforms can map to the same gene, we aggregate protein embeddings to the gene level.

1.  Create Gene-to-Protein Mapping (JSON):
    Generate a JSON file that maps each gene symbol to a list of its corresponding protein IDs from the processed FASTA file.
    python ../../scripts/map_gene_symbol_to_protein_ids.py --fasta_path human_hg38.pep.all_clean.fa --save_path hg38.gene_symbol_to_protein_ID.json
    

2.  Aggregate to Gene-Level Embeddings:
    Use the mapping file and the protein embeddings directory to create a single embedding vector per gene.
    python convert_protein_embeddings_to_gene_embeddings.py --embedding_dir hg38_pep.all_clean.fa_esm1b --gene_symbol_to_protein_ids_path hg38.gene_symbol_to_protein_ID.json --embedding_model ESM2 --save_path hg38_gene_symbol_to_embedding_ESM2.pt
    

Step 4: Prepare Data for SATURN Analysis

Repeat Steps 1 through 3 for all nine species. After completion, you should have for each species:
•   An h5ad file containing the gene expression data (prepared separately).

•   A .pt file containing the gene symbol to ESM-2 embedding mapping (e.g., hg38_gene_symbol_to_embedding_ESM2.pt).

Step 5: Run SATURN for Cell Type Annotation

Once all species-specific data files are ready, execute the SATURN analysis script.
sbatch -D `pwd` -N 1 --mem=60g -p gpu4 -n 1 -c 1 --export=ALL,input_data="your_input",nmac="6" run_SATURN.aggl.cell_type.sh

Note: Replace "your_input" with the path to your consolidated input data configuration.

Step 6: Visualize SATURN Results

After the SATURN analysis is complete, use the provided draw_embedding.py script to generate visualizations (e.g., 2D embeddings like UMAP/t-SNE plots) based on the output files.
python draw_embedding.py ./multiple_seeds_results/saturn_results/9species_batch_label_split_seed_0.h5ad ./multiple_seeds_results/saturn_results/9species_genes_to_macrogenes.pkl NA NA


Follow these steps sequentially for each species to prepare the necessary files for the cross-species gene embedding analysis and to generate the final result visualizations.