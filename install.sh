source /vast/xa73pav/tools/miniconda3/etc/profile.d/conda.sh
conda create -n mcfDNA_profiler_v0.1 python=3.12
## Note: installing one by one helped to solve the "Solving environment" issue. (Solving environment took too long when installing all at once)
conda install snakemake seqkit kraken2 bracken krakentools krona minimap2 samtools bedtools biopython circos matplotlib numpy pandas
conda activate mcfDNA_profiler_v0.1
# Optional SLURM executor:
# conda install -c conda-forge -c bioconda snakemake-executor-plugin-slurm

