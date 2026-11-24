# mcfDNA profiler

The workflow counts reads and read-length thresholds and classifies each
sample with Kraken2.

## Local execution

Activate the existing environment, adjust [config.yaml](config.yaml), and run
locally from this directory:

```bash
conda activate mcfDNA_profiler_v0.1
snakemake --configfile config.yaml --cores 24 --resources mem_mb=256000 --dry-run
snakemake --configfile config.yaml --cores 24 --resources mem_mb=256000
```

Each rule specifies its own `threads`. The `--cores` argument sets the total
limit for concurrent jobs and should match the cores allocated to the job.
A rule cannot use more threads than the total limit.

Kraken2 requests 24 threads and `mem_mb=256000` per sample by default. For local
execution, `--resources mem_mb=256000` limits scheduling to one Kraken2 job at
a time; adjust the limit to the memory allocated to your job. Snakemake's
resource limit controls scheduling, not actual memory use. Do not start a
job unless sufficient memory is allocated. Set `kraken2_mem_mb` according to
the database's requirements; the rule fails before starting Kraken2 if a global
memory limit would reduce its request below that setting.

### SLURM (optional)

To submit each executable rule as a separate SLURM job, install the optional
executor in the existing environment:

```bash
conda install -c conda-forge -c bioconda snakemake-executor-plugin-slurm
```

Then run from this directory, choosing a suitable maximum number of jobs:

```bash
snakemake --configfile config.yaml --profile profiles/slurm --jobs 10 --dry-run
snakemake --configfile config.yaml --profile profiles/slurm --jobs 10
```

The [SLURM profile](profiles/slurm/config.yaml) uses Snakemake's default
memory and disk estimates. Configure any required account, partition, time
limit and rule-specific resources for your cluster before submitting. The
Kraken2 rule explicitly requests `mem_mb=256000` per sample by default; the
SLURM executor translates this to a 250 GiB memory request. SLURM can schedule
multiple samples concurrently when the cluster has sufficient resources.
The Snakemake controller and compute jobs must be able to access the workflow,
inputs, outputs and active environment on a shared filesystem. Without the
profile, the local commands above are unchanged.

Paths in the configuration are relative to the working directory. `samples`
is a text file with one unique sample ID per line; each ID must correspond to
`fastq_dir/<ID>.fastq.gz`. IDs may contain letters, numbers, periods,
underscores and hyphens, and must start with a letter or number.
Set `kraken2_db` to a Kraken2 database containing `hash.k2d`, `opts.k2d`
and `taxo.k2d`. The configured database must be accessible to compute nodes.

The workflow writes `reads.tsv`, `reads_length_distribution.tsv`,
`stacked_barplot.png` and `stacked_barplot_log.png` under `output_dir`, plus
each sample's table at `output_dir/raw_files/<sample>/<sample>.reads.tsv`.
Kraken2 writes `<sample>.kraken.out`, `<sample>.kreport` and
`<sample>.kraken.log` in the same sample directory. It uses minimum hit groups
6, confidence 0.05 and minimizer data in its report, as in the original step.
Counts are calculated independently, so multiple samples can run concurrently.
Snakemake reruns affected samples when their FASTQ files change and rebuilds
the summary tables when the sample list changes. Use a new output directory
instead of pointing at existing results that should be preserved.

Rules are named `step_001_count_reads`, `step_002_summarize_reads`,
`step_003_plot_read_lengths` and `step_004_kraken2`. The numbers identify
steps; Kraken2 does not depend on the read-counting steps and can run in
parallel with them.

Resource benchmarks are written to `output_dir/benchmark/`: one report per
sample in `001_count_reads/` and `004_kraken2/`, plus
`002_summarize_reads.tsv` and `003_plot_read_lengths.tsv`. Existing benchmark
files in `benchmarks/` are not moved. The reports include elapsed time, peak
RAM (`max_rss`, in MiB) and CPU time (`cpu_time`, in seconds). Measurements are sampled, so short
jobs may report `NA` if they finish before a sample is collected. Existing
up-to-date outputs are not rerun just to create benchmark files. To benchmark
them, force the rules to rerun:

```bash
snakemake --configfile config.yaml --cores 24 --forcerun step_001_count_reads step_002_summarize_reads step_003_plot_read_lengths
```

This recalculates the outputs and plots.
