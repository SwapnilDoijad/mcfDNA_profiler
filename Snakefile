import os
import re


for key in ("samples", "fastq_dir", "output_dir", "kraken2_db"):
    if not isinstance(config.get(key), str) or not config[key]:
        raise ValueError(f"Configure a non-empty '{key}' in the config file")

if type(config.get("kraken2_mem_mb")) is not int or config["kraken2_mem_mb"] <= 0:
    raise ValueError("Configure a positive integer 'kraken2_mem_mb'")

with open(config["samples"]) as sample_file:
    SAMPLES = [line.strip() for line in sample_file if line.strip()]

if not SAMPLES:
    raise ValueError("The sample list is empty")
if len(SAMPLES) != len(set(SAMPLES)):
    raise ValueError("The sample list contains duplicate IDs")
if any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample) for sample in SAMPLES):
    raise ValueError("Sample IDs must contain only letters, numbers, '.', '_' or '-'")

OUTPUT_DIR = config["output_dir"]
SAMPLE_STATS = os.path.join(OUTPUT_DIR, "raw_files", "{sample}", "{sample}.reads.tsv")
READS = os.path.join(OUTPUT_DIR, "reads.tsv")
LENGTHS = os.path.join(OUTPUT_DIR, "reads_length_distribution.tsv")
PLOT = os.path.join(OUTPUT_DIR, "stacked_barplot.png")
PLOT_LOG = os.path.join(OUTPUT_DIR, "stacked_barplot_log.png")
KRAKEN_OUTPUT = os.path.join(OUTPUT_DIR, "raw_files", "{sample}", "{sample}.kraken.out")
KRAKEN_REPORT = os.path.join(OUTPUT_DIR, "raw_files", "{sample}", "{sample}.kreport")
BENCHMARK_DIR = os.path.join(OUTPUT_DIR, "benchmark")
COUNT_BENCHMARK = os.path.join(BENCHMARK_DIR, "001_count_reads", "{sample}.tsv")
KRAKEN_BENCHMARK = os.path.join(BENCHMARK_DIR, "004_kraken2", "{sample}.tsv")
BENCHMARK_SUMMARY = os.path.join(OUTPUT_DIR, "benchmarks.tsv")

rule all:
    input:
        READS,
        LENGTHS,
        PLOT,
        PLOT_LOG,
        expand(KRAKEN_OUTPUT, sample=SAMPLES),
        expand(KRAKEN_REPORT, sample=SAMPLES),


rule step_001_count_reads:
    input:
        os.path.join(config["fastq_dir"], "{sample}.fastq.gz"),
    output:
        SAMPLE_STATS,
    benchmark:
        os.path.join(OUTPUT_DIR, "benchmark", "001_count_reads", "{sample}.tsv")
    threads: 1
    shell:
        r"""
        seqkit fx2tab -j {threads} -n -l {input:q} |
        awk -F '\t' -v sample={wildcards.sample:q} '
            {{ n++; a += ($NF >= 35); b += ($NF >= 50); c += ($NF >= 100); d += ($NF >= 250) }}
            END {{ printf "%s\t%d\t%d\t%d\t%d\t%d\n", sample, n, a, b, c, d }}
        ' > {output:q}
        """

rule step_002_summarize_reads:
    input:
        samples=config["samples"],
        counts=expand(SAMPLE_STATS, sample=SAMPLES),
    output:
        reads=READS,
        lengths=LENGTHS,
    benchmark:
        os.path.join(OUTPUT_DIR, "benchmark", "002_summarize_reads.tsv")
    threads: 1
    shell:
        r"""
        python - {output.reads:q} {output.lengths:q} {input.counts:q} <<-'PY'
	import sys

	reads_path, lengths_path, *count_paths = sys.argv[1:]
	with open(reads_path, "w") as reads, open(lengths_path, "w") as lengths:
	    reads.write("ID\tReads\n")
	    lengths.write("sample\ttotal\tge35\tge50\tge100\tge250\n")
	    for path in count_paths:
	        with open(path) as stats:
	            row = stats.read().strip()
	        sample, total, _, _, _, _ = row.split("\t")
	        reads.write(sample + "\t" + total + "\n")
	        lengths.write(row + "\n")
	PY
        """


rule step_003_plot_read_lengths:
    input:
        lengths=LENGTHS,
        script=os.path.join(workflow.basedir, "scripts", "stacked_barplot.py"),
    output:
        PLOT,
        PLOT_LOG,
    benchmark:
        os.path.join(OUTPUT_DIR, "benchmark", "003_plot_read_lengths.tsv")
    params:
        output_dir=OUTPUT_DIR,
    threads: 1
    shell:
        "MPLBACKEND=Agg python {input.script:q} -i {input.lengths:q} -o {params.output_dir:q}"


rule step_004_kraken2:
    input:
        fastq=os.path.join(config["fastq_dir"], "{sample}.fastq.gz"),
        database=expand(
            os.path.join(config["kraken2_db"], "{file}.k2d"),
            file=["hash", "opts", "taxo"],
        ),
    output:
        classified=KRAKEN_OUTPUT,
        report=KRAKEN_REPORT,
    params:
        database=config["kraken2_db"],
        required_mem=config["kraken2_mem_mb"],
    log:
        os.path.join(OUTPUT_DIR, "raw_files", "{sample}", "{sample}.kraken.log"),
    benchmark:
        os.path.join(OUTPUT_DIR, "benchmark", "004_kraken2", "{sample}.tsv")
    threads: 24
    resources:
        mem_mb=config["kraken2_mem_mb"],
    shell:
        r"""
        if [ {resources.mem_mb} -lt {params.required_mem} ]; then
            echo "Kraken2 needs mem_mb={params.required_mem}; the global resource limit reduced it to {resources.mem_mb}." | tee {log:q} >&2
            exit 1
        fi
        kraken2 \
            --db {params.database:q} \
            --report-minimizer-data \
            --minimum-hit-groups 6 \
            --threads {threads} \
            --output {output.classified:q} \
            --report {output.report:q} \
            --confidence 0.05 \
            {input.fastq:q} 2> {log:q}
        """

rule step_100_collect_benchmarks:
    input:
        counts=expand(COUNT_BENCHMARK, sample=SAMPLES),
        kraken=expand(KRAKEN_BENCHMARK, sample=SAMPLES),
    output:
        BENCHMARK_SUMMARY,
    run:
        import csv

        sources = (
            ("step_001_count_reads", input.counts),
            ("step_004_kraken2", input.kraken),
        )

        with open(output[0], "w", newline="") as summary:
            writer = csv.writer(summary, delimiter="\t", lineterminator="\n")
            benchmark_header = None

            for rule_name, paths in sources:
                for sample_id, path in zip(SAMPLES, paths):
                    with open(path, newline="") as benchmark:
                        reader = csv.reader(benchmark, delimiter="\t")
                        header = next(reader, None)
                        if not header:
                            raise ValueError(f"Empty benchmark file: {path}")
                        if benchmark_header is None:
                            benchmark_header = header
                            writer.writerow(["sample_id", "rule", *header])
                        elif header != benchmark_header:
                            raise ValueError(f"Benchmark columns differ in {path}")

                        for row in reader:
                            writer.writerow([sample_id, rule_name, *row])

