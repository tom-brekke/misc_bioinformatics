# misc_bioinformatics

A collection of command-line utilities for FASTA/FASTQ files, GC analysis,
genomic regions, transcript annotations, and SAM/BAM processing. These scripts
are independent tools rather than a single pipeline.

## Requirements

- Python 3 for the Python scripts.
- [Biopython](https://biopython.org/) for scripts that read FASTA records:
  `CalcGC.py`, `CalcGC_v1.py`, `CalcGC_v2.py`, `Reformat_fasta.py`,
  `extract_chr_lengths_simple.py`, `extract_fasta_region.py`, and
  `extract_fasta_region_by_sam.py`.
- Bash and `samtools` for the SAM-to-BAM scripts. They also call
  `run_cmd.sh`, which should be available on your `PATH` or in the working
  directory.

Run Python scripts with `python3 script.py --help` where supported. Several
older utilities use positional arguments instead of an option parser.

## Sequence conversion and formatting

| Script | What it does | Example use |
| --- | --- | --- |
| `Convert_Fastq_to_Fasta.py` | Reads a four-line-per-record FASTQ file and writes FASTA records to standard output, keeping each read name and sequence. | Convert reads before using a FASTA-only utility: `python3 Convert_Fastq_to_Fasta.py reads.fastq > reads.fasta`. |
| `Convert_Tab_to_Fasta.py` | Converts a tab-delimited FASTA representation to regular FASTA. It expects the name, description, and sequence in the first three columns. | `python3 Convert_Tab_to_Fasta.py sequences.tsv sequences.fasta`. |
| `Reformat_fasta.py` | Reads FASTA and can write each sequence on one line, or create simple BED and/or GTF files describing each whole record. | Make whole-record BED entries: `python3 Reformat_fasta.py -F genome.fasta -b True`. To write one-line FASTA, add `-f1 True -o genome_oneline.fasta`. The generated BED/GTF contain one interval per FASTA record, not feature annotations. |

## GC content and coding-sequence analysis

| Script | What it does | Example use |
| --- | --- | --- |
| `CalcGC.py` | Calculates GC percentage for each FASTA record. `-c` adds GC percentages by codon position; `-u` adds codon counts when used with `-c` and is intended for in-frame sequences. | Summarize ORFs: `python3 CalcGC.py -f orfs.fasta -c -o gc.tsv`. |
| `CalcGC_v1.py` | Extends per-record GC analysis with optional sliding windows (`-w`), target sequence names (`-t`), GFF-based frame information (`-g`), codon-position percentages (`-c`), and codon usage (`-u`, together with `-c`). The window option does not work with several annotation/codon options; see `--help`. | Calculate windowed GC across a genome: `python3 CalcGC_v1.py -f genome.fasta -w 10000 -o windows.tsv`. |
| `CalcGC_v2.py` | Uses genome FASTA and optional GFF annotations to extract and join CDS features by transcript, then outputs those sequences as FASTA or reports their GC statistics. `-t` selects the output type (`fasta` or `GC_table`); `-c` adds codon-position GC. | Produce a GC table from annotated CDSs: `python3 CalcGC_v2.py -f genome.fasta -g genes.gff3 -t GC_table -o cds_gc.tsv`. |

GC percentages are based on A/C/G/T bases; ambiguous bases such as `N` are
excluded from the denominator. Codon-position and codon-usage results assume
the sequence is in the correct coding frame.

## FASTA and genomic-region utilities

| Script | What it does | Example use |
| --- | --- | --- |
| `extract_chr_lengths_simple.py` | Prints a CSV-like list of FASTA record descriptions and sequence lengths to standard output. | `python3 extract_chr_lengths_simple.py genome.fasta > lengths.csv`. |
| `extract_fasta_region.py` | Extracts one or more regions from FASTA, optionally reverse-complementing them. Coordinates are 1-based; a chromosome list can contain chromosome names, coordinates, and strand. Output goes to standard output unless `-o` is given. | Extract bases 100–500 of `chr1`: `python3 extract_fasta_region.py -f genome.fasta -chr chr1 -s 100 -e 500 -o region.fasta`. Use `--help` for chromosome-list column formats. |
| `extract_fasta_region_by_sam.py` | Experimental utility intended to compare aligned reads and variants with a reference, writing `GOOD.fasta` and `ODD.fasta` (tab-separated read name/sequence records despite the suffixes) and a hard-coded VCF output file. | It expects genome FASTA, SAM, and VCF paths as positional arguments. The VCF-writing portion is incomplete, so inspect and test the script for your data before relying on its results. |

## Transcript annotation conversion

| Script | What it does | Example use |
| --- | --- | --- |
| `fasta_to_gft.py` | Converts Trinity transcript FASTA headers containing a `length=` value into a GTF-style, transcript/exon-like record. It takes input and output filenames as positional arguments. | `python3 fasta_to_gft.py transcripts.fasta transcripts.gft`. |
| `tran_to_gtf.py` | Writes one transcript-level GTF record per matching Trinity FASTA header (`comp..._c..._seq... len=...`). | `python3 tran_to_gtf.py transcripts.fasta transcripts.gtf`. |

These annotation converters rely on the header formats described above; they
do not infer transcript structure from sequence data.

## SAM/BAM utilities

| Script | What it does | Example use |
| --- | --- | --- |
| `Convert_sam_to_sorted_bam.sh` | Uses `samtools` to filter a SAM by minimum mapping quality, convert and sort it as BAM, and create an index. An optional second argument sets the minimum mapping quality (default `0`). | `bash Convert_sam_to_sorted_bam.sh reads.sam 20`. |
| `Convert_sam_to_sorted_bam_PrimaryOnly.sh` | Performs the same conversion while also excluding secondary alignments (`samtools` flag `256`). | `bash Convert_sam_to_sorted_bam_PrimaryOnly.sh reads.sam 20`. |
| `run_cmd.sh` | Helper used by the conversion scripts: runs a command only when the named output file does not already exist. | Usually called by the SAM/BAM scripts; it accepts a command string and an output path. |

**Important:** Both SAM-to-BAM scripts delete the input SAM and intermediate
unsorted BAM after processing. Keep a backup if you need to retain the SAM.
