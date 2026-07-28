# Honeycomb - Graph Based Barcode Deconvolution

Barcode deconvolution and graph-based clustering pipeline for long-read single-cell RNA-seq data.

This repository extracts putative barcodes from FASTQ reads, matches them to a whitelist of anchor barcodes, builds a barcode graph, clusters the graph into communities, and writes per-community output tables.

## Overview

The main pipeline is driven by [main.py](main.py) and runs in four stages:

1. Barcode extraction from FASTQ reads.
2. Edge construction between putative barcodes and anchors.
3. Graph construction and Leiden clustering.
4. Read assignment to clustered communities.

## Installation

Create and activate a python virtual environment, then install the Python dependencies.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Currently, Python version 3.12 is recommended.

Finally clone the repository and it is ready to run.

```bash
git clone https://github.com/dwschaeper/LR-scRNA-seq-GBBD
```

## Input files

The pipeline expects:

- A FASTQ file containing long-read single-cell reads.
- A whitelist file containing one barcode per line.

## Run the pipeline

Example invocation:

```bash
python main.py \
  --cell_barcode_whitelist /path/to/whitelist \
  --input /path/to/reads.fastq.gz \
  --processes 10 \
  --expected_cells 1000 \
  --outdir /path/to/output
  --prefix [output_file_prefix]
```

Use `python main.py --help` to see the full list of options.

## Key outputs

The pipeline writes intermediate and final outputs into the configured output directory and prefix.

Typical outputs include:

- `*_anchor_probabilities.csv` - negative-binomial anchor probabilities.
- `*_community_sizes.csv` - per-community size summary.
- `*_cell_membership.tsv` - read-to-community membership table.
- `*_cell_sizes.tsv` - per-community read counts.
- `*_assigned_reads.fastq.gz` - assigned read shard merge.

## Repository layout

- [main.py](main.py) - command-line entry point.
- [config.py](config.py) - shared runtime configuration for output prefix and directory.
- [src/extract_barcodes.py](src/extract_barcodes.py) - barcode extraction and anchor estimation.
- [src/edge_construction.py](src/edge_construction.py) - barcode matching and edge construction.
- [src/make_graph.py](src/make_graph.py) - graph assembly and clustering.
- [src/assignment.py](src/assignment.py) - final read assignment to communities.
- [src/logging_config.py](src/logging_config.py) - logging helpers.

## Notes

- Output paths are controlled by the `--outdir` and `--prefix` arguments passed to [main.py](main.py).
- The pipeline uses multiprocessing, so large inputs benefit from tuning `--processes` and `--batch_size`.
