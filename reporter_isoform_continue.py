#!/usr/bin/env python3
"""Resume scNanoGPS reporter_isoform at isoform-matrix generation.

This script does not rerun IsoQuant or rebuild the merged BAM. It reads an
existing IsoQuant grouped transcript-count table and writes the scNanoGPS
isoform expression matrix.
"""
import time, os, sys, glob, re
import pandas as pd
import multiprocessing as mp
from functools import partial
from contextlib import contextmanager
from optparse import OptionParser
from scanner_core import misc
from curator_core import curator_io


def parse_txid2txname(gtf_path):
    """Map transcript_id to transcript_name from transcript records in a GTF."""
    tx_mapping = {}
    with open(gtf_path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = dict(re.findall(r'([^\s;]+)\s+"([^"]*)"', fields[8]))
            txid = attrs.get("transcript_id")
            if txid:
                tx_mapping[txid] = attrs.get("transcript_name", txid)
    return tx_mapping


def main():
    parser = OptionParser(
        description="Resume reporter_isoform at matrix generation; does not rerun IsoQuant."
    )

    parser.add_option("-d",           dest = "o_dir", default = "scNanoGPS_res",
                      help = "Output directory name. "
                             "Default: scNanoGPS_res")
    parser.add_option("--tmp_dir",    dest = "tmp_dir",  default = "tmp",
                      help = "Temporary folder name. "
                             "Default: tmp")
    parser.add_option("--CB_file",    dest = "CB_file",  default = "filtered_barcode_list.txt",
                      help = "File name for filtered barcode list. "
                             "Default: filtered_barcode_list.txt")
    parser.add_option("--ref_genome", dest = "ref_genome",  default = None,
                      help = "* Required ! "
                             "File for reference genome.")
    parser.add_option("--gtf",        dest = "gtf",      default = None,
                      help = "* Required ! "
                             "File for genome annotation.")
    parser.add_option("--bam_suf",    dest = "bam_suf",  default = ".curated.minimap2.bam",
                      help = "Suffix of the bam files. "
                             "Default: .curated.minimap2.bam")
    parser.add_option("-o",           dest = "o_name",      default = "isoform_exp_matrix.tsv.gz",
                      help = "Counting table name. "
                             "Default: isoform_exp_matrix.tsv.gz")
    parser.add_option("--log",        dest = "log_f_name",  default = 'logs/reporter_isoform.log.txt',
                      help = "Log file name. "
                             "Default: logs/reporter_isoform_continue.log.txt")
    parser.add_option("-t",           dest = "ncores",   default = 1,
                      help = "Number of cores for program running. "
                             "Default: 1", type = "int")
    parser.add_option("--batch_no",   dest = "batch_no",  default = 500,
                      help = "Batch number for merging bam files. "
                             "Default: 500", type = "int")
    parser.add_option("--samtools",   dest = "samtools",  default = "samtools",
                      help = "Path to samtools. "
                             "Default: samtools")
    parser.add_option("--isoquant",   dest = "isoquant",  default = "isoquant",
                      help = "Path to IsoQuant. "
                             "Default: isoquant")
    parser.add_option("--isoquant_d", dest = "isoquant_d", default = "scNanoGPS_res/IsoQuant_res",
                      help = "IsoQuant output directory. "
                             "Default: scNanoGPS_res/IsoQuant_res")
    parser.add_option("--isoquant_o", dest = "isoquant_o",  default = 'OUT/OUT.transcript_grouped_tag_RG_counts.tsv',
                      help = "Suffix of IsoQuant output file; \
                        was OUT/OUT.transcript_model_grouped_counts.tsv in isoQua)nt v3.7.0, \
                        now OUT/OUT.transcript_grouped_tag_RG_counts.tsv in isoQuant v4.0.0 "
                             "Default: OUT/OUT.transcript_grouped_tag_RG_counts.tsv")
    options, args = parser.parse_args()

    start_time = time.time()
    output_dir = os.path.abspath(options.o_dir)
    gtf_path = os.path.abspath(options.gtf)
    isoquant_dir = os.path.abspath(options.isoquant_d)
    input_path = os.path.join(isoquant_dir, options.isoquant_o)
    output_path = options.o_name if os.path.isabs(options.o_name) else os.path.join(output_dir, options.o_name)

    if not os.path.isdir(output_dir):
        parser.error(f"Output directory does not exist: {output_dir}")
    if not os.path.isfile(gtf_path):
        parser.error(f"GTF file does not exist: {gtf_path}")
    if not os.path.isfile(input_path):
        parser.error(f"IsoQuant count table does not exist: {input_path}")

    if options.log_f_name:
        log_path = os.path.abspath(options.log_f_name)
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        with open(log_path, "wt") as log:
            log.write("Command: " + " ".join(sys.argv) + "\n")
            log.write("Input: " + input_path + "\n")
            log.write("GTF: " + gtf_path + "\n")
            log.write("Output: " + output_path + "\n")

    print(f"Reading IsoQuant grouped counts: {input_path}", flush=True)
    res_df = pd.read_csv(input_path, header=0, sep="\t", compression="infer")
    if len(res_df.columns) == 0:
        raise RuntimeError(f"Input table has no columns: {input_path}")

    # scNanoGPS v2.0 expected '#feature_id'; accept IsoQuant 4.x 'feature_id' too.
    if "#feature_id" in res_df.columns:
        res_df = res_df.rename(columns={"#feature_id": "transcript_id"})
    elif "feature_id" in res_df.columns:
        res_df = res_df.rename(columns={"feature_id": "transcript_id"})
    elif "gene_id" in res_df.columns: ## Fixed: "gene_id" in IsoQuant 4.0.0 OUT.transcript_grouped_tag_RG_counts.tsv
        res_df = res_df.rename(columns={"gene_id": "transcript_id"})
    elif "transcript_id" not in res_df.columns:
        raise RuntimeError(
            "Expected a '#feature_id', 'feature_id', 'gene_id' or 'transcript_id' column. "
            f"Found columns: {res_df.columns.tolist()}"
        )

    # Do not silently treat IsoQuant's ungrouped feature_id/count output as a cell matrix.
    if len(res_df.columns) == 2 and str(res_df.columns[1]).strip().lower() == "count":
        raise RuntimeError(
            "Input appears to be an ungrouped feature_id/count table, not a "
            "feature-by-read-group/barcode matrix. Use the appropriate grouped counts output."
        )

    print("Parsing transcript names from GTF...", flush=True)
    tx_mapping = parse_txid2txname(gtf_path)
    if not tx_mapping:
        raise RuntimeError("No transcript_id values found on GTF transcript records.")

    res_df = res_df.loc[res_df["transcript_id"].isin(tx_mapping.keys())].copy()
    if res_df.empty:
        raise RuntimeError(
            "No input transcript IDs matched transcript IDs in the GTF. "
            "Check that the count table and GTF use matching IDs; no matrix was written."
        )

    res_df.insert(1, "transcript_name", res_df["transcript_id"].map(tx_mapping))
    count_columns = [c for c in res_df.columns if c not in ("transcript_id", "transcript_name")]
    res_df = res_df[["transcript_id", "transcript_name"] + count_columns]

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    compression = "gzip" if output_path.endswith(".gz") else None
    res_df.to_csv(output_path, sep="\t", header=True, index=False, compression=compression)

    print(f"Wrote matrix: {output_path}", flush=True)
    print(f"Matched transcripts: {len(res_df):,}", flush=True)
    print(f"Count/read-group columns: {len(count_columns):,}", flush=True)
    elapsed = time.time() - start_time
    print(f"Finished: {time.strftime('%a, %d %b %Y %H:%M:%S', time.localtime())}", flush=True)
    print(f"Elapsed: {elapsed / 3600:.0f} h {(elapsed % 3600) / 60:.0f} min {elapsed % 60:.1f} s", flush=True)


if __name__ == "__main__":
    main()
