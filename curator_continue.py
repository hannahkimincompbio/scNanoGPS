
#!/usr/bin/env python3

import os
import sys
import glob
import time
import subprocess
import multiprocessing as mp
from functools import partial

from curator_core import curator_io, preprocessing


def is_complete(bam_path, samtools):
    """Return True if the expected final BAM exists and passes quickcheck."""
    if not os.path.isfile(bam_path) or os.path.getsize(bam_path) == 0:
        return False

    result = subprocess.run(
        [samtools, "quickcheck", bam_path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return result.returncode == 0


def main():
    parser = preprocessing.getOptions()
    options, arguments = parser.parse_args()

    # Keep the continuation log separate from the original log.
    if options.log_f_name == "logs/curator.log.txt":
        options.log_f_name = "logs/curator_continue.log.txt"

    log_parent = os.path.dirname(os.path.abspath(options.log_f_name))
    os.makedirs(log_parent, exist_ok=True)

    # Validate inputs and normalize options as in the original Curator.
    parser, options = preprocessing.precheck(parser, options)

    if options.skip_curation:
        sys.exit(
            "ERROR: --skip_curation is enabled. "
            "This script expects curated.minimap2.bam outputs."
        )

    # Use a normalized path because curation_master derives the barcode
    # name from the FASTQ path and options.tmp_dir.
    options.tmp_dir = os.path.abspath(options.tmp_dir).rstrip(os.sep)

    # Match the thread-related environment settings in curator.py.
    for variable in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[variable] = str(options.ncores)

    fq_list = sorted(
        glob.glob(os.path.join(options.tmp_dir, "*.fastq.gz"))
    )

    if not fq_list:
        sys.exit(f"ERROR: No per-cell FASTQ files found in {options.tmp_dir}")

    pending = []

    print(f"Found {len(fq_list)} per-cell FASTQ files.", flush=True)

    for fq in fq_list:
        barcode = os.path.basename(fq)[:-len(".fastq.gz")]
        bam_path = os.path.join(
            options.tmp_dir, barcode + ".curated.minimap2.bam"
        )

        if is_complete(bam_path, options.samtools):
            print(f"SKIP: {barcode} (completed BAM found)", flush=True)
        else:
            print(f"QUEUE: {barcode} (missing/incomplete output)", flush=True)

            # The original Curator appends to this file. Remove it before
            # rerunning to avoid duplicate records from an interrupted run.
            stale_curated_sam = os.path.join(
                options.tmp_dir, barcode + ".curated.sam"
            )
            if os.path.exists(stale_curated_sam):
                os.remove(stale_curated_sam)

            # Remove a possibly incomplete final BAM before retrying.
            for suffix in (
                ".curated.minimap2.bam",
                ".curated.minimap2.bam.bai",
            ):
                stale = os.path.join(options.tmp_dir, barcode + suffix)
                if os.path.exists(stale):
                    os.remove(stale)

            pending.append(fq)

    print(
        f"\nCompleted: {len(fq_list) - len(pending)}"
        f"\nRemaining: {len(pending)}"
        f"\nWorkers:   {options.ncores}\n",
        flush=True,
    )

    if not pending:
        print("Nothing to resume; all final BAMs passed quickcheck.")
        return

    start_time = time.time()

    with open(options.log_f_name, "at") as logger:
        logger.write("\n=== curator_continue.py ===\n")
        logger.write(f"Time: {time.ctime()}\n")
        logger.write(f"FASTQ files to resume: {len(pending)}\n")

        with mp.Pool(processes=options.ncores) as pool:
            results = pool.map(
                partial(curator_io.curation_master, options=options),
                pending,
            )

        logger.write("\n".join(results) + "\n")
        logger.write(
            f"Continuation elapsed seconds: {time.time() - start_time:.1f}\n"
        )

    print("\nContinuation run finished.", flush=True)
    print(
        "Check the SLURM log and curator continuation log for errors.",
        flush=True,
    )


if __name__ == "__main__":
    main()

