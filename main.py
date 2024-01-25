"""This script checks and marks the duplicated reads in RNA-seq."""

import argparse
import logging
import numpy as np

import hailtop.batch as hb
import hailtop.fs as hfs

from utils import (DATA_PATHS_TABLE_ID,
                   DATA_PATHS_VIEW_ID,
                   DOWNSTREAM_ANALYSIS_TABLE_ID,
                   DOWNSTREAM_ANALYSIS_VIEW_ID,
                   RNA_SEQ_BASE_ID,
                   read_from_airtable,
                   switch_account_for_gtex)

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

REGION = ["us-central1"]
"""
GCP region.
"""
DEFAULT_MEMORY = "standard"
"""
Hail Batch default memory.
"""
DEFAULT_CPU = 2 ** 4
"""
Hail Batch default CPU.
"""
DEFAULT_STORAGE = "50G"
"""
Hail Batch default storage.
"""
DOCKER_IMAGE = "gcr.io/cmg-analysis/qc@sha256:612a6432e6365c7ea2156bab03d7e434bec9e5e7642399dc2270381c202cad55"
"""
Docker image.
"""
SAMPLE_ID_COL = "sample_id"
"""
Name of the column that has sample IDs.
"""
TISSUE_COL = "tissue"
"""
Name of the column that has tissue types.
"""
BAM_COL = "bam_path"
"""
Name of the column that has bam paths.
"""


def create_symbolic_links(batch, job, path, link_path, is_gtex):
    """
    Copy a file from GCP to the batch container.

    :param batch: Hail Batch.
    :param job: Hail Batch job.
    :param path: The GCP path to the file.
    :param link_path: The local path in the batch container.
    :param is_gtex: If True, run with GTEx samples; else, run with only RDG samples.
    """
    if is_gtex:
        job.command(f"gsutil -u {args.requester_pays_project} -m cp {path} {link_path}")
    else:
        localized_path = batch.read_input(path)
        job.command(f"ln -s {localized_path} {link_path}")


def get_mark_duplicates(job, sample_id, bam_file, is_gtex):
    """
    Mark duplicated reads in RNA-seq data using Picard.

    :param job: Hail Batch job.
    :param sample_id: ID of the current sample.
    :param bam_file: Local file path to the bam file of the current sample.
    :param is_gtex: If True, run with GTEx samples; else, run with only RDG samples.
    :return: Generate a txt file containing marked reads and copy the file to the GCP output path.
    """
    if is_gtex:
        switch_account_for_gtex(job)

    job.command("cd /io")
    job.command("ls -lh .")
    create_symbolic_links(batch, job, bam_file, f"{sample_id}.bam", is_gtex)
    job.command(f"java -jar /base/usr/picard/picard.jar MarkDuplicates "
                f"I={sample_id}.bam "
                f"O={sample_id}.mark_dup.bam "
                f"M=marked_dup_metrics.txt")
    job.command("ls -lh .")
    job.command(f"cp marked_dup_metrics.txt {job.ofile} ")
    batch.write_output(job.ofile, f"{args.out_dir}/{sample_id}_marked_dup_metrics.txt")


def get_fasqc(job, sample_id, bam_file, is_gtex):
    """
    Generate the FASTQC QC report.

    :param job: Hail Batch job.
    :param sample_id: ID of the current sample.
    :param bam_file: Local file path to the bam file of the current sample.
    :param is_gtex: If True, run with GTEx samples; else, run with only RDG samples.
    :return: Generate a HTML file and copy the file to the GCP output path.
    """
    if is_gtex:
        switch_account_for_gtex(job)

    job.command("cd /io")
    create_symbolic_links(batch, job, bam_file, f"{sample_id}.bam", is_gtex)
    job.command(f"fastqc {sample_id}.bam")
    job.command("ls -lh .")
    job.command(f"cp *.html {job.ofile} ")
    batch.write_output(job.ofile, f"{args.out_dir}/{sample_id}.html")


def get_read_count(job, sample_id, bam_file, is_gtex):
    """
    Generate read counts summary statistics using samtools.

    :param job: Hail Batch job.
    :param sample_id: ID of the current sample.
    :param bam_file: Local file path to the bam file of the current sample.
    :param is_gtex: If True, run with GTEx samples; else, run with only RDG samples.
    :return: Generate a read count txt file and copy the file to the GCP output path.
    """
    if is_gtex:
        switch_account_for_gtex(job)

    job.command("cd /io")
    create_symbolic_links(batch, job, bam_file, f"{sample_id}.bam", is_gtex)
    job.command(f"samtools view -c {sample_id}.bam > {job.ofile}")
    job.command(f"samtools view -c -F 260 {sample_id}.bam >> {job.ofile}")
    batch.write_output(job.ofile, f"{args.out_dir}/{sample_id}_read_counts.txt")


def run_samples(batch, id_list, bam_list, is_gtex):
    """
    Run the pipeline on all samples.

    :param batch: Hail Batch.
    :param id_list: A list containing sample IDs.
    :param bam_list: A list containing GCP paths to bam files of samples.
    :param is_gtex: If True, run with GTEx samples; else, run with only RDG samples.
    :return: Three files: marked duplicated reads, QC report, and read counts.
    """
    for i in range(len(id_list)):
        cur_id = id_list[i]
        cur_bam = bam_list[i]

        cur_fastqc_job = batch.new_job(f"{cur_id}_fastqc")
        cur_md_job = batch.new_job(f"{cur_id}_md")
        cur_read_count_job = batch.new_job(f"{cur_id}_read_count")

        cur_read_count_job.storage(hfs.stat(cur_bam).size)
        cur_read_count_job.cpu(1)

        get_fasqc(cur_fastqc_job, cur_id, cur_bam, is_gtex)
        get_mark_duplicates(cur_md_job, cur_id, cur_bam, is_gtex)
        get_read_count(cur_read_count_job, cur_id, cur_bam, is_gtex)


def main(args, batch):
    rd_samples = read_from_airtable(RNA_SEQ_BASE_ID,
                                    DATA_PATHS_TABLE_ID,
                                    DATA_PATHS_VIEW_ID)

    if args.s:
        ids_set = set(args.s)
        subset_samples = rd_samples[rd_samples[SAMPLE_ID_COL].isin(ids_set)]
        rd_sample_ids = np.array(subset_samples[SAMPLE_ID_COL])
        rd_bam_files = np.array(subset_samples[BAM_COL])

    else:
        new_rd_samples = rd_samples

        new_muscle_rd_samples = new_rd_samples[new_rd_samples[TISSUE_COL] == args.tissue]

        gtex_samples = read_from_airtable(RNA_SEQ_BASE_ID,
                                          DOWNSTREAM_ANALYSIS_TABLE_ID,
                                          DOWNSTREAM_ANALYSIS_VIEW_ID)
        gtex_samples = gtex_samples[gtex_samples[SAMPLE_ID_COL].str.startswith("GTEX")]
        gtex_muscle_samples = gtex_samples[gtex_samples[TISSUE_COL] == args.tissue]

        rd_sample_ids = np.array(new_muscle_rd_samples[SAMPLE_ID_COL])
        rd_bam_files = np.array(new_muscle_rd_samples[BAM_COL])

        gtex_sample_ids = np.array(gtex_muscle_samples[SAMPLE_ID_COL])
        gtex_bam_files = np.array(gtex_muscle_samples[BAM_COL])
        run_samples(batch, gtex_sample_ids, gtex_bam_files, True)

    run_samples(batch, rd_sample_ids, rd_bam_files, False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--billing-project",
        type=str,
        help="Project to bill under.",
    )
    parser.add_argument(
        "--requester-pays-project",
        type=str,
        help="Requester pays project to bill under.",
    )
    parser.add_argument(
        "--file-dir",
        type=str,
        help="The temporary directory for Hail",
    )
    parser.add_argument(
        "--tissue",
        type=str,
        help="Tissue type of the samples.",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        help="The GCP output directory.",
    )
    parser.add_argument(
        "-s",
        nargs="+",
        help="IDs of specific samples to process (non-GTEx for now).",
    )
    args = parser.parse_args()

    backend = hb.ServiceBackend(
        billing_project=args.billing_project,
        remote_tmpdir=args.file_dir,
        regions=REGION,
    )
    batch = hb.Batch(
        backend=backend,
        requester_pays_project=args.requester_pays_project,
        default_image=DOCKER_IMAGE,
        default_memory=DEFAULT_MEMORY,
        default_cpu=DEFAULT_CPU,
        default_storage=DEFAULT_STORAGE
    )

    main(args, batch)
    batch.run()
