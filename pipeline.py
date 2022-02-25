"""
===========================
Pipeline template
===========================

.. Replace the documentation below with your own description of the
   pipeline's purpose

Overview
========

This pipeline computes the word frequencies in the configuration
files :file:``config.yml` and :file:`conf.py`.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use cgat pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.yml` file.
cgatReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.yml` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_@template@.py config

Input files
-----------

None required except the pipeline configuration files.

Requirements
------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

Pipeline output
===============

.. Describe output files of the pipeline here

Glossary
========

.. glossary::


Code
====

"""

###########
# Imports #
###########

from ruffus import *
import sys
import os
import re
import pandas as pd
import cgatcore.experiment as E
from cgatcore import pipeline as P

#################
# Configuration #
#################


# Load parameters from config file, located in `./config/pipeline.yml`.
PARAMS = P.get_parameters(
    "%s/config/pipeline.yml" % os.path.dirname(os.path.realpath(__file__))
)


############
# Workflow #
############


@follows(mkdir("data/fastq"))
@subdivide(
    "config/fastq_links.tsv",
    formatter(),
    # Output parameter: Glob matches any number of output file names
    "data/fastq/*.fastq.gz"
)
def symlink_fastq(input_file, output_files):
    # the arguments 'input_file' and 'output_files' are not used here
    # input and output files are generated from 'config/fastq_files.tsv'

    fastq_links = pd.read_csv(input_file, sep="\t", names=["source", "link"])
    for i in fastq_links.index:
        source_file = os.path.abspath(fastq_links["source"][int(i)])
        destination_dir = os.path.dirname(os.path.realpath(__file__))
        destination_link = os.path.join(destination_dir, "data", "fastq", fastq_links["link"][int(i)])
        os.symlink(source_file, destination_link)


@follows(mkdir("results/qc/fastqc"))
@transform(
    symlink_fastq, regex(r"data/fastq/(.*).fastq.gz"), r"results/qc/fastqc/\1_fastqc.html"
)
def fastqc_on_fastq(infile, outfile):
    """
    Run `fastqc` on the FASTQ files.
    """

    statement = """
        fastqc 
            -o results/qc/fastqc
            --nogroup
            %(infile)s
            > %(outfile)s.log
            2>&1
    """

    P.run(statement, job_condaenv="pipeline_rnaseq_hisat2")


@follows(mkdir("results/reports/multiqc"))
@merge(fastqc_on_fastq, "results/reports/multiqc/fastq.html")
def multiqc_fastq(infiles, outfile):
    """
    Run `multiqc` on the files of quality metrics computed for the FASTQ files.
    """

    statement = """
        multiqc 
            -n fastq.html
            -o results/reports/multiqc
            results/qc/fastqc
            > %(outfile)s.log
            2>&1
    """

    P.run(statement, job_condaenv="pipeline_rnaseq_hisat2")


@follows(mkdir("results/hisat2"))
@collate(
    symlink_fastq,
    regex(r"data/fastq/(.+)-.+-[12].fastq.gz"),
    r"results/hisat2/\1.bam"
)
def hisat2_on_fastq(input_files, output_file):
    # the arguments 'input_file' and 'output_files' are not used here
    # input and output files are generated from 'config/fastq_files.tsv'

    hisat2_threads = PARAMS["hisat2"]["threads"]
    hisat2_genome = PARAMS["hisat2"]["genome"]
    hisat2_options = PARAMS["hisat2"]["options"]

    fastqs1 = [file for file in input_files if bool(re.search("-1.fastq.gz", file))]
    fastqs2 = [file for file in input_files if bool(re.search("-2.fastq.gz", file))]

    fastqs1 = ",".join(fastqs1)
    fastqs2 = ",".join(fastqs2)

    statement = """
        hisat2
            --threads %(hisat2_threads)s
            -x %(hisat2_genome)s
            -1 %(fastqs1)s
            -2 %(fastqs2)s
            %(hisat2_options)s
            --summary-file %(output_file)s.log
        | samtools sort
            -@ %(hisat2_threads)s
            -o %(output_file)s
            -
        && samtools index
            %(output_file)s
        """ % locals()

    P.run(
        statement,
        job_condaenv="pipeline_rnaseq_hisat2",
        job_threads=PARAMS["hisat2_threads"],
        job_memory=PARAMS["hisat2_memory_per_thread"]
    )


@follows(mkdir("results/qc/samtools/idxstats"))
@transform(
    hisat2_on_fastq,
    regex(r"results/hisat2/(.*).bam"),
    r"results/qc/samtools/idxstats/\1",
)
def idxstats_on_bam(infile, outfile):
    """
    Run `samtools idxstats` on the BAM files produced by HISAT2.
    """

    statement = """
        samtools idxstats
        %(infile)s
        > %(outfile)s
    """

    P.run(statement, job_condaenv="pipeline_rnaseq_hisat2")


@follows(mkdir("results/qc/samtools/flagstat"))
@transform(
    hisat2_on_fastq,
    regex(r"results/hisat2/(.*).bam"),
    r"results/qc/samtools/flagstat/\1",
)
def flagstat_on_bam(infile, outfile):
    """
    Run `samtools flagstat` on the BAM files produced by HISAT2.
    """

    statement = """
        samtools flagstat
        %(infile)s
        > %(outfile)s
    """

    P.run(statement, job_condaenv="pipeline_rnaseq_hisat2")


@follows(mkdir("results/qc/picard/CollectAlignmentSummaryMetrics"))
@transform(
    hisat2_on_fastq,
    regex(r"results/hisat2/(.*).bam"),
    r"results/qc/picard/CollectAlignmentSummaryMetrics/\1",
)
def picard_alignment_metrics_on_bam(infile, outfile):
    """
    Run `picard CollectAlignmentSummaryMetrics` on the BAM files produced by HISAT2.
    """

    statement = """
        picard CollectAlignmentSummaryMetrics
            -I %(infile)s
            -O %(outfile)s
            -R %(picard_genome)s
            2> %(outfile)s.log
        """

    P.run(statement, job_condaenv="pipeline_rnaseq_hisat2")


@follows(mkdir("results/qc/picard/CollectInsertSizeMetrics"))
@transform(
    hisat2_on_fastq,
    regex(r"results/hisat2/(.*).bam"),
    r"results/qc/picard/CollectInsertSizeMetrics/\1",
)
def picard_insert_size_on_bam(infile, outfile):
    """
    Run `picard CollectInsertSizeMetrics` on the BAM files produced by HISAT2.
    """

    statement = """
        picard CollectInsertSizeMetrics
        -I %(infile)s
        -O %(outfile)s
        -R %(picard_genome)s
        -H %(outfile)s.pdf
        2> %(outfile)s.log
    """

    P.run(statement, job_condaenv="pipeline_rnaseq_hisat2")


@follows(mkdir("results/featureCounts"))
@merge(hisat2_on_fastq, "results/featureCounts/counts")
def featurecounts_on_bam(infiles, outfile):
    """
    Run `featureCounts`
    """

    infile_list = " ".join(infiles)

    statement = """
        featureCounts
            -a %(feature_counts_gtf)s
            -o %(outfile)s
            -T %(feature_counts_threads)s
            %(feature_counts_options)s
            %(infile_list)s
            > %(outfile)s.log
            2>&1
        """

    P.run(
        statement,
        job_condaenv="pipeline_rnaseq_hisat2",
        job_threads=PARAMS["feature_counts_threads"],
    )


@merge(
    [
        idxstats_on_bam,
        flagstat_on_bam,
        picard_alignment_metrics_on_bam,
        picard_insert_size_on_bam,
        featurecounts_on_bam,
    ],
    "results/reports/multiqc/bam.html",
)
def multiqc_bam(infiles, outfile):
    """
    Run MultiQC on the files of quality metrics computed for the BAM files.
    """

    statement = """
        multiqc
            -n bam.html
            -o results/reports/multiqc
            results/hisat2
            results/qc/samtools/idxstats
            results/qc/samtools/flagstat
            results/qc/picard/CollectAlignmentSummaryMetrics
            results/qc/picard/CollectInsertSizeMetrics
            results/featureCounts
            > %(outfile)s.log
            2>&1
    """

    P.run(statement, job_condaenv="pipeline_rnaseq_hisat2")


@follows(multiqc_fastq, multiqc_bam, featurecounts_on_bam)
def full():
    pass


##################
# Main execution #
##################

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
