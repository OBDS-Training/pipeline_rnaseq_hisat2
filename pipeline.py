###########
# Imports #
###########

from ruffus import collate, follows, merge, mkdir, regex, transform
import sys
import os
import re
from cgatcore import pipeline as P

#################
# Configuration #
#################


# Load parameters from config file, located in `./config/pipeline.yml`.
PARAMS = P.get_parameters(
    "%s/pipeline.yml" % os.path.dirname(os.path.realpath(__file__))
)


############
# Workflow #
############


@follows(mkdir("results/qc/fastqc"))
@transform(
    "data/fastq/*.fastq.gz",
    regex(r"data/fastq/(.*).fastq.gz"),
    r"results/qc/fastqc/\1_fastqc.html"
)
def run_fastqc_on_fastq(infile, outfile):
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
@merge(run_fastqc_on_fastq, "results/reports/multiqc/fastq.html")
def run_multiqc_for_fastq(infiles, outfile):
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
    "data/fastq/*.fastq.gz",
    regex(r"data/fastq/(.+)_R[12].fastq.gz"),
    r"results/hisat2/\1.bam"
)
def run_hisat2_on_fastq(input_files, output_file):
    # the arguments 'input_file' and 'output_files' are not used here
    # input and output files are generated from 'config/fastq_files.tsv'

    hisat2_threads = PARAMS["hisat2"]["threads"]
    hisat2_genome = PARAMS["hisat2"]["genome"]
    hisat2_options = PARAMS["hisat2"]["options"]

    fastqs1 = [file for file in input_files if bool(re.search("_R1.fastq.gz", file))]
    fastqs2 = [file for file in input_files if bool(re.search("_R2.fastq.gz", file))]

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
    run_hisat2_on_fastq,
    regex(r"results/hisat2/(.*).bam"),
    r"results/qc/samtools/idxstats/\1",
)
def run_idxstats_on_bam(infile, outfile):
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
    run_hisat2_on_fastq,
    regex(r"results/hisat2/(.*).bam"),
    r"results/qc/samtools/flagstat/\1",
)
def run_flagstat_on_bam(infile, outfile):
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
    run_hisat2_on_fastq,
    regex(r"results/hisat2/(.*).bam"),
    r"results/qc/picard/CollectAlignmentSummaryMetrics/\1",
)
def run_picard_alignment_metrics_on_bam(infile, outfile):
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
    run_hisat2_on_fastq,
    regex(r"results/hisat2/(.*).bam"),
    r"results/qc/picard/CollectInsertSizeMetrics/\1",
)
def run_picard_insert_size_on_bam(infile, outfile):
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
@merge(run_hisat2_on_fastq, "results/featureCounts/counts")
def run_featurecounts_on_bam(infiles, outfile):
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
        job_threads=PARAMS["feature_counts"]["threads"],
    )


@merge(
    [
        run_idxstats_on_bam,
        run_flagstat_on_bam,
        run_picard_alignment_metrics_on_bam,
        run_picard_insert_size_on_bam,
        run_featurecounts_on_bam,
    ],
    "results/reports/multiqc/bam.html",
)
def run_multiqc_for_bam(infiles, outfile):
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


@follows(run_multiqc_for_fastq, run_multiqc_for_bam)
def full():
    pass


##################
# Main execution #
##################

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
