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


# Load parameters from the configuration file 'pipeline.yml'
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
def run_fastqc_on_fastq(input_file, output_file):
    """
    Run `fastqc` on the FASTQ files.
    """
    # Prepare a command from the input and output file names.
    statement = """
        fastqc
            -o results/qc/fastqc
            --nogroup
            %(input_file)s
            > %(output_file)s.log
            2>&1
    """
    # Run the command in the Conda environment.
    P.run(statement, job_condaenv="pipeline_rnaseq_hisat2")


@follows(mkdir("results/reports/multiqc"))
@merge(
    run_fastqc_on_fastq,
    "results/reports/multiqc/fastq.html"
)
def run_multiqc_for_fastq(input_files, output_file):
    """
    Run `multiqc` on the files of quality metrics computed for the FASTQ files.
    """
    # Prepare a command from the input folder and output file names.
    statement = """
        multiqc
            -n fastq.html
            -o results/reports/multiqc
            results/qc/fastqc
            > %(output_file)s.log
            2>&1
    """
    # Run the command in the Conda environment.
    P.run(statement, job_condaenv="pipeline_rnaseq_hisat2")


@follows(mkdir("results/hisat2"))
@collate(
    "data/fastq/*.fastq.gz",
    regex(r"data/fastq/(.+)_R[12].fastq.gz"),
    r"results/hisat2/\1.bam"
)
def run_hisat2_on_fastq(input_files, output_file):
    # Fetch relevant pipeline parameters.
    hisat2_threads = PARAMS['hisat2']['threads']
    hisat2_genome = PARAMS['hisat2']['genome']
    hisat2_options = PARAMS['hisat2']['options']
    # Separate the input file names.
    fastq1 = input_files[0]
    fastq2 = input_files[1]
    # Prepare a command from the input and output file names,
    # as well as the pipeline parameters.
    statement = """
        hisat2
            --threads %(hisat2_threads)s
            -x %(hisat2_genome)s
            -1 %(fastq1)s
            -2 %(fastq2)s
            %(hisat2_options)s
            --summary-file %(output_file)s.log
        | samtools sort
            -@ %(hisat2_threads)s
            -o %(output_file)s
            -
        && samtools index
            %(output_file)s
        """ % locals()
    # Run the command in the Conda environment,
    # using the appopriate number of threads.
    P.run(
        statement,
        job_condaenv='pipeline_rnaseq_hisat2',
        job_threads=PARAMS['hisat2_threads'],
        job_memory=PARAMS['hisat2_memory_per_thread']
    )


@follows(mkdir("results/qc/samtools/idxstats"))
@transform(
    run_hisat2_on_fastq,
    regex(r"results/hisat2/(.*).bam"),
    r"results/qc/samtools/idxstats/\1"
)
def run_idxstats_on_bam(input_file, output_file):
    """
    Run `samtools idxstats` on the BAM files produced by HISAT2.
    """
    # Prepare a command from the input and output file names.
    statement = """
        samtools idxstats
        %(input_file)s
        > %(output_file)s
    """
    # Run the command in the Conda environment.
    P.run(statement, job_condaenv="pipeline_rnaseq_hisat2")


@follows(mkdir("results/qc/samtools/flagstat"))
@transform(
    run_hisat2_on_fastq,
    regex(r"results/hisat2/(.*).bam"),
    r"results/qc/samtools/flagstat/\1"
)
def run_flagstat_on_bam(input_file, output_file):
    """
    Run `samtools flagstat` on the BAM files produced by HISAT2.
    """
    # Prepare a command from the input and output file names.
    statement = """
        samtools flagstat
        %(input_file)s
        > %(output_file)s
    """
    # Run the command in the Conda environment.
    P.run(statement, job_condaenv="pipeline_rnaseq_hisat2")


@follows(mkdir("results/qc/picard/CollectAlignmentSummaryMetrics"))
@transform(
    run_hisat2_on_fastq,
    regex(r"results/hisat2/(.*).bam"),
    r"results/qc/picard/CollectAlignmentSummaryMetrics/\1"
)
def run_picard_alignment_metrics_on_bam(input_file, output_file):
    """
    Run `picard CollectAlignmentSummaryMetrics` on the BAM files produced by
    HISAT2.
    """
    # Fetch relevant pipeline parameters.
    picard_genome = PARAMS['picard']['genome']
    # Prepare a command from the input and output file names,
    # as well as the pipeline parameters.
    statement = """
        picard CollectAlignmentSummaryMetrics
            -I %(input_file)s
            -O %(output_file)s
            -R %(picard_genome)s
            2> %(output_file)s.log
        """
    # Run the command in the Conda environment.
    P.run(statement, job_condaenv="pipeline_rnaseq_hisat2")


@follows(mkdir("results/qc/picard/CollectInsertSizeMetrics"))
@transform(
    run_hisat2_on_fastq,
    regex(r"results/hisat2/(.*).bam"),
    r"results/qc/picard/CollectInsertSizeMetrics/\1"
)
def run_picard_insert_size_on_bam(input_file, output_file):
    """
    Run `picard CollectInsertSizeMetrics` on the BAM files produced by HISAT2.
    """
    # Fetch relevant pipeline parameters.
    picard_genome = PARAMS['picard']['genome']
    # Prepare a command from the input and output file names,
    # as well as the pipeline parameters.
    statement = """
        picard CollectInsertSizeMetrics
        -I %(input_file)s
        -O %(output_file)s
        -R %(picard_genome)s
        -H %(output_file)s.pdf
        2> %(output_file)s.log
    """
    # Run the command in the Conda environment.
    P.run(statement, job_condaenv="pipeline_rnaseq_hisat2")


@follows(mkdir("results/featureCounts"))
@merge(
    run_hisat2_on_fastq,
    "results/featureCounts/counts"
)
def run_featurecounts_on_bam(input_files, output_file):
    """
    Run `featureCounts`
    """
    # Fetch relevant pipeline parameters.
    featurecounts_gtf = PARAMS['featurecounts']['gtf']
    featurecounts_threads = PARAMS['featurecounts']['threads']
    featurecounts_options = PARAMS['featurecounts']['options']
    # Combine the list of input file names.
    input_file_list = " ".join(input_files)
    # Prepare a command from the input and output file names,
    # as well as the pipeline parameters.
    statement = """
        featureCounts
            -a %(featurecounts_gtf)s
            -o %(output_file)s
            -T %(featurecounts_threads)s
            %(featurecounts_options)s
            %(input_file_list)s
            > %(output_file)s.log
            2>&1
        """
    # Run the command in the Conda environment,
    # using the appopriate number of threads.
    P.run(
        statement,
        job_condaenv="pipeline_rnaseq_hisat2",
        job_threads=PARAMS["featurecounts"]["threads"],
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
def run_multiqc_for_bam(input_files, output_file):
    """
    Run MultiQC on the files of quality metrics computed for the BAM files.
    """
    # Prepare a command from the input folder and output file names.
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
            > %(output_file)s.log
            2>&1
    """
    # Run the command in the Conda environment.
    P.run(statement, job_condaenv="pipeline_rnaseq_hisat2")


@follows(
    run_multiqc_for_fastq,
    run_multiqc_for_bam
)
def full():
    """
    Run all the pipeline.
    """
    pass


##################
# Main execution #
##################

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
