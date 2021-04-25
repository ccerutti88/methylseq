---
title: "User's guide of nf-core/methylSeq pipeline:"
date: "24/04/21"
output:
   html_document:
      toc: TRUE
      toc_float: TRUE
      code_folding: hide
---

**nf-core/methylseq** is a bioinformatics analysis pipeline used for Methylation (Bisulfite) sequencing data. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.

The pipeline is built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

[comment]: <All information about the *nf-core/methylseq pipeline* are available here: https://nf-co.re/methylseq/1.6>

## 1. Create a work directory:

    mkdir methylseq
    cd methylseq

## 2. Search the available *Nextflow* module and load the last version:

    search_module nextflow
    module load bioinfo/Nextflow-v20.10.0


## 3. Search the available *nf-core* module and load the last version:

    search_module nfcore
    module load bioinfo/nfcore-Nextflow-v20.10.0


## 4. Fetching the *nf-core/methylseq* pipeline code:

    nextflow pull nf-core/methylseq


## 5. Create a bash script to run the pipeline:
Create the bash script:

    touch script.sh


Copy and paste the following commands into *script.sh*:
```sh
#!/bin/bash
#SBATCH -J run20210424
#SBATCH -e /work2/genphyse/genepi/Chloe/methylseq/run20210424.err
#SBATCH -o /work2/genphyse/genepi/Chloe/methylseq/run20210424.log
#SBATCH -p workq
#SBATCH --workdir=/work2/genphyse/genepi/Chloe/methylseq
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

nextflow run nf-core/methylseq -profile genotoul \
                               -r 1.5 \
                               -name run20210424 \
                               --reads '/work2/genphyse/genepi/Chloe/methylseq/data/FRPL*.fastq.gz' \
                               --single_end \
                               --cytosine_report \
                               --fasta '/work2/genphyse/genepi/Chloe/methylseq/data/Coturnix_japonica.Coturnix_japonica_2.0.dna.toplevel.fa' \
                               --rrbs \
                               --email chloe.cerutti@inrae.fr\
                               --max_memory 64GB \
                               --max_cpus 16

```
[comment]: <Adapt the code with your work directories, paths, email,...>

Save and escape.

**SBATCH paramaters:**
- *-J*: specify a name for the job allocation
- *-e*: instruct Slurm to connect the batch script's standard error directly to the file name specified in the "filename pattern"
- *-o*: path to the file where the job (error) output is written to
- *-p*: request a specific partition for the resource allocation
- *--wordir*: path of the working directory
- *--mail-type*: turn on mail notification; type can be one of *BEGIN* (start of execution), *END*, *FAIL*, *REQUEUE* or *ALL*
- *--export*: identify which environment variables from the submission environment are propagated to the launched application. If *--export=ALL*, all of the users environment will be loaded
- *--cpus-per-task*: just one CPU core will be used.
- *--mem*: Memory (RAM) for the job. Number followed by unit prefix (4G)

[comment]: <Details about slurm parameters are available here: https://slurm.schedmd.com/sbatch.html>

**nf-core/methyseq pipeline paramaters:**
- *-profile*: specify the profile to download and use (example: profile *genotoul.config*)
- *-r*: specify the pipeline version number
- *-name*: specify job name
- *--reads*: path to fastq.gz reads files
- *--single_end*: specifies that the input is single-end reads
- *--cytosine_report*: output stranded cytosine report during Bismark's bismark_methylation_extractor step.
- *--fasta*: path to FASTA genome file
- *--rrbs*: turn on if dealing with MspI digested material
- *--email*: email address for completion summary
- *--max_memory*: maximum of memory allocated
- *--max_cpus*: maximum number of cpus allocated

[comment]: <Details about pipeline parameters are available here: https://nf-co.re/methylseq/1.6/parameters>

## 6. Run the script on a node cluster:

    sbatch script.sh


## 7. Results:
All output files are available into the *results* directory.

### FastQC: Read quality control:
FastQC gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences.

**Output files:**
- *fastqc/ (*_fastqc.html*):* FastQC report containing quality metrics for your untrimmed raw fastq files.
- *fastqc/zips/ (*_fastqc.zip*):* Zip archive containing the FastQC report, tab-delimited data file and plot images.
sample_fastqc.html
- *FastQC report*: containing quality metrics for your untrimmed raw fastq files
- *sample_fastqc.zip*: zip file containing the FastQC report, tab-delimited data file and plot images

### TrimGalore: Adapter trimming:
The nf-core/methylseq pipeline uses TrimGalore! for removal of adapter contamination and trimming of low quality regions. TrimGalore is a wrapper around Cutadapt and runs FastQC after it finishes.
MultiQC reports the percentage of bases removed by Cutadapt in the General Statistics table, along with a line plot showing where reads were trimmed.

**Output directory:** *results/trim_galore*
Contains FastQ files with quality and adapter trimmed reads for each sample, along with a log file describing the trimming.
- *sample_val_1.fq.gz*, *sample_val_2.fq.gz* : Trimmed FastQ data, reads 1 and 2
**NB:** Only saved if --save_trimmed has been specified.
- *logs/sample_val_1.fq.gz_trimming_report.txt*: Trimming report (describes which parameters that were used)
- *FastQC/sample_val_1_fastqc.zip*: FastQC report for trimmed reads
Single-end data will have slightly different file names and only one FastQ file per sample.

### Alignment: Aligning reads to reference genome:
Bismark and bwa-meth convert all Cytosines contained within the sequenced reads to Thymine in-silico and then align against a three-letter reference genome. This method avoids methylation-specific alignment bias. The alignment produces a BAM file of genomic alignments.

**Bismark output directory:** *results/bismark_alignments/*
- *sample.bam*: Aligned reads in BAM format.
**NB:** Only saved if --save_align_intermeds, --skip_deduplication or --rrbs is specified when running the pipeline.
- *logs/sample_PE_report.txt*: Log file giving summary statistics about alignment.
- *unmapped/unmapped_reads_1.fq.gz*, *unmapped/unmapped_reads_2.fq.gz*: Unmapped reads in FastQ format. Only saved if *--unmapped* specified when running the pipeline.

### Methylation Extraction: Calling cytosine methylation steps:
The methylation extractor step takes a BAM file with aligned reads and generates files containing cytosine methylation calls. It produces a few different output formats, described below.

**Bismark output directory:** *results/bismark_methylation_calls/*
- *methylation_calls/XXX_context_sample.txt.gz*: Individual methylation calls, sorted into files according to cytosine context.
- *methylation_coverage/sample.bismark.cov.gz*: Coverage text file summarising cytosine methylation values.
- *bedGraph/sample.bedGraph.gz*: Methylation statuses in bedGraph format, with 0-based genomic start and 1- based end coordinates.
- *m-bias/sample.M-bias.txt*: QC data showing methylation bias across read lengths. See the bismark documentation for more information.
- *logs/sample_splitting_report.txt*: Log file giving summary statistics about methylation extraction.

### Bismark Reports: Single-sample and summary analysis reports:
Bismark generates a HTML reports describing results for each sample, as well as a summary report for the whole run.

**Output directory:** *results/bismark_summary*

### Qualimap: Tool for genome alignments QC:
Qualimap BamQC is a general-use quality-control tool that generates a number of statistics about aligned BAM files. It's not specific to bisulfite data, but it produces several useful stats - for example, insert size and coverage statistics.

**Output directory:** *results/qualimap*
- *sample/qualimapReport.html*: Qualimap HTML report
- *sample/genome_results.txt*, *sample/raw_data_qualimapReport/.txt*: Text-based statistics that can be loaded into downstream programs

### Preseq: Tool for estimating sample complexity:
Preseq estimates the complexity of a library, showing how many additional unique reads are sequenced for increasing the total read count. A shallow curve indicates that the library has reached complexity saturation and further sequencing would likely not add further unique reads. The dashed line shows a perfectly complex library where total reads = unique reads.
Note that these are predictive numbers only, not absolute. The MultiQC plot can sometimes give extreme sequencing depth on the X axis - click and drag from the left side of the plot to zoom in on more realistic numbers.

**Output directory:** *results/preseq*
- *sample_ccurve.txt*: This file contains plot values for the complexity curve, plotted in the MultiQC report.

### MultiQC: Aggregate report describing results of the whole pipeline:
MultiQC is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability.

For more information about how to use MultiQC reports, see https://multiqc.info.

**Output files:** *results/MultiQC/*:
- *multiqc_report.html*: a standalone HTML file that can be viewed in your web browser.
- *multiqc_report_data/*: directory containing parsed statistics from the different tools used in the pipeline.
- *multiqc_plots/*: directory containing static images from the report in various formats.

### Pipeline information: Report metrics generated during the workflow execution:
Nextflow provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

**Output files:** *results/pipeline_info/*:
- Reports generated by Nextflow: *execution_report.html*, *execution_timeline.html*, *execution_trace.txt* and *pipeline_dag.dot/pipeline_dag.svg*.
- Reports generated by the pipeline: *pipeline_report.html*, *pipeline_report.txt* and *software_versions.csv*.
- Documentation for interpretation of results in HTML format: *results_description.html*.

[comment]: <Details about the output produced by the pipeline are available here: https://nf-co.re/methylseq/1.6/output>
