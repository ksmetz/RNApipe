# RNApipe
Phanstiel lab (temporary) pipeline for paired-end RNAseq data.

#### TABLE OF CONTENTS
- [Quickstart](#quickstart)
- [Detailed Usage](#detailed-usage)
	- [Input Parameters](#input-parameters)
	- [Output Parameters](#output-parameters)
	- [Run Parameters](#run-parameters)
	- [Command Parameters](#command-parameters)
- [Commands Run](#commands-run)
	- [CORE scripts](#core-scripts)
	- [MERGE scripts](#merge-scripts)
	- [FINAL scripts](#final-scripts)

------------------------------

## Quickstart

1. **Make sure your sample info is entered into the "Sequencing Data" spreadsheet on the lab Google Drive**
	- Edit *your copy* of the sheet (if you don't have one, ask Katie or [make your own](https://docs.google.com/spreadsheets/d/13dfdyM4HxgsEoFhiuDqrohGfLulVtxYHqSK1dA29Xvw/edit#gid=2031315238))
	- The [master sheet](https://docs.google.com/spreadsheets/d/14-j6QiyzX4oV378CgQhb6btfaaopbZTRXew1FxN1vag/edit#gid=673528436) should update with your run info automatically once your personal sheet is set up.
	- <span style="background-color: #ea9999">Pink columns</span> are required. See "Instructions" tab in spreadsheet for more details

2. **Make sure the most up-to-date sample sheet is on the cluster.**
	- The sample sheet lives on the cluster, here:\
	 <code>/proj/phanstiel_lab/software/RNApipe/config/SequencingData_YYYY-MM-DD.tsv</code>
	- Check that the date of the sample sheet is past when you entered your info
	- If not, run the command to upload it:
```
		module load launch_pipeline
		launch googleSync
```
3. **Log onto longleaf using your ONYEN**

	- <code>ssh onyen@longleaf.unc.edu</code>

4. **Navigate to your folder in the lab directory**

	- <code>cd /proj/phanstiel_lab/users/ONYEN</code>

5. **Load the python 2.7.12 module**

	- <code>module load python/2.7.12</code>

6. **Run RNApipe.py, selecting for the Project you want** \
	Note: See [Detailed Usage](#detailed-usage) for other selection/input options, output customization options, and pipeline options.
	- <code>python /proj/phanstiel_lab/software/RNApipe/RNApipe.py --p PROJ [options]</code>

7. **Review the run info printed to terminal screen**
	- This includes the samples selected, samples merged together, stages to be run, and full path to output directory
	- If the information looks right, answer <code>hell yeah</code> and monitor your run with <code>squeue -u onyen</code>! It's a good idea to double check that no errors appear while the jobs are submitted.
	- If the information looks wrong, answer <code>nah</code> and adjust your RNApipe.py options.

See below for more detailed usage, or use:\
<code>python /proj/phanstiel_lab/software/RNApipe/RNApipe.py --help</code>

-----------------

## Detailed Usage

### Input Parameters 
***Parameters that control which samples are run in the pipeline*** \
Many RNApipe.py options are used in combination to determine which files from the Sequencing Data sheet you would like to run. They correspond to the first several columns (C through M) of the data sheet. These options include:
- <code>-p</code> (<code>--p</code>, <code>--proj</code>, <code>--PROJ</code>, <code>--project</code>) \
**REQUIRED** (unless manual config file supplied)
- <code>--cell</code> (<code>--CELL</code>)
- <code>--geno</code> (<code>--g</code>, <code>--GENO</code>, <code>--genotype</code>)
- <code>--cond</code> (<code>--COND</code>, <code>--condition</code>)
- <code>--time</code> (<code>--TIME</code>)
- <code>--br</code> (<code>--bior</code>, <code>--BIOR</code>, <code>--biorep</code>, <code>--bio</code>)
- <code>--tr</code> (<code>--techr</code>, <code>--TECR</code>, <code>--techrep</code>, <code>--tech</code>)
- <code>--sr</code> (<code>--seqr</code>, <code>--SEQR</code>, <code>--seqrep</code>, <code>--seq</code>)
- <code>--tag</code> (<code>--TAG</code>) \
**DEFAULT:** S (ignores MiniSeq QC runs)
- <code>-u</code> (<code>--user</code>, <code>--USER</code>, <code>--onyen</code>, <code>--ONYEN</code>)

For all options not selected, with the exception of TAG, it will default to picking all available options for the project listed. The project option is required.

#### Manual samplesheet input
There is also a <code>--manual</code> (also <code>--MAN</code>, <code>--man</code>, <code>--MANUAL</code>) option that can override the Sequencing Data sheet. You use this command to point to a different config file, formatted the same way as the Master Sequencing Data sheet, and the program will run on all of the samples listed in that file. **Unforuntately,** for most reliable results (until Katie gets off her butt and rewrites this entire thing), it is best to format manual sheets after at least columns A-M of the Master Sheet. Additional columns can be added for metadata. 

#### Examples
To run all biological replicates for the project LIMA, you would enter: \
<code>python /proj/phanstiel_lab/software/RNApipe/RNApipe.py --proj LIMA </code>

To run biological reps 2 and 3 (but not 1) of the LIMA project, you could enter: \
<code>python /proj/phanstiel_lab/software/RNApipe/RNApipe.py --proj LIMA --br 2,3</code>

Alternatively, to run these samples with a manual sheet, you could make a sheet that looks as follows and run this command: \
<code>python /proj/phanstiel_lab/software/RNApipe/RNApipe.py --manual mySampleSheet.tsv</code>

| Unique_ID | Name | Project | Cell_Type | Genotype | Condition | Time | Bio_Rep | Tech_Rep | Seq_Rep | Tag |
|-----------|------|---------|-----------|----------|-----------|------|---------|----------|---------|-----|
| KR51EBE4	| LIMA_RNA_THP1_WT_LPIF_0000_S_3.1.1	| LIMA	| THP1	| WT	| LPIF	| 0000	| 3	| 1	| 1	| S |
| KR0A9DF4	| LIMA_RNA_THP1_WT_LPIF_0030_S_3.1.1	| LIMA	| THP1	| WT	| LPIF	| 0030	| 3	| 1	| 1	| S |
| KR751205	| LIMA_RNA_THP1_WT_LPIF_0060_S_3.1.1	| LIMA	| THP1	| WT	| LPIF	| 0060	| 3	| 1	| 1	| S |
| KR02C443	| LIMA_RNA_THP1_WT_LPIF_0090_S_3.1.1	| LIMA	| THP1	| WT	| LPIF	| 0090	| 3	| 1	| 1	| S |
| KR08B759	| LIMA_RNA_THP1_WT_LPIF_0120_S_3.1.1	| LIMA	| THP1	| WT	| LPIF	| 0120	| 3	| 1	| 1	| S |
| KR365BAB	| LIMA_RNA_THP1_WT_LPIF_0240_S_3.1.1	| LIMA	| THP1	| WT	| LPIF	| 0240	| 3	| 1	| 1	| S |
| KR8A7D54	| LIMA_RNA_THP1_WT_LPIF_0360_S_3.1.1	| LIMA	| THP1	| WT	| LPIF	| 0360	| 3	| 1	| 1	| S |
| KR0CE283	| LIMA_RNA_THP1_WT_LPIF_1440_S_3.1.1	| LIMA	| THP1	| WT	| LPIF	| 1440	| 3	| 1	| 1	| S |
| KR031F3D	| LIMA_RNA_THP1_WT_LPIF_0000_S_2.1.1	| LIMA	| THP1	| WT	| LPIF	| 0000	| 2	| 1	| 1	| S |
| KR13130A	| LIMA_RNA_THP1_WT_LPIF_0030_S_2.1.1	| LIMA	| THP1	| WT	| LPIF	| 0030	| 2	| 1	| 1	| S |
| KR772B9D	| LIMA_RNA_THP1_WT_LPIF_0060_S_2.1.1	| LIMA	| THP1	| WT	| LPIF	| 0060	| 2	| 1	| 1	| S |
| KR62E3AC	| LIMA_RNA_THP1_WT_LPIF_0090_S_2.1.1	| LIMA	| THP1	| WT	| LPIF	| 0090	| 2	| 1	| 1	| S |
| KR2BA90C	| LIMA_RNA_THP1_WT_LPIF_0120_S_2.1.1	| LIMA	| THP1	| WT	| LPIF	| 0120	| 2	| 1	| 1	| S |
| KR79BBD3	| LIMA_RNA_THP1_WT_LPIF_0240_S_2.1.1	| LIMA	| THP1	| WT	| LPIF	| 0240	| 2	| 1	| 1	| S |
| KR859E81	| LIMA_RNA_THP1_WT_LPIF_0360_S_2.1.1	| LIMA	| THP1	| WT	| LPIF	| 0360	| 2	| 1	| 1	| S |
| KR3136D3	| LIMA_RNA_THP1_WT_LPIF_1440_S_2.1.1	| LIMA	| THP1	| WT	| LPIF	| 1440	| 2	| 1	| 1	| S |

[Back to Top](#table-of-contents)

-----------------

### Output Parameters
***Parameters that control the name and location of the output files*** \
By default, the pipeline will create an output directory relative to where the program is run. Typically, you will move to the user directory at /proj/phanstiel_lab/users/ONYEN, and it will create the following output directory structure (based on the parts of the pipeline run):
```bash
. # PWD, or path set by --home
└── project
    └── {PROJ} # Based off 'Project' column in sample sheet
        └── rna
            └── {NAME} # Auto-generated OR set with --name, appended with --suffix
                └── proc
                    ├── config    # Subset samplesheet used; merged sequence info sheet
                    ├── debug     # Standard out/error log files
                    ├── fastq     # Links to fastq files used
                    ├── scripts   # SLURM scripts written and used
                    ├── (aligned) # BAM + BAI files
                    ├── (QC)      # MultiQC summary HTML report
                    ├── (quant)   # quant.sf files, tximport .RDS for use in DESeq2
                    ├── (signal)  # Signal tracks (bedgraphs and/or bigwigs)
                    └── (splice)  # LSV + other MAJIQ files
```

#### Defaults
The **{PROJ}** name is determined by the "Project" column of the samples you selected from the Master sheet, or provided in the manual sheet. 

The **{NAME}** of the final output directory is created automatically based on the "Project", "Cell_Type", "Genotype", "Condition", "Time", and "Tag" columns selected or provided. With the exception fo "Tag", it will only use columns that are shared between *all* slected samples.

For example, if I ran *all* LIMA samples, including MiniSeq Data (using <code>--tag Q,S</code>), the name would be <code>LIMA_THP1_WT_QS</code>. 

Within <code>proc</code>, directories will be created according to your run. The directories <code>config</code>, <code>debug</code>, <code>fastq</code>, and <code>scripts</code> will always be generatd, while the others will depend on which parts of the pipeline you are running (see [Detailed Usage: Run Parameters](#run-parameters)).


#### Options
- <code>--name</code> (<code>-n</code>, <code>--n</code>): Replace the auto-generated {NAME} with any name of your choosing 
- <code>--suffix</code> (<code>--suff</code>, <code>--SUFF</code>, <code>--SUFFIX</code>): Add a suffix to the auto-generated {NAME}\
**NOTE:** if you provide a manual sheet, it will by default append the suffix <code>MANUAL</code> 
- <code>--home</code> (<code>--o</code>, <code>--dir</code>, <code>--HOME</code>, <code>--OUT</code>, <code>--DIR</code>): Set the output directory to somewhere other than your PWD (top of diagram)

[Back to Top](#table-of-contents)

-----------------

### Run Parameters 
***Parameters that control what parts of the pipeline to run*** \
By default, the pipeline will trim low-quality reads, quantify transcripts for use in DESeq2, and align transcripts to create signal tracks and splicing info (see map below). There are several temporary files created in this process that are deleted by default during the pipeline (grey boxes below). If files already exist from previous runs, the default behavior is to overwrite them. These default settings can be changed by many run parameters, described below. 

![](img/RNApipe_flow-01.png)

- <code>--stage</code> (<code>--STAGE</code>, <code>--S</code>): Set what stages you want the pipeline to run, separated by commas **without spaces** \
**DEFAULT:** <code>QC,trim,align,quant,signal</code>\
Options include:
	- QC: Create individual and summary reports of QC metrics for sequencing and alignment, using FastQC + MultiQC
	- trim: Remove adapters and low-quality reads from fastq files, using Trim Galore!
	- align: Align reads to transcriptome (required for signal tracks or splicing info), using HISAT2 and Samtools
	- quant: Quantify reads and summarize at a gene level for use in DESeq2, using Salmon and txImport
	- splice: Identify splicing variants, using MAJIQ
	- signal: Create bedgraph or bigwig signal tracks
	- merge: Combine aligned files to make merged signal tracks
- <code>--merge</code> (<code>--mergeby</code>, <code>--MERGE</code>): Select how you would like to merge samples \
**DEFAULT:** If you include <code>--stage merge</code> it will combine any samples that have the same info in all columns except for Bio_Rep, Tech_Rep, and Seq_Rep (i.e. <code>LIMA_THP1_WT_0000_1.1.1_Q</code> would merge with <code>LIMA_THP1_WT_0000_2.3.5_Q</code> but not <code>LIMA_THP1_WT_0030_1.1.1_Q</code>.)
- <code>--temp</code> (<code>--TEMP</code>, <code>--keep</code>, <code>--KEEP</code>): True/False for whether to keep temporary files \
**DEFAULT:** False
- <code>--signalout</code> (<code>--sigout</code>, <code>--SIGNALOUT</code>, <code>--SIGOUT</code>): Select signal output file type. Options are "bigwig" and "bedgraph"\
**DEFAULT:** bigwig
- <code>--rerun</code> (<code>--RERUN</code>): True/False for whether to overwrite existing files that might exist from previous runs \
**DEFAULT:** True

[Back to Top](#table-of-contents)

-----------------

### Command Parameters 
***Parameters that are fed directly into the program commands*** \
The pipeline will run many programs. Some of the options for these programs are hard-coded, while others can be changed through the command parameters of RNApipe.py. These options include:

#### *txImport* (for --stage quant)
- <code>--ann</code> (<code>--ANN</code>, <code>--annotation</code>, <code>--gtf</code>, <code>--GTF</code>): Path to transcriptome annotation file\
**DEFAULT:** hg19 \
<code>/proj/phanstiel_lab/SHARE/genomes/GENCODE.v19/gencode.v19.annotation.gtf_withproteinids</code>

#### *Salmon* (for --stage quant)
- <code>--trns</code> (<code>--TRNS</code>, <code>--transcript</code>): Path to Salmon transcript index \
**DEFAULT:** hg19 \
<code>/proj/phanstiel_lab/SHARE/geonmes/GENCODE.v19/salmon_index</code>

#### *HISAT2* (for --stage align)
- <code>--hisatidx</code> (<code>--HSIDX</code>, <code>--hisat2idx</code>, <code>--hisat2index</code>): Path to HISAT2 index \
**DEFAULT:** hg19 \
<code>/proj/seq/data/HG19_UCSC/Sequence/HISAT2Index/genome</code>

#### *MAJIQ* (for --stage splice)
- <code>--majConfig</code> (<code>--mfqConfig</code>, <code>--mjqcon</code>, <code>--majcon</code>): Path to custom MAJIQ configuration file \
**DEFAULT:** automatically generated config based on sample config, other options \
<code>{HOME}/project/{PROJ}/rna/{NAME}/proc/config/MAJIQconfig_NAME_YYMMDD_HHMMSS.txt </code>
- <code>--length</code> (<code>--readlength</code>, <code>--readlen</code>, <code>--len</code>): Read length, for MAJIQ config file \
**DEFAULT:** 150
- <code>--genome</code>: Name of genome used, for MAJIQ config file \
**DEFAULT:** hg19
- <code>--strand</code> (<code>--strandness</code>): Strandness, chosen from ‘forward’, ‘reverse’ and ‘None’, for MAJIQ config file \
**DEFAULT:** forward
- <code>--mjqann</code> (<code>--MJQANN</code>, <code>--majiqannotation</code>, <code>-gff3</code>, <code>--GFF3</code>): GFF3 transcriptome annotation \
**DEFAULT:** hg19 \
<code>/proj/phanstiel_lab/SHARE/genomes/GENCODE.v19/gencode.v19.annotation.gtf_withproteinids</code>

[Back to Top](#table-of-contents)

-----------------

## Commands Run
***Overview of commands run by RNApipe.py SLURM scripts***\
This section will roughly review the commands (including built-in and user-provided settings) run by RNApipe.py, given each stage selected. Paths will be abbreviated with (...) replacing <code>{HOME}/project/{RPOJ}/rna/{NAME}/proc/</code>. “Name” represents the name of the sample (from the “Name” column of Config), while "{NAME}" represents the project name generated automatically or provided by the --name command (see Output Parameters). “MergeName” represents the combined name from merged samples (see [Detailed Usage: Output Parameters](#output-parameters)).

#### CORE SCRIPTS
These commands are written and run for every sample (line) in the samplesheet provided.

##### QC stage
1. <code>fastqc -o .../QC .../fastq/Name\_R1.fastq.gz .../fastq/Name\_R2.fastq.gz</code>

##### trim stage
2.  <code>trim\_galore -o .../fastq --paired .../fastq/Name\_R1.fastq.gz .../fastq/Name\_R2.fastq.gz</code>
3.  <code>mv .../fastq/Name\_R1\_val_1.\* .../fastq/Name\_1\_trimmed.fq.gz</code>
4.  <code>mv .../fastq/Name\_R1\*\_trimming\_report.txt ../fastq/Name\_1.fastq\_trimming\_report.txt</code>
5.  <code>mv .../fastq/Name\_R2\_val\_2.\* .../fastq/Name\_2\_trimmed.fq.gz</code>
6.  <code>mv .../fastq/Name\_R2\*\_trimming\_report.txt ../fastq/Name\_2.fastq\_trimming\_report.txt</code>

##### quant stage
7. <code>salmon quant --writeUnmappedNames --threads 1 -i {TRNS} -l A -1 .../fastq/Name\_1\_trimmed.fq.gz -2 .../fastqName\_2\_trimmed.fq.gz -o .../quant/Name</code>

##### align stage
8. <code>hisat2 -q -x HSIDX -1 .../fastq/Name\_1\_trimmed.fq.gz -2 .../fastq/Name\_2\_trimmed.fq.gz \| samtools view -u \| samtools sort -o .../align/Name\_sorted.bam</code>
9. <code>samtools flagstat .../align/Name\_sorted.bam > .../align/Name\_stats.txt</code>
10. <code>samtools index .../align/Name\_sorted.bam</code>

##### temp/keep set to False
11. <code>rm .../fastq/Name\_1\_trimmed.fq.gz</code>
12. <code>rm .../fastq/Name\_2\_trimmed.fq.gz</code>

##### signal stage, signalout set to bedgraph
13. <code>bedtools genomecov -bga -ibam .../align/Name\_sorted.bam > .../signal/Name.bedgraph</code>

##### signal stage, signalout set to bigwig
14. <code>bamCoverage -b .../align/Name\_sorted.bam -o .../signal/Name.bw</code>

#### MERGE SCRIPTS
These commands are run for every merged sample in <code>.../config/mergeList</code>, as deteremined by <code>--stage merge</code> and <code>--mergeby</code>.

##### merge stage
1. <code>samtools merge .../align/MERGE\_MergeName.bam .../align/Name1\_sorted.bam .../align/Name2\_sorted.bam</code>
2. <code>samtools index .../align/MERGE\_MergeName.bam</code>

##### signal stage, signalout set to bedgraph
3. <code>bedtools genomecov -bga -ibam .../align/MERGE\_MergeName.bam > .../signal/MERGE\_MergeName.bedgraph</code>

##### signal stage, signalout set to bigwig
4. <code>bamCoverage -b .../align/MERGE\_MergeName.bam -o .../signal/MERGE\_MergeName.bw</code>

##### temp/keep set to False
5. <code>rm .../align/MERGE\_MergeName.bam</code>
6. <code>rm .../align/MERGE\_MergeName.bam.bai</code>

#### FINAL SCRIPTS
Run only one time for the entire samplesheet

##### QC stage
1. <code>multiqc -f .../proc/* -o .../QC</code>
2. <code>mv .../QC/multiqc\_report.html .../QC/NAME\_multiqc\_report.html</code>

##### temp/keep set to False
3. <code>rm .../QC/\*fastqc.zip</code>
4. <code>rm .../QC/\*fastqc.html</code>
5. <code>rm .../fastq/\*fastq\_trimming\_report.txt</code>

##### quant stage
6. <code>Rscript /proj/phanstiel\_lab/software/RNApipe/txImporter.R .../config/config\_NAME\_YYMMDD\_HHMMSS.tsv {ANNT} .../quant {NAME}</code>

##### splice stage
7. <code>majiq build {MJQANN} -c {MJQCON} -j 8 -o .../splice</code>

[Back to Top](#table-of-contents)