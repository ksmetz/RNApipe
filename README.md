# RNApipe
Phanstiel lab temporary pipeline for paired-end RNAseq data

------------------------------

## Quickstart - RNApipe.py

1. **Make sure your sample info is entered into the "Sequencing Data" spreadsheet on the lab Google Drive**
	- Edit *your copy* of the sheet (if you don't have one, ask Katie or [make your own](https://docs.google.com/spreadsheets/d/13dfdyM4HxgsEoFhiuDqrohGfLulVtxYHqSK1dA29Xvw/edit#gid=2031315238))
	- The [master sheet](https://docs.google.com/spreadsheets/d/14-j6QiyzX4oV378CgQhb6btfaaopbZTRXew1FxN1vag/edit#gid=673528436) should update with your run info automatically once your personal sheet is set up.
	- <span style="background-color: #ea9999">Pink columns</span> are required. See "Instructions" tab in spreadsheet for more details

2. **Make sure the most up-to-date sample sheet is on the cluster.**
	- The sample sheet lives on the cluster, here:\
	 <code>/proj/phanstiel_lab/software/RNApipe/config/SequencingData_YYYY-MM-DD.tsv</code>
	- Check that the date of the sample sheet is past when you entered your info
	- If not, run the command to upload it:

		module load launch_pipeline
		launch googleSync

3. **Log onto longleaf using your ONYEN**

	- <code>ssh onyen@longleaf.unc.edu</code>

4. **Navigate to your folder in the lab directory**

	- <code>cd /proj/phanstiel_lab/users/ONYEN</code>

5. **Load the python 2.7.12 module**

	- <code>module load python/2.7.12</code>

6. **Run RNApipe.py, selecting for the Project you want** 
	Optional: select by Bio Rep, Condition, etc. (See [Detailed Usage: Input Selection Parameters](#input))
	- <code>python /proj/phanstiel_lab/software/RNApipe/RNApipe.py --p PROJ [options]</code>

7. **Review the run info printed to terminal screen**
	- This includes the samples selected, samples merged together, stages to be run, and full path to output directory
	- If the information looks right, answer "hell yeah" and monitor your run with <code>squeue -u onyen</code>! It's a good idea to double check that no errors appear while the jobs are submitted.
	- If the information looks wrong, answer "nah" and adjust your RNApipe.py options.

See below for more detailed usage, or use:\
<code>python /proj/phanstiel_lab/software/RNApipe/RNApipe.py --help</code>

-----------------

## Detailed Usage

### Input Selection Parameters {#input}
***Parameters that control which samples are run in the pipeline*** \
Many RNApipe.py options are used in combination to determine which files from the Sequencing Data sheet you would like to run. They correspond to the first several columns (C through M) of the data sheet. These options include:
- -p (--p, --proj, --PROJ, --project) \
**REQUIRED** (unless manual config file supplied)
- --cell (--CELL)
- --geno (--g, --GENO, --genotype)
- --cond (--COND, --condition)
- --time (--TIME)
- --br (--bior, --BIOR, --biorep, --bio)
- --tr (--techr, --TECR, --techrep, --tech)
- --sr (--seqr, --SEQR, --seqrep, --seq)
- --tag (--TAG) \
**DEFAULT:** S (ignores MiniSeq QC runs)
- -u (--user, --USER, --onyen, --ONYEN)

For all options not selected, with the exception of TAG, it will default to picking all available options for the project listed. The project option is required.

#### Manual samplesheet input
There is also a <code>--manual</code> (also <code>--MAN</code>, <code>--man</code>, <code>--MANUAL</code>) option that can override the Sequencing Data sheet. You use this command to point to a different config file, formatted the same way as the Master Sequencing Data sheet, and the program will run on all of the samples listed in that file. **Unforuntately,** for most reliable results (until Katie gets off her butt and rewrites this entire thing), it is best to format manual sheets after at least columns A-M of the Master Sheet. Additional columns can be added for metadata. 

#### Examples
To run all biological replicates for the project LIMA, you would enter:
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


### Output Parameters {#output}
***Parameters that control which samples are run in the pipeline*** \
By default, the pipeline will create an output directory relative to where the program is run. Typically, you will move to the user directory at /proj/phanstiel_lab/users/ONYEN, and it will create the following output directory structure (based on the parts of the pipeline run):
```bash
. # PWD, or path set by --home
└── project
   └── {PROJ} # Based off 'Project' column in sample sheet
   		└── rna
   			└── {PROJ}_{CELL}_{GENO}_{TRMT}_{TAG} # Auto-generated OR set with --name, appended with --suffix
   				└── proc
   					├── config		# Subset samplesheet used; merged sequence info sheet
   					├── debug		# Standard out/error log files
   					├── fastq		# Links to fastq files used
   					├── scripts		# SLURM scripts written and used
   					├── (aligned)	# BAM + BAI files
   					├── (QC)		# MultiQC summary HTML report
   					├── (quant)		# quant.sf files, tximport .RDS for use in DESeq2
   					├── (signal)	# Signal tracks (bedgraphs and/or bigwigs)
   					└── (splice)	# LSV + other MAJIQ files
```

#### Defaults
The **PROJ** name is determined by the "Project" column of the samples you selected from the Master sheet, or provided in the manual sheet. \

The **{PROJ}\_{CELL}\_{GENO}\_{TRMT}\_{TAG}** name of the final output directory is created automatically based on the "Project", "Cell_Type", "Genotype", "Condition", "Time", and "Tag" columns selected or provided. With the exception fo "Tag", it will only use columns that are shared between *all* slected samples.\

For example, if I ran *all* LIMA samples, including MiniSeq Data (using <code>--tag Q,S</code>), the name would be <code>LIMA_THP1_WT_QS</code>. 

Within <code>proc</code>, directories will be created according to your run. The directories <code>config</code>, <code>debug</code>, <code>fastq</code>, and <code>scripts</code> will always be generatd, while the others will depend on which parts of the pipeline you are running (see [Detailed Usage: Run Parameters](#running)).
- **config** - Sequencing Data sheet subet for the samples run; merged sequence info sheet
- **debug** - Standard out/error log files
- **fastq**	- Links to fastq files used
- **scripts** - SLURM scripts written and used


#### Other options
- --name (-n, --n): Replace the auto-generated name with any name of your choosing 
- --suffix (--suff, --SUFF, --SUFFIX): Add a suffix to the auto-generated name\
**NOTE:** if you provide a manual sheet, it will by default append the suffix <code>MANUAL</code> 
- --home (--o, --dir, --HOME, --OUT, --DIR): Set the output directory to somewhere other than your PWD (top of diagram)


### Run Parameters {#running}




