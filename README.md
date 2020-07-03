# RNApipe
Phanstiel lab temporary pipeline for paired-end RNAseq data

------------------------------

## Quickstart - RNApipe.py

1. **Make sure your sample info is entered into the "Sequencing Data" spreadsheet on the lab Google Drive**
	- Edit *your copy* of the sheet (if you don't have one, ask Katie or [make your own](https://docs.google.com/spreadsheets/d/13dfdyM4HxgsEoFhiuDqrohGfLulVtxYHqSK1dA29Xvw/edit#gid=2031315238))
	- The [master sheet](https://docs.google.com/spreadsheets/d/14-j6QiyzX4oV378CgQhb6btfaaopbZTRXew1FxN1vag/edit#gid=673528436) should update with your run info automatically once your personal sheet is set up.
	- <span style="background-color: #ea9999">Pink columns</span> are required. See "Instructions" tab in spreadsheet for more details

2. **Make sure the most up-to-date sample sheet is on the cluster.**
	- The sample sheet lives on the cluster, here: <code>/proj/phanstiel_lab/software/RNApipe/config/SequencingData_YYYY-MM-DD.tsv</code>
	- Check that the date of the sample sheet is past when you entered your info
	- If not, run the command to upload it:

		module load launch_pipeline
		launch googleSync

3. **Log onto longleaf using your ONYEN**

<code>ssh onyen@longleaf.unc.edu</code>

4. **Navigate to your folder in the lab directory**

<code>cd /proj/phanstiel_lab/users/ONYEN</code>

5. **Load the python 2.7.12 module**

<code>module load python/2.7.12</code>

6. **Run RNApipe.py, selecting for the Project you want** 

Optional: select by Bio Rep, Condition, etc. (See **Detailed Usage: Input Selection Parameters**)
<code>python /proj/phanstiel_lab/software/RNApipe/RNApipe.py --p PROJ [options]</code>

7. **Review the run info printed to terminal screen**
	- This includes the samples selected, samples merged together, stages to be run, and full path to output directory
	- If the information looks right, answer "hell yeah" and monitor your run with <code>squeue -u onyen</code>! It's a good idea to double check that no errors appear while the jobs are submitted.
	- If the information looks wrong, answer "nah" and adjust your RNApipe.py options.

See below for more detailed usage, or use:
<code>python /proj/phanstiel_lab/software/RNApipe/RNApipe.py --help</code>

-----------------

## Detailed Usage

### Input Selection Parameters
***Parameters that control which samples are run in the pipeline***
Many RNApipe.py options are used in combination to determine which files from the Sequencing Data sheet you would like to run. They correspond to the first several columns (C through M) of the data sheet. These options include:
- -p (--p, --proj, --PROJ, --project) 
**REQUIRED** (unless manual config file supplied)
- --cell (--CELL)
- --geno (--g, --GENO, --genotype)
- --cond (--COND, --condition)
- --time (--TIME)
- --br (--bior, --BIOR, --biorep, --bio)
- --tr (--techr, --TECR, --techrep, --tech)
- --sr (--seqr, --SEQR, --seqrep, --seq)
- --tag (--TAG)
**DEFAULT:** S (ignores MiniSeq QC runs)
- -u (--user, --USER, --onyen, --ONYEN)

For all options not selected, with the exception of TAG, it will default to picking all available options for the project listed. The project option is required.

#### Manual samplesheet input
There is also a <code>--manual</code> (also <code>--MAN</code>, <code>--man</code>, <code>--MANUAL</code>) option that can override the Sequencing Data sheet. You use this command to point to a different config file, formatted the same way as the Master Sequencing Data sheet, and the program will run on all of the samples listed in that file. For most reliable results (until Katie gets off her butt and rewrites this entire thing), it is best to keep columns A-M from the Master Sheet when providing a manual sheet. Additional columns can be added for metadata.

#### Examples
To run biological reps 2 and 3 (but not 1) of the LIMA project, you could enter:
<code>python /proj/phanstiel_lab/software/RNApipe/RNApipe.py --proj LIMA --br 2,3</code>

To run these samples with a manual sheet, you could provide a sheet that looks as follows:

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

To run all listed biological replicates, it would simply be:
<code>python /proj/phanstiel_lab/software/RNApipe/RNApipe.py --proj LIMA </code>



