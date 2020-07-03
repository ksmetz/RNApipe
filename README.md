# RNApipe
Phanstiel lab temporary pipeline for paired-end RNAseq data


## Quickstart - RNApipe.py
------------------------------

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
Optional: select by Bio Rep, Condition, etc. (See **RNApipe.py Detailed Usage: *Selection Parameters***)
<code>python /proj/phanstiel_lab/software/RNApipe/RNApipe.py --p PROJ [options]</code>

7. **Review the run info printed to terminal screen**
- This includes the samples selected, samples merged together, stages to be run, and full path to output directory
- If the information looks right, answer "hell yeah" and monitor your run with <code>squeue -u onyen</code>! It's a good idea to double check that no errors appear while the jobs are submitted.
- If the information looks wrong, answer "nah" and adjust your RNApipe.py options.

See below for more detailed usage, or use <code>python /proj/phanstiel_lab/software/RNApipe/RNApipe.py --help</code>


