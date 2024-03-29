#!/usr/bin/python

##########
# This is the RNA-Seq processing pipeline
#=========
# AUTHOR: Katie Metz Reed
# UPDATED: 2019-02-14
#=========
# This is what the pipeline does. (Update this at some point)
##########

# UPDATE: reconsider memory/node requirements?
# UPDATE: consider making MERGE jobs wait only on their respective samples
# UPDATE: change description!
# UPDATE: add example use to help, or "no project" error


#-------------------------------------- I N I T I A L I Z E --------------------------------------#

# Load any required libraries
import optparse
import pandas as pd
import os
import time
import sys
import commands


# Time stamp
print time.asctime(time.localtime())
stamp = time.strftime("%y%m%d_%H%M%S")
print "Run started ..."


# Functions
# OPTIONLISTER
#     Splits multi-item selection paramters into lists
def optionLister(optionString):
    if "," in optionString:
        optionList = optionString.split(",")
    else:
        optionList = [optionString]
    return(optionList) 

# NAMER 
#    Generates a project name based on the samples run by listing common column fields (1) and concatonating others (2)
def namer(columnNamesSingle, columnNamesCombine):

    nameList = []

    # For columns which you want to list ONLY if there is a single one
    for column in columnNamesSingle:
        columnSet = map(str, set(config[column]))       
        if len(columnSet) == 1:
            nameList.append(str(columnSet[0]))

    # For columns which you want to squish together all used
    for column in columnNamesCombine:
        columnSet = map(str, set(config[column])) 
        nameList.append("".join(columnSet))

    return("_".join(nameList))





#------------------------------------------- I N P U T -------------------------------------------#

# NON-ARG OPTIONS
# ===============

# Variables not passed through arguments (for now), but made variables to edit code easily in future

# version numbers for different programs used
FASTQCvers  = "0.11.5"
TrimVers    = "0.4.3"
SalmonVers  = "1.4.0"
rVers       = "3.3.1"
multiQCvers = "1.5" 
majiqVers   = "1.1.7"
hisatVers   = "2.1.0"
samVers     = "1.9"
bedVers     = "2.26"
deepVers    = "3.0.1"

# path to directory containing master config file
CDIR = "/proj/phanstiel_lab/software/RNApipe/config/" #cluster
# CDIR = "/Users/ksmetz/Phanstiel Lab Dropbox/Kathleen Metz/Work/Research/Code/Pipelines/local/GoogleSync/config/" #laptop
# CDIR = "/Users/phanstiel3/Dropbox/Work/Research/Code/Pipelines/local/GoogleSync/config/" #desktop



# READ IN OPTIONS 
# ===============

# Read in parameters with optparse
parser = optparse.OptionParser()


# SELECTION PARAMETERS
# Used to select rows from master config file

selGroup = optparse.OptionGroup(parser, "Selection Parameters",
                    "Use these options to select your desired samples from the master config file. "
                    "The master config lives in /proj/phanstiel_lab/software/RNApipe/config/. "
                    "Alternatively, provide your own config file with --manual.")

# Project (Column 3, C)
selGroup.add_option("-p", "--p", "--proj", "--PROJ", "--project",action = "store", type = "string", dest = "proj", 
                  help = "Name of project. Only one may be selected. Required unless supplying manual config file with --manual",default="all")
# Cell_Type (Column 4, D)
selGroup.add_option("--cell", "--CELL",action = "store", type = "string", dest = "cell", 
                  help = "Cell type selected.",default="all")
# Genotype (Column 5, E)
selGroup.add_option("-g", "--geno", "--GENO", "--genotype",action = "store", type = "string", dest = "geno", 
                  help = "Genotype selected.",default="all")
# Condition (Column 6, F)
selGroup.add_option("--cond", "--COND", "--condition",action = "store", type = "string", dest = "cond", 
                  help = "Conditions selected.",default="all")
# Time (Column 7, G)
selGroup.add_option("--time", "--TIME",action = "store", type = "string", dest = "time", 
                  help = "Time points selected.",default="all")
# Bio_Rep (Column 8, H)
selGroup.add_option("--br", "--bior", "--BIOR", "--biorep", "--bio",action = "store", type = "string", dest = "bior", 
                  help = "Biological replicates selected.",default="all")
# Tech_Rep (Column 9, I)
selGroup.add_option("--tr", "--techr", "--TECR", "--techrep", "--tech",action = "store", type = "string", dest = "tecr", 
                  help = "Technical replicates selected.",default="all")
# Seq_Rep (Column 10,J)
selGroup.add_option("--sr", "--seqr", "--SEQR", "--seqrep", "--seq",action = "store", type = "string", dest = "seqr", 
                  help = "Sequencing replicates selected.",default="all")
# Tag (Column 11, K)
selGroup.add_option("--tag", "--TAG",action = "store", type = "string", dest = "tagg", 
                  help = "Tags selected. 'S' for Sequencing runs, 'Q' for Quality control runs",default="S")
# User_ONYEN (Column 13, M)
selGroup.add_option("-u", "--user", "--USER", "--onyen", "--ONYEN",action = "store", type = "string", dest = "user", 
                  help = "Users selected.\n",default="all")

# Manual override: supply a path to a TSV file to use instead of master config in CDIR
selGroup.add_option("--manual", "--MAN", "--man", "--MANUAL",action = "store", type = "string", dest = "manual", 
                  help = "Path to a manually created config file. If supplied, 'project' not required.",default=None)

parser.add_option_group(selGroup)


# OUTPUT PARAMETERS 
# Parameters for naming the output directory

outGroup = optparse.OptionGroup(parser, "Output Parameters",
                    "These options help customize the output directory. "
                    "Structure is: HOMEDIR/PROJ/rna/NAME/<outputs>. "
                    "These options let you control HOMEDIR and NAME. "
                    "PROJ is generated based on input samples. ")

# Name: override the automatically generated name. Not compatible with Suffix
outGroup.add_option("-n", "--n", "--name",action = "store", type = "string", dest = "name", 
                  help = "Manually override the generated name for output directory.",default=None)
# Suffix: add a string to the end of the automatically generated name.
outGroup.add_option("--suff", "--suffix", "--SUFF", "--SUFFIX",action = "store", type = "string", dest = "suff", 
                  help = "Suffix to append to generated output directory NAME.",default=None)
# Home directory: directory where final output directory, files will be made (by default will build /PROJ/rna/outputdirectory).
outGroup.add_option("--home","--o", "--dir","--HOME","--OUT", "--DIR",action = "store", type = "string", dest = "home", 
                  help = "Path to home directory which will house the output directory in ./PROJ/rna/outputdirectory",default=None)

parser.add_option_group(outGroup)


# RUN PARAMETERS
# Parameters for the run itself; which parts to run, how to handle files

runGroup = optparse.OptionGroup(parser, "Run Parameters",
                    "These options control which parts of the pipeline to run. "
                    "Also control settings like sample merging, temp file deleting, overwriting.")

# Temporary: option to keep the temporary files that are normally deleted
runGroup.add_option("--temp", "--TEMP", "--keep", "--KEEP",action = "store", type = "string", dest = "temp", 
                  help = "Toggle (T/F) whether to keep temporary files (FastQC, trimmed FastQ, etc.)",default="F")
# Stage: pick the parts of the analysis you want to run (function to interpret???)
runGroup.add_option("--stage", "--STAGE", "--S",action = "store", type = "string", dest = "stage", 
                  help = "Select the stages of the code to run. Choose from 'QC' (create multiQC), trim' (create trimmed.fq), quant' (Salmon quant.sf files), 'align' (BAM files), 'signal' (bedgraphs/bigwigs, see --sigout), 'merge' (combine BAMs by Bio Rep for signal step), and 'splice' (MAJIQ).",
                  default="QC,trim,align,quant,signal") 
# Merge: choose what to match for BAM alignment merging (used to create merged signal tracks)
runGroup.add_option("--merge", "--mergeby", "--MERGE", action = "store", type = "string", dest = "merge", 
                  help = "Select which columns to merge when creating sample alignments for combined signal tracks. Please use exact column name.",
                  default="Bio_Rep,Tech_Rep,Seq_Rep") 
# Sigout: control the output of signal, if chosen
runGroup.add_option("--signalout", "--sigout", "--SIGNALOUT", "--SIGOUT",action = "store", type = "string", dest = "sigout", 
                  help = "Select the output signal file type. Choose from 'bigwig' and/or 'bedgraph",default="bigwig") 
# Rerun: option to overwrite existing files, rather than skip the steps that lead to them
runGroup.add_option("--rerun", "--RERUN", action = "store", type = "string", dest = "rerun", 
                  help = "Toggle (T/F) whether to overwrite existing files in output directory.",default="T")
# Mem: Select how much memory to run jobs with (default 32g)
runGroup.add_option("--mem", "--MEM", action = "store", type = "string", dest = "mem", 
                  help = "Memory for all SLURM jobs. (default=32g)",default="32g")

parser.add_option_group(runGroup)


# COMMAND PARAMETERS 
# Parameters passed to the actual commands

cmdGroup = optparse.OptionGroup(parser, "Command Parameters",
                    "These options are passed directly to the programs used. "
                    "Examples: annotation files for Salmon, indices, MAJIQ settings.")

# Annotation: GTF file for tx2gene
cmdGroup.add_option("--ann", "--ANN", "--annotation", "--gtf", "--GTF",action = "store", type = "string", dest = "annt", 
                  help = "Path to transcriptome annotation used to create tx2gene (default: GENCODE v.19)",
                  default="/proj/phanstiel_lab/SHARE/genomes/GENCODE.v19/gencode.v19.annotation.gtf_withproteinids")
# Kmer: kmer size for Salmon index
cmdGroup.add_option("--kmer", "--KMER", "--k",action = "store", type = "string", dest = "kmer", 
                  help = "kmer size; minimum acceptable length for a valid match during Salmon indexing (default: 31)",
                  default="31")
# Transcript: salmon index path
cmdGroup.add_option("--trns", "--TRNS", "--transcript",action = "store", type = "string", dest = "trns", 
                  help = "Path to Salmon-created transcript index (default: GENCODE v.19)",
                  default="/proj/phanstiel_lab/SHARE/genomes/GENCODE.v19/salmon_index")
# HiSat2Idx: HISAT2 index path
cmdGroup.add_option("--hisatidx", "--HSIDX", "--hisat2idx", "--hisat2index",action = "store", type = "string", dest = "hisat2idx", 
                  help = "Path to HISAT2-created index + prefix (default: HG19)",
                  default="/proj/seq/data/HG19_UCSC/Sequence/HISAT2Index/genome")
# MAJIQ config: MAJIQ config file path
cmdGroup.add_option("--majConfig", "--mjqConfig", "--mjqcon", "--majcon",action = "store", type = "string", dest = "mjqcon", 
                  help = "Path to custom MAJIQ configuration file. (Default: created using defaults, run config file)",
                  default=None)
# MAJIQ length: read length for MAJIQ
cmdGroup.add_option("--length", "--readlength", "--readlen", "--len",action = "store", type = "string", dest = "length", 
                  help = "Read length (required for MAJIQ; default = 150)",
                  default="150")
# MAJIQ genome: genome used for MAJIQ
cmdGroup.add_option("--genome",action = "store", type = "string", dest = "genome", 
                  help = "Genome used (required for MAJIQ; default = hg19)",
                  default="hg19")
# MAJIQ strandness: RNA strandness used for MAJIQ
cmdGroup.add_option("--strand", "--strandness",action = "store", type = "string", dest = "strand", 
                  help = "Strandness. Choose from 'forward', 'reverse', and 'None' (required for MAJIQ; default = forward)",
                  default="forward")
# MAJIQ annotation: transcriptome annotation used for MAJIQ, must be GFF3
cmdGroup.add_option("--mjqann", "--MJQANN", "--majiqannotation", "--gff3", "--GFF3",action = "store", type = "string", dest = "gff3", 
                  help = "Path to transcriptome annotation GFF3 file used for MAJIQ (default: hg19)",
                  default="/proj/phanstiel_lab/SHARE/genomes/GENCODE.v19/gencode.v19.annotation.gtf_withproteinids")

parser.add_option_group(cmdGroup)


(options, args) = parser.parse_args()



# ASSIGN TO VARIABLES 
# ===================

# Selection parameters
PROJ = optionLister(options.proj)
CELL = optionLister(options.cell)
GENO = optionLister(options.geno)
COND = optionLister(options.cond)
TIME = optionLister(options.time)
BIOR = optionLister(options.bior)
TECR = optionLister(options.tecr)
SEQR = optionLister(options.seqr)
TAGG = optionLister(options.tagg)
USER = optionLister(options.user)
MANUAL = options.manual

# Output parameters
HOME = options.home
NAME = options.name
SUFF = options.suff

# Run parameters
TEMP  = options.temp
STAGE = optionLister(options.stage)
SIGOUT = optionLister(options.sigout)
MERGE = optionLister(options.merge)
RERUN = options.rerun
MEM = options.mem

# Command parameters
ANNT = options.annt
KMER = options.kmer
TRNS = options.trns
HSIDX = options.hisat2idx
MJQCON = options.mjqcon
LENGTH = options.length
GENOME = options.genome
STRAND = options.strand
MJQANN = options.gff3



# EDIT/INTERPRET VARIABLES 
# ========================

# Interpret alternative TAGG settings
if "all" in TAGG or "A" in TAGG:
    TAGG = ["Q","S"]

# Interpret alternative STAGE settings
if "all" in STAGE:
    STAGE = ["qc", "trim", "quant", "align", "signal", "splice", "merge"]

# Interpret alternative SIGOUT settings
if "all" in SIGOUT:
    SIGOUT = ["bigwig", "bedgraph"]

# Interpret alternative T/F settings (TEMP, RERUN, MERGE)
truthResp = ["T", "t", "true", "True", "TRUE"]
falseResp = ["F", "f", "false", "False", "FALSE"]

if TEMP in truthResp:
    TEMP = True
elif TEMP in falseResp:
    TEMP = False
else:
    parser.error("Please only use true/false tags for the --temp option.")

if RERUN in truthResp:
    RERUN = True
elif RERUN in falseResp:
    RERUN = False
else:
    parser.error("Please only use true/false tags for the --rerun option.")

# Add default suffix if manual config supplied
if SUFF == None and MANUAL is not None:
    SUFF = "MANUAL"



# INPUT ERRORS 
# ============

# ...if no project listed (unless manual config supplied with MANUAL)
if "all" in PROJ and MANUAL == None:
    parser.error('Please select a project.')

# ...if more than one project listed
if len(PROJ) > 1:
    parser.error('Please select only one project.')

# ...if stage chosen is not viable
for item in STAGE:
    if item not in ['QC', 'trim', 'align', 'quant', 'signal', 'splice', 'merge']:
        parser.error("Please choose from 'QC', trim', quant', 'align', 'signal', 'merge', and 'splice' for --stage.")



# PREP FOR SELECTION 
# ==================

# Turn selection parameter inputs into dictionary
selectors = {
    'Project':PROJ,
    'Cell_Type':CELL,
    'Genotype':GENO,
    'Condition':COND,
    'Time':TIME,    
    'Bio_Rep':BIOR,
    'Tech_Rep':TECR,
    'Seq_Rep':SEQR,
    'Tag':TAGG,
    'User_ONYEN':USER
    }

# Limit dictionary to only non-'all' values (parameters to be selected by)
selectors = dict((k,v) for k,v in selectors.items() if 'all' not in v)


    


#------------------------------------------- S E L E C T -------------------------------------------#

# READ IN MASTER CONFIG 
# =====================

# Find current master config file synced from Google Drive, or provided directory if sampled manually
configFile = ""
if MANUAL == None:
    for file in os.listdir(CDIR):
        if file.endswith("rna.tsv"):
            configFile = CDIR + file
else:
    configFile = MANUAL

# ERROR if no master config file found
if os.path.isfile(configFile) == False and MANUAL == None:
    parser.error("No master config file found in" + CDIR + "\nTry downloading locally and re-syncing.")

# ERROR if manual config file not found
if os.path.isfile(configFile) == False and MANUAL is not None:
    parser.error("No config file found in manually defined config location " + configFile)

# Read in master config file
config = pd.read_table(configFile, dtype=str)



# SUBSET CONFIG 
# =============

# Build boolean dataframe
#     For each non-default input, build boolean series for matching config rows, build into dataframe
subsetter = pd.DataFrame(config.Exclude.isnull()) #    Always exclude rows where there is any input in "Exclude" column
for entry in selectors:
    boolKey = config[entry].isin(selectors[entry])
    subsetter[entry] = boolKey

# Find indices where all inputs match
allTrue = subsetter.all(axis='columns') # Creates series of length nrow, T/F for whether all columns are True
subsetter['Selected'] = allTrue

# Subset master config file for selected indices
config = config.loc[subsetter.Selected,:]

# ERROR if no indices where all inputs match
if subsetter.Selected.any() == False and MANUAL is None: # If it is false that there are any True values in "Selected" column
    parser.error("No samples were found matching input parameters. Try checking config file at " + configFile)



# SELECT MERGE PAIRS 
# ==================

if "merge" in STAGE:
	header = list(config)
	# ERROR if merge option is not viable
	for item in MERGE:
		if item not in header:
			parser.error(item + " not found in config file. Please choose from following column names for --merge. " + str(header))
	# Remove unique columns: Unique_ID (0), Name (1) 
	# Remove irrelevant columns for merging: everything after Exclude (12-24)
	header = header[2:11]
	# Remove columns being merged (set by MERGE)
	for item in MERGE:
		header.remove(item)
	# Build sets of samples to be merged into dataframe
	mergeDF = pd.DataFrame()
	for n in range(len(config)):               
		# for each row in config...
		row = config.iloc[n]
		# ...make a merged sample name based on everything that ISN'T merged, plus the "RNA" tag
		rowList = map(str, row[header])
		rowList.insert(1,'RNA')
		mergeName = "_".join(rowList)   
		# ...build a dataframe   
		merger = pd.DataFrame()                
		for column in header:                    #     for each column name that must match to merge...
			query = [row[column]]                #     ...turn the value for given row into a list
			boolKey = config[column].isin(query) #     ...create T/F series for other rows that match that row's entry for that column
			merger[column] = boolKey 
		# distill into one T/F series for rows that match ALL items; extract names from config rows that matched    
		matchNames = pd.DataFrame(config[merger.all(axis='columns')]['Name'].values, columns=[mergeName]) 
		# add extracted names to dataframe
		mergeDF = pd.concat([mergeDF,matchNames], axis=1) 
	# remove duplicate columns    
	mergeDF = mergeDF.T.drop_duplicates().T    
	# Transpose for printing purposes
	mergeDFprint = mergeDF.T
	mergeDFprint.index.name="MERGED NAME"
	mergeDFprint.columns = ["SAMPLE " + s for s in map(str, range(1,mergeDFprint.shape[1]+1))]



# GENERATE NAME 
# =============

# Name of final output directory
if NAME == None:
    NAME = namer(['Project', 'Cell_Type', 'Genotype', 'Condition', 'Time'],['Tag'])
    if SUFF is not None:
        NAME = "_".join([NAME, SUFF])

# Name of directory within homeDir/***** where output directories sit; named after project
PROJNAME = "".join(map(str, set(config["Project"])))



# FEEDBACK 
# ========

# Print out parameters to stdout to check before run

# Parameters used to select from master config
print "\n"
print "SELECTION PARAMETERS"
print "===================="
for key in selectors:
    print str(key) + "\t= " + str(selectors[key])

# Master config used
if MANUAL is not None:
    print "Manual config provided: " + configFile
else:
    print "Master config file used: " + configFile
print "\n"

# Files used
print "FASTQ FILES"
print "==========="
print str(config.shape[0]) + " samples found ..."
print config[['Name','Read1','Read2']].to_string(index=False)
print "\n"

# Merged pairs (if selected)
if "merge" in STAGE:
    print "MERGED PAIRS"
    print "============"
    print "Samples will be merged as follows:"
    print str(mergeDFprint.shape[0]) + " final merged files ..."
    print mergeDFprint.to_string()
    print "\n"

# Run parameters
print "RUN PARAMETERS"
print "=============="    
print "Keep temp files? = " + str(TEMP)
print "Stages run       = " + str(STAGE)
if "signal" in STAGE:
    print "Signal output    = " + str(SIGOUT)
if "merge" in STAGE:    
    print "Merge reps       = " + str(MERGE)
print "Rerun/overwrite? = " + str(RERUN)
print "\n"

# Command parameters
print "COMMAND PARAMETERS"
print "==================\n"    
if "quant" in STAGE:
    print "Salmon/TxImport"    
    print "---------------"
    print "Annotation file  = " + ANNT
    print "Kmer length      = " + KMER
    print "Transcript index = " + TRNS
    print "\n"
if "align" in STAGE:
    print "HISAT2"    
    print "------"
    print "Index            = " + HSIDX
    print "\n"
if "splice" in STAGE:
    print "MAJIQ"    
    print "-----"
    if MJQCON == None:
        print "Config file created in /config based on:"
        print "Read length      = " + LENGTH
        print "Genome           = " + GENOME
        print "Strandness       = " + STRAND
    else:
        print "Config file      = " + MJQCON
    print "\n"

# Output parameters
print "OUTPUT PARAMETERS"
print "================="
print "------------------------------------------------------------------------------"    
if HOME == None:
    HOME = os.getcwd()
    print "*****  Home directory not specified, running from current directory      *****"
print "*****         !!! DOUBLE CHECK OUTPUT DIRECTORY STRUCTURE !!!            *****"
print "------------------------------------------------------------------------------"
print "Home directory   = " + HOME
print "Output directory = " + HOME + "/" + PROJNAME + "/rna/" + NAME
print "\n\n\n"


# Check with user
UserCheck = raw_input("Does this all seem right? (y/n)")

positiveResp = ["y", "Y", "yes", "Yes", "YES", 'hell yeah', 'yaaas']
negativeResp = ["n", "N", "no", "No", "NO", 'naw', 'nah', 'nooo']

if UserCheck in negativeResp:
    parser.error("Please review input options and config file and try again!")
elif UserCheck in positiveResp:
    print "\nGreat, we'll keep going! :)"
else:
    parser.error("Sorry, I only understand yes or no :(")





#------------------------------------------- P R E P -------------------------------------------#

# BUILD BASE OUTPUT DIRECTORY
# ===========================

print "\n\n"
print "Making directories..."
print "\n\n"


# Build up to output directory: HOME/PROJNAME/rna/
procPath = HOME + "/" + PROJNAME + "/rna/" + NAME 

print "Making output directory:"
print "mkdir -p " + procPath
os.system("mkdir -p " + procPath)
print "\n"



# MAKE OTHER DIRECTORIES 
# ======================

# Make other directories, as needed
directories = {
    "fastq":procPath + "/fastq",
    'config':procPath + "/config",
    'debug':procPath + "/debug",
    'scripts':procPath + "/scripts",
    'QC':procPath + "/QC",    
    'align':procPath + "/aligned",
    'signal':procPath + "/signal",
    'quant':procPath + "/quant",
    'splice':procPath + "/splice"
    }

baseDirectories = ['fastq', 'config', 'debug', 'scripts']

print "Making subdirectories:"
for entry in directories:
    if entry in baseDirectories or entry in STAGE:
        os.system("mkdir -v " + directories[entry])



# MAKE RUN INFO LOG FILE  
# ======================  

# Make a log file containing 
#    (1) time run started, 
#    (2) command as entered in terminal, 
#    (3) all settings (prior to interpretations), and 
#    (4) versions of all software

# Create signal module list appropriate for settings
signalDict = {}

if "bedgraph" in SIGOUT:
    signalDict["Bedtools"] = bedVers

if "bigwig" in SIGOUT:
    signalDict["Deeptools"] = deepVers

# Create full module list
modules = {
    'QC':{"FASTQC":FASTQCvers, "MultiQC":multiQCvers},
    'trim':{"Trim Galore!":TrimVers},
    'quant':{"Salmon":SalmonVers, "R":rVers},
    'align':{"HISAT2":hisatVers, "Samtools":samVers},
    'signal':signalDict,    
    'splice':{"MAJIQ":majiqVers},
    'merge':{"Samtools":samVers}
    }

infoOutput = directories['debug'] + "/runinfo_" + NAME + "_" + stamp + ".out"
infoOut = open(infoOutput, "a")
infoOut.write(
    "Date run...\n" + 
    "===========\n" + 
    time.asctime(time.localtime()) + "\n\n" +
    "Command entered...\n" + 
    "==================\n" + 
    " ".join(sys.argv) + "\n\n" +
    "All settings...\n" + 
    "===============\n"
    )
infoOut.write(str(options) + "\n")
infoOut.write(
    "Versions...\n" + 
    "===========\n")

for entry in modules:
    if entry in STAGE:
        for k,v in modules[entry].items():
            infoOut.write(k + " : " + v)

infoOut.close()



# SAVE SUBSET CONFIG FILE 
# =======================

config.to_csv(directories['config'] + "/config_"+ NAME + "_" + stamp + ".tsv", sep="\t", index=False)



# SAVE MERGE PAIR CONFIG FILE 
# ===========================

if "merge" in STAGE:
    mergeDF.T.to_csv(directories['config'] + "/mergeList_"+ NAME + "_" + stamp + ".tsv", header=["SAMPLE " + s for s in map(str, range(1,mergeDF.T.shape[1]+1))], index_label='MERGED NAME')



# MAKE MAJIQ CONFIG FILE 
# ======================

# As accoriding to here: https://majiq.biociphers.org/docs/TUTORIAL%20-%20V1.5.pdf#conf_file
if "splice" in STAGE and MJQCON == None:
    MJQCON = directories['config'] + "/MAJIQconfig_" + NAME + "_" + stamp + ".txt"
    majiqConfig = open(MJQCON, "w")
    majiqConfig.write(
        "[info]\n" + 
        "readlen=" + LENGTH + "\n" +
        "samdir=" + directories['align'] + "\n" + 
        "genome=" + GENOME + "\n" +
        "strandness=" + STRAND + "\n" + 
        "\n" +
        "[experiments]\n" + 
        "Group=" + "_sorted,".join(config['Name']) + "_sorted")
    majiqConfig.close()



# POPULATE /FASTQ W LINKS 
# =======================

# for each row of config, make a link to the fastq files in /fastq, and rename according to the sample name
for n in range(len(config)):

    # Set variables for given row
    nName = config.iloc[n]["Name"]
    nProj = config.iloc[n]["Project"]
    nPath = config.iloc[n]["Sequencing_Directory"]
    nRead1 = config.iloc[n]["Read1"]
    nRead2 = config.iloc[n]["Read2"]
        
    # Path to the raw .fastq.gz in the original sequencing directory # UPDATE AFTER REWORKING CLUSTER???
    #rawRead1 = "/proj/phanstiel_lab/Data/Sequencing/fastq/project/" + nProj + "/rna/" + nPath + nRead1
    #rawRead2 = "/proj/phanstiel_lab/Data/Sequencing/fastq/project/" + nProj + "/rna/" + nPath + nRead2

    # For when whole path is defined in config file... UPDATE LATER????
    rawRead1 = nPath + "/" + nRead1
    rawRead2 = nPath + "/" + nRead2
    
    # Create link from raw file to /fastq
    link1Cmd = "ln -st " + directories['fastq'] + " " + rawRead1
    link2Cmd = "ln -st " + directories['fastq'] + " " + rawRead2
    os.system(link1Cmd)
    os.system(link2Cmd)

    # Rename according to the sample name, appending "R1" or "R2"
    mv1Cmd = "mv " + directories['fastq'] + "/" + nRead1 + " " + directories['fastq'] + "/" +  nName + "_R1.fastq.gz"
    mv2Cmd = "mv " + directories['fastq'] + "/" + nRead2 + " " + directories['fastq'] + "/" +  nName + "_R2.fastq.gz"
    os.system(mv1Cmd)
    os.system(mv2Cmd)





#---------------------------------------- S C R I P T S ----------------------------------------#

# CORE SCRIPTS 
# ============ 

# build empty job list; append with each job's job number at launch
jobList = list()

# Write the .SBATCH Slurm scripts for each sample in config
for n in range(len(config)):

    # Set variables for given row
    nName = config.iloc[n]["Name"]
    nProj = config.iloc[n]["Project"]

    # Make the script
    script = directories['scripts'] + "/CORE_" + nName + "_" + stamp + ".SBATCH"
    file = open(script, "w")

    # Write the commands according to run/command settings
    file.write(
        '#!/bin/bash -l\n' +
        '#SBATCH -J CORE_' + nName + '\n' +
        '#SBATCH -n 8\n' +
        '#SBATCH -N 1\n' +
        '#SBATCH --mem=' + MEM + '\n' +
        '#SBATCH -t 5760\n' +
        '#SBATCH -o ' + directories['debug'] + '/core_' + nName + '_' + stamp + '-%j.log.out\n' + 
        '#SBATCH -e ' + directories['debug'] + '/core_' + nName + '_' + stamp + '-%j.log.err\n' +
        '\n' +
        '#----------------------------------------------------------------------------#\n')

    if "QC" in STAGE:
        file.write(
            '\n\n' +
            '# STAGE: QC (FastQC)\n' +
            '# ==================\n' +
            '\n' +
            'printf "Creating FastQC reports (temp zip/html files)...\\n"\n' +
            'printf "================================================\\n"\n' + 
            '\n' +
            'date\n' +
            'module add fastqc/' + FASTQCvers + '\n' + 
            '\n')

        if RERUN == False and os.path.isfile(directories['QC'] + '/' + nName + '_R1_fastqc.zip') and os.path.isfile(directories['QC'] + '/' + nName + '_R2_fastqc.zip'):
            file.write('printf "FastQC files already exist for sample ' + nName + '\\n"\n' +
                '# COMMAND USED: \n#')

        if RERUN == False and os.path.isfile(directories['QC'] + '/multiqc_report.html'):
            file.write('printf "MultiQC report already exists\\n"\n'+
                '# COMMAND USED: \n#')

        file.write('fastqc -o ' + directories['QC'] + ' ' + directories['fastq']  + '/' + nName + '_R1.fastq.gz ' + directories['fastq'] + '/' + nName + '_R2.fastq.gz\n')

    if "trim" in STAGE:
        file.write(
            '\n\n' +
            '# STAGE: TRIM (Trim Galore!)\n' +
            '# ==========================\n' +
            '\n' +
            'printf "Creating trimmed fastq files (temp files)...\\n"\n' +
            'printf "============================================\\n"\n' + 
            '\n'+
            'date\n' +
            'module add trim_galore/' + TrimVers + '\n' + 
            '\n')

        if RERUN == False and os.path.isfile(directories['fastq'] + '/' + nName + '_1_trimmed.fq.gz') and os.path.isfile(directories['fastq'] + '/' + nName + '_2_trimmed.fq.gz'):
            file.write('printf "Trimmed fastq files already exist for sample ' + nName + '\\n"\n'+
                       '# COMMAND USED: \n#')

        file.write('trim_galore -o ' + directories['fastq'] + ' --paired ' + directories['fastq'] + '/' + nName + '_R1.fastq.gz ' + directories['fastq'] + '/' + nName + '_R2.fastq.gz\n')

        file.write(
            '\n' +
            'printf "Renaming trimmed fastq files...\\n"\n' +
            'printf "===============================\\n"\n' +
            '\n' +
            'date\n' +
            '\n' +
            'mv ' + directories['fastq'] + '/' + nName + '_R1_val_1.* ' + directories['fastq'] + '/' + nName + '_1_trimmed.fq.gz\n' +
            'mv ' + directories['fastq'] + '/' + nName + '_R1*_trimming_report.txt ' + directories['fastq'] + '/' + nName + '_1.fastq_trimming_report.txt\n' +
            'mv ' + directories['fastq'] + '/' + nName + '_R2_val_2.* ' + directories['fastq'] + '/' + nName + '_2_trimmed.fq.gz\n' +
            'mv ' + directories['fastq'] + '/' + nName + '_R2*_trimming_report.txt ' + directories['fastq'] + '/' + nName + '_2.fastq_trimming_report.txt\n')

    if "quant" in STAGE:
        file.write(
            '\n\n' +
            '# STAGE: QUANT (Salmon)\n' +
            '# =====================\n' +
            '\n' +
            'printf "Quantifying transcripts...\\n"\n' +
            'printf "==========================\\n"\n' +
            '\n' + 
            'date\n' +
            'module add salmon/' + SalmonVers + "\n" +
            '\n')

        if RERUN == False and os.path.isfile(directories['quant'] + '/' + nName + '/quant.sf'):
            file.write('printf "Salmon quant.sf files already exist for sample ' + nName + '\\n"\n' +
                '# COMMAND USED: \n#')

        file.write('salmon quant --writeUnmappedNames --threads 1 -i ' + TRNS + ' -l A -1 ' + directories['fastq'] + '/' + nName + '_1_trimmed.fq.gz'+ ' -2 ' + directories['fastq'] + '/' + nName + '_2_trimmed.fq.gz' + ' -o ' + directories['quant'] + '/' + nName + '\n')

    
    if "align" in STAGE:
        file.write(
            '\n\n' +
            '# STAGE: ALIGN (HISAT2, Samtools)\n' +
            '# ===============================\n' +
            '\n' +
            'printf "Creating alignment with HISAT2...\\n"\n' +
            'printf "=================================\\n"\n' +
            '\n' + 
            'date\n' +
            'module add hisat2/' + hisatVers + "\n" +
            'module add samtools/' + samVers + "\n" +
            '\n')

        if RERUN == False and os.path.isfile(directories['align'] + '/' + nName + '_sorted.bam'):
            file.write('printf "Sorted BAM file already exists for sample ' + nName + '\\n"\n'+
                       '# COMMAND USED: \n#')

        file.write('hisat2 -q -x ' + HSIDX + ' -1 ' + directories['fastq'] + '/' + nName + '_1_trimmed.fq.gz -2 ' + directories['fastq'] + '/' + nName + '_2_trimmed.fq.gz | samtools view -u | samtools sort -o ' + directories['align'] + '/' + nName +'_sorted.bam\n')
             
        file.write(
            '\n' + 
            'printf "Creating SAM stat files with Samtools...\\n"\n' +
            'printf "========================================\\n"\n' +
            '\n' +
            'date\n' +
            '\n')

        if RERUN == False and os.path.isfile(directories['align'] + '/' + nName + '_stats.txt'):
            file.write('printf "SAM stat file already exists for sample ' + nName + '\\n"\n'+
                   '# COMMAND USED: \n#')

        file.write('samtools flagstat ' + directories['align'] + '/' + nName + '_sorted.bam > ' + directories['align'] + '/' + nName + '_stats.txt\n')

        file.write(
            '\n' +
            'printf "Indexing BAM files with Samtools...\\n"\n' +
            'printf "===================================\\n"\n' +
            '\n'
            'date\n' +
            '\n')

        if RERUN == False and os.path.isfile(directories['align'] + '/' + nName + '_sorted.bam.bai'):
            file.write('printf "BAM indices already conducted for sample ' + nName + '\\n"\n'+
                       '# COMMAND USED: \n#')

        file.write('samtools index ' + directories['align'] + '/' + nName + '_sorted.bam \n')

    if TEMP == False:
        file.write(
            '\n\n' +
            '# REMOVING TEMP FILES\n' +
            '# ===================\n' +
            '\n' +
            'printf "Deleting temp trimmed.fq files...\\n"\n' +
            'printf "=================================\\n"\n' + 
            '\n' +
            'date\n' +
            'rm ' + directories['fastq'] + '/' + nName + '_1_trimmed.fq.gz \n' +
            'rm ' + directories['fastq'] + '/' + nName + '_2_trimmed.fq.gz \n')

    if "signal" in STAGE:
        file.write(
            '\n\n' +
            '# STAGE: SIGNAL (deeptools and/or bedtools) - unmerged\n' +
            '# ====================================================\n')

        if "bedgraph" in SIGOUT:
            file.write(
                '\n' +
                'printf "Creating bedgraphs with bedtools...\\n"\n' +
                'printf "===================================\\n"\n' +
                '\n' +
                'date\n' +
                'module add bedtools/' + bedVers + "\n" +
                '\n')

            if RERUN == False and os.path.isfile(directories['signal'] + '/' + nName + '.bedgraph'):
                file.write('printf "Bedgraph file already exists for sample ' + nName + '\\n"\n'+
                       '# COMMAND USED: \n#')

            file.write('bedtools genomecov -bga -ibam ' + directories['align'] + '/' + nName + '_sorted.bam > ' + directories['signal'] + '/' + nName + '.bedgraph\n')

        if "bigwig" in SIGOUT:
            file.write(
                '\n' +
                'printf "Creating bigwigs with deeptools...\\n"\n' +
                'printf "==================================\\n"\n' +
                '\n' +
                'date\n' +
                'module add deeptools/' + deepVers + "\n" +
                '\n')

            if RERUN == False and os.path.isfile(directories['signal'] + '/' + nName + '.bw'):
                file.write('printf "Bigwig file already exists for sample ' + nName + '\\n"\n'+
                       '# COMMAND USED: \n#')

            file.write('bamCoverage -b ' + directories['align'] + '/' + nName + '_sorted.bam -o ' + directories['signal'] + '/' + nName + '.bw\n' +
              'bamCoverage --filterRNAstrand forward -b ' + directories['align'] + '/' + nName + '_sorted.bam -o ' + directories['signal'] + '/' + nName + '_fwd.bw\n' + 
              'bamCoverage --filterRNAstrand reverse -b ' + directories['align'] + '/' + nName + '_sorted.bam -o ' + directories['signal'] + '/' + nName + '_rev.bw\n')

    file.close()

    coreCmd = 'sbatch < ' + script
    print "Submitting job CORE %s with command: %s" % (nName, coreCmd)

    status, jobStmt = commands.getstatusoutput(coreCmd)
    jobnum = jobStmt.split(" ")[-1]
    if (status == 0):
        print "Job number for CORE %s is %s" % (nName, jobnum) 
    else:
        print "*** ERROR SUBMITTING JOB ***"

    jobList.append(jobnum)
coreJobs = ":".join(jobList)



# MERGE SCRIPTS 
# ============ 

# Separate script for if "merge" in STAGE (merge):
#     Merge w/ Samtools
#    Signal tracks w/ bed/deeptools, if "signal"

# build empty job list; append with each job's job number at launch
jobList = list()

# Write the .SBATCH Slurm scripts for each sample in config
if "merge" in STAGE and mergeDF.shape[0] > 1:

    for mergeName in mergeDF:

        # Set variables for given mergeName 
        sampleList = mergeDF[mergeName]
        sampleList = sampleList[sampleList.notna()]

        if len(sampleList) > 1:

            # Make the script
            script = directories['scripts'] + "/MERGE_" + mergeName + "_" + stamp + ".SBATCH"
            file = open(script, "w")

            # Write the commands according to run/command settings
            file.write(
                '#!/bin/bash -l\n' +
                '#SBATCH -J MERGE_' + mergeName + '\n' +
                '#SBATCH -n 8\n' +
                '#SBATCH -N 1\n' +
                '#SBATCH --mem=' + MEM + '\n' +
                '#SBATCH -t 5760\n' +
                '#SBATCH -o ' + directories['debug'] + '/merge_' + mergeName + '_' + stamp + '-%j.log.out\n' + 
                '#SBATCH -e ' + directories['debug'] + '/merge_' + mergeName + '_' + stamp + '-%j.log.err\n' +
                '\n' +
                '#----------------------------------------------------------------------------#\n')
            
            file.write(
                '\n\n' +
                '# STAGE: MERGE (Samtools)\n' +
                '# =======================\n' +
                '\n' +
                'printf "Merging BAM files with Samtools (temp merge.BAM files)...\\n"\n' +
                'printf "=========================================================\\n"\n' +
                '\n' +
                'date\n' +
                'module add samtools/' + samVers + "\n" +
                '\n')

            if RERUN == False and os.path.isfile(directories['align'] + '/MERGE_' + mergeName + '.bam'):
                file.write('printf "Merged BAM file already exists for samples: ' + ", ".join(sampleList) + '\\n"\n'+
                             '# COMMAND USED: \n#')

            mergeStrList = directories['align'] + '/' + sampleList + '_sorted.bam'
            file.write('samtools merge ' + directories['align'] + '/MERGE_' + mergeName + '.bam ' + " ".join(mergeStrList) + '\n')
            
            file.write(
                '\n' +
                'printf "Indexing BAM files with Samtools...\\n"\n' +
                'printf "===================================\\n"\n' +
                '\n'
                'date\n' +
                '\n')

            if RERUN == False and os.path.isfile(directories['align'] + '/MERGE_' + mergeName + '.bam.bai'):
                file.write('printf "BAM indices already conducted for samples: ' + ", ".join(sampleList) + '\\n"\n'+
                           '# COMMAND USED: \n#')
    
            file.write('samtools index ' + directories['align'] + '/MERGE_' + mergeName + '.bam \n')

            if "signal" in STAGE:
                file.write(
                    '\n\n' +
                    '# STAGE: SIGNAL (deeptools and/or bedtools) - merged\n' +
                    '# ==================================================\n')

                if "bedgraph" in SIGOUT:
                    file.write(
                    '\n' +
                    'printf "Creating bedgraphs with bedtools...\\n"\n' +
                    'printf "===================================\\n"\n' +
                    '\n' +
                    'date\n' +
                    'module add bedtools/' + bedVers + "\n" +
                    '\n')

                    if RERUN == False and os.path.isfile(directories['signal'] + '/MERGE_' + mergeName + '.bedgraph'):
                        file.write('printf "Merged bedgraph file already exists for samples: ' + ", ".join(sampleList) + '\\n"\n'+
                           '# COMMAND USED: \n#')

                    file.write('bedtools genomecov -bga -ibam ' + directories['align'] + '/MERGE_' + mergeName + '.bam > ' + directories['signal'] + '/MERGE_' + mergeName + '.bedgraph\n')

                if "bigwig" in SIGOUT:
                    file.write(
                    '\n' +
                    'printf "Creating bigwigs with deeptools...\\n"\n' +
                    'printf "==================================\\n"\n' + 
                    '\n' +
                    'date\n' +
                    'module add deeptools/' + deepVers + "\n")

                    if RERUN == False and os.path.isfile(directories['signal'] + '/MERGE_' + mergeName + '.bw'):
                        file.write('printf "Merged bigwig file already exists for samples ' + ", ".join(sampleList) + '\\n"\n'+
                               '# COMMAND USED: \n#')

                    file.write('bamCoverage -b ' + directories['align'] + '/MERGE_' + mergeName + '.bam -o ' + directories['signal'] + '/MERGE_' + mergeName + '.bw\n' +
                      'bamCoverage --filterRNAstrand forward -b ' + directories['align'] + '/MERGE_' + mergeName + '.bam -o ' + directories['signal'] + '/MERGE_' + mergeName + '_fwd.bw\n' + 
                      'bamCoverage --filterRNAstrand reverse -b ' + directories['align'] + '/MERGE_' + mergeName + '.bam -o ' + directories['signal'] + '/MERGE_' + mergeName + '_rev.bw\n')

            if TEMP == False:
                file.write(
                    '\n\n' +
                    '# REMOVING TEMP FILES\n' +
                    '# ===================\n' +
                    '\n' +
                    'printf "Deleting temp merged BAM files...\\n"\n' +
                    'printf "=================================\\n"\n' + 
                    '\n' +
                    'date\n' +
                    'rm ' + directories['align'] + '/MERGE_' + mergeName + '.bam\n' +
                    'rm ' + directories['align'] + '/MERGE_' + mergeName + '.bam.bai\n')

            file.close()

            mergeCmd = 'sbatch --dependency=afterok:' + coreJobs + ' ' + script
            print "Submitting job MERGE %s with command: %s" % (mergeName, mergeCmd)

            status, jobStmt = commands.getstatusoutput(mergeCmd)
            jobnum = jobStmt.split(" ")[-1]
            if (status == 0):
                print "Job number for MERGE %s is %s" % (mergeName, jobnum)
            else:
                print "*** ERROR SUBMITTING JOB ***: status: %s message: %s " % (status, jobStmt)

            jobList.append(jobnum)
        mergeJobs = ":".join(jobList)



# FINAL SCRIPTS 
# =============

# Separate script for MultiQC if "QC", tximport if "quant" (final):
if "QC" in STAGE or "quant" in STAGE:

    # Make the script
    script = directories['scripts'] + "/FINAL_" + NAME + "_" + stamp + ".SBATCH"
    file = open(script, "w")

    # Write the commands according to run/command settings
    file.write(
        '#!/bin/bash -l\n' +
        '#SBATCH -J FINAL_' + NAME + '\n' +
        '#SBATCH -n 8\n' +
        '#SBATCH -N 1\n' +
        '#SBATCH --mem=' + MEM + '\n' +
        '#SBATCH -t 5760\n' +
        '#SBATCH -o ' + directories['debug'] + '/final_' + NAME + '_' + stamp + '-%j.log.out\n' + 
        '#SBATCH -e ' + directories['debug'] + '/final_' + NAME + '_' + stamp + '-%j.log.err\n' +
        '\n' +
        '#----------------------------------------------------------------------------#\n')

    file.write(
        '\n' +
        '#------------------------ Run Code --------------------------#\n')

    if "QC" in STAGE:
        file.write(
        '\n\n' +
        '# STAGE: QC SUMMARY (MultiQC)\n' +
        '# ===========================\n' +
        '\n' +
        'printf "Generating MultiQC report...\\n"\n' +
        'printf "============================\\n"\n' +
        '\n' +
        'date\n' +
        'module add multiqc/' + multiQCvers + "\n" +
        '\n' +
        'multiqc -f ' + procPath + '/* -o ' + directories['QC'] + '\n' +
        'mv ' + directories['QC'] + '/multiqc_report.html ' + directories['QC'] + '/' + NAME + '_multiqc_report.html\n')

        if TEMP == False:
            file.write(
                '\n\n' +
                '# REMOVING TEMP FILES\n' +
                '# ===================\n' +
                '\n' +
                'printf "Deleting temp FastQC files...\\n"\n' +
                'printf "=============================\\n"\n' +
                '\n' +
                'date\n' +
                'rm ' + directories['QC'] + '/*fastqc.zip \n' +
                'rm ' + directories['QC'] + '/*fastqc.html \n' +
                'rm ' + directories['fastq'] + '/*fastq_trimming_report.txt \n')

    if "quant" in STAGE:
        file.write(
        '\n\n' +
        '# STAGE: QUANT SUMMARY (txImport)\n' +
        '# ===============================\n' +
        '\n' +
        'printf "Generating tximport object...\\n"\n' +
        'printf "=============================\\n"\n' +
        '\n' +
        'date\n' +
        'module add r/' + rVers + "\n" +
        'Rscript /proj/phanstiel_lab/software/RNApipe/txImporter.R ' + directories['config'] + "/config_"+ NAME + "_" + stamp + ".tsv "  + ANNT + " " + directories['quant'] + " " + NAME)

    if "splice" in STAGE:
        file.write(
            '\n\n' +
            '# STAGE: SPLICE (MAJIQ)\n' +
            '# =====================\n' +
            '\n' +
            'printf "Running MAJIQ Builder (creating LSV)...\\n"\n' +
            'printf "=======================================\\n"\n' + 
            '\n' +
            'date\n' +
            'module add majiq/' + majiqVers + "\n" +
            '\n')

        if RERUN == False and os.path.isfile(directories['splice'] + '/*_sorted.majiq') and os.path.isfile(directories['splice'] + '/*_sortedsplicegraph.sql'):
            file.write('printf "MAJIQ files already exist for project.\\n"\n'+
                   '# COMMAND USED: \n#')

        file.write('majiq build ' + MJQANN + ' -c ' + MJQCON + ' -j 8 -o ' + directories['splice'] + '\n')

    file.close()

    if "merge" in STAGE and mergeDF.shape[0] > 1:    
        finalCmd = 'sbatch --dependency=afterok:' + mergeJobs + ' ' + script
    else:
        finalCmd = 'sbatch --dependency=afterok:' + coreJobs + ' ' + script
    print "Submitting job FINAL %s with command: %s" % (PROJNAME, finalCmd)

    status, jobStmt = commands.getstatusoutput(finalCmd)
    jobnum = jobStmt.split(" ")[-1]
    if (status == 0):
        print "Job number is %s" % jobnum
    else:
        print "*** ERROR SUBMITTING JOB *** status: %s message: %s " % (status, jobStmt)
