# MDMcleaner
## A pipeline for the assessment, classification and refinement of microbial dark matter SAGs and MAGs
MDMcleaner is a reference-DB contamination aware pipeline for reliable contig classification of metagenome assembled (MAG) and single-cell amplified (SAG) genomes.
It is based on the GTDB taxonomic system and uses GTDB representative genomes, as well as SILVA SSU and LSU and RefSeq eukaryotic and viral datasets as references. Classification is based on a "least common ancestor" (LCA) approach, that is implemented in a way that can recognize potential contaminants not only in the analyzed genome, but also in the underlying reference datasets. Furthermore each contig is classified only up to taxlevels that are actually supported by the corresponding alignment identities, thereby avoiding overclassification for organisms that are underrepresented in the reference database.

## Dependencies:
- python 3.7+
- biopython v.1.78+
- wget v.1.19+ (downloading reference datasets)
- ncbi-blast v.2.10.1+ (aligning nucleotide sequences)
- diamond v.2.0.6+ (aligning amino acid sequences)
- hmmer v.3.3.1+ (detecting conserved marker genes)
- barrnap v.0.9+ (detecting ribosomal RNA genes)
- aragorn v 1.2.38+ (detecting tRNA genes)
- prodigal v 2.6.3 (ORF/CDS-prediction)

This repository is hosted at [github]( https://github.com/KIT-IBG-5/mdmcleaner) and mirrored at [gitlab](https://git.scc.kit.edu/ww5070/mdmcleaner).

## Installation:
The MDMcleaner is now installable via **pip** (without dependencies), full availability via **Bioconda** recipe will follow shortly!
Until then, a convinient way to install all needed dependencies is to create a conda environment using [bioconda](http://www.ddocent.com//bioconda/) and the following command:
```
conda create -n mdmcleaner_env python=3.10 biopython blast=2.12.0 diamond hmmer barrnap aragorn prodigal
```
you can then always activate this environment with the command ```conda activate mdmcleaner.env```

After activating that environment you can then install mdmcleaner into it by using the command
```pip install mdmcleaner```

## Configuration
Several options can be passed directly as commandline arguments (see usage below), but basic settings, such as database location, should be provided in the form of ```mdmcleaner.config``` config files. The pipeline distinguishes between global (system/environment-wide settings) and local (individual) config files.

The hierarchy is as follows:
 - command line arguments override global and local setting
 - local settings override global settings
 - global settings are set in the ```mdmcleaner.config``` file within the mdmcleaner ```lib``` folder (best accessed via ```mdmcleaner.py set_configs -s global```

the settings that can be specified/adjusted in the config files are:
- location of blastn binaries (default is simply "blastn", which assumes it is present in PATH)
- location of blastp binaries (default is simply "blastp", which assumes it is present in PATH)
- location of barrnap binaries (default is simply "barrnap", which assumes it is present in PATH)
- location of aragorn binaries (default is simply "blastp", which assumes it is present in PATH)

To create a local config file in the current working directory, simply use ```mdmcleaner set_configs -s local [SETTING_ARGUMENTS]```
This file can be moved and copied and will be automatically recognized if present in the current working directory when running MDMcleaner. Alternative the path to a local config file can be passed to mdmcleaner via the "-c" argument of the "clean" and "makedb" workflows.

## overview of MDMcleaner commands
A list of mdmcleaner commands is returned when invoking the help function of MCMcleaner as follows: ```mdmcleaner -h```. Each command has it's own help function that can be invoked with ```mdmcleaner <COMMAND> -h```. The available commands are:
 - **set_configs** can be used to change global or local settings. Will modify or create 'mdmcleaner.config'-files
 - **show_configs** lists the currently applicable MDMcleaner settings/configurations
 - **makedb** downloads and processes reference data into a MDMclenaner reference database. May have a LONG run-time but can be aborted and resumed
 - **clean** the major MDMcleaner workflow for assessing and filtering genome contamination
 - **get_markers** an accessory command for extracting marker gene sequences from input genomes
 - **completeness** an accessory command for "quick-and-dirty" assessment of bin completeness based on universally required types or tRNAs
 - **refdb_contams** EXPERIMENTAL: evaluates refDBambiguity overviewfiles and adds obvious refDB contaminations to the blacklist
 - **acc2taxpath** Get full taxonomic path associated with a specific input accession. Currently only works for MDMcleaner/GTDB accessions, but support for NCBI accession-numbers will follow soon
 - **check_dependencies** Just check if all dependencies are being met
 - **version** show version info and quit 

#### usage of ```mdmcleaner set_configs```:
```
usage: mdmcleaner set_configs [-h] [-s {local,global}] [--blastp BLASTP] [--blastn BLASTN] [--diamond DIAMOND] [--barrnap BARRNAP] [--hmmsearch HMMSEARCH] [--aragorn ARAGORN] [--db_basedir DB_BASEDIR]
                                 [--threads THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -s {local,global}, --scope {local,global}
                        change settings in local or global config file. 'global' likely require admin privileges. 'local' will modify or create a mdmcleaner.config file in the current working directory.
                        default = 'local'
  --blastp BLASTP       path to blastp binaries (if not in PATH)
  --blastn BLASTN       path to blastn binaries (if not in PATH)
  --diamond DIAMOND     path to diamond binaries (if not in PATH)
  --barrnap BARRNAP     path to barrnap binaries (if not in PATH)
  --hmmsearch HMMSEARCH
                        path to hmmsearch binaries (if not in PATH)
  --aragorn ARAGORN     path to aragorn binaries (if not in PATH)
  --db_basedir DB_BASEDIR
                        path to basedirectory for reference database
  --threads THREADS     threads to use by default
```
## Downloading and creating the reference database
The reference database is created by running ```mdmcleaner.py makedb```. This downloads the most recent datasets from GTDB, RefSeq and Silva and processes them into the format used by MDMcleaner, ensuring that users can access the newest reference data independently of our ability to keep our database current. The download and processing may take a long time (>13h at 100 Mbit/s), but can be resumed from the last checkpoint, simply by running it again if aborted.

Alternatively you can use the ready made database used during the work on the mdmcleaner publication. This has been deposited at zenodo for reproducability reasons. You find and download it at [this link](https://zenodo.org/record/5698995#.YkylqjVCRhE), or use the ```--get_pub_data``` argument of ```mdmcleaner.py makedb```

**Remember to add the location of the database directory to a config file!** Preferrably to the global configuration file by running ```mdmcleaner.py set_configs --db_basedir <path/to/databasefolder>``` 
If you specify the target database directory in a config file, you do not need to use the ```-o``` option of makedb (instead the target folder will be automatically read from the configs file)

#### usage of ```mdmcleaner makedb```:
```
usage: mdmcleaner makedb [-h] [-o OUTDIR] [-c CONFIGFILE] [--verbose] [--quiet]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                        target base directory for reference-data. may not be the current working directory. Needs >100GB space! Default = './db/gtdb'
  -c CONFIGFILE, --config CONFIGFILE
                        provide a local config file with the target location to store database-files. default: looks for config files named 'mdmcleaner.config' in current working directory. settings in
                        the local config file will override settings in the global config file '/home/ww5070/temp_binsaga/sagabin_refiner/lib/mdmcleaner.config'
  --get_pub_data        get the ready made reference dataset used during the work on the mdmcleaner publication (WARNING: is outdated and only provided for reproducability and/or testing purposes)
  --verbose             verbose output (download progress etc)
  --quiet               quiet mode (suppress any status messages except Errors and Warnings)
```

## Processing MAGs and SAGs
To process MAGs and SAGs, simply use the "clean" function of MDMcleaner. If you have multiple genomes, specify them all at once instead of calling MDMcleaner individually for each genome. This apeeds up processing time by allowing to reuse of database objects without needing to reinitialize them for every query genome.
For safety reasons the number of threads used is set to 1 by default. However many steps of the MDMcleaner pipeline can profit from using multiple threads. Remember to specify as many threads as you can safely use on your system.

For each genome, intermediary results are stored in a seperate subfolder in "mdmcleaner_results". Delete them only when you are finished, as these allow resuming or rerunning analyses without having to redo all analyses from scratch...

By default, mdmcleaner will try to resolve any potential reference-database ambiguities it encounters. This will increase run time and may not always be successful. If speed is of essence, this can be skipped by using the ```--fast_run``` argument

If possible, we recommend to use at least 8 threads (```-t 8```) for mdmcleaner runs. With 8 threads it will take about 20-30 minutes per complete & average-sized bacterial genome

#### usage of ```mdmcleaner clean```:
```
usage: mdmcleaner clean [-h] -i INPUT_FASTAS [INPUT_FASTAS ...] [-o OUTPUT_FOLDER] [-v] [-c CONFIGFILE] [-t THREADS] [-f] [--overview_files_basename OVERVIEW_BASENAME] [-I IGNORELISTFILE]
                           [--no_filterfasta]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FASTAS [INPUT_FASTAS ...], --input_fastas INPUT_FASTAS [INPUT_FASTAS ...]
                        input fastas of genomes and/or bins
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        output-folder for MDMcleaner results. Default = 'mdmcleaner_output'
  -v, --version         show program's version number and exit
  -c CONFIGFILE, --config CONFIGFILE
                        provide a local config file with basic settings (such as the location of database-files). default: looks for config files named 'mdmcleaner.config' in current working directory.
                        settings in the local config file will override settings in the global config file '/home/ww5070/temp_binsaga/sagabin_refiner/lib/mdmcleaner.config'
  -t THREADS, --threads THREADS
                        Number of threads to use. Can also be set in the mdmcleaner.config file
  -f, --force           Force reclassification of pre-existing blast-results
  --overview_files_basename OVERVIEW_BASENAME
                        basename for overviewfiles (default="overview"
  -b BLACKLISTFILE, --blacklistfile BLACKLISTFILE
                        File listing reference-DB sequence-names that should be ignored during blast-analyses (e.g. known refDB-contaminations...
  --no_filterfasta      Do not write filtered contigs to final output fastas (Default = False)
  --fast_run            skips evaluation of reference-database-ambiguities. Runs will be faster but ambiguities will not be resolved on the go!
```

## OUTPUT FILES
Running the ```clean``` option of MDMcleaner.py will create a mdmcleaner_results folder in the current orking directory. Individual results for each genome are saved in individual subfolders of "mdmcleaner_results". This includes a detailed report on the individual contig classifications, a input-table for vizualization with KRONA and (optionally) filtered contig fastas, divided into the four categories "keep", "evaluate_low", "evaluate_high" and "delete":
#### individual output files per genome:
- fullcontiginfos_beforecleanup.tsv --> a tab seperated table with details for each contig
- \<genome-name\>_rRNA_lsu_rRNA.fasta --> 23s rRNA genes
- \<genome-name\>_rRNA_ssu_rRNA.fasta --> 16s rRNA genes
- \<genome-name\>_rRNA_tsu_rRNA.fasta --> 5s rRNA genes
- \<genome-name\>_tRNAs.fasta.gz --> tRNA genes
- \<genome-name\>_totalprots.faa --> total protein sequences
- \<genome-name\>_keep.fasta --> trusted contigs that are safe to submit
- \<genome-name\>_evaluate_low.fasta --> contigs that yielded some kind of refDB ambiguity and __may__ profit from some individual re-evaluation (cross-blast against RefSeq etc)
- \<genome-name\>_evaluate_high.fasta --> contigs that yielded indications of refDB contamiations and should not be submitted before detailed re-evaluation
- \<genome-name\>_delete.fasta --> untrusted contigs that should definitively not be included in the genome submission
- various intermediary files and progress markers...

Additionally, general overview files are written to the current working directory: "overview_allbeforecleanup.tsv", "overview_refdb_ambiguities.tsv" and "overview_errorlist.txt"
#### general overview files
- overview_allbeforecleanup.tsv --> tab seperated table listing majority classification and general metrics per analyzed bin
- overview_refdb_ambiguities.tsv --> tab seperated table listing detected refDB ambiguities together with evidence information and preliminary MDMcleaner assessments
- overview_errorlist.txt --> list of genomes that yielded errors during MDMcleaner assessments and may need to be rerun

## How to cite:
look out for our publication which has been submitted to Nucleic acids research and will hopefully be available soon.

Vollmers et al. ***How clear is our current view on microbial dark-matter? (Re-)assessing public MAG & SAG-datasets with “MDMcleaner”*** (submitted)

## For help or feedback:
please write an issue on the [github repository](https://github.com/KIT-IBG-5/mdmcleaner/issues) or send us an [email](mailto:ibg5-support@lists.kit.edu?subject=[MDMcleaner%20support])
