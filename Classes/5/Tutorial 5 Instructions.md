# Tutorial 5
Today we will explore differential gene expression analysis and also biologically motivated analysis through building gene lists. This tutorial highlights two previous examples and then asks you to think about a question answerable with single cell expression data and practice building a gene list to answer it.

We'll turn back to the HLCA dataset, but this time use the SmartSeq2 expression profiles from all 3 subjects.

## Login to FarmShare
```bash
ssh <SUNetID>@rice.stanford.edu
```
**Note**: Sherlock users, please use \<SUNetID\>@login.sherlock.stanford.edu

### Switch to the persistent terminal window
The resources we requested on Wednessday should be ready to go. To check specifics issue "squeue" command:
```bash
squeue -u <SUNetID>
```
You should be able to see that SLURM resource manager has assigned one of the compute nodes "wheat\<xx\>" for your use. If the resources have been allocated
, then proceed. If not, please let us know!
```bash
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           2107056    normal     srun   <SUNetID>  R      26:50      1 wheat<xx>
```
Login to wheat node assigned to you
```bash
ssh wheat<xx>
```
You can tell if your terminal has switched from something like "(base) \<SUNetID\>\@**rice**\<XX\>:\~\$" to "(base) \<SUNetID\>\@**wheat**\<XX\>:\~\$"

**Note:** Sherlock users will be connecting to a compute node instead of wheat.

Just like last week we will start a persistant session with "tmux". In case we lose connection, we can login again to the same assigned wheat node as above
and reattach to the same persistantly running "tmux" session. "tmux" session is specific to a given compute node not cluster-wide, therefore, the "tmux" ses
sion you just started is specific to you and is running only on the wheat\<xx\> (or rice\<xx\> node as we did in the first tutorial) where you started it
```bash
tmux
```
You should now see the green tmux bar along the bottom of your terminal.

### Login to jupyter in your browser, open a Terminal window (File > New... > Terminal), and activate the single cell environment
```bash
conda activate singlecell
```

### Pull the GitHub repository
```bash
cd $HOME/BIOC281
git pull
```

**Note:** On sherlock run
```bash
cd $GROUP_HOME/BIOC281
git pull
```

### Navigate to the current tutorial folder
```bash
cd $HOME/BIOC281/Classes/5
```

**Note:** On sherlock run
```bash
cd $GROUP_HOME/BIOC281/Classes/5
```

## Download the HLCA FACS data and metadata from all three patients as CSVs

```bash
wget --quiet http://hsc.stanford.edu/resources/DataObjects/krasnow_hlca_facs_counts.csv
wget --quiet http://hsc.stanford.edu/resources/DataObjects/krasnow_hlca_facs_metadata.csv
```

### Create folder to download gene data
```bash
mkdir dbs
cd dbs
```

### Download a list of transcription factors from AnimalTFDB
```bash
wget --quiet http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Homo_sapiens_TF -O AnimalTFDB.tsv
```

### Download a list of gene associated with mendelian diseases from OMIM
Visit https://www.omim.org/downloads and complete the download request form. Indicate "Use in a class on single cell expression analysis" as your description for using OMIM. Check your email and copy the link for genemap2.txt
```bash
wget --quiet <URL>/mgenemap2.txt -O OMIM.tsv
```

### Download a list of enzymes from ExPASy
```bash
wget --quiet ftp://ftp.expasy.org/databases/enzyme/enzyme.dat -O ExPASy.txt
```

### Download a list of receptors from Guide to Pharmacology
```bash
wget --quiet https://www.guidetopharmacology.org/DATA/targets_and_families.tsv -O GtP_receptors.tsv
```

### Download a list of ligands from Guide to Pharmacology
```bash
wget --quiet https://www.guidetopharmacology.org/DATA/ligands.tsv -O GtP_ligands.tsv
wget --quiet https://www.guidetopharmacology.org/DATA/GtP_to_HGNC_mapping.tsv
```

## Start jupyter
```bash
jupyter lab --no-browser
```
**Note:** Sherlock users will need to specify a port between 49152-65335 by adding --port=<Port#> as before.

## Create an SSH tunnel (on your system)
```bash
ssh -N -f -L <Port>:localhost:<Port> <SUNetID>@wheat<XX>.stanford.edu
```

**Note:** Sherlock users will need the more complex tunneling command they have used before. The last \<shXX-XXnXX\> part comes from the node they are on, it is the part of your command prompt bolded here: (singlecell) [\<SUNetID\>@**sh03-01n52** ~]$
```
ssh -L <Port#>:localhost:<Port#> <SUNetID>@login.sherlock.stanford.edu ssh -L <Port#>:localhost:<Port#> -N <shXX-XXnXX>
```

## Login to jupyter in your browser, open a Terminal window, and activate the single cell environment

## Open and complete BIOC281/Classes/5/Tutorial 5.ipypnb
When finished, please save the outputs of Tutorial 5.ipynb (File > Export Notebook As... > Export Notebook to HTML) when you have completed it and upload it to Canvas.

### Relinquish the current wheat allocation
After running the command below, your terminal should switch back to rice from wheat.
```bash
exit
```

### Run salloc
```bash
salloc --ntasks-per-node=1 --cpus-per-task=4 --mem=30G --time=0-5:00:00 --begin="13:30:00 10/30/20" --qos=interactive srun --pty bash -i -l
```
**Note:** Sherlock users should use --quit=normal as Sherlock does not have an interactive quality of service like FarmShare.

You can now detach from the "tmux" with control+b and then press "d". The green bar on the bottom should disappear. Note down the rice node (**rice\<xx\>**) where your "tmux" session is running (in case we need to get back to it) and the session ID (usually just a single session==0). After that, you can safely close your terminal window. See you next class!
