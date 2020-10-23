# Tutorial 3
Today we will be working through two notebooks, one written in R and the other in python, that use the Seurat (https://satijalab.org/seurat/) and scanpy (https://scanpy.readthedocs.io/en/stable/) packages, respectively, to cluster part of the human lung cell atlas single cell dataset (https://hlca.ds.czbiohub.org). The notebooks are designed to build intuition for how changing various paramters affects clustering results and illustrate the utility of iterative subclustering. In addition to the standard normalization, the notebooks also demo some of the more complex procedures we discussed during the lecture as well as methods for integrating datasets constructed with different technologies.

To begin the tutorials, we need to request dedicated resources on Stanford FarmShare, launch a jupyer notebbook, then download the data from Synapse.

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

## Activate the single cell environment
```bash
conda activate singlecell
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
```bash
conda activate singlecell
```

Install a missing dependency, louvain, for today's scanpy notebook and resolve a package conflict created by a new version of igraph (released between making the tutorials and the start of the class).
```bash
pip install louvain
conda install -c conda-forge python-igraph=0.8.2 seaborn=0.10
```

## Pull the GitHub repository
```bash
git pull
```

## Navigate to the current tutorial folder
```bash
cd ~/BIOC281/Classes/3
```
**Note:** On sherlock run
```bash
cd $GROUP_HOME/BIOC281/Classes/3
```

## Download Human Lung Cell Atlas data and metadata from hsc.stanford.edu (Travaglini et al (2020) _Nature_)
This dataset contains ~75,000 cells captured using 10x (90%) and SmartSeq2 (10%) from the lung of 3 human patients.

```bash
wget http://hsc.stanford.edu/resources/DataObjects/krasnow_hlca_10x_raw.rds

wget http://hsc.stanford.edu/resources/DataObjects/krasnow_hlca_10x_raw.h5ad

wget http://hsc.stanford.edu/resources/DataObjects/krasnow_hlca_facs_raw.rds

wget http://hsc.stanford.edu/resources/DataObjects/krasnow_hlca_facs_raw.h5ad
```

## Open and complete BIOC281/Classes/3/Tutorial 3 - Seurat.ipypnb
When finished, please save the outputs of Tutorial 3 - Seurat.ipynb (File > Export Notebook As... > Export Notebook to HTML) when you have completed it and upload it to Canvas.

## Open and complete BIOC281/Classes/3/Tutorial 3 - scanpy.ipynb
When finished, please save the outputs of Tutorial 3 - scanpy.ipynb (File > Export Notebook As... > Export Notebook to HTML) when you have completed it and upload it to Canvas.

## Request resources for the next class
After saving Tutorial 3 - Seurat.ipynb and Tutorial 3 - scanpy.ipynb, return to the original terminal window where you launched jupyter. Quit jupyter by pressing command+c twice.

### Relinquish the current wheat allocation
After running the command below, your terminal should switch back to rice from wheat.
```bash
exit
```
### Run salloc
Start "tmux" session on rice and then then request resources
```bash
tmux
## You should see the green tmux bar. Now run salloc requesting resources for the next class
salloc --ntasks-per-node=1 --cpus-per-task=4 --mem=30G --time=0-3:00:00 --begin="13:30:00 10/26/20" --qos=interactive srun --pty bash -i -l
```
**Note:** Sherlock users should use --quit=normal as Sherlock does not have an interactive quality of service like FarmShare.

You can now detach from the "tmux" with control+b and then press "d". The green bar on the bottom should disappear. Note down the rice node# (reice\<xx\>) w
here your "tmux" session is running (in case we need to get back to it) and the session ID (usually just a single session==0). After that, you can safely close your terminal window. See you next class!
