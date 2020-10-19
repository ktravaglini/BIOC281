# Tutorial 3
Today we will be working through two notebooks, one written in R and the other in python, that use the Seurat (https://satijalab.org/seurat/) and scanpy (https://scanpy.readthedocs.io/en/stable/) packages, respectively, to cluster part of the human lung cell atlas single cell dataset (https://hlca.ds.czbiohub.org). The notebooks are designed to build intuition for how changing various paramters affects clustering results and illustrate the utility of iterative subclustering. In addition to the standard normalization, the notebooks also demo some of the more complex procedures we discussed during the lecture as well as methods for integrating datasets constructed with different technologies.

To begin the tutorials, we need to request dedicated resources on Stanford FarmShare, launch a jupyer notebbook, then download the data from Synapse.

## Login to FarmShare
```bash
ssh <SUNetID>@rice.stanford.edu
```
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

## Create an SSH tunnel (on your system)
```bash
ssh -N -f -L <Port>:localhost:<Port> <SUNetID>@wheat<XX>.stanford.edu
```

## Login to jupyter in your browser, open a Terminal window, and activate the single cell environment
```bash
conda activate singlecell
```

Install a missing dependency, louvain, for today's  scanpy notebook
Louvain is largely replaced by leiden these days. However, louvain may still be usefull to recreate plots for some of published data for comparison
```bash
pip install louvain
```

## Pull the GitHub repository
```bash
git pull
```

## Navigate to the current tutorial folder
```bash
cd ~/BIOC281/Classes/3
```

## Login to synapse
```bash
synapse login -u <Username> --remember-me
```

## Download Human Lung Cell Atlas data and metadata from synapse (Travaglini et al (2020) _Nature_)
This dataset contains ~75,000 cells captured using 10x (90%) and SmartSeq2 (10%) from the lung of 3 human patients.

```bash
synapse get syn21560510

synapse get syn21560511

synapse get syn21560409

synapse get syn21560410
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
    
You can now detach from the "tmux" with control+b and then press "d". The green bar on the bottom should disappear. Note down the rice node# (reice\<xx\>) w
here your "tmux" session is running (in case we need to get back to it) and the session ID (usually just a single session==0). After that, you can safely cl
ose your terminal window. See you next class!
