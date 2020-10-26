# Tutorial 4
Today we will be working through two notebooks, one written in R and the other in python, that use the PAGA, velocyto, monocle3 and CytoTRACE, respectively, to map the dynamic process of hematopoietic stem cell (HSC) differentiation to granulocytes and monocytes (two of the main immune cell lineages). 

To begin the tutorials, we need to access dedicated resources on Stanford FarmShare, launch a jupyer notebbook, then download data from Tabula Muris Consortium et al (2018) _Nature_.

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

## Install ForceAtlas2
```bash
pip install fa2
```

## Pull the GitHub repository
```bash
git pull
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

## Navigate to the current tutorial folder
```bash
cd ~/BIOC281/Classes/4
```

**Note:** On sherlock run
```bash
cd $GROUP_HOME/BIOC281/Classes/4
```

## Download the Bone Marrow 10x dataset from Tabula Muris
This dataset contains ~3,000 cells captured using 10x from two mice. See Tabula Muris Consortium et al (2018) _Nature_ for details.

```bash
wget --quiet http://hsc.stanford.edu/resources/DataObjects/MACA_bonemarrow_10x_subsetted.loom -O MACA_bonemarrow_10x.loom 
wget --quiet http://hsc.stanford.edu/resources/DataObjects/MACA_bonemarrow_10x_metadata.txt
```

## Open and complete BIOC281/Classes/4/Tutorial 4 - PAGA and velocyto.ipypnb
When finished, please save the outputs of Tutorial 4 - PAGA and velocyto.ipynb (File > Export Notebook As... > Export Notebook to HTML) when you have completed it and upload it to Canvas.

## Open and complete BIOC281/Classes/4/Tutorial 4 - monocle3 and CytoTRACE.ipypnb
When finished, please save the outputs of Tutorial 4 - monocle3 and CytoTRACE.ipynb (File > Export Notebook As... > Export Notebook to HTML) when you have completed it and upload it to Canvas.

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
salloc --ntasks-per-node=1 --cpus-per-task=4 --mem=30G --time=0-5:00:00 --begin="13:30:00 10/28/20" --qos=interactive srun --pty bash -i -l
```
**Note:** Sherlock users should use --quit=normal as Sherlock does not have an interactive quality of service like FarmShare.

You can now detach from the "tmux" with control+b and then press "d". The green bar on the bottom should disappear. Note down the rice node (**rice\<xx\>**) where your "tmux" session is running (in case we need to get back to it) and the session ID (usually just a single session==0). After that, you can safely close your terminal window. See you next class!

