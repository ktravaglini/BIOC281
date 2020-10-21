# Tutorial 1 - For Sherlock

## Connect to the Sherlock cluster
Run this on Terminal on a Mac or linux computer and PowerShell on Windows 10
```bash
ssh <SUNetID>@login.sherlock.stanford.edu
```

## Backup your old \~/\.bashrc
```bash
cp ~/.bashrc ~/.bashrc.old
## make another backup just in case!!
cp ~/.bashrc ~/.bashrc.bak
```
## Create a clean ~/.bashrc
```bash
cat > ~/.bashrc << EOF
# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi

# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=

# User specific aliases and functions
EOF
```

**IMPORTANT:** Logout and login again to Sherlock

## Load gcc 10 for building R and python packages
```bash
ml load gcc/10.1.0
```

## Navigate to $GROUP_HOME
We're installing packages in $GROUP_HOME for two reasons: First, it has more space allocated per lab (1TB) compared to the $HOME folder (15GB), and second, the setup after install can be used by any lab member since this location is accessible to everyone in the group.
```bash
cd $GROUP_HOME
```

## Install conda
```bash
wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh
```

**Note:** Answer the prompts accordingly

```
Do you accept the license terms? [yes|no]
[no] >>> yes
```

```
Miniconda3 will now be installed into this location:
/home/users/<SUNetID>/miniconda3

  - Press ENTER to confirm the location
  - Press CTRL-C to abort the installation
  - Or specify a different location below

[/home/users/<SUNetID>/miniconda3] >>> /home/groups/<GroupName>/miniconda3
```

```
Do you wish the installer to initialize Miniconda3
by running conda init? [yes|no]
[no] >>> yes
```

## Source the ~/.bashrc that was modified by conda to activate "base" conda environment
```bash
source ~/.bashrc
rm Miniconda3-latest-Linux-x86_64.sh
```

## Similar to FarmShare check if "base" conda installed correctly
These should all point binaries in to /home/groups/\<GroupName\>/miniconda3/bin
```bash
which conda
which python
which pip
python --version
```

## Update conda
```bash
conda update -n base -c defaults conda
```

## Create "singlecell" environment and activate it
```bash
conda create --name singlecell python=3.8.5

conda activate singlecell
```

## Install python packages for "singlecell" environment
```bash
conda install -c conda-forge python-igraph leidenalg jupyterlab nodejs adjustText seaborn scikit-learn statsmodels numba pytables libcurl

conda install -c conda-forge cython cmake

conda install scvi -c bioconda -c conda-forge

conda install -c bioconda -c conda-forge samtools kallisto mygene

pip install scanpy==1.6 MulticoreTSNE bbknn velocyto

pip install seaborn==0.10

pip install synapseclient[pandas,pysftp]
```

## Fetch some binaries and scripts and install in "singlecell" environment
```bash
wget --quiet https://github.com/alexdobin/STAR/archive/2.6.1d.tar.gz
tar -xf 2.6.1d.tar.gz
cp ./STAR-2.6.1d/bin/Linux_x86_64/* ./miniconda3/envs/singlecell/bin/

wget --quiet http://downloads.sourceforge.net/project/skewer/Binaries/skewer-0.2.2-linux-x86_64 -O ./miniconda3/envs/singlecell/bin/skewer
chmod +x ./miniconda3/envs/singlecell/bin/skewer

wget --quiet https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary -O ./miniconda3/envs/singlecell/bin/bedtools
chmod a+x ~/miniconda3/envs/singlecell/bin/bedtools

wget --quiet https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz
tar -xf tophat-2.1.1.Linux_x86_64.tar.gz
cp ./tophat-2.1.1.Linux_x86_64/gtf_to_fasta ./miniconda3/envs/singlecell/bin/

rm -rf *2.6.1d* *tophat*
```

## Install Cell Ranger from 10X Genomics. The download page requires signup "https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest". Either signup and follow wget or curl download instruction or to avoid that simply download a pre-downloaded copy.
```bash
wget --quiet http://hsc.stanford.edu/resources/cellranger-4.0.0.tar.gz
tar xf ./cellranger-4.0.0.tar.gz
echo 'export PATH="$GROUP_HOME/cellranger-4.0.0:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

## Configure environment To access MariaDB server at UCSC and use their tools and utilities. Details: "https://genome.ucsc.edu/goldenPath/help/mysql.html"
```bash
echo -e "#US MariaDB server\ndb.host=genome-mysql.soe.ucsc.edu\ndb.user=genomep\ndb.password=password\ncentral.db=hgcentral\ncentral.host=genome-mysql.soe.ucsc.edu\ncentral.user=genomep\ncentral.password=password\ngbdbLoc1=http://hgdownload.soe.ucsc.edu/gbdb/\nforceTwoBit=on" > ~/.hg.conf

chmod 600 ~/.hg.conf

wget --quiet http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa -O ./miniconda3/envs/singlecell/bin/twoBitToFa
wget --quiet http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit -O ./miniconda3/envs/singlecell/bin/faToTwoBit
wget --quiet http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf -O ./miniconda3/envs/singlecell/bin/genePredToGtf

chmod +x ./miniconda3/envs/singlecell/bin/twoBitToFa \
./miniconda3/envs/singlecell/bin/faToTwoBit \
./miniconda3/envs/singlecell/bin/genePredToGtf
```

## Install SCILIAN, MAGIC, and SAM from source
```bash
git clone https://github.com/salzmanlab/SICILIAN.git
cd SICILIAN && pip install -r requirements.txt && cd ..

git clone git://github.com/atarashansky/self-assembling-manifold.git
pip install ./self-assembling-manifold/.

git clone git://github.com/KrishnaswamyLab/MAGIC.git
pip install ./MAGIC/python/.
```

## Install R packages

#### First deactivate "singlecell" conda environment and load Sherlock's R module
```bash
conda deactivate

ml load java R/4.0.2 physics gdal proj geos expat

which R
```

#### Set these modules so they are loaded everytime you login to sherlock
```bash
ml save
echo 'module load gcc/10.1.0 R/4.0.2 java physics gdal proj geos expat' >> ~/.bashrc
```

#### Setup R environment and profile
```bash
echo -e "AR=gcc-ar\nNM=gcc-nm\nRANLIB=gcc-ranlib" > ~/.Renviron
echo -e "options(Ncpus = 4)\nmessage(\"\\\033[38;05;208mHi <SUNetID>, welcome to R...\\\033[00m\")\nSys.setenv(RETICULATE_PYTHON = \"/home/groups/<PI_Name>/miniconda3/envs/singlecell/bin/python\")" >> ~/.Rprofile
```

#### Run R
```bash
R
```

#### Install IRKernel (in R) and quit R
```bash
install.packages("BiocManager")
```

**Note:** Answer the prompts accordingly

```
Warning in install.packages("BiocManager") :
  'lib = "/share/software/user/open/R/4.0.2/lib64/R/library"' is not writable
Would you like to use a personal library instead? (yes/No/cancel) yes
```

```
Would you like to create a personal library
‘~/R/x86_64-pc-linux-gnu-library/4.0’
to install packages into? (yes/No/cancel) yes
```

```
--- Please select a CRAN mirror for use in this session ---
Secure CRAN mirrors
1
```

#### Install additional packages with BiocManager
```bash
BiocManager::install("devtools")
BiocManager::install("remotes")
BiocManager::install("snoweye/pbdZMQ")
BiocManager::install("IRkernel")
q()
```

#### Activate "singlecell" environment and run R
```bash
conda activate singlecell

R
```

#### Run the following command in R to connect IRkernel to jupyter and then quit R
```bash
IRkernel::installspec(user = TRUE)
q()
```

#### Clone the class github repository
```bash
git clone https://github.com/ktravaglini/BIOC281.git
```

#### Request resources from Sherlock
```bash
salloc --ntasks-per-node=1 --cpus-per-task=4 --mem=30G --time=0-3:00:00 --qos=normal srun --pty bash -i -l
```

#### Once resources have been allocated to you, launch jupyter
**Note:** Pick a port anywhere in the range 49152-65335
```bash
conda activate singlecell
jupyter lab --no-browser --port=<Port#>
```
If the port you chose is not available (likely because it is in use by someone else), terminate jupyter lab and run it again with a different port.

#### Now on your local computer create an SSH tunnel to Sherlock
**Note:** that you need to input the specicific compute node name that was allocated to you.
It is the part of your command prompt bolded here: (singlecell) [\<SUNetID\>@**sh03-01n52** ~]$

On a Mac, run this command in Terminal, on a PC run it within PowerShell
```bash
ssh -L <Port#>:localhost:<Port#> <SUNetID>@login.sherlock.stanford.edu ssh -L <Port#>:localhost:<Port#> -N <shXX-XXnXX>
```

#### Copy the line with web address and open it in your browser
**Note:** The notebook will only be active for a little less than 3 hours, so finishes this part as quickly as possible

You can now access terminal through Jupyter to install the R packages that we need.

While you are working through the R installs (some of them can take a while), please open the Tutorial 1 excersize in your browser on jupyter (File > Open from Path... > "/home/groups/\<GroupName\>/BIOC281/Classes/1/Tutorial 1.ipynb") and complete it. Periodically, go back check the terminal to see if the package installations have completed.

**Note:** When finished, please save the outputs of Tutorial 1.ipynb (File > Export Notebook As... > Export Notebook to HTML) when you have completed it and upload it to Canvas.

#### Start R in jupyter (File > New... > Terminal)
```bash
R
```

#### Install Seurat (in R) and quit R
```bash
BiocManager::install(c("Seurat", "useful", "here", "RColorBrewer", "plotly", "genieclust"))
q()
``` 

#### Install a dependency for monocle3
```bash
wget --quiet ftp://ftp.unidata.ucar.edu/pub/udunits/udunits-2.2.26.tar.gz
tar -xf udunits-2.2.26.tar.gz
cd udunits-2.2.26
./configure --prefix=/home/groups/<GroupName>/udunitsv2.2.26
make -j8 && make install && cd ..
rm -rf udunits-2.2.26*
```

#### Add the udunit2 package's compiler and runtime PATH locations to your ~/.bashrc
```bash
echo 'export LIBRARY_PATH="$GROUP_HOME/udunitsv2.2.26/lib:$LIBRARY_PATH"' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH="$GROUP_HOME/udunitsv2.2.26/lib:$LD_LIBRARY_PATH"' >> ~/.bashrc
echo 'export PATH="$GROUP_HOME/udunitsv2.2.26/bin:$PATH"' >> ~/.bashrc
echo 'export CPATH="$GROUP_HOME/udunitsv2.2.26/include:$CPATH"' >> ~/.bashrc
source ~/.bashrc
which udunits2
```

#### Now start R
```bash
R
```

#### Install monocle3 from source
```bash
BiocManager::install("cole-trapnell-lab/monocle3")
```

#### Install some dependencies for CytoTRACE and quit R
```bash
BiocManager::install(c("sva", "HiClimR", "ccaPP", "nnls", "egg", "ggpubr"))
q()
```

#### Download CytoTRACE and install CytoTRACE
```bash
wget --quiet https://cytotrace.stanford.edu/CytoTRACE_0.3.3.tar.gz --no-check-certificate
R CMD INSTALL ./CytoTRACE_0.3.3.tar.gz 
rm -fr ./CytoTRACE_0.3.3.tar.gz
``` 
    
If you have finished the first exercise, saved the output and uploaded it to canvas you can now return to the terminal window you used to launch jupyter. Kill the program by pressing "ctrl+c" twice

#### Relinquish control of the "compute" node as we are done with computing
```bash
exit
```
You should be returned to the "login" node

## Request resources for the next class
    
#### Start a persistent terminal session
```bash
tmux
```
You should see a green bar at the bottom of the terminal
    
#### Run salloc
```bash
salloc --ntasks-per-node=1 --cpus-per-task=4 --mem=30G --time=0-3:00:00 --begin="13:30:00 10/21/20" --qos=normal srun --pty bash -i -l
```
This command will hang as its waiting until Wednesday to obtain resources. Detach from tmux by pressing "ctrl+b" followed by "d" key. This will return you from the "tmux" prompt to the main "login" prompt (no green bar). Write down which "login" node your persistent "tmux" session is running before logging out so we can get back to the same "tmux" session on Wednesday. You can then safely close terminal.
