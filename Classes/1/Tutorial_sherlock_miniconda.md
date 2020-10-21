## backup your old \~/\.bashrc
```bash
cp ~/.bashrc ~/.bashrc.old
## make another backup just in case!!
cp ~/.bashrc ~/.bashrc.bak
```
## create a clean new bashrc
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
# logout and login again to Sherlock

## After login change to $GROUP_HOME
```bash
cd $GROUP_HOME
```

## Install miniconda
```bash
wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh

## when asked specify the install location to $GROUP_HOME/miniconda3-- type full path as shown below. Add your groupname without "<>" and press Enter (or Return) key
Do you accept the license terms? [yes|no]
[no] >>> yes

Miniconda3 will now be installed into this location:
/home/users/SUNetID/miniconda3

  - Press ENTER to confirm the location
  - Press CTRL-C to abort the installation
  - Or specify a different location below

[/home/users/sinhar/miniconda3] >>> /home/groups/<yourgroupname>/miniconda3
```

## When asked if you want to activate conda-- type yes then Enter (Return)
```bash
Do you wish the installer to initialize Miniconda3
by running conda init? [yes|no]
[no] >>> yes
```

## Source the bashrc that was modified by conda to activate "base" conda environment
```bash
source ~/.bashrc
rm -fr Miniconda3-latest-Linux-x86_64.sh
```

## like farmshare check if "base" conda installed correctly
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

## Install essential packages for "singlecell" environment
```bash
conda install -c conda-forge python-igraph leidenalg jupyterlab nodejs adjustText seaborn scikit-learn statsmodels numba pytables libcurl
conda install -c conda-forge cython cmake
conda install scvi -c bioconda -c conda-forge
conda install -c bioconda -c conda-forge samtools kallisto mygene
pip install MulticoreTSNE bbknn velocyto
pip install seaborn==0.10
pip install synapseclient[pandas,pysftp]
```

## Fetch some binaries and scripts and install in "singlecell" environment
```bash
wget --quiet https://github.com/alexdobin/STAR/archive/2.6.1d.tar.gz
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

## Install SCILIAN, MAGIC, scVI from source
```bash
git clone https://github.com/salzmanlab/SICILIAN.git
cd SICILIAN && pip install -r requirements.txt && cd ..

git clone git://github.com/atarashansky/self-assembling-manifold.git
pip install ./self-assembling-manifold/.

git clone git://github.com/KrishnaswamyLab/MAGIC.git
pip install ./MAGIC/python/.
```

## Install R packages
#### First deactivate "singlecell" conda environment and check lates R module is loaded
```bash
conda deactivate
which R
ml load java R/4.0.2 physics gdal proj geos
ml
which R
## If you want these module to be loaded everytime you login to sherlock
ml save
## If you don't save you will need to load the modules everytime you login'
```
## stup R environment and profile
```bash
echo 'module load R/4.0.2 java physics gdal proj geos' >> ~/.bashrc
echo -e "AR=gcc-ar\nNM=gcc-nm\nRANLIB=gcc-ranlib" > ~/.Renviron

echo -e "options(Ncpus = 4)\nmessage(\"\\\033[38;05;208mHi <SUNetID>, welcome to R...\\\033[00m\")\nSys.setenv(RETICULATE_PYTHON = \"/home/<SUNetID>/miniconda3/envs/singlecell/bin/python\")" >> ~/.Rprofile
```

## Run R
```bash
R
```

## Install IRKernel (in R) and quit R
```bash
install.packages("BiocManager")
```

When prompted type yes (twice) and then Enter (Return) key. This will install R packages in you ~/R directory. You can select any mirror by typing its number and then Enter (Return) key.
```bash
Warning in install.packages("BiocManager") :
  'lib = "/share/software/user/open/R/4.0.2/lib64/R/library"' is not writable
Would you like to use a personal library instead? (yes/No/cancel) yes
Would you like to create a personal library
‘~/R/x86_64-pc-linux-gnu-library/4.0’
to install packages into? (yes/No/cancel) yes
--- Please select a CRAN mirror for use in this session ---
Secure CRAN mirrors 
```

#### now install all other current and future packages using BiocManager
```bash
BiocManager::install("devtools")
BiocManager::install("remotes")
BiocManager::install("snoweye/pbdZMQ")
BiocManager::install("IRkernel")
q()
```

## activate "singlecell" environment and run R
```bash
conda activate singlecell

R
```

## Run the following command in R to connect IRkernel to jupyter and then quit R
```bash
IRkernel::installspec(user = TRUE)
q()
```

## clone the class github repository
```bash
git clone https://github.com/ktravaglini/BIOC281.git
```
## Request resources
```bash
salloc --ntasks-per-node=1 --cpus-per-task=4 --mem=30G --time=0-3:00:00 --qos=normal srun --pty bash -i -l
```

Once resources have been allocated to you, run jupyter lab -- pick a port anywhere in the range 49152-65335. Example shown below
```bash
jupyter lab --no-browser --port=<Port#>
```

## Now on you Mac or PC local shell run the following command.
Note that you need to input the specicific compute node name that was allocated to you.
```bash
ssh -L <Port#>:localhost:<Port> sherlock ssh -L <Port>:localhost:<Port> -N sh-xxx-xx
```

At this you find that the port you have chosen is not available--probably used by someone else. If so, terminate jupyter lab and run it again with a different port.. the try tunnelling again with that port.

# IMPORTANT: visit https://vsoch.github.io/lessons/sherlock-jupyter/ and follow instructions to fully secure your broadcasted notebook
For now Launch browser on your system/laptop and login to jupyter running remotely on Sherlock -- we will be running these notebooks for only a few hours
Copy the line with web address "http://localhost:<Port>/?token=e968329a256d1264f643d2bf3fa72fc75292446d9d337b3a" from the terminal and paste it into your browser

You can now access terminal in the browser window connected to Jupyter running on rice, which we can use to start installing R packages that we need. While you are working through the R installs (some of them can take a while), please open BIOC281/Classes/1/Tutorial 1.ipynb through jupyter connected through your browser and begin working through the exercise. Periodically, go back check the terminal to see if the package installations have completed.

When finished, please save the outputs of Tutorial 1.ipynb (File > Export Notebook As... > Export Notebook to HTML) when you have completed it and upload it to Canvas. Aditotors if interested can email us the results.
## Start R in jupyter (File > New... > Terminal)
Note and confirm that terminal launches automatically in the "base" conda environment [(base) [sinhar\@sh03\-ln07 login\ ~]\$]
```bash
R
```

## Install Seurat (in R) and quit R
```bash
BiocManager::install(c("Seurat", "useful", "here", "RColorBrewer", "plotly", "genieclust"))
q()
```

## Install a dependency for monocle3
```bash
wget --quite ftp://ftp.unidata.ucar.edu/pub/udunits/udunits-2.2.26.tar.gz
tar -xf udunits-2.2.26.tar.gz
cd udunits-2.2.26
## sanity check
which gcc
ml
## If gcc 10.1.0 is lodaed and in system PATH, we are ready to configure and compile udunits
## please note that you need to replace your GroupName
./configure --prefix=/home/groups/<GroupName>/udunitsv2.2.26
make -j8 && make install && cd ..
rm -rf udunits-2.2.26*
ll udunitsv2.2.26/*

## add the newly installed location to your system by appending .bashrc
echo 'export LIBRARY_PATH="$GROUP_HOME/udunitsv2.2.26/lib:$LIBRARY_PATH"' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH="$GROUP_HOME/udunitsv2.2.26/lib:$LD_LIBRARY_PATH"' >> ~/.bashrc
echo 'export PATH="$GROUP_HOME/udunitsv2.2.26/bin:$PATH"' >> ~/.bashrc
echo 'export CPATH="$GROUP_HOME/udunitsv2.2.26/include:$CPATH"' >> ~/.bashrc
source ~/.bashrc
which udunits2
```

## Now start R
```bash
R
```

## Install monocle3 from source
```bash
BiocManager::install("cole-trapnell-lab/monocle3")
```

## Install some dependencies for CytoTRACE and quit R
```bash
BiocManager::install(c("sva", "HiClimR", "ccaPP", "nnls", "egg", "ggpubr"))
q()
```

## Download CytoTRACE and install CytoTRACE
```bash
wget --quiet https://cytotrace.stanford.edu/CytoTRACE_0.3.3.tar.gz --no-check-certificate
R CMD INSTALL ./CytoTRACE_0.3.3.tar.gz 
rm -fr ./CytoTRACE_0.3.3.tar.gz
``` 
## Relinquish control node as we are done with computing
```bash
exit
```
you should be back to login node
## Request resources for the next class
```bash
tmux
## you should see green bar
salloc --ntasks-per-node=1 --cpus-per-task=4 --mem=30G --time=0-3:00:00 --begin="13:30:00 10/21/20" --qos=normal srun --pty bash -i -l

## now detach tmux by pressing "ctrl+b" followed by "d" key. Note the login node you are before logging out in case we need to get back to the tmux session
