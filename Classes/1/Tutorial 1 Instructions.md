# Tutorial 1

## On Windows 10 launch PowerShell or install SecureCRT

#### Right click on the Windows "Start" icon on lower left and click "Windows PowerShell"

#### Alternatively, install Stanford provided terminal app SecureCRT "https://uit.stanford.edu/software/scrt_sfx" (recommended)
#### After installing SecureCRT click file-->Quick Connect
#### In the "Hostname" field type "rice.stanford.edu" without the quotes
#### Make sure that "Save session" and "Open in tab" options are checked
#### click "Connect" followed by "OK". When promted type your SUNetID as Username and SUNetpassword as password
#### You may want to save your username and password so its easier to connect in future
#### Next time simmply click on "Session Manager" and double click on the saved rice.stanford.edu session to connect

## On Mac launch Terminal app

#### Press command+space keys together to access search. Type "Terminal" in search field and click Terminal.app to launch terminal

## Once Terminal is running connect to FarmShare by typing the following command

#### Replace SUNetID with your SUNetID without the "<>" characters.
#### follow similar scheme for the rest of the tutorial: "<>" indicates a text string specicific to you that needs to be replaced
```bash
ssh <SUNetID>@rice.stanford.edu
```

## Do you have conda installed?

```bash
which conda
```

#### If not, install it (follow prompts)

```bash
wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh
rm -fr Miniconda3-latest-Linux-x86_64.sh
```

#### Activate conda
```bash
source ~/.bashrc
```

#### Check the installation (should point to Miniconda)

```bash
which conda
which python
which pip
```

#### Check python version (should be >3.8)
```bash
python --version
```

## Update conda
```bash
conda update -n base -c defaults conda
```

## Create "singlecell" environment
```bash
conda create --name singlecell python=3.8.5
```

## Activate the environment
```bash
conda activate singlecell
```
# Moving forward, when using pip, install all python packages within "singlecell" environment after activating it and not in the "base" environment

## Install scanpy, jupyter, MulticoreTSNE, bbknn, velocyto, samtools, and kallisto with conda and pip
```bash
conda install seaborn scikit-learn statsmodels numba pytables libcurl
conda install -c conda-forge python-igraph leidenalg jupyterlab nodejs
pip install scanpy
pip install MulticoreTSNE bbknn
pip install velocyto

conda install -c bioconda samtools kallisto
```

## Fetch some binaries and scripts from STAR, skewer, bedtools, and tophat packages and install them
```bash
wget --quiet https://github.com/alexdobin/STAR/archive/2.6.1d.tar.gz
tar -xf 2.6.1d.tar.gz
cp ~/STAR-2.6.1d/bin/Linux_x86_64/* ~/miniconda3/envs/singlecell/bin/
rm -fr ~/*2.6.1d*

wget --quiet http://downloads.sourceforge.net/project/skewer/Binaries/skewer-0.2.2-linux-x86_64 -O ~/miniconda3/envs/singlecell/bin/skewer
chmod +x ~/miniconda3/envs/singlecell/bin/skewer

wget --quiet https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary -O ~/miniconda3/envs/singlecell/bin/bedtools
chmod a+x ~/miniconda3/envs/singlecell/bin/bedtools

wget --quiet https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz
tar -xf tophat-2.1.1.Linux_x86_64.tar.gz
cp ~/tophat-2.1.1.Linux_x86_64/gtf_to_fasta ~/miniconda3/envs/singlecell/bin/
rm -fr ~/*tophat*
```

#### To access MariaDB server at UCSC to use their tools and utilities
#### Details: "https://genome.ucsc.edu/goldenPath/help/mysql.html"
```bash
echo -e "#US MariaDB server\ndb.host=genome-mysql.soe.ucsc.edu\ndb.user=genomep\ndb.password=password\ncentral.db=hgcentral\ncentral.host=genome-mysql.soe.ucsc.edu\ncentral.user=genomep\ncentral.password=password\ngbdbLoc1=http://hgdownload.soe.ucsc.edu/gbdb/\nforceTwoBit=on" > ~/.hg.conf

chmod 600 ~/.hg.conf

wget --quiet http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa -O ~/miniconda3/envs/singlecell/bin/twoBitToFa
wget --quiet http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit -O ~/miniconda3/envs/singlecell/bin/faToTwoBit
wget --quiet http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf -O ~/miniconda3/envs/singlecell/bin/genePredToGtf

chmod +x ~/miniconda3/envs/singlecell/bin/twoBitToFa \
~/miniconda3/envs/singlecell/bin/faToTwoBit \
~/miniconda3/envs/singlecell/bin/genePredToGtf
```

## Install SICIALIAN from GitHub source 
```bash
git clone https://github.com/salzmanlab/SICILIAN.git
cd SICILIAN && pip install -r requirements.txt && cd
```

# Install R 4.0.2
## Moving forwaed always install R-packages using BiocManager in "base" environment after decativating "singlecell"
```bash
conda deactivate
```

#### Notice your command prompt indicates now you have moved from "singlecell" to "base" conda environment
```bash
module load gcc/9.3.0

wget --quiet http://cran.rstudio.com/src/base/R-4/R-4.0.2.tar.gz
tar -xf R-4.0.2.tar.gz

cd R-4.0.2 && ./configure --prefix=$HOME/R AR=gcc-ar NM=gcc-nm RANLIB=gcc-ranlib LDFLAGS='-L/farmshare/software/free/gcc/9.3.0/lib -L/farmshare/software/free/gcc/9.3.0/lib64' CFLAGS=-I/farmshare/software/free/gcc/9.3.0/include

make -j32 && make install && cd

echo 'export PATH="$HOME/R/bin:$PATH"' >> ~/.bashrc

echo 'module load gcc/9.3.0' >> ~/.bashrc

echo -e "AR=gcc-ar\nNM=gcc-nm\nRANLIB=gcc-ranlib" > ~/.Renviron

echo -e "options(Ncpus = 4)\nmessage(\"\\\033[38;05;208mHi <Rahul>, welcome to R...\\\033[00m\")\nSys.setenv(RETICULATE_PYTHON = \"/home/<SUNetID>/miniconda3/envs/singlecell/bin/python\")" >> ~/.Rprofile

source ~/.bashrc

```

## Check that R is installed
```bash
which R
```

## Start R
```bash
R
```

## Install IRKernel (in R) and quit R
```R
install.packages("BiocManager")

BiocManager::install(c("remotes", "devtools"))

BiocManager::install("snoweye/pbdZMQ")

BiocManager::install("IRkernel")

q()
```

## activate "singlecell" conda environment
```bash
conda activate singlecell
```

## Start R
```bash
R
```

## Connect IRkernel to jupyter and quit R
```R
IRkernel::installspec(user = TRUE)

q()
```

## Clone the class GitHub repository
git clone https://github.com/ktravaglini/BIOC281.git

## Create an SSH tunnel (on your system/laptop) to connect to Jupyter lab running on rice

## On a Windows 10 PC
#### Right click on the Windows "Start" icon on lower left and click "Windows PowerShell" to lauch powershell

#### Alternatively, If you installed SecureCRT then click "File-->Connect Local Shell"

## On a Mac
#### Launch Terminal app as above 

## On Mac or PC local shell run the following command
```bash
ssh -N -f -L 9999:localhost:8888 <SUNetID>@rice<XX>.stanford.edu
```

## Now move back to the terminal connected to rice and Launch Jupyter
#### You must be in "singlecell" environment before you run the command below
```bash
jupyter lab --no-browser
```
#### output should look like:

```bash
[I 01:28:43.597 LabApp] JupyterLab extension loaded from /home/sinhar/miniconda3/envs/singlecell/lib/python3.8/site-packages/jupyterlab
[I 01:28:43.598 LabApp] JupyterLab application directory is /home/sinhar/miniconda3/envs/singlecell/share/jupyter/lab
[I 01:28:43.609 LabApp] Serving notebooks from local directory: /home/sinhar
[I 01:28:43.609 LabApp] Jupyter Notebook 6.1.4 is running at:
[I 01:28:43.609 LabApp] http://localhost:8888/?token=e968329a256d1264f643d2bf3fa72fc75292446d9d337b3a
[I 01:28:43.609 LabApp]  or http://127.0.0.1:8888/?token=e968329a256d1264f643d2bf3fa72fc75292446d9d337b3a
[I 01:28:43.609 LabApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 01:28:43.643 LabApp] 
    
    To access the notebook, open this file in a browser:
        file:///home/sinhar/.local/share/jupyter/runtime/nbserver-22427-open.html
    Or copy and paste one of these URLs:
        http://localhost:8888/?token=e968329a256d1264f643d2bf3fa72fc75292446d9d337b3a
     or http://127.0.0.1:8888/?token=e968329a256d1264f643d2bf3fa72fc75292446d9d337b3a
```

## Now Lauch browser on your system/laptop and login to jupyter running remotely on rice
#### copy the line with webaddress "http://localhost:8888/?token=e968329a256d1264f643d2bf3fa72fc75292446d9d337b3a" from the terminal
#### Paste it into your browser and replace port# 8888 with 9999 so the address now looks like:
```bash
http://localhost:9999/?token=e968329a256d1264f643d2bf3fa72fc75292446d9d337b3a
```

## You can now access terminal in the browser window connected to Jupyter running on rice
## Now we will start installing R packages in the terminal connected to rice
## While you are working through the R installs (some of them can take a while), please open BIOC281/Classes1/Tutorial.ipynb through jupyter connected through your browser and begin working through the exercise. Periodically go back check the terminal to see if the package installations have completed.

## Make sure that the terminal in Jupyter in the browser has started with the "base" environment


## Start R
```bash
R
```

## Install Seurat and Monocle3 (in R)

```R
BiocManager::install("Seurat")
BiocManager::install("cole-trapnell-lab/monocle3")
```

## Install some dependencies for CytoTRACE and quit R
```bash
BiocManager::install(c("sva", "HiClimR", "ccaPP", "nnls", "egg", "ggpubr"))
q()
```
## Download CytoTRACE
```bash
wget --quiet https://cytotrace.stanford.edu/CytoTRACE_0.3.3.tar.gz
```

## Install CytoTRACE
```R
R CND INSTALL ~/CytoTRACE_0.3.3.tar.gz    
```

## Cleanup CytoTRACE source
```bash
rm -fr ~/CytoTRACE_0.3.3.tar.gz
```
