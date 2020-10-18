# Tutorial 1
Today we will install the software we will be using throughout the course on Stanford FarmShare as well as completing a small data analysis exercise to familiarize you with jupyter. To begin, we need to connect to FarmShare using ssh (**s**ecure **sh**ell).

## On Windows 10 launch PowerShell or install SecureCRT

Right click on the Windows "Start" icon on lower left and click "Windows PowerShell"

Alternatively, install Stanford provided terminal app SecureCRT "https://uit.stanford.edu/software/scrt_sfx" (recommended)

After installing SecureCRT click file-->Quick Connect. In the "Hostname" field type "rice.stanford.edu" without the quotes. Make sure that "Save session" and "Open in tab" options are checked. Click "Connect" followed by "OK". When promted, type your SUNetID as Username and SUNetpassword as password. You may want to save your username and password so its easier to connect in future. Next time simply click on "Session Manager" and double click on the saved rice.stanford.edu session to connect.

## On Mac launch Terminal app

Press command+space keys together to access search. Type "Terminal" in search field and click Terminal.app to launch terminal

Once Terminal is running connect to FarmShare by typing the following command. Replace SUNetID with your SUNetID without the "<>" characters. Follow similar scheme for the rest of the tutorial: "<>" indicates a text string specific to you that needs to be replaced.

```bash
ssh <SUNetID>@rice.stanford.edu
```

## Switch to a persistent session
The current command line we opened on rice will exit when we disconnect from the server. Programs like tmux and screen allow for command lines to run on a remote server, even when your connection terminates. After running the command below, a green bar should appear at the bottom of your screen
```bash
tmux
```

## Do you have conda installed?

```bash
which conda
```

### If not, install it (follow prompts)

```bash
wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh
rm -fr Miniconda3-latest-Linux-x86_64.sh
```

### Activate conda
```bash
source ~/.bashrc
```

### Check the installation (should point to Miniconda)

```bash
which conda
which python
which pip
```

### Check python version (should be >3.8)
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

## Install python based packages
Moving forward, when using pip, install all python packages within "singlecell" environment after activating it and not in the "base" environment

### Install scanpy, jupyter, MulticoreTSNE, BBKNN, velocyto, samtools, kallisto, and synapse with conda and pip
Synapse is one of many resources where you can deposit processed single cell RNAseq data (others include the NIH's Gene Expression Omnibus or EMBL-EBI's BioStudies)
While these packages are installing, go to https://www.synapse.org and create a synapse account

```bash
conda install seaborn scikit-learn statsmodels numba pytables libcurl mygene adjustText

conda install -c conda-forge python-igraph leidenalg jupyterlab nodejs

pip install scanpy

pip install MulticoreTSNE bbknn

pip install velocyto

pip install seaborn==0.10

pip install synapseclient[pandas,pysftp]

conda install -c bioconda samtools kallisto

conda install scvi -c bioconda -c conda-forge
```

### Fetch some binaries and scripts from STAR, skewer, bedtools, and tophat packages and install them
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

### Install Cell Ranger from 10X Genomics. The download page requires signup "https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest". Either signup and follow wget or curl download instruction or to avoid that simply download a pre-downloaded copy.
```bash
wget --quiet http://hsc.stanford.edu/resources/cellranger-4.0.0.tar.gz
tar xf ~/cellranger-4.0.0.tar.gz
echo 'export PATH="$HOME/cellranger-4.0.0:$PATH"' >> ~/.bashrc
```

### Configure environment To access MariaDB server at UCSC and use their tools and utilities. Details: "https://genome.ucsc.edu/goldenPath/help/mysql.html"
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

### Install SICIALIAN, MAGIC, and scVI from GitHub source 
```bash
git clone https://github.com/salzmanlab/SICILIAN.git
cd SICILIAN && pip install -r requirements.txt && cd

git clone git://github.com/atarashansky/self-assembling-manifold.git
pip install ~/self-assembling-manifold/.

git clone git://github.com/KrishnaswamyLab/MAGIC.git
pip install ~/MAGIC/python/.
```

## Install R packages

### Deacticate the conda environment
```bash
conda deactivate
```

### Install R 4.0.2
**Note:** You need to replace <SUNetID> in the commands below.
```bash
module load gcc/9.3.0

wget --quiet http://cran.rstudio.com/src/base/R-4/R-4.0.2.tar.gz
tar -xf R-4.0.2.tar.gz

cd R-4.0.2 && ./configure --prefix=$HOME/R AR=gcc-ar NM=gcc-nm RANLIB=gcc-ranlib LDFLAGS='-L/farmshare/software/free/gcc/9.3.0/lib -L/farmshare/software/free/gcc/9.3.0/lib64' CFLAGS=-I/farmshare/software/free/gcc/9.3.0/include

make -j32 && make install && cd

echo 'export PATH="$HOME/R/bin:$PATH"' >> ~/.bashrc
echo 'module load gcc/9.3.0' >> ~/.bashrc
echo -e "AR=gcc-ar\nNM=gcc-nm\nRANLIB=gcc-ranlib" > ~/.Renviron

echo -e "options(Ncpus = 4)\nmessage(\"\\\033[38;05;208mHi <SUNetID>, welcome to R...\\\033[00m\")\nSys.setenv(RETICULATE_PYTHON = \"/home/<SUNetID>/miniconda3/envs/singlecell/bin/python\")" >> ~/.Rprofile

source ~/.bashrc
```

### Check that R is installed
```bash
which R
```

### Start R
```bash
R
```

### Install IRKernel (in R) and quit R
```R
install.packages("BiocManager")
BiocManager::install(c("remotes", "devtools"))
BiocManager::install("snoweye/pbdZMQ")
BiocManager::install("IRkernel")
q()
```

### Activate "singlecell" conda environment
```bash
conda activate singlecell
```

### Start R
```bash
R
```

### Connect IRkernel to jupyter and quit R
```R
IRkernel::installspec(user = TRUE)
q()
```

### Clone the class GitHub repository
```bash
git clone https://github.com/ktravaglini/BIOC281.git
```

### Run jupyter
```bash
jupyter lab --no-browser
```

### Output should look like:

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

### Create an SSH tunnel (on your system/laptop) to connect to Jupyter lab running on rice

**On a Windows 10 PC:** Right click on the Windows "Start" icon on lower left and click "Windows PowerShell" to lauch powershell.
Alternatively, If you installed SecureCRT then click "File-->Connect Local Shell"

**On a Mac:** Launch Terminal app as above 

On Mac or PC local shell, run the following command. The "\<Port\>" comes from the jupyter lab output (default is 8888, but changes if that's in use) and the specific rice system comes from your terminal window. For example, if I see "(base) ktrav@rice11:~$" I know I am on rice11.stanford.edu.
```bash
ssh -N -f -L <Port>:localhost:<Port> <SUNetID>@rice<XX>.stanford.edu
```

### Now Launch browser on your system/laptop and login to jupyter running remotely on rice
Copy the line with web address "http://localhost:8888/?token=e968329a256d1264f643d2bf3fa72fc75292446d9d337b3a" from the terminal and paste it into your browser

You can now access terminal in the browser window connected to Jupyter running on rice, which we can use to start installing R packages that we need. While you are working through the R installs (some of them can take a while), please open BIOC281/Classes/1/Tutorial 1.ipynb through jupyter connected through your browser and begin working through the exercise. Periodically, go back check the terminal to see if the package installations have completed.

When finished, please save the outputs of Tutorial 1.ipynb (File > Export Notebook As... > Export Notebook to HTML) when you have completed it and upload it to Canvas.

### Start R in jupyter (File > New... > Terminal)
Note and confirm thhat terminal launches automatically in the "base" conda environment ["(base) \<SUNetID\>@rice11:~$"]
```bash
R
```

### Install Seurat and Monocle3 (in R)

```R
BiocManager::install("Seurat", "useful", "here", "RColorBrewer", "plotly", "genieclust")
BiocManager::install("cole-trapnell-lab/monocle3")
```

### Install some dependencies for CytoTRACE and quit R
```bash
BiocManager::install(c("sva", "HiClimR", "ccaPP", "nnls", "egg", "ggpubr"))
q()
```
### Download CytoTRACE
```bash
wget --quiet https://cytotrace.stanford.edu/CytoTRACE_0.3.3.tar.gz
```

### Install CytoTRACE
```R
R CMD INSTALL ~/CytoTRACE_0.3.3.tar.gz    
```

### Cleanup CytoTRACE source
```bash
rm -fr ~/CytoTRACE_0.3.3.tar.gz
```

## Request resources for the next class
After saving Tutorial 1.ipynb and installing all the R packages, return to the original terminal window where you launched jupyter. Quit jupyter by pressing command+c twice (you should still see the green tmux bar). Now, request resources for the next class with salloc.
```bash
salloc --ntasks-per-node=1 --cpus-per-task=4 --mem=30G --time=0-3:00:00 --begin="13:30:00 10/21/20" --qos=interactive srun --pty bash -i -l
```
    
You can now detach from tmux with control+b and then press "d". If the green bar on the bottom disappears, you can safely close your terminal window. See you next class!
