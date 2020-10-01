# Tutorial 1

## Connect to FarmShare
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

## Create ‘singlecell’ environment
```bash
conda create --name singlecell python=3.8.5
```

## Activate the environment
```bash
source activate singlecell
```

## Install scanpy, jupyter, samtools, and kallisto with conda
```bash
conda install seaborn scikit-learn statsmodels numba pytables libcurl

conda install -c conda-forge python-igraph leidenalg jupyterlab nodejs
pip install scanpy

conda install -c bioconda samtools kallisto
```

## Fetch skewer, STAR, samtools, and kallisto, bedtools, ucsc-tools, tophat, SICILIAN binaries
```bash
wget --quiet https://github.com/alexdobin/STAR/archive/2.6.1d.tar.gz
tar -xzf 2.6.1d.tar.gz
cp ~/STAR-2.6.1d/bin/Linux_x86_64/* ~/miniconda3/envs/singlecell/bin/
rm -fr ~/*2.6.1d*

wget --quiet http://downloads.sourceforge.net/project/skewer/Binaries/skewer-0.2.2-linux-x86_64 -O ~/miniconda3/envs/singlecell/bin/skewer
chmod +x ~/miniconda3/envs/singlecell/bin/skewer

wget --quiet https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary -O ~/miniconda3/envs/singlecell/bin/bedtools
chmod a+x ~/miniconda3/envs/singlecell/bin/bedtools

wget --quiet http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa -O ~/miniconda3/envs/singlecell/bin/twoBitToFa
chmod +x ~/miniconda3/envs/singlecell/bin/twoBitToFa

wget --quiet http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit -O ~/miniconda3/envs/singlecell/bin/faToTwoBit
chmod +x ~/miniconda3/envs/singlecell/bin/faToTwoBit

wget --quiet http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf -O ~/miniconda3/envs/singlecell/bin/genePredToGtf
chmod +x ~/miniconda3/envs/singlecell/bin/genePredToGtf
echo -e "db.host=genome-mysql.soe.ucsc.edu\ndb.user=genomep\ndb.password=password\ncentral.db=hgcentral" > ~/.hg.conf
chmod 600 ~/.hg.conf

wget --quiet https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz
tar -xzf tophat-2.1.1.Linux_x86_64.tar.gz
cp ~/tophat-2.1.1.Linux_x86_64/gtf_to_fasta ~/miniconda3/envs/singlecell/bin/
rm -fr ~/*tophat*
```

## Install SICIALIAN and velocyto through PyPI
```bash
git clone https://github.com/salzmanlab/SICILIAN.git

cd SICILIAN && pip install -r requirements.txt && cd ~

pip install cython

pip install velocyto
```

## Install R 4.0.2
```bash
conda deactivate

module load gcc/9.3.0

wget --quiet http://cran.rstudio.com/src/base/R-4/R-4.0.2.tar.gz

tar -xzf R-4.0.2.tar.gz

cd R-4.0.2 && ./configure --prefix=$HOME/R AR=gcc-ar NM=gcc-nm RANLIB=gcc-ranlib LDFLAGS='-L/farmshare/software/free/gcc/9.3.0/lib -L/farmshare/software/free/gcc/9.3.0/lib64' CFLAGS=-I/farmshare/software/free/gcc/9.3.0/include

make -j32 && make install

cd
echo 'export PATH="/home/<SUNetID>/R/bin:$PATH"' >> ~/.bashrc
echo 'module load gcc/9.3.0' >> ~/.bashrc
echo -e "AR=gcc-ar\nNM=gcc-nm\nRANLIB=gcc-ranlib" > ~/.Renviron
echo -e "options(Ncpus = 4)\nmessage('\033[38;05;208mHi <Name>, welcome to R...\033[00m')\nSys.setenv(RETICULATE_PYTHON = '/home/<SUNetID>/miniconda3/envs/singlecell/bin/python')" >> ~/.Rprofile
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
install.packages("devtools")

devtools::install_github("snoweye/pbdZMQ")

install.packages("IRKernel")

q()
```

## activate 'singlecell' conda environment
```bash
source activate singlecell
```

## Start R
```bash
R
````

## Connect IRKernel to jupyter and quit R
```R
IRkernel::installspec(user = TRUE)

q()
```

## Clone the class GitHub repository
git clone https://github.com/ktravaglini/BIOC281.git


## Create an SSH tunnel (on your system)
```bash
# On a mac
ssh -N -f -L <Port>:localhost:<Port> <SUNetID>@riceXX.stanford.edu
```
On windows: PuTTY (Configuration -> Connection -> SSH -> Tunnel)

## Launch Jupyter
```bash
jupyter lab --no-browser
```

## Login to jupyter in your browser and open a Terminal window (continue there)
While you are working through the R installs (some of them can take time), please open BIOC281/Classes1/Tutorial.ipynb through jupyter and begin working through the exercise. Periodically check the terminal to see if the package installations have completed.


## Download CytoTRACE
```bash
wget --quiet https://cytotrace.stanford.edu/CytoTRACE_0.3.3.tar.gz
```

## Start R
```bash
R
```

## Install Seurat (in R)

```R
install.packages("Seurat")
```

## Install CytoTRACE
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("sva")

devtools::install_local("~/CytoTRACE_0.3.3.tar.gz")
```

## Install monocle3 and quit R
```R
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', 'limma', 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

devtools::install_github('cole-trapnell-lab/leidenbase')

devtools::install_github('cole-trapnell-lab/monocle3')

q()
```

## Cleanup CytoTRACE source
```bash
rm -fr ~/CytoTRACE_0.3.3.tar.gz
```