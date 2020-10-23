# Tutorial 1 - For SCG

## Connect to the SCG cluster
Run this on Terminal on a Mac or linux computer and PowerShell on Windows 10
```bash
ssh <SUNetID>@login.scg.stanford.edu
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

if [[ \$USER != "root" ]];
then
        export PS1='[\\[\\e[01;32m\\]\\u\\[\\e[m\\]@\\h \\W]\\\$ '
                else
                export PS1='[\\[\\e[01;31m\\]\\u\\[\\e[m\\]@\\h \\W]\\\$ '
fi

alias sbc='source ~/.bashrc'
alias la='ls -lA --color=auto'
alias lh='ls -lAh --color=auto'
alias lt='ls -ltrh --color=auto'

module load gcc/9.2.0-centos_7 pcre2/10.34 icu4c/67_1 texlive/2013 texinfo/6.7

SAPP=\$HOME/sysapps
APP=\$HOME/apps

export PKG_CONFIG_PATH=\\
\$SAPP/gdalv3.1.3/lib/pkgconfig:\\
\$SAPP/projv7.1.1/lib/pkgconfig:\\
\$SAPP/sqlitev3.33/lib/pkgconfig:\\
\$PKG_CONFIG_PATH

export LD_LIBRARY_PATH=\\
\$SAPP/gdalv3.1.3/lib:\\
\$SAPP/projv7.1.1/lib:\\
\$SAPP/sqlitev3.33/lib:\\
\$LD_LIBRARY_PATH

export  LIBRARY_PATH=\\
\$SAPP/gdalv3.1.3/lib:\\
\$SAPP/projv7.1.1/lib:\\
\$SAPP/sqlitev3.33/lib:\\
\$LIBRARY_PATH

export CPATH=\\
\$SAPP/gdalv3.1.3/include:\\
\$SAPP/projv7.1.1/include:\\
\$SAPP/sqlitev3.33/include:\\
\$CPATH

export PATH=\\
\$SAPP/gdalv3.1.3/bin:\\
\$SAPP/projv7.1.1/bin:\\
\$SAPP/sqlitev3.33/bin:\\
\$APP/cellranger-4.0.0:\\
\$APP/R/bin:\\
\$PATH

export CMAKE_PREFIX_PATH=\\
\$SAPP/gdalv3.1.3:\\
\$SAPP/projv7.1.1:\\
\$SAPP/sqlitev3.33:\\
\$CMAKE_PREFIX_PATH
EOF
```

**IMPORTANT:** Logout and login again to SCG

## Install sqlite3 PROJ, and GDAL
**NOTE:** You need to replace your SUNet ID in the commands below
```bash
mkdir sysapps && cd sysapps

wget --quiet https://sqlite.org/2020/sqlite-autoconf-3330000.tar.gz
tar xf sqlite-autoconf-3330000.tar.gz 
cd sqlite-autoconf-3330000
./configure --prefix=/home/<SUNetID>/sysapps/sqlitev3.33 --disable-static --enable-fts5 CFLAGS="-g -O2 -DSQLITE_ENABLE_FTS3=1 -DSQLITE_ENABLE_FTS4=1 -DSQLITE_ENABLE_COLUMN_METADATA=1 -DSQLITE_ENABLE_UNLOCK_NOTIFY=1 -DSQLITE_ENABLE_DBSTAT_VTAB=1 -DSQLITE_SECURE_DELETE=1 -DSQLITE_ENABLE_FTS3_TOKENIZER=1"
make -j64 && make install && cd ..

wget --quiet https://download.osgeo.org/proj/proj-7.1.1.tar.gz
tar xf proj-7.1.1.tar.gz
cd proj-7.1.1
./configure --prefix=/home/<SUNetID>/sysapps/projv7.1.1
make -j64 && make install && cd ..

wget --quiet https://github.com/OSGeo/gdal/releases/download/v3.1.3/gdal-3.1.3.tar.gz
tar xf gdal-3.1.3.tar.gz 
cd gdal-3.1.3
./configure --prefix=/home/<SUNetID>/sysapps/gdalv3.1.3 --with-proj=/home/<SUNetID>/sysapps/projv7.1.1
make -j64 && make install && cd ..

rm -rf sqlite-autoconf-3330000* proj-7.1.1* gdal-3.1.3*
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

On SCG it is ok to install conda in its default location. At this step just hit the Enter key
```
Miniconda3 will now be installed into this location:
/home/<SUNetID>/miniconda3

  - Press ENTER to confirm the location
  - Press CTRL-C to abort the installation
  - Or specify a different location below

[/home/<SUNetID>/miniconda3] >>> 
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
conda install -c conda-forge python-igraph leidenalg jupyterlab nodejs adjustText seaborn scikit-learn statsmodels numba pytables libcurl cython cmake

conda install -c bioconda -c conda-forge samtools kallisto mygene scvi

pip install scanpy==1.6 MulticoreTSNE bbknn velocyto

pip install seaborn==0.10

pip install synapseclient[pandas,pysftp]
```

## Fetch some binaries and scripts and install in "singlecell" environment
```bash
wget --quiet https://github.com/alexdobin/STAR/archive/2.6.1d.tar.gz
tar -xf 2.6.1d.tar.gz
cp ./STAR-2.6.1d/bin/Linux_x86_64/* ~/miniconda3/envs/singlecell/bin/

wget --quiet http://downloads.sourceforge.net/project/skewer/Binaries/skewer-0.2.2-linux-x86_64 -O ~/miniconda3/envs/singlecell/bin/skewer
chmod +x ~/miniconda3/envs/singlecell/bin/skewer

wget --quiet https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary -O ~/miniconda3/envs/singlecell/bin/bedtools
chmod a+x ~/miniconda3/envs/singlecell/bin/bedtools

wget --quiet https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz
tar -xf tophat-2.1.1.Linux_x86_64.tar.gz
cp ./tophat-2.1.1.Linux_x86_64/gtf_to_fasta ~/miniconda3/envs/singlecell/bin/

rm -rf *2.6.1d* *tophat*
```

## Install Cell Ranger from 10X Genomics. The download page requires signup "https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest". Either signup and follow wget or curl download instruction or to avoid that simply download a pre-downloaded copy.
```bash
cd && mkdir apps && cd apps
wget --quiet http://hsc.stanford.edu/resources/cellranger-4.0.0.tar.gz
tar xf ./cellranger-4.0.0.tar.gz
rm cellranger-4.0.0.tar.gz
which cellranger
```

## Configure environment To access MariaDB server at UCSC and use their tools and utilities. Details: "https://genome.ucsc.edu/goldenPath/help/mysql.html"
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

## Install SCILIAN, MAGIC, and SAM from source
```bash
git clone https://github.com/salzmanlab/SICILIAN.git
cd SICILIAN && pip install -r requirements.txt && cd ..

git clone git://github.com/atarashansky/self-assembling-manifold.git
pip install ./self-assembling-manifold/.

git clone git://github.com/KrishnaswamyLab/MAGIC.git
pip install ./MAGIC/python/.

rm -rf SICILIAN/ self-assembling-manifold/ MAGIC/
```

## Install R-4.0.3 and required R packages

#### First deactivate "singlecell" conda environment and load Sherlock's R module
```bash
conda deactivate

## download and compile latest R
wget --quiet https://stat.ethz.ch/R/daily/R-patched.tar.gz
tar xf R-patched.tar.gz 

cd R-patched/tools/ && ./rsync-recommended
./link-recommended && cd ..
./configure --prefix=$HOME/apps/R AR=gcc-ar NM=gcc-nm RANLIB=gcc-ranlib LDFLAGS='-L/scg/apps/software/gcc/9.2.0-centos_7/lib64 -L/scg/apps/software/gcc/9.2.0-centos_7/lib64' CFLAGS=-I/scg/apps/software/gcc/9.2.0-centos_7/include --enable-lto
make -j64 && make install && cd
rm -rf R-patched/

## setup R profile and environment
echo -e "AR=gcc-ar\nNM=gcc-nm\nRANLIB=gcc-ranlib" > ~/.Renviron
echo -e "options(Ncpus = 30)\nmessage(\"\\\033[38;05;208mHi <SUNetID>, welcome to R...\\\033[00m\")\nSys.setenv(RETICULATE_PYTHON = \"/home/<SUNetID>/miniconda3/envs/singlecell/bin/python\")" >> ~/.Rprofile

which R
```

#### Run R
```bash
R
```

#### Install IRKernel (in R) and quit R
```bash
install.packages("BiocManager")
```

Type 1 then press Enter when asked to select a CRAN mirror
```bash
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

#### Launch jupyter
```bash
jupyter lab --no-browser
```

Note down the \<Port\> in the jupyter lab web URL

#### Now on your local computer create an SSH tunnel to Sherlock

On a Mac, run this command in Terminal, on a PC run it within PowerShell
```bash
ssh -N -f -L <Port>:localhost:<Port> <SUNetID>@scg.stanford.edu
```

#### Copy the line with web address and open it in your browser

You can now access terminal through Jupyter to install the R packages that we need.

While you are working through the R installs (some of them can take a while), please open the Tutorial 1 excersize in your browser on jupyter (File > Open from Path... > "/home/groups/\<GroupName\>/BIOC281/Classes/1/Tutorial 1.ipynb") and complete it. Periodically, go back check the terminal to see if the package installations have completed.

**Note:** When finished, please save the outputs of Tutorial 1.ipynb (File > Export Notebook As... > Export Notebook to HTML) when you have completed it and upload it to Canvas.

#### Start R in jupyter (File > New... > Terminal)
```bash
R
```

#### Install Seurat (in R) 
```bash
BiocManager::install(c("Seurat", "useful", "here", "RColorBrewer", "plotly", "genieclust"))
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
rm ./CytoTRACE_0.3.3.tar.gz
``` 
    
If you have finished the first exercise, saved the output and uploaded it to canvas you can now return to the terminal window you used to launch jupyter. Kill the program by pressing "ctrl+c" twice


