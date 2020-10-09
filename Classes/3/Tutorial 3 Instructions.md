# Tutorial 3

## Launching jupyter

### Login to FarmShare
```bash
ssh <SUNetID>@rice.stanford.edu
```

### Request resources needed for alignments
```bash
srun --pty --qos=interactive --cpus-per-task=4 --mem=30G --time=240 SHELL -l
```

### Activate the single cell environment
```bash
source activate singlecell
```

### Start jupyter
```bash
jupyter lab --no-browser
```

### Create an SSH tunnel (on your system)
```bash
ssh -N -f -L <Port>:localhost:<Port> <SUNetID>@wheat<XX>.stanford.edu
```

### Login to jupyter in your browser, open a Terminal window, and activate the single cell environment
```bash
source activate singlecell
```

### Pull the GitHub repository
```bash
git pull
```

### Navigate to the current tutorial folder
```bash
cd ~/BIOC281/Classes/3
```

## Download Human Lung Cell Atlas data and metadata from synapse (Travaglini et al (2020) _Nature_)
This dataset contains ~75,000 cells captured using 10x (90%) and SmartSeq2 (10%) from the lung of 3 human patients.
We will use it to explore normalization, feature selection, dimensionality rediction, clustering, and data integration
in Seurat and scanpy

```bash
synapse get syn21560510

synapse get syn21560511

synapse get syn21560409

synapse get syn21560410
```

