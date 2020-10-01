# Tutorial 2

## Launching jupyter

### Login to FarmShare
```bash
ssh <SUNetID>@rice.stanford.edu
```

### Request resources needed for alignments
```bash
srun --pty --qos=interactive --cpus-per-task=4 --mem=30G SHELL -l
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
cd ~/BIOC281/Classes/2
```

## Build a STAR index without and without masking psuedoautosomal regions (PAR)

### Download human Y chromosome and ERCC fasta and annotation (GTF) from RefSeq (UCSC) and ThermoFisher
```bash
wget --quiet https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrY.fa.gz
gzip -d chrY.fa.gz
mv chrY.fa hg38.refGene.chrY.fa

genePredToGtf hg38 refGene hg38.refGene.gtf

wget --quiet https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip
unzip ERCC92.zip && rm -fr ERCC92.zip
```

### Check the genomes and annotations look OK
```bash
less hg38.refGene.chrY.fa

less hg38.refGene.gtf

less ERCC92.fa

less ERCC92.gtf
```

### Concatenate hg38 genome and annotation with ERCCs
```bash
cat hg38.refGene.chrY.fa ERCC92.fa > hg38.refGene.chrY.ERCC92.fa

cat hg38.refGene.gtf ERCC92.gtf > hg38.refGene.ERCC92.gtf
```

### Create a PAR mask from coordinates in Mangs et al (2007) Current Genomics
```bash
echo -e "chrY\t10001\t2781479\nchrY\t56887903\t57217415" > PAR1_2.bed

bedtools maskfasta -fi hg38.refGene.chrY.ERCC92.fa -bed PAR1_2.bed -fo hg38.refGene.chrY.PARhm.ERCC.fa
```

### Check if the bases are properly masked
```bash
faToTwoBit hg38.refGene.chrY.PARhm.ERCC.fa hg38.refGene.chrY.PARhm.ERCC.2bit

twoBitToFa hg38.refGene.chrY.PARhm.ERCC.2bit PAR1.fa -seq=chrY -start=10000 -end=2781480

twoBitToFa hg38.refGene.chrY.PARhm.ERCC.2bit PAR2.fa -seq=chrY -start=56887902 -end=57217416

less PAR1.fa

less PAR2.fa
```

### Build the indices in STAR
```bash

```

## Build a kallisto index
Open a new Terminal window while STAR is building the index and proceed below

### Build a transcriptome from hg38 genome and annotation
```bash
gtf_to_fasta hg38.refGene.ERCC92.gtf hg38.refGene.chrY.PARhm.ERCC.fa hg38.refGene.transcriptome.ERCC92.fa

cat hg38.refGene.transcriptome.chrY.PARhm.ERCC92.fa | perl -ne 'if ($_ =~/^\>\d+\s+\w+\s+(ERCC\S+)[\+\-]/){print ">$1\n"}elsif($_ =~ /\d+\s+([A-Z]+_\d+_?\d+?)/){print ">$1\n"}else{print $_}' > hg38.refGene.transcriptome.chrY.PARhm.ERCC92.clean.fa
```

### Check the transcriptome looks OK
```bash
less hg38.refSeq.transcriptome.fa
```

### Concatenate hg38 transcriptome with ERCCs
```bash
cat hg38.refSeq.transcriptome.fa ERCC92.fa > hg38.refSeq.transcriptome.ERCC92.fa
```

### Build the index in kallisto
```bash
kallisto index -i hg38.refSeq.transcriptome.ERCC92.idx hg38.refSeq.transcriptome.ERCC92.fa
```

