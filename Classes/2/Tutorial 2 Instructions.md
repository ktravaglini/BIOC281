# Tutorial 2

### Login to FarmShare
```bash
ssh <SUNetID>@rice.stanford.edu
```

### Request resources needed for alignments
```bash
salloc --ntasks-per-node=1 --cpus-per-task=4 --mem=30G --time=0-48:00:00 --qos=interactive srun --pty bash -i -l
```

### Activate the single cell environment
```bash
conda activate singlecell
```

### Start jupyter
```bash
jupyter lab --no-browser
```

### Initiate port forwarding through ssh tunnel on your system/laptop
```bash
ssh -N -f -L 9999:localhost:8888 <SUNetID>@wheat<XX>.stanford.edu
```

### Login to jupyter in your browser, open a Terminal window, and activate the single cell environment
```bash
conda activate singlecell
```

### Pull the GitHub repository
```bash
cd ~/BIOC281
git pull
```

### Navigate to the current tutorial folder
```bash
cd ~/BIOC281/Classes/2
```

## Build a STAR index with and without masking psuedoautosomal regions (PAR)

### Download NCBI Refseq sequences of human chromosomes 1, 2, X, Y, and M from UCSC
```bash
mkdir -p RefSeq_Oct2020/{Annotation,GenomeFasta,GenomeIndex/STARIndex/{100bp_PRMSK,100bp_YMSK}}

cd RefSeq_Oct2020/GenomeFasta && \
wget --quiet https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz && \
wget --quiet https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr2.fa.gz && \
wget --quiet https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrX.fa.gz && \
wget --quiet https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrY.fa.gz && \
wget --quiet https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz
```

### Concatenate the sequences to create a minigenome
```bash
cat chr1.fa.gz chr2.fa.gz chrX.fa.gz chrY.fa.gz chrM.fa.gz > hg38.RefSeq.mini.fa.gz
```
### Decompress the minigenome to get genome fasta file
```bash
gzip -d hg38.RefSeq.mini.fa.gz
```

### Many single cell sequencing experiments add ERCC spike-ins during library prep. Download the ERCC fasta and annotation from ThermoFisher
```bash
wget --quiet https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip
unzip ERCC92.zip
mv ERCC92.gtf ../Annotation/
```

### create a ERCC version of the minigenome
```bash
cat hg38.RefSeq.mini.fa ERCC92.fa > hg38.RefSeq.mini.ERCC.fa
```

### Create a hard PAR mask version of genome as per the coordinates described in Mangs et al (2007) Current Genomics paper
```bash
echo -e "chrY\t10001\t2781479\nchrY\t56887903\t57217415" > PAR1_2.bed
echo -e "chrY\t0\t57227415" > chrY.bed

bedtools maskfasta -fi hg38.RefSeq.mini.ERCC.fa -bed PAR1_2.bed -fo hg38.RefSeq.mini.PARYhMSK.ERCC.fa
bedtools maskfasta -fi hg38.RefSeq.mini.ERCC.fa -bed chrY.bed -fo hg38.RefSeq.mini.YhMSK.ERCC.fa
```

### To explore the newly created genomes efficintly convert the fasta to 2bit format
```bash
faToTwoBit hg38.RefSeq.mini.PARYhMSK.ERCC.fa hg38.RefSeq.mini.PARYhMSK.ERCC.2bit
faToTwoBit hg38.RefSeq.mini.YhMSK.ERCC.fa hg38.RefSeq.mini.YhMSK.ERCC.2bit
```

### Now explore the minigenome and also check that the PAR regions are correctly masked. We would try to extract the PAR region from the genome fasta such that the resulting sequence contains 20 base pairs upstream and 20 base pairs downstream of the PAR regions.
```bash
twoBitToFa hg38.RefSeq.mini.PARYhMSK.ERCC.2bit par1.fa -seq=chrY -start=9981 -end=2781499
twoBitToFa hg38.RefSeq.mini.PARYhMSK.ERCC.2bit par2.fa -seq=chrY -start=56887883 -end=57217435

less par1.fa
less par2.fa
```

### Test the original before masking
```bash
faToTwoBit hg38.RefSeq.mini.ERCC.fa hg38.RefSeq.mini.ERCC.2bit
twoBitToFa hg38.RefSeq.mini.ERCC.2bit par1_normal.fa -seq=chrY -start=9981 -end=2781499
twoBitToFa hg38.RefSeq.mini.ERCC.2bit par2_normal.fa -seq=chrY -start=56887883 -end=57217435

less par1_normal.fa
less par2_normal.fa
```

### Download the Refgene Annotation GTF file. Correct way to download refGene GTF is by using a tool called genePredToGtf from UCSC 
```bash
cd ../Annotation/
genePredToGtf hg38 refGene hg38.refGene.gtf
```

### Create ERCC version of the Annotation GTF
```bash
cat hg38.refGene.gtf ERCC92.gtf > hg38.refGene.ERCC.gtf
```

### You can use less to check and manually read and explore all the fasta and GTF file. Notice that less can also read the compresse gzipped (xxx.gz) files as well.
```bash
less chrM.fa.gz
less ERCC92.fa
less hg38.refGene.ERCC.gtf
```

### Build the STAR genome indices. We need both the genome fasta and the annotation GTF for this purpose
```bash
cd ../GenomeIndex/STARIndex/100bp_PRMSK/
ln -sv ../../../GenomeFasta/hg38.RefSeq.mini.PARYhMSK.ERCC.fa
ln -sv ../../../Annotation/hg38.refGene.ERCC.gtf
```
#### Create a STARIndex script. We have created text files using echo earlier now we will learn how we can use cat to do the same. This is also the same command we used to concatenate compressed or uncompressed files. You can inspect the script by less as shown above or edit using vim/vi or your favorite text editor.
```bash
cat > ./createPARYMSKSTARIndex_scr << EOF
#!/bin/bash

STAR --runThreadN 4 --runMode genomeGenerate \
--genomeDir ./ --genomeFastaFiles ./hg38.RefSeq.mini.PARYhMSK.ERCC.fa \
--sjdbGTFfile ./hg38.refGene.ERCC.gtf \
--sjdbOverhang 100
EOF

less createPARYMSKSTARIndex_scr
```

#### Now run the bash script to create STAR Index. 
```bash
bash createPARYMSKSTARIndex_scr
```

### While index generation is going on (approx ~10 min), open new rice session
```bash
ssh <SUNetID>@stanford.edu
```

### From the new rice session connect to the same whaeat node that was allocated to you earlier and is running STAR to generate Index and activate "singlecell". this way you can run multiple terminal sessions to work concurrently.
```bash
ssh wheat<XX>@stanford.edu
conda activate singlecell
```

### Change to the appropriate directory and create script for creating Index now for the Complete Y masked minigenome. Notice that the Annotation remains the same for different versions of the masked genomes
```bash
cd ~/BIOC281/Classes/2/RefSeq_Oct2020/GenomeIndex/STARIndex/100bp_YMSK

ln -sv ../../../GenomeFasta/hg38.RefSeq.mini.YhMSK.ERCC.fa
ln -sv ../../../Annotation/hg38.refGene.ERCC.gtf

cat > ./createYMSKSTARIndex_scr << EOF
STAR --runThreadN 4 --runMode genomeGenerate \
--genomeDir ./ --genomeFastaFiles ./hg38.RefSeq.mini.YhMSK.ERCC.fa \
--sjdbGTFfile ./hg38.refGene.ERCC.gtf \
--sjdbOverhang 100
EOF
```

### Wait for the first Index generation in the previos terminal to finish before running your new script. "STAR --runmode genomegenerate" takes 21-25G of RAM even with the minigenome and we requested only 30G as memory for our use.
```bash
bash createYMSKSTARIndex_scr
```

### Move back to your first terminal and prepare Index for 10x Cellranger, which also uses STAR. Create directories and softlinks for appropriate gennome fasta and annotation GTF
bash```
cd ~/BIOC281/Classes/2/RefSeq_Oct2020/GenomeIndex/
mkdir 10xIndex && cd 10xIndex
ln -sv ../../GenomeFasta/hg38.RefSeq.mini.PARYhMSK.ERCC.fa
ln -sv ../../Annotation/hg38.refGene.ERCC.gtf
```

#### First clean up annotation GTF to meet Cell Ranger requirements
```bash
cellranger mkgtf ./hg38.refGene.ERCC.gtf ./hg38.refGene.ERCC.filtered.gtf --attribute=gene_biotype:protein_coding
```

### Now we  use the filtered annotation to prepare 10X reference genome Index
```bash
cellranger mkref --genome=hg38.PARYhMSK.ERCC --fasta=hg38.RefSeq.mini.PARYhMSK.ERCC.fa --genes=hg38.refGene.ERCC.filtered.gtf
```
### You should see the following output indicating how to use reference index while mapping reads generated after a 10x sequencing run
```bash
>>> Reference successfully created! <<<

You can now specify this reference on the command line:
cellranger --transcriptome=/home/sinhar/BIOC281/Classes/2/RefSeq_Oct2020/GenomeIndex/10xIndex/hg38.PARYhMSK.ERCC ...
```

### Build a kallisto index
#### Either open a new Terminal window while 10X STAR Index is being built and proceed as below: first download the whole genome fasta then build transcriptome fasta using hg38 genome fasta and annotation and finally clean the messy transcript names with a perl oneliner "https://github.com/trinityrnaseq/Griffithlab_rnaseq_tutorial_wiki/blob/master/Kallisto.md"
```bash
cd ~/BIOC281/Classes/2/RefSeq_Oct2020/GenomeFasta
wget --quiet http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit
twoBitToFa ./hg38.2bit ./hg38.RefSeq.fa
cat ./hg38.RefSeq.fa ./ERCC92.fa > hg38.RefSeq.ERCC.fa

mkdir -p ~/BIOC281/Classes/2/RefSeq_Oct2020/GenomeIndex/kallistoIndex
cd ~/BIOC281/Classes/2/RefSeq_Oct2020/GenomeIndex/kallistoIndex

ln -sv ../../GenomeFasta/hg38.RefSeq.ERCC.fa
ln -sv ../../Annotation/hg38.refGene.ERCC.gtf

gtf_to_fasta ./hg38.refGene.ERCC.gtf ./hg38.RefSeq.ERCC.fa hg38.refGene.ERCC.transcriptome.fa

cat hg38.refGene.ERCC.transcriptome.fa | perl -ne 'if ($_ =~/^\>\d+\s+\w+\s+(ERCC\S+)[\+\-]/){print ">$1\n"}elsif($_ =~ /\d+\s+([A-Z]+_\d+_?\d+?)/){print ">$1\n"}else{print $_}' > hg38.refGene.ERCC.transcriptome.clean.fa

less hg38.refGene.ERCC.transcriptome.fa
less hg38.refGene.ERCC.transcriptome.clean.fa
```

### Now build kallisto transcriptome index
```bash
kallisto index -i hg38.refSeq.transcriptome.ERCC.idx ./hg38.refGene.ERCC.transcriptome.clean.fa --make-unique
```

## Common theme: all mappers/aligners of NGS data (fastq) require a step where you create a genome index which speeds up the mapping process. Some mappers like STAR also insert annotation information in the genome index so while mapping the program is aware of known (user provided) splice-junctions-- further speeding up the mapping/alignment process. We encourage you to look up other popular mappers and their genome index generation process like bowtie2, tophat, hisat2, bwa, all genome based mappers, and salmon, which is another transcriptome based mapper similar to kallisto.

# Mapping fastq reads using STAR
### The script below is adapted from ENCODE long-mRNA protocol
