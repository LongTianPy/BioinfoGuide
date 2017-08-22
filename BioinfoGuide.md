## Downloading genome assemblies from NCBI GenBank
Use a function that provided by the ANI calculation package [`pyani`](https://github.com/widdowquinn/pyani#script-genbank_get_genomes_by_taxonpy). Visit the link for detailed explanation of all arguments. But for lab use:

Go to [NCBI-taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy), enter the taxon of interest in the textbox, depends on the level or rank of the genome assemblies you want, you can enter family name, genus name, species name, etc. For this example, we want all genome assemblies of `Ralstonia solanacearum`.

Click on the link in bold.

Then end of the url of the following page contains the taxon ID of this taxon. For `Ralstonia solanacearum`, it is `305`.

For the script that we are going to execute, `-t` is the option for `taxon ID`, so, to download `Ralstonia solanacearum` genomes, use the command:
`python3 ~/Scripts/pyani/genbank_get_genomes_by_taxon.py -o Ralstonia_solanacearum_genomes -t 305 -v --email this@doesnt.matter`

where `-o` indicates the targeted folder to store downloaded genomes (it will be created if does not exist), `-t` the taxon ID, `-v` be the flag of information when accessing NCBI and processing genome sequences.

## Assembling genomes with SPAdes

Since our lab mainly deals with bacterial genomes, it seems that [SPAdes](http://bioinf.spbau.ru/spades) is a good assembler after I have compared it with other assemblers.

There is no one-click way of assembling genome from raw reads with the existence of quality control. For example, raw reads are in the folder `~/Data/raw_reads`.

Usually, paired-end reads are named with `R1` and `R2` to separate them, merging reads from the same end is the first step, so in the directory:
```bash
cd ~/Data/raw_reads
cat *R1* > merged_R1.fastq
cat *R2* > merged_R2.fastq
```

Check the quality of the raw reads by launching FastQC:
```bash
cd ~/Software/FastQC
./fastqc
```

This will open the GUI of FastQC. Open the
