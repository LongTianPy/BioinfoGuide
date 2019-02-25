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

This will open the GUI of FastQC. Open one of the merged reads, and check the report after waiting for a while.

Things to look at in the report:
`sequence length` in the basic statistics tab, `per base sequence quality` for the distribution of quality, `per base sequence content` for the approximate length of adapters.

Then trim the raw reads according to the quality report. There is no strict values of quality control, but adapter sequences have to be removed, quality score at least 30, and keep a reasonable length of raw reads so that it can be assembled.

[`Trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic) is used to trim the reads, for example:
```bash
cd ~/Data/raw_reads
java -jar \
~/Software//home/vinatzerlab/Software/Trimmomatic-0.32/trimmomatic-0.32.jar \
PE -phred33 merged_R1.fastq merged_R2.fastq \
trimmed_R1_paired.fq trimmed_R1_unpaired.fq \
trimmed_R2_paired.fq trimmed_R2_unpaired.fq \
HEADCROP:12 SLIDINGWINDOW:4:30 MINLEN:160
```

After trimmed the reads and only keep the high-quality ones, we can proceed to assemble them into genome.
SPAdes takes a range of k-mer lengths.
The current SPAdes version as of Feb 25, 2019 is 3.11.0, its scripts are located in `~/Software/SPAdes-3.13.0-Linux/bin/`
Detailed document can be viewed at (`https://github.com/ablab/spades#spades-command-line-options`).
But for paired-end Illumina MiSeq, the following should be enough:
```
~/Software/SPAdes-3.13.0-Linux/bin/spades.py -o <output_dir> --pe-1 trimmed_R1_paired.fq --pe-2 trimmed_R2_paired.fq -t 4 -m 30 -k 51,61,71,81,91,101,111
```
Then the software will assemble the genome with different k and decide the one with the optimal quality.


## BIGSdb installation information
BIGSdb is installed on the new server

BIGSdb documentation: (http://bigsdb.readthedocs.org/en/latest/ )

The library files are in `/usr/local/lib/BIGSdb/`   including the modified GenomeComparator.pm file (for calculating Tajima's D). Save this file separately if you are going to update BIGSdb.


The configuration files are in `/etc/bigsdb/`


The databases are in `/etc/bigsdb/dbases/`   
 and additional information needs to be added to the MySQL database for each new database created.

BIGSdb logs are in `/var/log/bigsdb/`

The website interface is installed in `/var/www/cgi-bin/bigsdb/`

Main website: http://128.173.74.68
The website addresses for the RAINS (Pantoea/Erwinia + NCBI relatives) database are:http://128.173.74.68/cgi-bin/bigsdb/bigsdb.pl?db=bigsdb_pantoea_isolates
http://128.173.74.68/cgi-bin/bigsdb/bigsdb.pl?db=bigsdb_pantoea_seqdef

the curator website addresses for the RAINS database (where you can modify/add to the databases)
http://128.173.74.68/cgi-bin/bigsdb/bigscurate.pl?db=bigsdb_pantoea_isolates
http://128.173.74.68/cgi-bin/bigsdb/bigscurate.pl?db=bigsdb_pantoea_seqdef  
Please ask Dr. Vinatzer for login info.

Adding genomes to BIGSdb
If you used the above script for downloading genomes from NCBI, tables containing most of the information needed to add genomes to BIGSdb have been automatically created. These tables must be modified before importing to BIGSdb, specifically several columns need to be added so that the table contains the following columns:
`isolate aliases references      strain  genus   species complete_name   location        sample  isolation_date  assembly_level  genome_representation   download_database       genbank_assembly_id     wgs_project     date_sequenced  total_sequence_length   number_of_contigs       number_of_plasmids      genbank_id_plasmids     comments`
If you assembled genomes manually or using the genome assembly script, you need to make a table with the above mentioned columns. The genome_assembly_pipeline.pl script creates the barebones of a table, however much of the information (such as genus species names) needs to be added manually.

You first have to add the isolate information, and then sequences may be added afterwards separately. The isolate information can be added all at once in a batch, however sequences must be added one isolate at time.

#Core and Pan genome analyses using GET_HOMOLOGUES
To compute the core and pan-genomes for a given set of genomes, the GET_HOMOLOGUES scripts can be used to automate a lot of the analyses.
Manual and examples: http://eead-csic-compbio.github.io/get_homologues/manual/

First annotate all of your genomes using RAST or prokka, and download them in genbank format. Place all of the genbank files for the genomes you wish to compare in a folder. The get_homologues script can use either orthoMCL or COGTriangles to calculate groups of genes. If you want to add more genomes later, just add them to the folder and get_homologues will reuse previous blast searches in order to save time. To run get homologues using orthoMCL, run the following in the directory containing your folder of genomes (genomefolder in this example):

`~/Software/get_homologues-x86_64-20140930/get_homologues.sh -d ./genomefolder -n 10 -M -t 0 > gethomologs.log`

In this case, `-n 10` means use 10 threads to speed up the analysis, `-M` means use orthoMCL, `-o` is what you want to call the output folder, and `-t` is the minimum size of clusters to report, in this case 0 means all clusters, which is necessary for pangenome analyses. If you want to also run the analyses using COGTriangles, replace `-M` with `-G`. This will reuse the blast searches and create a separate folder containing the COGTriangles analysis, which can be compared with the orthoMCL run later. The output files will be in the folder genomefolder_homologues. The next step is to run the `compare_clusters` script three times, which compares runs between orthoMCL and COGTriangles and also prepares matrices for the core and pan-genome analyses for individual and combined analyses:

```
~/Software/get_homologues-x86_64-20140930/compare_clusters.pl -d ./genomefolder_homologues/genomefolder_OMCL -o ./genomefolder_homologues/genomefolder_OMCL -m -T > compareclusters.out

~/Software/get_homologues-x86_64-20140930/compare_clusters.pl -d ./genomefolder/genomefolder_homologues/genomefolder_COG -o 
./genomefolder/genomefolder_homologues/genomefolder_COG -m -T > compareclusters.out
```

```
~/Software/get_homologues-x86_64-20140930/compare_clusters.pl -d ./genomefolder/genomefolder_homologues/genomefolder_OMCL,./genomefolder/genomefolder_homologues/genomefolder_COG -o ./genomefolder/genomefolder_homologues/genomefolder_OMCL-COG -m -T > compareclusters.out
```

In this case, -d is the input file folders containing the orthoMCL and COG analyses, separated by a comma with no spaces between. -o is the output directory, -m produces a pangenome matrix for the intersection of the orthoMCL and COG analyses, and -T makes a parsimony based phylogeny (based on presence or absence of genes, not the sequences themselves). 

The next step is to run the parse_pangenome_matrix script to calculate core and pan genomes for the individual analyses along with the combined orthoMCL/COGTriangles:

```
~/Software/get_homologues-x86_64-20140930/parse_pangenome_matrix.pl -m ./genomefolder/genomefolder_homologues/genomefolder_OMCL/pangenome_matrix_t0.tab -s 
```

```~/Software/get_homologues-x86_64-20140930/parse_pangenome_matrix.pl -m ./genomefolder/genomefolder_homologues/genomefolder_COG/pangenome_matrix_t0.tab -s 
```

```
~/Software/get_homologues-x86_64-20140930/parse_pangenome_matrix.pl -m ./genomefolder/genomefolder_homologues/genomefolder_OMCL-COG/pangenome_matrix_t0.tab -s 
```

You can also compare two sets of genomes from your analysis to identify if there are any differences in gene content between them. Just make two text files, A.txt and B.txt, each containing a list of the genome genbank files you want to compare. The script will find all genes present in A that are not in B:

```
~/Software/get_homologues-x86_64-20140930/parse_pangenome_matrix.pl -m ./genomefolder/genomefolder_homologues/genomefolder_OMCL-COG/pangenome_matrix_t0.tab -A A.txt -B B.txt -g 
```

To make the script a little more lenient in deciding what is included in A and not B, you can use the option -P, for example use -P 90 to consider genes present in 90% of A and not present in 90% of B. This is good when you have draft genomes.

```
~/Software/get_homologues-x86_64-20140930/parse_pangenome_matrix.pl -m ./genomefolder/genomefolder_homologues/genomefolder_OMCL-COG/pangenome_matrix_t0.tab -A A.txt -B B.txt -g -P 90
```
## LIN code assignment
Visit LINbase.org

## To copy files to the server

## List of software installed on the server

## Phylogenetic analysis with RAxML

Statistical model-based methods of phylogenetic tree construction (such as maximum likelihood or Bayesian methods) are generally better than parsimony or distance based ones. They rely on observed models of DNA or amino acid (AA) evolution rather than just the number of differences between sequences. The trees generated by these methods can also be statistically compared since they are based on models. RAxML and GARLI are two programs that perform maximum likelihood (ML) analyses and are well represented in the literature. MrBayes is a good program for Bayesian tree construction but it is very difficult to use and interpret the results.

To make a phylogeny with RAxML, several things have to be done first. 

### Aligning sequences

Use MAFFT (mafft dataset.fasta > dataset.aln) to align your sequences, then if you have a DNA alignment you need to run jModelTest, and if you have a protein alignment you need to run ProtTest. These will check for the “best” DNA or AA substitution rate model of evolution to use in the analysis. If you have a SNP alignment you need to use a model taking ascertainment bias into account, in RAxML these models start with ASC_ (example ASC_GTR is GTR model with ascertainment bias). The best model might also include gamma-distributed rate heterogeneity (GAMMA or model+G) or DNA/AA frequencies fixed or inferred from the dataset (model+F) rather than from the model; or finally an inferred proportion of invariant sites (model+I). 

For example, the WAG protein model with gamma distributed rate heterogeneity is WAG+G, or in RAxML PROTGAMMAWAG; with dataset frequences (+F) it would be PROTGAMMAWAGF.

Go with the best model under the AIC or BIC criterion, if the log-likelihood (lnL) of a model and that model with gamma correction is equal or close, go with the simpler model (ie GTR vs GTRGAMMA). The model with the smallest (ie least negative and closest to 0) lnL is usually the best one to use for that dataset.

To run jModelTest for DNA alignments on the server, run it as follows:

```
cd ~/Software/jmodeltest-2.1.6/
./runjmodeltest-gui.sh
```

In the window that pops up, click File->Load DNA alignment and load your alignment. Then click Analysis->Compute likelihood scores. When this is finished click Analysis->Do AIC calculations... and Analysis->Do BIC calculations...  

Finally click Results->Show results table to see the results. The model at the top of each list is usually the best one, but it may differ between analyses.

To run ProtTest for protein alignments on the server, run it as follows:

```
cd ~/Software/prottest-3.4.-20140123/
./runXProtTestHPC.sh
```

and do the same thing as for jModelTest.

Now that you have a model for your dataset, its time to run RAxML. You have to run it twice, first to find the best ML tree, and second to calculate bootstrap replicates.

To find the best ML tree out of 20 searches, run RAxML as follows:

```
~/Software/standard-RaxML-master/raxmlHPC-AVX -p 12345 -m GTRGAMMA -s dataset.aln -n datasetname -N 20 > dataset_raxml.out
```

To run 1000 bootstrap replicates, run RAxML as follows:

```
~/Software/standard-RaxML-master/raxmlHPC-AVX -p 12345 -b 12345 -m GTRGAMMA -s dataset.aln -n datasetname_boot -N 100 > datasetbootstrap_raxml.out
```

To map bootstrap support to branches of the best ML tree as a percentage of 1000 bootstrap replicates, first run the previous two commands and then run RAxML as follows:

```
~/Software/standard-RaxML-master/raxmlHPC-AVX -p 12345 -m GTRGAMMA -f b -t RAxML_bestTree.datasetname -z RAxML_bootstrap.datasetname_boot -n datasetname_raxml.tre
```

which produces the file:  
`RAxML_bipartitions.datasetname_raxml.tre`  
 which contains the best tree with bootstrap support mapped to branches. To collapse branches with less than 50% bootstrap support (or any other cutoff), open the tree with FigTree and save it. Then run the program TreeCollapseCL4:

```
java -jar ~/Software/TreeCollapseCL4.jar -t F -b 50 -f RAxML_bipartitions.datasetname_raxml.tre
```

which creates the file:`RAxML_bipartitions_50coll.datasetname_raxml.tre`

If you want to use multiple cpus to speed up the analysis (up to 16), instead use: raxml-PTHREADS-AVX and include -T 16 to run on 16 cpus.
Example:

```
~/Software/standard-RaxML-master/raxmlHPC-PTHREADS-AVX -T 16 -p 12345 -m GTRGAMMA -s dataset.aln -n datasetname -N 20 > dataset_raxml.out
```

If you want to run it on a cluster, use raxml-MPI-AVX with mpirun.


## Search and replace inside a file using a table of search and replace terms

