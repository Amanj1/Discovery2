# Discovery2
A pipeline that identifies common and uncommon variants of microorganisms and viruses. 
This is possible by translating nucleotide sequences into 6 frame protein sequences and matching the sequences against different databases. 
It also produces a classification of each sequence and an interactive table of all results. 

We also have other results that I have not worked with or developed further. 
In addition, a list of sequences that do not match anything in the databases used in the pipeline is produced for each sample.

This pipeline is divided into two parallel runs. One is for reads analysis and other are after assembly data is produced and analyzed.

![alt text](/UML_diagram/discovery2.png)


The UML (Unified Modeling Language) diagram displays the steps in the pipeline. "Preprocessing reads" is the input data and the data should be trimmed, quality checked and human DNA/reads should be filtered out before running the pipeline. Everything below labeled "tax reads" are reads based analysis for matching sequences using different methods. The rest are analysis based on assembly data. Where we include different detection method. In this pipeline I have focused mostly on extending the Diamond blastx run with a classifier and creating HTML interactive tables and tsv files. Everything below "tax contigs diamond" are processes that are extensions for the Diamond blastx analysis to make it more presentable for the user.

## Software requirements 
 Software required to run all processes in the pipeline.
 - [Nextflow DSL1](https://www.nextflow.io/)
 - [Python3](https://www.python.org/downloads/)
 - [megahit](https://github.com/voutcn/megahit)
 - [seqtk](https://github.com/lh3/seqtk)
 - [BBmap tools](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/)
    - bbwrap.sh
    - pileup.sh
 - [samtools](http://www.htslib.org/)
 - [metaphlan2](http://huttenhower.sph.harvard.edu/metaphlan2)
 - [kraken2](https://ccb.jhu.edu/software/kraken2/)
 - [FastViromeExplorer](https://fastviromeexplorer.readthedocs.io/en/latest/)
 - [kallisto](https://github.com/pachterlab/kallisto)
 - [diamond blastx](https://github.com/bbuchfink/diamond)
 - [R](https://www.r-project.org/)
 - [virfinder](https://github.com/jessieren/VirFinder)
 - [FragGeneScan](https://omics.informatics.indiana.edu/FragGeneScan/)
 - [BLAST blastdbcmd](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

#### Software versions currently tested on
| Software   | Version |
| --------   | ------- |
| Nextflow   | 19.07.0.5106 |
| Python3    | 3.6.9   |
| Seqtk      | 1.3-r106|
| Megahit    | v1.2.8  |
| BBmap      | 38.68   |
| Samtools   | 1.9     |
| MetaPhlAn2 | 2.7.7   |
| Kraken2    | 2.0.8-beta |
| FastViromeExplorer | 1.3 |
| Kallisto    | v0.43.1 |
| Diamond Blastx | 2.0.6 |
| R           | 3.6.1   |
| virfinder   | 1.1     |
| FragGeneScan| 1.31    |
| BLAST       | 2.5.0+  |

## Database requirements 

You will need to download several databses to be able to run the pipeline for five softwares mentioned below.

#### 1. Kraken2 
Kraken2 requires to download a database that consist of viruses, bacteria and fungi (based on what we worked on).  
How to build database for kraken2:
```
kraken2-build --download-taxonomy --db <database name>
kraken2-build --download-library bacteria --db <database name>
kraken2-build --download-library viral --db <database name>
kraken2-build --download-library fungi --db <database name>
kraken2-build --build --db <database name> --threads <number of threads>
```
The name you choose for the database should be exactly the same for all the commands above. If you want to build a different database you can open the kraken2 link in "software requirements".

#### 2. FastViromeExplorer
FastViromeExplorer requries a viruslist and a kallisto index.

#### 3. Diamond blastx
Diamond requries a nr database in dmnd format, a protein accession to taxid and taxonomy nodes. 

#### 4. blastdbcmd
Requires a BLAST databse for latest nr sequences (all sequences).

#### 5. contig classifer
Requries a complete taxonomic lineage to matched sequences in the Diamond blastx results










