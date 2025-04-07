# Discovery2
## Work in Progress
I'm currently developing a new and improved version of Discovery! This update will include:

    More advanced filtering options in the final results

    Enhanced result handling

    Improved performance and faster run times

    Multi-threading support for processes that previously ran single-threaded

    More complete and accurate taxonomic data

    DS1 -> DSL2 Nextflow

Once the new version is complete, I’ll share a link to the updated repository here.

A metagenomic pipeline that identifies common and uncommon variants of microorganisms and viruses. 
We find the most common by using different nucleotide-based searches/tools against databases and the uncommon by translating nucleotide sequences into 6 frame protein sequences and matching the sequences against the entire nr database. 
The pipeline also produces a classification of each sequence and an interactive table of all results. In addition, a list of sequences that do not match with anything in the databases used in the pipeline is produced for each sample.

This pipeline is divided into two parallel runs. One is for reads analysis and other are for analysis of assembly data.

![alt text](/UML_diagram/discovery2.png)


The UML (Unified Modeling Language) diagram displays the steps in the pipeline. "Preprocessing reads" is the input data and the data should be trimmed, quality checked and human DNA/reads should be filtered out before running the pipeline. Everything below labeled "tax reads" are reads based analysis for matching sequences using different methods. The rest are analysis based on assembly data. Where we include different detection method. In this pipeline I have focused mostly on extending the Diamond blastx run with a classifier and creating HTML interactive tables and tsv files. Everything below "tax contigs diamond" are processes that are extensions for the Diamond blastx analysis to make it more presentable for the user.

Below we can see one example of a table produced by one of the last processes. Please note that "aln_len" (alignment length) is based on aminoacid length and "contig_len" is based on the nucleotide sequence length.

![alt text](/UML_diagram/discovery2_html_table.png)

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
| FragGeneScan| 1.31    |
| BLAST       | 2.5.0+  |

| R-package  | Version |
| --------   | ------- |
| virfinder  | 1.1     |

## Database requirements 

You will need to download several databses to be able to run the pipeline for five softwares mentioned below. It may be a good idea to create a folder for each database.

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
You can read about how to install the database in the link provided in "Software requirements".

#### 3. Diamond blastx
Diamond requries a nr database in dmnd format, a protein accession to taxid and taxonomy nodes. 
First you will need to download a compressed file taxdump.gz.tar from this website.

```
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump.tar.gz
```
You will also need to protein accession to taxid in the link below.

```
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
```

When running Diamond blastx we are using the nr-protein (Non-redundant protein sequences) database. The link below is where you can download the fasta version of nr.
```
https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
```
or download it through terminal.
```
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
```
When you have the nr file you will need to run the following command below.
```
diamond makedb --in nr.gz -d nr
```
Final step to include taxnodes in the database. "pathwayToFile" should be the full path to the file.
```
diamond getseq --db diamond/nr.dmnd | diamond makedb --db nr.dmnd --taxonmap pathwayToFile/prot.accession2taxid --taxonnodes pathwayToFile/nodes.dmp
```

#### 4. blastdbcmd
Requires a BLAST databse for latest nr sequences (all sequences). In this step we will use the same nr.fasta we downloaded for 3.Diamond.
Run the following command below and replace <reference.fasta> with the nr fasta file.

```
makeblastdb -in <reference.fasta> -dbtype prot -parse_seqids -out nr -title "nr"
```

#### 5. contig classifer
Requries a complete taxonomic lineage to matched sequences in the Diamond blastx results. In the compressed file "new_taxdumps" that we downloaded for Diamond already contains the "fullnamelineage.dmp" for finding complete lineages for Diamond blastx results.

## Configuration file
In our configuration file "discovery2.config" we can add a profile. Each profile can be customized according to the available resources of your system. You can create your own profile through the following example below which you can find in "discovery2.config" file.

```
profiles {
  hamlet {
    includeConfig 'conf/hamlet.discovery.config'
  }

/*
  bianca {
    includeConfig 'conf/bianca.discovery.config'
  }
*/

 othello {
  includeConfig 'conf/othello.discovery.config'
 }

}
```
Above we can see three profiles. One of them is included in the repository created for Hamlet. You can follow the exemples below on how to modify hamlet and create youre own file and include this later in "discovery2.config" file.

In the configuration file "discovery2.config" we also have parameters that are used during the pipeline run.

```
/* Pipeline running parameters */
params{
  project_id='P12653'
  min_ctg_size=200
  publish_base_dir='discovery'
  html_dir='input_html'
  fastq_dir='preprocessing'
}
```
"min_ctg_size" is value for the minimum nucleotide sequence length after assemble run. The length filters out sequence size below the mimum value set in "min_ctg_size". fastq_dir is the input folder for your "fq.gz" data. The "project_id" should be changed to the label for your project.

### Database pathways 
In our configuration file found in conf/hamlet.discovery.config we will need to add all the paths for the databases.
Below you can see files and folders you will need to add full path to. 

```
 //1. FastViromeExplorer
 FVE_index='/PathToFile/FastViromeExplorer/ncbi-virus-kallisto-index-k31.idx'
 FVE_viruslist='/PathToFile/FastViromeExplorer/1.3/ncbi-viruses-list.txt'
 //2. Kraken2
 kraken2_db='/PathToFolder/kraken2/210205_bacvirfun'
 //3. Diamond databases
 diamond_db='/PathToFile/diamond/nr.dmnd'
 diamond_taxonmap='/PathToFile/ncbi_taxonomy/prot.accession2taxid
 diamond_taxonnodes='/PathToFile/ncbi_taxonomy/nodes.dmp'
 ncbi_full_lineage='/PathToFile/ncbi_taxonomy/fullnamelineage.dmp
```
If you scroll down you can find the process for contig classification here you will need to add the path for the nr database that you created in "4. blastdbcmd".
```
    withName: tax_contigs_classifier{
beforeScript='export BLASTDB=/PathToFolder/nr'
    }
}
```

### Software pathways 
In our configuration file found in conf/hamlet.discovery.config we will need to add all the paths for all softwares if they are not included in your PATH variable or conda environment.

Below you can see an example of a process which requires a path the BBmap tools. The default value for number of threads are 1 in each process if you want to increase you will need to add "cpus = 8" to increase it and in this case it is increased to 8 threads. You will need to do this for all processes if they are not part of your current PATH variable or conda environment. All processes are executed as a child process and adding the line "beforeScript='export PATH="/PathToFolder/tools/BBMap/38.68/bbmap:$PATH"'" does not save it to your PATH variable. It is only temprorary during the pipeline run. 
```
    withName: asm_map_reads_to_contigs{
		beforeScript='export PATH="/PathToFolder/tools/BBMap/38.68/bbmap:$PATH"'
        cpus = 8
    }
```

## Running pipeline
The user should create a folder one called 'preprocessing' and store all samples there. Each sample should have a folder labeled with sample name / samle ID. Each sample folder should contain paired-end fastq compressed gzip (GNU Zip) files and a file with unpaired reads in fastq gzip format. Each fastq file need to be labled with same sample names but different file extension names "_1.fq.gz", "_2.fq.gz" and "_unpaired.fq.gz" for the pipeline to recognize the sample.

If you prefer "fastq.gz" over "fq.gz" you can change in the line below and this line can be found in begning of "disovery2.nf" file.

```
fastq_files = Channel.fromFilePairs("${params.fastq_dir}/**/*_{1,2,unpaired}.fq.gz",size:3)
```

To run the pipeline in command line:
```
nextflow -C discovery2.config run discovery2.nr -profile hamlet
```
To run the pipeline in command line and resume from cache memory:
```
nextflow -C discovery2.config run discovery2.nr -profile hamlet -resume
```





