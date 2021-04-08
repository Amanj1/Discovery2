#!/usr/bin/env nextflow

/*
How to run:
nextflow -C discovery.nf.config run discovery.nf --fastq_files preprocessing -profile hamlet
*/

params.html_dir='input_html'
params.fastq_dir='preprocessing'
fastq_files = Channel.fromFilePairs("${params.fastq_dir}/**/*_{1,2,unpaired}.fq.gz",size:3)
html_files = Channel.fromFilePairs("${params.html_dir}/html_{start,end}.txt",size:2)
params.project_id='P12653'

fastq_files.into{
  asm_megahit_in;
  asm_metaspades_in;
  reads_to_map;
  tax_reads_metaphlan2_in;
  tax_reads_kraken2_in;
  tax_reads_FastViromeExplorer_in;
  print_reads
 }

html_files.into{html_each;
html_all;
html_simplified_all}

/*print_reads.println()*/

/**
ASSEMBLY Module

TODO: Decide if run megahit in meta-sensitive mode or in default?
TODO: Run also SPAdes standard pipeline or only MetaSpades?
**/

process asm_megahit{
  tag { "${sample_id}" }
  publishDir "${params.publish_base_dir}/${sample_id}/megahit", mode: 'copy', pattern: "1_assembly"

  input:
  set sample_id, reads from asm_megahit_in

  output:
  set sample_id,val('megahit'), "contigs.fa" optional true into asm_megahit_out
  file "1_assembly" optional true into asm_megahit_dir_out

  script:
  """
  megahit -t ${task.cpus} --presets meta-sensitive -1 ${reads[0]} -2 ${reads[1]} -r ${reads[2]} --cpu-only  -o 1_assembly
  if [ -s 1_assembly/final.contigs.fa ]; then cat 1_assembly/final.contigs.fa | sed -e 's/\\s\\+/,/g' > contigs.fa; fi
  """
}

/*
NOTE: MetaSPAdes does not support incorporating the single reads
TODO: figure out how to use a scratch folder? --tmp-dir opt
*/

process asm_metaspades{
  tag { "${sample_id}" }
  publishDir "${params.publish_base_dir}/${sample_id}/metaspades", mode: 'copy', pattern: "1_assembly"

  input:
  //set sample_id, reads from asm_metaspades_in
  set sample_id, reads from Channel.empty()

  output:
  set sample_id,val('metaspades'),"contigs.fa" optional true into asm_metaspades_out
  file "1_assembly" optional true into asm_metaspades_dir_out

  script:
  """
  spades.py --meta -t ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} -o 1_assembly -m ${params.max_spades_mem}
  if [ -s 1_assembly/contigs.fasta ]; then ln 1_assembly/contigs.fasta contigs.fa; fi
  """
}

/*
NOTE: Combine contigs into a single channels
*/
all_assemblies = asm_megahit_out.mix(asm_metaspades_out)

process asm_filter_contigs{
  tag { "${sample_id}/${assembler}" }
  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}/2_filt_contigs", mode:'link'

  input:
  set sample_id,assembler,"contigs.fa" from all_assemblies

  output:
  set sample_id,assembler,"contigs_filt.fa" into asm_filter_contigs_out

  script:
  """
  seqtk seq -L ${params.min_ctg_size} contigs.fa > contigs_filt.fa
  """
}

/*
Clone the contigs into different channels
*/

asm_filter_contigs_out.into{
  contigs_to_be_mapped;
  tax_contigs_diamond_in;
  tax_contigs_kraken2_in;
  tax_contigs_virfinder_in;
  tax_contigs_fgs_in;
  tax_contigs_virsorter_in
  /*TODO: COMPLETE THIS LIST */
}

/*
Combine reads and contigs for mapping
*/
asm_map_reads_to_contigs_in = reads_to_map.cross(contigs_to_be_mapped).map{ it[0]+it[1][1,2] }

/**
NOTE: bbwrap/bbmap parameters
  * kfilter is the minimum length of consecutive exact matches for an alignment (min kmer size of dBg assembler)
  * maxindel limits the indel size
  * subfilter limits the number of mismatches for an alignment
*/

process asm_map_reads_to_contigs{
  tag { "${sample_id}/${assembler}" }

  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}/2_filt_contigs", mode:'link'

  input:
  set sample_id, 'reads_*.fq.gz',assembler,"contigs_filt.fa" from asm_map_reads_to_contigs_in

  output:
  set sample_id, assembler,"reads_to_contigs.sam.gz" into asm_map_reads_to_contigs_out

  script:
  """
  bbwrap.sh ref=contigs_filt.fa in=reads_1.fq.gz,reads_3.fq.gz in2=reads_2.fq.gz,null out=reads_to_contigs.sam.gz usejni=t kfilter=22 subfilter=15 maxindel=80
  """
}

/* Send the mapped reads to
1) mapping stats with samtools flagstat
2) Calculate coverage with bbtools pileup.sh
*/

asm_map_reads_to_contigs_out.into{ asm_mapping_stats_in;
                                   asm_per_ctg_coverage_in }

process asm_mapping_stats{
  tag {"${sample_id}/${assembler}"}

  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}/2_filt_contigs", mode:'link'

  input:
  set sample_id, assembler,'reads_to_contigs.sam.gz' from asm_mapping_stats_in

  output:
  set sample_id, assembler, "${sample_id}_${assembler}_flagstat.txt" into asm_mapping_stats_out

  script:
  """
  samtools flagstat reads_to_contigs.sam.gz > ${sample_id}_${assembler}_flagstat.txt
  """
}

process asm_per_ctg_coverage{
  tag {"${sample_id}/${assembler}"}

  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}/2_filt_contigs", mode:'link'

  input:
  set sample_id, assembler,"reads_to_contigs.sam.gz" from asm_per_ctg_coverage_in

  output:
  set sample_id, assembler,"reads_to_contigs.cov.txt" into asm_per_ctg_coverage_out

  script:
  """
  pileup.sh in=reads_to_contigs.sam.gz out=reads_to_contigs.cov.txt
  """
}

/**
TAX ASSIGNMENT - READS
**/

process tax_reads_metaphlan2{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/reads", mode:'link'

  input:
  set sample_id, 'reads_*fq.gz' from tax_reads_metaphlan2_in

  output:
  set sample_id, "${sample_id}_metaphlan2.tsv" into tax_reads_metaphlan2_out

  script:
  """
  zcat *fq.gz > seq.fq
  metaphlan2.py --nproc ${task.cpus} --input_type multifastq --sample_id_key '#clade' \
  		--sample_id '${sample_id}' seq.fq ${sample_id}_metaphlan2.tsv
  rm seq.fq
  """
}

process tax_reads_kraken2{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/reads", mode:'link'

  input:
  set sample_id, 'reads_*.fq.gz' from tax_reads_kraken2_in

  output:
  set sample_id, "${sample_id}_kraken2.txt","${sample_id}_kraken2_report.txt" into tax_reads_kraken2_out
  set sample_id, "${sample_id}_kraken2_unmapped.fq" into tax_reads_unmapped_kraken2
  
  script:
  """
  kraken2 --db ${params.kraken2_db} --threads ${task.cpus} --output ${sample_id}_kraken2.txt \
    --unclassified-out ${sample_id}_kraken2_unmapped.fq \
    --report ${sample_id}_kraken2_report.txt --gzip-compressed reads_*.fq.gz
  """
}

tax_reads_unmapped_kraken2.into{unmapped_reads_kraken2_blastn; unmapped_reads_kraken2_blastx}

/**
TODO: Choose between  NCBI RefSeq or  IMG/VR database
**/

process tax_reads_FastViromeExplorer{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/reads", mode:'link'

  input:
  set sample_id, 'reads_*.fq.gz' from tax_reads_FastViromeExplorer_in

  output:
  file "${sample_id}_fastviromeexplorer.sam" into fve_out_1
  file "${sample_id}_fastviromeexplorer_abundance.tsv" into fve_out_2

  script:
  """
  java -cp \${FVE_PATH}/bin FastViromeExplorer -i ${params.FVE_index} -l ${params.FVE_viruslist} -1 reads_1.fq.gz -2 reads_2.fq.gz  -o ./
  mv FastViromeExplorer-reads-mapped-sorted.sam ${sample_id}_fastviromeexplorer.sam
  mv FastViromeExplorer-final-sorted-abundance.tsv ${sample_id}_fastviromeexplorer_abundance.tsv
  """
}

/**
TAX ASSIGNMENT - CONTIGS
**/

process tax_contigs_kraken2{
  tag {"${sample_id}_${assembler}"}

  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}", mode:'link'

  input:
  set sample_id, assembler, 'contigs.fa' from tax_contigs_kraken2_in

  output:
  set sample_id, assembler,"${sample_id}_${assembler}_kraken2.txt","${sample_id}_${assembler}_kraken2_report.txt" into tax_contigs_kraken2_out
  set sample_id, assembler,"${sample_id}_${assembler}_kraken2_unmapped.fa" into tax_contigs_unmapped_kraken2 /* Test unmapped reads */
  file "${sample_id}_${assembler}_kraken2_unmapped.fa" into tax_contigs_kraken2_unmapped

  script:
  """
  kraken2 --db ${params.kraken2_db} --threads ${task.cpus} --use-names \
    --unclassified-out ${sample_id}_${assembler}_kraken2_unmapped.fa \
    --output ${sample_id}_${assembler}_kraken2.txt \
    --report ${sample_id}_${assembler}_kraken2_report.txt contigs.fa
  """
}

process tax_contigs_diamond{
  tag {"${sample_id}_${assembler}"}

  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}", mode:'link'

  input:
  set sample_id, assembler, 'contigs.fa' from tax_contigs_diamond_in

  output:
  set sample_id, assembler,"${sample_id}_${assembler}_diamond.tsv" into tax_contigs_diamond_out
  set sample_id, assembler,"${sample_id}_${assembler}_diamond_unmapped.fa" into tax_contigs_unmapped_diamond /* Test unmapped reads */
  file "${sample_id}_${assembler}_diamond_unmapped.fa" into diamond_unmapped_out 						
  file 'diamond_blast_cols.txt' into tax_diamond_blast_cols
  
  

  script:
  blast_cols='qseqid sseqid qstart qend sstart send evalue score length pident nident mismatch positive ppos gapopen gaps staxids'
  """
  diamond blastx --sensitive --masking 1 -p ${task.cpus} --db ${params.diamond_db} \
    --taxonmap ${params.diamond_taxonmap} --taxonnodes ${params.diamond_taxonnodes} \
    --query contigs.fa --out ${sample_id}_${assembler}_diamond.tsv \
    --outfmt 6 ${blast_cols} \
    --un ${sample_id}_${assembler}_diamond_unmapped.fa
  echo '${blast_cols}' > diamond_blast_cols.txt 
  """
}

process tax_contigs_diamond_view{
  tag {"${sample_id}_${assembler}"}

  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}", mode:'link'

  input:
  set sample_id, assembler, 'diamond.daa' from Channel.empty()

  output:
  set sample_id, assembler,"${sample_id}_${assembler}_diamond.daa", "${sample_id}_${assembler}_diamond.tsv" into tax_contigs_diamond_view_out

  script:
  """
  diamond view -p ${task.cpus} --db ${params.diamond_db} \
    --outfmt 6 qseqid sseqid qstart qend sstart send evalue score length pident nident mismatch positive ppos gapopen gaps staxids \
    --daa diamond.daa \
    --out ${sample_id}_${assembler}_diamond.tsv
  """
}

/*
NOTE: https://github.com/EnvGen/toolbox/tree/master/scripts/assign_taxonomy_from_blast
*/

//process tax_diamond_lca{}

process tax_contigs_virfinder{
  tag {"${sample_id}_${assembler}"}

  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}", mode:'link'

  input:
  set sample_id, assembler, 'contigs.fa' from tax_contigs_virfinder_in

  output:
  set sample_id, assembler,"${sample_id}_${assembler}_virfinder.csv" into tax_contigs_virfinder_out

  script:
  """
  #!/usr/bin/env Rscript
  library(VirFinder)
  predResult <- VF.pred("contigs.fa")
  write.csv(predResult,file="${sample_id}_${assembler}_virfinder.csv")
  """
}

process tax_contigs_FragGeneScan{
  tag {"${sample_id}_${assembler}"}

  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}/orfs", mode:'link'

  input:
  set sample_id, assembler, 'contigs.fa' from tax_contigs_fgs_in

  output:
  set sample_id, assembler,"${sample_id}_${assembler}_fgs_orfs.*" into tax_contigs_fgs_out

  script:
  """
  run_FragGeneScan.pl -thread ${task.cpus} -complete 1 -train=illumina_10 \
    -genome=contigs.fa -out=${sample_id}_${assembler}_fgs_orfs
  """
}

unmapped_contigs = tax_contigs_unmapped_kraken2.combine(tax_contigs_unmapped_diamond, by: [0,1])

//unmapped_contigs.into{testPrint;testInput;unmapped_contigs1;}
/* testPrint.println() */

/* This process creates data for the extended discovery pipeline (discoveryExtended) */
process tax_contigs_unmapped_merged{
  tag {"${sample_id}_${assembler}"}

  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}", mode:'link'

  input:
  set sample_id, assembler, un_kraken2, un_diamond from unmapped_contigs
 
  output:
  set sample_id, assembler, "${sample_id}_${assembler}_unmapped_contigs.fa" into tax_contigs_merged_unmapped
  file "${sample_id}_${assembler}_unmapped_contigs.fa" into tax_contigs_merged_unmapped_out

  script:
   """
  #!/bin/bash
  cat ${un_kraken2} ${un_diamond} > tmp.fa 
  #This previous awk script was combining fasta reads from kraken2 and diamond and filtering duplicated reads
  #awk '/^>/{f=!d[\$1];d[\$1]=1}f' tmp.fa > ${sample_id}_${assembler}_unmapped_contigs.fa
  
  # New version of the script makes multiline fasta code into one line fasta code and keeps the fasta reads that exist in both kraken2 and diamond
  awk '/^>/ {printf("\\n%s\\n",\$0);next; } { printf("%s",\$0);}  END {printf("\\n");}' < tmp.fa > tmp1.fa
  rm tmp.fa
  awk 'seen[\$0]++ &&seen[\$0] > N' tmp1.fa > ${sample_id}_${assembler}_unmapped_contigs.fa
  rm tmp1.fa
  """
}

process tax_contigs_classifier{
 tag {"${sample_id}"}

 publishDir "${params.publish_base_dir}/${sample_id}/${assembler}/taxonomy", mode:'link'

 input:
 set sample_id, assembler, diamond from tax_contigs_diamond_out

 output:
 set sample_id, assembler, "${sample_id}_contigs_classified_all_sequences.tsv", "${sample_id}_contigs_classified_viruses.tsv", "${sample_id}_contigs_classified_bacteria.tsv","${sample_id}_contigs_classified_eukaryota.tsv", "${sample_id}_contigs_classified_other_sequences.tsv" into tax_contigs_taxonomy_out

 script:
 """
 #!/bin/bash
 touch "${sample_id}_contigs_classified.tsv"
 echo -e "sample_id" '\\t' "contig_id" '\\t' "accession_nr" '\\t' "scientific_name" '\\t' "title" '\\t' "percentage" '\\t' "evalue" '\\t' "mismatch" '\\t' "aln_len" '\\t' "contig_len" '\\t' "full_lineage" >> "${sample_id}_contigs_classified_all_sequences.tsv"
 echo -e "sample_id" '\\t' "contig_id" '\\t' "accession_nr" '\\t' "scientific_name" '\\t' "title" '\\t' "percentage" '\\t' "evalue" '\\t' "mismatch" '\\t' "aln_len" '\\t' "contig_len" '\\t' "full_lineage" >> "${sample_id}_contigs_classified_viruses.tsv"
 echo -e "sample_id" '\\t' "contig_id" '\\t' "accession_nr" '\\t' "scientific_name" '\\t' "title" '\\t' "percentage" '\\t' "evalue" '\\t' "mismatch" '\\t' "aln_len" '\\t' "contig_len" '\\t' "full_lineage" >> "${sample_id}_contigs_classified_bacteria.tsv"
 echo -e "sample_id" '\\t' "contig_id" '\\t' "accession_nr" '\\t' "scientific_name" '\\t' "title" '\\t' "percentage" '\\t' "evalue" '\\t' "mismatch" '\\t' "aln_len" '\\t' "contig_len" '\\t' "full_lineage" >> "${sample_id}_contigs_classified_eukaryota.tsv"
 echo -e "sample_id" '\\t' "contig_id" '\\t' "accession_nr" '\\t' "scientific_name" '\\t' "title" '\\t' "percentage" '\\t' "evalue" '\\t' "mismatch" '\\t' "aln_len" '\\t' "contig_len" '\\t' "full_lineage" >> "${sample_id}_contigs_classified_other_sequences.tsv"
 while read p; do
     full_lineage="-"
     taxid=\$(echo "\$p" | awk '{print \$17}')
     SUB=";"
     if [[ "\$taxid" == *"\$SUB"* ]]; then
     	IFS=';'
     	read -a taxid_arr <<< "\$taxid"
     	taxid=\${lineage_arr[0]}
     elif [[ -z \$taxid ]]; then
     	taxid="empty"
     fi
     accid=\$(echo "\$p" | awk '{print \$2}')
     contig_line=\$(echo "\$p" | awk '{print \$1}')
     IFS=','
     read -a contig_arr <<< "\$contig_line"     
     contig_id="\$(echo \${contig_arr[0]} | xargs)"
     contig_len="\$(echo \${contig_arr[3]} | xargs)"
     wordToRemove="len="
     contig_len=\${contig_len//\$wordToRemove/}
     full_lineage=\$(grep -w \$taxid ${params.ncbi_full_lineage})
     if [ -n "\$full_lineage" ]; then
     IFS='|'
     read -a lineage_arr <<< "\$full_lineage"
     	scientific_name="\$(echo \${lineage_arr[1]} | xargs | sed -e 's/ /_/g')"
     	full_lineage2="\$(echo \${lineage_arr[2]} | xargs | sed -e 's/ /_/g')"
     else
     	scientific_name="-"
     	full_lineage2="-"
     fi
     title="\$(blastdbcmd -db nr -entry "\$accid" -range 1 | head -n 1 | sed 's/>//g' | sed 's/:1-1//g' | sed "s/\$accid//g" | xargs | sed -e 's/ /_/g')"
     perc=\$(echo "\$p" | awk '{print \$10}')
     aln_len=\$(echo "\$p" | awk '{print \$9}')
     evalue=\$(echo "\$p" | awk '{print \$7}')
     mismatch=\$(echo "\$p" | awk '{print \$12}')
     if [[ "\$full_lineage2" == *"Virus"* ]]; then
        echo -e ${sample_id} '\\t' \$contig_id '\\t' \$accid '\\t' \$scientific_name '\\t' \$title '\\t' \$perc '\\t' \$evalue '\\t' \$mismatch '\\t' \$aln_len '\\t' \$contig_len '\\t' \$full_lineage2 >> "${sample_id}_contigs_classified_viruses.tsv"
     elif [[ "\$full_lineage2" == *"Bacteria"* ]]; then
        echo -e ${sample_id} '\\t' \$contig_id '\\t' \$accid '\\t' \$scientific_name '\\t' \$title '\\t' \$perc '\\t' \$evalue '\\t' \$mismatch '\\t' \$aln_len '\\t' \$contig_len '\\t' \$full_lineage2 >> "${sample_id}_contigs_classified_bacteria.tsv"
     elif [[ "\$full_lineage2" == *"Eukaryota"* ]]; then
        echo -e ${sample_id} '\\t' \$contig_id '\\t' \$accid '\\t' \$scientific_name '\\t' \$title '\\t' \$perc '\\t' \$evalue '\\t' \$mismatch '\\t' \$aln_len '\\t' \$contig_len '\\t' \$full_lineage2 >> "${sample_id}_contigs_classified_eukaryota.tsv"
     else
        echo -e ${sample_id} '\\t' \$contig_id '\\t' \$accid '\\t' \$scientific_name '\\t' \$title '\\t' \$perc '\\t' \$evalue '\\t' \$mismatch '\\t' \$aln_len '\\t' \$contig_len '\\t' \$full_lineage2 >> "${sample_id}_contigs_classified_other_sequences.tsv"
     fi
     echo -e ${sample_id} '\\t' \$contig_id '\\t' \$accid '\\t' \$scientific_name '\\t' \$title '\\t' \$perc '\\t' \$evalue '\\t' \$mismatch '\\t' \$aln_len '\\t' \$contig_len '\\t' \$full_lineage2 >> "${sample_id}_contigs_classified_all_sequences.tsv"
 done < ${diamond}

 """
}

tax_contigs_taxonomy_out.into{data_for_html_table;
collect_all_tsv;
collect_virus_tsv;
collect_bact_tsv;
collect_eukar_tsv;
collect_other_tsv}

html_data_each = data_for_html_table.combine(html_each)

process TSV_to_HTML_each_sample{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}/taxonomy", mode:'link'

  input:
  set sample_id, assembler, all, viral, bacteria, eukaryota, other, hvalue, html from html_data_each
 
  output:
  set sample_id, "${sample_id}_all_asm_results.html", "${sample_id}_viral_asm_results.html", "${sample_id}_bacteria_asm_results.html", "${sample_id}_eukaryota_asm_results.html", "${sample_id}_other_asm_results.html" into html_each_out
  script:
"""
#!/home/amanj/anaconda3/envs/amanjEnv/bin/python3
import string
import csv
import json
import collections

OrderedDict = collections.OrderedDict

def HTML_table(src, header, htmlEndStr, name):
    
    f = open(src)
    lines = f.readlines()
    f.close()
    f1 = open(header)
    header_lines = f1.readlines()
    f2 = open(htmlEndStr)
    end_lines = f2.readlines()
    f1.close()
    f2.close()
    
    with open(name, 'w') as f:
        for h_line in header_lines:
            f.write(h_line)
        for line in lines:
            f.write(line)
        for e_line in end_lines:
            f.write(e_line)
    f.close()
    return None

def TSV_file_into_JSON(src, dst, header):
    data = []
    with open(src, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter='\\t', quotechar='"')
        for row in reader:
            if row[0] == 'sample_id':
                print("\\n")
            else:
                if row[0].strip()[0] == '#':  #
                    continue
                row = filter(None, row)
                data.append(OrderedDict(zip(header, row)))

    with open(dst, 'w') as jsonfile:
        json.dump(data, jsonfile, indent=2)
    return None

header = ['sample_id', 'contig_id', 'accession_nr', 'scientific_name', 'title', 'percentage', 'evalue', 'mismatch', 'aln_len', 'contig_len', 'full_lineage']
src = "${all}"
dst = "all.json"  
TSV_file_into_JSON(src, dst, header)
HTML_table(dst, "${html[1]}", "${html[0]}", "${sample_id}_all_asm_results.html")

src = "${viral}"
dst = "viral.json"  
TSV_file_into_JSON(src, dst, header)
HTML_table(dst, "${html[1]}", "${html[0]}", "${sample_id}_viral_asm_results.html")

src = "${bacteria}"
dst = "bact.json"
TSV_file_into_JSON(src, dst, header)
HTML_table(dst, "${html[1]}", "${html[0]}", "${sample_id}_bacteria_asm_results.html")

src = "${eukaryota}"
dst = "eukar.json"
TSV_file_into_JSON(src, dst, header)
HTML_table(dst, "${html[1]}", "${html[0]}", "${sample_id}_eukaryota_asm_results.html")

src = "${other}"
dst = "other.json"
TSV_file_into_JSON(src, dst, header)
HTML_table(dst, "${html[1]}", "${html[0]}", "${sample_id}_other_asm_results.html")
"""
}

process collect_all_tsv_tables{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/tsv_tables", mode:'link'

  input:
  file "tsv_table" from collect_all_tsv.map{it[2]}.collect()
 
  output:
  file "${params.project_id}_all_asm_results.tsv" into all_tsv_collected
  script:
"""
#!/bin/bash
touch tmp.tsv
echo -e "sample_id\\tcontig_id\\taccession_nr\\tscientific_name\\ttitle\\tpercentage\\tevalue\\tmismatch\\taln_len\\tcontig_len\\tfull_lineage" > header.tsv

for sample_file in ${tsv_table}
do
	cat \$sample_file | sed '1d' >> tmp.tsv
done

cat tmp.tsv | sort -k1,1 -k2,2 > tmp2.tsv
cat header.tsv tmp2.tsv > "${params.project_id}_all_asm_results.tsv"
rm tmp.tsv tmp2.tsv header.tsv
"""
}

process collect_virus_tsv_tables{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/tsv_tables", mode:'link'

  input:
  file "tsv_table" from collect_virus_tsv.map{it[3]}.collect()

  output:
  file "${params.project_id}_virus_asm_results.tsv" into virus_tsv_collected
  script:
"""
#!/bin/bash
touch tmp.tsv
echo -e "sample_id\\tcontig_id\\taccession_nr\\tscientific_name\\ttitle\\tpercentage\\tevalue\\tmismatch\\taln_len\\tcontig_len\\tfull_lineage" > header.tsv

for sample_file in ${tsv_table}
do
        cat \$sample_file | sed '1d' >> tmp.tsv
done

cat tmp.tsv | sort -k1,1 -k2,2 > tmp2.tsv
cat header.tsv tmp2.tsv > "${params.project_id}_virus_asm_results.tsv"
rm tmp.tsv tmp2.tsv header.tsv
"""
}

process collect_bacteria_tsv_tables{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/tsv_tables", mode:'link'

  input:
  file "tsv_table" from collect_bact_tsv.map{it[4]}.collect()

  output:
  file "${params.project_id}_bacteria_asm_results.tsv" into bact_tsv_collected
  script:
"""
#!/bin/bash
touch tmp.tsv
echo -e "sample_id\\tcontig_id\\taccession_nr\\tscientific_name\\ttitle\\tpercentage\\tevalue\\tmismatch\\taln_len\\tcontig_len\\tfull_lineage" > header.tsv

for sample_file in ${tsv_table}
do
        cat \$sample_file | sed '1d' >> tmp.tsv
done

cat tmp.tsv | sort -k1,1 -k2,2 > tmp2.tsv
cat header.tsv tmp2.tsv > "${params.project_id}_bacteria_asm_results.tsv"
rm tmp.tsv tmp2.tsv header.tsv
"""
}

process collect_eukaryotic_tsv_tables{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/tsv_tables", mode:'link'

  input:
  file "tsv_table" from collect_eukar_tsv.map{it[5]}.collect()

  output:
  file "${params.project_id}_eukaryotic_asm_results.tsv" into eukar_tsv_collected
  script:
"""
#!/bin/bash
touch tmp.tsv
echo -e "sample_id\\tcontig_id\\taccession_nr\\tscientific_name\\ttitle\\tpercentage\\tevalue\\tmismatch\\taln_len\\tcontig_len\\tfull_lineage" > header.tsv

for sample_file in ${tsv_table}
do
        cat \$sample_file | sed '1d' >> tmp.tsv
done

cat tmp.tsv | sort -k1,1 -k2,2 > tmp2.tsv
cat header.tsv tmp2.tsv > "${params.project_id}_eukaryotic_asm_results.tsv"
rm tmp.tsv tmp2.tsv header.tsv
"""
}

process collect_otherSeq_tsv_tables{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/tsv_tables", mode:'link'

  input:
  file "tsv_table" from collect_other_tsv.map{it[6]}.collect()

  output:
  file "${params.project_id}_other_sequences_asm_results.tsv" into other_tsv_collected
  script:
"""
#!/bin/bash
touch tmp.tsv
echo -e "sample_id\\tcontig_id\\taccession_nr\\tscientific_name\\ttitle\\tpercentage\\tevalue\\tmismatch\\taln_len\\tcontig_len\\tfull_lineage" > header.tsv

for sample_file in ${tsv_table}
do
        cat \$sample_file | sed '1d' >> tmp.tsv
done

cat tmp.tsv | sort -k1,1 -k2,2 > tmp2.tsv
cat header.tsv tmp2.tsv > "${params.project_id}_other_sequences_asm_results.tsv"
rm tmp.tsv tmp2.tsv header.tsv
"""
}

collected_tsv = all_tsv_collected.combine(virus_tsv_collected).combine(bact_tsv_collected).combine(eukar_tsv_collected).combine(other_tsv_collected)

collected_tsv.into{collected_tsv_into_html;
collected_tsv_to_simplified}

process tsv_to_simplified_tsv{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/simplified_tsv_tables", mode:'link'

  input:
  set all, virus, bacteria, eukaryotic, other from collected_tsv_to_simplified

  output:
  set "${params.project_id}_all_simplified_asm_results.tsv", "${params.project_id}_virus_simplified_asm_results.tsv", "${params.project_id}_bacteria_simplified_asm_results.tsv", "${params.project_id}_eukaryota_simplified_asm_results.tsv", "${params.project_id}_other_simplified_asm_results.tsv" into tsv_simplifed_out
  script:
"""
#!/bin/bash
echo -e "sample_id\\tcontig_id\\taccession_nr\\tscientific_name\\ttitle\\tpercentage\\tevalue\\tmismatch\\taln_len\\tcontig_len\\tfull_lineage" > header.tsv
cat ${all} | sed 1d | sed 's/ //g' | datamash -s -g1,2 min 7 -f | awk '{print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" \$6 "\\t" \$7 "\\t" \$8 "\\t" \$9 "\\t" \$10 "\\t" \$11}' | sed 's/\\,/./g' | sort -k1,1 -k3,3 > all_tmp.tsv

cat ${virus} | sed 1d | sed 's/ //g' | datamash -s -g1,2 min 7 -f | awk '{print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" \$6 "\\t" \$7 "\\t" \$8 "\\t" \$9 "\\t" \$10 "\\t" \$11}' | sed 's/\\,/./g' | sort -k1,1 -k3,3 > virus_tmp.tsv

cat ${bacteria} | sed 1d | sed 's/ //g' | datamash -s -g1,2 min 7 -f | awk '{print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" \$6 "\\t" \$7 "\\t" \$8 "\\t" \$9 "\\t" \$10 "\\t" \$11}' | sed 's/\\,/./g' | sort -k1,1 -k3,3 > bact_tmp.tsv

cat ${eukaryotic} | sed 1d | sed 's/ //g' | datamash -s -g1,2 min 7 -f | awk '{print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" \$6 "\\t" \$7 "\\t" \$8 "\\t" \$9 "\\t" \$10 "\\t" \$11}' | sed 's/\\,/./g' | sort -k1,1 -k3,3 > eukar_tmp.tsv

cat ${other} | sed 1d | sed 's/ //g' | datamash -s -g1,2 min 7 -f | awk '{print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" \$6 "\\t" \$7 "\\t" \$8 "\\t" \$9 "\\t" \$10 "\\t" \$11}' | sed 's/\\,/./g' | sort -k1,1 -k3,3 > other_tmp.tsv

cat header.tsv all_tmp.tsv > "${params.project_id}_all_simplified_asm_results.tsv"
cat header.tsv virus_tmp.tsv > "${params.project_id}_virus_simplified_asm_results.tsv"
cat header.tsv bact_tmp.tsv > "${params.project_id}_bacteria_simplified_asm_results.tsv"
cat header.tsv eukar_tmp.tsv > "${params.project_id}_eukaryota_simplified_asm_results.tsv"
cat header.tsv other_tmp.tsv > "${params.project_id}_other_simplified_asm_results.tsv"
rm header.tsv all_tmp.tsv virus_tmp.tsv bact_tmp.tsv eukar_tmp.tsv other_tmp.tsv
"""
}


html_simplifed_in = tsv_simplifed_out.combine(html_simplified_all)
html_collection_tsv_in = collected_tsv_into_html.combine(html_all)

process simplified_tsv_to_html{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/simplified_html_tables", mode:'link'

  input:
  set all, virus, bacteria, eukaryotic, other, hvalue, html from html_simplifed_in

  output:
  set "${params.project_id}_all_simplified_asm_results.html", "${params.project_id}_virus_simplified_asm_results.html", "${params.project_id}_bacteria_simplified_asm_results.html", "${params.project_id}_eukaryota_simplified_asm_results.html", "${params.project_id}_other_simplified_asm_results.html" into html_simplifed_out
  script:
"""
#!/home/amanj/anaconda3/envs/amanjEnv/bin/python3
import string
import csv
import json
import collections

OrderedDict = collections.OrderedDict

def HTML_table(src, header, htmlEndStr, name):

    f = open(src)
    lines = f.readlines()
    f.close()
    f1 = open(header)
    header_lines = f1.readlines()
    f2 = open(htmlEndStr)
    end_lines = f2.readlines()
    f1.close()
    f2.close()

    with open(name, 'w') as f:
        for h_line in header_lines:
            f.write(h_line)
        for line in lines:
            f.write(line)
        for e_line in end_lines:
            f.write(e_line)
    f.close()
    return None

def TSV_file_into_JSON(src, dst, header):
    data = []
    with open(src, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter='\\t', quotechar='"')
        for row in reader:
            if row[0] == 'sample_id':
                print("\\n")
            else:
                if row[0].strip()[0] == '#':  #
                    continue
                row = filter(None, row)
                data.append(OrderedDict(zip(header, row)))

    with open(dst, 'w') as jsonfile:
        json.dump(data, jsonfile, indent=2)
    return None

header = ['sample_id', 'contig_id', 'accession_nr', 'scientific_name', 'title', 'percentage', 'evalue', 'mismatch', 'aln_len', 'contig_len', 'full_lineage']
src = "${all}"
dst = "all.json"
TSV_file_into_JSON(src, dst, header)
HTML_table(dst, "${html[1]}", "${html[0]}", "${params.project_id}_all_simplified_asm_results.html")

src = "${virus}"
dst = "viral.json"
TSV_file_into_JSON(src, dst, header)
HTML_table(dst, "${html[1]}", "${html[0]}", "${params.project_id}_virus_simplified_asm_results.html")

src = "${bacteria}"
dst = "bact.json"
TSV_file_into_JSON(src, dst, header)
HTML_table(dst, "${html[1]}", "${html[0]}", "${params.project_id}_bacteria_simplified_asm_results.html")

src = "${eukaryotic}"
dst = "eukar.json"
TSV_file_into_JSON(src, dst, header)
HTML_table(dst, "${html[1]}", "${html[0]}", "${params.project_id}_eukaryota_simplified_asm_results.html")

src = "${other}"
dst = "other.json"
TSV_file_into_JSON(src, dst, header)
HTML_table(dst, "${html[1]}", "${html[0]}", "${params.project_id}_other_simplified_asm_results.html")

"""
}

process collected_tsv_to_html{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/html_tables", mode:'link'

  input:
  set all, virus, bacteria, eukaryotic, other, hvalue, html from html_collection_tsv_in

  output:
  set "${params.project_id}_all_asm_results.html", "${params.project_id}_virus_asm_results.html", "${params.project_id}_bacteria_asm_results.html", "${params.project_id}_eukaryota_asm_results.html", "${params.project_id}_other_asm_results.html" into html_collection_tsv_out
  script:
"""
#!/home/amanj/anaconda3/envs/amanjEnv/bin/python3
import string
import csv
import json
import collections

OrderedDict = collections.OrderedDict

def HTML_table(src, header, htmlEndStr, name):

    f = open(src)
    lines = f.readlines()
    f.close()
    f1 = open(header)
    header_lines = f1.readlines()
    f2 = open(htmlEndStr)
    end_lines = f2.readlines()
    f1.close()
    f2.close()

    with open(name, 'w') as f:
        for h_line in header_lines:
            f.write(h_line)
        for line in lines:
            f.write(line)
        for e_line in end_lines:
            f.write(e_line)
    f.close()
    return None

def TSV_file_into_JSON(src, dst, header):
    data = []
    with open(src, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter='\\t', quotechar='"')
        for row in reader:
            if row[0] == 'sample_id':
                print("\\n")
            else:
                if row[0].strip()[0] == '#':  #
                    continue
                row = filter(None, row)
                data.append(OrderedDict(zip(header, row)))

    with open(dst, 'w') as jsonfile:
        json.dump(data, jsonfile, indent=2)
    return None

header = ['sample_id', 'contig_id', 'accession_nr', 'scientific_name', 'title', 'percentage', 'evalue', 'mismatch', 'aln_len', 'contig_len', 'full_lineage']
src = "${all}"
dst = "all.json"
TSV_file_into_JSON(src, dst, header)
HTML_table(dst, "${html[1]}", "${html[0]}", "${params.project_id}_all_asm_results.html")

src = "${virus}"
dst = "viral.json"
TSV_file_into_JSON(src, dst, header)
HTML_table(dst, "${html[1]}", "${html[0]}", "${params.project_id}_virus_asm_results.html")

src = "${bacteria}"
dst = "bact.json"
TSV_file_into_JSON(src, dst, header)
HTML_table(dst, "${html[1]}", "${html[0]}", "${params.project_id}_bacteria_asm_results.html")

src = "${eukaryotic}"
dst = "eukar.json"
TSV_file_into_JSON(src, dst, header)
HTML_table(dst, "${html[1]}", "${html[0]}", "${params.project_id}_eukaryota_asm_results.html")

src = "${other}"
dst = "other.json"
TSV_file_into_JSON(src, dst, header)
HTML_table(dst, "${html[1]}", "${html[0]}", "${params.project_id}_other_asm_results.html")
"""
}





