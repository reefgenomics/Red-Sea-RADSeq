#!/bin/bash

###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 00 Raw paired-end reads processing #####

####  00.01 Adapter removal ####
# Only reads at least 110 bp were kept

# Flow Cell 1 170213_D00658_0007_AC9N98ANXX
cd ~/RADseq/Raw_files/170213_D00658_0007_AC9N98ANXX
for i in {1..8}; do trimmomatic-0.36.jar PE -threads 30 -phred33 -trimlog trimlog-lib1.txt ./Lane${i}/M_17_*_${i}_NoIndex_L00${i}_R1_001.fastq.gz ./Lane${i}/M_17_*_${i}_NoIndex_L00${i}_R2_001.fastq.gz -baseout lib${i}-non-adapt.fq.gz ILLUMINACLIP:/Truseq2_adapters_trimmomatic/TruSeq2-PE.fa:2:30:10 MINLEN:110 &>> lib${i}-adapt-rem-trimmo.txt

# Flow Cell 2 170307_D00658_0013_AC9TH4ANXX
cd ~/RADseq/Raw_files/170307_D00658_0013_AC9TH4ANXX
for i in {1..8}; do trimmomatic-0.36.jar PE -threads 30 -phred33 -trimlog trimlog-lib1.txt ./Lane${i}/M_17_*_pool$((i+8))_NoIndex_L00${i}_R1_001.fastq.gz ./Lane${i}/M_17_*_pool$((i+8))_NoIndex_L00${i}_R2_001.fastq.gz -baseout lib$((i+8))-non-adapt.fq.gz ILLUMINACLIP:/Truseq2_adapters_trimmomatic/TruSeq2-PE.fa:2:30:10 MINLEN:110 &>> lib$((i+8))-adapt-rem-trimmo.txt

# Flow Cell 3 170307_D00658_0014_BC9T6WANXX
cd ~/RADseq/Raw_files/170307_D00658_0014_BC9T6WANXX
for i in {1..8}; do trimmomatic-0.36.jar PE -threads 30 -phred33 -trimlog trimlog-lib1.txt ./Lane${i}/M_17_*_pool$((i+16))_NoIndex_L00${i}_R1_001.fastq.gz ./Lane${i}/M_17_*_pool$((i+16))_NoIndex_L00${i}_R2_001.fastq.gz -baseout lib$((i+16))-non-adapt.fq.gz ILLUMINACLIP:/Truseq2_adapters_trimmomatic/TruSeq2-PE.fa:2:30:10 MINLEN:110 &>> lib$((i+16))-adapt-rem-trimmo.txt

# Flow Cell 4 170313_7001439_0202_BC9U0AANXX
cd ~/RADseq/Raw_files/170313_7001439_0202_BC9U0AANXX
for i in {1..8}; do trimmomatic-0.36.jar PE -threads 30 -phred33 -trimlog trimlog-lib1.txt ./Lane${i}/M_17_*_$((i+24))_NoIndex_L00${i}_R1_001.fastq.gz ./Lane${i}/M_17_*_$((i+24))_NoIndex_L00${i}_R2_001.fastq.gz -baseout lib$((i+24))-non-adapt.fq.gz ILLUMINACLIP:/Truseq2_adapters_trimmomatic/TruSeq2-PE.fa:2:30:10 MINLEN:110 &>> lib$((i+24))-adapt-rem-trimmo.txt


####  00.02 Libraries demultiplexing ####

# Flow Cell 1 170213_D00658_0007_AC9N98ANXX
cd ~/RADseq/Raw_files/170213_D00658_0007_AC9N98ANXX
for i in {1..8}; do stacks/process_radtags -i gzfastq -P -1 ./Lane${i}/lib${i}-non-adapt_1P.fq.gz -2 ./Lane${i}/lib${i}-non-adapt_2P.fq.gz -o ../../Samples -b ../process_RADtags_txt-files/barcodes_lib${i} --barcode_dist_1 1 -E phred33 --retain_header -D -e pstI -r -c -q &>> lib${i}-process-RADtag.log; done

# Flow Cell 2 170307_D00658_0013_AC9TH4ANXX
cd ~/RADseq/Raw_files/170307_D00658_0013_AC9TH4ANXX
for i in {1..8}; do stacks/process_radtags -i gzfastq -P -1 ./Lane${i}/lib$((i+8))-non-adapt_1P.fq.gz -2 ./Lane${i}/lib$((i+8))-non-adapt_2P.fq.gz -o ../../Samples -b ../process_RADtags_txt-files/barcodes_lib$((i+8)) --barcode_dist_1 1 -E phred33 --retain_header -D -e pstI -r -c -q &>> lib$((i+8))-process-RADtag.log; done

# Flow Cell 3 170307_D00658_0014_BC9T6WANXX
cd ~/RADseq/Raw_files/170307_D00658_0014_BC9T6WANXX
for i in {1..8}; do stacks/process_radtags -i gzfastq -P -1 ./Lane${i}/lib$((i+16))-non-adapt_1P.fq.gz -2 ./Lane${i}/lib$((i+16))-non-adapt_2P.fq.gz -o ../../Samples -b ../process_RADtags_txt-files/barcodes_lib$((i+16)) --barcode_dist_1 1 -E phred33 --retain_header -D -e pstI -r -c -q &>> lib$((i+16))-process-RADtag.log; done

# Flow Cell 4 170313_7001439_0202_BC9U0AANXX
cd ~/RADseq/Raw_files/170313_7001439_0202_BC9U0AANXX
for i in {1..8}; do stacks/process_radtags -i gzfastq -P -1 ./Lane${i}/lib$((i+24))-non-adapt_1P.fq.gz -2 ./Lane${i}/lib$((i+24))-non-adapt_2P.fq.gz -o ../../Samples -b ../process_RADtags_txt-files/barcodes_lib$((i+24)) --barcode_dist_1 1 -E phred33 --retain_header -D -e pstI -r -c -q &>> lib$((i+24))-process-RADtag.log; done


####  00.03 Quality check  ####
cd ~/RADseq/Samples
while IFS= read -r line || [[ -n "$line" ]]; do fastqc ./$line.1.fq.gz ./$line.2.fq.gz -o ./FastQC/ ; done < samples_id


#### 00.04 PCR duplicates removal ####
while IFS= read -r line || [[ -n "$line" ]]; do stacks/clone_filter -P -1 ./$line.1.fq.gz -2 ./$line.2.fq.gz -i gzfastq -o ./Clone-filtered-reads/ -D &>> libs_clone-filter_log.txt; done < samples_id


#### 00.05 Crop sequences to 110 bp ####
cd ~/RADseq/Samples/Clone-filtered-reads
while IFS= read -r line || [[ -n "$line" ]]; do trimmomatic-0.36.jar PE -phred33 -trimlog trimlog_$line.txt ./${i}.1.1.fq.gz ./${i}.2.2.fq.gz -baseout ${i}_qtrim.fq.gz SLIDINGWINDOW:10:20 MINLEN:110 CROP:110 &>>libs-trimmomatic-110_log.txt; done < samples_id


#### 00.06 Map sequences to reference genome ####
## Generate an Index for referene genome - bowtie2-build version 2.3.3##

# cd /home/RADseq/BOWTIE_INDEXES
# bowtie2-build Spis.genome.scaffold.final.fa Spis_Genome
# bowtie2-build Pver_genome_assembly_v1.0.fasta Pver_Genome

# Split each species samples
# grep "^S" samples_id > Spis_samples
# grep "^P" samples_id > Pver_samples 

while IFS= read -r line || [[ -n "$line" ]];do bowtie2 -x /home/RADseq/BOWTIE_INDEXES/Spis_Genome --local --very-sensitive-local --no-discordant --no-mixed -1 ./$line_qtrim_1P.fq.gz -2 ./$line_qtrim_2P.fq.gz -S ./$line.sam --al-conc-gz ../Alignments/al-conc/$line_%.conc.fq.gz --un-conc-gz ../Alignments/non-al/$line_%.un.fq.gz  -p 20 &>> ../Alignments/libs-Spis-bowtie2-log.txt; done < Spis_samples

while IFS= read -r line || [[ -n "$line" ]];do bowtie2 -x /home/RADseq/BOWTIE_INDEXES/Pver_Genome --local --very-sensitive-local --no-discordant --no-mixed -1 ./$line_qtrim_1P.fq.gz -2 ./$line_qtrim_2P.fq.gz -S ./$line.sam --al-conc-gz ../Alignments/al-conc/$line_%.conc.fq.gz --un-conc-gz ../Alignments/non-al/$line_%.un.fq.gz  -p 20 &>> ../Alignments/libs-Pver-bowtie2-log.txt; done < Pver_samples

# Convert SAM to BAM files
cd ~/RADseq/Samples/Alignments/al-conc
for i in *.sam; do samtools view -b -S -q 10 ${i} -@ 30 | samtools sort - -o ${i}.bam -@ 30; done
# rm -rf *.sam # remove SAM files as those occupy a lot of space

