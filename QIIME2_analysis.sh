#This script relates to pre-processing of demultiplexed paired-end 16S and ITS2 rDNA sequences (available at: https://www.ebi.ac.uk/ena/browser/view/PRJEB49004) 
#for Cullen et al. (2022) in Animal - Open Space. Demultiplexed paired-end sequences were imported (in Casava 1.8 demultiplexed paired-end format) 
#into QIIME2 v.2020.8.0 (Bolyen et al., 2019), which was installed on a virtual machine (VirtualBox 6.0). 
#Forward and reverse reads were quality assessed using the ‘qiime demux summarize’ command, FastQC v.0.11.5 and MultiQC v.1.9. 
#The 16S and ITS2 primers were removed from reads using the cutadapt plugin (Martin, 2011). 
#The QIIME2 DADA2 (Callahan et al., 2016) plugin was used for filtering and dereplication, chimera removal, merging paired-end reads 
#and to infer amplicon sequence variants (ASVs) in each sample after truncating reads to remove low quality bases. 
#Read 1 and read 2 of the 16S rDNA sequences were truncated at 267 and 183 bp, respectively, while ITS2 sequences were truncated at 266 and 187 bp, respectively. 
#For bacterial sequences, taxonomy was assigned to each ASV using a Naive Bayes classifier trained on 16S V3-V4 sequences from the Silva database (Version 138) 
#with the ‘q2-feature-classifier’ plugin, while taxonomy was assigned to fungal ASVs using a Naive Bayes classifier trained on full-length ITS sequences 
#from the UNITE v.8.3 database (Kõljalg et al., 2013).

###16S sequence pre-processing###

#Make a new directory to store analysis outputs
mkdir meth-dev-16S

#Move into the bead_beating_16S directory
cd meth-dev-16S

#Make a new directory for FastQC reports
mkdir fastqc_reports

#Run FastQC on raw data to assess quality
fastqc /home/qiime2/meth-dev-16S-files/* -O fastqc_reports

#Use MulitQC to compile a QC report of all samples
multiqc -n meth-dev-16S_report fastqc_reports

#Import demultiplexed data in Casava 1.8 demultiplexed (paired-end) format to generate .qza file
qiime tools import   --type 'SampleData[PairedEndSequencesWithQuality]' \  
--input-path /home/qiime2/meth-dev-16S-files \   
--input-format CasavaOneEightSingleLanePerSampleDirFmt \  
--output-path meth-dev-16S.qza 

#Create a .qzv file to summarise and visualise the demultiplexed data (view the .qzv file by dragging and dropping the file at https://view.qiime2.org)
#The number of reads per sample and an interactive quality plot can be viewed to assess quality
qiime demux summarize \
--i-data meth-dev-16S.qza \
--o-visualization meth-dev-16S-summary.qzv

#Use the cutadapt plugin to trim primer sequences from the reads and generate a new trimmed .qza file
qiime cutadapt trim-paired \
--i-demultiplexed-sequences meth-dev-16S.qza \
--p-front-f CCTACGGGNGGCWGCAG \
--p-front-r GACTACHVGGGTATCTAATCC \
--o-trimmed-sequences meth-dev-16S-trimmed.qza

#Create a new .qzv file to summarise and visualise the trimmed data (https://view.qiime2.org)
#Use the interactive quality plot to assess truncation parameters prior to DADA2 processing
qiime demux summarize  --i-data meth-dev-16S-trimmed.qza  --o-visualization meth-dev-16S-trimmed-demuz.qzv

#Use DADA2 plugin to denoise and dereplicate sequences, to remove chimeric sequences, merge reads and generate the amplicon sequence variant (ASV) 
#table, representaive sequences, and denoising stats 
#Set truncation parameters for read 1 and read 2 based on quality scores to remove low quality bases 
# The default minimum overlap required by DADA2 to merge read 1 and read 2 is 12 bases)
qiime dada2 denoise-paired \
--i-demultiplexed-seqs meth-dev-16S-trimmed.qza \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 267 \
--p-trunc-len-r 183 \
--o-table meth-dev-16S-trimmed-table.qza \
--p-n-threads 0 \
--o-representative-sequences meth-dev-16S-trimmed-rep-seqs.qza \
--o-denoising-stats meth-dev-16S-trimmed-denoising-stats.qza





