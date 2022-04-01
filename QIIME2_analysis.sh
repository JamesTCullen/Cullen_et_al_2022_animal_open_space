#This script relates to pre-processing of demultiplexed paired-end 16S and ITS2 rDNA sequences (available at: https://www.ebi.ac.uk/ena/browser/view/PRJEB49004) 
#for Cullen et al. (2022) published in Animal - Open Space. Demultiplexed paired-end sequences were imported (in Casava 1.8 demultiplexed paired-end format) 
#into QIIME2 v.2020.8.0 (Bolyen et al., 2019), which was installed on a virtual machine (VirtualBox 6.0). 
#Forward and reverse reads were quality assessed using the ‘qiime demux summarize’ command, FastQC v.0.11.5 and MultiQC v.1.9. 
#The 16S and ITS2 primers were removed from reads using the cutadapt plugin (Martin, 2011). 
#The QIIME2 DADA2 (Callahan et al., 2016) plugin was used for filtering and dereplication, chimera removal, merging paired-end reads 
#and to infer amplicon sequence variants (ASVs) in each sample after truncating reads to remove low quality bases. 
#Read 1 and read 2 of the 16S rDNA sequences were truncated at 267 and 183 bp, respectively, while ITS2 sequences were truncated at 266 and 187 bp, respectively. 
#For bacterial sequences, taxonomy was assigned to each ASV using a Naive Bayes classifier trained on 16S V3-V4 sequences from the Silva database (Version 138) 
#with the ‘q2-feature-classifier’ plugin, while taxonomy was assigned to fungal ASVs using a Naive Bayes classifier trained on full-length ITS sequences 
#from the UNITE v.8.3 database (Kõljalg et al., 2013).

###16S pre-processing###

#Make a new directory to store 16S analysis outputs
mkdir meth-dev-16S

#Move into the meth-dev-16S directory
cd meth-dev-16S

#Make a new directory for 16S FastQC reports
mkdir fastqc-reports

#Run FastQC on raw 16S data to assess quality
fastqc /home/qiime2/meth-dev-16S-files/* -O fastqc-reports

#Use MulitQC to compile a QC report of all 16S samples
multiqc -n meth-dev-16S-report fastqc-reports

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

#Use the cutadapt plugin to trim V3-V4 region primer sequences from the reads and generate a new trimmed .qza file
qiime cutadapt trim-paired \
--i-demultiplexed-sequences meth-dev-16S.qza \
--p-front-f CCTACGGGNGGCWGCAG \
--p-front-r GACTACHVGGGTATCTAATCC \
--o-trimmed-sequences meth-dev-16S-trimmed.qza

#Create a new .qzv file to summarise and visualise the trimmed data (https://view.qiime2.org)
#Use the interactive quality plot to assess truncation parameters prior to DADA2 processing
qiime demux summarize \
--i-data meth-dev-16S-trimmed.qza \
--o-visualization meth-dev-16S-trimmed-summary.qzv

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
--o-representative-sequences meth-dev-16S-trimmed-rep-seqs.qza \
--o-denoising-stats meth-dev-16S-trimmed-denoising-stats.qza

#Create .qzv files to summarise and visualise the denoising stats, representative sequences and the ASV table (https://view.qiime2.org)
qiime metadata tabulate \
--m-input-file meth-dev-16S-trimmed-denoising-stats.qza \
--o-visualization meth-dev-16S-trimmed-denoising-stats.qzv

qiime feature-table summarize \
--i-table meth-dev-16S-trimmed-table.qza \
--o-visualization meth-dev-16S-trimmed-table.qzv \
--m-sample-metadata-file meth-dev-metadata-16S.tsv

qiime feature-table tabulate-seqs \
--i-data meth-dev-16S-trimmed-rep-seqs.qza \
--o-visualization meth-dev-16S-trimmed-rep-seqs.qzv

#The next step is to assign taxonomy to the representative sequences using a Naive Bayes classifier trained on 16S V3-V4 sequences
#from the SILVA database (Version 138) with the ‘q2-feature-classifier’ 

#First download the full length reference sequences (clustered at 99% sequence similarity) and their associated taxonomy (available from the QIIME2 website)
wget --no-check-certificate https://data.qiime2.org/2020.8/common/silva-138-99-seqs.qza
wget --no-check-certificate https://data.qiime2.org/2020.8/common/silva-138-99-tax.qza

#Extract the V3-V4 region from the full length SILVA reference sequences using the V3-V4 primer sequences prior to training the classifier
qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GGACTACNVGGGTWTCTAAT \ 
  --o-reads V3V4-ref-seqs.qza

#Train the classifier on the V3-V4 sequences
qiime feature-classifier fit-classifier-naive-bayes  \
    --i-reference-reads V3V4-ref-seqs.qza \
    --i-reference-taxonomy silva-138-99-tax.qza \
    --o-classifier silva_138_99_nb_classifier.qza
    
#Assign taxonomy to our sequences using the classifier trained on the V3-V4 sequences from SILVA 
qiime feature-classifier classify-sklearn \
--i-classifier silva-138-99-nb-classifier.qza \
--i-reads meth-dev-16S-trimmed-rep-seqs.qza \
--o-classification meth-dev-16S-trimmed-taxonomy.qza 

#Create and visualise a .qzv file of the taxonomy (https://view.qiime2.org)
qiime metadata tabulate \
--m-input-file meth-dev-16S-trimmed-taxonomy.qza \
--o-visualization meth-dev-16S-trimmed-taxonomy.qzv

#Finally, create a phylogenetic tree (this is required for importing QIIME2 artifacts into R using the qiime2R package)
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences meth-dev-16S-trimmed-rep-seqs.qza \
--o-alignment meth-dev-16S-trimmed-rep-seqs-aligned.qza \
--o-masked-alignment meth-dev-16S-trimmed-rep-seqs-masked-aligned.qza \
--o-tree meth-dev-16S-trimmed-unrooted-tree.qza \
--o-rooted-tree meth-dev-16S-trimmed-rooted-tree.qza


###ITS pre-processing###

#Make a new directory to store ITS analysis outputs
mkdir meth-dev-its

#Move into the meth-dev-its directory
cd meth-dev-its

#Make a new directory for ITS FastQC reports
mkdir fastqc-reports-its

#Run FastQC on raw ITS data to assess quality
fastqc /home/qiime2/meth-dev-its-files/* -O fastqc-reports-its

#Use MulitQC to compile a QC report of all ITS samples
multiqc -n meth-dev-its-report fastqc-reports-its

#Import demultiplexed ITS data in Casava 1.8 demultiplexed (paired-end) format to generate .qza file
qiime tools import   --type 'SampleData[PairedEndSequencesWithQuality]' \  
--input-path /home/qiime2/meth-dev-its-files \   
--input-format CasavaOneEightSingleLanePerSampleDirFmt \  
--output-path meth-dev-its.qza 

#Create a .qzv file to summarise and visualise the demultiplexed data (view the .qzv file by dragging and dropping the file at https://view.qiime2.org)
#The number of reads per sample and an interactive quality plot can be viewed to assess quality
qiime demux summarize \
--i-data meth-dev-its.qza \
--o-visualization meth-dev-its-summary.qzv

#Use the cutadapt plugin to trim ITS primer sequences from the reads and generate a new trimmed .qza file
qiime cutadapt trim-paired \
--i-demultiplexed-sequences meth-dev-its.qza \
--p-adapter-f GCATATCAATAAGCGGAGGA \
--p-front-f GCATCGATGAAGAACGCAGC \
--p-adapter-r GCTGCGTTCTTCATCGATGC \
--p-front-r TCCTCCGCTTATTGATATGC \
--o-trimmed-sequences meth-dev-its-trimmed.qza

#Create a new .qzv file to summarise and visualise the trimmed ITS data (https://view.qiime2.org)
#Use the interactive quality plot to assess truncation parameters prior to DADA2 processing
qiime demux summarize \
--i-data meth-dev-its-trimmed.qza \
--o-visualization meth-dev-its-trimmed-summary.qzv

#Use DADA2 plugin to denoise and dereplicate sequences, to remove chimeric sequences, merge reads and generate the amplicon sequence variant (ASV) 
#table, representaive sequences, and denoising stats 
#Set truncation parameters for read 1 and read 2 based on quality scores to remove low quality bases 
# The default minimum overlap required by DADA2 to merge read 1 and read 2 is 12 bases)
qiime dada2 denoise-paired \
--i-demultiplexed-seqs meth-dev-its-trimmed.qza \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 266 \
--p-trunc-len-r 187 \
--o-table meth-dev-its-trimmed-table.qza \
--o-representative-sequences meth-dev-its-trimmed-rep-seqs.qza \
--o-denoising-stats meth-dev-its-trimmed-denoising-stats.qza

#Create .qzv files to summarise and visualise the ITS denoising stats, representative sequences and the ASV table (https://view.qiime2.org)
qiime metadata tabulate \
--m-input-file meth-dev-its-trimmed-denoising-stats.qza \
--o-visualization meth-dev-its-trimmed-denoising-stats.qzv

qiime feature-table summarize \
--i-table meth-dev-its-trimmed-table.qza \
--o-visualization meth-dev-its-trimmed-table.qzv \
--m-sample-metadata-file meth-dev-metadata-its.tsv

qiime feature-table tabulate-seqs \
--i-data meth-dev-its-trimmed-rep-seqs.qza \
--o-visualization meth-dev-its-trimmed-rep-seqs.qzv

#The next step is to assign taxonomy to the representative sequences using a Naive Bayes classifier trained on full lenggth ITS sequences
#from the UNITE database  with the ‘q2-feature-classifier’ 

#First download the full length reference sequences (clustered at 99% sequence similarity) and their associated taxonomy (available from the QIIME2 website)
wget --no-check-certificate https://data.qiime2.org/2020.8/common/silva-138-99-seqs.qza
wget --no-check-certificate https://data.qiime2.org/2020.8/common/silva-138-99-tax.qza

#Extract the V3-V4 region from the full length SILVA reference sequences using the V3-V4 primer sequences prior to training the classifier
qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GGACTACNVGGGTWTCTAAT \ 
  --o-reads V3V4-ref-seqs.qza

#Train the classifier on the V3-V4 sequences
qiime feature-classifier fit-classifier-naive-bayes  \
    --i-reference-reads V3V4-ref-seqs.qza \
    --i-reference-taxonomy silva-138-99-tax.qza \
    --o-classifier silva_138_99_nb_classifier.qza
    
#Assign taxonomy to our sequences using the classifier trained on the V3-V4 sequences from SILVA 
qiime feature-classifier classify-sklearn \
--i-classifier silva-138-99-nb-classifier.qza \
--i-reads meth-dev-16S-trimmed-rep-seqs.qza \
--o-classification meth-dev-16S-trimmed-taxonomy.qza 

#Create and visualise a .qzv file of the taxonomy (https://view.qiime2.org)
qiime metadata tabulate   --m-input-file meth-dev-16S-trimmed-taxonomy.qza   --o-visualization meth-dev-16S-trimmed-taxonomy.qzv

#Finally, create a phylogenetic tree (this is required for importing QIIME2 artifacts into R using the qiime2R package)
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences meth-dev-16S-trimmed-rep-seqs.qza \
--o-alignment meth-dev-16S-trimmed-rep-seqs-aligned.qza \
--o-masked-alignment meth-dev-16S-trimmed-rep-seqs-masked-aligned.qza \
--o-tree meth-dev-16S-trimmed-unrooted-tree.qza \
--o-rooted-tree meth-dev-16S-trimmed-rooted-tree.qza








