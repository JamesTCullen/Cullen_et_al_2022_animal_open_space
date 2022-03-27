#Demultiplexed paired-end 16S and ITS2 rDNA sequences (available at: https://www.ebi.ac.uk/ena/browser/view/PRJEB49004) were 
#imported (in Casava 1.8 demultiplexed paired-end format) into QIIME2 v.2020.8.0 (Bolyen et al., 2019), which was installed on 
#a virtual machine (VirtualBox 6.0). Forward and reverse reads were quality assessed using the ‘qiime demux summarize’ command, 
#FastQC v.0.11.5 and MultiQC v.1.9. The 16S and ITS2 primers were removed from reads using the cutadapt plugin (Martin, 2011). 
#The QIIME2 DADA2 (Callahan et al., 2016) plugin was used for filtering and dereplication, chimera removal, merging paired-end 
#reads and to infer amplicon sequence variants (ASVs) in each sample after truncating reads to remove low quality bases. 
#Read 1 and read 2 of the 16S rDNA sequences were truncated at 267 and 183 bp, respectively, while ITS2 sequences were truncated 
#at 266 and 187 bp, respectively. For bacterial sequences, taxonomy was assigned to each ASV using a Naive Bayes classifier 
#trained on 16S V3-V4 sequences from the Silva database (Version 138) with the ‘q2-feature-classifier’ plugin, while taxonomy 
#was assigned to fungal ASVs using a Naive Bayes classifier trained on full-length ITS sequences from the UNITE v.8.3 database (Kõljalg et al., 2013).


#Import demultiplexed data in Casava 1.8 demultiplexed (paired-end) format
qiime tools import   --type 'SampleData[PairedEndSequencesWithQuality]'   --input-path /home/qiime2/method-dev-files/   --input-format CasavaOneEightSingleLanePerSampleDirFmt   --output-path 16S-meth-dev.qza
