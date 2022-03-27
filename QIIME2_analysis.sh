#Import demultiplexed data in Casava 1.8 demultiplexed (paired-end) format
qiime tools import   --type 'SampleData[PairedEndSequencesWithQuality]'   --input-path /home/qiime2/method-dev-files/   --input-format CasavaOneEightSingleLanePerSampleDirFmt   --output-path 16S-meth-dev.qza
