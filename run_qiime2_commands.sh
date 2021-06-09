cd 

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path qiime/manifest.txt \
--output-path qiime2/lounsbery_seq2.qza \
--input-format PairedEndFastqManifestPhred64
