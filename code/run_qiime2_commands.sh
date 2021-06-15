cd 

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path /storage/thinc/data/projects/lounsbery/qiime/manifest.txt \
--output-path /storage/thinc/data/projects/lounsbery/qiime/lounsbery_seq2.qza \
--input-format PairedEndFastqManifestPhred33V2
