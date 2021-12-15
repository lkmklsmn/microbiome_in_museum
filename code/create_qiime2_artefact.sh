conda activate qiime2

cd /storage/thinc/data/projects/lounsbery/qiime_2

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest_seq2.txt \
  --output-path lounsberty_2.qza \
  --input-format PairedEndFastqManifestPhred33V2
