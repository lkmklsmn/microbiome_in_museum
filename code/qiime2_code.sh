conda activate qiime2

cd /storage/thinc/data/projects/lounsbery/qiime_2

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest_seq2.txt \
  --output-path lounsbery_2.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data lounsbery_2.qza \
  --o-visualization lounsbery_2_DEMUX_SEQ_SUMMARY_VIZ.qzv

qiime cutadapt trim-paired \
--i-demultiplexed-sequences lounsbery_2.qza \
--p-cores 8 \
--p-front-f CCTACGGGWGGCWGCAG \
--p-front-r CTACCRGGGTATCTAATCC \
--o-trimmed-sequences primer-trimmed-VL_16S_PE.qza \
--verbose \
&> primer_trimming.log

qiime demux summarize \
  --i-data primer-trimmed-VL_16S_PE.qza \
  --o-visualization lounsbery_2_DEMUX_SEQ_SUMMARY_VIZ_AFTER_TRIMMING.qzv

qiime dada2 denoise-paired \
--p-n-threads 8 \
--i-demultiplexed-seqs primer-trimmed-VL_16S_PE.qza \
--p-trunc-len-f 200 \
--p-trunc-len-r 200 \
--output-dir DADA2_denoising_output \
--verbose \
&> DADA2_denoising.log

qiime metadata tabulate \
--m-input-file DADA2_denoising_output/denoising_stats.qza \
--o-visualization DADA2_denoising_output/denoising_stats.qzv

qiime feature-table tabulate-seqs \
--i-data DADA2_denoising_output/representative_sequences.qza \
--o-visualization DADA2_denoising_output/rep_seqs.qzv

qiime feature-table summarize \
--i-table DADA2_denoising_output/table.qza \
--o-visualization DADA2_denoising_output/table.qzv

qiime diversity core-metrics \
--i-table DADA2_denoising_output/table.qza \
--m-metadata-file ../sample_sheet_all.tsv \
--p-sampling-depth 100 \
--output-dir test
