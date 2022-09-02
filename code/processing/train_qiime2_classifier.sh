#conda activate qiime2-2021.4
cd /storage/thinc/data/projects/lounsbery/genomes/rlappan/SILVA

qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna \
--output-path 99_otus_16S.qza

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path SILVA_132_QIIME_release/taxonomy/16S_only/99/majority_taxonomy_7_levels.txt \
--output-path 99_otus_16S_taxonomy.qza

qiime feature-classifier extract-reads \
--i-sequences 99_otus_16S.qza \
# 341f
--p-f-primer CCTACGGGWGGCWGCAG \
# 803r 
--p-r-primer CTACCRGGGTATCTAATCC \
--p-min-length 285 \
--p-max-length 600 \
--o-reads ref_seqs_16S_V3-V4.qza \
--verbose \
&> 16S_V3-V4_training.log

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads 99_otus_16S.qza \
--i-reference-taxonomy 99_otus_16S_taxonomy.qza \
--o-classifier classifier_16S_V3-V4.qza \
--verbose \
&> 16S_V3-V4_classifier.log
