conda activate qiime2

#wget https://data.qiime2.org/2021.8/common/silva-138-99-seqs.qza
#wget https://data.qiime2.org/2021.8/common/silva-138-99-tax.qza

cd /storage/thinc/data/projects/lounsbery/qiime/classifiers/silva_138

qiime feature-classifier extract-reads \
--i-sequences silva-138-99-seqs.qza \
--p-f-primer AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
--p-r-primer AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
--p-n-jobs 4 \
--p-min-length 100 \
--p-max-length 400 \
--o-reads ref-seqs.qza

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs.qza \
--i-reference-taxonomy silva-138-99-tax.qza \
--o-classifier classifier.qza
