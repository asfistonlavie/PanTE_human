#### Characterization of Transposable Elements Variations in the Human Pangenome
#### Shadi Shahatit, Master's Thesis, ISEM, UM 2023-2024

#### Content:
#### 1. Pangenomics data 
#### 	1.1 Pre-processing
#### 	1.2 Processing
#### 2. Frequency distribution, TE counts, PCA, and selection scans 
#### 3. Evolutionary age estimates
#### 4. Functional genomics analysis
#### 5. Genomic windows and statistical analysis
#### Appendix I - Tools and scripts

#### 1. Pangenomics data 

#### 	1.1 Pre-processing

## download decomposed VCF with the variants in the Minigraph-Cactus HPRC pangenome
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-grch38.vcf.gz
## pop 10M bubbles and records with keep SVs
vcfbub -l 0 -r 10000000 -i hprc-v1.0-mc-grch38.vcf.gz > bubed_hprc-v1.0-mc-grch38.vcf.gz
bcftools view -i 'MAX(STRLEN(ALT))>100 | STRLEN(REF)>100' bubed_hprc-v1.0-mc-grch38.vcf.gz -Oz -o hprc-v1.0-mc.grch38.vcfbub.r10M.svsite.vcf.gz
## split multiallelic records and keep allele with largest REF or ALT
bcftools norm -m -any hprc-v1.0-mc.grch38.vcfbub.r10M.svsite.vcf.gz | grep -v "#" | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,length($4),length($5),$4,$5,$6,$7,$8,$9, $10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45,$46,$47,$48,$49,$50,$51,$52,$53,$54}' | awk '{if ($5 <= 10000){print$0}}' | gzip > hprc-v1.0-mc.grch38.vcfbub.r10M.svsite.s10000.tsv.gz
## merge with wenwei's annotations 
## wenwei's files: hprc_marker_paper.fig3c_mc_data.bed.gz   hprc_marker_paper.fig3c_pggb_data.bed.gz / Columns 1-3: SV site in GRCh38 coordinates / Column 12: repeat type
Rscript match_svs.R
## remove the header
sed '1d' hprc_marker_paper.fig3c_mc_data.svinfo.tsv > tailed_hprc_marker_paper.fig3c_mc_data.svinfo.tsv
## filter for TE annotations only
grep -e 'DNA\/misc' -e 'LINE\/L1' -e 'LINE\/mis' -e 'LTR/ERV' -e 'LTR\/misc' -e 'RC\/Helitron' -e 'Retroposon\/SVA' -e 'SINE\/Alu' -e 'SINE\/MIR' tailed_hprc_marker_paper.fig3c_mc_data.svinfo.tsv > TE_hprc_marker_paper.fig3c_mc_data.svinfo.tsv
gzip TE_hprc_marker_paper.fig3c_mc_data.svinfo.tsv
## add the vcf header
(cat vcfheader.txt; zcat TE_hprc_marker_paper.fig3c_mc_data.svinfo.tsv.gz | grep -v ^#) | bgzip -c > TE_header_hprc_marker_paper.fig3c_mc_data.svinfo.tsv.gz

#### 	1.2 Processing

## filter based on site/inds missingness
## normalize/join multiallelic records
bcftools norm -m +any TE_header_hprc_marker_paper.fig3c_mc_data.svinfo.tsv.vcf.gz  -o TE_norm_header_hprc_marker_paper.fig3c_mc_data.svinfo.tsv.vcf.gz
## split autosomes and sex chromosomes and filter missing sites of 10% in autosomes:
vcftools --gzvcf TE_norm_header_hprc_marker_paper.fig3c_mc_data.svinfo.tsv.vcf.gz --max-missing 0.9  --not-chr X --not-chr Y --recode --recode-INFO-all --stdout | gzip -c > TE_norm_mm90_header_hprc_marker_paper.fig3c_mc_data.svinfo.tsv.vcf.gz
vcftools --gzvcf TE_norm_header_hprc_marker_paper.fig3c_mc_data.svinfo.tsv.vcf.gz --chr chrX --chr chrY --recode --recode-INFO-all --stdout | gzip -c > TE_norm_sexchr.vcf.gz
bcftools concat TE_norm_mm90_header_hprc_marker_paper.fig3c_mc_data.svinfo.tsv.vcf.gz TE_norm_sexchr.vcf.gz > TE_norm_mm90.vcf
## remove undefined freq sites from chrY abnd chrX: 
## nochrY/Xpositions.txt:
## chrY	22548519
## chrY	23706156
## chrY	23862043
## chrX	58521176
vcftools --vcf TE_norm_mm90.vcf --exclude-positions nochrXYpositions.txt --recode --recode-INFO-all --stdout | gzip -c > TE_norm_mm90.vcf.gz
## check missing-indv and site - supp fig
TE_norm_header_hprc_marker_paper.fig3c_mc_data.svinfo.tsv.vcf.gz
--missing-indv
vcftools --gzvcf TE_norm_header_hprc_marker_paper.fig3c_mc_data.svinfo.tsv.vcf.gz --missing-indv --out miss_indv
--missing-site
vcftools --gzvcf TE_norm_header_hprc_marker_paper.fig3c_mc_data.svinfo.tsv.vcf.gz --missing-site --out miss_site
Rscript miss.R
## remove duplicated entries and split into two datasets of biallelic records (1 alternative) and multiallelic records (multi alternative)
grep -v '#' TE_norm_mm90.vcf  | sort -u -k1,1 -k2,2 > TE_norm_mm90_nodup.vcf
## output TEvar numbers: 20 457 > 20 346
awk '{if ($5 ~ /,/) {print$0}}' TE_norm_mm90_nodup.vcf > TE_norm_mm90_nodup_multialt.vcf
awk '{if ($5 !~ /,/) {print$0}}' TE_norm_mm90_nodup.vcf > TE_norm_mm90_nodup_1alt.vcf
(cat headervcf.txt; cat TE_norm_mm90_nodup_1alt.vcf | grep -v ^#) | bgzip -c > TE_norm_mm90_nodup_1alt.vcf.gz
## output TEvar numbers: 20 346 = 13 939 + 6 407

#### 2. Frequency distribution, TE counts per individual, PCA, and selection scans

## frequency distribution

## calculate population specific allele frequencies: global population (AFR + outofAFR) 
bcftools view -S AFR.txt TE_norm_mm90_nodup_1alt.vcf > AFR_TE_norm_mm90_1alt.vcf
bcftools view -S outofAFR.txt TE_norm_mm90_nodup_1alt.vcf > outofAFR_TE_norm_mm90_1alt.vcf
vcftools --vcf AFR_TE_norm_mm90_1alt.vcf --freq --out freqAFR
vcftools --vcf outofAFR_TE_norm_mm90_1alt.vcf --freq --out freqoutofAFR
vcftools --vcf TE_norm_mm90_nodup_1alt.vcf --freq --out freqSuper
Rscript PanTE_human_Rscript.R

## TE counts per individual

bgzip TE_norm_mm90_nodup_1alt.vcf
bcftools sort TE_norm_mm90_nodup_1alt.vcf.gz -o sorted_TE_norm_mm90_nodup_1alt.vcf.gz
bcftools index --tbi --force sorted_TE_norm_mm90_nodup_1alt.vcf.gz 
## Filter for major
vcffilter -f "AF > 0.5" sorted_TE_norm_mm90_nodup_1alt.vcf.gz > major_TE_filtered.vcf 
bcftools stats -s- major_TE_filtered.vcf > major_TE_stat
## Filter for singlitons
vcftools --vcf TE_norm_mm90_nodup_1alt.vcf --singletons --out singletons_TE_filtered.vcf
cut -f5 singletons_TE_filtered.vcf.singletons | sort | uniq -c > R_singletons_TE_filtered.vcf.singletons
## Filter for shared
bcftools view -g ^miss TE_norm_mm90_1alt.vcf > notmiss_TE_norm_mm90_1alt.vcf 
bcftools view -e 'GT="RR"' notmiss_TE_norm_mm90_1alt.vcf > nonref_TE_filtered.vcf
bcftools view -H nonref_TE_filtered.vcf > nonref_TE_filtered_nohead.vcf
Rscript PanTE_human_Rscript.R

## PCA

## remove linked sites with plink
bcftools sort TE_norm_mm90_nodup_1alt.vcf > TE_sorted.vcf
plink --vcf TE_sorted.vcf --double-id --allow-extra-chr --vcf-half-call m --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out TE_linkage
## 50: a window of 50 Kb
## 10: a window step size - move 10 bp each time we linkage is calculated
## 0.1: r2 threshold, the threshold of linkage we are willing to tolerate removing what is greater than 0.1
## run PCA with prune in plink
plink --vcf TE_sorted.vcf --double-id --allow-extra-chr -vcf-half-call m --set-missing-var-ids @:# --extract TE_linkage.prune.in --make-bed --pca --out TE_pca ## --make-bed in case of admixture analysis
## plot the eigen vectors and values of the PCA wiht R
Rscript PanTE_human_Rscript.R

## frequency based selection scans

## calculate site specific fst
vcftools --vcf TE_norm_mm90_nodup_1alt.vcf --not-chr chrX --not-chr chrY --recode --recode-INFO-all --stdout | gzip -c > TE_norm_mm90_nodup_nosex_1alt.vcf.gz 
vcftools --vcf TE_norm_mm90_nodup_nosex_1alt.vcf --weir-fst-pop AFR.txt --weir-fst-pop outofAFR.txt	--keep AFR.txt --keep outofAFR.txt
Rscript PanTE_human_Rscript.R

## haplotype based selection scans

## get vcf file data for each pop with no sex chr
## AFR
vcftools --vcf AFR_TE_norm_mm90_1alt.vcf --max-missing 1.0 --not-chr chrX --not-chr chrY --recode --recode-INFO-all --stdout | gzip -c > AFR_TE_norm_mm90_1alt_filtered.vcf.gz 
gunzip AFR_TE_norm_mm90_1alt_filtered.vcf.gz  
awk '{ flag=0; for (i=1; i<=NF; i++) { if ($i == 0) { flag=1; break; } } } !flag' AFR_TE_norm_mm90_1alt_filtered.vcf > AFR_1.vcf
awk '{ flag=0; for (i=1; i<=NF; i++) { if ($i == 1) { flag=1; break; } } } !flag' AFR_1.vcf > AFR_2.vcf
bcftools sort -Oz AFR_2.vcf -o AFR_TE_norm_mm90_1alt.vcf.gz
## outofAFR
vcftools --vcf outofAFR_TE_norm_mm90_1alt.vcf --max-missing 1.0 --not-chr chrX --not-chr chrY --recode --recode-INFO-all --stdout | gzip -c > outofAFR_TE_norm_mm90_1alt_filtered.vcf.gz 
gunzip outofAFR_TE_norm_mm90_1alt_filtered.vcf.gz 
awk '{ flag=0; for (i=1; i<=NF; i++) { if ($i == 0) { flag=1; break; } } } !flag' outofAFR_TE_norm_mm90_1alt_filtered.vcf > outofAFR_1.vcf
awk '{ flag=0; for (i=1; i<=NF; i++) { if ($i == 1) { flag=1; break; } } } !flag' outofAFR_1.vcf > outofAFR_2.vcf
bcftools sort -Oz outofAFR_2.vcf -o outofAFR_TE_norm_mm90_1alt.vcf.gz
## seperate by chr
## AFR
for chr in chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr2 chr20 chr21 chr22 chr3 chr4 chr5 chr6 chr7 chr8 chr9;
do
zcat AFR_TE_norm_mm90_1alt.vcf.gz | grep -v "^#" | \
grep $chr | cut -f1-2 | \
awk '{print $1 "\t" $1"-"$2 "\t" $2 "\t" $2}' > $chr\_selscan.map;
done
## outofAFR
for chr in chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr2 chr20 chr21 chr22 chr3 chr4 chr5 chr6 chr7 chr8 chr9;
do
zcat outofAFR_TE_norm_mm90_1alt.vcf.gz | grep -v "^#" | \
grep $chr | cut -f1-2 | \
awk '{print $1 "\t" $1"-"$2 "\t" $2 "\t" $2}' > $chr\_selscan.map;
done
## create a vcf for each chr and pop
for chr in chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr2 chr20 chr21 chr22 chr3 chr4 chr5 chr6 chr7 chr8 chr9;
do
vcftools --gzvcf AFR_TE_norm_mm90_1alt.vcf.gz --chr ${chr} \
--remove-filtered-all --recode \
--stdout | gzip -c > AFR_1alt_phased_${chr}.vcf.gz;
vcftools --gzvcf outofAFR_TE_norm_mm90_1alt.vcf.gz --chr ${chr} \
--remove-filtered-all --recode \
--stdout | gzip -c > outofAFR_1alt_phased_${chr}.vcf.gz;
done
Rscript PanTE_human_Rscript.R

## ohana demography based selection scans

## remove chm13 ref and sex chr
vcftools --vcf TE_norm_mm90_nodup_1alt.vcf --remove-indv CHM13 --not-chr chrX --not-chr chrY --recode --recode-INFO-all --stdout | gzip -c > TE_norm_mm90_nochm13.vcf.gz

plink --vcf TE_norm_mm90_nochm13.vcf.gz --make-bed --vcf-half-call m --out MergedVariants
admixture -j8 --cv MergedVariants.bed 2

## genotype observations for TE_norm_mm90_1alt.vcf.gz 

plink --vcf TE_norm_mm90_nochm13.vcf.gz --recode12 --geno 0.0 --tab --vcf-half-call m --out TE
head -n 3 TE.ped | cut -f1-12
convert ped2dgm TE.ped geno.dgm

qpas geno.dgm -k 2 -qo q.matrix -fo f.matrix -mi 5
head -n 4 q.matrix
head -n 4 f.matrix | cut -c1-50
python plot-q.py q.matrix q-bar-chart.pdf

nemeco geno.dgm f.matrix -co c.matrix -mi 5
convert cov2nwk c.matrix tree.nwk
convert nwk2svg tree.nwk tree.svg

python sample-sites.py geno.dgm 5 geno_subset.dgm
qpas geno_subset.dgm -e 0.0001 -k 2 -qo Q_subset.matrix -fo F_subset.matrix
nemeco geno_subset.dgm F_subset.matrix -e 0.0 -co C_subset.matrix

qpas geno.dgm -qi Q_subset.matrix -fo F_subset.matrix -e 0.0001 -fq

selscan geno.dgm f.matrix C_subset.matrix > scansel.txt
Rscript PanTE_human_Rscript.R

#### 3. Evolutionary age estimates

## run repeamasker for hg38 vcf biallelec file:
sed -e 's/chr//' TE_norm_mm90_nodup_1alt.vcf | awk '{OFS="\t"; if (!/^#/){print "chromosome"$1"_"$2"_"$3,$5}}' | awk '{ gsub(/>/, "_") }1' | awk '{OFS="\t";{print ">"$1,$2}}' > vcf2fasta_hg38_1.txt
awk '{OFS=RS;$1=$1}1' vcf2fasta_hg38_1.txt > vcf2fasta_hg38.txt
RepeatMasker -a -gff –Species Human vcf2fasta_hg38.txt
## calcDivergenceFromAlign
/home/shadi/anaconda3/envs/projects/share/RepeatMasker/util/calcDivergenceFromAlign.pl -s genome.divsum vcf2fasta_hg38.txt.align 
## createRepeatLandscape
Gunzip hg38.fa.gz
faToTwoBit hg38.fa hg38.2bit
twoBitInfo hg38.2bit stdout | sort -k2rn > hg38.2bit.chrom.sizes
/home/shadi/anaconda3/envs/projects/share/RepeatMasker/util/createRepeatLandscape.pl -div genome.divsum -twoBit hg38.2bit > genome.html
## build Summary
/home/shadi/anaconda3/envs/projects/share/RepeatMasker/util/buildSummary.pl vcf2fasta_hg38.txt.out
## calcDivergenceFromAlign
/home/shadi/anaconda3/envs/projects/share/RepeatMasker/util/calcDivergenceFromAlign.pl -s vcf2fasta_h.txt.tbl  vcf2fasta_hg38.txt.out
## exctart K distance
grep -ioE "Kimura \(with divCpGMod\) = [0-9]+|chromosome[0-9XY]+_[0-9]+__[0-9]+_[0-9]+" vcf2fasta_hg38.txt.align > Kdistance.txt
awk '!seen[$0]++ || !/^chromosome/' Kdistance.txt > Kdistance_fil.txt
sed 's/ /_/g' Kdistance_fil.txt > Kdistance_filtered.txt
Rscript PanTE_human_Rscript.R

#### 4. Functional genomics analysis

## install annovar
wget  http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
tar -xvfz annovar.latest.tar.gz
## run annovar
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp30a humandb/
## gene based
perl convert2annovar.pl -format vcf4 -includeinfo -allsample -withfreq -outfile TEvar.annovar TE_norm_mm90_nodup_1alt.vcf
chmod a+x annotate_variation.pl
chmod a+x coding_change.pl
perl table_annovar.pl TEvar.annovar humandb/ --outfile anno_g --buildver hg38 --protocol refGene --operation g 
cut -f 4 anno_g.hg38_multianno.txt | sort | uniq -c 
cut -f 1,2,6 anno_g.hg38_multianno.txt > R_1alt_annovar_out.txt 
Rscript PanTE_human_Rscript.R

## extract candidate TE loci from selection scans and annotate with AnnotSV and VEP 

## download candi_all.bed from PanTE_human_Rscript.R
sed '1d' candi_all.bed | awk -F'\t' 'BEGIN {OFS="\t"} {$3=""; $1=$1; print}' > candi_all_filtered.bed
vcftools --vcf TE_norm_mm90_nodup_1alt.vcf --positions candi_all_filtered.bed --recode --recode-INFO-all --stdout | gzip -c > TE_candi_all.vcf.gz
gunzip TE_candi_all.vcf.gz
## upload in AnnotSV and VEP 

#### 5. Genomic windows, recombination rate, and statistical analysis

## background intervals for recombination rate posistion
## the data has 3713 entries of genomic windows
bedtools random -g hg38_chrom_sizes.tsv -n 3713 -l 10000 -seed 123 > random_3713_hg38.bed
bedtools shuffle -i random_3713_hg38.bed -g hg38_chrom_sizes.tsv -excl foreground_TE.bed > background_3713_hg38.bed
Rscript PanTE_human_Rscript.R

#### Appendix I - Tools and scripts

## conda v24.1.2
## vcfbub v0.1.0
## vcftools v0.1.17
## bcftools v1.19
## bedtools v2.25.0 
## plink v1.90b6.21
## RepeaMasker v4.1.5
## ohana v0.1
## ANNOVAR v3.3.9
## AnnotSV https://lbgi.fr/AnnotSV/
## VEP https://www.ensembl.org/info/docs/tools/vep/index.html 
## Public scripts; RepeaMasker: calcDivergenceFromAlign.pl, createRepeatLandscape.pl, buildSummary.pl; ANNOVAR: annotate_variation.pl, table_annovar.pl; Ohana: sample-sites.py


