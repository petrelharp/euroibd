#!/usr/bin/make
#
# Caveat utilitorum: This is intended as a guide. Making everything from this make file will at least result in syntax errors and long run times.
# Some of the below is pseudo-make code.  If you are serious about doing this and can't figure out what's supposed to happen, drop us a line.
# 
# Inputs not included here:
# Europeans.out -- sample info from POPRES
# POPRES data: transformed to BEAGLE input and named POPRES_chr%.beagle.gz
# HAPMAP trio-phased data
# europeans.imiss made with `plink --file data --missing `
# kin.genome.gz is also made with plink.

.PHONY : real-data false-pos-sims true-pos-sims chrom1 chrom2 chrom3 chrom4 chrom5 chrom6 chrom7 chrom8 chrom9 chrom10 chrom11 chrom12 chrom13 chrom14 chrom15 chrom16 chrom17 chrom18 chrom19 chrom20 chrom21 chrom22 

# Figures:

various-substructure.pdf sharing-rates-and-maps.pdf overlap-rates-all-chr-onepage.pdf substructure_summaries.pdf inversion-distributions.pdf inversion-distributions-coalrate.pdf inversion-boxplots-long.pdf inversion-boxplots.pdf inversion-boxplots-long-coal.pdf boxplotted-inversions.pdf selected-inversions.pdf example-inversion-bounds.pdf tinyinversions.pdf power-and-fp.pdf : inversion-writeup-plots.R eda-data-fine.Rdata all-blocks-winnowed-fine.Rdata overlaps.RData lentab.Rdata
	R CMD BATCH $<

spectra-comparisons.pdf full-inversions.pdf long-inversions.pdf complex-inversions.pdf age-distributions.pdf fp-inversions.pdf : sensitivity.R
	R CMD BATCH $<

pca-map.pdf Italy_PCA_vs_IBD.pdf UK_PCA_vs_IBD.pdf : PCA_vs_IBD.R eda-data-fine.Rdata
	R CMD BATCH $<

IBS_vs_dist.pdf IBS_vs_dist.svg : IBS_vs_dist.R Euro-samples-info.tsv kin.genome.gz
	R CMD BATCH $<

indiv-sharing-map-of-europe.pdf : ibd-pca.R all-blocks-winnowed-fine.Rdata eda-data-fine.Rdata
	R CMD BATCH $<

error-fit-plots.pdf : fit-error-model.R true-pos-everything.Rdata
	R CMD BATCH $<

false-pos-and-true-rate.pdf : beagle-false-positives.R false-pos-everything.Rdata true-pos-sims false-pos-sims
	# this one not so clean; beware.
	R CMD BATCH $<  

ibd-pop-info.csv ibd-blocklens.csv : anonymize.R eda-data-fine.Rdata all-blocks-winnowed-fine.Rdata
	# Anonymized data.
	R CMD BATCH $<

# Pre-processing of data

chrom% : POPRES_chr%.beagle.gz
	# WARNING TAKES A LOOOONG TIME:
	run-beagle $< 10

POPRES_chr%.combined.fibd.gz : *.POPRES_chr%.beagle.fibd.gz
	python plus-process-fibd.py 1e-7 5.0 marker.genetic%.gmap $+ | gzip -c > $@

POPRES_chr%.combined.winnowed.fibd.gz : POPRES_chr%.combined.fibd.gz 
	( zcat $< | head -n 1; zcat $< | awk '$6 <= 1E-9 && $5 >= 2 { print; }' ) | gzip -c > $@

real-data : POPRES_chr*.combined.winnowed.fibd.gz

true-pos-sims : recode-hapmap-legend-for-popres.R hapmap-phased-to-beagle.py make-hapmap-plus-popres.sh insert-ibd-segments.py $(gmaps)
	R CMD BATCH recode-hapmap-legend-for-popres.R
	POP=CEU
	for x in $(seq 22)
	do
		python hapmap-phased-to-beagle.py \
			-l genotypes_chr${x}_${POP}_r21_nr_fwd_legend_all_recoded_for_POPRES \
			-s genotypes_chr${x}_${POP}_r21_nr_fwd_sample.txt.gz \
			-p genotypes_chr${x}_${POP}_r21_nr_fwd_phased_all.gz \
			-o HAPMAP_PHASED_chr${x}_${POP}_r21_nr_all.beagle.gz \
			-m marker.genetic${x}.gmap
		PFIX=$RANDOM
		python insert-ibd-segments.py \
			-i HAPMAP/HAPMAP_PHASED_chr${x}_CEU_r21_nr_all.beagle.gz \
			-j temp.chr${x}.beagle.gz -m marker.genetic${x}.gmap \
			-o HAPMAP_POPRES_${PFIX}_chr${x}.beagle.gz \
			-g HAPMAP_POPRES_${PFIX}_chr${x}.TRUE.ibd.gz
		zcat HAPMAP_POPRES_${PFIX}_chr${x}.beagle.gz | tail -n +2 | cut -f 2 -d ' ' | sort >hapmap.popres.chr${x}.snplist
		(head -n 1 marker.genetic${x}.gmap; cat marker.genetic${x}.gmap | tail -n +2 | nl | sort -k 3b,3 \
			| join -1 1 -2 3 hapmap.popres.chr${x}.snplist - | sort -k 2n,2 \
			| awk '{ printf $3" "$1" "$4" "$5"\n" }' ) >hapmap.popres.genetic${x}.gmap
		rm hapmap.popres.chr${x}.snplist
		CHRBASENAME="HAPMAP_POPRES_${PFIX}_chr${x}"
		# remove gaps no longer than 5cM
		python /home/peter/projects/genome/plus-process-fibd.py 1e-7 5.0 hapmap.popres.genetic${x}.gmap \
			*${CHRBASENAME}.beagle*.fibd.gz | gzip -c > ${CHRBASENAME}.combined.fibd.gz
		( zcat ${CHRBASENAME}.combined.fibd.gz | head -n 1; \
			zcat ${CHRBASENAME}.combined.fibd.gz | awk '$6 <= 1E-9 && $5 >= 2 { print; }' ) \
			| gzip -c > ${CHRBASENAME}.combined.winnowed.fibd.gz
	done

lentab.Rdata : make-lentab.R
	R CMD BATCH $<

false-pos-sims : mixup-by-region.sh mixup-genomes.py Euro-samples-info.tsv $(gmaps)
	# see mixup-by-region.sh and do as in true-pos-sims

all-blocks-winnowed-fine.Rdata eda-data-fine.Rdata : blocks-eda.R real-data
	R CMD BATCH blocks-eda.R

overlaps.RData : endpoints-and-overlap.R all-blocks-winnowed.Rdata
	R CMD BATCH endpoints-and-overlap.R

true-pos-everything.Rdata : beagle-true-positives.R true-pos-sims
	R CMD BATCH beagle-true-positives.R

false-pos-everything.Rdata : beagle-false-positives.R false-pos-sims
	R CMD BATCH beagle-false-positives.R

# R script dependencies

blocks-eda.R : ibd-blocks-fns.R

PCA_vs_IBD.R : ibd-blocks-fns.R

ibd-pca.R : ibd-blocks-fns.R laplace-inversion-fns.R

fit-error-model.R : ibd-blocks-fns.R

beagle-true-positives.R : ibd-blocks-fns.R

beagle-false-positives.R : ibd-blocks-fns.R

inversion-writeup-plots.R : laplace-inversion-fns.R ibd-blocks-fns.R 

endpoints-and-overlap.R : ibd-blocks-fns.R laplace-inversion-fns.R 

sensitivity.R : sim-params.R parse-sims-fns.R ibd-blocks-fns.R laplace-inversion-fns.R

Euro-samples-info-fine.tsv citypair.dists.tsv : Europeans.out clean-sample-info.R dbgap_In_Europe.out kinship.above.cutoff euro_nooutlier.pcavec europeans.imiss
	R CMD BATCH clean-sample-info.R

# misc

gmaps = marker.genetic1.gmap marker.genetic2.gmap marker.genetic3.gmap marker.genetic4.gmap marker.genetic5.gmap marker.genetic6.gmap marker.genetic7.gmap marker.genetic8.gmap marker.genetic9.gmap marker.genetic10.gmap marker.genetic11.gmap marker.genetic12.gmap marker.genetic13.gmap marker.genetic14.gmap marker.genetic15.gmap marker.genetic16.gmap marker.genetic17.gmap marker.genetic18.gmap marker.genetic19.gmap marker.genetic20.gmap marker.genetic21.gmap marker.genetic22.gmap 

$(gmaps) : genetic.maps.R female.gmap.gz male.gmap.gz 
	R CMD BATCH genetic.maps.R

kinship.above.cutoff : find_close_rels_IBD.R kin.genome.gz
	R CMD BATCH find_close_rels_IBD.R
