# ## Recode the 0/1 alleles in the _legend_ file from HAPMAP
# ## so that it matches the POPRES data --
# 
# *[http://hapmap.ncbi.nlm.nih.gov/downloads/phasing/2006-07_phaseII/00README.txt phaseII phasing]: Phasing was done using the Phase I+II genotyping data files '''(non-redundant forward strand)''' / HapMap rel#21.
# *[http://hapmap.ncbi.nlm.nih.gov/downloads/genotypes/2006-07/00README.txt from genotypes/2006-07/rs_strand/non-redudant/] which is [ftp://ftp.hapmap.org/hapmap/00README.releasenotes_rel21 r21] get:
# Alleles are expressed in the same strand orientation as in dbSNP ("rs" sequence of refSNP cluster). 
# Thus, SNPs could be in the '''forward (+) or in the reverse (-) orientation relative to the reference 
# human genome''' (NCBI build 35 or UCSC hg17).
# To get from phased hapmap to popres:
# *# forward -> dbsnp: complement if strand = "-"
# *# popres agrees with dbsnp
# *the POPRES data gives '''sense relative to dbSNP''', so  can also check:
# *# dbsnp -> affy: complement if sense = "reverse"

pop <- commandArgs(TRUE)
if ( length(pop)==0 ) {
    stop("Enter a population -- CEU, YRI or JPT+CHB")
}
# pop <- "YRI"

complement <- function (x) {
    x <- factor(x, levels=c("A","C","G","T"))
    return( c("T","G","C","A")[x] )
}
dicomplement <- function (x) {
    x <- factor(x, levels=c("A/C", "A/G", "A/T", "C/G", "C/T", "G/T"))
    return( c("G/T", "C/T", "A/T", "C/G", "A/G", "A/C")[x] )
}
genotypify <- function (x,y) {
    z <- factor( paste(x,y,sep="/"), levels= c("A/A", "C/A", "G/A", "T/A", "A/C", "C/C", "G/C", "T/C", "A/G", "C/G", "G/G", "T/G", "A/T", "C/T", "G/T", "T/T") )
    return( c("A/A", "A/C", "A/G", "A/T", "A/C", "C/C", "C/G", "C/T", "A/G", "C/G", "G/G", "G/T", "A/T", "C/T", "G/T", "T/T")[z] )
}


popres <- read.table('/home/ibd/data/POPRES/original/phg000027.POPRES.genotype-calls.Affy500K.p1.MULTI.marker-info/POPRES_Snps_QC2.txt.gz', as.is=TRUE)
names(popres) <- c('snp.id', 'rs', 'alleles', 'sense')

phasedurl <- "http://hapmap.ncbi.nlm.nih.gov/downloads/phasing/2006-07_phaseII/all/"
phasedfile <- c("genotypes_chr",paste("_",pop,"_r21_nr_fwd_phased_all.gz",sep=""))
sampleurl <- "http://hapmap.ncbi.nlm.nih.gov/downloads/phasing/2006-07_phaseII/phased/"
samplefile <- c("genotypes_chr",paste("_",pop,"_r21_nr_fwd_sample.txt.gz",sep=""))
legendurl <- "http://hapmap.ncbi.nlm.nih.gov/downloads/phasing/2006-07_phaseII/all/"
legendfile <- c("genotypes_chr",paste("_",pop,"_r21_nr_fwd_legend_all.gz",sep=""))
genotypeurl <- "http://hapmap.ncbi.nlm.nih.gov/downloads/genotypes/2006-07/rs_strand/non-redundant/"
genotypefile <- c("genotypes_chr",paste("_",pop,"_r21_nr.txt.gz",sep=""))
snpinfofile <- c("snp_info_chr",paste("_",pop,"_r21.txt.gz",sep=""))
outputfile <- c("genotypes_chr",paste("_",pop,"_r21_nr_fwd_legend_all_recoded_for_POPRES",sep=""))

dl.if.needed <- function (furl,fname) {
    if (!file.exists(fname)) { download.file(furl,fname) }
}

for (chrom in 1:22) {
    print(chrom)

    samfile <- paste(samplefile[1],chrom,samplefile[2],sep="")
    phfile <- paste(phasedfile[1],chrom,phasedfile[2],sep="")
    lfile <- paste(legendfile[1],chrom,legendfile[2],sep="")
    gfile <- paste(genotypefile[1],chrom,genotypefile[2],sep="")
    sfile <- paste(snpinfofile[1],chrom,snpinfofile[2],sep="")
    ofile <- paste(outputfile[1],chrom,outputfile[2],sep="")

    # AVOID RAMPAGING ROBOT:
    Sys.sleep(60)
    dl.if.needed(paste(sampleurl,samfile,sep=""),samfile)
    dl.if.needed(paste(phasedurl,phfile,sep=""),phfile)
    dl.if.needed(paste(legendurl,lfile,sep=""),lfile)
    dl.if.needed(paste(genotypeurl,gfile,sep=""),gfile)
    # Extract just the marker information from the genotype file
    if (!file.exists(sfile)) { system(paste("zcat", gfile, "| cut -f 1-5 -d ' ' | sed -e 's/#//' | gzip -c > ", sfile)) }

    ## R script to make the map file
    hapmap.fwd <- read.table(lfile, as.is=TRUE, header=TRUE)
    hapmap.rs <- read.table(sfile, as.is=TRUE, header=TRUE)

    # careful: merge will reorder
    both <- merge(hapmap.fwd, popres, by="rs", all.x=TRUE)
    both <- merge(both, hapmap.rs, by="rs", all.x=TRUE)
    print( table(both$sense, both$strand) )
    # hapfwd and haprs should be related by strand; check this:
    # check: SNPalleles are rs, so (pop) alleles should be related by sense:
    both$rs.to.pop <- both$SNPalleles
    rs.switches <- !is.na(both$sense) & both$sense=="reverse" 
    both$rs.to.pop[ rs.switches ] <- dicomplement( both$rs.to.pop[ rs.switches ] )
    print(with(both, table(rs.to.pop, alleles)))  # yes!
    # now convert from fwd to dbsnp:
    both$flip.fwd.to.dbsnp <- !is.na(both$strand) & both$strand=="-"
    both$fwd.to.dbsnp.a0 <- both$a0
    both$fwd.to.dbsnp.a0[ both$flip.fwd.to.dbsnp ] <- complement( both$fwd.to.dbsnp.a0[ both$flip.fwd.to.dbsnp ] )
    both$fwd.to.dbsnp.a1 <- both$a1
    both$fwd.to.dbsnp.a1[ both$flip.fwd.to.dbsnp ] <- complement( both$fwd.to.dbsnp.a1[ both$flip.fwd.to.dbsnp ] )
    # # now for converting from fwd to affy:
    # both$flip.fwd.to.pop <- xor( both$strand=="-", both$sense=="reverse" )
    # both$flip.fwd.to.pop <- !is.na(both$flip.fwd.to.pop) & both$flip.fwd.to.pop 
    # both$fwd.to.pop.a0 <- both$a0
    # both$fwd.to.pop.a0[ both$flip.fwd.to.pop ] <- complement( both$fwd.to.pop.a0[ both$flip.fwd.to.pop ] )
    # both$fwd.to.pop.a1 <- both$a1
    # both$fwd.to.pop.a1[ both$flip.fwd.to.pop ] <- complement( both$fwd.to.pop.a1[ both$flip.fwd.to.pop ] )
    # both$fwd.to.pop.gen <- genotypify( both$fwd.to.pop.a0, both$fwd.to.pop.a1 )
    # print( with(both, table( alleles, fwd.to.pop.gen ) ) ) # yes!!!

    # reorder both to match hapmap.fwd 
    orig.order <- match(hapmap.fwd$position,both$position)
    if (any(is.na(orig.order))) { stop("Error! Missing some sites!") }
    both <- both[orig.order,]

    # check that ordering is consistent
    if ( any(diff(hapmap.fwd$position)<0) ) { stop("Something wrong in legend file -- sites not in order.") }
    if ( !all( hapmap.fwd$rs == both$rs ) ) { stop("Missing a site?") }
    # check that things make sense
    if (dim(both)[1] != system(paste('zcat', phfile, '| head -n 1 | wc -w'), intern=TRUE)) {
        stop( "Uh-oh: recoded legend file not the right length." )
    }

    both <- both[,c("snp.id","position","fwd.to.dbsnp.a0","fwd.to.dbsnp.a1")] 
    names(both) <- c("snp.id","position","a0","a1")
    write.table(both, file=ofile, quote=FALSE, row.names=FALSE)
}
