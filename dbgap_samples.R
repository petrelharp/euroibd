

write.table(file="Europeans.out",locations[locations$COUNTRY_SELF %in% c("Lausanne","LOLIPOP - European") & locations$GROUPING_COLLECTION=="European" & locations$STATUS_OVER=="Y",],col.names=F)



 euro<-read.table("Europeans.out")
 colnames(euro)<-c.names
 write.table(file="Europeans.out",euro)
 
 
 
 num.no.par<-apply(euro[,c("COUNTRY_MGF","COUNTRY_MGM","COUNTRY_PGF","COUNTRY_PGM")],1,function(grand){sum(grand=="")})
 
 in.Europe<-read.table("~/Dropbox/Ideas/relatedness/dbgap_In_Europe.out",as.is=T)
 europe.countries<-in.Europe[in.Europe$V1==1,2]
 
 reduced.euro<-euro[euro$COUNTRY_SELF %in% europe.countries,]
 
 
ancs.in.euro<- apply(reduced.euro[,c("COUNTRY_FATHER","COUNTRY_MOTHER","COUNTRY_MGF","COUNTRY_MGM","COUNTRY_PGF","COUNTRY_PGM")],1,function(kin){
 	all(kin %in% c(europe.countries,""))
 	
 })
 
 
  europe.count<-country.counts[names(country.counts) %in% europe.countries]
 
 
 write.table(file="keep.euro",cbind(reduced.euro$SUBJID,reduced.euro$SUBJID)[reduced.euro$STATUS_PASSED_QC2=="Y" & ancs.in.euro,],quote=F,row.nam=F,col.nam=F,sep="\t") 
 
 #including the X
 for(chr in 23:1){
 	
 	system("~/Documents/plink-1.07-mac-intel/plink --bfile POPRES_Genotypes_QC2_v2 --chr ",chr," --keep ../../keep.euro --recode --out POPRES_chr",chr)
 	
 }
 
 
 # map(regions = names(europe.count),col=heat.colors(as.numeric(europe.count)),fill=TRUE)
