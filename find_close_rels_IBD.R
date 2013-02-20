##find rels

kin<-read.table(paste(.pcadir,"kin.genome.gz",sep=""),as.is=T,head=T)

close.plink<-kin[kin$PI>=0.1,]
write.table(file=paste(.pcadir,"kinship.above.cutoff",sep=""),close.plink)

