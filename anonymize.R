#!/usr/bin/R

source("ibd-blocks-fns.R")
load("all-blocks-winnowed-fine.Rdata")
load("eda-data-fine.Rdata")

random.ids <- sample( max(indivinfo$SUBJID) )

anonblocks <- blocks[ c("id1","id2","chrom","maplen") ]
anonblocks$id1 <- random.ids[ anonblocks$id1 ]
anonblocks$id2 <- random.ids[ anonblocks$id2 ]
anonindivinfo <- indivinfo[ c("SUBJID", "COUNTRY_SELF") ]
anonindivinfo$SUBJID <- random.ids[ anonindivinfo$SUBJID ]

write.csv( anonblocks, file="ibd-blocklens.csv", row.names=FALSE )
write.csv( anonindivinfo, file="ibd-pop-info.csv", row.names=FALSE )
