library(dplyr)
outfolder <- "out_dada"
metadata.in <- read.csv("metadata.csv",head=T,stringsAsFactors=F)
metadata <- metadata.in[which(metadata.in$seqid!=""),]
rownames(metadata) <- metadata$seqid
metadata <- dplyr::select(metadata, -c(seqid))

saveRDS(metadata,file.path(outfolder,"metadata.rds"))
