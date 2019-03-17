# RENORMALIZE RAW READS BASED ON TAXONOMY

library(dada2); packageVersion("dada2")

#* 1) 16S/bacteria+archaea
# 2) 16S/chloroplasts
# 3) 16S/mitochondria --> possible to split further into protists/metazoa?
# 4) 16S/eukaryotes

# 5) 18S/metazoa
#* 6) 18S/protists

# 7) ITS/all

# combine multiple primers and multiple runs
# need to add quality flag, figure out how to deal with different naming schemes
# can't base plate well off of run name
source("microscape-s00-setup.r")
# 
# cleannames <- function(x) {
#   x <- gsub(patt="18sPlate",repl="18S-Plate",x=x)
#   x <- gsub(patt="16sPlate",repl="16S-Plate",x=x)
#   x <- gsub(patt="16SPlate",repl="16S-Plate",x=x)
#   x <- gsub(patt="ITSPlate",repl="ITS-Plate",x=x)
#   x <- gsub(patt="Plates",repl="Plate",x=x)
#   x <- gsub(patt="Plate-",repl="Plate",x=x)
#   x
# }
# seqtab <- readRDS(file.path(outfolder,"seqtab_final.rds")) # 
# rownames(seqtab) <- cleannames(rownames(seqtab))
# saveRDS(seqtab,file="seqtab_final.rds")

# read from disk
seqtab <- readRDS(file.path(outfolder,"seqtab_final.rds")) # 
rownames(seqtab) <- toupper(rownames(seqtab))

taxout <- readRDS(file.path(outfolder,"taxout.rds"))

metadata <- readRDS(file.path(outfolder,"metadata.rds"))
rownames(metadata) <- toupper(rownames(metadata))

#20160323/16S-Plate1Pool-Rep1_A01 seqtab
#add Rep1

#20160619/16S-Plate5_A01 seqtab
#20160619/16S-Plate6_A01 seqtab

#20160619/18S-Plate1_A01	seqtab
#20160619/18S-Plate3_A01	seqtab
#where's 18S-Plate2?
#where's 18S-Plate4?

#20170528/16S-Plate-L1ND_A02	seqtab
#20170528/18S-Plate-L2ND_A01	seqtab
#20170528/16S-Plate-L4ND_A09	seqtab
#20170528/18S-Plate-L5ND_A01	seqtab

# 20180616/ITSPlate-5Redo_A01 seqtab

# 20180611/16SPlate-17_A01 metadata
# 20180611/16sPlates-17_A01 seqtab

#20180510/ITSPlate-10_A01	seqtab
#20180521/ITSPlate-10_A01	metadata


# should be 0
rownames(metadata)[duplicated(rownames(metadata))]

# how many have valid sequences but not tags
length(setdiff(rownames(seqtab),rownames(metadata)))
# 247
# 421

# how many have valid tags but not sequences
length(setdiff(rownames(metadata),rownames(seqtab)))
# 1661
# 2031

length(intersect(rownames(seqtab),rownames(metadata)))
# 6055
# 5881

# define sample sets
# samples_ITS <- which(grepl("/ITS",rownames(seqtab),ignore.case=T))
# samples_16S <- which(grepl("/16S",rownames(seqtab),ignore.case=T))
# samples_18S <- which(grepl("/18S",rownames(seqtab),ignore.case=T))
samples_ITS <- intersect(rownames(seqtab),rownames(metadata)[metadata$primer=="ITS"])
samples_16S <- intersect(rownames(seqtab),rownames(metadata)[metadata$primer=="16S"])
samples_18S <- intersect(rownames(seqtab),rownames(metadata)[metadata$primer=="18S"])
intersect(samples_ITS,intersect(samples_16S,samples_18S)) #should be character(0)

# subset sequence table by primer
seqtab_ITS <- seqtab[samples_ITS,]
seqtab_16S <- seqtab[samples_16S,]
seqtab_18S <- seqtab[samples_18S,]

# subset sequence tables by taxa
cs_ITS <- colSums(seqtab_ITS>0)
cs_16S <- colSums(seqtab_16S>0)
cs_18S <- colSums(seqtab_18S>0)
cs_all <- unname(cbind(cs_ITS,cs_16S,cs_18S))

pick_ITS <- unname(which(cs_ITS > cs_16S+cs_18S))
pick_16S <- unname(which(cs_16S > cs_ITS+cs_18S))
pick_18S <- unname(which(cs_18S > cs_ITS+cs_16S))
pick_none <- unname(setdiff(1:ncol(seqtab),c(pick_ITS,pick_16S,pick_18S))) #should be integer(0) unless some missing from metadata

# re-subset sequence table by primer
# set all non-primer sequences to zero to keep order
# these are probably bad demultiplexing or index switching?
seqtab_ITS <- seqtab[samples_ITS,]
seqtab_ITS[,c(pick_16S,pick_18S,pick_none)] <- 0
rownames(seqtab_ITS) <- metadata[rownames(seqtab_ITS),"cellid"]
seqtab_16S <- seqtab[samples_16S,]
seqtab_16S[,c(pick_ITS,pick_18S,pick_none)] <- 0
rownames(seqtab_16S) <- metadata[rownames(seqtab_16S),"cellid"]
seqtab_18S <- seqtab[samples_18S,]
seqtab_18S[,c(pick_16S,pick_ITS,pick_none)] <- 0
rownames(seqtab_18S) <- metadata[rownames(seqtab_18S),"cellid"]

seqtab_none <- seqtab
seqtab_none[,c(pick_16S,pick_18S,pick_ITS)] <- 0
rownames(seqtab_none) <- metadata[rownames(seqtab),"cellid"]

# subset sequences by taxonomy

table_list <- (1:ncol(seqtab))*NA
proptab_all <- list()

# separate mets and organelles /16S

eukaryote_esv <- intersect(pick_16S,unname(which(taxout$ref_dada2_silva_nr_v132_train_set.fasta[,1]=="Eukaryota")))
seqtab_euk <- seqtab_16S[,eukaryote_esv]
proptab_all[["eukaryote_esv"]] <- norm_esv(seqtab_euk)
table_list[eukaryote_esv] <- "16S_eukaryote"

chloroplast_esv <- intersect(pick_16S,unname(which(taxout$ref_dada2_silva_nr_v132_train_set.fasta[,4]=="Chloroplast")))
seqtab_chloro <- seqtab_16S[,chloroplast_esv]
proptab_all[["chloroplast_esv"]] <- norm_esv(seqtab_chloro)
table_list[chloroplast_esv] <- "16S_chloroplast"

mitochondria_esv <- intersect(pick_16S,unname(which(taxout$ref_dada2_silva_nr_v132_train_set.fasta[,5]=="Mitochondria")))
seqtab_mito <- seqtab_16S[,mitochondria_esv]
proptab_all[["mitochondria_esv"]] <- norm_esv(seqtab_mito)
table_list[mitochondria_esv] <- "16S_mitochondria"

prokaryote_esv <- intersect(pick_16S,unname(which(taxout$ref_dada2_silva_nr_v132_train_set.fasta[,1]=="Bacteria" | taxout$ref_dada2_silva_nr_v132_train_set.fasta[,1]=="Archaea")))
prokaryote_esv <- setdiff(prokaryote_esv,c(eukaryote_esv,chloroplast_esv,mitochondria_esv))
seqtab_prok <- seqtab_16S[,prokaryote_esv]
proptab_all[["prokaryote_esv"]] <- norm_esv(seqtab_prok)
table_list[prokaryote_esv] <- "16S_prokaryote"

unknown_16S_esv <- intersect(pick_16S,setdiff(1:ncol(seqtab_16S),c(prokaryote_esv,eukaryote_esv,chloroplast_esv,mitochondria_esv)))
seqtab_16S_unknown <- seqtab_16S[,unknown_16S_esv]
#unname(taxout$ref_dada2_silva_nr_v132_train_set.fasta[colnames(seqtab_16S_unknown),])
proptab_all[["unknown_16S_esv"]] <- norm_esv(seqtab_16S_unknown)
table_list[unknown_16S_esv] <- "16S_unknown"

# separate metazoa and protists /18S, normalize to ppmillion
prokaryote_18S_esv <- intersect(pick_18S,unname(which(taxout$ref_dada2_silva_nr_v132_train_set.fasta[,1]=="Bacteria" | taxout$ref_dada2_silva_nr_v132_train_set.fasta[,1]=="Archaea")))
seqtab_prok_18S <- seqtab_18S[,prokaryote_18S_esv]
proptab_all[["prokaryote_18S_esv"]] <- norm_esv(seqtab_prok_18S)
table_list[prokaryote_18S_esv] <- "18S_prokaryote"

metazoa_esv <- intersect(pick_18S,unname(which(taxout$ref_dada2_silva_v132.fasta[,5]=="Metazoa_(Animalia)")))
seqtab_meta <- seqtab_18S[,metazoa_esv]
proptab_all[["metazoa_esv"]] <- norm_esv(seqtab_meta)
table_list[metazoa_esv] <- "18S_metazoa"

protist_esv <- intersect(pick_18S,setdiff(1:ncol(seqtab_18S),c(prokaryote_18S_esv,metazoa_esv)))
seqtab_prot <- seqtab_18S[,protist_esv]
proptab_all[["protist_esv"]] <- norm_esv(seqtab_prot)
table_list[protist_esv] <- "18S_protist"

unknown_18S_esv <- intersect(pick_18S,setdiff(1:ncol(seqtab_18S),c(prokaryote_18S_esv, metazoa_esv, protist_esv)))
seqtab_18S_unknown <- seqtab_18S[,unknown_18S_esv]
proptab_all[["unknown_18S_esv"]] <- norm_esv(seqtab_18S_unknown)
table_list[unknown_18S_esv] <- "18S_unknown"

# ITS databases suck so just keep them all together
proptab_all[["ITS_esv"]] <- norm_esv(seqtab_ITS[,pick_ITS])
table_list[pick_ITS] <- "ITS_all"

# # the remainder (have equal number of samples hit to 2 different primers?)
unknown_esv <- which(is.na(table_list))
seqtab_unknown <- seqtab_none[,pick_none]
proptab_all[["unknown_esv"]] <- norm_esv(seqtab_unknown)
# table_list[is.na(table_list)] <- "unknown"



# remerge them all back together
# this changes seqtab and removes/filters those without metadata
proptab <- proptab_all[[1]]
for(tab in 2:length(proptab_all)) {
    if(ncol(proptab_all[[tab]])==0) next
    print(names(proptab_all)[tab])
    proptab <- mergeSequenceTables(proptab,proptab_all[[tab]],repeats="sum")
}

#make sure they're the same size and in the same order as seqtab
all(dim(seqtab)==dim(proptab))
filtrows <- intersect(metadata[rownames(seqtab),"cellid"],rownames(proptab))
filtcols <- intersect(colnames(seqtab),colnames(proptab))
proptab_filt <- proptab[filtrows,filtcols]

# make integer/normalized version
normtab <- floor(1e6*proptab_all[[1]])
for(tab in 2:length(proptab_all)) {
    if(ncol(proptab_all[[tab]])==0) next
    print(names(proptab_all)[tab])
    normtab <- mergeSequenceTables(normtab, floor(1e6*proptab_all[[tab]]),repeats="sum")
}
normtab_filt <- normtab[filtrows,filtcols]


# merge renamed sequence tables
seqtab_filt <- seqtab_none[!is.na(rownames(seqtab_none)),]
seqtab_filt <- mergeSequenceTables(seqtab_filt,seqtab_16S,repeats="sum")
seqtab_filt <- mergeSequenceTables(seqtab_filt,seqtab_18S,repeats="sum")
seqtab_filt <- mergeSequenceTables(seqtab_filt,seqtab_ITS,repeats="sum")
seqtab_filt <- seqtab_filt[filtrows,filtcols]


# save the new tables
saveRDS(seqtab_filt, file.path(outfolder, "seqtab_filt.rds")) # filtered sequence table
saveRDS(proptab_filt, file.path(outfolder,"proptab_filt.rds")) # filtered proportional table 
saveRDS(normtab_filt, file.path(outfolder,"normtab_filt.rds")) # filtered normalized table (1e6)
saveRDS(table_list, file.path(outfolder,"table_list.rds")) # list of which sequence belongs with which subset
