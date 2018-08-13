library(dada2); packageVersion("dada2")
options(menu.graphics=FALSE)

#	MERGE SEQUENCING RUNS AND BASIC QUALITY CONTROL

#output folder
outfolder <- "out_data"

# Merge all runs
rdsfiles <- list.files(".","*seqtab.rds",recursive=T)

seqtab.all <- readRDS(rdsfiles[1])
print(sum(seqtab.all))
for(i in 2:length(rdsfiles)) {
	print(rdsfiles[i])
	seqtab.tmp <- readRDS(rdsfiles[i])
	print(sum(seqtab.tmp))
	if(sum(seqtab.tmp)>0) seqtab.all <- mergeSequenceTables(seqtab.tmp,seqtab.all)
}

dim(seqtab.all)

# Remove chimeras
seqtab <- removeBimeraDenovo(seqtab.all, method="consensus", multithread=TRUE, verbose=T)
sum(seqtab)/sum(seqtab.all)

# remove short sequences (<50 can't be assigned taxonomy)
# although this is maybe a bad idea ---- some ITS are actually <50bp...
seqtab.rm <- sapply(colnames(seqtab),nchar)<50
seqtab <- seqtab[,!seqtab.rm]
dim(seqtab)

#how many are retained?
sum(seqtab)/sum(seqtab.all)

par(mfrow=c(2,2))
plot(log10(sort(colSums(seqtab))), main="Reads per ESV")
plot(log10(sort(rowSums(seqtab))), main="Reads per Sample")

plot(log10(sort(colSums(seqtab>0))), main="Samples per ESV")
plot(log10(sort(rowSums(seqtab>0))), main="ESVs per Sample")

dev.off()

# remove orphans (ESVs present in only 1 sample)
seqtab.or <- colSums(seqtab>0) < 2
seqtab.orphans <- seqtab[,seqtab.or]
seqtab <- seqtab[,!seqtab.or]
dim(seqtab)
sum(seqtab)/sum(seqtab.all)

# remove singletons and doubletons
seqtab.sing <- colSums(seqtab) < 3
seqtab <- seqtab[,!seqtab.sing]
dim(seqtab)
sum(seqtab)/sum(seqtab.all)

# remove samples with less than 100 sequences
seqtab.sm <- rowSums(seqtab) < 101
seqtab.small <- seqtab[seqtab.sm,]
seqtab <- seqtab[!seqtab.sm,]
dim(seqtab)
sum(seqtab)/sum(seqtab.all)

# Write to disk

# final output:
# samples with more than 100 sequences each
# sequences present in more than 1 sample and with more than 2 sequences total
saveRDS(seqtab, file.path(outfolder,"seqtab_final.rds")) #
saveRDS(seqtab.orphans, file.path(outfolder,"seqtab_orphans.rds"))
saveRDS(seqtab.small, file.path(outfolder,"seqtab_small.rds"))

