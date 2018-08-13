library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")
library(ape)
library(cluster)

options(menu.graphics=FALSE)

####### merge taxa by phylogeny


library(msa)
library(ips) #for trimEnds()
library(ape) #for write.dna
library(txtplot)
library(DECIPHER) #for AlignSeqs
library(digest) #for digest

outfolder <- "out_dada"

# read from disk
seqtab <- readRDS(file.path(outfolder,"seqtab_final.rds")) # 
taxout <- readRDS(file.path(outfolder,"taxout.rds"))
bootout <- readRDS(file.path(outfolder,"bootout.rds"))

# assign primer sets
samples_ITS <- which(grepl("ITS",rownames(seqtab),ignore.case=T))
samples_16S <- which(grepl("16S",rownames(seqtab),ignore.case=T))
samples_18S <- which(grepl("18S",rownames(seqtab),ignore.case=T))

# for each taxonomy, make alignment and phylogenetic tree for each taxonomic group
for(thistax in names(taxout)) {
	
	nicetax <- strsplit(thistax,".",fixed=T)[[1]][1]
	dir.create(file.path(outfolder,"phylo",nicetax), recursive=T)

	print(paste0("################# ",nicetax," #################"))

	#remove unassigned samples and any sequences found only in them

	mytax <- taxout[[thistax]]
	myboot <- bootout[[thistax]]
	myseqtab <- seqtab
	
	if(thistax=="unite") {
		myseqtab <- myseqtab[samples_ITS,]
	} else if (thistax=="pr2") {
		myseqtab <- myseqtab[samples_18S,]
	} else if (any(grepl("ITS",thistax,ignore.case=T))) {
		myseqtab <- myseqtab[samples_ITS,]
	} else if (thistax=="phytoeuk") {
		myseqtab <- myseqtab[samples_18S,]
	} else if (thistax=="phytocyano") {
		myseqtab <- myseqtab[samples_16S,]
	} else if (any(grepl("18S",thistax,ignore.case=T))) {
		myseqtab <- myseqtab[samples_18S,]
	} else if (thistax=="silva_train") {
		myseqtab <- myseqtab[c(samples_16S,samples_18S),]	
	}
		
	#remove sequences that have no sequences after cutting samples
	toexclude <- colSums(myseqtab)>0
	mytax <- mytax[toexclude,]
	myboot <- myboot[toexclude,]
	myseqtab <- myseqtab[,toexclude]
	
	mytax <- cbind(mytax,"ESV"=paste0("esv",sprintf("%05.f", 1:ncol(myseqtab)), "_", colSums(myseqtab)) ) #add otu column with seq counts

	startingsize=sum(myseqtab)

	#remove NAs, replace with highest best taxa
	for(taxlevel in 2:ncol(mytax)) {
		gentax <- mytax[,taxlevel]
		nas <- unname(which(is.na(gentax)))
		gentax[nas] <- paste0(mytax[nas,taxlevel-1],"_",colnames(mytax)[taxlevel])
		mytax[,taxlevel] <- gentax
	}

	pdf(file=file.path(outfolder,"phylo",nicetax,"boots.pdf"),width=12,height=12)
	for(i in 1:ncol(myboot)) {
		if(i==1) plot(sort(myboot[,i]),ylim=c(0,100),col=i,type='l',xlab="number of sequences",ylab="bootstrap support")
	else {lines(sort(myboot[,i]),col=i)}
	}
	dev.off()

	#multiple sequence alignment to merge ESVs

	# find clusters of identical sequences after trimming / may be fixed after new cutadapt trimming? 
	# may still need to find non-rRNA sequences

	gen.save <- list()
	all.bad <- list()

	for(taxlevel in (ncol(mytax)-2):2) {
		gentax <- mytax[,taxlevel]
		uniqtax <- unique(gentax)

	###

	# how bad are the bootstraps?
	#    plot(0,xlim=c(0,100),ylim=c(0,100))
	
		pdf(file=file.path(outfolder,"phylo",nicetax,paste0(taxlevel,"-distmap.pdf")),width=24,height=24)
		for(gen in uniqtax) {
#			genseqs <- names(which(mytax[,taxlevel]==gen))
			genseqs <- rownames(mytax)[which(mytax[,taxlevel]==gen)]
			nicegen <- gsub("[/\\]","",gen)
			
			if(length(genseqs)<20) {next}
#			if(length(genseqs)>1000) {next}
			print(paste("Processing:", taxlevel,gen,length(genseqs),sep=" "))

	# how bad are the bootstraps?
	#        lines(x=(1:length(genseqs))*(100/length(genseqs)),y=sort(tax$boot[genseqs,6]),ylim=c(0,100))

			dnaseqs <- DNAStringSet(genseqs)
	#        names(dnaseqs) <- md5(mytax[genseqs,"OTU"])
			names(dnaseqs) <- sapply(dnaseqs,digest)
		
			#sometimes MUSCLE fails... might be because of truncated sample names? try just putting in otunums
			genmsa <- NULL
			try(genmsa <- msa(dnaseqs,type="dna","Muscle",order="input"))
			if(is.null(genmsa)) {
				genmsa <- AlignSeqs(dnaseqs)
				names(genmsa) <- mytax[genseqs,"ESV"]
				subseqs <- genmsa			
			} else {
				names(genmsa@unmasked) <- mytax[genseqs,"ESV"]
				subseqs <- genmsa@unmasked
			}

			#BrowseSeqs(genmsa)
			write.dna(x=as.character(genmsa),file=file.path(outfolder,"phylo",nicetax,paste0(taxlevel,"-",nicegen,".fasta")),format="fasta",nbcol=-1,colsep='', colw=1000, indent=0)
		#    print(genmsa,show="complete",halfNrow=-1)
	#        system(paste0("/Applications/ScienceApps/seaview.app/Contents/MacOS/seaview -printout -blocksize 1000 -landscape -svg 1200 -o ", outfolder,"/",taxlevel,"-",gen,".svg"," ", outfolder,"/",taxlevel,"-",gen,".fasta"))

			try({

				#observed problems:
					#1 big gaps at beginning (common) or end (less common)
					#2 tail on end (commonly primer) or beginning (commonly short)

	# maybe just filter out any sequences that have Ns after it's all over
	# 			#1: filter by missing start or end (+5 bases)
	# 			started <- apply(as.matrix(genmsa),2,function(x) sum(x!="-")/length(x))>0.6
	# 			firstcol <- which.max(started)
	# 			lastcol <- length(started) - which.max(rev(started))
	# 			msastart <- as.matrix(genmsa)[,1:(firstcol+5)]
	# 			msaend <- as.matrix(genmsa)[,(lastcol-5):length(started)]
	# 			whichbad1 <- which(apply(msastart,1,function(x) all(x=="-")) | apply(msaend,1,function(x) all(x=="-")))
	# 						
	# 			#2: find the tails: should be longer and have more at beginning/end
	# 			tails <- apply(msaend,1,function(x) sum(x!="-"))
	# 			longones <- nchar(dnaseqs)
	# 			
	#             # first filter for sequences that have more than median # of [>4 indels]
	#             genindel <- dist.dna(as.DNAbin(genmsa),model="indel")
	# 			gendist <- dist.dna(as.DNAbin(genmsa),model="JC69",pairwise=T)
	# 			txtplot(x=rowSums(as.matrix(genindel)),y=rowSums(as.matrix(gendist)))
	# 			
	#             indsums <- rowSums(as.matrix(genindel)>2)
	#             whichbad <- which(indsums>median(indsums))
	# 			all.bad[[gen]] <- genseqs[whichbad]
	# 
	# 			txtplot(sort(indsums))
	# 
	# #optimization of max indels
	#                savei <- NULL
	#                for(i in 1:10) {
	#             indsums <- rowSums(as.matrix(genindel)>i)
	#             whichbad <- which(indsums>median(indsums))
	#             txtplot(sort(indsums))
	#                savei[i] <- length(whichbad)
	#                }
	#                txtplot(savei)
	# 
	# 			#rownames(as.character(gentrim))
	# 			#how many are lost?
	# 			slost <- sum(unname(colSums(myseqtab[,which(unname(mytax[,7]) %in% as.character(names(whichbad))),drop=F])))
	# 			#what fraction?
	# 			stotal <- sum(unname(colSums(myseqtab[,which(unname(mytax[,7]) %in% as.character(names(indsums))),drop=F])))
	# 
	# 			print(paste0("removing ", length(whichbad), " of ",ncol(myseqtab)," ESVs (", round(100*length(whichbad)/ncol(myseqtab),digits=3)," %) while indel filtering ",gen))
	# 			print(paste0("removing ", slost, " of ",stotal," sequences (", round(100*slost/stotal,digits=3)," %) while indel filtering ",gen))
	# 
	# 			#update tables to remove bad sequences
	# 			if(length(whichbad)>0) {
	# 				myseqtab <- myseqtab[,-which(colnames(myseqtab) %in% genseqs[whichbad])]
	# 				mytax <- mytax[-which(rownames(mytax) %in% genseqs[whichbad]), ]
	# 				myboot <- myboot[-which(rownames(myboot) %in% genseqs[whichbad]), ]
	# 	            subseqs <- genmsa@unmasked[-whichbad]
	# 			} else {
	# 				subseqs <- genmsa@unmasked
	# 			}
	# 						

				#  trim ends if majority is shorter
				minn <- round(0.6*length(subseqs))
				if(minn<4) next
				gentrim <- ips::trimEnds(as.DNAbin(as.matrix(subseqs)), min.n.seq=minn)
				write.dna(x=toupper(as.character(gentrim)),file=paste0(outfolder,"/","phylo","/",nicetax,"/",taxlevel,"-",nicegen,".trim.fasta"),format="fasta",nbcol=-1,colsep='', colw=1000, indent=0)
	#            system(paste0("/Applications/ScienceApps/seaview.app/Contents/MacOS/seaview -printout -blocksize 1000 -landscape -svg 1200 -o ", outfolder,"/",taxlevel,"-",gen,".trim.svg"," ", outfolder,"/",taxlevel,"-",gen,".trim.fasta"))

	#dim(as.character(gentrim))[2]
	#length(dnaseqs)

				oldlength <- sum(nchar(genseqs))
				newlength <- sum(as.character(gentrim)!="-" & as.character(gentrim)!="n")
	#			newlength <- sum(unlist(lapply(as.list(gentrim),function(x) sum(as.character(x)!="-"))))

	#			print(paste0(round(100*newlength/oldlength,3)," % of original sequence bp retained (",newlength," from ",oldlength,") in ",gen))
				print(paste0("removing ", oldlength-newlength, " of ",oldlength," basepairs (", round(100*(1-newlength/oldlength),3)," %) after trimming ",gen))
	# [1] "removing 48 of 98455 basepairs (0.049 %) after trimming Mitochondria_Genus"
	# [1] NaN NaN

				if(newlength/oldlength < 0.9) print("WARNING LOW RETENTION")

				# calculate distances and print
				gendist <- dist.dna(gentrim, as.matrix=T, model="raw", pairwise.deletion=T) + dist.dna(gentrim, as.matrix=T, model="indel")
	#            print(range(gendist))
				print(range(gendist))
	 #           gendist[is.nan(gendist)] <- NA
				diag(gendist) <- NA
				if(length(genseqs)<400) {
					heatmap(gendist,scale="none",main=gen,na.rm=T)
				}
				gen.save[[colnames(mytax)[taxlevel]]][[gen]] <- gendist

				#cut sequences with Ns at highest taxonomic level
				#update myseqtab with trimmed sequences
				newseqs <- apply(as.matrix(as.character(gentrim)),1, function(x) toupper(paste0(x, collapse = '')))
				names(newseqs) <- genseqs #name them by old seqs

# 
# 		#THIS IS A BAD IDEA FOR ITS
# 				if(taxlevel==2) {
# 					ns <- which(grepl(x=newseqs,pattern="N"))
# 					newseqs <- gsub("-|N","",newseqs)
# 					if(length(ns)) {
# 						rmseqs <- newseqs[ns]
# 						print(paste0("removing ",ns," ESVs with Ns (",sum(myseqtab[,names(rmseqs)])," counts)"))
# 						myseqtab <- myseqtab[,-which(colnames(myseqtab) %in% names(rmseqs))]
# # [1] "removing 4 of 42343 basepairs (0.009 %) after trimming p__Mucoromycota"
# # [1]   0.0000 280.6242
# # Error in myseqtab[, rmseqs] : subscript out of bounds
# 					}
# 				} else {
# 					newseqs <- gsub("-|N","",newseqs)
# 				}
# 			

				newseqs <- gsub("-|N","",newseqs)

				##sh*tty hack to replace names
				tmpcol <- colnames(myseqtab)
				for(i in which(colnames(myseqtab) %in% names(newseqs))) {
	#				print(i)
					tmpcol[i] <- newseqs[tmpcol[i]]
				}
				colnames(myseqtab) <- tmpcol
				rownames(mytax) <- tmpcol
				rownames(myboot) <- tmpcol


				#duplicated only selects the second, not the first match
				replicated <- function(x) {x %in% unique(x[ duplicated(x)])}

				#find by duplicated rather than by genetic distance... 
				#this results in a lot less merging... the dist.dna function
				#removes gaps before calculating distance, so might have been merging seqs with indels
				
#				dups <- sum(gendist==0,na.rm=T)/2
				dups <- sum(replicated(colnames(myseqtab)))
				print(paste0("merging ",dups, " duplicates out of ",nrow(gendist)," ESVs"))

				#merge rows
				#might be faster if done gendist rowwise rather than pairwise?
				if(dups>0) {
#					zeros <- which(gendist==0,arr.ind=T)
					zeros <- which(replicated(colnames(myseqtab)))
#					zeroseqs <- unname(newseqs[which(unique(rownames(zeros)) %in% rownames(gentrim))])
					zeroseqs <- colnames(myseqtab)[zeros]
	#				zeroseqs <- as.character(gentrim[rownames(zeros)])
					#cut out matching columns, add them back one by one
#					newseqtab <- myseqtab[,-which(colnames(myseqtab) %in% zeroseqs)]
#					tmpseqtab <- myseqtab[,which(colnames(myseqtab) %in% zeroseqs)]
#					addseqtab <- tmpseqtab[,1,drop=F]
					newseqtab <- myseqtab[,-zeros]
					tmpseqtab <- myseqtab[,zeros]
					addseqtab <- tmpseqtab[,1,drop=F]
					if(length(zeros)>1) {
						for(i in 2:ncol(tmpseqtab)) {
							addseqtab <- mergeSequenceTables(addseqtab,tmpseqtab[,i,drop=F],repeats="sum")
						}
					}
				
					myseqtab <- mergeSequenceTables(newseqtab,addseqtab,repeats="sum")
	#				browser()
					#Q: why are any sequences being lost here?
					print(paste0(startingsize-sum(myseqtab)," hits removed so far from seqtab (",round(100-100*sum(myseqtab)/startingsize,digits=2),"%)"))

					mytax <- mytax[colnames(myseqtab),]
					myboot <- myboot[colnames(myseqtab),]
					mytax[,7] <- paste0("esv",sprintf("%05.f", 1:ncol(myseqtab)), "_", colSums(myseqtab))

				}

			})
#	         browser()
			}
dev.off()
	}


roworder <- rev(order(rowSums(myseqtab)))
#roworder <- roworder[rowSums(seqtab[roworder,])>=100]

# sort and cut out any sequences that have no reads
#colorder <- rev(order(colSums(seqtab)))
colorder <- 1:ncol(myseqtab)
colorder <- colorder[which(colSums(myseqtab[roworder,])>0)]
seqmat <- t(myseqtab[roworder,colorder])

seqmeta <- cbind("esv"=paste0("esv",sprintf("%05.f", 1:ncol(myseqtab))),"sequence"=colnames(myseqtab),mytax,myboot)
seqexp <- cbind(seqmeta[colorder,],"sum"=rowSums(seqmat),seqmat)
write.table(seqexp,file=paste0(outfolder,"/","phylo","/",nicetax,"/myseqexp.tsv"),sep="\t",quote=F,row.names=F)
saveRDS(list(myseqtab,mytax,myboot),file=paste0(outfolder,"/","phylo","/",nicetax,"/processed.rds"))


}


# browser()





#count losses at each step:
# demultiplex
# adapter trimming
# primer trimming
# quality trimming
# sequence alignment















