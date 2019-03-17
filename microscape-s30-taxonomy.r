library(dada2); packageVersion("dada2")
library(tidyr)

options(menu.graphics=FALSE)

#	DO TAXONOMIC ASSIGNMENT

# output folder
outfolder <- "out_dada"

# input file
seqtab <- readRDS(file.path(outfolder,"seqtab_final.rds")) #

# Actually assign taxonomies
started <- file.exists(file.path(outfolder,"taxout.rds"))
if(started) {
    taxout <- readRDS(file.path(outfolder,"taxout.rds"))
    bootout <- readRDS(file.path(outfolder,"bootout.rds"))
} else {
    taxout <- list()
    bootout <- list()
}

# get list of reference databases
refdbs <- list.files(".","ref_dada2.*fasta")

# assign taxonomy from each database
for(ref in refdbs) {

    # skip if already done
    if(!is.null(ncol(taxout[[ref]]))) next

    print(ref)
    taxlevels <- NULL
    # which taxa levels to use
    if (grepl("BOLD",ref)) { taxlevels = c("Superkingdom","Phylum","Class","Order","Family","Genus","Species","Accession")
    } else if (ref=="ref_dada2_ITSoneDB_rep_seq_1.131.fasta") { taxlevels = c("Root","Domain","Kingdom","Subkingdom","Phylum","Subphylum","Superclass","Class","Subclass","Infraclass","Superorder","Order","Suborder","Superfamily","Family","Subfamily","Genus","Species","Accession")
    } else if (ref=="ref_dada2_phytoref_cyano.fasta") { taxlevels = c("Domain","Phylum","Class","Order","Family","Genus","Species","Accession")
    } else if (ref=="ref_dada2_phytoref_euks.fasta") { taxlevels = c("Domain","Major_clade","Phylum","Class","Subclass","Order","Suborder","Family","Genus","Species","Accession")
    } else if (ref=="ref_dada2_pr2_version_4.10.0.fasta") {taxlevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species")
    } else if (ref=="ref_dada2_silva_nr_v132_train_set.fasta") { taxlevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    } else if (ref=="ref_dada2_silva_nr_v123_train_set.fasta") { taxlevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    } else if (ref=="ref_dada2_silva_v132.fasta") { taxlevels=c("Root","Domain","Major_clade","Superkingdom","Kingdom","Subkingdom","Infrakingdom","Superphylum","Phylum","Subphylum","Infraphylum","Superclass","Class","Subclass","Infraclass","Superorder","Order","Suborder","Superfamily","Family","Subfamily","Genus","Accession")
    } else { taxlevels = NULL}

    # UNITE is treated specially by DADA2
    if(ref=="ref_dada2_unite.fasta") {
        taxout.tmp <- assignTaxonomy(seqtab, ref, multithread=TRUE, minBoot=0, outputBootstraps=T, verbose=T)
    } else {
        taxout.tmp <- assignTaxonomy(seqtab, ref, multithread=TRUE, minBoot=0, outputBootstraps=T, verbose=T, taxLevels=taxlevels)
    }

    taxout[[ref]] <- taxout.tmp$tax
    bootout[[ref]] <- taxout.tmp$boot

    saveRDS(taxout,file=file.path(outfolder,"taxout.rds"))
    saveRDS(bootout,file=file.path(outfolder,"bootout.rds"))
}


# sets up names lists for shiny app

# cancel this for now since i need edited taxout/ can replace after re-running taxonomy or fixing taxout file
# not sure what ^^^ meant
# need to manually delete these to update
# if(file.exists(file.path(outfolder,"names_list.rds")) & file.exists(file.path(outfolder,"taxa_list.rds"))) {
# 
# 	taxa.list <- readRDS(file.path(outfolder,"taxa_list.rds"))
# 	names.list <- readRDS(file.path(outfolder,"names_list.rds"))
# 
# } else {
	#this sets up the taxa lists for each reference database
	# no need to do this at runtime?
	
	taxa.list <- list() # genus list for quick selection; edit to include all levels
	names.list <- list() # full concatenated name list
		
	for(ref in names(taxout)) {
	    print(ref)

		taxlevels <- NULL
		# which taxa levels to use
		if (grepl("BOLD",ref)) { taxlevels = c("Superkingdom","Phylum","Class","Order","Family","Genus","Species","Accession")
		} else if (ref=="ref_dada2_ITSoneDB_rep_seq_1.131.fasta") { taxlevels = c("Root","Domain","Kingdom","Subkingdom","Phylum","Subphylum","Superclass","Class","Subclass","Infraclass","Superorder","Order","Suborder","Superfamily","Family","Subfamily","Genus","Species","Accession")
		} else if (ref=="ref_dada2_phytoref_cyano.fasta") { taxlevels = c("Domain","Phylum","Class","Order","Family","Genus","Species","Accession")
		} else if (ref=="ref_dada2_phytoref_euks.fasta") { taxlevels = c("Domain","Major_clade","Phylum","Class","Subclass","Order","Suborder","Family","Genus","Species","Accession")
		} else if (ref=="ref_dada2_pr2_version_4.10.0.fasta") {taxlevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species")
		} else if (ref=="ref_dada2_silva_nr_v132_train_set.fasta") { taxlevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
		} else if (ref=="ref_dada2_silva_nr_v123_train_set.fasta") { taxlevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
		} else if (ref=="ref_dada2_silva_nr_v128_train_set.fasta") { taxlevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
		} else if (ref=="ref_dada2_silva_v132.fasta") { taxlevels=c("Root","Domain","Major_clade","Superkingdom","Kingdom","Subkingdom","Infrakingdom","Superphylum","Phylum","Subphylum","Infraphylum","Superclass","Class","Subclass","Infraclass","Superorder","Order","Suborder","Superfamily","Family","Subfamily","Genus","Accession")
		} else if (ref=="ref_dada2_unite.fasta") { taxlevels=c("Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species")
		} else if (ref=="ref_dada2_gg_13_8_train_set_97.fasta") { taxlevels=c("Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species")
		} else if (ref=="ref_dada2_rdp_train_set_16.fasta") { taxlevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
		} else { taxlevels = NULL}
	#    print(taxlevels)
		colnames(taxout[[ref]])[1:length(taxlevels)] <- taxlevels
		colnames(bootout[[ref]])[1:length(taxlevels)] <- taxlevels

		num.taxa <- ncol(seqtab)
		taxchoose <- intersect(taxlevels,c("Superkingdom","Kingdom","Phylum","Order","Family","Genus","Species","Accession"))
		names.list[[ref]] <- as.list(tidyr::unite(data.frame(taxout[[ref]][colnames(seqtab),]), taxchoose, sep=";"))
		esvs <- paste0("ESV_",1:num.taxa)
		names.list[[ref]] <- paste(names.list[[ref]][[1]],esvs,sep=";")
	
#		taxa.list[[ref]] <- unname(taxout[[ref]][colnames(seqtab),"Genus"]) #pick Genus, may want to update this to use taxonomic level selector
		taxa.list[[ref]] <- sort(unique(unlist(as.list(unname(taxout[[ref]])))))

		# replace NAs with lowest non-NA taxonomic name
		for(i in ncol(taxout[[ref]]):1) {
			taxa.list[[ref]][is.na(taxa.list[[ref]])] <- taxout[[ref]][colnames(seqtab)[is.na(taxa.list[[ref]])],i]
		}

	# 
	# 	#phyloseq
	# 	# do this outside and then subset it for plotting
	# 	# convert selected sequences to DNAbin
	# 	# but not yet merged
	# 	## seqtab <- readRDS(paste0(outfolder,"/seqtab_final.rds")) # 
	# 	## taxout <- readRDS(paste0(outfolder,"/taxout.rds"))
	# 	seqs.raw <- sapply(strsplit(colnames(proptab),""), tolower) 
	# 	names(seqs.raw) <- colnames(proptab)
	# 	seqs.bin <- as.DNAbin(seqs.raw)
	# 	# # read in metadata
	# 	# 
	# 	seqs.otu <- otu_table(t(proptab), taxa_are_rows = TRUE)
	# 	seqs.tax <- tax_table(taxout[[ref]][colnames(proptab),])
	# 	seqs.physeq <- phyloseq(seqs.otu, seqs.tax)
	# 	saveRDS(seqs.physeq,file=file.path(outfolder,"seqs.physeq.rds"))

	}
	
 	saveRDS(names.list, file=file.path(outfolder,"names_list.rds"))
 	saveRDS(taxa.list, file=file.path(outfolder,"taxa_list.rds"))
 	saveRDS(taxout, file=file.path(outfolder,"taxout_edit.rds"))
 	saveRDS(bootout, file=file.path(outfolder,"bootout_edit.rds"))
 	
# }







# Plot Classification Quality by taxonomic level

pdf(file=file.path(outfolder,"classification-taxlevel.pdf"),width=12,height=12)
for(taxl in c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {

    bhtaxlevels <- sapply(taxout, function(x) which(colnames(x)==taxl))

    plot(NULL,NULL,xlim=c(0,100),ylim=c(0,ncol(seqtab)),ylab="# ESVs classified",xlab="quality cutoff",main=taxl)

    for(i in 1:length(taxout)) {
        ref <- names(taxout)[i]
        qualsave <- NULL
        for(qual in 1:100) {
            bhtax <- taxout[[ref]]
            bhboot <- bootout[[ref]]
        
            coln <- unlist(bhtaxlevels[i])
            besthits <- cbind(as.data.frame(unclass(rle(sort(unname(bhtax[unname(bhboot[,coln]>qual),,drop=F])[,coln,drop=F])))))
            besthits <- besthits[rev(order(besthits[,1])),]
            qualsave <- c(qualsave,sum(besthits[,1]))
        }
        print(head(besthits,5))
        lines(qualsave,col=i)
        text(80,15000-500*i,ref,col=i)
    }
}
dev.off()



# Plot Classification Quality by database

pdf(file=file.path(outfolder,"classification-reference.pdf"),width=12,height=12)

for(i in 1:length(taxout)) {

    ref <- names(taxout)[i]
    bhtax <- taxout[[ref]]
    bhboot <- bootout[[ref]]

    bhtaxlevels <- colnames(bhtax)

    plot(NULL,NULL,xlim=c(0,100),ylim=c(0,ncol(seqtab)),ylab="# ESVs classified",xlab="quality cutoff",main=ref)

    for(coln in 1:length(bhtaxlevels)) {

        qualsave <- NULL
        for(qual in 1:100) {
            besthits <- cbind(as.data.frame(unclass(rle(sort(unname(bhtax[unname(bhboot[,coln]>qual),,drop=F])[,coln,drop=F])))))
            besthits <- besthits[rev(order(besthits[,1])),]
            qualsave <- c(qualsave,sum(besthits[,1]))
        }
        print(head(besthits,5))
        lines(qualsave,col=coln)
        text(80,15000-500*coln,colnames(bhtax)[coln],col=coln)
    }
}
dev.off()

