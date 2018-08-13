library(dada2); packageVersion("dada2")

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


