##### 
##### This program runs DADA2 on input directories of the form 20*/ which have demultiplexed files named as PLATE_WELL_READ.ctrimmed.fastq.gz
##### demultiplex using run-plate-prep.sh
##### trim primers using ./run-bbduk-remove-primers3.sh

# install requirements
# install.packages("ShortRead","ggplot2","ape","cluster")
# source("https://bioconductor.org/biocLite.R")
# biocLite("devtools")
# library("devtools")
# devtools::install_github("benjjneb/dada2")

library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")

options(menu.graphics=FALSE)

# if called with a directory name it will use only that directory, otherwise all subdirectories with matching file names
args = commandArgs(trailingOnly=TRUE)

if(length(args)>0) {
    fastqs <- list.files(args,recursive=T,pattern="*.ctrimmed.fastq.gz$")    
} else {
    fastqs <- list.files(".",recursive=T,pattern="*.ctrimmed.fastq.gz$")
}

# output.path
outfolder <- "out_dada"
# create directory
if(!file_test("-d", outfolder)) dir.create(outfolder, recursive=T)

fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fastqs <- fastqs[!grepl("Unassigned",fastqs)] #skip unassigned reads

fnF <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnR <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files

fnbig <- file.size(fnF)>20 & file.size(fnR)>20
fnF <- fnF[fnbig] #remove empty files
fnR <- fnR[fnbig] #remove empty files

# Get sample names from the first part of the forward read filenames
sample.names <- unname(sapply(fnF,function(x) sub("_R1[.].+$","",x,perl=T)))

fnlist <- c(fnF,fnR)

# separate into plates for error calling
# using plate instead of flow cell because PCR history can affect error rates 
plate.list <- unique(sapply(strsplit(fnF, "_"), `[`,1))

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(".", outfolder)
if(!file_test("-d", filt_path)) dir.create(filt_path)


for(plate in plate.list) {
    dada.save <- list()

    print(plate)
    platerun <- strsplit(plate,"/")[[1]][1]
    if( !file_test("-d", paste0(filt_path,"/",platerun))) dir.create(paste0(filt_path,"/",platerun))

	plate.samples <- grep(plate,sample.names,value=T)
	plate.paths <- grep(plate,fnlist,value=T)
	plate.f <- grep("_R1",plate.paths,value=T)
	plate.r <- grep("_R2",plate.paths,value=T)

	filtFs <- file.path(filt_path, platerun, paste0(basename(plate.samples), "_R1.filt.fastq.gz"))
	filtRs <- file.path(filt_path, platerun, paste0(basename(plate.samples), "_R2.filt.fastq.gz"))

	fto <- filterAndTrim(fwd=plate.f, filt=filtFs,
              rev=plate.r, filt.rev=filtRs,
              maxEE=2, truncQ=11, maxN=0, rm.phix=T,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

    print(fto)
	print(file.size(filtFs))
	
	#remove empty files
	plate.samples <- plate.samples[!is.na(file.size(filtFs))]
	filtFs <- filtFs[!is.na(file.size(filtFs))]
	filtRs <- filtRs[!is.na(file.size(filtRs))]
	

	#calculate error rates
	#First, it is not advised to pool samples that don’t share an “error history”, 
	#in particular samples that come from different sequencing runs or different PCR protocols.
	#Second, it is generally not necessary to estimate error rates across an entire sequencing run. 
	#Because error rate estimation requires multiple loops through the dada(...) algorithm, 
	#it increases compute time significantly. Therefore it is often desirable to estimate error rates 
	#on a subset of the samples, and then use those error rates to process all of the samples with selfConsist=FALSE.
	
    # Learn forward error rates
    errF <- learnErrors(filtFs, multithread=TRUE)
    # Learn reverse error rates
    errR <- learnErrors(filtRs, multithread=TRUE)

    mergers <- vector("list", length(plate.samples))
    names(mergers) <- plate.samples
    for(sam in 1:length(plate.samples)) {
      cat("Processing:", plate.samples[sam], "\n")
        derepF <- derepFastq(filtFs[[sam]], verbose=T)
        ddF <- dada(derepF, err=errF, multithread=TRUE)
        derepR <- derepFastq(filtRs[[sam]], verbose=T)
        ddR <- dada(derepR, err=errR, multithread=TRUE)
        #using trimOverhang=TRUE instead of other methods to clip ends 
        merger <- mergePairs(ddF, derepF, ddR, derepR, trimOverhang=TRUE, verbose=T, minOverlap=10)
        mergers[[plate.samples[sam]]] <- merger

		dada.save[[plate.samples[sam]]] <- list("dadaFs"=ddF,"dadaRs"=ddR,"derepFs"=derepF,"derepRs"=derepR)
    }

	seqtab <- makeSequenceTable(mergers)
    write.table(seqtab, paste0(outfolder,"/",plate,".tsv"), sep="\t", quote=F)
	saveRDS(seqtab, paste0(outfolder,"/",plate,".seqtab.rds"))
    saveRDS(dada.save,paste0(outfolder,"/",plate,".dada.save.rds"))

}

