############### Taxonomic Classification
#### not quite automated yet (3/14/2019)


install.packages(c("data.table","taxize","taxizedb","dplyr","parallel","taxonomizr","lme4"))
library(data.table)
library(taxize) #to reshape taxonomy files
library("taxizedb")
library("dplyr")
library(parallel)

# download and format taxonomic databases

taxlevels <- c("root","domain","major_clade","superkingdom","kingdom","subkingdom","infrakingdom","superphylum","phylum","subphylum","infraphylum","superclass","class","subclass","infraclass","superorder","order","suborder","superfamily","family","subfamily","genus","species")

x <- db_download_ncbi()
db_load_ncbi()
src_ncbi(x)


#SILVA
system("wget https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz; mv silva_nr_v132_train_set.fa.gz ref_dada2_silva_nr_v132_train_set.fa.gz")
system("wget https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz")
#10Gb system("wget https://www.mothur.org/w/images/3/32/Silva.nr_v132.tgz")
system("ln -s /scratch/cryomics/reference_dbs/SILVA/silva.full_v132.full.tax")
#AC201870.A6BReg23       root;Bacteria;Bacteria_mc;Bacteria_pk;Bacteria_ki;Bacteria_bk;Bacteria_ik;Bacteria_pp;Proteobacteria;Proteobacteria_bp;Proteobacteria_ip;Proteobacteria_pc;Gammaproteobacteria;Gammaproteo
system("ln -s /scratch/cryomics/reference_dbs/SILVA/silva.nr_v132.align")
# >AC201870.A6BReg23      94.48   Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Candidatus Regiella;
# ..................................................................................................................................................................................................................
#get NR into dada format
#paste -d '\n' <(awk '{print ">"$2 $1}' silva.full_v132.full.tax) <(grep -v '>' silva.nr_v132.align | tr -d '.-') > ref_dada2_silva_v132.fasta
#taxLevels = c("Root","Domain","Major_clade","Superkingdom","Kingdom","Subkingdom","Infrakingdom","Superphylum","Phylum","Subphylum","Infraphylum","Superclass","Class","Subclass","Infraclass","Superorder","Order","Suborder","Superfamily","Family","Subfamily","Genus")

#OLD SILVA via DADA2
system("wget https://zenodo.org/record/824551/files/silva_nr_v128_train_set.fa.gz")
system("gunzip silva_nr_v128_train_set.fa.gz")
system("mv silva_nr_v128_train_set.fa ref_dada2_silva_nr_v128_train_set.fasta")

#OLDER SILVA via DADA2
system("wget https://zenodo.org/record/158958/files/silva_nr_v123_train_set.fa.gz")
system("gunzip silva_nr_v123_train_set.fa.gz")
system("mv silva_nr_v123_train_set.fa ref_dada2_silva_nr_v123_train_set.fasta")

#RDP training set 16 via DADA2
system("wget https://zenodo.org/record/801828/files/rdp_train_set_16.fa.gz")
system("gunzip rdp_train_set_16.fa.gz")
system("mv rdp_train_set_16.fa ref_dada2_rdp_train_set_16.fasta")

#GG training set 13.8 via DADA2
system("wget https://zenodo.org/record/158955/files/gg_13_8_train_set_97.fa.gz")
system("gunzip gg_13_8_train_set_97.fa.gz")
system("mv gg_13_8_train_set_97.fa ref_dada2_gg_13_8_train_set_97.fasta")

#PR
system("wget https://github.com/vaulot/pr2_database/releases/download/4.10.0/pr2_version_4.10.0_dada2.fasta.gz")
system("gunzip pr2_version_4.10.0_dada2.fasta.gz")
system("mv pr2_version_4.10.0_dada2.fasta ref_dada2_pr2_version_4.10.0.fasta")
#>Eukaryota;Alveolata;Dinoflagellata;Dinophyceae;Peridiniales;Kryptoperidiniaceae;Unruhdinium;Unruhdinium_kevei;
#taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species")



#UNITE
system("wget -O unite.fasta.zip https://files.plutof.ut.ee/doi/B2/07/B2079372C79891519EF815160D4467BBF4AF1288A23E135E666BABF2C5779767.zip")
system("wget -O unite.fasta.zip https://files.plutof.ut.ee/doi/C8/E4/C8E4A8E6A7C4C00EACE3499C51E550744A259A98F8FE25993B1C7B9E7D2170B2.zip")
system("unzip *.zip")
system("cat sh_general_release_dynamic_01.12.2017.fasta sh_general_release_dynamic_s_01.12.2017.fasta > ref_dada2_unite.fasta")
#>Agaricus_chiangmaiensis|JF514531|SH174817.07FU|reps|k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__Agaricaceae;g__Agaricus;s__Agaricus_chiangmaiensis



#PhytoREF plastid sequences
system("wget http://phytoref.sb-roscoff.fr/static/downloads/PhytoRef_with_taxonomy.fasta")
#>356#|Eukaryota|Alveolata|Dinophyta|Dinophyceae|Dinophyceae_X|Suessiales|Suessiales_X|Suessiaceae|Symbiodinium|Symbiodinium cladeA
#>Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Coscinodiscophyceae;Thalassiosirales;Thalassiosirales;Thalassiosiraceae;Prorosira;Prorosira pseudodelicatula;FJ002167;
system("wget http://phytoref.sb-roscoff.fr/static/downloads/Cyanobacteria_p98.fasta")
#>AAVU01000008.771.2264|Bacteria|Cyanobacteria|Cyanobacteria|SubsectionIII|FamilyI|Lyngbya|Lyngbya|Lyngbya+sp.
#>Bacteria;Cyanobacteria;Cyanobacteria;SubsectionIII;FamilyI;Phormidium;Phormidium;Phormidium ambiguum;AB003167.1.1439;
system("wget http://phytoref.sb-roscoff.fr/static/downloads/PhytoRef_tb_sequence.tab") #table
system("wget -O phytoref.zip https://ndownloader.figshare.com/articles/4689826/versions/1; unzip phytoref.zip") #mothur formatted
#phytoref_mothur_version_1.0.fasta
#>202#HE610155
#phytoref_mothur_version_1.0.taxo
#202#HE610155    Eukaryota;Archaeplastida;Chlorophyta;Ulvophyceae;Ulvophyceae_X;Ulvophyceae_XX;Ulvophyceae_XXX;Ulvophyceae_XXXX;Desmochloris;Desmochloris_halophila;

euhead <- paste0(">",system("grep '>' PhytoRef_with_taxonomy.fasta | cut -f2 -d'#' | awk -F'|' 'BEGIN{OFS=\";\";} {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$1\";\"}' | tr '+' ' ' | sed -E 's/_X+;/;/g;s/_X+ / /g' ",intern=T))
euseqs <- system("grep -v '>' PhytoRef_with_taxonomy.fasta",intern=T)
write.table(file="ref_dada2_phytoref_euks.fasta",cbind(euhead,euseqs),sep="\n",quote=F,row.names=F,col.names=F)
#taxLevels = c("Domain","Major_clade","Phylum","Class","Subclass","Order","Suborder","Family","Genus","Species","Accession")

cyhead <- paste0(">",system("grep '>' Cyanobacteria_p98.fasta | cut -b2- | awk -F'|' 'BEGIN{OFS=\";\";} {print $2,$3,$4,$5,$6,$7,$8,$9,$1\";\"}' | tr '+' ' '",intern=T))
cyseqs <- system("grep -v '>' Cyanobacteria_p98.fasta",intern=T)
write.table(file="ref_dada2_phytoref_cyano.fasta",cbind(cyhead,cyseqs),sep="\n",quote=F,row.names=F,col.names=F)
#taxLevels = c("Domain","Phylum","Class","Order","Family","Genus","Species","Accession")

#giving up for now, using their taxonomy
#Prr_id  Original_id                                                                                    Longueur           P1       P2       type  s_Createur  s_Creation_date  s_Updateur      s_Update_Date  s_Accession
#1022    Prasinophyceae_cladeVIIA1_Genusnov_spnov_RCC1032                                               759                1        759      U     JD          13-12-17-10-12
#2417    NC_013088.94625.96116                                                                          1492               94625    96116    G     RC          14-01-07-16-22   NC_013088
# 

prin <- read.table("PhytoRef_tb_sequence.tab",sep="\t",head=T,row.names=1,strings=F)
taxid1 <- accessionToTaxa(paste0(prin$s_Accession,".1"),"accessionTaxa.sql")
taxid1[is.na(taxid1)] <- -1
taxy1 <- getTaxonomy(taxid1,taxaNodes,taxaNames,mc.cores=12,desiredTaxa=taxlevels)

intax <- system('grep ">" PhytoRef_with_taxonomy.fasta | cut -f1,11 -d"|" | cut -b2- | tr "#" "|" | sort -n',intern=T)
intab <- read.table(text = intax, sep = "|", colClasses = "character",comment="",row.names=1)
intab[,2] <- gsub(" sp.","",intab[,2])
intab[,2] <- gsub("_X+$","",intab[,2])
intab[is.na(intab)] <- -1

taxid2 <- getId(intab[,2],taxaNames) #ambiguous
taxy2 <- getTaxonomy(taxid2,taxaNodes,taxaNames,mc.cores=12,desiredTaxa=c("genus","species"))

less(cbind(intab,taxid2,taxy2))

allin <- merge(prin,intab,by="row.names",all=TRUE,sort=F)
allin <- allin[,-c(4:11)]

prtaxa <- getTaxonomy(taxaIds,taxaNodes,taxaNames, mc.cores=12, desiredTaxa=c("root","domain","major_clade","superkingdom","kingdom","subkingdom","infrakingdom","superphylum","phylum","subphylum","infraphylum","superclass","class","subclass","infraclass","superorder","order","suborder","superfamily","family","subfamily","genus"))
prtaxa[,2] <- prtaxa[,4]
prtaxa[,4] <- NA


infilltaxa <- function(inmat) {
    taxlevels <- c("root","domain","major_clade","superkingdom","kingdom","subkingdom","infrakingdom","superphylum","phylum","subphylum","infraphylum","superclass","class","subclass","infraclass","superorder","order","suborder","superfamily","family","subfamily","genus","species")
    output <- taxlevels %in% colnames(inmat)
    if(output[1]) { inmat[,"root"] <- "root" } else { inmat <- cbind("root"="root",inmat) }
    output[1] <- TRUE
    taxabb <- c("ro","do","mc","pk","ki","bk","ik","pp","ph","bp","ip","pc","cl","bc","ic","po","or","bo","pf","fa","bf","ge","sp")[output]
    cols <- ncol(inmat)
    tax.mat <- as.matrix(inmat)
    tmptax <- NULL
#    outlevels <- c("domain","phylum","class","order","family","genus")
    tax.mat.na <- is.na(tax.mat)
    for(i in 1:nrow(tax.mat)) {
        if (i%%1000==0) {print(i)}
        #this fills in the empty gaps by using the closest higher taxonomic level appended with an abbreviation for the current taxonomic level
        for(j in 1:cols) {
            ifelse(tax.mat.na[i,j], {tax.mat[i,j] <- paste0(tmptax,"_",taxabb[j])}, {tmptax <- tax.mat[i,j]} )
        }
    }
return(tax.mat)
}


#ITSoneDB
#have to click at http://itsonedb.cloud.ba.infn.it/index.jsp
#ITS1 located by both ENA and HMM
#Representative ITS1 Sequences with Flanking regions

# ITSoneDB_rep_seq_and_flanking_1.131.fasta
# >EU070644_ITS1_HMM|Eryngium fernandezianum|477861|ITS1 located by HMM annotation, 371bp
# tcgatgcctgcaaagcagaacgacccgcgaacacgtcaaaaataacgggcgagcggtccggggggcgcaagctccacgcgtccgcgaacccgcaggtcgagggcgtccctgggcgctcgacggccgcaaactcaccccggcgcggaatgcgccaaggaaatagaaccggactgaacgttctcgcccccgttcgcgggtggcgatggcgtc

#contains 18S sequence
#>root;Eukaryota;Fungi;Dikarya;Ascomycota;Ascomycota_bp;Ascomycota_pc;Ascomycota_cl;Ascomycota_bc;Ascomycota_ic;Ascomycota_po;Ascomycota_or;Ascomycota_bo;Ascomycota_pf;Ascomycota_fa;Ascomycota_bf;Tetrachaetum;Tetrachaetum_elegans;JX967531_ITS1_HMM;

itsfilename <- "ITSoneDB_rep_seq_1.131.fasta"
itsseq <- system(paste0('grep -v ">" ',itsfilename), intern=T)
itsseq <- toupper(itsseq)

itstax <- system(paste0('grep ">" ',itsfilename, ' | cut -b2-'),intern=T)
itstax <- gsub("['[&(){}/#?\"$%+*:;<>=^~!]]","",itstax) #remove weird characters
itstab <- as.data.frame(do.call(rbind, strsplit(itstax, split='|',fixed=T)),stringsAsFactors=F)
colnames(itstab) <- c("accession","query","taxid","marker")
itstab <- unique(cbind(itstab,itsseq))

#75% have hit from taxid // why don't the rest? hmmm now 99% have
#426 don't have a match and the ones I looked are at because they were very recently merged
itsclass3 <- unlist(mclapply(itstab$taxid,function(x) taxizedb::classification(x,db="ncbi",verbose=T,ask=F),mc.cores=12),recursive=F) #20k ids/min
itsclass4 <-taxize::classification(itstab$taxid[aa],db="ncbi",verbose=T)
itsclass3[names(bb)] <- itsclass4

# 
# #putting in one species id, getting out another :/
# itsclass <- unlist(mclapply(itstab$query,function(x) taxizedb::classification(x,db="ncbi",verbose=T,ask=F),mc.cores=12),recursive=F) #20k ids/min
# #itsclass <- unlist(itsclass,recursive=F)
# 
# #cut names down until the match // ok too much work, just match the genus sp, then just the genus
# query1 <- do.call(rbind,lapply(itstab$query,function(x) strsplit(x,split=' ',fixed=T)[[1]][1:2]))
# query2 <- query1[,1]
# query1 <- paste(query1[,1],query1[,2],sep=" ")
# 
# itsclass1 <- unlist(mclapply(query1[which(names(itsclass)=="")], function(x) taxizedb::classification(x,db="ncbi",verbose=T,ask=F),mc.cores=12),recursive=F) #20k ids/min
# itsclass2 <- unlist(mclapply(query2[which(names(itsclass)=="")], function(x) taxizedb::classification(x,db="ncbi",verbose=T,ask=F),mc.cores=12),recursive=F) #20k ids/min
# 
# istbak <- itsclass
# 
# itslist <- list()
# 
# itslist <- mclapply(itsclass, function(x) {
# 		if(is.null(x)) return(NA)
# 		n <- names(x)
# 		y <- x[[1]]
# 		if(is.na(y)) {print(n); return(NA)}
# #		if(length(grep("^Error",y))) {print(n); return(NA)} #include ask=T to get NAs
# 		print(n)
# 		
# 		y["label"] <- y[which(y["rank"]=="species"),"id"]
# 		return(y)
# 	},mc.cores=12)
# 
# 
# 

attr(itsclass3,"db") <- "ncbi"
attr(itsclass3,"class") <- "classification"
# 
# itsmat <- cbind(itslist)
# itmatbak <- itsmat

itsmat <- cbind(itsclass3)
dd<-itsmat
itsmat <- itsmat[,c(intersect(taxlevels,colnames(itsmat)))]
itsmat["taxid"] <- names(itsclass3)

itsmerge <- merge(itstab,itsmat,by.x="taxid",by.y="taxid",all=T)

itsmergetax <- infilltaxa(itsmerge[,colnames(itsmat)])
itsmerge[,colnames(itsmergetax)] <- itsmergetax

itsouthead <- apply(itsmerge,1,function(x) paste0(paste(x[c(intersect(taxlevels,colnames(itsmergetax)),"accession")],collapse=";"),";"))
x <- cbind(paste0(">",itsouthead),as.character(itsmerge$itsseq))
x <- x[order(x[,1]),]
write.table(x=x,file=paste0("ref_dada2_",itsfilename),sep="\n",quote=F,row.names=F,col.names=F)

#itsmat <- do.call(rbind,itslist)
# itstaxy1 <- getTaxonomy(itstab$taxid,taxaNodes,taxaNames,mc.cores=12,desiredTaxa=taxlevels) #rownames are not taxids
# itstaxy2 <- infilltaxa(itstaxy1,taxlevels) #rownames are not taxids
# itstaxy2 <- as.data.frame(itstaxy2)
# itstaxy2 <- cbind(itstaxy2,"taxid"=rownames(itstaxy2))
# itstaxy2 <- itstaxy2[order(itstaxy2$taxid),]
# itstaxtab <- merge(itstab,itstaxy2,by="row.names")




#BOLD
system("wget -O BOLD_12S.fasta http://v3.boldsystems.org/index.php/API_Public/sequence?marker=12S|g12S")
system("wget -O BOLD_16S.fasta http://v3.boldsystems.org/index.php/API_Public/sequence?marker=16S|g16S")
system("wget -O BOLD_18S.fasta http://v3.boldsystems.org/index.php/API_Public/sequence?marker=18S|g18S")
system("wget -O BOLD_28S.fasta http://v3.boldsystems.org/index.php/API_Public/sequence?marker=28S|g28S|28S-D2|28S-D2-D3|28S-D9-D10")
system("wget -O BOLD_ITS.fasta http://v3.boldsystems.org/index.php/API_Public/sequence?marker=ITS|ITS1|ITS2")
system("wget -O BOLD_58S.fasta http://v3.boldsystems.org/index.php/API_Public/sequence?marker=5.8S")
system("cat BOLD*.fasta > BOLD_all.fasta")

#>ASANO035-09|Camponotus MG005|12S|KJ427394
#TTTTATTAAGAT----AAAATTTAA-TTAATTATTAAATG-TGAATCATTAAAAATAAAAAATTTGGCGGTATTTATTATTTAGAGGAACTTGTTTATTAAATTGATAATCCACGAATGGAAAAACTTAA-T--TT--------TTTTTTATATATTGTTGTTGAAAAATTATATTAATAGGTAATGATTAATAAA-TTTAATTTTTAAA
#has Ns

intax <- system('grep ">" BOLD_all.fasta | cut -b2-',intern=T)
inseq <- system('grep -v ">" BOLD_all.fasta',intern=T)
#intab <- do.call(rbind, strsplit(intax, split='|',fixed=T)) #doesn't work because of missing delimiters
intab <- read.table(text=intax,colClasses="character",sep="|",fill=T)
colnames(intab) <- c("boldid","query","marker","accession")

intab <- unique(cbind(intab,inseq))
intab$inseq <- gsub("-","",intab$inseq)
intab$inseq <- gsub("\r","",intab$inseq)

#taxid1 <- accessionToTaxa(paste0(intab[,4],".1"),"accessionTaxa.sql")
#taxid3 <- get_ids(unique(taxtab[taxtab$taxid2=="NA",2]),db="ncbi")
#taxtab <- cbind(intab[,1:4],taxid1,taxid2,taxid1==taxid2))

markerkeep <- c("12S","16S","28S","18S","ITS","ITS1","ITS2","28S-D2","28S-D3-D5","28S-D1-D2","28S-D2-D3","18S-V4","5.8S")
keeptab <- intab[intab$marker %in% markerkeep,]

# testn <- ncol(ke
# keeptab <- keeptab[1:testn,]

#bold crashes on classification lookup on Phyla
bolds <- bold_search(sort(unique(keeptab$query[!grepl(pattern=" ",x=keeptab$query)])),verbose=T) #just search one-word taxa
boldphy <- bolds$input[bolds$tax_rank=="phylum"]
#taxid2 <- getId(keeptab$query,taxaNames) #ambiguous have commas
#taxy2 <- getTaxonomy(taxid2,taxaNodes,taxaNames,mc.cores=12,desiredTaxa="phylum")

# library(lme4) #for mclapply
# library(parallel) #for mclapply
#do.call(rbind,bold_search(unique(keeptab$query[grepl(pattern=" ",x=keeptab$query)])))
#bolds <- mclapply(unique(keeptab$query[!grepl(pattern=" ",x=keeptab$query)]),bold_search)

phytab <- keeptab[which(keeptab$query %in% boldphy),]
#phytab <- keeptab[which(taxy2==keeptab$query),]
classphy <- taxize::classification(unique(phytab$query),db="ncbi",return_id=T,rows=1)
mergephy <- merge(phytab,cbind(classphy))

#keeptab <- keeptab[-which(taxy2==keeptab$query),]
keeptab <- keeptab[-which(keeptab$query %in% boldphy),]

#get them 50 at a time with failures
queries <- sort(unique(keeptab$query))
numrounds <- (length(queries)%/%50)
firsts <- 50*0:numrounds+1
seconds <- 50*1:numrounds
classkeep <- list()
for(i in 1:numrounds) {
    trx <- try(error,silent=T)
    while(class(trx)=="try-error") { 
    print(i)
    trx = try({
        tosave <- taxize::classification(queries[firsts[i]:seconds[i]],db="bold",return_id=T,ask=F)
        },silent=TRUE)
    }
    classkeep <- c(classkeep,tosave)
}

ck <- names(classkeep) %in% names(which((unlist(lapply(classkeep,is.na)))))
boldlist <- classkeep[!ck]
attr(boldlist,"db") <- "bold"
attr(boldlist,"class") <- "classification"
boldmat <- cbind(boldlist)
#boldmat <- boldmat[,c(intersect(taxlevels,colnames(boldmat)),"query")]
#aa <- boldmat[,intersect(colnames(mergephy),colnames(boldmat))]
mergekeep <- merge(keeptab,boldmat,all=T)

mergeall <- merge(mergephy,mergekeep,all=T)
mergeall<- mergeall[,c(intersect(taxlevels,colnames(mergeall)),"boldid","marker","inseqs")]
mergefill <- as.data.frame(infilltaxa(mergeall))
mergefill$inseq <- gsub("N","",mergefill$inseq)


for (m in unique(mergefill$marker)) {
	print(m)
    mk <- mergeall$marker==m
    mergetaxout <- apply(mergefill[mk,],1,function(x) paste0(paste(x[c("phylum","class","order","family","genus","species","boldid","marker")],collapse=";"),";"))
    x <- cbind(paste0(">",mergetaxout),as.character(mergefill$inseq[mk]))
    x <- x[order(x[,1]),]
    write.table(x=x,file=paste0("ref_dada2_BOLD_",m,".fasta"),sep="\n",quote=F,row.names=F,col.names=F)
}





### GET TAXONOMY LEVELS FROM NCBI
### only kind of works
# function to get consensus taxonomy from NCBI
get_consensus_tax <- function(x) {
	aa <- unclass(rle(sort(x)))
	ab <- aa$values[rev(order(aa$lengths))][1:20]
	ab <- ab[!is.na(ab)]
	
	af <- as.data.frame(aa)
	rownames(af) <- af$values
	ag <- do.call(rbind,get_uid_(ab,db="ncbi"))
	if(is.null(ag)) return("unknown")
	ah <- merge(af,ag,by.x="row.names",by.y="scientificname")
	getmode(ah$rank)
}

# mode function
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}


# automatically assign taxonomic levels
refdbs <- list.files(".",pattern="ref_dada2_.*fasta")

for(ref in refdbs) {

    # UNITE is treated specially by DADA2
    if(ref=="ref_dada2_unite.fasta") {
        next
    }

    print(ref)
    
    somelabels <- system(paste0("grep '>' ",ref," | cut -b2-"),intern=T)
    labelmat <- as.data.frame(do.call(rbind, strsplit(somelabels, split=';',fixed=T)),stringsAsFactors=F) #last column is recyled

    taxlevellist <- list()
    labeltl <- list()
    for(i in 1:ncol(labelmat)) {
        print(i)
        labelcheck <- labelmat[,i]
        labelcheck <- labelcheck[!grepl("_",labelcheck)]
        labeltl <- c(labeltl,get_consensus_tax(labelcheck))
        print(labeltl)
    }
    taxlevellist[[ref]] <- unlist(labeltl)
}
saveRDS(taxlevellist, file.path(outfolder,"taxlevels.rds"))

