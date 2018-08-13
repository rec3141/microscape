library(shiny)
library(ape) #phylogenetics
#library(phytools) #
library(phyloseq)
library(ggplot2)
theme_set(theme_bw())
library(msa)
library(phangorn)
# source("microscape-s99-shiny.r"); shinyApp(ui = ui, server = server)

library(gdata) # for read.xls
library(marmap) # for fortify.bathy and plot.bathy
library(ncdf4) # for netCDF
library(vegan) #for metaMDS
library(Rtsne) #for Rtsne
library(spdep) #for Rotation
library(tidyr) 

# each tab filters the next? 
# no... should they affect each other?
# Step 1) Select samples
# Step 2) Select taxa
# Step 3) Select geography
# Step 4) Select environmental

#### load required data

outfolder <- "out_dada"
options(expressions=50000)
system("ulimit -s 16384")

taxout <- readRDS(file.path(outfolder,"taxout.rds"))
bootout <- readRDS(file.path(outfolder,"bootout.rds"))
proptab <- readRDS(file.path(outfolder,"proptab_final.rds"))
tt.cols <- ncol(proptab)
tt.rows <- nrow(proptab)
seqtab <- readRDS(file.path(outfolder,"seqtab_final.rds"))
rs.seqtab <- rowSums(seqtab)
sample.cex <- log10(rs.seqtab)/3

# old, may need to redo this but probs not
# sample.tsp.bray <- readRDS(file.path(outfolder,"sample_bray_tsp.rds"))
# convert_small <- proptab[sample.tsp.bray,]

#this sets up the taxa lists for each reference database
taxa.list <- list()
names.list <- list()
for(ref in names(taxout)) {
    print(ref)

    taxa.list[[ref]] <- unname(taxout[[ref]][colnames(proptab),6])
    for(i in ncol(taxout[[ref]]):1) {
        taxa.list[[ref]][is.na(taxa.list[[ref]])] <- taxout[[ref]][colnames(proptab)[is.na(taxa.list[[ref]])],i]
    }
    taxa.list[[ref]] <- c("",taxa.list[[ref]])

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
    print(taxlevels)
    colnames(taxout[[ref]]) <- taxlevels

    taxchoose <- intersect(taxlevels,c("Superkingdom","Kingdom","Phylum","Order","Family","Genus","Species","Accession"))
    names.list[[ref]] <- as.list(tidyr::unite(data.frame(taxout[[ref]][colnames(proptab),]), taxchoose, sep=";"))
    esvs <- paste0("ESV_",1:ncol(proptab))
    names.list[[ref]] <- paste(names.list[[ref]],esvs,sep=";")

    mapply(c, first, second, SIMPLIFY=FALSE)
    
}

# this sets up the sample list from the metadata sheet
metadata <- readRDS(file.path(outfolder,"metadata.rds"))
sample.list <- c("",metadata$desc)

# select which mappings to use for samples
sample.tsne <- readRDS(file.path(outfolder,"sample_bray_tsne.rds"))
# rotate for better view
sample.tsne.rot <- Rotation(sample.tsne,90*pi/180)
tsne.here <- data.frame(sample.tsne.rot)
colnames(tsne.here) <- c("x","y")

# select which mappings to use for sequences
seq.tsne <- readRDS(file.path(outfolder,"seq_bray_tsne.rds"))


leg.taxa <- unlist(lapply(newnames,function(x) paste0(strsplit(x,";")[[1]][6:7],collapse=" ")))

# load("save.weights.Rdata")
# load("save.graphs.Rdata")

#colors for biospatial plot
# color.biospatial <- rainbow(4*((ncol(proptab) %/% 4) +1)) #did this so that nearby taxa had different colors to easily distinguish them
# color.biospatial <- as.character(t(matrix(color.biospatial,nrow=4)))
# color.biospatial <- color.biospatial[order(newnames)]
color.biospatial <- rainbow(ncol(proptab))[order(newnames)]

#default colors for network plot
color.network <- rep("#427df4",times=ncol(proptab)) #blue default Bacteria
color.network[grep("Alphaproteobacteria",newnames)] <- "#0033ff" #blue
color.network[grep("Betaproteobacteria",newnames)] <- "#1897ff" #blue
color.network[grep("Gammaproteobacteria",newnames)] <- "#1dcdf9" #blue
color.network[grep("Bacteroidetes",newnames)] <- "#0ee0d2" #blue
color.network[grep("Archaea",newnames)] <- "#f9f03e" #yellow
color.network[grep("Eukaryota",newnames)] <- "#f44265" #red
color.network[grep("Metazoa",newnames)] <- "#f463f9" #purple
color.network[grep("Chloroplast",newnames)] <- "#0fe047" #green
color.network[grep("Ochrophyta",newnames)] <- "#0ee047" #green
color.network[grep("Chlorophyta",newnames)] <- "#0de047" #green
color.network[grep("Syndiniales",newnames)] <- "#f91bf5" #purple
color.network[grep("Mitochondria",newnames)] <- "#42e5f4" #orange

tax.levels <- colnames(taxout$ref_dada2_silva_nr_v132_train_set.fasta$tax)

#write to fasta format
dada2fasta <- function(esvs,taxanames,taxon) {
    which.taxa <- grepl(taxon,taxanames)
    saved <- paste0(">",taxanames[which.taxa],"\n",colnames(esvs)[which.taxa],"\n")
    cat(saved,sep="")
}

#phylo
# do this outside and then subset it for plotting
# convert selected sequences to DNAbin
# outfolder <- rev(dir("./",pattern="out_dada*"))[1]
# but not yet merged
## seqtab <- readRDS(paste0(outfolder,"/seqtab_final.rds")) # 
## taxout <- readRDS(paste0(outfolder,"/taxout.rds"))
# seqs.raw <- sapply(strsplit(colnames(proptab),""), tolower)
# 
# names(seqs.raw) <- colnames(proptab)
# seqs.bin <- as.DNAbin(seqs.raw)
# # read in metadata
# 
# seqs.otu <- otu_table(t(proptab), taxa_are_rows = TRUE)
# seqs.tax <- tax_table(taxout$ref_dada2_silva_nr_v132_train_set.fasta$tax[colnames(proptab),])
# seqs.physeq <- phyloseq(seqs.otu, seqs.tax)
# save(seqs.physeq,file="seqs.physeq.Rdata")
#load(file="seqs.physeq.Rdata")

#        seqs.meta <- sample_data()
#        seqs.physeq <- merge_phyloseq(seqs.physeq, sampledata, seqs.tree)


#mapping


# open the netCDF file
nc <- nc_open("/scratch/SKQ201813S/underway/ETOPO2v2g_f4.nc")
etopo.in <- ncvar_get(nc, "z") #get depths

# do some stuff I copied from the web
colnames(etopo.in) <- ncvar_get(nc, "y")
rownames(etopo.in) <- ncvar_get(nc, "x")
class(etopo.in) <- "bathy"

########### START WEBSITE
# Define UI for app that draws a histogram ----
ui <- fluidPage(

  # App title ----
  titlePanel("Microscape v0.00001"),

    fluidRow(
        column(
    
      # Input: Slider for the number of bins ----
      selectInput(inputId = "taxon", label = "Show single taxon", choices = taxa.list),
      textAreaInput(inputId = "plot.taxa.regex", label="Show all matching taxa (regex)", value=NULL),
      selectInput(inputId = "sample", label = "Show single sample", choices = sample.list),
      textAreaInput(inputId = "plot.sample.regex", label="Show all matching samples (regex)", value=NULL),
        #slider value set to one higher than maximum reported network correlation
      sliderInput("network.cor", "N: Minimum SparCC correlation", min = 0.01, max = 1, value = (0.01+length(save.weights)/100), step=0.01),
      sliderInput("min.taxa.hits", "N: Minimum samples per ESV (log2)", min = 0.0, max = 12.0, value = 0),
      sliderInput("min.reads", "S: Minimum sample reads (log2)", min = 0.0, max = 20.0, value = 17.0),
      sliderInput("min.richness", "S: Minimum sample richness (log2)", min = 0.0, max = 10.0, value = 1.0),

      width=3,class = "well"
    ),

    
    column(
     tabsetPanel(type = "tabs", id = "tabs",
        tabPanel("Samples", value="biospatial",
            span(textOutput("samplename"), style="font-size:120%; text-align:center"),
            div(style = "height:700px;", plotOutput("biospatial", hover = "biospatial_hover", click = "biospatial_click", dblclick="biospatial_dblclick"))
        ),
        tabPanel("Network", value="network",
            span(textOutput("taxaname"), style="font-size:95%; text-align:center"),
            div(style = "height:700px;", plotOutput("network", hover = "network_hover", click = "network_click", dblclick="network_dblclick"))
        ),
        tabPanel("Map", plotOutput("geographic")),
        tabPanel("Environmental", plotOutput("environmental")),
        tabPanel("Phylogenetic", value="phylogenetic", 
             div(style = "height:1000px;", plotOutput("phylogenetic", hover = "phylogenetic_hover", click = "phylogenetic_click", dblclick="phylogenetic_dblclick"))
        ),
        tabPanel("Functional", plotOutput("functional"))
        

     ),width=5, height=12, class="well"
    ),
    column(
      div(tableOutput("overunder"), style = "font-size:80%"),
      div(tableOutput("rhtable"), style = "font-size:80%"),
      sliderInput("plot.all.names", label="Show some random names", min=0, max=25, value = 0),
      textAreaInput(inputId = "plot.names", label="Show these names (regex)", value=NULL),
      sliderInput("point.scalar", "Point size scalar (log)", min = -20.0, max = 20.0, value = 0),
      br(),
      actionButton(inputId="reset.zoom", label="Reset Zoom"),
      br(),
      sliderInput("tree.height", "P: Tree height limits", min = -1000, max = 2000.0, value = c(0,20.0)),
      sliderInput("tree.width", "P: Tree width limits", min = -500, max = 2000, value = c(0,100)),
      sliderInput("tree.font.size", "P: Tree font size", min = 0, max = 20, value = 10),
      sliderInput("tree.colors", "P: Rotate tree colors", min = 0, max = 6, value = c(0,0), step=0.1),
      selectInput("tax.levels","P: Taxonomic levels", tax.levels, multiple=TRUE, selected="genus"),
      sliderInput("longitude", "M: Longitude", min = -200, max = 200, value = c(-170,-150)),
      sliderInput("latitude", "M: Latitude", min = -90, max = 90, value = c(55,75)),

      width=4, class="well"
    )
  )
)



# Define server logic
server <- function(input, output, session) {

    # set up reactive values
    ranges <- reactiveValues(biospatial_xlim=NULL, biospatial_ylim=NULL, network_xlim=c(-60,60), network_ylim=c(-60,60))
    oldxy <- reactiveValues(biospatial_x=0, network_x=0)
    oldzoom <- reactiveValues(z=0)
    picks <- reactiveValues(taxa=NULL,samples=NULL)
    point.size <- reactiveValues(biospatial=0, network=0, input=NULL)
    color <- reactiveValues(phylo=rep(NA,ncol(proptab)))
    last <- reactiveValues(seqs=NULL); isolate(last$seqs)

    # some functions
   xy_biospatial <- function(e) {
        if(is.null(e)) return(NULL)
        new.pos <- c(e$x, e$y) # New position
        # Compute distance to points and select nearest index
        tmp.dist <- colSums((t(sample.tsne.rot) - new.pos)^2)
        tmp.dist[!picks$samples] <- max(tmp.dist)
        nearest.idx <- which.min(tmp.dist)
        nearest.idx
    }
    
   xy_network <- function(e) {
        if(is.null(e)) return(NULL)
        new.pos <- c(e$x, e$y) # New position
        # Compute distance to points and select nearest index
        tmp.dist <- colSums((t(seq.tsne) - new.pos)^2)
        tmp.dist[!picks$taxa] <- max(tmp.dist)
        nearest.idx <- which.min(tmp.dist)
        nearest.idx
    }

    
    ###############
    # HEADERS
    ###############

    output$samplename <- renderText({
        if(is.null(input$biospatial_hover)) {
            as.character("​") #had to use zero width space (U+200B) to avoid wiping
        } else {
            #require it to be picked
            if(picks$samples[xy_biospatial(input$biospatial_hover)]) {
                metadata$desc[xy_biospatial(input$biospatial_hover)]
            } else {
                as.character("​") #had to use zero width space (U+200B) to avoid wiping
            }
        }
    })

    output$taxaname <- renderText({
        if(is.null(input$network_hover)) {
            as.character("​") #had to use zero width space (U+200B) to avoid wiping
        } else {
            if(picks$taxa[xy_network(input$network_hover)]) {
                newnames[xy_network(input$network_hover)]
            } else {
                as.character("​") #had to use zero width space (U+200B) to avoid wiping
            }
        }
    })

    ###############
    # biospatial plot
    ###############

      output$biospatial <- renderPlot({
      
      # MECHANICS AND ZOOM
        if(is.null(ranges$biospatial_xlim)) ranges$biospatial_xlim <- c(-50,50)
        if(is.null(ranges$biospatial_ylim)) ranges$biospatial_ylim <- c(-50,50)
        point.size$biospatial <- input$point.scalar
        
        if(input$reset.zoom > oldzoom$z & input$tabs=="biospatial") {
            ranges$biospatial_xlim=c(-50,50)
            ranges$biospatial_ylim=c(-50,50)
            oldzoom$z <- input$reset.zoom
        }

        # what to do if user clicks (zoom out) or double clicks (zoom in)
        if(!is.null(input$biospatial_dblclick$x)) {

            if(input$biospatial_dblclick$x != oldxy$biospatial_x) {

            e <- input$biospatial_dblclick
            xlim.diff = 0.4*(ranges$biospatial_xlim[2]-ranges$biospatial_xlim[1])
            ranges$biospatial_xlim <- c(e$x - xlim.diff/2, e$x + xlim.diff/2)
    
            ylim.diff = 0.4*(ranges$biospatial_ylim[2]-ranges$biospatial_ylim[1])
            ranges$biospatial_ylim <- c(e$y - ylim.diff/2, e$y + ylim.diff/2)
    
            oldxy$biospatial_x <- input$biospatial_dblclick$x
            }
        } else if (!is.null(input$biospatial_click$x)) {

            if(input$biospatial_click$x != oldxy$biospatial_x) {
            e <- input$biospatial_click
            xlim.diff = 1.25*(ranges$biospatial_xlim[2]-ranges$biospatial_xlim[1])
            ranges$biospatial_xlim <- c(e$x - xlim.diff/2, e$x + xlim.diff/2)
    
            ylim.diff = 1.25*(ranges$biospatial_ylim[2]-ranges$biospatial_ylim[1])
            ranges$biospatial_ylim <- c(e$y - ylim.diff/2, e$y + ylim.diff/2)
    
            oldxy$biospatial_x <- input$biospatial_click$x
            }
        }


        # PROCESS USER INPUT

        #### SAMPLES FIRST
        # Sample filters
        min.reads.samples <- rs.seqtab > 2^input$min.reads
        min.richness.samples <- rowSums(proptab>0) >= 2^input$min.richness

        # initialize 
        sample.user.regex <- rep(FALSE,tt.rows)
        # print(input$sample)
        # if user doesn't provide a lone sample or a regex, select all
        if(input$sample=="" & input$plot.sample.regex=="") {
            sample.user.regex <- rep(TRUE,tt.rows)
        # if user provides a lone sample and no regex, grep for it
        } else if (input$sample!="" & input$plot.sample.regex=="") {
            sample.user.regex <- metadata$desc==input$sample
        } else if(input$plot.sample.regex!="") {
        # if user provides a list of regex, split and search each one, ORing
        # override input$sample
            patterns <- strsplit(x=input$plot.sample.regex,split="\n")[[1]]
            for(patt in patterns) {
                # print(patt)
                if(patt=="\n") next
                sample.user.regex <- sample.user.regex | grepl(x=metadata$desc, pattern=patt)
            }
        } else {
            print("biospatial samples: how did you get here?")
        }


        # don't waste time calculating and plotting things that aren't in the view box
        samples.in.box <- sample.tsne.rot[,1] > ranges$biospatial_xlim[1] & 
        sample.tsne.rot[,1] < ranges$biospatial_xlim[2] & 
        sample.tsne.rot[,2] > ranges$biospatial_ylim[1] & 
        sample.tsne.rot[,2] < ranges$biospatial_ylim[2]

        # put all the sample filters together
        # QC and holdup

        picks$samples <- min.reads.samples & min.richness.samples & sample.user.regex & samples.in.box
        # print(paste(metadata$desc[picks$samples],sep="\n"))
        sps <- sum(picks$samples)

        if(sps==0) {print("biospatial stopped"); return(0)}

        map.samples <- min.reads.samples & min.richness.samples & samples.in.box
        # plot the basic outline of selected samples
        plot(sample.tsne.rot[map.samples,],pch=19,cex=sample.cex[picks$samples],col="grey",axes=FALSE,ann=FALSE,xlim=ranges$biospatial_xlim,ylim=ranges$biospatial_ylim)



        #### THEN TAXA

        # initialize
        taxa.user.regex <- rep(FALSE,tt.cols)

        # if user doesn't provide a taxon or a regex, select all
        # print(input$taxon)
        if(input$taxon=="" & input$plot.taxa.regex=="") {
            taxa.user.regex <- rep(TRUE,tt.cols)
        # else if user provides a lone taxon and no regex, grep it
        } else if (input$taxon!="" & input$plot.taxa.regex=="") {
            taxa.user.regex <- grepl(x=newnames,pattern=input$taxon)
        } else if(input$plot.taxa.regex!="") {
        # if user provides a list of regex, split and search each one, ORing
            patterns <- strsplit(x=input$plot.taxa.regex,split="\n")[[1]]
            for(patt in patterns) {
                # print(patt)
                if(patt=="\n") next
                taxa.user.regex <- taxa.user.regex | grepl(x=newnames, pattern=patt)
            }
        } else {
            print("biospatial taxa: how did you get here?")
        }

        # filter by what's shown in taxa box
        taxa.in.box <- seq.tsne[,1] > ranges$network_xlim[1] & 
        seq.tsne[,1] < ranges$network_xlim[2] & 
        seq.tsne[,2] > ranges$network_ylim[1] & 
        seq.tsne[,2] < ranges$network_ylim[2]

        #select taxa in multiple samples
        taxa_mult_samples <- colSums(proptab>0) >= 2^input$min.taxa.hits 
#         picks$taxa <- picks$taxa & taxa_mult_samples #infinite loop?

        # enforce slider limits over taxa choices
        picks$taxa <- taxa.user.regex & taxa_mult_samples & taxa.in.box
        # print(paste(newnames[picks$taxa],sep="\n"))

        # print(sum(picks$taxa))
        ptl <- length(picks$taxa)

        spt <- sum(picks$taxa)
        if(spt==0) { print("biospatial stopped");return(0)} #no valid taxa selected

        # --> diverge

        #only color selected taxa -- not sure if this is necessary
        color.taxa <- rep(NA,length(picks$taxa))
        color.taxa[which(picks$taxa)] <- color.biospatial[which(picks$taxa)]

        # if defined by phylogeny, add some color for plotting points
        color.taxa[which(!is.na(color$phylo))] <- color$phylo[which(!is.na(color$phylo))]

        # for each sample, sort pick.taxa from greatest to least, plot them
        for(n in which(picks$samples)) {
            if(n %% 100 == 0) print("working")
            seq.ra <- unname(proptab[n,,drop=F])
#             print(which(picks$taxa))
            seq.ra[!picks$taxa] <- 0
            seq.ra[seq.ra==0] <- NA
            seq.order <- order(seq.ra,decreasing=T,na.last=NA)
#             print(seq.order)

            vsize <- rep(NA,length(picks$taxa))
            vsize[seq.order] <- unname(sqrt(sqrt(proptab[n,seq.order,drop=F]/1000)))
            vsize <- vsize*1.1^point.size$biospatial

            # PLOT 4: Sample map with circles proportional to relative abundance in n
            tsne.plot <- cbind(rep(tsne.here[n,1],ptl),rep(tsne.here[n,2],ptl))
            tsne.plot <- data.frame(tsne.plot)
            colnames(tsne.plot) <- c("x","y")
    #        tsne.plot.taxa <- rep("",length(pick.taxa))
    #        tsne.plot.taxa[pick.taxa] <- newnames[which(pick.taxa)]
#             print(sort(vsize))
            points(tsne.plot,pch=19,cex=vsize,col=color.taxa)    
        }


#         print(input$plot.all.names)
#         print(input$plot.names)
        
         if(input$plot.all.names > 0 & input$tabs=="biospatial") {
                small.sample <- sample(which(picks$samples),input$plot.all.names, replace=T)
                if(length(small.sample)>0) text(sample.tsne.rot[small.sample,,drop=F],labels=metadata$desc[small.sample])
         }
#                   if(input$plot.all.names & input$tabs=="biospatial") text(sample.tsne.rot[picks$samples,],labels=metadata$desc[picks$samples])

         if(input$plot.names!="") {
            # print(input$plot.sample.names)
            patterns <- strsplit(x=input$plot.names,split="\n")[[1]]
            for(patt in patterns) {
                # print(patt)
                if(patt=="\n") next
                pick.user.samples <- grepl(pattern=patt, x=metadata$desc)
                if(sum(pick.user.samples)>0) text(sample.tsne.rot[pick.user.samples,,drop=F],labels=metadata$desc[pick.user.samples])
            }
        }
        print("done biospatial")
    }, width=700,height=700)


    ###############
    # network plot
    ###############
    output$network <- renderPlot({
    
    # pre-calculate weighted networks and load
    # to do? click on a taxon to show that taxon's network
    
        ### MECHANICS AND ZOOM
        if(input$reset.zoom > oldzoom$z & input$tabs=="network") {
            ranges$network_xlim=c(-60,60)
            ranges$network_ylim=c(-60,60)
            oldzoom$z <- input$reset.zoom
        }

        point.size$network <- input$point.scalar
        
        if(!is.null(input$network_dblclick$x)) {

            if(input$network_dblclick$x != oldxy$network_x) {

            e <- input$network_dblclick
            network_xlim.diff = 0.4*(ranges$network_xlim[2]-ranges$network_xlim[1])
            ranges$network_xlim <- c(e$x - network_xlim.diff/2, e$x + network_xlim.diff/2)

            network_ylim.diff = 0.4*(ranges$network_ylim[2]-ranges$network_ylim[1])
            ranges$network_ylim <- c(e$y - network_ylim.diff/2, e$y + network_ylim.diff/2)

            oldxy$network_x <- input$network_dblclick$x
            }
        } else if (!is.null(input$network_click$x)) {

            if(input$network_click$x != oldxy$network_x) {
            e <- input$network_click
            network_xlim.diff = 1.25*(ranges$network_xlim[2]-ranges$network_xlim[1])
            ranges$network_xlim <- c(e$x - network_xlim.diff/2, e$x + network_xlim.diff/2)

            network_ylim.diff = 1.25*(ranges$network_ylim[2]-ranges$network_ylim[1])
            ranges$network_ylim <- c(e$y - network_ylim.diff/2, e$y + network_ylim.diff/2)

            oldxy$network_x <- input$network_click$x
            }
        }

   
        ### PROCESS USER INPUT
        
        ### TAXA FIRST
      
        # initialize
        taxa.user.regex <- rep(FALSE,tt.cols)

        # if user doesn't provide a taxon or a regex, select all
        # print(input$taxon)
        if(input$taxon=="" & input$plot.taxa.regex=="") {
            taxa.user.regex <- rep(TRUE,tt.cols)
        # else if user provides a lone taxon and no regex, grep it
        } else if (input$taxon!="" & input$plot.taxa.regex=="") {
            taxa.user.regex <- grepl(x=newnames,pattern=input$taxon)
        } else if(input$plot.taxa.regex!="") {
        # if user provides a list of regex, split and search each one, ORing
            patterns <- strsplit(x=input$plot.taxa.regex,split="\n")[[1]]
            for(patt in patterns) {
                # print(patt)
                if(patt=="\n") next
                taxa.user.regex <- taxa.user.regex | grepl(x=newnames, pattern=patt)
            }
        } else {
            print("network taxa: how did you get here?")
        }
        
        #select taxa in multiple samples
        taxa_mult_samples <- colSums(proptab>0) >= 2^input$min.taxa.hits 

       ## Plot the basic chart using slider limits
        plot(seq.tsne[taxa_mult_samples,], pch=19, col="grey", cex=0.3, xlim=ranges$network_xlim, ylim=ranges$network_ylim, axes = FALSE, xlab = "", ylab = "")

        # enforce slider limits over taxa choices
        picks$taxa <- taxa.user.regex & taxa_mult_samples
        # print(paste(newnames[picks$taxa],sep="\n"))

        # print(sum(picks$taxa))
        ptl <- length(picks$taxa)


        spt <- sum(picks$taxa)
        if(spt==0) {print("network stopped"); return(0)} #no valid taxa selected



        ### SAMPLES NEXT

        # Sample filters
        min.reads.samples <- rs.seqtab > 2^input$min.reads
        min.richness.samples <- rowSums(proptab>0) >= 2^input$min.richness

        # initialize 
        sample.user.regex <- rep(FALSE,tt.rows)

        # print(input$sample)
        # if user doesn't provide a lone sample or a regex, select all
        if(input$sample=="" & input$plot.sample.regex=="") {
            sample.user.regex <- rep(TRUE,tt.rows)
        # if user provides a lone sample and no regex, grep for it
        } else if (input$sample!="" & input$plot.sample.regex=="") {
            sample.user.regex <- metadata$desc==input$sample
        } else if(input$plot.sample.regex!="") {
        # if user provides a list of regex, split and search each one, ORing
        # override input$sample
            patterns <- strsplit(x=input$plot.sample.regex,split="\n")[[1]]
            for(patt in patterns) {
                # print(patt)
                if(patt=="\n") next
                sample.user.regex <- sample.user.regex | grepl(x=metadata$desc, pattern=patt)
            }
        } else {
            print("network samples: how did you get here?")
        }

        
        # don't waste time calculating and plotting things that aren't in the view box
        samples.in.box <- sample.tsne.rot[,1] > ranges$biospatial_xlim[1] & 
        sample.tsne.rot[,1] < ranges$biospatial_xlim[2] & 
        sample.tsne.rot[,2] > ranges$biospatial_ylim[1] & 
        sample.tsne.rot[,2] < ranges$biospatial_ylim[2]

       # put all the sample filters together
        # QC and holdup

        picks$samples <- min.reads.samples & min.richness.samples & sample.user.regex & samples.in.box
        # print(paste(metadata$desc[picks$samples],sep="\n"))
        sps <- sum(picks$samples)

         if(sps==0) {print("network stopped");return(0)}

        # add some color for plotting points
        color.taxa <- rep(NA,ncol(proptab))
        color.taxa[which(picks$taxa)] <- color.network[which(picks$taxa)]
        color.taxa[which(!is.na(color$phylo))] <- color$phylo[which(!is.na(color$phylo))]
#         print(color$phylo[which(!is.na(color$phylo))])
        # --> diverge


        # select taxa in requested network space: only plot if SparCC correlation passes threshold
        
        # if input$network.cor is too high for selected taxa, lower it
        network_cor <- input$network.cor
        # bring in saved weights
        while(is.null(save.weights[network_cor*100][[1]])) {
            network_cor <- network_cor - 0.01
        }
        isolate({updateSliderInput(session, "network.cor", value = network_cor)})
        network_save_weights <- save.weights[[network_cor*100]]
        # select positive and negative correlations
        networked_taxa <- (network_save_weights[picks$taxa,,drop=F] > network_cor) | (network_save_weights[picks$taxa,,drop=F] < -1*network_cor)
 
        if(spt>1) {
            networked_taxa <- colSums(networked_taxa) > 0 # collapse samples
        }
        # re-enforce limits and merge here 
        networked_taxa <- networked_taxa & taxa_mult_samples
        networked_taxa <- networked_taxa | picks$taxa

        taxa.in.box <- seq.tsne[,1] > ranges$network_xlim[1] & 
        seq.tsne[,1] < ranges$network_xlim[2] & 
        seq.tsne[,2] > ranges$network_ylim[1] & 
        seq.tsne[,2] < ranges$network_ylim[2]


        # scale the point and edge sizes
        network_asize <- unname(proptab[picks$samples,])
        if(sps>1) network_asize <- colMeans(network_asize) #collapse
#        network_asize <- 250*sqrt(sqrt(network_asize * 1.8^point.size$network/1000)) # janky sizing
        network_asize <- 250*(sqrt(network_asize * 1.8^point.size$network/1000)) # janky sizing
        network_asize[network_asize==0] <- NA #don't plot zeros
        network_asize[!networked_taxa] <- NA #don't plot unnetworked
        network_asize[!taxa.in.box] <- NA #don't plot if not in box
    
        # get the network from the saved graphs file
        network_save_graph <- save.graphs[[network_cor*100]]
        # resize the edges
        network_edgesize <- rescale(E(network_save_graph)$weight,to=c(3,30))

        # plot edges and vertices
        if(!all(is.na(network_asize))) {
            plot(network_save_graph, layout=seq.tsne, xlim=ranges$network_xlim, ylim=ranges$network_ylim, vertex.shape="circle", vertex.size=network_asize, vertex.color=color.taxa, vertex.label=NA, vertex.label.cex=1.5, vertex.label.color="black", edge.label=NA, edge.width=network_edgesize, main=NULL, rescale=F, add=T)
        }

        # add random text annotations if requested
         if(input$plot.all.names > 0 & input$tabs=="network") {
            shown_taxa <- networked_taxa & taxa.in.box
                small.sample <- sample(which(shown_taxa),input$plot.all.names, replace=T)
                if(length(small.sample)>0) text(seq.tsne[small.sample,,drop=F],labels=newnames[small.sample])
         }

        # add specific text annotations if requested 
         if(input$plot.names!="") {
            patterns <- strsplit(x=input$plot.names,split="\n")[[1]]
            for(patt in patterns) {
                if(patt=="\n") next
                pick.user.taxa <- grepl(pattern=patt, x=newnames)
                if(sum(pick.user.taxa)>0) text(seq.tsne[pick.user.taxa,,drop=F],labels=sapply(which(pick.user.taxa), function(x) paste(strsplit(newnames[x],split=";")[[1]][c(6,7)],collapse=" ")))
            }
        }
        print("done network")
    }, width=700,height=700)







    ###############
    # geographic plot
    ###############
    output$geographic <- renderPlot({
    
        # set up the map range for ggplot
        range.lat <- input$latitude
        range.lon <- input$longitude

        print(range.lat)
        print(range.lon)
        # you can change the limits of the base map here
        lat.s <- which(as.numeric(colnames(etopo.in))==range.lat[1])
        lat.n <- which(as.numeric(colnames(etopo.in))==range.lat[2])

        lon.e <- which(as.numeric(rownames(etopo.in))==range.lon[1])
        lon.w <- which(as.numeric(rownames(etopo.in))==range.lon[2])

        etopo.chukchi <- etopo.in[lon.e:lon.w, lat.s:lat.n]
        class(etopo.chukchi) <- "bathy"

        # this reformats the data into a ggplot-friendly format
        chukchi.df <- fortify(etopo.chukchi)

        # we'll plot the map using ggplot2
        # ggplot is a layer-based plotting system, so we'll set up the layers one by one and then combine them at the end

        # prepare the base topographic map

        # each line adds a new layer of data or information, make sure to include the '+' at the end of each line before the final
        # the n=10 tells ggplot how many levels to use to interpolate the colors
        m.base <- ggplot(chukchi.df, aes(x=x, y=y)) +
            geom_raster(aes(fill=z), data=chukchi.df) + scale_fill_etopo(n=10) + # the raster image of the bathymetry
            scale_x_continuous(limits=c(range.lon)) + scale_y_continuous(limits=c(range.lat)) # limit the axes to the map range given previously

        # to more specifically focus on the ocean, we'll download a nice bathymetric color scheme from the shared drive, "scale_fill_etopo_bathy"
        source("http://share.sikuliaq.alaska.edu/public/Cruises/SKQ201813S/ObservationalOceanography/Eric/etopo_colors.r")

        # the limits tell ggplot which data to plot, here we're greying out everything below 3000 m and above 0 m
        m.base <- ggplot(chukchi.df, aes(x=x, y=y)) +
            geom_raster(aes(fill=z), data=chukchi.df) + scale_fill_etopo_bathy(n=10,limits=c(-3000,0)) + # the raster image of the bathymetry
            scale_x_continuous(limits=c(range.lon)) + scale_y_continuous(limits=c(range.lat)) # limit the axes to the map range given previously

        # you can add more layers as necessary, here we'll add some contour lines
        m.fancy <- m.base +
            geom_contour(aes(z=z),breaks=c(-10, -20, -30, -40, -50, -60, -70, -80, -90, -100),colour="white", size=0.2) + # contour lines for shallow depths
            geom_contour(aes(z=z),breaks=c(-200, -250, -300, -350, -400, -450),colour="grey", size=0.2) + # contour lines for medium depths
            geom_contour(aes(z=z),breaks=c(-500, -1000, -1500, -2000, -2500, -3000, -3500, -4000),colour="blue", size=0.2) + # contour lines for deep depths
            geom_contour(aes(z=z),breaks=c(0),colour="black", size=0.4) # contour lines for the shoreline

        # we can also read in some waypoints from previous cruises
#         waypoint.in <- read.xls(xls="http://share.sikuliaq.alaska.edu/public/Cruises/SKQ201813S/ObservationalOceanography/Eric/Station_Plan.xls",sheet="Plan")
        waypoint.data <- waypoint.in[,1:4]
        waypoint.data <- waypoint.data[complete.cases(waypoint.data),]

        # give them some styling
        # for points: small black circles
        l.waypoint.point <- geom_point(data=waypoint.data, aes(x=longitude, y=latitude), size=0.2, color="black", fill="white")

        # and add them to the map
        m.fancy <- m.fancy + l.waypoint.point

        # you can even add labels, but it gets a bit messy
        # styling for labels: small text, avoiding overlaps, opaque, and nudged to the right of the points 
        l.waypoint.label <-	geom_text(data=waypoint.data, aes(x=longitude, y=latitude, label=site), size=2, check_overlap=T, alpha="1", nudge_x=0.01*(abs(range.lat[1]-range.lat[2])))

        # and add them to the map
        m.fancy <- m.fancy + l.waypoint.label

        # if we want, we can remove the depth legend as follows
        m.fancy <- m.fancy + theme(legend.position="none")

        # and add a title and labels
        m.fancy <- m.fancy + ggtitle("Temperature") + xlab("Longitude") + ylab("Latitude")

        # finally, read in underway data from the Sikuliaq LDS by its URL
        data.in <- read.table("http://data.sikuliaq.alaska.edu/archive/SKQ201812T/lds/proc/20180601Z/SKQ201812T_underway_20180601Z.txt",comment="%")

        # name the columns
        data.names <- c(
            "Year",
            "Month",
            "Day",
            "Hour",
            "Minute",
            "Second",
            "ins_seapath_position Latitude",
            "ins_seapath_position Longitude",
            "ins_seapath_position Speed over ground",
            "ins_seapath_position Course over Ground",
            "ins_seapath_position Heading",
            "ins_seapath_position Roll",
            "ins_seapath_position Pitch",
            "ins_seapath_position Heading",
            "ins_seapath_position Heave",
            "mb_em302_centerbeam Latitude",
            "mb_em302_centerbeam Longitude",
            "mb_em302_centerbeam Depth",
            "met_ptu307 Atmospheric pressure",
            "met_ptu307 Air Temperature",
            "met_ptu307 Relative Humidity",
            "fluoro_turner-c6 Phycoeryth",
            "fluoro_turner-c6 CDOM",
            "fluoro_turner-c6 Chlorophyll_a",
            "fluoro_turner-c6 Crude Oil",
            "fluoro_turner-c6 Turbidity",
            "fluoro_turner-c6 Depth",
            "fluoro_turner-c6 Temp",
            "wind_gill_fwdmast_true Wind direction",
            "wind_gill_fwdmast_true Wind speed N",
            "rad_psp-pir LW",
            "rad_psp-pir SW",
            "rad_qsr2150a PAR",
            "flow_krohne_fwd flow speed",
            "flow_krohne_fwd volume flow",
            "flow_krohne_fwd coil temperature",
            "flow_krohne_fwd conductivity",
            "tsg_sbe45_fwd Temperature",
            "tsg_sbe45_fwd Conductivity",
            "tsg_sbe45_fwd Salinity",
            "tsg_sbe45_fwd Speed of sound",
            "thermo_sbe38_fwd Intake Temperature",
            "oxygen_sbe43_fwd_calc Dissolved oxygen",
            "oxygen_sbe43_fwd_calc Temperature",
            "oxygen_sbe43_fwd_calc Salinity",
            "oxygen_sbe43_fwd_calc Oxygen Solubility",
            "oxygen_sbe43_fwd_calc Dissolved Oxygen"
        )

        #make unique names with no spaces
        colnames(data.in) <- make.names(data.names)
        colnames(data.in)

        # add temperature track to fancy map
        # the way we'll do this is to plot the lat/long of our cruise track, and use colors to show the temperature
        # we'll get the lat/lon from the EM302
        # first we need to get some colors, let's do a quick and dirty version
        # since we know all of the temperature data is going to fall between 0 and 20 C, we'll use a built-in color ramp from 1:20
        twenty.colors <- heat.colors(20) #or topo.colors(), terrain.colors(), rainbow(), etc.
        simple.colors <- twenty.colors[round(data.in$tsg_sbe45_fwd.Temperature)]

        m.temperature <- m.fancy +
            geom_point(data=data.in, aes(x=mb_em302_centerbeam.Longitude, y=mb_em302_centerbeam.Latitude), color=simple.colors, size=2)

        # here's one way to make a nicer color gradient
        col.data <- data.in$tsg_sbe45_fwd.Temperature
        range.temperature <- range(col.data)
        col.row <- round(scales::rescale(col.data,to=c(1,nrow(data.in)), from=range.temperature))
        col.row[col.data < range.temperature[1]] <- 1
        col.row[col.data > range.temperature[2]] <- nrow(data.in)
        mapcolors <- rev(topo.colors(nrow(data.in)))[col.row]

        m.temperature <- m.fancy +
            geom_point(data=data.in, aes(x=mb_em302_centerbeam.Longitude, y=mb_em302_centerbeam.Latitude), color=mapcolors, size=2)

        print(m.temperature)

        # and if you want to plot the temperature against time in ggplot
        # first convert the time
#         data.in$gmt_time <- as.POSIXct(strptime(paste(data.in$Year,data.in$Month,data.in$Day,data.in$Hour,data.in$Minute,data.in$Second),format="%Y %m %d %H %M %S"))
#         data.in$local_time <- data.in$gmt_time -8*60*60
# 
#         p.temperature <- ggplot() +
#             geom_point(data=data.in, aes(x=local_time, y=tsg_sbe45_fwd.Temperature), col="grey") + #biochem lab
#             geom_point(data=data.in, aes(x=local_time, y=tsg_sbe45_fwd.Temperature), col=mapcolors) + #portside 
#             xlab("local time") + ylab("Temperature") + ylim(range.temperature)
    
print("done geographic")
    
    }, width=700, height=700)

    ###############
    # environmental plot
    ###############
    output$environmental <- renderPlot({})

    ###############
    # phylogenetic plot
    ###############
    output$phylogenetic <- renderPlot({
    
        ### tree annotated with sample heatmap
        ### PROCESS USER INPUT
        
        #### TAXA FIRST

        taxa.user.regex <- rep(FALSE,tt.cols)

        # if user doesn't provide a taxon or a regex, select all
        # print(input$taxon)
        if(input$taxon=="" & input$plot.taxa.regex=="") {
            taxa.user.regex <- rep(TRUE,tt.cols)
        # else if user provides a lone taxon and no regex, grep it
        } else if (input$taxon!="" & input$plot.taxa.regex=="") {
            taxa.user.regex <- grepl(x=newnames,pattern=input$taxon)
        } else if(input$plot.taxa.regex!="") {
        # if user provides a list of regex, split and search each one, ORing
            patterns <- strsplit(x=input$plot.taxa.regex,split="\n")[[1]]
            for(patt in patterns) {
                # print(patt)
                if(patt=="\n") next
                taxa.user.regex <- taxa.user.regex | grepl(x=newnames, pattern=patt)
            }
        } else {
            print("phylogenetic taxa: how did you get here?")
        }

        # filter by what's shown in taxa box
        taxa.in.box <- seq.tsne[,1] > ranges$network_xlim[1] & 
        seq.tsne[,1] < ranges$network_xlim[2] & 
        seq.tsne[,2] > ranges$network_ylim[1] & 
        seq.tsne[,2] < ranges$network_ylim[2]

        #select taxa in multiple samples
        taxa_mult_samples <- colSums(proptab>0) >= 2^input$min.taxa.hits 
#         picks$taxa <- picks$taxa & taxa_mult_samples #infinite loop?

        # enforce slider limits over taxa choices
        picks$taxa <- taxa.user.regex & taxa_mult_samples & taxa.in.box
        # print(paste(newnames[picks$taxa],sep="\n"))

        # print(sum(picks$taxa))
        ptl <- length(picks$taxa)

        spt <- sum(picks$taxa)
        if(spt<3) {print("phylogenetic stopped"); plot(0,0,main="not enough sequences",axes = FALSE, xlab = "", ylab = "", pch=""); return(0)} #can't do phylogeny on <3 sequences
        


        #### THEN SAMPLES
        # Sample filters
        min.reads.samples <- rs.seqtab > 2^input$min.reads
        min.richness.samples <- rowSums(proptab>0) >= 2^input$min.richness

        # initialize 
        sample.user.regex <- rep(FALSE,tt.rows)
        # print(input$sample)
        # if user doesn't provide a lone sample or a regex, select all
        if(input$sample=="" & input$plot.sample.regex=="") {
            sample.user.regex <- rep(TRUE,tt.rows)
        # if user provides a lone sample and no regex, grep for it
        } else if (input$sample!="" & input$plot.sample.regex=="") {
            sample.user.regex <- metadata$desc==input$sample
        } else if(input$plot.sample.regex!="") {
        # if user provides a list of regex, split and search each one, ORing
        # override input$sample
            patterns <- strsplit(x=input$plot.sample.regex,split="\n")[[1]]
            for(patt in patterns) {
                # print(patt)
                if(patt=="\n") next
                sample.user.regex <- sample.user.regex | grepl(x=metadata$desc, pattern=patt)
            }
        } else {
            print("phylogenetic samples: how did you get here?")
        }


        # don't waste time calculating and plotting things that aren't in the view box
        samples.in.box <- sample.tsne.rot[,1] > ranges$biospatial_xlim[1] & 
        sample.tsne.rot[,1] < ranges$biospatial_xlim[2] & 
        sample.tsne.rot[,2] > ranges$biospatial_ylim[1] & 
        sample.tsne.rot[,2] < ranges$biospatial_ylim[2]

        # put all the sample filters together
        # QC and holdup

        picks$samples <- min.reads.samples & min.richness.samples & sample.user.regex & samples.in.box
        # print(paste(metadata$desc[picks$samples],sep="\n"))
        sps <- sum(picks$samples)
        if(sps==0) {print("phylogenetic stopped"); return(0)}


        # keep it small for now
        if(spt > 100) {
            #want selected taxa in the top 100 in abundance for the selected 
            seqs.taxa <- rep(FALSE,length(picks$taxa)) # F F F F
            order.taxa <- order(colSums(proptab[picks$samples,picks$taxa]),decreasing=T) # 3 1 2 4
            seqs.taxa[picks$taxa] <- order.taxa #probably FALSE got converted to zeros
          
            seqs.taxa <- seqs.taxa < 100 & seqs.taxa > 0
#             order.taxa.pick <- order.taxa[picks$taxa] # 1 4 
#             order.taxa.pick.order <- order(order.taxa.pick)  # 1 2 
#             order.taxa.pick.order.which <- which(order.taxa.pick.order <=100) # 1 2
#             order.taxa.pick.order.which.which <- which(
#             order.max <- max(order.taxa[picks$taxa])         
#             seqs.taxa[order.taxa <= order.max & picks$taxa] <- TRUE
#             seqs.taxa[sample(which(picks$taxa),100)] <- TRUE
        } else {
            seqs.taxa <- picks$taxa
        }
#         sub.physeq <- subset_taxa(seqs.physeq, seqs.taxa)
#         sub.physeq <- subset_samples(sub.physeq, picks$samples)
#        sub.bin <- seqs.bin[seqs.taxa]
#        checkAlignment(sub.bin)

        seqs.ss <- DNAStringSet(colnames(proptab))
        sub.ss <- seqs.ss[seqs.taxa,]
        # align sequences / will fail miserably with ITS
#        seqs.aln <- muscle(sub.bin, quiet=FALSE, MoreArgs="-maxiters 1 -diags")
		sub.aln <- msa(sub.ss,type="dna","ClustalOmega",order="input")
        sub.aln.c <- as.character(sub.aln@unmasked)
        sub.aln.raw <- t(sapply(strsplit(sub.aln.c,""), tolower))
        rownames(sub.aln.raw) <- colnames(proptab)[seqs.taxa]
        sub.bin <- as.DNAbin(sub.aln.raw)
        # get distances
        sub.dist <- dist.dna(sub.bin,as.matrix=T)
        sub.dist[is.na(sub.dist)] <- 1
        sub.dist[is.infinite(sub.dist)] <- 1
        sub.dist[is.nan(sub.dist)] <- 1
        # make tree
        sub.tree <- bionjs(sub.dist)
        sub.tree <- midpoint(sub.tree)
        
        #assign taxonomic labels
        user.levels <- input$tax.levels
        print(user.levels)
        if(length(user.levels)==0) { 
            sub.tree$tip.label <- paste0("ESV_",which(seqs.taxa)) 
        } else {
            user.taxlevels <- unname(taxout$ref_dada2_silva_nr_v132_train_set.fasta$tax[colnames(proptab)[seqs.taxa],user.levels])
            user.taxlevels <- cbind(user.taxlevels,paste0("ESV_",which(seqs.taxa)))
            sub.tree$tip.label <- apply(user.taxlevels,1,paste,collapse=";")
        }
#         sub.physeq <- subset_taxa(seqs.physeq, seqs.taxa)
#         sub.physeq <- subset_samples(sub.physeq, picks$samples)
#         sub.otu <- otu_table(t(proptab[,seqs.taxa]), taxa_are_rows = TRUE)
#         sub.tax <- tax_table(taxout$ref_dada2_silva_nr_v132_train_set.fasta$tax[colnames(proptab)[seqs.taxa],])
#         sub.physeq <- phyloseq(sub.otu, sub.tax, sub.tree)
    
        #scale to relative abundances
        if(sps>1) {
            sub.edge.weight <- unname(colSums(proptab[picks$samples,seqs.taxa,drop=F]))
        } else {
            sub.edge.weight <- unname(proptab[picks$samples,seqs.taxa,drop=F])
        }

        # now need to calculate for internal nodes
        sub.node.weight <- rep(0,nrow(sub.tree$edge)+1)
        for(tip in 1:length(sub.tree$tip.label)) {
            tip.anc <- c(tip,Ancestors(sub.tree,type="all",node=tip))
#             sub.node.weight[tip] <- sub.edge.weight[tip]
#             sub.node.weight[tip.anc] <- sub.node.weight[tip.anc] + sub.edge.weight[tip]
            for(t in 2:length(tip.anc)) {
                tt <- which(sub.tree$edge[,1]==tip.anc[t] & sub.tree$edge[,2]==tip.anc[t-1])
#                 print(tt)
                sub.node.weight[tt] <- sub.node.weight[tt] + sub.edge.weight[tip]
            }
        }
        sub.node.weight <- rescale(sqrt(sqrt(sub.node.weight)),to=c(1,10))
#        sub.node.weight <- sub.node.weight[match(1:length(sub.node.weight),sub.tree$edge[,1])]
        # new distance matrix from branch lengths
        sub.dist.tree <- cophenetic.phylo(sub.tree)
        # jitter to avoid Rtsne complaining
        sub.dist.jit <- apply(sub.dist.tree,1,jitter)
        # sometimes Rtsne fails "perplexity too high"
        tsne.complete <- 0
        perp <- round(ncol(sub.dist.jit)/3)
        counter <- 0
        while(tsne.complete==0) {
            counter <- counter+1
            if(counter>100) return(NULL)
            try({
                perp <- perp-1
               set.seed(42) #I don't actually understand why this works to prevent infinite loop
                sub.tsne <- unname(metaMDS(comm=sub.dist, k=3, distance="euclidean")$points)
#                sub.tsne <- Rtsne(X=sub.dist.jit,dims=3,perplexity=perp)$Y
                tsne.complete <- 1
            })
        }

#         print(color$phylo[seqs.taxa])

        # if colors for any member of this set of taxa are not set, and no rotation is requested, make default colors
        if(any(is.na(color$phylo[seqs.taxa])) | all(input$tree.colors==0)) {
            sub.colors <- cbind(sub.tsne[,1],sub.tsne[,2],sub.tsne[,3])
            sub.colors.rot <- rgb(rescale(sub.colors,to=c(0,0.8)))

#             print("default colors")

        # else use rotated color
        } else if (any(input$tree.colors!=0)) {
            sub.colors <- cbind(sub.tsne[,1],sub.tsne[,2],sub.tsne[,3])
            sub.colors.rot <- cbind(Rotation(sub.colors[,1:2],input$tree.colors[1]),Rotation(sub.colors[,c(2,3)],input$tree.colors[2]))
            sub.colors.rot <- rgb(rescale(sub.colors.rot,to=c(0,0.8)))
#             print("new colors")
        }

        #kill infinite loop
        isolate({color$phylo[seqs.taxa] <- sub.colors.rot})

        # plot tree, auto-calculate scaling only if first time
        if(is.null(last$seqs) | !identical(last$seqs,picks$taxa)) { #if defaults
            sub.plot <- plot(sub.tree, no.margin=TRUE, plot=FALSE, edge.width=1, label.offset=input$tree.width[2]/50000)
            updateSliderInput(session, "tree.height", value = sub.plot$y.lim*10)
            updateSliderInput(session, "tree.width", value = sub.plot$x.lim*1000)
        }
       
       plot(sub.tree, no.margin=TRUE, x.lim=input$tree.width/1000,y.lim=input$tree.height/10, cex=input$tree.font.size/10, tip.color=sub.colors.rot, edge.width=sub.node.weight, show.node.label = FALSE, label.offset=input$tree.width[2]/50000)
        #tiplabels() for relative abundance data
       isolate({last$seqs <- picks$taxa})

        print("done phylogenetic")
    }, width=700, height=900)


    ###############
    # functional plot
    ###############
    output$functional <- renderPlot({})

    ###############
    # taxa table
    ###############
      output$rhtable <- renderTable({
        if(input$tabs == "biospatial") {
            #we want a table of selected taxa in the sample being hovered over

            xy <- xy_biospatial(input$biospatial_hover)
            if(is.null(xy)) return()
            # print(xy)
            if(picks$samples[xy]) {
                sample.count <- unname(seqtab[xy,picks$taxa,drop=F])
                sample.sum <- sum(unname(seqtab[xy,,drop=F]))
                rhtable.num <- round(1000*sample.count/sample.sum,1)[order(sample.count,decreasing=T)] #bc of sprintf
                rhtable.ppt <- sprintf('%.2f',rhtable.num)
                rhtable.taxa <- newnames[picks$taxa][order(sample.count,decreasing=T)]
                  data.frame(
                    cbind(
                        "per mille"=rhtable.ppt[rhtable.num>0],                  
                        "taxa"=rhtable.taxa[rhtable.num>0]
                        )
                    )
            }
        } else if(input$tabs == "network") {
            #we want a table of networked taxa connected to the taxa being hovered over
            # or do we want a table of samples that the hovered taxa is in
            xy <- xy_network(input$network_hover) #gives closest point
            if(is.null(xy)) return()
            if(picks$taxa[xy]) {
                taxa.count <- unname(proptab[picks$samples,xy,drop=F]) #raw reads of taxon in selected samples
                taxa.sums <- rowSums(unname(proptab[picks$samples,,drop=F])) #raw reads of all taxa in selected samples
                rhtable.num <- round(taxa.count/1000,1)[order(taxa.count,decreasing=T)] #bc of sprintf
                rhtable.ppt <- sprintf('%.2f',rhtable.num)
                rhtable.taxa <- metadata$desc[picks$samples][order(taxa.count,decreasing=T)]
                  data.frame(
                    cbind(
                        "per mille"=rhtable.ppt[rhtable.num>0],                  
                        "samples"=rhtable.taxa[rhtable.num>0]
                        )
                    )
            }
        } else {
            print("nothing to see here")
            data.frame()
        }
      }, align='rl',digits=2)



    ###############
    # overunder table
    ###############

      output$overunder <- renderTable({
        if(input$tabs == "biospatial") {
            #we want a table of over and underrepresented taxa in the sample being hovered over

            xy <- xy_biospatial(input$biospatial_hover)
            # if not hovering, then table is for whole viewport
            if(is.null(xy)) {
                
            
            } else {
            # if hovering then table is just for that sample
                if(picks$samples[xy]) {
                    sample.count <- unname(seqtab[xy,picks$taxa,drop=F])
                    sample.sum <- sum(unname(seqtab[xy,,drop=F]))
                    overunder.num <- round(1000*sample.count/sample.sum,1)[order(sample.count,decreasing=T)] #bc of sprintf
                    overunder.ppt <- sprintf('%.2f',overunder.num)
                    overunder.taxa <- newnames[picks$taxa][order(sample.count,decreasing=T)]
                      data.frame(
                        cbind(
                            "per mille"=overunder.ppt[overunder.num>0],                  
                            "taxa"=overunder.taxa[overunder.num>0]
                            )
                        )
                }
            }
        } else if(input$tabs == "network") {
            # we want a table of samples which have enriched abundances of the taxa in the FOV
            xy <- xy_network(input$network_hover) #gives closest point
            if(is.null(xy)) return()
            if(picks$taxa[xy]) {
                taxa.count <- unname(proptab[picks$samples,xy,drop=F]) #raw reads of taxon in selected samples
                taxa.sums <- rowSums(unname(proptab[picks$samples,,drop=F])) #raw reads of all taxa in selected samples
                overunder.num <- round(taxa.count/1000,1)[order(taxa.count,decreasing=T)] #bc of sprintf
                overunder.ppt <- sprintf('%.2f',overunder.num)
                overunder.taxa <- metadata$desc[picks$samples][order(taxa.count,decreasing=T)]
                  data.frame(
                    cbind(
                        "per mille"=overunder.ppt[overunder.num>0],                  
                        "samples"=overunder.taxa[overunder.num>0]
                        )
                    )
            }
        } else {
            print("nothing to see here")
            data.frame()
        }
      }, align='rl',digits=2)
    
    


}

# Create Shiny app ----
shinyApp(ui = ui, server = server)



#network
# do this once outside of shiny app and import it


# can this be vastly simplified by just saving vectors of weights and colors instead of the whole graph?
#network_correlations <- props.sparcc$Cor #does it make sense to try to choose networked taxa here?
#save.graphs <- list()
#save.weights <- list()

# for(i in 100:1) {
#     cutoff <- i/100
#     print(cutoff)
#     # define graph
#     network_matrix <- abs(network_correlations) > cutoff
#     diag(network_matrix) <- 0
#     if(sum(network_matrix)==0) { save.graphs[[i]] <- NULL; save.weights[[i]] <- NULL; next }
#     network_matrix <- as.matrix(tril(network_matrix) + t(tril(network_matrix))) #to make it symmetricso igraph doesn't crash
# #    network_matrix <- Matrix(network_matrix, sparse=TRUE) #does this actually speed anything up?
#     network_graph <- adj2igraph(network_matrix)
#     #define weights and colors
#     network_weights <- (network_matrix & upper.tri(network_matrix))*network_correlations
#     #scale the weights to x^4
#     E(network_graph)$weight <- network_weights[abs(network_weights) > cutoff]^3 #preserve sign
#     #recenter the weights
#     E(network_graph)$weight <- E(network_graph)$weight/max(abs(E(network_graph)$weight))
#     #color the rescaled weights
#     E(network_graph)$color <- rgb(0,0,0,0)
#     network_neg <- E(network_graph)$weight < 0
#     network_pos <- E(network_graph)$weight > 0
#     E(network_graph)$color[network_neg] <- rgb(abs(E(network_graph)$weight[network_neg]),0,1,abs(E(network_graph)$weight[network_neg]))
#     E(network_graph)$color[network_pos] <- rgb(0,0,abs(E(network_graph)$weight[network_pos]),abs(E(network_graph)$weight[network_pos]))
#     network_weights[network_weights==1] <- 0
#     save.graphs[[i]] <- network_graph
#     save.weights[[i]] <- network_weights
# }
# save(file="save.graphs.Rdata",save.graphs) #730Mb each... must be a better way but for now...
# save(file="save.weights.Rdata",save.weights)



#         print(input$point.scalar)
#         print(point.size$network)
#         if(input$point.scalar != 0 & input$point.scalar != point.size$network) {
#             point.size$network <- point.size$network + input$point.scalar
#             updateSliderInput(session, "point.scalar", value = 0)
#         }

#  old
# 
#         min.reads.samples <- rs.seqtab > 2^input$min.reads
#         min.richness.samples <- rowSums(proptab>0) > 2^input$min.richness
# 
#         sample.user.regex <- metadata$desc==input$sample
# #        sample.user.regex <- rep(FALSE,tt.rows)
# 
#         if(input$sample=="") {
#             sample.user.regex <- rep(TRUE,tt.rows)
#         } else {
#         sample.user.regex <- metadata$desc==input$sample
# 
#             if(input$plot.sample.regex != "") {
#                 patterns <- strsplit(x=input$plot.sample.regex,split="\n")[[1]]
#                 for(patt in patterns) {
#                     # print(patt)
#                     if(patt=="\n") next
#                     sample.user.regex <- sample.user.regex | grepl(x=metadata$desc, pattern=patt)
#                 }
#             }
#         }
# 
# 
#         picks$samples <- min.reads.samples & min.richness.samples & sample.user.regex
#         sps <- sum(picks$samples)
#         
#         

#
# old trash if works 
#     #choose networked taxa to plot
#     taxon <- input$taxon
#     
#     picks$taxa <- rep(FALSE,tt.cols)
#     taxa.user.regex <- rep(FALSE,tt.cols)
#     if(input$plot.taxa.regex != "") {
#         patterns <- strsplit(x=input$plot.taxa.regex,split="\n")[[1]]
#         for(patt in patterns) {
#             # print(patt)
#             if(patt=="\n") next
#             taxa.user.regex <- taxa.user.regex | grepl(x=newnames, pattern=patt)
#         }
#         picks$taxa <- taxa.user.regex #logical
#     } else {
#         picks$taxa <- grepl(pattern=taxon,x=newnames,fixed=T) #logical
#     }
#     


# can't figure out how to not make it scale both graphs
#         print(input$point.scalar)
#         print(input$tabs)
#         if(input$point.scalar != 0 & input$point.scalar != point.size$biospatial) {
#             point.size$biospatial <- point.size$biospatial + input$point.scalar
#             updateSliderInput(session, "point.scalar", value = 0)
#         }

# forget autoscaling again
# 
#         print(networked_taxa)
#         if(is.null(ranges$network_xlim)) {
# #             xlim <- c(-60,60)
#             print(networked_taxa)
#             print(seq.tsne[networked_taxa,])
#             xlim2 <- range(seq.tsne[networked_taxa,1,drop=F],na.rm=T,finite=T)
# #             ox <- mean(seq.tsne[networked_taxa,1,drop=F],na.rm=T,finite=T)
# #             xd <- max(abs(xlim2-ox))
# #             zd <- c(xd,yd)[which.max(abs(c(xd,yd)))]
# #             xlim2[1] <- ox - xd*1.1
# #             xlim2[2] <- ox + xd*1.1
# #             if(xlim2[1] < xlim[1]) xlim2[1] <- xlim[1]
# #             if(xlim2[2] > xlim[2]) xlim2[2] <- xlim[2]
#             ranges$network_xlim <- xlim2
#         }
# 
#         if(is.null(ranges$network_ylim)) {
# #            ylim <- c(-60,60)
#             ylim2 <- range(seq.tsne[networked_taxa,2,drop=F])
# #             oy <- mean(seq.tsne[networked_taxa,2,drop=F])
# #             yd <- max(abs(ylim2-oy))
# #             zd <- c(xd,yd)[which.max(abs(c(xd,yd)))]
# #             ylim2[1] <- oy - yd*1.1
# #             ylim2[2] <- oy + yd*1.1
# #             if(ylim2[1] < ylim[1]) ylim2[1] <- ylim[1]
# #             if(ylim2[2] > ylim[2]) ylim2[2] <- ylim[2]
#             ranges$network_ylim <- ylim2
#         }
# 
#         print(paste0("ranges:", ranges$network_xlim))
# 
# 
# library(shiny)
# runApp( list(ui = bootstrapPage(
#   verbatimTextOutput("results"),
#   tags$script('
#     $(document).on("keypress", function (e) {
#        Shiny.onInputChange("mydata", e.which);
#     });
#   ') 
# )
# , server = function(input, output, session) {
# 
#   output$results = renderPrint({
#     input$mydata
#   })
# }
# ))
# for keydown events you can substitute:
# 
#   tags$script('
#     $(document).on("keydown", function (e) {
#        Shiny.onInputChange("mydata", e.which);
#     });
#   ') 

