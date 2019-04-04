# Microscape: An Interactive Web Application for Microbiome Analysis
Visualizing Microbial Landscapes

You can try out the R Shiny app at https://cryomics.shinyapps.io/microscape/

If you get a 503 Error try hitting refresh a couple of times

But it will be much faster if you install it yourself locally, and you can use your own data! It is not currently nicely packaged, you have to just download the scripts from github.

This program uses DADA2 as the workhorse: https://benjjneb.github.io/dada2/

The basic output of the DADA2 pipeline is 2 files:

1) seqtab: a count table with samples in the rows and ESV sequences in the columns. The columns are named by the complete sequence.
2) taxtab: a taxonomic classification table with taxonomic levels in the columns and ESV sequences in the rows (again named by sequence).


These are all quality controls steps that I use to go from the MiSeq to R, and may not be applicable to anyone else's setup
- microscape-s00-demultiplex.sh
- microscape-s00-remove-primers.sh
- microscape-s00-countseqs.sh
- microscape-s00-plot-seq-map.R

microscape-s00-setup.r
- This provides some functions to be used later

microscape-s01-setup-tax-dbs.r
- This ought to set up the taxonomic databases I use, but is not quite automated at the moment, but I can provide these to anyone who needs them

microscape-s10-dada.r
- This runs dada2 on each sequencing plate/run + primer combination and produces a sequence table for each one because plates need to be run separately due to the error profiles being different for each one -- DADA2 uses those to correct errors. 

microscape-s20-merge.r
- This merges the sequence tables and and does some quality control including chimera removal

microscape-s25-phylogeny.r
- This optional step is still in development, but it uses multiple sequence alignment and phylogenetics to correct additional errors, for example ESVs that differ only by an overhang due to bad trimming or other artifacts. It also outputs some diversity plots

microscape-s30-taxonomy.r
- This does the taxonomic classification using the available databases

microscape-s40-metadata.r
- This loads a metadata file, not currently used except for renaming samples, but I have plans to provide more plotting capabilities in the future using phyloseq

microscape-s50-renormalize.r
- I don't believe in throwing away perfectly good data, so this uses the taxonomic classification to separate out ESVs based on primer set + taxonomic affiliation into separate tables for:
  - 16S/bacteria+archaea
  - 16S/chloroplasts
  - 16S/mitochondria  
  - 16S/eukaryotes (e.g. mispriming)
  - 18S/metazoa
  - 18S/protists
  - ITS/all

- It produces a master count table with all of these put together, and another one where they are normalized to a counts per million sequences (for input into next steps).

microscape-s60-clusters.r
- This does the clustering to make two t-SNE maps, one where each point is a sample, and one where each point is an ESV

microscape-s70-network.r
- This runs the SparCC network correlation analysis

microscape-s99-shiny.r
- This produces the shiny app that uses outputs from the previous steps to make the interactive plots
