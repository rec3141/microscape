# Microscape Shiny App — Feature TODO

## Sample Explorer
- [x] t-SNE scatter plot with plotly WebGL
- [x] Point size proportional to read depth
- [x] Color by metadata field
- [ ] **Taxa composition overlay** — colored circles per-taxon on each sample, size = relative abundance (the signature "microscape" view)
- [ ] Single sample dropdown selector
- [ ] Sample regex filter
- [ ] Min reads slider filter
- [ ] Min richness slider filter
- [ ] Primer/target pair filter (16S_prokaryote, 18S_protist, ITS_all)
- [ ] Taxonomy group filter (affects which taxa are shown in overlay)
- [ ] Taxonomy regex filter
- [ ] Click sample → show top taxa table in sidebar
- [ ] Hover → show sample name + metadata
- [ ] Random name labels toggle

## ASV Network
- [x] t-SNE scatter plot with plotly WebGL
- [x] Points colored by taxonomic group
- [ ] SparCC correlation edges between ASVs
- [ ] Correlation threshold slider
- [ ] Min prevalence (samples per ASV) slider
- [ ] Taxonomy regex filter
- [ ] Primer/target filter
- [ ] Color by detailed taxonomy (Alphaproteo=blue, Bacteroidetes=teal, Archaea=yellow, etc.)
- [ ] Point size = mean abundance across selected samples
- [ ] When sample selected in Sample Explorer, highlight that sample's ASVs
- [ ] Hover → full taxonomy string
- [ ] Click ASV → show samples containing it in sidebar

## Sidebar / Cross-tab
- [ ] Taxonomy database selector (switch between SILVA, PR2, etc.)
- [ ] Primer/target group checkboxes (shared across tabs)
- [ ] Selected sample display with top taxa table
- [ ] Selected ASV display with sample list
- [ ] Bootstrap confidence display

## Summary Tables
- [x] Searchable/sortable ASV info table
- [x] Searchable/sortable sample info table
- [ ] CSV download buttons

## Pipeline Features
- [ ] SRA import (download FASTQ from accession numbers)
- [ ] MIMARKS template validation on metadata import
- [ ] Phylogenetic tree tab (using pre-built tree from pipeline)
- [ ] Geographic map tab (when lat/lon in metadata)

## Performance
- [x] plotly WebGL for 100K+ points
- [x] data.table for all filtering
- [x] Pre-computed data bundle (app_data.rds)
- [ ] Debounced slider inputs
- [ ] Lazy-load network edges (only when tab is active)
