#!/usr/bin/env Rscript
# hello world!
# Load required libraries
library(Gviz)
library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

# load biomaRt object for ch37
# note: doesn't seem to work offline
# load("scripts/biomart_grch37.RData")


# load ucsc cytoband for hg19
cytoBands <- readRDS("scripts/hg19_cytobands.rds")

# load txdb
txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Get command line arguments (optional)
args <- commandArgs(trailingOnly = TRUE)
symbol <- args[1] # First argument as gene of interest
input_file <- args[2]  # Second argument as input file
output_file <- args[3]  # Third argument as output file
#bigwig_files <- args[4] # Fourth argument as a list of bw files

files_arg <- args[which(grepl("--files=", args))]

bigwig_files <- strsplit(gsub("--files=", "", files_arg), ",")[[1]]
bigwig_files <- as.list(bigwig_files)
bigwig_files

# If no arguments provided, use defaults
if(length(args) == 0) {
  symbol <- "CA9"
  input_file <- "input_data.csv"
  output_file <- "results.csv"
}

# Print some information
print(paste("Starting analysis at", Sys.time()))
print(paste("Gene of interest:", symbol))
print(paste("Input file:", input_file))
print(paste("Output file:", output_file))
print(paste("bigwig files:", bigwig_files))

# ---------- #
#  Functions #
# ---------- #


# Function to Extract coordinates from a gene
gene_coord <- function(symbol, Txdb = TxDb.Hsapiens.UCSC.hg19.knownGene){
  symbol <- toupper(symbol) # make gene symbol upper case
  
  # find the entrez id for that symbol
  entrez_id <- select(org.Hs.eg.db, keys=symbol, keytype="SYMBOL", columns="ENTREZID")
  
  # find the coordinates of that gene based on entrez id
  coords <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene, filter = list(gene_id=entrez_id$ENTREZID))
  
  print(paste("The ENTREZID for the following gene:",symbol,"is:",entrez_id$ENTREZID))
  
  print(paste(
      "The coordinates for the following gene:",
      symbol,
      "is:",
      as.character(seqnames(coords)),start(coords),"-",end(coords)
  ))
  return(coords)
}

# Function to get max signal from one or more bigWig files in a given range
getMaxSignal <- function(bigwig_files, chr, start, end) {
  # Handle single file input
  if (is.character(bigwig_files) && length(bigwig_files) == 1) {
    bigwig_files <- list(bigwig_files)
  }
  
  # Get maximum signal across all files
  max_signals <- sapply(bigwig_files, function(bw) {
    # Import specific genomic range from bigWig
    data <- import(bw,
                   format = "bigWig",
                   which = GRanges(
                     seqnames = chr,
                     ranges = IRanges(start = start,
                                      end = end)
                   ))
    
    if (length(data) == 0)
      return(0)
    return(max(data$score))
  })
  
  return(max(max_signals))
}

# function that uses TxDb to make a GeneRegionTrack with proper Gene Symbol display
MakeGeneRegionTrack_txdb <- function(txdb, chr, start, end, name = "gene"){
  gr <- GeneRegionTrack(txdb,
                        chromosome = chr,
                        start = start,
                        end = end,
                        name = name)
  # Map transcript IDs to gene symbols (TxDb doesn't store gene symbols)
  symbols <-
    unlist(mapIds(org.Hs.eg.db, gene(gr), "SYMBOL", "ENTREZID", multiVals = "first"))
  symbol(gr) <- symbols[gene(gr)]
  return(gr)
}

# Function to Create and plot tracks with dynamic y-axis
createScaledTracks <- function(bigwig_files,
                               chr,
                               start,
                               end,
                               track_names = NULL,
                               colors = NULL,
                               padding_factor = 1.1,
                               overlay = "no") {
  # Handle input parameters
  if (is.character(bigwig_files) && length(bigwig_files) == 1) {
    bigwig_files <- list(bigwig_files)
  }
  
  # Generate default track names if not provided
  if (is.null(track_names)) {
    track_names <- paste("Signal", seq_along(bigwig_files))
  } else if (length(track_names) != length(bigwig_files)) {
    stop("Number of track names must match number of bigWig files")
  }
  
  # Generate default colors if not provided
  if (is.null(colors)) {
    colors <- rainbow(length(bigwig_files))
  } else if (length(colors) != length(bigwig_files)) {
    stop("Number of colors must match number of bigWig files")
  }
  
  # Get maximum signal across all files
  max_signal <- getMaxSignal(bigwig_files, chr, start, end)
  
  # Create DataTrack for each bigWig file
  tracks <- list()
  for (i in seq_along(bigwig_files)) {
    tracks[[i]] <- DataTrack(
      range = bigwig_files[[i]],
      chromosome = chr,
      from = start,
      to = end,
      name = track_names[i],
      type = "l",
      col = colors[i],
      # Explicitly assign track color by index
      baseline = 0,
      # Ensure histogram starts at 0
      col.baseline = "black",
      # set baseline color
      groups = factor(track_names[i],
                      levels = track_names),
     legend = TRUE,
      ylim = c(0, max_signal * padding_factor)
    )
  }
  
  # return scaled tracks:
  # option to overlay tracks
  if (tolower(overlay) == "yes") {
    tracks[[1]]@name <- "signal" # set shared y-axis title
    overlayTrack <- OverlayTrack(trackList = tracks)
    return(overlayTrack)
  } else{
    # default to stack tracks
    return(tracks)
  }
  
}

# ------------ #
#  Plotting    #
# ------------ #

# extract coordinates for a given gene symbol
coords <- gene_coord(toupper(symbol))

# assign coordinate info from symbol
chr <- as.character(seqnames(coords))
start <- start(coords)
end <- end(coords)

# add genome coordinates track
print("Create GenomeAxisTrack")
gtrack <- GenomeAxisTrack()

# add ideogram track (visual representation of chromosome)
print("Create IdeogramTrack")
itrack <- IdeogramTrack(genome = "hg19", bands = cytoBands, chromosome = chr)

# Custom track names and colors
track_names <- c("normoxia", "hypoxia")
track_colors <- c("red", "blue")

print("Create scaledtracks")
# create scaled tracks
scaledTracks <- createScaledTracks(
  bigwig_files = bigwig_files,
  chr = chr,
  start = start,
  end = end,
  track_names = track_names,
  colors = track_colors,
  padding_factor = 1.2  # 20% padding above max signal
)

scaledTracks
# create Gene region track
txTr <- MakeGeneRegionTrack_txdb(txdb_hg19, chr = chr, start = start, end = end)

# biomTrack <-
#   BiomartGeneRegionTrack(
#     genome = "hg19",
#     name = "ENSEMBL",
#     symbol = symbol,
#     biomart = bm,
#     col = "transparent",
#     filter = list(with_refseq_mrna = TRUE),
#     # filtering to RefSeq models only
#     collapseTranscripts = "longest"
#   )

# single data track

bgFile <- system.file("extdata/test.bedGraph", package = "Gviz")
bamFile <- system.file("extdata/test.bam", package = "Gviz")

dTrack.ATAC1 <- DataTrack(range = bamFile, genome = "hg19", type = "l", 
                     chromosome = "chr1", window = -1,
                     name = "bedGraph")

dTrack4 <- DataTrack(range = bgFile, genome = "hg19", type = "l", 
                     name = "Coverage", window = -1, 
                     chromosome = "chr1", stream = TRUE)

print("Attempt to Plot Tracks")
# plot results
#png(paste0("results/","test",".png"))
pdf(paste0("results/","test.pdf"))
 plotTracks(c(list(itrack,
                  gtrack),
             scaledTracks,
             txTr),
           from = start,
           to = end,
           col = c("red", "blue"),
           transcriptAnnotation = "symbol",
           stream = TRUE)

plotTracks(dTrack.ATAC1,  from = 189990000, to = 190000000, background.panel = "#FFFEDB")
displayPars(dTrack4) <- list(fontcolour = "black")
plotTracks(dTrack4,  from = 189990000, to = 190000000, background.panel = "#FFFEDB")

#plotTracks(scaledTracks, from = start, to = end, transcriptAnnotation = "symbol")
#plotTracks(dTrack, from = start, to = end, transcriptAnnoation = "symbol")
#plotTracks(dTrack.ATAC1, chromosome = "chr6", from = start, to = end, transcriptAnnoation = "symbol")

dev.off()
# save results as PDF

# Example script use:
# Rscript scripts/genome-viewer.R EGFR --files="/cluster/projects/wouterslab/RW555/deeptools/heatmap/RW555-1_treat_pileup.bw,/cluster/projects/wouterslab/RW555/deeptools/heatmap/RW555-2_treat_pileup.bw"

# show displayPars
displayPars(dTrack.ATAC1)
print("next")

