exonToCDS <- function(exons, cdsStart, cdsEnd) {
  if (is.na(cdsStart) || is.na(cdsEnd) || cdsEnd <= cdsStart) {
    return(NULL); } 
  cds <- exons[(
    (GenomicRanges::start(exons) < cdsEnd) &
      (GenomicRanges::end(exons) > cdsStart))]
  GenomicRanges::start(cds) <- pmax(
    GenomicRanges::start(cds), cdsStart)
  GenomicRanges::end(cds) <- pmin(
    GenomicRanges::end(cds), cdsEnd)
  cds
}

plotTranscript1 <- 
  function(txn, reg, i, beg, end,
           isoformHeight, padHeight, txn.font.size) {
    
    txn_name <- names(txn)[1]; exons <- txn[[1]]
    meta <- as.data.frame(GenomicRanges::mcols(txn))
    plt.width <- end - beg
    txn.beg <- max(beg, min(GenomicRanges::start(exons)))
    txn.end <- min(end, max(GenomicRanges::end(exons)))
    exons <- subsetByOverlaps(exons, reg)
    txn.strand <- as.character(GenomicRanges::strand(exons[1]))
    lined <- (c(txn.beg, txn.end)-beg) / plt.width # direction is in arrow ends
    
    y.bot <- (i-1) * isoformHeight + padHeight
    y.bot.exon <- y.bot + padHeight
    y.hei <- isoformHeight - 2 * padHeight
    
    ## transcript name
    g <- gList(grid.text(sprintf('%s (%s) [%s : %s]', meta$gene_name, txn_name, meta$transcript_name, meta$transcript_type),
                         x=mean(lined), y=y.bot + y.hei + padHeight * 0,
                         just=c('center','bottom'),
                         gp = gpar(fontfamily="sans", fontsize = txn.font.size), draw=FALSE))
    
    ## plot transcript line
    g <- gList(g, gList(grid.lines(
      x=lined, y=y.bot+y.hei/2, arrow=arrow(length=unit(0.06, "inches"),
                                            ends=ifelse(txn.strand == "+", "last", "first")), draw=FALSE)))
    
    g <- gList(g, gList(grid.lines(x=c(0,1), y=y.bot+y.hei/2,
                                   gp=gpar(lty='dotted'), draw=FALSE)))
    
    ## plot exons
    g <- gList(g, gList(
      grid.rect((GenomicRanges::start(exons)-beg)/plt.width,
                y.bot + y.hei/2 - y.hei/3, GenomicRanges::width(exons)/plt.width,
                y.hei/3*2, gp=gpar(fill='grey30', lwd=0),
                just=c('left','bottom'), draw=FALSE)
    ))
    
    #plot TSS (transcription start site aka first base of 1st exon for a transcript)
    g <- gList(g, gList(grid.segments(
      x0 = ifelse(txn.strand == "+", ((GenomicRanges::start(exons)[1]-beg)/plt.width), ((GenomicRanges::start(exons)[1]+(GenomicRanges::width(exons)[1]-1))-beg)/plt.width), 
      y0=(y.bot + y.hei/2 - y.hei/3), 
      x1= ifelse(txn.strand == "+", ((GenomicRanges::start(exons)[1]-beg)/plt.width), ((GenomicRanges::start(exons)[1]+(GenomicRanges::width(exons)[1]-1))-beg)/plt.width),
      y1=((y.bot + y.hei/2 - y.hei/3)+y.hei/3*2),
      gp = gpar(col="magenta", lwd=2),
      draw=FALSE)
    ))
    
    
    ## plot cds
    cds <- exonToCDS(exons, as.integer(meta$cdsStart), as.integer(meta$cdsEnd))
    if (length(cds) > 0) {
      g <- gList(g, gList(
        grid.rect((GenomicRanges::start(cds)-beg)/plt.width,
                  y.bot + y.hei/2 - y.hei/6, GenomicRanges::width(cds)/plt.width,
                  y.hei/6*2, gp=gpar(fill='red', lwd=0),
                  just=c('left','bottom'), draw=FALSE)))
    }
    g
  }

plotTranscripts <- function(
    txns, reg, beg, end,
    txn.types = NULL, txn.font.size = txn.font.size) {
  
  if (is.null(txn.types)) {
    txn.types <- unique(GenomicRanges::mcols(txns)$transcript_type)
    txns <- txns[
      GenomicRanges::mcols(txns)$transcript_type %in% txn.types]
  }
  
  if (!is.null(txn.types)) {
    txns <- txns[
      GenomicRanges::mcols(txns)$transcript_type %in% txn.types] 
  }
  
  if (length(txns) == 0) {
    return(gList(
      grid.rect(0,0.1,1,0.8, just = c('left','bottom'), draw=FALSE),
      grid.text('No transcript found', x=0.5, y=0.5, gp = gpar(fontfamily="sans", fontsize = txn.font.size), draw=FALSE)))
  }
  
  isoformHeight <- 1/length(txns)
  padHeight <- isoformHeight*0.2
  
  do.call(gList, lapply(seq_along(txns), function(i) {
    plotTranscript1(txns[i], reg, i, beg, end,
                    isoformHeight, padHeight, txn.font.size)
  }))
  
}





plotMapLines <- function(probes, beg, end) {
  nprobes <- length(probes)
    x00 <- ((GenomicRanges::start(probes) - beg) / (end - beg))
    y0 <- rep(0.5, length.out=length(probes))
    x1 <- ((seq_len(nprobes) - 0.5)/nprobes)
    y1 <- rep(0, length.out=nprobes)
    x0 <- c(x00, x00)
    x1 <- c(x1, x00)
    y0 <- c(y0, rep(0.5, length.out=length(probes)))
    y1 <- c(y1, rep(1, length.out=length(probes)))
    grid.segments(x0, y0, x1, y1, draw=FALSE)
}

##to add lines that intersect at the CpG Loci coordinate on the all of the transcripts
plotLociTxnIntersectionLines <- function(probes, beg, end, txns) {
  nprobes <- length(probes)
  
#making height of 1 txn section
  i=as.numeric(length(seq_along(txns)))
  
  isoformHeight <- 1/length(txns)
  padHeight <- isoformHeight*0.2
  
  y.bot <- (i-1) * isoformHeight + padHeight
  y.bot.exon <- y.bot + padHeight
  y.hei <- isoformHeight - 2 * padHeight
  
  one_txn_section_height =  as.numeric(isoformHeight*i)
  
  grid.segments(x0 = ((GenomicRanges::start(probes) - beg) / (end - beg)), y0 =  0,
                x1 = ((GenomicRanges::start(probes) - beg) / (end - beg)), y1 =  as.numeric(one_txn_section_height*2)
                , 
                gp = gpar(lwd=0.5),
                draw=FALSE)
  
}

plotCytoBand <- function(
    chrom, beg, end, genomeInfo, txn.font.size) {
  
  cytoBand <- genomeInfo$cytoBand
  
  ## set cytoband color
  requireNamespace("pals")
  cytoBand2col <- setNames(
    pals::ocean.gray(10)[seq(9,3)],
    c('stalk', 'gneg', 'gpos25', 'gpos50', 'gpos75', 'gpos100'))
  cytoBand2col['acen'] <- 'red'
  cytoBand2col['gvar'] <- cytoBand2col['gpos75']
  
  ## chromosome range
  cytoBand.target <- cytoBand[cytoBand$chrom == chrom,]
  chromEnd <- max(cytoBand.target$chromEnd)
  chromBeg <- min(cytoBand.target$chromStart)
  chromWid <- chromEnd - chromBeg
  bandColor <- cytoBand2col[as.character(cytoBand.target$gieStain)]
  
  pltx0 <- (c(beg, end)-chromBeg)/chromWid
  gList(
    grid.text( # coordinate name
      sprintf("%s:%d-%d", chrom, beg, end), 0, 0.9,
      just = c('left','bottom'), gp = gpar(fontfamily="sans", fontsize = as.numeric(1.1*txn.font.size)), draw = FALSE),
    ## cytoband box
    grid.rect(0, 0.35, 1, 0.35, just = c("left", "bottom"),
        gp = gpar(col = "black", lwd=2, lty="solid"), draw = FALSE),
    grid.rect( # cytoband
      vapply(cytoBand.target$chromStart,
             function(x) (x-chromBeg)/chromWid, 1),
      0.35,
      (cytoBand.target$chromEnd - cytoBand.target$chromStart)/chromWid,
      0.35, gp = gpar(fill = bandColor, col = bandColor),
      just = c('left','bottom'), draw = FALSE),
  #  grid.rect(0, 0.35, 1, 0.2, just = c("left", "bottom"), draw = FALSE),
  #  grid.segments( # projection lines
  #    x0 = pltx0, y0 = 0.3, x1 = c(0,1), y1 = 0.1, draw = FALSE),
    grid.segments( # sentinel bar
      x0 = pltx0, y0 = 0.1, x1 = pltx0, y1 = 0.9,
      gp = gpar(col = "red"), draw = FALSE))
}


rename_probes_to_loci_and_de_duplicate_if_needed <- function(probes.list.or.betas, probes) {
  
  
  probes <- GenomicRanges::as.data.frame(probes)
  probes["genomic.loci"] <-  paste0(probes$seqnames, ":", probes$start, "-", probes$end)
  probes["probes"] <-  rownames(probes)
  
  #converting probes to genomic loci in probes.list.or.betas 'rownames' //START
  probes.list.or.betas <- GenomicRanges::as.data.frame(probes.list.or.betas)
  probes.list.or.betas["probes"] <- rownames(probes.list.or.betas)
  
  probes.list.or.betas  <- merge(probes[c("genomic.loci", "probes")], probes.list.or.betas, by = "probes")
  
  #below 'if/else statement' will determine if any genomic loci coordinates are repeated, if so, de-duplicate by paste0'ing together probe_ID and genomic loci by an '_' in between and setting that to the row.name, instead of just the genomic.loci
  #for duplicates only, the paste0'd (combined probe_ID and genomic loci coordinates) will be their heatmap "x-axis", otherwise only genomic loci coordinate will be used
  if (any(duplicated(probes.list.or.betas$genomic.loci))) {
    
    #identifies *all* duplicate loci - by marking TRUE in the column "redundant_loci" for the corresponding row
    probes.list.or.betas[probes.list.or.betas$genomic.loci %in% unique(probes.list.or.betas[duplicated(probes.list.or.betas$genomic.loci),"genomic.loci"]), "redundant_loci"] <- "TRUE"
    
    #wherever there is TRUE in redundant loci, probe and genomic loci and combined via paste0 with a '_' delimiter in between - and then the rowname is changed to the paste0 output
    rownames(probes.list.or.betas)[which(probes.list.or.betas$redundant_loci == TRUE)] <- paste0(probes.list.or.betas[which(probes.list.or.betas$redundant_loci == TRUE),"probes"], "_", probes.list.or.betas[which(probes.list.or.betas$redundant_loci == TRUE),"genomic.loci"])
    
    #wherevere there is an NA in "redunant_loci" column change (where no duplicates exist) - change the rowname to the genomic loci
    rownames(probes.list.or.betas)[which(is.na(probes.list.or.betas$redundant_loci))] <- probes.list.or.betas[which(is.na(probes.list.or.betas$redundant_loci)), "genomic.loci"]
    
    #NULL/clear out redundant loci column
    probes.list.or.betas$redundant_loci <- NULL
    probes.list.or.betas
    
  } else {
    
    row.names(probes.list.or.betas) <- probes.list.or.betas$genomic.loci
  }
  probes.list.or.betas$probes <- NULL
  probes.list.or.betas$genomic.loci <- NULL
  probes.list.or.betas <- as.matrix(probes.list.or.betas)
  #converting probes to genomic loci in probes.list.or.betas 'rownames' //END
  probes.list.or.betas
}
        
#this is a function to convert probes to loci and also de-duplicate loci (by combining probeID with the genomic loci)
#this is made for the assemble_plots function (the same code was more or less being used twice, so to reduce lines of code, it was made into one function)


#' assemble plots
#'
#' @param betas beta value
#' @param txns transcripts GRanges
#' @param probes probe GRanges
#' @param plt.txns transcripts plot objects
#' @param plt.mapLines map line plot objects
#' @param plt.cytoband cytoband plot objects
#' @param heat.height heatmap height (auto inferred based on rows)
#' @param show.probeNames whether to show probe names
#' @param show.samples.n number of samples to show (default: all)
#' @param show.sampleNames whether to show sample names
#' @param sample.name.fontsize sample name font size
#' @param dmin data min
#' @param dmax data max
#' @return a grid object
assemble_plots <- function(
    betas, txns, probes, plt.txns, plt.mapLines, plt.cytoband, plt.loci_txn_intersectionLines, font.size.scaling.factor,
    heat.height = NULL, mapLine.height = 0.2,
    show.probeNames = TRUE, show.samples.n = NULL, hypothesis_generation = TRUE,
    show.sampleNames = TRUE, sample.name.fontsize = 10,
    dmin = 0, dmax = 1) {
  
  #hypothesis_generation = TRUE is only really meant to be used with (mouse) mm39, since the regulatory features were obtained from Ensembl 108 (mouse) mm39
  #it will show loci, instead of probes, genomic regulatory features, and cpg island
  if (hypothesis_generation) {
    #START - Addition by Pratik - bring in regulatory features from mft
    mft <- sesameDataGet('MM285.mm39.manifest')
    probe.list <- data.frame(mft[which(names(mft) %in% rownames(betas))])
    rownames(probe.list) <- probe.list$probeID
    probe.list <- probe.list[row.names(betas),]
    
    cpg.and.regulatory.feature.probe.list <- rename_probes_to_loci_and_de_duplicate_if_needed(probe.list, probes)
    #END - Addition by Pratik - bring in regulatory features from mft
    
    #START - PRATIK - add fixed colors for annotations
    #tol.rainbow color pallete is color-blind friendly
    regulatory_features_4_colors <- pals::tol.rainbow(n=12)[1:4] 
    regulatory_features_ctcf_colors <- pals::tol.rainbow(n=12)[5]
    
    #from unique(mft$regulatory_feature)[1:4]
    names(regulatory_features_4_colors) <- c("predicted enhancer", "predicted promoter", "open chromatin region", "TF binding site")
    #from unique(mft$CTCF_binding_site)[2]
    names(regulatory_features_ctcf_colors) <- "CTCF binding site"
    
    cpg_island_color <- c("#009e73", "#d55e00", "#000000", "#e69f00", "#0072B2")
    #below from unique(mft$CpG_Island)[2]
    names(cpg_island_color) <- c("CpG Island",  "CpG Shore", "CpG Shore-Shelf Transition", "CpG Shelf",  "Open Sea")
    #START - PRATIK - add fixed colors for annotations
    
    betas <- rename_probes_to_loci_and_de_duplicate_if_needed(betas, probes)
  }
  
  
  if (is.null(show.samples.n)) { show.samples.n <- ncol(betas); }
  if (is.null(heat.height) && length(txns) > 0) {
    heat.height <- 15 / length(txns); }
  w <- WGrob(plt.txns, name = 'txn')
  w <- w + WGrob(plt.mapLines, name = 'mpl', Beneath(pad=0, height=mapLine.height))
  w <- w + WGrob(plt.loci_txn_intersectionLines, TopOf('mpl'))
  w <- w + WHeatmap(
    t(betas), Beneath('mpl', height = heat.height),
    name = 'betas',
    cmp = CMPar(dmin=dmin, dmax=dmax),
    xticklabels = show.probeNames,
    xticklabel.rotat = 45,
    yticklabels = show.sampleNames,
    xticklabel.pad = .4,
    yticklabel.fontsize = as.numeric(sample.name.fontsize*font.size.scaling.factor),
    xticklabel.fontsize = as.numeric(sample.name.fontsize*font.size.scaling.factor),
    yticklabels.n = show.samples.n,
    xticklabels.n = length(probes))
  
  w <- w + WGrob(plt.cytoband, TopOf('txn', height=0.15))
  
  if (hypothesis_generation) {    
    w <- w + WLegendV(x = "betas", RightOf('betas', h.scale = 'betas', h.scale.proportional = TRUE, pad = .02), n.text = 2, "betalegend", decreasing = TRUE)  +
      WLabel(x= "Beta Values", TopOf("betalegend"), name = 'betaleglab', fontsize = as.numeric(11*font.size.scaling.factor)) +
      WColorBarH(probe.list$CpG_Location_Type, Beneath('betas'), cmp=CMPar(label2color = cpg_island_color), label = "CpG Island", label.fontsize = as.numeric(9*font.size.scaling.factor), label.side = 'l', "CpG") +
      WColorBarH(probe.list$regulatory_feature, Beneath(v.scale.proportional = TRUE), cmp=CMPar(label2color = regulatory_features_4_colors), label.fontsize = as.numeric(9*font.size.scaling.factor), label = "Genomic Feature Type", label.side = 'l', 'GFT') +
      WColorBarH(probe.list$CTCF_binding_site, Beneath(v.scale.proportional = TRUE), cmp=CMPar(label2color = regulatory_features_ctcf_colors),label.fontsize = as.numeric(9*font.size.scaling.factor), label = "CTCF Binding", label.side = 'l', 'CTCF') +        
      WLabel(x= "CpG Island Annotations", BottomRightOf(x = "betalegend", just = c('center', 'top'), v.pad = -.2), name = 'cpgleglab', fontsize = as.numeric(11*font.size.scaling.factor)) +
      WLegendV(x= "CpG",name = "cpgleg", Beneath("cpgleglab"),  height = as.numeric(heat.height/show.samples.n), cmp=CMPar(label2color = cpg_island_color), label.fontsize = as.numeric(8.5*font.size.scaling.factor)) +
      WLabel(x= "Ensembl Regulatory Build Annotations", name = "ensregleg", Beneath("cpgleg", pad = 0.1), fontsize = as.numeric(11*font.size.scaling.factor)) +
      WLegendV(x= "GFT", name = "ensreg", Beneath("ensregleg", pad = 0), label.fontsize = as.numeric(8.5*font.size.scaling.factor), height = as.numeric(heat.height/show.samples.n)) +
      WLegendV(x= "CTCF",Beneath("ensreg", pad = 0), label.fontsize = as.numeric(8.5*font.size.scaling.factor), height = as.numeric(heat.height/show.samples.n))
    
    #changing order to cpg island, shore, shore-shelf transition, shelf, open sea [  rather than an order that is automatically-generated and doesn't go with the "island" theme ;)  ]
    cpg_island_annotation_reordered <- na.omit(w$children[[12]]$cm$colors[c("CpG Island",  "CpG Shore", "CpG Shore-Shelf Transition", "CpG Shelf",  "Open Sea")])
    attr(cpg_island_annotation_reordered,"na.action") <- NULL
    
    w$children[[12]]$cm$colors <- cpg_island_annotation_reordered
    w$children[[12]]$data <- t(t(w$children[[12]]$data[c(names(cpg_island_annotation_reordered)),]))
  }
  w +  WCustomize(mar.bottom = ( 1/ ( .3* (length(probes)+length(txns)) ) ) )
}

