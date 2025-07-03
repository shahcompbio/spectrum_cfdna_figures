sv_plot <- function(hscn, clustering, chrfilt = NULL, maxCN = 30, svwidth = 0.4, geneannot = NULL, svs = NULL, posticks = F, pointsize = 1, homolog = FALSE, font_size = 8){
  
  sample_clusters <- clustering
  hscn_new <- hscn
  
  glist <- list()
  for (myclone in sort(unique(sample_clusters$clone_id))){
    
    if (is.null(svs)){
      svsforplot <- NULL
    } else{
      svsforplot <- svs %>% 
        filter( clone_id == myclone) %>% 
        filter(VAF > 0)
      if (nrow(svsforplot) == 0){
        svsforplot <- NULL
      }
    }
    
    mycells <- sample_clusters %>% filter(clone_id == myclone) %>% pull(cell_id)
    if (homolog){
      glist[[myclone]] <- hscn_new[cell_id %in% mycells] %>% 
        consensuscopynumber() %>% 
        plotCNprofileBAF(., svwidth = svwidth,returnlist = F,font_size = gg_font_size,
                      annotateregions = signals:::get_gene_idx(geneannot, chr = chrfilt), raster = raster,
                      SV = svsforplot,positionticks = posticks, tickwidth = 50,
                      chrfilt = chrfilt, pointsize = pointsize, homolog = homolog,
                      y_axis_trans = "squashy", maxCN = maxCN, legend.position = "none") +
        ggtitle(paste0("Clone ",myclone, " (", length(mycells), " cells)")) +
        theme_cowplot(font_size = font_size) +
        theme(legend.position = "none") +
        theme(plot.title = element_text(size=font_size))
    } else{
      glist[[myclone]] <- hscn_new[cell_id %in% mycells] %>% 
        consensuscopynumber() %>% 
        plotCNprofile(., svwidth = svwidth,returnlist = F,font_size = gg_font_size,
                      annotateregions = signals:::get_gene_idx(geneannot, chr = chrfilt), 
                      SV = svsforplot, positionticks = posticks, tickwidth = 25,
                      chrfilt = chrfilt, pointsize = pointsize,raster = raster,
                      y_axis_trans = "squashy", maxCN = maxCN, legend.position = "none") +
        theme_cowplot(font_size = font_size) +
        theme(legend.position = "none") +
        ggtitle(paste0("Clone ",myclone, " (", length(mycells), " cells)")) +
        theme(plot.title = element_text(size=font_size))
    }
    
  }
  
  gall <- plot_grid(plotlist = glist)
  return(list(plot = gall, plotlist = glist))
}

sv_plot_clone <- function(hscn,clustering, svline_size = 0.2, chrfilt = NULL, maxCN = 30, svwidth = 0.4, geneannot = NULL, svs = NULL, posticks = F, pointsize = 1, homolog = FALSE, line_size = 0.25,  font_size = 7, ...){
  
  sample_clusters <- clustering
  hscn_new <- hscn
  myclones <- unique(hscn$clone_id) 
  
  glist <- list()
  glistsv <- list()
  for (myclone in myclones){
    print(myclone)
    mycells <- sample_clusters %>% filter(clone_id == myclone) %>% pull(cell_id)
    if (is.null(svs)){
      svsforplot <- NULL
    } else{
      svsforplot <- svs %>% 
        filter( clone_id == myclone & (chromosome_1 %in% chrfilt | chromosome_2 %in% chrfilt)) %>% 
        filter(VAF > 0)
      if (nrow(svsforplot) == 0){
        svsforplot <- NULL
      }
    }
    
    if (homolog){
      glist[[myclone]] <- hscn_new %>% 
        filter(clone_id == myclone) %>% 
        plotCNprofileBAF(., svwidth = svwidth,returnlist = F,font_size = gg_font_size,
                         ideogram = T,
                         annotateregions = signals:::get_gene_idx(geneannot, chr = chrfilt), 
                         SV = svsforplot,positionticks = posticks, tickwidth = 50,
                         chrfilt = chrfilt, pointsize = pointsize, homolog = homolog,
                         y_axis_trans = "squashy", maxCN = maxCN, legend.position = "none", ...) +
        ggtitle(paste0("Clone ",myclone, " (", length(mycells), " cells)")) +
        theme_cowplot(font_size = font_size, line_size = line_size) +
        theme(legend.position = "none") +
        theme(plot.title = element_text(size=font_size))
      
      glistsv[[myclone]] <- svsforplot %>% 
        plotSVlines(chrfilt = chrfilt) +
        theme_cowplot(font_size = font_size, line_size = line_size) +
        ggtitle(paste0("Clone ",myclone, " (", length(mycells), " cells)")) +
        theme(plot.title = element_text(size=font_size))
    } else{
      glist[[myclone]] <- hscn_new %>% 
        filter(clone_id == myclone) %>% 
        plotCNprofile(., svwidth = svwidth,returnlist = F,font_size = gg_font_size,
                      ideogram = T,
                      annotateregions = signals:::get_gene_idx(geneannot, chr = chrfilt), 
                      SV = svsforplot, positionticks = posticks, tickwidth = 50,
                      chrfilt = chrfilt, pointsize = pointsize,
                      y_axis_trans = "squashy", maxCN = maxCN, legend.position = "none", ...) +
        theme_cowplot(font_size = font_size, line_size = line_size) +
        theme(legend.position = "none") +
        ggtitle(paste0("Clone ",myclone, " (", length(mycells), " cells)")) +
        theme(plot.title = element_text(size=font_size))
      
      if (is.null(svsforplot)){
        glistsv[[myclone]] <- NULL
      } else {
        glistsv[[myclone]] <- svsforplot %>% 
          plotSVlines(chrfilt = chrfilt, line_width = svline_size) +
          theme_cowplot(font_size = font_size, line_size = line_size) +
          ggtitle(paste0("Clone ",myclone, " (", length(mycells), " cells)")) +
          theme(plot.title = element_text(size=font_size))
      }
    }
    
  }
  
  glistcnsv <- list()
  myclones <- myclones[myclones != "Pseudobulk"]
  for (myclone in myclones){
    glistcnsv[[myclone]] <- plot_grid(glistsv[[myclone]] + removexaxis +  
                theme(legend.position = "none",
                      #axis.line.y = ggplot2::element_blank(),
                      axis.text.y = ggplot2::element_blank(),
                      axis.ticks.y = ggplot2::element_blank()),
              NULL,
              glist[[myclone]] + ggtitle(NULL), ncol = 1,
              rel_heights = c(0.3, -0.3, 1),
              align = "hv", axis = "tblr")
  }
  
  
  gall <- plot_grid(plotlist = glist)
  gallcnsv <- plot_grid(plotlist = glistcnsv)
  return(list(plot = gall,plotsv = gallcnsv, plotlist = glist, plotlistcnsv = glistcnsv, plotlistsv = glistsv))
}
