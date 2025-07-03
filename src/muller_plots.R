normalize_vafs <- function(svdf_clone, small = 0.001234){
  x1 <- svdf_clone %>% 
    filter(nclones == 1) %>% 
    group_by(days) %>% 
    mutate(cloneVAFn = cloneVAFn / sum(cloneVAFn)) %>% 
    mutate(cloneVAF = cloneVAF / sum(cloneVAF)) %>% 
    #mutate(cloneVAFn = ifelse(cloneVAFn == 0, 0.0001, cloneVAFn)) %>% 
    ungroup() %>% 
    mutate(cloneVAFn = ifelse(is.nan(cloneVAFn), small, cloneVAFn)) %>% 
    mutate(cloneVAF = ifelse(is.nan(cloneVAF), small, cloneVAF))
  return(x1)
}

smooth_vafs <- function(x1){
  clones <- unique(x1$clone_id)
  timeline <- seq(min(x1$days),max(x1$days), 2)
  dfnew <- data.frame()
  for (myclone in clones){
    clonedat <- x1 %>% filter(clone_id == myclone)
    fun <- splinefun(x = clonedat$days, y = clonedat$cloneVAFn, method = "monoH.FC")
    dfnew <- data.frame(days = timeline, cloneVAFn = fun(timeline), clone_id = myclone) %>% 
      bind_rows(dfnew, .)
  }
  
  dfnew <- dfnew %>% 
    mutate(cloneVAFn = ifelse(cloneVAFn < 0, 0, cloneVAFn)) %>%
    mutate(cloneVAFn = ifelse(cloneVAFn > 1, 1, cloneVAFn)) %>% 
    #mutate(cloneVAFn = ifelse(cloneVAFn == 0, 0.0001, cloneVAFn)) %>% 
    group_by(days) %>% 
    mutate(cloneVAFn = cloneVAFn / sum(cloneVAFn)) %>% 
    ungroup() %>% 
    mutate(cloneVAFn = ifelse(is.nan(cloneVAFn), 0, cloneVAFn))
  return(dfnew)
}


smooth_vafs2 <- function(x1){
  clones <- unique(x1$clone_id)
  timeline <- seq(min(x1$days),max(x1$days), 2)
  dfnew <- data.frame()
  for (myclone in clones){
    clonedat <- x1 %>% filter(clone_id == myclone)
    fun <- splinefun(x = clonedat$days, y = clonedat$cloneVAF, method = "monoH.FC")
    dfnew <- data.frame(days = timeline, cloneVAF = fun(timeline), clone_id = myclone) %>% 
      bind_rows(dfnew, .)
  }
  
  dfnew <- dfnew %>% 
    mutate(cloneVAF = ifelse(cloneVAF < 0, 0, cloneVAF)) %>%
    mutate(cloneVAF = ifelse(cloneVAF > 1, 1, cloneVAF)) %>% 
    #mutate(cloneVAFn = ifelse(cloneVAFn == 0, 0.0001, cloneVAFn)) %>% 
    mutate(cloneVAF = ifelse(is.nan(cloneVAF), 0, cloneVAF)) %>% 
    mutate(cloneVAFn = cloneVAF)
  return(dfnew)
}

clone_tree_plot <- function(mysample){
  print(mysample)
  tree <- ape::read.tree(find_allfiles(config$root, pattern = "newick", patient = mysample))
  tree$node.label <- tree$node.label %>% str_remove(., "_")
  clones <- read_allfiles("clones.csv", patients = mysample)
  if (length(unique(clones$clone_id)) == 1){
    return(ggplot())
  }
  
  clonecols <- config$clonecols
  names(clonecols) <- LETTERS[1:length(clonecols)]
  
  clone_tree <- get_clone_tree(tree, clones, removeinternal = T) %>% compute.brlen(.,1)
  tree_plot <- ggtree(clone_tree, size = 0.75) + 
    #geom_tiplab(aes(col = label), size = 2.5, hjust = -1.5) + 
    geom_tippoint(aes(col = label), size = 1.5, shape = 15) +
    scale_color_manual(values = clonecols) +
    theme(legend.position = "none") #+geom_treescale(x = 0)
  tree_plot
  return(tree_plot)
}

plot_timeline_newv <- function(mypatient, 
                               recurrence,
                               svdf_clone,
                               surgeries,
                               treatcols,
                               ca125,
                               purity,
                               meta_long,
                               zerotimepoints = NULL,
                               scale = NULL,
                               expandval = 0.01,
                               gap = -0.05,
                               gap2 = -0.05,
                               add_rec = TRUE,
                               includeca125 = TRUE,
                               removelegend = TRUE,
                               add_surgery = FALSE,
                               plotclonalSV = FALSE,
                               plotTP53 = FALSE,
                               onlyTP53 = FALSE,
                               onlySV = FALSE,
                               muller_smooth_cutoff = 1.0,
                               overwritelastT = FALSE,
                               plot_title = NULL,
                               plotlist = FALSE,
                               alphaval = 1.0,
                               rename_cols = TRUE,
                               type = "muller1", ...){
  
  #treatments_patient <- treatments %>% filter(str_detect(patient_id,mypatient))
  recurrence_patient <- recurrence %>% filter(str_detect(patient_id,mypatient))
  surgeries_patient <- surgeries %>% filter(str_detect(patient_id,mypatient))
  svdf_clone_patient <- svdf_clone %>% filter(str_detect(patient,mypatient))
  ca125_patient <- ca125 %>% filter(str_detect(patient_id,mypatient))
  purity_patient <- purity %>% filter(str_detect(patient,mypatient))
  meta_long_patient <- meta_long %>% filter(str_detect(ptid,mypatient))
  
  if (rename_cols == TRUE){
    svdf_clone_patient <- rename(svdf_clone_patient, 
                                 cloneVAF = clone_frequency,
                                 cloneVAFn = clone_frequency_normalized)
    recurrence_patient <- rename(recurrence_patient, d = date)
    purity_patient <- rename(purity_patient, TP53 = tumorfraction_TP53,
                             Clonal_sv = tumorfraction_ClonalSV)
  }
  
  if (!is.null(zerotimepoints)){
    df_to_add <- expand.grid(zerotimepoints, unique(svdf_clone_patient$clone_id)) %>% 
      dplyr::rename(days = Var1, clone_id = Var2) %>% 
      mutate(cloneVAFn = 0, patient = paste0("SPECTRUM-OV-", mypatient), nclones = 1)
    svdf_clone2 <- bind_rows(svdf_clone_patient, df_to_add) %>% 
      replace_na(list(cloneVAF = 0, fracSV = 0, nSVs = 0, T = 0))
  } else{
    svdf_clone2 <- svdf_clone_patient %>% 
      filter(str_detect(patient, mypatient))
  }
  
  svdf_clone2 <- svdf_clone2 %>% 
    filter(clone_id %in% LETTERS)
  
  clonal2 <- svdf_clone2 %>% 
    #filter(nclones == 1) %>% 
    #normalize_vafs() %>% 
    smooth_vafs2() %>% 
    group_by(days) %>% 
    summarize(cloneVAF = sum(cloneVAF))
  
  days_data <- svdf_clone2 %>% distinct(days, timepoint) %>% 
    filter(!days %in% zerotimepoints)
  
  #nclones <- length(unique(svdf_clone2$clone_id))
  
  cp <- MetBrewer::met.brewer("VanGogh2", override.order = TRUE)
  cp <- cp[c(1,3,5,7,8,2,4,6)]
  names(cp) <- LETTERS[1:length(cp)]
  
  maxx <- ceiling(max(c(svdf_clone2$days)) / 10) * 10
  minx <- min(c(svdf_clone2$days, surgeries_patient$surgery_date, 0))
  if (overwritelastT == TRUE){
    maxx <- max(ca125_patient$days)
  } else if (is.numeric(overwritelastT)){
    maxx <- overwritelastT
  }
  
  options(scipen = 999)
  if (type == "raw"){
    gclone <- svdf_clone2 %>% 
      #filter(nclones == 1) %>% 
      filter(!is.na(patient)) %>% 
      ggplot(aes(x = days, y = cloneVAF, col = clone_id)) +
      geom_vline(data = days_data, aes(xintercept = days, col = NULL, x = NULL, y = NULL), lty = 2, alpha = 0.5, lwd = 0.2) +
      geom_point(shape = 16, size = 0.5, position = position_jitter(height = 0, width = 0)) +
      geom_line(alpha = 0.75) +
      scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.00001, base = 10), 
                         breaks = c(0.0, 0.0001, 0.001, 0.01, 0.1, 1.0),
                         labels = c("0.0", "0.0001", "0.001", "0.01", "0.1", "1.0"),
                         limits = c(0, NA)) +
      # geom_segment(data = days_data, aes(x = days, xend = days, 
      #                                    y = -0.0003, yend = 0, col = NULL), 
      #              arrow = arrow(angle = 21, length = unit(0.04, "npc"))) +
      #scale_fill_manual(values = cp) +
      scale_color_manual(values = cp) +
      scale_x_continuous(limits = c(minx, maxx), expand = c(expandval, expandval)) +
      xlab("Days") +
      ylab("Clone\nfraction") +
      theme_cowplot(font_size = font_size, ...) +
      theme(legend.position = "top") +
      guides(col = guide_legend(nrow = 1, title = element_text("Clone")))
  } else if (type == "muller1") { 
    gclone <- svdf_clone2 %>%
      arrange(days) %>% 
      #filter(nclones == 1) %>% 
      #normalize_vafs() %>% 
      smooth_vafs2() %>% 
      ggplot( aes(x = days, fill = clone_id)) +
      geom_area(aes(y = cloneVAF), size = 1) +
      geom_segment(data = days_data, aes(x = days, xend = days, 
                                         y = -max(clonal2$cloneVAF) * 0.1, yend = 0, fill = NULL), size = 0.2,
                   arrow = arrow(angle = 21, type = "closed", length = unit(0.03, "npc"))) +
      xlab("Days") +
      ylab("Clone\nfraction") +
      scale_x_continuous(limits = c(minx, maxx), expand = c(expandval, expandval)) +
      theme_cowplot(font_size = font_size, ...) +
      #jcolors::scale_color_jcolors(palette = "default") +
      scale_fill_manual(values = cp) +
      theme(legend.position = "top") +
      guides(fill = guide_legend(nrow = 1, title = element_text("Clone")))
  } else if (type == "muller2") { 
    gclone <- svdf_clone2 %>%
      #filter(nclones == 1) %>% 
      normalize_vafs() %>% 
      smooth_vafs2() %>% 
      group_by(days) %>% 
      mutate(sumvaf = sum(cloneVAFn)) %>% 
      mutate(cloneVAFn = ifelse(sumvaf > muller_smooth_cutoff, cloneVAFn / sumvaf, cloneVAFn))  %>% 
      ggplot( aes(x = days, fill = clone_id)) +
      geom_area(aes(y = cloneVAFn), size = 1, alpha = alphaval) +
      geom_segment(data = days_data, aes(x = days, xend = days, 
                                         y = -0.1, yend = 0, fill = NULL), size = 0.2,
                   arrow = arrow(angle = 21, type = "closed", length = unit(0.03, "npc"))) +
      xlab("Days") +
      ylab("Clone\nfraction") +
      scale_x_continuous(limits = c(minx, maxx), expand = c(expandval, expandval)) +
      theme_cowplot(font_size = font_size, ...) +
      #jcolors::scale_color_jcolors(palette = "default") +
      scale_fill_manual(values = cp) +
      theme(legend.position = "top") +
      guides(fill = guide_legend(nrow = 1, title = element_text("Clone")))
  }
  
  if (plotTP53 == TRUE){
    if (is.null(scale)){
      scale <- max(clonal2$cloneVAF) / max(purity_patient$TP53vaf)
    }
    if (all(is.na(purity_patient$TP53vaf))){
      purity_patient$TP53vaf <- purity_patient$purity
    }
    gclone <- gclone +
      geom_point(data = purity_patient, 
                 aes(y = scale * TP53vaf, fill = NULL), 
                 shape = 4) +
      geom_line(data = purity_patient, 
                aes(y = scale * TP53vaf, fill = NULL), 
                lty = 2) +
      scale_y_continuous(sec.axis = sec_axis(~.*scale, name="TP53 VAF"))
  }
  
  if (plotclonalSV == TRUE){
    if (is.null(scale)){
      scale <- max(clonal2$cloneVAF) / max(purity_patient$Clonal_sv)
    }
    if (all(is.na(purity_patient$TP53vaf))){
      purity_patient$TP53vaf <- purity_patient$purity
    }
    if (type == "raw"){
      gclone <- gclone +
        geom_point(data = purity_patient, 
                   aes(y = scale * Clonal_sv, fill = NULL), 
                   col = "black",
                   shape = 16, size = 0.5) +
        geom_line(data = purity_patient, 
                  aes(y = scale * Clonal_sv, fill = NULL), 
                  col = "black",
                  lty = 2) +
        scale_y_continuous(trans = 
                             scales::pseudo_log_trans(sigma = 0.00001, base = 10), 
                           breaks = c(0.0, 0.0001, 0.001, 0.01, 0.1, 1.0),
                           labels = c("0.0", "0.0001", "0.001", "0.01", "0.1", "1.0"),
                           limits = c(0, NA))
    } else{
      gclone <- gclone +
        geom_point(data = purity_patient, 
                   aes(y = scale * Clonal_sv, fill = NULL), 
                   col = "black",
                   shape = 16, size = 0.5) +
        geom_line(data = purity_patient, 
                  aes(y = scale * Clonal_sv, fill = NULL), 
                  col = "black",
                  lty = 2) +
        scale_y_continuous(sec.axis = sec_axis(~.*scale, name="Clonal SV VAF"))
    }
  }
  
  if (onlyTP53 == TRUE){
    gclone <- purity_patient %>% 
      ggplot() +
      geom_point(
        aes(y = TP53vaf, x = days, fill = NULL), 
        shape = 4) +
      geom_line(
        aes(y = TP53vaf, x = days, fill = NULL), 
        lty = 2) +
      scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.0001, base = 10), breaks = c(0.0, 0.001, 0.01, 0.1, 1.0), limits = c(0, NA)) +
      ylab("TP53 VAF") +
      theme_cowplot(font_size = font_size, ...) +
      scale_x_continuous(limits = c(minx, maxx), expand = c(expandval, expandval))
  }
  
  if (onlySV == TRUE){
    gclone <- purity_patient %>% 
      ggplot() +
      geom_point(
        aes(y = Clonal_sv, x = days, fill = NULL), 
        shape = 4) +
      geom_line(
        aes(y = Clonal_sv, x = days, fill = NULL), 
        lty = 2) +
      scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.0001, base = 10), breaks = c(0.0, 0.001, 0.01, 0.1, 1.0), limits = c(0, NA)) +
      ylab("SV VAF") +
      theme_cowplot(font_size = font_size, ...) +
      scale_x_continuous(limits = c(minx, maxx), expand = c(expandval, expandval))
  }
  
  if (removelegend == TRUE){
    gclone <- gclone + theme(legend.position = "none")
  }
  
  gSV_TP53 <- purity_patient %>% 
    ggplot() +
    geom_point(
      aes(y = Clonal_sv, x = days, fill = NULL), 
      shape = 17, size = 0.5, alpha = 0.75) +
    geom_line(
      aes(y = Clonal_sv, x = days, fill = NULL), 
      linetype = "dotted", size = 0.25, alpha = 0.75) +
    geom_point(
      aes(y = TP53, x = days, fill = NULL), 
      shape = 18, col = "firebrick4", size = 0.75, alpha = 0.75) +
    geom_line(
      aes(y = TP53, x = days, fill = NULL), 
      col = "firebrick4", size = 0.25, alpha = 0.75) +
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.00001, base = 10), 
                       breaks = c(0.0,0.0001, 0.001, 0.01, 0.1, 1.0), limits = c(0, NA),
                       labels = scales::trans_format("log10", scales::math_format(10^.x)),
                       guide = guide_axis(check.overlap = T)) +
    ylab("VAF") +
    theme_cowplot(font_size = font_size, ...) +
    scale_x_continuous(limits = c(minx, maxx), expand = c(expandval, expandval))
  
  gca125 <- ca125_patient %>% 
    filter(days < maxx) %>% 
    ggplot(aes(x = days, y = CA125)) +
    geom_line(size = 0.2) +
    scale_x_continuous(limits = c(minx, maxx), expand = c(expandval, expandval)) +
    scale_y_continuous(guide = guide_axis(check.overlap = T)) +
    geom_point(size = 0.2) +
    xlab("Days") +
    theme_cowplot(font_size = font_size, ...) +
    ylab("CA-125")
  
  # gtreatment <- treatments_patient %>% 
  #   filter(days < maxx) %>% 
  #   #filter(!str_detect(treatment_generic_name, "TRAS")) %>% 
  #   ggplot(aes(x = days, y = treatment_generic_name, col = treatment_generic_name)) +
  #   geom_point(size = 1, shape = 15) +
  #   theme(legend.position = "none") +
  #   scale_x_continuous(breaks = scales::breaks_width(365), limits = c(minx, maxx)) +
  #   xlab("Time from first surgery (days)") +
  #   ylab("") +
  #   #scale_color_brewer(palette = "Dark2") +
  #   scale_color_manual(values = treatcols) +
  #   theme_cowplot(font_size = font_size, ...) +
  #   theme(panel.grid.major.y = element_line(size=.1, color="grey50", linetype = "dashed")) +
  #   theme(legend.position = "none", axis.text.y = element_text(size = 5))
  
  gtreatment <-  meta_long_patient %>% 
    ggplot(aes(x = time_from_dx, y = fct_reorder(ptid, Final_tp))) +
    geom_blank() +
    ### treatment lines/layers
    geom_linerange(data = meta_long_patient[meta_long_patient$event_type == "treatments", ], 
                   aes(xmin = time_from_dx - 15, xmax = time_from_dx + 15, # change to modulate overlap of   treatment datapoints
                       y = as.numeric(fct_reorder(ptid, Final_tp)) + 3 * treatment_level, 
                       color = event_name), 
                   lwd = 2.0, alpha = 1) +
    scale_color_manual(values = treatcols) +
    theme_cowplot(font_size = font_size, ...) +
    theme(legend.position = "none") +
    ylab("") +
    scale_y_discrete(labels = c("Treatment")) +
    xlab("Time from first surgery (days)") +
    scale_x_continuous(breaks = scales::breaks_width(365), limits = c(minx, maxx), 
                       expand = c(expandval, expandval))
  
  if (add_rec == TRUE){
    gclone <- gclone + geom_rect(data = recurrence_patient,
                                 aes(xmin = d, x = NULL, y = NULL, fill = NULL, col = NULL), 
                                 xmax = Inf, ymin = -Inf, ymax = Inf, 
                                 fill = "grey60", color = NA, alpha = .2) #+
    #guides(col = guide_legend(nrow = 1, title = element_text("Clone")), fill = "none")
    gtreatment <- gtreatment + geom_rect(data = recurrence_patient,
                                         aes(xmin = d, x = NULL, y = NULL, fill = NULL, col = NULL), 
                                         xmax = Inf, ymin = -Inf, ymax = Inf, 
                                         fill = "grey60", color = NA, alpha = .2)
    gca125 <- gca125 + geom_rect(data = recurrence_patient,
                                 aes(xmin = d, x = NULL, y = NULL, fill = NULL, col = NULL), 
                                 xmax = Inf, ymin = -Inf, ymax = Inf, 
                                 fill = "grey60", color = NA, alpha = .2)
    gSV_TP53 <- gSV_TP53 + geom_rect(data = recurrence_patient,
                                     aes(xmin = d, x = NULL, y = NULL, fill = NULL, col = NULL), 
                                     xmax = Inf, ymin = -Inf, ymax = Inf, 
                                     fill = "grey60", color = NA, alpha = .2)
  }
  
  if (add_surgery == TRUE){
    gclone <- gclone + geom_vline(xintercept = surgeries_patient$surgery_date, lty = 2)
    gtreatment <- gtreatment + geom_vline(xintercept = surgeries_patient$surgery_date, lty = 2)
  }
  
  if (!is.null(plot_title)){
    gclone <- gclone + 
      ggtitle(plot_title) +
      theme(plot.title = element_text(face = "plain"))
  }
  
  if (includeca125 == TRUE){
    final1 <- plot_grid(gclone + removexaxis,# + ggtitle(mypatient),
                        NULL,
                        gca125 + removexaxis,
                        NULL,
                        gtreatment, 
                        ncol = 1,
                        rel_heights = c(0.8,gap,0.65,gap2, 0.65), 
                        align = "v", axis = "lr")
  } else {
    final1 <- cowplot::plot_grid(gclone + removexaxis,# + ggtitle(mypatient),
                                 NULL,
                                 gtreatment, 
                                 ncol = 1,
                                 rel_heights = c(1.5,gap,0.8), 
                                 align = "v", axis = "lr")
  }
  
  if (type == "muller2"){
    final1 <- cowplot::plot_grid(gclone + removexaxis,# + ggtitle(mypatient),
                                 NULL,
                                 gSV_TP53 + removexaxis,
                                 NULL,
                                 gca125 + removexaxis,
                                 NULL,
                                 gtreatment, 
                                 ncol = 1,
                                 rel_heights = c(0.8,gap,0.65,gap,0.65,gap2, 0.85), 
                                 align = "v", axis = "lr")
  }
  
  if (plotlist == TRUE){
    final1 <- list(gclone = gclone, gca125 = gca125, gtreatment = gtreatment, gSV_TP53 = gSV_TP53)
  }
  
  return(final1)
}

