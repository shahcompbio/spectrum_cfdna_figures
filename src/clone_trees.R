
get_clone_tree_old <- function(dat, loci_nodes = NULL, plottree = FALSE){
  
  mypatient <- dat$sample
  #get clones
  if ("clones" %in% names(dat) == FALSE){
    dat <- get_clones(dat, loci_nodes = loci_nodes, plottree = plottree)
  }
  
  dat$clones <- dat$clones %>% filter(clone_id != "0")
  
  #set clone_id to minor clone if there is one
  clonedata <- dat$clones %>% 
    mutate(clone_id = ifelse(is.na(clone_id.minor), clone_id, clone_id.minor), 
           node_id = ifelse(is.na(node_id.minor), node_id, node_id.minor), 
           node_label = ifelse(is.na(node_label.minor), node_label,
                               node_label.minor)) %>% 
    distinct(clone_id, node_id, node_label) %>% 
    as.data.frame() 
  row.names(clonedata) <- clonedata$node_label
  nodes_to_retain <- clonedata$node_label
  clone_names <- clonedata$clone_id
  
  #replace node label with clone_id 
  mytree <- dat$tree
  mytree$node.label <- str_remove(mytree$node.label, "locus_")
  for (i in 1:nrow(clonedata)){
    mytree$node.label <- str_replace(mytree$node.label, clonedata$node_label[i], clonedata$clone_id[i])
  }
  
  #remove all other labels
  mytree$node.label[!mytree$node.label %in% clone_names] <- ""
  
  if (plottree == TRUE){
    g4 <- ggtree(mytree %>% compute.brlen(.,1)) + geom_nodelab(size = 5)
  }
  
  #iterate down the tree removing tips if not clone_id
  stopping <- all(mytree$tip.label %in% clone_names)
  while (stopping == FALSE){
    mytree <- drop.tip(mytree, setdiff(mytree$tip.label, clone_names),
                       trim.internal = F, 
                       collapse.singles = F)
    stopping <- all(mytree$tip.label %in% clone_names)
  }
  
  clonedatafortree <- dat$clones %>% 
    mutate(clone_id = ifelse(is.na(clone_id.minor), clone_id, clone_id.minor), 
           node_id = ifelse(is.na(node_id.minor), node_id, node_id.minor), 
           node_label = ifelse(is.na(node_label.minor), node_label, node_label.minor)) %>% 
    select(cell_id, clone_id) %>% as.data.frame()
  
  if (plottree == TRUE){
    g1 <- ggtree(mytree) + geom_tiplab() + geom_nodelab()
    g2 <- plot_grid(plotlist = dat$cloneplots)
    g3 <- signals:::make_tree_ggplot(dat$tree %>% compute.brlen(.,1), clones = clonedatafortree)
    plot_grid(g1, g2, g3, g4) %>% print()
  }
  
  return(mytree)
}


get_clone_tree <- function(tree, clones, loci_nodes = NULL, removeinternal = TRUE){
  source(here("src/tree_utils.R"))
  
  
  clones <-clones %>% filter(clone_id != "0")
  
  #set clone_id to minor clone if there is one
  clonedata <- clones %>% 
    # mutate(clone_id = ifelse(is.na(clone_id.minor), clone_id, clone_id.minor), 
    #        node_id = ifelse(is.na(node_id.minor), node_id, node_id.minor), 
    #        node_label = ifelse(is.na(node_label.minor), node_label,
    #                            node_label.minor)) %>% 
    distinct(clone_id, node_id, node_label) %>% 
    as.data.frame() 
  row.names(clonedata) <- clonedata$node_label
  nodes_to_retain <- clonedata$node_label
  clone_names <- clonedata$clone_id
  
  #replace node label with clone_id 
  mytree <- tree
  mytree$node.label <- str_remove(mytree$node.label, "locus_")
  for (i in 1:nrow(clonedata)){
    mytree$node.label <- str_replace(mytree$node.label, clonedata$node_label[i], clonedata$clone_id[i])
  }
  
  #remove all other labels
  #mytree$node.label[!mytree$node.label %in% clone_names] <- ""
  
  #iterate down the tree removing tips if not clone_id
  stopping <- all(mytree$tip.label %in% clone_names)
  while (stopping == FALSE){
    mytree <- drop.tip(mytree, setdiff(mytree$tip.label, clone_names),
                       trim.internal = F, 
                       collapse.singles = F)
    stopping <- all(mytree$tip.label %in% clone_names)
  }
  
  if (removeinternal == TRUE){
    mytree <- collapse.singles(mytree)
  }
  
  clonedatafortree <- clones %>% 
    mutate(clone_id = ifelse(is.na(clone_id.minor), clone_id, clone_id.minor), 
           node_id = ifelse(is.na(node_id.minor), node_id, node_id.minor), 
           node_label = ifelse(is.na(node_label.minor), node_label, node_label.minor)) %>% 
    select(cell_id, clone_id) %>% as.data.frame()
  
  return(mytree)
}
