create_edge_table <- function(my_tree){
  
  if (!"node.label" %in% names(my_tree)){
    my_tree$node.label <- paste0("node_", seq(1:my_tree$Nnode))
  }
  
  node_labels_in_edge <- my_tree$node.label[my_tree$edge[,1]-ape::Ntip(my_tree)]
  tips_nodes <- my_tree$edge[,2]
  
  select.tip.or.node <- function(element, tree) {
    ifelse(element < ape::Ntip(tree)+1, tree$tip.label[element], tree$node.label[element-ape::Ntip(tree)])
  }
  
  edge_table <- data.frame(
    "parent" = my_tree$edge[,1],
    "par.name" = sapply(my_tree$edge[,1], select.tip.or.node, tree = my_tree),
    "child" = my_tree$edge[,2],
    "chi.name" = sapply(my_tree$edge[,2], select.tip.or.node, tree = my_tree)
  )
}
