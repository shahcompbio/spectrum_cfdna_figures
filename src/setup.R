options(java.parameters = "-Xmx8000m")
library(tidyverse)
library(devtools)
library(data.table)
library(cowplot)
library(yaml)
library(glue)
library(signals)
library(ape)
library(ggtree)
library(cowplot)
library(ggpubr)
library(xlsx)
#source(here("src/tree_utils.R"))

config <- read_yaml(here("config.yaml"))

font_size <- config$font_size
gg_font_size <- font_size / 0.8 
gg_theme <- theme_cowplot(font_size = 7, font_family = "sans", line_size = 0.25)
theme_set(gg_theme)
#theme_update(plot.margin = unit(c(5.5,5.5,5.5,5.5), "points"))

## read i some files

purity <- readxl::read_xlsx(config$tables, sheet = "S2 - Tumor Fractions")

svs <- readxl::read_xlsx(config$tables, sheet = "S7 - SV read counts")
svs_md <- readxl::read_xlsx(config$tables, sheet = "S9 - SV metadata")
svs <- svs %>% left_join(svs_md,by = join_by(patient, target_idx))

snvs <- readxl::read_xlsx(config$tables, sheet = "S8 - SNV read counts")
snvs_md <- readxl::read_xlsx(config$tables, sheet = "S10 - SNV metadata")
snvs <- snvs %>% left_join(snvs_md,by = join_by(target_idx, clonality, Hugo_Symbol))

recurrence <- readxl::read_xlsx(config$tables, sheet = "S12 - recurrence data") %>% 
  mutate(date = as.numeric(date))

svdf_clone <- readxl::read_xlsx(config$tables, sheet = "S3 - Clone frequencies (SVs)") %>% 
  mutate(nclones = 1) %>% 
  left_join(purity %>% select(patient, timepoint, tumorfraction_ClonalSV)) %>% 
  group_by(patient, timepoint) %>% 
  mutate(max_clone_frequency = max(clone_frequency)) %>% 
  ungroup() %>% 
  left_join(recurrence %>% filter(recurrence_number == 1) %>% rename(patient = patient_id, recurrence_date = date))

svs$Type <- lapply(svs$type, stringr::str_to_title) %>% unlist()
svs_md$Type <- lapply(svs_md$type, stringr::str_to_title) %>% unlist()
svs$Bamtype <- lapply(svs$bamtype, stringr::str_to_title) %>% unlist()
snvs$Bamtype <- lapply(snvs$bamtype, stringr::str_to_title) %>% unlist()

removexaxis <- theme(axis.line.x=element_blank(),
                     axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())

removeyaxis <- theme(axis.line.y=element_blank(),
                     axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())

format_tree <- function(tree, removeloci = TRUE){
  if (removeloci){
    #tree <- removeleafloci(tree)
    tip.loci <- grep('locus', tree$tip.label, value = T)
    while (length(tip.loci) > 0) {
      tree <- ape::drop.tip(tree, tip.loci, trim.internal = FALSE, collapse.singles = FALSE)
      tip.loci <- grep('locus', tree$tip.label, value = T)
    }
  }
  tree$tip.label <- str_remove(tree$tip.label, "cell_")
  tree$node.label <- str_remove(tree$node.label, "locus_")
  #tree <- collapse.singles(tree)
  #tree <- ape::compute.brlen(tree, 1)
  return(tree)
}

read_file <- function(file, sep = "_", verbose = FALSE){
  if (verbose == TRUE){
    print(file)
  }
  df <- fread(file)
  df$patient <- strsplit(basename(file), sep)[[1]][1]
  return(df)
}

read_allfiles <- function(pattern, patients = NULL, sep = "_", ncores = 4, verbose = FALSE){

  files <- find_allfiles(pattern, patients, sep = sep)

  if (ncores > 1){
  df <- files %>%
    parallel::mclapply(., function(x) read_file(x, sep = sep, verbose = verbose), mc.cores = ncores) %>%
    rbindlist(fill = TRUE)
  } else{
    df <- files %>%
      lapply(., function(x) read_file(x, sep = sep, verbose = verbose)) %>%
      rbindlist(fill = TRUE)
  }

  print(df$patient %>% unique %>% sort)

  return(df)
}

find_allfiles <- function(pattern, patients = NULL, sep = "_", ncores = 4, verbose = FALSE){
  files <- list.files(config$root, pattern = pattern, recursive = T, full.names = T)
  if (is.null(patients)){
    patients <- config$patients
  }
  patient_pattern <- paste(patients, collapse = "|")
  files <- files[grepl(patient_pattern, files)]
  print(files)
  return(files)
}

logyscale <- function(sigma = 0.0001, limits = c(0.0, NA), scientific = TRUE) {
  mybreaks <-  unlist(lapply(1:10, function(x) sigma*10^x))
  mybreaks <- c(0, mybreaks[mybreaks < 2])
  if(scientific == FALSE){
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = sigma, base = 10), breaks = mybreaks, limits = limits)
  } else{
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = sigma, base = 10), breaks = mybreaks, limits = limits, labels = scales::trans_format("log10", scales::math_format(10^.x)))
  }
}

logxscale <- function(sigma = 0.0001, limits = c(0.0, NA), scientific = TRUE) {
  mybreaks <-  unlist(lapply(1:10, function(x) sigma*10^x))
  mybreaks <- c(0, mybreaks[mybreaks < 2])
  if(scientific == FALSE){
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma = sigma, base = 10), breaks = mybreaks, limits = limits)
  } else{
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma = sigma, base = 10), breaks = mybreaks, limits = limits, labels = scales::trans_format("log10", scales::math_format(10^.x)))
  }
}