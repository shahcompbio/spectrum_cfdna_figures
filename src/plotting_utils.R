get_clone_bar_plot <- function(mypt){
  
  
  cp <- config$clonecols
  names(cp) <- LETTERS[1:length(cp)]
  
  svdf_clone %>% 
    filter(tumorfraction_ClonalSV > 0) %>% 
    filter(str_detect(patient, mypt)) %>% 
    filter(nclones == 1) %>% 
    filter(max_clone_frequency > 0) %>% 
    mutate(recurrence = days >= recurrence_date) %>% 
    group_by(patient, timepoint, days, recurrence) %>% 
    mutate(f= clone_frequency / sum(clone_frequency)) %>% 
    ungroup() %>% 
    dplyr::select(recurrence, patient, clone_id, timepoint, days, f) %>% 
    group_by(patient) %>% 
    filter(timepoint == 1 | timepoint == max(timepoint)) %>% 
    mutate(has_recurrence = any(recurrence)) %>% 
    filter(has_recurrence == TRUE) %>% 
    ungroup() %>% 
    mutate(pt = str_remove(patient, "SPECTRUM-OV-")) %>% 
    ggplot(aes(x = recurrence, y = f, fill = clone_id)) +
    geom_col(width = 0.5) +
    scale_fill_manual(values = cp) +
    scale_x_discrete(labels = c("B", "R"),
                     guide = guide_axis(angle = 0)) +
    xlab("") +
    labs(fill = "") +
    ylab("Clone frequency")
}