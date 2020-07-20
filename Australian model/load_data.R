load_data <- function(site, year){
  
  # loading observed data table from chosen year, filtering by site 
  
  file_name <- paste("prev_data_by_sites_", year, ".txt", sep = "")
  data_2019<-read.table(as.character(file_name))
  data_2019 %>% 
    filter(site_abbrev == site) %>% 
    select("min_date", "number_of_tests") -> samplesize_dat
  colnames(samplesize_dat) <- c("time","samplesize")
  
  data_2019 %>% 
    filter(site_abbrev == "CLU") %>% 
    select("min_date", "hen_num_pos") -> pos_dat
  colnames(pos_dat) <- c("time","pos")
  
  return(list(samplesize_dat, pos_dat))
}