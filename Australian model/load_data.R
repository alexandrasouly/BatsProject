load_data <- function(site, year){
  
  # loading observed data table from chosen year, filtering by site 
  
if (as.character(year) =="all"){
  file_name <- paste("prev_data_by_sites_2017.txt")
  data_2017<-read.table(as.character(file_name))
  
  file_name <- paste("prev_data_by_sites_2018.txt")
  data_2018<-read.table(as.character(file_name))
  
  file_name <- paste("prev_data_by_sites_2019.txt")
  data_2019<-read.table(as.character(file_name))
  
  merge(merge(data_2017, data_2018, all = TRUE), data_2019, all = TRUE) -> data_all_year
  
  
  # file_name <- paste("prev_data_by_sites.txt")
  # data_all_year<-read.table(as.character(file_name))
  # strptime(data_all_year[['min_date']],format = "%Y-%m-%d") -> data_all_year[['min_date']]
  # format(data_all_year[['min_date']],format = "%j") -> data_all_year[['min_date']]
  # 
  # data_all_year[["min_date"]] = as.numeric(data_all_year[["min_date"]])
  # 
  data_all_year %>%
    filter(site_abbrev == site) -> data_all_year
    
  data_all_year %>%
    select("min_date", "number_of_tests") -> samplesize_dat
  colnames(samplesize_dat) <- c("time","samplesize")

  data_all_year %>%
    filter(site_abbrev == site) %>%
    select("min_date", "hen_num_pos") -> pos_dat
  colnames(pos_dat) <- c("time","pos")
  
}else {
  
  file_name <- paste("prev_data_by_sites_", year, ".txt", sep = "")
  data_2019<-read.table(as.character(file_name))
  data_2019 %>% 
    filter(site_abbrev == site) %>% 
    select("min_date", "number_of_tests") -> samplesize_dat
  colnames(samplesize_dat) <- c("time","samplesize")
  
  data_2019 %>% 
    filter(site_abbrev == site) %>% 
    select("min_date", "hen_num_pos") -> pos_dat
  colnames(pos_dat) <- c("time","pos")
}
  
  return(list(samplesize_dat, pos_dat))
}