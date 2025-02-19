
library(lipdR)
library(tidyverse)
library(readxl)
library(lubridate)

dat1 <- read_excel("data/cside_sst_final_forNick.xlsx",sheet = 1)
dat2 <- read_excel("data/cside_sst_final_forNick.xlsx",sheet = 2)


lipdNames <- c("dataSetName","paleoData_core","geo_latitude","geo_longitude","geo_accZone","geo_location","geo_elevation","values_depth","values_age","values_temperature","paleoData_proxy","pub1_citation","pub1_doi")


names(dat1) <- names(dat2) <- lipdNames


allDsn <- unique(dat1$dataSetName)

createLipd <- function(thisDataSetName,dataSource){

  thisDsn <- filter(dataSource,dataSetName == thisDataSetName)

  #get the values data
  valOnly <- select(thisDsn,starts_with("values"))

  variableNames <- names(valOnly) %>% str_remove_all("values_")

  thisTs <- select(thisDsn,-starts_with("values")) %>% distinct()

  #make metadata fixes as needed
  thisTs$geo_elevation <- -abs(thisTs$geo_elevation ) #make it negative (meters above sea level)
  thisTs$dataSetName <- str_remove_all(thisTs$dataSetName,"[^A-Za-z0-9]")
  thisTs$datasetId <- paste0("cside2lipd",thisTs$dataSetName)

  if(nrow(thisTs) != 1){
    stop("there should be exactly 1 row")
  }
  thisRow <- thisTs
  while(nrow(thisTs) < length(variableNames)){
    thisTs <- bind_rows(thisTs,thisRow)
  }

  thisTs$paleoData_variableName <- variableNames
  thisTs$paleoData_values <- NA
  thisTs$paleoData_units <- NA
  thisTs$paleoData_TSid <- NA
  thisTs$paleoData_number <- NA




  for(v in 1:length(variableNames)){
    thisTs$paleoData_values[v] <- list(as.matrix(valOnly[,v]))
    thisTs$paleoData_units[v] <- case_when(thisTs$paleoData_variableName[v] == "depth" ~ "m",
                                           thisTs$paleoData_variableName[v] == "age" ~ "ka",
                                           thisTs$paleoData_variableName[v] == "temperature" ~ "deg C")
    thisTs$paleoData_TSid[v] <- paste0("cside2lipd",thisTs$dataSetName[v],thisTs$paleoData_units[v]) %>%
      str_remove_all("[^A-Za-z0-9]")
    thisTs$paleoData_number[v] <- v
  }


  L <- purrr::transpose(thisTs) %>%
    collapseTs(force = TRUE)

  L$archiveType <- "MarineSediment"
  L$lipdVersion <- 1.3
  L$createdBy <- "https://github.com/nickmckay/cside2lipd"
  L <- initializeChangelog(L)

  return(L)
}


L <- createLipd(allDsn[3],dat1)

validLipd(L)

