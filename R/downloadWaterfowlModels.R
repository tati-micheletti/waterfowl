downloadWaterfowlModels <- function(speciesURL, dataFolder){
  speciesPaths <- rbindlist(lapply(1:NROW(speciesURL), function(index){
    sp <- speciesURL[index, species]
    URL <- speciesURL[index, URL]
    if (!dir.exists(file.path(dataFolder, 
                              sp))){
      filePath <- file.path(dataFolder,
                            paste0(sp, ".zip"))
      
      drive_download(file = URL, path = filePath)
      unzip(zipfile = filePath, exdir = dataFolder)
      speciesPath <- data.table(species = sp,
                                spInputPath = file.path(dataFolder, 
                                                        sp))
    } else {
      speciesPath <- data.table(species = sp,
                                spInputPath = file.path(dataFolder, 
                                                        sp))
    } 
    return(speciesPath)
  }))
  return(speciesPath)
}