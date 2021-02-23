prepareBioCov <- function(Decade, 
                          climateModel, 
                          pathInputs,
                          rasterToMatch = NULL,
                          studyArea = NULL){
  browser()
  bioCovs <- data.table(decade = c(2050, 
                                   2070),
                        URL = c("http://biogeo.ucdavis.edu/data/climate/cmip5/5m/cc85bi50.zip", 
                                "http://biogeo.ucdavis.edu/data/climate/cmip5/5m/cc85bi70.zip"),
                        climmod = c("CCSM4", 
                                    "CCSM4"))
  climPath <- checkPath(file.path(pathInputs,
                                  climateModel),
                        create = TRUE)
  prepInputs(url = bioCovs[decade == Decade, URL],
             destinationPath = climPath)
  stk <- raster::stack(lapply(X = list.files(climPath, pattern = ".tif", 
                                             full.names = TRUE), 
                              FUN = raster))
  currNames <- usefulFuns::substrBoth(strng = names(stk), 
                                      howManyCharacters = 2, 
                                      fromEnd = TRUE)
  newNames <- paste0("bio", currNames)
  if (any(!is.null(studyArea), !is.null(rasterToMatch))){
    stk[] <- stk[]
    stk <- postProcess(stk, studyArea = studyArea,
                       rasterToMatch = rasterToMatch)
  }
  names(stk) <- newNames
  return(stk)
}