prepareBioCov <- function(Decade, 
                          climateModel, 
                          pathInputs,
                          climateURLTable,
                          rasterToMatch = NULL,
                          studyArea = NULL){
  
  climPath <- checkPath(file.path(pathInputs,
                                  climateModel),
                        create = TRUE)
  message(crayon::yellow("Downloding climate data for model ", climateModel, 
                         " for decade ", Decade, " from \n", 
                         climateURLTable[decade == Decade, URL]))
  suppressMessages(prepInputs(url = climateURLTable[decade == Decade, URL],
                              destinationPath = climPath))
  stk <- raster::stack(lapply(X = list.files(climPath, pattern = ".tif", 
                                             full.names = TRUE), 
                              FUN = raster))
  currNames <- usefulFuns::substrBoth(strng = names(stk), 
                                      howManyCharacters = 2, 
                                      fromEnd = TRUE)
  # Need to fix the stack names... Remove the leading zeros
  currNames <- sub("^0+", "", currNames)
  newNames <- paste0("bio", currNames)
  if (any(!is.null(studyArea), !is.null(rasterToMatch))){
    stk <- stack(lapply(names(stk), function(lay){
      message(crayon::magenta(paste0("Post-processing ", 
                                     lay," layer")))
      ras <- stk[[lay]]
      ras[] <- ras[]
      rasProc <- postProcess(ras, studyArea = studyArea,
                         rasterToMatch = rasterToMatch)
      return(rasProc)
    }))
  }
  names(stk) <- newNames
  return(stk)
}