prepareTreeCov <- function(cohortData,
                           pixelGroupMap,
                           biomassMap,
                           rasterToMatch){
  
  # Individual species biomass
  speciesNames <- as.character(unique(cohortData[["speciesCode"]]))
  speciesRasters <- stack(lapply(X = speciesNames, FUN = function(sp){
    subsCohort <- cohortData[speciesCode == sp, ]
    zeroedMap <- raster(pixelGroupMap)
    zeroedMap[] <- getValues(pixelGroupMap)
    vals <- getValues(x = zeroedMap)
    vals[!is.na(vals)] <- 0
    zeroedMap <- setValues(x = zeroedMap, values = vals)
    
    if (NROW(subsCohort) == 0){
      assign(x = sp, value = zeroedMap)
      names(zeroedMap) <- paste0("Species_", sp)
      return(zeroedMap)
    } else {
      valsCoho <- data.table(pixelID = 1:ncell(pixelGroupMap), 
                             pixelGroup = getValues(x = pixelGroupMap))
      joinOn <- c("speciesCode", "pixelGroup")
      newCohoVals <- valsCoho[subsCohort[, list(sumBiomass = sum(B)), by = joinOn], 
                              on = "pixelGroup"]
      zeroedMap[newCohoVals$pixelID] <- newCohoVals$sumBiomass
      assign(x = as.character(sp), value = zeroedMap)
      names(zeroedMap) <- paste0("Species_", sp)
      return(zeroedMap)
    }
  }))

  # "SpeciesGroups_Needleleaf_Spp" and "SpeciesGroups_Broadleaf_Spp"  
  decidousSp <- c("Betu_Pap", "Popu_Tre", "Popu_Bal")
  treeSpecies <- data.table(speciesCode = c(decidousSp,
                                            speciesNames[!speciesNames %in% decidousSp]),
                            treeType = c(rep("broadleaf", 
                                            times = length(decidousSp)),
                                        rep("conifer", 
                                            times = length(speciesNames[!speciesNames %in% decidousSp]))))
  cohortData <- merge(cohortData, treeSpecies)
  # Calculate the sum of treeType biomass per pixel group
  cohortData[, treeTypeBiomass := sum(B), by = c("pixelGroup", "treeType")]
  cohortDataReduced <- unique(cohortData[, c("pixelGroup", "treeType", "treeTypeBiomass")])
  
  # Conifers
  conifers <- cohortDataReduced[treeType == "conifer"]
  SpeciesGroups_Needleleaf_Spp <- rasterizeReduced(reduced = conifers, 
                                                   fullRaster = pixelGroupMap, 
                                                   newRasterCols = "treeTypeBiomass", 
                                                   mapcode = "pixelGroup")
  names(SpeciesGroups_Needleleaf_Spp) <- "SpeciesGroups_Needleleaf_Spp"
  
  # SpeciesGroups_Broadleaf_Spp
  broadleaf <- cohortDataReduced[treeType == "broadleaf"]
  SpeciesGroups_Broadleaf_Spp <- rasterizeReduced(reduced = broadleaf, 
                                                   fullRaster = pixelGroupMap, 
                                                   newRasterCols = "treeTypeBiomass", 
                                                   mapcode = "pixelGroup")
  names(SpeciesGroups_Broadleaf_Spp) <- "SpeciesGroups_Broadleaf_Spp"
  
  SpeciesGroups <- stack(SpeciesGroups_Broadleaf_Spp, SpeciesGroups_Needleleaf_Spp)
  
  treeCov <- raster::stack(speciesRasters, SpeciesGroups)
  # Need to rescale all tree layers using raster::rescale(lays) 
  # (accepting all defaults of rescaling and centering) from Antoine Adde 
  # [22FEB21 per email]
  treeCov <- raster::stack(lapply(names(treeCov), FUN = function(lay){
    message(paste0("Rescaling and centering variable ", lay))
    ras <- treeCov[[lay]]
    ras[] <- ras[]
    rasRescaled <- suppressWarnings(raster::scale(ras))
    # And also fill in the holes
    rasRescaled[is.na(rasRescaled[]) & 
                  rasterToMatch[] == 1] <- 0
    return(rasRescaled)
  }))
  
  return(treeCov)
}