## For each species, make new projections
predictBiomod <- function(sp, 
                          inputPath,
                          bioCov,
                          treeCov,
                          currentTime,
                          overwritePred = FALSE){
  message(crayon::yellow(paste0("Predictions starting for ", sp, "...")))
  originalWD <- getwd()
  on.exit(setwd(originalWD))
  setwd(dirname(inputPath))
  predRas <- file.path(getwd(), paste0(sp, "/proj_future/individual_projections/", 
                                       sp,"_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.grd"))
  if (any(!file.exists(predRas), 
          overwritePred)){
    fl <- list.files(path = inputPath, 
                     pattern = "FirstModeling.models.out",
                     full.names = TRUE, 
                     recursive = TRUE)
    modInd <- get(load(file = fl))
    # Needs to decompress
    fls <- list.files(path = inputPath, pattern = "gz",
                      full.names = TRUE, 
                      recursive = TRUE)
    if (any(unique(tools::file_ext(fls)) %in% "gz")){
      lapply(fls, R.utils::gunzip)
    }
    biomod2::BIOMOD_LoadModels(modInd)
    sel.cov <- modInd@expl.var.names
    currentStack <- raster::stack(bioCov, treeCov)
    # Make sure that any missing layers are created and added
    missingLays <- setdiff(sel.cov, names(currentStack))
    if (length(missingLays) != 0){
      zeroLays <- stack(lapply(missingLays, function(lay){
        template <- currentStack[[1]] # template
        template[] <- template[]
        template[!is.na(template[])] <- 0
        names(template) <- lay
        return(template)
      })) 
    } else zeroLays <- NA
    # Now we need to select which layers will be used
    current.env.stack <- stack(currentStack[[which(names(currentStack) %in% sel.cov)]])
    # And add the missing layers
    current.env.stack <- stack(current.env.stack, zeroLays)
    # Check if all layers are there
    testthat::expect_true(all(names(current.env.stack) %in% sel.cov))
    models_proj_future <- BIOMOD_Projection(modeling.output = modInd,
                                            new.env = current.env.stack,
                                            proj.name = "future",
                                            binary.meth = "TSS", # or ROC ...
                                            output.format = ".grd",
                                            do.stack = TRUE)
    fl <- list.files(path = inputPath, 
                     pattern = "FirstModelingensemble.models.out",
                     full.names = TRUE, 
                     recursive = TRUE)
    modEns <- get(load(file = fl))
    ## Ensemble projections
    ensemble_models_proj_future <- BIOMOD_EnsembleForecasting(EM.output = modEns,
                                                              projection.output = models_proj_future,
                                                              binary.meth = "TSS", # or ROC ...
                                                              output.format = ".grd",
                                                              do.stack = FALSE)
    
    testthat::expect_true(file.exists(predRas))
  }
  # Now the file exists. Load, correct and return
  predRas <- raster::raster(predRas)
  # Correcting for integer forcing (x1000) we need to divide the results
  # by 1000
  predRas <- predRas/1000
  names(predRas) <- paste0(sp, "_Year", currentTime)
  message(crayon::green(paste0("Predictions done for ", 
                               sp, " for year ", currentTime)))
  return(predRas)
}
