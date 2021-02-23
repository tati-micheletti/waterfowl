## For each species, make new projections
predictBiomod <- function(sp, 
                          inputPath,
                          bioCov,
                          treeCov){
  message(crayon::yellow(paste0("Predictions starting for ", sp, "...")))
  originalWD <- getwd()
  on.exit(setwd(originalWD))
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
  setwd(Paths$inputPath)
  biomod2::BIOMOD_LoadModels(modInd)
  sel.cov <- modInd@expl.var.names
  currentStack <- raster::stack(bioCov, treeCov)
  current.env.stack <- stack(currentStack[[which(names(currentStack) %in% sel.cov)]])
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
  ensemble_models_proj_future <-BIOMOD_EnsembleForecasting(EM.output = modEns,
                                                           projection.output = models_proj_future,
                                                           binary.meth = "TSS", # or ROC ...
                                                           output.format = ".grd",
                                                           do.stack = FALSE)/1000 # Correcting for integer forcing
  
  message(crayon::green(paste0("Predictions done for ", sp)))
  return(ensemble_models_proj_future)
}
