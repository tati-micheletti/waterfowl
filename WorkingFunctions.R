library("Require")
Require("biomod2")
Require("SpaDES")
Require("raster")
Require("data.table")
Require("googledrive")

googledrive::drive_auth(email = "tati.micheletti@gmail.com")

speciesURL <- data.table(species = c("REDH", "BWTE", "BUFF"),
                         URL = c("https://drive.google.com/file/d/1JZ6ihDdn_WCSK05SIYrVqDq-YhRoGirj/view?usp=sharing",
                                 "https://drive.google.com/file/d/1BruuNFQlEw1t9IrYgN_OP3j92rF3ObcZ/view?usp=sharing",
                                 "https://drive.google.com/file/d/1escVgr4I15iizMYa47M_vBZIKJTunnFt/view?usp=sharing"))

covariates <- "https://drive.google.com/file/d/1T729j8czLKugK3yUmJ3SZVLplISQQ6n-/view?usp=sharing"

rasterOptions(tmpdir="C:/Users/Tati/Google Drive/Postdoc PFC-UBC/WBI/waterfowl/rasterTemp")

setPaths(cachePath = "C:/Users/Tati/Google Drive/Postdoc PFC-UBC/WBI/waterfowl/cache", 
         inputPath = "C:/Users/Tati/Google Drive/Postdoc PFC-UBC/WBI/waterfowl/inputs",
         outputPath = "C:/Users/Tati/Google Drive/Postdoc PFC-UBC/WBI/waterfowl/outputs")

# Download the models and data for each species

  speciesPaths <- rbindlist(lapply(1:NROW(speciesURL), function(index){
    sp <- speciesURL[index, species]
    URL <- speciesURL[index, URL]
    if (!dir.exists(file.path(Paths$inputPath, 
                   sp))){
      filePath <- file.path(Paths$inputPath,
                            paste0(sp, ".zip"))
      
      drive_download(file = URL, path = filePath)
      unzip(zipfile = filePath, exdir = Paths$inputPath)
      speciesPath <- data.table(species = sp,
                                spInputPath = file.path(Paths$inputPath, 
                                                        sp))
      } else {
      speciesPath <- data.table(species = sp,
                                spInputPath = file.path(Paths$inputPath, 
                                                        sp))
    } 
    return(speciesPath)
  }))

# Retrieve current environmental covariate layers
## load all covariates for the specific period of time
prepareBioCov <- function(Decade, climateModel, 
                          variables = NULL, 
                          rasterToMatch = NULL,
                          studyArea = NULL){
  bioCovs <- data.table(decade = c(2050, 
                                   2070),
                        URL = c("http://biogeo.ucdavis.edu/data/climate/cmip5/5m/cc85bi50.zip", 
                                "http://biogeo.ucdavis.edu/data/climate/cmip5/5m/cc85bi70.zip"),
                        climmod = c("CCSM4", 
                                    "CCSM4"))
  climPath <- checkPath(file.path(Paths$inputPath,
                      climateModel),
                      create = TRUE)
  prepInputs(url = bioCovs[decade == Decade, URL],
             destinationPath = climPath)
  stk <- raster::stack(lapply(X = list.files(climPath, pattern = ".tif", 
                                             full.names = TRUE), 
                              FUN = raster))
  currNames <- usefulFuns::substrBoth(strng = names(stk), howManyCharacters = 2, 
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

# bioCov <- prepareBioCov(Decade = 2050, 
#                         climateModel = "CCSM4") # NOT WORKING YET
#                         
# treeCov <- rasterizeReducedForspeciesURL() # TODO

# NEED TO MAKE THE ONES BELOW
# "SpeciesGroups_Broadleaf_Spp" --> Check exactly what they are!! 
# "SpeciesGroups_Needleleaf_Spp" --> Check exactly what they are!!

# JUST TO TEST. SHOULD BE REMOVED WHEN prepareBioCov() is working
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
allCov <- prepInputs(targetFile = "covariates.rds",
                                    url = covariates, fun = "readRDS", 
                                    destinationPath = Paths$inputPath) # For now use the ready covars
bioCov <- raster::subset(allCov, names(allCov)[1:19])
treeCov <- raster::subset(allCov, names(allCov)[20:nlayers(allCov)])
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## For each species, make new projections
predictForSpecies <- function(sp, 
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

allPredictions <- lapply(1:NROW(speciesURL), function(index){
  sp <- speciesURL[index, species]
  filePath <- file.path(Paths$inputPath, sp)
  pred <- predictBiomod(sp = sp, 
                            inputPath = filePath,
                            bioCov = bioCov,
                            treeCov = treeCov)
})
names(allPredictions) <- speciesURL[["species"]]


