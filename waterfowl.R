defineModule(sim, list(
  name = "waterfowl",
  description = paste0("This module forecasts waterfowl habitat suitability ",
                       "index for 3 species of waterfowl for a given study ",
                       "aerea within the Western Boreal Region"),
  keywords = "WBI",
  authors = structure(list(list(given = "Tati", family = "Micheletti", 
                                role = c("aut", "cre"), 
                                email = "tati.micheletti@gmail.com", 
                                comment = NULL)), class = "person"),
  childModules = character(0),
  version = list(SpaDES.core = "1.0.5", waterfowl = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = deparse(list("README.txt", "waterfowl.Rmd")),
  reqdPkgs = list("biomod2", "raster", "data.table", "googledrive", "crayon"),
  parameters = rbind(
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first plot event should occur."),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "Describes the simulation time interval between plot events."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first save event should occur."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    paste("Should this entire module be run with caching activated?",
                          "This is generally intended for data-type modules, where stochasticity",
                          "and time are not relevant")),
    defineParameter("predictionsDecades", "numeric", c(2050, 2070), NA, NA,
                    paste0("Describes the simulation time interval between prediction events. ",
                           "Defaults to 2050 and 2070 --> These are ensambles ",
                           "but uses simulated data of the closest years ending in 1 to match the 10 year ",
                           "interval of data saving")),
    defineParameter("decadesMatchingYears", "numeric", c(2051, 2071), NA, NA,
                    paste0("Years for which this module should run (i.e. ",
                           "years that should match the predictionsDecades for which data ",
                           "(cohort data and pixelGroup map) are available",
                           "interval of data saving")),
    defineParameter("overwritePred", "logical", FALSE, NA, NA,
                    paste0("If predictions exist, should be overwritten?"))
  ),
  inputObjects = bindrows(
    expectsInput(objectName = "speciesURL", objectClass = "data.table", 
                 desc = paste0("Data.table with 2 columns, species (4-letter code)",
                               " and URL (url of google drive FOLDER where the ",
                               "fitted models from biomod2 were saved). It is ",
                               "important that the structure of the models is ",
                               "kept intact as the predict function depends on ",
                               "it."), 
                 sourceURL = NA),
    expectsInput(objectName = "species", objectClass = "character", 
                 desc = paste0("Character vector with species (4-letter code). ",
                               "These need to be available in speciesURL table."), 
                 sourceURL = NA),
    expectsInput(objectName = "studyArea", objectClass = "SpatialPolygonDataFrame", 
                 desc = "Study area for the prediction. Currently only available for NWT", 
                 sourceURL = "https://drive.google.com/open?id=1P4grDYDffVyVXvMjM-RwzpuH1deZuvL3"),
    expectsInput(objectName = "rasterToMatch", objectClass = "RasterLayer",
                 desc = "All spatial outputs will be reprojected and resampled to it", 
                 sourceURL = "https://drive.google.com/open?id=1P4grDYDffVyVXvMjM-RwzpuH1deZuvL3"),
    expectsInput(objectName = "climateURLTable", objectClass = "data.table",
                 desc = "Data.table with information on decade, url and climate model", 
                 sourceURL = NA)
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "biomassMap", objectClass = "RasterLayer", 
                  desc = "Total biomass map"),
    createsOutput(objectName = "cohortData", objectClass = "data.table", 
                  desc = "Table with cohort information (biomass per species per pixelGroup)"),
    createsOutput(objectName = "pixelGroupMap", objectClass = "RasterLayer", 
                  desc = "Mapping raster to pixelGroup"),   
    createsOutput(objectName = "waterfowlPredictions", objectClass = "RasterLayer", 
                  desc = "Mapping raster to pixelGroup"),
    createsOutput(objectName = "bioCov", objectClass = "RasterStack", 
                  desc = "Mapping raster to pixelGroup") ,
    createsOutput(objectName = "treeCov", objectClass = "RasterStack", 
                  desc = "Mapping raster to pixelGroup")
  )
))

doEvent.waterfowl = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      sim$.schedulingCounter <- 1
      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$decadesMatchingYears[sim$.schedulingCounter], "waterfowl", "modelsPrep")
      sim <- scheduleEvent(sim, P(sim)$decadesMatchingYears[sim$.schedulingCounter], "waterfowl", "climateDataPrep")
      sim <- scheduleEvent(sim, P(sim)$decadesMatchingYears[sim$.schedulingCounter], "waterfowl", "simulatedLayersPrep")
      sim <- scheduleEvent(sim, P(sim)$decadesMatchingYears[sim$.schedulingCounter], "waterfowl", "predictions")
      sim <- scheduleEvent(sim, P(sim)$decadesMatchingYears[sim$.schedulingCounter], "waterfowl", "updateScheduler")
    },
    modelsPrep = {
      # Prepare models: download data
      downloadWaterfowlModels(speciesURL = sim$speciesURL,
                              dataFolder = dataPath(sim))
    },
    climateDataPrep = {
      # Prepare climate data
      sim$bioCov <- prepareBioCov(Decade = P(sim)$predictionsDecades[which(
                        P(sim)$decadesMatchingYears[sim$.schedulingCounter] == time(sim))],
                        climateModel = "CCSM4",
                        pathInputs = dataPath(sim),
                        studyArea = sim$studyArea,
                        climateURLTable = sim$climateURLTable,
                        rasterToMatch = sim$rasterToMatch) #TODO Add more climate scenarios
      
      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$decadesMatchingYears[sim$.schedulingCounter + 1], "birdsNWT", "climateDataPrep")
    },
    simulatedLayersPrep = {
      # Prepare simulated layers
      if (!is.null(sim$cohortData)){
        mod$cohortData <- sim$cohortData
      } else {
        mod$cohortData <- createModObject(data = "cohortData", sim = sim,
                                          pathInput = inputPath(sim), currentTime = time(sim))
      }
      
      if (!is.null(sim$pixelGroupMap)){
        mod$pixelGroupMap <- sim$pixelGroupMap
      } else {
        mod$pixelGroupMap <- createModObject(data = "pixelGroupMap", sim = sim, 
                                             pathInput = inputPath(sim), currentTime = time(sim))
      }
      
      if (!is.null(sim$simulatedBiomassMap)){
        mod$simulatedBiomassMap <- sim$simulatedBiomassMap
      } else {
        mod$simulatedBiomassMap <- createModObject(data = "simulatedBiomassMap", sim = sim, 
                                                   pathInput = inputPath(sim), currentTime = time(sim))
      }
      if (any(is.null(mod$pixelGroupMap), is.null(mod$cohortData), 
              is.null(mod$simulatedBiomassMap))) {
        stop(paste0("This module can't run without data. ",
                    "Please make sure you are providing the correct path to ",
                    "the cohortData and pixelGroupMap saved objects, and the ",
                    "correct years for those in the parameter ",
                    "'decadesMatchingYears'"))
      }
      sim$treeCov <-  prepareTreeCov(cohortData = mod$cohortData,
                                     pixelGroupMap = mod$pixelGroupMap,
                                     biomassMap = mod$simulatedBiomassMap,
                                     rasterToMatch = sim$rasterToMatch)
      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$decadesMatchingYears[sim$.schedulingCounter + 1], "birdsNWT", "simulatedLayersPrep")
    },
    predictions = {
      sim$waterfowlPredictions <- lapply(sim$species, function(sp){ 
        # TODO Use future for parallelizing
        filePath <- file.path(dataPath(sim), sp)
        pred <- predictBiomod(sp = sp, 
                              inputPath = filePath,
                              bioCov = sim$bioCov,
                              currentTime = time(sim),
                              overwritePred = P(sim)$overwritePred,
                              treeCov = sim$treeCov)
        return(pred)
      })
      names(sim$waterfowlPredictions) <- sim$species
      
      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$decadesMatchingYears[sim$.schedulingCounter + 1], "birdsNWT", "predictions")
    },
    updateScheduler = {
      # First you schedule the event as the previous ones
      sim <- scheduleEvent(sim, P(sim)$decadesMatchingYears[sim$.schedulingCounter + 1], "waterfowl", "updateScheduler")
      # Then you update the counter
      sim$.schedulingCounter <- sim$.schedulingCounter +1
},
  warning(paste("Undefined event type: \'", current(sim)[1, "eventType", with = FALSE],
                  "\' in module \'", current(sim)[1, "moduleName", with = FALSE], "\'", sep = ""))
  )
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")
  if (!suppliedElsewhere("species", sim = sim)){
    warning("species was not provided, defaulting to c('REDH','BWTE','BUFF')", 
            immediate. = TRUE)
    sim$species <- c("REDH", "BWTE", "BUFF")
  }
  if (!suppliedElsewhere("speciesURL", sim = sim)){
    sim$speciesURL <- data.table(species = c("REDH", "BWTE", "BUFF"),
                             URL = c("https://drive.google.com/file/d/1JZ6ihDdn_WCSK05SIYrVqDq-YhRoGirj/view?usp=sharing",
                                     "https://drive.google.com/file/d/1BruuNFQlEw1t9IrYgN_OP3j92rF3ObcZ/view?usp=sharing",
                                     "https://drive.google.com/file/d/1escVgr4I15iizMYa47M_vBZIKJTunnFt/view?usp=sharing"))
    sim$speciesURL <- sim$speciesURL[species %in% sim$species, ]
    if (NROW(sim$speciesURL) == 0) stop("None of the provided species have models")
    unavail <- setdiff(sim$speciesURL[["species"]], sim$species)
    if (!length(unavail) == 0) warning(paste0("The following species are not available: ", 
                                              unavail), 
                                       immediate. = TRUE)
  }
  if (!suppliedElsewhere("studyArea", sim = sim, where = "sim")){
    if (quickPlot::isRstudioServer()) options(httr_oob_default = TRUE)
    
    message("No specific study area was provided. Croping to the Edehzhie Indigenous Protected Area (Southern NWT)")
    Edehzhie.url <- "https://drive.google.com/open?id=1klq0nhtFJZv47iZVG8_NwcVebbimP8yT"
    sim$studyArea <- Cache(prepInputs,
                           url = Edehzhie.url,
                           destinationPath = inputPath(sim),
                           omitArgs = c("destinationPath"))
  }
  
  if (!suppliedElsewhere("rasterToMatch", sim = sim, where = "sim")){
    sim$rasterToMatch <- Cache(prepInputs, url = "https://drive.google.com/open?id=1fo08FMACr_aTV03lteQ7KsaoN9xGx1Df", 
                               studyArea = sim$studyArea,
                               targetFile = "RTM.tif", destinationPath = inputPath(sim),
                               filename2 = NULL,
                               omitArgs = c("destinationPath", "filename2"))
  }
  if (!suppliedElsewhere("climateURLTable", sim = sim, where = "sim")){
    sim$climateURLTable <- data.table(decade = c(2050, 2070),
                                               URL = c("http://biogeo.ucdavis.edu/data/climate/cmip5/5m/cc85bi50.zip", 
                                                       "http://biogeo.ucdavis.edu/data/climate/cmip5/5m/cc85bi70.zip"),
                                               climmod = c("CCSM4", 
                                                           "CCSM4"))
  }
  
  return(invisible(sim))
}

### add additional events as needed by copy/pasting from above
