addPredictor = function(model, newdata, predictor) {
  xpredictors = predict(model, newdata=newdata, type="terms")
  colnames(xpredictors) = gsub("[\\(\\) \\*]*", "", gsub(",.*", "", colnames(xpredictors)))
  return(xpredictors[, predictor])
}


writeConfigFile = function(model, scenario, fishing, period, path, 
                           multiplier=NULL, nltl=4, ndt=24, test=FALSE) {
  
  tmp0 = "osmose.configuration.calibration;../../../LIB/osmose/parameters/calibration-parameters.csv\n
  osmose.configuration.main;../../../LIB/osmose/parameters/main-parameters_%s.csv\n
  osmose.configuration.output;../../../LIB/osmose/parameters/output-parameters_fishmip.csv\n
  osmose.configuration.fishing;../../../LIB/osmose/parameters/%s-parameters%s.csv\n
  osmose.configuration.simulation;../../../LIB/osmose/parameters/simulation-parameters.csv\n
  simulation.restart.file;../../../LIB/osmose/initial_conditions/%s.nc\n
  ltl.netcdf.file;../../../LIB/osmose/%s_%s_ltl-osmose_humboldt-n_15days_%s.nc\n
  plankton.multiplier.plk0;%s\n
  plankton.multiplier.plk1;%s\n
  plankton.multiplier.plk2;%s\n
  plankton.multiplier.plk3;%s\n
  simulation.time.nyear;%s\n
  ltl.nstep;%s\n
  simulation.restart.spinup;%s\n"
  
  tmp1 = "osmose.configuration.calibration;../../../LIB/osmose/parameters/calibration-parameters.csv\n
  osmose.configuration.main;../../../LIB/osmose/parameters/main-parameters_%s.csv\n
  osmose.configuration.output;../../../LIB/osmose/parameters/output-parameters_test.csv\n
  osmose.configuration.fishing;../../../LIB/osmose/parameters/%s-parameters%s.csv\n
  osmose.configuration.simulation;../../../LIB/osmose/parameters/simulation-parameters_short.csv\n
  simulation.restart.file;../../../LIB/osmose/initial_conditions/%s.nc\n
  ltl.netcdf.file;../../../LIB/osmose/%s_%s_ltl-osmose_humboldt-n_15days_%s.nc\n
  plankton.multiplier.plk0;%s\n
  plankton.multiplier.plk1;%s\n
  plankton.multiplier.plk2;%s\n
  plankton.multiplier.plk3;%s\n
  simulation.time.nyear;%s\n
  ltl.nstep;%s\n
  simulation.restart.spinup;%s\n"
  
  
  tmp = ifelse(!isTRUE(test), tmp0, tmp1) 
  
  if(is.null(multiplier)) multiplier = rep(1, nltl)
  
  F = ifelse(scenario=="historical", "historical", "future")
  F = ifelse(fishing=="fishing", paste0("_",F), "")
  start = as.numeric(unlist(strsplit(period, "_"))[1]) 
  restart = sprintf("%s_%s_%s_restart", model, fishing, start)
  restart = ifelse(scenario=="historical", "hindcast_restart", restart)
  
  years = as.list(setNames(as.numeric(unlist(strsplit(period, "_"))), nm=c("start", "end")))
  years$T = years$end - years$start + 1
  years$nstep = ndt*ifelse(scenario=="climatological", 1, years$T)
  
  periodX = ifelse(scenario=="historical", "historical", "future")
  
  txt = sprintf(tmp, periodX, fishing, F, restart, model, scenario, period,
                multiplier[1],multiplier[2],multiplier[3],multiplier[4], 
                years$T, years$nstep, years$T)
  
  fileName0 = sprintf("%s_%s_%s_%s.csv", model, scenario, fishing, period)
  fileName1 = sprintf("%s_%s_%s_%s-test.csv", model, scenario, fishing, period)
  
  fileName = ifelse(!isTRUE(test), fileName0, fileName1)
  fileName = file.path(path, tolower(fileName))
  
  cat(txt, file = fileName)
  return(invisible(fileName))
}


createRestart = function(input, output) {
  confName = basename(input)
  restartName = .restartName(confName)
  input = file.path(input, "restart")
  files = dir(path=input)
  file.copy(from=file.path(input, files), to=output)
  invisible(lapply(file.path(output, files), FUN=.initialRestart))
  newFiles = paste0(restartName, ".", seq_along(files) - 1)
  out = file.rename(from=file.path(output, files), 
                    to=file.path(output, newFiles))
  return(invisible(out))
}


.restartName = function(confName) {
  temp = unlist(strsplit(confName, "_"))
  model = temp[1]
  fishing = temp[3]
  start = as.numeric(temp[5]) + 1
  text = sprintf("%s_%s_%s_restart.nc", model, fishing, start) 
  return(text)
}

.initialRestart = function(file) {
  nc = open.ncdf(file, write=TRUE)
  att.put.ncdf(nc, varid=0, attname="step", attval="-1", prec="text")
  close.ncdf(nc)  
  return(invisible())
}

monthlyMean = function(x) {
  out = 0.5*(x[,,,c(TRUE, FALSE)] + x[,,,c(FALSE, TRUE)])
  return(out)
}


reduceSpatializeOutputs = function(path, output, prefix=NULL, clean=FALSE, 
                                   compression=9, prec="float", 
                                   missval=1e+20) {
  gc(verbose=FALSE)
  file = "run_spatialized_Simu0.nc"
  fileNew = paste0(basename(path), ".nc4")
  if(!is.null(fileNew)) fileNew = paste(prefix, fileNew, sep="_")
  
  spatial = nc_open(file.path(path, file))
  biomass = monthlyMean(ncvar_get(spatial, "biomass"))
  ltl     = monthlyMean(ncvar_get(spatial, "ltl_biomass"))
  yield   = monthlyMean(ncvar_get(spatial, "yield"))
  lat     = ncvar_get(spatial, "latitude")[1,]
  lon     = ncvar_get(spatial, "longitude")[,1]
  time    = spatial$dim$time$vals
  species = spatial$dim$species$vals
  ltl.sp  = spatial$dim$ltl$vals
  
  nc_close(spatial)
  
  dimLon  = ncdim_def("lon", "degrees", vals=lon)
  dimLat  = ncdim_def("lat", "degrees", vals=lat)
  dimTime = ncdim_def("time", "month", vals=seq_len(dim(biomass)[4]))
  dimSp   = ncdim_def("species", "", vals=species)
  dimLtl  = ncdim_def("ltl", "", vals=ltl.sp)
  
  BIO = ncvar_def(name="biomass", units="tonnes", dim=list(dimLon, dimLat, dimSp, dimTime), 
                  missval=missval, longname="Biomass", prec=prec,
                  compression=compression)
  
  LTL   = ncvar_def(name="ltl_biomass", units="tonnes", dim=list(dimLon, dimLat, dimLtl, dimTime), 
                    missval=missval, longname="LTL biomass", prec=prec,
                    compression=compression)
  
  YLD   = ncvar_def(name="yield", units="tonnes", dim=list(dimLon, dimLat, dimSp, dimTime), 
                    missval=missval, longname="Yield", prec=prec,
                    compression=compression)
  
  if(!file.exists(output)) dir.create(output, recursive=TRUE)
  
  ncNew = nc_create(filename=file.path(output, fileNew), vars=list(BIO, LTL, YLD))
  ncvar_put(ncNew, BIO, biomass) 
  ncvar_put(ncNew, LTL, ltl) 
  ncvar_put(ncNew, YLD, yield) 
  nc_close(ncNew)
  
  out = NULL
  if(isTRUE(clean)) out=file.remove(file.path(path, file))
  
  gc(verbose=FALSE)
  return(invisible(out))
}


#' @param file Name of the ncdf file with osmose outputs. The variables
#' must be named 'biomass', ltl_biomass' and 'yield' and follow
#' OSMOSE conventions.
#' @param output Output folder to save FISH-MIP processed outputs.
#' @param pattern C-style string formating for FISH-MIP outputs with
#' just one '%s' for variable name. 
#' @return Write FISH-MIP formatted outputs in the 'output' folder.
fish_mip = function(file, info, output, mass2C=0.1, 
                    compression=9, prec="float", missval=1e+20) {
  
  
  xnames = unlist(strsplit(gsub(pattern="\\.nc[0-9]$", "", basename(file)), split="_"))
  names(xnames) = c("model", "gcm", "scenario", "fishing", "start", "end")
  info = c(info, as.list(xnames))
  info$period = paste(info$start, info$end, sep="_")
  
  pattern = with(info, paste(model, gcm, scenario, diaz, fishing, oa, "%s",
                             region, temp, period, sep="_"))
  
  DateStamp("Processing file", sQuote(file))
  path = file.path(output, toupper(info$gcm))
  
  if(!file.exists(path)) dir.create(path, recursive = TRUE)
  
  spatial = nc_open(file)
  
  biomass = ncvar_get(spatial, "biomass")
  ltl     = ncvar_get(spatial, "ltl_biomass")
  yield   = ncvar_get(spatial, "yield")
  lat     = ncvar_get(spatial, "lat")
  lon     = ncvar_get(spatial, "lon")
  
  grid = createGridAxes(lat = range(lat), lon=range(lon), 
                        dx=abs(mean(diff(lon))), 
                        dy=abs(mean(diff(lat))))
  
  area = as.numeric(grid$area) # area in km^2
  
  #mandatory
  tsb = mass2C*(apply(biomass, c(1,2,4), sum) + apply(ltl, c(1,2,4), sum))/area
  tcb = mass2C*(apply(biomass, c(1,2,4), sum) + apply(ltl[,,info$ltlTL>1,], c(1,2,4), sum))/area
  b10cm = mass2C*(apply(biomass[,,info$Linf>=10,], c(1,2,4), sum))/area
  b30cm = mass2C*(apply(biomass[,,info$Linf>=30,], c(1,2,4), sum))/area
  tc = apply(yield, c(1,2,4), sum)/area
  tla = apply(yield, c(1,2,4), sum)/area
  #optional
  bcom = mass2C*(apply(biomass[,, info$isCommercial,], c(1,2,4), sum))/area
  
  # create dimensions
  dimLon  = ncdim_def("lon", "degrees", vals=lon)
  dimLat  = ncdim_def("lat", "degrees", vals=lat)
  dimTime = ncdim_def("time", "month", seq_len(dim(tsb)[3]))
  
  # create variables
  units1 = "gC/ m^2"
  units2 = "g wet biomass/m^2"
  
  TSB   = ncvar_def(name="tsb", units=units1, dim=list(dimLon, dimLat, dimTime), 
                    missval=missval, longname="Total system biomass density", prec=prec,
                    compression=compression)
  TCB   = ncvar_def(name="tcb", units=units1, dim=list(dimLon, dimLat, dimTime), 
                    missval=missval, longname="Total consumer biomass density", prec=prec,
                    compression=compression)
  B10CM = ncvar_def(name="b10cm", units=units1, dim=list(dimLon, dimLat, dimTime), 
                    missval=missval, longname="Biomass density of consumers >10cm", prec=prec,
                    compression=compression)
  B30CM = ncvar_def(name="b30cm", units=units1, dim=list(dimLon, dimLat, dimTime), 
                    missval=missval, longname="Biomass density of consumers >30cm", prec=prec,
                    compression=compression)
  TC    = ncvar_def(name="tc", units=units2, dim=list(dimLon, dimLat, dimTime), 
                    missval=missval, longname="Total catch", prec=prec,
                    compression=compression)
  TLA   = ncvar_def(name="tla", units=units2, dim=list(dimLon, dimLat, dimTime), 
                    missval=missval, longname="Total landings", prec=prec,
                    compression=compression)
  BCOM  = ncvar_def(name="bcom", units=units1, dim=list(dimLon, dimLat, dimTime), 
                    missval=missval, longname="Biomass density of commercial species", prec=prec,
                    compression=compression)
  
  
  vars = c("tsb", "tcb", "b10cm", "b30cm", "tc", "tla", "bcom")
  
  csv = list()
  for(var in vars) {
    csv[[var]] = apply(get(var), 3, mean, na.rm=TRUE)
  }
  names(csv) = vars
  csv = as.data.frame(csv)
  fileNew = file.path(output, sprintf(paste0(pattern, ".csv"), "all"))
  write.csv(csv, file=fileNew)
  
  for(var in vars) {
    fileNew = file.path(path, sprintf(paste0(pattern, ".nc4"), var))
    ncNew = nc_create(filename=fileNew, vars=get(toupper(var)), force_v4=TRUE)
    ncvar_put(ncNew, var, get(var)) 
    nc_close(ncNew)
  }
  
  nc_close(spatial)
  
  return(invisible())
}


regridOsmose = function(file, output, newGrid, 
                        compression=9, prec="float", missval=1e+20) {
  
  nc = nc_open(file)
  lat = ncvar_get(nc, "lat")
  lon = ncvar_get(nc, "lon")
  old = createGridAxes(lat=range(lat), lon=range(lon),
                       dx=mean(abs(diff(lon))),
                       dy=mean(abs(diff(lat))))
  
  time    = ncvar_get(nc, "time")
  species = ncvar_get(nc, "species")
  ltl.sp  = ncvar_get(nc, "ltl")
  
  bio = ncvar_get(nc, "biomass")/as.numeric(old$area)
  ltl = ncvar_get(nc, "ltl_biomass")/as.numeric(old$area)
  yld = ncvar_get(nc, "yield")/as.numeric(old$area)
  
  newBio = regrid(bio, new=newGrid, old=old)*as.numeric(newGrid$area)
  newLtl = regrid(ltl, new=newGrid, old=old)*as.numeric(newGrid$area)
  newYld = regrid(yld, new=newGrid, old=old)*as.numeric(newGrid$area)
  
  dimLon  = ncdim_def("lon", "degrees", vals=newGrid$lon)
  dimLat  = ncdim_def("lat", "degrees", vals=newGrid$lat)
  dimTime = ncdim_def("time", "month", vals=seq_len(dim(bio)[4]))
  dimSp   = ncdim_def("species", "", vals=species)
  dimLtl  = ncdim_def("ltl", "", vals=ltl.sp)
  
  BIO = ncvar_def(name="biomass", units="tonnes", dim=list(dimLon, dimLat, dimSp, dimTime), 
                  missval=missval, longname="Biomass", prec=prec,
                  compression=compression)
  
  LTL   = ncvar_def(name="ltl_biomass", units="tonnes", dim=list(dimLon, dimLat, dimLtl, dimTime), 
                    missval=missval, longname="LTL biomass", prec=prec,
                    compression=compression)
  
  YLD   = ncvar_def(name="yield", units="tonnes", dim=list(dimLon, dimLat, dimSp, dimTime), 
                    missval=missval, longname="Yield", prec=prec,
                    compression=compression)
  
  ncNew = nc_create(filename=output, vars=list(BIO, LTL, YLD))
  ncvar_put(ncNew, BIO, newBio) 
  ncvar_put(ncNew, LTL, newLtl) 
  ncvar_put(ncNew, YLD, newYld) 
  
  nc_close(ncNew)
  nc_close(nc)
  
  return(invisible())
}

readMask = function(maskFile) {
  nc = nc_open(maskFile)
  mask = ncvar_get(nc, "mask")
  lat = ncvar_get(nc, "lat")
  lon = ncvar_get(nc, "lon")
  old = createGridAxes(lat=range(lat), lon=range(lon),
                       dx=mean(abs(diff(lon))),
                       dy=mean(abs(diff(lat))))
  old$mask = mask
  nc_close(nc)  
  return(old)
}



.getVarSize = function(file) {
  if(length(file)==0) return(NULL)
  nc = nc_open(file)
  out = nc$var[[1]]$varsize 
  nc_close(nc)
  return(out)
}


.checkVarSize = function(files) {
  sizes = lapply(files, .getVarSize)
  sizes = sizes[!sapply(sizes, FUN=is.null)]
  OK = all(sapply(sizes, FUN=identical, y=sizes[[1]]))
  if(!isTRUE(OK)) stop("LTL dimensions don't match")
  return(sizes[[1]])
}


.readLTLData = function(file, dim) {
  if(length(file)==0) return(array(dim=dim))
  xtemp = nc_open(file)
  out = ncvar_get(xtemp, names(xtemp$var))
  return(out)
}

.readLTLDimensions = function(files) {
  files = files[sapply(files, FUN=length)>0]
  xtemp = nc_open(files[[1]])
  dimNames = names(xtemp$dim) 
  dims = list()
  for(dimName in dimNames) {
    dims[[dimName]] = ncvar_get(xtemp, dimName)    
  }
  return(dims)
}


.getLTLFiles = function(input, pattern, ltlNames) {
  
  .getLTLFiles1 = function(input, pattern, ltlNames, res) {
    files = dir(path=input, pattern=pattern)
    files = grep(patt=res, x=files, value = TRUE)
    files = file.path(input, files)
    
    files = list(sphy=grep(patt=ltlNames[1], x=files, value = TRUE),
                 lphy=grep(patt=ltlNames[2], x=files, value = TRUE),
                 szoo=grep(patt=ltlNames[3], x=files, value = TRUE),
                 lzoo=grep(patt=ltlNames[4], x=files, value = TRUE))
    
    return(files)
  }
  
  byMonth = .getLTLFiles1(input, pattern, ltlNames, res="monthly")
  byYear  = .getLTLFiles1(input, pattern, ltlNames, res="annual")
  bm = sum(sapply(byMonth, FUN=length)>0)
  by = sum(sapply(byYear, FUN=length)>0)
  out = if(bm>=by) byMonth else byYear
  return(out)
  
}


mergeLTLData = function(input, output, varName="ltl", cf=NULL) {
  
  ltlNames = c("sphy", "lphy", "szoo", "lzoo")
  pattern = "_zint_"
  
  if(is.null(cf)) cf = as.list(setNames(rep(1,length(ltlNames)), nm=ltlNames))
  if(length(cf)==1) cf = as.list(setNames(rep(as.numeric(cf),length(ltlNames)), nm=ltlNames))
  
  # get files
  files = .getLTLFiles(input, pattern, ltlNames)
  
  size = .checkVarSize(files)
  dims = .readLTLDimensions(files)
  
  ind = which(sapply(files, FUN=length)>0)[1]
  outputFile = basename(gsub(pattern=ltlNames[ind], replacement=varName,
                             x=files[[ind]])) 
  outputFile = file.path(output, outputFile)
  
  LTL = list()
  for(iFile in seq_along(files)) {
    plkName = names(files)[iFile]
    LTL[[plkName]] = .readLTLData(file=files[[iFile]], dim=size)
  }
  
  dimLon  = ncdim_def("lon", "degrees", vals=dims[[1]])
  dimLat  = ncdim_def("lat", "degrees", vals=dims[[2]])
  dimTime = ncdim_def("time", "month", vals=dims[[3]])
  
  SPHY   = ncvar_def(name="sphy", units="mol C/m^2", dim=list(dimLon, dimLat, dimTime), 
                     missval=1e+20, longname="Small phytoplankton carbon concentration", prec="float",
                     compression=9)
  
  LPHY   = ncvar_def(name="lphy", units="mol C/m^2", dim=list(dimLon, dimLat, dimTime), 
                     missval=1e+20, longname="Large phytoplankton (diatoms) carbon concentration", 
                     prec="float", compression=9)
  
  SZOO   = ncvar_def(name="szoo", units="mol C/m^2", dim=list(dimLon, dimLat, dimTime), 
                     missval=1e+20, longname="small zooplankton carbon concentration", prec="float",
                     compression=9)
  
  LZOO   = ncvar_def(name="lzoo", units="mol C/m^2", dim=list(dimLon, dimLat, dimTime), 
                     missval=1e+20, longname="meso zooplankton carbon concentration", prec="float",
                     compression=9)
  
  if(!file.exists(output)) dir.create(output, recursive=TRUE)
  
  ncNew = nc_create(filename=outputFile, vars=list(SPHY, LPHY, SZOO, LZOO))
  ncvar_put(ncNew, SPHY, cf$sphy*LTL$sphy) 
  ncvar_put(ncNew, LPHY, cf$lphy*LTL$lphy) 
  ncvar_put(ncNew, SZOO, cf$szoo*LTL$szoo) 
  ncvar_put(ncNew, LZOO, cf$lzoo*LTL$lzoo) 
  nc_close(ncNew)
  
  return(invisible())
}

.capitalize = function(x) {
  n = nchar(x)
  first = substr(x, 1, 1)
  remainder = substr(x, 2, n)
  out = paste0(toupper(first), remainder)
  return(out)
}

# downscaling -------------------------------------------------------------


createDownscalingData = function(mask, model) {
  
  mbase = data.frame(lat=as.numeric(mask$LAT),
                     lon=as.numeric(mask$LON),
                     mask=as.numeric(mask$mask),
                     area=as.numeric(mask$area))
  
  coast = model$input$coast
  shelf = model$input$shelf
  
  xdist = getSignedDistance(data=mbase, ref=shelf, abs=coast)
  mbase$dc    = xdist$abs
  mbase$shelf = xdist$dist
  
  models = model$models
  
  for(iModel in names(models)) {
    mod = models[[iModel]]
    proxyName = paste0("proxy", .capitalize(iModel))
    mbase[, proxyName] = addPredictor(mod, mbase, "teshelf")
  }
  
  return(mbase)
}

loadShelfModel = function(file) {
  cat("Loading", file, "...\n")
  load(file)
  objs = ls()
  ind = sapply(mget(objs), inherits, what = "ShelfModel")
  n = length(which(ind))
  msg1 = "The file %s does not contain a valid shelf model."
  msg2 = "The file %s contains multiple shelf models, only one allowed."
  if(n==0) stop(sprintf(msg1, file))
  if(n>1)  stop(sprintf(msg2, file))
  output = get(x = objs[which(ind)])
  return(output)
}

loadMultipliers = function(file) {
  cat("Loading", file, "...\n")
  load(file)
  objs = ls()
  ind = sapply(mget(objs), inherits, what = "plk.multiplier")
  n = length(which(ind))
  msg1 = "The file %s does not contain valid plankton multiplier."
  msg2 = "The file %s contains multiple set of plankton multipliers, only one allowed."
  if(n==0) stop(sprintf(msg1, file))
  if(n>1)  stop(sprintf(msg2, file))
  output = get(x = objs[which(ind)])
  return(output)
}

downscale = function(object, old, new, proxy, thr, ...) {
  UseMethod("downscale")
}

downscale.matrix = function(object, old, new, proxy, ...) {
  
  stopifnot(exists("lat", where=old), exists("lon", where=old))
  stopifnot(exists("lat", where=new), exists("lon", where=new)) 
  stopifnot(exists(proxy, where=old), exists(proxy, where=new)) 
  
  nlat = length(unique(new$lat))
  nlon = length(unique(new$lon))
  
  newDim = c(nlon, nlat)
  
  if(all(is.na(object))) return(array(dim=newDim))
  
  old$ivar = as.numeric(object)
  old$iproxy = old[, proxy]
  
  new$iproxy = new[, proxy]
  
  iMod = gam(ivar ~ te(lat,lon) + iproxy, data=old, family=gaussian(log))
  
  out = predict(iMod, newdata=new, type="response")*new$mask 
  
  newmap = matrix(out, ncol=nlat, nrow=nlon)
  
  return(newmap)
  
}


downscale.array = function(object, old, new, proxy, thr=NULL, ...) {
  
  stopifnot(exists("lat", where=old), exists("lon", where=old))
  stopifnot(exists("lat", where=new), exists("lon", where=new)) 
  stopifnot(exists(proxy, where=old), exists(proxy, where=new)) 
  
  nlat = length(unique(new$lat))
  nlon = length(unique(new$lon))
  
  newDim = c(nlon, nlat, dim(object)[-c(1,2)])
  
  if(all(is.na(object))) return(array(dim=newDim))
  
  ndim = seq_along(dim(object))[-c(1,2)]
  newmap = apply(object, ndim, FUN=downscale, old=old, new=new, proxy=proxy, ...)
  
  dim(newmap) = newDim
  
  if(!is.null(thr)) {
    newmap[newmap>thr] = NA
    newmap = interpolateMap(newmap, anomalies=TRUE)    
  }
  
  return(newmap)
  
}


downscaleLTLData = function(file, old, new, model, outputFile,
                            thr=NULL) {
  
  vars = names(model$models)
  
  nc   = nc_open(file)
  ncvars = names(nc$var)
  
  size = .checkVarSize(file)
  stopifnot(identical(size[1:2], dim(old$mask)))
  
  dims = .readLTLDimensions(file)
  newDims = c(dim(new$mask), length(dims[[3]]))
  
  newBase = createDownscalingData(mask=new, model=model)
  oldBase = createDownscalingData(mask=old, model=model)
  
  LTL = list()
  
  for(ivar in vars) {
    if(ivar %in% ncvars) {
      DateStamp("\t\tDownscaling for", ivar)
      iproxy = paste0("proxy", .capitalize(ivar))
      xvar = ncvar_get(nc, ivar)    
      LTL[[ivar]] = downscale(xvar, old=oldBase, new=newBase, proxy=iproxy, thr=thr[[ivar]])      
    } else LTL[[ivar]] = array(dim=newDims)
  }
  
  dimLon  = ncdim_def("lon", "degrees", vals=new$lon)
  dimLat  = ncdim_def("lat", "degrees", vals=new$lat)
  dimTime = ncdim_def("time", "month", vals=dims[[3]])
  
  # generalize
  SPHY   = ncvar_def(name="sphy", units="mol C/m^2", dim=list(dimLon, dimLat, dimTime), 
                     missval=1e+20, longname="Small phytoplankton carbon concentration", prec="float",
                     compression=9)
  
  LPHY   = ncvar_def(name="lphy", units="mol C/m^2", dim=list(dimLon, dimLat, dimTime), 
                     missval=1e+20, longname="Large phytoplankton (diatoms) carbon concentration", 
                     prec="float", compression=9)
  
  SZOO   = ncvar_def(name="szoo", units="mol C/m^2", dim=list(dimLon, dimLat, dimTime), 
                     missval=1e+20, longname="small zooplankton carbon concentration", prec="float",
                     compression=9)
  
  LZOO   = ncvar_def(name="lzoo", units="mol C/m^2", dim=list(dimLon, dimLat, dimTime), 
                     missval=1e+20, longname="meso zooplankton carbon concentration", prec="float",
                     compression=9)
  
  if(!file.exists(dirname(outputFile))) 
    dir.create(dirname(outputFile), recursive=TRUE)
  
  ncNew = nc_create(filename=outputFile, vars=list(SPHY, LPHY, SZOO, LZOO))
  ncvar_put(ncNew, SPHY, LTL$sphy) 
  ncvar_put(ncNew, LPHY, LTL$lphy) 
  ncvar_put(ncNew, SZOO, LTL$szoo) 
  ncvar_put(ncNew, LZOO, LTL$lzoo) 
  nc_close(ncNew)
  
  return(invisible())
}

readLTL = function(file) {
  nc = nc_open(file)
  vars = names(nc$var)
  output = list()
  for(var in vars) {
    output[[var]] = ncvar_get(nc, var)
  }
  nc_close(nc)
  return(output)
}

readLTL2 = function(input, var="ltl") {
  scen = list()
  hrFile = dir(path=input, patt=sprintf("_%s-hr_", var))
  lrFile = dir(path=input, patt=sprintf("_%s_", var))
  scen[["hr"]] = readLTL(file.path(input, hrFile))
  scen[["lr"]] = readLTL(file.path(input, lrFile))
  return(scen)
}


plankton2Biomass = function(file, cf, pattern="plk", replacement="ltl") {
  
  opatt = paste0("_", pattern)
  npatt = paste0("_", replacement)
  
  outputFile = gsub(patt=opatt, rep=npatt, x=file)
  
  dims = .readLTLDimensions(file)
  LTL  = readLTL(file)
  
  stopifnot(all(names(LTL) %in% names(cf)), is.list(cf))
  
  grid = createGridAxes(lat = range(dims$lat), lon=range(dims$lon), 
                        dx=abs(mean(diff(dims$lon))), 
                        dy=abs(mean(diff(dims$lat))))
  
  area = grid$area # area in km^2
  
  .byArea = function(x) x*as.numeric(area)
  
  LTL = rapply(LTL, f=.byArea, classes="array", how="list")
  
  dimLon  = ncdim_def("lon", "degrees", vals=dims$lon)
  dimLat  = ncdim_def("lat", "degrees", vals=dims$lat)
  dimTime = ncdim_def("time", "month", vals=dims$time)
  
  # generalize
  SPHY   = ncvar_def(name="sphy", units="tonnes", dim=list(dimLon, dimLat, dimTime), 
                     missval=1e+20, longname="Small phytoplankton biomass", prec="float",
                     compression=9)
  
  LPHY   = ncvar_def(name="lphy", units="tonnes", dim=list(dimLon, dimLat, dimTime), 
                     missval=1e+20, longname="Large phytoplankton biomass", 
                     prec="float", compression=9)
  
  SZOO   = ncvar_def(name="szoo", units="tonnes", dim=list(dimLon, dimLat, dimTime), 
                     missval=1e+20, longname="small zooplankton biomass", prec="float",
                     compression=9)
  
  LZOO   = ncvar_def(name="lzoo", units="tonnes", dim=list(dimLon, dimLat, dimTime), 
                     missval=1e+20, longname="meso zooplankton biomass", prec="float",
                     compression=9)
  
  if(!file.exists(dirname(outputFile))) 
    dir.create(dirname(outputFile), recursive=TRUE)
  
  ncNew = nc_create(filename=outputFile, vars=list(SPHY, LPHY, SZOO, LZOO))
  ncvar_put(ncNew, SPHY, cf$sphy*LTL$sphy) 
  ncvar_put(ncNew, LPHY, cf$lphy*LTL$lphy) 
  ncvar_put(ncNew, SZOO, cf$szoo*LTL$szoo) 
  ncvar_put(ncNew, LZOO, cf$lzoo*LTL$lzoo) 
  nc_close(ncNew)
  
}

.guessNdt = function(file) {
  isMonthly = grepl(patt="monthly", x=file)
  isAnnual  = grepl(patt="annual", x=file)
  if(isMonthly) return(12)
  if(isAnnual) return(1)
  stop("Could not guess the number of time steps per year, add the ndt value")
}

.defaultCF = function(plk) {
  x = rep(1, length(plk))
  names(x) = plk
  return(as.list(x))
}

.getTLabel = function(freq) {
  stopifnot(is.numeric(freq))
  out = switch(as.character(freq), 
               "1"  = "years",
               "4"  = "quarter",
               "12" = "months",
               "24" = "15days",
               "52" = "1week",
               sprintf("around %d days :)", floor(365/freq))
  )
  return(out)
}

# coef for tonnes/km2 if density=TRUE
createOsmoseLTL = function(file, outputFile, freq, ndt=NULL, 
                           density=FALSE, cf=NULL, skip=NULL,
                           length=NULL) {
  
  dims = .readLTLDimensions(file)
  LTL  = readLTL(file)
  
  plk = names(LTL)
  
  tlab = .getTLabel(freq)
  
  if(is.null(cf)) cf = .defaultCF(plk)
  if(is.null(ndt)) ndt = .guessNdt(file)
  
  msg =  sprintf("freq must be a multiple of ndt=%d", ndt)
  if((freq/ndt)%%1 != 0) stop(msg)
  
  
  if(isTRUE(density)) {
    grid = createGridAxes(lat = range(dims$lat), lon=range(dims$lon), 
                          dx=abs(mean(diff(dims$lon))), 
                          dy=abs(mean(diff(dims$lat))))
    area = as.numeric(grid$area)
  } else area = 1
  
  ltlNew = array(dim=append(dim(LTL[[1]]), length(plk), 2))
  
  for(i in seq_along(plk)) 
    ltlNew[, , i, ] = cf[[plk[i]]]*LTL[[plk[i]]]*area
  
  ind = .getIndex(n=dim(ltlNew)[4], skip=skip, length=length)
  
  if(!is.null(ind)) ltlNew = ltlNew[, , , ind]
  dts    = rep(seq_len(dim(ltlNew)[4]), each=freq/ndt)
  ltlNew = ltlNew[, , , dts]
  
  dimLon  = dim.def.ncdf("nx", "degrees", seq_along(dims$lon))
  dimLat  = dim.def.ncdf("ny", "degrees", seq_along(dims$lat))
  dimLtl  = dim.def.ncdf("ltl", "group", seq_along(plk))
  dimTime = dim.def.ncdf("time", tlab, seq_along(dts))
  
  ltl_biomass = var.def.ncdf(name="ltl_biomass", units="tonnes", 
                             dim=list(dimLon, dimLat, dimLtl, dimTime), 
                             missval=1.00000001504747e+30, 
                             longname="LTL biomass integrated", 
                             prec="float")
  
  if(!file.exists(dirname(outputFile))) 
    dir.create(dirname(outputFile), recursive=TRUE)
  
  nc.new = create.ncdf(outputFile, list(ltl_biomass))
  put.var.ncdf(nc.new, "ltl_biomass", ltlNew)
  att.put.ncdf(nc.new, varid=0, attname="LTL names", attval=plk, prec="text")
  close.ncdf(nc.new)
  return(plk)
  
}

.getIndex = function(n, skip=NULL, length=NULL) {
  ind = seq_len(n)
  if(is.null(skip)&is.null(length)) return(NULL)
  if(!is.null(skip)) ind = tail(ind, -skip)
  if(is.null(length)) return(ind)
  if(length(ind)<length) stop("length argument greater than available steps")
  ind = ind[seq_len(length)]
  return(ind)
}

.doRunTest = function(logFile, multiplier) {
  if(!file.exists(logFile)) return(TRUE)
  lastMultiplier = as.numeric(unlist(
    strsplit(tail(readLines(logFile),1), " ")))
  return(!identical(lastMultiplier, multiplier))
}


runOsmoseTest = function(model, scenario, fishing, period, multiplier, input, output, jar) {
  
  config  = writeConfigFile(model=model, scenario=scenario, 
                            fishing=fishing, period=period, path=input, 
                            multiplier=multiplier, test=TRUE)
  
  DateStamp()
  confName = gsub("^[0-9]*_", "", gsub("\\.csv", "", basename(config)))
  confName = tolower(confName)
  cat(sprintf("\n\n\nRunning configuration %s\n\n\n", confName))
  
  logFile = file.path(output, "test", paste0(confName, ".log"))
  
  if(is.null(multiplier)) multiplier = rep(1, 4)
  
  doRun = .doRunTest(logFile, multiplier)
  
  if(!file.exists(dirname(logFile))) dir.create(dirname(logFile))
  if(!file.exists(logFile)) file.create(logFile)
  
  if(doRun) {
    runOsmose(osmose = jar, 
              input = file.path(input, config),
              output = file.path(output, "test", confName))
    cat(date(), "\n", file = logFile, append = TRUE)
    write(x = multiplier, file = logFile, append = TRUE)
    
  }
  
  out = osmose2R(file.path(output, "test", confName))
  out$multiplier = multiplier
  
  file.remove(config)
  
  return(out)
}


compareOsmoseBiomass = function(...) {
  
  old.par = par(no.readonly = TRUE)
  on.exit(par(old.par))
  
  models = list(...)
  stopifnot(all(sapply(models, inherits, what="osmose")))
  modelNames = as.character(match.call())[-1]
  names(models) = modelNames
  bio = lapply(models, getVar, var="biomass")
  dims = sapply(bio, dim)
  allMatch = all(dims[2,-1]==dims[2,1])
  stopifnot(allMatch)
  out = array(dim=c(apply(dims, 1, max), length(bio)))
  for(i in seq_along(bio)) {
    ind = seq_len(nrow(bio[[i]]))
    out[ind,,i] = bio[[i]]
  }
  
  dimnames(out)[[2]] = colnames(models[[1]])
  dimnames(out)[[3]] = modelNames
  par(mfrow=kali::getmfrow(dim(out)[2]), mar=c(3,3,1,1),
      oma=c(0,0,2,0))
  apply(out, 2, matplot, type="l", lty=1)
  title(main=paste(modelNames, collapse=" - "), outer=TRUE)
  return(invisible(out))
}

getOsmoseLtl = function(file) {
  nc = nc_open(file)
  ltl = ncvar_get(nc, "ltl_biomass")
  out = apply(ltl, 4:3, sum, na.rm=TRUE)
  gc(verbose=FALSE)
  return(out)
}
