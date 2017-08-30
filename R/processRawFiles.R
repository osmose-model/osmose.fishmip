
# Merge environmental data ------------------------------------------------

mergeEnvData = function(input, output, varName="env", 
                        pattern = "_zall_", cf=NULL) {
  
  varNames = c("to", "so", "o2", "ph")
  
  # get files
  files = osmose.fishmip:::.getEnvFiles(input, pattern, varNames)
  
  size = osmose.fishmip:::.checkVarSize(files)
  dims = osmose.fishmip:::.readNcdfDimensions(files)
  
  ind = which(sapply(files, FUN=length)>0)[1]
  outputFile = basename(gsub(pattern=paste("_", varNames[ind], "_", sep=""), 
                             replacement=paste("_", varName, "_", sep=""), 
                             x=files[[ind]])) 
  outputFile = gsub(pattern = pattern, replacement = "_", x=outputFile)
  outputFile = file.path(output, outputFile)
  
  ENV = list()
  for(iFile in seq_along(files)) {
    envName = names(files)[iFile]
    ENV[[envName]] = osmose.fishmip:::.readEnvData(file=files[[iFile]], dim=size)
  }
  
  dimLon  = ncdim_def("lon", "degrees", vals=dims[[1]])
  dimLat  = ncdim_def("lat", "degrees", vals=dims[[2]])
  dimTime = ncdim_def("time", "month", vals=dims[[4]])
  
  SST   = ncvar_def(name="sst", units="ºC", dim=list(dimLon, dimLat, dimTime), 
                    missval=1e+20, longname="Sea Surface Temperature", prec="float",
                    compression=9)
  
  SBT   = ncvar_def(name="sbt", units="ºC", dim=list(dimLon, dimLat, dimTime), 
                    missval=1e+20, longname="Sea Bottom Temperature", prec="float",
                    compression=9)
  
  DT15  = ncvar_def(name="dt15", units="meters", dim=list(dimLon, dimLat, dimTime), 
                    missval=1e+20, longname="Depth of the isotherm of 15ºC", prec="float",
                    compression=9)
  
  TC   = ncvar_def(name="tc", units="meters", dim=list(dimLon, dimLat, dimTime), 
                   missval=1e+20, longname="Depth of the thermocline", prec="float",
                   compression=9)
  
  SSS   = ncvar_def(name="sss", units="psu", dim=list(dimLon, dimLat, dimTime), 
                    missval=1e+20, longname="Sea Surface Salinity", 
                    prec="float", compression=9)
  
  DO2   = ncvar_def(name="do2", units="meters", dim=list(dimLon, dimLat, dimTime), 
                    missval=1e+20, longname="Depth of the 2 mL/L oxygen", prec="float",
                    compression=9)
  
  OC    = ncvar_def(name="oc", units="meters", dim=list(dimLon, dimLat, dimTime), 
                    missval=1e+20, longname="Depth of the oxycline", prec="float",
                    compression=9)
  
  PH    = ncvar_def(name="ph", units="", dim=list(dimLon, dimLat, dimTime), 
                    missval=1e+20, longname="pH", prec="float",
                    compression=9)
  
  if(!file.exists(output)) dir.create(output, recursive=TRUE)
  
  ncNew = nc_create(filename=outputFile, 
                    vars=list(SST, SBT, DT15, TC, SSS, DO2, OC, PH))
  
  ENV$to = ENV$to - 273.15 # K -> ºC
  ENV$o2 = 1000*16*1e-3*(5.3/7.6)*ENV$o2 # mol/m3 -> mL/L
  
  DateStamp("Adding SST...\n")
  ncvar_put(ncNew, SST, ENV$to[,,1,]) 
  # ncvar_put(ncNew, SBT, osmose.fishmip:::.getLast(ENV$to))
  # ncvar_put(ncNew, DT15, calculateIsoline(ENV$to, dims$depth, ref=15))
  # ncvar_put(ncNew, TC, calculateCline(ENV$to, dims$depth))
  DateStamp("Adding SSS...\n")
  ncvar_put(ncNew, SSS, ENV$so[,,1,])
  DateStamp("Adding oxycline...\n")
  ncvar_put(ncNew, DO2, calculateIsoline(ENV$o2, dims$depth, 
                                         ref=2, k=2, xtr=c(+1, -1))) 
  # ncvar_put(ncNew, OC, calculateOxycline(ENV$o2, dims$depth))
  DateStamp("Adding PH...\n")
  ncvar_put(ncNew, PH, ENV$ph[,,1,])
  
  nc_close(ncNew)
  
  return(invisible())
}

mergeEnvData2 = function(input, output, varName="env", 
                        pattern = "_zs_", cf=NULL) {
  
  varNames = c("to", "so", "o2", "ph")
  
  # get files
  files = osmose.fishmip:::.getEnvFiles(input, pattern, varNames)
  
  size = osmose.fishmip:::.checkVarSize(files)
  dims = osmose.fishmip:::.readNcdfDimensions(files)
  
  ind = which(sapply(files, FUN=length)>0)[1]
  outputFile = basename(gsub(pattern=paste("_", varNames[ind], "_", sep=""), 
                             replacement=paste("_", varName, "_", sep=""), 
                             x=files[[ind]])) 
  outputFile = gsub(pattern = pattern, replacement = "_", x=outputFile)
  outputFile = file.path(output, outputFile)
  
  ENV = list()
  for(iFile in seq_along(files)) {
    envName = names(files)[iFile]
    ENV[[envName]] = osmose.fishmip:::.readEnvData(file=files[[iFile]], dim=size)
  }
  
  dimLon  = ncdim_def("lon", "degrees", vals=dims[[1]])
  dimLat  = ncdim_def("lat", "degrees", vals=dims[[2]])
  dimTime = ncdim_def("time", "month", vals=dims[[4]])
  
  SST   = ncvar_def(name="sst", units="ºC", dim=list(dimLon, dimLat, dimTime), 
                    missval=1e+20, longname="Sea Surface Temperature", prec="float",
                    compression=9)
  
  # SBT   = ncvar_def(name="sbt", units="ºC", dim=list(dimLon, dimLat, dimTime), 
  #                   missval=1e+20, longname="Sea Bottom Temperature", prec="float",
  #                   compression=9)
  

  SSS   = ncvar_def(name="sss", units="psu", dim=list(dimLon, dimLat, dimTime), 
                    missval=1e+20, longname="Sea Surface Salinity", 
                    prec="float", compression=9)
  
  O2   = ncvar_def(name="o2", units="mL/L", dim=list(dimLon, dimLat, dimTime), 
                    missval=1e+20, longname="Depth of the 2 mL/L oxygen", prec="float",
                    compression=9)
  
  PH    = ncvar_def(name="ph", units="", dim=list(dimLon, dimLat, dimTime), 
                    missval=1e+20, longname="pH", prec="float",
                    compression=9)
  
  if(!file.exists(output)) dir.create(output, recursive=TRUE)
  
  ncNew = nc_create(filename=outputFile, 
                    vars=list(SST, SSS, O2, PH))
  
  ENV$o2 = 1000*16*1e-3*(5.3/7.6)*ENV$o2 # mol/m3 -> mL/L
  
  DateStamp("Adding SST...\n")
  ncvar_put(ncNew, SST, ENV$to) 
  DateStamp("Adding SSS...\n")
  ncvar_put(ncNew, SSS, ENV$so)
  DateStamp("Adding oxycline...\n")
  ncvar_put(ncNew, O2, ENV$o2) 
  DateStamp("Adding PH...\n")
  ncvar_put(ncNew, PH, ENV$ph)
  
  nc_close(ncNew)
  
  return(invisible())
}

mergePPData = function(input, output, varName="npp", 
                        pattern = "_zint_", cf=NULL,  varNames = c("lpp", "spp")) {
  
  # get files
  files = osmose.fishmip:::.getEnvFiles(input, pattern, varNames)
  
  size = osmose.fishmip:::.checkVarSize(files)
  dims = osmose.fishmip:::.readNcdfDimensions(files)
  
  ind = which(sapply(files, FUN=length)>0)[1]
  outputFile = basename(gsub(pattern=paste("_", varNames[ind], "_", sep=""), 
                             replacement=paste("_", varName, "_", sep=""), 
                             x=files[[ind]])) 
  outputFile = gsub(pattern = pattern, replacement = "_", x=outputFile)
  outputFile = file.path(output, outputFile)
  
  ENV = list()
  for(iFile in seq_along(files)) {
    envName = names(files)[iFile]
    ENV[[envName]] = osmose.fishmip:::.readEnvData(file=files[[iFile]], dim=size)
  }
  
  dimLon  = ncdim_def("lon", "degrees", vals=dims[[1]])
  dimLat  = ncdim_def("lat", "degrees", vals=dims[[2]])
  dimTime = ncdim_def("time", "month", vals=dims[[3]])
  
  NPP   = ncvar_def(name="npp", units="mgC/m2/day", dim=list(dimLon, dimLat, dimTime), 
                    missval=1e+20, longname="Net primary production", prec="float",
                    compression=9)
  
  if(!file.exists(output)) dir.create(output, recursive=TRUE)
  
  ncNew = nc_create(filename=outputFile, 
                    vars=list(NPP))
  
  cf = (12*1000)*(60*60*24) # mol/m2/s -> mcG/m2/day
  
  DateStamp("Adding NPP...\n")
  ncvar_put(ncNew, NPP, cf*(ENV[[1]] + ENV[[2]])) 
  
  nc_close(ncNew)
  
  return(invisible())
}

