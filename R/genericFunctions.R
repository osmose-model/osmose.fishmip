
# Clines ------------------------------------------------------------------

calculateCline = function(obj, depth, k=1) {
  
  out = apply(obj, -3, .getCline, depth=depth, k=k)
  dim(out) = dim(obj)[-3]
  return(out)
  
}

.getCline = function(x, depth, k=1) {
  
  # check depths, convert to negative if necessary
  sig = 1
  if(all(is.na(x))) return(NA)
  if(any(is.na(depth))) 
    stop("NAs are not allowed for depths.")
  if(all(depth>=0)) {
    depth = -depth
    sig   = -1
  }
  if(any(depth>0)) stop("Not consistent depths provided.")
  yind = sort(depth, index.return=TRUE)$ix
  y = depth[yind]
  x = x[yind]
  
  yf = .fitProfile(x=x, y=y, k=k)
  ry = range(y, na.rm=TRUE)
  
  yyf = seq(from=ry[1], to=ry[2], by=1)
  dxf = yf(yyf, deriv=1)
  
  clina = yyf[which.max(dxf)] # clina
  
  return(sig*clina)
}

.getOxycline = function(x, depth, profile=FALSE) {
  
  # check depths, convert to negative if necessary
  sig = 1
  if(all(is.na(x))) return(NA)
  if(any(is.na(depth))) 
    stop("NAs are not allowed for depths.")
  if(all(depth>=0)) {
    depth = -depth
    sig   = -1
  }
  if(any(depth>0)) stop("Not consistent depths provided.")
  yind = sort(depth, index.return=TRUE)$ix
  y = depth[yind]
  x = x[yind]
  
  yf = .fitProfile(x=x, y=y, k=2, xtr=c(+1, -1))
  ry = range(y, na.rm=TRUE)
  
  yyf = seq(from=ry[1], to=ry[2], by=1)
  dxf = yf(yyf, deriv=1)
  
  cline = yyf[which.max(dxf)] # clina
  
  return(sig*cline)
}

fitOxygenProfile = function(x, depth) {
  
  # check depths, convert to negative if necessary
  sig = 1
  if(all(is.na(x))) return(NA)
  if(any(is.na(depth))) 
    stop("NAs are not allowed for depths.")
  if(all(depth>=0)) {
    depth = -depth
    sig   = -1
  }
  if(any(depth>0)) stop("Not consistent depths provided.")
  yind = sort(depth, index.return=TRUE)$ix
  y = depth[yind]
  x = x[yind]
  
  yf = .fitProfile(x=x, y=y, k=2, xtr=c(+1, -1))
  ry = range(y, na.rm=TRUE)
  
  yyf = seq(from=ry[1], to=ry[2], by=1)
  dxf = yf(yyf, deriv=1)
  xxf = yf(yyf)
  
  cline = yyf[which.max(dxf)] # clina
  oxy   = xxf[which.max(dxf)] # oxygen at clina
  
  out = list(x=x, y=y, xfit=xxf, yfit=yyf, 
             cline=cline, oxy=oxy)
  class(out) = c("oxygenProfile", class(out))
  return(out)
}

plot.oxygenProfile = function(object, cex=1, ...) {
  plot(object$x, object$y, pch=19, 
       xlab="Oxygen concentration (mL/L)",
       ylab="Depth (m)", axes=FALSE)
  axis(1)
  axis(2, at=axTicks(2), labels=abs(axTicks(2)))
  lines(object$xfit, object$yfit, col="blue", lwd=2)
  abline(h=object$cline, lty=2, col="red")
  abline(v=object$oxy, lty=2, col="red")
  txt1 = sprintf("Oxycline depth = %0.0f m", abs(object$cline))
  txt2 = sprintf("Oxygen concentration at oxycline = %0.1f mL/L",
                 object$oxy)
  mtext(txt1, 1, adj=0.95, line=-3, cex=cex)
  mtext(txt2, 1, adj=0.95, line=-2, cex=cex)
  box()
  return(invisible())
}


.fitProfile = function(x, y, k=1, xtr=1) {
  
  if(k==1 & length(xtr)==1) xtr = rep(xtr, 2)
  if(length(xtr)!=2) stop("Incorrect xtr specification.")
  
  # remove NAs from data
  ind = !is.na(x) 
  n = sum(ind)
  y = y[ind]
  x = x[ind]
  
  # checking monotonicity
  mp = makeMonotonic(x, y, k=k, xtr=xtr)
  x = mp$x
  y = mp$y
  
  if(mp$n == 1) {
    yf = splinefun(y, x, method="monoH.FC")
  } else {
    x1 = NULL
    y1 = NULL
    for(i in seq_len(mp$n)) {
      x0 = x[mp$index[[i]]]
      y0 = y[mp$index[[i]]]
      yf = splinefun(y0, x0, method="monoH.FC")
      ys = 0.5*(head(y0, -1) + tail(y0, -1))
      xs = 0.5*(head(x0, -1) + tail(x0, -1))
      xs = 0.5*(yf(ys) + xs)
      x0 = c(x0, xs)
      y0 = c(y0, ys)
      ind = sort(y0, index.return=TRUE)$ix
      yf = splinefun(y0[ind], x0[ind], method="monoH.FC")
      ry = range(y0, na.rm=TRUE)
      yyf = seq(from=ry[1], to=ry[2], by=10)
      x1 = c(x1, yf(yyf))
      y1 = c(y1, yyf)
    }
    yf = splinefun(y1, x1)
    # yf = smooth.spline(y1, x1, df = 0.6*length(y1))
  }
  
  return(yf)
  
}

.removeOne = function(x, y, xtr) {
  
  xs = rle(sign(diff(x)))
  hh = head(xs$lengths, 1)
  hv = head(xs$values, 1)
  tt = tail(xs$lengths, 1)
  tv = tail(xs$values, 1)
  n  = length(xs$lengths)
  
  if(hh==1 & hv==xtr[1]) xs$lengths[1] = 2
  if(tt==1 & tv==xtr[2]) xs$lengths[n] = 2
  
  l = rle(xs$length == 1)
  nrow   = cumsum(xs$lengths) 
  ind = cumsum(l$lengths)[which(l$value)]
  # index = cbind(begin=c(1, head(nrow, -1)+1), end=nrow+1)
  if(length(ind)==0) 
    return(list(x=x, y=y))
  out = nrow[ind]
  x0 = 0.5*(x[out] + x[out+1])
  y0 = 0.5*(y[out] + y[out+1])
  x1 = x[unique(-c(out, out+1))]
  y1 = y[unique(-c(out, out+1))]
  
  x = c(x0, x1)
  y = c(y0, y1)
  ind = sort(y, index.return=TRUE)$ix
  x = x[ind]
  y = y[ind]
  # rind = setdiff(ind, nrow(index))
  # index[rind+1, 1] = index[rind+1, 1] + 1
  # index = index[-ind,]
  return(list(x=x, y=y))
  
}

.removeOne0 = function(x, y, xtr) {
  
  xs = rle(sign(diff(x)))
  hh = head(xs$lengths, 1)
  hv = head(xs$values, 1)
  tt = tail(xs$lengths, 1)
  tv = tail(xs$values, 1)
  n  = length(xs$lengths)
  
  if(hh==1 & hv==xtr[1]) xs$lengths[1] = 2
  if(tt==1 & tv==xtr[2]) xs$lengths[n] = 2
  
  l = rle(xs$length == 1)
  nrow   = cumsum(xs$lengths) 
  ind = cumsum(l$lengths)[which(l$value)]
  # index = cbind(begin=c(1, head(nrow, -1)+1), end=nrow+1)
  if(length(ind)==0) 
    return(list(x=x, y=y))
  out = nrow[ind] + 1
  x1 = x[-out]
  y1 = y[-out]
  # rind = setdiff(ind, nrow(index))
  # index[rind+1, 1] = index[rind+1, 1] + 1
  # index = index[-ind,]
  return(list(x=x1, y=y1))
  
}

makeMonotonic = function(x, y, k=1, iter.max=15, xtr=1) {
  
  if(k==1 & length(xtr)==1) xtr = rep(xtr, 2)
  if(length(xtr)!=2) stop("Incorrect xtr specification.")
  
  k0 = length(rle(sign(diff(x)))$values)
  .FUN = function(index) {
    index = unlist(index)
    return(index[1]:index[2])
  }
  if(k0<k) warning("Less monotonic pieces than requested.")
  if(k0<=k) {
    xs = rle(sign(diff(x)))
    nrow   = cumsum(xs$lengths) 
    index = cbind(begin=c(1, head(nrow, -1)+1), end=nrow+1)
    index = lapply(apply(index, 1, list), .FUN)
    return(list(x=x, y=y, n=k0, index=index))
  }
  
  iter = 0
  while(k0>k & iter<iter.max) {
    p = .removeOne(x=x, y=y, xtr=xtr)
    x = p$x
    y = p$y
    k0 = length(rle(sign(diff(x)))$values)
    iter = iter + 1
  }
  if(k0>k) warning("Maximum of iterations reached.")
  
  xs = rle(sign(diff(x)))
  nrow   = cumsum(xs$lengths) 
  index = cbind(begin=c(1, head(nrow, -1)+1), end=nrow+1)
  index = lapply(apply(index, 1, list), .FUN)
  
  return(list(x=x, y=y, n=k0, index=index))
}


# Isolines ----------------------------------------------------------------


calculateIsoline = function(obj, depth, ref, k=1) {
  
  out = apply(obj, -3, .getIsoline, depth=depth, ref=ref, k=k)
  dim(out) = dim(obj)[-3]
  return(out)
  
}

.getIsoline = function(x, depth, ref, k=1) {
  
  rx = range(x, na.rm=TRUE)
  test = ref >= rx[1] & ref <= rx[2]
  if(!test) return(NA)
  
  sig = 1
  if(all(is.na(x))) return(NA)
  if(any(is.na(depth))) 
    stop("NAs are not allowed for depths.")
  if(all(depth>=0)) {
    depth = -depth
    sig   = -1
  }
  if(any(depth>0)) stop("Not consistent depths provided.")
  yind = sort(depth, index.return=TRUE)$ix
  y = depth[yind]
  x = x[yind]
  
  yf = .fitProfile(x=x, y=y, k=k)
  ry = range(y, na.rm=TRUE)
  
  yyf = seq(from=ry[1], to=ry[2], by=1)
  xxf = yf(yyf)
  fit = (xxf - ref)^2 
  
  iso = yyf[which.min(fit)] # isoline
  
  return(sig*iso)
  
}


# Var at bottom -----------------------------------------------------------


.getLast = function(obj) {
  .last = function(x) {
    if(all(is.na(x))) return(NA)
    out = tail(na.omit(x), 1)
    return(out)
  }
  out = apply(obj, -3, .last)
  dim(out) = dim(obj)[-3]
  return(out)
}

