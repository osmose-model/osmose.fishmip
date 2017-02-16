
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
  
  # remove NAs from data
  ind = !is.na(x) 
  n = sum(ind)
  y = y[ind]
  x = x[ind]
  
  # checking monotonicity
  mp = makeMonotonic(x, y, k=k)
  x = mp$x
  y = mp$y
  
  if(mp$n == 1) {
    yf = splinefun(y, x, method="monoH.FC")
    ry = range(y, na.rm=TRUE)
  } else {
    x1 = NULL
    y1 = NULL
    for(i in seq_len(mp$n)) {
      x0 = x[mp$index[[i]]]
      y0 = y[mp$index[[i]]]
      yf = splinefun(y0, x0, method="monoH.FC")
      ry = range(y0, na.rm=TRUE)
      yyf = seq(from=ry[1], to=ry[2], by=10)
      x1 = c(x1, yf(yyf))
      y1 = c(y1, yyf)
    }
    yf = splinefun(y1, x1)
    ry = range(y1, na.rm=TRUE)
  }
  
  yyf = seq(from=ry[1], to=ry[2], by=1)
  dxf = yf(yyf, deriv=1)
  
  clina = yyf[which.max(dxf)] # clina
  
  return(sig*clina)
}

.removeOne = function(x, y) {
  
  xs = rle(sign(diff(x)))
  
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

makeMonotonic = function(x, y, k=1, iter.max=15) {
  
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
    p = .removeOne(x, y)
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


calculateIsoline = function(obj, depth, ref, thr=0.15) {
  
  out = apply(obj, -3, .getIsoline, depth=depth, ref=ref, thr=thr)
  dim(out) = dim(obj)[-3]
  return(out)
  
}




.getIsoline = function(x, depth, ref) {
  
# in range, non-monotonicity

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

