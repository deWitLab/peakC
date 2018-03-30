
################################################
#4C specific functions
################################################


significant.fragments <- function( p.value, pos, window = 21, FDR = 0.01 ){
  #correct the nominal p-value for multiple hypothesis testing
  p.combined <- p.adjust(p.value, method="fdr")
  #determine the significant windows and select the fragments therein
  sig.i <- which(p.combined < FDR)
  if(length(sig.i)>0) {
		sig.i.start <- sig.i-floor(window/2); sig.i.end <- sig.i+floor(window/2)
		sig.i <- unique(multi.seq(sig.i.start,sig.i.end))
		sig.i <- sig.i[sig.i >= 1 & sig.i <= length(pos)]
		sigFrags <- pos[sig.i]
		return(sigFrags)
  } else {
    return(NULL)
  }
}

righttailgamma = function(r,k,n) 1 - pgamma(-log(r/(n+1)^k),k,scale=1)

rank.product.p <- function( data, num.exp,method="diff"){
  if(method=="diff") {
  stats <- data[,2:(num.exp+1)]-data[,(2:(num.exp+1))+num.exp]
  } else {
  stats <- data[,2:(num.exp+1)]/data[,(2:(num.exp+1))+num.exp]
  }
  rp <- nrow(data)-apply(stats,2,rank)+1
  rp <- apply(rp,1,prod)
  p <- righttailgamma(rp,num.exp,length(rp))
}


getWindowedFrags <- function(x,frags,wSize=21) {

	outFrags <- ((match(x,frags)-floor(wSize/2))):((match(x,frags)+floor(wSize/2)))
	outFrags <- outFrags[outFrags>=1&outFrags<=length(frags)]

	return(frags[outFrags])

}

#set a dynamic threshold for the residuals
getThreshold <- function(resids,qW=5) {
	q75 <- quantile(resids,probs=0.75) #75% quantile of the residuals
	qd50 <- diff(quantile(resids,probs=c(0.25,0.75))) #the range between the 25% and 75% quantiles
	threshold <- q75 + qW*qd50
	return(threshold)
}


#these are not "significant" frags just above an arbitrary threshold
thresholdFrags <- function(resids,frags,wSize=21,qW=5) {

	qMax <- getThreshold(resids=resids,qW=qW)

	sel.i <- which(resids > qMax)
  if(length(sel.i)>0) {
		sel.i.start <- sel.i-floor(wSize/2); sel.i.end <- sel.i+floor(wSize/2)
		sel.i <- unique(multi.seq(sel.i.start,sel.i.end))
		sel.i <- sel.i[sel.i >= 1 & sel.i <= length(frags)]
		selFrags <- frags[sel.i]
		return(selFrags)
	}else{
		return(NULL)
	}


}

#' Single experiment 4C/Capture C analysis
#'
#' @param data list containing the 4C/CapC data in two column format, an additional element num.exp describes the number of experiments
#' @param vp.pos viewpoint position, this can be a single value or a two values to analyse a viewpoint region
#' @param wSize number of fragments in a window
#' @param alphaFDR false-discovery rate threshold
#' @param qWd threshold for difference from the background
#' @param qWr threshold for ratio over the background
#' @param minDist minimal region around the viewpoint to exclude for the significance analysis
#'
#' @description Function for identifying interaction peaks above a background distribution. A list containing a 4C/Capture-C dataset is required as input. The viewpoint position is given in the vp.pos argument.
#'
#' @return a list containing a matrix with the data and the background model and a vector with the significant fragments
#' @export
#'
#' @examples
#' data <- readMultiple(f[1:3], vp.pos = 65923803)
#' res <- combined.analysis(data, num.exp=3, vp.pos = 65923803)
#'
#'
#'
single.analysis <- function(data, vp.pos, wSize = 21, qWd = 1.5, qWr = 1, minDist = 15e3) {

  #create two element vector containing the viewpoint position
  #if only one viewpoint is given
  if(length(vp.pos) == 1){
    vp.pos <- c(vp.pos,vp.pos)
  }
  vp.pos <- sort(vp.pos)


  db <- get.single.background(data=data, num.exp = 1, vp.pos=vp.pos)
  #running mean over the data
  db[,2] <- caTools::runmean(x=db[,2],k=wSize,endrule="mean")
  #running mean over the isotonic regression line
  db[,3] <- caTools::runmean(x=db[,3],k=5,endrule="mean")

  #add a pseudocount to improve the calculations
  pseudoCount <- non.zero.quantile(x=db[,2],probs=0.05)
  ratios <- cbind(db[, 1], (db[, 2] + pseudoCount)/(db[,3] + pseudoCount))
  deltas <- cbind(db[, 1], db[,2]-db[,3])

  #frags <- db[, 1]
  #distFrags <- frags[abs(frags-vp.pos)>=minDist]

  #select the fragments that are more than minDist from the viewpoint
  sel.frag <- db[which( (db[,1] < vp.pos[1] & vp.pos[1]-db[,1] > minDist) | (db[,1] > vp.pos[2] & db[,1]-vp.pos[2] > minDist) ),1]

  rTFrags <- thresholdFrags(resids=ratios[,2],frags=ratios[,1],wSize=wSize,qW=qWr)
  dTFrags <- thresholdFrags(resids=deltas[,2],frags=deltas[,1],wSize=wSize,qW=qWd)
  peakFrags <- intersect(intersect(rTFrags,dTFrags), sel.frag)

  return(list(dbR=db,ratios=ratios,deltas=deltas,peak=peakFrags, num.exp=1))

}

#' Combined 4C/Capture C analysis
#'
#' @param data list containing the 4C/CapC data in two column format, an additional element num.exp describes the number of experiments
#' @param num.exp number of experiments, used in conjuction with "data" and "multi" in type, default is 0 which means the number in the data list is used, a different number overwrites the default number
#' @param vp.pos viewpoint position, this can be a single value or a two values to analyse a viewpoint region
#' @param wSize number of fragments in a window
#' @param alphaFDR false-discovery rate threshold
#' @param qW threshold for absolute difference
#' @param minDist minimal region around the viewpoint to exclude for the significance analysis
#'
#' @description Function for identifying interaction peaks above a background distribution. A list of 4C/Capture-C datasets are required as input. The viewpoint position is given in the vp.pos argument.
#'
#' @return a list containing a matrix with the data and the background model and a vector with the significant fragments
#' @export
#'
#' @examples
#' data <- readMultiple(f[1:3], vp.pos = 65923803)
#' res <- combined.analysis(data, num.exp=3, vp.pos = 65923803)
#'
#'
#'
combined.analysis <- function( data, num.exp = 0, vp.pos, wSize = 21, alphaFDR = 0.1, qWr = 1, minDist = 15e3 ){
  #set the number of experiments
  if(num.exp == 0){
    num.exp = data$num.exp

  }

  #create two element vector containing the viewpoint position
  #if only one viewpoint is given
  if(length(vp.pos) == 1){
    vp.pos <- c(vp.pos,vp.pos)
  }
  vp.pos <- sort(vp.pos)


	db <- combine.experiments(data,num.exp, vp.pos)

	# make a data.frame where a running mean is already applied to apply all statistics to those data -> stronger (positive) dependency between statistics but less variance

	dbR <- db
	#running mean over the data
	dbR[,2:(num.exp+1)] <- apply(db[,2:(num.exp+1)],2,caTools::runmean,k=wSize,endrule="mean")
	#running mean over the isotonic regression line (window is 5)
	dbR[,2:(num.exp+1)+num.exp] <- apply(db[,2:(num.exp+1)+num.exp],2,caTools::runmean,k=5,endrule="mean")

	# for the ratios add a small pseudocount to avoid dividing by 0. Calculate the ratios and deltas (diff) on the runmean data

	#add a small pseudo count so that there will be no divide/0 errors
	pseudoCount <- apply(db[,2:(num.exp+1)], 2, non.zero.quantile, probs=0.05)
	pseudoCount <- sum(pseudoCount)/num.exp
	#calculate the ratio of the data with the regression line
	ratio <- cbind(db[,1],(dbR[,2:(num.exp+1)]+pseudoCount)/(dbR[,(2:(num.exp+1))+num.exp]+pseudoCount))
	#calculate the differene between the data and the regression line
	delta <- cbind(db[,1],dbR[,2:(num.exp+1)]-dbR[,(2:(num.exp+1))+num.exp])

	#determine the per-window p-value using rank products based on the ratio
	p.val <- rank.product.p(data = dbR, num.exp = num.exp,method="diff")
	#select the significant fragments
	sfr <- significant.fragments(p.value = p.val, pos = db[, 1], window = wSize, FDR = alphaFDR)

	#select the fragments that are more than minDist from the viewpoint
	sel.frag <- db[which( (db[,1] < vp.pos[1] & vp.pos[1]-db[,1] > minDist) | (db[,1] > vp.pos[2] & db[,1]-vp.pos[2] > minDist) ),1]
	idx <- delta[,1]%in%sel.frag

	#set a threshold on the minimal delta threshold, this threshold is defined empirically
	tfr <- thresholdFrags(resids=apply(ratio[idx,2:(num.exp+1)],1,mean),frags=ratio[idx,1],wSize=wSize,qW=qWr)

	sfr <- intersect(sfr,tfr)
	list(dbR=dbR, peak=sfr, num.exp = num.exp, p.value=p.val, ratio = apply(ratio[,2:(num.exp+1)],1,mean), sel=sel.frag )
}

#' Take a result from the combined.analysis function and generate a chromosomal map of the result
#'
#' @param data list containing the output of combined.analysis (i.e. 4C data and significant fragments)
#' @param num.exp number of experiments
#' @param y.min bottom limit of the plot
#' @param y.max top limit of the plot
#'
#' @return Nothing, a plot is drawn
#' @export
#'
#' @examples
plot_C <- function(data, num.exp = 0, y.min=0, y.max=3000, ...){
  if(num.exp == 0){
    num.exp = data$num.exp
  }
  pos <- data$dbR[,1]
  if(num.exp == 1){
    y.ave <- data$dbR[,2]
  }else{
    y.ave <- apply(data$dbR[,2:(num.exp+1)], 1, median)
  }
  plot(pos, y.ave, type='h', col=ifelse(pos%in%data$peak, "red", "grey"), axes=F, xlab="chromosomal position", ylab="4C signal", ylim=c(y.min,y.max), ... )
  axis(2, at=c(0,y.max), las=2)
  at <- seq(200e3*floor(min(pos)/200e3), ceiling(max(pos)/200e3)*200e3, by=200e3)
  axis(1, at=seq(0,1e9,by=200e3), lab=sprintf("%.1f", seq(0,1e9,by=200e3)/1e6), cex.axis=1.5)
}


#' Combine list of experiments into a matrix with the background model
#'
#' @param data list containing the 4C/CapC data in two column format
#' @param num.exp number of experiments, default is 0, which means the the number in the data list is taken, other values allow for choosing a subset of the experiments
#' @param vp.pos position of the viewpoint
#'
#' @return merged matrix containing: 1. the position of the fragments, 2:(n+1) the data, (n+1):(n+n+1) background models corresponding to the respective datasets
#' @export
#'
#' @examples
combine.experiments <- function( data, num.exp = 0, vp.pos ){
  if(num.exp == 0){
    num.exp = data$num.exp
  }
  #create two element vector containing the viewpoint position
  #if only one viewpoint is given
  if(length(vp.pos) == 1){
    vp.pos <- c(vp.pos,vp.pos)
  }
  vp.pos <- sort(vp.pos)
  data.m <- data[[1]]
  for( i in 2:num.exp ){
    data.m <- merge(data.m, data[[i]], by=1)
  }
  #create the background model for the upstream regions
  data.bg <- data.m
  for( i in 1:num.exp ){
    data.bg[data.m[,1] < vp.pos[1],i+1] <- get.background(data.m[data.m[,1] < vp.pos[1],c(1,i+1)], vp.pos[1] )
  }
  #and for the downstream regions
  for( i in 1:num.exp ){
    data.bg[data.m[,1] > vp.pos[2],i+1] <- get.background(data.m[data.m[,1] > vp.pos[2],c(1,i+1)], vp.pos[2] )
  }
  #if two viewpoint fragments are given set the intervening fragments
  #to zero
  #set background to 1 to prevent NaN in the ratio
  if(vp.pos[1] != vp.pos[2]){
    for( i in 1:num.exp){
      data.m[data.m[,1] >= vp.pos[1] & data.m[,1] <= vp.pos[2],i+1] <- 0
      data.bg[data.m[,1] >= vp.pos[1] & data.m[,1] <= vp.pos[2],i] <- 1
    }
  }
  cbind(data.m, data.bg[,-1])
}

#perform pava regression and return the background regression line
get.background <- function( data, vp.pos, weight.factor=0, fractile=F){
  require(isotone)
  switched = FALSE
  weights <- (1:nrow(data))**weight.factor
  if(data[1,1] > vp.pos){
    data[,1] <- -data[,1] #reverse the sign to make the trend increasing
    switched = TRUE
    weights <- rev(weights)
  }
  #create the isotonic regression
  if(fractile){
    lm <- gpava(data[,1], data[,2], solver=weighted.fractile, weights=NULL, p=0.75)
  }else{
    lm <- gpava(data[,1], data[,2], solver=weighted.mean)
  }

  if(switched)
    pred.data <- data.frame( -lm$z, lm$x )
  else
    pred.data <- data.frame( lm$z, lm$x )

  pred.data[order(pred.data[,1]),2]
}

get.single.background <- function(data, num.exp = 1, vp.pos) {

  if (length(vp.pos) == 1) {
    vp.pos <- c(vp.pos, vp.pos)
  }
  vp.pos <- sort(vp.pos)

  data.bg <- data
  data.bg[data[, 1] < vp.pos[1], 2] <- get.background(data[data[, 1] < vp.pos[1], c(1, 2)], vp.pos[1])
  data.bg[data[, 1] > vp.pos[2], 2] <- get.background(data[data[, 1] > vp.pos[2], c(1, 2)], vp.pos[2])

  return(cbind(data, data.bg[, -1]))

}

#' Single experiment Jaccard similarity calculation
#'
#' @param data list containing the 4C/CapC data in two column format, an additional element num.exp describes the number of experiments
#' @param vp.pos viewpoint position, this can be a single value or a two values to analyse a viewpoint region
#' @param wSize number of fragments in a window
#' @param alphaFDR false-discovery rate threshold
#' @param qWd threshold for difference from the background
#' @param qWr threshold for ratio over the background
#' @param minDist minimal region around the viewpoint to exclude for the significance analysis
#' @description Function for identifying interaction peaks above a background distribution. A list of 4C/Capture-C datasets are required as input. The viewpoint position is given in the vp.pos argument.
#'
#' @return a matrix containing the pairwise Jaccard similarity scores
#' @export
#'
#' @examples
#' data <- readMultiple(f[1:3], vp.pos = 65923803)
#' res <- combined.analysis(data, num.exp=3, vp.pos = 65923803)
#'
#'
#'
pairwise.jaccard <- function(data, vp.pos, wSize = 21, qWd = 1.5, qWr = 1, minDist = 15e3) {
  num.exp <- data$num.exp
  js.mat <- matrix(0, nrow=num.exp, ncol=num.exp)
  diag(js.mat) <- 1
  res.list <- list()
  #first calculate the single experiment peak calling
  for( i in 1:num.exp){
    res.list[[i]] <- single.analysis(data[[i]], vp.pos = vp.pos, wSize = wSize, qWd = qWd, qWr = qWr, minDist = minDist)
  }
  for(i in 1:(num.exp-1)){
    for(j in (i+1):num.exp){
      js <- jaccardSim(res.list[[i]]$peak, res.list[[j]]$peak)
      js.mat[i,j] <- js; js.mat[j,i] <- js
    }
  }
  js.mat
}

#calculate the jaccard index for two vectors
jaccardSim <- function( p1, p2 ){
  intersection <- sum(p1 %in% p2)
  union <- length(unique(c(p1,p2)))
  intersection/union
}

##################################################################
#General functions
##################################################################

#running mean function
#' Title
#'
#' @param x numeric vector
#' @param n window size for the running window
#'
#' @return a numeric vector with the windowed means
#' @export
#'
#' @examples
running<-function(x,n=20){
  cumsum(x)->sum.v
  sum.v<-c(0,sum.v)
  #(sum.v[(n+1):length(x)]-sum.v[1:(length(x)-n)])/n
  diff(sum.v,n)/n
}

#running sum
runsum<-function(x,n=20){
  cumsum(x)->sum.v
  sum.v<-c(0,sum.v)
  diff(sum.v,n)
}

#make running function compatible vector
#' Remove leading and trailing values
#'
#' @param a vector
#' @param n window size of the corresponding running function
#' @description remove n/2 elements from the front and the end
#' @return vector, shortened by (n - 1) elements
#' @export
#'
#' @examples
rem <- function(a, n ){
  half.window <- floor(n/2)
  head(tail(a, -half.window),-half.window)
}

#quick way of generating a vector with the required indexes
multi.seq <- function( start, end ){
  x <- rep(start, end-start+1)->x
  df <- diff(x)
  df <- df + 1
  low <- which(df > 1)
  df[low] <- -diff(c(0,low))+1
  add <- c(0,cumsum(df))
  x + add
}

#wrapper function to calculate the quantile distribution of all
#non-zero values
non.zero.quantile <- function( x, probs ){
  quantile(x[x > 0], probs)
}


