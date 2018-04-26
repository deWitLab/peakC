

#' Read in a list of wig files
#'
#' @param files vector containg paths to wig files
#' @param vp.pos position of the viewpoint
#' @param window flanking sequence to be considered in the analysis
#'
#' @return A list containing two column matrices with the position and the normalized coverage score
#' @export
#'
#' @examples
readMultipleWig <- function( files, vp.pos, window = 700e3 ){
  data.list <- list()
  quality <- list()
  i <- 1
  for(f in files){
    d <- readqWig(f, vp.pos=vp.pos, window=window)
    if(i == 1){
      quality <- d$quality
    }else{
      quality$percentage.capture.100kb[i] <- d$quality$percentage.capture.100kb
      quality$percentage.capture.1Mb[i] <- d$quality$percentage.capture.1Mb
      quality$percentage.capture.cis[i] <- d$quality$percentage.capture.cis
      quality$total.read.cis[i] <- d$quality$total.read.cis
    }
    data.list[[i]] <- d$data
    i <- i+1
  }
  data.list$num.exp = length(data.list)
  data.list$quality = quality
  data.list
}

#' Read in a list of matrix files
#'
#' @param files vector containg paths to two column matrix files
#' @param vp.pos position of the viewpoint
#' @param window flanking sequence to be considered in the analysis
#'
#' @return A list containing two column matrices with the position and the normalized coverage score
#' @export
#'
#' @examples
readMultiple <- function( files, vp.pos, window = 700e3, normalize=T ){
  data.list <- list()
  quality <- list()
  i <- 1
  for(f in files){
    d <- readMatrix(f, vp.pos=vp.pos, window=window, normalize=normalize)
    data.list[[i]] <- d$data
    if(i == 1){
        quality <- d$quality
    }else{
      quality$percentage.capture.100kb[i] <- d$quality$percentage.capture.100kb
      quality$percentage.capture.1Mb[i] <- d$quality$percentage.capture.1Mb
      quality$percentage.capture.cis[i] <- d$quality$percentage.capture.cis
      quality$total.read.cis[i] <- d$quality$total.read.cis
    }
    i <- i+1
  }
  data.list$num.exp = length(data.list)
  data.list$quality = quality
  data.list
}

#' Read in a multiple experiment file
#'
#' @param file path to the experiment file
#' @param vp.pos position of the viewpoint
#' @param window flanking sequence to be considered in the analysis
#' @param num.exp Number of experiments
#'
#' @return A list containing two column matrices with the position and the normalized coverage score
#' @export
#'
#' @examples
readMultiColumnFile <- function(file, vp.pos, window=700e3, num.exp=4){
  d <- read.delim(file, h=F, stringsAsFactors=F)
  d <- d[,1:(num.exp+1)] #in case there are weird empty columns

  quality <- list()

  #normalize to 1M reads
  for( i in 1:num.exp){
    #calculate and store the quality characteristics
    q.temp <- quality.metrics(d[,c(1,i+1)], vp.pos)
    quality$percentage.capture.100kb[i] <- q.temp$percentage.capture.100kb
    quality$percentage.capture.1Mb[i] <- q.temp$percentage.capture.1Mb
    quality$percentage.capture.cis[i] <- q.temp$percentage.capture.cis
    quality$total.read.cis[i] <- q.temp$total.read.cis
    num.reads <- sum(d[,i+1], na.rm=T)
    d[,i+1] <- 1e6*d[,i+1]/num.reads
  }
  d <- d[d[,1] > vp.pos-window & d[,1] < vp.pos+window,]
  #make a list out of the matrix to make it compatible with the peak caller
  data.list <- list()
  for( i in 2:(num.exp+1)){
    d.sub <- d[,c(1,i)]
    colnames(d.sub) <- c("frag_pos", "frag_score")
    data.list[[i-1]] <- d.sub
  }
  data.list$num.exp = length(data.list)
  data.list$quality = quality
  data.list
}

#' Transform a data frame to a list that can be used as input for peakC
#'
#' @param df data.frame containing the experimental data
#' @param vp.pos position of the viewpoint
#' @param window flanking sequence to be considered in the analysis
#' @param num.exp Number of experiments
#'
#' @return A list containing two column matrices with the position and the normalized coverage score
#' @export
#'
#' @examples
data.frame.to.peakC <- function( df, vp.pos, window, num.exp ){
  df <- df[,1:(num.exp+1)] #in case there are weird empty columns

  quality <- list()

  #normalize to 1M reads
  for( i in 1:num.exp){
    #calculate and store the quality characteristics
    q.temp <- quality.metrics(df[,c(1,i+1)], vp.pos)
    quality$percentage.capture.100kb[i] <- q.temp$percentage.capture.100kb
    quality$percentage.capture.1Mb[i] <- q.temp$percentage.capture.1Mb
    quality$percentage.capture.cis[i] <- q.temp$percentage.capture.cis
    quality$total.read.cis[i] <- q.temp$total.read.cis
    num.reads <- sum(df[,i+1], na.rm=T)
    df[,i+1] <- 1e6*df[,i+1]/num.reads
  }
  df <- df[df[,1] > vp.pos-window & df[,1] < vp.pos+window,]
  #make a list out of the matrix to make it compatible with the peak caller
  data.list <- list()
  for( i in 2:(num.exp+1)){
    d.sub <- df[,c(1,i)]
    colnames(d.sub) <- c("frag_pos", "frag_score")
    data.list[[i-1]] <- d.sub
  }
  data.list$num.exp = length(data.list)
  data.list$quality = quality
  data.list
}



#qwigly read a wig file (with only one chromosome)
#' Quickly read and normalize a wig formatted file
#'
#' @param file path to a wiggle file
#' @param window genomic window around the viewpoint to read in
#' @param vp.pos postion of the viewpoint in the genome
#'
#' @return a matrix with two columns, position and score
#' @export
#' @description Wrapper function for reading files that are formatted as wig files. The data is also normalized to 1 million sequencing reads.
#' @examples
#' data <- readqWig(file="alpha.wig", window = 700e3, vp.pos = 32224333 )
readqWig <- function( file, window, vp.pos ){
  wig <- scan(file, skip = 2, quiet = T)
  d <- matrix(wig, ncol=2, byrow=T)
  d <- d[-which.max(d[,2]),]
  d <- d[d[,1]!=vp.pos,]

  #calculate the quality metrics for the experiment
  quality <- quality.metrics(d, vp.pos)

  #check if the file contains any data
  if(sum(d[,2]) > 0){
    d[,2] <- 1e6*d[,2]/sum(d[,2])
  }else{
    stop("Data file does not contain any data")
  }

  if(window > 0){
    d <- d[d[,1] > vp.pos - window & d[,1] < vp.pos + window,]
  }else{
    #select a genomic region around the viewpoint with a given amount of coverage
    i <- range(which( running(d[,2]>0,2001) > 0.2))+1000
    print(i)
    if(any(is.infinite(i))){
      window = 100e3
      d <- d[d[,1] > vp.pos - window & d[,1] < vp.pos + window,]
    }else{
      d <- d[i[1]:i[2],]
    }

  }
  colnames(d) <- c("frag_pos", "frag_score")
  #d
  list(data=d, quality=quality)
}

#qwigly read a 4C file (with only one chromosome)
#' Quickly read and normalize a matrix formatted file
#'
#' @param file path to a two column file
#' @param window genomic window around the viewpoint to read in
#' @param vp.pos postion of the viewpoint in the genome
#'
#' @return a matrix with two columns, position and score
#' @export
#' @description Wrapper function for reading files that are formatted as wig files. The data is also normalized to 1 million sequencing reads.
#' @examples
#'
readMatrix <- function( file, window, vp.pos, normalize = T ){
  vec <- scan(file, quiet = T)
  d <- matrix(vec, ncol=2, byrow=T)
  d <- d[-which.max(d[,2]),]
  d <- d[d[,1]!=vp.pos,]

  #calculate the quality metrics for the experiment
  quality <- quality.metrics(d, vp.pos)

  if(normalize){
    if(sum(d[,2]) > 0){
      d[,2] <- 1e6*d[,2]/sum(d[,2])
    }else{
      stop("Data file does not contain any data")
    }
  }
  if(window > 0){
    d <- d[d[,1] > vp.pos - window & d[,1] < vp.pos + window,]
  }else{
    #select a genomic region around the viewpoint with a given amount of coverage
    i <- range(which( running(d[,2]>0,2001) > 0.2))+1000
    print(i)
    if(any(is.infinite(i))){
      window = 100e3
      d <- d[d[,1] > vp.pos - window & d[,1] < vp.pos + window,]
    }else{
      d <- d[i[1]:i[2],]
    }

  }
  #add names to the columns
  colnames(d) <- c("frag_pos", "frag_score")
  #d
  list(data=d, quality=quality)
}

###############################################
#Quality metrics function
###############################################

#internal function for calculating the quality metrics for a specific
#experiment
#the quality metrics that are assessed are:
#% captured fragments within the first 100kb
#% captured fragments within the first 1Mb
#% captured fragments total chromosome
#total number of reads in cis
quality.metrics <- function(data, vp.pos){
  quality <- list()
  quality$percentage.capture.100kb <- 100*mean( data[data[,1] > vp.pos-100e3 & data[,1] < vp.pos+100e3,2] > 0)
  quality$percentage.capture.1Mb   <- 100*mean( data[data[,1] > vp.pos-1e6   & data[,1] < vp.pos+1e6,2] > 0)
  quality$percentage.capture.cis   <- 100*mean(data[,2] > 0)
  quality$total.read.cis <- sum(data[,2])
  quality
}
