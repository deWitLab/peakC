

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
  i <- 1
  for(f in files){
    d <- readqWig(f, vp.pos=vp.pos, window=window)
    data.list[[i]] <- d
    i <- i+1
  }
  data.list$num.exp = length(data.list)
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
  i <- 1
  for(f in files){
    d <- readMatrix(f, vp.pos=vp.pos, window=window, normalize=normalize)
    data.list[[i]] <- d
    i <- i+1
  }
  data.list$num.exp = length(data.list)
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
  #normalize to 1M reads
  for( i in 1:num.exp){
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
  #select only the non-blind fragments
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
  d
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
  d
}
