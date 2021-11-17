#########################################
# Adapted from - Splinter et al. 2012 Methods
#Supplemental R script
#
#Within code is contained that identifies significant
#contacts in a 4C experiment. The function analyzeAndWriteBatchWig
#requires a list (vector) of wiggle filenames (which should also end
#with wig). The significant contacts on the cis and trans chromosomes are
#calculated seperatedly as detailed in the main text.


#written by Elzo de Wit (2012) ; edited by Felipe S D

#########################################

#quickly calculate the running sum
running.sum <- function( x, n ){
  sum.v <- c(0,cumsum(x))
  diff(sum.v,n)
}	

#calculate the Z-score as described by Splinter et al. (2011, 2012)
#sigma is estimated from the binomial distribution in a window size (large.window)
#that greatly exceed the regular window size (window)
z.score <- function( x, window=20, large.window=3000){
  p.large <- running.sum(x > 0, large.window)/large.window
  first <- p.large[1]
  last  <- tail(p.large,1)
  p.large <- c(rep(first,large.window/2), p.large, rep(last,large.window/2-1))
  p.large <- tail(p.large, n=-window/2)
  p.large <- head(p.large, n=-window/2+1)
  
  
  #calculate the Z score
  X <- running.sum(x > 0, window)
  mu    <- p.large*window
  sigma <- window*p.large*(1-p.large)
  
  Z <- (X - mu)/(sigma**0.5)
  
  Z
}

fdr.calc <- function(x, window=20, large.window=3000, iterations=1, fdr=0.01){
  z.real <- z.score(x, window=window, large.window=large.window)
  
  #now iteratively calculate the random Z-score
  z.random <- c()
  for(i in 1:iterations){
    z.random <- c(z.random, z.score(sample(x), window=window, large.window=large.window))
  }
  id.vec <- c(rep(1,length(z.real)), rep(-1,length(z.random)))
  z.vec <- c(z.real, z.random)
  
  id.vec <- id.vec[rev(order(z.vec))]
  z.vec <- rev(sort(z.vec))
  
  n.z.real   <- cumsum(ifelse(id.vec > 0, 1, 0))
  n.z.random <- cumsum(ifelse(id.vec < 0, 1/iterations, 0))
  
  fdr.i <- max(which(n.z.random/n.z.real < fdr))
  z.vec[fdr.i]
}

fast.fdr.trans <- function(gff, window=20, iterations=1000, fdr=0.01){
  x <- gff[,6] > 0
  true <-	table( factor(running.sum(x,window), levels=0:window) )
  for(i in 1:iterations){
    shuffled <- factor(running.sum(sample(x),window), levels=0:window)
    if(!exists("random")){
      random <- table(shuffled)
    }else{
      random <- random + table(shuffled)
    }
  }
  sum.data <- cbind(true,random)
  true <- cumsum(rev(sum.data[,1]))
  random <- cumsum(rev(sum.data[,2])/iterations)
  fractions <- random[true>0]/true[true>0]
  sig <- fractions[fractions < fdr]
  cut.off <- min(as.numeric(names(sig)))
  cut.off
}	

#qwiqly read a wig file (with only one chromosome)
readqWig <- function( file, skip = 2 ){
  wig <- scan(file, skip = skip)
  matrix(wig, ncol=2, byrow=T)
}

#function for reading in wiggle tracks into an R object
#(list) that can be easily processed
readWig <- function( file ){
  wig <- scan(file, what=list(""), sep="\t", quote="")
  chrom.i <- grep("variable", wig[[1]])
  chrom.i <- c(chrom.i, length(wig[[1]]))
  chrom.vec <- c()
  pos.data <- NULL
  for(i in 1:(length(chrom.i)-1)){
    start <- chrom.i[i]+1
    if(i == length(chrom.i)-1 )
      end   <- chrom.i[i+1]
    else
      end   <- chrom.i[i+1]-1
    chrom <- sub("variableStep  chrom=", "", wig[[1]][start-1])
    #skip if there are no rows in the dataset
    if(start==end+1)
      next
    wig.sub <- as.numeric(unlist(strsplit(wig[[1]][start:end]," ")))
    chrom.data <- matrix(wig.sub,ncol=2, byrow=T)
    pos.data   <- rbind(pos.data, chrom.data)
    chrom.vec <- c(chrom.vec, rep(chrom, nrow(chrom.data)))
  }
  wig.out <- data.frame(chr = chrom.vec, pos.data)
  wig.out
}


getVPpos <- function( file ){
  pos <- scan(file, nlines=1, what="string")[4]
  as.numeric(sub(".*:([0-9]+).*","\\1", pos))
}	

getChrom <- function( file ){
  pos <- scan(file, nlines=1, what="string")[4]
  sub(".*(chr.*):[0-9]+.*","\\1", pos)
}


#function for exactly deterimining the blocks
smallBlocks <- function(chrom, start, end, filename ){
  len <- length(start)   #should return 1
  diff <- which(start[-1]-end[-len] > 1e5)
  new.start <- start[c(1,diff+1)]
  new.end   <- end[c(diff,len)]
  cbind( new.start, new.end)
}	


trans.interactions <- function( gff, window=500, iterations=1000, filename="MCF7_L3MBTL3-B_basic4cseq.wig"){
  
  domains <- NULL
  chromosomes <- unique(gff[,1])
  for( chrom in chromosomes ){
    #skip chromosomes smaller than 1000 fragments
    if(nrow(gff[gff[,1]==chrom,]) < 1000)
      next
    print(paste("Doing", chrom))
    cut.off <- fast.fdr.trans(gff=gff[gff[,1]==chrom,],window=window, iterations=iterations)
    if(!is.finite(cut.off))
      next
    
    window.data <- running.sum(gff[gff[,1]==chrom,6] > 0, window)
    pos <- gff[gff[,1]==chrom,4]
    x <- which(window.data >= cut.off)
    start <- pos[x]
    end   <- pos[x+window]
    end[is.na(end)] <- max(pos)
    
    chrom.dom <- smallBlocks(chrom = 6, start = 130339274, end = 130339558, filename = "MCF7_L3MBTL3-B_basic4cseq.wig")
    chrom.dom <- data.frame(chrom = rep(chrom,nrow(chrom.dom)), d_start = chrom.dom[,1], d_end = chrom.dom[,2])
    domains <- rbind(domains, chrom.dom)
  }
  domains
}	

#from the indexes create a domain
getDomain <- function(index, gff, max.diff=100){
  start <- c(T, diff(index) > max.diff)
  end <-   c(diff(index) > max.diff, T)
  
  start.pos <- gff[index[start],4]
  end.pos   <- gff[index[end],4]
  data.frame(d_start=start.pos, d_end=end.pos)
}

#determine the indices that belong go over the FDR Z-score
get.indices.z <- function(gff, fdr, window=20, large.window=3000, dist.remove=2e6, vp.pos=vp.pos){
  rem.index <- which(gff[,4] > vp.pos - dist.remove & gff[,4] < vp.pos + dist.remove)
  index <- which(z.score(gff[,6],window=window, large.window=large.window)>fdr)
  full.index <- c(rep(index,each=window)+rep(0:(window-1),times=length(index)))
  full.index <- unique(full.index)
  full.index <- full.index[!full.index%in%rem.index]
  full.index
}	

calculateIndex <- function( gff, large.window=3000, window=100, fdr=0.01, iterations=100, vp.pos=vp.pos){
  #calculate the z-scores
  z <- z.score(gff[,6],large.window=large.window, window=window)
  fdr <- fdr.calc(gff[,6], large.window=large.window, window=window, fdr=fdr, iterations=iterations)
  index <- get.indices.z(gff, fdr=fdr, large.window=large.window, window=window, vp.pos=vp.pos)
  index
}



analyzeAndWrite <- function( filename, vp.chrom, vp.pos ){ 
  data <- readWig(filename)
  #create a gff structure
  gff <- data.frame(data[,1],"4C_seq", ".", data[,2], data[,2]+1, data[,3], ".")
  gff[,1] <- sub("chr","",gff[,1])
  vp.chrom <- sub("chr", "", vp.chrom)
  
  cis.gff <- gff[gff[,1]==vp.chrom,]
  trans.gff <- gff[gff[,1]!=vp.chrom,]
  index.4c <- calculateIndex(cis.gff, vp.pos=vp.pos)
  dom <- getDomain(index.4c, cis.gff)
  
  
  trans.data <- trans.interactions(trans.gff, iterations=100)
  print(dom)
  print(trans.data)
  all.dom <- rbind(data.frame(chrom = vp.chrom,dom), trans.data)
  cis <- ifelse(all.dom[,1] == vp.chrom, 1, 0)
  upload.df <- data.frame(all.dom, cis=cis)
  
  #manage the writing of the data to a file
  out.file <- filename
  out.file <- sub(".wig", "_targets.txt", out.file)
  write.table(upload.df, out.file, quote=F, col=F, row=F, sep="\t")
  
}


#example
trans_data <- analyzeAndWrite("MCF7_L3MBTL3-B_basic4cseq.wig", vp.chrom = "chr6", vp.pos = NA)
