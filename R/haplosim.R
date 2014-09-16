#################################################################################
### HAPLOSIM - METHOD FOR FINDING SIMILAR HAPLOTYPES (version dated 07/12/09) ###
#################################################################################
# main input data is a matrix of 0s and 1s, where each   
# row codes for a chromosome, and each column codes for a 
# unique biallelic SNP, with missing values coded as 9 by default.
# position = vector containing either genetic or physical distance
# miss.code = code for missing value, default as 9
# snp.miss = vector of the same length as the number of SNPs with elements denoting the amount of missingness at each SNP.
# focal = TRUE/FALSE statement about whether there is a focal position; if FALSE, focal.flag and focal.weight will be ignored even if specified. 
# focal.flag = index of SNP in file to be the focal SNP
# focal.weight = enforced weighting of focal SNP 
# pos.flag = which column contains information about position, only necessary if focal = T
# distance.measure = either "genetic" or "physical" to represent whether genetic (cM) or physical distance (bp) is used
# n.hap = the number of canonical haplotypes to find, default to "auto" for automatic detection
# tolerance = scale representing degree of similarity required for clustering haplotypes, with 1 at most stringent, with a recommended maximum of 5 for fuzzy clustering.
# sticky = an odd integer, representing the number of SNPs to consider jointly when considering mosaic of unmapped haplotypes. 

haplosim <- function(input, position, miss.code = 9, snp.miss = NULL, focal = TRUE, focal.flag = NA, focal.weight = -1, distance.measure = "genetic", n.hap = "auto", tolerance = 1, sticky = 5){
   require(lattice)   
   n.chr <- dim(input)[1]
   n.snp <- dim(input)[2]
   miss.matrix <- matrix(0, nr = n.chr, nc = n.snp)
   if (n.chr >= n.snp){
      for (i in 1:n.snp){
         na.flag <- which(input[,i] == miss.code)
         if (length(na.flag) > 0) {
            input[na.flag,i] <- NA
            miss.matrix[na.flag,i] <- 1
         }
      }
   }

   if (n.chr < n.snp){
      for (i in 1:n.chr){
         na.flag <- which(input[i,] == miss.code)
         if (length(na.flag) > 0) {
            input[i,na.flag] <- NA
            miss.matrix[i,na.flag] <- 1
         }
      }
   }

   print(paste("Analysing ", n.chr, " chromosomes across ", n.snp, " SNPs.", sep=""))
   if (focal == T){
      ifelse(is.na(focal.flag) == T, focal.pos <- mean(range(position)), focal.pos <- position[focal.flag])
      if (is.na(focal.flag)) print("Focal SNP not specified, center of genomic region used as focal point.")
      if (is.na(focal.flag) == F & focal.weight == -1) print(paste("SNP with index ", focal.flag, " chosen as focal SNP at position ", focal.pos, ".", sep=""))
      if (is.na(focal.flag) == F & focal.weight != -1) print(paste("SNP with index ", focal.flag, " chosen as focal SNP at position ", focal.pos, " with a weighting of ", focal.weight, ".", sep=""))
   }
   if (focal == F){
      print("No focal position defined, every SNP contributes equally to haplotype clustering.")
      focal.pos <- NA
   }
   mu.snp <- apply(input, 2, mean, na.rm = T)
   p.snp <- mu.snp
   m.mat <- matrix(, nr = n.chr, nc = n.snp)
   for (i in 1:n.snp){
      m.mat[,i] <- (input[,i] - mu.snp[i])/sqrt(p.snp[i] * (1 - p.snp[i]))
   }
   flag.monomorphic <- which(p.snp == 0 | p.snp == 1)
   print(paste("Input haplotypes include ", length(flag.monomorphic), " monomorphic SNPs.", sep=""))
   if (length(flag.monomorphic) > 0) m.mat[,flag.monomorphic] <- 0
   c.array <- array(, dim = c(n.chr, n.chr, n.snp))
   for (i in 1:n.chr){
      for (j in i:n.chr){
         c.array[i,j,] <- (m.mat[i,] - m.mat[j,])^2
         c.array[j,i,] <- (m.mat[j,] - m.mat[i,])^2
         na.flag <- which(is.na(c.array[i,j,]))
         c.array[i,j,na.flag] <- NA
         c.array[j,i,na.flag] <- NA
      }
   }
   print("Finished standardizing the input matrix.")
   
   miss.snp <- snp.miss
   # Missingness weight #
   if (is.null(snp.miss)){
      miss.snp <- double(n.snp)
      for (i in 1:n.snp){
         miss.snp[i] <- mean(is.na(input[,i]))
      }
   }
   
   # Distance from focal/centre #
   if (focal == T){
      if (distance.measure == "physical") dis.snp <- exp(- abs(position - focal.pos)/100000) # for physical distance in bp
      if (distance.measure == "genetic") dis.snp <- exp(- abs(position - focal.pos)/0.15)   # for genetic distance in cM
      if (focal.weight != -1 & is.na(focal.flag) == F){
         temp.1 <- sum(dis.snp[-focal.flag])
         temp.2 <- focal.weight * temp.1 / (temp.1 - focal.weight * temp.1)
         dis.snp[focal.flag] <- temp.2 * temp.1 
      }
      if (focal.weight != -1 & is.na(focal.flag)) print("Index of focal SNP not provided, so center of region used instead. Enforced weighting for focal SNP is thus ignored.")
   }
   if (focal == F){
      dis.snp <- 1/n.snp
   }   

   # Weight for each SNP #
   w.snp <- dis.snp * (1 - miss.snp)
   w.snp <- w.snp / sum(w.snp) * n.snp
   print("Finished calculating the weights for the SNPs.")
   
   # Calculating the X matrix #   
   x.mat.temp <- matrix(, nr = n.chr, nc = n.chr)
   penalty.mat <- matrix(0, nr = n.chr, nc = n.chr)
   for (i in 1:n.chr){
      for (j in i:n.chr){
         flag.miss <- which(is.na(c.array[i,j,]))
         n.miss.snp <- length(flag.miss)
         if (n.miss.snp == 0){
            x.mat.temp[i,j] <- n.snp - sqrt(sum(w.snp * c.array[i,j,]))
            x.mat.temp[j,i] <- x.mat.temp[i,j]
         }
         if (n.miss.snp > 0){
            x.mat.temp[i,j] <- (n.snp - n.miss.snp) - sqrt(sum(w.snp[-flag.miss] * c.array[i,j,-flag.miss]))
            x.mat.temp[j,i] <- x.mat.temp[i,j]
            penalty.mat[i,j] <- (sum(w.snp[flag.miss]))^(-1) * sum(w.snp[flag.miss] * p.snp[flag.miss] * (1 - p.snp[flag.miss]))/(n.snp - n.miss.snp)
            penalty.mat[j,i] <- penalty.mat[i,j]
         }
      }
   }
   print("Finished calculating the covariance and penalty matrices")
   diag.x.mat.temp <- matrix(0, nr = n.chr, nc = n.chr)
   diag(diag.x.mat.temp) <- diag(x.mat.temp)
   chol.x.mat.temp <- chol(diag.x.mat.temp)
   x.mat.corr <- solve(chol.x.mat.temp) %*% x.mat.temp %*% solve(chol.x.mat.temp)
   x.mat <- x.mat.corr - penalty.mat
   diag(x.mat) <- 1
   
   print("Penalized correlation matrix obtained, begin eigen-decomposition.")
   eigen.out <- eigen(x.mat)
   e.values <- eigen.out$values
   valid.evalues <- rev(sort(e.values))[-1]
   n.iter <- 1 + length(which(valid.evalues > mean(valid.evalues) + (2 + tolerance) * sd(valid.evalues)))
   print(paste("Automatic detection found ", n.iter, " possible canonical haplotypes", sep=""))
   if (n.hap == "auto") n.hap <- n.iter
   print(paste("Finding ", n.hap, " canonical haplotypes...", sep=""))
      
   ### First canonical haplotype ###
   hap.group <- double(n.chr)
   hap.ranking <- double(n.chr)
   eigen.score <- eigen.out$vectors[,2]
   rank.eigen.score <- rank(order(eigen.score))
   diff.sort.eigen.score <- diff(sort(eigen.score))
   cutoff <- median(diff.sort.eigen.score) + tolerance * diff(quantile(diff.sort.eigen.score, prob = c(0.25, 0.75)))
   flag.temp <- which(diff.sort.eigen.score > cutoff)
   max.common <- which(diff(c(0, flag.temp, n.chr)) == max(diff(c(0, flag.temp, n.chr))))[1]
   if (max.common == 1){ 
      flag <- 1:flag.temp[1] 
      flag.1 <- rank.eigen.score[1:max(flag)]
   }
   if (max.common != 1){ 
      flag <- (flag.temp[max.common-1]+1):min(c(flag.temp[max.common], n.chr), na.rm = T) 
      flag.1 <- rank.eigen.score[min(flag):max(flag)]
   }
   hap.group[flag.1] <- 1   
   hap.ranking[flag.1] <- 1:length(flag)
   x.remain <- x.mat[-flag.1, -flag.1]
   index.remain <- (1:n.chr)[-flag.1]
   end.1 <- length(flag)
   k <- 2
   n.chr.k <- dim(x.remain)[1]
   print(paste("Canonical haplotype 1 found with ", length(flag), " entries.", sep=""))

   ### Stepwise haplotype sorting ### 
   sorted.input <- input
   sorted.input[1:end.1,] <- input[flag.1,]
   stepwise.ranking <- double(n.chr)
   stepwise.ranking[flag.1] <- 1:length(flag)
   s <- 3
   n.chr.s <- n.chr - length(flag.1)
   index.remain.s <- (1:n.chr)[-flag.1]
   end.2 <- end.1
   eigen.scores <- eigen.out$vectors[index.remain.s,s] 
   rank.eigen.scores <- rank(order(eigen.scores))
   diff.sort.eigen.scores <- diff(sort(eigen.scores))
   cutoff <- median(diff.sort.eigen.scores) + (tolerance + 1) * diff(quantile(diff.sort.eigen.scores, prob = c(0.25, 0.75)))
   flag.temp <- which(diff.sort.eigen.scores > cutoff)            
   while (n.chr.s > 10 & length(flag.temp) > 0){
      if (flag.temp[1] >= (n.chr.s - max(flag.temp))){
         flag <- 1:flag.temp[1]
         flag.1 <- rank.eigen.scores[min(flag):max(flag)]
      }
      if (flag.temp[1] < (n.chr.s - max(flag.temp))){
         flag <- (max(flag.temp)+1):n.chr.s
         flag.1 <- rank.eigen.scores[min(flag):max(flag)]
      }
      sorted.input[(end.2+1):(end.2+length(flag.1)),] <- input[index.remain.s[flag.1],]
      n.chr.s <- n.chr.s - length(flag.1)
      stepwise.ranking[index.remain.s[flag.1]] <- (end.2+1):(end.2+length(flag))
      index.remain.s <- index.remain.s[-flag.1]
      end.2 <- end.2 + length(flag.1)
      s <- s + 1
      eigen.scores <- eigen.out$vectors[index.remain.s,s] 
      rank.eigen.scores <- rank(order(eigen.scores))
      diff.sort.eigen.scores <- diff(sort(eigen.scores))
      cutoff <- median(diff.sort.eigen.scores) + (tolerance + 1) * diff(quantile(diff.sort.eigen.scores, prob = c(0.25, 0.75)))
      flag.temp <- which(diff.sort.eigen.scores > cutoff)            
   }
   stepwise.ranking[index.remain.s] <- (end.2+1):n.chr

   ### Subsequent canonical haplotypes ###
   while (k <= n.hap & n.chr.k > 5){
      eigen.outk <- eigen(x.remain)
      eigen.scorek <- eigen.outk$vectors[,2]   
      rank.eigen.scorek <- rank(order(eigen.scorek))
      diff.sort.eigen.scorek <- diff(sort(eigen.scorek))
      cutoff <- median(diff.sort.eigen.scorek) + tolerance * diff(quantile(diff.sort.eigen.scorek, prob = c(0.25, 0.75)))
      flag.temp <- which(diff.sort.eigen.scorek > cutoff)
      max.common <- which(diff(c(0, flag.temp, n.chr.k)) == max(diff(c(0, flag.temp, n.chr.k))))[1]
      if (max.common == 1){ 
         flag <- 1:flag.temp[1] 
         flag.1 <- rank.eigen.scorek[1:max(flag)]
      }
      if (max.common != 1){ 
         flag <- (flag.temp[max.common-1]+1):min(c(flag.temp[max.common], n.chr.k), na.rm = T) 
         flag.1 <- rank.eigen.scorek[min(flag):max(flag)]
      }
      hap.group[index.remain[flag.1]] <- k
      hap.ranking[index.remain[flag.1]] <- (end.1+1):(end.1+length(flag))
      x.remain <- x.remain[-flag.1, -flag.1]
      index.remain <- index.remain[-flag.1]
      end.1 <- end.1 + length(flag)
      n.chr.k <- dim(x.remain)[1]
      print(paste("Canonical haplotype ", k, " found with ", length(flag), " entries.", sep=""))
      k <- k + 1
   }
   flag.0 <- which(hap.group == 0)
   if (length(flag.0) > 0){
      hap.group[flag.0] <- -1
      last.ranking <- rank(order(eigen.score[flag.0]))
      hap.ranking[flag.0[last.ranking]] <- (end.1+1):n.chr
   }
   
   ### Simple sorting ### 
   print("Performing simple sorting...")
   ranking <- rev(rank(order(eigen.out$vectors[,2])))   
   simple.sort <- input[ranking,]
   
   ### k-means sorting ### 
   print("Performing k-means sorting...")
   kmeans.out <- kmeans(eigen.out$vectors[,2:n.hap], n.hap)
   kmeans.ranking <- double(n.chr)
   end.k <- 0
   kmeans.unsorted <- input      
   kmeans.sort <- {}
   for (k in 1:n.hap){
      flag.k <- which(kmeans.out$cluster == k)
      kmeans.unsorted[flag.k,] <- k
      x.temp <- x.mat[flag.k, flag.k] 
      eigen.kmeans <- eigen(x.temp)
      ranking.k.temp <- rank(order(eigen.kmeans$vectors[,2]))
      kmeans.ranking[(end.k+1):(end.k + length(flag.k))] <- flag.k[ranking.k.temp]
      end.k <- end.k + length(flag.k)
   }

   ### Matching to canonical haplotypes ### $
   print("Creating mosaics of unmapped haplotypes...")
   haplosim.sort <- input
   canonical.hap <- matrix(, nr = n.hap, nc = n.snp)
   for (k in 1:n.hap){
      flag.hap <- which(hap.group == k)
      canonical.hap[k,] <- apply(input[flag.hap,], 2, mean)
      haplosim.sort[flag.hap,] <- k   
   }   
   canonical.m <- matrix(, nr = n.hap, nc = n.snp)
   for (i in 1:n.snp){
      canonical.m[,i] <- (canonical.hap[,i] - mu.snp[i])/(sqrt(p.snp[i] * (1 - p.snp[i])))
   }
   unmapped.chr <- which(hap.group == -1)
   n.unmap.chr <- length(unmapped.chr)
   c.unmap <- array(, dim = c(n.hap, n.unmap.chr, n.snp))
   for (k in 1:n.hap){
      for (j in 1:n.unmap.chr){
         c.unmap[k,j,] <- (canonical.hap[k,] - m.mat[unmapped.chr[j],])^2
         na.flag <- which(is.na(c.unmap[k,j,]))
         c.unmap[k,j,na.flag] <- NA
      }
   }   
   mosaic.hap <- matrix(, nr = n.unmap.chr, nc = n.snp)
   mosaic.window.length <- n.snp - sticky + 1
   sticky.center <- ceiling(sticky/2)
   for (j in 1:n.unmap.chr){
      mosaic.window <- double(mosaic.window.length)
      for (i in 1:mosaic.window.length){
         temp.mosaic <- apply(w.snp[i:(i+sticky-1)]*c.unmap[,j,i:(i+sticky-1)], 1, mean, na.rm = T) 
         mosaic.window[i] <- which(temp.mosaic == min(temp.mosaic))[1]
      }
      smooth.window.length <- mosaic.window.length - sticky + 1
      for (i in 1:smooth.window.length){
         temp.smooth <- mosaic.window[i:(i+sticky-1)]
         if (temp.smooth[sticky.center] != temp.smooth[sticky.center-1] & temp.smooth[sticky.center] != temp.smooth[sticky.center+1]){
            tab.temp <- tabulate(temp.smooth)
            if (temp.smooth[sticky.center-1] == temp.smooth[sticky.center+1]){ 
               mosaic.window[(i:(i+sticky-1))[sticky.center]] <- temp.smooth[sticky.center-1]
            }
            if (temp.smooth[sticky.center-1] != temp.smooth[sticky.center+1]){
               mode.hap <- which(tab.temp == max(tab.temp))
               if (length(mode.hap) == 1) mosaic.window[(i:(i+sticky-1))[sticky.center]] <- mode.hap
               if (length(mode.hap) > 1 & is.element(temp.smooth[sticky.center-1], mode.hap)) mosaic.window[(i:(i+sticky-1))[sticky.center]] <- temp.smooth[sticky.center-1]
               if (length(mode.hap) > 1 & is.element(temp.smooth[sticky.center-1], mode.hap) == F & is.element(temp.smooth[sticky.center+1], mode.hap)) mosaic.window[(i:(i+sticky-1))[sticky.center]] <- temp.smooth[sticky.center+1]
               if (length(mode.hap) > 1 & is.element(temp.smooth[sticky.center-1], mode.hap) == F & is.element(temp.smooth[sticky.center+1], mode.hap) == F){
                  print("No clear mosaic membership, unable to smooth effectively.")
               }
            }
         }
      }
      mosaic.hap.j <- c(rep(mosaic.window[1], floor(sticky/2)), mosaic.window, rep(mosaic.window[mosaic.window.length], floor(sticky/2)))
      mosaic.hap[j, ] <- mosaic.hap.j
   }
   for (j in 1:n.unmap.chr){
      haplosim.sort[unmapped.chr[j],] <- mosaic.hap[j,]
   }
   haplosim.unsorted <- haplosim.sort
   if (n.chr >= n.snp){
      for (i in 1:n.snp){
         na.flag <- which(miss.matrix[,i] == 1)
         if (length(na.flag) > 0) {
            kmeans.unsorted[na.flag,i] <- 0
            haplosim.unsorted[na.flag,i] <- 0
         }
      }
   }
   if (n.chr < n.snp){
      for (i in 1:n.chr){
         na.flag <- which(miss.matrix[i,] == 1)
         if (length(na.flag) > 0) {
            kmeans.unsorted[i,na.flag] <- 0
            haplosim.unsorted[i,na.flag] <- 0
         }
      }
   }
   kmeans.sort <- kmeans.unsorted[kmeans.ranking,]   
   haplosim.sort[hap.ranking,] <- haplosim.sort      
   
   output <- list(unsorted = input, simple.sort = simple.sort, stepwise.sort = sorted.input, kmeans.unsorted = kmeans.unsorted, kmeans.sort = kmeans.sort, haplosim.unsorted = haplosim.unsorted, haplosim.sort = haplosim.sort, hap.group = hap.group, hap.ranking = hap.ranking, simple.ranking = ranking, stepwise.ranking = stepwise.ranking, kmeans.ranking = kmeans.ranking, focal.flag = focal.flag, focal.pos = focal.pos, n.chr = n.chr, n.snp = n.snp, eigenvalues = eigen.out$values, eigenvectors = eigen.out$vectors, n.hap = n.hap)
   print("Done! Recommend the use of HAPVISUAL to visualize the results.")
   output
}


