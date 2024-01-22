### PCA analysis

# libraries
library(devtools); library(bigutilsr)

#Use bam file
# - 2: PCA
# - genotype likelihoods (from .bam file)
# - pipe to PCANGSD
# - PCA admixture
# - R plotting


#PCAngsd
#####################
# https://github.com/Rosemeis/pcangsd/tree/master

# run in python to get these files
# pcangsd -h
# Genotype likelihood file in Beagle format with 2 eigenvectors using 64 threads
# pcangsd -b input.beagle.gz -e 2 -t 64 -o pcangsd
# Outputs by default log-file (pcangsd.log) and covariance matrix (pcangsd.cov)
# PLINK files (using file-prefix, *.bed, *.bim, *.fam)
# pcangsd -p input.plink -e 2 -t 64 -o pcangsd
# Perform PC-based selection scan and estimate admixture proportions
# pcangsd -b input.beagle.gz -e 2 -t 64 -o pcangsd --selection --admix
# Outputs the following files:
# log-file (pcangsd.log)
# covariance matrix (pcangsd.cov)
# selection statistics (pcangsd.selection)
# admixture proportions (pcangsd.admix.3.Q)
# ancestral allele frequencies (pcangsd.admix.3.F)

C <- as.matrix(read.table("pcangsd.cov")) # Reads estimated covariance matrix
D <- as.matrix(read.table("pcangsd.selection")) # Reads PC based selection statistics

# Plot PCA plot
e <- eigen(C)
plot(e$vectors[,1:2], xlab="PC1", ylab="PC2", main="PCAngsd")

# Obtain p-values from PC-based selection scan
p <- pchisq(D, 1, lower.tail=FALSE)











#####################

####################
# Run pcadapt scan #
####################
#
#	args[1]: zscores matrix
#	args[2]: prefix for output files
#

args = commandArgs(trailingOnly=TRUE)

zscores <- read.table(args[1])[,1]
K <- ncol(zscores)

# For one component only
if (K == 1) {
  d2 <- (zscores - median(zscores))^2
} else {
  d2 <- dist_ogk(zscores)
}

write.table(d2, file=paste0(args[2], ".pcadapt.test.txt"), quote=F, row.names=F, col.names=F)
write.table(pchisq(d2, df=K, lower.tail=F), file=paste0(args[2], ".pcadapt.pval.txt"), quote=F, row.names=F, col.names=F)



# angsd
###################
# angsd
# https://github.com/ANGSD/angsd
# https://www.popgen.dk/angsd/index.php/PCA_MDS#MDS/PCA_using_R

## MDS
name <- "angsdput.ibsMat"
m <- as.matrix(read.table(name))
mds <- cmdscale(as.dist(m))
plot(mds,lwd=2,ylab="Dist",xlab="Dist",main="multidimensional scaling",col=rep(1:3,each=10))

name <- "angsdput.covMat"
m <- as.matrix(read.table(name))
e <- eigen(m)
plot(e$vectors[,1:2],lwd=2,ylab="PC 2",xlab="PC 2",main="Principal components",col=rep(1:3,each=10),pch=16)

## heatmap / clustering / trees
name <- "angsdput.ibsMat" # or covMat
m <- as.matrix(read.table(name))
#heat map
heatmap(m)
#neighbour joining
plot(ape::nj(m))
plot(hclust(dist(m), "ave")
     
#####################
