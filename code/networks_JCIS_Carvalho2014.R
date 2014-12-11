# Copyleft (or the one to blame): Carvalho, LMF (2013)
# created in April 29th, 2013  
# last update: July 7th, 2014
require(igraph)
source("tera.R")
source("tera_parallel.R")
#################
## Load distance matrices and metadata
M.NT <- as.matrix(read.table("../data//input//MAT_NT.txt", sep = ";"))
M.AA <- as.matrix(read.table("../data/input//MAT_AA.txt", sep = ";"))
M.ANT <- as.matrix(read.table("../data/input//MAT_ANT.txt", sep = ";"))
country.vec <- read.table("../data//input//countries_O.txt")$V1 # countries
###
s_min <- seq(0, 1, 0.01) # build 101 thresholds

exporta <- TRUE
output.NT <- tera(M.NT, s_min, exporta, path = "../data//output/NT/")
output.AA <- tera(M.AA, s_min, exporta, path = "../data//output/AA/")
output.ANT <- tera(M.ANT, s_min, exporta, path = "../data//output/ANT/")
#
adj.mats.NT <- output.NT[[1]]
bet.NT <- output.NT[[2]]
indexes.NT <- output.NT[[3]]
#
adj.mats.AA <- output.AA[[1]]
bet.AA <- output.AA[[2]]
indexes.AA <- output.AA[[3]]
#
adj.mats.ANT <- output.ANT[[1]]
bet.ANT <- output.ANT[[2]]
indexes.ANT <- output.ANT[[3]]

######################
binary.entropy <- function(p) -p*log2(p) - (1-p)*log2(1-p)
pNT <- apply(adj.mats.NT, 3, mean)
pAA <- apply(adj.mats.AA, 3, mean)
pANT <-apply(adj.mats.ANT, 3, mean)   
ENT <- binary.entropy(pNT)
EAA <- binary.entropy(pAA)
EANT <- binary.entropy(pANT)
plot(s_min, ENT, ylab = expression(H(X)),
     xlab = expression(sigma), type = "l", lwd = 2, xlim = c(0, 1.2))
lines(s_min, EAA, col = 2)
lines(s_min, EANT, col = 3)
legend(x = "topright", legend = c("Nucleotides", "Aminoacids", "Antigenic"),
       col = c(1, 2, 3), lwd = 3, cex = .7, bty  =  "n" )
######################
#exporting indexes table
write.table(indexes.NT, file = "../data//output/tabela_du_pain_NT.txt",
            row.names = FALSE, sep = "\t")
write.table(indexes.AA, file = "../data//output/tabela_du_pain_AA.txt",
            row.names = FALSE, sep = "\t")
write.table(indexes.ANT, file = "../data//output/tabela_du_pain_ANT.txt",
            row.names = FALSE, sep = "\t")
######################
# Exporting vertex betweenness tables for plotting purposes
namiz <- paste("s_min=", s_min, sep = "")
dbNT <- data.frame(matrix(unlist(bet.NT), 167, 101)); names(dbNT) <- namiz
dbAA <- data.frame(matrix(unlist(bet.AA), 167, 101)); names(dbAA) <- namiz
dbANT <- data.frame(matrix(unlist(bet.ANT), 167, 101)); names(dbANT) <- namiz
## Creating the betweenness tables
write.table(dbNT, file =  "../data//output/NT/Betweenness_NT.txt", sep = "\t")
write.table(dbAA, file = "../data//output/AA/Betweenness_AA.txt", sep = "\t")
write.table(dbANT, file = "../data//output/ANT/Betweenness_ANT.txt", sep = "\t")
######################
## Plotting
plot(ecdf(M.NT/max(M.NT)), main = "Cumulative distribution of normalized distances",
     xlab = expression(D/max(D)))
lines(ecdf(M.AA/max(M.AA)), col = 2)
lines(ecdf(M.ANT/max(M.ANT)), col = 3)
legend(x = "topleft", legend = c("Nucleotides", "Aminoacids", "Antigenic"),
       col = c(1, 2, 3), lwd = 3, cex = 1, bty = "n" )
#
pdf("../figures/diameter.pdf")
plot(s_min, indexes.NT$D, type = "l", lwd = 3, xlab = expression(sigma),
     ylab = "L", main = "Diameter", ylim = c(1, 9))
lines(s_min, indexes.AA$D, lwd = 3, col = 2)
lines(s_min, indexes.ANT$D, lwd = 3, col = 3)
legend(x = "topright", legend = c("Nucleotides", "Aminoacids", "Antigenic"),
       col = c(1, 2, 3), lwd = 3, cex = 1, bty  =  "n" )
dev.off()
#
pdf("../figures/degree.pdf")
plot(s_min, indexes.NT$K, type = "l", lwd = 3, xlab = expression(sigma),
     ylab = "<k>", main = "Degree")
lines(s_min, indexes.AA$K, lwd = 3, col = 2)
lines(s_min, indexes.ANT$K, lwd = 3, col = 3)
legend(x = "topleft", legend = c("Nucleotides", "Aminoacids", "Antigenic"),
       col = c(1, 2, 3), lwd = 3, cex = 1, bty  =  "n" )
dev.off()
#
pdf("../figures/clustering.pdf")
plot(s_min, indexes.NT$Aclu, type = "l", lwd = 3, xlab = expression(sigma),
     ylab = "<c>", main = "Clustering coefficient", ylim = c(.82, 1))
lines(s_min, indexes.AA$Aclu, lwd = 3, col = 2)
lines(s_min, indexes.ANT$Aclu, lwd = 3, col = 3)
legend(x = "bottomright", legend = c("Nucleotides", "Aminoacids", "Antigenic"),
       col = c(1, 2, 3), lwd = 3, cex = 1, bty  =  "n" )
dev.off()
#
pdf("../figures/entropy.pdf")
plot(s_min, indexes.NT$ent, type = "l", lwd = 3, xlab = expression(sigma),
     ylab = "Graph entropy", main = "Entropy")
lines(s_min, indexes.AA$ent, lwd = 3, col = 2)
lines(s_min, indexes.ANT$ent, lwd = 3, col = 3)
legend(x = "topleft", legend = c("Nucleotides", "Aminoacids", "Antigenic"),
       col = c(1, 2, 3), lwd = 3, cex = 1, bty  =  "n" )
dev.off()
#
pdf("../figures/spath.pdf")
plot(s_min, indexes.NT$L, type = "l", lwd = 3, xlab = expression(sigma), ylim = c(1, 3),
     ylab = "<l>", main = "Shortest path")
lines(s_min, indexes.AA$L, lwd = 3, col = 2)
lines(s_min, indexes.ANT$L, lwd = 3, col = 3)
legend(x = "topright", legend = c("Nucleotides", "Aminoacids", "Antigenic"),
       col = c(1, 2, 3), lwd = 3, cex = 1, bty  =  "n" )
dev.off()
#
pdf("../figures/degree_assort.pdf")
plot(s_min, indexes.NT$AS, type = "l", lwd = 3, xlab = expression(sigma),
     ylab = "Assortativity", main = "Degree assortativity")
lines(s_min, indexes.AA$AS, lwd = 3, col = 2)
lines(s_min, indexes.ANT$AS, lwd = 3, col = 3)
legend(x = "topright", legend = c("Nucleotides", "Aminoacids", "Antigenic"),
       col = c(1, 2, 3), lwd = 3, cex = 1, bty  =  "n")
dev.off()
#
pdf("../figures/country_assort.pdf")
plot(s_min, indexes.NT$ASG, type = "l", lwd = 3, xlab = expression(sigma),
     ylab = "Assortativity", main = "Country nominal assortativity")
lines(s_min, indexes.AA$ASG, lwd = 3, col = 2)
lines(s_min, indexes.ANT$ASG, lwd = 3, col = 3)
legend(x = "topright", legend = c("Nucleotides", "Aminoacids", "Antigenic"),
       col = c(1, 2, 3), lwd = 3, cex = 1, bty  =  "n")
dev.off()
#
pdf("../figures/betweenness.pdf")
plot(s_min, lapply(bet.NT, mean), type = "l", lwd = 3, xlab = expression(sigma),
     ylab = "Average betweenness", main = "Betweenness", ylim = c(1, 120))
lines(s_min, lapply(bet.AA, mean), lwd = 3, col = 2)
lines(s_min, lapply(bet.ANT, mean), lwd = 3, col = 3)
legend(x = "topright", legend = c("Nucleotides", "Aminoacids", "Antigenic"),
       col = c(1, 2, 3), lwd = 3, cex = 1, bty  =  "n")
dev.off()
#########
# Small-world analysis
#########
plot(s_min, indexes.NT$cK, type = "l", lwd = 3, lty = 2, xlab = expression(sigma),
     ylab = "<k'>", main = "Small-world -- <k'>")
lines(s_min, indexes.NT$K, lwd = 3)
lines(s_min, indexes.AA$cK, lwd = 3, lty = 2, col = 2)
lines(s_min, indexes.AA$K, lwd = 3, col = 2)
lines(s_min, indexes.ANT$cK, lwd = 3, lty = 2, col = 3)
lines(s_min, indexes.ANT$K, lwd = 3, col = 3)
legend(x = "bottomright", legend = c("Nucleotides", "Aminoacids", "Antigenic"),
       col = c(1, 2, 3), lwd = 3, cex = 1, bty  =  "n")
#
plot(s_min, indexes.NT$cC, type = "l", lwd = 3, lty = 2, xlab = expression(sigma),
     ylab = "<c'>", main = "Small-world -- <c'>")
lines(s_min, indexes.AA$cC, lwd = 3, lty = 2, col = 2)
lines(s_min, indexes.ANT$cC, lwd = 3, lty = 2, col = 3)
legend(x = "bottomright", legend = c("Nucleotides", "Aminoacids", "Antigenic"),
       col = c(1, 2, 3), lwd = 3, cex = 1, bty  =  "n")
#
plot(s_min, indexes.NT$cL, type = "l", lwd = 3, lty = 2, xlab = expression(sigma),
     ylab = "<l'>", main = "Small-world -- <l'>")
lines(s_min, indexes.AA$cL, lwd = 3, lty = 2, col = 2)
lines(s_min, indexes.ANT$cL, lwd = 3, lty = 2, col = 3)
legend(x = "topright", legend = c("Nucleotides", "Aminoacids", "Antigenic"),
       col = c(1, 2, 3), lwd = 3, cex = 1, bty  =  "n")