### Implementation of the Threshold Edge Removal Algorithm (TERA) 
#   described in Carvalho, Silveira & Santos, 2014 ###
# Copyleft (or the one to blame): Carvalho,  LMF (2014)
# Last update: 30/05/2014  
# TODO: make it more general and create documentation
##################################################
tera <- function(M, s_min, export = "FALSE", path){
  # M is the distance matrix of the original weighted graph
  # s_min is the vector with thresholds (\sigma \in [0, 1]) 
  # export is boolean for whether to create files with output,  which will be placed in 'path'
  source("graph_entropy.R")
  require(igraph)
  M_norm <- M/max(M)
  V <- ncol(M) # number of vertices
  # empty vectors of the appropriate size
  K <- cK <- Gclu <- cC <- Aclu <- D <- L <- E <- cL <- ent <- AS <- ASG <- rep(NA, length(s_min)) 
  TA <- array(NA, dim = c(dim(M)[1], dim(M)[2], length(s_min)))
  B <- vector(length(s_min), mode = "list")# empty objects to store the results
  for (t in s_min){
    i <- match(t, s_min) 
    cat("S-min =", t, "\n")
    tim <- matrix(NA, ncol = V, nrow = V )
    tim[.Internal(which(M_norm <= t))] <- 1
    tim[.Internal(which(M_norm > t))] <- 0
    diag(tim) <- 0
    TA[, , i] <- tim
    G <- graph.adjacency(tim, mode = c("undirected")) # transforms each matrix into a graph
    if(export == "TRUE"){
      write.table(tim, file = paste(path, "matrix_", s_min[i], ".txt"), sep = "\t")
      write.graph(G, file = paste(path, "pajek/graph", s_min[i], ".net", sep = ""),
                  format = "pajek")      
    }
    B[[i]] <- betweenness(G) #calculates the betweenness for every connected edge in the graph
    K[i] <- mean(igraph::degree(G))
    Gclu[i] <- transitivity(G, type = c("global"))
    Aclu[i] <- transitivity(G, type = c("average"))
    D[i] <- diameter(G)
    L[i] <- average.path.length(G) 
    ent[i] <- graph.entropy(tim)
    AS[i] <- assortativity.degree(G)    
    ASG[i] <- assortativity.nominal(G, types = country.vec)
    ## 'small-world-corrected' indexes
     E[i] <- ecount(G)
    cC[i] <- 2*E[i]/(V*(V-1)) 
    cK[i] <- cC[i]*(V-1)
    cL[i] <- log(V)/log(cK[i])    
  }
  return(list(matrices = TA,
              betweenness = B,
              indexes = data.frame(nat = max(M)*s_min, t = s_min, K, cK, Gclu, # nat = threshold in the original scale
                                   Aclu, cC,D, L, E, cL,
                                   ent, AS, ASG)))
}
#######################################################
## Function to extract influent nodes (vertices),  with high betweenness
extract.influent <- function(B, cp, alpha) {
  # B is the list with betweeness values for all sigmas
  # cp is a vector with all desired critical points
  # alpha is the quantile from which "influent" points are drawn
 influent <- vector(length(cp), mode = "list") # empty object to store the results
 subB <- B[cp] # subsampling the whole list
 whichinfs <- lapply (subB,
          function(x, alpha)  which(x >= quantile(x, probs = alpha)),
                                    alpha = alpha)
 # which positions are greater the F(alpha)?
 for (i in 1:length(cp)){
   influent [[i]] <- subB[[i]][ whichinfs[[i]] ]
 }
  return(influent)}