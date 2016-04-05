### PARALLEL implementation of the Threshold Edge Removal Algorithm (TERA) 
#   described in Carvalho,  Santos & Silveira,  2014 ###
# Copyleft (or the one to blame): Carvalho,  LMF (2014)
# Last update: 30/05/2014  
##################################################
source("graph_entropy.R")
tera_step <- function(t,  M,  export = "FALSE",  path){
  # M is the distance matrix of the original weighted graph
  # s_min is the vector with thresholds (\sigma)
  # export is boolean for whether to create files with output,  which will be placed in 'path'
  require(igraph)
  M_norm <- M/max(M)
  V <-  ncol(M) # number of vertices
  cat("S-min = ", t, "\n")
  tim <- M
  tim[.Internal(which(M_norm <= t))] <- 1
  tim[.Internal(which(M_norm > t))] <- 0
  diag(tim) <- 0
  G <- graph.adjacency(tim, mode = c("undirected"))#transforms each matrix into a graph
  if(export == "TRUE"){
    write.table(tim, file = paste(path, "matrix_", t, ".txt"), sep = "\t")
    write.graph(G, file = paste(path, "pajek/graph", t, ".net", sep = ""), format = "pajek")      
  }
  Bet  <- betweenness(G) #calculates the betweenness for every connected edge in the graph
  K <- mean(igraph::degree(G))
  Gclu <- transitivity(G, type = c("global"))
  Aclu <- transitivity(G, type = c("average"))
  D <- diameter(G)
  L <- average.path.length(G) 
  ent <- graph.entropy(tim)
  AS <- assortativity.degree(G)    
  ASG <- assortativity.nominal(G, types = country.vec)
  ## small-world 'corrected' indexes
  E <- ecount(G)
  cC <- 2*E/(V*(V-1)) 
  cK <- cC*(V-1)
  cL <- log(V)/log(cK)    
  return(list(Matrix = tim, betweenness = Bet, indexes = c(K = K, cK = cK, Gclu = Gclu,
                                             Aclu = Aclu, cC = cC, D = D, L = L,
                                             E = E, cL = cL, ent = ent, AS = AS,
                                             ASG = ASG)))
}
#######################################################
tera.parallel <- function(M, s_min, export = "FALSE", path, cores = 4){
  require(doMC)
  registerDoMC()
  options(cores = cores)
  out <-  foreach(i = s_min) %dopar% tera_step(i, M, export, path)
  return(out)
}
