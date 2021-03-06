graph.entropy<-function(adj){
  graph <- graph.adjacency(adj, mode = c("undirected"))
  Ks <- (igraph::degree(graph))
  mean.k <- mean(Ks)
  bins <- sort(unique(Ks))-1
  q.k <- table(Ks)/nrow(adj) 
  entropy <- -sum(q.k*log2(q.k))
  return(entropy)}
###
graphEntropy <- function(adj, type = "SoleValverde") {
  if (type == "SoleValverde") {
    return(graphEntropySoleValverde(adj))
  }
  else {
    return(graphEntropyWang(adj))
  }
}
#
graphEntropySoleValverde <- function(adj) {
  # Calculate Sole & Valverde, 2004 graph entropy
  # Uses Equations 1 and 4
  # First we need the denominator of q(k)
  # To get it we need the probability of each degree
  # First get the number of nodes with each degree
  # Special thanks to Dr. Jason J. Jones, who wrote the code was taken from http://stackoverflow.com/questions/6950791/how-do-i-calculate-the-entropy-of-a-graph
  graph<-graph.adjacency(adj,mode=c("undirected"))
  existingDegrees = (igraph::degree(graph))
  maxDegree = nrow(adj) - 1
  allDegrees = 0:maxDegree
  degreeDist = matrix(0, 3, length(allDegrees)+1) # Need an extra zero prob degree for later calculations
  degreeDist[1, ] = 0:(maxDegree+1)
  for(aDegree in allDegrees) {
    degreeDist[2, aDegree+1] = sum(existingDegrees == aDegree)
  }
  # Calculate probability of each degree
  for(aDegree in allDegrees) {
    degreeDist[3, aDegree+1] = degreeDist[2, aDegree+1]/sum(degreeDist[2, ])
  }
  # Sum of all degrees mult by their probability
  sumkPk = 0
  for(aDegree in allDegrees) {
    sumkPk = sumkPk + degreeDist[2, aDegree+1] * degreeDist[3, aDegree+1]
  }
  # Equivalent is sum(degreeDist[2,] * degreeDist[3,])
  # Now we have all the pieces we need to calculate graph entropy
  graphEntropy = 0
  for(aDegree in 1:maxDegree) {
    q.of.k = ((aDegree + 1)*degreeDist[3, aDegree+2])/sumkPk
    # 0 log2(0) is defined as zero
    if (q.of.k != 0) {
      graphEntropy = graphEntropy + -1 * q.of.k * log2(q.of.k)
    }
  }
  return(graphEntropy)
}
graphEntropyWang <- function(adj) {
  # Calculate Wang, 2008 graph entropy
  # Uses Equation 14
  # bigN is simply the number of nodes
  # littleP is the link probability.  That is the same as graph density calculated by sna with gden().
  bigN = nrow(adj)
  littleP = gden(adj)
  graphEntropy = 0
  if (littleP != 1 && littleP != 0) {
    graphEntropy = -1 * .5 * bigN * (bigN - 1) * (littleP * log2(littleP) + (1-littleP) * log2(1-littleP))
  }
  return(graphEntropy)
}