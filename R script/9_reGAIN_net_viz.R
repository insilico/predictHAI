#
#
# reGAIN network visualization
# 
# By Saeid Parvandeh
#
# ----------------------------------------------------------
rm(list=ls())

# load("baylor_regain_resid.RData")
rownames(baylor_regain_resid) <- colnames(baylor_regain_resid)

# Selected top 200 genes reGAIN
Bay_reGAIN_t200 <- baylor_regain_resid[bay_tg, bay_tg]


# Function to estimate the threshold and create adjacency matrix 
library(igraph)
density_fun <- function(thresh.points,reGAIN)
{
  gAdjBool <- ifelse(abs(reGAIN) > thresh.points, 1, 0)
  diag(gAdjBool) <- 0
  g.DB <- graph.adjacency(gAdjBool, "undirected")
  edges <- ecount(g.DB)
  full <- (ncol(gAdjBool)*(ncol(gAdjBool)-1))/2
  return (edges/full) 
}

# # Here, we estimate the threshold using reGAIN matrix
thresh.points <- seq(0, max(abs(Bay_reGAIN_t200)), .05)
thresh.vec <- sapply(thresh.points, density_fun, Bay_reGAIN_t200)
# plot(thresh.points, thresh.vec)

# We are using estimated threshold here to create adjacency matrices 
# Baylor reGAIN to adjacency
Baylor_reGAIN_Adj <- ifelse(abs(Bay_reGAIN_t200) > 4.5, 1, 0)
diag(Baylor_reGAIN_Adj) <- 0

# Visulaize reGAIN network
bay_graph <- graph.adjacency(Baylor_reGAIN_Adj, "undirected")
plot(bay_graph, vertex.size=2)

