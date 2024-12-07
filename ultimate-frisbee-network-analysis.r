# Ultimate Frisbee Team Network Analysis
# Author: Tyler WW
# Date: 2023-09-19
# Purpose: Analyze network characteristics of ultimate frisbee teams across regions

# Load required libraries
library(rmarkdown)
library(knitr)
library(igraph)
library(moments)
library(scatterplot3d)
library(corrplot)
library(pso)
library(psych)
library(GPArotation)
library(lavaan)

# Data Import and Initial Graph Creation ------------------------------------

# Read node and edge data
# Note: Replace with your local file path or use a relative path
nodes <- read.csv("STAT 414 Nodes.csv")
edges <- read.csv("STAT 414 Edges.csv")

# Create undirected graph from data
graph <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)

# Basic Network Characteristics -------------------------------------------

# Count vertices and edges
num_vertices <- vcount(graph)
num_edges <- ecount(graph)

# Simplify graph (remove multiple edges)
graph <- simplify(graph)

# Create adjacency matrix
adjacency_matrix <- as.matrix(as_adjacency_matrix(graph))

# Visualization of Adjacency Matrix
png("adjacency_matrix_plot.png", width = 800, height = 800)
corrplot(adjacency_matrix, method = "color", 
         title = "Ultimate Frisbee Team Connection Matrix")
dev.off()

# Network Centrality Analysis --------------------------------------------

# Calculate network centrality measures
degree_centrality <- degree(graph)
knn_values <- knn(graph)$knn

# Degree-Connectivity Relationship
png("degree_connectivity_plot.png", width = 800, height = 600)
plot(x = jitter(log(degree_centrality)), 
     y = jitter(log(knn_values)), 
     xlab = "Log Degree", 
     ylab = "Log Connectivity",
     main = "Degree-Connectivity Relationship")
dev.off()

# Fit linear model for degree-connectivity
degree_connectivity_model <- lm(log(knn_values) ~ log(degree_centrality))

# Network Topology Metrics -----------------------------------------------

# Calculate key network metrics
avg_path_length <- average.path.length(graph)
network_diameter <- diameter(graph)
avg_clustering_coefficient <- transitivity(graph, type = "average")

# Centrality Measures ----------------------------------------------------

# Calculate centrality metrics
betweenness_centrality <- betweenness(graph)
closeness_centrality <- closeness(graph)
eigenvector_centrality <- eigen_centrality(graph)$vector

# Identify top and bottom performers in centrality
top_betweenness <- head(sort(betweenness_centrality, decreasing = TRUE))
top_closeness <- head(sort(closeness_centrality, decreasing = TRUE))
top_eigenvector <- head(sort(eigenvector_centrality, decreasing = TRUE))

# Assortativity Analysis ------------------------------------------------

# Calculate assortativity metrics
nominal_assortativity <- assortativity_nominal(graph, types = as.factor(V(graph)$Region))
degree_assortativity <- assortativity_degree(graph)

# Network Robustness Simulation -----------------------------------------

# Attack simulation function (as defined in original script)
attack_sim <- function(g, m){
  G <- g
  
  p_inf <- rep(NA, vcount(G))
  k <- (sum(degree(G)^2)/vcount(G))/(sum(degree(G))/vcount(G))
  if(is.na(k)){
      k <- 0
  }
  if(k > 2){
    p_inf[1] <- max(components(G)$csize)/vcount(G)
  }else{
    p_inf <- rep(0, vcount(G))
    return(p_inf)
  }
  
  for(i in 2:vcount(graph)){
    if(m == 1){
      node <- names(which.max(degree(G))[1])
    }else{
      node <- sample(V(G)$name,1)
    }
    
    G <- G-vertex(node)
    
    k <- (sum(degree(G)^2)/vcount(G))/(sum(degree(G))/vcount(G))
    if(is.na(k)){
      k <- 0
    }
    if(k > 2){
      p_inf[i] <- max(components(G)$csize)/vcount(G)
    }else{
      p_inf[which(is.na(p_inf))] <- 0
      return(p_inf)
    }
  }
  return(p_inf)
}

# Run attack simulations
attack_targeted <- attack_sim(graph, 1)
attack_random <- rowMeans(replicate(n = 250, attack_sim(graph, 2)))

# Visualize attack tolerance
png("network_attack_tolerance.png", width = 800, height = 600)
plot((1:vcount(graph))/vcount(graph), attack_targeted, 
     type = "l", col = "orchid", lwd = 3,
     xlab = "Fraction of Nodes Removed", 
     ylab = "Proportion of Largest Component",
     main = "Network Robustness under Attack")
lines((1:vcount(graph))/vcount(graph), attack_random, 
      type = "l", col = "forestgreen", lwd = 3)
legend("topright", 
       legend = c("Targeted Attack", "Random Removal"), 
       col = c("orchid", "forestgreen"), 
       lwd = 3)
dev.off()

# Community Detection --------------------------------------------------

# Detect communities using edge betweenness
communities <- cluster_edge_betweenness(graph)

# Generate dendrogram of communities
png("community_dendrogram.png", width = 800, height = 600)
dendPlot(communities, main = "Community Structure Dendrogram")
dev.off()

# Community membership
community_membership <- membership(communities)

# Save key results -----------------------------------------------------
save(list = c("graph", "communities", "community_membership", 
              "degree_centrality", "betweenness_centrality", 
              "closeness_centrality", "eigenvector_centrality"),
     file = "ultimate_frisbee_network_analysis.RData")
