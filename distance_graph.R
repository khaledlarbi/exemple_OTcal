library(igraph)

compute_graph_distance <- function(data, normalisation = TRUE, k = NULL) {
  # data : data.frame avec colonnes "latitude" et "longitude"
  coords <- as.matrix(data[, c("latitude", "longitude")])
  N <- nrow(coords)

  # Matrice de distances euclidiennes
  dist_matrix <- as.matrix(dist(coords))


  # C_max
  C_max <- max(dist_matrix)

  if(normalisation){
   dist_matrix <- dist_matrix/C_max
   C_max <- 1
  }
  # Matrice des poids w_ij = C_max - d_ij
  W <- C_max - dist_matrix
  diag(W) <- 0  # w_ii = 0

  # Optionnel : ne garder que les k plus proches voisins pour chaque point
  if (!is.null(k)) {
    for (i in 1:N) {
      neighbors <- order(dist_matrix[i, ])[1:(N-k-1)]  # garder seulement les k plus proches
      W[i, neighbors] <- 0
    }
  }

  # Symétriser pour graphe non dirigé
  W <- pmax(W, t(W))

  # Création du graphe pondéré
  g <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE, diag = FALSE)

  # Distance de plus court chemin sur le graphe
  D_graph <- distances(g, weights = E(g)$weight)

  return(list(
    weights = W,
    graph_distance = D_graph
  ))
}


dist <- compute_graph_distance(simulated_data)$graph_distance
