library("spsurvey")
library("dplyr")
library("ggplot2")
library("spdep")
library("geoR")


devtools::load_all() #Chargement du package OTcalibration

generate_clusters_var <- function(n_points = 500, n_clusters = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  centers <- matrix(runif(2 * n_clusters), ncol = 2)

  # Taille et poids aléatoires
  cluster_sd <- runif(n_clusters, 0.02, 0.08)
  weights <- runif(n_clusters)
  weights <- weights / sum(weights)
  points_per_cluster <- rmultinom(1, n_points, weights)

  X <- matrix(NA, nrow = n_points, ncol = 2)
  start <- 1
  for (i in 1:n_clusters) {
    n_i <- points_per_cluster[i]
    cluster_points <- cbind(
      rnorm(n_i, mean = centers[i, 1], sd = cluster_sd[i]),
      rnorm(n_i, mean = centers[i, 2], sd = cluster_sd[i])
    )
    end <- start + n_i - 1
    X[start:end, ] <- cluster_points
    start <- end + 1
  }

  X <- pmin(pmax(X, 0), 1)

  list(points = X, centers = centers)
}


generate_spatial_data <- function(N = 1000, rho = 0.5, beta = 1, sigma2 = 1, seed = NULL,
                                  position_uniforme = TRUE, params_tirage_non_unif = NULL) {
  # Packages nécessaires
  if (!requireNamespace("spdep", quietly = TRUE)) {
    stop("Le package 'spdep' est nécessaire. Veuillez l'installer avec install.packages('spdep').")
  }

  # Fixer la graine pour la reproductibilité si précisée
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Étape 1 : Génération des coordonnées spatiales
  if(position_uniforme){
    latitude <- runif(N, 0, 1)
    longitude <- runif(N, 0, 1)
    coords <- cbind(latitude, longitude)
  } else {
    if(is.null(params_tirage_non_unif)){
      stop("`params_tirage_non_unif ne peut être NULL si position_uniforme est à TRUE.")
    }
    coords <- generate_clusters_var(n_points = N, n_clusters = params_tirage_non_unif$n_clusters)
    coords <- coords$points
    colnames(coords) <- c("latitude", "longitude")
    latitude <- coords[, 1]
    longitude <- coords[, 2]
  }

  # Étape 2 : Génération de la variable explicative X
  X <- rnorm(N)

  # Étape 3 : Calcul de la matrice de voisinage W
  dist_matrix <- as.matrix(dist(coords))

  # Distance minimale pour qu'un point ait au moins un voisin
  threshold <- max(apply(dist_matrix, 1, function(x) min(x[x > 0])))

  # Matrice binaire de voisinage (1 si voisin, 0 sinon)
  W <- (dist_matrix <= threshold) * 1
  diag(W) <- 0  # Pas de boucle sur soi-même

  # Normalisation des lignes pour que chaque ligne somme à 1
  W <- W / rowSums(W)

  # Étape 4 : Génération de la variable d'intérêt Y selon un modèle SAR
  I <- diag(N)
  inverse_matrix <- solve(I - rho * W)

  # Bruit gaussien avec variance contrôlée
  epsilon <- rnorm(N, mean = 0, sd = sqrt(sigma2))

  # Calcul de Y
  Y <- inverse_matrix %*% (X * beta + epsilon)

  # Transformation pour assurer que Y est positive
  Y <- Y - min(Y)

  # Retourner les données sous forme de DataFrame et la matrice de voisinage
  data <- data.frame(latitude, longitude, X, Y = as.numeric(Y))
  return(list("data" = data, "adjacence" = inverse_matrix))
}


#Génération en utilisant un processus gaussien
generate_spatial_data_gp <- function(n, sigma2 = 1, phi = 1){
  sim1 <- grf(n, grid = "irreg",
              cov.model = "matern", cov.pars = c(sigma2, phi))
  sim <- cbind(sim1$coord, sim1$data)
  sim <- as.data.frame(sim)
  colnames(sim) <- c("longitude","latitude","X")
  sim$Y <- 2 + 2*(sim$X)^2 + rnorm(nrow(sim),10)
  compute_moran_index(sim, "Y", c("longitude","latitude"))
  plot_spatial_distribution2(sim)
  return(sim)
}



convert_to_sf_projected <- function(data, lon_col, lat_col, crs_projected) {
  # Vérification des prérequis
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Le package 'sf' est nécessaire. Veuillez l'installer avec install.packages('sf').")
  }

  if (!all(c(lon_col, lat_col) %in% names(data))) {
    stop("Les colonnes spécifiées pour longitude et latitude n'existent pas dans le DataFrame.")
  }

  # Création de l'objet sf en coordonnées géographiques (EPSG:4326)
  sf_object <- sf::st_as_sf(
    data,
    coords = c(lon_col, lat_col),
    crs = crs_projected
  )
  # Transformation en système de coordonnées projetées
  sf_projected <- sf::st_transform(sf_object, crs = crs_projected)

  return(sf_projected)
}


#Tirage en utilisant la méthode GRTS.
tirage_grts <- function(data_coord_pop, n,
                        lon_name = "longitude", lat_name = "latitude", proj = 4326){
  data_coord_pop$id <- 1:nrow(data_coord_pop)
  population_sf <- convert_to_sf_projected(data_coord_pop, lon_name, lat_name, proj)
  sample <- grts(
    sframe = population_sf,
    n_base = n, projcrs_check = FALSE
  )

  sf_sample <- sample$sites_base$id

  return(sf_sample)
}


#Affichage des données
plot_spatial_data <- function(data, color_var = "Y", point_size = 3) {
  # Vérifier que les colonnes nécessaires sont présentes
  if (!all(c("latitude", "longitude", color_var) %in% colnames(data))) {
    stop("Les colonnes 'latitude', 'longitude', et la variable spécifiée par 'color_var' doivent exister dans le DataFrame.")
  }

  # Charger ggplot2 si nécessaire
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Le package 'ggplot2' est nécessaire. Veuillez l'installer avec install.packages('ggplot2').")
  }

  # Chargement du package ggplot2
  library(ggplot2)

  # Création du graphique
  p <- ggplot(data, aes(x = latitude, y = longitude, color = !!sym(color_var))) +
    geom_point(size = point_size) +
    scale_color_viridis_c(option = "plasma") +  # Palette de couleurs agréable
    labs(title = "Scatterplot des données spatiales",
         x = "Latitude",
         y = "Longitude",
         color = paste("Valeurs de", color_var)) +
    theme_minimal() +
    coord_fixed()  # Maintient l'échelle entre x et y

  # Affichage du graphique
  print(p)
}



plot_spatial_distribution <- function(data) {
  # Vérification des colonnes nécessaires
  if (!all(c("latitude", "longitude", "Y") %in% colnames(data))) {
    stop("Le DataFrame doit contenir les colonnes 'latitude', 'longitude', et 'Y'.")
  }

  # Création du graphique
  p <- ggplot(data, aes(x = longitude, y = latitude, color = Y)) +
    geom_point(size = 3, alpha = 0.7) +  # Points avec transparence
    scale_color_gradient(low = "lightpink", high = "red") +  # Gradient de couleur
    labs(
      x = "longitude",
      y = "latitude",
      color = "Y",
      title = "Répartition de la variable Y"
    ) +
    theme_minimal() +  # Thème minimaliste
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Centrage du titre
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 10)
    )

  # Afficher le graphique
  print(p)
}

compute_moran_index <- function(data, var_name, coords) {
  # Vérifications des colonnes nécessaires
  if (!all(c(coords[1], coords[2], var_name) %in% colnames(data))) {
    stop(paste(
      "Le DataFrame doit contenir les colonnes suivantes :",
      paste(c(coords[1], coords[2], var_name), collapse = ", ")
    ))
  }

  # Extraction des coordonnées et de la variable
  coords_matrix <- as.matrix(data[, coords])
  variable <- data[[var_name]]

  # Étape 1 : Création de la matrice de voisinage (k plus proches voisins)
  knn <- knearneigh(coords_matrix, k = 1)  # Garantir au moins 1 voisin
  nb <- knn2nb(knn)
  listw <- nb2listw(nb, style = "W")  # Conversion en objet "listw"

  # Étape 2 : Calcul de l'indice de Moran
  moran <- moran.test(variable, listw)

  # Retourner les résultats
  return(list(
    moran_index = moran$estimate["Moran I statistic"],
    p_value = moran$p.value,
    z_score = moran$statistic
  ))
}

assign_named_zones <- function(data, coords = c("latitude", "longitude"), d = 10) {
  # Vérifications des colonnes nécessaires
  if (!all(coords %in% colnames(data))) {
    stop(paste(
      "Le DataFrame doit contenir les colonnes suivantes :",
      paste(coords, collapse = ", ")
    ))
  }

  # Calcul des indices de zones
  data$zone <- with(data, {
    x_zone <- ceiling(get(coords[2]) * d)  # Indice colonne (longitude)
    y_zone <- ceiling(get(coords[1]) * d)  # Indice ligne (latitude)
    paste(y_zone, x_zone, sep = "_")  # Nom de la zone sous forme "i_j"
  })

  # Retourne les données avec la colonne "zone"
  return(data)
}

modif_couts <- function(liste_params, n_ind){
  C <- liste_params$params_dist$C
  liste_params$params_dist$C <- C[n_ind,][, n_ind]
  return(liste_params)
}


echantillon_et_preparation <- function(donnees,
                                       n = 15,
                                       var_calage = c("copper"),
                                       liste_params = list(),
                                       spatialement = FALSE){
  #Tirage de l'échantillon à partir du data.frame
  #`donnees`
  if(!spatialement){
    n_ech <- sample(1:nrow(donnees), size = n, replace = FALSE)
  } else {
    n_ech <- tirage_grts(donnees, n)
  }

  #Calcul des poids
  donnees$poids <- nrow(donnees)/n #SRS
  terme_normalisation <- sum(donnees[n_ech,]$poids)
  donnees$proba <- donnees$poids/terme_normalisation


  #Calcul des marges
  marges <- colSums(donnees[, var_calage, drop = FALSE])
  marges <- as.matrix(marges, ncol = 1)

  #Calcul de la table des variables de calage
  var_x <- (donnees[n_ech, ])[, match(var_calage, colnames(donnees)), drop = FALSE]


  #Normalisation des marges
  marges <- marges/terme_normalisation

  #Ajout de la constante
  var_x <- cbind(var_x, "cst" = 1)
  marges <- rbind(marges, 1)
  d0_prob <- matrix(donnees[n_ech, "proba"], ncol = 1)

  res <- list(
    "d" = d0_prob,
    "X" = as.matrix(var_x),
    "M" = marges,
    "indice" = n_ech
  )

  res <- lapply(X = liste_params,
                FUN = function(element){element <- modif_couts(element, n_ech);
                return(c(element, res))})
  return(res)
}


application_calage <- function(liste_param){
  res <- lapply(X = liste_param, FUN = function(param){
    param <- param[names(param) %in% names(formals(calibration))];
    do.call(what = calibration, args = param)})
  ind <- unique(lapply(X = liste_param, `[[`, "indice"))
  if(length(ind) != 1){
    stop("Il y a un problème : différents indices de tirage dans la liste liste_param.")
  } else {
    ind <- ind[[1]]
  }
  return(list("ponderation" = res, "ind" = ind))
}

calcul_estimation_y <- function(res_calage, donnees_pop, donnees_calage, var_interet = c("cadmium","soil")){
  if(any(duplicated(var_interet))){
    stop("Il y a des doublons dans le vecteur renseigné dans var_interet.")
  }
  ind_ech <- res_calage$ind
  poids_cales <- Reduce(cbind, lapply(X = res_calage[[1]], FUN = `[[`, 1))
  poids_initiaux <- donnees_calage[[1]]$d
  poids_cales <- cbind(poids_cales, poids_initiaux)
  var_calage <- as.matrix(donnees_pop[ind_ech, match(var_interet, colnames(donnees_pop)), drop = FALSE])
  result_matrix <- lapply(1:ncol(var_calage), function(j) poids_cales * var_calage[, j])
  estim_y <- lapply(X = result_matrix, colSums)
  estim_y <- as.data.frame(matrix(Reduce(f = rbind, estim_y), nrow = length(var_interet)))
  estim_y <- cbind(var_interet, estim_y)
  colnames(estim_y) <- c("var_interet", paste0("methode_", 1:(ncol(poids_cales) - 1)), "HT")
  return(estim_y)
}

simulation_one_step <- function(donnees, n, var_calage, var_interet, liste_params, spatialement = FALSE){
  donnees_calage <- echantillon_et_preparation(donnees = donnees,
                                               n = n,
                                               var_calage = var_calage,
                                               liste_params = liste_params,
                                               spatialement = spatialement)
  pond_cales <- application_calage(donnees_calage)
  res_y <- calcul_estimation_y(pond_cales, donnees, donnees_calage, var_interet)
  #Ajout HT estimateur
  return(res_y)
}

nb_to_adj_matrix <- function(nb_object) {
  # Nombre de nœuds
  n <- length(nb_object)

  # Initialiser une matrice n x n avec des zéros
  adj_matrix <- matrix(0, nrow = n, ncol = n)

  # Remplir la matrice avec 1 là où il y a une connexion
  for (i in seq_along(nb_object)) {
    neighbors <- nb_object[[i]]
    adj_matrix[i, neighbors] <- 1
  }

  diag(adj_matrix) <- 1
  return(adj_matrix)
}

plot_spatial_distribution2 <- function(data) {
  # Vérification des colonnes nécessaires
  if (!all(c("latitude", "longitude", "Y") %in% colnames(data))) {
    stop("Le DataFrame doit contenir les colonnes 'latitude', 'longitude', et 'Y'.")
  }

  # Création du graphique
  p <- ggplot(data, aes(x = longitude, y = latitude, color = Y, size = Y)) +
    geom_point(alpha = 0.7) +  # Points avec transparence
    scale_color_gradient(low = "lightpink", high = "red") +  # Gradient de couleur
    scale_size_continuous(range = c(1, 10)) +  # Mise à l'échelle des tailles des points
    labs(
      x = "longitude",
      y = "latitude",
      color = "Y",
      size = "Y (taille)",
      title = "Répartition de la variable Y"
    ) +
    theme_minimal() +  # Thème minimaliste
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Centrage du titre
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 10)
    )

  # Afficher le graphique
  print(p)
}


k_nearest_matrix <- function(distance_matrix, k = 3) {
  # Vérifier que la matrice est carrée
  if (!is.matrix(distance_matrix) || nrow(distance_matrix) != ncol(distance_matrix)) {
    stop("distance_matrix doit être une matrice carrée.")
  }

  # Vérifier que k est un entier valide
  if (!is.numeric(k) || k <= 0 || k > nrow(distance_matrix)) {
    stop("k doit être un entier positif et inférieur ou égal au nombre de lignes de la matrice.")
  }

  # Initialiser la matrice de résultat
  n <- nrow(distance_matrix)
  result_matrix <- matrix(1, n, n)

  # Remplir la matrice avec des 0 pour les k plus proches voisins
  for (i in 1:n) {
    # Trouver les indices des k plus petites distances (hors l'individu i lui-même)
    nearest_indices <- order(distance_matrix[i, ], decreasing = FALSE)[1:(k + 1)]
    nearest_indices <- setdiff(nearest_indices, i) # Exclure i de ses propres voisins
    nearest_indices <- head(nearest_indices, k)   # Garder uniquement les k plus proches

    # Marquer les voisins avec 0
    result_matrix[i, nearest_indices] <- 0
  }
  diag(result_matrix) <- 0
  return(result_matrix)
}



transformation_df_res <- function(res){
  id <- sapply(X = res, nrow)
  id <- rep(1:length(id), id)
  res <- Reduce(f = rbind, res)
  res$id <- id
  print(colSums(is.na(res)))
  return(res)
}

