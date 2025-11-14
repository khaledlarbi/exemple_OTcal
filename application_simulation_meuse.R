
#Ajout des matrices de distance
library("igraph")
set.seed(1995)
library("gstat")
library("dplyr")
devtools::load_all()
library(spsurvey)
library(ggplot2)

data("meuse.all")


simulated_data <- meuse.all
simulated_data[is.na(simulated_data$om), "om"] <- median(simulated_data$om,na.rm = TRUE)
simulated_data$log_zinc <- log(simulated_data$zinc)
colnames(simulated_data)[2:3] <- c("longitude","latitude")

d_mat  <- as.matrix(dist(simulated_data[, c("longitude","latitude")]))       # matrice des distances euclidiennes
d_mat <- d_mat/max(d_mat)

d_prime <- compute_graph_distance(simulated_data)$graph_distance


plot_meuse_variable <- function(data, varname) {
  if (!all(c("longitude", "latitude") %in% names(data))) {
    stop("Le jeu de données doit contenir les colonnes 'longitude' et 'latitude'.")
  }

  if (!varname %in% names(data)) {
    stop(paste("La variable", varname, "n'existe pas dans le jeu de données."))
  }

  ggplot(data, aes_string(x = "longitude", y = "latitude", color = varname)) +
    geom_point(size = 3) +
    scale_color_viridis_c(option = "plasma", name = varname) +
    coord_equal() +
    labs(title = paste("Carte de la variable", varname),
         x = "Longitude",
         y = "Latitude") +
    theme_minimal()
}


plot_meuse_variable(simulated_data, "cadmium")
plot_meuse_variable(simulated_data, "log_zinc")
plot_meuse_variable(simulated_data, "soil")


add_geographical_zones <- function(data, k = 4) {
  if (!all(c("longitude", "latitude") %in% names(data))) {
    stop("Le jeu de données doit contenir les colonnes 'longitude' et 'latitude'.")
  }

  coords <- data[, c("longitude", "latitude")]
  kmeans_result <- kmeans(coords, centers = k, nstart = 25)

  data$zone <- (kmeans_result$cluster)
  return(data)
}

one_hot_encode_zone <- function(data, zone_col = "zone") {
  if (!zone_col %in% names(data)) {
    stop(paste("La colonne", zone_col, "n'existe pas dans le jeu de données."))
  }

  # S'assurer que la colonne est un facteur
  data[[zone_col]] <- factor(data[[zone_col]])

  # Créer le one-hot encoding
  dummies <- model.matrix(~ . - 1, data = data[zone_col])

  # Renommer les colonnes proprement
  colnames(dummies) <- gsub(paste0("^", zone_col), "", colnames(dummies))
  colnames(dummies) <- paste0(zone_col, "_", colnames(dummies))

  # Ajouter au dataset original
  data <- cbind(data, dummies)
  data[[zone_col]] <- as.numeric(data[[zone_col]])
  return(data)
}


nb_sim <- 1000
taille_ech <- 30
var_int_init <- c("cadmium", "soil", "log_zinc")
var_calage <- c("elev", "ffreq", "om")
var_calage %in% colnames(simulated_data)

simulated_data <- add_geographical_zones(simulated_data)
simulated_data <- one_hot_encode_zone(simulated_data)




zone_columns <- grep("^zone_", names(simulated_data), value = TRUE)  # Sélection des colonnes 'zone_...'
# Calcul des produits
for (col in zone_columns) {
  for(int in var_int_init){
    simulated_data[[paste0(col, int, collapse = "_")]] <- simulated_data[[col]] * simulated_data[[int]]
  }
}

var_int <- c(var_int_init, as.vector(outer(zone_columns, var_int_init, paste, sep = "")))




liste_params <- list(
  list("method" = "wass_ent",
       "regularisation_param" = 1e-4,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_mat, "epsilon" = c(0.1))),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-4,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_prime, "epsilon" = c(0.1))),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-4,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_mat, "epsilon" = c(0.3))),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-4,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_prime, "epsilon" = c(0.3))),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-4,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_prime, "epsilon" = c(0.2))),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-4,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_prime, "epsilon" = c(0.4))),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-4,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_mat, "epsilon" = c(0.5))),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-4,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_prime, "epsilon" = c(0.5))),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-4,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_mat, "epsilon" = c(0.8))),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-4,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_prime, "epsilon" = c(0.8))),
  list("method" = "linear",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE),
  list("method" = "raking",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE)
)


res <- lapply(1:nb_sim, FUN = function(i){
  print(i)
  simulation_one_step(simulated_data,
                      taille_ech,
                      var_calage,
                      var_int,
                      liste_params,
                      spatialement = TRUE)
})



res <- transformation_df_res(res)

totals <- function(donnees, var_interet){
  #We assumed that sum of initial weights on the sample =
  #size of population
  initial_tot <- colSums(donnees[, match(var_interet,colnames(donnees)), drop = FALSE])
  N <- nrow(donnees)
  return(initial_tot/N)
}

biais_relatif <- function(res, var_interet, relatif = FALSE){
  nb_methode <- sum(startsWith(colnames(res), "methode"))
  resultat_par_methode <- res[, paste0("methode_", 1:nb_methode)]
  tot_par_methode <- matrix(rep(tot[res$var_interet], nb_methode), ncol = nb_methode)
  resultat_par_methode <- resultat_par_methode - tot_par_methode
  if(relatif){
    resultat_par_methode <- resultat_par_methode/tot_par_methode
  }

  res_ag <- aggregate.data.frame(x = resultat_par_methode,
                                 by = list(res$var_interet),
                                 FUN = "mean",
                                 na.rm = TRUE)
  return(res_ag)
}

variance_MC <- function(res, var_interet, relatif = FALSE){
  nb_methode <- ncol(res) - 3 #La première colonne est le nom de la variable
  resultat_par_methode <- res[, c(paste0("methode_", 1:nb_methode), "HT")]
  tot_par_methode <- matrix(rep(tot[res$var_interet], nb_methode), ncol = nb_methode)
  resultat_par_methode <- (resultat_par_methode - tot_par_methode)^2

  res_ag <- aggregate.data.frame(x = resultat_par_methode,
                                 by = list(res$var_interet),
                                 FUN = "mean", na.rm = TRUE)
  return(res_ag)
}

tot <- totals(as.matrix(simulated_data[,var_int, drop = FALSE]),
              var_int)


taille_pop <- nrow(simulated_data)
f <- taille_ech/taille_pop

res_mse <- variance_MC(res,
                       c("Y"))



res_mse <- res_mse %>%
  mutate(across(starts_with(c("methode_","HT")), ~ .x / HT))






weight_mat <- 1 / d_mat
diag(weight_mat) <- 0  # On ne veut pas se pondérer soi-même
nb <- spdep::mat2listw(weight_mat, style = "W", zero.policy = TRUE)
spdep::moran.test(simulated_data$log_zinc, nb)
spdep::moran.test(simulated_data$cadmium, nb)
spdep::moran.test(simulated_data$soil, nb)

