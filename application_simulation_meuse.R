
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

plot_meuse_variable(simulated_data, "soil")
plot_meuse_variable(simulated_data, "copper")
plot_meuse_variable(simulated_data, "om")
plot_meuse_variable(simulated_data, "cadmium")
plot_meuse_variable(simulated_data, "lead")



#plot_meuse_variable(simulated_data, "zinc")
summary(lm(data = simulated_data, cadmium ~ elev  - 1))
summary(lm(data = simulated_data, lead ~ elev - 1))
summary(lm(data = simulated_data, copper ~ elev - 1))
summary(lm(data = simulated_data, om ~ elev - 1))
summary(lm(data = simulated_data, soil ~ elev - 1))

#j'ai l'impression que ça marche mieux quabd mauvais lien
summary(lm(data = simulated_data, cadmium ~  dist.m ))
summary(lm(data = simulated_data, lead ~ dist.m ))
summary(lm(data = simulated_data, copper ~ dist.m ))
summary(lm(data = simulated_data, om ~ dist.m ))
summary(lm(data = simulated_data, soil ~ dist.m))


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
var_int_init <- c("cadmium", "lead", "copper","om","soil")

simulated_data <- add_geographical_zones(simulated_data)
simulated_data <- one_hot_encode_zone(simulated_data)
simulated_data$zone <- as.factor(simulated_data$zone)

summary(lm(data = simulated_data, cadmium ~  zone - 1))
summary(lm(data = simulated_data, lead ~ zone - 1))
summary(lm(data = simulated_data, copper ~ zone - 1))
summary(lm(data = simulated_data, om ~ zone - 1))
summary(lm(data = simulated_data, soil ~ zone - 1))



zone_columns <- grep("^zone_", names(simulated_data), value = TRUE)  # Sélection des colonnes 'zone_...'
# Calcul des produits
for (col in zone_columns) {
  for(int in var_int_init){
    simulated_data[[paste0(col, int, collapse = "_")]] <- simulated_data[[col]] * simulated_data[[int]]
  }
}

var_int <- c(var_int_init, as.vector(outer(zone_columns, var_int_init, paste, sep = "")))

dist_k1 <- k_nearest_matrix(d_prime, 1)
dist_k4 <- k_nearest_matrix(d_prime, 4)
dist_k10 <- k_nearest_matrix(d_prime, 10)


liste_params <- list(
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = dist_k1, "epsilon" = 0.3)),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = dist_k4, "epsilon" = 0.3)),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = dist_k4, "epsilon" = 0.5)),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = dist_k10, "epsilon" = rev(seq(0.1, 0.5,0.05)))),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_prime, "epsilon" = 2)),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_prime, "epsilon" = 1.5)),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_prime, "epsilon" = 1)),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_prime, "epsilon" = 0.75)),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_prime, "epsilon" = 0.5)),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_prime, "epsilon" = 0.4)),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_prime, "epsilon" = 0.35)),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_prime, "epsilon" = 0.3)),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_prime, "epsilon" = 0.25)),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_prime, "epsilon" = 0.2)),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_prime, "epsilon" = 0.15)),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = d_prime, "epsilon" = c(0.2, 0.15, 0.14, 0.13, 0.12, 0.11,0.1))),
  list("method" = "linear",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE),
  list("method" = "raking",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE))

#var_calage <- c("dist.m")
simulated_data$const <- 1

#var_calage <- c("const")
#var_calage <- c("zone_1", "zone_2", "zone_3", "zone_4")
var_calage %in% colnames(simulated_data)

var_calage <- c("const", "dist.m")

res <- lapply(1:nb_sim, FUN = function(i){
  print(i)
  simulation_one_step(simulated_data,
                      taille_ech,
                      var_calage,
                      var_int,
                      liste_params,
                      spatialement = FALSE)
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


#OM Pas efficace

moran_I <- function(y, d_mat, method = "inverse", threshold = NULL) {
  n <- length(y)

  if (method == "inverse") {
    W <- 1 / d_mat
    diag(W) <- 0
    if (!is.null(threshold)) {
      W[d_mat > threshold] <- 0
    }
  }

  # Normalisation par ligne
  W <- W / rowSums(W)

  y_bar <- mean(y)
  y_dev <- y - y_bar
  S0 <- sum(W)

  numerator <- sum(W * (outer(y_dev, y_dev)))
  denominator <- sum(y_dev^2)

  I <- (n / S0) * numerator / denominator
  return(I)
}


weight_mat <- 1 / d_mat
diag(weight_mat) <- 0  # On ne veut pas se pondérer soi-même
nb <- spdep::mat2listw(weight_mat, style = "W", zero.policy = TRUE)
spdep::moran.test(simulated_data$cadmium, nb)
moran_I(simulated_data$cadmium, d_mat)
spdep::moran.test(simulated_data$copper, nb)
spdep::moran.test(simulated_data$lead, nb)
spdep::moran.test(simulated_data$om, nb)
spdep::moran.test(simulated_data$soil, nb)


#Tirage
i <- 1
liste_param <- list(list("method" = "wass_ent",
                         "regularisation_param" = 1e-6,
                         "normalisation"= TRUE,
                         "params_dist" = list("C" = d_mat,
                                              "epsilon" = rev(seq(1, 10, 0.1)))))
prep_donnees_calage <- echantillon_et_preparation(donnees = simulated_data,
                                                  n = taille_ech,
                                                  var_calage = var_calage,
                                                  liste_params = liste_param,
                                                  spatialement = FALSE)
pond_cales <- application_calage(prep_donnees_calage)


plan_opti <- sinkhorn(a = prep_donnees_calage[[i]]$d,
                      b = pond_cales$ponderation[[i]]$weights,
                      C = prep_donnees_calage[[i]]$params_dist$C,
                      epsilon = min(prep_donnees_calage[[i]]$params_dist$epsilon))
plan_opti$convergence
plot_optimal_transport(P = plan_opti$P, axis_digits = 1)

#Fonctionne mieux sur cadiaum et lead

