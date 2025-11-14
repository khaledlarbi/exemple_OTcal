
#####################
#### APPLICATION ####
#####################
library("MASS")   # Pour mvrnorm
library("fields") # Pour rdist
library("spsurvey")
library("dplyr")
library("ggplot2")
library("spdep")
library("geoR")
library("OTcalibration")

totals <- function(donnees, var_interet){
  #We assumed that sum of initial weights on the sample =
  #size of population
  initial_tot <- colSums(donnees[, match(var_interet,colnames(donnees)), drop = FALSE])
  N <- nrow(donnees)
  return(initial_tot/N)
}



biais_relatif <- function(res, var_interet, relatif = FALSE){
  nb_methode <- ncol(res) - 2 #La première colonne est le nom de la variable
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

# Fonction pour générer une matrice de corrélation spatiale
spatial_correlation_matrix <- function(locations, range, nugget, corr_fn = "exponential") {
  dists <- rdist(locations)  # Distances euclidiennes
  if (corr_fn == "exponential") {
    corr <- exp(-dists / range)
  } else if (corr_fn == "matern") {
    nu <- 2.1
    corr <- fields::Matern(d = dists, range = range, nu = nu)
  } else {
    stop("Unsupported correlation function")
  }
  corr <- (1 - nugget) * corr
  diag(corr) <- 1
  return(corr)
}




simulated_data <- generate_spatial_data(N = 1000, rho = 0.95,
                                        beta = 2, sigma2 = 1,
                                        position_uniforme = FALSE,
                                        params_tirage_non_unif = list("n_clusters" = 20))
summary(lm(data = simulated_data$data, Y ~ X))
plot_spatial_distribution(simulated_data$data)
compute_moran_index(simulated_data$data, "Y", c("latitude","longitude"))


simulated_data <- assign_named_zones(data = simulated_data$data, d = 2)
simulated_data <- fastDummies::dummy_cols(simulated_data, "zone")


#Ajout croisement Y et zone
# Calcul des produits
zone_columns <- grep("^zone_", names(simulated_data), value = TRUE)  # Sélection des colonnes 'zone_...'
# Calcul des produits
for (col in zone_columns) {
  simulated_data[[paste0("Y_", col)]] <- simulated_data$Y * simulated_data[[col]]  # Création d'une nouvelle colonne
}

dist <- compute_graph_distance(simulated_data)$graph_distance
dist_k1 <- k_nearest_matrix(dist, 1)
dist_k4 <- k_nearest_matrix(dist, 4)
dist_k10 <- k_nearest_matrix(dist, 10)


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
       "params_dist" = list("C" = dist, "epsilon" = 0.35)),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = dist, "epsilon" = 0.25)),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = dist, "epsilon" = 0.2)),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = dist, "epsilon" = c(0.2, 0.15, 0.14, 0.13, 0.12, 0.11,0.1))),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = dist, "epsilon" = c(0.2, 0.15, 0.14, 0.13, 0.12, 0.11, 0.1, 0.09, 0.085,
                                                      0.080, 0.075,
                                                      0.07, 0.065, 0.06, 0.055, 0.05))),
  list("method" = "wass_ent",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE,
       "params_dist" = list("C" = dist, "epsilon" = c(0.5, 0.3, 0.2, rev(seq(0.01, 0.1, 0.01))))),
  list("method" = "linear",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE),
  list("method" = "raking",
       "regularisation_param" = 1e-6,
       "normalisation"= TRUE))


nb_sim <- 250
taille_ech <- 50
var_int <- colnames(simulated_data)[startsWith(colnames(simulated_data),"Y")]

var_calage <- c("X")

res <- lapply(1:nb_sim, FUN = function(i){
  print(i)
  simulation_one_step(donnees = simulated_data,
                      n = taille_ech,
                      var_calage = var_calage,
                      var_interet = var_int,
                      liste_params = liste_params,
                      spatialement = TRUE)
})


res_srs <- lapply(1:nb_sim, FUN = function(i){
  print(i)
  simulation_one_step(donnees = simulated_data,
                      n = taille_ech,
                      var_calage = var_calage,
                      var_interet = var_int,
                      liste_params = liste_params,
                      spatialement = TRUE)
})


res <- transformation_df_res(res)
res_srs <- transformation_df_res(res_srs)


tot <- totals(as.matrix(simulated_data[,var_int, drop = FALSE]),
              var_int)


taille_pop <- nrow(simulated_data)
f <- taille_ech/taille_pop

res_mse <- variance_MC(res, c("Y"))

res_mse <- res_mse %>%
  mutate(across(starts_with(c("methode_","HT")), ~ .x / HT))


res_mse_srs <- variance_MC(res_srs, c("Y"))

res_mse_srs <- res_mse_srs %>%
  mutate(across(starts_with(c("methode_","HT")), ~ .x / HT))








#Visualisation

#Tirage
i <- 1
liste_param <- list(list("method" = "wass_ent",
                    "regularisation_param" = 1e-6,
                    "normalisation"= TRUE,
                    "params_dist" = list("C" = dist,
                                         "epsilon" = rev(seq(0.35, 0.5, 0.05)))))
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

plot_optimal_transport(P = plan_opti$P, axis_digits = 1)

