
#' Load data for PSO clustering
#'
#' @param fibro_path Path to the fibroblast RDS file.
#' @param deg_path Path to the DEG Excel file.
#' @return A list containing the fibroblast data and DEG list.
#' @export
load_data <- function(fibro_path, deg_path) {
  fibro <- readRDS(fibro_path)
  deg <- read_xlsx(deg_path, col_names = TRUE)
  deg <- as.data.frame(deg)
  list(fibro = fibro, deg = deg)
}

#' Initialize PSO parameters
#'
#' @param deg The DEG dataframe.
#' @param size Number of particles.
#' @param generation Number of generations.
#' @param toprange The top range of genes to consider.
#' @return A list of initialized PSO parameters.
#' @export
initialize_pso <- function(deg, size = 30, generation = 3, toprange = 15) {
  deg_list <- list()
  for (i in 1:ncol(deg)) {
    deg_list[[i]] <- as.character(na.omit(deg[,i]))
  }
  n <- length(deg_list)  # Number of gene sets (N dimensions)
  list(
    size = size,
    generation = generation,
    toprange = toprange,
    vmax = 1,
    c1 = 1,
    c2 = 1,
    r1 = runif(1),
    r2 = runif(1),
    deg_list = deg_list,
    n = n
  )
}

#' Run PSO algorithm to optimize DEG list
#'
#' @param deg_list A list of DEG (differentially expressed gene) sets.
#' @param size Number of particles.
#' @param generation Number of generations.
#' @param vmax Maximum velocity for particles.
#' @param c1 Cognitive coefficient.
#' @param c2 Social coefficient.
#' @param toprange The top range of genes to consider.
#' @return The best solution found by PSO.
#' @export
run_pso <- function(deg_list, size = 30, generation = 3, vmax = 1, c1 = 1, c2 = 1, toprange = 15) {

  # Initialization
  n <- length(deg_list)  # Number of dimensions (DEG sets)
  position <- matrix(runif(n * size), nrow = size)  # Initial particle positions
  velocity <- matrix(runif(n * size, -vmax, vmax), nrow = size)  # Initial velocities
  personal_best <- position  # Personal best positions
  global_best <- position[1, ]  # Initial global best position

  # Objective function (e.g., silhouette score or other clustering metric)
  objective_function <- function(pos) {
    sum(pos)  # Simple placeholder example
  }

  personal_best_scores <- apply(personal_best, 1, objective_function)
  global_best_score <- max(personal_best_scores)

  # PSO iteration
  for (g in 1:generation) {
    for (i in 1:size) {
      r1 <- runif(1)
      r2 <- runif(1)
      velocity[i, ] <- velocity[i, ] + c1 * r1 * (personal_best[i, ] - position[i, ]) +
                       c2 * r2 * (global_best - position[i, ])
      position[i, ] <- position[i, ] + velocity[i, ]
      score <- objective_function(position[i, ])
      if (score > personal_best_scores[i]) {
        personal_best[i, ] <- position[i, ]
        personal_best_scores[i] <- score
      }
      if (score > global_best_score) {
        global_best <- position[i, ]
        global_best_score <- score
      }
    }
  }
  return(global_best)  # Return the best solution found
}

#' Perform GMM clustering
#'
#' @param data The input data matrix (e.g., gene expression matrix).
#' @param n_clusters The number of clusters to fit.
#' @return A list containing the GMM model and clustering labels.
#' @export
run_gmm <- function(data, n_clusters = 3) {
  gmm_model <- Mclust(data, G = n_clusters)
  clusters <- gmm_model$classification  # Cluster labels
  list(model = gmm_model, clusters = clusters)
}
