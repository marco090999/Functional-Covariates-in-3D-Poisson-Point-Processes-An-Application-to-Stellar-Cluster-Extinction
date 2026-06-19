# ============================================================================
# Functional spatial covariates in three-dimensional Poisson point processes
# Stellar-cluster application code
# ============================================================================
#
# This script reproduces the application to stellar clusters. It constructs
# cluster-level extinction quantile curves, computes FPCA scores, fits the
# functional and scalar-summary Poisson point process models, runs the 3D
# inhomogeneous K-function diagnostic, and performs the parametric bootstrap
# for the reconstructed coefficient function.
#
# Before running this script, source `01_simulation_study.R`. The input file
# `data/processed/stars_with_extinction.rds` must contain a data frame with
# ID_CL_PAPER, X, Y, Z, and A0_ML. An optional cluster metadata file can contain
# ID_CL_PAPER and group_final.
# ============================================================================

if (!exists("project_root", mode = "function")) {
  project_root <- function() {
    root <- Sys.getenv("PROJECT_ROOT", unset = getwd())
    normalizePath(root, winslash = "/", mustWork = FALSE)
  }
}

application_script_path <- function() {
  candidate <- Sys.getenv("APPLICATION_SCRIPT_PATH", unset = "")
  if (nzchar(candidate)) return(normalizePath(candidate, winslash = "/", mustWork = FALSE))
  file.path(project_root(), "R", "02_stellar_cluster_application.R")
}

load_simulation_functions <- function(root = project_root()) {
  required <- c("FunPois.fit_v2", "fpca_grid", "s3dpm", "trapz_weights")
  missing <- required[!vapply(required, exists, logical(1), mode = "function")]
  if (length(missing)) {
    source(file.path(root, "R", "01_simulation_study.R"))
  }
  invisible(TRUE)
}

check_application_packages <- function() {
  required <- c("mgcv", "GET", "plot3D", "ggplot2", "dplyr", "tidyr", "future", "future.apply")
  unavailable <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(unavailable)) stop("Install required packages: ", paste(unavailable, collapse = ", "))
  invisible(required)
}

make_application_directories <- function(root = project_root()) {
  directories <- c(
    file.path(root, "results", "application"),
    file.path(root, "figures", "application"),
    file.path(root, "tables", "application")
  )
  invisible(lapply(directories, dir.create, recursive = TRUE, showWarnings = FALSE))
  list(results = directories[1], figures = directories[2], tables = directories[3])
}

read_application_inputs <- function(root = project_root(),
                                    stars_file = file.path(root, "data", "processed", "stars_with_extinction.rds"),
                                    metadata_file = file.path(root, "data", "processed", "cluster_metadata.rds")) {
  if (!file.exists(stars_file)) {
    stop("Missing stellar input file: ", stars_file,
         ". Save a data frame with ID_CL_PAPER, X, Y, Z, and A0_ML as an RDS file.")
  }
  stars <- readRDS(stars_file)
  required <- c("ID_CL_PAPER", "X", "Y", "Z", "A0_ML")
  absent <- setdiff(required, names(stars))
  if (length(absent)) stop("The stellar data are missing: ", paste(absent, collapse = ", "))
  metadata <- if (file.exists(metadata_file)) readRDS(metadata_file) else NULL
  list(stars = as.data.frame(stars), metadata = metadata)
}

# -----------------------------------------------------------------------------
# Cluster-level functional covariate construction
# -----------------------------------------------------------------------------

select_fpca_dimension <- function(curve_matrix, u_grid, minimum_components = 2L,
                                  maximum_components = 25L, pve_target = 0.95,
                                  plateau_threshold = 0.0005, plateau_runs = 1L) {
  stopifnot(is.matrix(curve_matrix), length(u_grid) == ncol(curve_matrix))
  maximum_components <- min(maximum_components, nrow(curve_matrix) - 1L, ncol(curve_matrix))
  minimum_components <- max(1L, min(minimum_components, maximum_components))
  weighted_fpca <- fpca_grid(curve_matrix, u_grid, K = maximum_components)
  cumulative_pve <- cumsum(weighted_fpca$d^2) / weighted_fpca$totalSS
  gains <- c(cumulative_pve[1], diff(cumulative_pve))
  selected <- maximum_components
  small_gain_run <- 0L
  for (k in seq_along(cumulative_pve)) {
    if (k < minimum_components) next
    if (cumulative_pve[k] >= pve_target) {
      selected <- k
      break
    }
    if (cumulative_pve[k] >= 0.90) {
      small_gain_run <- if (gains[k] < plateau_threshold) small_gain_run + 1L else 0L
      if (small_gain_run >= plateau_runs) {
        selected <- k
        break
      }
    }
  }
  list(K = selected, pve = cumulative_pve, gain = gains, fpca = weighted_fpca)
}

prepare_cluster_application_data <- function(stars, metadata = NULL,
                                             minimum_cluster_size = 30L,
                                             u_grid = seq(0, 1, length.out = 101),
                                             K = 4L) {
  stars <- dplyr::as_tibble(stars) |>
    dplyr::filter(!is.na(A0_ML)) |>
    dplyr::mutate(ID_CL_PAPER = as.character(ID_CL_PAPER))
  retained_ids <- stars |>
    dplyr::count(ID_CL_PAPER, name = "n_stars") |>
    dplyr::filter(n_stars >= minimum_cluster_size) |>
    dplyr::pull(ID_CL_PAPER)
  stars <- stars |> dplyr::filter(ID_CL_PAPER %in% retained_ids)

  clusters <- stars |>
    dplyr::group_by(ID_CL_PAPER) |>
    dplyr::summarise(
      n_stars = dplyr::n(),
      x = stats::median(X, na.rm = TRUE),
      y = stats::median(Y, na.rm = TRUE),
      z = stats::median(Z, na.rm = TRUE),
      .groups = "drop"
    )
  if (!is.null(metadata) && all(c("ID_CL_PAPER", "group_final") %in% names(metadata))) {
    clusters <- clusters |> dplyr::left_join(
      metadata |> dplyr::mutate(ID_CL_PAPER = as.character(ID_CL_PAPER)) |>
        dplyr::select(ID_CL_PAPER, group_final),
      by = "ID_CL_PAPER"
    )
  }

  curve_matrix <- vapply(clusters$ID_CL_PAPER, function(cluster_id) {
    values <- stars$A0_ML[stars$ID_CL_PAPER == cluster_id]
    stats::quantile(values, probs = u_grid, type = 8, names = FALSE, na.rm = TRUE)
  }, numeric(length(u_grid)))
  curve_matrix <- t(curve_matrix)
  rownames(curve_matrix) <- clusters$ID_CL_PAPER
  colnames(curve_matrix) <- paste0("u", seq_along(u_grid))

  fpca <- fpca_grid(curve_matrix, u_grid, K = K)
  scores <- fpca$scores
  colnames(scores) <- paste0("ext_pc", seq_len(K))
  cluster_data <- dplyr::bind_cols(clusters, as.data.frame(scores))
  functional_pattern <- s3dpm(cluster_data, names = c("n_stars", colnames(scores)))
  attr(functional_pattern, "cluster_id") <- clusters$ID_CL_PAPER

  weights <- fpca$w / sum(fpca$w)
  scalar_mean <- as.numeric(curve_matrix %*% weights)
  scalar_pattern <- functional_pattern
  scalar_pattern$df$ext_mean <- scalar_mean

  list(
    stars = stars,
    clusters = clusters,
    curve_matrix = curve_matrix,
    u_grid = u_grid,
    fpca = fpca,
    functional_pattern = functional_pattern,
    scalar_pattern = scalar_pattern,
    scalar_mean = scalar_mean,
    K = K
  )
}

idw_predict_multivariate <- function(new_coordinates, observed_coordinates, observed_values,
                                    power = 4, epsilon = 1e-8, chunk_size = 10000L) {
  new_coordinates <- as.matrix(new_coordinates)
  observed_coordinates <- as.matrix(observed_coordinates)
  observed_values <- as.matrix(observed_values)
  output <- matrix(NA_real_, nrow(new_coordinates), ncol(observed_values))
  colnames(output) <- colnames(observed_values)
  observed_norm <- rowSums(observed_coordinates^2)
  for (start in seq(1L, nrow(new_coordinates), by = chunk_size)) {
    end <- min(start + chunk_size - 1L, nrow(new_coordinates))
    current <- new_coordinates[start:end, , drop = FALSE]
    distance_squared <- outer(rowSums(current^2), observed_norm, "+") -
      2 * tcrossprod(current, observed_coordinates)
    distances <- sqrt(pmax(distance_squared, 0))
    weights <- (distances + epsilon)^(-power)
    predictions <- weights %*% observed_values / rowSums(weights)
    exact <- which(apply(distances, 1L, min) < epsilon * 10)
    if (length(exact)) {
      for (row in exact) predictions[row, ] <- observed_values[which.min(distances[row, ]), ]
    }
    output[start:end, ] <- predictions
  }
  output
}

make_idw_field_function <- function(pattern, value_columns, power = 4) {
  observed_coordinates <- as.matrix(pattern$df[, c("x", "y", "z")])
  observed_values <- as.matrix(pattern$df[, value_columns, drop = FALSE])
  storage.mode(observed_coordinates) <- "double"
  storage.mode(observed_values) <- "double"
  force(observed_coordinates); force(observed_values); force(value_columns); force(power)
  function(coordinates) {
    coordinates <- as.matrix(coordinates)
    colnames(coordinates) <- c("x", "y", "z")
    output <- idw_predict_multivariate(coordinates, observed_coordinates, observed_values, power = power)
    colnames(output) <- value_columns
    output
  }
}

build_application_formula <- function(covariate_columns, xy_basis = 30L, z_basis = 10L) {
  rhs <- c(sprintf("s(x, y, k = %d)", xy_basis), sprintf("s(z, k = %d)", z_basis), covariate_columns)
  stats::as.formula(paste("~", paste(rhs, collapse = " + ")))
}

fit_application_models <- function(application_data, idw_power = 4L,
                                   mult = 30L, ncube = 20L,
                                   xy_basis_functional = 30L, xy_basis_scalar = 20L,
                                   z_basis = 10L, seed = 123L) {
  functional_columns <- grep("^ext_pc[0-9]+$", names(application_data$functional_pattern$df), value = TRUE)
  functional_field <- make_idw_field_function(application_data$functional_pattern, functional_columns, idw_power)
  scalar_field <- make_idw_field_function(application_data$scalar_pattern, "ext_mean", idw_power)

  functional_fit <- FunPois.fit_v2(
    X = application_data$functional_pattern,
    formula = build_application_formula(functional_columns, xy_basis_functional, z_basis),
    covs = NULL, marked = FALSE, spatial.cov = FALSE, verbose = TRUE,
    mult = mult, ncube = ncube, seed = seed, grid = FALSE,
    mark.c = TRUE, mark.c.mode = "oracle", mark.oracle_fun = functional_field,
    interp.p = idw_power, process.type = "s3d", z_weight = FALSE,
    sphere = FALSE
  )
  scalar_fit <- FunPois.fit_v2(
    X = application_data$scalar_pattern,
    formula = build_application_formula("ext_mean", xy_basis_scalar, z_basis),
    covs = NULL, marked = FALSE, spatial.cov = FALSE, verbose = TRUE,
    mult = mult, ncube = ncube, seed = seed, grid = FALSE,
    mark.c = TRUE, mark.c.mode = "oracle", mark.oracle_fun = scalar_field,
    interp.p = idw_power, process.type = "s3d", z_weight = FALSE,
    sphere = FALSE
  )
  list(functional = functional_fit, scalar = scalar_fit,
       functional_field = functional_field, scalar_field = scalar_field,
       functional_columns = functional_columns)
}

extract_score_coefficients <- function(fit, prefix = "^ext_pc[0-9]+$") {
  coefficients <- stats::coef(fit$mod_global)
  names_found <- grep(prefix, names(coefficients), value = TRUE)
  if (!length(names_found)) stop("No FPCA-score coefficients were found in the fitted model.")
  names_found <- names_found[order(as.integer(sub("ext_pc", "", names_found)))]
  estimates <- as.numeric(coefficients[names_found])
  names(estimates) <- names_found
  estimates
}

reconstruct_functional_coefficient <- function(score_coefficients, fpca, u_grid, center = TRUE) {
  phi <- fpca$phi[, seq_along(score_coefficients), drop = FALSE]
  theta <- as.numeric(phi %*% score_coefficients)
  if (center) theta <- theta - sum(theta * fpca$w) / sum(fpca$w)
  data.frame(u = u_grid, theta = theta)
}

make_model_summary_table <- function(functional_fit, scalar_fit, functional_columns) {
  scalar_coefficients <- summary(scalar_fit$mod_global)$p.table
  functional_coefficients <- summary(functional_fit$mod_global)$p.table
  collect <- function(tab, terms) {
    out <- tab[intersect(terms, rownames(tab)), , drop = FALSE]
    data.frame(
      term = rownames(out), estimate = out[, 1], standard_error = out[, 2],
      statistic = out[, 3], p_value = out[, 4], row.names = NULL
    )
  }
  list(
    functional_parametric = collect(functional_coefficients, c("(Intercept)", functional_columns)),
    scalar_parametric = collect(scalar_coefficients, c("(Intercept)", "ext_mean")),
    functional_smooth = as.data.frame(summary(functional_fit$mod_global)$s.table),
    scalar_smooth = as.data.frame(summary(scalar_fit$mod_global)$s.table)
  )
}

# -----------------------------------------------------------------------------
# Three-dimensional inhomogeneous K-function diagnostic
# -----------------------------------------------------------------------------

runif_box3 <- function(n, box) {
  cbind(
    x = stats::runif(n, box$xrange[1], box$xrange[2]),
    y = stats::runif(n, box$yrange[1], box$yrange[2]),
    z = stats::runif(n, box$zrange[1], box$zrange[2])
  )
}

Kinhom3D_box <- function(coordinates, intensity, r, box, correction = TRUE,
                         intensity_floor = NULL) {
  coordinates <- as.matrix(coordinates)
  intensity <- as.numeric(intensity)
  if (!is.null(intensity_floor)) intensity <- pmax(intensity, intensity_floor)
  n <- nrow(coordinates)
  volume <- diff(box$xrange) * diff(box$yrange) * diff(box$zrange)
  if (n < 2L) return(rep(0, length(r)))
  pairs <- utils::combn(n, 2L)
  dx <- abs(coordinates[pairs[1, ], 1] - coordinates[pairs[2, ], 1])
  dy <- abs(coordinates[pairs[1, ], 2] - coordinates[pairs[2, ], 2])
  dz <- abs(coordinates[pairs[1, ], 3] - coordinates[pairs[2, ], 3])
  distance <- sqrt(dx^2 + dy^2 + dz^2)
  intensity_product <- intensity[pairs[1, ]] * intensity[pairs[2, ]]
  keep <- is.finite(distance) & is.finite(intensity_product) & intensity_product > 0
  dx <- dx[keep]; dy <- dy[keep]; dz <- dz[keep]
  distance <- distance[keep]; intensity_product <- intensity_product[keep]
  if (!length(distance)) return(rep(0, length(r)))
  edge_weight <- if (correction) volume / ((diff(box$xrange) - dx) * (diff(box$yrange) - dy) * (diff(box$zrange) - dz)) else rep(1, length(distance))
  contribution <- 2 * edge_weight / intensity_product
  order_index <- order(distance)
  cumulative <- cumsum(contribution[order_index])
  last_index <- findInterval(r, distance[order_index])
  output <- numeric(length(r))
  output[last_index > 0] <- cumulative[last_index[last_index > 0]] / volume
  output
}

make_lambda_predictor <- function(fit, field_function, covariate_columns) {
  force(fit); force(field_function); force(covariate_columns)
  function(coordinates) {
    coordinates <- as.matrix(coordinates)
    colnames(coordinates) <- c("x", "y", "z")
    covariates <- field_function(coordinates)
    newdata <- data.frame(x = coordinates[, 1], y = coordinates[, 2], z = coordinates[, 3])
    newdata <- cbind(newdata, as.data.frame(covariates[, covariate_columns, drop = FALSE]))
    intensity <- exp(stats::predict(fit$mod_global, newdata = newdata, type = "link"))
    intensity[!is.finite(intensity)] <- 0
    intensity
  }
}

simulate_inhomogeneous_box <- function(lambda_function, box, n_mc = 100000L,
                                       safety_factor = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  volume <- diff(box$xrange) * diff(box$yrange) * diff(box$zrange)
  mc_points <- runif_box3(n_mc, box)
  lambda_max <- max(lambda_function(mc_points), na.rm = TRUE) * safety_factor
  if (!is.finite(lambda_max) || lambda_max <= 0) stop("The thinning bound is invalid.")
  candidates <- runif_box3(stats::rpois(1L, lambda_max * volume), box)
  if (!nrow(candidates)) return(list(points = candidates, lambda = numeric(0), lambda_max = lambda_max))
  lambda_candidates <- lambda_function(candidates)
  kept <- stats::runif(nrow(candidates)) <= pmin(lambda_candidates / lambda_max, 1)
  list(points = candidates[kept, , drop = FALSE], lambda = lambda_candidates[kept], lambda_max = lambda_max)
}

run_get_diagnostic <- function(fit, observed_pattern, field_function, covariate_columns,
                               B = 200L, n_mc = 100000L, n_r = 80L,
                               seed = 123L) {
  observed_coordinates <- as.matrix(observed_pattern$df[, c("x", "y", "z")])
  box <- list(
    xrange = range(observed_coordinates[, 1]),
    yrange = range(observed_coordinates[, 2]),
    zrange = range(observed_coordinates[, 3])
  )
  lambda_function <- make_lambda_predictor(fit, field_function, covariate_columns)
  lambda_observed <- lambda_function(observed_coordinates)
  r_max <- min(diff(box$xrange), diff(box$yrange), diff(box$zrange)) / 4
  r_grid <- seq(0, r_max, length.out = n_r + 1L)[-1L]
  observed_k <- Kinhom3D_box(observed_coordinates, lambda_observed, r_grid, box)
  simulations <- vector("list", B)
  K_sim <- matrix(NA_real_, B, length(r_grid))
  for (b in seq_len(B)) {
    simulated <- simulate_inhomogeneous_box(lambda_function, box, n_mc = n_mc, seed = seed + b)
    simulations[[b]] <- simulated
    K_sim[b, ] <- Kinhom3D_box(simulated$points, simulated$lambda, r_grid, box)
  }
  curve_set <- GET::create_curve_set(list(r = r_grid, obs = observed_k, sim_m = t(K_sim)))
  global_test <- GET::global_envelope_test(curve_set, type = "erl", alternative = "two.sided")
  list(
    box = box, r = r_grid, K_obs = observed_k, K_sim = K_sim,
    K_low = apply(K_sim, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
    K_med = apply(K_sim, 2, stats::median, na.rm = TRUE),
    K_high = apply(K_sim, 2, stats::quantile, probs = 0.975, na.rm = TRUE),
    K0 = (4 / 3) * pi * r_grid^3, get = global_test,
    lambda_observed = lambda_observed, lambda_function = lambda_function
  )
}

plot_get_diagnostic <- function(result, title = NULL) {
  data <- data.frame(r = result$r, observed = result$K_obs, median = result$K_med,
                     lower = result$K_low, upper = result$K_high, K0 = result$K0)
  ggplot2::ggplot(data, ggplot2::aes(x = r)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "grey80") +
    ggplot2::geom_line(ggplot2::aes(y = observed), linewidth = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = median), linetype = "dashed", colour = "steelblue", linewidth = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = K0), linetype = "dotted", colour = "firebrick", linewidth = 0.8) +
    ggplot2::labs(x = expression(rho), y = expression(K[inhom](rho)), title = title) +
    ggplot2::theme_classic(base_size = 12)
}

# -----------------------------------------------------------------------------
# Parametric bootstrap for the reconstructed coefficient function
# -----------------------------------------------------------------------------

bootstrap_functional_coefficient <- function(functional_fit, application_data,
                                             field_function, functional_columns,
                                             mult = 30L, ncube = 20L,
                                             idw_power = 4L, B = 199L,
                                             confidence = 0.95, n_mc = 100000L,
                                             safety_factor = 2, workers = 1L,
                                             seed = 31001L) {
  coordinates <- as.matrix(application_data$functional_pattern$df[, c("x", "y", "z")])
  box <- list(xrange = range(coordinates[, 1]), yrange = range(coordinates[, 2]), zrange = range(coordinates[, 3]))
  lambda_function <- make_lambda_predictor(functional_fit, field_function, functional_columns)
  mc_points <- runif_box3(n_mc, box)
  lambda_max <- max(lambda_function(mc_points), na.rm = TRUE) * safety_factor
  if (!is.finite(lambda_max) || lambda_max <= 0) stop("Invalid bootstrap thinning bound.")
  base_gamma <- extract_score_coefficients(functional_fit)
  base_theta <- reconstruct_functional_coefficient(base_gamma, application_data$fpca, application_data$u_grid)$theta
  formula <- formula(functional_fit$mod_global)

  one_bootstrap <- function(b) {
    set.seed(seed + b)
    volume <- diff(box$xrange) * diff(box$yrange) * diff(box$zrange)
    candidate_coordinates <- runif_box3(stats::rpois(1L, lambda_max * volume), box)
    if (!nrow(candidate_coordinates)) return(NULL)
    candidate_lambda <- lambda_function(candidate_coordinates)
    points <- candidate_coordinates[stats::runif(nrow(candidate_coordinates)) <= pmin(candidate_lambda / lambda_max, 1), , drop = FALSE]
    if (nrow(points) < 20L) return(NULL)
    scores <- field_function(points)
    bootstrap_pattern <- s3dpm(data.frame(points, as.data.frame(scores)), names = functional_columns)
    fit <- try(
      FunPois.fit_v2(
        X = bootstrap_pattern, formula = formula, covs = NULL, marked = FALSE,
        spatial.cov = FALSE, verbose = FALSE, mult = mult, ncube = ncube,
        seed = seed + 100000L + b, grid = FALSE, mark.c = TRUE,
        mark.c.mode = "oracle", mark.oracle_fun = field_function,
        interp.p = idw_power, process.type = "s3d", z_weight = FALSE,
        sphere = FALSE, box3 = box
      ), silent = TRUE
    )
    if (inherits(fit, "try-error")) return(NULL)
    gamma <- try(extract_score_coefficients(fit), silent = TRUE)
    if (inherits(gamma, "try-error")) return(NULL)
    theta <- reconstruct_functional_coefficient(gamma, application_data$fpca, application_data$u_grid)$theta
    list(gamma = gamma, theta = theta, n = nrow(points))
  }

  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  replicates <- future.apply::future_lapply(seq_len(B), one_bootstrap, future.seed = TRUE,
                                             future.packages = c("mgcv"))
  replicates <- replicates[!vapply(replicates, is.null, logical(1))]
  if (!length(replicates)) stop("No bootstrap replicate completed successfully.")
  theta_matrix <- do.call(rbind, lapply(replicates, `[[`, "theta"))
  gamma_matrix <- do.call(rbind, lapply(replicates, `[[`, "gamma"))
  alpha <- 1 - confidence
  theta_sd <- apply(theta_matrix, 2, stats::sd, na.rm = TRUE)
  standardized <- sweep(abs(sweep(theta_matrix, 2, base_theta, "-")), 2, theta_sd, "/")
  simultaneous_critical <- stats::quantile(apply(standardized, 1, max, na.rm = TRUE), confidence, na.rm = TRUE)
  data <- data.frame(
    u = application_data$u_grid,
    theta_hat = base_theta,
    theta_median = apply(theta_matrix, 2, stats::median, na.rm = TRUE),
    theta_low = apply(theta_matrix, 2, stats::quantile, probs = alpha / 2, na.rm = TRUE),
    theta_high = apply(theta_matrix, 2, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE),
    theta_sim_low = base_theta - simultaneous_critical * theta_sd,
    theta_sim_high = base_theta + simultaneous_critical * theta_sd
  )
  list(B = B, B_ok = length(replicates), data = data, theta_matrix = theta_matrix,
       gamma_matrix = gamma_matrix, lambda_max = lambda_max, box = box)
}

extract_significant_intervals <- function(bootstrap_data, lower_column, upper_column,
                                          minimum_grid_points = 2L) {
  state <- ifelse(bootstrap_data[[lower_column]] > 0, "positive",
                  ifelse(bootstrap_data[[upper_column]] < 0, "negative", "crosses_zero"))
  runs <- rle(state)
  end <- cumsum(runs$lengths)
  start <- c(1L, head(end, -1L) + 1L)
  output <- data.frame(state = runs$values, start = start, end = end, n_grid = runs$lengths)
  output$u_from <- bootstrap_data$u[output$start]
  output$u_to <- bootstrap_data$u[output$end]
  output$mean_theta <- mapply(function(a, b) mean(bootstrap_data$theta_hat[a:b]), output$start, output$end)
  output |> dplyr::filter(state != "crosses_zero", n_grid >= minimum_grid_points)
}

plot_bootstrap_bands <- function(bootstrap_result, simultaneous = FALSE) {
  data <- bootstrap_result$data
  lower <- if (simultaneous) "theta_sim_low" else "theta_low"
  upper <- if (simultaneous) "theta_sim_high" else "theta_high"
  median_line <- if (simultaneous) NULL else data$theta_median
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = u)) +
    ggplot2::geom_ribbon(ggplot2::aes_string(ymin = lower, ymax = upper), fill = "grey75") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_line(ggplot2::aes(y = theta_hat), linewidth = 0.9) +
    ggplot2::labs(x = expression(tau), y = expression(hat(theta)[F](tau))) +
    ggplot2::theme_classic(base_size = 12)
  if (!is.null(median_line)) plot <- plot + ggplot2::geom_line(ggplot2::aes(y = theta_median), linetype = "dashed", colour = "steelblue")
  plot
}

# -----------------------------------------------------------------------------
# Figures for the paper
# -----------------------------------------------------------------------------

plot_stars_and_clusters <- function(application_data) {
  stars_plot <- ggplot2::ggplot(application_data$stars, ggplot2::aes(X, Y, colour = A0_ML)) +
    ggplot2::geom_point(size = 0.25, alpha = 0.7) + ggplot2::coord_equal() +
    ggplot2::labs(x = "X [pc]", y = "Y [pc]", colour = expression(A[0])) +
    ggplot2::theme_classic(base_size = 12)
  clusters_plot <- ggplot2::ggplot(application_data$clusters, ggplot2::aes(x, y, colour = application_data$scalar_mean)) +
    ggplot2::geom_point(size = 1.6, alpha = 0.85) + ggplot2::coord_equal() +
    ggplot2::labs(x = "X [pc]", y = "Y [pc]", colour = expression(bar(A)[0](s))) +
    ggplot2::theme_classic(base_size = 12)
  list(stars = stars_plot, clusters = clusters_plot)
}

plot_cluster_curves <- function(application_data, example_cluster = 1L) {
  curve_data <- data.frame(
    u = rep(application_data$u_grid, nrow(application_data$curve_matrix)),
    value = as.vector(t(application_data$curve_matrix)),
    cluster = rep(seq_len(nrow(application_data$curve_matrix)), each = length(application_data$u_grid))
  )
  all_curves <- ggplot2::ggplot(curve_data, ggplot2::aes(u, value, group = cluster)) +
    ggplot2::geom_line(alpha = 0.12, colour = "steelblue") +
    ggplot2::labs(x = expression(tau), y = expression(Q[A[0], c](tau))) +
    ggplot2::theme_classic(base_size = 12)
  observed <- application_data$curve_matrix[example_cluster, ]
  scores <- application_data$fpca$scores[example_cluster, ]
  reconstructed <- application_data$fpca$mu + as.numeric(application_data$fpca$phi %*% scores)
  reconstruction <- data.frame(u = application_data$u_grid, observed = observed, reconstructed = reconstructed)
  example_plot <- ggplot2::ggplot(reconstruction, ggplot2::aes(u)) +
    ggplot2::geom_line(ggplot2::aes(y = observed), linewidth = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = reconstructed), linetype = "dashed", colour = "firebrick", linewidth = 0.8) +
    ggplot2::labs(x = expression(tau), y = expression(Q[A[0], c](tau))) +
    ggplot2::theme_classic(base_size = 12)
  list(all_curves = all_curves, reconstruction = example_plot)
}

plot_application_intensity <- function(functional_fit, application_data, field_function, functional_columns) {
  coordinates <- as.matrix(application_data$functional_pattern$df[, c("x", "y", "z")])
  lambda <- make_lambda_predictor(functional_fit, field_function, functional_columns)(coordinates)
  data <- data.frame(x = coordinates[, 1], y = coordinates[, 2], z = coordinates[, 3], lambda = lambda)
  plot3D::scatter3D(data$x, data$y, data$z, colvar = data$lambda, pch = 19,
                    xlab = "X [pc]", ylab = "Y [pc]", zlab = "Z [pc]")
  invisible(data)
}


plot_cluster_extinction_3d <- function(application_data) {
  data <- application_data$scalar_pattern$df
  plot3D::scatter3D(
    data$x, data$y, data$z, colvar = data$ext_mean, pch = 19,
    xlab = "X [pc]", ylab = "Y [pc]", zlab = "Z [pc]",
    clab = expression(bar(A)[0](s))
  )
  invisible(data)
}

save_application_3d_figures <- function(fits, application_data, directories) {
  save_base_plot_pdf(
    file.path(directories$figures, "Figure_09_cluster_extinction_3D.pdf"), 8, 6,
    function() plot_cluster_extinction_3d(application_data)
  )
  save_base_plot_pdf(
    file.path(directories$figures, "Figure_13_fitted_intensity_3D.pdf"), 8, 6,
    function() plot_application_intensity(
      fits$functional, application_data, fits$functional_field, fits$functional_columns
    )
  )
  invisible(TRUE)
}

APPLICATION_SETTINGS <- list(
  minimum_cluster_size = 30L,
  u_grid = seq(0, 1, length.out = 101),
  K = 4L,
  idw_power = 4L,
  mult = 30L,
  ncube = 20L,
  xy_basis_functional = 30L,
  xy_basis_scalar = 20L,
  z_basis = 10L,
  get_B = 200L,
  bootstrap_B = 199L,
  workers = 1L,
  seed = 123L
)

run_stellar_cluster_application <- function(settings = APPLICATION_SETTINGS,
                                            root = project_root(),
                                            rerun_bootstrap = FALSE) {
  load_simulation_functions(root)
  check_application_packages()
  directories <- make_application_directories(root)
  inputs <- read_application_inputs(root)
  application_data <- prepare_cluster_application_data(
    inputs$stars, inputs$metadata, settings$minimum_cluster_size,
    settings$u_grid, settings$K
  )
  fits <- fit_application_models(
    application_data, idw_power = settings$idw_power, mult = settings$mult,
    ncube = settings$ncube, xy_basis_functional = settings$xy_basis_functional,
    xy_basis_scalar = settings$xy_basis_scalar, z_basis = settings$z_basis, seed = settings$seed
  )
  model_table <- make_model_summary_table(fits$functional, fits$scalar, fits$functional_columns)
  utils::write.csv(model_table$functional_parametric, file.path(directories$tables, "functional_parametric_coefficients.csv"), row.names = FALSE)
  utils::write.csv(model_table$scalar_parametric, file.path(directories$tables, "scalar_parametric_coefficients.csv"), row.names = FALSE)

  get_functional <- run_get_diagnostic(fits$functional, application_data$functional_pattern,
                                       fits$functional_field, fits$functional_columns,
                                       B = settings$get_B, seed = settings$seed + 10L)
  get_scalar <- run_get_diagnostic(fits$scalar, application_data$scalar_pattern,
                                   fits$scalar_field, "ext_mean",
                                   B = settings$get_B, seed = settings$seed + 20L)
  saveRDS(get_functional, file.path(directories$results, "get_functional_model.rds"))
  saveRDS(get_scalar, file.path(directories$results, "get_scalar_model.rds"))

  bootstrap_file <- file.path(directories$results, "functional_bootstrap.rds")
  bootstrap <- if (file.exists(bootstrap_file) && !rerun_bootstrap) readRDS(bootstrap_file) else {
    output <- bootstrap_functional_coefficient(
      fits$functional, application_data, fits$functional_field, fits$functional_columns,
      mult = settings$mult, ncube = settings$ncube, idw_power = settings$idw_power,
      B = settings$bootstrap_B, workers = settings$workers, seed = settings$seed + 30L
    )
    saveRDS(output, bootstrap_file)
    output
  }

  stars_clusters <- plot_stars_and_clusters(application_data)
  curve_plots <- plot_cluster_curves(application_data)
  save_ggplot_pdf(stars_clusters$stars, file.path(directories$figures, "Figure_08_stars.pdf"), 7, 6)
  save_ggplot_pdf(stars_clusters$clusters, file.path(directories$figures, "Figure_08_clusters.pdf"), 7, 6)
  save_ggplot_pdf(curve_plots$all_curves, file.path(directories$figures, "Figure_09_quantile_curves.pdf"), 7, 5)
  save_ggplot_pdf(curve_plots$reconstruction, file.path(directories$figures, "Figure_10_example_reconstruction.pdf"), 7, 5)
  save_ggplot_pdf(plot_get_diagnostic(get_functional), file.path(directories$figures, "Figure_11_GET_functional.pdf"), 7, 5)
  save_ggplot_pdf(plot_get_diagnostic(get_scalar), file.path(directories$figures, "Figure_11_GET_scalar.pdf"), 7, 5)
  save_ggplot_pdf(plot_bootstrap_bands(bootstrap, simultaneous = FALSE), file.path(directories$figures, "Figure_12_pointwise_band.pdf"), 7, 5)
  save_ggplot_pdf(plot_bootstrap_bands(bootstrap, simultaneous = TRUE), file.path(directories$figures, "Figure_12_simultaneous_band.pdf"), 7, 5)
  save_application_3d_figures(fits, application_data, directories)

  list(application_data = application_data, fits = fits, model_table = model_table,
       get_functional = get_functional, get_scalar = get_scalar,
       bootstrap = bootstrap, directories = directories)
}
