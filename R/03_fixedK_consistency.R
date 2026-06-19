# ============================================================================
# Functional spatial covariates in three-dimensional Poisson point processes
# Fixed-K consistency study
# ============================================================================
#
# This script reproduces the fixed-K Monte Carlo study supporting the consistency
# theorem. It uses checkpointed blocks so an interrupted run can resume without
# discarding completed replications.
#
# The script sources `01_simulation_study.R` when the required simulation
# functions are not already available in the current R session.
# ============================================================================

load_consistency_dependencies <- function(root = Sys.getenv("PROJECT_ROOT", unset = getwd())) {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  required <- c("trapz_weights", "w_center", "make_true_phi_fourier",
                "make_beta_scenarios", "run_one_rep_fpca")
  missing <- required[!vapply(required, exists, logical(1), mode = "function")]
  if (length(missing)) source(file.path(root, "R", "01_simulation_study.R"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install package 'ggplot2' before running the consistency workflow.")
  }
  invisible(TRUE)
}

load_consistency_dependencies()
.has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)

CONSISTENCY_SETTINGS <- list(
  base_seed = 1234L, nsim = 100L, target_N_seq = c(100L, 250L, 500L, 1000L, 1500L),
  u_grid = seq(0, 1, length.out = 101), beta_name = "single", K0 = 4L, K_fixed = 4L,
  range = 0.25, nugget = 0.05, sigma_eps = 0.05, trunc_eps = 4, trunc_xi = 4, L_rff = 128L,
  fpca_basis_mode = "event", interp.p = 16, mult_fit = 30L, mark_modes = c("oracle", "idw"),
  mc_n = 20000L, n_eval = 5000L, beta_xyz = c(0, 0, 0)
)

# -----------------------------------------------------------------------------
# Checkpointed consistency simulation functions
# -----------------------------------------------------------------------------

safe_save_rds <- function(object, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(file, ".tmp")
  saveRDS(object, tmp)
  if (file.exists(file)) unlink(file)
  ok <- file.rename(tmp, file)
  if (!ok) {
    stop("Could not move temporary file into place: ", file)
  }
  invisible(file)
}

project_curve_on_true_basis <- function(beta_curve, Phi_true, w_u) {
  as.numeric(crossprod(Phi_true * w_u, beta_curve))
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x)
}

safe_rmse <- function(err) {
  err <- err[is.finite(err)]
  if (!length(err)) return(NA_real_)
  sqrt(mean(err^2))
}

block_id_consistency <- function(target_N, mark_mode) {
  sprintf("N%04d_%s", as.integer(target_N), as.character(mark_mode))
}

block_file_consistency <- function(save_dir, target_N, mark_mode) {
  file.path(save_dir, paste0(block_id_consistency(target_N, mark_mode), ".rds"))
}

build_consistency_truth <- function(settings = CONSISTENCY_SETTINGS) {
  u_grid <- settings$u_grid
  w_u    <- trapz_weights(u_grid)

  Phi_true <- make_true_phi_fourier(u_grid, K0 = settings$K0)

  beta_list <- make_beta_scenarios(u_grid)
  if (!(settings$beta_name %in% names(beta_list))) {
    stop("beta_name '", settings$beta_name, "' not found in make_beta_scenarios().")
  }

  beta_true <- beta_list[[settings$beta_name]]
  beta_true <- w_center(beta_true, w_u)

  gamma_true_basis <- project_curve_on_true_basis(beta_true, Phi_true, w_u)

  list(
    u_grid = u_grid,
    w_u = w_u,
    Phi_true = Phi_true,
    beta_true = beta_true,
    gamma_true_basis = gamma_true_basis
  )
}

run_one_rep_consistency_fixedK <- function(target_N,
                                           rep_id,
                                           mark_mode = c("oracle", "idw"),
                                           settings = CONSISTENCY_SETTINGS,
                                           truth = NULL) {
  mark_mode <- match.arg(mark_mode)

  if (is.null(truth)) truth <- build_consistency_truth(settings)

  seed_rep <- settings$base_seed +
    100000L * match(mark_mode, settings$mark_modes) +
    1000L   * as.integer(target_N) +
    as.integer(rep_id)

  out <- run_one_rep_fpca(
    beta_u = truth$beta_true,
    u_grid = truth$u_grid,
    Phi_true = truth$Phi_true,
    range = settings$range,
    nugget = settings$nugget,
    target_N = as.integer(target_N),
    beta_xyz = settings$beta_xyz,

    K_fpca = settings$K_fixed,
    interp.p = settings$interp.p,
    mult_fit = settings$mult_fit,

    fpca_basis_mode = settings$fpca_basis_mode,
    mark_mode = mark_mode,

    K0 = settings$K0,
    sigma_eps = settings$sigma_eps,
    trunc_eps = settings$trunc_eps,
    trunc_xi  = settings$trunc_xi,
    L_rff = settings$L_rff,

    n_eval = settings$n_eval,
    mc_n   = settings$mc_n,
    seed   = seed_rep
  )

  if (isTRUE(out$ok)) {
    gamma_hat_true_basis <- project_curve_on_true_basis(
      beta_curve = out$beta_hat,
      Phi_true   = truth$Phi_true,
      w_u        = truth$w_u
    )

    gamma_err_true_basis <- gamma_hat_true_basis - truth$gamma_true_basis

    out$gamma_hat_true_basis <- gamma_hat_true_basis
    out$gamma_true_basis     <- truth$gamma_true_basis
    out$gamma_err_true_basis <- gamma_err_true_basis

    out$scaled_gamma_err_true_basis <- sqrt(target_N) * gamma_err_true_basis
    out$scaled_beta_rmse            <- sqrt(target_N) * out$rmse_beta
    out$scaled_lambda_rmse          <- sqrt(target_N) * out$rmse_lambda
  } else {
    out$gamma_hat_true_basis        <- rep(NA_real_, settings$K_fixed)
    out$gamma_true_basis            <- truth$gamma_true_basis
    out$gamma_err_true_basis        <- rep(NA_real_, settings$K_fixed)
    out$scaled_gamma_err_true_basis <- rep(NA_real_, settings$K_fixed)
    out$scaled_beta_rmse            <- NA_real_
    out$scaled_lambda_rmse          <- NA_real_
  }

  out$target_N  <- as.integer(target_N)
  out$rep_id    <- as.integer(rep_id)
  out$mark_mode <- mark_mode

  out
}

run_consistency_block_fixedK <- function(target_N,
                                         mark_mode,
                                         settings = CONSISTENCY_SETTINGS,
                                         truth = NULL,
                                         save_dir = "consistency_fixedK_blocks",
                                         overwrite = FALSE,
                                         verbose = TRUE) {
  if (is.null(truth)) truth <- build_consistency_truth(settings)

  block_file <- block_file_consistency(save_dir, target_N, mark_mode)
  block_name <- block_id_consistency(target_N, mark_mode)

  if (file.exists(block_file) && !overwrite) {
    block_obj <- readRDS(block_file)


    same_nsim <- identical(as.integer(block_obj$settings$nsim), as.integer(settings$nsim))
    same_K    <- identical(as.integer(block_obj$settings$K_fixed), as.integer(settings$K_fixed))
    same_N    <- identical(as.integer(block_obj$target_N), as.integer(target_N))
    same_mode <- identical(as.character(block_obj$mark_mode), as.character(mark_mode))

    if (!(same_nsim && same_K && same_N && same_mode)) {
      stop("Existing block file is incompatible with current settings: ", block_file)
    }

    results <- block_obj$results
    done <- which(vapply(results, function(x) !is.null(x), logical(1)))
    start_rep <- if (length(done)) max(done) + 1L else 1L

    if (isTRUE(verbose)) {
      cat(sprintf("[RESUME BLOCK %s] completed %d / %d reps\n",
                  block_name, length(done), settings$nsim))
    }
  } else {
    results <- vector("list", settings$nsim)
    start_rep <- 1L

    block_obj <- list(
      block_id   = block_name,
      target_N   = as.integer(target_N),
      mark_mode  = as.character(mark_mode),
      settings   = settings,
      started_at = Sys.time(),
      last_update = Sys.time(),
      completed_reps = integer(0),
      results = results
    )

    safe_save_rds(block_obj, block_file)

    if (isTRUE(verbose)) {
      cat(sprintf("[NEW BLOCK %s] starting %d reps\n", block_name, settings$nsim))
    }
  }

  if (start_rep > settings$nsim) {
    if (isTRUE(verbose)) {
      cat(sprintf("[SKIP BLOCK %s] already complete\n", block_name))
    }
    return(invisible(readRDS(block_file)))
  }

  block_start_time <- Sys.time()

  for (rep_id in seq.int(start_rep, settings$nsim)) {
    rep_start <- Sys.time()

    if (isTRUE(verbose)) {
      cat(sprintf("  -> [%s] rep %d / %d ... ",
                  block_name, rep_id, settings$nsim))
      flush.console()
    }

    res <- tryCatch(
      run_one_rep_consistency_fixedK(
        target_N = target_N,
        rep_id   = rep_id,
        mark_mode = mark_mode,
        settings = settings,
        truth = truth
      ),
      error = function(e) {
        list(
          ok = FALSE,
          error_message = conditionMessage(e),
          target_N = as.integer(target_N),
          rep_id = as.integer(rep_id),
          mark_mode = as.character(mark_mode),
          gamma_hat_true_basis = rep(NA_real_, settings$K_fixed),
          gamma_true_basis = truth$gamma_true_basis,
          gamma_err_true_basis = rep(NA_real_, settings$K_fixed),
          scaled_gamma_err_true_basis = rep(NA_real_, settings$K_fixed),
          scaled_beta_rmse = NA_real_,
          scaled_lambda_rmse = NA_real_
        )
      }
    )

    block_obj$results[[rep_id]] <- res
    block_obj$completed_reps <- which(vapply(block_obj$results, function(x) !is.null(x), logical(1)))
    block_obj$last_update <- Sys.time()

    safe_save_rds(block_obj, block_file)

    if (isTRUE(verbose)) {
      rep_elapsed <- difftime(Sys.time(), rep_start, units = "secs")
      total_elapsed <- difftime(Sys.time(), block_start_time, units = "mins")
      cat(sprintf("done (%.1f sec | block elapsed %.2f min)\n",
                  as.numeric(rep_elapsed), as.numeric(total_elapsed)))
    }
  }

  if (isTRUE(verbose)) {
    cat(sprintf("[DONE BLOCK %s] saved in %s\n", block_name, block_file))
  }

  invisible(readRDS(block_file))
}

run_consistency_simulation_fixedK <- function(settings = CONSISTENCY_SETTINGS,
                                              save_dir = "consistency_fixedK_blocks",
                                              overwrite = FALSE,
                                              verbose = TRUE) {
  truth <- build_consistency_truth(settings)

  task_grid <- expand.grid(
    target_N = settings$target_N_seq,
    mark_mode = settings$mark_modes,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(verbose)) {
    cat("============================================================\n")
    cat("Consistency simulation with block checkpointing\n")
    cat("Save directory:", normalizePath(save_dir, winslash = "/", mustWork = FALSE), "\n")
    cat("Blocks to run:", nrow(task_grid), "\n")
    cat("Replicates per block:", settings$nsim, "\n")
    cat("============================================================\n")
  }

  for (i in seq_len(nrow(task_grid))) {
    tt <- task_grid[i, ]
    run_consistency_block_fixedK(
      target_N  = tt$target_N,
      mark_mode = tt$mark_mode,
      settings  = settings,
      truth     = truth,
      save_dir  = save_dir,
      overwrite = overwrite,
      verbose   = verbose
    )
  }

  invisible(TRUE)
}

load_consistency_blocks_fixedK <- function(settings = CONSISTENCY_SETTINGS,
                                           save_dir = "consistency_fixedK_blocks",
                                           require_complete = TRUE) {
  truth <- build_consistency_truth(settings)

  task_grid <- expand.grid(
    target_N = settings$target_N_seq,
    mark_mode = settings$mark_modes,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  raw_list <- list()
  task_rows <- list()
  idx <- 0L

  for (i in seq_len(nrow(task_grid))) {
    tt <- task_grid[i, ]
    f <- block_file_consistency(save_dir, tt$target_N, tt$mark_mode)

    if (!file.exists(f)) {
      if (require_complete) {
        stop("Missing block file: ", f)
      } else {
        next
      }
    }

    block_obj <- readRDS(f)
    results <- block_obj$results

    done <- which(vapply(results, function(x) !is.null(x), logical(1)))

    if (require_complete && length(done) < settings$nsim) {
      stop("Incomplete block file: ", f, " (", length(done), "/", settings$nsim, " reps done)")
    }

    for (j in done) {
      idx <- idx + 1L
      nm <- sprintf("N%s_%s_rep%03d", tt$target_N, tt$mark_mode, j)
      raw_list[[idx]] <- results[[j]]
      names(raw_list)[idx] <- nm

      task_rows[[idx]] <- data.frame(
        target_N = tt$target_N,
        mark_mode = tt$mark_mode,
        rep_id = j,
        stringsAsFactors = FALSE
      )
    }
  }

  tasks <- if (length(task_rows)) do.call(rbind, task_rows) else data.frame()

  list(
    settings = settings,
    truth = truth,
    tasks = tasks,
    raw = raw_list
  )
}

extract_consistency_tables <- function(sim_obj) {
  settings <- sim_obj$settings
  truth    <- sim_obj$truth
  raw      <- sim_obj$raw

  metrics_rows <- vector("list", length(raw))
  gamma_rows   <- vector("list", length(raw) * settings$K_fixed)
  kk <- 0L

  for (i in seq_along(raw)) {
    rr <- raw[[i]]

    metrics_rows[[i]] <- data.frame(
      target_N      = rr$target_N,
      mark_mode     = rr$mark_mode,
      rep_id        = rr$rep_id,
      ok            = isTRUE(rr$ok),
      n_obs         = rr$n %||% NA_integer_,
      rmse_beta     = rr$rmse_beta %||% NA_real_,
      ise_beta      = rr$ise_beta %||% NA_real_,
      rmse_lambda   = rr$rmse_lambda %||% NA_real_,
      mise_lambda   = rr$mise_lambda %||% NA_real_,
      logrmse_lambda = rr$logrmse_lambda %||% NA_real_,
      I_true        = rr$I_true %||% NA_real_,
      I_hat         = rr$I_hat %||% NA_real_,
      relerr_I      = rr$relerr_I %||% NA_real_,
      pve_cum_last  = rr$pve_cum_last %||% NA_real_,
      scaled_beta_rmse   = rr$scaled_beta_rmse %||% NA_real_,
      scaled_lambda_rmse = rr$scaled_lambda_rmse %||% NA_real_,
      error_message = rr$error_message %||% NA_character_,
      stringsAsFactors = FALSE
    )

    for (k in seq_len(settings$K_fixed)) {
      kk <- kk + 1L
      gamma_rows[[kk]] <- data.frame(
        target_N   = rr$target_N,
        mark_mode  = rr$mark_mode,
        rep_id     = rr$rep_id,
        ok         = isTRUE(rr$ok),
        component  = paste0("k", k),
        gamma_true = rr$gamma_true_basis[k],
        gamma_hat  = rr$gamma_hat_true_basis[k],
        gamma_err  = rr$gamma_err_true_basis[k],
        scaled_gamma_err = rr$scaled_gamma_err_true_basis[k],
        stringsAsFactors = FALSE
      )
    }
  }

  metrics_df <- do.call(rbind, metrics_rows)
  gamma_df   <- do.call(rbind, gamma_rows)

  beta_curve_list <- lapply(raw, function(rr) {
    list(
      target_N  = rr$target_N,
      mark_mode = rr$mark_mode,
      rep_id    = rr$rep_id,
      ok        = isTRUE(rr$ok),
      beta_hat  = rr$beta_hat %||% NULL,
      beta_true = truth$beta_true
    )
  })

  list(
    metrics_df = metrics_df,
    gamma_df = gamma_df,
    beta_curve_list = beta_curve_list
  )
}

summarise_consistency_metrics <- function(metrics_df) {
  split_key <- interaction(metrics_df$target_N, metrics_df$mark_mode, drop = TRUE)
  sp <- split(metrics_df, split_key)

  out <- lapply(sp, function(df) {
    ok <- df$ok

    data.frame(
      target_N     = unique(df$target_N),
      mark_mode    = unique(df$mark_mode),
      n_rep        = nrow(df),
      success_rate = mean(ok),
      mean_n_obs   = safe_mean(df$n_obs[ok]),
      sd_n_obs     = safe_sd(df$n_obs[ok]),

      mean_rmse_beta   = safe_mean(df$rmse_beta[ok]),
      median_rmse_beta = safe_median(df$rmse_beta[ok]),
      sd_rmse_beta     = safe_sd(df$rmse_beta[ok]),

      mean_ise_beta    = safe_mean(df$ise_beta[ok]),
      mean_rmse_lambda = safe_mean(df$rmse_lambda[ok]),
      median_rmse_lambda = safe_median(df$rmse_lambda[ok]),
      sd_rmse_lambda   = safe_sd(df$rmse_lambda[ok]),

      mean_mise_lambda = safe_mean(df$mise_lambda[ok]),
      mean_logrmse_lambda = safe_mean(df$logrmse_lambda[ok]),

      mean_relerr_I    = safe_mean(df$relerr_I[ok]),
      median_relerr_I  = safe_median(df$relerr_I[ok]),
      sd_relerr_I      = safe_sd(df$relerr_I[ok]),

      mean_scaled_beta_rmse   = safe_mean(df$scaled_beta_rmse[ok]),
      mean_scaled_lambda_rmse = safe_mean(df$scaled_lambda_rmse[ok]),

      mean_pve_cum_last = safe_mean(df$pve_cum_last[ok]),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$mark_mode, out$target_N), ]
}

summarise_consistency_gamma <- function(gamma_df) {
  split_key <- interaction(gamma_df$target_N, gamma_df$mark_mode, gamma_df$component, drop = TRUE)
  sp <- split(gamma_df, split_key)

  out <- lapply(sp, function(df) {
    ok <- df$ok
    err <- df$gamma_err[ok]
    serr <- df$scaled_gamma_err[ok]

    data.frame(
      target_N    = unique(df$target_N),
      mark_mode   = unique(df$mark_mode),
      component   = unique(df$component),
      n_rep       = nrow(df),
      success_rate = mean(ok),

      bias_gamma  = safe_mean(err),
      sd_gamma    = safe_sd(df$gamma_hat[ok]),
      rmse_gamma  = safe_rmse(err),

      mean_scaled_gamma_err = safe_mean(serr),
      sd_scaled_gamma_err   = safe_sd(serr),

      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$component, out$mark_mode, out$target_N), ]
}

build_oracle_idw_gap <- function(summary_metrics, summary_gamma = NULL) {
  oracle <- summary_metrics[summary_metrics$mark_mode == "oracle", ]
  idw    <- summary_metrics[summary_metrics$mark_mode == "idw", ]

  commonN <- intersect(oracle$target_N, idw$target_N)

  gap_metrics <- merge(
    oracle[oracle$target_N %in% commonN, ],
    idw[idw$target_N %in% commonN, ],
    by = "target_N",
    suffixes = c("_oracle", "_idw")
  )

  gap_metrics$gap_rmse_beta   <- gap_metrics$mean_rmse_beta_idw   - gap_metrics$mean_rmse_beta_oracle
  gap_metrics$gap_rmse_lambda <- gap_metrics$mean_rmse_lambda_idw - gap_metrics$mean_rmse_lambda_oracle
  gap_metrics$gap_relerr_I    <- gap_metrics$mean_relerr_I_idw    - gap_metrics$mean_relerr_I_oracle

  out <- list(metrics_gap = gap_metrics)

  if (!is.null(summary_gamma)) {
    sg_or <- summary_gamma[summary_gamma$mark_mode == "oracle", ]
    sg_id <- summary_gamma[summary_gamma$mark_mode == "idw", ]

    gap_gamma <- merge(
      sg_or,
      sg_id,
      by = c("target_N", "component"),
      suffixes = c("_oracle", "_idw")
    )

    gap_gamma$gap_bias_gamma <- gap_gamma$bias_gamma_idw - gap_gamma$bias_gamma_oracle
    gap_gamma$gap_rmse_gamma <- gap_gamma$rmse_gamma_idw - gap_gamma$rmse_gamma_oracle

    out$gamma_gap <- gap_gamma
  }

  out
}

plot_consistency_metric <- function(summary_metrics,
                                    metric = c("mean_rmse_beta",
                                               "mean_rmse_lambda",
                                               "mean_relerr_I",
                                               "success_rate"),
                                    title = NULL) {
  metric <- match.arg(metric)
  if (!.has_ggplot2) stop("ggplot2 is required for plotting.")

  lab_map <- c(
    mean_rmse_beta   = "Mean RMSE of reconstructed beta(u)",
    mean_rmse_lambda = "Mean RMSE of intensity",
    mean_relerr_I    = "Mean relative error on integral intensity",
    success_rate     = "Success rate"
  )

  ggplot2::ggplot(summary_metrics,
                  ggplot2::aes(x = target_N, y = .data[[metric]], color = mark_mode)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_x_continuous(breaks = sort(unique(summary_metrics$target_N))) +
    ggplot2::labs(
      x = expression(N[target]),
      y = lab_map[[metric]],
      color = "Mode",
      title = title %||% lab_map[[metric]]
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "top"
    )
}

plot_consistency_gamma <- function(summary_gamma,
                                   value = c("bias_gamma", "rmse_gamma"),
                                   title = NULL) {
  value <- match.arg(value)
  if (!.has_ggplot2) stop("ggplot2 is required for plotting.")

  lab_map <- c(
    bias_gamma = "Bias of coefficient on true basis",
    rmse_gamma = "RMSE of coefficient on true basis"
  )

  ggplot2::ggplot(summary_gamma,
                  ggplot2::aes(x = target_N, y = .data[[value]], color = mark_mode)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.3) +
    ggplot2::facet_wrap(~ component, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = sort(unique(summary_gamma$target_N))) +
    ggplot2::labs(
      x = expression(N[target]),
      y = lab_map[[value]],
      color = "Mode",
      title = title %||% lab_map[[value]]
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "top"
    )
}

plot_oracle_idw_gap <- function(gap_metrics,
                                metric = c("gap_rmse_beta", "gap_rmse_lambda", "gap_relerr_I"),
                                title = NULL) {
  metric <- match.arg(metric)
  if (!.has_ggplot2) stop("ggplot2 is required for plotting.")

  lab_map <- c(
    gap_rmse_beta   = "IDW - Oracle gap: beta RMSE",
    gap_rmse_lambda = "IDW - Oracle gap: intensity RMSE",
    gap_relerr_I    = "IDW - Oracle gap: relative integral error"
  )

  ggplot2::ggplot(gap_metrics,
                  ggplot2::aes(x = target_N, y = .data[[metric]])) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
    ggplot2::geom_line(linewidth = 1, color = "steelblue4") +
    ggplot2::geom_point(size = 2.5, color = "steelblue4") +
    ggplot2::scale_x_continuous(breaks = sort(unique(gap_metrics$target_N))) +
    ggplot2::labs(
      x = expression(N[target]),
      y = lab_map[[metric]],
      title = title %||% lab_map[[metric]]
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold")
    )
}

plot_scaled_gamma_boxplot <- function(gamma_df) {
  if (!.has_ggplot2) stop("ggplot2 is required for plotting.")

  ggplot2::ggplot(
    gamma_df[gamma_df$ok, ],
    ggplot2::aes(x = factor(target_N), y = scaled_gamma_err, fill = mark_mode)
  ) +
    ggplot2::geom_boxplot(outlier.size = 0.6) +
    ggplot2::facet_wrap(~ component, scales = "free_y") +
    ggplot2::labs(
      x = expression(N[target]),
      y = expression(sqrt(N[target]) * (hat(gamma) - gamma)),
      fill = "Mode",
      title = "Scaled coefficient errors on the true basis"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "top"
    )
}

run_consistency_workflow_fixedK <- function(settings = CONSISTENCY_SETTINGS,
                                            save_dir = "consistency_fixedK_blocks",
                                            overwrite = FALSE,
                                            verbose = TRUE,
                                            require_complete = TRUE,
                                            save_final_rds = TRUE,
                                            final_rds_name = "consistency_fixedK_workflow.rds") {

  run_consistency_simulation_fixedK(
    settings = settings,
    save_dir = save_dir,
    overwrite = overwrite,
    verbose = verbose
  )


  sim_obj <- load_consistency_blocks_fixedK(
    settings = settings,
    save_dir = save_dir,
    require_complete = require_complete
  )


  tables <- extract_consistency_tables(sim_obj)

  summary_metrics <- summarise_consistency_metrics(tables$metrics_df)
  summary_gamma   <- summarise_consistency_gamma(tables$gamma_df)
  gap_obj         <- build_oracle_idw_gap(summary_metrics, summary_gamma)

  plots <- list()
  if (.has_ggplot2) {
    plots$rmse_beta    <- plot_consistency_metric(summary_metrics, "mean_rmse_beta")
    plots$rmse_lambda  <- plot_consistency_metric(summary_metrics, "mean_rmse_lambda")
    plots$relerr_I     <- plot_consistency_metric(summary_metrics, "mean_relerr_I")
    plots$success_rate <- plot_consistency_metric(summary_metrics, "success_rate")
    plots$gamma_bias   <- plot_consistency_gamma(summary_gamma, "bias_gamma")
    plots$gamma_rmse   <- plot_consistency_gamma(summary_gamma, "rmse_gamma")
    plots$gap_beta     <- plot_oracle_idw_gap(gap_obj$metrics_gap, "gap_rmse_beta")
    plots$gap_lambda   <- plot_oracle_idw_gap(gap_obj$metrics_gap, "gap_rmse_lambda")
    plots$scaled_gamma <- plot_scaled_gamma_boxplot(tables$gamma_df)
  }

  out <- list(
    settings = settings,
    save_dir = save_dir,
    truth = sim_obj$truth,
    tasks = sim_obj$tasks,
    raw = sim_obj$raw,
    metrics_df = tables$metrics_df,
    gamma_df = tables$gamma_df,
    beta_curve_list = tables$beta_curve_list,
    summary_metrics = summary_metrics,
    summary_gamma = summary_gamma,
    gap = gap_obj,
    plots = plots
  )

  if (isTRUE(save_final_rds)) {
    safe_save_rds(out, file.path(save_dir, final_rds_name))
  }

  out
}


# -----------------------------------------------------------------------------
# Settings used in the manuscript
# -----------------------------------------------------------------------------


save_consistency_figures <- function(consistency_output, root = project_root()) {
  output_directory <- file.path(root, "figures", "supplement")
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  save_ggplot_pdf(consistency_output$plots$rmse_beta,
                  file.path(output_directory, "Figure_14_consistency_beta_rmse.pdf"), 7, 5)
  save_ggplot_pdf(consistency_output$plots$gamma_rmse,
                  file.path(output_directory, "Figure_15_consistency_gamma_rmse.pdf"), 8, 6)
  save_ggplot_pdf(consistency_output$plots$gap_beta,
                  file.path(output_directory, "Figure_16_oracle_idw_gap.pdf"), 7, 5)
  invisible(output_directory)
}

run_fixedK_consistency_study <- function(settings = CONSISTENCY_SETTINGS,
                                         root = project_root(),
                                         overwrite = FALSE) {
  output_directory <- file.path(root, "results", "consistency")
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  output <- run_consistency_workflow_fixedK(
    settings = settings,
    save_dir = output_directory,
    overwrite = overwrite,
    verbose = TRUE,
    require_complete = TRUE,
    save_final_rds = TRUE,
    final_rds_name = "fixedK_consistency_workflow.rds"
  )
  save_consistency_figures(output, root)
  utils::write.csv(output$summary_metrics,
                   file.path(root, "tables", "consistency_summary_metrics.csv"), row.names = FALSE)
  utils::write.csv(output$summary_gamma,
                   file.path(root, "tables", "consistency_summary_gamma.csv"), row.names = FALSE)
  output
}
