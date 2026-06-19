# ============================================================================
# Functional spatial covariates in three-dimensional Poisson point processes
# Simulation study code
# ============================================================================
#
# This script contains the complete simulation workflow used for the manuscript:
#   (1) the main recovery experiment;
#   (2) the FPCA-truncation and IDW-completion sensitivity analysis;
#   (3) the comparison with scalar-summary and coordinate-only models; and
#   (4) the robustness analysis for FPCA construction and score completion.
#
# The script defines functions only. It does not launch computationally intensive
# simulations when sourced. Run `run_simulation_study()` explicitly after setting
# the paths and computational settings near the end of this file.
# ============================================================================

check_simulation_packages <- function() {
  required <- c(
    "mgcv", "splines", "spatstat.geom", "splancs",
    "future", "future.apply", "ggplot2", "dplyr", "tidyr", "plot3D"
  )
  unavailable <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(unavailable)) {
    stop("Install the following R packages before running this script: ",
         paste(unavailable, collapse = ", "))
  }
  invisible(required)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

project_root <- function() {
  root <- Sys.getenv("PROJECT_ROOT", unset = getwd())
  normalizePath(root, winslash = "/", mustWork = FALSE)
}

make_output_directories <- function(root = project_root()) {
  directories <- c(
    file.path(root, "results", "simulation"),
    file.path(root, "figures", "simulation"),
    file.path(root, "tables", "simulation")
  )
  invisible(lapply(directories, dir.create, recursive = TRUE, showWarnings = FALSE))
  list(
    results = directories[1],
    figures = directories[2],
    tables = directories[3]
  )
}

save_ggplot_pdf <- function(plot_object, filename, width, height) {
  ggplot2::ggsave(filename, plot_object, width = width, height = height,
                  units = "in", device = grDevices::cairo_pdf)
  invisible(filename)
}

save_base_plot_pdf <- function(filename, width, height, plotting_function) {
  grDevices::pdf(filename, width = width, height = height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  plotting_function()
  invisible(filename)
}

# -----------------------------------------------------------------------------
# Cubature-based three-dimensional Poisson estimator
# -----------------------------------------------------------------------------

build_points_covs <- function(quad_p, covs, formula, type.cov.values,
                              interp.p = 81, process.type = c("s2d","st","s3d")) {
  points.covs <- as.data.frame(quad_p)


  if (missing(process.type) || is.null(process.type)) {
    nms <- names(points.covs)
    if ("t" %in% nms) {
      process.type <- "st"
    } else if ("z" %in% nms) {
      process.type <- "s3d"
    } else {
      process.type <- "s2d"
    }
  } else {
    process.type <- match.arg(tolower(process.type), c("s2d","st","s3d"))
  }


  allowed_coords <- switch(
    process.type,
    "s2d" = c("x","y"),
    "st"  = c("x","y","t"),
    "s3d" = c("x","y","z")
  )

  cov_names <- intersect(names(covs), all.vars(formula))
  if (!length(cov_names)) {
    stop("No external covariates found in 'covs' following the model formula.")
  }

  allowed_types <- c("min","max","interp")

  coord_order <- function(nms) intersect(allowed_coords, nms)

  coerce_time_numeric <- function(df, cols) {
    if ("t" %in% cols && inherits(df$t, c("POSIXct","POSIXt","Date"))) {
      df$t <- as.numeric(df$t)
    }
    df
  }

  for (nm in cov_names) {
    choice <- type.cov.values[[nm]]
    if (is.null(choice) || !(choice %in% allowed_types)) {
      stop("type.cov.values for covariate '", nm, "' must be one of: ",
           paste(allowed_types, collapse = ", "))
    }

    if (choice %in% c("min","max")) {

      dfc0 <- covs[[nm]]
      if (is.null(dfc0)) {
        stop("Covariate '", nm, "' not found in 'covs'.")
      }
      present <- coord_order(names(dfc0))
      if (length(present) == 0L || length(present) > length(allowed_coords)) {
        stop(sprintf(
          "Cov '%s': 1-%d coordinate columns needed in external covariate dataset with colnames in %s.",
          nm, length(allowed_coords), paste(allowed_coords, collapse = ", ")
        ))
      }
      if (!all(present %in% names(as.data.frame(quad_p)))) {
        miss <- setdiff(present, names(as.data.frame(quad_p)))
        stop(sprintf("Cov '%s': columns missing in data points: %s",
                     nm, paste(miss, collapse = ", ")))
      }
      dfp <- quad_p[, present, drop = FALSE]
      dfc <- dfc0[, present, drop = FALSE]
      dfp <- coerce_time_numeric(dfp, present)
      dfc <- coerce_time_numeric(dfc, present)

      D <- distN_fast(dfp, dfc, d = length(present), coord_names = present)
      agg_fun <- if (choice == "min") min else max
      vals <- apply(D, 1L, agg_fun, na.rm = TRUE)
      points.covs[[nm]] <- vals

    } else if (choice == "interp") {

      dfc0 <- covs[[nm]]
      if (is.null(dfc0)) {
        stop("Covariate '", nm, "' not found in 'covs'.")
      }

      coord_p <- coord_order(names(as.data.frame(quad_p)))
      if (!length(coord_p)) {
        stop("Missing coordinates ", paste(allowed_coords, collapse = "/"),
             " in data points for covariate '", nm, "'.")
      }
      req <- c(coord_p, nm)
      if (!all(req %in% names(dfc0))) {
        miss <- setdiff(req, names(dfc0))
        stop(sprintf("Cov '%s' (interp): missing columns in covariate dataset: %s",
                     nm, paste(miss, collapse = ", ")))
      }

      pts <- quad_p[, coord_p, drop = FALSE]
      dfc <- dfc0[, req, drop = FALSE]
      pts <- coerce_time_numeric(pts, coord_p)
      dfc <- coerce_time_numeric(dfc, coord_p)

      cov.pred <- idw_interp_fast(
        points = pts,
        covs   = dfc,
        p      = interp.p,
        d      = length(coord_p)
      )
      points.covs[[nm]] <- cov.pred

    } else {
      stop("No managed type for covariate '", nm, "': ", choice)
    }
  }

  points.covs
}

.select_coords <- function(df, d = NULL, coord_names = c("x","y","z","t"), require_all = FALSE) {
  if (!is.null(d)) {
    if (ncol(df) < d) stop("Argument 'd' exceeds the number of available columns.")
    return(as.matrix(df[, seq_len(d), drop = FALSE]))
  }
  if (!is.null(colnames(df))) {
    if (require_all) {
      if (!all(coord_names %in% colnames(df))) stop("Missing required coordinates: ", paste(coord_names, collapse = ", "))
      return(as.matrix(df[, coord_names, drop = FALSE]))
    } else {
      cn <- intersect(coord_names, colnames(df))
      if (length(cn) == 0L) stop("Unable to infer coordinate columns. Provide 'd' or name columns among {x,y,z,t}.")
      return(as.matrix(df[, cn, drop = FALSE]))
    }
  }
  stop("Provide 'd' or column names to infer coordinates.")
}

distN_fast <- function(A, B, d = NULL, coord_names = c("x","y","z","t"), require_all = FALSE) {
  if (is.null(d)) {
    A <- .select_coords(A, coord_names = coord_names, require_all = require_all)
    B <- .select_coords(B, coord_names = colnames(A), require_all = TRUE)
  } else {
    A <- .select_coords(A, d = d, coord_names = coord_names, require_all = require_all)
    B <- .select_coords(B, d = d, coord_names = coord_names, require_all = require_all)
  }
  if (!is.null(colnames(A)) && "t" %in% colnames(A) && !is.numeric(A[, "t"])) A[, "t"] <- as.numeric(A[, "t"])
  if (!is.null(colnames(B)) && "t" %in% colnames(B) && !is.numeric(B[, "t"])) B[, "t"] <- as.numeric(B[, "t"])
  storage.mode(A) <- "double"
  storage.mode(B) <- "double"
  A2 <- rowSums(A^2)
  B2 <- rowSums(B^2)
  D2 <- outer(A2, B2, "+") - 2 * (A %*% t(B))
  D2[D2 < 0] <- 0
  sqrt(D2)
}

idw_interp_fast <- function(points, covs, p = 2, d = NULL,
                            coord_names = c("x","y","z","t"),
                            eps = 1e-10, exact_on_zero = TRUE, na.rm = TRUE,
                            require_all = FALSE, rescale01 = FALSE) {
  if (ncol(covs) < 2) stop("'covs' must have at least coordinates plus one value column.")
  value <- covs[[ncol(covs)]]
  if (!is.numeric(value)) stop("The last column of 'covs' must be numeric (values to interpolate).")
  if (na.rm) {
    keep <- !is.na(value)
    covs  <- covs[keep, , drop = FALSE]
    value <- value[keep]
  }
  if (nrow(covs) == 0) stop("No observations available in 'covs' after removing NAs.")

  if (is.null(d)) {
    P <- .select_coords(points, coord_names = coord_names, require_all = require_all)
    C <- .select_coords(covs,   coord_names = colnames(P), require_all = TRUE)
    d_use <- ncol(P)
  } else {
    P <- .select_coords(points, d = d, coord_names = coord_names, require_all = require_all)
    C <- .select_coords(covs,   d = d, coord_names = coord_names, require_all = require_all)
    d_use <- d
  }

  if (!is.null(colnames(P)) && "t" %in% colnames(P) && !is.numeric(P[, "t"])) P[, "t"] <- as.numeric(P[, "t"])
  if (!is.null(colnames(C)) && "t" %in% colnames(C) && !is.numeric(C[, "t"])) C[, "t"] <- as.numeric(C[, "t"])

  if (rescale01) {
    for (j in seq_len(ncol(P))) {
      combo <- c(P[, j], C[, j])
      combo_sc <- scale_to_range(combo, 0, 1)
      nP <- nrow(P)
      P[, j] <- combo_sc[seq_len(nP)]
      C[, j] <- combo_sc[(nP + 1L):length(combo_sc)]
    }
  }

  storage.mode(P) <- "double"
  storage.mode(C) <- "double"

  D <- distN_fast(P, C, d = d_use, coord_names = coord_names, require_all = FALSE)

  res <- numeric(nrow(P))
  if (exact_on_zero) {
    zero_mask <- (D == 0)
    rows_zero <- which(rowSums(zero_mask) > 0)
    if (length(rows_zero)) {
      for (i in rows_zero) res[i] <- mean(value[zero_mask[i, ]], na.rm = TRUE)
    }
    rows_rest <- setdiff(seq_len(nrow(P)), rows_zero)
  } else {
    rows_rest <- seq_len(nrow(P))
  }
  if (length(rows_rest)) {
    D_sub <- D[rows_rest, , drop = FALSE]
    D_sub[D_sub == 0] <- eps
    W <- 1 / (D_sub^p)
    denom <- rowSums(W)
    numer <- W %*% matrix(value, ncol = 1)
    res[rows_rest] <- as.vector(numer / denom)
  }
  res
}

build_poisson_formula <- function(formula,
                                  add_offset = FALSE,
                                  offset_var = "str_Klocal",
                                  data_cols = NULL) {

  if (inherits(formula, "formula")) {
    if (length(formula) %in% c(2L, 3L)) {
      form_full <- update(formula, y_resp ~ .)
    } else {
      stop("Invalid 'formula': unexpected structure.")
    }
  } else if (is.character(formula)) {
    rhs_txt <- trimws(sub("^~", "", formula[1]))
    form_full <- as.formula(paste0("y_resp ~ ", rhs_txt))
  } else {
    stop("'formula' must be a formula or a character string.")
  }


  if (isTRUE(add_offset)) {
    form_full <- update(form_full,
                        as.formula(paste0(". ~ . + offset(", offset_var, ")")))
  }


  if (!is.null(data_cols)) {
    vars <- setdiff(all.vars(form_full), "y_resp")
    miss <- setdiff(vars, data_cols)
    if (length(miss)) {
      stop("Variables not found in data: ", paste(miss, collapse = ", "))
    }
  }

  form_full
}

counting.weights <- function(id, volumes) {
  id <- as.integer(id)
  fid <- factor(id, levels = seq_along(volumes))
  counts <- table(fid)
  w <- volumes[id] / counts[id]
  w <- as.vector(w)
  names(w) <- NULL
  return(w)
}

grid1.index <- function(x, xrange, nx) {
  i <- ceiling(nx * (x - xrange[1]) / diff(xrange))
  i <- pmax.int(1, i)
  i <- pmin.int(i, nx)
  i
}

rstpp.2D <- function (lambda = 500, nsim = 1, verbose = FALSE, par = NULL,
                      minX = 0, maxX = 1, minY = 0, maxY = 1) {
  if (is.numeric(lambda)) {
    par <- log(lambda)
    lambda <- function(x, y, a) exp(a[1])
  }
  if (nsim != 1) pp0 <- list(l = nsim)
  for (i in 1:nsim) {
    if (isTRUE(verbose)) progressreport(i, nsim)
    lam  <- lambda(1, 1, par)
    candn <- rpois(1, lam)
    candx <- runif(candn, minX, maxX)
    candy <- runif(candn, minY, maxY)
    d     <- runif(candn)
    lam2  <- lambda(candx, candy, par)
    lmax  <- max(lam2)
    keep  <- (d < lam2/lmax)
    lon   <- candx[keep]; lat <- candy[keep]
    if (nsim != 1) {
      pp0[[i]] <- spatstat.geom::ppp(x = lon, y = lat,
                                     window = spatstat.geom::owin(xrange = range(lon), yrange = range(lat)))
    } else {
      pp0 <- spatstat.geom::ppp(x = lon, y = lat,
                                window = spatstat.geom::owin(xrange = range(lon), yrange = range(lat)))
    }
  }
  return(pp0)
}

default.ncube.2D <- function(X){
  guess.ngrid <- floor((splancs::npts(X) / 2) ^ (1 / 3))
  max(5, guess.ngrid)
}

grid.index.2D <- function(x, y, xrange, yrange, nx, ny) {
  ix <- grid1.index(x, xrange, nx)
  iy <- grid1.index(y, yrange, ny)
  list(ix = ix, iy = iy, index = as.integer((iy - 1) * nx + ix))
}

FunPois.fit_v2 <- function(X, formula, covs = NULL,
                           marked = FALSE, spatial.cov = FALSE,
                           verbose = TRUE, mult = 4, seed = NULL, ncube = NULL,
                           grid = FALSE, mark.c = FALSE,
                           mark.c.mode = c("idw","oracle"),
                           mark.oracle_fun = NULL,
                           process.type = NULL, type.cov.values = NULL,
                           Klocal = FALSE, interp.p = 81, prob = NULL,
                           alpha = 0.5, z_weight = TRUE,
                           n_stars_col = "n_stars",
                           sphere = TRUE,
                           sphere_center = c(0,0,0),
                           mc_volume_points = 200000L,

                           box3 = NULL
){
  if (!inherits(X, c("ppp","stp","stpm","s3dp","s3dpm")))
    stop("X should be one of classes: 'ppp','stp','stpm','s3dp','s3dpm'")

  if (is.null(process.type)) {
    process.type <- ifelse(inherits(X, c("s3dp","s3dpm")), "s3d", ifelse(inherits(X, "ppp"), "s2d", "st"))
  } else {
    process.type <- match.arg(tolower(process.type), c("s2d","st","s3d"))
  }
  third_name <- ifelse(process.type == "st", "t", ifelse(process.type == "s3d", "z", NA))

  time1 <- Sys.time()
  if(verbose) cat("Starting stppm.WPI (process.type = ", process.type, ")", "\n\n")

  if (!is.numeric(mult) || mult <= 0) stop("'mult' must be a positive numeric value")
  if(!is.null(ncube)){
    if(!is.numeric(ncube) || ncube <= 0) stop("'ncube' must be a positive numeric value")
  }

  X0 <- X




  if (process.type == "s2d") {

    if (inherits(X0, "ppp")) {
      X <- data.frame(x = X0$x, y = X0$y)
    } else if (inherits(X0, c("stp","stpm","s3dp","s3dpm"))) {

      X <- data.frame(x = X0$df[,1], y = X0$df[,2])
    } else {
      stop("For process.type = 's2d', X must be 'ppp' or includes at least x, y in $df")
    }
    nX <- nrow(X)
    x  <- X[,1]
    y  <- X[,2]

    s.region <- splancs::sbox(cbind(x, y), xfrac = 0, yfrac = 0)
    if (verbose) {
      cat("Built 2D spatial window (s.region)\n\n")
    }

  } else {
    X <- X0$df
    nX <- nrow(X)
    x  <- X[,1]; y <- X[,2]; d3 <- X[,3]

    if (process.type == "s3d" && isTRUE(sphere)) {

      R <- max(sqrt((x - sphere_center[1])^2 + (y - sphere_center[2])^2 + (d3 - sphere_center[3])^2), na.rm = TRUE)


      s.region <- matrix(c(-R, -R,  R, -R,  -R, R,  R, R), ncol = 2, byrow = TRUE)
      third.region <- c(-R, R)

      if (verbose) {
        cat("Built SPHERICAL domain for s3d: center = (",
            paste(sphere_center, collapse = ","), "), R = ", R, "\n\n")
      }
    } else {

      if (!is.null(box3)) {
        if (!is.list(box3)) stop("'box3' must be a list like list(xrange=..., yrange=..., zrange=...) (or trange for st).")

        if (!all(c("xrange","yrange") %in% names(box3))) {
          stop("'box3' must contain at least 'xrange' and 'yrange'.")
        }
        xr <- as.numeric(box3$xrange); yr <- as.numeric(box3$yrange)
        if (length(xr) != 2 || length(yr) != 2) stop("'box3$xrange' and 'box3$yrange' must be length-2 numeric.")


        if (process.type == "st") {
          if ("trange" %in% names(box3)) {
            tr <- as.numeric(box3$trange)
          } else if ("zrange" %in% names(box3)) {
            tr <- as.numeric(box3$zrange)
          } else {
            stop("For process.type='st', 'box3' must contain 'trange' (or 'zrange' as fallback).")
          }
          if (length(tr) != 2) stop("'box3$trange' must be length-2 numeric.")
          third.region <- tr
        } else {
          if (!("zrange" %in% names(box3))) stop("For process.type='s3d', 'box3' must contain 'zrange'.")
          zr <- as.numeric(box3$zrange)
          if (length(zr) != 2) stop("'box3$zrange' must be length-2 numeric.")
          third.region <- zr
        }

        s.region <- matrix(c(xr[1], yr[1],
                             xr[2], yr[1],
                             xr[1], yr[2],
                             xr[2], yr[2]), ncol = 2, byrow = TRUE)

        if (verbose) {
          cat("Using FIXED box domain from 'box3': ",
              "x in [", s.region[1,1], ",", s.region[2,1], "], ",
              "y in [", s.region[1,2], ",", s.region[3,2], "], ",
              third_name, " in [", third.region[1], ",", third.region[2], "]\n\n", sep = "")
        }

      } else {

        s.region <- splancs::sbox(cbind(x, y), xfrac = 0, yfrac = 0)
        r3 <- range(d3, na.rm = TRUE)
        w3 <- diff(r3)
        third.region <- c(r3[1] - 0.0 * w3, r3[2] + 0.0 * w3)

        if (verbose) {
          cat("Built spatial window (s.region) and third-dimension range [",
              third.region[1], ",", third.region[2], "]", "\n\n")
        }
      }
    }
  }

  if (verbose) cat("Observed points: nX = ", nX, "\n\n")

  HomLambda <- nX
  rho <- mult * HomLambda




  if (grid) {
    if (verbose) cat("Generating dummy points on a grid...", "\n\n")

    if (process.type == "s2d") {
      ff <- max(1L, floor(rho^(1/2)))
      x0 <- y0 <- seq_len(ff)

      x0 <- scale_to_range(x0, s.region[1, 1], s.region[2, 1])
      y0 <- scale_to_range(y0, s.region[1, 2], s.region[3, 2])

      df0 <- expand.grid(x0, y0)
      colnames(df0) <- c("x","y")

      dummy_points <- df0

    } else {
      ff <- max(1L, floor(rho^(1/3)))

      if (process.type == "s3d" && isTRUE(sphere)) {
        df0 <- grid_in_sphere3(ff = ff, R = R, center = sphere_center)
        colnames(df0) <- c("x","y","z")
        dummy_points <- df0
      } else {

        x0 <- y0 <- d30 <- seq_len(ff)
        x0  <- scale_to_range(x0,  s.region[1, 1], s.region[2, 1])
        y0  <- scale_to_range(y0,  s.region[1, 2], s.region[3, 2])
        d30 <- scale_to_range(d30, third.region[1], third.region[2])
        df0 <- expand.grid(x0, y0, d30)
        colnames(df0) <- c("x","y", third_name)
        dummy_points <- s3dp(cbind(df0$x, df0$y, df0[[third_name]]))$df
      }
    }

    if(verbose) cat("Dummy points (grid): ", nrow(dummy_points), "\n\n")

  } else {
    if(verbose) cat("Generating dummy points at random...", "\n\n")
    if (!is.null(seed)) { set.seed(seed); if(verbose) cat("Random seed = ", as.character(seed), "\n") }

    if (process.type == "s2d") {
      dummy_points <- rstpp.2D(
        lambda = rho, nsim = 1, verbose = FALSE,
        minX = s.region[1,1], maxX = s.region[2,1],
        minY = s.region[1,2], maxY = s.region[3,2]
      )
      dummy_points <- data.frame(x = dummy_points$x, y = dummy_points$y)

    } else {
      if (process.type == "s3d" && isTRUE(sphere)) {
        dummy_points <- runif_sphere3(n = rho, R = R, center = sphere_center)
      } else {
        dummy_points <- r_st_s3dpp(
          lambda = rho, nsim = 1, verbose = FALSE, process.type = process.type,
          minX = s.region[1, 1], maxX = s.region[2, 1],
          minY = s.region[1, 2], maxY = s.region[3, 2],
          minT = third.region[1],  maxT = third.region[2]
        )$points
      }
    }

    if(verbose) cat("Dummy points (random): ", nrow(dummy_points), "\n\n")
  }


  if (process.type == "s2d") {

    quad_p <- rbind(
      as.matrix(X[, 1:2]),
      as.matrix(dummy_points[, 1:2])
    )
    colnames(quad_p) <- c("x", "y")
    if (verbose) cat("Quadrature points assembled (2D): total (data + dummy) = ",
                     nrow(quad_p), "\n\n")

    xx <- quad_p[, 1]
    xy <- quad_p[, 2]


    win <- spatstat.geom::owin(
      xrange = range(xx, na.rm = TRUE),
      yrange = range(xy, na.rm = TRUE)
    )


    if (is.null(ncube)) ncube <- default.ncube.2D(quad_p)
    ncube <- rep.int(ncube, 2)
    nx <- ncube[1]; ny <- ncube[2]

    nxy <- nx * ny


    cubearea <- spatstat.geom::area.owin(win) / nxy
    volumes <- rep.int(cubearea, nxy)


    id <- grid.index.2D(
      xx, xy,
      win$xrange, win$yrange,
      nx, ny
    )$index

    w <- counting.weights(id, volumes)

  } else {

    quad_p <- rbind(
      as.matrix(X[, 1:3]),
      as.matrix(dummy_points[, 1:3])
    )
    colnames(quad_p) <- c("x", "y", third_name)

    if (verbose) cat("Quadrature points assembled (3D): total (data + dummy) = ",
                     nrow(quad_p), "\n\n")

    xx  <- quad_p[, 1]
    xy  <- quad_p[, 2]
    xd3 <- quad_p[, 3]

    if (process.type == "s3d" && isTRUE(sphere)) {
      win <- spatstat.geom::box3(
        xrange = c(-R, R) + sphere_center[1],
        yrange = c(-R, R) + sphere_center[2],
        zrange = c(-R, R) + sphere_center[3]
      )
    } else {

      if (!is.null(box3) && is.list(box3) && all(c("xrange","yrange") %in% names(box3))) {
        xr <- as.numeric(box3$xrange)
        yr <- as.numeric(box3$yrange)
        zr <- if (process.type == "st") {
          if ("trange" %in% names(box3)) as.numeric(box3$trange) else as.numeric(box3$zrange)
        } else {
          as.numeric(box3$zrange)
        }
        win <- spatstat.geom::box3(
          xrange = xr,
          yrange = yr,
          zrange = zr
        )
      } else {
        win <- spatstat.geom::box3(
          xrange = range(xx,  na.rm = TRUE),
          yrange = range(xy,  na.rm = TRUE),
          zrange = range(xd3, na.rm = TRUE)
        )
      }
    }

    if (is.null(ncube)) ncube <- default.ncube(quad_p)
    ncube <- rep.int(ncube, 3)
    nx <- ncube[1]; ny <- ncube[2]; nt <- ncube[3]

    nxyt <- nx * ny * nt

    if (process.type == "s3d" && isTRUE(sphere)) {
      volumes <- sphere_cell_volumes_mc(
        nx = nx, ny = ny, nz = nt,
        R = R,
        center = sphere_center,
        mc_points = mc_volume_points,
        xrange = win$xrange, yrange = win$yrange, zrange = win$zrange
      )
      if (verbose) {
        cat("Cell volumes (cell ∩ sphere) computed by MC.",
            " Sum(volumes)=", signif(sum(volumes),4),
            " True sphere volume=", signif((4/3)*pi*R^3,4), "\n\n")
      }
    } else {
      cubevolume <- spatstat.geom::volume(win) / nxyt
      volumes <- rep.int(cubevolume, nxyt)
    }

    id <- grid.index(
      xx, xy, xd3,
      win$xrange, win$yrange, win$zrange,
      nx, ny, nt
    )$index

    w <- counting.weights(id, volumes)
  }


  ndata  <- nrow(X)
  ndummy <- nrow(dummy_points)

  Wdat <- w[1:ndata]
  Wdum <- w[(ndata + 1):(ndata + ndummy)]


  if (!is.null(prob)) {
    base_mult <- prob
    if (length(base_mult) != ndata) stop("'prob' must have length ndata.")
  } else {
    base_mult <- rep(1, ndata)
  }

  if(z_weight) {
    if (!n_stars_col %in% names(X)) stop("Column '", n_stars_col, "' not found in X.")
    z_data <- base_mult * (X[[n_stars_col]] ^ alpha)
  }

  if (marked == TRUE && spatial.cov == TRUE) {

    if (verbose) cat("Case: marked = TRUE, spatial.cov = TRUE. Building external covariates on quadrature...", "\n\n")


    points.covs <- build_points_covs(
      quad_p          = quad_p,
      covs            = covs,
      formula         = formula,
      type.cov.values = type.cov.values,
      interp.p        = interp.p,
      process.type    = process.type
    )

    coord_cols <- if (process.type == "s2d") 1:2 else 1:3
    if (verbose) cat(
      "External covariates ready: ",
      max(0, ncol(points.covs) - length(coord_cols)),
      "\n\n"
    )

    if (verbose) cat("Building replicated quadrature for categorical marks...", "\n\n")

    marked.process <- replicated.cubature.dummy(
      X           = X,
      formula     = formula,
      dummy_points = dummy_points,
      Wdum        = Wdum,
      Wdat        = Wdat,
      ndata       = ndata,
      ndummy      = ndummy,
      process.type = process.type
    )


    df_dumb <- as.data.frame(marked.process$dumb)
    n_dumb  <- nrow(df_dumb)
    n_comb  <- marked.process$n_comb_levels

    if (verbose) cat(
      "Replicated quadrature: combinations = ",
      n_comb, " replicated dummy rows = ", n_dumb, "\n\n"
    )


    if(z_weight) {
      z <- c(z_data, rep(0, nrow(df_dumb)))
    } else {
      z <- c(rep(1, ndata), rep(0, nrow(df_dumb)))
    }
    w_final <- c(w[1:ndata],   marked.process$Wdumb)
    y_resp  <- z / w_final


    if (process.type == "s2d") {

      dati.modello <- data.frame(
        y_resp = y_resp,
        w      = w_final,
        x      = c(X$x, df_dumb$x),
        y      = c(X$y, df_dumb$y),
        check.names = FALSE
      )
    } else {

      third_name <- if (process.type == "st") "t" else "z"
      dati.modello <- data.frame(
        y_resp = y_resp,
        w      = w_final,
        x      = c(X$x,            df_dumb$x),
        y      = c(X$y,            df_dumb$y),
        third  = c(X[[third_name]], df_dumb$z),
        check.names = FALSE
      )
      names(dati.modello)[names(dati.modello) == "third"] <- third_name
    }


    dati.modello <- cbind(
      dati.modello,
      rbind(marked.process$df_marks, marked.process$total_dummy_marks)
    )


    stopifnot(nrow(dati.modello) == (ndata + ndummy) * n_comb)
    if (verbose) cat(
      "Model frame (with replicated marks) assembled: ",
      nrow(dati.modello), " rows.", "\n\n"
    )


    cov_names_only <- names(points.covs)[-coord_cols]

    if (length(cov_names_only)) {
      if (verbose) cat("Replicating external covariates across mark combinations...", "\n\n")

      dati.interpolati.rep <- as.data.frame(
        matrix(NA_real_, nrow = nrow(dati.modello), ncol = length(cov_names_only)),
        check.names = FALSE
      )
      colnames(dati.interpolati.rep) <- cov_names_only

      stopifnot(nrow(points.covs) == (ndata + ndummy))

      for (nm in cov_names_only) {
        v_data  <- points.covs[1:ndata, nm]
        v_dummy <- points.covs[(ndata+1):(ndata+ndummy), nm]

        dati.interpolati.rep[[nm]] <- c(
          v_data,
          rep(v_dummy, each = n_comb),
          if (n_comb > 1) rep(v_data, each = n_comb - 1) else NULL
        )
      }
      dati.cov.marks <- cbind(dati.modello, dati.interpolati.rep)
    } else {
      dati.cov.marks <- dati.modello
    }

    if (verbose) cat(
      "Design matrix ready (rows: ",
      nrow(dati.cov.marks), ", cols: ", ncol(dati.cov.marks), ")\n\n"
    )
  } else if (marked == FALSE & spatial.cov == TRUE) {

    if (verbose) cat("Case: marked = FALSE, spatial.cov = TRUE. Building external covariates on quadrature...", "\n\n")


    points.covs <- build_points_covs(
      quad_p          = quad_p,
      covs            = covs,
      formula         = formula,
      type.cov.values = type.cov.values,
      interp.p        = interp.p,
      process.type    = process.type
    )


    if(z_weight) {
      z <- c(z_data, rep(0, ndummy))
    } else {
      z <- c(rep(1, ndata), rep(0, ndummy))
    }
    y_resp  <- z / w


    dati.cov.marks <- cbind(
      y_resp = y_resp,
      w      = w,
      points.covs
    )

    if (verbose) cat(
      "Design matrix ready (rows: ",
      nrow(dati.cov.marks), ", cols: ", ncol(dati.cov.marks), ")\n\n"
    )

  } else if (marked == FALSE & spatial.cov == FALSE) {

    if (verbose) cat("Case: marked = FALSE, spatial.cov = FALSE. Using only coordinates.\n")


    if (process.type == "s2d") {
      colnames(quad_p) <- c("x", "y")
    } else {
      third_name <- if (process.type == "st") "t" else "z"
      colnames(quad_p) <- c("x", "y", third_name)
    }


    if(z_weight) {
      z <- c(z_data, rep(0, ndummy))
    } else {
      z <- c(rep(1, ndata), rep(0, ndummy))
    }
    y_resp <- z / w


    dati.cov.marks <- cbind(
      y_resp = y_resp,
      w      = w,
      quad_p
    )

    if (verbose) cat(
      "Design matrix ready (rows: ",
      nrow(dati.cov.marks), ", cols: ", ncol(dati.cov.marks), ")\n\n"
    )

  } else if (marked == TRUE && spatial.cov == FALSE) {

    if (verbose) cat("Case: marked = TRUE, spatial.cov = FALSE. Building replicated quadrature for categorical marks...\n\n")

    marked.process <- replicated.cubature.dummy(
      X            = X,
      formula      = formula,
      dummy_points = dummy_points,
      Wdum         = Wdum,
      Wdat         = Wdat,
      ndata        = ndata,
      ndummy       = ndummy,
      process.type = process.type
    )


    df_dumb <- as.data.frame(marked.process$dumb)
    n_dumb  <- nrow(df_dumb)


    if(z_weight) {
      z <- c(z_data, rep(0, nrow(df_dumb)))
    } else {
      z <- c(rep(1, ndata), rep(0, nrow(df_dumb)))
    }
    w_final <- c(w[1:ndata],   marked.process$Wdumb)
    y_resp  <- z / w_final


    if (process.type == "s2d") {

      dati.modello <- data.frame(
        y_resp = y_resp,
        w      = w_final,
        x      = c(X$x, df_dumb$x),
        y      = c(X$y, df_dumb$y),
        check.names = FALSE
      )
    } else {

      third_name <- if (process.type == "st") "t" else "z"
      dati.modello <- data.frame(
        y_resp = y_resp,
        w      = w_final,
        x      = c(X$x,            df_dumb$x),
        y      = c(X$y,            df_dumb$y),
        third  = c(X[[third_name]], df_dumb$z),
        check.names = FALSE
      )
      names(dati.modello)[names(dati.modello) == "third"] <- third_name
    }


    if (ncol(marked.process$df_marks) > 0) {
      dati.modello <- cbind(
        dati.modello,
        rbind(marked.process$df_marks, marked.process$total_dummy_marks)
      )
    }

    dati.cov.marks <- dati.modello

    if (verbose) cat(
      "Design matrix ready (rows: ",
      nrow(dati.cov.marks), ", cols: ", ncol(dati.cov.marks), ")\n\n"
    )
  }




  if (mark.c) {

    if (is.null(mark.c.mode)) {
      mark.c.mode <- "idw"
    } else {
      mark.c.mode <- match.arg(mark.c.mode, c("idw","oracle"))
    }

    if (isTRUE(verbose)) cat("Handling continuous marks (mode = ", mark.c.mode, ")...\n\n", sep = "")

    dati.cov.marks <- as.data.frame(dati.cov.marks)


    if (process.type == "s2d") {
      coord_cols <- c("x", "y")
    } else {
      third_name <- if (process.type == "st") "t" else "z"
      coord_cols <- c("x", "y", third_name)
    }


    form_vars       <- all.vars(formula)
    mark_candidates <- setdiff(colnames(X), coord_cols)
    numeric_marks   <- if (length(mark_candidates)) {
      mark_candidates[vapply(X[, mark_candidates, drop = FALSE], is.numeric, logical(1))]
    } else character(0)

    cont_mark_names <- intersect(form_vars, numeric_marks)

    if (!length(cont_mark_names)) {
      if (isTRUE(verbose)) cat("No continuous marks found in formula; skipping.\n\n")
    } else {


      n_base  <- ndata + ndummy
      n_total <- nrow(dati.cov.marks)
      if (n_total %% n_base != 0)
        stop("Quadrature replication mismatch: nrow(dati.cov.marks) is not a multiple of (ndata + ndummy).")
      n_comb <- as.integer(n_total / n_base)

      if (!all(coord_cols %in% names(dummy_points))) {
        miss <- setdiff(coord_cols, names(dummy_points))
        stop("dummy_points missing coordinate column(s): ", paste(miss, collapse = ", "))
      }
      if (!all(coord_cols %in% names(X))) {
        miss <- setdiff(coord_cols, names(X))
        stop("X missing coordinate column(s): ", paste(miss, collapse = ", "))
      }


      M_obs <- as.matrix(X[, cont_mark_names, drop = FALSE])
      storage.mode(M_obs) <- "double"
      if (any(!is.finite(M_obs))) {
        bad <- which(!is.finite(M_obs), arr.ind = TRUE)
        stop("Non-finite values in continuous marks. Example col(s): ",
             paste(unique(colnames(M_obs)[bad[,2]]), collapse = ", "))
      }
      K <- ncol(M_obs)

      coords_data  <- as.matrix(X[, coord_cols, drop = FALSE]); storage.mode(coords_data) <- "double"
      coords_dummy <- as.matrix(dummy_points[, coord_cols, drop = FALSE]); storage.mode(coords_dummy) <- "double"


      M_dummy <- matrix(NA_real_, nrow = ndummy, ncol = K)
      colnames(M_dummy) <- cont_mark_names

      if (identical(mark.c.mode, "idw")) {

        for (j in seq_along(cont_mark_names)) {
          nm <- cont_mark_names[j]
          covs_mark <- X[, c(coord_cols, nm), drop = FALSE]
          pred_dummy <- idw_interp_fast(
            points = dummy_points[, coord_cols, drop = FALSE],
            covs   = covs_mark,
            p      = interp.p,
            d      = length(coord_cols)
          )
          M_dummy[, j] <- pred_dummy
        }

      } else if (identical(mark.c.mode, "oracle")) {

        if (is.null(mark.oracle_fun) || !is.function(mark.oracle_fun)) {
          stop("mark.c.mode='oracle' requires 'mark.oracle_fun' (a function) that returns an ndummy x K matrix of marks.")
        }
        M_or <- mark.oracle_fun(dummy_points[, coord_cols, drop = FALSE])
        M_or <- as.matrix(M_or)
        if (nrow(M_or) != ndummy) stop("oracle marks: nrow does not match ndummy.")
        if (ncol(M_or) != K) stop("oracle marks: ncol does not match number of continuous marks in formula.")
        colnames(M_or) <- cont_mark_names
        M_dummy[,] <- M_or
      }

      M_dumdum <- M_dummy[rep(seq_len(ndummy), each = n_comb), , drop = FALSE]
      if (n_comb > 1L) {
        M_dumdat <- M_obs[rep(seq_len(ndata), each = n_comb - 1L), , drop = FALSE]
        M_all <- rbind(M_obs, M_dumdum, M_dumdat)
      } else {
        M_all <- rbind(M_obs, M_dumdum)
      }

      stopifnot(nrow(M_all) == n_total)
      dati.cov.marks <- cbind(dati.cov.marks, as.data.frame(M_all, check.names = FALSE))

      if (isTRUE(verbose)) cat("Continuous marks appended (", mark.c.mode,
                               "). cols added: ", K, "\n\n", sep = "")
    }
  }

  if (isTRUE(Klocal)) {
    if(verbose) cat("Computing local K-function offset...", "\n\n")

    third_name <- if (process.type == "st") "t" else "z"

    if (process.type == "st") {
      obj_st <- stp(as.matrix(X[, c("x","y","t")]))
      kout   <- Khat_st(obj_st, correction = "translate", Klocal = TRUE)

      Kloc  <- kout$Khat
      Ktheo <- kout$Ktheo

      npt <- dim(Kloc)[3]
      w_obs <- vapply(seq_len(npt), function(i) {
        sum((Kloc[ , , i] - Ktheo) / Ktheo, na.rm = TRUE)
      }, numeric(1))

    } else {
      obj_s3d <- s3dp(as.matrix(X[, c("x","y","z")]))
      kout    <- Khat_s3d(obj_s3d, correction = "translate", Klocal = TRUE)

      Kloc  <- kout$Khat
      Ktheo <- kout$Ktheo

      rel   <- sweep(Kloc, 1, Ktheo, FUN = function(a, b) (a - b)/b)
      w_obs <- colSums(rel, na.rm = TRUE)
    }

    X_Klocal <- X[, c("x","y", third_name), drop = FALSE]
    X_Klocal$w_klocal <- w_obs

    kpred <- idw_interp_fast(
      points = dummy_points[, c("x","y", third_name), drop = FALSE],
      covs   = X_Klocal[, c("x","y", third_name, "w_klocal"), drop = FALSE],
      p      = interp.p,
      d      = 3
    )

    if (isTRUE(marked) && exists("marked.process") && !is.null(marked.process$n_comb_levels)) {
      m <- marked.process$n_comb_levels
      offset_vec <- c(
        w_obs,
        rep(kpred, m),
        rep(w_obs,  m - 1L)
      )
    } else {
      offset_vec <- c(w_obs, kpred)
    }

    offset_vec <- offset_vec - mean(offset_vec, na.rm = TRUE)
    dati.cov.marks$str_Klocal <- offset_vec
    if(verbose) cat("Local K offset computed and appended.", "\n\n")
  }

  if(verbose) print(summary(dati.cov.marks))
  dati.cov.marks <- as.data.frame(dati.cov.marks)

  add_off <- isTRUE(Klocal) && "str_Klocal" %in% names(dati.cov.marks)
  form_full <- build_poisson_formula(formula, add_offset = add_off)

  suppressWarnings(
    mod_global <- try(
      mgcv::gam(
        formula  = form_full,
        family   = poisson(link = "log"),
        data     = dati.cov.marks,
        weights  = w
      ),
      silent = TRUE
    )
  )

  if (inherits(mod_global, "try-error")) {
    stop("Model fit failed: ", as.character(mod_global))
  }
  if(verbose) cat("Model fit complete.", "\n\n")


  pred_data <- predict(mod_global, newdata = dati.cov.marks[1:ndata, , drop = FALSE], type = "response")


  idx_dummy_start <- ndata + 1L
  idx_dummy_end   <- nrow(dati.cov.marks)
  pred_dummy_all  <- predict(mod_global,
                             newdata = dati.cov.marks[idx_dummy_start:idx_dummy_end, , drop = FALSE],
                             type = "response")


  pred_dummy_mean <- pred_dummy_all

  res_global  <- coef(mod_global)

  time2 <- Sys.time()
  elapsed <- round(as.numeric(difftime(time2, time1, units = "sec")), 3)
  if(verbose) cat("Done in ", elapsed, " secs.", "\n\n")

  list.obj <- list(IntCoefs = coef(mod_global), X = X0, nX = ndata, formula = formula, mod_global = mod_global,

                   newdata = dati.cov.marks[1:ndata, , drop = FALSE],
                   dummy.data = dati.cov.marks[(ndata+1):nrow(dati.cov.marks), , drop = FALSE],

                   l = as.vector(pred_data), l_dummy = as.vector(pred_dummy_mean),

                   l_dummy_all = as.vector(pred_dummy_all),

                   w_all = w, y_resp = y_resp, ncube = ncube,

                   alpha = alpha, sphere = sphere, sphere_center = sphere_center,
                   sphere_radius = if (exists("R")) R else NA_real_, mult = mult, seed = seed,
                   mark_c = mark.c, mark_c_mode = mark.c.mode,
                   box3 = box3,
                   time = paste0(elapsed, " sec"))

  if(process.type == "st") {
    class(list.obj) <- "stppm"
  } else {
    class(list.obj) <- "s3dppm"
  }
  return(list.obj)
}


# Minimal helper functions required by the cubature estimator.
clamp_int <- function(x, lower, upper) {
  as.integer(max(lower, min(upper, as.integer(x))))
}

rmse_vec <- function(x, y) {
  sqrt(mean((x - y)^2))
}

scale_to_range <- function(x, lower = 0, upper = 1) {
  x <- as.numeric(x)
  (x - min(x)) / (max(x) - min(x) + 1e-12) * (upper - lower) + lower
}

r_st_s3dpp <- function(lambda, nsim = 1, verbose = FALSE, process.type = "s3d",
                        minX, maxX, minY, maxY, minT, maxT) {
  n <- as.integer(round(lambda))
  points <- data.frame(
    x = stats::runif(n, minX, maxX),
    y = stats::runif(n, minY, maxY),
    z = stats::runif(n, minT, maxT)
  )
  list(points = points)
}

default.ncube <- function(quad_p, target_ppc = 8) {
  n <- nrow(quad_p)
  max(3L, as.integer(floor((n / target_ppc)^(1 / 3))))
}

grid.index <- function(x, y, z, xrange, yrange, zrange, nx, ny, nz) {
  sx <- (x - xrange[1]) / (xrange[2] - xrange[1] + 1e-12)
  sy <- (y - yrange[1]) / (yrange[2] - yrange[1] + 1e-12)
  sz <- (z - zrange[1]) / (zrange[2] - zrange[1] + 1e-12)
  ix <- pmin(nx, pmax(1L, floor(sx * nx) + 1L))
  iy <- pmin(ny, pmax(1L, floor(sy * ny) + 1L))
  iz <- pmin(nz, pmax(1L, floor(sz * nz) + 1L))
  list(index = as.integer(ix + (iy - 1L) * nx + (iz - 1L) * nx * ny))
}

s3dpm <- function(df, names = NULL) {
  object <- list(df = as.data.frame(df), names = names)
  class(object) <- "s3dpm"
  object
}

# -----------------------------------------------------------------------------
# Functional data-generating mechanism and Monte Carlo experiments
# -----------------------------------------------------------------------------

trapz_weights <- function(u) {
  dx <- diff(u)
  c(dx[1], (dx[-1] + dx[-length(dx)]) / 2, dx[length(dx)])
}

w_center <- function(x, w) x - sum(x * w) / sum(w)

orthonormalize_weighted <- function(Phi, w) {
  sw <- sqrt(w)
  A  <- Phi * sw
  Q  <- qr.Q(qr(A))
  Q / sw
}

make_true_phi_fourier <- function(u, K0 = 6) {
  stopifnot(K0 %% 2 == 0)
  M <- length(u)
  Khalf <- K0 / 2
  Phi_raw <- matrix(NA_real_, nrow = M, ncol = K0)
  col <- 1
  for (k in seq_len(Khalf)) {
    Phi_raw[, col] <- sin(2*pi*k*u); col <- col + 1
    Phi_raw[, col] <- cos(2*pi*k*u); col <- col + 1
  }
  w <- trapz_weights(u)
  Phi <- orthonormalize_weighted(Phi_raw, w)
  colnames(Phi) <- paste0("phi", seq_len(K0))
  Phi
}

make_beta_scenarios <- function(u, beta_sd_target = 0.30) {
  beta1 <- u
  beta_single <- exp(-0.5 * ((u - 0.55) / 0.12)^2)
  beta2 <- exp(-0.5 * ((u - 0.30) / 0.08)^2) + 0.8 * exp(-0.5 * ((u - 0.72) / 0.10)^2)
  beta3 <- tanh((u - 0.50) / 0.08)
  beta4 <- sin(2 * pi * u)
  beta5 <- exp(-0.5 * ((u - 0.60) / 0.04)^2)

  betas <- list(
    monotone  = beta1,
    single    = beta_single,
    double    = beta2,
    signflip  = beta3,
    oscillate = beta4,
    narrow    = beta5
  )

  betas <- lapply(betas, function(b) {
    b <- b - mean(b)
    b <- b / sd(b) * beta_sd_target
    b
  })
  betas
}

make_rff_field <- function(K0, range, nugget, L = 128L, seed = NULL, trunc_xi = 4) {
  if (!is.null(seed)) set.seed(seed)
  d <- 3
  omega <- matrix(rnorm(L * d, sd = 1 / range), nrow = L, ncol = d)
  b <- runif(L, 0, 2*pi)

  nug2 <- nugget^2
  sig2 <- max(1 - nug2, 1e-8)
  W <- matrix(rnorm(L * K0), nrow = L, ncol = K0) * sqrt(sig2)

  eval_xi <- function(coords) {
    coords <- as.matrix(coords)
    PhiL <- coords %*% t(omega)
    PhiL <- sweep(PhiL, 2, b, `+`)
    PhiL <- sqrt(2 / L) * cos(PhiL)
    xi_sp <- PhiL %*% W
    xi_ng <- matrix(rnorm(nrow(coords) * K0), nrow(coords), K0) * sqrt(nug2)
    xi <- xi_sp + xi_ng
    xi <- pmax(pmin(xi, trunc_xi), -trunc_xi)
    xi
  }

  list(eval_xi = eval_xi, range = range, nugget = nugget, L = L, trunc_xi = trunc_xi)
}

simulate_X_curves <- function(coords, field, Phi_true, sigma_eps = 0.05, trunc_eps = 4) {
  coords <- as.matrix(coords)
  xi <- field$eval_xi(coords)
  Xsig <- xi %*% t(Phi_true)
  n <- nrow(Xsig); M <- ncol(Xsig)
  eps <- matrix(rnorm(n * M, sd = sigma_eps), n, M)
  if (is.finite(trunc_eps) && trunc_eps > 0) {
    lo <- -trunc_eps * sigma_eps
    hi <-  trunc_eps * sigma_eps
    eps <- pmax(pmin(eps, hi), lo)
  }
  Xsig + eps
}

calibrate_a0_functional <- function(beta_w, Phi_true, target_N,
                                    beta_xyz = c(0,0,0),
                                    K0 = 6,
                                    sigma_eps = 0.05, trunc_eps = 4,
                                    trunc_xi = 4,
                                    mc_n = 20000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  M <- length(beta_w)
  beta_xyz <- as.numeric(beta_xyz)
  if (length(beta_xyz) != 3) stop("beta_xyz must be length 3 (betax,betay,betaz).")


  xi <- matrix(rnorm(mc_n * K0), mc_n, K0)
  xi <- pmax(pmin(xi, trunc_xi), -trunc_xi)


  Xsig <- xi %*% t(Phi_true)

  eps <- matrix(rnorm(mc_n * M, sd = sigma_eps), mc_n, M)
  if (is.finite(trunc_eps) && trunc_eps > 0) {
    lo <- -trunc_eps * sigma_eps
    hi <-  trunc_eps * sigma_eps
    eps <- pmax(pmin(eps, hi), lo)
  }
  X <- Xsig + eps

  integ <- as.vector(X %*% beta_w)


  xx <- runif(mc_n); yy <- runif(mc_n); zz <- runif(mc_n)
  lin <- beta_xyz[1]*xx + beta_xyz[2]*yy + beta_xyz[3]*zz

  Eexp <- mean(exp(integ + lin))

  log(target_N) - log(Eexp)
}

lmax_bound_functional <- function(a0, beta_w, Phi_true,
                                  beta_xyz = c(0,0,0),
                                  trunc_xi = 4,
                                  sigma_eps = 0.05, trunc_eps = 4) {
  beta_xyz <- as.numeric(beta_xyz)
  if (length(beta_xyz) != 3) stop("beta_xyz must be length 3.")

  max_abs_signal <- rowSums(abs(Phi_true) * trunc_xi)
  max_abs_eps <- trunc_eps * sigma_eps
  max_abs_X <- max_abs_signal + max_abs_eps
  integral_max <- sum(abs(beta_w) * max_abs_X)

  max_abs_lin <- sum(abs(beta_xyz))

  exp(a0 + integral_max + max_abs_lin)
}

simulate_one_pattern_box3 <- function(beta_w, a0, lmax,
                                      field, Phi_true,
                                      beta_xyz = c(0,0,0),
                                      sigma_eps = 0.05, trunc_eps = 4,
                                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  candn <- rpois(1, lmax)
  if (candn == 0L) {
    return(list(coords = matrix(numeric(0), 0, 3),
                Xcurves = matrix(numeric(0), 0, ncol(Phi_true))))
  }

  coords <- cbind(runif(candn), runif(candn), runif(candn))
  colnames(coords) <- c("x","y","z")

  Xcand <- simulate_X_curves(coords, field, Phi_true, sigma_eps = sigma_eps, trunc_eps = trunc_eps)

  integ <- as.vector(Xcand %*% beta_w)
  lin   <- beta_xyz[1]*coords[,1] + beta_xyz[2]*coords[,2] + beta_xyz[3]*coords[,3]
  lam   <- exp(a0 + integ + lin)

  keep <- runif(candn) < (lam / lmax)
  if (!any(keep)) {
    return(list(coords = matrix(numeric(0), 0, 3),
                Xcurves = matrix(numeric(0), 0, ncol(Phi_true))))
  }

  list(coords = coords[keep, , drop = FALSE],
       Xcurves = Xcand[keep, , drop = FALSE])
}

fpca_grid <- function(F, grid, K = 5) {
  stopifnot(is.matrix(F), length(grid) == ncol(F), K >= 1)

  ids <- rownames(F); if (is.null(ids)) ids <- as.character(seq_len(nrow(F)))

  dx <- diff(grid)
  w  <- c(dx[1], (dx[-1] + dx[-length(dx)]) / 2, dx[length(dx)])
  sw <- sqrt(w)

  mu <- colMeans(F, na.rm = TRUE)
  Fc <- sweep(F, 2, mu, "-")
  Fc[!is.finite(Fc)] <- 0

  Xw <- sweep(Fc, 2, sw, "*")

  totalSS <- sum(Xw^2)
  sv <- svd(Xw, nu = K, nv = K)

  scores <- sv$u %*% diag(sv$d[1:K], K)
  rownames(scores) <- ids
  colnames(scores) <- paste0("age_pc", seq_len(K))

  phi_w <- sv$v
  phi   <- sweep(phi_w, 1, sw, "/")
  colnames(phi) <- paste0("age_pc", seq_len(K))

  pve <- (sv$d[1:K]^2) / totalSS
  pve_cum <- cumsum(pve)

  list(mu = mu, phi = phi, scores = scores, grid = grid, w = w,
       d = sv$d[1:K], totalSS = totalSS, pve = pve, pve_cum = pve_cum)
}

fpca_project_scores <- function(F_new, fp) {

  stopifnot(is.matrix(F_new), length(fp$mu) == ncol(F_new))
  Fc <- sweep(F_new, 2, fp$mu, "-")
  S  <- Fc %*% (fp$phi * fp$w)
  colnames(S) <- colnames(fp$phi)
  S
}

beta_hat_from_gamma_fpca <- function(gamma, fp, w_u) {
  b <- as.numeric(fp$phi %*% gamma)
  w_center(b, w_u)
}

build_s3dpm_fpca <- function(coords, scores) {
  coords <- as.matrix(coords)
  scores <- as.matrix(scores)
  df <- data.frame(
    x = coords[,1], y = coords[,2], z = coords[,3],
    n_stars = rep(1, nrow(coords)),
    as.data.frame(scores),
    check.names = FALSE
  )
  mark_cols <- colnames(scores)
  X_s3d <- s3dpm(df = df, names = c("n_stars", mark_cols))
  X_s3d
}

fit_funpois_fpca <- function(X_s3d,
                             rhs_terms,
                             interp.p,
                             mult_fit,
                             seed_fit,
                             mark_mode = c("idw","oracle"),
                             mark_oracle_fun = NULL,
                             box3 = list(xrange=c(0,1), yrange=c(0,1), zrange=c(0,1))) {

  mark_mode <- match.arg(mark_mode)

  if (!exists("FunPois.fit_v2")) stop("FunPois.fit_v2 not found in environment.")

  form <- as.formula(paste0("y_resp ~ ", rhs_terms))

  FunPois.fit_v2(
    X = X_s3d,
    formula = form,
    covs = NULL,
    marked = FALSE,
    spatial.cov = FALSE,
    verbose = FALSE,
    mult = mult_fit,
    seed = seed_fit,
    ncube = NULL,
    grid = FALSE,
    mark.c = TRUE,
    mark.c.mode = mark_mode,
    mark.oracle_fun = mark_oracle_fun,
    process.type = "s3d",
    type.cov.values = NULL,
    Klocal = FALSE,
    interp.p = interp.p,
    prob = NULL,
    alpha = 0.5,
    z_weight = FALSE,
    sphere = FALSE,
    box3 = box3
  )
}

make_oracle_fun_scores <- function(field, Phi_true, u_grid, fp, sigma_eps, trunc_eps) {
  force(field); force(Phi_true); force(u_grid); force(fp); force(sigma_eps); force(trunc_eps)
  function(dummy_coords_df) {
    coords <- as.matrix(dummy_coords_df)

    Xdum <- simulate_X_curves(coords, field, Phi_true, sigma_eps = sigma_eps, trunc_eps = trunc_eps)
    Sdum <- fpca_project_scores(Xdum, fp)
    Sdum
  }
}

run_one_rep_fpca <- function(beta_u, u_grid, Phi_true,
                             range, nugget, target_N,
                             beta_xyz = c(0,0,0),

                             K_fpca = 5,
                             interp.p = 4,
                             mult_fit = 30,

                             fpca_basis_mode = c("event","background"),
                             bg_mult = 50, N_bg_min = 3000, N_bg_max = 10000,

                             mark_mode = c("idw","oracle"),

                             K0 = 6,
                             sigma_eps = 0.05, trunc_eps = 4,
                             trunc_xi = 4, L_rff = 128L,

                             n_eval = 5000,

                             mc_n = 20000L,
                             seed = NULL) {

  fpca_basis_mode <- match.arg(fpca_basis_mode)
  mark_mode <- match.arg(mark_mode)

  if (!is.null(seed)) set.seed(seed)


  w_u <- trapz_weights(u_grid)
  beta_u <- w_center(beta_u, w_u)
  beta_w <- beta_u * w_u


  a0 <- calibrate_a0_functional(
    beta_w = beta_w,
    Phi_true = Phi_true,
    target_N = target_N,
    beta_xyz = beta_xyz,
    K0 = K0,
    sigma_eps = sigma_eps, trunc_eps = trunc_eps,
    trunc_xi = trunc_xi,
    mc_n = mc_n,
    seed = if (!is.null(seed)) seed + 11L else NULL
  )

  lmax <- lmax_bound_functional(
    a0 = a0, beta_w = beta_w, Phi_true = Phi_true,
    beta_xyz = beta_xyz,
    trunc_xi = trunc_xi,
    sigma_eps = sigma_eps, trunc_eps = trunc_eps
  )


  field <- make_rff_field(
    K0 = K0, range = range, nugget = nugget,
    L = L_rff,
    seed = if (!is.null(seed)) seed + 101L else NULL,
    trunc_xi = trunc_xi
  )


  sim <- simulate_one_pattern_box3(
    beta_w = beta_w, a0 = a0, lmax = lmax,
    field = field, Phi_true = Phi_true,
    beta_xyz = beta_xyz,
    sigma_eps = sigma_eps, trunc_eps = trunc_eps,
    seed = if (!is.null(seed)) seed + 1001L else NULL
  )

  coords  <- sim$coords
  Xevent  <- sim$Xcurves
  n <- nrow(coords)


  if (n < max(15, K_fpca + 3)) {
    return(list(
      ok = FALSE, n = n,
      beta_hat = rep(NA_real_, length(u_grid)),
      gamma_hat = rep(NA_real_, K_fpca),
      coef_xyz = c(betax = NA_real_, betay = NA_real_, betaz = NA_real_),
      rmse_beta = NA_real_, ise_beta = NA_real_,
      logrmse_lambda = NA_real_, rmse_lambda = NA_real_,
      mise_lambda = NA_real_,
      I_true = NA_real_, I_hat = NA_real_, relerr_I = NA_real_,
      fpca_basis_mode = fpca_basis_mode,
      mark_mode = mark_mode,
      N_bg = NA_integer_,

      pve = NA_real_,
      pve_cum = NA_real_,
      pve_cum_last = NA_real_
    ))
  }


  N_bg_used <- NA_integer_
  if (fpca_basis_mode == "event") {
    rownames(Xevent) <- as.character(seq_len(n))
    fp <- fpca_grid(Xevent, u_grid, K = K_fpca)
  } else {
    N_bg_used <- clamp_int(bg_mult * n, N_bg_min, N_bg_max)
    coords_bg <- cbind(runif(N_bg_used), runif(N_bg_used), runif(N_bg_used))
    Xbg <- simulate_X_curves(coords_bg, field, Phi_true, sigma_eps = sigma_eps, trunc_eps = trunc_eps)
    rownames(Xbg) <- as.character(seq_len(N_bg_used))
    fp <- fpca_grid(Xbg, u_grid, K = K_fpca)
  }


  pve     <- fp$pve
  pve_cum <- fp$pve_cum
  pve_cum_last <- if (length(pve_cum)) pve_cum[length(pve_cum)] else NA_real_


  scores_event <- fpca_project_scores(Xevent, fp)
  colnames(scores_event) <- paste0("age_pc", seq_len(K_fpca))


  X_s3d <- build_s3dpm_fpca(coords, scores_event)
  mark_cols <- paste0("age_pc", seq_len(K_fpca))


  rhs <- paste(mark_cols, collapse = " + ")
  if (!isTRUE(all(beta_xyz == 0))) {
    if (beta_xyz[1] != 0) rhs <- paste("x +", rhs)
    if (beta_xyz[2] != 0) rhs <- paste("y +", rhs)
    if (beta_xyz[3] != 0) rhs <- paste("z +", rhs)
  }


  oracle_fun <- NULL
  if (mark_mode == "oracle") {
    oracle_fun <- make_oracle_fun_scores(field, Phi_true, u_grid, fp, sigma_eps, trunc_eps)
  }


  fit <- try(
    fit_funpois_fpca(
      X_s3d = X_s3d,
      rhs_terms = rhs,
      interp.p = interp.p,
      mult_fit = mult_fit,
      seed_fit = if (!is.null(seed)) seed + 2001L else NULL,
      mark_mode = mark_mode,
      mark_oracle_fun = oracle_fun,
      box3 = list(xrange=c(0,1), yrange=c(0,1), zrange=c(0,1))
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    return(list(
      ok = FALSE, n = n,
      beta_hat = rep(NA_real_, length(u_grid)),
      gamma_hat = rep(NA_real_, K_fpca),
      coef_xyz = c(betax = NA_real_, betay = NA_real_, betaz = NA_real_),
      rmse_beta = NA_real_, ise_beta = NA_real_,
      logrmse_lambda = NA_real_, rmse_lambda = NA_real_,
      mise_lambda = NA_real_,
      I_true = NA_real_, I_hat = NA_real_, relerr_I = NA_real_,
      fpca_basis_mode = fpca_basis_mode,
      mark_mode = mark_mode,
      N_bg = N_bg_used,

      pve = pve,
      pve_cum = pve_cum,
      pve_cum_last = pve_cum_last
    ))
  }

  cc <- coef(fit$mod_global)


  if (!all(mark_cols %in% names(cc))) {
    return(list(
      ok = FALSE, n = n,
      beta_hat = rep(NA_real_, length(u_grid)),
      gamma_hat = rep(NA_real_, K_fpca),
      coef_xyz = c(betax = NA_real_, betay = NA_real_, betaz = NA_real_),
      rmse_beta = NA_real_, ise_beta = NA_real_,
      logrmse_lambda = NA_real_, rmse_lambda = NA_real_,
      mise_lambda = NA_real_,
      I_true = NA_real_, I_hat = NA_real_, relerr_I = NA_real_,
      fpca_basis_mode = fpca_basis_mode,
      mark_mode = mark_mode,
      N_bg = N_bg_used,

      pve = pve,
      pve_cum = pve_cum,
      pve_cum_last = pve_cum_last
    ))
  }
  gamma_hat <- as.numeric(cc[mark_cols])


  bx_hat <- if ("x" %in% names(cc)) as.numeric(cc["x"]) else 0
  by_hat <- if ("y" %in% names(cc)) as.numeric(cc["y"]) else 0
  bz_hat <- if ("z" %in% names(cc)) as.numeric(cc["z"]) else 0


  beta_hat <- beta_hat_from_gamma_fpca(gamma_hat, fp, w_u)


  rmse_b <- sqrt(mean((beta_hat - beta_u)^2))
  ise_b  <- sum(((beta_hat - beta_u)^2) * w_u) / sum(w_u)


  coords_eval <- cbind(runif(n_eval), runif(n_eval), runif(n_eval))
  Xeval <- simulate_X_curves(coords_eval, field, Phi_true, sigma_eps = sigma_eps, trunc_eps = trunc_eps)
  integ_eval <- as.vector(Xeval %*% beta_w)
  lin_eval <- beta_xyz[1]*coords_eval[,1] + beta_xyz[2]*coords_eval[,2] + beta_xyz[3]*coords_eval[,3]
  lam_true <- exp(a0 + integ_eval + lin_eval)

  scores_eval <- fpca_project_scores(Xeval, fp)
  colnames(scores_eval) <- mark_cols
  df_eval <- data.frame(
    x = coords_eval[,1], y = coords_eval[,2], z = coords_eval[,3],
    as.data.frame(scores_eval),
    check.names = FALSE
  )
  lam_hat <- as.numeric(predict(fit$mod_global, newdata = df_eval, type = "response"))

  logrmse_l <- sqrt(mean((log(lam_hat) - log(lam_true))^2))
  rmse_l    <- sqrt(mean((lam_hat - lam_true)^2))
  mise_l    <- mean((lam_hat - lam_true)^2)
  I_true <- mean(lam_true)
  I_hat  <- mean(lam_hat)
  relI   <- (I_hat - I_true) / I_true

  list(
    ok = TRUE,
    n = n,
    beta_hat = beta_hat,
    gamma_hat = gamma_hat,
    coef_xyz = c(betax = bx_hat, betay = by_hat, betaz = bz_hat),
    rmse_beta = rmse_b,
    ise_beta  = ise_b,
    logrmse_lambda = logrmse_l,
    rmse_lambda = rmse_l,
    mise_lambda = mise_l,
    I_true = I_true,
    I_hat  = I_hat,
    relerr_I = relI,
    fpca_basis_mode = fpca_basis_mode,
    mark_mode = mark_mode,
    N_bg = N_bg_used,

    pve = pve,
    pve_cum = pve_cum,
    pve_cum_last = pve_cum_last
  )
}

scenario_grid_main <- function(ranges = c(0.10, 0.25, 0.50),
                               nuggets = c(0.05, 0.15, 0.30),
                               target_Ns = c(100, 500, 1000)) {
  grid <- expand.grid(
    grf_range = ranges,
    grf_nugget = nuggets,
    target_N = target_Ns,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$scenario_key <- sprintf("R=%.2f__NG=%.2f__N=%d",
                               grid$grf_range, grid$grf_nugget, grid$target_N)
  grid
}

scenario_grid_joint_sens <- function(interp_ps = c(2,4,8,16),
                                     K_fpcas  = c(2,4,6,8,10)) {
  base <- data.frame(
    scenario_class = c("worst","best"),
    grf_range  = c(0.50, 0.25),
    grf_nugget = c(0.05, 0.05),
    target_N   = c(100,  500),
    stringsAsFactors = FALSE
  )

  combos <- expand.grid(
    interp.p = interp_ps,
    K_fpca   = K_fpcas,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  grid <- merge(base, combos, by = NULL)
  grid$scenario_key <- sprintf(
    "%s__R=%.2f__NG=%.2f__N=%d__p=%g__Kfpca=%d",
    grid$scenario_class,
    grid$grf_range, grid$grf_nugget, grid$target_N,
    grid$interp.p, grid$K_fpca
  )
  grid
}

scenario_grid_stress <- function() {
  base <- data.frame(
    scenario_class = c("worst","best"),
    grf_range  = c(0.50, 0.25),
    grf_nugget = c(0.05, 0.05),
    target_N   = c(100,  1000),
    stringsAsFactors = FALSE
  )

  combos <- expand.grid(
    fpca_basis_mode = c("event","background"),
    mark_mode = c("idw","oracle"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  grid <- merge(base, combos, by = NULL)





  grid$intensity_case <- "fun_only"

  grid$scenario_key <- sprintf(
    "%s__R=%.2f__NG=%.2f__N=%d__basis=%s__mark=%s",
    grid$scenario_class,
    grid$grf_range, grid$grf_nugget, grid$target_N,
    grid$fpca_basis_mode, grid$mark_mode
  )
  grid
}

run_plan <- function(grid,
                     nsim = 100,
                     base_seed = 1234,
                     u_grid = seq(0,1,length.out=101),
                     K0 = 6,

                     sigma_eps = 0.05,
                     trunc_eps = 4,
                     trunc_xi  = 4,
                     L_rff = 128L,
                     mc_n = 20000L,

                     mult_fit = 30,

                     bg_mult = 50,
                     N_bg_min = 3000,
                     N_bg_max = 10000,

                     n_eval = 5000,

                     default_interp_p = 4,
                     default_K_fpca   = 5,
                     default_fpca_basis_mode = "event",
                     default_mark_mode = "idw",
                     default_intensity_case = "fun_only",

                     save_dir = NULL,
                     save_name = "sim_results.rds",
                     verbose = TRUE) {

  betas <- make_beta_scenarios(u_grid)
  beta_names <- names(betas)
  if (!length(beta_names)) stop("make_beta_scenarios returned empty list.")

  Phi_true <- make_true_phi_fourier(u_grid, K0 = K0)


  w_u <- trapz_weights(u_grid)
  betas_centered <- lapply(betas, function(b) w_center(b, w_u))


  results <- setNames(vector("list", length(beta_names)), beta_names)
  for (bn in beta_names) results[[bn]] <- setNames(vector("list", nrow(grid)), grid$scenario_key)

  for (i in seq_len(nrow(grid))) {
    sc <- grid[i, , drop = FALSE]
    key <- sc$scenario_key

    rr <- as.numeric(sc$grf_range)
    ng <- as.numeric(sc$grf_nugget)
    Nt <- as.integer(sc$target_N)

    ip <- as.numeric(sc$interp.p %||% default_interp_p)
    Kf <- as.integer(sc$K_fpca   %||% default_K_fpca)

    basis_mode <- as.character(sc$fpca_basis_mode %||% default_fpca_basis_mode)
    mark_mode  <- as.character(sc$mark_mode %||% default_mark_mode)
    case_int   <- as.character(sc$intensity_case %||% default_intensity_case)


    beta_xyz <- switch(case_int,
                       fun_only = c(0,0,0),
                       x_fun    = c(1,0,0),
                       xyz_fun  = c(1,1,1),
                       stop("Unknown intensity_case: ", case_int)
    )

    if (isTRUE(verbose)) {
      cat("\n[SCENARIO ", i, "/", nrow(grid), "] ", key,
          " | R=", rr, " NG=", ng, " N=", Nt,
          " | interp.p=", ip, " K_fpca=", Kf,
          " | basis=", basis_mode, " mark=", mark_mode,
          " | case=", case_int, "\n", sep = "")
    }

    for (bn in beta_names) {

      beta_u <- betas_centered[[bn]]

      beta_hat_mat <- matrix(NA_real_, nrow = length(u_grid), ncol = nsim)
      rmse_beta <- rep(NA_real_, nsim)
      ise_beta  <- rep(NA_real_, nsim)

      logrmse_l <- rep(NA_real_, nsim)
      rmse_l    <- rep(NA_real_, nsim)
      mise_l    <- rep(NA_real_, nsim)
      I_true    <- rep(NA_real_, nsim)
      I_hat     <- rep(NA_real_, nsim)
      relI      <- rep(NA_real_, nsim)

      n_obs     <- integer(nsim)
      ok_rep    <- logical(nsim)

      coef_xyz_mat <- matrix(NA_real_, nrow = nsim, ncol = 3,
                             dimnames = list(NULL, c("betax","betay","betaz")))

      N_bg_vec <- rep(NA_integer_, nsim)

      pve_cum_last <- rep(NA_real_, nsim)

      beta_idx <- match(bn, beta_names)
      if (is.na(beta_idx)) beta_idx <- 0L

      for (s in seq_len(nsim)) {

        seed_s <- base_seed +
          1000000L * as.integer(beta_idx) +
          100000L  * as.integer(i) +
          1000L    * as.integer(s)

        out <- run_one_rep_fpca(
          beta_u = beta_u,
          u_grid = u_grid,
          Phi_true = Phi_true,
          range = rr,
          nugget = ng,
          target_N = Nt,
          beta_xyz = beta_xyz,
          K0 = K0,
          K_fpca = Kf,
          interp.p = ip,
          mult_fit = mult_fit,
          fpca_basis_mode = basis_mode,
          bg_mult = bg_mult, N_bg_min = N_bg_min, N_bg_max = N_bg_max,
          mark_mode = mark_mode,
          sigma_eps = sigma_eps, trunc_eps = trunc_eps,
          trunc_xi = trunc_xi, L_rff = L_rff,
          n_eval = n_eval,
          mc_n = mc_n,
          seed = seed_s
        )

        ok_rep[s] <- isTRUE(out$ok)
        n_obs[s]  <- out$n
        N_bg_vec[s] <- out$N_bg

        beta_hat_mat[, s] <- out$beta_hat
        rmse_beta[s] <- out$rmse_beta
        ise_beta[s]  <- out$ise_beta

        logrmse_l[s] <- out$logrmse_lambda
        rmse_l[s]    <- out$rmse_lambda
        mise_l[s]    <- out$mise_lambda
        I_true[s]    <- out$I_true
        I_hat[s]     <- out$I_hat
        relI[s]      <- out$relerr_I

        coef_xyz_mat[s, ] <- out$coef_xyz

        pve_cum_last[s] <- out$pve_cum_last
      }


      mean_beta <- rowMeans(beta_hat_mat, na.rm = TRUE)
      sd_beta   <- apply(beta_hat_mat, 1, sd, na.rm = TRUE)

      results[[bn]][[key]] <- list(
        beta_name = bn,
        scenario_key = key,
        params = list(
          range = rr, nugget = ng, target_N = Nt,
          interp.p = ip, K_fpca = Kf,
          fpca_basis_mode = basis_mode,
          mark_mode = mark_mode,
          intensity_case = case_int,
          beta_xyz_true = beta_xyz,
          mult_fit = mult_fit,
          bg_mult = bg_mult, N_bg_min = N_bg_min, N_bg_max = N_bg_max
        ),
        u_grid = u_grid,
        beta_true = beta_u,
        beta_hat = beta_hat_mat,
        mean_beta = mean_beta,
        sd_beta = sd_beta,
        ok = ok_rep,
        n_obs = n_obs,
        N_bg = N_bg_vec,
        rmse_beta = rmse_beta,
        ise_beta  = ise_beta,
        logrmse_lambda = logrmse_l,
        rmse_lambda = rmse_l,
        mise_lambda = mise_l,
        I_true = I_true,
        I_hat  = I_hat,
        relerr_I = relI,
        coef_xyz_hat = coef_xyz_mat,
        pve_cum_last = pve_cum_last
      )
    }
  }

  out <- list(
    grid = grid,
    nsim = nsim,
    base_seed = base_seed,
    u_grid = u_grid,
    K0 = K0,
    sigma_eps = sigma_eps,
    trunc_eps = trunc_eps,
    trunc_xi  = trunc_xi,
    L_rff = L_rff,
    mult_fit = mult_fit,
    bg_mult = bg_mult,
    N_bg_min = N_bg_min,
    N_bg_max = N_bg_max,
    n_eval = n_eval,
    mc_n = mc_n,
    results = results
  )

  if (!is.null(save_dir)) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(out, file = file.path(save_dir, save_name))
  }

  out
}

summarise_table <- function(out_obj) {
  res <- out_obj$results
  beta_names <- names(res)

  rows <- list()
  k <- 0L

  for (bn in beta_names) {
    scn_list <- res[[bn]]
    for (key in names(scn_list)) {
      sc <- scn_list[[key]]
      if (is.null(sc)) next

      ok_idx <- which(as.logical(sc$ok) &
                        is.finite(sc$rmse_beta) &
                        is.finite(sc$mise_lambda))

      n_ok <- length(ok_idx)

      p <- sc$params
      beta_xyz_true <- p$beta_xyz_true

      rmse_b_m <- if (n_ok) mean(sc$rmse_beta[ok_idx]) else NA_real_
      rmse_b_s <- if (n_ok > 1) sd(sc$rmse_beta[ok_idx]) else NA_real_
      ise_b_m  <- if (n_ok) mean(sc$ise_beta[ok_idx]) else NA_real_

      mise_m <- if (n_ok) mean(sc$mise_lambda[ok_idx]) else NA_real_
      logrmse_m <- if (n_ok) mean(sc$logrmse_lambda[ok_idx]) else NA_real_
      rmse_l_m  <- if (n_ok) mean(sc$rmse_lambda[ok_idx]) else NA_real_

      relI_m <- if (n_ok) mean(sc$relerr_I[ok_idx]) else NA_real_

      pve_last_m <- if (n_ok) mean(sc$pve_cum_last[ok_idx], na.rm = TRUE) else NA_real_

      bx_hat_m <- if (n_ok) mean(sc$coef_xyz_hat[ok_idx, "betax"], na.rm = TRUE) else NA_real_
      by_hat_m <- if (n_ok) mean(sc$coef_xyz_hat[ok_idx, "betay"], na.rm = TRUE) else NA_real_
      bz_hat_m <- if (n_ok) mean(sc$coef_xyz_hat[ok_idx, "betaz"], na.rm = TRUE) else NA_real_

      k <- k + 1L
      rows[[k]] <- data.frame(
        beta_name = bn,
        scenario_key = key,
        grf_range = as.numeric(p$range),
        grf_nugget = as.numeric(p$nugget),
        target_N = as.integer(p$target_N),
        interp.p = as.numeric(p$interp.p),
        K_fpca   = as.integer(p$K_fpca),
        fpca_basis_mode = as.character(p$fpca_basis_mode),
        mark_mode = as.character(p$mark_mode),
        intensity_case = as.character(p$intensity_case),
        betax_true = beta_xyz_true[1],
        betay_true = beta_xyz_true[2],
        betaz_true = beta_xyz_true[3],
        n_ok = n_ok,
        n_obs_mean = mean(sc$n_obs, na.rm = TRUE),
        rmse_beta_mean = rmse_b_m,
        rmse_beta_sd   = rmse_b_s,
        ise_beta_mean  = ise_b_m,
        mise_lambda_mean = mise_m,
        logrmse_lambda_mean = logrmse_m,
        rmse_lambda_mean    = rmse_l_m,
        relerr_I_mean       = relI_m,
        pve_cum_last_mean   = pve_last_m,
        betax_hat_mean = bx_hat_m,
        betay_hat_mean = by_hat_m,
        betaz_hat_mean = bz_hat_m,
        stringsAsFactors = FALSE
      )
    }
  }

  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  tab
}

plot_beta_bands <- function(scn_res, main = NULL, lwd = 2) {
  u <- scn_res$u_grid
  btrue <- scn_res$beta_true
  m <- scn_res$mean_beta
  s <- scn_res$sd_beta

  if (is.null(main)) main <- paste(scn_res$beta_name, "-", scn_res$scenario_key)

  plot(u, btrue, type="l", xlab="u", ylab="beta(u)", main=main, lwd=lwd)
  lines(u, m, lty=2, lwd=lwd)
  lines(u, m + s, lty=3, lwd=lwd)
  lines(u, m - s, lty=3, lwd=lwd)
}

scenario_score_S <- function(out_main) {
  res <- out_main$results
  beta_names <- names(res)


  keys <- names(res[[beta_names[1]]])

  S <- numeric(length(keys))
  MISE <- numeric(length(keys))
  NOK <- integer(length(keys))

  names(S) <- keys
  names(MISE) <- keys
  names(NOK) <- keys

  for (k in seq_along(keys)) {
    key <- keys[k]
    vals_rmse <- c()
    vals_mise <- c()
    n_ok_tot <- 0L

    for (bn in beta_names) {
      sc <- res[[bn]][[key]]
      if (is.null(sc)) next
      ok_idx <- which(as.logical(sc$ok) & is.finite(sc$rmse_beta) & is.finite(sc$mise_lambda))
      n_ok_tot <- n_ok_tot + length(ok_idx)
      if (length(ok_idx)) {
        vals_rmse <- c(vals_rmse, sc$rmse_beta[ok_idx])
        vals_mise <- c(vals_mise, sc$mise_lambda[ok_idx])
      }
    }

    S[k]    <- if (length(vals_rmse)) median(vals_rmse) else NA_real_
    MISE[k] <- if (length(vals_mise)) median(vals_mise) else NA_real_
    NOK[k]  <- n_ok_tot
  }

  data.frame(
    scenario_key = keys,
    S_med_rmse_beta = as.numeric(S),
    MISE_med = as.numeric(MISE),
    n_ok_total = as.integer(NOK),
    stringsAsFactors = FALSE
  )
}

pick_best_mid_worst <- function(out_main) {
  tabS <- scenario_score_S(out_main)
  tabS <- tabS[is.finite(tabS$S_med_rmse_beta), , drop = FALSE]
  tabS <- tabS[order(tabS$S_med_rmse_beta), , drop = FALSE]
  if (nrow(tabS) < 3) stop("Not enough scenarios with finite S to pick best/mid/worst.")

  best <- tabS[1, , drop = FALSE]
  worst <- tabS[nrow(tabS), , drop = FALSE]
  mid <- tabS[ceiling(nrow(tabS)/2), , drop = FALSE]

  list(best = best, intermediate = mid, worst = worst, ranking = tabS)
}

recommend_tuning <- function(out_sens) {
  tab <- summarise_table(out_sens)




  ok <- is.finite(tab$rmse_beta_mean)
  tab <- tab[ok, , drop = FALSE]

  agg <- aggregate(rmse_beta_mean ~ interp.p + K_fpca, data = tab, FUN = mean)
  agg <- agg[order(agg$rmse_beta_mean), , drop = FALSE]
  if (nrow(agg) == 0) stop("No finite RMSE in sensitivity results.")

  list(interp.p = agg$interp.p[1], K_fpca = agg$K_fpca[1], table = agg)
}

run_one_rep_fpca_compare <- function(beta_u, u_grid, Phi_true,
                                     range, nugget, target_N,
                                     intensity_case = c("x_fun","xyz_fun"),
                                     K_fpca = 6,
                                     interp.p = 16,
                                     mult_fit = 30,

                                     fpca_basis_mode = "event",

                                     mark_mode = "idw",

                                     K0 = 6,
                                     sigma_eps = 0.05, trunc_eps = 4,
                                     trunc_xi = 4, L_rff = 128L,
                                     n_eval = 5000,
                                     mc_n = 20000L,
                                     seed = NULL) {

  intensity_case <- match.arg(intensity_case, c("x_fun","xyz_fun"))

  if (!is.null(seed)) set.seed(seed)


  beta_xyz <- switch(intensity_case,
                     x_fun   = c(1,0,0),
                     xyz_fun = c(1,1,1))


  w_u <- trapz_weights(u_grid)
  beta_u <- w_center(beta_u, w_u)
  beta_w <- beta_u * w_u


  a0 <- calibrate_a0_functional(
    beta_w = beta_w,
    Phi_true = Phi_true,
    target_N = target_N,
    beta_xyz = beta_xyz,
    K0 = K0,
    sigma_eps = sigma_eps, trunc_eps = trunc_eps,
    trunc_xi = trunc_xi,
    mc_n = mc_n,
    seed = if (!is.null(seed)) seed + 11L else NULL
  )

  lmax <- lmax_bound_functional(
    a0 = a0, beta_w = beta_w, Phi_true = Phi_true,
    beta_xyz = beta_xyz,
    trunc_xi = trunc_xi,
    sigma_eps = sigma_eps, trunc_eps = trunc_eps
  )


  field <- make_rff_field(
    K0 = K0, range = range, nugget = nugget,
    L = L_rff,
    seed = if (!is.null(seed)) seed + 101L else NULL,
    trunc_xi = trunc_xi
  )


  sim <- simulate_one_pattern_box3(
    beta_w = beta_w, a0 = a0, lmax = lmax,
    field = field, Phi_true = Phi_true,
    beta_xyz = beta_xyz,
    sigma_eps = sigma_eps, trunc_eps = trunc_eps,
    seed = if (!is.null(seed)) seed + 1001L else NULL
  )

  coords  <- sim$coords
  Xevent  <- sim$Xcurves
  n <- nrow(coords)

  if (n < max(15, K_fpca + 3)) {
    return(list(ok = FALSE, n = n, intensity_case = intensity_case))
  }


  rownames(Xevent) <- as.character(seq_len(n))
  fp <- fpca_grid(Xevent, u_grid, K = K_fpca)

  pve_cum_last <- if (length(fp$pve_cum)) fp$pve_cum[length(fp$pve_cum)] else NA_real_

  scores_event <- fpca_project_scores(Xevent, fp)
  colnames(scores_event) <- paste0("age_pc", seq_len(K_fpca))

  X_s3d <- build_s3dpm_fpca(coords, scores_event)
  mark_cols <- paste0("age_pc", seq_len(K_fpca))


  oracle_fun <- NULL
  if (mark_mode == "oracle") {
    oracle_fun <- make_oracle_fun_scores(field, Phi_true, u_grid, fp, sigma_eps, trunc_eps)
  }


  rhs_correct <- paste(mark_cols, collapse = " + ")
  if (beta_xyz[1] != 0) rhs_correct <- paste("x +", rhs_correct)
  if (beta_xyz[2] != 0) rhs_correct <- paste("y +", rhs_correct)
  if (beta_xyz[3] != 0) rhs_correct <- paste("z +", rhs_correct)

  fit_correct <- try(
    fit_funpois_fpca(
      X_s3d = X_s3d,
      rhs_terms = rhs_correct,
      interp.p = interp.p,
      mult_fit = mult_fit,
      seed_fit = if (!is.null(seed)) seed + 2001L else NULL,
      mark_mode = mark_mode,
      mark_oracle_fun = oracle_fun
    ),
    silent = TRUE
  )


  rhs_miss <- if (intensity_case == "x_fun") "x" else "x + y + z"

  fit_miss <- try(
    fit_funpois_fpca(
      X_s3d = X_s3d,
      rhs_terms = rhs_miss,
      interp.p = interp.p,
      mult_fit = mult_fit,
      seed_fit = if (!is.null(seed)) seed + 3001L else NULL,
      mark_mode = mark_mode,
      mark_oracle_fun = oracle_fun
    ),
    silent = TRUE
  )

  ok_correct <- !inherits(fit_correct, "try-error")
  ok_miss    <- !inherits(fit_miss, "try-error")


  coords_eval <- cbind(runif(n_eval), runif(n_eval), runif(n_eval))
  Xeval <- simulate_X_curves(coords_eval, field, Phi_true, sigma_eps = sigma_eps, trunc_eps = trunc_eps)
  integ_eval <- as.vector(Xeval %*% beta_w)
  lin_eval <- beta_xyz[1]*coords_eval[,1] + beta_xyz[2]*coords_eval[,2] + beta_xyz[3]*coords_eval[,3]
  lam_true <- exp(a0 + integ_eval + lin_eval)
  I_true <- mean(lam_true)


  scores_eval <- fpca_project_scores(Xeval, fp)
  colnames(scores_eval) <- mark_cols
  df_eval <- data.frame(
    x = coords_eval[,1], y = coords_eval[,2], z = coords_eval[,3],
    as.data.frame(scores_eval),
    check.names = FALSE
  )

  out <- list(
    ok = TRUE,
    n = n,
    intensity_case = intensity_case,
    pve_cum_last = pve_cum_last,
    I_true = I_true
  )


  if (ok_correct) {
    cc <- coef(fit_correct$mod_global)
    gamma_hat <- as.numeric(cc[mark_cols])
    beta_hat <- beta_hat_from_gamma_fpca(gamma_hat, fp, w_u)

    rmse_b <- sqrt(mean((beta_hat - beta_u)^2))
    ise_b  <- sum(((beta_hat - beta_u)^2) * w_u) / sum(w_u)

    lam_hat <- as.numeric(predict(fit_correct$mod_global, newdata = df_eval, type = "response"))

    out$correct <- list(
      rmse_beta = rmse_b,
      ise_beta  = ise_b,
      mise_lambda = mean((lam_hat - lam_true)^2),
      logrmse_lambda = sqrt(mean((log(lam_hat) - log(lam_true))^2)),
      relerr_I = (mean(lam_hat) - I_true) / I_true,
      coef_xyz = c(
        betax = if ("x" %in% names(cc)) as.numeric(cc["x"]) else 0,
        betay = if ("y" %in% names(cc)) as.numeric(cc["y"]) else 0,
        betaz = if ("z" %in% names(cc)) as.numeric(cc["z"]) else 0
      )
    )
  } else {
    out$correct <- list(rmse_beta = NA_real_, ise_beta = NA_real_,
                        mise_lambda = NA_real_, logrmse_lambda = NA_real_,
                        relerr_I = NA_real_,
                        coef_xyz = c(betax=NA_real_, betay=NA_real_, betaz=NA_real_))
  }


  if (ok_miss) {
    cc2 <- coef(fit_miss$mod_global)
    lam_hat2 <- as.numeric(predict(fit_miss$mod_global, newdata = df_eval, type = "response"))

    out$misspec <- list(
      mise_lambda = mean((lam_hat2 - lam_true)^2),
      logrmse_lambda = sqrt(mean((log(lam_hat2) - log(lam_true))^2)),
      relerr_I = (mean(lam_hat2) - I_true) / I_true,
      coef_xyz = c(
        betax = if ("x" %in% names(cc2)) as.numeric(cc2["x"]) else 0,
        betay = if ("y" %in% names(cc2)) as.numeric(cc2["y"]) else 0,
        betaz = if ("z" %in% names(cc2)) as.numeric(cc2["z"]) else 0
      )
    )
  } else {
    out$misspec <- list(mise_lambda = NA_real_, logrmse_lambda = NA_real_,
                        relerr_I = NA_real_,
                        coef_xyz = c(betax=NA_real_, betay=NA_real_, betaz=NA_real_))
  }

  out
}

run_plan_coords_misspec <- function(selected_triplet,
                                    nsim = 100,
                                    base_seed = 1234,
                                    u_grid = seq(0,1,length.out=101),
                                    K0 = 6,
                                    sigma_eps = 0.05,
                                    trunc_eps = 4,
                                    trunc_xi  = 4,
                                    L_rff = 128L,
                                    mc_n = 20000L,
                                    mult_fit = 30,
                                    interp.p = 16,
                                    K_fpca = 6,
                                    n_eval = 5000,
                                    mark_mode = "idw",
                                    verbose = TRUE) {

  betas <- make_beta_scenarios(u_grid)
  beta_names <- names(betas)
  w_u <- trapz_weights(u_grid)
  betas_centered <- lapply(betas, function(b) w_center(b, w_u))
  Phi_true <- make_true_phi_fourier(u_grid, K0 = K0)


  sc_tab <- rbind(
    transform(selected_triplet$best, scenario_class="best"),
    transform(selected_triplet$intermediate, scenario_class="intermediate"),
    transform(selected_triplet$worst, scenario_class="worst")
  )


  parse_key <- function(key) {

    R  <- as.numeric(sub(".*R=([0-9.]+).*", "\\1", key))
    NG <- as.numeric(sub(".*NG=([0-9.]+).*", "\\1", key))
    N  <- as.integer(sub(".*__N=([0-9]+)$", "\\1", key))
    c(range=R, nugget=NG, target_N=N)
  }

  grid <- do.call(rbind, lapply(seq_len(nrow(sc_tab)), function(i) {
    key <- sc_tab$scenario_key[i]
    pars <- parse_key(key)
    data.frame(
      scenario_class = sc_tab$scenario_class[i],
      grf_range = pars["range"],
      grf_nugget = pars["nugget"],
      target_N = pars["target_N"],
      intensity_case = c("x_fun","xyz_fun"),
      stringsAsFactors = FALSE
    )
  }))

  grid$interp.p <- interp.p
  grid$K_fpca <- K_fpca
  grid$scenario_key <- sprintf("%s__R=%.2f__NG=%.2f__N=%d__case=%s",
                               grid$scenario_class, grid$grf_range, grid$grf_nugget,
                               grid$target_N, grid$intensity_case)


  results <- setNames(vector("list", length(beta_names)), beta_names)
  for (bn in beta_names) results[[bn]] <- setNames(vector("list", nrow(grid)), grid$scenario_key)

  for (i in seq_len(nrow(grid))) {
    sc <- grid[i, , drop = FALSE]
    key <- sc$scenario_key

    rr <- as.numeric(sc$grf_range)
    ng <- as.numeric(sc$grf_nugget)
    Nt <- as.integer(sc$target_N)
    case_int <- as.character(sc$intensity_case)

    if (isTRUE(verbose)) {
      cat("\n[COORDS+MISSPEC ", i, "/", nrow(grid), "] ", key,
          " | R=", rr, " NG=", ng, " N=", Nt,
          " | interp.p=", interp.p, " K_fpca=", K_fpca,
          " | case=", case_int, "\n", sep = "")
    }

    for (bn in beta_names) {
      beta_u <- betas_centered[[bn]]

      ok_rep <- logical(nsim)
      n_obs  <- integer(nsim)
      pve_last <- rep(NA_real_, nsim)


      rmse_beta_c <- rep(NA_real_, nsim)
      ise_beta_c  <- rep(NA_real_, nsim)
      mise_c      <- rep(NA_real_, nsim)
      logrmse_c   <- rep(NA_real_, nsim)
      relI_c      <- rep(NA_real_, nsim)
      coef_c      <- matrix(NA_real_, nsim, 3, dimnames=list(NULL, c("betax","betay","betaz")))


      mise_m      <- rep(NA_real_, nsim)
      logrmse_m   <- rep(NA_real_, nsim)
      relI_m      <- rep(NA_real_, nsim)
      coef_m      <- matrix(NA_real_, nsim, 3, dimnames=list(NULL, c("betax","betay","betaz")))

      beta_idx <- match(bn, beta_names)
      for (s in seq_len(nsim)) {
        seed_s <- base_seed +
          1000000L * as.integer(beta_idx) +
          100000L  * as.integer(i) +
          1000L    * as.integer(s)

        out <- run_one_rep_fpca_compare(
          beta_u = beta_u, u_grid = u_grid, Phi_true = Phi_true,
          range = rr, nugget = ng, target_N = Nt,
          intensity_case = case_int,
          K_fpca = K_fpca,
          interp.p = interp.p,
          mult_fit = mult_fit,
          fpca_basis_mode = "event",
          mark_mode = mark_mode,
          K0 = K0,
          sigma_eps = sigma_eps, trunc_eps = trunc_eps,
          trunc_xi = trunc_xi, L_rff = L_rff,
          n_eval = n_eval,
          mc_n = mc_n,
          seed = seed_s
        )

        ok_rep[s] <- isTRUE(out$ok) && is.finite(out$n)
        n_obs[s]  <- out$n %||% NA_integer_
        pve_last[s] <- out$pve_cum_last %||% NA_real_

        if (isTRUE(out$ok)) {

          rmse_beta_c[s] <- out$correct$rmse_beta
          ise_beta_c[s]  <- out$correct$ise_beta
          mise_c[s]      <- out$correct$mise_lambda
          logrmse_c[s]   <- out$correct$logrmse_lambda
          relI_c[s]      <- out$correct$relerr_I
          coef_c[s, ]    <- out$correct$coef_xyz


          mise_m[s]      <- out$misspec$mise_lambda
          logrmse_m[s]   <- out$misspec$logrmse_lambda
          relI_m[s]      <- out$misspec$relerr_I
          coef_m[s, ]    <- out$misspec$coef_xyz
        }
      }

      results[[bn]][[key]] <- list(
        beta_name = bn,
        scenario_key = key,
        params = list(range=rr, nugget=ng, target_N=Nt,
                      interp.p=interp.p, K_fpca=K_fpca,
                      intensity_case=case_int,
                      mark_mode=mark_mode),
        ok = ok_rep,
        n_obs = n_obs,
        pve_cum_last = pve_last,
        correct = list(
          rmse_beta = rmse_beta_c,
          ise_beta  = ise_beta_c,
          mise_lambda = mise_c,
          logrmse_lambda = logrmse_c,
          relerr_I = relI_c,
          coef_xyz_hat = coef_c
        ),
        misspec = list(
          mise_lambda = mise_m,
          logrmse_lambda = logrmse_m,
          relerr_I = relI_m,
          coef_xyz_hat = coef_m
        )
      )
    }
  }

  list(grid = grid, nsim = nsim, base_seed = base_seed, results = results,
       interp.p = interp.p, K_fpca = K_fpca)
}

run_plan_parallel <- function(grid,
                              nsim = 100,
                              base_seed = 1234,
                              u_grid = seq(0,1,length.out=101),
                              K0 = 6,
                              sigma_eps = 0.05,
                              trunc_eps = 4,
                              trunc_xi  = 4,
                              L_rff = 128L,
                              mc_n = 20000L,
                              mult_fit = 30,
                              bg_mult = 50,
                              N_bg_min = 3000,
                              N_bg_max = 10000,
                              n_eval = 5000,
                              default_interp_p = 4,
                              default_K_fpca   = 5,
                              default_fpca_basis_mode = "event",
                              default_mark_mode = "idw",
                              default_intensity_case = "fun_only",
                              save_dir = NULL,
                              save_name = "sim_results.rds",
                              verbose = TRUE,
                              future_workers = NULL,
                              future_packages = c("mgcv","splines","spatstat.geom")) {

  if (!is.null(future_workers)) {
    future::plan(future::multisession, workers = future_workers)
  }

  betas <- make_beta_scenarios(u_grid)
  beta_names <- names(betas)
  if (!length(beta_names)) stop("make_beta_scenarios returned empty list.")
  Phi_true <- make_true_phi_fourier(u_grid, K0 = K0)

  w_u <- trapz_weights(u_grid)
  betas_centered <- lapply(betas, function(b) w_center(b, w_u))

  tasks <- expand.grid(
    i = seq_len(nrow(grid)),
    bn = beta_names,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  worker <- function(row_idx) {
    i  <- tasks$i[row_idx]
    bn <- tasks$bn[row_idx]

    sc  <- grid[i, , drop = FALSE]
    key <- sc$scenario_key

    rr <- as.numeric(sc$grf_range)
    ng <- as.numeric(sc$grf_nugget)
    Nt <- as.integer(sc$target_N)

    ip <- as.numeric(sc$interp.p %||% default_interp_p)
    Kf <- as.integer(sc$K_fpca   %||% default_K_fpca)

    basis_mode <- as.character(sc$fpca_basis_mode %||% default_fpca_basis_mode)
    mark_mode  <- as.character(sc$mark_mode %||% default_mark_mode)
    case_int   <- as.character(sc$intensity_case %||% default_intensity_case)

    beta_xyz <- switch(case_int,
                       fun_only = c(0,0,0),
                       x_fun    = c(1,0,0),
                       xyz_fun  = c(1,1,1),
                       stop("Unknown intensity_case: ", case_int))

    beta_u <- betas_centered[[bn]]
    beta_hat_mat <- matrix(NA_real_, nrow = length(u_grid), ncol = nsim)

    rmse_beta <- rep(NA_real_, nsim)
    ise_beta  <- rep(NA_real_, nsim)

    logrmse_l <- rep(NA_real_, nsim)
    rmse_l    <- rep(NA_real_, nsim)
    mise_l    <- rep(NA_real_, nsim)

    I_true    <- rep(NA_real_, nsim)
    I_hat     <- rep(NA_real_, nsim)
    relI      <- rep(NA_real_, nsim)

    n_obs     <- integer(nsim)
    ok_rep    <- logical(nsim)

    coef_xyz_mat <- matrix(NA_real_, nrow = nsim, ncol = 3,
                           dimnames = list(NULL, c("betax","betay","betaz")))

    N_bg_vec <- rep(NA_integer_, nsim)
    pve_cum_last <- rep(NA_real_, nsim)

    beta_idx <- match(bn, beta_names)
    if (is.na(beta_idx)) beta_idx <- 0L

    for (s in seq_len(nsim)) {
      seed_s <- base_seed +
        1000000L * as.integer(beta_idx) +
        100000L  * as.integer(i) +
        1000L    * as.integer(s)

      out <- run_one_rep_fpca(
        beta_u = beta_u,
        u_grid = u_grid,
        Phi_true = Phi_true,
        range = rr,
        nugget = ng,
        target_N = Nt,
        beta_xyz = beta_xyz,
        K0 = K0,
        K_fpca = Kf,
        interp.p = ip,
        mult_fit = mult_fit,
        fpca_basis_mode = basis_mode,
        bg_mult = bg_mult, N_bg_min = N_bg_min, N_bg_max = N_bg_max,
        mark_mode = mark_mode,
        sigma_eps = sigma_eps, trunc_eps = trunc_eps,
        trunc_xi = trunc_xi, L_rff = L_rff,
        n_eval = n_eval,
        mc_n = mc_n,
        seed = seed_s
      )

      ok_rep[s] <- isTRUE(out$ok)
      n_obs[s]  <- out$n
      N_bg_vec[s] <- out$N_bg

      beta_hat_mat[, s] <- out$beta_hat
      rmse_beta[s] <- out$rmse_beta
      ise_beta[s]  <- out$ise_beta

      logrmse_l[s] <- out$logrmse_lambda
      rmse_l[s]    <- out$rmse_lambda
      mise_l[s]    <- out$mise_lambda

      I_true[s]    <- out$I_true
      I_hat[s]     <- out$I_hat
      relI[s]      <- out$relerr_I

      coef_xyz_mat[s, ] <- out$coef_xyz
      pve_cum_last[s] <- out$pve_cum_last
    }

    mean_beta <- rowMeans(beta_hat_mat, na.rm = TRUE)
    sd_beta   <- apply(beta_hat_mat, 1, sd, na.rm = TRUE)

    scenario_obj <- list(
      beta_name = bn,
      scenario_key = key,
      params = list(
        range = rr, nugget = ng, target_N = Nt,
        interp.p = ip, K_fpca = Kf,
        fpca_basis_mode = basis_mode,
        mark_mode = mark_mode,
        intensity_case = case_int,
        beta_xyz_true = beta_xyz,
        mult_fit = mult_fit,
        bg_mult = bg_mult, N_bg_min = N_bg_min, N_bg_max = N_bg_max
      ),
      u_grid = u_grid,
      beta_true = beta_u,
      beta_hat = beta_hat_mat,
      mean_beta = mean_beta,
      sd_beta = sd_beta,
      ok = ok_rep,
      n_obs = n_obs,
      N_bg = N_bg_vec,
      rmse_beta = rmse_beta,
      ise_beta  = ise_beta,
      logrmse_lambda = logrmse_l,
      rmse_lambda = rmse_l,
      mise_lambda = mise_l,
      I_true = I_true,
      I_hat  = I_hat,
      relerr_I = relI,
      coef_xyz_hat = coef_xyz_mat,
      pve_cum_last = pve_cum_last
    )

    list(bn = bn, key = key, scenario_obj = scenario_obj)
  }

  if (isTRUE(verbose)) {
    cat("\n[run_plan_parallel] Parallel tasks (scenario x beta): ", nrow(tasks), "\n", sep = "")
  }

  outs <- future.apply::future_lapply(
    X = seq_len(nrow(tasks)),
    FUN = worker,
    future.seed = TRUE,
    future.packages = future_packages
  )


  results <- setNames(vector("list", length(beta_names)), beta_names)
  for (bn in beta_names) results[[bn]] <- setNames(vector("list", nrow(grid)), grid$scenario_key)

  for (o in outs) {
    results[[o$bn]][[o$key]] <- o$scenario_obj
  }

  out <- list(
    grid = grid,
    nsim = nsim,
    base_seed = base_seed,
    u_grid = u_grid,
    K0 = K0,
    sigma_eps = sigma_eps,
    trunc_eps = trunc_eps,
    trunc_xi  = trunc_xi,
    L_rff = L_rff,
    mult_fit = mult_fit,
    bg_mult = bg_mult,
    N_bg_min = N_bg_min,
    N_bg_max = N_bg_max,
    n_eval = n_eval,
    mc_n = mc_n,
    results = results
  )

  if (!is.null(save_dir)) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(out, file = file.path(save_dir, save_name))
  }

  out
}

run_plan_coords_misspec_parallel <- function(selected_triplet,
                                             nsim = 100,
                                             base_seed = 1234,
                                             u_grid = seq(0,1,length.out=101),
                                             K0 = 6,
                                             sigma_eps = 0.05,
                                             trunc_eps = 4,
                                             trunc_xi  = 4,
                                             L_rff = 128L,
                                             mc_n = 20000L,
                                             mult_fit = 30,
                                             interp.p = 16,
                                             K_fpca = 6,
                                             n_eval = 5000,
                                             mark_mode = "idw",
                                             verbose = TRUE,
                                             future_packages = c("mgcv","splines","spatstat.geom")) {

  betas <- make_beta_scenarios(u_grid)
  beta_names <- names(betas)
  w_u <- trapz_weights(u_grid)
  betas_centered <- lapply(betas, function(b) w_center(b, w_u))
  Phi_true <- make_true_phi_fourier(u_grid, K0 = K0)


  sc_tab <- rbind(
    transform(selected_triplet$best, scenario_class="best"),
    transform(selected_triplet$intermediate, scenario_class="intermediate"),
    transform(selected_triplet$worst, scenario_class="worst")
  )

  parse_key <- function(key) {
    R  <- as.numeric(sub(".*R=([0-9.]+).*", "\\1", key))
    NG <- as.numeric(sub(".*NG=([0-9.]+).*", "\\1", key))
    N  <- as.integer(sub(".*__N=([0-9]+)$", "\\1", key))
    c(range=R, nugget=NG, target_N=N)
  }

  grid <- do.call(rbind, lapply(seq_len(nrow(sc_tab)), function(i) {
    key <- sc_tab$scenario_key[i]
    pars <- parse_key(key)
    data.frame(
      scenario_class = sc_tab$scenario_class[i],
      grf_range = pars["range"],
      grf_nugget = pars["nugget"],
      target_N = pars["target_N"],
      intensity_case = c("x_fun","xyz_fun"),
      stringsAsFactors = FALSE
    )
  }))

  grid$interp.p <- interp.p
  grid$K_fpca <- K_fpca
  grid$scenario_key <- sprintf("%s__R=%.2f__NG=%.2f__N=%d__case=%s",
                               grid$scenario_class, grid$grf_range, grid$grf_nugget,
                               grid$target_N, grid$intensity_case)


  tasks <- expand.grid(
    i = seq_len(nrow(grid)),
    bn = beta_names,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  worker <- function(row_idx) {
    i  <- tasks$i[row_idx]
    bn <- tasks$bn[row_idx]

    sc <- grid[i, , drop = FALSE]
    key <- sc$scenario_key

    rr <- as.numeric(sc$grf_range)
    ng <- as.numeric(sc$grf_nugget)
    Nt <- as.integer(sc$target_N)
    case_int <- as.character(sc$intensity_case)

    beta_u <- betas_centered[[bn]]


    ok_rep <- logical(nsim)
    n_obs  <- integer(nsim)
    pve_last <- rep(NA_real_, nsim)


    rmse_beta_c <- rep(NA_real_, nsim)
    ise_beta_c  <- rep(NA_real_, nsim)
    mise_c      <- rep(NA_real_, nsim)
    logrmse_c   <- rep(NA_real_, nsim)
    relI_c      <- rep(NA_real_, nsim)
    coef_c      <- matrix(NA_real_, nsim, 3, dimnames=list(NULL, c("betax","betay","betaz")))


    mise_m      <- rep(NA_real_, nsim)
    logrmse_m   <- rep(NA_real_, nsim)
    relI_m      <- rep(NA_real_, nsim)
    coef_m      <- matrix(NA_real_, nsim, 3, dimnames=list(NULL, c("betax","betay","betaz")))

    beta_idx <- match(bn, beta_names)
    for (s in seq_len(nsim)) {
      seed_s <- base_seed +
        1000000L * as.integer(beta_idx) +
        100000L  * as.integer(i) +
        1000L    * as.integer(s)

      out <- run_one_rep_fpca_compare(
        beta_u = beta_u, u_grid = u_grid, Phi_true = Phi_true,
        range = rr, nugget = ng, target_N = Nt,
        intensity_case = case_int,
        K_fpca = K_fpca,
        interp.p = interp.p,
        mult_fit = mult_fit,
        fpca_basis_mode = "event",
        mark_mode = mark_mode,
        K0 = K0,
        sigma_eps = sigma_eps, trunc_eps = trunc_eps,
        trunc_xi = trunc_xi, L_rff = L_rff,
        n_eval = n_eval,
        mc_n = mc_n,
        seed = seed_s
      )

      ok_rep[s] <- isTRUE(out$ok)
      n_obs[s]  <- out$n %||% NA_integer_
      pve_last[s] <- out$pve_cum_last %||% NA_real_

      if (isTRUE(out$ok)) {

        rmse_beta_c[s] <- out$correct$rmse_beta
        ise_beta_c[s]  <- out$correct$ise_beta
        mise_c[s]      <- out$correct$mise_lambda
        logrmse_c[s]   <- out$correct$logrmse_lambda
        relI_c[s]      <- out$correct$relerr_I
        coef_c[s, ]    <- out$correct$coef_xyz


        mise_m[s]      <- out$misspec$mise_lambda
        logrmse_m[s]   <- out$misspec$logrmse_lambda
        relI_m[s]      <- out$misspec$relerr_I
        coef_m[s, ]    <- out$misspec$coef_xyz
      }
    }

    scenario_obj <- list(
      beta_name = bn,
      scenario_key = key,
      params = list(range=rr, nugget=ng, target_N=Nt,
                    interp.p=interp.p, K_fpca=K_fpca,
                    intensity_case=case_int,
                    mark_mode=mark_mode),
      ok = ok_rep,
      n_obs = n_obs,
      pve_cum_last = pve_last,
      correct = list(
        rmse_beta = rmse_beta_c,
        ise_beta  = ise_beta_c,
        mise_lambda = mise_c,
        logrmse_lambda = logrmse_c,
        relerr_I = relI_c,
        coef_xyz_hat = coef_c
      ),
      misspec = list(
        mise_lambda = mise_m,
        logrmse_lambda = logrmse_m,
        relerr_I = relI_m,
        coef_xyz_hat = coef_m
      )
    )

    list(bn = bn, key = key, scenario_obj = scenario_obj)
  }

  if (isTRUE(verbose)) {
    cat("\n[run_plan_coords_misspec_parallel] Parallel tasks (scenario x beta): ",
        nrow(tasks), "\n", sep = "")
  }

  outs <- future.apply::future_lapply(
    X = seq_len(nrow(tasks)),
    FUN = worker,
    future.seed = TRUE,
    future.packages = future_packages
  )

  results <- setNames(vector("list", length(beta_names)), beta_names)
  for (bn in beta_names) results[[bn]] <- setNames(vector("list", nrow(grid)), grid$scenario_key)

  for (o in outs) {
    results[[o$bn]][[o$key]] <- o$scenario_obj
  }

  list(grid = grid, nsim = nsim, base_seed = base_seed,
       interp.p = interp.p, K_fpca = K_fpca,
       results = results)
}

.parse_key_main <- function(key) {
  R  <- as.numeric(sub(".*R=([0-9.]+).*", "\\1", key))
  NG <- as.numeric(sub(".*NG=([0-9.]+).*", "\\1", key))
  N  <- as.integer(sub(".*__N=([0-9]+)$", "\\1", key))
  list(R=R, NG=NG, N=N)
}

curve_mean_weighted <- function(F, u_grid) {
  w_u <- trapz_weights(u_grid)
  as.vector(F %*% (w_u / sum(w_u)))
}

build_s3dpm_mean <- function(coords, xmean, var_name = "xmean_fun") {
  coords <- as.matrix(coords)

  df <- data.frame(
    x = coords[,1],
    y = coords[,2],
    z = coords[,3],
    n_stars = rep(1, nrow(coords)),
    xmean_fun = as.numeric(xmean),
    check.names = FALSE
  )

  names(df)[5] <- var_name

  s3dpm(df = df, names = c("n_stars", var_name))
}

run_one_rep_mean_scalar <- function(beta_u, u_grid, Phi_true,
                                    range, nugget, target_N,
                                    beta_xyz = c(0,0,0),
                                    interp.p = 16,
                                    mult_fit = 30,
                                    mark_mode = c("idw","oracle"),
                                    mean_var_name = "xmean_fun",
                                    K0 = 6,
                                    sigma_eps = 0.05,
                                    trunc_eps = 4,
                                    trunc_xi = 4,
                                    L_rff = 128L,
                                    n_eval = 5000,
                                    mc_n = 20000L,
                                    seed = NULL) {

  mark_mode <- match.arg(mark_mode)

  if (!is.null(seed)) set.seed(seed)


  w_u <- trapz_weights(u_grid)
  beta_u <- w_center(beta_u, w_u)
  beta_w <- beta_u * w_u


  a0 <- calibrate_a0_functional(
    beta_w = beta_w,
    Phi_true = Phi_true,
    target_N = target_N,
    beta_xyz = beta_xyz,
    K0 = K0,
    sigma_eps = sigma_eps,
    trunc_eps = trunc_eps,
    trunc_xi = trunc_xi,
    mc_n = mc_n,
    seed = if (!is.null(seed)) seed + 11L else NULL
  )

  lmax <- lmax_bound_functional(
    a0 = a0,
    beta_w = beta_w,
    Phi_true = Phi_true,
    beta_xyz = beta_xyz,
    trunc_xi = trunc_xi,
    sigma_eps = sigma_eps,
    trunc_eps = trunc_eps
  )


  field <- make_rff_field(
    K0 = K0,
    range = range,
    nugget = nugget,
    L = L_rff,
    seed = if (!is.null(seed)) seed + 101L else NULL,
    trunc_xi = trunc_xi
  )


  sim <- simulate_one_pattern_box3(
    beta_w = beta_w,
    a0 = a0,
    lmax = lmax,
    field = field,
    Phi_true = Phi_true,
    beta_xyz = beta_xyz,
    sigma_eps = sigma_eps,
    trunc_eps = trunc_eps,
    seed = if (!is.null(seed)) seed + 1001L else NULL
  )

  coords <- sim$coords
  Xevent <- sim$Xcurves
  n <- nrow(coords)

  if (n < 15) {
    return(list(
      ok = FALSE,
      n = n,
      coef_xyz = c(betax = NA_real_, betay = NA_real_, betaz = NA_real_),
      coef_mean = NA_real_,
      mise_lambda = NA_real_,
      logrmse_lambda = NA_real_,
      relerr_I = NA_real_
    ))
  }


  xmean <- curve_mean_weighted(Xevent, u_grid)


  X_s3d_mean <- build_s3dpm_mean(coords, xmean, var_name = mean_var_name)


  rhs <- mean_var_name
  if (beta_xyz[1] != 0) rhs <- paste("x +", rhs)
  if (beta_xyz[2] != 0) rhs <- paste("y +", rhs)
  if (beta_xyz[3] != 0) rhs <- paste("z +", rhs)

  fit <- try(
    fit_funpois_fpca(
      X_s3d = X_s3d_mean,
      rhs_terms = rhs,
      interp.p = interp.p,
      mult_fit = mult_fit,
      seed_fit = if (!is.null(seed)) seed + 2001L else NULL,
      mark_mode = mark_mode,
      mark_oracle_fun = NULL,
      box3 = list(xrange = c(0,1), yrange = c(0,1), zrange = c(0,1))
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    return(list(
      ok = FALSE,
      n = n,
      coef_xyz = c(betax = NA_real_, betay = NA_real_, betaz = NA_real_),
      coef_mean = NA_real_,
      mise_lambda = NA_real_,
      logrmse_lambda = NA_real_,
      relerr_I = NA_real_
    ))
  }

  cc <- coef(fit$mod_global)

  bx_hat <- if ("x" %in% names(cc)) as.numeric(cc["x"]) else 0
  by_hat <- if ("y" %in% names(cc)) as.numeric(cc["y"]) else 0
  bz_hat <- if ("z" %in% names(cc)) as.numeric(cc["z"]) else 0
  bm_hat <- if (mean_var_name %in% names(cc)) as.numeric(cc[mean_var_name]) else NA_real_



  if (!is.null(seed)) set.seed(seed + 5001L)

  coords_eval <- cbind(runif(n_eval), runif(n_eval), runif(n_eval))
  Xeval <- simulate_X_curves(
    coords_eval, field, Phi_true,
    sigma_eps = sigma_eps,
    trunc_eps = trunc_eps
  )

  integ_eval <- as.vector(Xeval %*% beta_w)
  lin_eval <- beta_xyz[1]*coords_eval[,1] +
    beta_xyz[2]*coords_eval[,2] +
    beta_xyz[3]*coords_eval[,3]
  lam_true <- exp(a0 + integ_eval + lin_eval)

  xmean_eval <- curve_mean_weighted(Xeval, u_grid)

  df_eval <- data.frame(
    x = coords_eval[,1],
    y = coords_eval[,2],
    z = coords_eval[,3],
    xmean_fun = xmean_eval,
    check.names = FALSE
  )
  names(df_eval)[4] <- mean_var_name

  lam_hat <- as.numeric(predict(fit$mod_global, newdata = df_eval, type = "response"))

  logrmse_l <- sqrt(mean((log(lam_hat) - log(lam_true))^2))
  mise_l    <- mean((lam_hat - lam_true)^2)
  I_true    <- mean(lam_true)
  I_hat     <- mean(lam_hat)
  relI      <- (I_hat - I_true) / I_true

  list(
    ok = TRUE,
    n = n,
    coef_xyz = c(betax = bx_hat, betay = by_hat, betaz = bz_hat),
    coef_mean = bm_hat,
    mise_lambda = mise_l,
    logrmse_lambda = logrmse_l,
    relerr_I = relI
  )
}

run_plan4_mean_scalar_parallel <- function(out_old,
                                           u_grid = seq(0,1,length.out=101),
                                           K0 = 6,
                                           sigma_eps = 0.05,
                                           trunc_eps = 4,
                                           trunc_xi  = 4,
                                           L_rff = 128L,
                                           mc_n = 20000L,
                                           mult_fit = 30,
                                           n_eval = 5000,
                                           mark_mode = "idw",
                                           verbose = TRUE,
                                           future_packages = c("mgcv","splines","spatstat.geom")) {

  grid <- out_old$grid
  nsim <- out_old$nsim
  base_seed <- out_old$base_seed
  interp.p <- out_old$interp.p
  K_fpca <- out_old$K_fpca

  betas <- make_beta_scenarios(u_grid)
  beta_names <- names(betas)
  w_u <- trapz_weights(u_grid)
  betas_centered <- lapply(betas, function(b) w_center(b, w_u))
  Phi_true <- make_true_phi_fourier(u_grid, K0 = K0)

  tasks <- expand.grid(
    i = seq_len(nrow(grid)),
    bn = beta_names,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  worker <- function(row_idx) {
    i  <- tasks$i[row_idx]
    bn <- tasks$bn[row_idx]

    sc <- grid[i, , drop = FALSE]
    key <- sc$scenario_key

    rr <- as.numeric(sc$grf_range)
    ng <- as.numeric(sc$grf_nugget)
    Nt <- as.integer(sc$target_N)
    case_int <- as.character(sc$intensity_case)

    beta_u <- betas_centered[[bn]]
    beta_xyz <- beta_xyz_from_case(case_int)

    ok_rep <- logical(nsim)
    n_obs  <- integer(nsim)
    mise_v <- rep(NA_real_, nsim)
    logrmse_v <- rep(NA_real_, nsim)
    relI_v <- rep(NA_real_, nsim)
    coef_xyz_v <- matrix(NA_real_, nsim, 3, dimnames = list(NULL, c("betax","betay","betaz")))
    coef_mean_v <- rep(NA_real_, nsim)

    beta_idx <- match(bn, beta_names)

    for (s in seq_len(nsim)) {
      seed_s <- base_seed +
        1000000L * as.integer(beta_idx) +
        100000L  * as.integer(i) +
        1000L    * as.integer(s)

      out <- run_one_rep_mean_scalar(
        beta_u = beta_u,
        u_grid = u_grid,
        Phi_true = Phi_true,
        range = rr,
        nugget = ng,
        target_N = Nt,
        beta_xyz = beta_xyz,
        interp.p = interp.p,
        mult_fit = mult_fit,
        mark_mode = mark_mode,
        K0 = K0,
        sigma_eps = sigma_eps,
        trunc_eps = trunc_eps,
        trunc_xi = trunc_xi,
        L_rff = L_rff,
        n_eval = n_eval,
        mc_n = mc_n,
        seed = seed_s
      )

      ok_rep[s] <- isTRUE(out$ok)
      n_obs[s]  <- out$n %||% NA_integer_

      if (isTRUE(out$ok)) {
        mise_v[s]      <- out$mise_lambda
        logrmse_v[s]   <- out$logrmse_lambda
        relI_v[s]      <- out$relerr_I
        coef_xyz_v[s,] <- out$coef_xyz
        coef_mean_v[s] <- out$coef_mean
      }
    }

    scenario_obj <- list(
      beta_name = bn,
      scenario_key = key,
      params = list(
        range = rr,
        nugget = ng,
        target_N = Nt,
        interp.p = interp.p,
        K_fpca = K_fpca,
        intensity_case = case_int,
        mark_mode = mark_mode
      ),
      ok = ok_rep,
      n_obs = n_obs,
      mean_scalar = list(
        mise_lambda = mise_v,
        logrmse_lambda = logrmse_v,
        relerr_I = relI_v,
        coef_xyz_hat = coef_xyz_v,
        coef_mean_hat = coef_mean_v
      )
    )

    list(bn = bn, key = key, scenario_obj = scenario_obj)
  }

  if (isTRUE(verbose)) {
    cat("\n[run_plan4_mean_scalar_parallel] Parallel tasks (scenario x beta): ",
        nrow(tasks), "\n", sep = "")
  }

  outs <- future.apply::future_lapply(
    X = seq_len(nrow(tasks)),
    FUN = worker,
    future.seed = TRUE,
    future.packages = future_packages
  )

  results <- setNames(vector("list", length(beta_names)), beta_names)
  for (bn in beta_names) {
    results[[bn]] <- setNames(vector("list", nrow(grid)), grid$scenario_key)
  }

  for (o in outs) {
    results[[o$bn]][[o$key]] <- o$scenario_obj
  }

  list(
    grid = grid,
    nsim = nsim,
    base_seed = base_seed,
    interp.p = interp.p,
    K_fpca = K_fpca,
    results = results
  )
}

merge_plan4_three_models <- function(out_old, out_mean) {
  out_new <- out_old

  for (bn in names(out_old$results)) {
    for (key in names(out_old$results[[bn]])) {
      out_new$results[[bn]][[key]]$mean_scalar <- out_mean$results[[bn]][[key]]$mean_scalar
      out_new$results[[bn]][[key]]$ok_mean_scalar <- out_mean$results[[bn]][[key]]$ok
      out_new$results[[bn]][[key]]$n_obs_mean_scalar <- out_mean$results[[bn]][[key]]$n_obs
    }
  }

  out_new
}

summarise_plan4_three_models <- function(out_3m) {
  rows <- list()
  kk <- 1L

  for (bn in names(out_3m$results)) {
    for (key in names(out_3m$results[[bn]])) {
      obj <- out_3m$results[[bn]][[key]]

      ok_old  <- obj$ok
      ok_mean <- obj$ok_mean_scalar

      ok_common <- ok_old & ok_mean
      n_common <- sum(ok_common, na.rm = TRUE)

      rows[[kk]] <- data.frame(
        beta_name = bn,
        scenario_key = key,
        intensity_case = obj$params$intensity_case %||% NA_character_,
        n_common = n_common,

        mise_correct = mean(obj$correct$mise_lambda[ok_common], na.rm = TRUE),
        mise_coords_only = mean(obj$misspec$mise_lambda[ok_common], na.rm = TRUE),
        mise_coords_mean = mean(obj$mean_scalar$mise_lambda[ok_common], na.rm = TRUE),

        logrmse_correct = mean(obj$correct$logrmse_lambda[ok_common], na.rm = TRUE),
        logrmse_coords_only = mean(obj$misspec$logrmse_lambda[ok_common], na.rm = TRUE),
        logrmse_coords_mean = mean(obj$mean_scalar$logrmse_lambda[ok_common], na.rm = TRUE),

        relI_correct = mean(obj$correct$relerr_I[ok_common], na.rm = TRUE),
        relI_coords_only = mean(obj$misspec$relerr_I[ok_common], na.rm = TRUE),
        relI_coords_mean = mean(obj$mean_scalar$relerr_I[ok_common], na.rm = TRUE),

        rmse_beta_correct = mean(obj$correct$rmse_beta[ok_common], na.rm = TRUE),
        ise_beta_correct = mean(obj$correct$ise_beta[ok_common], na.rm = TRUE),

        coef_x_correct = mean(obj$correct$coef_xyz_hat[ok_common, "betax"], na.rm = TRUE),
        coef_y_correct = mean(obj$correct$coef_xyz_hat[ok_common, "betay"], na.rm = TRUE),
        coef_z_correct = mean(obj$correct$coef_xyz_hat[ok_common, "betaz"], na.rm = TRUE),

        coef_x_coords_only = mean(obj$misspec$coef_xyz_hat[ok_common, "betax"], na.rm = TRUE),
        coef_y_coords_only = mean(obj$misspec$coef_xyz_hat[ok_common, "betay"], na.rm = TRUE),
        coef_z_coords_only = mean(obj$misspec$coef_xyz_hat[ok_common, "betaz"], na.rm = TRUE),

        coef_x_coords_mean = mean(obj$mean_scalar$coef_xyz_hat[ok_common, "betax"], na.rm = TRUE),
        coef_y_coords_mean = mean(obj$mean_scalar$coef_xyz_hat[ok_common, "betay"], na.rm = TRUE),
        coef_z_coords_mean = mean(obj$mean_scalar$coef_xyz_hat[ok_common, "betaz"], na.rm = TRUE),

        coef_mean_scalar = mean(obj$mean_scalar$coef_mean_hat[ok_common], na.rm = TRUE),

        stringsAsFactors = FALSE
      )

      kk <- kk + 1L
    }
  }

  do.call(rbind, rows)
}

build_plan4_relmise_targetN <- function(out_old, out_mean,
                                        beta_agg = c("median", "mean")) {
  beta_agg <- match.arg(beta_agg)
  eps <- 1e-12
  rows <- list()
  k <- 0L

  beta_names <- names(out_old$results)

  for (bn in beta_names) {
    keys <- names(out_old$results[[bn]])

    for (key in keys) {
      obj_old  <- out_old$results[[bn]][[key]]
      obj_mean <- out_mean$results[[bn]][[key]]

      ok_common <- as.logical(obj_old$ok) & as.logical(obj_mean$ok)
      if (!any(ok_common)) next

      target_N <- as.numeric(obj_old$params$target_N)
      denom <- max(target_N^2, eps)

      scenario_class <- sub("__.*$", "", key)
      intensity_case <- obj_old$params$intensity_case %||% NA_character_

      reps <- which(ok_common)


      k <- k + 1L
      rows[[k]] <- data.frame(
        beta_name = bn,
        scenario_key = key,
        scenario_class = scenario_class,
        intensity_case = intensity_case,
        model = "Solo coordinate",
        rep = reps,
        target_N = target_N,
        rMISE = obj_old$misspec$mise_lambda[ok_common] / denom,
        stringsAsFactors = FALSE
      )


      k <- k + 1L
      rows[[k]] <- data.frame(
        beta_name = bn,
        scenario_key = key,
        scenario_class = scenario_class,
        intensity_case = intensity_case,
        model = "Coordinate + media covariata",
        rep = reps,
        target_N = target_N,
        rMISE = obj_mean$mean_scalar$mise_lambda[ok_common] / denom,
        stringsAsFactors = FALSE
      )


      k <- k + 1L
      rows[[k]] <- data.frame(
        beta_name = bn,
        scenario_key = key,
        scenario_class = scenario_class,
        intensity_case = intensity_case,
        model = "Coordinate + covariata funzionale",
        rep = reps,
        target_N = target_N,
        rMISE = obj_old$correct$mise_lambda[ok_common] / denom,
        stringsAsFactors = FALSE
      )
    }
  }

  df_rel_rep <- bind_rows(rows)

  df_rel_beta <- df_rel_rep %>%
    group_by(beta_name, scenario_class, intensity_case, model) %>%
    summarise(
      rMISE_beta = mean(rMISE, na.rm = TRUE),
      .groups = "drop"
    )

  if (beta_agg == "median") {
    df_rel_heat <- df_rel_beta %>%
      group_by(scenario_class, intensity_case, model) %>%
      summarise(
        rMISE_value = median(rMISE_beta, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    df_rel_heat <- df_rel_beta %>%
      group_by(scenario_class, intensity_case, model) %>%
      summarise(
        rMISE_value = mean(rMISE_beta, na.rm = TRUE),
        .groups = "drop"
      )
  }

  df_rel_heat %>%
    mutate(
      scenario_class = factor(scenario_class,
                              levels = c("worst", "intermediate", "best")),
      intensity_case = factor(intensity_case,
                              levels = c("x_fun", "xyz_fun"),
                              labels = c("x + funzione", "x, y, z + funzione")),
      model = factor(model,
                     levels = c("Solo coordinate",
                                "Coordinate + media covariata",
                                "Coordinate + covariata funzionale"))
    )
}

plot_plan4_relmise_heatmap_targetN <- function(df_rel_heat) {

  facet_labs <- c(
    "x + funzione" = "lambda(s) * ',' ~~ s == x",
    "x, y, z + funzione" = "lambda(s) * ',' ~~ s == '(x,y,z)'"
  )

  ggplot(df_rel_heat, aes(x = scenario_class, y = model, fill = rMISE_value)) +
    geom_tile(color = "grey80", linewidth = 0.3) +
    geom_text(
      aes(label = sprintf("%.3f", rMISE_value)),
      color = "black",
      fontface = "bold",
      size = 10
    ) +
    facet_wrap(
      ~ intensity_case,
      labeller = as_labeller(facet_labs, label_parsed)
    ) +
    scale_y_discrete(
      labels = c(
        expression(hat(lambda)[coord](s)),
        expression(hat(lambda)[mean](s)),
        expression(hat(lambda)[fun](s))
      )
    ) +
    scale_fill_gradient(
      low = "white",
      high = "red",
      name = expression(rMISE[lambda(s)])
    ) +
    scale_x_discrete(
      labels = c(
        "worst" = "Worst",
        "intermediate" = "Intermediate",
        "best" = "Best"
      )
    ) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 28, face = "bold"),
      axis.title = element_text(size = 32),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 28),
      legend.title = element_text(size = 32),
      legend.text = element_text(size = 30),
      legend.key.height = unit(35, "mm"),
      legend.key.width = unit(6, "mm")
    ) +
    guides(
      fill = guide_colorbar(
        barheight = unit(200, "mm"),
        barwidth = unit(10, "mm"),
        title.position = "top",
        title.hjust = 0.5
      )
    )
}

build_plan4_coord_rmse <- function(out_old, out_mean,
                                   beta_agg = c("median", "mean")) {
  beta_agg <- match.arg(beta_agg)

  rows <- list()
  k <- 0L

  beta_names <- names(out_old$results)

  for (bn in beta_names) {
    keys <- names(out_old$results[[bn]])

    for (key in keys) {
      obj_old  <- out_old$results[[bn]][[key]]
      obj_mean <- out_mean$results[[bn]][[key]]

      ok_common <- as.logical(obj_old$ok) & as.logical(obj_mean$ok)
      if (!any(ok_common)) next

      intensity_case <- obj_old$params$intensity_case %||% NA_character_
      scenario_class <- sub("__.*$", "", key)

      beta_true <- beta_xyz_from_case(intensity_case)
      names(beta_true) <- c("betax", "betay", "betaz")


      params_use <- if (identical(intensity_case, "x_fun")) {
        "betax"
      } else if (identical(intensity_case, "xyz_fun")) {
        c("betax", "betay", "betaz")
      } else {
        stop("Unknown intensity_case: ", intensity_case)
      }


      coef_coords_only <- obj_old$misspec$coef_xyz_hat[ok_common, , drop = FALSE]
      coef_coords_mean <- obj_mean$mean_scalar$coef_xyz_hat[ok_common, , drop = FALSE]
      coef_functional  <- obj_old$correct$coef_xyz_hat[ok_common, , drop = FALSE]

      for (pp in params_use) {
        true_val <- as.numeric(beta_true[pp])


        rmse_coords_only <- sqrt(mean((coef_coords_only[, pp] - true_val)^2, na.rm = TRUE))
        rmse_coords_mean <- sqrt(mean((coef_coords_mean[, pp] - true_val)^2, na.rm = TRUE))
        rmse_functional  <- sqrt(mean((coef_functional[, pp]  - true_val)^2, na.rm = TRUE))

        k <- k + 1L
        rows[[k]] <- data.frame(
          beta_name = bn,
          scenario_key = key,
          scenario_class = scenario_class,
          intensity_case = intensity_case,
          parameter = pp,
          model = "Solo coordinate",
          RMSE = rmse_coords_only,
          stringsAsFactors = FALSE
        )

        k <- k + 1L
        rows[[k]] <- data.frame(
          beta_name = bn,
          scenario_key = key,
          scenario_class = scenario_class,
          intensity_case = intensity_case,
          parameter = pp,
          model = "Coordinate + media covariata",
          RMSE = rmse_coords_mean,
          stringsAsFactors = FALSE
        )

        k <- k + 1L
        rows[[k]] <- data.frame(
          beta_name = bn,
          scenario_key = key,
          scenario_class = scenario_class,
          intensity_case = intensity_case,
          parameter = pp,
          model = "Coordinate + covariata funzionale",
          RMSE = rmse_functional,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  df_rmse_beta <- bind_rows(rows)


  if (beta_agg == "median") {
    df_rmse_heat <- df_rmse_beta %>%
      group_by(scenario_class, intensity_case, parameter, model) %>%
      summarise(
        RMSE_value = median(RMSE, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    df_rmse_heat <- df_rmse_beta %>%
      group_by(scenario_class, intensity_case, parameter, model) %>%
      summarise(
        RMSE_value = mean(RMSE, na.rm = TRUE),
        .groups = "drop"
      )
  }

  df_rmse_heat <- df_rmse_heat %>%
    mutate(
      scenario_class = factor(
        scenario_class,
        levels = c("worst", "intermediate", "best")
      ),
      intensity_case = factor(
        intensity_case,
        levels = c("x_fun", "xyz_fun"),
        labels = c("x + funzione", "x, y, z + funzione")
      ),
      parameter = factor(
        parameter,
        levels = c("betax", "betay", "betaz"),
        labels = c("x", "y", "z")
      ),
      model = factor(
        model,
        levels = c(
          "Solo coordinate",
          "Coordinate + media covariata",
          "Coordinate + covariata funzionale"
        )
      )
    )

  list(
    df_rmse_beta = df_rmse_beta,
    df_rmse_heat = df_rmse_heat
  )
}

plot_plan4_coord_rmse_heatmap_compact <- function(df_rmse_heat) {

  df_plot <- df_rmse_heat %>%
    filter(
      (intensity_case == "x + funzione" & parameter == "x") |
        (intensity_case == "x, y, z + funzione" & parameter %in% c("x", "y", "z"))
    )

  facet_labs_intensity <- c(
    "x + funzione" = "lambda(s) * ',' ~~ s == x",
    "x, y, z + funzione" = "lambda(s) * ',' ~~ s == '(x,y,z)'"
  )

  facet_labs_param <- c(
    "x" = "theta[x]",
    "y" = "theta[y]",
    "z" = "theta[z]"
  )

  ggplot(df_plot, aes(x = scenario_class, y = model, fill = RMSE_value)) +
    geom_tile(color = "grey80", linewidth = 0.3) +
    geom_text(
      aes(label = sprintf("%.3f", RMSE_value)),
      color = "black",
      fontface = "bold",
      size = 10
    ) +
    facet_grid(
      parameter ~ intensity_case,
      drop = TRUE,
      labeller = labeller(
        parameter = as_labeller(facet_labs_param, label_parsed),
        intensity_case = as_labeller(facet_labs_intensity, label_parsed)
      )
    ) +
    scale_y_discrete(
      labels = c(
        expression(hat(lambda)[coord](s)),
        expression(hat(lambda)[mean](s)),
        expression(hat(lambda)[fun](s))
      )
    ) +
    scale_fill_gradient(
      low = "white",
      high = "red",
      name = expression(RMSE[bold(theta)[Z]])
    ) +
    scale_x_discrete(
      labels = c(
        "worst" = "Worst",
        "intermediate" = "Intermediate",
        "best" = "Best"
      )
    ) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 28, face = "bold"),
      axis.title = element_text(size = 32),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 28),
      legend.title = element_text(size = 32),
      legend.text = element_text(size = 30),
      legend.key.height = unit(25, "mm"),
      legend.key.width = unit(6, "mm")
    ) +
    guides(
      fill = guide_colorbar(
        barheight = unit(200, "mm"),
        barwidth = unit(10, "mm"),
        title.position = "top",
        title.hjust = 0.5
      )
    )
}

get_scn <- function(out, beta_name, scenario_key) {
  out$results[[beta_name]][[scenario_key]]
}

collect_rmse_rep <- function(out, scenario_key) {
  do.call(rbind, lapply(beta_levels, function(bn) {
    sc <- get_scn(out, bn, scenario_key)
    if (is.null(sc)) return(NULL)
    ok <- as.logical(sc$ok)
    data.frame(
      beta_name = bn,
      rep = seq_along(sc$rmse_beta),
      ok = ok,
      rmse_beta = sc$rmse_beta,
      ise_beta  = sc$ise_beta,
      mise_lambda = sc$mise_lambda,
      n_obs = sc$n_obs,
      stringsAsFactors = FALSE
    )
  }))
}

build_relmise_rep <- function(out) {
  rows <- list()
  k <- 0L
  eps <- 1e-12

  for (bn in beta_levels) {
    scn_list <- out$results[[bn]]
    for (key in names(scn_list)) {
      sc <- scn_list[[key]]
      if (is.null(sc)) next

      ok <- as.logical(sc$ok) &
        is.finite(sc$mise_lambda) &
        is.finite(sc$I_true) &
        (sc$I_true > 0)

      if (!any(ok)) next

      rel_mise <- sc$mise_lambda[ok] / pmax(sc$I_true[ok]^2, eps)
      rel_rmse <- sc$rmse_lambda[ok] / pmax(sc$I_true[ok], eps)
      k <- k + 1L

      rows[[k]] <- data.frame(
        beta_name = bn,
        scenario_key = key,
        grf_range  = as.numeric(sc$params$range),
        grf_nugget = as.numeric(sc$params$nugget),
        target_N   = as.integer(sc$params$target_N),
        rep = which(ok),
        n_obs = sc$n_obs[ok],
        rmse_beta = sc$rmse_beta[ok],
        mise_lambda = sc$mise_lambda[ok],
        I_true = sc$I_true[ok],
        nMISE_lambda = rel_mise,
        rRMSE_lambda = rel_rmse,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}

build_rmsebeta_rep_main <- function(out) {
  rows <- list(); k <- 0L
  for (bn in beta_levels) {
    scn_list <- out$results[[bn]]
    for (key in names(scn_list)) {
      sc <- scn_list[[key]]
      if (is.null(sc)) next

      ok <- as.logical(sc$ok) & is.finite(sc$rmse_beta)
      if (!any(ok)) next

      k <- k + 1L
      rows[[k]] <- data.frame(
        beta_name = bn,
        scenario_key = key,
        grf_range  = as.numeric(sc$params$range),
        grf_nugget = as.numeric(sc$params$nugget),
        target_N   = as.integer(sc$params$target_N),
        rep = which(ok),
        rmse_beta = sc$rmse_beta[ok],
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

build_rmsebeta_rep_sens <- function(out) {
  rows <- list(); k <- 0L
  for (bn in beta_levels) {
    scn_list <- out$results[[bn]]
    for (key in names(scn_list)) {
      sc <- scn_list[[key]]
      if (is.null(sc)) next

      sc_class <- if (grepl("^best__", key)) "best" else if (grepl("^worst__", key)) "worst" else NA_character_
      if (is.na(sc_class)) next

      ok <- as.logical(sc$ok) & is.finite(sc$rmse_beta)
      if (!any(ok)) next

      k <- k + 1L
      rows[[k]] <- data.frame(
        beta_name = bn,
        scenario_class = sc_class,
        interp.p = as.numeric(sc$params$interp.p),
        K_fpca   = as.integer(sc$params$K_fpca),
        rep = which(ok),
        rmse_beta = sc$rmse_beta[ok],
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

build_relmise_rep_sens <- function(out) {
  rows <- list(); k <- 0L; eps <- 1e-12
  for (bn in beta_levels) {
    scn_list <- out$results[[bn]]
    for (key in names(scn_list)) {
      sc <- scn_list[[key]]
      if (is.null(sc)) next
      ok <- as.logical(sc$ok) & is.finite(sc$mise_lambda) & is.finite(sc$I_true) & (sc$I_true > 0)
      if (!any(ok)) next


      sc_class <- if (grepl("^best__", key)) "best" else if (grepl("^worst__", key)) "worst" else NA
      if (is.na(sc_class)) next

      k <- k + 1L
      rows[[k]] <- data.frame(
        beta_name = bn,
        scenario_class = sc_class,
        interp.p = as.numeric(sc$params$interp.p),
        K_fpca   = as.integer(sc$params$K_fpca),
        nMISE_lambda = sc$mise_lambda[ok] / pmax(sc$I_true[ok]^2, eps),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

build_nMISE_rep_stress <- function(out) {
  rows <- list(); k <- 0L
  eps <- 1e-12

  for (bn in beta_levels) {
    scn_list <- out$results[[bn]]
    for (key in names(scn_list)) {
      sc <- scn_list[[key]]
      if (is.null(sc)) next

      ok <- as.logical(sc$ok) &
        is.finite(sc$mise_lambda) &
        is.finite(sc$I_true) &
        (sc$I_true > 0)

      if (!any(ok)) next

      scenario_class <- if (grepl("^best__", key)) "best" else if (grepl("^worst__", key)) "worst" else NA_character_
      if (is.na(scenario_class)) next

      k <- k + 1L
      rows[[k]] <- data.frame(
        beta_name = bn,
        scenario_key = key,
        scenario_class = scenario_class,
        fpca_basis_mode = as.character(sc$params$fpca_basis_mode),
        mark_mode = as.character(sc$params$mark_mode),
        rep = which(ok),
        rmse_beta = sc$rmse_beta[ok],
        logrmse_lambda = sc$logrmse_lambda[ok],
        nMISE_lambda = sc$mise_lambda[ok] / pmax(sc$I_true[ok]^2, eps),
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}

make_delta_from_baseline_nMISE <- function(df) {
  out <- df
  for (sc in unique(df$scenario_class)) {
    b <- df %>% filter(scenario_class == sc,
                       fpca_basis_mode == "event",
                       mark_mode == "idw")
    stopifnot(nrow(b) == 1)

    out <- out %>%
      group_by(scenario_class) %>%
      mutate(
        d_rmse_beta   = rmse_beta_mean - b$rmse_beta_mean,
        d_logrmse     = logrmse_lambda_mean - b$logrmse_lambda_mean,
        d_nMISE       = nMISE_mean - b$nMISE_mean
      ) %>%
      ungroup()
  }
  out
}

parse_cm_key <- function(key) {

  sc_class <- sub("__.*", "", key)
  case <- sub(".*__case=", "", key)
  list(scenario_class = sc_class, intensity_case = case)
}

build_cm_rep_df <- function(out_cm) {
  res <- out_cm$results
  rows <- list()
  k <- 0L
  eps <- 1e-12

  for (bn in names(res)) {
    scn_list <- res[[bn]]
    for (key in names(scn_list)) {
      sc <- scn_list[[key]]
      if (is.null(sc)) next

      meta <- parse_cm_key(key)
      sc_class <- meta$scenario_class
      case_int <- meta$intensity_case


      ok <- as.logical(sc$ok)
      nrep <- length(ok)
      if (!nrep) next



      Nt <- as.numeric(sc$params$target_N)


      rmse_beta_c <- sc$correct$rmse_beta
      ise_beta_c  <- sc$correct$ise_beta
      mise_c      <- sc$correct$mise_lambda
      logrmse_c   <- sc$correct$logrmse_lambda
      relI_c      <- sc$correct$relerr_I
      coef_c      <- sc$correct$coef_xyz_hat


      mise_m      <- sc$misspec$mise_lambda
      logrmse_m   <- sc$misspec$logrmse_lambda
      relI_m      <- sc$misspec$relerr_I
      coef_m      <- sc$misspec$coef_xyz_hat



      nMISE_c <- mise_c / pmax(Nt^2, eps)
      nMISE_m <- mise_m / pmax(Nt^2, eps)


      for (r in seq_len(nrep)) {
        if (!ok[r]) next

        k <- k + 1L
        rows[[k]] <- data.frame(
          beta_name = bn,
          scenario_key = key,
          scenario_class = sc_class,
          intensity_case = case_int,
          target_N = Nt,
          rep = r,
          model = "correct",
          rmse_beta = rmse_beta_c[r],
          ise_beta  = ise_beta_c[r],
          mise_lambda = mise_c[r],
          nMISE_lambda = nMISE_c[r],
          logrmse_lambda = logrmse_c[r],
          relerr_I = relI_c[r],
          betax_hat = coef_c[r, "betax"],
          betay_hat = coef_c[r, "betay"],
          betaz_hat = coef_c[r, "betaz"],
          stringsAsFactors = FALSE
        )

        k <- k + 1L
        rows[[k]] <- data.frame(
          beta_name = bn,
          scenario_key = key,
          scenario_class = sc_class,
          intensity_case = case_int,
          target_N = Nt,
          rep = r,
          model = "misspec",
          rmse_beta = NA_real_,
          ise_beta  = NA_real_,
          mise_lambda = mise_m[r],
          nMISE_lambda = nMISE_m[r],
          logrmse_lambda = logrmse_m[r],
          relerr_I = relI_m[r],
          betax_hat = coef_m[r, "betax"],
          betay_hat = coef_m[r, "betay"],
          betaz_hat = coef_m[r, "betaz"],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

true_coef <- function(case) {
  if (case == "x_fun") c(betax=1, betay=0, betaz=0) else c(betax=1, betay=1, betaz=1)
}

# -----------------------------------------------------------------------------
# Figure and table builders
# -----------------------------------------------------------------------------

beta_shape_names <- c("monotone", "single", "double", "signflip", "oscillate", "narrow")

plot_beta_panel_six <- function(simulation_output, scenario_key, band_multiplier = 1,
                                main_suffix = NULL) {
  previous <- par(mfrow = c(3, 2), mar = c(4.0, 4.2, 2.2, 0.8))
  on.exit(par(previous), add = TRUE)
  for (beta_name in beta_shape_names) {
    scenario <- simulation_output$results[[beta_name]][[scenario_key]]
    if (is.null(scenario)) {
      plot.new()
      next
    }
    title_text <- beta_name
    if (!is.null(main_suffix)) title_text <- paste(beta_name, main_suffix)
    plot_beta_bands(scenario, main = title_text, lwd = 2)
  }
}

build_main_heatmap_data <- function(main_output) {
  tab <- summarise_table(main_output)
  dplyr::summarise(
    dplyr::group_by(tab, grf_range, grf_nugget, target_N),
    rmse_beta = mean(rmse_beta_mean, na.rm = TRUE),
    relative_mise = mean(mise_lambda_mean / pmax(target_N^2, 1e-12), na.rm = TRUE),
    .groups = "drop"
  )
}

plot_main_heatmaps <- function(main_output) {
  data <- build_main_heatmap_data(main_output)
  left <- ggplot2::ggplot(data, ggplot2::aes(
    x = factor(grf_nugget), y = factor(grf_range), fill = rmse_beta
  )) +
    ggplot2::geom_tile(colour = "grey85") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", rmse_beta)), size = 3.1) +
    ggplot2::facet_wrap(~target_N, labeller = ggplot2::label_bquote(rows = n == .(target_N))) +
    ggplot2::labs(x = expression(N[g]), y = expression(R[g]), fill = expression(RMSE[theta[F]])) +
    ggplot2::theme_classic(base_size = 12)
  right <- ggplot2::ggplot(data, ggplot2::aes(
    x = factor(grf_nugget), y = factor(grf_range), fill = relative_mise
  )) +
    ggplot2::geom_tile(colour = "grey85") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", relative_mise)), size = 3.1) +
    ggplot2::facet_wrap(~target_N, labeller = ggplot2::label_bquote(rows = n == .(target_N))) +
    ggplot2::labs(x = expression(N[g]), y = expression(R[g]), fill = "rMISE") +
    ggplot2::theme_classic(base_size = 12)
  list(beta_rmse = left, intensity_rmise = right, data = data)
}

build_sensitivity_pve_data <- function(sensitivity_output) {
  tab <- summarise_table(sensitivity_output)
  dplyr::summarise(
    dplyr::group_by(tab, scenario_key, interp.p, K_fpca),
    pve_cum = mean(pve_cum_last_mean, na.rm = TRUE),
    .groups = "drop"
  )
}

plot_sensitivity_pve_heatmap <- function(sensitivity_output) {
  data <- build_sensitivity_pve_data(sensitivity_output)
  ggplot2::ggplot(data, ggplot2::aes(
    x = factor(interp.p), y = factor(K_fpca), fill = pve_cum
  )) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", pve_cum)), size = 3.4) +
    ggplot2::facet_wrap(~scenario_key) +
    ggplot2::scale_fill_gradient(low = "white", high = "firebrick") +
    ggplot2::labs(x = expression(p[IDW]), y = "Retained FPCA components", fill = expression(PVE[cum](K))) +
    ggplot2::theme_classic(base_size = 12)
}

plot_sensitivity_heatmaps <- function(sensitivity_output) {
  tab <- summarise_table(sensitivity_output)
  beta_plot <- ggplot2::ggplot(tab, ggplot2::aes(
    x = factor(interp.p), y = factor(K_fpca), fill = rmse_beta_mean
  )) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", rmse_beta_mean)), size = 3.0) +
    ggplot2::facet_wrap(~scenario_key) +
    ggplot2::labs(x = expression(p[IDW]), y = "Retained FPCA components", fill = expression(RMSE[theta[F]])) +
    ggplot2::theme_classic(base_size = 12)
  intensity_plot <- ggplot2::ggplot(tab, ggplot2::aes(
    x = factor(interp.p), y = factor(K_fpca), fill = mise_lambda_mean / pmax(target_N^2, 1e-12)
  )) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", mise_lambda_mean / pmax(target_N^2, 1e-12))), size = 3.0) +
    ggplot2::facet_wrap(~scenario_key) +
    ggplot2::labs(x = expression(p[IDW]), y = "Retained FPCA components", fill = "rMISE") +
    ggplot2::theme_classic(base_size = 12)
  list(beta_rmse = beta_plot, intensity_rmise = intensity_plot)
}

plot_illustrative_simulation <- function(seed = 2026L, target_N = 200L, K0 = 6L,
                                         beta_name = "single") {
  set.seed(seed)
  u_grid <- seq(0, 1, length.out = 101)
  w_u <- trapz_weights(u_grid)
  phi_true <- make_true_phi_fourier(u_grid, K0 = K0)
  beta_true <- w_center(make_beta_scenarios(u_grid)[[beta_name]], w_u)
  beta_w <- beta_true * w_u
  field <- make_rff_field(K0 = K0, range = 0.25, nugget = 0.05, seed = seed + 1L)
  intercept <- calibrate_a0_functional(beta_w, phi_true, target_N = target_N,
                                        K0 = K0, seed = seed + 2L)
  lambda_max <- lmax_bound_functional(intercept, beta_w, phi_true, K0 = K0)
  simulated <- simulate_one_pattern_box3(beta_w, intercept, lambda_max, field, phi_true,
                                         seed = seed + 3L)
  curve_data <- data.frame(
    u = rep(u_grid, nrow(simulated$Xcurves)),
    value = as.vector(t(simulated$Xcurves)),
    curve = rep(seq_len(nrow(simulated$Xcurves)), each = length(u_grid))
  )
  mean_curve <- data.frame(u = u_grid, value = colMeans(simulated$Xcurves))
  spatial_plot <- plot3D::scatter3D(
    x = simulated$coords[, 1], y = simulated$coords[, 2], z = simulated$coords[, 3],
    pch = 19, cex = 0.5, xlab = "x", ylab = "y", zlab = "z", main = ""
  )
  curve_plot <- ggplot2::ggplot(curve_data, ggplot2::aes(u, value, group = curve)) +
    ggplot2::geom_line(alpha = 0.15, colour = "steelblue") +
    ggplot2::geom_line(data = mean_curve, linewidth = 0.8) +
    ggplot2::labs(x = expression(tau), y = expression(C(s, tau))) +
    ggplot2::theme_classic(base_size = 12)
  list(points = simulated$coords, curves = simulated$Xcurves, curve_plot = curve_plot,
       spatial_plot = spatial_plot)
}

write_simulation_outputs <- function(main_output, sensitivity_output, plan4_output,
                                     stress_output, directories = make_output_directories()) {
  selected <- pick_best_mid_worst(main_output)
  main_plots <- plot_main_heatmaps(main_output)
  sensitivity_plots <- plot_sensitivity_heatmaps(sensitivity_output)
  pve_plot <- plot_sensitivity_pve_heatmap(sensitivity_output)
  plan4_rmise <- plot_plan4_relmise_heatmap_targetN(
    build_plan4_relmise_targetN(plan4_output$functional_coordinate, plan4_output$scalar_mean)
  )
  plan4_coef <- plot_plan4_coord_rmse_heatmap_compact(
    build_plan4_coord_rmse(plan4_output$functional_coordinate, plan4_output$scalar_mean)
  )

  save_ggplot_pdf(main_plots$beta_rmse, file.path(directories$figures, "Figure_03_beta_rmse_heatmap.pdf"), 10, 4.5)
  save_ggplot_pdf(main_plots$intensity_rmise, file.path(directories$figures, "Figure_03_intensity_rmise_heatmap.pdf"), 10, 4.5)
  save_ggplot_pdf(sensitivity_plots$beta_rmse, file.path(directories$figures, "Figure_04_beta_sensitivity.pdf"), 9, 4.5)
  save_ggplot_pdf(sensitivity_plots$intensity_rmise, file.path(directories$figures, "Figure_04_intensity_sensitivity.pdf"), 9, 4.5)
  save_ggplot_pdf(pve_plot, file.path(directories$figures, "Figure_05_pve_heatmap.pdf"), 9, 4.5)
  save_ggplot_pdf(plan4_rmise, file.path(directories$figures, "Figure_06_rmise_comparison.pdf"), 9, 5)
  save_ggplot_pdf(plan4_coef, file.path(directories$figures, "Figure_07_coordinate_rmse.pdf"), 10, 6)
  save_base_plot_pdf(file.path(directories$figures, "Figure_02_beta_recovery_favourable.pdf"), 8, 8,
                     function() plot_beta_panel_six(main_output, selected$best$scenario_key[1], main_suffix = "Favourable"))
  save_base_plot_pdf(file.path(directories$figures, "Figure_02_beta_recovery_challenging.pdf"), 8, 8,
                     function() plot_beta_panel_six(main_output, selected$worst$scenario_key[1], main_suffix = "Challenging"))

  utils::write.csv(summarise_table(main_output), file.path(directories$tables, "main_simulation_summary.csv"), row.names = FALSE)
  utils::write.csv(summarise_table(sensitivity_output), file.path(directories$tables, "sensitivity_summary.csv"), row.names = FALSE)
  utils::write.csv(summarise_table(stress_output), file.path(directories$tables, "robustness_summary.csv"), row.names = FALSE)
  invisible(selected)
}

# -----------------------------------------------------------------------------
# Reproducible workflow
# -----------------------------------------------------------------------------

SIMULATION_SETTINGS <- list(
  seed = 1234L,
  nsim = 100L,
  u_grid = seq(0, 1, length.out = 101),
  K0 = 6L,
  K_fpca_main = 6L,
  K_fpca_sensitivity = c(2L, 4L, 6L, 8L, 10L),
  idw_powers = c(2, 4, 8, 16),
  mult_fit = 30L,
  n_eval = 5000L,
  mc_n = 20000L,
  workers = 1L
)

run_simulation_study <- function(settings = SIMULATION_SETTINGS,
                                 root = project_root(),
                                 overwrite = FALSE) {
  check_simulation_packages()
  directories <- make_output_directories(root)
  future::plan(future::multisession, workers = settings$workers)
  on.exit(future::plan(future::sequential), add = TRUE)

  main_file <- file.path(directories$results, "main_simulation.rds")
  sensitivity_file <- file.path(directories$results, "sensitivity_simulation.rds")
  comparison_file <- file.path(directories$results, "functional_coordinate_comparison.rds")
  scalar_file <- file.path(directories$results, "scalar_summary_comparison.rds")
  robustness_file <- file.path(directories$results, "robustness_simulation.rds")

  main_output <- if (file.exists(main_file) && !overwrite) {
    readRDS(main_file)
  } else {
    run_plan_parallel(
      grid = scenario_grid_main(), nsim = settings$nsim, base_seed = settings$seed,
      u_grid = settings$u_grid, K0 = settings$K0, mult_fit = settings$mult_fit,
      n_eval = settings$n_eval, mc_n = settings$mc_n, default_interp_p = 4,
      default_K_fpca = settings$K_fpca_main, save_dir = directories$results,
      save_name = basename(main_file)
    )
  }

  sensitivity_output <- if (file.exists(sensitivity_file) && !overwrite) {
    readRDS(sensitivity_file)
  } else {
    run_plan_parallel(
      grid = scenario_grid_joint_sens(interp_ps = settings$idw_powers,
                                      K_fpcas = settings$K_fpca_sensitivity),
      nsim = settings$nsim, base_seed = settings$seed + 100L,
      u_grid = settings$u_grid, K0 = settings$K0, mult_fit = settings$mult_fit,
      n_eval = settings$n_eval, mc_n = settings$mc_n, default_interp_p = 4,
      default_K_fpca = settings$K_fpca_main, save_dir = directories$results,
      save_name = basename(sensitivity_file)
    )
  }

  selected <- pick_best_mid_worst(main_output)
  functional_coordinate <- if (file.exists(comparison_file) && !overwrite) {
    readRDS(comparison_file)
  } else {
    result <- run_plan_coords_misspec_parallel(
      selected_triplet = selected, nsim = settings$nsim,
      base_seed = settings$seed + 200L, u_grid = settings$u_grid, K0 = settings$K0,
      mult_fit = settings$mult_fit, interp.p = 16, K_fpca = 4,
      n_eval = settings$n_eval, mc_n = settings$mc_n
    )
    saveRDS(result, comparison_file)
    result
  }

  scalar_mean <- if (file.exists(scalar_file) && !overwrite) {
    readRDS(scalar_file)
  } else {
    result <- run_plan4_mean_scalar_parallel(
      out_old = functional_coordinate, u_grid = settings$u_grid, K0 = settings$K0,
      mult_fit = settings$mult_fit, n_eval = settings$n_eval, mc_n = settings$mc_n,
      interp.p = 16, mark_mode = "idw"
    )
    saveRDS(result, scalar_file)
    result
  }

  robustness_output <- if (file.exists(robustness_file) && !overwrite) {
    readRDS(robustness_file)
  } else {
    result <- run_plan_parallel(
      grid = scenario_grid_stress(), nsim = settings$nsim, base_seed = settings$seed + 300L,
      u_grid = settings$u_grid, K0 = settings$K0, mult_fit = settings$mult_fit,
      n_eval = settings$n_eval, mc_n = settings$mc_n, default_interp_p = 16,
      default_K_fpca = 4, save_dir = directories$results,
      save_name = basename(robustness_file)
    )
    result
  }

  plan4_output <- list(functional_coordinate = functional_coordinate, scalar_mean = scalar_mean)
  write_simulation_outputs(main_output, sensitivity_output, plan4_output, robustness_output, directories)

  list(
    main = main_output,
    sensitivity = sensitivity_output,
    functional_coordinate = functional_coordinate,
    scalar_mean = scalar_mean,
    robustness = robustness_output,
    selected = selected,
    directories = directories
  )
}
