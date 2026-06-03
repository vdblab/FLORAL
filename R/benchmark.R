# ── Benchmarking utilities for FLORAL knockoff feature selection ──────────────
#
# Public functions
#   benchmark_knockoff()   – one simulation replicate, all methods
#   run_benchmark_sweep()  – iterate configs × replicates, collect results
#   plot_benchmark()       – FDR / power visualisation
#
# Internal helpers
#   .bk_fdr_power()        – compute FDR and power from selected + true sets
#   .bk_selected_from_names() – extract taxa indices from name vectors
# ─────────────────────────────────────────────────────────────────────────────


# ── Internal helpers ──────────────────────────────────────────────────────────

#' @noRd
.bk_fdr_power <- function(selected_idx, true_idx) {
  n_sel <- length(selected_idx)
  tp    <- length(intersect(selected_idx, true_idx))
  fp    <- n_sel - tp
  list(
    fdr       = if (n_sel == 0L) 0 else fp / n_sel,
    power     = tp / max(length(true_idx), 1L),
    n_sel     = n_sel,
    n_tp      = tp
  )
}

#' Convert "taxa3" style names → integer indices
#' @noRd
.bk_idx_from_names <- function(nms, prefix = "taxa") {
  if (is.null(nms) || length(nms) == 0L) return(integer(0L))
  as.integer(gsub(prefix, "", nms, fixed = TRUE))
}


# ── benchmark_knockoff() ──────────────────────────────────────────────────────

#' Run one benchmarking replicate for knockoff feature selection
#'
#' Simulates a Gaussian log-ratio dataset using \code{\link{simu}} and evaluates
#' FDR and power for several feature-selection methods.
#'
#' @section Methods evaluated:
#' \describe{
#'   \item{\code{"FLORAL-KN-VAE"}}{Log-ratio lasso with VAE-generated knockoffs
#'     (the FLORAL knockoff method). Controlled at \code{fdr_target}.}
#'   \item{\code{"FLORAL-KN-2order"}}{Log-ratio lasso with second-order model-X
#'     knockoffs (\code{knockoff::create.second_order}).}
#'   \item{\code{"FLORAL-lmin"}}{Log-ratio lasso with \code{ncv}-fold CV;
#'     features selected at \code{lambda.min}.}
#'   \item{\code{"FLORAL-l1se"}}{Same as above but at \code{lambda.1se}
#'     (more conservative).}
#'   \item{\code{"KN-fixed"}}{Standard lasso knockoff filter with fixed-X
#'     knockoffs (\code{knockoff::create.fixed}). Skipped unless \code{n > p}.}
#'   \item{\code{"KN-2order"}}{Standard lasso knockoff filter with second-order
#'     model-X knockoffs applied directly to log1p features.}
#' }
#'
#' @param n Integer. Sample size.
#' @param p Integer. Number of taxa (features).
#' @param rho Numeric in [0, 1). Correlation parameter for the taxa covariance
#'   matrix in \code{\link{simu}}.
#' @param weak Integer. Number of weak-signal taxa.
#' @param strong Integer. Number of strong-signal taxa.
#' @param weaksize Numeric. Effect size for weak signals.
#' @param strongsize Numeric. Effect size for strong signals.
#' @param pct.sparsity Numeric in [0, 1). Fraction of zero counts per sample.
#' @param fdr_target Numeric. Nominal FDR level for knockoff methods (default
#'   \code{0.2}).
#' @param ncv Integer. Number of CV folds for \code{FLORAL-lmin} and
#'   \code{FLORAL-l1se} (default \code{5}).
#' @param vae_params Named list of VAE hyperparameters forwarded to the
#'   \code{kn_method} list inside \code{\link{FLORAL}}.  Any element not
#'   provided falls back to the FLORAL default.  Useful elements:
#'   \code{latent_dim}, \code{lambda_kl}, \code{lambda_abun},
#'   \code{lambda_pres}, \code{gamma_full}, \code{gamma_swap},
#'   \code{lambda_moments}, \code{delta_corr}, \code{epochs},
#'   \code{batch_size}, \code{lr}, \code{weight_decay}.
#' @param methods Character vector.  Subset of methods to run.  Defaults to
#'   all six (see Details).
#' @param seed Integer or \code{NULL}.  Random seed passed to \code{set.seed}
#'   before data simulation.
#' @param progress Logical.  If \code{TRUE} (default \code{FALSE}), forwards
#'   progress output from model fitting.
#'
#' @return A \code{data.frame} with one row per method and columns:
#' \describe{
#'   \item{\code{method}}{Method name.}
#'   \item{\code{fdr_empirical}}{Realised FDP (false discovery proportion).}
#'   \item{\code{power}}{Proportion of true signals selected.}
#'   \item{\code{n_selected}}{Total features selected.}
#'   \item{\code{n_tp}}{True positives selected.}
#'   \item{\code{n_signal}}{Number of true signals in the simulation.}
#'   \item{\code{n}}{Sample size.}
#'   \item{\code{p}}{Number of taxa.}
#'   \item{\code{rho}}{Correlation parameter.}
#'   \item{\code{fdr_target}}{Nominal FDR (knockoff methods).}
#'   \item{\code{seed}}{Seed used.}
#'   \item{\code{error}}{Error message if a method failed; \code{NA} otherwise.}
#' }
#'
#' @examples
#' \dontrun{
#' # Single replicate with default settings
#' res <- benchmark_knockoff(n = 150, p = 50, seed = 1)
#'
#' # Vary VAE loss weights
#' res <- benchmark_knockoff(
#'   n = 150, p = 50, seed = 1,
#'   vae_params = list(gamma_full = 2, gamma_swap = 2, lambda_moments = 0.5),
#'   methods = c("FLORAL-KN-VAE", "FLORAL-KN-2order", "FLORAL-lmin")
#' )
#' }
#'
#' @seealso \code{\link{run_benchmark_sweep}}, \code{\link{plot_benchmark}}
#' @importFrom stats setNames
#' @export
benchmark_knockoff <- function(n            = 150,
                               p            = 50,
                               rho          = 0,
                               weak         = 2,
                               strong       = 4,
                               weaksize     = 0.125,
                               strongsize   = 0.25,
                               pct.sparsity = 0.5,
                               fdr_target   = 0.2,
                               ncv          = 5,
                               vae_params   = list(),
                               methods      = c("FLORAL-KN-VAE",
                                                "FLORAL-KN-2order",
                                                "FLORAL-lmin",
                                                "FLORAL-l1se",
                                                "KN-fixed",
                                                "KN-2order"),
                               seed         = NULL,
                               progress     = FALSE) {

  # ── 0. Setup ────────────────────────────────────────────────────────────────
  if (!is.null(seed)) set.seed(seed)

  # ── 1. Simulate data ────────────────────────────────────────────────────────
  dat <- simu(n            = n,
              p            = p,
              model        = "linear",
              weak         = weak,
              strong       = strong,
              weaksize     = weaksize,
              strongsize   = strongsize,
              pct.sparsity = pct.sparsity,
              rho          = rho)

  true_idx  <- dat$idx                     # integer vector, 1-based
  n_signal  <- length(true_idx)
  xcount    <- dat$xcount                  # raw counts  (n × p)
  xlog1p    <- dat$x                       # log1p counts (n × p)
  y         <- dat$y

  # ── 2. Build default VAE kn_method, then override with user vae_params ──────
  default_vae <- list(
    model          = "VAE",
    latent_dim     = 32L,
    lambda_kl      = 0,
    lambda_abun    = 0,
    lambda_pres    = 0,
    gamma_full     = 1,
    gamma_swap     = 1,
    lambda_moments = 0,
    delta_corr     = 0,
    sigma_list     = c(1, 2, 4, 8, 16, 32, 64, 128),
    epochs         = 100L,
    batch_size     = 50L,
    lr             = 1e-3,
    weight_decay   = 1e-2,
    seed           = if (!is.null(seed)) seed else 123L
  )
  kn_vae <- utils::modifyList(default_vae, vae_params)

  kn_2order <- list(model = "2order")

  # ── 3. Helper: run one method safely, return metrics row ───────────────────
  run_method <- function(method_name, expr) {
    result <- tryCatch(
      expr,
      error = function(e) structure(list(error = conditionMessage(e)),
                                    class = "bk_error")
    )

    if (inherits(result, "bk_error")) {
      return(data.frame(
        method        = method_name,
        fdr_empirical = NA_real_,
        power         = NA_real_,
        n_selected    = NA_integer_,
        n_tp          = NA_integer_,
        n_signal      = n_signal,
        n             = n,
        p             = p,
        rho           = rho,
        fdr_target    = fdr_target,
        seed          = if (is.null(seed)) NA_integer_ else seed,
        error         = result$error,
        stringsAsFactors = FALSE
      ))
    }

    metrics <- .bk_fdr_power(result, true_idx)
    data.frame(
      method        = method_name,
      fdr_empirical = metrics$fdr,
      power         = metrics$power,
      n_selected    = metrics$n_sel,
      n_tp          = metrics$n_tp,
      n_signal      = n_signal,
      n             = n,
      p             = p,
      rho           = rho,
      fdr_target    = fdr_target,
      seed          = if (is.null(seed)) NA_integer_ else seed,
      error         = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  # ── 4. Run each requested method ────────────────────────────────────────────
  rows <- list()

  # --- FLORAL-KN-VAE ---
  if ("FLORAL-KN-VAE" %in% methods) {
    rows[["FLORAL-KN-VAE"]] <- run_method("FLORAL-KN-VAE", {
      fit <- FLORAL(xcount, y,
                    family    = "gaussian",
                    ncov      = 0,
                    ncv       = "knockoff",
                    kn_method = kn_vae,
                    fdr       = fdr_target,
                    progress  = progress,
                    plot      = FALSE,
                    step2     = FALSE)
      .bk_idx_from_names(fit$selected.features)
    })
  }

  # --- FLORAL-KN-2order ---
  if ("FLORAL-KN-2order" %in% methods) {
    rows[["FLORAL-KN-2order"]] <- run_method("FLORAL-KN-2order", {
      fit <- FLORAL(xcount, y,
                    family    = "gaussian",
                    ncov      = 0,
                    ncv       = "knockoff",
                    kn_method = kn_2order,
                    fdr       = fdr_target,
                    progress  = progress,
                    plot      = FALSE,
                    step2     = FALSE)
      .bk_idx_from_names(fit$selected.features)
    })
  }

  # --- FLORAL-lmin and FLORAL-l1se (shared fit) ---
  if (any(c("FLORAL-lmin", "FLORAL-l1se") %in% methods)) {
    cv_fit_result <- tryCatch(
      FLORAL(xcount, y,
             family   = "gaussian",
             ncov     = 0,
             ncv      = ncv,
             progress = progress,
             plot     = FALSE,
             step2    = FALSE),
      error = function(e) structure(list(error = conditionMessage(e)),
                                    class = "bk_error")
    )

    if ("FLORAL-lmin" %in% methods) {
      rows[["FLORAL-lmin"]] <- if (inherits(cv_fit_result, "bk_error")) {
        run_method("FLORAL-lmin", stop(cv_fit_result$error))
      } else {
        run_method("FLORAL-lmin", {
          sel_nms <- names(which(cv_fit_result$best.beta$min != 0))
          .bk_idx_from_names(sel_nms)
        })
      }
    }

    if ("FLORAL-l1se" %in% methods) {
      rows[["FLORAL-l1se"]] <- if (inherits(cv_fit_result, "bk_error")) {
        run_method("FLORAL-l1se", stop(cv_fit_result$error))
      } else {
        run_method("FLORAL-l1se", {
          sel_nms <- names(which(cv_fit_result$best.beta$`1se` != 0))
          .bk_idx_from_names(sel_nms)
        })
      }
    }
  }

  # --- KN-fixed (requires n > p) ---
  if ("KN-fixed" %in% methods) {
    rows[["KN-fixed"]] <- run_method("KN-fixed", {
      if (n <= p) {
        stop(sprintf(
          "Fixed-X knockoffs require n > p (n=%d, p=%d). Skipping.", n, p
        ))
      }
      fit_kn <- knockoff::knockoff.filter(
        X          = xlog1p,
        y          = y,
        knockoffs  = knockoff::create.fixed,
        statistic  = knockoff::stat.glmnet_coefdiff,
        fdr        = fdr_target,
        offset     = 1
      )
      as.integer(fit_kn$selected)
    })
  }

  # --- KN-2order (model-X, standard lasso) ---
  if ("KN-2order" %in% methods) {
    rows[["KN-2order"]] <- run_method("KN-2order", {
      fit_kn <- knockoff::knockoff.filter(
        X          = xlog1p,
        y          = y,
        knockoffs  = knockoff::create.second_order,
        statistic  = knockoff::stat.glmnet_coefdiff,
        fdr        = fdr_target,
        offset     = 1
      )
      as.integer(fit_kn$selected)
    })
  }

  # ── 5. Combine and return ────────────────────────────────────────────────────
  do.call(rbind, rows[methods[methods %in% names(rows)]])
}


# ── run_benchmark_sweep() ─────────────────────────────────────────────────────

#' Sweep over hyperparameter configurations and simulation scenarios
#'
#' Calls \code{\link{benchmark_knockoff}} across a grid of VAE hyperparameter
#' configurations and/or simulation settings, repeating each combination for
#' \code{n_rep} independent replicates.
#'
#' @param param_grid A \code{data.frame} where each row defines one
#'   configuration.  Columns whose names match VAE hyperparameter names
#'   (e.g. \code{gamma_full}, \code{lambda_moments}) are forwarded to
#'   \code{vae_params}.  Columns matching simulation arguments (\code{n},
#'   \code{p}, \code{rho}, \code{weak}, \code{strong}, \code{weaksize},
#'   \code{strongsize}, \code{pct.sparsity}) override the corresponding
#'   defaults.  Use \code{\link{default_vae_grid}} to get a sensible starting
#'   grid.
#' @param n_rep Integer.  Number of independent simulation replicates per
#'   configuration row (default \code{20}).
#' @param n Integer. Default sample size (overridden by \code{param_grid$n}).
#' @param p Integer. Default number of taxa (overridden by \code{param_grid$p}).
#' @param rho Numeric. Default correlation (overridden by \code{param_grid$rho}).
#' @param weak Integer. Default number of weak signals.
#' @param strong Integer. Default number of strong signals.
#' @param weaksize Numeric. Default weak effect size.
#' @param strongsize Numeric. Default strong effect size.
#' @param pct.sparsity Numeric. Default sparsity.
#' @param fdr_target Numeric. Nominal FDR for knockoff methods.
#' @param ncv Integer. CV folds for FLORAL-lmin / FLORAL-l1se.
#' @param methods Character vector of methods to include (see
#'   \code{\link{benchmark_knockoff}}).
#' @param seed_start Integer. First random seed; subsequent replicates use
#'   \code{seed_start + rep - 1}.
#' @param parallel Logical. If \code{TRUE}, uses \code{foreach} with the
#'   currently registered parallel backend (e.g. via \code{doParallel}).
#'   Default \code{FALSE}.
#' @param progress Logical. Print a progress message per configuration row.
#'
#' @return A \code{data.frame} with all results (one row per method × replicate
#'   × configuration), plus a \code{config_id} column identifying the row of
#'   \code{param_grid}.
#'
#' @examples
#' \dontrun{
#' library(doParallel)
#' registerDoParallel(4)
#'
#' grid <- default_vae_grid()
#' results <- run_benchmark_sweep(grid, n_rep = 10, n = 150, p = 50)
#' plot_benchmark(results)
#' }
#'
#' @seealso \code{\link{benchmark_knockoff}}, \code{\link{default_vae_grid}},
#'   \code{\link{plot_benchmark}}
#' @export
run_benchmark_sweep <- function(param_grid,
                                n_rep       = 20,
                                n           = 150,
                                p           = 50,
                                rho         = 0,
                                weak        = 2,
                                strong      = 4,
                                weaksize    = 0.125,
                                strongsize  = 0.25,
                                pct.sparsity = 0.5,
                                fdr_target  = 0.2,
                                ncv         = 5,
                                methods     = c("FLORAL-KN-VAE",
                                                "FLORAL-KN-2order",
                                                "FLORAL-lmin",
                                                "FLORAL-l1se",
                                                "KN-2order"),
                                seed_start  = 1L,
                                parallel    = FALSE,
                                progress    = TRUE) {

  # Column names that are simulation scenario arguments
  simu_args <- c("n", "p", "rho", "weak", "strong",
                 "weaksize", "strongsize", "pct.sparsity")

  # Column names that are VAE hyperparameters
  vae_arg_names <- c("latent_dim", "lambda_kl", "lambda_abun", "lambda_pres",
                     "gamma_full", "gamma_swap", "lambda_moments", "delta_corr",
                     "epochs", "batch_size", "lr", "weight_decay")

  n_cfg <- nrow(param_grid)
  if (n_cfg == 0L) stop("`param_grid` has no rows.")

  run_one_cfg <- function(cfg_idx) {
    cfg <- param_grid[cfg_idx, , drop = FALSE]

    # Build vae_params from this row
    vp_cols <- intersect(names(cfg), vae_arg_names)
    vae_params <- as.list(cfg[, vp_cols, drop = FALSE])

    # Override simulation defaults from this row
    get_val <- function(col, default) {
      if (col %in% names(cfg)) cfg[[col]] else default
    }
    n_cfg_val   <- get_val("n",            n)
    p_cfg_val   <- get_val("p",            p)
    rho_cfg_val <- get_val("rho",          rho)
    wk_val      <- get_val("weak",         weak)
    st_val      <- get_val("strong",       strong)
    wks_val     <- get_val("weaksize",     weaksize)
    sts_val     <- get_val("strongsize",   strongsize)
    spar_val    <- get_val("pct.sparsity", pct.sparsity)

    if (progress)
      message(sprintf("[config %d/%d] n=%d, p=%d, rho=%.2f | VAE: %s",
                      cfg_idx, n_cfg, n_cfg_val, p_cfg_val, rho_cfg_val,
                      paste(names(vae_params),
                            round(unlist(vae_params), 3),
                            sep = "=", collapse = ", ")))

    reps <- lapply(seq_len(n_rep), function(rep_i) {
      seed_i <- seed_start + (rep_i - 1L)  # shared across configs: paired design
      res <- tryCatch(
        benchmark_knockoff(
          n            = n_cfg_val,
          p            = p_cfg_val,
          rho          = rho_cfg_val,
          weak         = wk_val,
          strong       = st_val,
          weaksize     = wks_val,
          strongsize   = sts_val,
          pct.sparsity = spar_val,
          fdr_target   = fdr_target,
          ncv          = ncv,
          vae_params   = vae_params,
          methods      = methods,
          seed         = seed_i,
          progress     = FALSE
        ),
        error = function(e) {
          data.frame(method = methods,
                     fdr_empirical = NA_real_,
                     power = NA_real_,
                     n_selected = NA_integer_,
                     n_tp = NA_integer_,
                     n_signal = NA_integer_,
                     n = n_cfg_val, p = p_cfg_val, rho = rho_cfg_val,
                     fdr_target = fdr_target,
                     seed = seed_i,
                     error = conditionMessage(e),
                     stringsAsFactors = FALSE)
        }
      )
      res$rep       <- rep_i
      res$config_id <- cfg_idx
      # Attach config columns for traceability
      for (col in names(cfg)) res[[col]] <- cfg[[col]]
      res
    })

    do.call(rbind, reps)
  }

  if (parallel) {
    i <- NULL  # satisfy R CMD CHECK for foreach variable
    all_results <- foreach::foreach(
      i          = seq_len(n_cfg),
      .combine   = rbind,
      .packages  = "FLORAL",
      .export    = c("run_one_cfg")
    ) %dopar% run_one_cfg(i)
  } else {
    all_results <- do.call(rbind, lapply(seq_len(n_cfg), run_one_cfg))
  }

  rownames(all_results) <- NULL
  all_results
}


# ── default_vae_grid() ────────────────────────────────────────────────────────

#' Default hyperparameter grid for VAE loss-weight benchmarking
#'
#' Constructs a focused grid for the seven VAE loss weights.  The strategy is
#' a two-stage sweep:
#'
#' \enumerate{
#'   \item Fix the deep-knockoff MMD weights (\code{gamma_full},
#'         \code{gamma_swap}) at their defaults and sweep the
#'         VAE-specific weights (\code{lambda_kl}, \code{lambda_abun},
#'         \code{lambda_pres}, \code{lambda_moments}, \code{delta_corr}).
#'   \item Fix the best VAE-specific weights from stage 1 and sweep the MMD
#'         weights.
#' }
#'
#' By default the function returns the stage-1 grid (the more impactful
#' sweep), which has 32 rows.  Set \code{include_mmd_sweep = TRUE} to append
#' the stage-2 rows.
#'
#' @param include_mmd_sweep Logical.  If \code{TRUE}, also includes a sweep
#'   of \code{gamma_full} and \code{gamma_swap} (stage 2).  Default
#'   \code{FALSE}.
#'
#' @return A \code{data.frame} where each row is one hyperparameter
#'   configuration, ready to pass to \code{\link{run_benchmark_sweep}}.
#'
#' @examples
#' grid <- default_vae_grid()
#' nrow(grid)  # 32 rows (stage 1 only)
#'
#' grid_full <- default_vae_grid(include_mmd_sweep = TRUE)
#' nrow(grid_full)
#'
#' @seealso \code{\link{run_benchmark_sweep}}
#' @export
default_vae_grid <- function(include_mmd_sweep = FALSE) {

  # ── Stage 1: sweep VAE-specific weights, fix MMD at defaults ────────────
  # lambda_kl: KL divergence weight (0 = pure AE, >0 = regularised latent)
  # lambda_abun: abundance reconstruction weight
  # lambda_pres: presence (Bernoulli) reconstruction weight
  # lambda_moments: second-order moment-matching weight
  # delta_corr: decorrelation penalty weight

  stage1 <- expand.grid(
    gamma_full      = 1,
    gamma_swap      = 1,
    lambda_kl       = c(0, 0.1),
    lambda_abun     = c(0, 0.5),
    lambda_pres     = c(0, 0.5),
    lambda_moments  = c(0, 0.5),
    delta_corr      = c(0, 1),
    stringsAsFactors = FALSE
  )

  stage1$sweep_stage <- "vae_weights"
  stage1$config_label <- apply(stage1[, 1:7], 1, function(r)
    paste(names(r), round(r, 3), sep = "=", collapse = "|"))

  if (!include_mmd_sweep) return(stage1)

  # ── Stage 2: sweep MMD weights, fix VAE weights at plausible mid-point ──
  stage2 <- expand.grid(
    gamma_full      = c(0.5, 1, 2, 4),
    gamma_swap      = c(0.5, 1, 2, 4),
    lambda_kl       = 0.1,
    lambda_abun     = 0.5,
    lambda_pres     = 0.5,
    lambda_moments  = 0.5,
    delta_corr      = 1,
    stringsAsFactors = FALSE
  )

  stage2$sweep_stage <- "mmd_weights"
  stage2$config_label <- apply(stage2[, 1:7], 1, function(r)
    paste(names(r), round(r, 3), sep = "=", collapse = "|"))

  rbind(stage1, stage2)
}


# ── plot_benchmark() ──────────────────────────────────────────────────────────

#' Visualise FDR and power from a benchmarking sweep
#'
#' Produces a two-panel ggplot2 figure (FDR and power) comparing methods
#' across simulation replicates.  If \code{results} contains multiple
#' configurations (e.g. from a hyperparameter sweep), the plot facets or
#' colours by \code{facet_by} / \code{colour_by}.
#'
#' @param results A \code{data.frame} returned by \code{\link{run_benchmark_sweep}}
#'   or \code{\link{benchmark_knockoff}}.
#' @param facet_by Character.  Name of a column in \code{results} to use for
#'   faceting (e.g. \code{"sweep_stage"} or \code{"rho"}).  \code{NULL} for no
#'   faceting (default).
#' @param colour_by Character.  Column to map to colour.  Defaults to
#'   \code{"method"}.
#' @param fdr_ref Numeric.  Reference FDR level drawn as a horizontal dashed
#'   line (defaults to the first unique value of \code{results$fdr_target}, or
#'   \code{0.2} if absent).
#' @param ncol Integer.  Number of columns for facets (default \code{2}).
#' @param point_size Numeric.  Size of summary points (default \code{2}).
#' @param errorbar_width Numeric.  Width of error bars (default \code{0.3}).
#'
#' @return A \code{ggplot} object (two panels arranged with
#'   \code{patchwork::wrap_plots}).
#'
#' @examples
#' \dontrun{
#' results <- run_benchmark_sweep(default_vae_grid(), n_rep = 5, n = 100, p = 40)
#' plot_benchmark(results)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_hline
#'   facet_wrap scale_colour_brewer theme_bw theme labs element_text
#'   position_dodge
#' @importFrom stats sd
#' @export
plot_benchmark <- function(results,
                           facet_by      = NULL,
                           colour_by     = "method",
                           fdr_ref       = NULL,
                           ncol          = 2L,
                           point_size    = 2,
                           errorbar_width = 0.3) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plot_benchmark().")

  # Determine FDR reference line
  if (is.null(fdr_ref)) {
    fdr_ref <- if ("fdr_target" %in% names(results))
      unique(results$fdr_target)[1L] else 0.2
  }

  # Summarise (mean ± SE) per grouping variables
  group_vars <- unique(c(colour_by, facet_by, "method"))
  group_vars <- intersect(group_vars, names(results))

  summarise_metric <- function(metric) {
    agg <- do.call(
      rbind,
      lapply(split(results, results[, group_vars, drop = FALSE]), function(grp) {
        vals <- grp[[metric]]
        vals <- vals[!is.na(vals)]
        data.frame(
          mean_val = mean(vals),
          se_val   = if (length(vals) > 1L) sd(vals) / sqrt(length(vals)) else 0,
          n_reps   = length(vals),
          grp[1L, group_vars, drop = FALSE],
          stringsAsFactors = FALSE
        )
      })
    )
    rownames(agg) <- NULL
    agg
  }

  fdr_sum   <- summarise_metric("fdr_empirical")
  power_sum <- summarise_metric("power")

  # Shared aesthetics
  aes_base <- ggplot2::aes(
    x      = .data[[colour_by]],
    y      = .data$mean_val,
    colour = .data[[colour_by]],
    ymin   = pmax(.data$mean_val - .data$se_val, 0),
    ymax   = pmin(.data$mean_val + .data$se_val, 1)
  )

  make_panel <- function(df, y_label, ref_line = NULL) {
    p <- ggplot2::ggplot(df, aes_base) +
      ggplot2::geom_point(size = point_size,
                          position = ggplot2::position_dodge(width = 0.4)) +
      ggplot2::geom_errorbar(width = errorbar_width,
                             position = ggplot2::position_dodge(width = 0.4)) +
      ggplot2::scale_colour_brewer(palette = "Set1") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_text(angle = 30, hjust = 1),
        legend.position = "bottom"
      ) +
      ggplot2::labs(x = NULL, y = y_label, colour = colour_by)

    if (!is.null(ref_line))
      p <- p + ggplot2::geom_hline(yintercept = ref_line,
                                   linetype = "dashed", colour = "grey40")

    if (!is.null(facet_by) && facet_by %in% names(df))
      p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)),
                                   ncol = ncol)
    p
  }

  p_fdr   <- make_panel(fdr_sum,   "Empirical FDR",   ref_line = fdr_ref)
  p_power <- make_panel(power_sum, "Power (sensitivity)")

  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(p_fdr, p_power, ncol = 2L)
  } else {
    message("Install patchwork for a combined two-panel layout.",
            "\nReturning FDR panel only; power panel is the second element.")
    list(fdr = p_fdr, power = p_power)
  }
}
