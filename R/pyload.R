# ── Internal helper ───────────────────────────────────────────────────────────

#' Source knockoff.py into the current Python session
#' @noRd
.floral_source_knockoff <- function(pkgname = "FLORAL") {
  py_file <- system.file("python", "knockoff.py", package = pkgname)
  if (py_file == "") {
    warning("Could not find knockoff.py in inst/python. ",
            "Python functions will not be available.", call. = FALSE)
    return(invisible(FALSE))
  }
  tryCatch({
    reticulate::source_python(py_file)
    invisible(TRUE)
  }, error = function(e) {
    warning(paste0(
      "Failed to load knockoff.py: ", e$message,
      "\nPython functions may not be available.",
      "\nEnsure torch and numpy are installed in your Python environment:",
      "\n  pip install torch numpy"
    ), call. = FALSE)
    invisible(FALSE)
  })
}


# ── .onLoad ───────────────────────────────────────────────────────────────────

.onLoad <- function(libname, pkgname) {

  # Reticulate is required
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    warning(
      "The reticulate package is required for Python-based knockoff generation ",
      "but is not installed.\n",
      "Install it with: install.packages('reticulate')\n",
      "Python functions will not be available.",
      call. = FALSE
    )
    return()
  }

  # ── Step 1: honour explicit user preference (set via use_floral_python()) ──
  #
  # We read from package options first.  reticulate::use_*() must be called
  # before Python is initialised; if Python is already running we skip the
  # configuration step (use_floral_python() will have already done it) but
  # still try to source knockoff.py at the end.

  py_initialized <- reticulate::py_available(initialize = FALSE)

  if (!py_initialized) {
    pref_python     <- getOption("FLORAL.python")
    pref_virtualenv <- getOption("FLORAL.virtualenv")
    pref_condaenv   <- getOption("FLORAL.condaenv")

    if (!is.null(pref_python)) {
      tryCatch(
        reticulate::use_python(pref_python, required = FALSE),
        error = function(e) NULL
      )
    } else if (!is.null(pref_virtualenv)) {
      tryCatch(
        reticulate::use_virtualenv(pref_virtualenv, required = FALSE),
        error = function(e) NULL
      )
    } else if (!is.null(pref_condaenv)) {
      tryCatch(
        reticulate::use_condaenv(pref_condaenv, required = FALSE),
        error = function(e) NULL
      )
    } else {
      # ── Step 2: auto-discovery fallback ─────────────────────────────────
      tryCatch({
        python_paths <- c(
          Sys.which("python3"),
          Sys.which("python"),
          "/usr/bin/python3",
          "/usr/local/bin/python3",
          "/opt/homebrew/bin/python3"
        )
        python_paths <- unique(python_paths[nzchar(python_paths)])

        for (py_path in python_paths) {
          if (file.exists(py_path)) {
            tryCatch({
              reticulate::use_python(py_path, required = FALSE)
              if (reticulate::py_available()) break
            }, error = function(e) NULL)
          }
        }

        if (!reticulate::py_available()) {
          cfg <- reticulate::py_discover_config()
          if (!is.null(cfg$python))
            reticulate::use_python(cfg$python, required = FALSE)
        }
      }, error = function(e) NULL)
    }
  }

  # ── Step 3: verify Python is usable ──────────────────────────────────────
  if (!reticulate::py_available(initialize = TRUE)) {
    warning(
      "Python is not available.  FLORAL knockoff functions will not work.\n\n",
      "Quick fix — call one of the following BEFORE library(FLORAL):\n\n",
      "  # Option A – bare Python interpreter\n",
      "  use_floral_python(python = '/usr/bin/python3')\n\n",
      "  # Option B – virtual environment (recommended)\n",
      "  use_floral_python(virtualenv = 'floral-env')\n\n",
      "  # Option C – conda environment\n",
      "  use_floral_python(condaenv = 'floral')\n\n",
      "Then restart R and reload the package.\n",
      "See ?use_floral_python for full details.",
      call. = FALSE
    )
    return()
  }

  # ── Step 4: source the knockoff script ───────────────────────────────────
  .floral_source_knockoff(pkgname)
}


# ── Public API ────────────────────────────────────────────────────────────────

#' Configure the Python environment used by FLORAL
#'
#' FLORAL's knockoff generation relies on PyTorch via
#' \href{https://rstudio.github.io/reticulate/}{reticulate}.  Python must be
#' configured \emph{before} reticulate initialises it — which happens the first
#' time any Python function is called (including during \code{library(FLORAL)}).
#'
#' The recommended workflow is therefore:
#'
#' \enumerate{
#'   \item Call \code{use_floral_python()} \strong{once} in your
#'         \code{.Rprofile} or at the very top of your script, before
#'         \code{library(FLORAL)}.
#'   \item Load the package normally: \code{library(FLORAL)}.
#' }
#'
#' Calling this function \emph{after} Python has already been initialised in
#' the current session will print an informative message but cannot change the
#' active interpreter; restart R to apply a new setting.
#'
#' @param python Character string. Full path to a Python binary, e.g.
#'   \code{"/usr/bin/python3"} or \code{"/opt/homebrew/bin/python3"}.  Mutually
#'   exclusive with \code{virtualenv} and \code{condaenv}.
#' @param virtualenv Character string. Name or path of a Python virtual
#'   environment created with \code{python -m venv} or
#'   \code{reticulate::virtualenv_create()}.  Mutually exclusive with
#'   \code{python} and \code{condaenv}.
#' @param condaenv Character string. Name of a conda environment.  Mutually
#'   exclusive with \code{python} and \code{virtualenv}.
#' @param required Logical.  If \code{TRUE} (default), an error is raised when
#'   the requested environment cannot be found.  Set to \code{FALSE} to allow
#'   silent fallback.
#'
#' @return Invisibly returns a list with fields \code{python_available},
#'   \code{torch_available}, and \code{knockoff_loaded} so callers can check
#'   the outcome programmatically.
#'
#' @section Installing dependencies:
#' The knockoff script requires \code{torch} and \code{numpy}.  After creating
#' your environment, install them with:
#'
#' \preformatted{
#' # virtual environment
#' reticulate::virtualenv_create("floral-env")
#' reticulate::virtualenv_install("floral-env", c("torch", "numpy"))
#'
#' # conda
#' reticulate::conda_create("floral")
#' reticulate::conda_install("floral", c("pytorch", "numpy"), channel = "pytorch")
#' }
#'
#' @examples
#' \dontrun{
#' # --- Put this in .Rprofile or at the top of your script ---
#'
#' # Option A: point to a specific Python binary
#' use_floral_python(python = "/usr/bin/python3")
#'
#' # Option B: use a virtual environment (recommended for isolation)
#' use_floral_python(virtualenv = "floral-env")
#'
#' # Option C: use a conda environment
#' use_floral_python(condaenv = "floral")
#'
#' # --- Then load the package as usual ---
#' library(FLORAL)
#' }
#'
#' @seealso \code{\link{check_python_setup}} for diagnosing the current setup.
#' @export
use_floral_python <- function(python     = NULL,
                              virtualenv = NULL,
                              condaenv   = NULL,
                              required   = TRUE) {

  # ── Validate arguments ─────────────────────────────────────────────────────
  n_specified <- sum(!is.null(python), !is.null(virtualenv), !is.null(condaenv))
  if (n_specified == 0L) {
    stop(
      "Please specify exactly one of `python`, `virtualenv`, or `condaenv`.\n\n",
      "Examples:\n",
      "  use_floral_python(python     = '/usr/bin/python3')\n",
      "  use_floral_python(virtualenv = 'floral-env')\n",
      "  use_floral_python(condaenv   = 'floral')",
      call. = FALSE
    )
  }
  if (n_specified > 1L) {
    stop(
      "Specify only one of `python`, `virtualenv`, or `condaenv`, not multiple.",
      call. = FALSE
    )
  }

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop(
      "The reticulate package is not installed.\n",
      "Install it with: install.packages('reticulate')",
      call. = FALSE
    )
  }

  # ── Store preference in package options ───────────────────────────────────
  #    .onLoad reads these when the package (re)loads, so the preference
  #    persists across detach/library cycles within the same session.
  options(FLORAL.python     = python)
  options(FLORAL.virtualenv = virtualenv)
  options(FLORAL.condaenv   = condaenv)

  # ── Check if Python has already been initialised ──────────────────────────
  py_initialized <- reticulate::py_available(initialize = FALSE)

  result <- list(
    python_available = FALSE,
    torch_available  = FALSE,
    knockoff_loaded  = FALSE
  )

  if (py_initialized) {
    message(
      "Python is already initialised in this R session.  The active ",
      "interpreter cannot be changed without restarting R.\n\n",
      "Your preference has been saved and will take effect the next time ",
      "you start R and call library(FLORAL).\n\n",
      "To apply it now:\n",
      "  1. Restart R  (Session > Restart R in RStudio)\n",
      "  2. Call use_floral_python() again before library(FLORAL)\n",
      "  3. Then library(FLORAL)"
    )
    result$python_available <- reticulate::py_available()
    result$torch_available  <- reticulate::py_module_available("torch")
    # Try to re-source knockoff.py in the already-running interpreter anyway
    result$knockoff_loaded  <- .floral_source_knockoff()
    return(invisible(result))
  }

  # ── Configure Python ───────────────────────────────────────────────────────
  tryCatch({
    if (!is.null(python)) {
      if (!file.exists(python))
        warning("Python binary not found at: ", python,
                "\nDouble-check the path.", call. = FALSE)
      reticulate::use_python(python, required = required)
      msg_what <- paste0("Python binary: ", python)

    } else if (!is.null(virtualenv)) {
      reticulate::use_virtualenv(virtualenv, required = required)
      msg_what <- paste0("virtual environment: ", virtualenv)

    } else {
      reticulate::use_condaenv(condaenv, required = required)
      msg_what <- paste0("conda environment: ", condaenv)
    }
  }, error = function(e) {
    stop(
      "Failed to configure Python environment.\n",
      "Details: ", e$message, "\n\n",
      "Run check_python_setup() for diagnostics.",
      call. = FALSE
    )
  })

  # ── Verify Python is reachable ────────────────────────────────────────────
  if (!reticulate::py_available(initialize = TRUE)) {
    stop(
      "Python environment was configured but Python could not be initialised.\n",
      "Run check_python_setup() for diagnostics.",
      call. = FALSE
    )
  }
  result$python_available <- TRUE

  # ── Check torch ───────────────────────────────────────────────────────────
  result$torch_available <- reticulate::py_module_available("torch")
  if (!result$torch_available) {
    message(
      "Python configured successfully using ", msg_what, ".\n",
      "However, 'torch' was not found in this environment.\n\n",
      "Install it with one of:\n",
      if (!is.null(virtualenv))
        paste0("  reticulate::virtualenv_install('", virtualenv,
               "', c('torch', 'numpy'))\n")
      else if (!is.null(condaenv))
        paste0("  reticulate::conda_install('", condaenv,
               "', c('pytorch', 'numpy'), channel = 'pytorch')\n")
      else
        "  pip install torch numpy   (in your terminal)\n",
      "\nKnockoff generation will not work until torch is installed."
    )
  } else {
    message("Python configured successfully using ", msg_what, ".")
  }

  # ── Source knockoff.py ────────────────────────────────────────────────────
  result$knockoff_loaded <- .floral_source_knockoff()

  invisible(result)
}


#' Check Python and PyTorch setup
#'
#' Diagnostic function to check if Python, PyTorch, and required dependencies
#' are properly configured.  Call this when \code{train_vae()} or other Python
#' functions fail unexpectedly.
#'
#' @param auto_configure Logical.  If \code{TRUE} (default), attempts to
#'   automatically configure Python if it is installed but not yet available to
#'   reticulate.
#'
#' @return A named list with diagnostic information.  Key fields:
#' \describe{
#'   \item{\code{python_available}}{Logical — is Python reachable?}
#'   \item{\code{python_version}}{Version string of the active interpreter.}
#'   \item{\code{torch_available}}{Logical — is the \code{torch} module importable?}
#'   \item{\code{torch_version}}{Version string of PyTorch.}
#'   \item{\code{numpy_available}}{Logical — is \code{numpy} importable?}
#'   \item{\code{vae_func_loaded}}{Logical — was \code{knockoff.py} sourced correctly?}
#'   \item{\code{fix_suggestion}}{Character — human-readable fix advice (if any).}
#' }
#'
#' @examples
#' \dontrun{
#' check_python_setup()
#' }
#'
#' @seealso \code{\link{use_floral_python}} to configure the environment.
#' @export
check_python_setup <- function(auto_configure = TRUE) {
  result <- list()

  # reticulate present?
  result$reticulate_available <- requireNamespace("reticulate", quietly = TRUE)
  if (!result$reticulate_available) {
    result$error          <- "reticulate package is not installed"
    result$fix_suggestion <- "Install with: install.packages('reticulate')"
    return(result)
  }

  # Python available?
  result$python_available <- reticulate::py_available()

  if (!result$python_available && auto_configure) {
    result$auto_config_attempted <- TRUE

    python_paths <- unique(c(
      Sys.which("python3"), Sys.which("python"),
      "/usr/bin/python3", "/usr/local/bin/python3",
      "/opt/homebrew/bin/python3"
    ))
    python_paths <- python_paths[nzchar(python_paths) & file.exists(python_paths)]
    result$python_paths_found <- python_paths

    for (py_path in python_paths) {
      tryCatch({
        reticulate::use_python(py_path, required = FALSE)
        if (reticulate::py_available()) {
          result$python_available   <- TRUE
          result$python_path_used   <- py_path
          result$auto_config_success <- TRUE
          break
        }
      }, error = function(e) NULL)
    }

    if (!result$python_available) {
      tryCatch({
        cfg <- reticulate::py_discover_config()
        if (!is.null(cfg$python) && file.exists(cfg$python)) {
          reticulate::use_python(cfg$python, required = FALSE)
          if (reticulate::py_available()) {
            result$python_available   <- TRUE
            result$python_path_used   <- cfg$python
            result$auto_config_success <- TRUE
          }
        }
      }, error = function(e) {
        result$discover_config_error <- e$message
      })
    }
  }

  if (!result$python_available) {
    result$error <- "Python is not available"
    result$fix_suggestion <- paste0(
      "Call use_floral_python() BEFORE library(FLORAL), e.g.:\n\n",
      "  use_floral_python(python     = '/usr/bin/python3')\n",
      "  use_floral_python(virtualenv = 'floral-env')\n",
      "  use_floral_python(condaenv   = 'floral')\n\n",
      "Then restart R and reload the package."
    )
    return(result)
  }

  # Python version
  tryCatch(
    result$python_version <- reticulate::py_config()$version,
    error = function(e) result$python_version <<- "unknown"
  )

  # VAE function loaded?
  result$vae_func_loaded <- tryCatch({
    if (exists("VAE_func_DK", envir = .GlobalEnv, inherits = TRUE)) {
      TRUE
    } else {
      main <- reticulate::import_main()
      reticulate::py_has_attr(main, "VAE_func_DK")
    }
  }, error = function(e) {
    result$vae_func_error <<- e$message
    FALSE
  })

  # torch
  result$torch_available <- FALSE
  result$torch_version   <- "not installed"
  tryCatch({
    if (reticulate::py_module_available("torch")) {
      result$torch_available <- TRUE
      torch <- reticulate::import("torch")
      result$torch_version <- torch$`__version__`
    }
  }, error = function(e) result$torch_error <<- e$message)

  # numpy
  result$numpy_available <- tryCatch(
    reticulate::py_module_available("numpy"),
    error = function(e) FALSE
  )

  # Overall status message
  if (result$python_available && result$torch_available && result$vae_func_loaded) {
    result$status <- "OK — FLORAL Python environment is fully configured."
  } else {
    issues <- c(
      if (!result$torch_available) "torch not installed",
      if (!result$numpy_available) "numpy not installed",
      if (!result$vae_func_loaded) "knockoff.py not loaded"
    )
    result$status <- paste0(
      "Issues detected: ", paste(issues, collapse = ", "), ".\n",
      "Run use_floral_python() to reconfigure."
    )
  }

  result
}
