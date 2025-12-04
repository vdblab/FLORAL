#' Train VAE using PyTorch (replacement for R torch version)
#'
#' This function wraps the Python PyTorch implementation of VAE_func_DK,
#' providing the same interface as the original R torch version.
#'
#' @param x Input data matrix (n_samples x n_features)
#' @param latent_dim Dimension of latent space (default: 32)
#' @param lambda_kl Weight for KL divergence loss (default: 0)
#' @param lambda_abun Weight for abundance loss (default: 0)
#' @param lambda_pres Weight for presence loss (default: 0)
#' @param gamma_full Weight for MMD full loss (default: 1)
#' @param gamma_swap Weight for MMD swap loss (default: 1)
#' @param lambda_moments Weight for moment matching loss (default: 0)
#' @param delta_corr Weight for correlation loss (default: 0)
#' @param sigma_list List of sigma values for MMD kernel (default: c(1,2,4,8,16,32,64,128))
#' @param epochs Number of training epochs (default: 100)
#' @param batch_size Batch size for training (default: 50)
#' @param lr Learning rate (default: 1e-3)
#' @param weight_decay Weight decay for optimizer (default: 1e-2)
#' @param seed Random seed (default: 123)
#'
#' @return A list with components:
#'   \item{knockoff_x}{Generated knockoff data matrix}
#'   \item{recon_x}{Reconstructed data matrix}
#'
#' @export
train_vae <- function(x,
                      latent_dim      = 32,
                      lambda_kl       = 0,
                      lambda_abun     = 0,
                      lambda_pres     = 0,
                      gamma_full       = 1,
                      gamma_swap       = 1,
                      lambda_moments  = 0,
                      delta_corr      = 0,
                      sigma_list      = c(1,2,4,8,16,32,64,128),
                      epochs          = 100,
                      batch_size      = 50,
                      lr              = 1e-3,
                      weight_decay    = 1e-2,
                      seed            = 123) {
  
  # Check if Python is available
  if (!reticulate::py_available()) {
    stop("Python is not available. Please install Python and reticulate package.", call. = FALSE)
  }
  
  # Configure PyTorch to avoid threading issues with R
  # This must be done before importing torch
  tryCatch({
    # Set environment variables to control PyTorch threading
    # This prevents conflicts with R's threading
    Sys.setenv(OMP_NUM_THREADS = "1")
    Sys.setenv(MKL_NUM_THREADS = "1")
    Sys.setenv(NUMEXPR_NUM_THREADS = "1")
    
    # Try to configure torch if it's already loaded
    if (reticulate::py_module_available("torch")) {
      torch <- reticulate::import("torch")
      # Set threads to 1 to avoid conflicts with R
      tryCatch({
        torch$set_num_threads(1L)
        torch$set_num_interop_threads(1L)
      }, error = function(e) {
        # Ignore if threads can't be set (might be set already)
      })
    }
  }, error = function(e) {
    warning("Could not configure PyTorch threading. This may cause issues.", call. = FALSE)
  })
  
  # Get the Python function (source_python makes it available)
  # Try direct access first, then fallback to main module
  VAE_func_DK_py <- NULL
  tryCatch({
    if (exists("VAE_func_DK", envir = .GlobalEnv, inherits = TRUE)) {
      VAE_func_DK_py <- get("VAE_func_DK", envir = .GlobalEnv)
    } else {
      # Access from Python main module
      main <- reticulate::import_main()
      if (!reticulate::py_has_attr(main, "VAE_func_DK")) {
        stop("Python function VAE_func_DK is not loaded. Did .onLoad fail? Please ensure knockoff.py is loaded.", call. = FALSE)
      }
      VAE_func_DK_py <- reticulate::py_get_attr(main, "VAE_func_DK")
    }
  }, error = function(e) {
    stop(paste("Failed to access Python function VAE_func_DK:", e$message, 
               "\nTry reloading the package or check if Python dependencies are installed."), 
         call. = FALSE)
  })
  
  # Validate input
  if (missing(x) || is.null(x)) {
    stop("Input data 'x' is required.", call. = FALSE)
  }
  
  # Convert R matrix/array to appropriate format for Python
  # reticulate will handle the conversion automatically, but we ensure it's a matrix
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  
  # Check for valid dimensions
  if (nrow(x) < 2) {
    stop("Input data must have at least 2 rows.", call. = FALSE)
  }
  if (ncol(x) < 1) {
    stop("Input data must have at least 1 column.", call. = FALSE)
  }
  
  # Check for NaN or Inf values
  if (any(!is.finite(x))) {
    warning("Input data contains non-finite values (NaN or Inf). These will be handled by Python but may cause issues.", call. = FALSE)
  }
  
  # Call Python function with error handling
  result <- tryCatch({
    # Call Python function
    # reticulate automatically converts R types to Python types
    VAE_func_DK_py(
      x = x,
      latent_dim = as.integer(latent_dim),
      lambda_kl = as.numeric(lambda_kl),
      lambda_abun = as.numeric(lambda_abun),
      lambda_pres = as.numeric(lambda_pres),
      gamma_full = as.numeric(gamma_full),
      gamma_swap = as.numeric(gamma_swap),
      lambda_moments = as.numeric(lambda_moments),
      delta_corr = as.numeric(delta_corr),
      sigma_list = as.numeric(sigma_list),
      epochs = as.integer(epochs),
      batch_size = as.integer(batch_size),
      lr = as.numeric(lr),
      weight_decay = as.numeric(weight_decay),
      seed = as.integer(seed)
    )
  }, error = function(e) {
    # Provide more informative error messages
    error_msg <- e$message
    if (grepl("torch", error_msg, ignore.case = TRUE)) {
      stop(paste("PyTorch error:", error_msg, 
                 "\n\nPossible solutions:",
                 "1. Ensure PyTorch is installed: pip install torch",
                 "2. Try restarting R session",
                 "3. Check if there are memory issues"), 
           call. = FALSE)
    } else if (grepl("memory", error_msg, ignore.case = TRUE)) {
      stop(paste("Memory error:", error_msg,
                 "\n\nTry reducing batch_size or using smaller data."), 
           call. = FALSE)
    } else {
      stop(paste("Python function error:", error_msg), call. = FALSE)
    }
  })
  
  # Validate result
  if (is.null(result)) {
    stop("Python function returned NULL. Check Python console for errors.", call. = FALSE)
  }
  
  # Convert Python dict to R list
  # reticulate converts Python dicts to R lists automatically
  # Ensure the structure matches the original R function output
  tryCatch({
    return(list(
      knockoff_x = result$knockoff_x,
      recon_x = result$recon_x
    ))
  }, error = function(e) {
    stop(paste("Failed to extract results from Python output:", e$message,
               "\nResult structure:", paste(names(result), collapse = ", ")), 
         call. = FALSE)
  })
}

#' Check Python and PyTorch setup
#'
#' Diagnostic function to check if Python, PyTorch, and required dependencies are properly configured.
#' This function also attempts to auto-configure Python if it's installed but not yet configured.
#'
#' @param auto_configure Logical. If TRUE (default), attempts to automatically configure Python
#'   if it's installed but not yet available to reticulate.
#'
#' @return A list with diagnostic information about Python setup
#' @export
check_python_setup <- function(auto_configure = TRUE) {
  result <- list()
  
  # Check reticulate
  result$reticulate_available <- requireNamespace("reticulate", quietly = TRUE)
  if (!result$reticulate_available) {
    result$error <- "reticulate package is not installed"
    result$fix_suggestion <- "Install with: install.packages('reticulate')"
    return(result)
  }
  
  # Check Python
  result$python_available <- reticulate::py_available()
  
  # Try to auto-configure if Python is not available
  if (!result$python_available && auto_configure) {
    result$auto_config_attempted <- TRUE
    
    # Try common Python paths
    python_paths <- c(
      Sys.which("python3"),
      Sys.which("python"),
      "/usr/bin/python3",
      "/usr/local/bin/python3",
      "/opt/homebrew/bin/python3",
      "C:/Python*/python.exe",
      "C:/Program Files/Python*/python.exe"
    )
    
    # Filter out empty paths and check if files exist
    python_paths <- python_paths[python_paths != ""]
    python_paths <- python_paths[file.exists(python_paths)]
    
    result$python_paths_found <- python_paths
    
    # Try each path
    for (py_path in python_paths) {
      tryCatch({
        reticulate::use_python(py_path, required = FALSE)
        if (reticulate::py_available()) {
          result$python_available <- TRUE
          result$python_path_used <- py_path
          result$auto_config_success <- TRUE
          break
        }
      }, error = function(e) {
        # Try next path
      })
    }
    
    # If still not available, try py_discover_config
    if (!result$python_available) {
      tryCatch({
        config <- reticulate::py_discover_config()
        if (!is.null(config$python) && file.exists(config$python)) {
          reticulate::use_python(config$python, required = FALSE)
          if (reticulate::py_available()) {
            result$python_available <- TRUE
            result$python_path_used <- config$python
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
    result$fix_suggestion <- paste(
      "To fix this:\n",
      "1. Install Python (>= 3.8) from https://www.python.org/downloads/\n",
      "2. Configure reticulate: reticulate::use_python('/path/to/python')\n",
      "3. Or use a virtual environment: reticulate::use_virtualenv('envname')\n",
      "4. Or use conda: reticulate::use_condaenv('envname')"
    )
    return(result)
  }
  
  # Get Python version
  tryCatch({
    result$python_version <- reticulate::py_config()$version
  }, error = function(e) {
    result$python_version <- "unknown"
  })
  
  # Check if VAE_func_DK is loaded
  result$vae_func_loaded <- FALSE
  tryCatch({
    if (exists("VAE_func_DK", envir = .GlobalEnv, inherits = TRUE)) {
      result$vae_func_loaded <- TRUE
    } else {
      main <- reticulate::import_main()
      result$vae_func_loaded <- reticulate::py_has_attr(main, "VAE_func_DK")
    }
  }, error = function(e) {
    result$vae_func_loaded <- FALSE
    result$vae_func_error <- e$message
  })
  
  # Check PyTorch
  result$torch_available <- FALSE
  result$torch_version <- "unknown"
  tryCatch({
    if (reticulate::py_module_available("torch")) {
      result$torch_available <- TRUE
      torch <- reticulate::import("torch")
      result$torch_version <- torch$`__version__`
    }
  }, error = function(e) {
    result$torch_error <- e$message
  })
  
  # Check numpy
  result$numpy_available <- FALSE
  tryCatch({
    result$numpy_available <- reticulate::py_module_available("numpy")
  }, error = function(e) {
    result$numpy_error <- e$message
  })
  
  return(result)
}