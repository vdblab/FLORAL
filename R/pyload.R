.onLoad <- function(libname, pkgname) {
  # Check if reticulate is available
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    warning("reticulate package is required but not installed. ",
            "Install it with: install.packages('reticulate')",
            "\nPython functions will not be available.", call. = FALSE)
    return()
  }
  
  # Try to configure Python if not already available
  if (!reticulate::py_available(initialize = TRUE)) {
    # Try to find and use Python automatically
    tryCatch({
      # Check if Python is installed but not configured
      python_paths <- c(
        Sys.which("python3"),
        Sys.which("python"),
        "/usr/bin/python3",
        "/usr/local/bin/python3",
        "/opt/homebrew/bin/python3",
        "C:/Python*/python.exe",
        "C:/Program Files/Python*/python.exe"
      )
      
      # Filter out empty paths
      python_paths <- python_paths[python_paths != ""]
      
      # Try to use the first available Python
      if (length(python_paths) > 0) {
        for (py_path in python_paths) {
          if (file.exists(py_path)) {
            tryCatch({
              reticulate::use_python(py_path, required = FALSE)
              if (reticulate::py_available()) {
                break
              }
            }, error = function(e) {
              # Try next path
            })
          }
        }
      }
      
      # If still not available, try py_discover_config
      if (!reticulate::py_available()) {
        tryCatch({
          config <- reticulate::py_discover_config()
          if (!is.null(config$python)) {
            reticulate::use_python(config$python, required = FALSE)
          }
        }, error = function(e) {
          # Ignore errors
        })
      }
    }, error = function(e) {
      # Ignore configuration errors
    })
  }
  
  # Check if Python is available after configuration attempts
  if (!reticulate::py_available()) {
    warning("Python is not available. ",
            "To fix this:\n",
            "1. Install Python (>= 3.8) from https://www.python.org/downloads/\n",
            "2. Configure reticulate: reticulate::use_python('/path/to/python')\n",
            "3. Or use a virtual environment: reticulate::use_virtualenv('envname')\n",
            "4. Or use conda: reticulate::use_condaenv('envname')\n",
            "\nPython functions will not be available until Python is configured.",
            call. = FALSE)
    return()
  }
  
  # Find the path to the Python script inside the installed package
  py_file <- system.file("python", "knockoff.py", package = pkgname)
  
  if (py_file == "") {
    warning("Could not find knockoff.py in inst/python. ",
            "Python functions will not be available.", call. = FALSE)
    return()
  }
  
  # Try to source the Python file with error handling
  tryCatch({
    reticulate::source_python(py_file)
  }, error = function(e) {
    warning(paste("Failed to load knockoff.py:", e$message, 
                  "\nPython functions may not be available. ",
                  "Check Python dependencies (torch, numpy).",
                  "\nInstall with: pip install torch numpy"),
            call. = FALSE)
  })
}