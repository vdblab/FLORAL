concordance_index <- function(x, surv.time, surv.event) {
  # Remove NA values
  cc_ix <- stats::complete.cases(x, surv.time, surv.event)
  if (sum(cc_ix) < 3) {
    return(list(c.index = NA))
  }
  x <- x[cc_ix]
  surv.time <- surv.time[cc_ix]
  surv.event <- surv.event[cc_ix]

  # No events - cannot compute concordance
  if (sum(surv.event) == 0) {
    return(list(c.index = NA))
  }

  result <- tryCatch(
    survival::concordance(
      survival::Surv(surv.time, surv.event) ~ x
    ),
    error = function(e) NULL,
    warning = function(w) NULL
  )

  if (is.null(result)) {
    return(list(c.index = NA))
  }

  c_index <- as.numeric(result$concordance)
  if (is.na(c_index) || !is.finite(c_index)) {
    return(list(c.index = NA))
  }

  list(c.index = c_index)
}
