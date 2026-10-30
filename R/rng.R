## ---------------------------------------------------------------------------
##  Random-number helpers. Functions with a `seed` argument must not leave
##  the caller's stream altered, and bootstrap replicates must give the same
##  result serially and in parallel, so each replicate gets its own seed.
## ---------------------------------------------------------------------------

.rng_get <- function() {
  get0(".Random.seed", envir = globalenv(), inherits = FALSE)
}

.rng_set <- function(state) {
  if (is.null(state)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  } else {
    assign(".Random.seed", state, envir = globalenv())
  }
  invisible(NULL)
}

# Seed the stream for the rest of the calling function and restore the
# caller's state when it exits. Does nothing when `seed` is NULL.
.local_seed <- function(seed, envir = parent.frame()) {
  if (is.null(seed)) return(invisible(NULL))
  old <- .rng_get()
  do.call(on.exit, list(substitute(.rng_set(old), list(old = old)),
                        add = TRUE), envir = envir)
  set.seed(seed)
  invisible(NULL)
}
