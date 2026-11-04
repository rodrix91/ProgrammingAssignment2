# Base-R regression checks for the original caching assignment.
# Run from the repository root: Rscript tests/cachematrix-regressions.R
source("cachematrix.R")

a <- matrix(c(2, 1, 1, 3), nrow = 2)
cached <- makeCacheMatrix(a)
expected <- solve(a)
stopifnot(isTRUE(all.equal(cacheSolve(cached), expected)))
stopifnot(identical(cacheSolve(cached), expected))
stopifnot(identical(cached$getInverse(), expected))

# Replacing the source matrix must invalidate the previous inverse.
b <- matrix(c(4, 1, 2, 3), nrow = 2)
cached$set(b)
stopifnot(is.null(cached$getInverse()))
stopifnot(isTRUE(all.equal(cacheSolve(cached), solve(b))))

# A solve(A, b) result is not an inverse and must not leak through
# the original no-argument getInverse() accessor.
rhs <- c(1, 2)
solution <- cacheSolve(cached, b = rhs)
stopifnot(isTRUE(all.equal(solution, solve(b, rhs))))
stopifnot(is.null(cached$getInverse()))
stopifnot(identical(cacheSolve(cached, b = rhs), solution))

# Returning to the inverse must recompute it, not reuse the RHS solution.
stopifnot(isTRUE(all.equal(cacheSolve(cached), solve(b))))
stopifnot(isTRUE(all.equal(cached$getInverse(), solve(b))))

# A failed solve must not overwrite a previously valid cached inverse.
singular <- matrix(c(1, 2, 2, 4), nrow = 2)
cached$set(singular)
stopifnot(inherits(try(cacheSolve(cached), silent = TRUE), "try-error"))
stopifnot(is.null(cached$getInverse()))

cat("cachematrix regression checks passed\n")
