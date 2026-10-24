# Base-R regression tests for the preserved Coursera matrix-cache assignment.
# Run with: Rscript --vanilla tests/run_tests.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) {
        stop("Run this test file using Rscript", call. = FALSE)
}
script <- sub("^--file=", "", file_arg)
root <- normalizePath(file.path(dirname(normalizePath(script)), ".."))
source(file.path(root, "cachematrix.R"), local = TRUE)

assert_equal_matrix <- function(actual, expected) {
        stopifnot(
                is.matrix(actual),
                identical(dim(actual), dim(expected)),
                isTRUE(all.equal(actual, expected, tolerance = 1e-12))
        )
}

# Fresh matrix inverse and the published accessor contract.
first <- matrix(c(2, 1, 1, 3), nrow = 2)
cache <- makeCacheMatrix(first)
stopifnot(
        is.list(cache),
        identical(sort(names(cache)), sort(c("set", "get", "setInverse", "getInverse"))),
        identical(cache$get(), first),
        is.null(cache$getInverse())
)
computed <- cacheSolve(cache, tol = 1e-12)
assert_equal_matrix(computed, solve(first))
assert_equal_matrix(cache$getInverse(), solve(first))
assert_equal_matrix(cacheSolve(cache), computed)

# A stored inverse must be reused, not recomputed. A singular matrix is useful:
# if solve() ran again instead of using the cache, it would throw an error.
singular <- matrix(c(1, 2, 2, 4), nrow = 2)
singular_cache <- makeCacheMatrix(singular)
sentinel <- diag(2)
singular_cache$setInverse(sentinel)
assert_equal_matrix(suppressMessages(cacheSolve(singular_cache)), sentinel)

# Changing the matrix invalidates the cached inverse and recomputes it.
second <- matrix(c(4, 1, 2, 3), nrow = 2)
cache$set(second)
stopifnot(identical(cache$get(), second), is.null(cache$getInverse()))
second_inverse <- cacheSolve(cache)
assert_equal_matrix(second_inverse, solve(second))
stopifnot(!isTRUE(all.equal(computed, second_inverse)))
assert_equal_matrix(suppressMessages(cacheSolve(cache)), second_inverse)

# A failed inverse must not populate the cache.
failed <- makeCacheMatrix(singular)
error <- tryCatch(cacheSolve(failed), error = function(e) e)
stopifnot(inherits(error, "error"), is.null(failed$getInverse()))

cat("PASS: matrix inverse caching, reuse, invalidation, accessors, and errors\n")
