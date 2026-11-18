# Base-R regression coverage for complex-valued linear algebra.
# No external packages, no network, no modifications to the assignment source.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) stop("Run with Rscript", call. = FALSE)
script <- sub("^--file=", "", file_arg)
root <- normalizePath(file.path(dirname(normalizePath(script)), ".."))
source(file.path(root, "cachematrix.R"), local = TRUE)

same <- function(actual, expected) {
    stopifnot(
        identical(dim(actual), dim(expected)),
        isTRUE(all.equal(actual, expected, tolerance = 1e-12))
    )
}
cache_messages <- function(expression) {
    messages <- character()
    value <- withCallingHandlers(
        expression,
        message = function(m) {
            messages <<- c(messages, conditionMessage(m))
            invokeRestart("muffleMessage")
        }
    )
    list(value = value, messages = messages)
}

# Non-real, invertible matrix: retaining the imaginary component is essential.
original <- matrix(c(2 + 1i, 1 - 2i, 3 + 0i, 4 + 1i), nrow = 2)
cache <- makeCacheMatrix(original)
inverse <- cacheSolve(cache)
stopifnot(is.matrix(inverse), is.complex(inverse))
same(inverse, solve(original))
same(cache$getInverse(), inverse)
hit <- cache_messages(cacheSolve(cache))
same(hit$value, inverse)
stopifnot(length(hit$messages) == 1L,
          grepl("getting cached data", hit$messages[[1]]))

# Complex RHS vectors and matrices are solutions, not matrix inverses.
rhs <- c(1 + 2i, 2 - 1i)
solution <- cacheSolve(cache, b = rhs)
same(solution, solve(original, b = rhs))
stopifnot(is.complex(solution), is.null(cache$getInverse()))
rhs_hit <- cache_messages(cacheSolve(cache, b = rhs))
same(rhs_hit$value, solution)
stopifnot(length(rhs_hit$messages) == 1L)
other_rhs <- c(2 - 1i, 1 + 3i)
other_solution <- cacheSolve(cache, b = other_rhs)
same(other_solution, solve(original, b = other_rhs))
stopifnot(!isTRUE(all.equal(solution, other_solution)))

rhs_matrix <- matrix(c(1 + 0i, 2 - 1i, 3 + 2i, 4 - 3i), nrow = 2)
matrix_solution <- cacheSolve(cache, b = rhs_matrix)
same(matrix_solution, solve(original, b = rhs_matrix))
stopifnot(is.null(cache$getInverse()))

# Switching back to inverse form must not reuse any RHS solution.
same(cacheSolve(cache), solve(original))
same(cache$getInverse(), inverse)

# Copy-on-modify semantics preserve the original closure's source matrix.
original[1, 1] <- 99 + 4i
same(cache$get(), matrix(c(2 + 1i, 1 - 2i, 3 + 0i, 4 + 1i), nrow = 2))
same(cacheSolve(cache), inverse)

# A genuinely singular complex matrix must fail without leaving a cached value.
singular <- matrix(c(1 + 1i, 2 + 2i, 2 + 2i, 4 + 4i), nrow = 2)
cache$set(singular)
stopifnot(is.null(cache$getInverse()))
failed <- tryCatch(cacheSolve(cache), error = function(e) e)
stopifnot(inherits(failed, "error"), is.null(cache$getInverse()))

# A fresh invertible complex matrix can recover the same closure.
recovered <- diag(c(2 + 1i, 3 - 2i))
cache$set(recovered)
same(cacheSolve(cache), solve(recovered))
same(cache$getInverse(), solve(recovered))

cat("PASS: complex matrix inverses, RHS isolation, cache hits and recovery\n")
