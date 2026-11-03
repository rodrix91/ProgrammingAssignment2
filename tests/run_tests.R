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

# Distinct closures must never share cached values or invalidation state.
left <- makeCacheMatrix(diag(c(2, 3)))
right <- makeCacheMatrix(diag(c(5, 7)))
left_inverse <- cacheSolve(left)
right_inverse <- cacheSolve(right)
assert_equal_matrix(left_inverse, solve(left$get()))
assert_equal_matrix(right_inverse, solve(right$get()))
left$set(matrix(c(6, 1, 2, 5), nrow = 2))
stopifnot(is.null(left$getInverse()))
assert_equal_matrix(right$getInverse(), right_inverse)
assert_equal_matrix(suppressMessages(cacheSolve(right)), right_inverse)
assert_equal_matrix(cacheSolve(left), solve(left$get()))
stopifnot(!isTRUE(all.equal(left$getInverse(), right$getInverse())))

# The assignment's closure works with square matrices of varying dimensions,
# including a 1x1 matrix and a non-symmetric 3x3 matrix.
single <- makeCacheMatrix(matrix(8, nrow = 1))
assert_equal_matrix(cacheSolve(single), matrix(1 / 8, nrow = 1))
assert_equal_matrix(suppressMessages(cacheSolve(single)), single$getInverse())

third_matrix <- matrix(c(3, 0, 0, 1, 4, 0, 2, 1, 5), nrow = 3)
third_cache <- makeCacheMatrix(third_matrix)
third_inverse <- cacheSolve(third_cache)
assert_equal_matrix(third_inverse, solve(third_matrix))
stopifnot(
        identical(dim(third_inverse), c(3L, 3L)),
        isTRUE(all.equal(third_matrix %*% third_inverse, diag(3), tolerance = 1e-12))
)

# cacheSolve must forward solve() arguments on cache misses. A deliberately
# strict tolerance rejects an ill-conditioned but invertible diagonal matrix.
# Failed numerical checks must leave the cached value unset.
ill_conditioned <- diag(c(1, 1000))
tolerance_cache <- makeCacheMatrix(ill_conditioned)
strict_error <- tryCatch(cacheSolve(tolerance_cache, tol = 0.5), error = function(e) e)
stopifnot(inherits(strict_error, "error"), is.null(tolerance_cache$getInverse()))
assert_equal_matrix(cacheSolve(tolerance_cache, tol = 0), solve(ill_conditioned, tol = 0))
assert_equal_matrix(suppressMessages(cacheSolve(tolerance_cache)), tolerance_cache$getInverse())

# Cache hits are only valid when solve() is called with the SAME options.
# A lenient tolerance must never bypass a later, stricter matrix check.
option_cache <- makeCacheMatrix(ill_conditioned)
lenient_inverse <- cacheSolve(option_cache, tol = 0)
assert_equal_matrix(lenient_inverse, solve(ill_conditioned, tol = 0))
messages <- character()
cached_lenient <- withCallingHandlers(
        cacheSolve(option_cache, tol = 0),
        message = function(m) {
                messages <<- c(messages, conditionMessage(m))
                invokeRestart("muffleMessage")
        }
)
assert_equal_matrix(cached_lenient, lenient_inverse)
stopifnot(length(messages) == 1L, grepl("getting cached data", messages[[1]]))
strict_after_lenient <- tryCatch(
        cacheSolve(option_cache, tol = 0.5),
        error = function(e) e
)
stopifnot(inherits(strict_after_lenient, "error"))
# A failed stricter solve leaves the successful lenient cached entry intact
# without ever returning it to callers using the stricter option.
assert_equal_matrix(option_cache$getInverse(), lenient_inverse)
assert_equal_matrix(suppressMessages(cacheSolve(option_cache, tol = 0)), lenient_inverse)

# solve() also accepts a right-hand-side vector via b=..., yielding a
# different type of result. Switching b or omitting it must never reuse
# the cached result for a different problem.
rhs_cache <- makeCacheMatrix(first)
rhs_one <- cacheSolve(rhs_cache, b = c(3, 5))
stopifnot(isTRUE(all.equal(rhs_one, solve(first, b = c(3, 5)))))
# A solution vector is not the inverse matrix exposed by getInverse().
stopifnot(is.null(rhs_cache$getInverse()))
rhs_messages <- character()
cached_rhs <- withCallingHandlers(
        cacheSolve(rhs_cache, b = c(3, 5)),
        message = function(m) {
                rhs_messages <<- c(rhs_messages, conditionMessage(m))
                invokeRestart("muffleMessage")
        }
)
stopifnot(isTRUE(all.equal(cached_rhs, rhs_one)))
stopifnot(length(rhs_messages) == 1L, grepl("getting cached data", rhs_messages[[1]]))
rhs_two <- cacheSolve(rhs_cache, b = c(7, 11))
stopifnot(isTRUE(all.equal(rhs_two, solve(first, b = c(7, 11)))))
stopifnot(!isTRUE(all.equal(rhs_one, rhs_two)), is.null(rhs_cache$getInverse()))

# Positional RHS vectors and matrix right-hand sides are solutions too.
positional_rhs <- cacheSolve(rhs_cache, c(3, 5))
stopifnot(isTRUE(all.equal(positional_rhs, rhs_one)), is.null(rhs_cache$getInverse()))
matrix_rhs <- 2 * diag(2)
assert_equal_matrix(cacheSolve(rhs_cache, b = matrix_rhs), solve(first, b = matrix_rhs))
stopifnot(is.null(rhs_cache$getInverse()))

# An explicit NULL right-hand side is an inverse in base R, not a
# solution to solve(A, b) with an actual vector/matrix b.
null_rhs_cache <- makeCacheMatrix(first)
null_rhs_inverse <- cacheSolve(null_rhs_cache, b = NULL)
assert_equal_matrix(null_rhs_inverse, solve(first))
assert_equal_matrix(null_rhs_cache$getInverse(), null_rhs_inverse)
assert_equal_matrix(suppressMessages(cacheSolve(null_rhs_cache, b = NULL)), null_rhs_inverse)
null_positional_cache <- makeCacheMatrix(first)
assert_equal_matrix(cacheSolve(null_positional_cache, NULL), solve(first))
assert_equal_matrix(null_positional_cache$getInverse(), solve(first))
# A subsequent genuine RHS result must still be hidden from getInverse().
null_rhs_solution <- cacheSolve(null_rhs_cache, b = c(3, 5))
stopifnot(isTRUE(all.equal(null_rhs_solution, solve(first, b = c(3, 5)))))
stopifnot(is.null(null_rhs_cache$getInverse()))

# Named tolerance settings still compute genuine inverses.
rhs_inverse <- cacheSolve(rhs_cache, tol = 0)
assert_equal_matrix(rhs_inverse, solve(first, tol = 0))
assert_equal_matrix(rhs_cache$getInverse(), rhs_inverse)
assert_equal_matrix(suppressMessages(cacheSolve(rhs_cache, tol = 0)), rhs_inverse)

# An explicit manual cache entry retains the original no-options behavior.
stopifnot(identical(sort(names(singular_cache)),
                    sort(c("set", "get", "setInverse", "getInverse"))))
assert_equal_matrix(suppressMessages(cacheSolve(singular_cache)), sentinel)
strict_manual <- tryCatch(cacheSolve(singular_cache, tol = 0.5), error = function(e) e)
stopifnot(inherits(strict_manual, "error"))

# A valid cached inverse must not survive a transition to a singular matrix.
# A failed recomputation must leave the same closure reusable once data recover.
recovering <- makeCacheMatrix(diag(c(2, 3)))
before_failure <- cacheSolve(recovering)
assert_equal_matrix(before_failure, diag(c(0.5, 1 / 3)))
recovering$set(singular)
stopifnot(is.null(recovering$getInverse()))
singular_after_cached <- tryCatch(cacheSolve(recovering), error = function(e) e)
stopifnot(inherits(singular_after_cached, "error"), is.null(recovering$getInverse()))
recovered_matrix <- matrix(c(5, 2, 1, 3), nrow = 2)
recovering$set(recovered_matrix)
recovered_inverse <- cacheSolve(recovering)
assert_equal_matrix(recovered_inverse, solve(recovered_matrix))
stopifnot(!isTRUE(all.equal(before_failure, recovered_inverse)))
assert_equal_matrix(suppressMessages(cacheSolve(recovering)), recovered_inverse)

# R arguments are lazy: the closure must capture the initial matrix even
# when the caller changes its variable before the first get() or cacheSolve().
source_before_first_read <- diag(c(3, 4))
unread_cache <- makeCacheMatrix(source_before_first_read)
source_before_first_read[1, 1] <- 99
assert_equal_matrix(unread_cache$get(), diag(c(3, 4)))
assert_equal_matrix(cacheSolve(unread_cache), diag(c(1 / 3, 1 / 4)))
assert_equal_matrix(suppressMessages(cacheSolve(unread_cache)), unread_cache$getInverse())

# Rebinding the caller's variable must also leave the original matrix intact.
source_before_reassignment <- diag(c(5, 6))
unread_reassigned <- makeCacheMatrix(source_before_reassignment)
source_before_reassignment <- diag(c(10, 12))
assert_equal_matrix(unread_reassigned$get(), diag(c(5, 6)))
assert_equal_matrix(cacheSolve(unread_reassigned), diag(c(1 / 5, 1 / 6)))

# R uses copy-on-modify for matrices: mutating the caller's original variable
# after the closure forced its argument must not change its stored matrix.
source_matrix <- diag(c(3, 4))
snapshot <- makeCacheMatrix(source_matrix)
snapshot_inverse <- cacheSolve(snapshot)
source_matrix[1, 1] <- 99
assert_equal_matrix(snapshot$get(), diag(c(3, 4)))
assert_equal_matrix(suppressMessages(cacheSolve(snapshot)), snapshot_inverse)

# Setter arguments also need eager evaluation: otherwise the stored value
# can follow a later caller-side reassignment before its first read.
replacement <- diag(c(2, 5))
setter_cache <- makeCacheMatrix(diag(2))
setter_cache$set(replacement)
replacement <- diag(c(11, 13))
assert_equal_matrix(setter_cache$get(), diag(c(2, 5)))
assert_equal_matrix(cacheSolve(setter_cache), diag(c(1 / 2, 1 / 5)))

# Mutation of the caller's original matrix after set() must likewise be safe.
replacement <- diag(c(3, 7))
setter_cache$set(replacement)
replacement[1, 1] <- 99
assert_equal_matrix(setter_cache$get(), diag(c(3, 7)))
assert_equal_matrix(cacheSolve(setter_cache), diag(c(1 / 3, 1 / 7)))

# The public manual-cache setter must freeze values as well as the matrix.
manual <- makeCacheMatrix(diag(c(2, 3)))
supplied_inverse <- diag(c(1 / 2, 1 / 3))
manual$setInverse(supplied_inverse)
supplied_inverse <- diag(c(9, 9))
assert_equal_matrix(manual$getInverse(), diag(c(1 / 2, 1 / 3)))
assert_equal_matrix(suppressMessages(cacheSolve(manual)), diag(c(1 / 2, 1 / 3)))

supplied_inverse <- diag(c(1 / 4, 1 / 5))
manual$setInverse(supplied_inverse)
supplied_inverse[1, 1] <- 999
assert_equal_matrix(manual$getInverse(), diag(c(1 / 4, 1 / 5)))

# Manual options also must be frozen, so cache hits match the original solve call.
manual_args <- list(tol = 0)
manual$setInverse(diag(c(1 / 4, 1 / 5)), manual_args)
manual_args$tol <- 0.5
assert_equal_matrix(manual$getInverse(list(tol = 0)), diag(c(1 / 4, 1 / 5)))
stopifnot(is.null(manual$getInverse(list(tol = 0.5))))

# A failed evaluation of manual solve options must not partially overwrite
# the previously valid inverse while leaving the old options in place.
atomic <- makeCacheMatrix(diag(c(2, 3)))
known_inverse <- diag(c(1 / 2, 1 / 3))
atomic$setInverse(known_inverse, list(tol = 0))
bad_manual <- tryCatch(
        atomic$setInverse(diag(c(9, 9)), stop("bad manual options")),
        error = function(e) e
)
stopifnot(inherits(bad_manual, "error"),
          grepl("bad manual options", conditionMessage(bad_manual)))
assert_equal_matrix(atomic$getInverse(list(tol = 0)), known_inverse)
stopifnot(is.null(atomic$getInverse(list(tol = 0.5))))

cat("PASS: eager setter snapshots preserve matrix, inverse, and options\\n")

cat("PASS: singular transition, failed inverse recovery, and caller-side mutation\\n")

cat("PASS: matrix inverse caching, reuse, invalidation, accessors, and errors\n")
cat("PASS: independent closures, matrix dimensions, and forwarded solve tolerance\n")
