## Matrix-inverse caching via lexical scoping.
## The cache is cleared when the underlying matrix changes.

## Return a matrix-like object with accessors for the matrix and its cached inverse.
makeCacheMatrix <- function(x = matrix()) {
        ## Evaluate the caller's matrix now, before its variable can change.
        ## Otherwise R's lazy promise can resolve only on the first get().
        force(x)
        inverse <- NULL
        inverseArgs <- list()

        set <- function(y) {
                x <<- y
                inverse <<- NULL
                inverseArgs <<- list()
        }
        get <- function() x
        ## Keep the original four accessors. An optional argument lets
        ## cacheSolve check whether the stored inverse used the same solve()
        ## options, while getInverse() without arguments still returns the
        ## ordinary cached value for assignment callers.
        setInverse <- function(value, solveArgs = list()) {
                inverse <<- value
                inverseArgs <<- solveArgs
        }
        getInverse <- function(solveArgs) {
                if (missing(solveArgs) || identical(solveArgs, inverseArgs)) {
                        return(inverse)
                }
                NULL
        }

        list(
                set = set,
                get = get,
                setInverse = setInverse,
                getInverse = getInverse
        )
}


## Return a cached inverse if available; otherwise compute and store it.
cacheSolve <- function(x, ...) {
        solveArgs <- list(...)
        inverse <- x$getInverse(solveArgs)
        if (!is.null(inverse)) {
                message("getting cached data")
                return(inverse)
        }

        ## A changed tolerance or other solve() option must be checked by
        ## solve() rather than silently returning a previous calculation.
        ## Store the options only after solve() has succeeded.
        inverse <- do.call(solve, c(list(x$get()), solveArgs))
        x$setInverse(inverse, solveArgs)
        inverse
}
