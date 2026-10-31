## Matrix-inverse caching via lexical scoping.
## The cache is cleared when the underlying matrix changes.

## Return a matrix-like object with accessors for the matrix and its cached inverse.
makeCacheMatrix <- function(x = matrix()) {
        ## Evaluate the caller's matrix now, before its variable can change.
        ## Otherwise R's lazy promise can resolve only on the first get().
        force(x)
        inverse <- NULL

        set <- function(y) {
                x <<- y
                inverse <<- NULL
        }
        get <- function() x
        setInverse <- function(value) {
                inverse <<- value
        }
        getInverse <- function() inverse

        list(
                set = set,
                get = get,
                setInverse = setInverse,
                getInverse = getInverse
        )
}


## Return a cached inverse if available; otherwise compute and store it.
cacheSolve <- function(x, ...) {
        inverse <- x$getInverse()
        if (!is.null(inverse)) {
                message("getting cached data")
                return(inverse)
        }

        inverse <- solve(x$get(), ...)
        x$setInverse(inverse)
        inverse
}
