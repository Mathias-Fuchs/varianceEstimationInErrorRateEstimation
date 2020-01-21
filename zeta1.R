require(RcppArmadillo)
draw <- function(n, p) cbind(rnorm(n), matrix(rnorm(n * p, 0, p), nrow=n))
boundedGamma <- function(g, data) {
     stopifnot(nrow(data) == g + 1)
     coef <- fastLmPure(as.matrix(data[1:g, -1]), data[1:g, 1])$coefficients
     atan(as.numeric(data[g + 1, 1] - (sum(coef * data[g + 1, -1])))^2)*2/pi
 }
shifter <- function(x, n) if (n == 0) x else c(tail(x, -n), head(x, n))
GammaPairOverlapOne <- function(g, data, kernelFunction, overlapsize=1) {
    n <- nrow(data)
    stopifnot(n >= 2 * g + 2)
    indices <- sample(n, 2 * (g + 1))
    indicesOne  <- indices[1:(g+1)]
    indicesTwo  <- indices[(g+1 - overlapsize + 1):(2 * g + 1  - overlapsize + 1)]    
    evaluationOne <- mean(sapply(0:g, function(i) 
        kernelFunction(
            g=g,
            data=data[shifter(indicesOne, i), ]
        )
        )
        )
    evaluationTwo <- mean(sapply(0:g, function(i) 
        kernelFunction(
            g=g,
            data=data[shifter(indicesTwo, i), ]
        )
        )
        )
    c(evaluationOne, evaluationTwo)
}
g <- 13
p <- 3
for (n in c(30, 35, 40, 50)) {
    data <- draw(n, p)
    d <- replicate(1e5, GammaPairOverlapOne(g, data, boundedGamma))
    print(cov(d[1, ], d[2,]))
}
