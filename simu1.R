                                        # Set working directory:

setwd("~/")


## number of covariables
p <- 3

## learning set size
g <- 13

                                        # true model: X is N(0, 9)
                                        # Y is 0*X + e = e is N(0, 1) - null hypothesis
                                        # data Z = (X, Y)
                                        # it is known that the intercept coefficient is zero
                                        # so that the regression only needs to estimate p regression coefficients.

                                        # returns a sample of size n in p independent covariable dimensions, and a response
                                        # first column contains the Y observations
                                        # other columns contain the X observations
draw <- function(n, p)
    cbind(rnorm(n), matrix(rnorm(n * p, 0, p), nrow=n))

## mean(replicate(1e6,boundedGamma(g=5,draw(n=6,p=3))))
## [1] 0.4837652

cIJiang <- function(data, conf.level=.95, kernelFunction, B=1e3) {
    n <- nrow(data)
    bootstrapSamples <- replicate(B, sample(n, replace=TRUE))
                                        # counts how often each observation is contained in each bootstrap sample
    m <- apply(bootstrapSamples, 2, function(x) sapply(1:n, function (j) sum(j == x)))
                                        # bccb is a vector of length B
    bccb <- apply(
                                        # returns a matrix of dim B times n
        sapply(
            1:n,
                                        # this function returns a vector of length B, containing the value of the loss function on each bootstrap sample, for the leftOutIndex as the test index
            function(leftOutIndex) {
                trainingSets <- apply(bootstrapSamples, 2, function(x) x[x != leftOutIndex])
                m[leftOutIndex, ] * unlist(
                                        lapply(
                                            trainingSets,
                                            function(ts)
                                                kernelFunction(
                                                    g = length(ts),
                                                    data=rbind(data[ts, ], data[leftOutIndex, ])
                                                )
                                        )
                                    )
            }
        ),
        1,
        mean
    )
                                        # uncorrected bootstrap point estimate of the error
    bccv <- mean(bccb)
                                        # upper limit
    bccvpUpper <- sort(bccb)[floor(B*(1 - (1-conf.level)/2))]
                                        # point estimate
    point <- sort(bccb)[floor(B*.5)]
                                        # lower limit
    bccvpLower <- sort(bccb)[max(1, floor(B*(1-conf.level)/2))]


                                        # Jiang's correction
                                        # usual loo point estimate of the error
    loocv <- mean(
        sapply(
            1:n,
            function(i)
                kernelFunction(g=n-1, data=data[c((1:n)[-i], i), ])
        )
    )

                                        # corrected upper bound
    cub <- bccvpUpper - bccv + loocv
                                        # corrected point estimate
    cpe <- point - bccv + loocv
                                        # corrected lower bound
    clb <- bccvpLower - bccv + loocv
    return(c(clb, cpe, cub))
}


boundedGamma <- function(g, data) {
    require("RcppArmadillo")
    stopifnot(nrow(data) == g + 1)
    coef <- fastLmPure(as.matrix(data[1:g, -1]), data[1:g, 1])$coefficients
    atan(as.numeric(data[g + 1, 1] - (sum(coef * data[g + 1, -1])))^2)*2/pi
}

                                        # this function takes a sample of size n and returns two elementary estimators of the error rate, on two non-overlapping random subsamples
GammaPair <- function(g, data, kernelFunction) {
    n <- nrow(data)

    stopifnot(n >= 2 * g + 2)
    indices <- sample(n, 2 * (g + 1))
    c(
        kernelFunction(
            g=g,
            data=data[indices[1:(g+1)], ]
        ),
        kernelFunction(
            g=g,
            data=data[indices[(g+2):(2*(g + 1))], ]
        )
    )
}



shifter <- function(x, n) {
    if (n == 0) x else c(tail(x, -n), head(x, n))
}


                                        # this function takes a sample of size n and returns two elementary estimators of the error rate, on two random subsamples overlapping on one but not more observations
GammaPairOverlapOne <- function(g, data, kernelFunction, overlapsize=1) {
    n <- nrow(data)

    stopifnot(n >= 2 * g + 2)
                                        # replace is FALSE by default
    indices <- sample(n, 2 * (g + 1))
    indicesOne  <- indices[1:(g+1)]
                                        # here, we start at g+1, so we use exactly one observation from evaluationOne again
    indicesTwo  <- indices[(g+1 - overlapsize + 1):(2 * g + 1  - overlapsize + 1)]

                                        # the symmetrized kernel on the indices indicesOne
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

                                        # the confidence interval for the mean
cI <- function(g, data, N, conf.level=.95, kernelFunction) {
    n <- nrow(data)

    stopifnot(g <= (n + 2)/2)
    alpha <- 1 - conf.level


                                        # qt(1-alpha/2, df=n-1)
    q <- qnorm(1 - alpha/2)

    mcresults <- replicate(
        N,
        GammaPair(
            g=g,
            data=data,
            kernelFunction=kernelFunction
        )
    )

    meanhat <- mean(as.vector(mcresults))


                                        # note that since we use the "immediate trick", we do not have to divide by n here
    vhat <- meanhat^2 - mean(mcresults[1, ] * mcresults[2, ])

                                        #    if (kernelFunction == Gamma) cat(paste("True error is", 1 + p/(g - p - 1)), "\n")
    return(
        c(
            meanhat - q * sqrt(vhat),
            meanhat,
            meanhat + q * sqrt(vhat)
        )
    )
}

                                        # the confidence interval for the mean, using the u(n) part ii from the paper, as suggested by the referee
cIunii <- function(g, data, N, conf.level=.95, kernelFunction) {
    n <- nrow(data)

    stopifnot(g <= (n + 2)/2)
    alpha <- 1 - conf.level
                                        # qt(1-alpha/2, df=n-1)
    q <- qnorm(1 - alpha/2)

    ## mcresults <- replicate(
    ##     N,
    ##     GammaPair(
    ##         g=g,
    ##         data=data,
    ##         kernelFunction=kernelFunction
    ##     )
    ## )

    mcresultsOverlapOne <- replicate(
        N,
        GammaPairOverlapOne(
            g=g,
            data=data,
            kernelFunction=kernelFunction,
            overlapsize=1
        )
    )

                                        # the U statistic itself is just computed by averaging as many Gamma evaluations as possible
    meanhat <- mean(c(as.vector(mcresultsOverlapOne)))

                                        # to estimate kappa_1, we need to average out a lot of evaluations of the product of two symmetrized kernels with overlap one
                                        # the difference term is the estimated covariance between two evaluations of the symmetrized kernel on sets with overlap one

                                        # vhat <- (g + 1)^2 / n * (mean(mcresultsOverlapOne[1, ] * mcresultsOverlapOne[2, ]) - mean(mcresults[1, ] * mcresults[2, ]))
                                        # an alternative
    vhat <- (g + 1)^2 / n * cov(mcresultsOverlapOne[1, ], mcresultsOverlapOne[2, ])

                                        #    if (kernelFunction == Gamma) cat(paste("True error is", 1 + p/(g - p - 1)), "\n")
    return(
        c(
            meanhat - q * sqrt(vhat),
            meanhat,
            meanhat + q * sqrt(vhat)
        )
    )
}

                                        # using bootstrapping
cIboot <- function(g, data, N, conf.level = 0.95, kernelFunction) {
    require("boot")
    n <- nrow(data)

    stopifnot(g <= (n + 2)/2)
    alpha <- 1 - conf.level

    computeStatistics <- function(data, index, g, N, conf.level=.95, kernelFunction)
    {
        mean(as.vector(replicate(
            N,
            GammaPair(
                g=g,
                data=data[index, ],
                kernelFunction=kernelFunction
            )
        )))
    }
    tmp <- boot(data, computeStatistics, R = 999, sim = 'ordinary',
                g = g, N = N, conf.level = .95, kernelFunction = kernelFunction)

    interv <- boot.ci(tmp, type = "basic")
    return(
        c(
            interv$basic[4],
            interv$t0,
            interv$basic[5]
        )
    )

}





## to approximate the true error
                                        # returns 0.48 plus/minus 0.001 for g=5
                                        # returns 0.3835935 plus/minus 0.01 for g=13
tE <- function (g, p) {
    require("RcppArmadillo")
    k <- rep(NA, 1e6)
                                        # "burn in"
    for (i in 1:10) {
        learnData <- draw(n=g, p=p)
        coef <- fastLmPure(as.matrix(learnData[, -1]), learnData[, 1])$coefficients
        testData <- draw(n=100000, p=p)
        k[i] <- mean(atan(as.numeric(testData[, 1] - testData[, -1] %*% coef)^2)*2/pi)
    }

    certainty <- 1/0
    currentIndex <- 11
    while (certainty > 1e-2) {
        print(certainty)

        learnData <- draw(n=g, p=p)
        coef <- fastLmPure(as.matrix(learnData[, -1]), learnData[, 1])$coefficients
        testData <- draw(n=100000, p=p)
        k[currentIndex] <- mean(atan(as.numeric(testData[, 1] - testData[, -1] %*% coef)^2)*2/pi)
        ci <- t.test(k, na.rm=TRUE)$conf.int
        certainty <- ci[2] - ci[1]
        currentIndex <- currentIndex + 1
    }
    mean(k, na.rm=TRUE)
}



####################
## the simulation ##
####################

set.seed(1234)
seeds <- sample(1e4:1e6, size=7000)


                                        # for N = 5e5
## number of resampling iterations for the confidence interval, should be at least 1e5
N <- 5e5




                                        # Start the cluster:

cl <- makeCluster()


                                        # Export the objects in the workspace to the
                                        # parallel jobs:

clusterExport(cl, list=ls())




if (file.exists("./SimulationRevision/log1_N5e5.txt")) system("./SimulationRevision/rm log1_N5e5.txt")
time1 <- system.time(simulationResults <- t(
                         parSapply(cl, 1:7,
                                   function(i) {

                                       set.seed(seeds[i])

                                       n <- floor((i-1)/1000) * 5 + 30

                                       write(paste("repetition ", i), file="./SimulationRevision/log1_N5e5.txt", append=TRUE)
                                       print(i)

                                       sple <- draw(n=n, p=p)
                                       interval <- cI(g=g, data=sple, N=N, kernelFunction=boundedGamma)

                                       if (any(is.na(interval)))
                                           warning("NAs produced. You should choose a higher N!")

                                       line <-  paste(i, n, paste(interval, collapse=" "))

                                       write(line, file="./SimulationRevision/aggregateOurIntg13_N5e5.csv", append=TRUE)
                                       interval
                                   }
                                   )
                     )

                     )

save(time1, file = './SimulationRevision/time1_N5e5.R')


if (file.exists("./SimulationRevision/log2_N5e5.txt")) system("./SimulationRevision/rm log2_N5e5.txt")
time2 <- system.time(simulationResultsJiangWithnEqualg <- t(
                         parSapply(cl, 1:5,
                                   function(i) {

                                       write(paste("repetition ", i), file="./SimulationRevision/log2_N5e5.txt", append=TRUE)
                                       print(i)

                                       sple2 <- draw(n=g, p=p)
                                       intervalJiang2 <- cIJiang(data=sple2, kernelFunction=boundedGamma, B=1e4)

                                       line <- paste(i, paste(intervalJiang2, collapse=" "))
                                       write(line, file="./SimulationRevision/aggregateJiangg13_N5e5.csv", append=TRUE)
                                       intervalJiang2
                                   }
                                   )
                     )

                     )

save(time2, file = './SimulationRevision/time2_N5e5.R')


if (file.exists("./SimulationRevision/log3_N5e5.txt")) system("./SimulationRevision/rm log3_N5e5.txt")
time3 <- system.time(simulationResults <- t(
                         parSapply(cl, 1:7,
                                   function(i) {

                                       set.seed(seeds[i])

                                       n <- floor((i-1)/1000) * 5 + 30

                                       write(paste("repetition ", i), file="./SimulationRevision/log3_N5e5.txt", append=TRUE)
                                       print(i)

                                       sple <- draw(n=n, p=p)
                                       interval <- cIunii(g=g, data=sple, N=N, kernelFunction=boundedGamma)

                                       if (any(is.na(interval)))
                                           warning("NAs produced. You should choose a higher N!")

                                       line <-  paste(i, n, paste(interval, collapse=" "))

                                       write(line, file="./SimulationRevision/aggregateIunii_N5e5.csv", append=TRUE)
                                       interval
                                   }
                                   )
                     )

                     )

save(time3, file = './SimulationRevision/time3_N5e5.R')


if (file.exists("./SimulationRevision/log4_N5e5.txt")) system("./SimulationRevision/rm log4_N5e5.txt")
time4 <- system.time(simulationResults <- t(
                         parSapply(cl, 1:7,
                                   function(i) {

                                       set.seed(seeds[i])

                                       n <- floor((i-1)/1000) * 5 + 30

                                       write(paste("repetition ", i), file="./SimulationRevision/log4_N5e5.txt", append=TRUE)
                                       print(i)

                                       sple <- draw(n=n, p=p)
                                       interval <- cIboot(g=g, data=sple, N=N, kernelFunction=boundedGamma)

                                       if (any(is.na(interval)))
                                           warning("NAs produced. You should choose a higher N!")

                                       line <-  paste(i, n, paste(interval, collapse=" "))

                                       write(line, file="./SimulationRevision/aggregatecIboot_N5e5.csv", append=TRUE)
                                       interval
                                   }
                                   )
                     )

                     )

save(time4, file = './SimulationRevision/time4_N5e5.R')



                                        # for N = 5e4
## sacrificing accuracy, as asked by a reviewer
N <- 5e4


sfExportAll()

if (file.exists("./SimulationRevision/log1_N5e4.txt")) system("./SimulationRevision/rm log1_N5e4.txt")
time1b <- system.time(simulationResults <- t(
                          parSapply(cl, 1:7,
                                    function(i) {

                                        set.seed(seeds[i])

                                        n <- floor((i-1)/1000) * 5 + 30

                                        write(paste("repetition ", i), file="./SimulationRevision/log1_N5e4.txt", append=TRUE)
                                        print(i)

                                        sple <- draw(n=n, p=p)
                                        interval <- cI(g=g, data=sple, N=N, kernelFunction=boundedGamma)

                                        if (any(is.na(interval)))
                                            warning("NAs produced. You should choose a higher N!")

                                        line <-  paste(i, n, paste(interval, collapse=" "))

                                        write(line, file="./SimulationRevision/aggregateOurIntg13_N5e4.csv", append=TRUE)
                                        interval
                                    }
                                    )
                      )

                      )

save(time1b, file = './SimulationRevision/time1_N4e5.R')


if (file.exists("./SimulationRevision/log2_N5e4.txt")) system("./SimulationRevision/rm log2_N5e4.txt")
time2b <- system.time(simulationResultsJiangWithnEqualg <- t(
                          parSapply(cl, 1:5,
                                    function(i) {

                                        write(paste("repetition ", i), file="./SimulationRevision/log2_N5e4.txt", append=TRUE)
                                        print(i)

                                        sple2 <- draw(n=g, p=p)
                                        intervalJiang2 <- cIJiang(data=sple2, kernelFunction=boundedGamma, B=1e4)

                                        line <- paste(i, paste(intervalJiang2, collapse=" "))
                                        write(line, file="./SimulationRevision/aggregateJiangg13_N5e4.csv", append=TRUE)
                                        intervalJiang2
                                    }
                                    )
                      )

                      )

save(time2b, file = './SimulationRevision/time2_N4e5.R')


if (file.exists("./SimulationRevision/log3_N5e4.txt")) system("./SimulationRevision/rm log3_N5e4.txt")
time3b <- system.time(simulationResults <- t(
                          parSapply(cl, 1:7,
                                    function(i) {

                                        set.seed(seeds[i])

                                        n <- floor((i-1)/1000) * 5 + 30

                                        write(paste("repetition ", i), file="./SimulationRevision/log3_N5e4.txt", append=TRUE)
                                        print(i)

                                        sple <- draw(n=n, p=p)
                                        interval <- cIunii(g=g, data=sple, N=N, kernelFunction=boundedGamma)

                                        if (any(is.na(interval)))
                                            warning("NAs produced. You should choose a higher N!")

                                        line <-  paste(i, n, paste(interval, collapse=" "))

                                        write(line, file="./SimulationRevision/aggregateIunii_N5e4.csv", append=TRUE)
                                        interval
                                    }
                                    )
                      )

                      )

save(time3b, file = './SimulationRevision/time3_N4e5.R')

if (file.exists("./SimulationRevision/log4_N5e4.txt")) system("./SimulationRevision/rm log4_N5e4.txt")
time4b <- system.time(simulationResults <- t(
                          parSapply(cl, 1:7,
                                    function(i) {

                                        set.seed(seeds[i])

                                        n <- floor((i-1)/1000) * 5 + 30

                                        write(paste("repetition ", i), file="./SimulationRevision/log4_N5e4.txt", append=TRUE)
                                        print(i)

                                        sple <- draw(n=n, p=p)
                                        interval <- cIboot(g=g, data=sple, N=N, kernelFunction=boundedGamma)

                                        if (any(is.na(interval)))
                                            warning("NAs produced. You should choose a higher N!")

                                        line <-  paste(i, n, paste(interval, collapse=" "))

                                        write(line, file="./SimulationRevision/aggregatecIboot_N5e4.csv", append=TRUE)
                                        interval
                                    }
                                    )
                      )

                      )

save(time4b, file = './SimulationRevision/time4_N4e5.R')


save(time1, time2, time3, time4, time1b, time2b, time3b, time4b, file='./SimulationRevision/computationalTimes.Rdata')
