
library(ggplot2)
library(reshape)
library(scales)

# as assessed by function tE in supplementary code, has been determined up to a confidence interval of width 1e-3
trueError <- .3836

simulationResultsJiang <- read.csv("aggregateJiangg13.csv", sep=" ", header=FALSE)
simulationResultsJiang <- simulationResultsJiang[, -1]
colnames(simulationResultsJiang) <- c("lj", "pj", "uj")
simulationResultsJiang <- cbind(simulationResultsJiang, simulationResultsJiang$lj <= trueError & simulationResultsJiang$uj >= trueError)
colnames(simulationResultsJiang)[4] <- "success"
simulationResultsJiang <- simulationResultsJiang[complete.cases(simulationResultsJiang), ]


cat(paste("average length of Jiang's CI:", mean(simulationResultsJiang[, 3] - simulationResultsJiang[, 1]), "\n"))


k <- binom.test(
    x=sum(simulationResultsJiang$success),
    n=length(simulationResultsJiang$success)
)

simulationResults <- read.csv("aggregateOurIntg13.csv", sep=" ", header=FALSE)
simulationResults <- simulationResults[complete.cases(simulationResults), ]


                                        # sort by runId, just to keep it nice
simulationResults <- simulationResults[order(simulationResults$V1), ]

                                        # remove the runId
simulationResults <- simulationResults[, -1]


colnames(simulationResults) <- c("sample_size", "l", "p", "u")

cat(paste("average length of our CI:", mean(simulationResults[, 4] - simulationResults[, 2]), "\n"))



simulationResults <- cbind(simulationResults, simulationResults$l <= trueError & simulationResults$u >= trueError)

colnames(simulationResults)[5] <- "successAsymptExactMethod"

                                        # now, we need to manufacture data in the format of a dataframe with columns c("sample_size", "whose_method", "lower_conf_limit_for_coverage_probability", "point", "upper_conf_limit_for_coverage_probability")

ourConfInt <- cbind(method="asymptotically.exact", simulationResults[, names(simulationResults) %in% c("sample_size", "successAsymptExactMethod")])
colnames(ourConfInt) <- c("method", "sample_size", "coverage_yes_no")


                                        
                                        # we need a confidence interval function

cI2 <- function(x) {
    x <- x[!is.na(x)]
    cd <- try(binom.test(x=sum(x), n=length(x)), TRUE)
    if (!(inherits(cd, "try-error")))  return(c(cd$conf.int[1], as.numeric(cd$estimate), cd$conf.int[2]))
    else return(c(0, .5, 1))
}


c <- cast(ourConfInt, formula=sample_size~., fun.aggregate=cI2, value="coverage_yes_no")
c <- as.data.frame(c)



names(c) <-  c("sample_size", "lower_conf_limit_for_coverage_probability", "point", "upper_conf_limit_for_coverage_probability")
c <- subset(c, !is.na(sample_size))



c$Jianglower <- k$conf.int[1]
c$Jiangpoint <- k$estimate
c$Jiangupper <- k$conf.int[2]

require(scales)

ggplot(
    data = c,
    aes(
        x = sample_size,
        y = point
                                        #        group = factor(method),
                                        #        colour = factor(method)
    )
) +
    geom_ribbon(aes(ymin=lower_conf_limit_for_coverage_probability, ymax=upper_conf_limit_for_coverage_probability), fill="grey") +
        geom_line() +
            geom_ribbon(aes(ymin=Jianglower, ymax=Jiangupper), fill="lightgrey") +
                geom_line(aes(x=sample_size, y=Jiangpoint)) +
                    xlab("sample size n") +
                        ylab("coverage probability") +
                            theme(aspect.ratio=1, legend.position="none") +
                                    geom_line(aes(x=sample_size, y=.95, cex=17, legend=NULL)) +
                                        scale_y_continuous(labels=percent, breaks=seq(.9, 1.01, length.out=12)) +
                                      #  scale_y_continuous(labels=percent,limits=c(.9, 1), breaks = seq(.9, 1, length.out=11)) +
                                           scale_x_continuous(breaks=c$sample_size)

        
        
ggsave("results.eps")
