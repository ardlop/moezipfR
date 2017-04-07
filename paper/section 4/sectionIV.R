###########################################
# 
# Usage examples of the moezipfR package.
#
# @author  A. Duarte-López <aduarte@ac.upc.edu>
# @version 0.0.1, 07/04/2017
#
###########################################

library(moezipfR)
set.seed(2882)

# Figure 1. 
N <- 100
alpha <- 1.8
betas <- c(0.2, 0.5, 9, 15)
options(scipen=999)
pdf('density_beta.pdf')
plot(dmoezipf(1:N, alpha, beta = 1), log='xy', type='l', lwd = 1.3, lty=1, col=1,
     xlab='x', ylab='Prob', yaxt='n')
axis(2,at=c(0.0001, 0.0005, 0.0050, 0.0500, 0.500),labels=T)

title(expression(paste("MOEZipf(", alpha,' = 1.8, ', beta,")")))
for(b in 1:length(betas)){
  lines(dmoezipf(1:N, alpha, beta = betas[b]), type='l', lwd = 1.3, lty=b+1, col=b+1)
}
graphics::legend("topright",  legend = c(expression(paste(beta,' = 0.2')),
                                         expression(paste(beta,' = 0.5')),
                                         expression(paste(beta,' = 1')),
                                         expression(paste(beta,' = 9')),
                                         expression(paste(beta,' = 15'))),
                 col=c(2,3,1,4,5),
                 lwd= rep(2, times=length(betas)),
                 lty = c(2,3,1,4,5),
                 bty = "n")
dev.off()

#Distribution functions
smoezipf(1:10, alpha = 3.5, beta = 0.75, show.plot = F)
pmoezipf(1:10, alpha = 3.5, beta = 0.75, show.plot = F)
moezipfR.moments(2, alpha = 3.5, beta = 0.75, tolerance = 10^(-6))


# Figure 2. The expected value as function of alpha and beta.
alphas <- seq(2.0, 5, by = 0.1)
betas <- c(0.5, 1, 1.5, 3)
alphaExpctValues <- matrix(ncol = length(alphas), nrow=length(betas))

alphas1 <- c(2.5, 3, 5, 10)
betas1 <- seq(0.2, 6, by = 0.2)
betasExpctValues <- matrix(ncol = length(betas1), nrow=length(alphas1))

for(i in 1:length(betas)){
  print(betas[i])
  alphaExpctValues[i, ] <- sapply(alphas, moezipfR.mean, beta = betas[i], tolerance = 10^(-4))
  betasExpctValues[i, ] <- sapply(betas1, moezipfR.mean, alpha = alphas1[i], tolerance = 10^(-4))
}

pdf('mean_beta.pdf')
plot(alphas, alphaExpctValues[1, ],main=bquote("Expected value as function of "~beta),
     xlab=bquote(alpha),ylab="Expected Value", type = 'l', ylim = c(1, 9))
for(i in 2:length(betas)){
  lines(alphas, alphaExpctValues[i,], col=i, lty=i)
}
legend("topright", legend=c(as.expression(bquote(beta == .(betas[1]))),
                            as.expression(bquote(beta == .(betas[2]))),
                            as.expression(bquote(beta == .(betas[3]))),
                            as.expression(bquote(beta == .(betas[4])))),
       col=1:4,lty=1:4, cex=1, y.intersp=0.8, bty='n')

dev.off()

pdf('mean_alpha.pdf')
plot(betas1, betasExpctValues[1, ],main=bquote("Expected value as function of "~alpha),
     xlab=bquote(beta),ylab="Expected Value", type = 'l', ylim = c(1, 6))
for(i in 2:length(alphas1)){
  lines(betas1, betasExpctValues[i,], col=i, lty=i)
}
legend("topleft", legend=c(as.expression(bquote(alpha == .(alphas1[1]))),
                           as.expression(bquote(alpha == .(alphas1[2]))),
                           as.expression(bquote(alpha == .(alphas1[3]))),
                           as.expression(bquote(alpha == .(alphas1[4])))),
       col=1:4,lty=1:4, cex=1, y.intersp=0.8, bty='n')
dev.off()

# Functions used to fit the MOEZipf function and to calculate the maximum log-likelihood.
data <- rmoezipf(n = 10000, alpha = 2.5, beta = 1.3)
dataMatrix <- moezipfR.utils.getDataMatrix(data)
(mloglik <- moezipfR.loglikelihood(dataMatrix, alpha = 2.5, beta = 1.3))

(initial <- moezipfR.utils.getInitialValues(dataMatrix))
fit <- moezipfR.fit(dataMatrix, init_alpha = initial$init_alpha, init_beta = initial$init_beta)
print(fit)

# Zipf Log-likelihood fn
get_AIC <- function(loglike, K) {
  -2*loglike + 2*K
}

zeta_mle <- function(alpha, nSize, freq, values){
  -( -alpha*sum(freq * log(values)) - nSize*log(VGAM::zeta(alpha)))
}

zipf_pmf <- function(k, alpha){
  (k^(-alpha))/VGAM::zeta(alpha)
}

# Application 1: Out-degree of the Patents network.
data <- read.table(file = '/home/ariel/Documents/ac/recerca/trunk/temporal/degreeAnalysis/datasets/outPatentsDegrees.txt', header=F)

#MOEZipf fit.
dataMatrix <-moezipfR.utils.getDataMatrix(data$V1)
initVal <- moezipfR.utils.getInitialValues(dataMatrix)
estimation <- moezipfR.fit(dataMatrix, init_alpha = initVal$init_alpha, init_beta = initVal$init_beta)

#Zipf fit.
M <- sum(data$V1)
N <- length(data$V1)
M_prime <-sum(log(data$V1))

zipf.est  <- stats::optim(initVal$init_alpha, zeta_mle,
                          nSize =N, freq = dataMatrix[, 2],
                          values = dataMatrix[, 1], hessian = TRUE)

zipf.sd <- sqrt(solve(zipf.est$hessian))
offset <-round(stats::qnorm(1-(0.05/2)), 2)*zipf.sd
zipf.lwAlpha <- zipf.est$par - offset*zipf.sd
zipf.upAlpha <- zipf.est$par + offset*zipf.sd

zipf.llike<- - zipf.est$value
zipf.AIC <- get_AIC(zipf.llike, 1)

zipf.fitted <- N*sapply(dataMatrix[, 1], zipf_pmf, alpha = as.numeric(zipf.est$par[1]))

pdf('outPatentsFit.pdf')
graphics::plot(dataMatrix[, 1], dataMatrix[, 2], log="xy",
               xlab="ln(Observation)", ylab="ln(Frequency)",
               main="Patents network (out-degree)", type = 'l', lwd = 2)

graphics::lines(dataMatrix[, 1], fitted(estimation), col="blue", lty = 2, lwd = 2)
graphics::lines(dataMatrix[, 1], zipf.fitted, col="red", lty = 3, lwd = 2)

graphics::legend("topright",  legend = c('Observations', 'Zipf Distribution', 'MOEZipf Distribution'),
                 col=c('black', 'red', 'blue'), lty=c(1, 3, 2), lwd=c(2, 2, 2),
                 bty = 'n')

dev.off()

## Application 2: Out-degree of the email network.
data <- read.table(file ='/home/ariel/Documents/ac/recerca/trunk/ZPSS/datasets/degrees_Final/outEmailEuAll_zeros.txt', header=F)
data <- subset(data, V1 != 0)

#MOEZipf fit.
dataMatrix <-moezipfR.utils.getDataMatrix(data$V1)
initVal <- moezipfR.utils.getInitialValues(dataMatrix)
estimation <- moezipfR.fit(dataMatrix, init_alpha = initVal$init_alpha, init_beta = initVal$init_beta)

#Zipf fit.
M <- sum(data$V1)
N <- length(data$V1)
M_prime <-sum(log(data$V1))

zipf.est  <- stats::optim(initVal$init_alpha, zeta_mle,
                       nSize =N, freq = dataMatrix[, 2],
                       values = dataMatrix[, 1], hessian = TRUE)

zipf.sd <- sqrt(solve(zipf.est$hessian))
offset <-round(stats::qnorm(1-(0.05/2)), 2)*zipf.sd
zipf.lwAlpha <- zipf.est$par - offset*zipf.sd
zipf.upAlpha <- zipf.est$par + offset*zipf.sd

zipf.llike<- - zipf.est$value
zipf.AIC <- get_AIC(zipf.llike, 1)

zipf.fitted <- N*sapply(dataMatrix[, 1], zipf_pmf, alpha = as.numeric(zipf.est$par[1]))

pdf('outEmailEuAll.pdf')
graphics::plot(dataMatrix[, 1], dataMatrix[,2], log="xy",
               xlab="ln(Observation)", ylab="ln(Frequency)",
               main="Email network (out-degree)", type ='l', lwd = 2)

graphics::lines(dataMatrix[, 1], fitted(estimation), col="blue", lty = 2, lwd = 2)
graphics::lines(dataMatrix[, 1], zipf.fitted, col="red", lty = 3, lwd = 2)

graphics::legend("topright",  legend = c('Observations', 'Zipf Distribution', 'MOEZipf Distribution'),
                 col=c('black', 'red', 'blue'), lty=c(1, 3, 2), lwd=c(2, 2, 2),
                 bty = 'n')

dev.off()