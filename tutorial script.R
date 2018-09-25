##################### Everyday Bayes Script #####################
#	A brief tutorial on Bayesian statistics	first		
#			given at Northern Arizona University on
#			September 28, 2018			
#							
#	Created: 	9/15/2018				
#	Modified:	9/23/2018				
#							
#	By:			MJ Meyer		
#							
#################################################################

### install libraries ###
# ONLY NEED TO INSTALL ONCE # 
install.packages('MCMCpack')
install.packages('ISwR')
install.packages('mcmcplots')
install.packages('mvtnorm')

## load libraries ###
library(MCMCpack)
library(ISwR)
library(mcmcplots)
library(mvtnorm)

## Each of the following examples assumes the data has been cleaned,
##	i.e. there are no missing values in any of the datasets

### One Proportion ###

## Example I.1 ##

# original data #
?graft.vs.host
X			<- with(graft.vs.host, sum(gvhd))	# number of successes
N			<- with(graft.vs.host, length(gvhd))	# number of trials
B			<- 20000				# number of posterior samples

postSamp	<- rbeta(B, X + 1, N - X + 1)	# we add 1 to X and N - X to fully determine
						#	the posterior distribution

# summaries #
quantile(postSamp, probs = c(0.5, 0.025, 0.975))	# posterior median and credible intervals

# slide 23 graph
thetaMat			<- matrix(postSamp, ncol = 1)
colnames(thetaMat)	<- 'theta' # or 'Proportion with GVHD'
denplot(thetaMat, greek = TRUE, lwd = 3)

# more informative graph
thetaMat			<- matrix(postSamp, ncol = 1)
colnames(thetaMat)	<- 'Posterior Distribution of Proportion with GVHD'
denplot(thetaMat, lwd = 3)

# add credible interval---not the most elegant graph
thetaMat			<- matrix(postSamp, ncol = 1)
colnames(thetaMat)	<- 'Posterior Distribution of Proportion with GVHD'
denplot(thetaMat, lwd = 3, ci = 0.95)

# using original data for prior on "new" data #
X			<- 83						# number of successes from new study
X0			<- with(graft.vs.host, sum(gvhd))		# number of successes from prior study
N			<- 200						# number of trials from new study
N0			<- with(graft.vs.host, length(gvhd))		# number of trials from prior study
B			<- 20000					# number of posterior samples

infoPost	<- rbeta(B, X + X0 + 1, (N + N0) - (X + X0) + 1)

# summaries #
quantile(infoPost, probs = c(0.5, 0.025, 0.975))	# posterior median and credible intervals

# slide 28 graph
thetaMat			<- cbind(postSamp, infoPost)
colnames(thetaMat)	<- c('Informative Prior', 'Resulting Posterior')
denplot(thetaMat, lwd = 3, xlim = c(0.1, 0.8))


## Exercise ##

# 	Check the influence of the prior by re-running the informative model setting
#		X0 and N0 to 0. Compare your results from the two runs. How much does the prior
#		influence the results?



### Paired-t and t-test Equivalent ###

?sleep # data description
X	<- with(sleep, extra[group == 1] - extra[group == 2]) 	# differences for paired, for
								# 	t, use the vector of data
B	<- 40000						# number of posterior samples

burnin	<- B/2							# typically drop first half of
								#	samples when using Gibbs
N		<- length(X)					# sample size

mu		<- vector('numeric', length = B)	# an empty vector to store mu
sigma	<- vector('numeric', length = B)	# an empty vector to store sigma

mu[1]		<- mean(X)
sigma[1]	<- var(X)

for(b in 2:B){
	# draw mu #
	mu[b]		<- rnorm(1, mean(X), sqrt(sigma[b-1]/N))
	
	# draw sigma #
	sigma[b]	<- rinvgamma(1, N/2, (1/2)*sum((X - mu[b-1])^2))
}

# drop first half #
mub		<- mu[-c(1:burnin)]
sigmab	<- sigma[-c(1:burnin)]

quantile(mub, probs = c(0.5, 0.025, 0.975))	# posterior median and credible intervals
quantile(sigmab, probs = c(0.5, 0.025, 0.975))	# posterior median and credible intervals

mean(mu > 0)	# posterior probability mean difference is greater than 0

postMat	<- cbind(mu, sigma)
colnames(postMat)	<- c('Mean Difference', 'Variance of Differences')
denplot(postMat, lwd = 3)


## Exercise ##

# Use your own data to conduct a one sample t-equivalent. Test a null hypothesis
#	of your chosing. How would you conduct this test using Bayes?

# If you're having trouble thinking of a one sample data problem, use the dataset
#	bp.obese which contains data from a random sample of adults in a small California
#	town. Test to see if the mean systolic blood pressure from this group is equal to
#	the normal blood pressure level of 120. (The variable name is bp. use
#	X	<- with(bp.obese, bp) to define X.)



### Two sample t Equivalent ###

?CO2
X1	<- with(CO2, uptake[Treatment == 'nonchilled'])		# data for group 1
X2	<- with(CO2, uptake[Treatment == 'chilled'])		# data for group 2
B	<- 40000						# number of posterior samples

burnin	<- B/2
N1		<- length(X1)					# number of subjects in group 1
N2		<- length(X2)					# number of subjects in group 2

mun		<- vector('numeric', length = B)
sigman	<- vector('numeric', length = B)
muc		<- vector('numeric', length = B)
sigmac	<- vector('numeric', length = B)

mun[1]		<- mean(X1)
sigman[1]	<- var(X1)
muc[1]		<- mean(X2)
sigmac[1]	<- var(X2)

for(b in 2:B){
	## update mus ##
	mun[b] <- rnorm(1, mean(X1), sqrt(sigman[b-1]/N1))
	muc[b] <- rnorm(1, mean(X2), sqrt(sigmac[b-1]/N2))

	## update sigmas ##
	sigman[b]	<- rinvgamma(1, N1/2, (1/2)*sum((X1 - mun[b-1])^2))
	sigmac[b]	<- rinvgamma(1, N2/2, (1/2)*sum((X2 - muc[b-1])^2))
}

munb	<- mun[-c(1:burnin)]
mucb	<- muc[-c(1:burnin)]

sigmanb	<- sigman[-c(1:burnin)]
sigmacb	<- sigmac[-c(1:burnin)]

# separate summaries #
quantile(munb, probs = c(0.5, 0.025, 0.975))
quantile(mucb, probs = c(0.5, 0.025, 0.975))

muMat	<- cbind(munb, mucb)
colnames(muMat)	<- c('Non-chilled', 'Chilled')
denplot(muMat, lwd = 3, xlim = c(min(mucb), max(munb)))

# summaries of the difference #
diff	<- munb - mucb

quantile(diff, probs = c(0.5, 0.025, 0.975))
var(diff)
mean(diff > 0)

diffMat	<- matrix(diff, ncol = 1)
colnames(diffMat)	<- 'Difference in Mean CO2 Uptake (non-chilled - chilled)'
denplot(diffMat, lwd = 3)



### Two proportion Equivalent ###

X1			<- 225		# number of successes in group 1
N1			<- 411		# number of trails in group 1
X2			<- 253		# number of successes in group 2
N2			<- 410		# number of trials in group 2

B			<- 20000		# number of posterior samples

theta1	<- rbeta(B, X1 + 1, N1 - X1 + 1)	# we add 1 to X and N - X to fully determine
						#	the posterior distribution
theta2	<- rbeta(B, X2 + 1, N2 - X2 + 1)	# we add 1 to X and N - X to fully determine
						#	the posterior distribution

quantile(theta1, probs = c(0.5, 0.025, 0.975))
quantile(theta2, probs = c(0.5, 0.025, 0.975))

diff	<- theta1 - theta2
quantile(diff, probs = c(0.5, 0.025, 0.975))

mean(diff > 0)

diffMat				<- matrix(diff, ncol = 1)
colnames(diffMat)	<- 'Posterior Distribution of Risk Difference'
denplot(diffMat, lwd = 3)



### Regression Equivalent ###
?longley
X	<- model.matrix(~ GNP + Armed.Forces + Population, data = longley)	# Design matrix, contains X variables
Y	<- with(longley, Employed)						# Outcome, Y
B	<- 40000

burnin	<- B/2
k		<- ncol(X)
n		<- nrow(X)

sig		<- rep(0, B)
beta	<- matrix(0, nrow = B, ncol = k)
bhat	<- c(solve(t(X)%*%X)%*%(t(X)%*%Y))
v		<- solve(t(X)%*%X) 
m		<- v%*%(t(X)%*%Y)

sig[1]		<- (1/(n-k))*sum((Y-X%*%bhat)^2)
beta[1,]	<- bhat

# may take a few moments to run, esp if k is large
for(b in 2:B) { 
	# Sample Beta #
	beta[b,]	<- c(rmvnorm(1, m, sig[b-1]*v)) 
	# beta[b,]	<- c(rmvn(1, m, sig[b-1]*v)) # in package mvnfast, can be slightly faster

	# Sample Sigma^2 #
	sig[b]		<- rinvgamma(1, n/2, sum((Y-X%*%beta[b-1,])^2)/2)
}
colnames(beta)	<- colnames(X)

betab	<- beta[-c(1:burnin),]
sigb	<- sig[-c(1:burnin)]

t(apply(betab, 2, quantile, probs = c(0.5, 0.025, 0.975)))
quantile(sigb, probs = c(0.5, 0.025, 0.975))

denplot(betab, lwd = 3)

sigMat	<- matrix(sigb, ncol = 1)
colnames(sigMat)	<- 'Variance'
denplot(sigMat, lwd = 3)


## Exercise ##

# Return to the CO2 data from the two sample t example. Using uptake as Y and
#	Treatment as X fit a linear regression model. How do your results compare
#	between the two approaches?

