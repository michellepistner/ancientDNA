#####Part 2: Selecting the Best Model for Time#####

gamm.select = 1 ##Baseline signal to noise ratio, will be tuned later

###Prepping the OTU table
otu.closed = miniclo_array(otu.filtered, part = 1)
Y = otu_table(otu.filtered, taxa_are_rows=TRUE)

###Setting priors based on Y
upsilon = ntaxa(Y)+3
Omega = diag(ntaxa(Y))
G = cbind(diag(ntaxa(Y)-1), -1)
Xi = (upsilon-ntaxa(Y))*G%*%Omega%*%t(G)

###Intercept only model
X = t(model.matrix(~1, data = metadata_prep))

###Setting priors based on X
Theta = matrix(0, ntaxa(Y)-1, nrow(X))
Gamma = diag(nrow(X))*gamm.select

posterior = pibble(Y, X, upsilon, Theta, Gamma, Xi, multDirichletBoot=.5)
#posterior = to_clr(posterior)
###Checking the posterior predictive checks
ppc(posterior) 
ppc_summary(posterior)

ppc(posterior, from_scratch=TRUE) 
ppc_summary(posterior, from_scratch=TRUE)

var.exp.Y = var.explained.Y(posterior,otu.closed)

###Percentage of variance explained

print("Percent of Variation Explained by Intercept-Only Model (Counts):")
mean(var.exp.Y)

###Date model
X = t(model.matrix(~BlackDeath_1346_1353 - 1, data = metadata_prep))

###Setting priors based on X
Theta = matrix(0, ntaxa(Y)-1, nrow(X))
Gamma = diag(nrow(X))*gamm.select

posterior = pibble(Y, X, upsilon, Theta, Gamma, Xi, multDirichletBoot = 1)
posterior = to_clr(posterior)
###Checking the posterior predictive checks
ppc(posterior) 
ppc_summary(posterior)

ppc(posterior, from_scratch=TRUE)
ppc_summary(posterior, from_scratch=TRUE)


var.exp.Y = var.explained.Y(posterior,otu.closed)

###Percentage of variance explained

print("Percent of Variation Explained by Date Model (Counts):")
mean(var.exp.Y)

###Cemetery model
X = t(model.matrix(~Cemetry - 1, data = metadata_prep))

###Setting priors based on X
Theta = matrix(0, ntaxa(Y)-1, nrow(X))
Gamma = diag(nrow(X))*gamm.select

posterior = pibble(Y, X, upsilon, Theta, Gamma, Xi, multDirichletBoot = 1)
posterior = to_clr(posterior)

var.exp.Y = var.explained.Y(posterior,otu.closed)

###Percentage of variance explained

print("Percent of Variation Explained by Cemetery Model (Counts):")
mean(var.exp.Y)


#####Part 3: Cross-Validation to Determine Signal to Noise Ratio#####

###Going to maximize the percentage of variance explained
gamm.opt = seq(1,5, by = 1)
avg.var = rep(NA, length(gamm.opt))

for(i in 1:length(gamm.opt)){
  ###Date model
  X = t(model.matrix(~BlackDeath_1346_1353 - 1, data = metadata_prep))
  
  ###Setting priors based on X
  Theta = matrix(0, ntaxa(Y)-1, nrow(X))
  Gamma = diag(nrow(X))*gamm.opt[i]
  
  posterior = pibble(Y, X, upsilon, Theta, Gamma, Xi, multDirichletBoot = 1)
  posterior = to_clr(posterior)
  
  var.exp = var.explained.Y(posterior, otu.closed)
  ###Percentage of variance explained
  avg.var[i] =  mean(var.exp)
  print(i)
}

gamm.select = which.max(avg.var)


####Fitting optimal model
###Date_100 model
X = t(model.matrix(~BlackDeath_1346_1353 - 1, data = metadata_prep))

###Setting priors based on X
Theta = matrix(0, ntaxa(Y)-1, nrow(X))
Gamma = diag(nrow(X))*gamm.select

posterior = pibble(Y, X, upsilon, Theta, Gamma, Xi, multDirichletBoot = 1)
posterior = to_clr(posterior)
###Checking the posterior predictive checks
ppc(posterior) 
ppc_summary(posterior)

ppc(posterior, from_scratch=TRUE)
ppc_summary(posterior, from_scratch=TRUE)


var.exp = var.explained.Y(posterior, otu.closed)
###Percentage of variance explained

print("Percent of Variation Explained by Date_100 Model (Counts):")
mean(var.exp)
quantile(var.exp, c(0.025,0.975))

#####Part 4: Percent of Variability Explained by Other Sources#####


###First, getting the part "explained by" Date_100
Y.pred = predict(posterior, response = "Y", from_scratch = TRUE)

###Variation explained by other sources
Y.samp = apply(Y.pred, MARGIN=c(1,2), FUN = "mean")

otu.closed.Ysamp = miniclo_array(Y.samp, parts = 1)

rownames(Y.samp) = rownames(Y)
names(Y.samp) = names(Y)
Y.samp = otu_table(Y.samp, taxa_are_rows=TRUE)


##Sanity check, should be close to zero
X = t(model.matrix(~1, data = metadata_prep))

Theta = matrix(0, ntaxa(Y.samp)-1, nrow(X))
Gamma = diag(nrow(X))*gamm.select

posterior <- pibble(Y.samp, X, upsilon, Theta, Gamma, Xi, multDirichletBoot = 1)
posterior = to_clr(posterior)


var.exp.int = var.explained.Y(posterior, otu.closed.Ysamp)
###Percentage of variance explained

print("Percent of Variation in Time Explained by Intercept-Only Model (Counts):")

mean(var.exp.int)
quantile(rowMeans(var.exp.int), c(0.025,0.975))

##Another sanity check, should be close to 100%
X = t(model.matrix(~BlackDeath_1346_1353 - 1, data = metadata_prep))

Theta = matrix(0, ntaxa(Y.samp)-1, nrow(X))
Gamma = diag(nrow(X))*gamm.select

posterior <- pibble(Y.samp, X, upsilon, Theta, Gamma, Xi, multDirichletBoot = 1)
posterior = to_clr(posterior)

var.exp.date = var.explained.Y(posterior, otu.closed.Ysamp)
###Percentage of variance explained

print("Percent of Variation in Time Explained by Date_100 Model (Counts):")
mean(var.exp.date)
quantile((var.exp.date), c(0.025,0.975))

X = t(model.matrix(~Cemetry - 1, data = metadata_prep))

Theta = matrix(0, ntaxa(Y.samp)-1, nrow(X))
Gamma = diag(nrow(X))*gamm.select

posterior <- pibble(Y.samp, X, upsilon, Theta, Gamma, Xi, multDirichletBoot = 1)
posterior = to_clr(posterior)

var.exp.cem = var.explained.Y(posterior, otu.closed.Ysamp)

###Percentage of variance explained

print("Percent of Variation in Time Explained by Cemetry Model (Counts):")
mean(var.exp.cem)
quantile((var.exp.cem), c(0.025,0.975))


X = t(model.matrix(~Tooth - 1, data = metadata_prep))

Theta = matrix(0, ntaxa(Y.samp)-1, nrow(X))
Gamma = diag(nrow(X))*gamm.select

posterior <- pibble(Y.samp, X, upsilon, Theta, Gamma, Xi, multDirichletBoot = 1)
posterior = to_clr(posterior)

var.exp.tooth = var.explained.Y(posterior, otu.closed.Ysamp)

###Percentage of variance explained

print("Percent of Variation in Time Explained by Tooth Model (Counts):")
mean(var.exp.tooth)
quantile((var.exp.tooth), c(0.025,0.975))

