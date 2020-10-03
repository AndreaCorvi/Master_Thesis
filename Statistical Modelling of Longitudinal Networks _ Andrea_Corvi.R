##### Complete Analysis #####

#### Packages required ####
library(readr)
library(corrr)
library(dplyr)
library(tidyr)
library(network)
library(networkDynamic)
library(statnet)
library(htmlwidgets)
library(latticeExtra)
library(tergm)
library(btergm)
library(texreg)
library(lattice)
library(PerformanceAnalytics)
library(xts)
library(GGally)

# Loading of data which has already been cleaned from missing values
# "Month" column has been added with excel software based on the "Date" column 
# and the special character "^" from the indexes names has been removed
DataG = read_csv("DataG.csv") 

#### Covariates #### 
# Firstly I create the covariates that will be used in the build of the
# network as nodal attributes

cov_data = DataG[,2:22]
Mydata = DataG %>% dplyr::select(-IXIC,-NYA,-GDAXI,-AEX,-N225,-KS11)

cov_data_list = split(cov_data, as.factor(cov_data$Month))

for (i in 1:188){
  a = cov_data_list[[i]]
  cov_data_list[[i]] = a[,-2]
}

betas = matrix(nrow = 13,ncol = 188 )
betas = as.data.frame(betas)

compute_beta = function(x){
  x = xts(x = x[,-1], order.by = as.Date(x$Date))
  a = CAPM.beta(x$SMSG, x$KS11)
  b = CAPM.beta(x$TOSH, x$N225)
  c = CAPM.beta(x$FUJI, x$N225)
  d = CAPM.beta(x$SONY, x$N225)
  e = CAPM.beta(x$AAPL, x$IXIC)
  f = CAPM.beta(x$AMZN, x$IXIC)
  g = CAPM.beta(x$CSCO, x$IXIC)
  h = CAPM.beta(x$GOOG, x$IXIC)
  i = CAPM.beta(x$MSFT, x$IXIC)
  j = CAPM.beta(x$PHG, x$AEX)
  k = CAPM.beta(x$SAP, x$NYA)
  l = CAPM.beta(x$SIEGY, x$GDAXI)
  m = CAPM.beta(x$STM, x$GDAXI) 
  out = c(a,b,c,d,e,f,g,h,i,j,k,l,m)
}

for (i in 1:188){
  temp = cov_data_list[[i]]
  betas[,i] = compute_beta(temp)
}

rnam = colnames(Mydata)
rnam = rnam[4:16]
rownames(betas) = rnam


# Now I create a dataframe for the sharpe ratios to use it as nodal attributes during the building of the network
# Yield_data = DataG %>% dplyr::select(-IXIC,-NYA,-GDAXI,-AEX,-N225,-KS11) %>% dplyr::select(-1)
# write.csv(Yield_data, "C:/Users/Yield_data.csv")
# Added year column through excel +  Yield 
# Yield values source:
# https://www.macrotrends.net/2016/10-year-treasury-bond-rate-yield-chart
Yield_data <- read_csv("Yield_data.csv") %>% dplyr::select(-1)
yield_list <- split(Yield_data, as.factor(Yield_data$Month))

SR = matrix(nrow = 13,ncol = 188 )
SR = as.data.frame(betas)

compute_sharpe_ratio = function(x){
  x = xts(x = x[,4:17], order.by = as.Date(x$Date))
  a = SharpeRatio(x$SMSG, Rf = x$Yield, p = 0.95, FUN = "StdDev")
  b = SharpeRatio(x$TOSH, Rf = x$Yield, p = 0.95, FUN = "StdDev")
  c = SharpeRatio(x$FUJI, Rf = x$Yield, p = 0.95, FUN = "StdDev")
  d = SharpeRatio(x$SONY, Rf = x$Yield, p = 0.95, FUN = "StdDev")
  e = SharpeRatio(x$AAPL, Rf = x$Yield, p = 0.95, FUN = "StdDev")
  f = SharpeRatio(x$AMZN, Rf = x$Yield, p = 0.95, FUN = "StdDev")
  g = SharpeRatio(x$CSCO, Rf = x$Yield, p = 0.95, FUN = "StdDev")
  h = SharpeRatio(x$GOOG, Rf = x$Yield, p = 0.95, FUN = "StdDev")
  i = SharpeRatio(x$MSFT, Rf = x$Yield, p = 0.95, FUN = "StdDev")
  j = SharpeRatio(x$PHG, Rf = x$Yield, p = 0.95, FUN = "StdDev")
  k = SharpeRatio(x$SAP, Rf = x$Yield, p = 0.95, FUN = "StdDev")
  l = SharpeRatio(x$SIEGY, Rf = x$Yield, p = 0.95, FUN = "StdDev")
  m = SharpeRatio(x$STM, Rf = x$Yield, p = 0.95, FUN = "StdDev")
  out = c(a,b,c,d,e,f,g,h,i,j,k,l,m)
}

for (i in 1:188){
  temp = yield_list[[i]]
  SR[,i] = compute_sharpe_ratio(temp)
}

rnam = colnames(Mydata)
rnam = rnam[4:16]
rownames(SR) = rnam

# Now I set factors for the covariates based on quantiles
library(gtools)


c = matrix(nrow = 13, ncol = 188)
c = as.data.frame(c)

for (i in 1:188){
  a = SR[,i]
  b = quantcut(a, q = 4)
  b = (as.numeric(b))
  c[,i] = b
}


# for the betas
d = matrix(nrow = 13, ncol = 188)
d = as.data.frame(d)

for (i in 1:188){
  a = betas[,i]
  b = quantcut(a, q = 4)
  b = (as.numeric(b))
  d[,i] = b
}


#### Correlation matrix filtering and network building ####

# Means that there are 188 time points
month_names = as.vector(unique((DataG[,3]))) 

# All subdataframe list which will be my list of networks 
# each rappresented by an adjacency matrix

# I remove the indexes
Mydata = DataG %>% dplyr::select(-IXIC,-NYA,-GDAXI,-AEX,-N225,-KS11)
net_list1 <- split(Mydata, as.factor(Mydata$Month))

# User define function that takes a dataframe as an input a dataframe
# Takes only the columns with the returns and correlates them one vs one
# then the correlation is filtered with respect to a treshold. If it is upper then the cell is set as 1
# meaning that there is a relation, otherwise it is set = 0
extract_correlation1 = function(x){
  x = x[,4:16]
  x = mutate_all(x, function(x) as.numeric(as.character(x)))
  temp = correlate(x,diagonal = 0)  
  temp1 = temp[,-1]
  temp1[temp1 < 0.6] <- 0
  temp1[temp1 >= 0.6] <- 1
  out = temp1
  # out = cbind(temp$rowname, temp1) adds the rowname columns -> with it I don't have a square matrix
  return(out)
}

# Now the function is applied to each subdataframe (one per month) and eache dataframe is firstly
# converted as a matrix and then as a network object. In the end we got a list of adjacency matrices
# the network is directed
for (i in 1:188){
  net_list1[[i]] = extract_correlation1(net_list1[[i]])
  net_list1[[i]] = as.matrix(net_list1[[i]])
  net_list1[[i]] = network(net_list1[[i]],matrix.type="adjacency",directed=FALSE)
  network::set.vertex.attribute(net_list1[[i]], "Beta", as.vector(betas[,i])) # Add the betas as nodal
  # attributes
  network::set.vertex.attribute(net_list1[[i]],"SR", as.vector(SR[,i]))
  network::set.vertex.attribute(net_list1[[i]],"SR_factorized", as.vector(c[,i]))
  network::set.vertex.attribute(net_list1[[i]],"Beta_factorized", as.vector(d[,i]))
  #nam <- paste("CorMat", i, sep = ".")   useful to change the name of the single object but a
  #assign(nam, a)                         different loop has to be built in order to exploit it
}

# Create now second network with correlation filtering = 0.7 #
net_list2 <- split(DataG, as.factor(DataG$Month))

extract_correlation2 = function(x){
  x = x[,4:16]
  x = mutate_all(x, function(x) as.numeric(as.character(x)))
  temp = correlate(x,diagonal = 1)
  temp1 = temp[,-1]
  temp1[temp1 < 0.7] <- 0
  temp1[temp1 >= 0.7] <- 1
  out = temp1
  # out = cbind(temp$rowname, temp1) adds the rowname columns -> with it I don't have a square matrix
  return(out)
}

for (i in 1:188){
  net_list2[[i]] = extract_correlation2(net_list2[[i]])
  net_list2[[i]] = as.matrix(net_list2[[i]])
  net_list2[[i]] = network(net_list2[[i]],matrix.type="adjacency",directed=FALSE)
  network::set.vertex.attribute(net_list2[[i]], "Beta", as.vector(betas[,i]))
  network::set.vertex.attribute(net_list2[[i]],"SR", as.vector(SR[,i]))
  network::set.vertex.attribute(net_list2[[i]],"SR_factorized", as.vector(c[,i]))
  network::set.vertex.attribute(net_list2[[i]],"Beta_factorized", as.vector(d[,i]))
}

# Create third network, correlation filtering 0.65 #
net_list3 <- split(DataG, as.factor(DataG$Month))

extract_correlation3 = function(x){
  x = x[,4:16]
  x = mutate_all(x, function(x) as.numeric(as.character(x)))
  temp = correlate(x,diagonal = 1)
  temp1 = temp[,-1]
  temp1[temp1 < 0.65] <- 0
  temp1[temp1 >= 0.65] <- 1
  out = temp1
  # out = cbind(temp$rowname, temp1) adds the rowname columns -> with it I don't have a square matrix
  return(out)
}

for (i in 1:188){
  net_list3[[i]] = extract_correlation3(net_list3[[i]])
  net_list3[[i]] = as.matrix(net_list3[[i]])
  net_list3[[i]] = network(net_list3[[i]],matrix.type="adjacency",directed=FALSE)
  network::set.vertex.attribute(net_list3[[i]], "Beta", as.vector(betas[,i]))
  network::set.vertex.attribute(net_list3[[i]],"SR", as.vector(SR[,i]))
  network::set.vertex.attribute(net_list3[[i]],"SR_factorized", as.vector(c[,i]))
  network::set.vertex.attribute(net_list3[[i]],"Beta_factorized", as.vector(d[,i]))
}

#Create fourth network, correlation filtered at .75 #
net_list4 <- split(DataG, as.factor(DataG$Month))

extract_correlation4 = function(x){
  x = x[,4:16]
  x = mutate_all(x, function(x) as.numeric(as.character(x)))
  temp = correlate(x,diagonal = 1)
  temp1 = temp[,-1]
  temp1[temp1 < 0.75] <- 0
  temp1[temp1 >= 0.75] <- 1
  out = temp1
  # out = cbind(temp$rowname, temp1) adds the rowname columns -> with it I don't have a square matrix
  return(out)
}

for (i in 1:188){
  net_list4[[i]] = extract_correlation4(net_list4[[i]])
  net_list4[[i]] = as.matrix(net_list4[[i]])
  net_list4[[i]] = network(net_list4[[i]],matrix.type="adjacency",directed=FALSE)
  network::set.vertex.attribute(net_list4[[i]], "Beta", as.vector(betas[,i]))
  network::set.vertex.attribute(net_list4[[i]],"SR", as.vector(SR[,i]))
  network::set.vertex.attribute(net_list4[[i]],"SR_factorized", as.vector(c[,i]))
  network::set.vertex.attribute(net_list4[[i]],"Beta_factorized", as.vector(d[,i]))
}

# Four different network have been built. That will allow to study if a different correlation
# threshold has an impact on the modeling and goodness of fit of the model


#### Exploratory and Graphical Analysis  ####

## I create a networkDynamic object to perform some visual analysis (plot will be added)

tech_net1 = networkDynamic(network.list = net_list)
tech_net2 = networkDynamic(network.list = net_list2)
# tech_net3 = networkDynamic(network.list = net_list3)
# tech_net4 = networkDynamic(network.list = net_list4)


# Plotting the evolution of the second network
render.d3movie(tech_net2,
plot.par=list(displaylabels=T),
output.mode = 'htmlWidget') # using htmlwidgets package here

# Networks density
n1_dens = list()
n2_dens = list()
n3_dens = list()
n4_dens = list()

for (i in 1:188){
  n1_dens[[i]] = network.density(net_list1[[i]])
  n2_dens[[i]] = network.density(net_list2[[i]])
  n3_dens[[i]] = network.density(net_list3[[i]])
  n4_dens[[i]] = network.density(net_list4[[i]])
}

mean(unlist(n1_dens))
mean(unlist(n2_dens))
mean(unlist(n3_dens))
mean(unlist(n4_dens))

EDA_net1 = net_list4[[176]] 
EDA_net2 = net_list4[[188]]

# network plot

ggnet2(EDA_net1, color = "steelblue1", edge.color = "darkorange",
        edge.size = 0.8, alpha = 0.6, size = 20,label.alpha = 0.75, edge.alpha = 0.4, label = TRUE, label.size = 3)

ggnet2(EDA_net2, color = "steelblue1", edge.color = "darkorange",
        edge.size = 0.8, alpha = 0.6, size = 20,label.alpha = 0.75, edge.alpha = 0.4, label = TRUE, label.size = 3)

#### Analysis with TERGMs ####
### Prediction for last month given last year data.  ###

set.seed(123456)

# THe model is fitted on 11 months data.
net_fit = btergm(net_list4[176:187] ~edges+nodecov("SR")+
                   absdiff("SR")+nodecov("Beta")+nodecov("Beta_factorized")+
                   nodecov("SR_factorized")+absdiff("Beta")+
                   absdiffcat("SR_factorized")+absdiffcat("Beta_factorized")+
                   timecov(transform = function(t)t)+balance+
                   gwesp(decay=0, fixed=TRUE, cutoff=30)
                 +gwdsp(decay=0, fixed=TRUE, cutoff=30)
                 +gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                   memory(type = "autoregression"),
                 R = 200, parallel = "multicore", ncpus = 6, 
                 verbose = TRUE)


# The goodness of fit is computed using 200 simulation with prediction on the twelfth month.
gof_net = gof(net_fit,nsim = 200, target = net_list4[[188]],
              statistics =c(rocpr,esp, dsp,deg,kstar, geodesic, triad.undirected
              ),
              coef = coef(net_fit) ,seed = 12, formula = net_list4[187:188] ~ edges+
                nodecov("SR")+absdiff("SR")+nodecov("Beta")+
                nodecov("Beta_factorized")+nodecov("SR_factorized")+
                absdiff("Beta")+absdiffcat("SR_factorized")+
                absdiffcat("Beta_factorized")+timecov(transform =  function(t)t)+balance
              +gwesp(decay=0, fixed=TRUE, cutoff=30)
              +gwdsp(decay=0, fixed=TRUE, cutoff=30)+
                gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                memory(type = "autoregression"))

plot(gof_net[[7]], main = "Network 1", col = "dodgerblue", lcol = "orange", outline = TRUE, median.col = "orange")
Net_auc = gof_net[[1]]$auc.roc
Net_auc

# The model fitting has been then repeated for each network, alongside the plot which has to be made
# also for each statistics on which the gof has to evaluated.


# The results need to be validated even through the rolling origin methodology.

#### Rolling Origin ####
# In the following loops, for each network the model are fitted through the rolling origin
# methodoloy  using 1,2 and 3 years data to predicting the network state at one specific time point.
### Network 1 ###

### Annual Basis
Net1_Roc_auc1 = list()
Net1_Pr_auc1 = list()

i = 1
while (i < 178){
  set.seed(123456)
  a = i
  b = i+10
  c = i+11
  net1_fit = btergm(net_list1[a:b] ~edges+nodecov("SR")+
                      absdiff("SR")+nodecov("Beta")+nodecov("Beta_factorized")+
                      nodecov("SR_factorized")+absdiff("Beta")+
                      absdiffcat("SR_factorized")+absdiffcat("Beta_factorized")+
                      timecov(transform = function(t)t)+balance+
                      gwesp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdsp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                      memory(type = "stability"),
                    R = 100, parallel = "multicore", ncpus = 6, 
                    verbose = TRUE)
  
  gof_net1_fit = gof(net1_fit,nsim = 100, target = net_list1[[c]],
                     statistics =c(rocpr),
                     coef = coef(net1_fit) ,seed = 12, formula = net_list1[b:c] ~ edges+
                       nodecov("SR")+absdiff("SR")+nodecov("Beta")+
                       nodecov("Beta_factorized")+nodecov("SR_factorized")+
                       absdiff("Beta")+absdiffcat("SR_factorized")+
                       absdiffcat("Beta_factorized")+timecov(transform =  function(t)t)+balance
                     +gwesp(decay=0, fixed=TRUE, cutoff=30)
                     +gwdsp(decay=0, fixed=TRUE, cutoff=30)+
                       gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                       memory(type = "stability"))
  
  Net1_Roc_auc1[[i]] =gof_net1_fit[[1]]$auc.roc
  Net1_Pr_auc1[[i]] = gof_net1_fit[[1]]$auc.pr
  i = i+1
}



### Two years basis
Net1_Roc_auc2 = list()
Net1_Pr_auc2 = list()

i = 1
while (i < 166){
  set.seed(123456)
  a = i
  b = i+22
  c = i+23
  net1_fit = btergm(net_list1[a:b] ~edges+nodecov("SR")+
                      absdiff("SR")+nodecov("Beta")+nodecov("Beta_factorized")+
                      nodecov("SR_factorized")+absdiff("Beta")+
                      absdiffcat("SR_factorized")+absdiffcat("Beta_factorized")+
                      timecov(transform = function(t)t)+balance+
                      gwesp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdsp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                      memory(type = "stability"),
                    R = 100, parallel = "multicore", ncpus = 6, 
                    verbose = TRUE)
  
  gof_net1_fit = gof(net1_fit,nsim = 100, target = net_list1[[c]],
                     statistics =c(rocpr),
                     coef = coef(net1_fit) ,seed = 12, formula = net_list1[b:c] ~ edges+
                       nodecov("SR")+absdiff("SR")+nodecov("Beta")+
                       nodecov("Beta_factorized")+nodecov("SR_factorized")+
                       absdiff("Beta")+absdiffcat("SR_factorized")+
                       absdiffcat("Beta_factorized")+timecov(transform =  function(t)t)+balance
                     +gwesp(decay=0, fixed=TRUE, cutoff=30)
                     +gwdsp(decay=0, fixed=TRUE, cutoff=30)+
                       gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                       memory(type = "stability"))
  
  Net1_Roc_auc2[[i]] =gof_net1_fit[[1]]$auc.roc
  Net1_Pr_auc2[[i]] = gof_net1_fit[[1]]$auc.pr
  i = i+1
}



### Three years basis
Net1_Roc_auc3 = list()
Net1_Pr_auc3 = list()

i = 1
while (i < 154){
  set.seed(123456)
  a = i
  b = i+34
  c = i+35
  net1_fit = btergm(net_list1[a:b] ~edges+nodecov("SR")+
                      absdiff("SR")+nodecov("Beta")+nodecov("Beta_factorized")+
                      nodecov("SR_factorized")+absdiff("Beta")+
                      absdiffcat("SR_factorized")+absdiffcat("Beta_factorized")+
                      timecov(transform = function(t)t)+balance+
                      gwesp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdsp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                      memory(type = "stability"),
                    R = 100, parallel = "multicore", ncpus = 6, 
                    verbose = TRUE)
  
  gof_net1_fit = gof(net1_fit,nsim = 100, target = net_list1[[c]],
                     statistics =c(rocpr),
                     coef = coef(net1_fit) ,seed = 12, formula = net_list1[b:c] ~ edges+
                       nodecov("SR")+absdiff("SR")+nodecov("Beta")+
                       nodecov("Beta_factorized")+nodecov("SR_factorized")+
                       absdiff("Beta")+absdiffcat("SR_factorized")+
                       absdiffcat("Beta_factorized")+timecov(transform =  function(t)t)+balance
                     +gwesp(decay=0, fixed=TRUE, cutoff=30)
                     +gwdsp(decay=0, fixed=TRUE, cutoff=30)+
                       gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                       memory(type = "stability"))
  
  Net1_Roc_auc3[[i]] =gof_net1_fit[[1]]$auc.roc
  Net1_Pr_auc3[[i]] = gof_net1_fit[[1]]$auc.pr
  i = i+1
}



### Network 2 ###
Net2_Roc_auc1 = list()
Net2_Pr_auc1 = list()

i = 1
while (i < 178){
  set.seed(123456)
  a = i
  b = i+10
  c = i+11
  net2_fit = btergm(net_list2[a:b] ~edges+nodecov("SR")+
                      absdiff("SR")+nodecov("Beta")+nodecov("Beta_factorized")+
                      nodecov("SR_factorized")+absdiff("Beta")+
                      absdiffcat("SR_factorized")+absdiffcat("Beta_factorized")+
                      timecov(transform = sqrt)+balance+
                      gwesp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdsp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                      memory(type = "stability"),
                    R = 100, parallel = "multicore", ncpus = 6, 
                    verbose = TRUE)
  
  gof_net2_fit = gof(net2_fit,nsim = 100, target = net_list2[[c]],
                     statistics =c(rocpr),
                     coef = coef(net2_fit) ,seed = 12, formula = net_list2[b:c] ~ edges+
                       nodecov("SR")+absdiff("SR")+nodecov("Beta")+
                       nodecov("Beta_factorized")+nodecov("SR_factorized")+
                       absdiff("Beta")+absdiffcat("SR_factorized")+
                       absdiffcat("Beta_factorized")+timecov(transform =  sqrt)+balance
                     +gwesp(decay=0, fixed=TRUE, cutoff=30)
                     +gwdsp(decay=0, fixed=TRUE, cutoff=30)+
                       gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                       memory(type = "stability"))
  
  Net2_Roc_auc1[[i]] =gof_net2_fit[[1]]$auc.roc
  Net2_Pr_auc1[[i]] = gof_net2_fit[[1]]$auc.pr
  i = i+1
}


### Two years basis
Net2_Roc_auc2 = list()
Net2_Pr_auc2 = list()

i = 1
while (i < 166){
  set.seed(123456)
  a = i
  b = i+22
  c = i+23
  net2_fit = btergm(net_list2[a:b] ~edges+nodecov("SR")+
                      absdiff("SR")+nodecov("Beta")+nodecov("Beta_factorized")+
                      nodecov("SR_factorized")+absdiff("Beta")+
                      absdiffcat("SR_factorized")+absdiffcat("Beta_factorized")+
                      timecov(transform = function(t)t)+balance+
                      gwesp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdsp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                      memory(type = "stability"),
                    R = 100, parallel = "multicore", ncpus = 6, 
                    verbose = TRUE)
  
  gof_net2_fit = gof(net2_fit,nsim = 100, target = net_list2[[c]],
                     statistics =c(rocpr),
                     coef = coef(net2_fit) ,seed = 12, formula = net_list2[b:c] ~ edges+
                       nodecov("SR")+absdiff("SR")+nodecov("Beta")+
                       nodecov("Beta_factorized")+nodecov("SR_factorized")+
                       absdiff("Beta")+absdiffcat("SR_factorized")+
                       absdiffcat("Beta_factorized")+timecov(transform =  function(t)t)+balance
                     +gwesp(decay=0, fixed=TRUE, cutoff=30)
                     +gwdsp(decay=0, fixed=TRUE, cutoff=30)+
                       gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                       memory(type = "stability"))
  
  Net2_Roc_auc2[[i]] =gof_net2_fit[[1]]$auc.roc
  Net2_Pr_auc2[[i]] = gof_net2_fit[[1]]$auc.pr
  i = i+1
}


### Three years basis
Net2_Roc_auc3 = list()
Net2_Pr_auc3 = list()

i = 1
while (i < 154){
  set.seed(123456)
  a = i
  b = i+34
  c = i+35
  net2_fit = btergm(net_list2[a:b] ~edges+nodecov("SR")+
                      absdiff("SR")+nodecov("Beta")+nodecov("Beta_factorized")+
                      nodecov("SR_factorized")+absdiff("Beta")+
                      absdiffcat("SR_factorized")+absdiffcat("Beta_factorized")+
                      timecov(transform = function(t)t)+balance+
                      gwesp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdsp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                      memory(type = "stability"),
                    R = 100, parallel = "multicore", ncpus = 6, 
                    verbose = TRUE)
  
  gof_net2_fit = gof(net2_fit,nsim = 100, target = net_list2[[c]],
                     statistics =c(rocpr),
                     coef = coef(net2_fit) ,seed = 12, formula = net_list2[b:c] ~ edges+
                       nodecov("SR")+absdiff("SR")+nodecov("Beta")+
                       nodecov("Beta_factorized")+nodecov("SR_factorized")+
                       absdiff("Beta")+absdiffcat("SR_factorized")+
                       absdiffcat("Beta_factorized")+timecov(transform =  function(t)t)+balance
                     +gwesp(decay=0, fixed=TRUE, cutoff=30)
                     +gwdsp(decay=0, fixed=TRUE, cutoff=30)+
                       gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                       memory(type = "stability"))
  
  Net2_Roc_auc3[[i]] =gof_net2_fit[[1]]$auc.roc
  Net2_Pr_auc3[[i]] = gof_net2_fit[[1]]$auc.pr
  i = i+1
}

### Network 3 ###

### Annual basis
Net3_Roc_auc1 = list()
Net3_Pr_auc1 = list()

i = 1
while (i < 178){
  set.seed(123456)
  a = i
  b = i+10
  c = i+11
  net3_fit = btergm(net_list3[a:b] ~edges+nodecov("SR")+
                      absdiff("SR")+nodecov("Beta")+nodecov("Beta_factorized")+
                      nodecov("SR_factorized")+absdiff("Beta")+
                      absdiffcat("SR_factorized")+absdiffcat("Beta_factorized")+
                      timecov(transform = function(t)t)+balance+
                      gwesp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdsp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                      memory(type = "stability"),
                    R = 100, parallel = "multicore", ncpus = 6, 
                    verbose = TRUE)
  
  gof_net3_fit = gof(net3_fit,nsim = 100, target = net_list3[[c]],
                     statistics =c(rocpr),
                     coef = coef(net3_fit) ,seed = 12, formula = net_list3[b:c] ~ edges+
                       nodecov("SR")+absdiff("SR")+nodecov("Beta")+
                       nodecov("Beta_factorized")+nodecov("SR_factorized")+
                       absdiff("Beta")+absdiffcat("SR_factorized")+
                       absdiffcat("Beta_factorized")+timecov(transform =  function(t)t)+balance
                     +gwesp(decay=0, fixed=TRUE, cutoff=30)
                     +gwdsp(decay=0, fixed=TRUE, cutoff=30)+
                       gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                       memory(type = "stability"))
  
  Net3_Roc_auc1[[i]] =gof_net3_fit[[1]]$auc.roc
  Net3_Pr_auc1[[i]] = gof_net3_fit[[1]]$auc.pr
  i = i+1
}



### Two years basis
Net3_Roc_auc2 = list()
Net3_Pr_auc2 = list()

i = 1
while (i < 166){
  set.seed(123456)
  a = i
  b = i+22
  c = i+23
  net3_fit = btergm(net_list3[a:b] ~edges+nodecov("SR")+
                      absdiff("SR")+nodecov("Beta")+nodecov("Beta_factorized")+
                      nodecov("SR_factorized")+absdiff("Beta")+
                      absdiffcat("SR_factorized")+absdiffcat("Beta_factorized")+
                      timecov(transform = function(t)t)+balance+
                      gwesp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdsp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                      memory(type = "stability"),
                    R = 100, parallel = "multicore", ncpus = 6, 
                    verbose = TRUE)
  
  gof_net3_fit = gof(net1_fit,nsim = 100, target = net_list3[[c]],
                     statistics =c(rocpr),
                     coef = coef(net3_fit) ,seed = 12, formula = net_list3[b:c] ~ edges+
                       nodecov("SR")+absdiff("SR")+nodecov("Beta")+
                       nodecov("Beta_factorized")+nodecov("SR_factorized")+
                       absdiff("Beta")+absdiffcat("SR_factorized")+
                       absdiffcat("Beta_factorized")+timecov(transform =  function(t)t)+balance
                     +gwesp(decay=0, fixed=TRUE, cutoff=30)
                     +gwdsp(decay=0, fixed=TRUE, cutoff=30)+
                       gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                       memory(type = "stability"))
  
  Net3_Roc_auc2[[i]] =gof_net3_fit[[1]]$auc.roc
  Net3_Pr_auc2[[i]] = gof_net3_fit[[1]]$auc.pr
  i = i+1
}



### Three years basis
Net3_Roc_auc3 = list()
Net3_Pr_auc3 = list()

i = 1
while (i < 154){
  set.seed(123456)
  a = i
  b = i+34
  c = i+35
  net3_fit = btergm(net_list3[a:b] ~edges+nodecov("SR")+
                      absdiff("SR")+nodecov("Beta")+nodecov("Beta_factorized")+
                      nodecov("SR_factorized")+absdiff("Beta")+
                      absdiffcat("SR_factorized")+absdiffcat("Beta_factorized")+
                      timecov(transform = function(t)t)+balance+
                      gwesp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdsp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                      memory(type = "stability"),
                    R = 100, parallel = "multicore", ncpus = 6, 
                    verbose = TRUE)
  
  gof_net3_fit = gof(net3_fit,nsim = 100, target = net_list3[[c]],
                     statistics =c(rocpr),
                     coef = coef(net3_fit) ,seed = 12, formula = net_list3[b:c] ~ edges+
                       nodecov("SR")+absdiff("SR")+nodecov("Beta")+
                       nodecov("Beta_factorized")+nodecov("SR_factorized")+
                       absdiff("Beta")+absdiffcat("SR_factorized")+
                       absdiffcat("Beta_factorized")+timecov(transform =  function(t)t)+balance
                     +gwesp(decay=0, fixed=TRUE, cutoff=30)
                     +gwdsp(decay=0, fixed=TRUE, cutoff=30)+
                       gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                       memory(type = "stability"))
  
  Net3_Roc_auc3[[i]] =gof_net3_fit[[1]]$auc.roc
  Net3_Pr_auc3[[i]] = gof_net3_fit[[1]]$auc.pr
  i = i+1
}



### Network 4 ###

### Annual basis
Net4_Roc_auc1 = list()
Net4_Pr_auc1 = list()

i = 1
while (i < 178){
  set.seed(123456)
  a = i
  b = i+10
  c = i+11
  net4_fit = btergm(net_list4[a:b] ~edges+nodecov("SR")+
                      absdiff("SR")+nodecov("Beta")+nodecov("Beta_factorized")+
                      nodecov("SR_factorized")+absdiff("Beta")+
                      absdiffcat("SR_factorized")+absdiffcat("Beta_factorized")+
                      timecov(transform = function(t)t)+balance+
                      gwesp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdsp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                      memory(type = "stability"),
                    R = 100, parallel = "multicore", ncpus = 6, 
                    verbose = TRUE)
  
  gof_net4_fit = gof(net4_fit,nsim = 100, target = net_list4[[c]],
                     statistics =c(rocpr),
                     coef = coef(net4_fit) ,seed = 12, formula = net_list4[b:c] ~ edges+
                       nodecov("SR")+absdiff("SR")+nodecov("Beta")+
                       nodecov("Beta_factorized")+nodecov("SR_factorized")+
                       absdiff("Beta")+absdiffcat("SR_factorized")+
                       absdiffcat("Beta_factorized")+timecov(transform =  function(t)t)+balance
                     +gwesp(decay=0, fixed=TRUE, cutoff=30)
                     +gwdsp(decay=0, fixed=TRUE, cutoff=30)+
                       gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                       memory(type = "stability"))
  
  Net4_Roc_auc1[[i]] =gof_net4_fit[[1]]$auc.roc
  Net4_Pr_auc1[[i]] = gof_net4_fit[[1]]$auc.pr
  i = i+1
}



### Two years basis
Net4_Roc_auc2 = list()
Net4_Pr_auc2 = list()

i = 1
while (i < 166){
  set.seed(123456)
  a = i
  b = i+22
  c = i+23
  net4_fit = btergm(net_list4[a:b] ~edges+nodecov("SR")+
                      absdiff("SR")+nodecov("Beta")+nodecov("Beta_factorized")+
                      nodecov("SR_factorized")+absdiff("Beta")+
                      absdiffcat("SR_factorized")+absdiffcat("Beta_factorized")+
                      timecov(transform = function(t)t)+balance+
                      gwesp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdsp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                      memory(type = "stability"),
                    R = 100, parallel = "multicore", ncpus = 6, 
                    verbose = TRUE)
  
  gof_net4_fit = gof(net4_fit,nsim = 100, target = net_list4[[c]],
                     statistics =c(rocpr),
                     coef = coef(net4_fit) ,seed = 12, formula = net_list4[b:c] ~ edges+
                       nodecov("SR")+absdiff("SR")+nodecov("Beta")+
                       nodecov("Beta_factorized")+nodecov("SR_factorized")+
                       absdiff("Beta")+absdiffcat("SR_factorized")+
                       absdiffcat("Beta_factorized")+timecov(transform =  function(t)t)+balance
                     +gwesp(decay=0, fixed=TRUE, cutoff=30)
                     +gwdsp(decay=0, fixed=TRUE, cutoff=30)+
                       gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                       memory(type = "stability"))
  
  Net4_Roc_auc2[[i]] =gof_net4_fit[[1]]$auc.roc
  Net4_Pr_auc2[[i]] = gof_net4_fit[[1]]$auc.pr
  i = i+1
}



### Three years basis
Net4_Roc_auc3 = list()
Net4_Pr_auc3 = list()

i = 1
while (i < 154){
  set.seed(123456)
  a = i
  b = i+34
  c = i+35
  net4_fit = btergm(net_list4[a:b] ~edges+nodecov("SR")+
                      absdiff("SR")+nodecov("Beta")+nodecov("Beta_factorized")+
                      nodecov("SR_factorized")+absdiff("Beta")+
                      absdiffcat("SR_factorized")+absdiffcat("Beta_factorized")+
                      timecov(transform = function(t)t)+balance+
                      gwesp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdsp(decay=0, fixed=TRUE, cutoff=30)
                    +gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                      memory(type = "stability"),
                    R = 100, parallel = "multicore", ncpus = 6, 
                    verbose = TRUE)
  
  gof_net4_fit = gof(net4_fit,nsim = 100, target = net_list4[[c]],
                     statistics =c(rocpr),
                     coef = coef(net4_fit) ,seed = 12, formula = net_list4[b:c] ~ edges+
                       nodecov("SR")+absdiff("SR")+nodecov("Beta")+
                       nodecov("Beta_factorized")+nodecov("SR_factorized")+
                       absdiff("Beta")+absdiffcat("SR_factorized")+
                       absdiffcat("Beta_factorized")+timecov(transform =  function(t)t)+balance
                     +gwesp(decay=0, fixed=TRUE, cutoff=30)
                     +gwdsp(decay=0, fixed=TRUE, cutoff=30)+
                       gwdegree(decay=0, fixed=TRUE, cutoff=30, levels=NULL)+
                       memory(type = "stability"))
  
  Net4_Roc_auc3[[i]] =gof_net4_fit[[1]]$auc.roc
  Net4_Pr_auc3[[i]] = gof_net4_fit[[1]]$auc.pr
  i = i+1
}



mean(unlist(Net1_Roc_auc1))
mean(unlist(Net1_Roc_auc2))
mean(unlist(Net1_Roc_auc3))

mean(unlist(Net2_Roc_auc1))
mean(unlist(Net2_Roc_auc2))
mean(unlist(Net2_Roc_auc3))

mean(unlist(Net3_Roc_auc1))
mean(unlist(Net3_Roc_auc2))
mean(unlist(Net3_Roc_auc3))

mean(unlist(Net4_Roc_auc1))
mean(unlist(Net4_Roc_auc2))
hist(unlist(Net4_Roc_auc3))

### Histogram of the distributions of the AUC for the Roc curve
# of the fourth network based on three years data
library(ggplot2)
a = as.data.frame(cbind(xaxis, unlist(Net4_Roc_auc3)))
xaxis = seq(1,153,1)
a = as.data.frame(cbind(xaxis, unlist(Net4_Roc_auc3)))

ggplot(a, aes(x= V2)) + geom_histogram(bins = 25, alpha = 0.6, color = "dodgerblue", fill = "dodgerblue")+
  geom_vline(aes(xintercept=mean(V2)), size = 2, color = "orange", alpha = 0.7)+
  labs(x="Roc AUC Network 4")+
  xlim(c(0,1))
