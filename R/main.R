library(kernlab)
library(drf)
library(Matrix)
library(dplyr)
library(doParallel)
library(doRNG)
library(parallel)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(drf)
library(ggplot2)
library(fastDummies)
library(Hmisc)
library(gridExtra)
library(designmatch)

source("compute_causaldrf_vimp.R")

set.seed(10)


#run <- c("toy_examples","simulation_illustration" , "wage", "lalonde", "pension")
run<-"toy_examples"








if ("toy_examples" %in% run){

######################
#####Toy Examples #######
######################
n<-500
p<-5


Wtrain<- matrix(sample(c(0,1), size=n, prob=c(0.5,0.5), replace=T)) #matrix(c(rep(0,n/2), rep(1,n/2)))

Xtrain<-matrix(  runif(n=n*p) +2, nrow=n, ncol=p    )


Ytrain0 <- matrix(rnorm(n=sum(Wtrain==0))) #matrix(rnorm(n=n/2))

Ytrain1<-matrix(rnorm(n=sum(Wtrain==1))) #matrix(rnorm(n=n/2))
Ytrain2<-matrix(rnorm(n=sum(Wtrain==1), mean=Xtrain[,1])) #matrix(rnorm(n=n/2, mean=Xtrain[,1]))
Ytrain3<-matrix(rnorm(n=sum(Wtrain==1), sd=sqrt(1/Xtrain[,1]^2))) #matrix(rnorm(n=n/2, sd=sqrt(1/Xtrain[,1]^2)))
Ytrain4<-matrix(rnorm(n=sum(Wtrain==1), mean=Xtrain, sd=sqrt(Xtrain[,1]^2)))#matrix(rnorm(n=n/2, mean=Xtrain, sd=sqrt(Xtrain[,1]^2)))

Wtrain<-matrix(sort(Wtrain))


fit1 <- drf(X=Xtrain, Y=rbind(Ytrain0,Ytrain1), W=Wtrain, num.trees = 3000, ci.group.size = 100)
fit2 <- drf(X=Xtrain, Y=rbind(Ytrain0,Ytrain2), W=Wtrain, num.trees = 3000, ci.group.size = 100)
fit3 <- drf(X=Xtrain, Y=rbind(Ytrain0,Ytrain3), W=Wtrain, num.trees = 3000, ci.group.size = 100)
fit4 <- drf(X=Xtrain, Y=rbind(Ytrain0,Ytrain4), W=Wtrain, num.trees = 3000, ci.group.size = 100)

x<-matrix(rep(2.5,p), nrow=1, ncol=p)


witness1 <- predict_witness(fit1, alpha = 0.05, newdata = x, newtreatment = matrix(1) )
witness2 <- predict_witness(fit2, alpha = 0.05, newdata = x, newtreatment = matrix(1))
witness3 <- predict_witness(fit3, alpha = 0.05, newdata = x, newtreatment = matrix(1))
witness4 <- predict_witness(fit4, alpha = 0.05, newdata = x, newtreatment = matrix(1))

par(mfrow=c(2,4))


### Plot 1: No effect ######

# Generate data points for the densities
grid <- seq(min(c(Ytrain0,Ytrain1))-1, max(c(Ytrain0,Ytrain1))+1, length.out = 100)
density1 <- dnorm(grid, mean = 0, sd = 1)
density2 <- dnorm(grid, mean = 0, sd = 1)

# Create a data frame for ggplot
data1 <- data.frame(x = grid, density = density1, group = "Mean 0")
data2 <- data.frame(x = grid, density = density2, group = "Mean 1.5")

# Combine the data frames
data <- rbind(data1, data2)

# Plot the Gaussian densities
p1<-ggplot() +
  geom_line(data = data1, aes(x = x, y = density), color = "blue", size = 1.5) +
  geom_line(data = data2, aes(x = x, y = density), color = "red", size = 1.5) +
  labs(title = "", x = "X", y = "Density") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text = element_text(size = 12)    # Increase axis text size
  )




p2<-tibble(
  Y     = c(Ytrain0,Ytrain1),
  w     = witness1[1, ],
  lower = witness1[2, ],
  upper = witness1[3, ]
) %>%
  ggplot(aes(x = Y, y = w)) +
  geom_line() +
  geom_hline(yintercept = 0, size = 1.5)+
  geom_line(aes(y = lower), lty = 2) +
  geom_line(aes(y = upper), lty = 2)

grid.arrange(p1, p2, ncol = 2)

### Plot 1: No effect ######


### Plot 2: Mean effect ######

# Generate data points for the densities
grid <- seq(min(c(Ytrain0,Ytrain2))-1, max(c(Ytrain0,Ytrain2))+1, length.out = 100)
density1 <- dnorm(grid, mean = 0, sd = 1)
density2 <- dnorm(grid, mean = x, sd = 1)

# Create a data frame for ggplot
data1 <- data.frame(x = grid, density = density1, group = "Mean 0")
data2 <- data.frame(x = grid, density = density2, group = "Mean 1.5")

# Combine the data frames
data <- rbind(data1, data2)

# Plot the Gaussian densities
p1<-ggplot() +
  geom_line(data = data1, aes(x = x, y = density), color = "blue", size = 1.5) +
  geom_line(data = data2, aes(x = x, y = density), color = "red", size = 1.5) +
  labs(title = "", x = "X", y = "Density") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text = element_text(size = 12)    # Increase axis text size
  )




p2<-tibble(
  Y     = c(Ytrain0,Ytrain2),
  w     = witness2[1, ],
  lower = witness2[2, ],
  upper = witness2[3, ]
) %>%
  ggplot(aes(x = Y, y = w)) +
  geom_line() +
  geom_hline(yintercept = 0, size = 1.5)+
  geom_line(aes(y = lower), lty = 2) +
  geom_line(aes(y = upper), lty = 2)

grid.arrange(p1, p2, ncol = 2)

### Plot 2: Mean effect ######



### Plot 3: Variance effect ######

# Generate data points for the densities
grid <- seq(min(c(Ytrain0,Ytrain3))-1, max(c(Ytrain0,Ytrain3))+1, length.out = 100)
density1 <- dnorm(grid, mean = 0, sd = 1)
density2 <- dnorm(grid, mean = 0, sd = sqrt(1/x^2))

# Create a data frame for ggplot
data1 <- data.frame(x = grid, density = density1, group = "Mean 0")
data2 <- data.frame(x = grid, density = density2, group = "Mean 1.5")

# Combine the data frames
data <- rbind(data1, data2)

# Plot the Gaussian densities
p1<-ggplot() +
  geom_line(data = data1, aes(x = x, y = density), color = "blue", size = 1.5) +
  geom_line(data = data2, aes(x = x, y = density), color = "red", size = 1.5) +
  labs(title = "", x = "X", y = "Density") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text = element_text(size = 12)    # Increase axis text size
  )




p2<-tibble(
  Y     = c(Ytrain0,Ytrain3),
  w     = witness3[1, ],
  lower = witness3[2, ],
  upper = witness3[3, ]
) %>%
  ggplot(aes(x = Y, y = w)) +
  geom_line() +
  geom_hline(yintercept = 0, size = 1.5)+
  geom_line(aes(y = lower), lty = 2) +
  geom_line(aes(y = upper), lty = 2)


grid.arrange(p1, p2, ncol = 2)

### Plot 3: Variance effect ######


### Plot 4: Mean + Variance effect ######

# Generate data points for the densities
grid <- seq(min(c(Ytrain0,Ytrain4))-1, max(c(Ytrain0,Ytrain4))+1, length.out = 100)
density1 <- dnorm(grid, mean = 0, sd = 1)
density2 <- dnorm(grid, mean = x, sd = x)

# Create a data frame for ggplot
data1 <- data.frame(x = grid, density = density1, group = "Mean 0")
data2 <- data.frame(x = grid, density = density2, group = "Mean 1.5")

# Combine the data frames
data <- rbind(data1, data2)

# Plot the Gaussian densities
p1<-ggplot() +
  geom_line(data = data1, aes(x = x, y = density), color = "blue", size = 1.5) +
  geom_line(data = data2, aes(x = x, y = density), color = "red", size = 1.5) +
  labs(title = "", x = "X", y = "Density") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text = element_text(size = 12)    # Increase axis text size
  )




p2<-tibble(
  Y     = c(Ytrain0,Ytrain4),
  w     = witness4[1, ],
  lower = witness4[2, ],
  upper = witness4[3, ]
) %>%
  ggplot(aes(x = Y, y = w)) +
  geom_line() +
  geom_hline(yintercept = 0, size = 1.5)+
  geom_line(aes(y = lower), lty = 2) +
  geom_line(aes(y = upper), lty = 2)

grid.arrange(p1, p2, ncol = 2)

### Plot 4: Mean+ Variance effect ######


#####Study behavior of two separate DRFs###
Y1<-Ytrain1

fitDRF0 <- drf(X=Xtrain[Wtrain==0,], Y=Ytrain0, num.trees = 3000, ci.group.size = 100)
fitDRF1 <- drf(X=Xtrain[Wtrain==1,], Y=Y1, num.trees = 3000, ci.group.size = 100)

W0<-predict(fitDRF0, newdata=x)

data<-data.frame(W=Wtrain, X=Xtrain, Y=rbind(Ytrain0,Y1))

res<-drf_separate(data, num_trees=3000, Xtest=x, ci_group_size=100)


psep<-tibble(
  Y     = c(Ytrain0,Y1),
  w     = res$witness,
  lower = res$lower,
  upper = res$upper
) %>%
  ggplot(aes(x = Y, y = w)) +
  geom_line() +
  geom_hline(yintercept = 0, size = 1.5)+
  geom_line(aes(y = lower), lty = 2) +
  geom_line(aes(y = upper), lty = 2)


grid.arrange(psep,pjoint , ncol = 2)



}


if ("simulation_illustration" %in% run){
  
  source("simulation_study_setup.R")
  require(Matrix)
  N<-1000
  d<-5
  seed<-1
  num_trees<-2500
  
  d = c(5)
  num_trees = c(50 * 50)
  ci_group_size<-50
  
  dat <- simulate_both(N, d, seed = seed)
  truth <- both_truth
  dat$index <- 1:nrow(dat)
  
  Xtest <- matrix(c(0.7, 0.3, 0.5, 0.68, 0.43), nrow = 1)
  
  # Our method
  combined <- drf_combined(dat, num_trees, Xtest, ci_group_size, witness_type = "original")
  combined$truth <- truth(dat$Y)
  
  # Separate fits
  separate <- drf_separate(dat, num_trees, Xtest, ci_group_size)
  separate$truth <- truth(dat$Y)
  
  
  combined$Y<-dat$Y
  separate$Y<-dat$Y
  

  
  
p1<-  combined %>%
    ggplot(aes(x = Y, y = witness)) +
    geom_line(size = 1.5) +  # Thicker main line
    geom_line(aes(y = lower), lty = 2, size = 1.5) +  # Thicker dashed lines
    geom_line(aes(y = upper), lty = 2, size = 1.5) +  # Thicker dashed lines
    geom_line(aes(y = truth), size = 1.5, col="gray") +  # Thicker dashed lines
    geom_hline(yintercept = 0, size = 1.5)+
    scale_x_continuous() +
    scale_y_continuous(limits = c(-1,1)) +
    labs(y = NULL) +  # Remove y-axis label+
    theme(
      axis.title = element_text(size = 14),  # Increase axis title size
      axis.text = element_text(size = 12)    # Increase axis text size
    )
  
  
p2<-  separate %>%
    ggplot(aes(x = Y, y = witness)) +
    geom_line(size = 1.5) +  # Thicker main line
    geom_line(aes(y = lower), lty = 2, size = 1.5) +  # Thicker dashed lines
    geom_line(aes(y = upper), lty = 2, size = 1.5) +  # Thicker dashed lines
    geom_line(aes(y = truth), size = 1.5, col="gray") +  # Thicker dashed lines
    geom_hline(yintercept = 0, size = 1.5)+
    scale_x_continuous() +
    scale_y_continuous(limits = c(-1,1)) +
    labs(y = NULL) +  # Remove y-axis label+
    theme(
      axis.title = element_text(size = 14),  # Increase axis title size
      axis.text = element_text(size = 12)    # Increase axis text size
    )
  

grid.arrange(p1, p2, ncol = 2)

  
}


if ("wage" %in% run){


######################
#####Wage Data #######
######################

load("wage_data.Rdata")



which = rep(TRUE, nrow(wage))
which = which & (wage$age >= 17)
which = which & (wage$weeks_worked > 48)
which = which & (wage$hours_worked > 16)
which = which & (wage$employment_status == 'employed')
which = which & (wage$employer != 'self-employed')
which[is.na(which)] = FALSE

data = wage[which, ]
sum(is.na(data))
colSums(is.na(data))
rownames(data) = 1:nrow(data)

data$log_wage = log(data$salary / (data$weeks_worked * data$hours_worked))

## Prepare data and fit drf

X = data[,c(
  'age',
  'race',
  'hispanic_origin',
  'citizenship',
  'nativity',

  'marital',
  'family_size',
  'children',

  'education_level',
  'english_level',

  'economic_region'
)]
X$occupation = unlist(lapply(as.character(data$occupation), function(s){return(substr(s, 1, 2))}))
X$occupation = as.factor(X$occupation)
X$industry = unlist(lapply(as.character(data$industry), function(s){return(substr(s, 1, 2))}))
X$industry[X$industry %in% c('32', '33', '3M')] = '31'
X$industry[X$industry %in% c('42')] = '41'
X$industry[X$industry %in% c('45', '4M')] = '44'
X$industry[X$industry %in% c('49')] = '48'
X$industry[X$industry %in% c('92')] = '91'
X$industry = as.factor(X$industry)
X=dummy_cols(X, remove_selected_columns=TRUE)
X = as.matrix(X)

Y = data[,c('log_wage')]
W = data[,c('sex')]


x<-X[1+10000,,drop=F]


point_description = function(test_point){
  out = ''

  out = paste(out, 'job: ', test_point$occupation_description[1], sep='')
  out = paste(out, '\nindustry: ', test_point$industry_description[1], sep='')

  out = paste(out, '\neducation: ', test_point$education[1], sep='')
  out = paste(out, '\nemployer: ', test_point$employer[1], sep='')
  out = paste(out, '\nregion: ', test_point$economic_region[1], sep='')

  out = paste(out, '\nmarital: ', test_point$marital[1], sep='')
  out = paste(out, '\nfamily_size: ', test_point$family_size[1], sep='')
  out = paste(out, '\nchildren: ', test_point$children[1], sep='')

  out = paste(out, '\nnativity: ', test_point$nativity[1], sep='')
  out = paste(out, '\nhispanic: ', test_point$hispanic_origin[1], sep='')
  out = paste(out, '\nrace: ', test_point$race[1], sep='')
  out = paste(out, '\nage: ', test_point$age[1], sep='')

  return(out)
}

point_description(data[1+10000,])


i<-0
witnesslist<-list()
Ylist<-list()
for (n in c(500,1000,2000)){
i<-i+1
train_idx = sample(1:10000, n, replace=FALSE)

## Focus on training data
Ytrain=Y[train_idx]

#Ytrain<-scale(Ytrain)

Wtrain<-W[train_idx]
Wtrain<- matrix(ifelse(Wtrain == 'male', 1, 0), ncol=1)
#as.numeric(W[train_idx])
Xtrain=X[train_idx,]

fit <- drf(X=Xtrain, Y=Ytrain, W=Wtrain, num.trees = 3000, ci.group.size = 150)


witnesslist[[i]]<- predict_witness(fit, alpha = 0.05, newdata = x, newtreatment = matrix(1))
Ylist[[i]] <- Ytrain
}

for (i in 1:length(witnesslist)){

witness <- witnesslist[[i]]
Ytrain<- Ylist[[i]]
data <- tibble(
  Y     = Ytrain,
  w     = witness[1, ],
  lower = witness[2, ],
  upper = witness[3, ]
)


p <- data %>%
  ggplot(aes(x = Y, y = w)) +
  geom_line(size = 1.5) +  # Thicker main line
  geom_line(aes(y = lower), lty = 2, size = 1.5) +  # Thicker dashed lines
  geom_line(aes(y = upper), lty = 2, size = 1.5) +  # Thicker dashed lines
  geom_hline(yintercept = 0, size = 1.5)+
  scale_x_continuous(limits =  range(Ylist[[length(witnesslist)]])) +
  scale_y_continuous(limits = c(-0.4,0.4)) +
  labs(y = NULL) +  # Remove y-axis label+
  theme(
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text = element_text(size = 12)    # Increase axis text size
  )


# # Create the plot
# p <- data %>%
#   ggplot(aes(x = Y, y = w)) +
#   geom_line() +
#   geom_line(aes(y = lower), lty = 2) +
#   geom_line(aes(y = upper), lty = 2)

# Print the plot
print(p)
}

}


if ("lalonde" %in% run){

########################
#Lalonde data
##########################

set.seed(10)
###https://search.r-project.org/CRAN/refmans/designmatch/html/lalonde.html
data(lalonde)

lalonde$employed74<-ifelse(lalonde$re74 > 0, 1, 0)
lalonde$employed75<-ifelse(lalonde$re75 > 0, 1, 0)



#lalonde[lalonde$re74==0,"re74"] <- NA
#lalonde[lalonde$re75==0,"re75"] <- NA

XYW<-lalonde[lalonde$re78 > 0,]


W<-as.matrix(XYW[,"treatment"])
Y<-as.matrix(XYW[, "re78"])
X<-as.matrix(XYW[,!(colnames(XYW)%in% c("treatment", "re78"))])


fit <- drf(X=X, Y=Y, W=W, num.trees = 3000, ci.group.size = 150)

x<-as.matrix(lalonde[lalonde$re78 ==0,][1,!(colnames(lalonde)%in% c("treatment", "re78"))])

witness <- predict_witness(fit, alpha = 0.05, newdata = x, newtreatment = matrix(1))
data <- tibble(
  Y     = Y,
  w     = witness[1, ],
  lower = witness[2, ],
  upper = witness[3, ]
)


p <- data %>%
  ggplot(aes(x = Y, y = w)) +
  geom_line(size = 1.5) +  # Thicker main line
  geom_line(aes(y = lower), lty = 2, size = 1.5) +  # Thicker dashed lines
  geom_line(aes(y = upper), lty = 2, size = 1.5) +  # Thicker dashed lines
  geom_hline(yintercept = 0, size = 1.5)+
  scale_x_continuous() +
  scale_y_continuous(limits = c(-0.4,0.4)) +
  labs(y = NULL) +  # Remove y-axis label+
  theme(
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text = element_text(size = 12)    # Increase axis text size
  )

print(p)

}

if ("pension" %in% run){

######################
#####Pension 401(k) data set #######
######################

library(hdm)
data(pension)

pension<-pension[complete.cases(pension),]

X<-pension[1:nrow(pension), c("age", "inc", "fsize", "educ", "marr", "twoearn", "hown", "db", "pira")]

W<-pension[1:nrow(pension),"e401", drop=F]

# Y=(Net Financial Assets, Net Non-401(k) Financial Assets, Total Wealth)
Y<-pension[1:nrow(pension),c("net_tfa", "net_nifa", "tw")]
#Y<-pension[2:nrow(pension),c( "tw"), drop=F]
Y<-scale(Y)

## As in cite(...) we consider the free wealth measures, Net Financial Assets, Net Non-401(k) Financial
## Assets and Total Wealth as the dependent variables.
## While citet(Hansen) use separate quantile regressions and cite(...) combine
## the wealth measures into one variable beforehand, we use the ability of Causal-DRF
## to estimate the CKTE of the three wealth measures in total.

testid<-c(1:10)


#
# i<-0
# witnesslist<-list()
# Ylist<-list()
# for (n in 8000){
#   i<-i+1
#   train_idx = sample(1:nrow(X)-1, n, replace=FALSE)
#
#   ## Focus on training data
#   Ytrain=as.matrix(Y[train_idx,, drop=F])
#   Wtrain<-as.matrix(W[train_idx,,drop=F])
#   #as.numeric(W[train_idx])
#   Xtrain=as.matrix(X[train_idx,])
#
#   fit <- drf(X=Xtrain, Y=Ytrain[,1,drop=F], W=Wtrain, num.trees = 4000, ci.group.size = 100)
#
#
#   witnesslist[[i]]<- predict_witness(fit, alpha = 0.05, newdata = x, newtreatment = matrix(1))
#   Ylist[[i]] <- Ytrain[,1,drop=F]
# }
#
# for (i in 1:length(witnesslist)){
#
#   witness <- witnesslist[[i]]
#   data <- tibble(
#     Y     = Ylist[[i]],
#     w     = witness[1, ],
#     lower = witness[2, ],
#     upper = witness[3, ]
#   )
#
#
#   p <- data %>%
#     ggplot(aes(x = rowSums(Y), y = w)) +
#     geom_line(size = 1.5) +  # Thicker main line
#     geom_line(aes(y = lower), lty = 2, size = 1.5) +  # Thicker dashed lines
#     geom_line(aes(y = upper), lty = 2, size = 1.5) +  # Thicker dashed lines
#     geom_hline(yintercept = 0, size = 1.5)+
#     scale_x_continuous(limits =  range(Ylist[[length(witnesslist)]])) +
#     scale_y_continuous(limits = c(-0.4,0.4)) +
#     labs(y = NULL) +  # Remove y-axis label+
#     theme(
#       axis.title = element_text(size = 14),  # Increase axis title size
#       axis.text = element_text(size = 12)    # Increase axis text size
#     )
#   print(p)
# }


train_idx =setdiff(1:nrow(X),testid)#sample(2:nrow(X), 8000, replace=FALSE)


## Focus on training data
Ytrain=as.matrix(Y[train_idx,, drop=F])
Wtrain<-as.matrix(W[train_idx,,drop=F])
#as.numeric(W[train_idx])
Xtrain=as.matrix(X[train_idx,])


fit <- drf(X=Xtrain, Y=Ytrain, W=Wtrain, num.trees = 4000, ci.group.size = 100)

###Change testpoint here: We use 1 and 10
x<-as.matrix(pension[testid[1], c("age", "inc", "fsize", "educ", "marr", "twoearn", "hown", "db", "pira")])
#x<-as.matrix(pension[testid, c("age", "inc", "fsize", "educ", "marr", "twoearn", "hown", "db", "pira")])

#wx0<- predict(fit, newdata=x, newtreatment=0, bootstrap=F)$weights
#wx1<- predict(fit, newdata=x, newtreatment=1, bootstrap=F)$weights
#wx<-wx1 - wx0

wx<-predict(fit, newdata=x, newtreatment=NULL, bootstrap=F)$weights


#wxSb0<- predict(fit, newdata=x, newtreatment=0, bootstrap=T)$weights
#wxSb1<- predict(fit, newdata=x, newtreatment=1, bootstrap=T)$weights
#wxSb<-lapply(1:length(wxSb0), function(j) wxSb1[[j]] - wxSb0[[j]]  )

wxSb<-predict(fit, newdata=x, newtreatment=NULL, bootstrap=T)$weights



bandwidth_Y <- fit$bandwidth
k_Y <- rbfdot(sigma = 1/(2*bandwidth_Y^2) )


##Do Test of equality
Ky <- t(kernlab::kernelMatrix(k_Y, Ytrain, y = Ytrain))
H0list<-do.call(c, lapply(wxSb, function(w) as.numeric((w-wx)%*%Ky%*%t(w-wx)) ))
q<-quantile(H0list, 1-0.05)

teststat<-as.numeric(wx%*%Ky%*%t(wx))

print(ifelse(teststat >= q, "We reject at 5%", "We do not reject at 5%" )) ##We reject a t the 5 % level!

p<-list()

##Only plot first dimension
for (j in 1:3){


bandwidth_Yj <- drf:::medianHeuristic(Ytrain[,j])
k_Yj <- rbfdot(sigma = 1/(2*bandwidth_Yj^2) )


Kyj <- kernlab::kernelMatrix(k_Y, Ytrain[,j], y =  Ytrain[,j])
H0listj<-do.call(c, lapply(wxSb, function(w) as.numeric((w-wx)%*%Kyj%*%t(w-wx)) ))
qj<-quantile(H0listj, 1-0.05)


##Make this into an actual function!!! and then a grid from min to max of Y.
fvals<- function(y){

  Kyjf <- kernlab::kernelMatrix(k_Yj, Ytrain[,j], y =  y)

  return(as.numeric(unlist(wx%*%Kyjf)))
}

yseq<-seq(min(Ytrain[,j]), max(Ytrain[,j]), by=0.005    )

data <- tibble(
  Y     = yseq,
  w     = fvals(yseq),
  lower = fvals(yseq)-sqrt(qj),
  upper = fvals(yseq)+sqrt(qj)
)


p[[j]] <- data %>%
  ggplot(aes(x = Y, y = w)) +
  geom_line(size = 1.5) +  # Thicker main line
  geom_line(aes(y = lower), lty = 2, size = 1.5) +  # Thicker dashed lines
  geom_line(aes(y = upper), lty = 2, size = 1.5) +  # Thicker dashed lines
  geom_hline(yintercept = 0, size = 1.5)+
  scale_x_continuous(limits = c(-5,5) ) + # range(yseq)
  scale_y_continuous(limits = c(-0.2,0.2)) +
  labs(y = NULL) +  # Remove y-axis label+
  theme(
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text = element_text(size = 12)    # Increase axis text size
  )

}



grid.arrange(p[[1]], p[[2]], p[[3]], ncol = 3)


vimp<-compute_causaldrf_vimp(X, Y, W, X_test = NULL, num.trees = 300, silent = FALSE)

}







