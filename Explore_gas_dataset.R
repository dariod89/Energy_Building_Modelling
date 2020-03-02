# A computer model simulating monthly gas consumption in a private house has been run N=1000 times, for N different
# choices of 8 model parameters. These are, for example, boiler efficiency, floor thickness, wall insulation, etc.
# Actual observations of monthly consumptions for the modelled house are also available. 

# This script explores the dataset of simulated gas consumptions, and identifies the main factors (and interactions
# thereof) which explain the observed variability in the simulated gas consumptions. An emulator of the simulated 
# consumptions is then built: the emulator predicts the simulator response at any input in the form of a random 
# variable.


# LOAD LIBRARIES AND CCUSTOM FUNCTIONS

library(xlsx)
library(leaps)
source("../Various_Functions.R")
source("../emulation.R")

####################################################################
##  LOAD THE DATA: SIMULATION INPUTS AND OUTPUTS, AND OBSERVATIONS  
####################################################################

month.names <- c("Jan", "Feb", "Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec")

# Inputs (1000 x 8)
Temp <- read.csv("Data/inputs-batch2.csv")
Inputs <- Temp[, -1]            # delete first column with numbers from 1 to 1000
var.names <- colnames(Inputs)   # store names of input variables

# Outputs (1000 x 12)
Temp <- read.csv("Data/results-batch2.csv", skip = 3008, nrow=1001)
Outputs_Gas_2 <- Temp[-1,]            # remove baseline simulation
Outputs_Gas_2[,"File.Name"] <- NULL   # remove first column with file name of runs
rownames(Outputs_Gas) <- 1:1000       # name the simulations by numbers
colnames(Outputs_Gas) <- month.names  # assign month names
rm(Temp)

# Observations (1 x 12)
Obs_Gas <- read.xlsx("Data/Observation vs Simulated - Hailiang.xlsx", sheetIndex =1, 
                      startRow = 19, endRow = 31, colIndex = 3) # dim: 12x1
Obs_Gas <- t(Obs_Gas)                 # transpose to have months in columns: 1x12
colnames(Obs_Gas) <- month.names      # name the column by months


###########################################
##  RESCALE INPUTS LINEARLY WITHIN [-1,1]  
###########################################

# Set the ranges...
N_var<-8
input_ranges<-matrix(NA, nrow=2, ncol=N_var)

input_ranges[1,1]<-17.5     # Heating setpoint [Celcius degrees]
input_ranges[2,1]<-20.5

input_ranges[1,2]<-0.6      # Gas boiler seasonal efficiency
input_ranges[2,2]<-0.75

input_ranges[1,3]<-0.04     # External wall thickness [m]
input_ranges[2,3]<-0.063

input_ranges[1,4]<-0.15     # ~Roof thickness [m]
input_ranges[2,4]<-0.21

input_ranges[1,5]<-0.045    # ~Floor thickness [m] 
input_ranges[2,5]<-0.055

input_ranges[1,6]<-0.2      # Infiltration [ac/h]
input_ranges[2,6]<-0.95 

input_ranges[1,7]<-6.15e-06 # DHW consumption [litre/day]
input_ranges[2,7]<-2.20e-05

input_ranges[1,8]<-1.05     # Cooking [W/m2]
input_ranges[2,8]<-6.3

# ... and rescale the input variables within these
for (i in 1:N_var){
  a <- input_ranges[1,i]
  b <- input_ranges[2,i]
  x <- Inputs_2[,i]
  Inputs_2[,i] <- 2*(x-a)/(b-a) -1
}
rm(input_ranges)


###################################################
##   ANALYSE THE SIMULATED INPUT-OUTPUT RELATION
###################################################

# Build a matrix of all 2-way interactions of factors (orthogonal polynomials of order 2)
Interactions <- poly(data.matrix(Inputs), degree=2)

month<-1
y <- Outputs_Gas[, month]
# The next line considers linear models with different number of covariates (from 1 to nvmax), 
# and selects best model of each size
# (for fixed size, best model is the same regardless of criterion used - adjr2, Cp etc)
L <- summary(regsubsets(y~., data=Interactions, method = "exhaustive", nvmax = 8))
ind <- which.max(L$adjr2)  # index (between 1 and nvmax) of model with max adj-R2
regr <- L$which[ind,-1]    # logical vector with regressors corresponding to selected model
names(regr[regr==T])       # just prints selected variables
fit <- lm(y~ ., data = as.data.frame(Interactions[,regr]))
summary(fit)

colnames(Inputs) <- c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8')

###########################################################################
# CARRY OUT EMULATION BASED ON REGRESSORS AND INTERACTIONS SELECTED ABOVE
###########################################################################

beta <- fit$coefficients
Cov.beta <- vcov(fit)

# Number of inputs at which the emulator will provide a prediction of simulated response 
N <- 1.e6

# N random points in the cube [-1,1]^Nvar (Nvar=8 in our case)
Test.points <- matrix(2*runif(N*N_var)-1, nrow = N) 
colnames(New.points) <- colnames(Inputs)

# Regressors associated to test points
Test.regr <- predict(Interactions, newdata = Test.points)
Test.regr <- cbind(1, Test.regr[, regr, drop=F])      # add a column of 1 for intercept

Design.points <- Inputs                               # design points used to train the emulator
Design.regr <- cbind(1, Interactions[,regr, drop=F])  # regressors: add a column of 1s for intercept

d <- 0.4* replicate(N_var, 1)     # Correlation lengths
sigma2.tot <- var(fit$residuals)  # Prior cumulative variance of homoschedastic Gaussian process
sigma2.tot <- sigma2.tot/4        # Specific value to be set via, eg, cross-validation
nugget <- 0.02                    # Fraction of residual variability not explained by regressors
nu <- nugget*sigma2.tot           # Variance of nugget term (Gaussian noise)
sig2 <- (1-nugget)*sigma2.tot     # Variance of Gaussian process

# Call the function which carries out emulation, on the selected regressors and on output y
res <- emul(Test.points, Test.regr, Design.points, Design.regr, y, beta, Cov.beta, 1*d, nu, sig2, 'exp2')




# All the following is yet to come

#################################################################################
###  VALIDATE THE BUILT EMULATOR (GRAPHICALLY AND THROUGH STATISTICAL MEASURES)

####################################################################################################
###   SELECT REGIONS OF THE INPUT SPACE WHICH MATCH THE OBSERVED DATA, VIA IMPLAUSIBILITY MEASURES

###########################################################################################
###   CHECK WHETHER NON-EMPTY REGIONS CAN BE FOUND WHICH MATCH ALL MONTHS SIMULTANEOUSLY

##################################################
###   REFOCUS OF SMALLER REGIONS







  
