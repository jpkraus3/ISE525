## Start work
# Save the matrix in the datasets folder as a comma-delimited .csv
data <- read.csv(file = "C:/Users/kraus_000/Documents/R/Code/ISE525/Data/HW1.csv", header = T, sep = ",", dec = ".")
data

# Data frames are easy to work with. Place the data into a dataframe
HW1 <- as.data.frame(x = data)
names(HW1) #displays the names of the columns, see which is variables and response
# Prints '[1] "y"  "x1" "x2" "x3"'

# Time for matrix algebra, define X and Y as matrices
X <- data.matrix(HW1[2:4])  # Columns 2 through 4 make up the X matrix
X <- cbind(I=1, X)  # Add column of 1's to the front of X for constant
Y <- data.matrix(HW1[1])  # Column 1 is the Y vector
n <- nrow(HW1)
k <- ncol(X) - 1


### 1 Estimate the parameters of a regression model that relates x1, x2, x3 to y
## Find Beta parameters from linear model Y = X %*% Beta + e, solve for Beta = (X'X)^-1 * X' * Y

# Solve (X'X) coefficients by hand
i <- 2  # Remember column 1 is constant coefficient, 2 is x1, ...
j <- 1
t(X[,i]) %*% X[,j]  # (X'X) coefficient for i, j
t(X) %*% X  # Entire X'X matrix
XpX <- t(X) %*% X 

# Solve the (X'Y) coefficients by hand
i <- 1  # Remember 1 is constant coefficient, and just sums Y's
t(X[,i]) %*% Y
t(X) %*% Y  # Entire X'Y vector
XpY <- t(X) %*% Y

# Solve the Beta vector by hand 
solve(t(X) %*% X) %*% t(X) %*% Y  # Beta vector
Beta <- solve(t(X) %*% X) %*% t(X) %*% Y


### 2 Test for significance of the regression model using the Analysis of Variance. Interpret at 5% significance level
## Fit ANOVA table SSr, SSe, SSt, compute regression F-statistic, compare to F-critical
# 

# ANOVA Table for Full Model
# Regression Model Terms
# SSr
SSr <- t(Beta) %*% t(X) %*% Y - sum(Y)^2/n
# dofR = k
# MSr
MSr <- (t(Beta) %*% t(X) %*% Y - sum(Y)^2/n) / k

# Error terms
# SSe
SSe <- t(Y) %*% Y - t(Beta) %*% t(X) %*% Y
# dofE = n-k-1 = n-p
# MSe
MSe <- (t(Y) %*% Y - t(Beta) %*% t(X) %*% Y) / (n-k-1)

# Total Terms
SSt <- t(Y) %*% Y - sum(Y)^2/n
# dofT = n-1, there is no MSt

# F-observed = MSr/MSe
# F-test H0: all regression terms equal to zero, H1: at least one is not zero
# Reject H0 at F-observed > Fcritical ==> F(alpha, dof MSr, dof MSe)
((t(Beta) %*% t(X) %*% Y - sum(Y)^2/n) / k) / ((t(Y) %*% Y - t(Beta) %*% t(X) %*% Y) / (n-k-1))
Fobs <- ((t(Beta) %*% t(X) %*% Y - sum(Y)^2/n) / k) / ((t(Y) %*% Y - t(Beta) %*% t(X) %*% Y) / (n-k-1))
Fobs  # equals 147.36
# Fcritical = F(.05, 3, 12) = 3.49 < Fobs, reject null hypothesis, model is significant

### 3 Test for significance on individual components using t-Tests
# t-test H0: Beta[j] = 0, H1: Beta[j] =/= 0
# Reject H0 if abs(t_j) > tcritical ==> t(alpha/2, n-k-1)
# t_j = Beta[j]/sqrt(var(Beta[j])) = Beta[j]/sqrt(MSE*C[j,j])

# t-Tests for individual coefficients Beta[j], t
C <- solve(t(X)%*%X)

# Remember 1 is for constant, 2 is for x1, etc.
# tcritical = t(.025, 12) = 2.179
t.x1 <- Beta[2] / sqrt(MSe*C[2,2])  # abs(t.x1) = 14.41 > 2.179, reject H0, x1 is significant
t.x2 <- Beta[3] / sqrt(MSe*C[3,3])  # abs(t.x2) = 11.21 > 2.179, reject H0, x2 is significant
t.x3 <- Beta[4] / sqrt(MSe*C[4,4])  # abs(t.x3) = 10.42 > 2.179, reject H0, x3 is significant


## NOT NECESSARY FOR SOLUTION
# Solve the Y.hat matrix by hand
# Solve residual vector by hand
# Solve for R-square by hand
# Solve for adjusted R-square by hand
#####
# Solve the Y.hat matrix by hand
j <- 1  # Row corresponding to Y.hat treatment combination
X[j,] %*% Beta  # Y.hat for tc j
X %*% Beta  # Entire Y.hat vector
Y.hat <- X %*% Beta

# Solve residual vector by hand
j <- 1  # Row corresponding to residual treatment combination  
Y[j] - Y.hat[j]
Y - Y.hat  # Entire residual vector
residual <- Y - Y.hat

# R.square = SSr/SSt = 1 - SSe/SSt
(t(Beta) %*% t(X) %*% Y - sum(Y)^2/n) / (t(Y) %*% Y - sum(Y)^2/n)
R.square <- (t(Beta) %*% t(X) %*% Y - sum(Y)^2/n) / (t(Y) %*% Y - sum(Y)^2/n)

# Adjusted R.square = 1 - (SSe/(n-p)) / (SSt/(n-1)) = 1 - (n-1)/(n-p) *(1-R.square), p <- terms in model
1 - ((n-1)/(n-k))*(1-R.square)
adj.R.square <- 1 - ((n-1)/(n-k))*(1-R.square)