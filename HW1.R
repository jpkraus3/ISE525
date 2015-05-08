### ISE525 HW 1 in R

## Basic Equations for Multiple Linear Regression
# Y = X*Beta + error
#
# Y - response column vector, length 'n' responses
#
# X - design matrix, size 'n x k' responses by variables
# X takes the form [1, x(1,1), x(1,2), ..., x(1,k)
#                   1, x(2,1), x(2,2), ..., x(2,k)
#                   ...                     
#                   1, x(n,1), x(n,2), ..., x(n,k)]
#
# Beta - coefficient column vector, length 'k' variables
#
# error - error column vector, length 'n'
# Multiple regression modelling is the application of the 'least squares' estimator to find the vector Beta that minimizes L = sum((error_i)^2, i=[1:n]) = t(error) %*% error = t(Y - X %*% Beta) %*% (Y - X %*% Beta)
#
# The least squares estimator is also the partial derivative solved for 0, for each coefficient
# d(L)/d(Beta)|b = -2*(t(X) %*% Y) + 2*(t(X) %*% X %*% Beta)
# this equation simplifies to: Beta = solve(t(X) %*% X) %*% t(X) %*% Y
#
# Beta = X'X
# Beta = [                 n,  t(X[,1]) %*% X[,2],  ...,  t(X[,1]) %*% X[,k],
#         t(X[,2]) %*% X[,1],  t(X[,2]) %*% X[,2],  ...,  t(X[,2]) %*% X[,k],
#         ...
#         t(X[,k]) %*% X[,1],  t(X[,k]) %*% X[,1],  ...,  t(X[,k]) %*% X[,k]  ]
#
# Y.hat = X %*% Beta
#
# residual = Y - Y.hat
#
#
## ANOVA Table Calculations
# Source      SS                      dof         MS
# Regression  Beta'*X'*Y - sum(Y)^2/n k           SSR/dofR
# Error       Y'*Y - Beta'*X'*Y       n-k-1       SSE/dofE
# Total       Y'*Y - sum(Y)^2/n       n-1         
#
# Remember: SSR + SSE = SST
#
# F-Test: is there a linear relation between some of the variables and the response?
# H0: Beta[1] = Beta[2] = ... = 0, H1: at least one Beta[k] =/= 0
# Fobs = MSR/MSE = (SSr/k) / (SSe/(n-k-1))
# Fcritical = F(alpha, k, n-k-1) = F(alpha, dofR, dofE)
# Reject H0 if Fobs > Fcritical
# 
# R.square = SSr/SSt = 1 - SSe/SSt
# adj.R.square = 1 - (SSe/(n-p)) / (SSt/(n-1)) = 1 - (n-1)/(n-p) *(1-R.square), p <- terms in model
#
# t-Test: which design variables relate to the response?
# H0: Beta[j] = 0, H1: Beta[j] =/= 0
# t_j = Beta[j]/sqrt(var(Beta[j])) = Beta[j]/sqrt(MSE*C[j,j])
# Reject H0 if |t_j| > t(alpha/2, n-k-1) 
#
# C[j,j] is variance-covariance matrix, Sigma
#
# Sigma = MSE*(X'X)^-1
# Sigma = [ var((X'X)^-1[0,0]),  cov((X'X)^-1[0,1]),  ...,  cov((X'X)^-1[0,k])
#           cov((X'X)^-1[0,1]),  var((X'X)^-1[1,1]),  ...,  cov((X'X)^-1[1,k])
#           cov((X'X)^-1[k,0])                     ,  ..., var((X'X)^-1[k,k]) ]
# 
# Also, sqrt(MSE*C[j,j]) is the standard error of the regression coefficient Beta[j]




## Start work
# Save the matrix in the datasets folder as a comma-delimited .csv
data <- read.csv(file = "C:/Users/kraus_000/Documents/R/Datasets/HW1.csv", header = T, sep = ",", dec = ".")
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

# Solve the Beta matrix by hand 
solve(t(X) %*% X) %*% t(X) %*% Y  # Beta vector
Beta <- solve(t(X) %*% X) %*% t(X) %*% Y

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


## ANOVA Table for Full Model

# Regression Model Terms
# SSr
t(Beta) %*% t(X) %*% Y - sum(Y)^2/n
SSr <- t(Beta) %*% t(X) %*% Y - sum(Y)^2/n
# dofR
k
# MSr
(t(Beta) %*% t(X) %*% Y - sum(Y)^2/n) / k
MSr <- (t(Beta) %*% t(X) %*% Y - sum(Y)^2/n) / k

# Error terms
# SSe
t(Y) %*% Y - t(Beta) %*% t(X) %*% Y
SSe <- t(Y) %*% Y - t(Beta) %*% t(X) %*% Y
# dofE
n-k-1
# MSe
(t(Y) %*% Y - t(Beta) %*% t(X) %*% Y) / (n-k-1)
MSe <- (t(Y) %*% Y - t(Beta) %*% t(X) %*% Y) / (n-k-1)

# Total Terms
t(Y) %*% Y - sum(Y)^2/n
SSt <- t(Y) %*% Y - sum(Y)^2/n
# dofT
n - 1

# F-observed = MSr/MSe
((t(Beta) %*% t(X) %*% Y - sum(Y)^2/n) / k) / ((t(Y) %*% Y - t(Beta) %*% t(X) %*% Y) / (n-k-1))
Fobs <- ((t(Beta) %*% t(X) %*% Y - sum(Y)^2/n) / k) / ((t(Y) %*% Y - t(Beta) %*% t(X) %*% Y) / (n-k-1))

# R.square = SSr/SSt = 1 - SSe/SSt
(t(Beta) %*% t(X) %*% Y - sum(Y)^2/n) / (t(Y) %*% Y - sum(Y)^2/n)
R.square <- (t(Beta) %*% t(X) %*% Y - sum(Y)^2/n) / (t(Y) %*% Y - sum(Y)^2/n)

# Adjusted R.square = 1 - (SSe/(n-p)) / (SSt/(n-1)) = 1 - (n-1)/(n-p) *(1-R.square), p <- terms in model
1 - ((n-1)/(n-k))*(1-R.square)
adj.R.square <- 1 - ((n-1)/(n-k))*(1-R.square)


# t-Tests for individual coefficients Beta[j], t
C <- solve(t(X)%*%X)
j <- 2  # Remember 1 is for constant, 2 is for x1, etc.
Beta[j] / sqrt(MSe*C[j,j])

