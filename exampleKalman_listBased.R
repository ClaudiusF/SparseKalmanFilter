library(Matrix)
library(KFAS)
library(FKF)

###################### time-invariant model ####################
mod2 <- SSModel(
  cbind(numeric(1000), numeric(1000))~
    SSMtrend(2, list(diag(2), diag(2)*0.01), a1=rnorm(4), P1=diag(4)*5, P1inf=diag(4)*0)+
    SSMseasonal(100, diag(2), a1=rnorm(198, sd=100), P1=diag(198), P1inf=diag(198)*0),
  H = diag(2)
)

# simulate y from ss form
set.seed(202303)
Y <- simulateSSM(mod2, "signals", conditional = FALSE)[,,1]
plot(ts(Y))

mod2$y = Y

#lists
T <- list()
for (i in 1:dim(mod2$T)[3]){
  T[i] <- list(Matrix(mod2$T[,,i], sparse = TRUE))  
}

Q <- list()
for (i in 1:dim(mod2$Q)[3]){
  Q[i] <- list(Matrix(mod2$Q[,,i], sparse = TRUE))
}

Z <- list()
for (i in 1:dim(mod2$Z)[3]){
  Z[i] <- list(Matrix(mod2$Z[,,i], sparse = TRUE))
}

H <- list()
for (i in 1:dim(mod2$H)[3]){
  H[i] <- list(Matrix(mod2$H[,,i], sparse = TRUE))
}

R <- list()
for (i in 1:dim(mod2$R)[3]){
  R[i] <- list(Matrix(mod2$R[,,i], sparse = TRUE))
}

V <- list()
for (i in 1:dim(mod2$Q)[3]){
  V[i] <- list(Matrix(R[[i]] %*% Q[[i]] %*% t(R[[i]]), sparse = TRUE))
}

a1 <- Matrix(mod2$a1, sparse = TRUE)
P1 <- Matrix(mod2$P1, sparse = TRUE)
Yt <- Matrix(t(Y), sparse = TRUE)
Ytd <- t(Y)

#benchmark
timeViaListe <- microbenchmark::microbenchmark(
  kfs  <- KFS(mod2, filtering = "state", smoothing = c("state", "disturbance"), simplify = FALSE),
  kff  <- kalman(T, V, Q, Z, H, R, a1, P1, Ytd, c("Smoother", "Disturbance")), # Armadillo sparse full filter
  times = 5
)

plot(timeViaListe)


#
#
#
#
#
#
########################### time-varying model ############
set.seed(1)
x1 <- rnorm(1000)
x1
x2 <- rnorm(1000)
x2

mod2 <- SSModel(
  cbind(numeric(1000), numeric(1000))~
    SSMregression(~ x1 + x2, list(diag(2), diag(2)*0.01), a1=rnorm(4), P1=diag(4)*1, P1inf=diag(4)*0) +
    SSMtrend(2, list(diag(2), diag(2)*0.01), a1=rnorm(4), P1=diag(4)*5, P1inf=diag(4)*0)+
    SSMseasonal(100, diag(2), a1=rnorm(198, sd=100), P1=diag(198), P1inf=diag(198)*0),
  H = diag(2)
)

# simulate y from ss form
set.seed(202303)
Y <- simulateSSM(mod2, "signals", conditional = FALSE)[,,1]
plot(ts(Y))

mod2$y = Y

#lists
T <- list()
for (i in 1:dim(mod2$T)[3]){
  T[i] <- list(Matrix(mod2$T[,,i], sparse = TRUE))  
}

Q <- list()
for (i in 1:dim(mod2$Q)[3]){
  Q[i] <- list(Matrix(mod2$Q[,,i], sparse = TRUE))
}

Z <- list()
for (i in 1:dim(mod2$Z)[3]){
  Z[i] <- list(Matrix(mod2$Z[,,i], sparse = TRUE))
}

H <- list()
for (i in 1:dim(mod2$H)[3]){
  H[i] <- list(Matrix(mod2$H[,,i], sparse = TRUE))
}

R <- list()
for (i in 1:dim(mod2$R)[3]){
  R[i] <- list(Matrix(mod2$R[,,i], sparse = TRUE))
}

V <- list()
for (i in 1:dim(mod2$Q)[3]){
  V[i] <- list(Matrix(R[[i]] %*% Q[[i]] %*% t(R[[i]]), sparse = TRUE))
}

a1 <- Matrix(mod2$a1, sparse = TRUE)
P1 <- Matrix(mod2$P1, sparse = TRUE)
Yt <- Matrix(t(Y), sparse = TRUE)
Ytd <- t(Y)

#benchmark
timeViaListe <- microbenchmark::microbenchmark(
  kfs  <- KFS(mod2, filtering = "state", smoothing = c("state", "disturbance"), simplify = FALSE),
  kff  <- kalman(T, V, Q, Z, H, R, a1, P1, Ytd, c("Smoother", "Disturbance")), # Armadillo sparse full filter
  times = 5
)

plot(timeViaListe)
