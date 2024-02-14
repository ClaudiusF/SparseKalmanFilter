library(Matrix)
library(KFAS)
library(FKF)
library(spatstat.sparse)

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

# matrix 3d array
T <- as.sparse3Darray(mod2$T[,,])
Z <- as.sparse3Darray(mod2$Z[,,])
H <- as.sparse3Darray(mod2$H[,,])
Q <- as.sparse3Darray(mod2$Q[,,])
R <- as.sparse3Darray(mod2$R[,,])

for (i in 1:dim(mod2$Q)[3]){
  V <- as.sparse3Darray(R[,,i] %*% Q[,,i] %*% t(R[,,i]))
}

a1 <- Matrix(mod2$a1, sparse = TRUE)
P1 <- Matrix(mod2$P1, sparse = TRUE)
Yt <- Matrix(t(Y), sparse = TRUE)

Ytd <- t(Y)

X <- data.frame(NULL, NULL, NULL)
iX <- data.frame(NULL, NULL, NULL)
iT <- data.frame(NULL, NULL, NULL)
iV <- data.frame(NULL, NULL, NULL)
iZ <- data.frame(NULL, NULL, NULL)
iH <- data.frame(NULL, NULL, NULL)
iQ <- data.frame(NULL, NULL, NULL)
iR <- data.frame(NULL, NULL, NULL)

p1 = 1
p2 = 1
p3 = 1
p4 = 1
p5 = 1
p6 = 1

fZ <- 1
fT <- 1
fV <- 1
fH <- 1
fQ <- 1
fR <- 1
################### matrix generation (X, iZ, iT, ...)

for (k in 1:(dim(Y)[1] + 1)){
  i=0
  if (k != (dim(Y)[1] + 1)){
    if (fZ == 1 && dim(Z)[3]>1){
      for (s in 1:dim(Z)[1]){
        for (j in 1:NCOL(t(Z[s,,k]))){
          if (Z[s,j,k] != 0) {
            i <- i + 1
            p1 <- p1 + 1
            X[i,k] <- Z[s,j,k]
            iZ[p1,1] <- k - 1
            iZ[p1,2] <- j - 1
            iZ[p1,3] <- s - 1
          }
        } 
      }
    }  else { fZ <- 0
    Z <- Matrix(Z[,,1], sparse = TRUE)
    iZ <- Matrix(NA, nrow=2, ncol=2)
    }
  }    else { p1 <- p1 + 1
  iZ[p1, ] <- NA
  Z <- Matrix(Z[,,1], sparse = TRUE)
  }
  
  if (k != (dim(Y)[1] + 1)){
    if (fT == 1 && dim(T)[3]>1){
      for (s in 1:dim(T)[1]){
        for(j in 1:NCOL(t(T[s,,k]))){
          if (T[s,j,k] != 0) {
            i <- i + 1
            p2 <- p2 + 1
            X[i,k] <- T[s,j,k]
            iT[p2,1] <- k - 1
            iT[p2,2] <- j - 1
            iT[p2,3] <- s - 1
          }
        }
      }
    }  else { fT <- 0
    T <- Matrix(T[,,1], sparse = TRUE)
    iT <- Matrix(NA, nrow=2, ncol=2)
    }
  } else { p2 <- p2 + 1
  iT[p2, ] <- NA
  T <- Matrix(T[,,1], sparse = TRUE)
  }
  
  if (k != (dim(Y)[1] + 1)){
    if (fV == 1 && dim(V)[3]>1){
      for (s in 1:dim(V)[1]){
        for(j in 1:NCOL(t(V[s,,k]))){
          if (V[s,j,k] != 0) {
            i <- i + 1
            p3 <- p3 + 1
            X[i,k] <- V[s,j,k]
            iV[p3,1] <- k - 1
            iV[p3,2] <- j - 1
            iV[p3,3] <- s - 1
          }
        }
      }
    } else { fV <- 0
    V <- Matrix(V[,,1], sparse = TRUE)
    iV <- Matrix(NA, nrow=2, ncol=2)
    }
  } else { p3 <- p3 + 1
  iV[p3, ] <- NA
  V <- Matrix(V[,,1], sparse = TRUE)
  }
  
  if (k != (dim(Y)[1] + 1)){
    if (fH == 1 && dim(H)[3]>1){
      for (s in 1:dim(H)[1]){
        for(j in 1:NCOL(t(H[s,,k]))){
          if (H[s,j,k] != 0) {
            i <- i + 1
            p4 <- p4 + 1
            X[i,k] <- H[s,j,k]
            iH[p4,1] <- k - 1
            iH[p4,2] <- j - 1
            iH[p4,3] <- s - 1
          }
        }
      }
    } else {fH <- 0
    H<- Matrix(H[,,1], sparse = TRUE)
    iH <- Matrix(NA, nrow=2, ncol=2)
    }
  } else { p4 <- p4 + 1
  iH[p4, ] <- NA
  H<- Matrix(H[,,1], sparse = TRUE)
  }
  
  if (k != (dim(Y)[1] + 1)){
    if (fQ == 1 && dim(Q)[3]>1){
      for (s in 1:dim(Q)[1]){
        for(j in 1:NCOL(t(Q[s,,k]))){
          if (Q[s,j,k] != 0) {
            i <- i + 1
            p5 <- p5 + 1
            X[i,k] <- Q[s,j,k]
            iQ[p5,1] <- k - 1
            iQ[p5,2] <- j - 1
            iQ[p5,3] <- s - 1
          }
        }
      }
    } else {fQ <- 0
    Q <- Matrix(Q[,,1], sparse = TRUE)
    iQ <- Matrix(NA, nrow=2, ncol=2)
    }
  } else { p5 <- p5 + 1
  iQ[p5, ] <- NA
  Q <- Matrix(Q[,,1], sparse = TRUE)
  }
  
  if (k != (dim(Y)[1] + 1)){
    if (fR == 1 && dim(R)[3]>1){
      for (s in 1:dim(R)[1]){
        for(j in 1:NCOL(t(R[s,,k]))){
          if (R[s,j,k] != 0) {
            i <- i + 1
            p6 <- p6 + 1
            X[i,k] <- R[s,j,k]
            iR[p6,1] <- k - 1
            iR[p6,2] <- j - 1
            iR[p6,3] <- s - 1
          }
        }
      }
    } else {fR <- 0
    R <- Matrix(R[,,1], sparse = TRUE)
    iR <- Matrix(NA, nrow=2, ncol=2)
    }
  } else { p6 <- p6 + 1
  iR[p6, ] <- NA
  R <- Matrix(R[,,1], sparse = TRUE)
  }
}
##################################

################### transform dataframe into matrix
X <- as.matrix(X)
iZ <- as.matrix(iZ)
iV <- as.matrix(iV)
iH <- as.matrix(iH)
iT <- as.matrix(iT)
iQ <- as.matrix(iQ)
iR <- as.matrix(iR)

#benchmark
timeViaIndici <- microbenchmark::microbenchmark(
  kfs  <- KFS(mod2, filtering = "state", smoothing = c("state", "disturbance"), simplify = FALSE),
  kff  <- kalman(T, V, Z, H, Q, R, a1, P1, X, iZ, iH, iV, iT, iQ, iR, Ytd, c("Smoother","Disturbance")),
  times = 3
)

plot(timeViaIndici)

#
#
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

# matrix 3d array
T <- as.sparse3Darray(mod2$T[,,])
Z <- as.sparse3Darray(mod2$Z[,,])
H <- as.sparse3Darray(mod2$H[,,])
Q <- as.sparse3Darray(mod2$Q[,,])
R <- as.sparse3Darray(mod2$R[,,])

for (i in 1:dim(mod2$Q)[3]){
  V <- as.sparse3Darray(R[,,i] %*% Q[,,i] %*% t(R[,,i]))
}

a1 <- Matrix(mod2$a1, sparse = TRUE)
P1 <- Matrix(mod2$P1, sparse = TRUE)
Yt <- Matrix(t(Y), sparse = TRUE)

Ytd <- t(Y)

X <- data.frame(NULL, NULL, NULL)
iX <- data.frame(NULL, NULL, NULL)
iT <- data.frame(NULL, NULL, NULL)
iV <- data.frame(NULL, NULL, NULL)
iZ <- data.frame(NULL, NULL, NULL)
iH <- data.frame(NULL, NULL, NULL)
iQ <- data.frame(NULL, NULL, NULL)
iR <- data.frame(NULL, NULL, NULL)

p1 = 1
p2 = 1
p3 = 1
p4 = 1
p5 = 1
p6 = 1

fZ <- 1
fT <- 1
fV <- 1
fH <- 1
fQ <- 1
fR <- 1
################### matrix generation (X, iZ, iT, ...)

for (k in 1:(dim(Y)[1] + 1)){
  i=0
  if (k != (dim(Y)[1] + 1)){
    if (fZ == 1 && dim(Z)[3]>1){
      for (s in 1:dim(Z)[1]){
        for (j in 1:NCOL(t(Z[s,,k]))){
          if (Z[s,j,k] != 0) {
            i <- i + 1
            p1 <- p1 + 1
            X[i,k] <- Z[s,j,k]
            iZ[p1,1] <- k - 1
            iZ[p1,2] <- j - 1
            iZ[p1,3] <- s - 1
          }
        } 
      }
    }  else { fZ <- 0
    Z <- Matrix(Z[,,1], sparse = TRUE)
    iZ <- Matrix(NA, nrow=2, ncol=2)
    }
  }    else { p1 <- p1 + 1
  iZ[p1, ] <- NA
  Z <- Matrix(Z[,,1], sparse = TRUE)
  }
  
  if (k != (dim(Y)[1] + 1)){
    if (fT == 1 && dim(T)[3]>1){
      for (s in 1:dim(T)[1]){
        for(j in 1:NCOL(t(T[s,,k]))){
          if (T[s,j,k] != 0) {
            i <- i + 1
            p2 <- p2 + 1
            X[i,k] <- T[s,j,k]
            iT[p2,1] <- k - 1
            iT[p2,2] <- j - 1
            iT[p2,3] <- s - 1
          }
        }
      }
    }  else { fT <- 0
    T <- Matrix(T[,,1], sparse = TRUE)
    iT <- Matrix(NA, nrow=2, ncol=2)
    }
  } else { p2 <- p2 + 1
  iT[p2, ] <- NA
  T <- Matrix(T[,,1], sparse = TRUE)
  }
  
  if (k != (dim(Y)[1] + 1)){
    if (fV == 1 && dim(V)[3]>1){
      for (s in 1:dim(V)[1]){
        for(j in 1:NCOL(t(V[s,,k]))){
          if (V[s,j,k] != 0) {
            i <- i + 1
            p3 <- p3 + 1
            X[i,k] <- V[s,j,k]
            iV[p3,1] <- k - 1
            iV[p3,2] <- j - 1
            iV[p3,3] <- s - 1
          }
        }
      }
    } else { fV <- 0
    V <- Matrix(V[,,1], sparse = TRUE)
    iV <- Matrix(NA, nrow=2, ncol=2)
    }
  } else { p3 <- p3 + 1
  iV[p3, ] <- NA
  V <- Matrix(V[,,1], sparse = TRUE)
  }
  
  if (k != (dim(Y)[1] + 1)){
    if (fH == 1 && dim(H)[3]>1){
      for (s in 1:dim(H)[1]){
        for(j in 1:NCOL(t(H[s,,k]))){
          if (H[s,j,k] != 0) {
            i <- i + 1
            p4 <- p4 + 1
            X[i,k] <- H[s,j,k]
            iH[p4,1] <- k - 1
            iH[p4,2] <- j - 1
            iH[p4,3] <- s - 1
          }
        }
      }
    } else {fH <- 0
    H<- Matrix(H[,,1], sparse = TRUE)
    iH <- Matrix(NA, nrow=2, ncol=2)
    }
  } else { p4 <- p4 + 1
  iH[p4, ] <- NA
  H<- Matrix(H[,,1], sparse = TRUE)
  }
  
  if (k != (dim(Y)[1] + 1)){
    if (fQ == 1 && dim(Q)[3]>1){
      for (s in 1:dim(Q)[1]){
        for(j in 1:NCOL(t(Q[s,,k]))){
          if (Q[s,j,k] != 0) {
            i <- i + 1
            p5 <- p5 + 1
            X[i,k] <- Q[s,j,k]
            iQ[p5,1] <- k - 1
            iQ[p5,2] <- j - 1
            iQ[p5,3] <- s - 1
          }
        }
      }
    } else {fQ <- 0
    Q <- Matrix(Q[,,1], sparse = TRUE)
    iQ <- Matrix(NA, nrow=2, ncol=2)
    }
  } else { p5 <- p5 + 1
  iQ[p5, ] <- NA
  Q <- Matrix(Q[,,1], sparse = TRUE)
  }
  
  if (k != (dim(Y)[1] + 1)){
    if (fR == 1 && dim(R)[3]>1){
      for (s in 1:dim(R)[1]){
        for(j in 1:NCOL(t(R[s,,k]))){
          if (R[s,j,k] != 0) {
            i <- i + 1
            p6 <- p6 + 1
            X[i,k] <- R[s,j,k]
            iR[p6,1] <- k - 1
            iR[p6,2] <- j - 1
            iR[p6,3] <- s - 1
          }
        }
      }
    } else {fR <- 0
    R <- Matrix(R[,,1], sparse = TRUE)
    iR <- Matrix(NA, nrow=2, ncol=2)
    }
  } else { p6 <- p6 + 1
  iR[p6, ] <- NA
  R <- Matrix(R[,,1], sparse = TRUE)
  }
}
##################################

################### transform dataframe into matrix
X <- as.matrix(X)
iZ <- as.matrix(iZ)
iV <- as.matrix(iV)
iH <- as.matrix(iH)
iT <- as.matrix(iT)
iQ <- as.matrix(iQ)
iR <- as.matrix(iR)

#benchmark
timeViaIndici <- microbenchmark::microbenchmark(
  kfs  <- KFS(mod2, filtering = "state", smoothing = c("state", "disturbance"), simplify = FALSE),
  kff  <- kalman(T, V, Z, H, Q, R, a1, P1, X, iZ, iH, iV, iT, iQ, iR, Ytd, c("Smoother","Disturbance")),
  times = 3
)

plot(timeViaIndici)
