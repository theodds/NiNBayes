FitDPMN <- function(Y, R, K = NULL,
                    n.adapt = 1000, n.save = 5000, thin = 1, n.chains = 1) {



  ## Fits the nonparametric model that does not share across
  ## treatment groups. Works by stratifying; we do not pass in
  ## treatment indicators.
  
  get.R.1m <- function(R) {
    stopifnot(is.matrix(R))
    J <- dim(R)[2]
    R.1m <- 1 - R
    for(j in J:2) R.1m[,j] <- ifelse(R[,j-1] == 1, R.1m[,j], NA)
    return(R.1m)
  }

  load.module("glm")
  load.module("dic")
  J <- dim(Y)[2]
  N <- dim(Y)[1]
  n.iter <- n.save * thin
  if(is.null(K)) {
    K <- 20
    warning("\nUsing Default K = 20\n")
  }
  Y[R == 0] <- NA
  R.1m <- get.R.1m(R)
  datum <- list(
    K = K,
    Y = Y,
    J = J, N = N,
    ## R_1M = R.1m,
    R  = R,
    gtt = get.c.vars(Y)
  )
  
  
  collect <- c("Y.new", "R.new", "mass", "a", "phi",
               "tau", "xi", "zeta", "gamma", "C")

  out <- run.jags(model = noshare_str, monitor = collect, data = datum,
                  adapt = n.adapt, burnin = n.adapt, sample = n.save, n.chains = n.chains)
  
  return(out)
}

get.c.vars <- function(Y) {
  
  J <- dim(Y)[2]
  ## Get the parameters used in the GARP prior
  
  s <- prelim.norm(Y)
  e <- em.norm(s)
  params <- getparam.norm(s,e)
  
  G <- params$sigma
  gtt <- diag(G)
  gtt.tm1 <- gtt
  for(j in 2:J) {
    gtt.tm1[j] <- gtt[j] - G[j,1:(j-1),drop=FALSE] %*%
      solve(G[1:(j-1),1:(j-1)]) %*% G[1:(j-1),j,drop=FALSE]
  }
  return(gtt.tm1)
}


noshare_str <-
  "

model {

## Hyperparameters
tau.mass <- 1 / 4
tau.mu.j <- 1 / 5 / 5
tau.mu.phi <- 1 / 2.5 / 2.5


mass <- 1
for(k in 1:K) {
mass.vec[k] <- mass / K
}
xi ~ ddirich(mass.vec)


# Prior for a
for(j in 1:J) {
for(k in 1:K) {
a[k,j] ~ dnorm(mu[j], tau.a[j])
}
}


## Prior on phi[k,j,l]
for(k in 1:K) {
for(j in 1:J) {
for(l in 1:J) {
mm[k,j,l] <- ifelse(l==(j-1), mu.phi, 0)
phi[k,j,l] ~ dnorm(mm[k,j,l], tau.phi[j])
}
}
}




lambda ~ dgamma(4,4)
shape ~ dgamma(5,1)
for(j in 1:J) {
muu[j] <- lambda / gtt[j] * 2
rate[j] <- shape / muu[j]
}
for(k in 1:K) {
tau[k,1] ~ dgamma(shape, rate[1])
for(j in 2:J) {
tau[k,j] ~ dgamma(shape + 1, rate[j])
}
}

## Prior for mu[j]
for(j in 1:J) {
mu[j] ~ dnorm(0, .1)
}
## Prior for sigma.a[j] and tau.a[j]
shape.a ~ dgamma(J,1)
lambda.a ~ dgamma(4, 4)
for(j in 1:J) {
sigma.a.raw[j] ~ dgamma(shape.a,1)
sigma.a[j] <- sigma.a.raw[j] * sqrt(gtt[j] / 2) / shape.a * lambda.a
tau.a[j] <- 1 / sigma.a[j] / sigma.a[j]
}
## Prior for mu.phi
mu.phi ~ dt(0, tau.mu.phi, 1)
## Prior for tau.phi
scale.phi ~ dgamma(1,2 / gtt[1])
shape.phi ~ dgamma(J, 1)
for(j in 1:J) {
sigma.phi.raw[j] ~ dgamma(shape.phi, 1)
sigma.phi[j] <- sigma.phi.raw[j] * scale.phi / shape.phi
tau.phi[j] <- 1/pow(sigma.phi[j],2)
}
## Prior for psi
psi ~ dt(0, tau.mu.j, 1)
## Prior for omega.inv
omega.inv ~ dgamma(0.1, 0.1)

## Prior on zeta and gamma
for(k in 1:K) {
for(j in 1:J) {
gamma[k,j] ~ dnorm(mu.gamma, tau.gamma)
zeta[k,j] ~ dnorm(mu.zeta, tau.zeta)
}
}
## And associated hypers
mu.zeta ~ dt(0, tau.mu.j, 1)
tau.zeta ~ dgamma(0.1, 0.1)
mu.gamma ~ dt(0, tau.mu.phi, 1)
tau.gamma ~ dgamma(0.1, 0.1)

## Finally, model the data
for(n in 1:N) {
C[n] ~ dcat(xi[])
Y[n,1] ~ dnorm(a[C[n],1], tau[C[n],1])
for(j in 2:J) {
m[n,j] <- a[C[n], j] + inprod(phi[C[n],j,1:(j-1)], Y[n,1:(j-1)])
Y[n,j] ~ dnorm(m[n,j], tau[C[n],j])
}
## LEAVE OFF WITH DROPOUT SITUATION
for(j in 1:(J-1)) {
prob[n,j] <- (1 - ilogit(zeta[C[n], j] 
+ gamma[C[n],j] * Y[n,j])) *
R[n,j]
R[n,j+1] ~ dbern(prob[n,j])
}

}

## Predict New
C.new ~ dcat(xi[])
tmp <- 1
R.new[1] ~ dbern(tmp)
Y.new[1] ~ dnorm(a[C.new,1], tau[C.new, 1])
for(j in 2:J) {
m.new[j] <- a[C.new, j] + inprod(phi[C.new,j,1:(j-1)], Y.new[1:(j-1)])
Y.new[j] ~ dnorm(m.new[j], tau[C.new, j])
}
for(j in 1:(J-1)) {
prob.new[j] <- (1 - ilogit(zeta[C.new, j] +
gamma[C.new,j] * Y.new[j])) * R.new[j]
R.new[j+1] ~ dbern(prob.new[j])
}
}


"

GCompDPMN <- function(mcmc, Y, K, sens.param) {

  noshare.impute <- function(a, sigma, phi,
                             zeta, gamma, xi, sens.param,
                             num.sim = 10000) {

    ## Wrapper for the compiled code to calculate the MNAR means when
    ## NOT sharing the baseline distribution
    
    J <- dim(a)[2]
    K <- dim(a)[1]
    num.iter <- dim(a)[3]

    if(missing(sens.param))
      sens.param <- array(0, c(num.iter, J))

    .Call("noshare_impute", a, sigma, phi, zeta, gamma, xi, num.sim,
          num.iter, J, K, sens.param)
  }

  get.means <- function(cs, num.save, sens.param, K, J) {

    ## USES THE PACKAGE to calculate the means. Does this by
    ## extracting the relevant paramters calling the compiled code.

    ## Returns matrix of mnar means
    
    a <- jags.extract.2(cs[[1]], "a", c(K,J), num.save)
    num.iter <- dim(a)[3]
    sigma <- 1 / sqrt(jags.extract.2(cs[[1]], "tau", c(K,J), num.save))
    phi <- jags.extract.2(cs[[1]], "phi", c(K,J,J), num.save)
    zeta <- jags.extract.2(cs[[1]], "zeta", c(K,J), num.save)
    gamma <- jags.extract.2(cs[[1]], "gamma", c(K,J), num.save)
    xi <- jags.extract.2(cs[[1]], "xi", K, num.save)
    if(missing(sens.param))
      sens.param <- array(0, c(num.iter, J))
    
    noshare.impute(a,sigma,phi,zeta,gamma,xi,sens.param)
  }

  jags.extract.2 <- function(js, str, d, n.iter) {

    ## Same as jags.extract but BETTER because it puts the iteration
    ## in the last component instead of the first, i.e. a[i,j,k]
    ## returns a[i,j] for iteration k
    
    n.dim <- length(d)
    if(n.dim==0) {
      return(js[,str])
    }
    out <- array(0, c(d, n.iter))
    if(n.dim==1) {
      for(j in 1:d) {
        str2 <- paste(str, "[", j, "]", sep = '')
        out[j,] <- js[,str2]
      }
    }
    if(n.dim==2) {
      for(j1 in 1:d[1]) {
        for(j2 in 1:d[2]) {
          str2 <- paste(str, "[", j1, ",", j2, "]", sep = '')
          out[j1,j2,] <- js[,str2]
        }
      }
    }
    if(n.dim==3) {
      for(j1 in 1:d[1]) {
        for(j2 in 1:d[2]) {
          for(j3 in 1:d[3]) {
            str2 <- paste(str, "[", j1, ",", j2, ",", j3, "]", sep = '')
            out[j1,j2,j3,] <- js[,str2]
          }
        }
      }
    }
    return(out)
  }

  cs <- mcmc
  num.save <- nrow(cs[[1]])
  J <- ncol(Y)

  return(get.means(cs, num.save, sens.param, K, J))


}

