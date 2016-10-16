########################################################################
# General Purpose Param Transformations for Constrained Optimization
#
# Example:
#
#   # Set up bounds on the parameter values
#   upper     = c( Inf, Inf, Inf, Inf, Inf,Inf,Inf, Inf, Inf)
#   lower     = c(-Inf,-Inf,-Inf,-Inf,-Inf,  0,  0,-Inf,-Inf)
#   bd        = parse_bounds(lower,upper)
#
#   # Set up shorter function that encodes upper and lower bounds
#   paratransf_short = function(para, to_unconstr)
#                       paratransf(para, lower, upper, bd, to_unconstr)
#
#   # Convert initial x0 to unconstrained domain
#   tx0 <- paratransf_short(x0,TRUE)
#
#   # Convert transformed init params back again
#   x0 <- paratransf_short(tx0,FALSE)
#
########################################################################

  # Given upper and lower bound vectors (u,l) for parameters, return
  # bool arrays that say which elements bounded above, below, or both
  # Entries Inf or -Inf in (l,u) denote NOT bounded below or above resp
  parse_bounds <- function(l,u) {
    both  = !is.infinite(u) & !is.infinite(l)
    below = !both & !is.infinite(l)
    above = !both & !is.infinite(u)

    return(data.frame(both,above,below))
  }

  # Fcns to convert y in (-Inf,Inf) to x in [l,u]. Args are vectors
  to_constr_both  <- function(y,l,u) 0.5*(l+u) + 0.5*((u-l)*y) / sqrt(1+y^2)
  to_constr_below <- function(y,l)   l + exp(y)
  to_constr_above <- function(y,u)   u - exp(y)

  # Fcns to convert x in [l,u] to y in [-Inf,Inf]. Args are vectors
  to_unconstr_both  <- function(x,l,u) (2*(x-(l+u)/2) / (u-l)) / sqrt(1-(2*(x-(l+u)/2) / (u-l))^2)
  to_unconstr_below <- function(x,l)   log(x-l)
  to_unconstr_above <- function(x,u)   log(u-x)


  # Transform parameters to unconstrained domain or constrained domain
  paratransf <- function(para, l, u, bd, to_unconstr) {
    tpara <- para

    if (to_unconstr) {
      tpara[bd$both]  <- to_unconstr_both( para[bd$both],  l[bd$both], u[bd$both])
      tpara[bd$below] <- to_unconstr_below(para[bd$below], l[bd$below])
      tpara[bd$above] <- to_unconstr_above(para[bd$above], u[bd$above])
    }
    else {
      tpara[bd$both]  <- to_constr_both( para[bd$both],  l[bd$both], u[bd$both])
      tpara[bd$below] <- to_constr_below(para[bd$below], l[bd$below])
      tpara[bd$above] <- to_constr_above(para[bd$above], u[bd$above])
    }

    return(tpara)
  }


########################################################################
## Setting Initial Parameter Values for Question 2 #####################
########################################################################

  # Function to compute AR coefficient
  ar1coeff <- function(x) cor(x[2:length(x)], x[1:(length(x)-1)])

  # Put data in dataframe with lags
  regdata = data.frame(ydata[2:Tfull,], ydata[1:(Tfull-1),])

  # Estimate AR(1) for GDP and UNRATE
  regGDP    <- lm(GDP    ~ GDP.1,    data=regdata)
  regUNRATE <- lm(UNRATE ~ UNRATE.1, data=regdata)

  # Take resids from AR(1) res above, do PCA on them to extract common
  # component, i.e. the *common* source of variation that is unexplained
  # by own lags in each regression.
  #
  # Note: PCA done after *standardizing* residuals to have same variance
  # so one variable's residuals does not dominate
  X <- cbind(regGDP$residuals/sd(regGDP$residuals),
             regUNRATE$residuals/sd(regUNRATE$residuals))
  eigdecomp <- eigen(t(X) %*% X, TRUE)
  pcs       <- X %*% eigdecomp$vectors

  # Proxy z process with first PC. Standardize this process so that the
  # variance of the shock in an AR(1) estimated on this process has
  # variance=1, as we will impose in the estimated state space model
  z         <- pcs[,1]                    # z is first PC
  zcorr     <- ar1coeff(z)                # Get AR(1) coeff
  z         <- (z/sd(z))/sqrt(1-zcorr^2)  # Standardize
  regdata$z  <- z[2:Tfull]                # Add to data matrix

  # Estimate naive proxy of state space model
  reg2GDP    <- lm(GDP    ~ GDP.1    + z, data=regdata)
  reg2UNRATE <- lm(UNRATE ~ UNRATE.1 + z, data=regdata)


########################################################################
## Question 1 ##########################################################
########################################################################

# Estimate VAR
Q1var   = mgnldnsty(ydata, Nlags, vprior=list(sig=c(0.005,0.005),w=1))

# Part a: Forecast eight quarters ahead
Q1fcast = fcast(ydata[(Tfull-Nlags+1):Tfull,],
                Q1var$var$By, Q1var$var$Bx,
                xdata=NULL, const=TRUE, NQfcast)

# Part b, c: Estimated residuals and covariance matrix of innovations
Q1resids <- Q1var$var$u[1:(Tfull-Nlags),]
Q1covmat <- print(crossprod(Q1resids)/dim(Q1resids)[1])

# Impulse responses: Nvar x Nshocks x Nstep array
Q1irfs = impulsdtrf(Q1var$var)


########################################################################
## Question 2 ##########################################################
########################################################################

# Initialize Parameter values, set using naive estimates above
rho      = c(reg2GDP$coefficients["GDP.1"],
              reg2UNRATE$coefficients["UNRATE.1"])
gamma    = zcorr
alpha    = c(reg2GDP$coefficients["z"],
            reg2UNRATE$coefficients["z"])
stdy     = c(sd(reg2GDP$residuals),
             sd(reg2UNRATE$residuals))
delta    = c(reg2GDP$coefficients["(Intercept)"],
             reg2UNRATE$coefficients["(Intercept)"])
params0  = c(rho, gamma, alpha, stdy, delta)


# Get parameter estimates and state estimates given those params
# Note: Must transform parameters back to constrained domain, since
# csminwel estimates are on the transformed unconstrained domain
params   <- paratransf_short(minimized$xh, FALSE)
filtered <- likelihood(params, TRUE)
filtered$yactual <- ydata[2:Tfull,]   # Attach sample

# Construct system matrices from estimated parameters
rho   <- diag(params[1:2])
gamma <- params[3]
alph  <- params[4:5]
stdy  <- params[6:7]
delta <- params[8:9]
G <- diag(c(gamma, 1)) # State Transition matrix
M <- diag(c(1, 0))      # Standard deviations in state transition equation
N <- diag(stdy)         # Standard deviations in measurement equation
H <- matrix(c(alph, delta), nrow = 2) # State to observable mapping matrix


# Forecast next Nquarters using estimated model
Tkf            <- Tfull-1
Q2fcast_states <- rbind(filtered$shat[(Tkf-1):Tkf,],    matrix(0,NQfcast,2))
Q2fcast        <- rbind(filtered$yactual[(Tkf-1):Tkf,], matrix(0,NQfcast,2))
for (n in 2:(NQfcast+1)) {
  spred <- G %*% t(t(Q2fcast_states[n,]))
  yprev <- t(t(Q2fcast[n,]))

  Q2fcast_states[n+1,] <- spred
  Q2fcast[n+1,]        <- rho %*% yprev + H %*% spred
}

# Calculate covariance matrix of innovations
Q2covmat <- G %*% filtered$sig[Tkf,,] %*% t(G) + N %*% N

# Estimate forecast errors
filtered$ypred  <- matrix(0, nrow(filtered$yactual), 2)  # Predicted y
filtered$errors <- matrix(0, nrow(filtered$yactual), 2)
for (n in 1:(Tkf-1)) {
  # Filtered state at time n
  s <- t(t(filtered$shat[n,]))

  # Predicted state and observable for time n+1
  spred <- (G %*% s)
  ypred <- rho %*% t(t(filtered$yactual[n,])) + H %*% spred

  # Store predicted y at n+1 and compute forecast error
  filtered$ypred[n+1,]  <- ypred
  filtered$errors[n+1,] <- filtered$yactual[n+1,] - ypred
}


# Calculate the IRFs

  # Initial Conditions
  #y0 <- filtered$yactual[Tkf,]
  #s0 <- filtered$shat[Tkf,]
  y0 <- c(0,0)
  s0 <- c(0,1)

  # Shocks to hit system with: GDP, UNRATE, z, nothing (for constant)
  shocks = diag(c(stdy,1,0))

  # Initialize everything
  Nstate  <- length(s0)
  Nshocks <- nrow(shocks)
  Nsteps  <- dim(Q1irfs)[3]
  Q2baseline_s <- array(0, dim=c(Nstate, Nshocks, Nsteps)) # Evolution of no shock baseline
  Q2baseline_y <- array(0, dim=c(Nvar,   Nshocks, Nsteps)) # Evolution of no shock baseline
  Q2shocked_s  <- array(0, dim=c(Nstate, Nshocks, Nsteps)) # Evolution of shocked system
  Q2shocked_y  <- array(0, dim=c(Nvar,   Nshocks, Nsteps)) # Evolution of shocked system
  Q2irfs       <- array(0, dim=c(Nvar,   Nshocks, Nsteps))
  colnames(Q2irfs) <- c("GDP", "UNRATE", "z", "Constant")
  rownames(Q2irfs) <- c("GDP", "UNRATE")
  for (sh in 1:Nshocks) {
    Q2baseline_y[,sh,1] <- y0
    Q2baseline_s[,sh,1] <- s0
    Q2shocked_y[,sh,1]  <- y0+shocks[sh,1:2]
    Q2shocked_s[,sh,1]  <- s0+shocks[sh,3:4]
    Q2irfs[,sh,1]       <- Q2shocked_y[,sh,1] - Q2baseline_y[,sh,1]

    for (n in 1:(Nsteps-1)) {

      # Predict next-period s in shocked and baseline case
      spred_baseline <- (G %*% t(t(Q2baseline_s[,sh,n])))
      spred_shocked  <- (G %*% t(t(Q2shocked_s[,sh,n])))

      # Grab current y in shocked and baseline case
      y_baseline  <- t(t(Q2baseline_y[,sh,n]))
      y_shocked   <- t(t(Q2shocked_y[,sh,n]))

      # Next-period state and observable in baseline case
      Q2baseline_s[,sh,n+1] <- spred_baseline
      Q2baseline_y[,sh,n+1] <- (rho %*% y_baseline) + (H %*% spred_baseline)

      # Next-period state and observable in shocked case
      Q2shocked_s[,sh,n+1] <- spred_shocked
      Q2shocked_y[,sh,n+1] <- (rho %*% y_shocked) + (H %*% spred_shocked)

      # Compute IRF
      Q2irfs[,sh,n+1] <- Q2shocked_y[,sh,n+1] - Q2baseline_y[,sh,n+1]
    }
  }



########################################################################
## Generate Plots ######################################################
########################################################################

  tsQ1fcast <- ts(rbind(ydata[(Tfull-33):Tfull,], Q1fcast[3:10,]), frequency=4, start=c(2008,1))
  tsQ2fcast <- ts(rbind(ydata[(Tfull-33):Tfull,], Q2fcast[3:10,]), frequency=4, start=c(2008,1))

  # Functions for plotting and saving
  lastdata <- function() abline(v=c(2016.25), lty=2)
  myplot <- function(name, x, ylab, lastdata=FALSE) {
    pdf(name)
    ts.plot(x, ylab=ylab, lwd=2)
    if (lastdata) {
      lastdata()
    }
    dev.off()
  }
  myqq <- function(name, x, Tbreak) {
    N = length(x)
    d <- qqnorm(x, plot.it=FALSE)
    pch   = rep(4, N)
    pch[Tbreak:N] = 1
    pdf(name)
    plot(d$x, d$y, pch=pch, ylab="Theoretical Quantiles", xlab="Sample Quantiles")
    legend(-2, 3, c("Early Sample", "Later Sample"), pch=c(4,1))
    dev.off()
  }

  # Plot VAR forecasts
  myplot("VAR_GDP.pdf", tsQ1fcast[,1], "Log Real GDP", lastdata=TRUE)
  myplot("VAR_UNRATE.pdf", tsQ1fcast[,2], "Unemployment Rate", lastdata=TRUE)

  # Plot factor model forecasts
  myplot("Factor_GDP.pdf", tsQ2fcast[,1], "Log Real GDP", lastdata=TRUE)
  myplot("Factor_UNRATE.pdf", tsQ2fcast[,2], "Unemployment Rate", lastdata=TRUE)

  # Plot growth forecasts
  myplot("VAR_GDPgrowth.pdf", 100*diff(tsQ1fcast[,1]), "Log Real GDP", lastdata=TRUE)
  myplot("Factor_GDPgrowth.pdf", 100*diff(tsQ2fcast[,1]), "Log Real GDP", lastdata=TRUE)

  # Covariance matrix of innovations
  print(Q1covmat)
  print(Q2covmat)


  # Plot standardized forecast errors
  tskfsample   <- function(x) ts(x, frequency=4, start=c(1948,2))
  tsfullsample <- function(x) ts(x, frequency=4, start=c(1948,3))
  Q1resids_GDP    <- tsfullsample(Q1resids[,1]/sqrt(Q1covmat[1,1]))
  Q1resids_UNRATE <- tsfullsample(Q1resids[,2]/sqrt(Q1covmat[2,2]))
  Q2resids_GDP    <- tskfsample(filtered$fcsterr[,1]/sqrt(Q2covmat[1,1]))
  Q2resids_UNRATE <- tskfsample(filtered$fcsterr[,2]/sd(filtered$fcsterr[,2]))

  myplot("VAR_resids_GDP.pdf", Q1resids_GDP, "Standardized Log GDP Residuals")
  myplot("VAR_resids_UNRATE.pdf", Q1resids_UNRATE, "Standardized UNRATE Residuals")
  myplot("Factor_resids_GDP.pdf", Q2resids_GDP, "Standardized Log GDP Residuals")
  myplot("Factor_resids_UNRATE.pdf", Q2resids_UNRATE, "Standardized UNRATE Residuals")

  Tbreak = 146
  myqq("VAR_qqresids_GDP.pdf", Q1resids_GDP, Tbreak)
  myqq("VAR_qqresids_UNRATE.pdf", Q1resids_UNRATE, Tbreak)
  myqq("Factor_qqresids_GDP.pdf", Q2resids_GDP, Tbreak)
  myqq("Factor_qqresids_UNRATE.pdf", Q2resids_UNRATE, Tbreak)


  # Plot the IRFs
  plotir(Q1irfs, file="Q1irfs.pdf")
  plotir(Q2irfs, shock=c(1,2,3), file="Q2irfs.pdf")
  plotir(Q2irfs, shock=c(2), file="Q2irfsb.pdf")



