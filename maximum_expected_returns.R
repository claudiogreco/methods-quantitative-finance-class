#
# Coursera's Mathematical Methods for Quantitative Finance.
# Solution of the Maximum Expected Returns optimization problem
# using Lagrangian multipliers plus the Newton method.
#

script.dir <- dirname(sys.frame(1)$ofile)

source(paste(script.dir, "/newton_method.R", sep=""))

maximum_expected_returns_optimization <- function(x, mu, Sigma, sigmaP2)
{
  num_iterations <- 25
  
  # We want to maximize expected return subject to constraint of risk.
  # This can be solved by finding the critical point of the Lagrangian
  # relative to this optimization problem. Finding the critical point
  # of the Lagrangian can be done by looking for the point in which
  # its gradient is equal to zero. For this reason, we define
  # the gradient of the Lagrangian G(x, mu, Sigma, sigmaP2),
  # which we want to solve for G(x, mu, Sigma, sigmaP2) = 0
  # and the gradient of the gradient of the Lagrangian, or
  # the gradient of G, which we need to use the Newton method.
  G <- function(x, mu, Sigma, sigmaP2)
  {
    n <- length(mu)
    
    c(mu + rep(x[n+1], n) + 2 * x[n+2] * (Sigma %*% x[1:n]),
      sum(x[1:n]) - 1, t(x[1:n]) %*% Sigma %*% x[1:n] - sigmaP2)
  }
  
  # Gradient of the gradient of the Lagrangian G(x, mu, Sigma, sigmaP2),
  # or gradient of G, which we need to use the Newton method.
  DG <- function(x, mu, Sigma, sigmaP2)
  {
    n <- length(mu)
    grad <- matrix(0.0, n + 2, n + 2)
    grad[1:n, 1:n] <- 2 * x[n+2] * Sigma
    grad[1:n, n+1] <- 1
    grad[1:n, n+2] <- 2 * (Sigma %*% x[1:n])
    grad[n+1, 1:n] <- 1
    grad[n+2, 1:n] <- 2 * t(x[1:n]) %*% Sigma
    grad
  }
  
  # Executes the Newton method update iterations.
  for(i in 1:num_iterations)
  {
    x <- x - solve(DG(x, mu, Sigma, sigmaP2), G(x, mu, Sigma, sigmaP2))
  }

  # Computes the eigenvalues of the upper-left n x n block
  # of the matrix evaluated by DG(x, mu, Sigma, sigmaP2).
  eigenvalues <- eigen(DG(x, mu, Sigma, sigmaP2)[1:5, 1:5])$values

  is_maximum <- TRUE

  # Checks the second order condition to be sure that
  # the computed solution is a constrained maximum
  # of the given optimization problem.
  for(i in eigenvalues)
  {
    if(i > 0)
    {
      is_maximum <- FALSE 
    }
  }

  # Returns the zero of the function G(x, mu, Sigma, sigmaP2)
  # and the flag which specifies if the computed solution
  # is a constrained maximum of the given optimization problem.
  c(x, is_maximum)
}

maximum_expected_returns_optimization_demo <- function()
{
  # Portfolio weights.
  x <- c(rep(0.5, 5), 1, 1)

  # Asset expected returns.
  mu <- c(0.08, 0.10, 0.13, 0.15, 0.20)

  # Asset returns covariance matrix.
  Sigma <- matrix(c(0.019600, -0.007560, 0.012880, 0.008750, -0.009800,
                    -0.007560, 0.032400, -0.004140, -0.009000, 0.009450,
                    0.012880, -0.004140, 0.052900, 0.020125, 0.020125,
                    0.008750, -0.009000, 0.020125, 0.062500, -0.013125,
                    -0.009800, 0.009450, 0.020125, -0.013125, 0.122500),
                  nrow=5,
                  ncol=5)

  # Portfolio target risk.
  sigmaP2 <- 0.25^2

  # Calls the Maximum Expected Returns optimization method.
  result <- maximum_expected_returns_optimization(x, mu, Sigma, sigmaP2)

  if(result[8] == TRUE)
  {
    print("The solution is a constrained maximum.")
  }

  else
  {
    print("The solution is not a constrained maximum.")
  }

  # Returns the result.
  result[1:7]
}
