#
# Coursera's Mathematical Methods for Quantitative Finance.
# Solution of the Implied Volatility problem using Newton method.
# The method looks for the zero of the Black-Scholes equation
# for the price of a European call option equation using:
# - The bisection method
# - The Newton method

script.dir <- dirname(sys.frame(1)$ofile)

source(paste(script.dir, "/bisection_method.R", sep=""))
source(paste(script.dir, "/newton_method.R", sep=""))
source(paste(script.dir, "/black_scholes.R", sep=""))

implied_volatility_bisection_solver <- function(S, T, t, K, r, q, option_sold_price)
{
  # Black-Scholes equation respect to sigma.
  fsigma <- function(sigma)
  {
    black_scholes_call_option(S, T, t, K, r, sigma, q) - option_sold_price
  }

  # Calls the Bisection univariate method.
  bisection_univariate_method(fsigma, 0.1, 0.3)
}

implied_volatility_newton_solver <- function(S, T, t, K, r, q, option_sold_price)
{
  # Black-Scholes equation respect to sigma.
  fsigma <- function(sigma)
  {
    black_scholes_call_option(S, T, t, K, r, sigma, q) - option_sold_price
  }

  # Vega is usually referred as the derivative of the
  # Black-Scholes equation respect to sigma.
  vega <- function(sigma)
  {
    d_plus <- (log(S/K) + (r - q + 0.5 * sigma^2) * (T - t)) / (sigma * sqrt(T - t))
    S * sqrt(T - t) * exp(-q * (T - t)) * pnorm(d_plus)
  }

  # Calls the Newton univariate method.
  newton_univariate_method(fsigma, vega, 15)
}
