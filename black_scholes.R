#
# Coursera's Mathematical Methods for Quantitative Finance.
# Black-Scholes formulas for the price of European call and put options.
#

black_scholes_call_option <- function(S, T, t, K, r, s, q)
{
  d1 <- (log(S/K) + (r - q + 0.5 * s^2) * (T - t)) / (s * sqrt(T - t))
  d2 <- d1 - s * sqrt(T - t)
  S * exp(-q * (T - t)) * pnorm(d1) - K * exp(-r * (T - t)) * pnorm(d2)
}

black_scholes_put_option <- function(S, T, t, K, r, s, q)
{
  d1 <- (log(S/K) + (r - q + 0.5 * s^2) * (T - t)) / (s * sqrt(T - t))
  d2 <- d1 - s * sqrt(T - t)
  K * exp(-r * (T - t)) * pnorm(-d2) - S * pnorm(-d1)
}
