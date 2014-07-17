#
# Coursera's Mathematical Methods for Quantitative Finance.
# Bisection method to find the zero of a univariate function.
#

bisection_univariate_method <- function(f, a, b, tolerance=10e-3)
{
  if(sign(a) == sign(b))
  {
    stop("The interval does not contain a zero for the given function.")
  }

  if(a >= b)
  {
    stop("The first endpoint of the interval must be less than the second endpoint.")
  }

  if(tolerance <= 0)
  {
    stop("The tolerance must be greater than 0.")
  }

  while(b - a > tolerance)
  {
    c <- (a + b) / 2

    if(sign(f(c)) == sign(f(a)))
    {
      a <- c
    }

    else
    {
      b <- c
    }
  }

  (a + b) / 2
}

bisection_univariate_method_demo <- function()
{
  # We are looking for the square root of the number 2.
  # This is equivalent to look for the zero of the following function:
  # f(x) = x^2 - 2.
  f <- function(x)
  {
    x^2 - 2
  }

  # First derivative of the function f(x).
  f_derivative <- function(x)
  {
    2 * x
  }

  # Calls the univariate bisection method.
  bisection_univariate_method(f, 0, 2)
}
