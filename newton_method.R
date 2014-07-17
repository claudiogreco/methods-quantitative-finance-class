#
# Coursera's Mathematical Methods for Quantitative Finance.
# Newton methods to find the zero of univariate and multivariate functions.
#

newton_univariate_method <- function(f, f_derivative, x0, max_iterations=100)
{
  if(max_iterations <= 0)
  {
    stop("The number of iterations must be greater than 0.")
  }

  x <- x0

  for(i in 1:max_iterations)
  {
    x <- x - f(x) / f_derivative(x)
  }

  x
}

newton_univariate_method_demo <- function()
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

  # Calls the univariate Newton method.
  newton_univariate_method(f, f_derivative, 10)
}

newton_multivariate_method <- function(f, f_gradient, x0, max_iterations=100)
{
  if(max_iterations <= 0)
  {
    stop("The number of iterations must be greater than 0.")
  }

  x <- x0

  for(i in 1:max_iterations)
  {
    x <- x - solve(f_gradient(x), f(x))
  }

  x
}

newton_multivariate_method_demo <- function()
{
  # We are looking for the maximum of the following function:
  # g(x, y) = 1 - (x - 1)^4 - (y - 1)^4.
  # This is equivalent to look for the zero of the following function:
  # f(x, y) = transpose of [4(x - 1)^3 4(y - 1)^3],
  # which is the gradient of the function g(x, y).
  f <- function(x)
  {
    x1 <- x[1]
    x2 <- x[2]
    c(4 * (x1 - 1)^3, 4 * (x2 - 1)^3)
  }

  # Gradient of the function f(x, y), which we want to solve for f(x, y) = 0.
  f_gradient <- function(x)
  {
    x1 <- x[1]
    x2 <- x[2]
    diag(c(12 * (x1 - 1)^2, 12 * (x2 - 1)^2))
  }

  # Calls the multivariate Newton method.
  newton_multivariate_method(f, f_gradient, c(0, 0))
}
