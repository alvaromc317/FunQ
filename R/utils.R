quantile_function_delete = function(x, quantile.value=0.5)
{
  # Quantile check function
  return(0.5 * base::abs(x) + (quantile.value - 0.5) * x)
}
