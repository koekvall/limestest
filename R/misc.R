# get the row and column index from the vectorization of the
# lower triangular part of an n x n matrix
get_row_col_ltri <- function(idx, n)
{
  last_idx <- as.integer(n * (n + 1) / 2)
  stopifnot(idx <= last_idx)
  # Count backwards from the end to our idx
  back_idx <- last_idx - idx + 1
  # Column counting from the right
  back_col <- ceiling(0.5 * (-1 + sqrt(1 + 8 * back_idx)))
  # Column counting from the left
  col <- n - back_col + 1

  # Get the row
  row <- 0.5 * back_col * (back_col + 1)
  row <- (n - back_col) + (row - back_idx) + 1
  c(row, col)
}

# Get the index in vectorization of lower triangular
get_idx_ltri <- function(row, col, n)
{
  stopifnot(row >= col && row <= n && col <= n)
  back_col <- n - col + 1
  # Indexing counting backwards from the end
  back_idx <- 0.5 * back_col * (back_col + 1) + col - row
  0.5 * n * (n + 1) - back_idx + 1
}
