
# TENSOR DESIGN MATRICES -------------------------------------------------------

#' @title Build Tensor Design Matrix
#' @description Inner function to construct the tensor product design matrix required for functional data modeling. It efficiently scales the B-spline basis evaluations by the subject-specific functional principal component scores across all observed time points.
#' @param scores A matrix of functional Principal Component Analysis (fPCA) scores for each subject.
#' @param Y.mask A logical mask matrix of dimensions (N x T) indicating which time points are observed (`TRUE`) or missing (`FALSE`) for each subject.
#' @param spline.basis The evaluation matrix of the B-spline basis functions at the time points.
#' @param intercept Boolean. If TRUE (default) the design starts with an intercept block built from \code{spline.basis[, -1]}; if FALSE only the score-scaled spline blocks are produced.
#' @return A matrix containing the stacked tensor products of the (optional) intercept bases and the scaled spline bases for all valid observations.
build_tensor_matrix <- function(scores, Y.mask, spline.basis, intercept = TRUE)
{
  # Define dimensions
  n.obs <- nrow(scores)
  npc <- ncol(scores)
  n.basis <- ncol(spline.basis)
  intercept.spline.basis <- if(intercept) spline.basis[, -1, drop = FALSE] else NULL
  n.int <- if(intercept) ncol(intercept.spline.basis) else 0L

  # Identify the specific time-point indices where each subject has valid observations
  idx_list <- lapply(seq_len(n.obs), function(i) which(Y.mask[i, ]))
  n_i <- lengths(idx_list)
  n_rows <- sum(n_i)

  # Pre-allocate tensor matrix
  out <- matrix(0, nrow = n_rows, ncol = n.int + npc * n.basis)

  # Pre-calculate the exact column blocks for each score component
  int_cols <- seq_len(n.int)
  col_idx  <- lapply(seq_len(npc), function(c) n.int + (c - 1L) * n.basis + seq_len(n.basis))

  pos <- 1L
  for (i in seq_len(n.obs))
  {
    # Retrieve the valid time indices and the observation count for the current subject
    idx <- idx_list[[i]]
    ni  <- n_i[i]
    if (ni == 0L) next

    # Determine the exact row range in the pre-allocated matrix where this subject's data belongs
    rows <- seq.int(pos, length.out = ni)

    # Subset the standard spline basis for the specific time points where this subject was observed
    spl <- spline.basis[idx, , drop = FALSE]
    if (intercept) {
      out[rows, int_cols] <- intercept.spline.basis[idx, , drop = FALSE]
    }

    # Scale the subject's splines by their scores and slot them into the correct pre-calculated column blocks
    sc <- scores[i, ]
    for (c in seq_len(npc))
    {
      out[rows, col_idx[[c]]] <- spl * sc[c]
    }

    pos <- pos + ni
  }
  return(out)
}

#' @title Build Sparse Tensor Design Matrix
#' @description Inner function to construct the tensor product design matrix as a sparse \code{SparseM::matrix.csr}, exploiting the local support of the B-spline basis (each row has only a handful of non-zero entries). The column layout is identical to \code{cbind(1, build_tensor_matrix(...))}, i.e. it includes the leading global intercept column, so the result can be passed directly to \code{quantreg::rq.fit.sfn}.
#' @param scores A matrix of functional Principal Component Analysis (fPCA) scores for each subject.
#' @param Y.mask A logical mask matrix of dimensions (N x T) indicating which time points are observed (`TRUE`) or missing (`FALSE`) for each subject.
#' @param spline.basis The evaluation matrix of the B-spline basis functions at the time points.
#' @return A \code{matrix.csr} sparse matrix containing the intercept column, the intercept bases and the scaled spline bases for all valid observations.
#' @importClassesFrom SparseM matrix.csr
build_tensor_matrix_sparse <- function(scores, Y.mask, spline.basis)
{
  n.obs <- nrow(scores)
  npc <- ncol(scores)
  n.basis <- ncol(spline.basis)
  intercept.spline.basis <- spline.basis[, -1, drop = FALSE]
  n.int <- ncol(intercept.spline.basis)
  n.time <- nrow(spline.basis)
  ncols <- 1L + n.int + npc * n.basis

  # Non-zero structure per time point (B-splines are locally supported)
  nz.spl <- lapply(seq_len(n.time), function(t) which(spline.basis[t, ] != 0))
  nz.int <- lapply(seq_len(n.time), function(t) which(intercept.spline.basis[t, ] != 0))
  nzs.len <- lengths(nz.spl)
  nzi.len <- lengths(nz.int)

  idx_list <- lapply(seq_len(n.obs), function(i) which(Y.mask[i, ]))
  n_i <- lengths(idx_list)
  n_rows <- sum(n_i)
  tt <- unlist(idx_list, use.names = FALSE)   # time index per design row
  ii <- rep(seq_len(n.obs), n_i)              # subject index per design row

  # Global intercept entries (one per row)
  one.rows <- seq_len(n_rows)

  # Intercept basis entries
  int.rows <- rep(one.rows, nzi.len[tt])
  int.nz <- unlist(nz.int[tt], use.names = FALSE)
  int.vals <- intercept.spline.basis[cbind(rep(tt, nzi.len[tt]), int.nz)]

  # Spline entries, replicated per component with score scaling
  spl.rows.base <- rep(one.rows, nzs.len[tt])
  spl.nz <- unlist(nz.spl[tt], use.names = FALSE)
  spl.vals.base <- spline.basis[cbind(rep(tt, nzs.len[tt]), spl.nz)]
  spl.subj <- ii[spl.rows.base]

  rows <- c(one.rows, int.rows, rep(spl.rows.base, npc))
  cols <- c(rep(1L, n_rows), 1L + int.nz, unlist(lapply(seq_len(npc), function(c) 1L + n.int + (c - 1L) * n.basis + spl.nz), use.names = FALSE))
  vals <- c(rep(1, n_rows), int.vals, unlist(lapply(seq_len(npc), function(c) spl.vals.base * scores[spl.subj, c]), use.names = FALSE))

  # Assemble in compressed sparse row format (entries sorted by row, then column)
  ord <- order(rows, cols)
  ia <- c(1L, cumsum(tabulate(rows, nbins = n_rows)) + 1L)
  methods::new(
    "matrix.csr",
    ra = vals[ord],
    ja = as.integer(cols[ord]),
    ia = as.integer(ia),
    dimension = as.integer(c(n_rows, ncols)))
}
