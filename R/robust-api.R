#' Principal representations from consecutive moments
#'
#' @param moments E(X), ..., E(X^v), excluding m0=1. Accepts numeric vectors,
#'   decimal strings, or an Rmpfr vector. Strings and MPFR values preserve input
#'   precision; converting rounded doubles cannot recover lost information.
#' @param support Finite increasing support endpoints. Precision input currently
#'   requires `[0,1]`; numerical support transformation is available for doubles.
#' @param tol Absolute moment reproduction tolerance, separate from rank tests.
#' @param precision Working precision in bits; NULL selects the double backend
#'   for numeric inputs and the MPFR backend otherwise. Set explicitly to request
#'   MPFR arithmetic for numeric input. Precision refinement is automatic.
#' @param accuracy Required absolute stability of nodes and weights on refinement.
#' @param max_precision Maximum working precision in bits.
#' @return A `principal_representations` object with lower/upper node-weight data
#'   frames, metadata and diagnostics. MPFR output is retained in
#'   `high_precision` for evaluation; data frames are rounded display values.
#' @details Small raw eigenvalues do not establish exact boundary membership.
#'   Ambiguous long double-input sequences fail explicitly. The MPFR backend
#'   constructs orthogonal-polynomial recurrences and refines Gaussian nodes
#'   and Christoffel weights at increasing precision. Returned numerical results
#'   are precision-checked, not interval-arithmetic certificates. Diagnostics
#'   distinguish input precision and working precision. Arbitrary unrecognised
#'   singular boundary sequences may be reported as unresolved.
#' @importFrom Rmpfr mpfr
#' @export
#' @examples
#' reps <- principal_representations(1 / (2:4))
#' reps$lower
#' # Preserve precision from the start:
#' m <- 1 / Rmpfr::mpfr(2:7, 256)
#' reps <- principal_representations(m)
principal_representations <- function(moments, support = c(0, 1), tol = 1e-12,
                                      precision = NULL, accuracy = 1e-8,
                                      max_precision = 1024L) {
  validate_precision_controls(tol, precision, accuracy, max_precision)
  precise <- is.character(moments) || is_mp(moments) || !is.null(precision)
  if (!precise) {
    result <- principal_double(moments, support, tol)
    simple <- all(moments == moments[1]) ||
      all(abs(moments - moments[1]^seq_along(moments)) <=
        16 * .Machine$double.eps * pmax(abs(moments), 1e-300))
    result$diagnostics$status <- if (!result$boundary) "checked_double" else if (simple) "boundary_at_input_precision" else "numerical_boundary"
    result$diagnostics$input_precision <- 53L
    result$diagnostics$working_precision <- 53L
    return(result)
  }
  if (!identical(as.numeric(support), c(0, 1))) {
    stop("Precision input currently requires support [0,1]; rescale moments before rounding.",
      call. = FALSE
    )
  }
  precision_representations(moments, seq_along(moments), tol, precision, accuracy,
    max_precision,
    consecutive = TRUE
  )
}

#' Principal representations from selected power moments
#'
#' @inheritParams principal_representations
#' @param orders Strictly increasing finite positive real exponents, excluding 0.
#' @param maxit Maximum Newton iterations per correction; total corrections are also
#'   bounded by a budget of four times this value at each working precision.
#' @return An object inheriting from `generalised_principal_representations` and
#'   `principal_representations`, usable by [evaluate_representations()] and
#'   [power_moment_bounds()].
#' @details Support is `[0,1]`. Equally spaced powers reduce to consecutive moments.
#'   Other systems use MPFR Newton corrections with ordered nodes and positive
#'   weights. Warm starts are refined against relative moment equations; on
#'   failure, continuation first deforms a uniform quadrature's exponents and
#'   then its moments to the requested system. No clipping or missing-moment
#'   imputation is performed. Nonconvergence is not a proof of infeasibility.
#'   Numerical precision checks and local input-sensitivity estimates are not
#'   rigorous global uncertainty enclosures. Singular nontrivial boundaries can
#'   be unresolved; point masses and endpoint mixtures have analytic paths.
#' @export
#' @examples
#' reps <- generalised_principal_representations(1 / c(3, 4, 6), c(2, 3, 5))
#' power_moment_bounds(reps, 1)
generalised_principal_representations <- function(moments, orders, tol = 1e-9,
                                                  maxit = 100L, precision = NULL,
                                                  accuracy = 1e-8,
                                                  max_precision = 1024L) {
  validate_precision_controls(tol, precision, accuracy, max_precision)
  valid <- function(x) {
    is.numeric(x) && !is.complex(x) && is.null(dim(x)) &&
      length(x) > 0 && all(is.finite(x))
  }
  if (!valid(orders) || any(orders <= 0) || any(diff(orders) <= 0) ||
    length(moments) != length(orders)) {
    stop("Supply finite moments and equally long, strictly increasing positive orders.", call. = FALSE)
  }
  if (!valid(maxit) || length(maxit) != 1L || maxit < 1 || maxit != floor(maxit)) {
    stop("Invalid maxit.", call. = FALSE)
  }
  arithmetic <- all(orders / orders[1] == seq_along(orders))
  if (arithmetic) {
    result <- principal_representations(moments,
      tol = tol, precision = precision,
      accuracy = accuracy, max_precision = max_precision
    )
    if (is.null(result$high_precision)) {
      transformed <- lapply(result[c("lower", "upper")], function(rule) {
        rule$nodes <- rule$nodes^(1 / orders[1])
        rule
      })
      loss <- vapply(transformed, function(rule) {
        max(abs(drop(outer(orders, rule$nodes, function(k, x) x^k) %*%
          rule$weights) - moments))
      }, numeric(1))
      if (any(loss > tol)) {
        return(generalised_principal_representations(moments, orders, tol,
          maxit,
          precision = 128L, accuracy = accuracy,
          max_precision = max_precision
        ))
      }
      result$lower <- transformed$lower
      result$upper <- transformed$upper
    } else {
      for (side in c("lower", "upper")) {
        if (!is.null(result$precision_reference)) {
          result$precision_reference[[side]]$nodes <- result$precision_reference[[side]]$nodes^
            (1 / mp_value(orders[1], result$diagnostics$working_precision))
        }
        result$high_precision[[side]]$nodes <- result$high_precision[[side]]$nodes^(1 / mp_value(orders[1], result$diagnostics$working_precision))
        result[[side]]$nodes <- mp_num(result$high_precision[[side]]$nodes)
      }
    }
    result$orders <- orders
    if (!is.null(result$precision_reference)) {
      changes <- vapply(c("lower", "upper"), function(side) {
        a <- result$high_precision[[side]]
        b <- result$precision_reference[[side]]
        mp_num(max(abs(a$nodes - b$nodes), abs(a$weights - b$weights)))
      }, numeric(1))
      if (max(changes) > accuracy / 10) {
        stop("Unresolved power transformation: increase precision.", call. = FALSE)
      }
      result$diagnostics$refinement_change <- changes
    }
  } else {
    result <- precision_representations(moments, orders, tol, precision, accuracy,
      max_precision,
      maxit = maxit
    )
  }
  class(result) <- c("generalised_principal_representations", "principal_representations")
  result
}

precision_representations <- function(moments, orders, tol, precision, accuracy,
                                      max_precision, consecutive = FALSE,
                                      maxit = 100L) {
  numeric_input <- is.numeric(moments) && !is_mp(moments)
  if ((!numeric_input && !is.character(moments) && !is_mp(moments)) ||
    is.complex(moments) || !is.null(dim(moments)) || !length(moments)) {
    stop("Supply a non-empty vector of finite moments.", call. = FALSE)
  }
  for (x in list(tol, accuracy, max_precision)) {
    if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
      stop("Invalid tol, accuracy or max_precision.", call. = FALSE)
    }
  }
  if (tol >= 1 || accuracy >= 1 || max_precision < 128 || max_precision != floor(max_precision)) {
    stop("Invalid tol, accuracy or max_precision.", call. = FALSE)
  }
  if (!is.null(precision) && (!is.numeric(precision) || length(precision) != 1L ||
    !is.finite(precision) || precision < 80 || precision != floor(precision) || precision > max_precision)) {
    stop("Invalid precision.", call. = FALSE)
  }
  input_bits <- if (is_mp(moments)) min(Rmpfr::getPrec(moments)) else if (numeric_input) 53L else NA_integer_
  supplied_bits <- if (is_mp(moments)) {
    input_bits
  } else if (is.character(moments)) {
    digits <- nchar(gsub("[^0-9]", "", sub("[eE].*$", "", moments)))
    ceiling(max(digits) * log2(10))
  } else {
    53L
  }
  bits <- if (is.null(precision)) {
    # Avoid needlessly rounding accurate input down to a precision at which
    # its nearly-boundary moment vector may appear infeasible. Reserve room
    # for the independent refinement calculation.
    max(80L, min(floor(max_precision / 2L), max(128L, 4L * length(orders) + 64L, supplied_bits)))
  } else {
    precision
  }
  bits <- min(bits, max_precision)
  previous <- NULL
  warm <- if (!consecutive) {
    tryCatch(suppressWarnings(generalised_double(
      mp_num(mp_value(moments, 128)),
      orders, tol, maxit
    )), error = function(e) NULL)
  } else {
    NULL
  }
  if (numeric_input && !is.null(warm) && isTRUE(warm$boundary)) {
    warm$diagnostics$status <- "boundary_at_input_precision"
    warm$diagnostics$input_precision <- 53L
    warm$diagnostics$working_precision <- 53L
    return(warm)
  }
  last_error <- "precision refinement did not converge"
  repeat {
    m <- tryCatch(mp_value(moments, bits), error = function(e) NULL)
    if (is.null(m) || any(!is.finite(m)) || any(m < 0 | m > 1) || any(diff(m) > 0)) {
      stop("Infeasible or invalid moments on [0,1].", call. = FALSE)
    }
    candidate <- tryCatch(
      {
        # Exact endpoint mixtures and point masses (verified against every power).
        p <- m[1]^(1 / mp_value(orders[1], bits))
        predicted <- p^orders
        eps <- mp_value(2, bits)^(-bits + 12L)
        point <- length(m) > 1L && all(abs(predicted - m) <= eps * Rmpfr::pmax(abs(m), abs(predicted)))
        if (all(m == m[1]) && (length(m) > 1L || m[1] == 0 || m[1] == 1)) {
          rule <- list(nodes = mp_value(c(0, 1), bits), weights = c(1 - m[1], m[1]))
          keep <- rule$weights > 0
          rule <- lapply(rule, function(x) x[keep])
          list(lower = rule, upper = rule)
        } else if (point) {
          rule <- list(nodes = p, weights = mp_value(1, bits))
          list(lower = rule, upper = rule)
        } else if (consecutive) {
          mp_consecutive(c(mp_value(1, bits), m), bits)
        } else {
          mp_selected(m, orders, bits, maxit, warm)
        }
      },
      error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(candidate)) {
      relative_residuals <- numeric()
      residuals <- vapply(candidate, function(rule) {
        f <- methods::as(mp_basis(mp_value(c(0, orders), bits), rule$nodes) %*% rule$weights, "mpfr") - c(mp_value(1, bits), m)
        if (any(rule$weights < 0) || any(rule$nodes < 0 | rule$nodes > 1)) {
          return(Inf)
        }
        target_scale <- abs(c(mp_value(1, bits), m))
        target_scale[target_scale == 0] <- 1
        relative_residuals <<- c(relative_residuals, mp_num(max(abs(f) / target_scale)))
        mp_num(max(abs(f)))
      }, numeric(1))
      if (max(residuals) <= tol && max(relative_residuals) < accuracy / 100 && !is.null(previous)) {
        changes <- vapply(c("lower", "upper"), function(side) {
          a <- previous[[side]]
          b <- candidate[[side]]
          if (length(a$nodes) != length(b$nodes)) {
            return(Inf)
          }
          mp_num(max(abs(a$nodes - b$nodes), abs(a$weights - b$weights)))
        }, numeric(1))
        if (max(changes) < accuracy / 10) {
          result <- list(
            moments = moments, orders = orders, support = c(0, 1), order = length(orders),
            high_precision = candidate, precision_reference = previous,
            boundary = length(candidate$lower$nodes) == length(candidate$upper$nodes) &&
              all(candidate$lower$nodes == candidate$upper$nodes) &&
              all(candidate$lower$weights == candidate$upper$weights)
          )
          for (side in c("lower", "upper")) {
            result[[side]] <- data.frame(
              nodes = mp_num(candidate[[side]]$nodes),
              weights = mp_num(candidate[[side]]$weights)
            )
          }
          sensitivity <- if (numeric_input) mp_sensitivity(candidate, orders, m, bits) else NA_real_
          if (numeric_input && sensitivity > accuracy) {
            stop("Unresolved input precision: rounded double moments cannot support the requested node/weight accuracy. Supply accurate decimal strings or MPFR moments.",
              call. = FALSE
            )
          }
          result$diagnostics <- list(
            status = "precision_checked", residuals = residuals,
            relative_residuals = relative_residuals,
            input_sensitivity = sensitivity,
            working_precision = bits, input_precision = input_bits,
            refinement_change = changes, accuracy = accuracy,
            input_interpretation = "Supplied values are constraints; input uncertainty is not a certified enclosure."
          )
          class(result) <- "principal_representations"
          return(result)
        }
      }
      previous <- candidate
      if (!consecutive) warm <- candidate
    }
    if (bits >= max_precision) break
    bits <- min(max_precision, bits * 2L)
  }
  stop(paste(
    "Generalised quadrature did not converge: unresolved precision calculation:", last_error,
    "Supply more accurate moments or increase max_precision."
  ), call. = FALSE)
}


validate_precision_controls <- function(tol, precision, accuracy, max_precision) {
  for (x in list(tol, accuracy, max_precision)) {
    if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
      stop("Invalid tol, accuracy or max_precision.", call. = FALSE)
    }
  }
  if (tol >= 1 || accuracy >= 1 || max_precision < 128 || max_precision != floor(max_precision)) {
    stop("Invalid tol, accuracy or max_precision.", call. = FALSE)
  }
  if (!is.null(precision) && (!is.numeric(precision) || is.complex(precision) ||
    length(precision) != 1L || !is.finite(precision) || precision < 80 ||
    precision != floor(precision) || precision > max_precision)) {
    stop("Invalid precision.", call. = FALSE)
  }
}
