#' Calculate weights for a model
#'
#' Calculates weighting for a model. Each record receives a numeric weight.
#'
#' @param w Either logical in which case \code{TRUE} (default) causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'}) \emph{or} a numeric list of weights, one per row in \code{data} \emph{or} the name of the column in \code{data} that contains site weights. If \code{FALSE}, then each datum gets a weight of 1.
#' @param Data frame
#' @param resp Name of response column
#'
#' @returns A numeric vector.
#' @keywords internal

.calcWeights <- function(w, data, resp, family) {

	if (inherits(w, 'logical')) {
		if (w & (family %in% c('binomial', 'quasibinomial'))) {
			posCases <- sum(data[ , resp, drop=TRUE] == 1)
			negCases <- sum(data[ , resp, drop=TRUE] == 0)
			w <- c(rep(1, posCases), rep(posCases / negCases, negCases))
		} else {
			w <- rep(1, nrow(data))
		}
	} else if (inherits(w, 'character')) {
		w <- data[ , w, drop=TRUE]
	}
	w <- w / max(w)
	w
	
}

