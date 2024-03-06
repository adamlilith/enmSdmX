# Author: Robert J. Hijmans, r.hijmans@gmail.com
# Date: 2009-2021
# Version 0.1
# Licence GPL v3

# NB These functions copied directly from predicts package because cannot import/depends/etc. on predicts::predict() function. Would love another solution so I can just Imports on predict().

methods::setGeneric(name = "predictME", def = function(object, ...) standardGeneric("predictME"))

.maxent_predict <- function(object, mxe, args, x) {
	lambdas <- paste(object@lambdas, collapse='\n')
	variables <- colnames(object@presence)
	x <- x[,variables,drop=FALSE]
	if (inherits(x, "data.frame")) {
		for (i in 1:ncol(x)) {
			if (inherits(x[,i], "factor")) {
				x[,i] <- as.numeric(as.character(x[,i]))
			} else if (inherits(x[,i], "character")) {
				x[,i] <- as.numeric(x[,i])
			}
		}
	} else {
		x[] <- as.numeric(x)
	}
	
	out <- rep(NA, times=nrow(x))
	ok <- rowSums(is.na(x)) == 0
	if (sum(ok) > 0) {
		x <- as.matrix(x[ok, ,drop=FALSE])
		p <- rJava::.jcall(mxe, "[D", "predict", lambdas, rJava::.jarray(colnames(x)), rJava::.jarray(x, dispatch=TRUE), args) 
		p[p == -9999] <- NA
		out[ok] <- p
	}
	out
}


setMethod("predictME", signature(object="MaxEnt_model"), 
	function(object, x, ext=NULL, args="", filename="", ...) {

		stopifnot(predicts::MaxEnt(silent=TRUE))

		args <- c(args, "")
		lambdas <- paste(object@lambdas, collapse="\n")
		variables <- colnames(object@presence)
			
		mxe <- rJava::.jnew("mebridge") 		
		args <- c("-z", args)
		tst <- rJava::.jcall(mxe, "S", "testPredictArgs", lambdas, args) 
		if (!is.null(tst)) {
			stop("args not understood:\n", tst)
		}

		if (!inherits(x, "SpatRaster")) {
			if (! all(variables %in% colnames(x))) {
				stop("missing predictor variables in x")
			}
			me <- .maxent_predict(object, mxe, args, x)
			return(me)
		}

		filename <- trimws(filename)
			
		if (! all(variables  %in%  names(x) )) {
			stop("missing layers (or wrong names)")
		}
		if (terra::nlyr(x) > length(variables)) {
			x <- x[[variables]]
		}
		if (!is.null(ext)) {
			x <- terra::crop(x, ext)
		}
		out <- terra::rast(x, nlyr=1)
		names(out)  <- "maxent"
		ncols <- terra::ncol(out)
		terra::readStart(x)
		on.exit(terra::readStop(x))
		b <- terra::writeStart(out, filename, ...)
		for (i in 1:b$n) {
			rowvals <- terra::readValues(x, b$row[i], b$nrows[i], 1, ncol(x), TRUE, FALSE)
			p <- .maxent_predict(object, mxe, args, rowvals)
			terra::writeValues(out, p, b$row[i], b$nrows[i])
		}
		terra::writeStop(out)
		return(out)
	}
)
