#' @name troubleshooting_parallel_operations
#'
#' @title Troubleshooting parallel operations
#'
#' @description This is a guide to solving issues with running functions that can use more than one core. This includes the \code{train}\emph{XYZ} functions, \code{\link{bioticVelocity}} function, and the \code{\link{predictEnmSdm}}. Each of these function has the argument \code{cores}. By default, the value of \code{cores} is 1, so the function will use only one core. By setting this higher, you can use more cores on your machine.  However, occasionally you will run into the error:\cr
#' \code{Error in checkForRemoteErrors(lapply(cl, recvResult)) :}\cr
#' \code{  2 nodes produced errors; first error: object '.doSnowGlobals' not found}\cr
#' This means that the worker "nodes" (different instances of \code{R} started by the function to run in parallel) cannot find the \pkg{doParallel} package, even if it is installed on your system.
#'
#' There are several solutions to this issue. One of them may work for you, and none are inherent to \pkg{enmSdmX}, as far as I can tell.
#'
#' @section Anti-virus is blocking R:
#' Strangely enough, running R in parallel sometimes looks like you are accessing the internet to anti-virus software. So, it may block access to other instances of R. You will have to do some surgery on your anti-virus software settings to find where to change this.
#'
#' @section Your R packages are not stored in the "traditional" place:
#' R has a default directory where packages are stored on any system. If your packages are stored in a different place, worker nodes may not be able to find them \emph{if} you use \code{\link{setwd}} to change the working directory. I do not know if you have to set the working directory back to the default for your system, or if you have to change it to the folder that \emph{contains} the folder where your R packages reside (for me, they are the same directory). You can see what your current working directory is using \code{\link{getwd}}. RStudio will often change this directory automatically.
#'
#' So, if you get this error, try using \code{\link{setwd}} to set your working directory to the default one for your system, or to the folder that contains the folder that contains your packages.
#' 
#' @section Let me know:
#' I'm always game to help you track down your problems (with this package, not necessarily in general). The best way is to create an issue on \href{https://github.com/adamlilith/enmSdmX/issues}{GitHub}.
#'
#' @section Exorcise your computer:
#' Not responsible for damage to your computer.
#'
#' @keywords tutorial
NULL
