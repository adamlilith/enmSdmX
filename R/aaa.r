# .onLoad <- function(lib, pkg) {

# }

.onAttach <- function(lib, pkg) {

	ver <- read.dcf(file=system.file('DESCRIPTION', package=pkg), fields='Version')
	packageStartupMessage(paste0('This is ', pkg, ' ', ver, '.'))

	packageStartupMessage(paste0('Back-incompatible changes starting starting with enmSdmX version 1.1.1 (2023-06-11):'))
	packageStartupMessage(paste0('* trainRF() replaces use of the randomForest package with the faster ranger package, which produces random forests that are statistically equivalent.'))
	packageStartupMessage(paste0('* geoFold() uses ', dQuote('complete'), ' clustering by default (was ', dQuote('single'), ').'))
	
	packageStartupMessage(paste0('Back-incompatible changes starting starting with enmSdmX version 1.2.10 (2024-12-07):'))
	packageStartupMessage(paste0('* By default, trainGLM() and trainNS() remove unconverged models and models with boundary issues.'))
	
}
