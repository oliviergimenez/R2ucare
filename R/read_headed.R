#' Read capture-recapture data with Headed format used by program E-SURGE
#'
#' This function reads in capture-recapture dataset with the Headed format.
#' It ignores all forms of censorship for now, and drops continuous covariates because no goodness-of-fit test exists for such models
#' @param file text file with Headed format
#' @return list with first component the matrix of encounter histories, second components the vector of number of individuals with corresponding histories and, if relevant, third component vector/matrix with group(s)
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>
#' @keywords package
#' @export
#' @examples
#' # read in Dipper dataset
#' dipper = system.file("extdata", "ed.txt", package = "R2ucare")
#' read_headed(dipper)
#' # read in Geese dataset
#' geese = system.file("extdata", "geese.txt", package = "R2ucare")
#' read_headed(geese)

read_headed <- function(file){

# read in data, all columns as character, and ignore comments
data = utils::read.table(file,header=T,colClasses='character', comment.char = '/')

# get columns names
headers = names(data)

# get the columns with pre-defined E-SURGE labels (note that ':' has become '.' and '$' is now 'X' after reading the dataset)
where_occasions = which(as.logical(stringr::str_count(headers, "H.")))
where_sample_size = which(as.logical(stringr::str_count(headers, "S.")))
where_covariates = which(as.logical(stringr::str_count(headers, "X.COV.")))

# initialize the group covariate(s)
group = NULL

# if there is a single column with the headed 'H:', occasions in columns are not separated by spaces
if (length(where_occasions) == 1){
	enc_hist = matrix(as.numeric(unlist(strsplit(data[, where_occasions], ''))),nrow = length(strsplit(data[, where_occasions], '')),byrow=T) # get encounter histories
	counts = as.numeric(data[, where_sample_size]) # counts
	if (length(where_covariates)!=0) group = data[, where_covariates] # groups if needed
} else {
# if there is more than one column with the headed H:, it means that occasions in columns are separated by spaces
	enc_hist = data.matrix(data[, where_occasions]) # get encounter histories
	counts = as.numeric(data[, where_sample_size]) # get counts
	if (length(where_covariates)!=0) group = data[, where_covariates] # groups if needed
}

# return list of results
if (is.null(group)){
	return(list(encounter_histories=enc_hist,sample_size=counts))
} else {
	return(list(encounter_histories=enc_hist,sample_size=counts,groups=group))
}

} # end of function
