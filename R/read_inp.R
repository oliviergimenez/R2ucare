#' Read capture-recapture data with Input (.inp) format used by program MARK
#'
#' This function reads in capture-recapture dataset with the Input format.
#' It is a wrapper for the function convert.inp from package RMark. It drops continuous covariates because no goodness-of-fit test exists for such models
#' @param file text file with Input format (extension .inp)
#' @param group.df dataframe with grouping variables; contains a row for each group defined in the input file row1=group1, row2=group2 etc. Names and number of columns in the dataframe is set by user to define grouping variables in RMark dataframe
#' @return list with first component the matrix of encounter histories, second components the vector of number of individuals with corresponding histories and, if relevant, third component vector/matrix with group(s)
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>
#' @keywords package
#' @export
#' @examples
#' # read in Dipper dataset
#' dipper = system.file("extdata", "ed.inp", package = "R2ucare")
#' read_inp(dipper,group.df=data.frame(sex=c('Male','Female')))
#' # read in Geese dataset
#' geese = system.file("extdata", "geese.inp", package = "R2ucare")
#' read_inp(geese)

read_inp <- function(file,group.df=NULL){

# read in data, all columns as character, and ignore comments
data = RMark::convert.inp(file,group.df=group.df)

# add spaces between columns:
enc_hist = matrix(as.numeric(unlist(strsplit(data$ch, ''))),nrow=nrow(data),byrow=T)
counts = data$freq

# return list of results
if (is.null(group.df)){
	return(list(encounter_histories=enc_hist,sample_size=counts))
} else {
	return(list(encounter_histories=enc_hist,sample_size=counts,groups=data[,names(group.df)]))
}

} # end of function
