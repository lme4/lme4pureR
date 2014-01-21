##' Make template for relative covariance factor
##'
##' @param nc Number of columns in a dense model matrix for a particular
##' random effects term
##' @examples
##' mkTemplate(5)
mkTemplate <- function(nc){
                                        # generate row (i) and column (j) indices
                                        # of the upper triangular template matrix
    i <- sequence(1:nc); j <- rep(1:nc,1:nc)
                                        # generate theta:
                                        # 1) 1's along the diagonal (i.e. when i==j)
                                        # 2) 0's above the diagonal (i.e. when i!=j)
    theta <- 1*(i==j)
                                        # return the template using triplet (i,j,theta)
    sparseMatrix(i=i,j=j,x=theta)
}

##' Make templates for relative covariance factor
mkTemplates <- function(nc) lapply(nc, mkTemplate)

##' Make vector of indices giving the mapping from theta to Lambdat
mkLind <- function(nl, nc){
                                        # number of thetas per term (i.e. number
                                        # of elements in the upper triangle of
                                        # each of the template matrices)
    nTheta <- choose(nc+1, 2)
                                        # Lind per template
    templateLind <- lapply(nTheta, seq_len)
                                        # 0-based pointers to where each term
                                        # begins in theta
    offsetTheta <- c(0,cumsum(nTheta[-length(nTheta)]))
                                        # add offsets (i.e. pointers)
    templateLindWithOffset <- mapply("+", templateLind, offsetTheta)
                                        # repeat template-specific Lind vectors
                                        # once for each term and return Lind
    unlist(rep(templateLindWithOffset, nl))
}

##' Make initial relative covariance factor from list of templates
mkLambdat <- function(templates, nl){
                                        # repeat templates once for each level
    repTemplates <- rep(templates, nl)
                                        # return Lambdat by putting blocks
                                        # together along the diagonal
    .bdiag(repTemplates)
}    

##' Make initial theta from list of templates
mkTheta <- function(templates){
    thetas <- lapply(templates, slot, "x")
    unlist(thetas)
}
