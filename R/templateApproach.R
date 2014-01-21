##' Create a section of a transposed random effects model matrix
##'
##' @param grp Grouping factor for a particular random effects term.
##' @param mm Dense model matrix for a particular random effects term.
##' @return Section of a random effects design matrix corresponding to a
##' particular term.
##' @examples
##' ## consider a term (x | g) with:
##' ## number of observations, n = 6
##' ## number of levels, nl = 3
##' ## number of columns ('predictors'), nc = 2
##' (X <- cbind("(Intercept)"=1,x=1:6)) # an intercept in the first column
##'                                     # and 1:6 predictor in the other
##' (g <- as.factor(letters[rep(1:3,2)])) # grouping factor
##' nrow(X) # n = 6
##' nrow(X) == length(g) # and consistent n between X and g
##' ncol(X) # nc = 2
##' nlevels(g) # nl = 3
##' Zsection(g, X)
mkZtSection <- function(grp,mm) {
    Jt <- as(as.factor(grp), Class="sparseMatrix")
    KhatriRao(Jt,t(mm))
}

##' Make transposed random-effects model matrix
mkZt <- function(grp,mm){
    ZtSections <- mapply(mkZtSection, grp, mm)
    rBind(ZtSections)
}


##' Make a single template for a relative covariance factor
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

##' Make list of templates for relative covariance factor
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

