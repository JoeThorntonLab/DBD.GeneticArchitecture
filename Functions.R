########################################
#Functions for analyzing epistasis data#
########################################
###Options for connecting genotypes in sequence space include:
#code=S Standard genetic code
#code=Z Standard genetic code that splits serine codons into two disconnected groups
#code=N No genetic code, i.e. Hamming distance

##Identify neighboring genotypes
#Get all possible genetic neighbors of a genotype 
get.neighbors <- function(AAseq,code="S"){
  if(code == "S") {
    AA.NEIGHBORS <- AA.NEIGHBORS.S
  } else if (code == "Z") {
    AA.NEIGHBORS <- AA.NEIGHBORS.Z
  } else if (code == "N") {
    AA.NEIGHBORS <- AA.NEIGHBORS.N  
  }  
  seq <- unlist(strsplit(as.character(AAseq),split=""))
  mutant <- vector()
  J <- 1:length(seq)
  for(j in J) {
    for(i in AAs.extend[AA.NEIGHBORS[seq[j],]==1]){
      AAmut <- AAseq
      substr(AAmut,j,j) <- i
      mutant <- c(mutant,AAmut)
    }
  }
  return(mutant)
}
#Get adjacent genotypes from a list
index.adjaceny <- function(list,code="S") {
  lapply(list,FUN = function(y){
    x <- which(list %in% get.neighbors(y,code=code))
    return(x)})
}
#Make adjaceny matrix
make.adjaceny <- function(list,index) {
  ADJM <- as(sparseMatrix(
    i = rep(1:length(list),times=sapply(index,length)),
    j = unlist(index)),"dgCMatrix")
  colnames(ADJM) <- list; rownames(ADJM) <- list
  return(ADJM)
}

##Genetic distance and routes between genotypes
#Get genetic distance between two sequences
get.distance <- function(SOURCE,TARGET,code="S",vectorized=FALSE) {
  if(code == "S") {
    AA.DISTANCES <- AA.DISTANCES.S
  } else if (code == "Z") {
    AA.DISTANCES <- AA.DISTANCES.Z
  } else if (code == "N") {
    AA.DISTANCES <- AA.DISTANCES.N
  }
  I <- 1:nchar(SOURCE)
  if(vectorized == TRUE) {
    DIST <- diag(as.matrix(AA.DISTANCES[substring(SOURCE,I,I),substring(TARGET,I,I)]))
  } else {
  DIST <- sum(diag(as.matrix(AA.DISTANCES[substring(SOURCE,I,I),substring(TARGET,I,I)])))
  }
  return(DIST)
}
#Faster version of genetic distance with hamming distance
get.HAM.distance <- function(SOURCE,TARGET) {
  SOURCE.SEQ <- unlist(strsplit(as.character(SOURCE),split=""))
  TARGET.SEQ <- unlist(strsplit(as.character(TARGET),split=""))
  
  if(length(SOURCE.SEQ) != length(TARGET.SEQ)) {
    stop("Sequences of different lengths")
  }
  else {
    DIST <- sum(SOURCE.SEQ != TARGET.SEQ)
  }
  return(DIST)
}
#Get all possible shortest paths between two genotypes
get.routes <- function(SOURCE,TARGET,code="S",vectorized=FALSE) {
  if(code == "S") {
    NUM.ROUTES <- GC.ROUTES.S
  } else if (code == "Z") {
    NUM.ROUTES <- GC.ROUTES.Z
  }
  I <- 1:nchar(SOURCE)
  if(vectorized == TRUE) {
    ROUTES <- diag(as.matrix(NUM.ROUTES[substring(SOURCE,I,I),substring(TARGET,I,I)]))
  } else {
    ROUTES <- sum(diag(as.matrix(NUM.ROUTES[substring(SOURCE,I,I),substring(TARGET,I,I)])))
  }
  return(ROUTES)
}

##Functions for approximate maximum likelihood fit of proportional odds model using glmnet
#Modifies code from glmnetcr package to include a forward proportional odds model
#Modifies code from glmnet package to allow for long vectors in sparse matrices

#Check proportional odds assumption
test.PO <- function(y){
  c("Y>=1"=qlogis(mean( y >= 1)), "Y>=2"=qlogis(mean( y >= 2)))
}
#Fitting of logistic models that can handle long vectors and sparse matrices
lognet64 <- function (x, is.sparse, ix, jx, y, weights, offset, alpha, nobs, 
                      nvars, jd, vp, cl, ne, nx, nlam, flmin, ulam, thresh, isd, 
                      intr, vnames, maxit, kopt, family) {
  nc = dim(y)
  maxit = as.integer(maxit)
  if (is.null(nc)) {
    y = as.factor(y)
    ntab = table(y)
    minclass = min(ntab)
    if (minclass <= 1) 
      stop("one multinomial or binomial class has 1 or 0 observations; not allowed")
    if (minclass < 8) 
      warning("one multinomial or binomial class has fewer than 8  observations; dangerous ground")
    classnames = names(ntab)
    nc = as.integer(length(ntab))
    y = diag(nc)[as.numeric(y), ]
  }
  else {
    noo = nc[1]
    if (noo != nobs) 
      stop("x and y have different number of rows in call to glmnet", call. = FALSE)
    nc = as.integer(nc[2])
    classnames = colnames(y)
  }
  maxvars = .Machine$integer.max/(nlam * nc)
  if (nx > maxvars) 
    stop(paste("Integer overflow; num_classes*num_lambda*pmax should not exceed .Machine$integer.max. Reduce pmax to be below", trunc(maxvars)))
  if (!missing(weights)) 
    y = y * weights
  weights = drop(y %*% rep(1, nc))
  o = weights > 0
  if (!all(o)) {
    y = y[o, ]
    if (is.sparse) {
      x = sparseMatrix(i = jx, p = ix - 1, x = x, dims = c(nobs,nvars))[o, , drop = FALSE]
      ix = as.integer(x@p + 1)
      jx = as.integer(x@i + 1)
      x = as.double(x@x)
    }
    else x = x[o, , drop = FALSE]
    nobs = sum(o)
  }
  else o = NULL
  if (family == "binomial") {
    if (nc > 2) 
      stop("More than two classes; use multinomial family instead in call to glmnet", call. = FALSE)
    nc = as.integer(1)
    y = y[, c(2, 1)]
  }
  storage.mode(y) = "double"
  if (is.null(offset)) {
    offset = y * 0
    is.offset = FALSE
  }
  else {
    offset = as.matrix(offset)
    if (!is.null(o)) 
      offset = offset[o, , drop = FALSE]
    do = dim(offset)
    if (do[[1]] != nobs) 
      stop("offset should have the same number of values as observations in binomial/multinomial call to glmnet",  call. = FALSE)
    if ((do[[2]] == 1) & (nc == 1)) 
      offset = cbind(offset, -offset)
    if ((family == "multinomial") & (do[[2]] != nc)) 
      stop("offset should have same shape as y in multinomial call to glmnet", call. = FALSE)
    storage.mode(offset) = "double"
    is.offset = TRUE
  }
  fit = if (is.sparse) 
    .C64("splognet", SIGNATURE = c("double", "integer","integer","integer","double","integer","integer",
    							                 "double", "double", "integer","double", "double",
                                   "integer","integer","integer","double", "double",
                                   "double", "integer","integer","integer","integer",
                                   "integer","double", "double", "integer","integer",
                                   "double", "double", "double", "integer","integer"),
         parm = alpha, nobs, nvars, nc, x, ix, jx, 
         y, offset, jd, vp, cl, 
         ne = ne, nx, nlam,flmin, ulam, 
         thresh, isd, intr, maxit, kopt, 
         lmu = integer(1), a0 = double(nlam * nc), ca = double(nx * nlam * nc), ia = integer(nx), nin = integer(nlam), 
         nulldev = double(1),dev = double(nlam), alm = double(nlam), nlp = integer(1),jerr = integer(1), 
         INTENT = c(rep("rw",4),"r",rep("rw",17),rep("w",10)), 
         PACKAGE = "glmnet",
         NAOK = TRUE)
   else .C64("lognet", SIGNATURE =   c("double", "integer","integer","integer","double",
                                       "double", "double", "integer","double", "double",
                                       "integer","integer","integer","double", "double",
                                       "double", "integer","integer","integer","integer",
                                       "integer","double", "double", "integer","integer",
                                       "double", "double", "double", "integer","integer"),
          parm = alpha,nobs,           nvars,               nc,      as.double(x), 
          y,           offset,         jd,                  vp,      cl, 
          ne,          nx,             nlam,                flmin,   ulam, 
          thresh,      isd,            intr,                maxit,   kopt, 
          lmu = integer(1),    a0 = double(nlam * nc), ca = double(nx * nlam * nc), ia = integer(nx), nin = integer(nlam), 
          nulldev = double(1), dev = double(nlam),     alm = double(nlam),          nlp = integer(1), jerr = integer(1), 
          INTENT = c(rep("rw",4),"r",rep("rw",15),rep("w",10)),     
          PACKAGE = "glmnet",
          NAOK = TRUE)
  if (fit$jerr != 0) {
    errmsg = jerr(fit$jerr, maxit, pmax = nx, family)
    if (errmsg$fatal) 
      stop(errmsg$msg, call. = FALSE)
    else warning(errmsg$msg, call. = FALSE)
  }
  if (family == "binomial") {
    outlist = getcoef(fit, nvars, nx, vnames)
  }
  else outlist = getcoef.multinomial(fit, nvars, nx, vnames, 
                                     nc, classnames)
  dev = fit$dev[seq(fit$lmu)]
  outlist = c(outlist, list(dev.ratio = dev, nulldev = fit$nulldev, 
                            npasses = fit$nlp, jerr = fit$jerr, offset = is.offset, 
                            classnames = classnames))
  if (family == "multinomial") {
    if (kopt == 2) 
      grouped = TRUE
    else grouped = FALSE
    outlist$grouped = grouped
  }
  class(outlist) = switch(family, binomial = "lognet", multinomial = "multnet")
  outlist
}
#Same as glmnet, but uses dotCall64 package to handle passing long vectors to fortran
glmnet64 <- function (x, y, family = c("gaussian", "binomial", "poisson","multinomial", "cox", "mgaussian"), weights, offset = NULL, 
                      alpha = 1, nlambda = 100, lambda.min.ratio = ifelse(nobs <  nvars, 0.01, 1e-04), lambda = NULL, standardize = TRUE, 
                      intercept = TRUE, thresh = 1e-07, dfmax = nvars + 1, pmax = min(dfmax *2 + 20, nvars), exclude, penalty.factor = rep(1, nvars), 
                      lower.limits = -Inf, upper.limits = Inf, maxit = 1e+05, type.gaussian = ifelse(nvars < 500, "covariance", "naive"), 
                      type.logistic = c("Newton", "modified.Newton"), standardize.response = FALSE, type.multinomial = c("ungrouped",  "grouped")) {
  family = match.arg(family)
  if (alpha > 1) {
    warning("alpha >1; set to 1")
    alpha = 1
  }
  if (alpha < 0) {
    warning("alpha<0; set to 0")
    alpha = 0
  }
  alpha = as.double(alpha)
  this.call = match.call()
  nlam = as.integer(nlambda)
  y = drop(y)
  np = dim(x)
  if (is.null(np) | (np[2] <= 1)) 
    stop("x should be a matrix with 2 or more columns")
  nobs = as.integer(np[1])
  if (missing(weights)) 
    weights = rep(1, nobs)
  else if (length(weights) != nobs) 
    stop(paste("number of elements in weights (", length(weights), ") not equal to the number of rows of x (", nobs, ")", sep = ""))
  nvars = as.integer(np[2])
  dimy = dim(y)
  nrowy = ifelse(is.null(dimy), length(y), dimy[1])
  if (nrowy != nobs) 
    stop(paste("number of observations in y (", nrowy, ") not equal to the number of rows of x (", nobs, ")", sep = ""))
  vnames = colnames(x)
  if (is.null(vnames)) 
    vnames = paste("V", seq(nvars), sep = "")
  ne = as.integer(dfmax)
  nx = as.integer(pmax)
  if (missing(exclude)) 
    exclude = integer(0)
  if (any(penalty.factor == Inf)) {
    exclude = c(exclude, seq(nvars)[penalty.factor == Inf])
    exclude = sort(unique(exclude))
  }
  if (length(exclude) > 0) {
    jd = match(exclude, seq(nvars), 0)
    if (!all(jd > 0)) 
      stop("Some excluded variables out of range")
    penalty.factor[jd] = 1
    jd = as.integer(c(length(jd), jd))
  }
  else jd = as.integer(0)
  vp = as.double(penalty.factor)
  internal.parms = glmnet.control()
  if (any(lower.limits > 0)) {
    stop("Lower limits should be non-positive")
  }
  if (any(upper.limits < 0)) {
    stop("Upper limits should be non-negative")
  }
  lower.limits[lower.limits == -Inf] = -internal.parms$big
  upper.limits[upper.limits == Inf] = internal.parms$big
  if (length(lower.limits) < nvars) {
    if (length(lower.limits) == 1) 
      lower.limits = rep(lower.limits, nvars)
    else stop("Require length 1 or nvars lower.limits")
  }
  else lower.limits = lower.limits[seq(nvars)]
  if (length(upper.limits) < nvars) {
    if (length(upper.limits) == 1) 
      upper.limits = rep(upper.limits, nvars)
    else stop("Require length 1 or nvars upper.limits")
  }
  else upper.limits = upper.limits[seq(nvars)]
  cl = rbind(lower.limits, upper.limits)
  if (any(cl == 0)) {
    fdev = glmnet.control()$fdev
    if (fdev != 0) {
      glmnet.control(fdev = 0)
      on.exit(glmnet.control(fdev = fdev))
    }
  }
  storage.mode(cl) = "double"
  isd = as.integer(standardize)
  intr = as.integer(intercept)
  if (!missing(intercept) && family == "cox") 
    warning("Cox model has no intercept")
  jsd = as.integer(standardize.response)
  thresh = as.double(thresh)
  if (is.null(lambda)) {
    if (lambda.min.ratio >= 1) 
      stop("lambda.min.ratio should be less than 1")
    flmin = as.double(lambda.min.ratio)
    ulam = double(1)
  }
  else {
    flmin = as.double(1)
    if (any(lambda < 0)) 
      stop("lambdas should be non-negative")
    ulam = as.double(rev(sort(lambda)))
    nlam = as.integer(length(lambda))
  }
  is.sparse = FALSE
  ix = jx = NULL
  if (inherits(x, "sparseMatrix")) {
    is.sparse = TRUE
    x = as(x, "CsparseMatrix")
    x = as(x, "dgCMatrix")
    ix = as.integer(x@p + 1)
    jx = as.integer(x@i + 1)
    x = as.double(x@x)
  }
  kopt = switch(match.arg(type.logistic), Newton = 0, modified.Newton = 1)
  if (family == "multinomial") {
    type.multinomial = match.arg(type.multinomial)
    if (type.multinomial == "grouped") 
      kopt = 2
  }
  kopt = as.integer(kopt)
  fit = switch(family, 
               gaussian    = elnet(   x, is.sparse, ix, jx, y, weights, offset, type.gaussian, alpha, nobs, nvars, jd, vp, cl, ne, nx, nlam, flmin, ulam, thresh, isd,      intr, vnames, maxit), 
               poisson     = fishnet( x, is.sparse, ix, jx, y, weights, offset,                alpha, nobs, nvars, jd, vp, cl, ne, nx, nlam, flmin, ulam, thresh, isd,      intr, vnames, maxit), 
               binomial    = lognet64(x, is.sparse, ix, jx, y, weights, offset,                alpha, nobs, nvars, jd, vp, cl, ne, nx, nlam, flmin, ulam, thresh, isd,      intr, vnames, maxit, kopt, family), 
               multinomial = lognet(  x, is.sparse, ix, jx, y, weights, offset,                alpha, nobs, nvars, jd, vp, cl, ne, nx, nlam, flmin, ulam, thresh, isd,      intr, vnames, maxit, kopt, family), 
               cox         = coxnet(  x, is.sparse, ix, jx, y, weights, offset,                alpha, nobs, nvars, jd, vp, cl, ne, nx, nlam, flmin, ulam, thresh, isd,            vnames, maxit), 
               mgaussian   = mrelnet( x, is.sparse, ix, jx, y, weights, offset,                alpha, nobs, nvars, jd, vp, cl, ne, nx, nlam, flmin, ulam, thresh, isd, jsd, intr, vnames, maxit))
  #if (is.null(lambda)) 
    #fit$lambda = fix.lam(fit$lambda)
  fit$call = this.call
  fit$nobs = nobs
  class(fit) = c(class(fit), "glmnet")
  fit
}

#Restructures a model matrix to implement a cumulative odds model with the proportional odds assumption
po.forward <- function (x, y, weights) {
  yname <- as.character(substitute(y))
  if (!is.factor(y)) 
    y <- factor(y, exclude = NA)
  ylevels <- levels(y)
  kint <- length(ylevels) - 1
  y <- as.integer(y)
  names <- dimnames(x)[[2]]
  if (length(names) == 0) 
    names <- paste("V", 1:dim(x)[2], sep = "")
  expand <- list()
  for (k in 1:kint) {
    expand[[k]] <- cbind(y,weights, x)
    expand[[k]][, 1] <- ifelse(expand[[k]][, 1] <= k, 1, 0)
    cp <- matrix(rep(0, dim(expand[[k]])[1] * kint), ncol = kint)
    cp[, k] <- 1
    dimnames(cp)[[2]] <- paste("cp", 1:kint, sep = "")
    expand[[k]] <- cbind(expand[[k]], cp)
    dimnames(expand[[k]])[[2]] <- c("y", "weights", names,paste("cp", 1:kint, sep = ""))
  }
  newx <- expand[[1]]
  for (k in 2:kint) newx <- rbind(newx, expand[[k]])
  newx
}

#Same as glmnetcr, but includes approximate proportional odds model and sparse matrices
glmnetpo <- function (x, y, method = "backward", weights = NULL, offset = NULL, 
                      alpha = 1, nlambda = 100, lambda.min.ratio = NULL, lambda = NULL, 
                      standardize = TRUE, thresh = 1e-04, exclude, penalty.factor = NULL, 
                      maxit = 100) {
  if (length(unique(y)) == 2) 
    stop("Binary response: Use glmnet with family='binomial' parameter")
  n <- nobs <- dim(x)[1]
  p <- m <- nvars <- dim(x)[2]
  k <- length(unique(y))
  #x <- as.matrix(x)
  x <- as(x, "dgCMatrix")
  if (is.null(penalty.factor)) 
    penalty.factor <- rep(1, nvars)
  else penalty.factor <- penalty.factor
  if (is.null(lambda.min.ratio)) 
    lambda.min.ratio <- ifelse(nobs < nvars, 0.01, 1e-04)
  if (is.null(weights)) 
    weights <- rep(1, length(y))
  if (c("backward", "forward", "cumulative")[charmatch(method, c("backward", "forward", "cumulative"))] == "backward") {
    restructure <- cr.backward(x = x, y = y, weights = weights)
  }
  if (c("backward", "forward", "cumulative")[charmatch(method, c("backward", "forward", "cumulative"))] == "forward") {
    restructure <- cr.forward(x = x, y = y, weights = weights)
  }
  if (c("backward", "forward", "cumulative")[charmatch(method, c("backward", "forward", "cumulative"))] == "cumulative") {
    restructure <- po.forward(x = x, y = y, weights = weights)
  }
  glmnet.data <- list(x = restructure[, -c(1, 2)], y = restructure[,"y"], weights = restructure[, "weights"])
  object <- glmnet(glmnet.data$x, glmnet.data$y, family = "binomial", 
                   weights = glmnet.data$weights, offset = offset, alpha = alpha, 
                   nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, 
                   lambda = lambda, standardize = standardize, thresh = thresh, 
                   exclude = exclude, penalty.factor = c(penalty.factor, rep(0, k - 1)), 
                   maxit = maxit, type.gaussian = ifelse(nvars < 500, "covariance", "naive"))
  object$x <- x
  object$y <- y
  object$method <- method
  class(object) <- "glmnetcr"
  object
}

#Same as glmnetpo, but calls glmnet64 for handling long vectors
glmnetpo64 <- function (x, y, method = "backward", weights = NULL, offset = NULL, 
                        alpha = 1, nlambda = 100, lambda.min.ratio = NULL, lambda = NULL, 
                        standardize = TRUE, thresh = 1e-04, exclude, penalty.factor = NULL, 
                        maxit = 100) {
  if (length(unique(y)) == 2) 
    stop("Binary response: Use glmnet with family='binomial' parameter")
  n <- nobs <- dim(x)[1]
  p <- m <- nvars <- dim(x)[2]
  k <- length(unique(y))
  #x <- as.matrix(x)
  x <- as(x, "dgCMatrix")
  if (is.null(penalty.factor)) 
    penalty.factor <- rep(1, nvars)
  else penalty.factor <- penalty.factor
  if (is.null(lambda.min.ratio)) 
    lambda.min.ratio <- ifelse(nobs < nvars, 0.01, 1e-04)
  if (is.null(weights)) 
    weights <- rep(1, length(y))
  if (c("backward", "forward", "cumulative")[charmatch(method, c("backward", "forward", "cumulative"))] == "backward") {
    restructure <- cr.backward(x = x, y = y, weights = weights)
  }
  if (c("backward", "forward", "cumulative")[charmatch(method, c("backward", "forward", "cumulative"))] == "forward") {
    restructure <- cr.forward(x = x, y = y, weights = weights)
  }
  if (c("backward", "forward", "cumulative")[charmatch(method, c("backward", "forward", "cumulative"))] == "cumulative") {
    restructure <- po.forward(x = x, y = y, weights = weights)
  }
  glmnet.data <- list(x = restructure[, -c(1, 2)], y = restructure[,"y"], weights = restructure[, "weights"])
  object <- glmnet64(glmnet.data$x, glmnet.data$y, family = "binomial", 
                   weights = glmnet.data$weights, offset = offset, alpha = alpha, 
                   nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, 
                   lambda = lambda, standardize = standardize, thresh = thresh, 
                   exclude = exclude, penalty.factor = c(penalty.factor, rep(0, k - 1)), 
                   maxit = maxit, type.gaussian = ifelse(nvars < 500, "covariance", "naive"))
  object$x <- x
  object$y <- y
  object$method <- method
  class(object) <- "glmnetcr"
  object
}

#Same as cv.glmnet function, but can take a glmnetcr object
cv.glmnetcr <- function (x, y, weights, offset = NULL, lambda = NULL, 
                         type.measure = c("mse","deviance", "class", "auc", "mae"), nfolds = 10, foldid, 
                         grouped = TRUE, keep = FALSE, parallel = FALSE, ...) {
    if (missing(type.measure)) 
      type.measure = "default"
    else type.measure = match.arg(type.measure)
    if (!is.null(lambda) && length(lambda) < 2) 
      stop("Need more than one value of lambda for cv.glmnet")
    N = nrow(x)
    if (missing(weights)) 
      weights = rep(1, N)
    else weights = as.double(weights)
    y = drop(y)
    glmnet.call = match.call(expand.dots = TRUE)
    which = match(c("type.measure", "nfolds", "foldid", "grouped", "keep"), names(glmnet.call), F)
    if (any(which)) 
      glmnet.call = glmnet.call[-which]
    glmnet.call[[1]] = as.name("glmnet")
    glmnet.object = glmnetcr(x, y, weights = weights, offset = offset, lambda = lambda, ...)
    glmnet.object$call = glmnet.call
    subclass = class(glmnet.object)[[1]]
    subclass = "lognet"
    type.measure = cvtype(type.measure, subclass)
    is.offset = glmnet.object$offset
    if (inherits(glmnet.object, "multnet") && !glmnet.object$grouped) {
      nz = predict(glmnet.object, type = "nonzero")
      nz = sapply(nz, function(x) sapply(x, length))
      nz = ceiling(apply(nz, 1, median))
    } 
    #else nz = sapply(predict(glmnet.object, type = "nonzero"), length)
    nz = glmnet.object$df
    if (missing(foldid)) 
      foldid = sample(rep(seq(nfolds), length = N))
    else nfolds = max(foldid)
    if (nfolds < 3) 
      stop("nfolds must be bigger than 3; nfolds=10 recommended")
    outlist = as.list(seq(nfolds))
    if (parallel) {
      outlist = foreach(i = seq(nfolds), .packages = c("glmnet")) %dopar% 
      {
        which = foldid == i
        if (length(dim(y)) > 1) 
          y_sub = y[!which, ]
        else y_sub = y[!which]
        if (is.offset) 
          offset_sub = as.matrix(offset)[!which, ]
        else offset_sub = NULL
        glmnetcr(x[!which, , drop = FALSE], y_sub, lambda = lambda,  offset = offset_sub, weights = weights[!which],...)
      }
    }
    else {
      for (i in seq(nfolds)) {
        which = foldid == i
        if (is.matrix(y)) 
          y_sub = y[!which, ]
        else y_sub = y[!which]
        if (is.offset) 
          offset_sub = as.matrix(offset)[!which, ]
        else offset_sub = NULL
        outlist[[i]] = glmnetcr(x[!which, , drop = FALSE],y_sub, lambda = lambda, offset = offset_sub,weights = weights[!which], ...)
      }
    }
    fun = paste("cv.cr", subclass, sep = ".")
    lambda = glmnet.object$lambda
	cvstuff = do.call(fun, list(outlist, lambda, x, y, weights, offset, foldid, type.measure, grouped, keep))
    cvm = cvstuff$cvm
    cvsd = cvstuff$cvsd
    nas = is.na(cvsd)
    if (any(nas)) {
      lambda = lambda[!nas]
      cvm = cvm[!nas]
      cvsd = cvsd[!nas]
      nz = nz[!nas]
    }
    cvname = names(cvstuff$type.measure)
    names(cvname) = cvstuff$type.measure
    out = list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm + cvsd, cvlo = cvm - cvsd, nzero = nz, name = cvname, glmnet.fit = glmnet.object)
    if (keep) 
      out = c(out, list(fit.preval = cvstuff$fit.preval, foldid = foldid))
    lamin = if (cvname == "AUC") 
      getmin(lambda, -cvm, cvsd)
    else getmin(lambda, cvm, cvsd)
    obj = c(out, as.list(lamin))
    class(obj) = "cv.glmnet"
    obj
  }

#Same as cv.glmnetcr function, but can take a glmnetpo object
cv.glmnetpo <- function (x, y, weights, offset = NULL, lambda = NULL, 
                         type.measure = c("mse","deviance", "class", "auc", "mae"), nfolds = 10, foldid, 
                         grouped = TRUE, keep = FALSE, parallel = FALSE, ...) {
    if (missing(type.measure)) 
      type.measure = "default"
    else type.measure = match.arg(type.measure)
    if (!is.null(lambda) && length(lambda) < 2) 
      stop("Need more than one value of lambda for cv.glmnet")
    N = nrow(x)
    if (missing(weights)) 
      weights = rep(1, N)
    else weights = as.double(weights)
    y = drop(y)
    glmnet.call = match.call(expand.dots = TRUE)
    which = match(c("type.measure", "nfolds", "foldid", "grouped", "keep"), names(glmnet.call), F)
    if (any(which)) 
      glmnet.call = glmnet.call[-which]
    glmnet.call[[1]] = as.name("glmnet")
    glmnet.object = glmnetpo(x, y, weights = weights, offset = offset, lambda = lambda, ...)
    glmnet.object$call = glmnet.call
    #subclass = class(glmnet.object)[[1]]
    subclass = "lognet"
    type.measure = cvtype(type.measure, subclass)
    is.offset = glmnet.object$offset
    if (inherits(glmnet.object, "multnet") && !glmnet.object$grouped) {
      nz = predict.glmnetpo(glmnet.object, type = "nonzero")
      nz = sapply(nz, function(x) sapply(x, length))
      nz = ceiling(apply(nz, 1, median))
    } 
    #else nz = sapply(predict.glmnetpo(glmnet.object, type = "nonzero"), length)
    nz = glmnet.object$df
    if (missing(foldid)) 
      foldid = sample(rep(seq(nfolds), length = N))
    else nfolds = max(foldid)
    if (nfolds < 3) 
      stop("nfolds must be bigger than 3; nfolds=10 recommended")
    outlist = as.list(seq(nfolds))
    if (parallel) {
      outlist = foreach(i = seq(nfolds), .packages = c("glmnet")) %dopar% 
      {
        which = foldid == i
        if (length(dim(y)) > 1) 
          y_sub = y[!which, ]
        else y_sub = y[!which]
        if (is.offset) 
          offset_sub = as.matrix(offset)[!which, ]
        else offset_sub = NULL
        glmnetpo(x[!which, , drop = FALSE], y_sub, lambda = lambda,  offset = offset_sub, weights = weights[!which],...)
      }
    }
    else {
      for (i in seq(nfolds)) {
        which = foldid == i
        if (is.matrix(y)) 
          y_sub = y[!which, ]
        else y_sub = y[!which]
        if (is.offset) 
          offset_sub = as.matrix(offset)[!which, ]
        else offset_sub = NULL
        outlist[[i]] = glmnetpo(x[!which, , drop = FALSE],y_sub, lambda = lambda, offset = offset_sub,weights = weights[!which], ...)
      }
    }
    fun = paste("cv.po", subclass, sep = ".")
    lambda = glmnet.object$lambda
	cvstuff = do.call(fun, list(outlist, lambda, x, y, weights, offset, foldid, type.measure, grouped, keep))
    cvm = cvstuff$cvm
    cvsd = cvstuff$cvsd
    nas = is.na(cvsd)
    if (any(nas)) {
      lambda = lambda[!nas]
      cvm = cvm[!nas]
      cvsd = cvsd[!nas]
      nz = nz[!nas]
    }
    cvname = names(cvstuff$type.measure)
    names(cvname) = cvstuff$type.measure
    out = list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm + cvsd, cvlo = cvm - cvsd, nzero = nz, name = cvname, glmnet.fit = glmnet.object)
    if (keep) 
      out = c(out, list(fit.preval = cvstuff$fit.preval, foldid = foldid))
    lamin = if (cvname == "AUC") 
      getmin(lambda, -cvm, cvsd)
    else getmin(lambda, cvm, cvsd)
    obj = c(out, as.list(lamin))
    class(obj) = "cv.glmnet"
    obj
}

#Same as cv.glmnetcr function, but can take a glmnetpo object
cv.glmnetpo64 <- function (x, y, weights, offset = NULL, lambda = NULL, 
                           type.measure = c("mse","deviance", "class", "auc", "mae"), nfolds = 10, foldid, 
                           grouped = TRUE, keep = FALSE, parallel = FALSE, ...) {
  if (missing(type.measure)) 
    type.measure = "default"
  else type.measure = match.arg(type.measure)
  if (!is.null(lambda) && length(lambda) < 2) 
    stop("Need more than one value of lambda for cv.glmnet")
  N = nrow(x)
  if (missing(weights)) 
    weights = rep(1, N)
  else weights = as.double(weights)
  y = drop(y)
  glmnet.call = match.call(expand.dots = TRUE)
  which = match(c("type.measure", "nfolds", "foldid", "grouped", "keep"), names(glmnet.call), F)
  if (any(which)) 
    glmnet.call = glmnet.call[-which]
  glmnet.call[[1]] = as.name("glmnet")
  glmnet.object = glmnetpo64(x, y, weights = weights, offset = offset, lambda = lambda, ...)
  glmnet.object$call = glmnet.call
  #subclass = class(glmnet.object)[[1]]
  subclass = "lognet"
  type.measure = cvtype(type.measure, subclass)
  is.offset = glmnet.object$offset
  if (inherits(glmnet.object, "multnet") && !glmnet.object$grouped) {
    nz = predict.glmnetpo64(glmnet.object, type = "nonzero")
    nz = sapply(nz, function(x) sapply(x, length))
    nz = ceiling(apply(nz, 1, median))
  } 
  #else nz = sapply(predict.glmnetpo64(glmnet.object, type = "nonzero"), length)
  nz = glmnet.object$df
  if (missing(foldid)) 
    foldid = sample(rep(seq(nfolds), length = N))
  else nfolds = max(foldid)
  if (nfolds < 3) 
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  outlist = as.list(seq(nfolds))
  if (parallel) {
    outlist = foreach(i = seq(nfolds), .packages = c("glmnet")) %dopar% 
    {
      which = foldid == i
      if (length(dim(y)) > 1) 
        y_sub = y[!which, ]
      else y_sub = y[!which]
      if (is.offset) 
        offset_sub = as.matrix(offset)[!which, ]
      else offset_sub = NULL
      glmnetpo(x[!which, , drop = FALSE], y_sub, lambda = lambda,  offset = offset_sub, weights = weights[!which],...)
    }
  }
  else {
    for (i in seq(nfolds)) {
      which = foldid == i
      if (is.matrix(y)) 
        y_sub = y[!which, ]
      else y_sub = y[!which]
      if (is.offset) 
        offset_sub = as.matrix(offset)[!which, ]
      else offset_sub = NULL
      outlist[[i]] = glmnetpo64(x[!which, , drop = FALSE],y_sub, lambda = lambda, offset = offset_sub,weights = weights[!which], ...)
    }
  }
  outlist
  fun = paste("cv.po64", subclass, sep = ".")
  lambda = glmnet.object$lambda
  #cvstuff = do.call(fun, list(outlist, lambda, x, y, weights, offset, foldid, type.measure, grouped, keep))
  cvstuff = cv.po64.lognet(outlist, lambda, x, y, weights, offset, foldid, type.measure, grouped, keep)
  cvm = cvstuff$cvm
  cvsd = cvstuff$cvsd
  nas = is.na(cvsd)
  if (any(nas)) {
    lambda = lambda[!nas]
    cvm = cvm[!nas]
    cvsd = cvsd[!nas]
    nz = nz[!nas]
  }
  cvname = names(cvstuff$type.measure)
  names(cvname) = cvstuff$type.measure
  out = list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm + cvsd, cvlo = cvm - cvsd, nzero = nz, name = cvname, glmnet.fit = glmnet.object)
  if (keep) 
    out = c(out, list(fit.preval = cvstuff$fit.preval, foldid = foldid))
  lamin = if (cvname == "AUC") 
    getmin(lambda, -cvm, cvsd)
  else getmin(lambda, cvm, cvsd)
  obj = c(out, as.list(lamin))
  class(obj) = "cv.glmnet"
  obj
}

#Modified cv.lognet function to handle glmnetcr objects
cv.cr.lognet <- function (outlist, lambda, x, y, weights, offset, foldid, type.measure, 
    grouped, keep = FALSE) {
    prob_min = 1e-05
    prob_max = 1 - prob_min
    nc = dim(y)
    if (is.null(nc)) {
        y = as.factor(y)
        ntab = table(y)
        nc = as.integer(length(ntab))
        y = diag(nc)[as.numeric(y), ]
    }
    N = nrow(y)
    nfolds = max(foldid)
    if ((N/nfolds < 10) && type.measure == "auc") {
        warning("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds", 
            call. = FALSE)
        type.measure = cvtype("deviance", "lognet")
    }
    if ((N/nfolds < 3) && grouped) {
        warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold", 
            call. = FALSE)
        grouped = FALSE
    }
    if (!is.null(offset)) {
        is.offset = TRUE
        offset = drop(offset)
    }
    else is.offset = FALSE
    mlami = max(sapply(outlist, function(obj) min(obj$lambda)))
    which_lam = lambda >= mlami
    #predmat = matrix(NA, nrow(y), length(lambda))
    predmat1 = matrix(NA, nrow(y), length(lambda))
    predmat2 = matrix(NA, nrow(y), length(lambda))
    predmat3 = matrix(NA, nrow(y), length(lambda))
    nlams = double(nfolds)
    for (i in seq(nfolds)) {
        which = foldid == i
        fitobj = outlist[[i]]
        if (is.offset) 
            off_sub = offset[which]
        preds = predict.glmnetpo(fitobj, x[which, , drop = FALSE], s = lambda[which_lam], newoffset = off_sub, type = "response")$probs
        nlami = sum(which_lam)
        #predmat[which, seq(nlami)] = preds
        
        predmat1[which, seq(nlami)] = preds[,1,seq(nlami)]
		predmat2[which, seq(nlami)] = preds[,2,seq(nlami)]
		predmat3[which, seq(nlami)] = preds[,3,seq(nlami)]
        nlams[i] = nlami
    }
    if (type.measure == "auc") {
        cvraw = matrix(NA, nfolds, length(lambda))
        good = matrix(0, nfolds, length(lambda))
        for (i in seq(nfolds)) {
            good[i, seq(nlams[i])] = 1
            which = foldid == i
            for (j in seq(nlams[i])) {
                cvraw[i, j] = auc.mat(y[which, ], predmat[which, 
                  j], weights[which])
            }
        }
        N = apply(good, 2, sum)
        weights = tapply(weights, foldid, sum)
    }
    else {
        ywt = apply(y, 1, sum)
        y = y/ywt
        weights = weights * ywt
        #N = nrow(y) - apply(is.na(predmat), 2, sum)
        N = nrow(y) - apply(is.na(predmat1), 2, sum)
        cvraw = switch(type.measure, mse = (y[, 1] - (1 - predmat))^2 + 
            (y[, 2] - predmat)^2, mae = abs(y[, 1] - (1 - predmat)) + 
            abs(y[, 2] - predmat), deviance = {
            predmat = pmin(pmax(predmat, prob_min), prob_max)
            lp = y[, 1] * log(1 - predmat) + y[, 2] * log(predmat)
            ly = log(y)
            ly[y == 0] = 0
            ly = drop((y * ly) %*% c(1, 1))
            2 * (ly - lp)
        }, class = y[, 1] * (predmat1 < predmat2) +
        		   y[, 1] * (predmat1 < predmat3) +
        		   y[, 2] * (predmat2 < predmat1) +
        		   y[, 2] * (predmat2 < predmat3) +
        		   y[, 3] * (predmat3 < predmat1) +
        		   y[, 3] * (predmat3 < predmat2))
        #class = y[, 1] * (predmat > 0.5) + y[, 2] * (predmat <= 0.5))
        if (grouped) {
            cvob = cvcompute(cvraw, weights, foldid, nlams)
            cvraw = cvob$cvraw
            weights = cvob$weights
            N = cvob$N
        }
    }
    cvm = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
    cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean, 
        w = weights, na.rm = TRUE)/(N - 1))
    out = list(cvm = cvm, cvsd = cvsd, type.measure = type.measure)
    if (keep) 
        out$fit.preval = predmat
    out
}

#Modified cv.lognet function to handle glmnetpo objects
cv.po.lognet <- function (outlist, lambda, x, y, weights, offset, foldid, type.measure, grouped, keep = FALSE) {
    prob_min = 1e-05
    prob_max = 1 - prob_min
    nc = dim(y)
    if (is.null(nc)) {
        y = as.factor(y)
        ntab = table(y)
        nc = as.integer(length(ntab))
        y = diag(nc)[as.numeric(y), ]
    }
    N = nrow(y)
    nfolds = max(foldid)
    if ((N/nfolds < 10) && type.measure == "auc") {
        warning("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds", call. = FALSE)
        type.measure = cvtype("deviance", "lognet")
    }
    if ((N/nfolds < 3) && grouped) {
        warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold", call. = FALSE)
        grouped = FALSE
    }
    if (!is.null(offset)) {
        is.offset = TRUE
        offset = drop(offset)
    }
    else is.offset = FALSE
    mlami = max(sapply(outlist, function(obj) min(obj$lambda)))
    which_lam = lambda >= mlami
    #predmat = matrix(NA, nrow(y), length(lambda))
    predmat1 = matrix(NA, nrow(y), length(lambda))
    predmat2 = matrix(NA, nrow(y), length(lambda))
    predmat3 = matrix(NA, nrow(y), length(lambda))
    nlams = double(nfolds)
    for (i in seq(nfolds)) {
        which = foldid == i
        fitobj = outlist[[i]]
        if (is.offset) 
            off_sub = offset[which]
        preds = predict.glmnetpo(fitobj, x[which, , drop = FALSE], s = lambda[which_lam], newoffset = off_sub)$probs
        nlami = sum(which_lam)
        #predmat[which, seq(nlami)] = preds
        predmat1[which, seq(nlami)] = preds[,1,seq(nlami)]
		predmat2[which, seq(nlami)] = preds[,2,seq(nlami)]
		predmat3[which, seq(nlami)] = preds[,3,seq(nlami)]
        nlams[i] = nlami
    }
    if (type.measure == "auc") {
        cvraw = matrix(NA, nfolds, length(lambda))
        good = matrix(0, nfolds, length(lambda))
        for (i in seq(nfolds)) {
            good[i, seq(nlams[i])] = 1
            which = foldid == i
            for (j in seq(nlams[i])) {
                cvraw[i, j] = auc.mat(y[which, ], predmat[which, j], weights[which])
            }
        }
        N = apply(good, 2, sum)
        weights = tapply(weights, foldid, sum)
    }
    else {
        ywt = apply(y, 1, sum)
        y = y/ywt
        weights = weights * ywt
        #N = nrow(y) - apply(is.na(predmat), 2, sum)
        N = nrow(y) - apply(is.na(predmat1), 2, sum)
        cvraw = switch(type.measure, 
        	mse = (y[, 1] - (1 - predmat))^2 + (y[, 2] - predmat)^2,
            mae = abs(y[, 1] - (1 - predmat)) + abs(y[, 2] - predmat),
            deviance = {
            	predmat = pmin(pmax(predmat, prob_min), prob_max)
            	lp = y[, 1] * log(1 - predmat) + y[, 2] * log(predmat)
           		ly = log(y)
            	ly[y == 0] = 0
            	ly = drop((y * ly) %*% c(1, 1))
            	2 * (ly - lp)
            	}, 
            class = y[, 1] * (predmat1 < predmat2) +
        		    y[, 1] * (predmat1 < predmat3) +
        		    y[, 2] * (predmat2 < predmat1) +
        		    y[, 2] * (predmat2 < predmat3) +
        		    y[, 3] * (predmat3 < predmat1) +
        		    y[, 3] * (predmat3 < predmat2))
        #class = y[, 1] * (predmat > 0.5) + y[, 2] * (predmat <= 0.5))
        if (grouped) {
            cvob = cvcompute(cvraw, weights, foldid, nlams)
            cvraw = cvob$cvraw
            weights = cvob$weights
            N = cvob$N
        }
    }
    cvm = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
    cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean, w = weights, na.rm = TRUE)/(N - 1))
    out = list(cvm = cvm, cvsd = cvsd, type.measure = type.measure)
    if (keep) 
        out$fit.preval = predmat
    out
}

#Modified cv.lognet function to handle glmnetpo objects
cv.po64.lognet <- function (outlist, lambda, x, y, weights, offset, foldid, type.measure, grouped, keep = FALSE) {
    prob_min = 1e-05
    prob_max = 1 - prob_min
    nc = dim(y)
    if (is.null(nc)) {
        y = as.factor(y)
        ntab = table(y)
        nc = as.integer(length(ntab))
        y = diag(nc)[as.numeric(y), ]
    }
    N = nrow(y)
    nfolds = max(foldid)
    if ((N/nfolds < 10) && type.measure == "auc") {
        warning("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds", call. = FALSE)
        type.measure = cvtype("deviance", "lognet")
    }
    if ((N/nfolds < 3) && grouped) {
        warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold", call. = FALSE)
        grouped = FALSE
    }
    if (!is.null(offset)) {
        is.offset = TRUE
        offset = drop(offset)
    }
    else is.offset = FALSE
    mlami = max(sapply(outlist, function(obj) min(obj$lambda)))
    which_lam = lambda >= mlami
    #predmat = matrix(NA, nrow(y), length(lambda))
    predmat1 = matrix(NA, nrow(y), length(lambda))
    predmat2 = matrix(NA, nrow(y), length(lambda))
    predmat3 = matrix(NA, nrow(y), length(lambda))
    nlams = double(nfolds)
    for (i in seq(nfolds)) {
        which = foldid == i
        fitobj = outlist[[i]]
        if (is.offset) 
            off_sub = offset[which]
        preds = predict.glmnetpo64(fitobj, x[which, , drop = FALSE], s = lambda[which_lam], newoffset = off_sub)$probs
        nlami = sum(which_lam)
        #predmat[which, seq(nlami)] = preds
        predmat1[which, seq(nlami)] = preds[,1,seq(nlami)]
		predmat2[which, seq(nlami)] = preds[,2,seq(nlami)]
		predmat3[which, seq(nlami)] = preds[,3,seq(nlami)]
        nlams[i] = nlami
    }
    if (type.measure == "auc") {
        cvraw = matrix(NA, nfolds, length(lambda))
        good = matrix(0, nfolds, length(lambda))
        for (i in seq(nfolds)) {
            good[i, seq(nlams[i])] = 1
            which = foldid == i
            for (j in seq(nlams[i])) {
                cvraw[i, j] = auc.mat(y[which, ], predmat[which, j], weights[which])
            }
        }
        N = apply(good, 2, sum)
        weights = tapply(weights, foldid, sum)
    }
    else {
        ywt = apply(y, 1, sum)
        y = y/ywt
        weights = weights * ywt
        #N = nrow(y) - apply(is.na(predmat), 2, sum)
        N = nrow(y) - apply(is.na(predmat1), 2, sum)
        cvraw = switch(type.measure, 
        	mse = (y[, 1] - (1 - predmat))^2 + (y[, 2] - predmat)^2,
            mae = abs(y[, 1] - (1 - predmat)) + abs(y[, 2] - predmat),
            deviance = {
            	predmat = pmin(pmax(predmat, prob_min), prob_max)
            	lp = y[, 1] * log(1 - predmat) + y[, 2] * log(predmat)
           		ly = log(y)
            	ly[y == 0] = 0
            	ly = drop((y * ly) %*% c(1, 1))
            	2 * (ly - lp)
            	}, 
            class = y[, 1] * (predmat1 < predmat2) +
        		    y[, 1] * (predmat1 < predmat3) +
        		    y[, 2] * (predmat2 < predmat1) +
        		    y[, 2] * (predmat2 < predmat3) +
        		    y[, 3] * (predmat3 < predmat1) +
        		    y[, 3] * (predmat3 < predmat2))
        #class = y[, 1] * (predmat > 0.5) + y[, 2] * (predmat <= 0.5))
        if (grouped) {
            cvob = cvcompute(cvraw, weights, foldid, nlams)
            cvraw = cvob$cvraw
            weights = cvob$weights
            N = cvob$N
        }
    }
    cvm = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
    cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean, w = weights, na.rm = TRUE)/(N - 1))
    out = list(cvm = cvm, cvsd = cvsd, type.measure = type.measure)
    if (keep) 
        out$fit.preval = predmat
    out
}

#Same as predict.glmnetcr but for glmnetpo 
predict.glmnetpo <- function (object, newx = NULL, s = NULL,...) {
    y <- object$y
    if (is.null(newx)) 
        newx <- object$x
    method <- object$method
    j <- length(unique(y))
    if (class(newx) == "numeric") 
        newx <- matrix(newx, ncol = dim(object$x)[2])
    n <- dim(newx)[1]
    p <- dim(newx)[2]
    y.mat <- matrix(0, nrow = n, ncol = j)
    for (i in 1:n) y.mat[i, y[i]] <- 1
    if(is.null(s) | length(s) > 1) {
		beta.est <- object$beta
		alpha.est <- object$a0
	} else {
		beta.est <- matrix(object$beta[,s])
		row.names(beta.est) <- dimnames(object$beta)[[1]]
		alpha.est <- object$a0[s]
	}	
    k <- apply(beta.est, 2, function(x) sum(x != 0))
    glmnet.BIC <- glmnet.AIC <- numeric()
    pi <- array(NA, dim = c(n, j, dim(beta.est)[2]))
    p.class <- matrix(NA, nrow = n, ncol = dim(beta.est)[2])
    for (i in 1:dim(beta.est)[2]) {
        beta <- beta.est[, i]
        logit <- matrix(rep(0, n * (j - 1)), ncol = j - 1)
        for (h in 1:(j - 1)) {
            cp <- paste("cp", h, sep = "")
            logit[, h] <- alpha.est[i] + beta[names(beta) == cp] + beta[1:p] %*% t(as.matrix(newx))
        }
        delta <- matrix(rep(0, n * (j - 1)), ncol = j - 1)
        for (h in 1:(j - 1)) {
            delta[, h] <- exp(logit[, h])/(1 + exp(logit[, h]))
        }
        minus.delta <- 1 - delta
        if (c("backward", "forward", "cumulative")[charmatch(method, c("backward", "forward", "cumulative"))] == "backward") {
            for (h in j:2) {
                if (h == j) {
                  pi[, h, i] <- delta[, j - 1]
                }
                else if (class(minus.delta[, h:(j - 1)]) == "numeric") {
                  pi[, h, i] <- delta[, h - 1] * minus.delta[, h]
                }
                else if (dim(minus.delta[, h:(j - 1)])[2] >= 2) {
                  pi[, h, i] <- delta[, h - 1] * apply(minus.delta[, h:(j - 1)], 1, prod)
                }
            }
            if (n == 1) 
                pi[, 1, i] <- 1 - sum(pi[, 2:j, i])
            else pi[, 1, i] <- 1 - apply(pi[, 2:j, i], 1, sum)
        }
        if (c("backward", "forward", "cumulative")[charmatch(method, c("backward", "forward", "cumulative"))] == "forward") {
            for (h in 1:(j - 1)) {
                if (h == 1) {
                  pi[, h, i] <- delta[, h]
                }
                else if (h == 2) {
                  pi[, h, i] <- delta[, h] * minus.delta[, h - 1]
                }
                else if (h > 2 && h < j) {
                  pi[, h, i] <- delta[, h] * apply(minus.delta[, 1:(h - 1)], 1, prod)
                }
            }
            if (n == 1) 
                pi[, j, i] <- 1 - sum(pi[, 1:(j - 1), i])
            else pi[, j, i] <- 1 - apply(pi[, 1:(j - 1), i], 1, sum)
        }
        if (c("backward", "forward", "cumulative")[charmatch(method, c("backward", "forward", "cumulative"))] == "cumulative") {
            for (h in 1:(j - 1)) {
                if (h == 1) {
                  pi[, h, i] <- delta[, h]
                }
                else if (h >= 2) {
                  pi[, h, i] <- delta[, h] - delta[,h-1]
                }
            }
            if (n == 1) 
                pi[, j, i] <- 1 - sum(pi[, 1:(j - 1), i])
            else pi[, j, i] <- 1 - apply(pi[, 1:(j - 1), i], 1, sum)
        } 
        LL <- 0
        if (c("backward", "forward", "cumulative")[charmatch(method, c("backward", "forward", "cumulative"))] == "backward") {
            for (h in 1:(j - 1)) {
                if (class(y.mat[, 1:h]) == "matrix") 
                  ylth <- apply(y.mat[, 1:h], 1, sum)
                else ylth <- y.mat[, 1]
                LL <- LL + log(delta[, h]) * y.mat[, h + 1] + log(1 - delta[, h]) * ylth
            }
        }
        if (c("backward", "forward", "cumulative")[charmatch(method, c("backward", "forward", "cumulative"))] == "forward") {
            for (h in 1:(j - 1)) {
                if (class(y.mat[, h:j]) == "matrix") 
                  ygeh <- apply(y.mat[, h:j], 1, sum)
                else ygeh <- y.mat[, j]
                LL <- LL + log(delta[, h]) * y.mat[, h] + log(1 - delta[, h]) * ygeh
            }
        }
        if (c("backward", "forward", "cumulative")[charmatch(method, c("backward", "forward", "cumulative"))] == "cumulative") {
            for (h in 1:j) {
                #if (class(y.mat[, h:j]) == "matrix") 
                #ygeh <- apply(y.mat[, h:j], 1, sum)
                #else ygeh <- y.mat[, j]
                #LL <- LL + log(delta[, h]) * y.mat[, h] + log(1 - delta[, h]) * ygeh
                LL <- LL + log(pi[,h,])*y.mat[,h]
            }
        }
        LL <- sum(LL)
        glmnet.BIC[i] <- -2 * LL + k[i] * log(n)
        glmnet.AIC[i] <- -2 * LL + 2 * k[i]
        if (n == 1) 
            p.class[, i] <- which.max(pi[, , i])
        else p.class[, i] <- apply(pi[, , i], 1, which.max)
    }
    class <- matrix(levels(object$y)[p.class], ncol = ncol(p.class))
    if(is.null(s) | length(s) > 1) {
    	names(glmnet.BIC) <- names(glmnet.AIC) <- dimnames(p.class)[[2]] <- dimnames(pi)[[3]] <- names(object$a0)
    } else {
    	names(glmnet.BIC) <- names(glmnet.AIC) <- dimnames(p.class)[[2]] <- dimnames(pi)[[3]] <- list(names(object$a0)[s])
    }
    dimnames(pi)[[2]] <- levels(object$y)
    list(BIC = glmnet.BIC, AIC = glmnet.AIC, class = class, probs = pi)
}

#Same as predict.glmnetcr but for glmnetpo 
predict.glmnetpo64 <- function (object, newx = NULL, s = NULL,...) {
    y <- object$y
    if (is.null(newx)) 
        newx <- object$x
    method <- object$method
    j <- length(unique(y))
    if (class(newx)[1] == "numeric") {
        #newx <- matrix(newx, ncol = dim(object$x)[2])
        newx <- as(newx, "DgCMatrix")
    }   
    n <- dim(newx)[1]
    p <- dim(newx)[2]
    y.mat <- matrix(0, nrow = n, ncol = j)
    for (i in 1:n) y.mat[i, y[i]] <- 1
    if(is.null(s) | length(s) > 1) {
		beta.est <- object$beta
		alpha.est <- object$a0
	} else {
		beta.est <- matrix(object$beta[,s])
		row.names(beta.est) <- dimnames(object$beta)[[1]]
		alpha.est <- object$a0[s]
	}	
    k <- apply(beta.est, 2, function(x) sum(x != 0))
    glmnet.BIC <- glmnet.AIC <- numeric()
    pi <- array(NA, dim = c(n, j, dim(beta.est)[2]))
    p.class <- matrix(NA, nrow = n, ncol = dim(beta.est)[2])
    for (i in 1:dim(beta.est)[2]) {
        beta <- beta.est[, i]
        logit <- matrix(rep(0, n * (j - 1)), ncol = j - 1)
        for (h in 1:(j - 1)) {
            cp <- paste("cp", h, sep = "")
            #logit[, h] <- alpha.est[i] + beta[names(beta) == cp] + beta[1:p] %*% t(as.matrix(newx))
			      logit[, h] <- as.matrix(alpha.est[i] + beta[names(beta) == cp] + tcrossprod(beta[1:p],newx))     
        }
        delta <- matrix(rep(0, n * (j - 1)), ncol = j - 1)
        for (h in 1:(j - 1)) {
            delta[, h] <- exp(logit[, h])/(1 + exp(logit[, h]))
        }
        minus.delta <- 1 - delta
        if (c("backward", "forward", "cumulative")[charmatch(method, c("backward", "forward", "cumulative"))] == "backward") {
            for (h in j:2) {
                if (h == j) {
                  pi[, h, i] <- delta[, j - 1]
                }
                else if (class(minus.delta[, h:(j - 1)]) == "numeric") {
                  pi[, h, i] <- delta[, h - 1] * minus.delta[, h]
                }
                else if (dim(minus.delta[, h:(j - 1)])[2] >= 2) {
                  pi[, h, i] <- delta[, h - 1] * apply(minus.delta[, h:(j - 1)], 1, prod)
                }
            }
            if (n == 1) 
                pi[, 1, i] <- 1 - sum(pi[, 2:j, i])
            else pi[, 1, i] <- 1 - apply(pi[, 2:j, i], 1, sum)
        }
        if (c("backward", "forward", "cumulative")[charmatch(method, c("backward", "forward", "cumulative"))] == "forward") {
            for (h in 1:(j - 1)) {
                if (h == 1) {
                  pi[, h, i] <- delta[, h]
                }
                else if (h == 2) {
                  pi[, h, i] <- delta[, h] * minus.delta[, h - 1]
                }
                else if (h > 2 && h < j) {
                  pi[, h, i] <- delta[, h] * apply(minus.delta[, 1:(h - 1)], 1, prod)
                }
            }
            if (n == 1) 
                pi[, j, i] <- 1 - sum(pi[, 1:(j - 1), i])
            else pi[, j, i] <- 1 - apply(pi[, 1:(j - 1), i], 1, sum)
        }
        if (c("backward", "forward", "cumulative")[charmatch(method, c("backward", "forward", "cumulative"))] == "cumulative") {
            for (h in 1:(j - 1)) {
                if (h == 1) {
                  pi[, h, i] <- delta[, h]
                }
                else if (h >= 2) {
                  pi[, h, i] <- delta[, h] - delta[,h-1]
                }
            }
            if (n == 1) 
                pi[, j, i] <- 1 - sum(pi[, 1:(j - 1), i])
            else pi[, j, i] <- 1 - apply(pi[, 1:(j - 1), i], 1, sum)
        } 
        LL <- 0
        if (c("backward", "forward", "cumulative")[charmatch(method, c("backward", "forward", "cumulative"))] == "backward") {
            for (h in 1:(j - 1)) {
                if (class(y.mat[, 1:h]) == "matrix") 
                  ylth <- apply(y.mat[, 1:h], 1, sum)
                else ylth <- y.mat[, 1]
                LL <- LL + log(delta[, h]) * y.mat[, h + 1] + log(1 - delta[, h]) * ylth
            }
        }
        if (c("backward", "forward", "cumulative")[charmatch(method, c("backward", "forward", "cumulative"))] == "forward") {
            for (h in 1:(j - 1)) {
                if (class(y.mat[, h:j]) == "matrix") 
                  ygeh <- apply(y.mat[, h:j], 1, sum)
                else ygeh <- y.mat[, j]
                LL <- LL + log(delta[, h]) * y.mat[, h] + log(1 - delta[, h]) * ygeh
            }
        }
        if (c("backward", "forward", "cumulative")[charmatch(method, c("backward", "forward", "cumulative"))] == "cumulative") {
            for (h in 1:j) {
                #if (class(y.mat[, h:j]) == "matrix") 
                #ygeh <- apply(y.mat[, h:j], 1, sum)
                #else ygeh <- y.mat[, j]
                #LL <- LL + log(delta[, h]) * y.mat[, h] + log(1 - delta[, h]) * ygeh
                LL <- LL + log(pi[,h,])*y.mat[,h]
            }
        }
        LL <- sum(LL)
        glmnet.BIC[i] <- -2 * LL + k[i] * log(n)
        glmnet.AIC[i] <- -2 * LL + 2 * k[i]
        if (n == 1) 
            p.class[, i] <- which.max(pi[, , i])
        else p.class[, i] <- apply(pi[, , i], 1, which.max)
    }
    class <- matrix(levels(object$y)[p.class], ncol = ncol(p.class))
    if(is.null(s) | length(s) > 1) {
    	names(glmnet.BIC) <- names(glmnet.AIC) <- dimnames(p.class)[[2]] <- dimnames(pi)[[3]] <- names(object$a0)
    } else {
    	names(glmnet.BIC) <- names(glmnet.AIC) <- dimnames(p.class)[[2]] <- dimnames(pi)[[3]] <- list(names(object$a0)[s])
    }
    dimnames(pi)[[2]] <- levels(object$y)
    list(BIC = glmnet.BIC, AIC = glmnet.AIC, class = class, probs = pi)
}

#Error bar function from glmnet for plotting cv results
error.bars <-	function(x, upper, lower, width = 0.02, ...) {
	xlim <- range(x)
	barw <- diff(xlim) * width
	segments(x, upper, x, lower, ...)
	segments(x - barw, upper, x + barw, upper, ...)
	segments(x - barw, lower, x + barw, lower, ...)
	range(upper, lower)
}

#Cross validation plotting with log10 base
plot.cv.glmnet2 <- function (x, sign.lambda = -1, ...) {
    cvobj = x
    xlab = "log10(Lambda)"
    if (sign.lambda < 0) 
        xlab = paste("-", xlab, sep = "")
    plot.args = list(x = sign.lambda * log10(cvobj$lambda), y = cvobj$cvm, 
        ylim = range(cvobj$cvup, cvobj$cvlo), xlab = xlab, ylab = cvobj$name, 
        type = "n")
    new.args = list(...)
    if (length(new.args)) 
        plot.args[names(new.args)] = new.args
    do.call("plot", plot.args)
    error.bars(sign.lambda * log10(cvobj$lambda), cvobj$cvup, cvobj$cvlo, width = 0.01, col = "darkgrey")
    points(sign.lambda * log10(cvobj$lambda), cvobj$cvm, pch = 20, col = "red")
    axis(side = 3, at = sign.lambda * log10(cvobj$lambda), labels = paste(cvobj$nz), tick = FALSE, line = 0)
    abline(v = sign.lambda * log10(cvobj$lambda.min), lty = 3)
    abline(v = sign.lambda * log10(cvobj$lambda.1se), lty = 3)
    invisible()
}




#Shortcut function for quickly calculating likelihoods using subsets of model parameters
#Set up a NEWX
JOINT.FULL.LL.CALC <- function() {
  for (h in 1:(JOINT.PE.FULL.LEVELS - 1)) {
    cp <- paste("cp", h, sep = "")
    JOINT.PE.FULL.LOGIT[, h] <- JOINT.PE.FULL.ALPHA + JOINT.PE.FULL.THRESH[names(JOINT.PE.FULL.THRESH) == cp] + JOINT.PE.FULL.BETA %*% t(NEWX)
  }
  
  JOINT.PE.FULL.DELTA <- matrix(rep(0, N.OBS * (JOINT.PE.FULL.LEVELS - 1)), ncol = JOINT.PE.FULL.LEVELS - 1)
  for (h in 1:(JOINT.PE.FULL.LEVELS - 1)) {
    JOINT.PE.FULL.DELTA[, h] <- exp(JOINT.PE.FULL.LOGIT[, h])/(1 + exp(JOINT.PE.FULL.LOGIT[, h]))
  }
  
  JOINT.PE.FULL.PI[, 1] <- JOINT.PE.FULL.DELTA[, 1]
  JOINT.PE.FULL.PI[, 2] <- JOINT.PE.FULL.DELTA[, 2] - JOINT.PE.FULL.DELTA[,1]
  JOINT.PE.FULL.PI[, 3] <- 1 - JOINT.PE.FULL.DELTA[,2]
  
  LL.TEST <- 0
  for (h in 1:JOINT.PE.FULL.LEVELS) {
    LL.TEST <- LL.TEST + log(JOINT.PE.FULL.PI[,h])*JOINT.PE.FULL.Y.MAT[,h]
  }
  return(sum(LL.TEST))
}

###Counts the number of identical shortest paths between two nodes on different networks
path.identity <- function(NETWORK.1,NETWORK.2,FROM,TO) {
  PATHS.1 <- lapply(FROM,FUN=function(a) {
    lapply(TO,FUN=function(b) {
      SHORT.PATHS.1 <- all_shortest_paths(NETWORK.1,from=a,to=b)$res
      sapply(SHORT.PATHS.1,FUN=names)
    })
  })
  PATHS.2 <- lapply(FROM,FUN=function(c) {
    lapply(TO,FUN=function(d) {
      SHORT.PATHS.2 <- all_shortest_paths(NETWORK.2,from=c,to=d)$res
      sapply(SHORT.PATHS.2,FUN=names)
    })
  })
  COMBINE <- mapply(FUN=function(x,y) {
    if(nrow(x) != nrow(y)) {
      return(0)
    }else {
      cbind(x,y)
    }
  },unlist(PATHS.1,recursive=FALSE),unlist(PATHS.2,recursive = FALSE),SIMPLIFY = FALSE)

  OVERLAP <- sapply(COMBINE,FUN=function(x) {
    sum(duplicated(t(x))) 
  })
  return(OVERLAP)
}

###Count number of shortest paths between two nodes in a network
num.paths <- function(NETWORK,FROM,TO) {
  sapply(FROM,FUN=function(a) {
    sapply(TO,FUN=function(b) {
      NUM.SHORT.PATHS <- length(all_shortest_paths(NETWORK,from=a,to=b)$res)
      return(NUM.SHORT.PATHS)
    })
  })
}

##Use for number of hamming distance paths
num.paths.HD <- function(FROM,TO) {
  sapply(FROM,FUN=function(x) {
    sapply(TO,FUN=function(y) {
      factorial(get.HAM.distance(x,y))
    })
  })
}

##Use for number of genetic code paths
num.paths.GC <- function(FROM,TO,CODE="S") {
  sapply(FROM,FUN=function(x) {
    sapply(TO,FUN=function(y) {
      DIST <- get.distance(x,y,code = CODE,vectorized = TRUE)
      ROUTES <- get.routes(x,y,code = CODE,vectorized = TRUE)
      prod(ROUTES)*factorial(sum(DIST))/prod(factorial(DIST))
    })
  })
}

##Functions for simulating across network
get.genotype.class <- function(GENOTYPE,NETWORKTYPE) {
  if(NETWORKTYPE == "T") {
    DT.JOINT[AAseq == paste0("E",GENOTYPE),TE.JOINT.CLASS]
  }else if(NETWORKTYPE == "P") {
    DT.JOINT[AAseq == paste0("E",GENOTYPE),PE.JOINT.CLASS]
  }else if(NETWORKTYPE == "M") {
    DT.JOINT[AAseq == paste0("E",GENOTYPE),ME.JOINT.CLASS]
  }
}

ms.simulate <- function(STARTING.GENOS,NETWORK,NETWORKTYPE,REPS) {
  I <- 1:REPS
  J <- 1:length(STARTING.GENOS)
  OUT.N <- matrix(0,nrow=length(I),ncol=length(J))
  OUT.GENO <- matrix(0,nrow=length(I),ncol=length(J))
  colnames(OUT.N) <- STARTING.GENOS
  colnames(OUT.GENO) <- STARTING.GENOS
  G.LIST <- list()
  for(j in J) {
    START.G <- STARTING.GENOS[j]
    G.LIST[[j]] <- list()
    for(i in I) {
      N <- 0
      G.LIST[[j]][[i]] <- START.G
      CURRENT.G <- START.G
      while(-1 < N) {
        if(get.genotype.class(CURRENT.G,NETWORKTYPE) == "SRE-specific" | length(neighbors(NETWORK,CURRENT.G)) == 0) {
          OUT.N[i,j] <- N
          OUT.GENO[i,j] <- CURRENT.G
          break()
        } else {
          N <- N + 1
          NEIGHBORS.G <- names(neighbors(NETWORK,CURRENT.G))
          CURRENT.G <- sample(NEIGHBORS.G,1)
          G.LIST[[j]][[i]] <- c(G.LIST[[j]][[i]],CURRENT.G) 
        }
      }
    }
  }
  names(G.LIST) <- STARTING.GENOS
  BIND <- list("Number.Steps"=OUT.N,"End.Genotype"=OUT.GENO,"Path"=G.LIST)
  return(BIND)
}

sel.simulate <- function(STARTING.GENOS,NETWORK,NETWORKTYPE,REPS,POP.SIZE=1000,CUTOFF=10000) {
  POP.SIZE <- POP.SIZE
  I <- 1:REPS
  J <- 1:length(STARTING.GENOS)
  OUT.N <- matrix(0,nrow=length(I),ncol=length(J))
  OUT.GENO <- matrix(0,nrow=length(I),ncol=length(J))
  colnames(OUT.N) <- STARTING.GENOS
  colnames(OUT.GENO) <- STARTING.GENOS
  G.LIST <- list()
  for(j in J) {
    START.G <- STARTING.GENOS[j]
    G.LIST[[j]] <- list()
    for(i in I) {
      N <- 0
      G.LIST[[j]][[i]] <- START.G
      CURRENT.G <- START.G
      while(-1 < N) {
        if(get.genotype.class(CURRENT.G,NETWORKTYPE) == "SRE-specific" | length(neighbors(NETWORK,CURRENT.G)) == 0 | N > CUTOFF) {
          OUT.N[i,j] <- N
          OUT.GENO[i,j] <- CURRENT.G
          break()
        } else {
          N <- N + 1
          NEIGHBORS.G <- names(neighbors(NETWORK,CURRENT.G))
          
          #Get phenotypes
          if(NETWORKTYPE == "T") {
            CURRENT.G.ERE.PHENOTYPE <- DT.JOINT[which(DT.JOINT$AAseq %in% paste0("E",CURRENT.G)),"PRED.TE.LINK"]
            CURRENT.G.SRE.PHENOTYPE <- DT.JOINT[which(DT.JOINT$AAseq %in% paste0("S",CURRENT.G)),"PRED.TE.LINK"]
            NEIGHBORS.ERE.PHENOTYPE <- DT.JOINT[which(DT.JOINT$AAseq %in% paste0("E",NEIGHBORS.G)),"PRED.TE.LINK"]
            NEIGHBORS.SRE.PHENOTYPE <- DT.JOINT[which(DT.JOINT$AAseq %in% paste0("S",NEIGHBORS.G)),"PRED.TE.LINK"]
          }else if (NETWORKTYPE == "P") {
            CURRENT.G.ERE.PHENOTYPE <- DT.JOINT[which(DT.JOINT$AAseq %in% paste0("E",CURRENT.G)),"PRED.PE.LINK"]
            CURRENT.G.SRE.PHENOTYPE <- DT.JOINT[which(DT.JOINT$AAseq %in% paste0("S",CURRENT.G)),"PRED.PE.LINK"]
            NEIGHBORS.ERE.PHENOTYPE <- DT.JOINT[which(DT.JOINT$AAseq %in% paste0("E",NEIGHBORS.G)),"PRED.PE.LINK"]
            NEIGHBORS.SRE.PHENOTYPE <- DT.JOINT[which(DT.JOINT$AAseq %in% paste0("S",NEIGHBORS.G)),"PRED.PE.LINK"]
          }else if (NETWORKTYPE == "M") {
            CURRENT.G.ERE.PHENOTYPE <- DT.JOINT[which(DT.JOINT$AAseq %in% paste0("E",CURRENT.G)),"PRED.ME.LINK"]
            CURRENT.G.SRE.PHENOTYPE <- DT.JOINT[which(DT.JOINT$AAseq %in% paste0("S",CURRENT.G)),"PRED.ME.LINK"]
            NEIGHBORS.ERE.PHENOTYPE <- DT.JOINT[which(DT.JOINT$AAseq %in% paste0("E",NEIGHBORS.G)),"PRED.ME.LINK"]
            NEIGHBORS.SRE.PHENOTYPE <- DT.JOINT[which(DT.JOINT$AAseq %in% paste0("S",NEIGHBORS.G)),"PRED.ME.LINK"]
          }
          #Get specificity
          CURRENT.SPEC <- CURRENT.G.ERE.PHENOTYPE - CURRENT.G.SRE.PHENOTYPE
          NEIGHBORS.SPEC <- NEIGHBORS.ERE.PHENOTYPE - NEIGHBORS.SRE.PHENOTYPE
          
          #Convert to selection coefficient
          S <- unlist(NEIGHBORS.SPEC) - unlist(CURRENT.SPEC)
          #Scale set such that a population size of 100 gives NS ~ 1 for a specificity difference of 3 (difference between null and weak thresholds)
          S.SCALED <- S/300
          
          #Get fixation probabilities
          FIX.PROBABILITY <- (1 - exp(-2*S.SCALED))/(1 - exp(-4*POP.SIZE*S.SCALED))
          RELATIVE.PROBABILITY <- FIX.PROBABILITY/sum(FIX.PROBABILITY)
          
          CURRENT.G <- sample(NEIGHBORS.G,1,prob=RELATIVE.PROBABILITY)
          G.LIST[[j]][[i]] <- c(G.LIST[[j]][[i]],CURRENT.G) 
        }
      }
    }
  }
  names(G.LIST) <- STARTING.GENOS
  BIND <- list("Number.Steps"=OUT.N,"End.Genotype"=OUT.GENO,"Path"=G.LIST)
  return(BIND)
}

##Generic error bar function
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

##Permutation tests
#Matrices of different sizes
permutation.test.matrix <- function(M1,M2,REP) {
  I <- 1:REP
  TEST.STAT <- numeric(length(I))
  for(i in I) {
    NR1 <- nrow(M1); NR2 <- nrow(M2); NC1 <- ncol(M1); NC2 <- ncol(M2)
    NROW <- NR1 + NR2
    NCOL <- NC1 + NC2
    
    ROW.SAMPLE <- sample(1:NROW,NR1,replace=FALSE)
    COL.SAMPLE <- sample(1:NCOL,NC1,replace=FALSE)
    
    S1.ROW.INDEX <- ROW.SAMPLE[which(ROW.SAMPLE <= NR1)]
    S1.COL.INDEX <- COL.SAMPLE[which(COL.SAMPLE <= NC1)]
    S2.ROW.INDEX <- ROW.SAMPLE[which((ROW.SAMPLE - NR1) > 0)]-NR1
    S2.COL.INDEX <- COL.SAMPLE[which((COL.SAMPLE - NC1) > 0)]-NC1
    
    M1.1 <- M1[ S1.ROW.INDEX, S1.COL.INDEX]
    M1.2 <- M2[ S2.ROW.INDEX, S2.COL.INDEX]
    M2.1 <- M1[-S1.ROW.INDEX,-S1.COL.INDEX]
    M2.2 <- M2[-S2.ROW.INDEX,-S2.COL.INDEX]
    
    TEST.STAT[i] <-   (sum(M1.1,na.rm=TRUE) + sum(M1.2,na.rm=TRUE))/(NR1*NC1)  - (sum(M2.1,na.rm=TRUE) + sum(M2.2,na.rm=TRUE))/(NR2*NC2)
  }
  return(TEST.STAT)
}  

#vectors of different lengths
permutation.test.vector <- function(V1,V2,REP) {
  I <- 1:REP
  TEST.STAT <- numeric(length(I))
  MASTER.VECTOR <- c(V1[is.finite(V1)],V2[is.finite(V2)])
  
  for(i in I) {
    SAMPLE <- sample(1:length(MASTER.VECTOR),size=length(V1),replace=FALSE)
    V1 <- MASTER.VECTOR[SAMPLE]
    V2 <- MASTER.VECTOR[-SAMPLE]
    
    TEST.STAT[i] <- mean(V1,na.rm=TRUE) - mean(V2,na.rm=TRUE)
  }
  return(TEST.STAT)
}


###Calculate how distinct sets of paths between genotypes are 
distinct.paths <- function(NETWORK,START,END) {
  I <- 1:length(START)
  J <- 1:length(END)
  DISTINCT.MAT <- matrix(0,nrow=length(I),ncol=length(J)); rownames(DISTINCT.MAT) <- START; colnames(DISTINCT.MAT) <- END
  
  for(i in I) {
    for(j in J) {
      SHORTEST.PATHS <- all_shortest_paths(NETWORK,from=START[i],to=END[j])$res
      if(length(SHORTEST.PATHS) == 0) {
        DISTINCT.MAT[i,j] <- NA
      } else if (length(SHORTEST.PATHS) == 1) {
        DISTINCT.MAT[i,j] <- 1
      } else {
        NAMES <- t(sapply(SHORTEST.PATHS,FUN=names))
        N.STEPS <- ncol(NAMES) - 1
        K <- 1:(ncol(NAMES)-1)
        NI.EDGES <- numeric(length(K))
        for(k in K) {
          EDGES <- unique(NAMES[,k:(k+1)])
          NI.EDGES[k] <- nrow(EDGES)
        }  
        N.GENOS <- length(unique(c(NAMES)))
        N.EDGES <- sum(NI.EDGES) 
        DISTINCT.MAT[i,j] <- (((N.GENOS-2)/(N.STEPS-1))^2)/(N.EDGES/N.STEPS)
      }  
    }
  }
  return(DISTINCT.MAT)
}

#Distinct paths calculation for NG network
distinct.paths.NG <- function(START,END) {
  I <- 1:length(START)
  J <- 1:length(END)
  DISTINCT.MAT <- matrix(0,nrow=length(I),ncol=length(J)); rownames(DISTINCT.MAT) <- START; colnames(DISTINCT.MAT) <- END
  
  for(i in I) {
    for(j in J) {
      #Extract subgraph with relevant genotypes to reduce computational time
      #Split start and ending genotype into individual amino acids
      START.SPLIT <- unlist(strsplit(START[i],""))
      END.SPLIT   <- unlist(strsplit(END[j],""))
      
      #Get intermediate amino acids between start and end amino acid for each site
      K <- 1:4
      SPLIT <- list()
      for(k in K) {
        SPLIT[k] <- AA.INT.Z[which(rownames(AA.INT.Z)==START.SPLIT[k]),which(colnames(AA.INT.Z)==END.SPLIT[k])]
      }  
      SPLIT <- sapply(SPLIT,strsplit,"")
      
      #Make combinations of amino acids at different sites
      GENO.COMBOS <- as.matrix(expand.grid(SPLIT[[1]],SPLIT[[2]],SPLIT[[3]],SPLIT[[4]],stringsAsFactors = FALSE))
      GENO.COMBOS <- apply(GENO.COMBOS,1,paste0,collapse="")
      
      #Get subgraph of only visited genotypes
      SUB <- induced_subgraph(NG.ACT.NET,GENO.COMBOS)
      
      DISTINCT.MAT[i,j] <- distinct.paths(SUB,START[i],END[j])
    }
  }
  return(DISTINCT.MAT)
}  
