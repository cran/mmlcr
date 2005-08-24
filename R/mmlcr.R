###################################################
### chunk number 1: 
###################################################
"mmlcr"<- function(object = NULL, ...) UseMethod("mmlcr")





###################################################
### chunk number 6: mmlcrEval
###################################################
"mmlcr.mmlcr"<-
function(object, max.iter = 50, trace = TRUE, tol = 0.005, ...)
{
	require(nnet)
	require(nlme)
	require(survival) # needed for censored normal functions
	
	# Put object components in more convenient form:
	components <- object$components
	c.length <- length(components)
	prior.prob <- object$prior.prob
	post.prob <- object$post.prob
	gamma.matrix <- object$gamma.matrix
	outer <- object$outer
	loglikelihood <- object$loglikelihood


	n.classes <- dim(post.prob)[2]
	n.id <- dim(post.prob)[1]
	
	# Turn down warnings, but restore them on exiting:
#	warn.level <- options()$warn
#	options(warn = -1)
#	on.exit(options(warn = warn.level))
	
	# Build a big version of the outer dataframe, 
	# with each observation replacated by the number of latent classes:
	latent.class <- as.factor(matrix(rep(1:n.classes, dim(object$outer.df)[1]), 
						ncol = n.classes, byrow = TRUE))
	big.df <- object$outer.df
	for(k in 2:n.classes)
		big.df <- rbind(big.df, object$outer.df)
	big.df <- data.frame(latent.class, big.df)
	
		
		
	# Get counters and convergence criteria initialized:	
	loglikegoal <- NA
	loglike2 <- NA
	iter <- 0
	attr(object, "warn") <- NULL
	multinom.result <- NULL
	# Begin main loop:
	while(ifelse(is.na(loglikegoal - loglikelihood), TRUE, !(abs(loglikegoal - loglikelihood) < tol))) {
		if(iter >= max.iter) {
#			options(warn = warn.level)
			warning("Maximum iterations reached")
			attr(object, "warn") <- "Maximum iterations reached"
			break
		}
		iter <- iter + 1
		loglike3 <- loglike2
		loglike2 <- loglikelihood
		
		## Fit each component ##
		for(k in 1:n.classes) {
			for(i in 1:c.length)
				components[[i]] <- mmlcrfit(components[[i]], 
				  weights = post.prob, classnumber = k)
		}#
		## Calculate likelihoods from each component and then calculate posterior probability ##
		temp.like <- matrix(rep(1, n.classes * n.id), ncol = n.classes)
		for(i in 1:c.length)
			for(k in 1:n.classes)
				temp.like[, k] <- temp.like[, k] * 
					mmlcrlike(components[[i]],
					weights = post.prob, 
					classnumber = k)
		post.prob <- prior.prob * temp.like

		# Make sure rows add up to 1:
		post.prob <- as.matrix(post.prob/(apply(post.prob, 1, sum)))
			
		# If the likelihood steps failed for some observations, default to
		# prior prob for those.
		post.prob[is.na(post.prob)] <- prior.prob[is.na(post.prob)]
		post.prob <- data.frame(post.prob, 
			row.names = row.names(object$outer.df))

		
		## Fit covariate part using multinom ##
		big.df$.posterior2 <- unlist(post.prob, use.names = FALSE)
		if (n.classes > 1) {
			if (is.null(multinom.result))
				multinom.result <- multinom(outer, data = big.df, weights = 
				.posterior2, trace = FALSE, maxit = 250)
			else {
				multinom.result <- multinom(outer, data = big.df, weights = 
					.posterior2, trace = FALSE, maxit = 10, Wts = multinom.result$wts)
			}
			multinom.result$weights[is.na(multinom.result$weights)] <- 0
		}#
		
		## Finish up loop ##
		if (n.classes > 1) {
			gamma.matrix <- rbind(0, coef(multinom.result))
			prior.prob <- data.frame(exp(model.matrix(outer, big.df)[1:n.id,
				] %*% t(gamma.matrix)))
			prior.prob <- prior.prob/apply(prior.prob, 1, sum)
			prior.prob <- apply(prior.prob, 2, pmax, .Machine$double.xmin)
			temp.like <- apply(temp.like, 2, pmax, .Machine$double.xmin)
		} 
		else{ 
			prior.prob <- matrix(1, ncol = 1, nrow = n.id)
		}
		
		loglikelihood <- sum(post.prob * (log(prior.prob) + log(temp.like)))
		convergence.index <- (loglikelihood - loglike2)/(loglike2 - loglike3)
		loglikegoal <- loglike2 + 1/(1 - convergence.index) * (loglikelihood - loglike2)
		if(trace) {
			if(iter == 1)
				cat("loglike\t conv. index\t ll goal\t Class Percentages:\n"
				  )
			group.pct <- apply(post.prob, 2, sum) /  dim(post.prob)[1]
			cat(paste(format(round(loglikelihood, 2), nsmall=2 ), 
				"\t", format(round(convergence.index, 2), nsmall=2), 
				"\t", format(round(loglikegoal, 2),nsmall=2), 
				"\t"), format(round(100 * unlist(group.pct, use.names = FALSE), 1), nsmall=1), 
				"\n")
		}
		if(is.na(loglikelihood) || is.infinite(loglikelihood) 
			|| is.nan(loglikelihood) || is.nan(loglikegoal))
			break
	} # End main loop
	for(i in 1:c.length)
		for(k in 1:n.classes)
			components[[i]] <- mmlcrlike(components[[i]],
				weights = post.prob, classnumber = k,
				final = TRUE)
	# n is the number of observations, m the number of parameters
	n <- 0
	m <- length(gamma.matrix[-1, ])
	for(i in 1:c.length) {
		n <- n + components[[i]]$n
		m <- m + components[[i]]$m
	}
	if (n.classes > 1) {
		dimnames(gamma.matrix) <- list(multinom.result$lev, multinom.result$coefnames)
	} else {

			temp.model.matrix <- model.matrix(outer, big.df)
			gamma.matrix <- matrix(0, ncol = ncol(temp.model.matrix), nrow = 1)
			dimnames(gamma.matrix) <- list("1", dimnames(temp.model.matrix)[[2]])
	}
	object$prior.prob <- data.frame(prior.prob, row.names = row.names(object$outer.df))
	names(object$prior.prob) <- paste("PriorProb", 1:n.classes, sep="")
	object$post.prob <- post.prob
	names(object$post.prob) <- paste("PostProb", 1:n.classes, sep="")
	object$loglikelihood <- loglikelihood
	object$components <- components
	object$gamma.matrix <- gamma.matrix
	object$BIC <- -2 * loglikelihood + m * log(n)
	object$resid.df <- n - m
	object$df <- m
	object
}


###################################################
### chunk number 7: 
###################################################
"mmlcr.default"<-
function(object = NULL, outer, components, 
	data = error("data must be given (as a data.frame)"), subset, 
	n.groups = 2, prior.prob = NULL, post.prob = NULL, no.fit = FALSE, 
	max.iter = 50, trace = TRUE, tol = 0.005, ...)
{
	require(nnet)
	require(nlme)
	require(survival) # needed for censored normal functions

	call <- match.call()
	m <- match.call(expand.dots = FALSE)
	m$outer <- m$components <- m$n.groups <- m$prior.prob <- m$post.prob <- 
		m$no.fit <- m$max.iter <- m$trace <- m$tol <- NULL
	m$drop.unused.levels <- TRUE
	m.components <- m
	m$formula <- asOneFormula(outer)
	m[[1]] <- as.name("model.frame")
	m <- eval(m, data)
	groups <- getGroups(m, outer)
	original.group.order <- order(groups)
	m <- m[original.group.order,  , drop = FALSE]
	groups <- as.vector(groups[original.group.order])
	gunique <- unique(groups)
	firstInGroup <- match(gunique, groups)
	m <- m[firstInGroup,  , drop = FALSE]
	
	# Make data.frame to be used as components:
	component.vars <- all.vars(outer)
	for (j in 1:length(components))
		component.vars <- c( component.vars, all.vars(components[[j]]$formula))
	m.components$formula <-  eval(parse(text = 
		paste("~", paste(component.vars, collapse = "+")))[[1]])
	if (is.null(m.components$na.action)) m.components$na.action <- na.pass
	m.components[[1]] <- as.name("model.frame")
	m.components <- eval(m.components, data)
	m.components <- m.components[original.group.order,  ]
	
	x <- list()
	x$call <- call
	x$components <- list()
	x$outer <- formula(paste("latent.class", paste(splitFormula(outer, sep
		 = "|")[[1]], collapse = " ")))
	x$outer.df <- data.frame(m, row.names = groups[firstInGroup])
	x$loglikelihood <- NA
	if(is.null(prior.prob))
		prior.prob <- 
			matrix(1/n.groups, ncol = n.groups, nrow = dim(x$outer.df)[1])
	x$prior.prob <- data.frame(prior.prob, row.names = groups[firstInGroup])
	names(x$prior.prob) <- paste("PriorProb", 1:n.groups, sep="")
	
	if(is.null(post.prob)) {
		post.prob <- matrix(nrow = dim(x$outer.df)[1], ncol = n.groups)
		for(k in 1:n.groups)
			post.prob[, k] <- sample(c(0.001, 1), dim(post.prob)[1],
				replace = TRUE, prob = c(0.95, 0.5))
	}
	x$post.prob <- data.frame(post.prob, row.names = groups[firstInGroup])
	names(x$post.prob) <- paste("PostProb", 1:n.groups, sep="")
	
	.grouping <- deparse(splitFormula(call$outer, sep = "|")[[2]][[2]])
	
	
	# Here we initialize each component. It's not necessary
	# to fully fit each one.
	for(j in 1:length(components)) {
		attr(components[[j]], "class") <- components[[j]]$class
		x$components[[j]] <- mmlcrcomponentinit(object = components[[j]], 
			n.groups = n.groups,  prob = x$post.prob, data = m.components,
			grouping = .grouping)
	}
	attr(x, "class") <- "mmlcr"
	if(no.fit)
		x
	else {
		mmlcr.mmlcr(x, max.iter = max.iter, tol = tol, trace = trace)
	}
}


###################################################
### chunk number 8: 
###################################################



"mmlcrcomponentinit"<-
function(object, n.groups, prob, data, grouping)
UseMethod("mmlcrcomponentinit")

"mmlcrcomponentinit.default" <- 
function(object, n.groups, prob, data, grouping)
{
	stop(paste("No mmlcrcomponentinit method implemented for class",
		class(object)))
}


"mmlcrfit"<-
function(object, weights, classnumber)
UseMethod("mmlcrfit")

"mmlcrfit.default" <- 
function(object, weights, classnumber)
{
	warning("mmlcrfit.nofit being used for mmlcrfit.default\n")
	object
}

"mmlcrlike"<-
function(object, weights, classnumber, final = FALSE)
UseMethod("mmlcrlike")

"mmlcrlike.default" <- 
function(object, weights, classnumber, final)
{
	stop("Missing class for mmlcrlike\n")
}



###################################################
### chunk number 9: 
###################################################
"mmlcrcomponentinit.normlong"<-
function(object, n.groups, prob, data, grouping)
{
	m <- match.call(expand = FALSE)
	m$object <- m$n.groups <- m$prob <- m$... <- NULL
	m$drop.unused.levels <- TRUE
	m$grouping <- NULL
	longformula <- formula(paste(object$formula[[2]], paste(
		unclass(object$formula), grouping, sep = "|")[[3]],
		sep = "~"))
	m$formula <- asOneFormula(longformula)
	m$na.action <- na.exclude
	m[[1]] <- as.name("model.frame")
	m <- eval.parent(m)
	attr(m, "terms") <- terms(longformula)
	attr(m, "formula") <- NULL
	.form <- longformula
	x <- list(data = m, coef = list(), formula = .form)
	x$shortform <- object$formula
	lm.result <- lm(formula(x$shortform), m, na.action = na.exclude)
	for(j in 1:n.groups) {
		x$coef[[j]] <- list(beta = coef(lm.result),
			sigma2 = sum(lm.result$residuals^2) / lm.result$df.resid)
	}
	x$n <- sum(!is.na(getResponse(m, .form)))
	x$m <- (length(x$coef[[1]]$beta) + 1) * n.groups
	attr(x, "class") <- c(attr(object, "class"), "mmlcrlong")
	x
}

"mmlcrfit.normlong"<-
function(object, weights, classnumber)
{
	data <- object$data
	form.full <- formula(attr(data, "terms"))
	form <- object$shortform
	mf <- model.frame(form, data, na.action = na.exclude)
	x <- model.matrix(attr(mf,"terms"), mf)
	y <- model.response(mf, "numeric")
	offset <- model.offset(mf)
	.wts <- weights[match(getGroups(data, form.full, 1), row.names(
		weights)), classnumber]
	var.min <- max(max(sapply(
		object$coef, function(x) x$sigma2)[ - classnumber])/100, 
		.Machine$double.xmin, na.rm=TRUE)
	var.max <- min(min(sapply(
		object$coef, function(x) x$sigma2)[ - classnumber]) * 100,
		.Machine$double.xmax, na.rm=TRUE)
	fit <- lm.wfit(x, y, .wts, offset)
#	fit <- lm(form, data = data, weights = .wts, na.action = na.exclude)
	temp.var <- (sum(.wts * fit$residuals^2) 
		/ (sum(.wts) - length(coef(fit))))
	object$coef[[classnumber]]$beta <- coef(fit)
	object$coef[[classnumber]]$sigma2 <- ifelse(temp.var < var.min, var.min,
		ifelse(temp.var > var.max, var.max, temp.var))
	object
}

"mmlcrlike.normlong" <- 
function(object, weights, classnumber, final = FALSE)
{
	ids <- row.names(weights)
	form.full <- formula(attr(object$data, "terms"))
	groups <- getGroups(object$data, form.full)
	include.vec <- is.element(groups, ids)
	data <- object$data[include.vec,  ]
	form <- object$shortform
	coef <- object$coef[[classnumber]]
	beta <- coef$beta
	sigma <- sqrt(coef$sigma2)
	model.mat <- model.matrix.lm(form, data)
	mu <- model.mat %*% as.matrix(beta)
	if(!final) {
		temp.like <- apply(matrix(getResponse(data), ncol = 1),
			2, dnorm, mean = mu, sd = sigma)
		temp.like[is.na(temp.like)] <- 1
		prod.like <- tapply(temp.like, groups, prod)
		prod.like[match(ids, names(prod.like))]
	}
	else {
		
		form.full <- formula(object$data)
		.wts <- weights[match(groups, row.names(weights)), classnumber]
		if(classnumber == 1) {
			object$fitted <- rep(NA, nrow(object$data))
			object$fitted[include.vec] <- mu * .wts
		}
		else object$fitted[include.vec] <- object$fitted[
				include.vec] + mu * .wts
		object$residuals <- rep(NA, nrow(object$data))
		object$residuals[include.vec] <- getResponse(data) -
			object$fitted[include.vec]
		object
	}
}

"summary.normlong" <- 
function(object, ...)
{
	cobject <- list()
	cobject$response <- getResponseFormula(object$data)[[2]]
	cobject$ResponseClass <- attr(object, "class")[1]
	cobject$coefficients <- matrix(
		sapply(object$coef, function(x) x$beta),
		nrow = length(object$coef), byrow = TRUE)
	dimnames(cobject$coefficients) <- 
		list(1:nrow(cobject$coefficients), 
		names(object$coef[[1]]$beta))
	cobject$sigma2 <- sapply(object$coef, function(x) x$sigma2)
	names(cobject$sigma2) <- 1:length(cobject$sigma2)
	cobject
}


###################################################
### chunk number 10: 
###################################################

"BIC.mmlcr" <- 
function (object, ...) 
{
    if ((rt <- nargs()) > 1) {
        object <- list(object, ...)
        val <- as.data.frame(t(sapply(object, function(el) c(el$df,
        	 BIC(el)))))
        names(val) <- c("df", "BIC")
        row.names(val) <- as.character(match.call()[-1])
        val
    }
    else {
        object$BIC
    }
}

"fitted.mmlcr" <- 
function(object, which = 1:length(object$components), ...)
{
	fit.object <- vector("list", length(which))
	for(i in which)
		fit.object[[i]] <- object$components[[i]]$fitted
	fit.object
}

"formula.mmlcr" <- 
function(x, ...)
x$outer


"logLik.mmlcr" <-
function (object, ...) 
{
    if (length(list(...))) 
        warning("extra arguments discarded")
    val <- object$loglikelihood
    attr(val, "df") <- object$df
    class(val) <- "logLik"
    val
}



"print.mmlcr" <- 
function(x, ...)
{
	if(!is.null(cl <- x$call)) {
		cat("Call:\n")
		dput(cl)
	}
	cat("\nCoefficients:\n")
	print(x$gamma.matrix)
	cat("\nClass Percentages:\n")
	group.pct <- numeric(dim(x$post.prob)[2])
	names(group.pct) <- as.character(1:length(group.pct))
	for(k in 1:length(group.pct)) {
		temp.divisor <- 0
		temp.divisor <- temp.divisor + sum(x$post.prob[, k])
		group.pct[k] <- temp.divisor/dim(x$post.prob)[1]
	}
	print(round(100 * group.pct, 1))
	cat("\nAIC:", format(AIC(x), nsmall = 2), "\n")
	cat("BIC:", format(x$BIC, nsmall = 2), "\n")
	cat("loglikelihood:", format(x$loglikelihood, nsmall = 2),
		"\n")
	if(!is.null(attr(x, "warn")))
		cat("Warning: ", attr(x, "warn"), "\n")
	invisible(x)
}

"print.summary.mmlcr" <- 
function(x, digits = x$digits, ...)
{
	if(!is.null(cl <- x$call)) {
		cat("Call:\n")
		dput(cl)
	}
	cat("\nCoefficients:\n")
	print(x$coefficients, digits = digits)
	cat("\nStd. Errors:\n")
	print(x$standard.errors, digits = digits)
	if(!is.null(x$Wald.ratios)) {
		cat("\nValue/SE (Wald statistics):\n")
		print(x$coefficients/x$standard.errors, digits = digits
			)
	}
	cat("\nResidual Deviance:", format(x$deviance, nsmall = 2),
		"\n")
	cat("AIC:", format(x$AIC, nsmall = 2), "\n")
	cat("BIC:", format(x$BIC, nsmall = 2), "\n")

	if(!is.null(correl <- x$correlation)) {
		p <- dim(correl)[2]
		if(p > 1) {
			cat("\nCorrelation of Coefficients:\n")
			ll <- lower.tri(correl)
			correl[ll] <- format(round(correl[ll], digits))
			correl[!ll] <- ""
			print(correl[-1,  - p], quote = FALSE, ...)
		}
	}
	for(j in 1:length(x$components.summary)) {
		cat("\nResponse Component ", j, ":\n", sep = "")
		for(k in 1:length(x$components.summary[[j]])) {
			cat("\t", names(x$components.summary[[j]])[
				k], ":\n", sep = "")
			print(x$components.summary[[j]][[k]], digits = 
				digits)
		}
	}
	invisible(x)
}

"residuals.mmlcr" <- 
function(object, which = 1:length(object$components), ...)
{
	resid.object <- vector("list", length(which))
	for(i in which)
		resid.object[[i]] <- object$components[[i]]$residuals
	resid.object
}

"summary.mmlcr" <- 
function(object, correlation = TRUE, digits = options()$digits, 
	Wald.ratios = FALSE, vc = vcov(object), ...)
{
	r <- dim(object$gamma.matrix)[2]
	if (ncol(object$post.prob) > 1){
		se <- sqrt(diag(vc))
		coef <- object$gamma.matrix[-1,  , drop = FALSE]
		stderr <- matrix(se, nrow = dim(coef)[1], byrow = TRUE)
		dimnames(stderr) <- dimnames(coef)
		object$digits <- digits
		object$coefficients <- coef
		object$standard.errors <- stderr
		if(Wald.ratios)
			object$Wald.ratios <- coef/stderr
		if(correlation)
			object$correlation <- vc/outer(se, se)
	}
	object$AIC <- AIC(object)
	class(object) <- "summary.mmlcr"
	object$components.summary <- list()
	for(i in 1:length(object$components))
		object$components.summary[[i]] <- summary(object$
			components[[i]])
	object
}

"vcov.mmlcr" <- 
function(object, info = FALSE, ...)
{
	n.groups <- dim(object$gamma.matrix)[1]
	n.X <- dim(object$gamma.matrix)[2]
	post <- object$post.prob
	prior <- object$prior.prob
	if(all.equal(getCovariateFormula(object$outer), formula( ~ 1)) ==
		TRUE)
		mm <- matrix(rep(1, nrow(object$outer.df)), ncol = 1,
			dimnames = list(NULL, "(Intercept)"))
	else mm <- model.matrix(getCovariateFormula(object$outer),
			object$outer.df)
	cov <- matrix(0, ncol = length(object$gamma.matrix[-1,  ]),
		nrow = length(object$gamma.matrix[-1,  ]))
	cf <- t(outer(2:n.groups, dimnames(object$gamma.matrix)[[2]],
		function(x, y)
	paste(x, y, sep = ":")))
	for(k in 2:n.groups) {
		pk <- post[, k] - prior[, k]
		for(kk in 2:n.groups) {
			pkk <- post[, kk] - prior[, kk]
			cov[(k - 2) * n.X + 1:n.X, (kk - 2) * n.X +
				1:n.X] <- crossprod(mm * pk, mm * pkk)
		}
	}
	if(!info)
	    require(MASS)
		cov <- ginv(cov)
	structure(cov, dimnames = list(cf, cf))
}

"anova.mmlcr" <- 
function(object, ..., test = c("Chisq", "BIC", "AIC", "none"), 
	model.names = FALSE)
{
	test <- match.arg(test)
	dots <- list(...)
	if(length(dots) == 0)
		stop("anova is not implemented for a single mmlcr object"
			)
	mlist <- list(object, ...)
	nt <- length(mlist)
	temp.call <- match.call()
	temp.call$model.names <- temp.call$test <- NULL
	model.names.vec <- rep(as.character(temp.call[[2]]), nt)
	for(j in 2:nt)
		model.names.vec[j] <- as.character(temp.call[[j + 1]])
	dflis <- sapply(mlist, function(x) x$resid.df)
	s <- sort.list( - dflis)
	dflis <- dflis
	mlist <- mlist[s]
	if(any(!sapply(mlist, inherits, "mmlcr")))
		stop("not all objects are of class `mmlcr'")
	rsp <- unique(sapply(mlist, function(x)
	paste(formula(x)[2])))
	mds <- sapply(mlist, function(x)
	paste(formula(x)[3]))
	if(model.names)
		mds <- model.names.vec
	dfs <- dflis[s]
	lls <- sapply(mlist, function(x)
	x$loglikelihood)
	ll2s <- sapply(mlist, function(x)
	-2 * x$loglikelihood)
	n.groups <- sapply(mlist, function(x)
	dim(x$post.prob)[2])
	tss <- c("", paste(1:(nt - 1), 2:nt, sep = " vs "))
	aics <- sapply(mlist, AIC)
	bics <- sapply(mlist, BIC)
	df <- c(NA,  - diff(dfs))
	x2 <- c(NA,  - diff(ll2s))
	pr <- c(NA, 1 - pchisq(x2[-1], df[-1]))
	table <- data.frame(Model = mds, N.Groups = n.groups, DFs = dfs,
		LogLikelihood = lls, Test = tss, Df = df, LRtest = x2,
		Prob = pr, AIC = aics, BIC = bics)
	names(table) <-  c("Model", "No. of Classes", "DFs", 
		"Log Likelihood", "Test", "   Df", "LR stat.", 
		"Pr(Chi)", "AIC", "BIC")
	if(test == "none")
		table <- table[, 1:7]
	if(test == "AIC")
		table <- table[, c(1:4, 9)]
	if(test == "BIC")
		table <- table[, c(1:4, 10)]
	if(test == "Chisq")
		table <- table[, 1:8]
	structure(table, heading = c("Mixed-Mode Latent Class Regression Models\n", 
		paste("Response:", rsp)), class = c("Anova", "data.frame"))

}

"mmlcrclassify" <- 
function(object)
{
	temp.fn <- function(x)
	{
		x.max <- max(x)
		which(x == x.max)
	}
	temp.class <- apply(object$post.prob, 1, temp.fn)
	data.frame(class = temp.class, row.names = row.names(object$
		post.prob))
}

"mmlcrfit.nofit" <- 
function(object, weights, classnumber)
{
	object
}



###################################################
### chunk number 11: 
###################################################
"mmlcrcomponentinit.cnormlong" <- 
function(object, n.groups, prob, data, grouping)
{
	m <- match.call(expand = FALSE)
	m$object <- m$n.groups <- m$prob <- m$... <- NULL
	m$drop.unused.levels <- TRUE
	m$grouping <- NULL
	longformula <- formula(paste(object$formula[[2]], paste(
		unclass(object$formula), grouping, sep = "|")[[3]],
		sep = "~"))
	m$formula <- asOneFormula(longformula)
	m$na.action <- na.exclude
	m[[1]] <- as.name("model.frame")
	m <- eval.parent(m)
	attr(m, "terms") <- terms(longformula)
	attr(m, "formula") <- NULL
	.form <- longformula
	
	.temp.response <- data.frame(getResponse(m, .form))
	.status <- rep(1, length(.temp.response[, 1]))
	if (is.null(object$min)) object$min <- .Machine$double.xmin
	if (is.null(object$max)) object$max <- .Machine$double.xmax
	.status[.temp.response <= object$min] <- 2
	.status[.temp.response >= object$max] <- 0
	.data2 <- data.frame(m, surv.y = Surv(getResponse(m, .form),
		.temp.response, .status, type = "interval"))
	attr(.data2, "terms") <- terms(.form)
	
	
	x <- list(data = .data2, coef = list(), formula = .form)
	x$shortform <- formula(paste("surv.y", paste(splitFormula(
		.form, sep = "|")[[1]], collapse = " ")))
	survreg.result <- survreg(x$shortform, .data2, dist = 
		"gaussian", na.action = na.exclude)
	for(j in 1:n.groups) {
		x$coef[[j]] <- list(beta = coef(survreg.result),
			sigma2 = survreg.result$scale^2)
	}
	x$n <- sum(!is.na(getResponse(.data2, .form)))
	x$m <- (length(x$coef[[1]]$beta) + 1) * n.groups
	x$min <- object$min
	x$max <- object$max
	attr(x, "class") <- c(attr(object, "class"), "mmlcrlong")
	x
}

"mmlcrfit.cnormlong" <- 
function(object, weights, classnumber)
{
	data <- object$data
	form.full <- formula(attr(data, "terms"))
	form <- object$shortform
	data$.wts <- .wts <- weights[match(getGroups(data, form.full, 1),
		row.names(weights)), classnumber]
	scale.min <- max(sqrt(max(sapply(object$coef, function(x)
	x$sigma2)[ - classnumber]))/10, .Machine$double.xmin, na.rm = TRUE
		)
	scale.max <- min(sqrt(min(sapply(object$coef, function(x)
	x$sigma2)[ - classnumber])) * 10, .Machine$double.xmax, na.rm
		 = TRUE)
	survreg.result <- survreg(form, data, dist = "gaussian", 
		weights = .wts, subset = .wts > 0, na.action = 
		na.exclude)
	if((survreg.result$scale < scale.min) | (survreg.result$scale >
		scale.max)) {
		.minimum <- object$min
		.maximum <- object$max
		model.mat <- model.matrix.lm(form, data)
		.y.mean <- model.mat %*% as.matrix(object$coef[[
			classnumber]]$beta)
		.y <- getResponse(data)
		minimizeme.fn <- function(sigma)
		{
 - sum(.wts * log(pmax(ifelse(is.na(.y | is.na(.y.mean)), 1, ifelse(
				.y <= .minimum, pnorm(.minimum, mean = 
				.y.mean, sd = sigma), ifelse(.y >= 
				.maximum, 1 - pnorm(.maximum, mean = 
				.y.mean, sd = sigma), dnorm(.y, mean = 
				.y.mean, sd = sigma)))), .Machine$
				double.xmin)))
		}
		scale <- optim(sqrt(object$coef[[classnumber]]$sigma2),
			minimizeme.fn, lower = scale.min, upper = 
			scale.max, method = "L-BFGS-B")$par
		survreg.result <- survreg(form, data, dist = "gaussian",
			weights = .wts, subset = .wts > 0, na.action = 
			na.exclude, scale = scale)
		survreg.result$scale <- scale
	}
	object$coef[[classnumber]]$beta <- coef(survreg.result)
	object$coef[[classnumber]]$sigma2 <- survreg.result$scale^
		2
	object
}

"mmlcrlike.cnormlong" <- 
function(object, weights, classnumber, final = FALSE)
{
	ids <- row.names(weights)
	form.full <- formula(attr(object$data, "terms"))
	groups <- getGroups(object$data, form.full)
	include.vec <- is.element(groups, ids)
	data <- object$data[include.vec,  ]
	form <- object$shortform
	coef <- object$coef[[classnumber]]
	beta <- coef$beta
	model.mat <- model.matrix.lm(form, data)
	mu <- model.mat %*% as.matrix(beta)
	if(!final) {
		temp.like <- 
		dcnorm(getResponse(data, form)[,
			1], mean = mu, sigma = rep(sqrt(coef$sigma2),
			length(mu)), min = rep(object$min, length(
			mu)), max = rep(object$max, length(mu)))
		temp.like[is.na(temp.like)] <- 1
		prod.like <- tapply(temp.like, groups, prod)
		prod.like[match(ids, names(prod.like))]
	}
	else {
		.wts <- weights[match(groups, row.names(weights)), classnumber]
		mu <- pmin(mu, object$max)
		mu <- pmax(mu, object$min)
		if(classnumber == 1) {
			object$fitted <- rep(NA, nrow(object$data))
			object$fitted[include.vec] <- mu * .wts
		}
		else object$fitted[include.vec] <- object$fitted[
				include.vec] + mu * .wts
		object$residuals <- rep(NA, nrow(object$data))
		object$residuals[include.vec] <- as.numeric(getResponse(
			data)) - object$fitted[include.vec]
		object
	}
}

"summary.cnormlong" <- 
function(object, ...)
{
	cobject <- list()
	cobject$response <- getResponseFormula(object$data)[[2]]
	cobject$ResponseClass <- attr(object, "class")[1]
	cobject$minimum <- object$min
	cobject$maximum <- object$max
	cobject$coefficients <- matrix(sapply(object$coef, function(x)
	x$beta), nrow = length(object$coef), byrow = TRUE)
	dimnames(cobject$coefficients) <- list(1:nrow(cobject$
		coefficients), names(object$coef[[1]]$beta))
	cobject$sigma2 <- sapply(object$coef, function(
		x)
	x$sigma2)
	names(cobject$sigma2) <- 1:length(cobject$sigma2)
	cobject
}

"dcnorm" <- 
function(x, mean = rep(0, length(x)), sigma = rep(1, length(x)), min = 
	rep(-(.Machine$double.xmax), length(x)), max = rep(.Machine$double.xmax, length(
	x)))
{
	y <- dnorm(x, mean = mean, sd = sigma)
	if(any(x <= min))
		y[x <= min] <- pnorm(min[x <= min], mean = mean[x <=
			min], sd = sigma[x <= min])
	if(any(x >= max))
		y[x >= max] <- 1 - pnorm(max[x >= max], mean = mean[
			x >= max], sd = sigma[x >= max])
	y[is.na(y)] <- 1
	y
}



###################################################
### chunk number 12: 
###################################################

"mmlcrcomponentinit.normonce" <- 
function(object, n.groups, prob, data, grouping)
{
	m <- match.call(expand = FALSE)
	m$object <- m$n.groups <- m$prob <- m$... <- NULL
	m$drop.unused.levels <- TRUE
	m$grouping <- NULL
	longformula <- formula(paste(object$formula[[2]], paste(
		unclass(object$formula), grouping, sep = "|")[[3]],
		sep = "~"))
	m$formula <- asOneFormula(longformula)
	m$na.action <- na.exclude
	m[[1]] <- as.name("model.frame")
	m <- eval.parent(m)
	attr(m, "terms") <- terms(longformula)
	attr(m, "formula") <- NULL
	.form <- longformula
	# drop all but the first observation for each individual
	m <- m[match(row.names(prob), getGroups(m, .form)),]
		x <- list(data = m, coef = list(), formula = .form)
	x$shortform <- object$formula
	lm.result <- lm(formula(x$shortform), m, na.action = na.exclude)
	for(j in 1:n.groups) {
		x$coef[[j]] <- list(beta = coef(lm.result),
			sigma2 = sum(lm.result$residuals^2) / lm.result$df.resid)
	}
	x$n <- sum(!is.na(getResponse(m, .form)))
	x$m <- (length(x$coef[[1]]$beta) + 1) * n.groups
	attr(x, "class") <- c(attr(object, "class"), "normlong", "mmlcronce")
	x
}

"mmlcrcomponentinit.cnormonce" <- 
function(object, n.groups, prob, data, grouping)
{
	m <- match.call(expand = FALSE)
	m$object <- m$n.groups <- m$prob <- m$... <- NULL
	m$drop.unused.levels <- TRUE
	m$grouping <- NULL
	longformula <- formula(paste(object$formula[[2]], paste(
		unclass(object$formula), grouping, sep = "|")[[3]],
		sep = "~"))
	m$formula <- asOneFormula(longformula)
	m$na.action <- na.exclude
	m[[1]] <- as.name("model.frame")
	m <- eval.parent(m)
	attr(m, "terms") <- terms(longformula)
	attr(m, "formula") <- NULL
	.form <- longformula
	# drop all but the first observation for each individual
	m <- m[match(row.names(prob), getGroups(m, .form)),]
		x <- list(data = m, coef = list(), formula = .form)
	.temp.response <- data.frame(getResponse(m, .form))
	.status <- rep(1, length(.temp.response[, 1]))
	if (is.null(object$min)) object$min <- .Machine$double.xmin
	if (is.null(object$max)) object$max <- .Machine$double.xmax
	.status[.temp.response <= object$min] <- 2
	.status[.temp.response >= object$max] <- 0
	.data2 <- data.frame(m, surv.y = Surv(getResponse(m, .form),
		.temp.response, .status, type = "interval"))
	attr(.data2, "terms") <- terms(.form)
	
	
	x <- list(data = .data2, coef = list(), formula = .form)
	x$shortform <- formula(paste("surv.y", paste(splitFormula(
		.form, sep = "|")[[1]], collapse = " ")))
	survreg.result <- survreg(x$shortform, .data2, dist = 
		"gaussian", na.action = na.exclude)
	for(j in 1:n.groups) {
		x$coef[[j]] <- list(beta = coef(survreg.result),
			sigma2 = survreg.result$scale^2)
	}
	x$n <- sum(!is.na(getResponse(.data2, .form)))
	x$m <- (length(x$coef[[1]]$beta) + 1) * n.groups
	x$min <- object$min
	x$max <- object$max
	attr(x, "class") <- c(attr(object, "class"), "cnormlong", "mmlcronce")
	x
}


###################################################
### chunk number 13: 
###################################################

"mmlcrcomponentinit.multinomlong" <- 
function(object, n.groups, prob, data, grouping)
{
	m <- match.call(expand = FALSE)
	m$object <- m$n.groups <- m$prob <- m$... <- NULL
	m$drop.unused.levels <- TRUE
	m$grouping <- NULL
	longformula <- formula(paste(object$formula[[2]], paste(
		unclass(object$formula), grouping, sep = "|")[[3]],
		sep = "~"))
	m$formula <- asOneFormula(longformula)
	m$na.action <- na.exclude
	m[[1]] <- as.name("model.frame")
	m <- eval.parent(m)
	attr(m, "terms") <- terms(longformula)
	attr(m, "formula") <- NULL
	.form <- longformula
	x <- list(data = m, coef = list(), formula = .form)
	x$shortform <- object$formula
	x$coef <- mmlcrfit.multinomlong(x, prob, 1)$coef
	x$n <- sum(!is.na(getResponse(m)))
	x$m <- length(x$coef[[1]]) * n.groups
	attr(x, "class") <- c(attr(object, "class"), "mmlcrlong")
	x
}

"mmlcrcomponentinit.multinomonce" <- 
function(object, n.groups, prob, data, grouping)
{
	m <- match.call(expand = FALSE)
	m$object <- m$n.groups <- m$prob <- m$... <- NULL
	m$drop.unused.levels <- TRUE
	m$grouping <- NULL
	longformula <- formula(paste(object$formula[[2]], paste(
		unclass(object$formula), grouping, sep = "|")[[3]],
		sep = "~"))
	m$formula <- asOneFormula(longformula)
	m$na.action <- na.exclude
	m[[1]] <- as.name("model.frame")
	m <- eval.parent(m)
	attr(m, "terms") <- terms(longformula)
	attr(m, "formula") <- NULL
	.form <- longformula
	# drop all but the first observation for each individual
	m <- m[match(row.names(prob), getGroups(m, .form)),]
		x <- list(data = m, coef = list(), formula = .form)
	x <- list(data = m, coef = list(), formula = .form)
	x$shortform <- object$formula
	x$coef <- mmlcrfit.multinomlong(x, prob, 1)$coef
	x$n <- sum(!is.na(getResponse(m)))
	x$m <- length(x$coef[[1]]) * n.groups
	attr(x, "class") <- c(attr(object, "class"), "multinomlong", "mmlcronce")
	x
}

"mmlcrfit.multinomlong" <- 
function(object, weights, classnumber)
{
	data <- object$data
	form.full <- formula(attr(data, "terms"))
	form <- object$shortform
	data$.wts <- weights[match(getGroups(data, form.full, 1),
		row.names(weights)), classnumber]
	component.multinom.result <- multinom(form, data = data, weights = .wts,
		na.action = na.exclude, subset = (.wts > 0), maxit = 
		150, trace = FALSE)
	object$coef[[classnumber]] <- 
		pmin(pmax(coef(component.multinom.result), -30), 30)
	object
}

"mmlcrlike.multinomlong" <- 
function(object, weights, classnumber, final = FALSE)
{
	ids <- row.names(weights)
	form.full <- formula(attr(object$data, "terms"))
	groups <- getGroups(object$data, form.full)
	include.vec <- is.element(groups, ids)
	data <- object$data[include.vec,  ]
	form <- object$shortform
	coef <- rbind(0, object$coef[[classnumber]])
	model.mat <- model.matrix.lm(form,data)
	like <- data.frame(exp(model.mat %*% t(coef)))
	like <- like/apply(like, 1, sum)
	like <- apply(like, 2, pmax, .Machine$double.xmin)
	if(!final) {
		temp.like <- like[cbind(1:dim(like)[1], as.numeric(
			getResponse(data)))]
		temp.like[is.na(temp.like)] <- 1
		prod.like <- tapply(temp.like, groups, prod)
		prod.like[match(ids, names(prod.like))]
	}
	else {
		.wts <- weights[match(groups, row.names(weights)), classnumber]
		all.responses <- matrix(1:ncol(like), ncol = ncol(
			like), nrow = nrow(like), byrow = TRUE)
		mu <- apply(like * all.responses, 1, sum)
		if(classnumber == 1) {
			object$fitted <- rep(NA, nrow(object$data))
			object$fitted[include.vec] <- mu * .wts
		}
		else object$fitted[include.vec] <- object$fitted[
				include.vec] + mu * .wts
		object$residuals <- rep(NA, nrow(object$data))
		object$residuals[include.vec] <- as.numeric(getResponse(
			data)) - object$fitted[include.vec]
		object
	}
}

"summary.multinomlong" <- 
function(object, ...)
{
	cobject <- list()
	cobject$response <- getResponseFormula(object$data)[[2]]
	cobject$ResponseClass <- attr(object, "class")[1]
	cobject$coefficients <- object$coef
	names(cobject$coefficients) <- paste("LatentClass", 1:length(
		cobject$coefficients), sep = "")
	cobject
}


###################################################
### chunk number 14: 
###################################################

"mmlcrcomponentinit.poislong" <- 
function(object, n.groups, prob, data, grouping)
{
	m <- match.call(expand = FALSE)
	m$object <- m$n.groups <- m$prob <- m$... <- NULL
	m$drop.unused.levels <- TRUE
	m$grouping <- NULL
	longformula <- formula(paste(object$formula[[2]], paste(
		unclass(object$formula), grouping, sep = "|")[[3]],
		sep = "~"))
	m$formula <- asOneFormula(longformula)
	m$na.action <- na.exclude
	m[[1]] <- as.name("model.frame")
	m <- eval.parent(m)
	attr(m, "terms") <- terms(longformula)
	attr(m, "formula") <- NULL
	.form <- longformula
	x <- list(data = m, coef = list(), formula = .form)
	x$shortform <- object$formula
	x$coef <- mmlcrfit.poislong(x, prob, 1)$coef
	x$n <- sum(!is.na(getResponse(m)))
	x$m <- length(x$coef[[1]]) * n.groups
	attr(x, "class") <- c(attr(object, "class"), "mmlcrlong")
	x
}

"mmlcrcomponentinit.poisonce" <- 
function(object, n.groups, prob, data, grouping)
{
	m <- match.call(expand = FALSE)
	m$object <- m$n.groups <- m$prob <- m$... <- NULL
	m$drop.unused.levels <- TRUE
	m$grouping <- NULL
	longformula <- formula(paste(object$formula[[2]], paste(
		unclass(object$formula), grouping, sep = "|")[[3]],
		sep = "~"))
	m$formula <- asOneFormula(longformula)
	m$na.action <- na.exclude
	m[[1]] <- as.name("model.frame")
	m <- eval.parent(m)
	attr(m, "terms") <- terms(longformula)
	attr(m, "formula") <- NULL
	# drop all but the first observation for each individual
	m <- m[match(row.names(prob), getGroups(m, .form)),]
	.form <- longformula
	x <- list(data = m, coef = list(), formula = .form)
	x$shortform <- object$formula
	x$coef <- mmlcrfit.poislong(x, prob, 1)$coef
	x$n <- sum(!is.na(getResponse(m)))
	x$m <- length(x$coef[[1]]) * n.groups
	attr(x, "class") <- c(attr(object, "class"), "poislong", "mmlcronce")
	x
}

"mmlcrfit.poislong" <- 
function(object, weights, classnumber)
{
	data <- object$data
	form.full <- formula(attr(data, "terms"))
	form <- object$shortform
	data$.wts <- weights[match(getGroups(data, form.full, 1), row.names(
		weights)), classnumber]
	glm.result <- glm(form, data = data, weights = .wts, quasi(
		link = "log", var = "mu"), na.action = na.exclude,
		subset = (.wts > 0))
	object$coef[[classnumber]] <- coef(glm.result)
	object$over[[classnumber]] <- sum(data$.wts * (resid(glm.result, "pearson")^2))/
			(sum(data$.wts) - glm.result$rank)# Was:
#	object$over[[classnumber]] <- sum(.wts * (resid(glm.result,
#		"response")^2/predict(glm.result, type = "response")))/
#		(sum(.wts) - glm.result$rank)
	object
}

"mmlcrlike.poislong" <- 
function(object, weights, classnumber, final = FALSE)
{
	ids <- row.names(weights)
	form.full <- formula(attr(object$data, "terms"))
	groups <- getGroups(object$data, form.full)
	include.vec <- is.element(groups, ids)
	data <- object$data[include.vec,  ]
	form <- object$shortform
	coef <- object$coef[[classnumber]]
	model.mat <- model.matrix.lm(form, data)
	log.lambda <- model.mat %*% as.matrix(coef)
	if(!final) {
		temp.like <- exp(getResponse(data) * log.lambda - exp(
			log.lambda) -lgamma(getResponse(data) + 1))
		temp.like[is.na(temp.like)] <- 1
		prod.like <- tapply(temp.like, groups, prod)
		prod.like[match(ids, names(prod.like))]
	}
	else {
		.wts <- weights[match(groups, row.names(weights)), classnumber]
		if(classnumber == 1) {
			object$fitted <- rep(NA, nrow(object$data))
			object$fitted[include.vec] <- exp(log.lambda) *
				.wts
		}
		else object$fitted[include.vec] <- object$fitted[
				include.vec] + exp(log.lambda) * .wts
		object$residuals <- rep(NA, nrow(object$data))
		object$residuals[include.vec] <- getResponse(data) -
			object$fitted[include.vec]
		object
	}
}

"summary.poislong" <- 
function(object, ...)
{
	cobject <- list()
	cobject$response <- getResponseFormula(object$data)[[2]]
	cobject$ResponseClass <- attr(object, "class")[1]
	cobject$coefficients <- matrix(unlist(object$coef), nrow = 
		length(object$coef), byrow = TRUE)
	dimnames(cobject$coefficients) <- list(1:nrow(cobject$
		coefficients), names(object$coef[[1]]))
	cobject$overdispersion <- object$over
	names(cobject$overdispersion) <- 1:length(object$over)
	cobject
}


###################################################
### chunk number 15: 
###################################################
#"dnb1" <- 
#function(y, mu, alpha)
#{
#	newalpha <- mu/alpha
#	exp(lgamma(y + newalpha) - (lgamma(y + 1) + lgamma(newalpha)) -
#		newalpha * log(1 + alpha) + y * log(alpha/(1 + alpha)))
#}

"dnb2" <- 
function(y, mu, th)
{
	exp(lgamma(th + y) - lgamma(th) + th * log(th/(th + mu)) + y *
		log((mu + (y == 0))/(th + mu)) - lgamma(y + 1))
}

"mmlcrcomponentinit.nblong" <- 
function(object, n.groups, prob, data, grouping)
{
	m <- match.call(expand = FALSE)
	m$object <- m$n.groups <- m$prob <- m$... <- NULL
	m$drop.unused.levels <- TRUE
	m$grouping <- NULL
	longformula <- formula(paste(object$formula[[2]], paste(
		unclass(object$formula), grouping, sep = "|")[[3]],
		sep = "~"))
	m$formula <- asOneFormula(longformula)
	m$na.action <- na.exclude
	m[[1]] <- as.name("model.frame")
	m <- eval.parent(m)
	attr(m, "terms") <- terms(longformula)
	attr(m, "formula") <- NULL
	.form <- longformula
	theta <- vector("list", n.groups)
	for(k in 1:n.groups)
		theta[[k]] <- 4
	x <- list(data = m, coef = vector("list", n.groups), theta
		 = theta, formula = .form)
	x$shortform <- object$formula
	x$coef[[1]] <- mmlcrfit.nblong(x, prob, 1)$coef[[1]]
	for(k in 2:n.groups)
		x$coef[[k]] <- x$coef[[1]]
	x$n <- sum(!is.na(getResponse(m)))
	x$m <- (length(x$coef[[1]]) + 1) * n.groups
	attr(x, "class") <- c(attr(object, "class"), "mmlcrlong")
	x
}

"mmlcrcomponentinit.nbonce" <- 
function(object, n.groups, prob, data, grouping)
{
	m <- match.call(expand = FALSE)
	m$object <- m$n.groups <- m$prob <- m$... <- NULL
	m$drop.unused.levels <- TRUE
	m$grouping <- NULL
	longformula <- formula(paste(object$formula[[2]], paste(
		unclass(object$formula), grouping, sep = "|")[[3]],
		sep = "~"))
	m$formula <- asOneFormula(longformula)
	m$na.action <- na.exclude
	m[[1]] <- as.name("model.frame")
	m <- eval.parent(m)
	attr(m, "terms") <- terms(longformula)
	attr(m, "formula") <- NULL
	.form <- longformula
	# drop all but the first observation for each individual
	m <- m[match(row.names(prob), getGroups(m, .form)),]
		x <- list(data = m, coef = list(), formula = .form)
	theta <- vector("list", n.groups)
	for(k in 1:n.groups)
		theta[[k]] <- 4
	x <- list(data = m, coef = vector("list", n.groups), theta
		 = theta, formula = .form)
	x$shortform <- object$formula
	x$coef[[1]] <- mmlcrfit.nblong(x, prob, 1)$coef[[1]]
	for(k in 2:n.groups)
		x$coef[[k]] <- x$coef[[1]]
	x$n <- sum(!is.na(getResponse(m)))
	x$m <- (length(x$coef[[1]]) + 1) * n.groups
	attr(x, "class") <- c(attr(object, "class"), "nblong", "mmlcronce")
	x
}

"mmlcrfit.nblong" <- 
function(object, weights, classnumber)
{
	data <- object$data
	form.full <- formula(attr(data, "terms"))
	form <- object$shortform
	data$.wts <- weights[match(getGroups(data, form.full, 1), row.names(
		weights)), classnumber]
	glm.result <- glm(form, data = data, weights = .wts, na.action
		= na.exclude, subset = (.wts > 1e-012), 
		neg.bin(object$theta[[classnumber]]))
	object$coef[[classnumber]] <- coef(glm.result)
	if(coef(glm.result)[1] > 1000 || coef(glm.result)[1] < -1000)
		glm.result <- glm(form, data = data, weights = .wts,
			na.action = na.exclude, subset = (.wts > 0),
			family = poisson)
	.log.lambda <- model.matrix.lm(form, data) %*% as.matrix(
		object$coef[[classnumber]])
	object$theta[[classnumber]] <- as.vector(theta.mmmod(
		getResponse(data), exp(.log.lambda), sum(data$.wts) - 
		length(coef(glm.result)), data$.wts))
	object
}



"mmlcrlike.nblong" <- 
function(object, weights, classnumber, final = FALSE)
{
	ids <- row.names(weights)
	form.full <- formula(attr(object$data, "terms"))
	groups <- getGroups(object$data, form.full)
	include.vec <- is.element(groups, ids)
	data <- object$data[include.vec,  ]
	form <- object$shortform
	theta <- object$theta[[classnumber]]
	coef <- object$coef[[classnumber]]
	model.mat <- model.matrix.lm(form, data)
	log.lambda <- model.mat %*% as.matrix(coef)
	if(!final) {
		temp.like <- dnb2(getResponse(data), mu = exp(
			log.lambda), th = theta)
		temp.like[is.na(temp.like)] <- 1
		prod.like <- tapply(temp.like, groups, prod)
		prod.like[match(ids, names(prod.like))]
	}
	else {
		.wts <- weights[match(groups, row.names(weights)), classnumber]
		if(classnumber == 1) {
			object$fitted <- rep(NA, nrow(object$data))
			object$fitted[include.vec] <- exp(log.lambda) *
				.wts
		}
		else object$fitted[include.vec] <- object$fitted[
				include.vec] + exp(log.lambda) * .wts
		object$residuals <- rep(NA, nrow(object$data))
		object$residuals[include.vec] <- getResponse(data) -
			object$fitted[include.vec]
		object
	}
}


"summary.nblong" <- 
function(object, ...)
{
	cobject <- list()
	cobject$response <- getResponseFormula(object$data)[[2]]
	cobject$ResponseClass <- attr(object, "class")[1]
	cobject$coefficients <- matrix(unlist(object$coef), nrow = 
		length(object$coef), byrow = TRUE)
	dimnames(cobject$coefficients) <- list(1:nrow(cobject$
		coefficients), names(object$coef[[1]]))
	cobject$theta <- unlist(object$theta)
	names(cobject$theta) <- 1:length(cobject$theta)
	cobject
}

"theta.mmmod" <- 
function(y, u, dfr, wts, limit = 10, eps =  .Machine$double.eps^0.25)
{
	n <- sum(wts)
	t0 <- n/sum(wts * (y/u - 1)^2)
	it <- 0
	del <- 1
	while((it <- it + 1) < limit && abs(del) > eps) {
		t0 <- abs(t0)
		del <- (sum((wts * (y - u)^2)/(u + u^2/t0)) - dfr)/
			sum((wts * (y - u)^2)/(u + t0)^2)
		t0 <- t0 - del
	}
	if(t0 < 0) {
		t0 <- 1e-006
		warning("estimator truncated at zero")
		attr(t0, "warn") <- "estimate truncated at zero"
	}
	if(is.infinite(t0) | t0 > 100000000)
		t0 <- 100000000
	t0
}


###################################################
### chunk number 16: 
###################################################

"plot.mmlcr" <- 
function(x, which = 1:length(x$components), ...)
{
	for(k in which)
		plot(x$components[[k]], x$post.prob, ...)
}

"plot.mmlcrlong" <- 
function(x, post.prob, smooth = 0, xlab = names(data)[2], ylab = 
	names(data)[1], ylim = c(min(yy, na.rm = TRUE), max(yy, na.rm = TRUE)
	), cols = rep(1, dim(post.prob)[2]), pch = as.character(1:
	(dim(post.prob)[2])), ...)
{
	data <- x$data
	x <- sort(unique(data[, 2]))
	yy <- matrix(ncol = length(x), nrow = (dim(post.prob)[2]))
	form.full <- formula(attr(data, "terms"))
	if(smooth) {
		require(modreg)
		for(k in 1:(dim(post.prob)[2])) {
			wts <- post.prob[match(getGroups(data, form.full, 1), 
				row.names(post.prob)), k]
			yy[k,  ] <- supsmu(data[, 2], data[, 1], wt = 
				wts, bass = smooth)$y
		}
	}
	else {
		for(k in 1:(dim(post.prob)[2])) {
			wts <- post.prob[match(getGroups(data, form.full, 1), 
				row.names(post.prob)), k]
				for(i in 1:length(x))
				yy[k, i] <- weighted.mean(data[data[
					, 2] == x[i], 1], wts[data[
					, 2] == x[i]], na.rm = TRUE)
		}
	}
	y <- yy[1,  ]
	plot(x, y, type = "b", pch = pch[1], xlab = xlab, ylab = ylab,
		ylim = ylim, col = cols[1], ...)
	if (dim(post.prob)[2] > 1){
		for(k in 2:(dim(post.prob)[2])) {
			lines(x, yy[k,  ], type = "b", pch = pch[k], col = cols[
				k])
		}
	}
}

"plot.mmlcronce" <- 
function(x, post.prob, ...)
{
	y <- getResponse(x$data, form = x$formula)
	y.min <- min(y)
	y.max <- max(y)
	plot(range(y), c(1 - 0.5, ncol(post.prob) + 0.5), type = "n",
		axes = FALSE, xlab = deparse(x$formula[[2]]), ylab = 
		"Class")
	box()
	axis(2, at = 1:ncol(post.prob))
	axis(1)
	for(i in 1:ncol(post.prob)) {
		abline(h = i, lty = 2)
		y.mean <- weighted.mean(y, post.prob[, i], na.rm = TRUE)
		y.sd <- sqrt(cov.wt(as.data.frame(y), post.prob[, i])$cov)
		lines(c(max(y.mean - y.sd, y.min), y.mean, min(y.mean +
			y.sd, y.max)), rep(i, 3), type = "b", pch = 18)
	}
	invisible()
}


 
 "plot.multinomonce" <- 
 function(x, post.prob, ...)
	warning("Sorry, no plot method for multinomonce.\n")
	
 "plot.multinomlong" <- 
 function(x, post.prob, ...)
	warning("Sorry, no plot method for multinomlong.\n")

.First.lib <- function(lib, pkg)
{
    require(nlme)
    require(nnet)
    require(survival)
}

 "postprob" <- 
 function(object)
 	object$post.prob



