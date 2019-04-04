perturb <- function(mod,pvars=NULL,prange=NULL,ptrans=NULL,postrun=NULL,pfac=NULL,uniform=FALSE,niter=100) {
	cutsp<-function(indx,tbl) {
		findInterval(runif(1),tbl[indx,],rightmost.closed=TRUE)
	}

	if (class(try(getCall(mod,silent=TRUE)))=="try-error") stop("Not a valid model call")

	stopifnot(is.list(pfac)||is.null(pfac))

	# Method used to extract fixed and random effects from an lme object
	nms<-GetNames(mod)
	stopifnot(all(pvars %in% nms))


	# nlme::lme requires a different method for extracting the data
	modData<-GetData(mod)
	modData<-as.data.frame(modData)

	result<-NULL
	ncases<-nrow(modData)
	ModCall<-deparse(getCall(mod),width.cutoff = 500)
	result$ModCall<-ModCall
	# Use the coefs method defined below
	# Defaults to coef but with special features for lme4 and nlme
	allb<-coefs(mod)

	# Use data frame "tmpData" for re-estimating the model
	ModCall<-gsub("data = (\\w+)","data = tmpData",ModCall)


	# Check the pvars specification
	if (length(pvars) > 0) {
		stopifnot(is.vector(modData[,pvars])||is.data.frame(modData[,pvars]))
		stopifnot(length(pvars)==length(prange))
		result$pvars<-pvars
		result$prange<-prange
		if (length(ptrans)>0) result$ptrans<-ptrans
	}

	if (length(pfac[[1]]) > 0) {
	  rcls.tbl<-NULL
	  if (is.list(pfac[[1]])) n<-length(pfac)
	  else n<-1
	  for (i in 1:n) {
	    if (n == 1) lstnm<-pfac
	    else lstnm<-pfac[[i]]
	    stopifnot(all(lstnm[[1]] %in% nms))
	    lstnm$data <- quote(modData)
	    rcls<-do.call("reclassify",lstnm)
	    rcls.tbl<-c(rcls.tbl,list(rcls))
	  }
	  result$reclassify.tables<-rcls.tbl
	}
  result$ModCall2<-ModCall

	if (uniform) {
		result$distribution<-"uniform"
	}
	else {
		result$distribution<-"normal"
	}

  if (!is.null(ptrans)){
    # Add "tmpData$" prefix to the variables in "ptrans"
    # First, get the variable names
    tvars<-all.vars(parse(text=ptrans))
    # Use "\<" and "\>" in the regex, beginning or end of word
    for (str in tvars) {
      ptrans<-gsub(paste0("\\<(",str,")\\>"),"tmpData$\\1",ptrans)
    }
  }

  # result$iter$perturbed is a list of data frames containing
  # modified variables at each iteration, for recreating the results at a specific iteration
  result$perturbed <- vector("list", niter+1)
  result$perturbed[[1]] <- as.data.frame(modData[,c(pvars,pfac[[1]])])
  result$messages[1] <- NULL

  if (!is.null(postrun)){
  	# Create a list of dataframes with the results from postrun
  	result$postrun <- vector("list",niter+1)
  	# The first dataframe will be from the postrun specification on the input model
  	result$postrun[[1]] <- as.data.frame(eval(parse(text=postrun)))
  	# Modify the postrun specification so it will use the "mod2" model object
  	result$postrun$call <- postrun
  	postrun <- sub(deparse(substitute(mod)),"mod2",postrun)
  }

	for (k in 1:niter) {
	  # Make a copy of the data
	  tmpData<-modData
		# add random perturbances to pvars using values in prange
		if (length(prange)>0) {
			for (i in 1:length(prange)) {
			  if (uniform) {
			    updtStr=paste0("tmpData$",pvars[i],"<-tmpData$",pvars[i],"+","runif(ncases,",-prange[i],",",prange[i],")")
			  }
			  else {
			    updtStr=paste0("tmpData$",pvars[i],"<-tmpData$",pvars[i],"+","rnorm(ncases,0,",prange[i],")")
			  }
			  eval(parse(text=updtStr))
			}
		}
		# re-execute the transformations
		for (trans in ptrans) eval(parse(text=trans))
		# reclassify factors here
		if (length(pfac[[1]])>0) {
			for (i in 1:length(rcls.tbl)) {
				tbl<-rcls.tbl[[i]]$cum.reclass.prob
				tbl<-cbind(0,tbl)
				tmpvar<-eval(parse(text=paste0("tmpData$",rcls.tbl[[i]]$variable)))
				eval(parse(text=paste0("tmpData$",pfac[[i]],"<-sapply(tmpvar,cutsp,tbl)")))
				eval(parse(text=paste0("tmpData$",pfac[[i]],"<-as.factor(tmpData$",pfac[[i]],")")))
			}
		}
		# re-estimate the model using the perturbed variables
		mod2<-eval(parse(text=ModCall))
		# collect the coefficients
		allb<-rbind(allb,coefs(mod2))
		# result$perturbed is a list of data frames containing modified variables
		result$perturbed[[k+1]] <- as.data.frame(tmpData[,c(pvars,pfac[[1]])])
		result$messages[k+1] <- GetWarnings(mod2)
		# Re-run the postrun specification, add to the results
		if (!is.null(postrun)){ result$postrun[[k+1]] <- as.data.frame(eval(parse(text=postrun))) }
	}
	# "allb" is the rowname value for the first row of allb; remove
	rownames(allb)<-NULL
	result$coeff.table<-allb
	class(result)<-"perturb"
	result
}

summary.perturb <-function(object,dec.places=3,full=FALSE,...) {
	coeffs<-object$coeff.table
	mysumm<-cbind(apply(coeffs,2,mean,na.rm=TRUE),apply(coeffs,2,sd,na.rm=TRUE),apply(coeffs,2,min,na.rm=TRUE),apply(coeffs,2,max,na.rm=TRUE))
	colnames(mysumm)<-c("mean","s.d.","min","max")
	object$coeff.table<-NULL
	object$summ<-mysumm
	object$dec.places<-dec.places
	object$full<-full
	dots<-substitute(expression(...))
	dots<-sub("^expression\\((.*)\\)$","\\1", deparse(dots))
	object$dots<-dots
	class(object)<-"summary.perturb"
	object
}

print.summary.perturb <-function(x,...) {
	if (x$full) {
		cat("formula:\n",x$formula,"\nformula2:\n",x$formula2,"\n\n")
	}
	if (length(x$pvars)>0) {
		cat("Perturb variables:\n")
		if (x$distribution=="uniform") {
			for (i in 1:length(x$pvars)) {
				prnt<-paste("uniform[",-round(x$prange[i],1),",",round(x$prange[i],1),"]",sep="")
				cat(x$pvars[i],"\t\t",prnt,"\n")
			}
		}
		else {
			for (i in 1:length(x$pvars)) {
				prnt<-paste("normal(0,",round(x$prange[i],1),")",sep="")
				cat(x$pvars[i],"\t\t",prnt,"\n")
			}
		}
		cat("\n")
	}
	if (length(x$ptrans)>0) {
		cat("Transformations:\n")
		for (trans in x$ptrans) {
			cat(trans,"\n")
		}
		if (x$full) {
			cat("\nTransformations2:\n")
			for (trans in x$ptrans2) {
				cat(trans,"\n")
			}
		}
		cat("\n")
	}
	if (!is.null(x$reclassify.tables)) {
		for (i in 1:length(x$reclassify.tables)) {
			if (x$dots=="") print(x$reclassify.tables[[i]],dec.places=x$dec.places,full=x$full,...)
			else eval(parse(text=paste("print(x$reclassify.tables[[i]],dec.places=x$dec.places,full=x$full,",x$dots,",...)")))
		}
	}
	cat("Impact of perturbations on coefficients:\n")
	#print(round(x$summ,x$dec.places),...)
	eval(parse(text=paste("print(round(x$summ,x$dec.places),",x$dots,",...)")))
}

# A function for returning a Variance-Covariance matrix for a mixed model,
# Similar to lme4:VarCorr but should return a vector of numbers rather than a list
# as VarCorr does. Currently only defined for lme4 models (class merMod)
VarCorr2 <- function (x, ...) {
  UseMethod("VarCorr2")
}

# Returns the covariance structure of an lme4 model (class merMod) as a named vector
# VarCorr returns a list, which can't be used as such in perturb
VarCorr2.merMod <- function(ModelObject) {
  vc<-as.data.frame(lme4::VarCorr(ModelObject))
  result<-vc$sdcor
  nms<-paste(vc$grp,vc$var1,vc$var2,sep="|")
  nms<-gsub("\\|NA","",nms)
  names(result)<-nms
  result
}

VarCorr2.lme<- function(ModelObject) {
  vc<-nlme::VarCorr(ModelObject)
  # vc is an x by 3 matrix
  # The first column is for variances, the second for standard deviations,the third for correlations
  # The LME output shows standard deviations, use these
  v<-vc[,2]
  # The values are empty for variables but have values for label "(Intercept)"
  # Copy the variables from previous entries
  nms<-names(v)
  x=which("(Intercept)"==nms)
  nms[x]=nms[x-1]
  names(v)<-nms
  # v contains quoted numeric values and empty strings
  # Force numeric but suppress warnings
  suppressWarnings(class(v)<-"numeric")
  # Empty cells are now NA, drop these
  v<-v[!is.na(v)]
  # Names contain equal signs, remove these
  names(v) <- gsub("^(.*)\\s[=\\s]*=\\s*$","Intercept(\\1)",names(v))
  v
}

# A function for returning a Variance-Covariance matrix for a mixed model,
# Similar to lme4:VarCorr but should return a vector of numbers rather than a list
# as VarCorr does. Currently only defined for lme4 models (class merMod)
Repeatedcoefs <- function (x, ...) {
  UseMethod("Repeatedcoefs")
}

# Doesn't work yet
# Classes are "corAR1"    "corStruct","varIdent" "varFunc"
# Can't simply extract values
Repeatedcoefs.lme <- function(ModelObject) {
  c(
    coef(ModelObject$modelStruct$corStruct,unconstrained=FALSE,allCoef=TRUE),
    coef(ModelObject$modelStruct$varStruct,unconstrained=FALSE,allCoef=TRUE)
  )
}

coefs <- function(x, ...) {
  UseMethod("coefs")
}

coefs.default <- function(obj) {
  coef(obj)
}

coefs.multinom <- function(obj) {
  gdata::unmatrix(coefficients(obj),byrow=TRUE)
}

# Returns the fixed and random effects of an lme4 model (class merMod)
# coef.default returns fixed effect coefficients but random intercepts and slopes for each case
coefs.merMod <- function(obj) {
  c(lme4::fixef(obj),VarCorr2(obj))
}

coefs.lme <- function(obj) {
  c(nlme::fixef(obj),VarCorr2(obj),Repeatedcoefs(obj))
}

# Get any warning messages from the models
GetWarnings <- function(x, ...) {
  UseMethod("GetWarnings")
}

GetWarnings.default <- function(obj) {
  NULL
}

GetWarnings.merMod <- function(obj) {
  paste(obj@optinfo$conv$lme4$messages,sep = "",collapse="; ")
}

# Get the names of variables used in the model, including random effects
# GetNames.default works for lme4::lmer, an alternative method is needed
# for nlme:lme
GetNames <- function (x, ...) {
  UseMethod("GetNames")
}
GetNames.default <- function(obj) {
  all.vars(terms(obj))
}
GetNames.lme <- function(obj) {
  c(all.vars(obj$call$fixed),all.vars(obj$call$random))
}

# Get the data frame used in the model
# Initial capital to distinguish from nlme::getData
GetData <- function (x, ...) {
  UseMethod("GetData")
}
GetData.default <- function(obj) {
  model.frame(obj)
}
GetData.lme <- function(obj) {
  obj$data
}
