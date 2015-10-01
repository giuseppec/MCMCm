##[8()
############################################################################
## sample from the posterior distribution of a factor analysis model
## model in R using linked C++ code in Scythe.
##
## The model is:
##
## Mimic model
##
## Alexander Raach, Sven Steinert
## Ludwig-Maximilians Universität München
## 25/07/2006
##
############################################################################

"MCMCm" <-
	function(ind.form, fixed.form = NA, covar.form = NA, 			###	3 Formeln für Response, Fixed Eff., Covariates
	         offset.form=NA,
 					data=parent.environment(), factors=1, 
 					lambda.constraints= NULL, burnin = 500, mcmc = 1000,
          thin=1, verbose = TRUE, seed = NA,
          beta_stern = 0, Beta_stern = 0, 						## Priori für Param.vektor beta
          gamma_stern = 0, Gamma_stern = 0, 					## Priori für Param.vektor gamma
          v_stern = 0, s_stern = 0,	## Priori für Param.vektor phi* (metrische Var.)
          beta.start = 0, gamma.start = 0,						## Startwerte für Param.vektoren beta und gamma
          store.scores = FALSE, store.z = FALSE,
          mh=FALSE, tune=0.10, gm = TRUE,						## Angabe, ob Grouped Move Step 1 durchgeführt werden soll
          DIC=FALSE, sim = 0, ... ) {
		call <- match.call()

#############################################################################
######  2	-	Daten validieren und Datenmatrizen aufbereiten						 ######
#############################################################################

cat ("Daten validieren und Datenmatrizen aufbereiten \n")
# All needed column_names are stored in the vector COLUMN_NAMES

##################################################################

# Namen der Indikatoren in COLUMN-NAMES speichern
  if (!is.na(as.list(ind.form)[1])){
	ind.terms <- terms(ind.form, data=data)
  y_sel_columns <- attr(ind.terms, "term.labels")
	del <- NULL
	for (i in 1:length(y_sel_columns)){
	     	r <- regexpr("[:\\*]",y_sel_columns[i])
		if(r[1] != -1){
			name1 <- substr(y_sel_columns[i],0,r[1]-1)
			name2 <- substr(y_sel_colomns[i],r[1]+1,100)
			y_sel_columns <- c(y_sel_columns, name1, name2)
			del <- c(del, -i)
		}
	}
	if (length(del > 0)){
		y_sel_columns <- y_sel_columns[del]
	}
	y_sel_columns <- unique(y_sel_columns)
	y_sel_columns <- gsub("\\([^\\]*\\)", "", y_sel_columns)
	COLUMN_NAMES <- c(y_sel_columns)
}

## Namen der direkten Kovariaten in COLUMN-NAMES speichern
#    if(!is.na(as.list(fixed.form)[1])) { 
#    	fixed.terms <- terms(fixed.form)  
#    	COLUMN_NAMES <- c(COLUMN_NAMES, attr(fixed.terms, "term.labels"))
#    }

################################################################
#
#     Änderung um semiparametrische Schätzung der direkten Effekte zu ermöglichen
#
    if(!is.na(as.list(fixed.form)[1])) {
    	fixed.terms <- terms(fixed.form)
    	xdir_sel_columns <- attr(fixed.terms, "term.labels")
    	DELdir <- NULL
    	for (i in 1:length(xdir_sel_columns)) {
    		r <- regexpr("[:\\*]", xdir_sel_columns[i])
    		if (r[1] != -1)
				{
					NAMEdir1 <- substr(xdir_sel_columns[i], 0, r[1]-1)
					NAMEdir2 <- substr(xdir_sel_columns[i], r[1]+1, 100)
					xdir_sel_columns <- c(xdir_sel_columns, NAMEdir1, NAMEdir2)
					DELdir <- c(DELdir, -i)
    		}
   	}
    	if (length(DELdir >0))
    		xdir_sel_columns <- xdir_sel_columns[DEL]
    	xdir_sel_columns <- unique(xdir_sel_columns)

    	xdir_sel_columns <- gsub("\\([^\\)]*\\)", "", xdir_sel_columns)
    	COLUMN_NAMES <- c(COLUMN_NAMES, xdir_sel_columns)
    }
##################################################################
        
# Namen der indirekten Kovariaten in COLUMN-NAMES speichern
    if(!is.na(as.list(covar.form)[1])){
    	covar.terms <- terms(covar.form) 
    	x_sel_columns <- attr(covar.terms, "term.labels")
    	DEL <- NULL
    	for (i in 1:length(x_sel_columns)) {
    		r <- regexpr("[:\\*]", x_sel_columns[i])
    		if (r[1] != -1)
				{
					NAME1 <- substr(x_sel_columns[i], 0, r[1]-1)
					NAME2 <- substr(x_sel_columns[i], r[1]+1, 100)
					x_sel_columns <- c(x_sel_columns, NAME1, NAME2)
					DEL <- c(DEL, -i)
    		}
    	}
    	if (length(DEL >0))    	
    		x_sel_columns <- x_sel_columns[DEL]     	
    	x_sel_columns <- unique(x_sel_columns)
			
    	x_sel_columns <- gsub("\\([^\\)]*\\)", "", x_sel_columns)
    	COLUMN_NAMES <- c(COLUMN_NAMES, x_sel_columns)
    }

###########################################################################
## Namen des Offset in COLUMN-NAMES speichern
    if(!is.na(as.list(offset.form)[1])) {
    	offset.terms <- terms(offset.form)
    	COLUMN_NAMES <- c(COLUMN_NAMES, attr(offset.terms, "term.labels"))
    }
###########################################################################


# Löschen aller Spalten in "data", die nicht benötigt werden und aller 
# Beobachtungen in "data", die ein "NA" enthalten 
    cat("Information about datamatrix \"data\":\n")
    if(missing(data)) data <- sys.frame(sys.parent())
		if(dim(data)[2]<=2)
			stop("data matrix doesn´t contain enough data columns (<=1).\n")
    data <- data[,COLUMN_NAMES]					# All not-needed columns are deleted
    # All observations with NAs are deleted
    cat("\tTotal number of observations: ", dim(data)[1],"\n")
    NotNA <- !is.na(data[,1])
    for (i in 2:dim(data)[2]) {
    	NotNA <- NotNA & !is.na(data[,i])
    }
    data <- data[NotNA,]
    rownames(data) <- 1:dim(data)[1]
    cat("\tNumber of observations without NAs: ", dim(data)[1],"\n\n")
 
#############################################################################
##
## Verarbeitung der Formel ind.form
##
if (!is.na(as.list(ind.form)[1])){
	ind.vec <- ind.form
	ind.vec <- strsplit(paste(ind.vec, sep ="")[[2]],"[\\+~]")[[1]]
	ind.vec <- gsub(" *", "" , ind.vec)
	ind.names <- 0
	ind.types <- 0
	for (i in 1:length(ind.vec)){
		r <- regexpr('^[^\\(]*', ind.vec[i], perl =TRUE)
		ind.names[i] <- substr(ind.vec[i],r[1],attr(r, "match.length") + r[1]-1)
		r <- regexpr('\\([^)]*\\)', ind.vec[i], perl=TRUE)
		r <- regexpr('\\([a-zA-Z]*[\\)]', ind.vec[i], perl = TRUE)
		ind.types[i] <- substr(ind.vec[i], r[1]+1, attr(r,"match.length")+r[1]-2)
	}
}
if (attr(ind.terms, "response") > 0){
	stop("Response not allowed in \"ind.form\" in MCMCm(). \n")
}
attributes(ind.terms)$intercept <- 0
Yterm.length <- length(ind.names)
Y <- subset(data, select=as.character(ind.names)[1:Yterm.length])
Yrow <- nrow(Y)
p <- ncol(Y)
ncat <- matrix(NA, p , 1)
p1 <-  0
tempflagc <- 0
tempflagp <- 0
p2 <- p3 <- 0
for (i in 1:p){
 	if (ind.types[i] == "c" & is.numeric(Y[,i])) {
		ncat[i] <- -999
		tempflagc <- 1
		p3 <- p3 +1
 	}
 else if (ind.types[i] == "p" & is.integer(Y[,i])) {
		if (tempflagc == 1){
			stop("Indicator in \"ind.form\" not sorted! First ordinal, second poisson, then matric variables!\n")
		}
		ncat[i] <- -999
		tempflagp <- 1
		p2 <- p2 +1
	}
 else if (ind.types[i] == "o" & is.ordered(Y[,i])) {
		if (tempflagc == 1 || tempflagp == 1){
			stop("Indicator in \"ind.form\" not sorted! First ordinal, second poisson, then matric variables!\n")
		}
		ncat[i] <- nlevels(Y[,i])
		Y[,i] <- as.integer(Y[,i])
		p1 <- p1 +1
	}
 else {
	stop("Manifest variable ", dimnames(Y)[[2]][i],
	     " neither ordered factor, nor poisson, nor numeric variable. \n")
	}
}
if (p != p1 + p2 + p3 ){
	stop("error: number of manifest variables !=  number of ordinal + number of poisson + number of continous. \n")
}

    tot_ncut <- 0
    for (i in 1:p1)
    	tot_ncut <- tot_ncut + ncat[i] - 2
    Y <- as.matrix(Y)
    Yvars <- dimnames(Y)[[2]] # Y variable names
    Yobs <- dimnames(Y)[[1]]  # observation names
   	if (is.null(Yobs)){
      Yobs <- 1:Yrow
    }
  cat ("Information about indicators y:\n\tNumber of observations: \t\t",Yrow, "\n")
  cat ("\tTotal number of indicators p: \t\t", p, "\n")
  cat ("\tNumber of ordinal indicators p1: \t", p1, "\n")
  cat ("\tNumber of Poisson-distributed indicators p2: \t", p2, "\n")
  cat ("\tNumber of metrical indicators p3: \t", p3, "\n")
  cat ("\tIndicator variable names: \t\t", Yvars, "\n")
  cat ("\tNumber of categories for each j: \t", ncat, "\n")
  cat ("\tTotal number of cutpoints: \t\t", tot_ncut, "\n\n")
  
#######################################################################

## Verarbeitung der Formel "offset.form" für die festen Effekte 1...d
    if(!is.na(as.list(offset.form)[1])) {
	    if (attr(offset.terms, "response") > 0)
 	     stop("Response not allowed in \"offset.terms\" in MCMCm().\n")
 	   attributes(offset.terms)$intercept <- TRUE
 	   Offs <- model.matrix(offset.terms, data, contrasts)
 	   Offs <- as.matrix(Offs[,-1])						# Intercept der ersten Column löschen
 	   if(is.null(dimnames(Offs)[[2]])) {
 	   		dimnames(Offs)[[2]] <- as.vector(attr(offset.terms, "term.labels"))
 	   }
 	   Offsvars <- dimnames(Offs)[[2]]
 	   Offsobs <- dimnames(Offs)[[1]]
 	   Offsrow <- nrow(Offs)	      # number of observations
 	   Offscol <- ncol(Offs)

 		} 	# endif !is.na(fixed.form)
 		else {
 			Offs <- matrix(-999, nrow=Yrow, ncol=1)
 			Offsvars <- Offsobs <- NULL
			Offsrow <- Offscol <- 0
 		}
 		if (Offsrow>1) offset <- 1 else offset <- 0
 	 cat ("Information about offset :\n\tNumber of observations: \t\t",Offsrow, "\n")
    cat ("\toffset  included: \t", offset,"\n")
   cat ("\toffset effects variable names: \t\t", Offsvars,"\n\n")
   
###########################################################################


################################################################################
## Verarbeitung der Formel "fixed.form" für die festen Effekte 1...d
#    if(!is.na(as.list(fixed.form)[1])) {
#	    if (attr(fixed.terms, "response") > 0) 
# 	     stop("Response not allowed in \"fixed.terms\" in MCMCm().\n")
# 	   nrFixedVar <- length(attr(fixed.terms, "variables"))-1
# 	   attributes(fixed.terms)$intercept <- TRUE
# 	   Omega <- model.matrix(fixed.terms, data, contrasts)
# 	   Omega <- as.matrix(Omega[,-1])						# Intercept der ersten Column löschen
# 	   if(is.null(dimnames(Omega)[[2]])) { 
# 	   		dimnames(Omega)[[2]] <- as.vector(attr(fixed.terms, "term.labels")) 
# 	   }
# 	   Omegavars <- dimnames(Omega)[[2]] 
# 	   Omegaobs <- dimnames(Omega)[[1]]
# 	   Omegarow <- nrow(Omega)	      # number of observations      
# 	   d <- ncol(Omega)        # number of fixed effects d  
# 		} 	# endif !is.na(fixed.form)   
# 		else {
# 			Omega <- matrix(-999, nrow=Yrow, ncol=1)
# 			Omegavars <- Omegaobs <- NULL
#			Omegarow <- d <- nrFixedVar <- 0
# 		}	
# 		if (Omegarow>1) Omega_exists <- 1 else Omega_exists <- 0
## 	 cat ("Information about fixed effects omega:\n\tNumber of observations: \t\t",Omegarow, "\n")
##    cat ("\tNumber of fixed effects variables: \t", nrFixedVar, "\n" )  
##    cat ("\tTotal number of fixed effects d: \t", d,"\n")
##    cat ("\tFixed effects variable names: \t\t", Omegavars,"\n\n")

################################################################
#
#     Änderung um semiparametrische Schätzung der direkten Effekte zu ermöglichen
#
# Verarbeitung der Formel "fixed.form"
		d <-0
		sdir <-0
		if (!is.na(as.list(fixed.form)[1])) {
			fixed.vec <- fixed.form
			fixed.vec <- strsplit(paste(fixed.vec, sep="")[[2]], "[\\+~]")[[1]]
			fixed.vec <- gsub(" *", "", fixed.vec)
			fixed.names 	<- fixed.names_nonpar  	<- 0				# speichert alle Covariates-Namen im Vektor
			fixed.types   <- fixed.types_nonpar 	<- 0				# speichert Typ der Regression
			fixed.options <- fixed.options_nonpar <- 0				# speichert ggf. Optionen

			for (i in 1:length(fixed.vec)) {
				r <- regexpr('^[^\\(]*', fixed.vec[i], perl=TRUE)
				fixed.names[i] <- substr(fixed.vec[i], r[1], attr(r, "match.length")+r[1]-1)
				r <- regexpr('\\([^)]*\\)', fixed.vec[i], perl=TRUE)
				if (r[1] == -1) {
					fixed.types[i] <- "std"
					fixed.options[i] <- "none"
				}
				else {
					r <- regexpr('\\([a-zA-Z0-9]*[\\),]', fixed.vec[i], perl=TRUE)
					fixed.types[i] <- substr(fixed.vec[i], r[1]+1, attr(r, "match.length")+r[1]-2)
					r <- regexpr(',[^\\)]*\\)', fixed.vec[i], perl=TRUE)
					if (r[1] == -1)
						fixed.options[i] <- "none"
					else
						fixed.options[i] <- substr(fixed.vec[i], r[1]+1, attr(r, "match.length")+r[1]-2)
				}
			}

			# Vektoren für nonpar-Variablen erstellen und Prädiktor für die fixen Effekte fixed.form_fixed herstellen
      
      fixed.form_fixed <- "~"
			for (i in 1:length(fixed.vec))
 	   		if (fixed.types[i]=="std") {
 	   			fixed.form_fixed <- paste(fixed.form_fixed, fixed.names[i],"+",sep="")
 	   			d <- d+1
    		}
  	  	else {
    			sdir <- sdir+1
    			fixed.names_nonpar[sdir] 	<- fixed.names[i]
    			fixed.types_nonpar[sdir] 	<- fixed.types[i]
    			fixed.options_nonpar[sdir] <- fixed.options[i]
    		}
    	r <- regexpr("[^=]*", fixed.form_fixed)
    	fixed.form_fixed <- substr(fixed.form_fixed, 1, attr(r, "match.length")-1)
		}
#		cat("covar.names_nonpar\n")
#		print(covar.names_nonpar)
#		cat("\ncovar.types_nonpar\n")
#  print(covar.types_nonpar)
#cat("\ncovar.options_nonpar\n")
#	print(covar.options_nonpar)

		bsplinedir <- list()
		degrdir		<- NULL  # Holder for degree of P-splines
		var_cofdir			<- NULL	 # Indicator if a var_cof is to be estimated
		# Bearbeitung der nonparametrischen Covariates
		if (sdir>0) {
			Xdir_par <- list()			# Holder for design matrices for each 1<=t<=sdir
			Kdir_par <- list()			# Holder for penalty matrices for each 1<=t<=sdir
			xdir_sorted <- list()	# Holder for different, sorted values of x for each 1<=t<=s
			intervalsdir <- NULL	# Holder for number of intervals for P-splines
			adir_par	<- NULL			# Holder for tau-prior a for each 1<=t<=sdir
			bdir_par <- NULL			# Holder for tau-prior b for each 1<=t<=sdir
			pdir_nr  <- NULL			# Holder for different parameters for design matrix column and parameter row numbers
			rank_Kdir <- NULL
			tdir <- 1
			while (tdir <= sdir) {
				if (fixed.types_nonpar[tdir] != "r1" && fixed.types_nonpar[tdir] != "r2"
				 && fixed.types_nonpar[tdir] != "p"  && fixed.types_nonpar[tdir] != "s"
				 && fixed.types_nonpar[tdir] != "v")
					stop("Unknown type of nonparametric effect ", fixed.types_nonpar[t], "!\n\n")
				xdir_par_orig	<- data[,fixed.names_nonpar[tdir]]

				if (fixed.types_nonpar[tdir] == "r1") {
					xdir_sorted[[tdir]]				<- sort(unique(xdir_par_orig))
					pdir_nr[tdir] 	<- length(xdir_sorted[[tdir]])
					Xdir_par[[tdir]]	<- matrix(0, nrow=Yrow	 , ncol=pdir_nr[tdir])
					### Creation of design matrix Xdir_par[[tdir]] for 1- and 2-order random walk
					for (i in 1:length(xdir_sorted[[tdir]]))
						Xdir_par[[tdir]][,i] <- xdir_par_orig == xdir_sorted[[tdir]][i]

					### Creation of penalty matrix Kdir_par[[tdir]] for 1-order random walk
					Kdir_par[[tdir]] <- create_K_r1(xdir_sorted[[tdir]])

					degrdir[tdir] <- -1
					var_cofdir[tdir] <- 0
				} # End random walk first order
				else if (fixed.types_nonpar[tdir] == "r2") {
					xdir_sorted[[tdir]]				<- sort(unique(xdir_par_orig))
					pdir_nr[tdir] 	<- length(xdir_sorted[[tdir]])
					Xdir_par[[tdir]]	<- matrix(0, nrow=Yrow	 , ncol=pdir_nr[tdir])
					Kdir_par[[tdir]]	<- matrix(0, nrow=pdir_nr[tdir], ncol=pdir_nr[tdir])
					### Creation of design matrix Xdir_par[[tdir]] for 1- and 2-order random walk
					for (i in 1:length(xdir_sorted[[tdir]]))
						Xdir_par[[tdir]][,i] <- xdir_par_orig == xdir_sorted[[tdir]][i]

					### Creation of penalty matrix K_par[[t]] for 2-order random walk
					Kdir_par[[tdir]] <- create_K_r2(xdir_sorted[[tdir]])

					degrdir[tdir] <- -1
					var_cofdir[tdir] <- 0
				}	# End random walk second order
				else if (fixed.types_nonpar[tdir] == "p") {
					#### Hole Werte für Anzahl Knoten (intervalsdir), Degree des Splines (degrdir) und random walk
					r <- regexpr('i=[0-9]*', fixed.options_nonpar[[tdir]], perl=TRUE)
					if (r != -1)
						intervalsdir[tdir] <- as.numeric(substr(fixed.options_nonpar[[tdir]], r[1]+2, attr(r, "match.length")+r[1]-1))
					else intervalsdir[tdir] <- 20
					r <- regexpr('d=[0-9]*', fixed.options_nonpar[[tdir]], perl=TRUE)
					if (r != -1)
						degrdir[tdir] <- as.numeric(substr(fixed.options_nonpar[[tdir]], r[1]+2, attr(r, "match.length")+r[1]-1))
					else degrdir[tdir] <- 3
					var_cofdir[tdir] <- 0
					r <- regexpr('r=[12]', fixed.options_nonpar[[tdir]], perl=TRUE)
					if (r != -1)
						psplinesdir_rw <- as.numeric(substr(fixed.options_nonpar[[tdir]], r[1]+2, attr(r, "match.length")+r[1]-1))
					else
						psplinesdir_rw <- 2

					pdir_nr[tdir] <- intervalsdir[tdir] + degrdir[tdir]
					Xdir_par[[tdir]]	<- matrix(0, nrow=Yrow	 , ncol=pdir_nr[tdir])

					xdir_min <- min(xdir_par_orig)
					xdir_max <- max(xdir_par_orig)
					#cat("Minimum value for P-splines (dir): ", xdir_min, "\n")
					#cat("Maximum value for P-splines (dir): ", xdir_max, "\n")

					cat("Start calculating B-Spline base...\n")
					### Creation of design matrix Xdir_par[[tdir]]
					for (i in 1:Yrow)
						Xdir_par[[tdir]][i,] <- bbase(xdir_par_orig[i], xdir_min, xdir_max, intervalsdir[tdir], degrdir[tdir])
					cat("For-Schleife beendet !\n")
					Xdir_par[[tdir]][abs(Xdir_par[[tdir]])<1e-8] <- 0 # Numerische Korrektur für Xdir_par[[tdir]]
					cat("Design matrix calculated (dir).....\n")
					maxidir <- 0
    			for (i in 1:dim(Xdir_par[[tdir]])[1]) {
    				row_elementsdir <- sum(Xdir_par[[tdir]][i,]>0)
    				maxidir <- max(maxidir, row_elementsdir)
    				if (row_elementsdir > degrdir[tdir]+1) {
    					cat("Bandbreite von Xdir_par[[", tdir,"]][", i, ",] = ",sum(Xdir_par[[tdir]][i,]>0),
    							 " ist größer degree+1 = ", degrdir[tdir]+1, " !!!\n")
    					return(Xdir_par[[tdir]][i,])
    				}
    			}
					#if (maxidir < degrdir[tdir]+1)
					#	stop("Maximale Anzahl Zeilenelemente ist ", maxidir, ", was kleiner ist als degrdir[tdir] + 1 = ", degrdir[tdir]+1, " !!!\n")

					### Create K-matrix with first- or second-order random walk
					dxdir <- (xdir_max-xdir_min) / intervalsdir[tdir]

					knotsdir <- 1:pdir_nr[tdir]
					cat("Bestimmte RW für P-Splines (dir) .... psplines_rw = ", psplinesdir_rw,"\n")
					if (psplinesdir_rw == 1)
						Kdir_par[[tdir]] <- create_K_r1(knotsdir)
					else
						Kdir_par[[tdir]] <- create_K_r2(knotsdir)
					cat("Number of knots(dir) = ", intervalsdir[tdir],"\n")
					cat("Degree (dir) = ", degrdir[tdir],"\n")
					xdir_sorted[[tdir]] <- 1:pdir_nr[tdir]

					b_lengthdir <- 50      ## Dieser Parameter steuert die Anzahl von function evaluations für Graphen
					bsplinedir[[tdir]] <- matrix(0, nrow=b_lengthdir, ncol=1+intervalsdir[tdir]+degrdir[tdir])
					xdir_dx <- seq(xdir_min, xdir_max, (xdir_max-xdir_min)/(b_lengthdir-1))
					for (i in 1:b_lengthdir)
						bsplinedir[[tdir]][i,] <- cbind(xdir_dx[i], bbase(xdir_dx[i], xdir_min, xdir_max, intervalsdir[tdir], degrdir[tdir]))
					bsplinedir[[tdir]][abs(bsplinedir[[tdir]])<1e-8] <- 0 # Numerische Korrektur für bsplinedir[[tdir]]
					cat("P-Splines calculation finished (dir) ...\n")
				}	# End p-splines
				else if (fixed.types_nonpar[tdir] == "s") {		## Treatment of spatial effects
					r <- regexpr('map=[^,^)]*', fixed.options_nonpar[[tdir]], perl=TRUE)
					if (r != -1)
						mapdir <- substr(fixed.options_nonpar[[tdir]], r[1]+4, attr(r, "match.length")+r[1]-1)
					else stop("No map indicated for (direct) spatial effect !\n")

					mapdir <- get(mapdir, parent.frame())
					mapdir_vec <- mapdir[[1]]
					pdir_nr[tdir] <- length(mapdir_vec)
					Xdir_par[[tdir]] <- matrix(0, nrow=Yrow, ncol=pdir_nr[tdir])
					Kdir_par[[tdir]] <- mapdir[[2]]

					### Construction of Design-Matrix Xdir_par[[tdir]] for spatial effect
					for (i in 1:pdir_nr[tdir])
						Xdir_par[[tdir]][,i] <- xdir_par_orig == mapdir_vec[i]

					xdir_sorted[[tdir]] <- mapdir_vec

					degrdir[tdir] <- -1
					var_cofdir[tdir] <- 0
				}
				else if (fixed.types_nonpar[tdir] == "v") {
					sdir_old <- sdir
					kategoriendir	<- sort(unique(xdir_par_orig))
					nr_kategoriendir <- length(kategoriendir)
					prevdir_index <- tdir - 1
					if (prevdir_index == 0)
						stop ("Only varying coefficient is not valid !\n")
					var_namedir <- fixed.names_nonpar[tdir]
					if (tdir==sdir) lastdir <- TRUE
					  else lastdir <- FALSE

					for (vdir in 2:nr_kategoriendir) {

						Kdir_par[[tdir]] <- Kdir_par[[prevdir_index]]
						degrdir[tdir]  <- degrdir[prevdir_index]
						adir_par[tdir] <- adir_par[prevdir_index]
						bdir_par[tdir] <- bdir_par[prevdir_index]
						rank_Kdir[tdir] <- rank_Kdir[prevdir_index]
						xdir_sorted[[tdir]] <- xdir_sorted[[prevdir_index]]
						pdir_nr[[tdir]] <- pdir_nr[[prevdir_index]]
						var_cofdir[tdir] <- 1														# Indicator for C++-Code not to centre this function
						zdir_vektor <- as.numeric(xdir_par_orig == kategoriendir[vdir])
						Xdir_par[[tdir]] <- zdir_vektor * Xdir_par[[prevdir_index]]

						tdir <- tdir + 1
						if (vdir > 2) sdir <- sdir + 1
					}
					rm(zdir_vektor)
					tdir <- tdir - 1

					# Variablennamen updaten im File
					if (lastdir == TRUE) {
							fixed.names_nonpar	 <- c(fixed.names_nonpar[1:prevdir_index], paste(var_namedir, "cat", 2:nr_kategoriendir, sep=""))
							fixed.types_nonpar	 <- c(fixed.types_nonpar[1:prevdir_index], rep("v", nr_kategoriendir-1) )
							fixed.options_nonpar <- c(fixed.options_nonpar[1:prevdir_index], rep("none", nr_kategoriendir-1) )
					} else {
							fixed.names_nonpar	 <- c(fixed.names_nonpar[1:prevdir_index], paste(var_namedir,"cat", 2:nr_kategoriendir, sep="-"),
																				fixed.names_nonpar[(prevdir_index+2):sdir_old])
							fixed.types_nonpar	 <- c(fixed.types_nonpar[1:prevdir_index], rep("v", nr_kategoriendir-1),
																				fixed.types_nonpar[(prevdir_index+2):sdir_old])
							fixed.options_nonpar <- c(fixed.options_nonpar[1:prevdir_index], rep("none", nr_kategoriendir-1),
																				fixed.options_nonpar[(prevdir_index+2):sdir_old])
					}


				} ## End of varying coefficients

				### Create Prioris for smoothing parameter tau
				if (fixed.types_nonpar[tdir] != "v") {
					r <- regexpr('a=[-0-9.]*', fixed.options_nonpar[[tdir]], perl=TRUE)
					if (r != -1)
						adir_par[tdir] <- as.numeric(substr(fixed.options_nonpar[[tdir]], r[1]+2, attr(r, "match.length")+r[1]-1))
					else
						adir_par[tdir] <- 0.001
					r <- regexpr('b=[-0-9.]*', fixed.options_nonpar[[tdir]], perl=TRUE)
					if (r != -1)
						bdir_par[tdir] <- as.numeric(substr(fixed.options_nonpar[[tdir]], r[1]+2, attr(r, "match.length")+r[1]-1))
					else
						bdir_par[tdir] <- 0.001

					# Create ranks of K
					rank_Kdir[tdir]		<- qr(Kdir_par[[tdir]])$rank
				}
				tdir <- tdir + 1
			} # end while-loop for tdir <= sdir

		}

#		cat("\n\ncovar.names_nonpar\n")
#		print(covar.names_nonpar)
#		cat("\ncovar.types_nonpar\n")
#		print(covar.types_nonpar)
#		cat("\ncovar.options_nonpar\n")
#		print(covar.options_nonpar)

    if (d>0) {
	    # Matrix X für fixen Effekte mit Anzahl d erstellen
 	   	fixed.terms <- terms(as.formula(fixed.form_fixed))
		  if (attr(fixed.terms, "response") > 0)
 		 	  stop("Response not allowed in \"fixed.form\" in MCMCm().\n")
 	  	nrCovardir <- length(attr(fixed.terms, "variables"))-1
	 	 	attributes(fixed.terms)$intercept <- TRUE
 		 	Omega <- model.matrix(fixed.terms, data, contrasts)
 			Omega <- as.matrix(Omega[,-1])						# Intercept der ersten Column rausschmeissen
#X <- as.matrix(X[,-4])
#return(X)
	 	 	if(is.null(dimnames(Omega)[[2]]))
 		 		dimnames(Omega)[[2]] <- as.vector(attr(fixed.terms, "term.labels"))

 		 	Omegavars <- dimnames(Omega)[[2]]
 		 	Omegaobs <- dimnames(Omega)[[1]]
 		 	Omegarow <- nrow(Omega)	      # number of observations
	 	 	d <- ncol(Omega)        # number of covariates d
		} else {   ## d=0
 			Omega <- matrix(-999, nrow=Yrow, ncol=1)
 			Omegavars <- Omegaobs <- NULL
			Omegarow <- d <- nrCovardir <- 0
 		}
 		if ( d>0) Omega_exists 		 <- 1 else 		 Omega_exists <- 0
 		if ( sdir>0) Xdir_par_exists <- 1 else Xdir_par_exists <- 0

# 		return(X)
 		##### Hier temporäre Effektkodierung eingebaut !!! #####
# 		for (xz in 1:dim(X)[1]) {
# 			if (sum(X[xz,])==0) {
# 				X[xz,] <- -1
 			#cat("Umgewandelt für i = ", i)
# 			}
# 		}
 		########################################################


 		# Output wesentlicher Größen des Covariate-Prädiktors für latente Variablen
 		cat ("Information about direct covariates xdir:\n\tNumber of observations: \t\t",Omegarow, "\n")
    cat ("\tTotal number of direct covariates d+sdir: \t\t", d+sdir,"\n")
    cat ("\tNumber of direct nonparametric covariates sdir: \t\t", sdir,"\n")
    cat ("\tNumber of direct fixed effects covariates d: \t\t", d,"\n")



################################################################################
# Verarbeitung der Formel "covar.form" für die Covariates der latenten Variablen 1...q
		q <-0
		s <-0
		if (!is.na(as.list(covar.form)[1])) {
			covar.vec <- covar.form
			covar.vec <- strsplit(paste(covar.vec, sep="")[[2]], "[\\+~]")[[1]]
			covar.vec <- gsub(" *", "", covar.vec)
			covar.names 	<- covar.names_nonpar  	<- 0				# speichert alle Covariates-Namen im Vektor
			covar.types   <- covar.types_nonpar 	<- 0				# speichert Typ der Regression
			covar.options <- covar.options_nonpar <- 0				# speichert ggf. Optionen
			for (i in 1:length(covar.vec)) {
				r <- regexpr('^[^\\(]*', covar.vec[i], perl=TRUE)
				covar.names[i] <- substr(covar.vec[i], r[1], attr(r, "match.length")+r[1]-1)
				r <- regexpr('\\([^)]*\\)', covar.vec[i], perl=TRUE)
				if (r[1] == -1) {
					covar.types[i] <- "std"
					covar.options[i] <- "none"
				}
				else {
					r <- regexpr('\\([a-zA-Z0-9]*[\\),]', covar.vec[i], perl=TRUE)
					covar.types[i] <- substr(covar.vec[i], r[1]+1, attr(r, "match.length")+r[1]-2)
					r <- regexpr(',[^\\)]*\\)', covar.vec[i], perl=TRUE)
					if (r[1] == -1)
						covar.options[i] <- "none"
					else
						covar.options[i] <- substr(covar.vec[i], r[1]+1, attr(r, "match.length")+r[1]-2)
				}
			}
			# Vektoren für nonpar-Variablen erstellen und Prädiktor für die fixen Effekte covar.form_fixed herstellen 
			covar.form_fixed <- "~"
			for (i in 1:length(covar.vec)) 
 	   		if (covar.types[i]=="std") {
 	   			covar.form_fixed <- paste(covar.form_fixed, covar.names[i],"+",sep="")
 	   			q <- q+1
    		}
  	  	else {
    			s <- s+1
    			covar.names_nonpar[s] 	<- covar.names[i]	
    			covar.types_nonpar[s] 	<- covar.types[i]	
    			covar.options_nonpar[s] <- covar.options[i]	
    		}
    	r <- regexpr("[^=]*", covar.form_fixed)
    	covar.form_fixed <- substr(covar.form_fixed, 1, attr(r, "match.length")-1)
		}
#		cat("covar.names_nonpar\n")
#		print(covar.names_nonpar)
#		cat("\ncovar.types_nonpar\n")
#   print(covar.types_nonpar)
#   cat("\ncovar.options_nonpar\n")
#	  print(covar.options_nonpar)

###############################################################################
#
#

# Vergleich: direkte und indirekte Kovariaten
if((d+sdir)>0 & (q+s)>0){
vglindikator1 <- 0
for (i in 1:length(fixed.names)) {
  for ( j in 1: length(covar.names)){
    if (fixed.names[i] == covar.names[j]){
      vglindikator1 <- 1
      stop("Variable nicht zugleich direkte und indirekte Kovariate!")
    }
  }
}
}

# Vergleich: direkte Kovariaten und Indikatoren
if(p>0 && (d+sdir)>0){
vglindikator2 <- 0
for (i in 1:length(fixed.names)) {
  for ( j in 1: length(ind.names)){
    if (fixed.names[i] == ind.names[j]){
      vglindikator2 <- 1
      stop("Variable nicht zugleich direkte Kovariable und Indikator!")
    }
  }
}
}
# Vergleich: indirekte Kovariaten und Indikatoren
if(p>0 && (q+s)>0){
vglindikator2 <- 0
for (i in 1:length(covar.names)) {
  for ( j in 1: length(ind.names)){
    if (covar.names[i] == ind.names[j]){
      vglindikator2 <- 1
      stop("Variable nicht zugleich indirekte Kovariable und Indikator!")
    }
  }
}
}

#
#
##############################################################################


		bspline <- list()
		degr		<- NULL  # Holder for degree of P-splines
		var_cof			<- NULL	 # Indicator if a var_cof is to be estimated
		# Bearbeitung der nonparametrischen Covariates
		if (s>0) {
			X_par <- list()			# Holder for design matrices for each 1<=t<=s
			K_par <- list()			# Holder for penalty matrices for each 1<=t<=s
			x_sorted <- list()	# Holder for different, sorted values of x for each 1<=t<=s
			intervals <- NULL	# Holder for number of intervals for P-splines
			a_par	<- NULL			# Holder for tau-prior a for each 1<=t<=s
			b_par <- NULL			# Holder for tau-prior b for each 1<=t<=s
			p_nr  <- NULL			# Holder for different parameters for design matrix column and parameter row numbers
			rank_K <- NULL
			t <- 1
			while (t <= s) {
				if (covar.types_nonpar[t] != "r1" && covar.types_nonpar[t] != "r2" 
				 && covar.types_nonpar[t] != "p"  && covar.types_nonpar[t] != "s"
				 && covar.types_nonpar[t] != "r"
				 && covar.types_nonpar[t] != "v")
					stop("Unknown type of nonparametric effect ", covar.types_nonpar[t], "!\n\n")
				x_par_orig	<- data[,covar.names_nonpar[t]]

				if (covar.types_nonpar[t] == "r1") {
					x_sorted[[t]]				<- sort(unique(x_par_orig))
					p_nr[t] 	<- length(x_sorted[[t]]) 
					X_par[[t]]	<- matrix(0, nrow=Yrow	 , ncol=p_nr[t])
					### Creation of design matrix X_par[[t]] for 1- and 2-order random walk
					for (i in 1:length(x_sorted[[t]]))
						X_par[[t]][,i] <- x_par_orig == x_sorted[[t]][i]

					### Creation of penalty matrix K_par[[t]] for 1-order random walk
					K_par[[t]] <- create_K_r1(x_sorted[[t]])
					
					degr[t] <- -1
					var_cof[t] <- 0
				} # End random walk first order
				else if (covar.types_nonpar[t] == "r2") {
					x_sorted[[t]]				<- sort(unique(x_par_orig))
					p_nr[t] 	<- length(x_sorted[[t]]) 
					X_par[[t]]	<- matrix(0, nrow=Yrow	 , ncol=p_nr[t])
					K_par[[t]]	<- matrix(0, nrow=p_nr[t], ncol=p_nr[t])
					### Creation of design matrix X_par[[t]] for 1- and 2-order random walk
					for (i in 1:length(x_sorted[[t]]))
						X_par[[t]][,i] <- x_par_orig == x_sorted[[t]][i]
					
					### Creation of penalty matrix K_par[[t]] for 2-order random walk				
					K_par[[t]] <- create_K_r2(x_sorted[[t]])
					
					degr[t] <- -1				
					var_cof[t] <- 0						
				}	# End random walk second order
				else if (covar.types_nonpar[t] == "p") {
					#### Hole Werte für Anzahl Knoten (intervals), Degree des Splines (degr) und random walk
					r <- regexpr('i=[0-9]*', covar.options_nonpar[[t]], perl=TRUE)
					if (r != -1)
						intervals[t] <- as.numeric(substr(covar.options_nonpar[[t]], r[1]+2, attr(r, "match.length")+r[1]-1))
					else intervals[t] <- 20
					r <- regexpr('d=[0-9]*', covar.options_nonpar[[t]], perl=TRUE)
					if (r != -1)
						degr[t] <- as.numeric(substr(covar.options_nonpar[[t]], r[1]+2, attr(r, "match.length")+r[1]-1))
					else degr[t] <- 3
					var_cof[t] <- 0
					r <- regexpr('r=[12]', covar.options_nonpar[[t]], perl=TRUE)					
					if (r != -1)
						psplines_rw <- as.numeric(substr(covar.options_nonpar[[t]], r[1]+2, attr(r, "match.length")+r[1]-1))
					else 
						psplines_rw <- 2
						
					p_nr[t] <- intervals[t] + degr[t]
					X_par[[t]]	<- matrix(0, nrow=Yrow	 , ncol=p_nr[t])																
										
					x_min <- min(x_par_orig)
					x_max <- max(x_par_orig)
					#cat("Minimum value for P-splines: ", x_min, "\n")
					#cat("Maximum value for P-splines: ", x_max, "\n")
					
					cat("Start calculating B-Spline base...\n")
					### Creation of design matrix X_par[[t]]
					for (i in 1:Yrow)
						X_par[[t]][i,] <- bbase(x_par_orig[i], x_min, x_max, intervals[t], degr[t])
					cat("For-Schleife beendet !\n")
					X_par[[t]][abs(X_par[[t]])<1e-8] <- 0 # Numerische Korrektur für X_par[[t]] 
					cat("Design matrix calculated.....\n")
					maxi <- 0
    			for (i in 1:dim(X_par[[t]])[1]) {
    				row_elements <- sum(X_par[[t]][i,]>0)
    				maxi <- max(maxi, row_elements)
    				if (row_elements > degr[t]+1) {
    					cat("Bandbreite von X_par[[", t,"]][", i, ",] = ",sum(X_par[[t]][i,]>0), 
    							 " ist größer degree+1 = ", degr[t]+1, " !!!\n") 
    					return(X_par[[t]][i,])
    				}
    			}
					#if (maxi < degr[t]+1) 
					#	stop("Maximale Anzahl Zeilenelemente ist ", maxi, ", was kleiner ist als degr[t] + 1 = ", degr[t]+1, " !!!\n")
					
					### Create K-matrix with first- or second-order random walk
					dx <- (x_max-x_min) / intervals[t]
					#knots <- seq(0, (p_nr[t]-1) * dx, by = dx)
					knots <- 1:p_nr[t]
					cat("Bestimmte RW für P-Splines.... psplines_rw = ", psplines_rw,"\n")
					if (psplines_rw == 1)
						K_par[[t]] <- create_K_r1(knots)
					else
						K_par[[t]] <- create_K_r2(knots)							
					cat("Number of knots = ", intervals[t],"\n")
					cat("Degree = ", degr[t],"\n")
					x_sorted[[t]] <- 1:p_nr[t]
					
					b_length <- 50      ## Dieser Parameter steuert die Anzahl von function evaluations für Graphen
					bspline[[t]] <- matrix(0, nrow=b_length, ncol=1+intervals[t]+degr[t])
					x_dx <- seq(x_min, x_max, (x_max-x_min)/(b_length-1)) 
					for (i in 1:b_length) 
						bspline[[t]][i,] <- cbind(x_dx[i], bbase(x_dx[i], x_min, x_max, intervals[t], degr[t]))
					bspline[[t]][abs(bspline[[t]])<1e-8] <- 0 # Numerische Korrektur für bspline[[t]]					
					cat("P-Splines calculation finished...\n")
				}	# End p-splines
				else if (covar.types_nonpar[t] == "s") {		## Treatment of spatial effects
					r <- regexpr('map=[^,^)]*', covar.options_nonpar[[t]], perl=TRUE)
					if (r != -1)
						map <- substr(covar.options_nonpar[[t]], r[1]+4, attr(r, "match.length")+r[1]-1)
					else stop("No map indicated for spatial effect !\n")

					map <- get(map, parent.frame())
					map_vec <- map[[1]]
					p_nr[t] <- length(map_vec)
					X_par[[t]] <- matrix(0, nrow=Yrow, ncol=p_nr[t]) 
					K_par[[t]] <- map[[2]]
					
					### Construction of Design-Matrix X_par[[t]] for spatial effect
					for (i in 1:p_nr[t])
						X_par[[t]][,i] <- x_par_orig == map_vec[i]
					
					x_sorted[[t]] <- map_vec
					
					degr[t] <- -1
					var_cof[t] <- 0
				}
				else if (covar.types_nonpar[t]=="r"){
				  if(t==1){
				    stop("First unstructured not allowed \n")
		      }
		      if(covar.types_nonpar[(t-1)] != "s"){
		        stop("stuctured spatial effect must be before unstructured \n")
          }
          if (covar.names_nonpar[(t-1)] != covar.names_nonpar[t]){
            stop("Different variables for structured and unstructured spatial effect not allowed. \n")
          }
 					# map <- get(map, parent.frame())
					# map_vec <- map[[1]]
					p_nr[t] <- length(map_vec)
					X_par[[t]] <- matrix(0, nrow=Yrow, ncol=p_nr[t]) 
					K_par[[t]] <- diag(length(map_vec))
					
					### Construction of Design-Matrix X_par[[t]] for unstructured spatial effect
					for (i in 1:p_nr[t])
						X_par[[t]][,i] <- x_par_orig == map_vec[i]
					
					x_sorted[[t]] <- map_vec
					
					degr[t] <- -1
					var_cof[t] <- 0
				}
		      
				else if (covar.types_nonpar[t] == "v") {
					s_old <- s
					kategorien	<- sort(unique(x_par_orig))
					nr_kategorien <- length(kategorien)
					prev_index <- t - 1
					if (prev_index == 0) 
						stop ("Only varying coefficient is not valid !\n")
					var_name <- covar.names_nonpar[t]
					if (t==s) last <- TRUE
					  else last <- FALSE
					for (v in 2:nr_kategorien) {
						K_par[[t]] <- K_par[[prev_index]]
						degr[t]  <- degr[prev_index]
						a_par[t] <- a_par[prev_index]
						b_par[t] <- b_par[prev_index]
						rank_K[t] <- rank_K[prev_index]
						x_sorted[[t]] <- x_sorted[[prev_index]]
						p_nr[[t]] <- p_nr[[prev_index]]							
						var_cof[t] <- 1														# Indicator for C++-Code not to centre this function
						z_vektor <- as.numeric(x_par_orig == kategorien[v])
						X_par[[t]] <- z_vektor * X_par[[prev_index]]
						t <- t + 1
						if (v > 2) s <- s + 1
					}
					rm(z_vektor)
					t <- t - 1
					# Variablennamen updaten im File
					if (last == TRUE) {
							covar.names_nonpar	 <- c(covar.names_nonpar[1:prev_index], paste(var_name, "cat", 2:nr_kategorien, sep=""))
							covar.types_nonpar	 <- c(covar.types_nonpar[1:prev_index], rep("v", nr_kategorien-1) )
							covar.options_nonpar <- c(covar.options_nonpar[1:prev_index], rep("none", nr_kategorien-1) )
					} else {
							covar.names_nonpar	 <- c(covar.names_nonpar[1:prev_index], paste(var_name,"cat", 2:nr_kategorien, sep="-"),
																				covar.names_nonpar[(prev_index+2):s_old])
							covar.types_nonpar	 <- c(covar.types_nonpar[1:prev_index], rep("v", nr_kategorien-1),
																				covar.types_nonpar[(prev_index+2):s_old])
							covar.options_nonpar <- c(covar.options_nonpar[1:prev_index], rep("none", nr_kategorien-1),
																				covar.options_nonpar[(prev_index+2):s_old])																		
					}
				} ## End of varying coefficients
				
				### Create Prioris for smoothing parameter tau
				if (covar.types_nonpar[t] != "v") {
					r <- regexpr('a=[-0-9.]*', covar.options_nonpar[[t]], perl=TRUE)
					if (r != -1)
						a_par[t] <- as.numeric(substr(covar.options_nonpar[[t]], r[1]+2, attr(r, "match.length")+r[1]-1))
					else 
						a_par[t] <- 0.001					
					r <- regexpr('b=[-0-9.]*', covar.options_nonpar[[t]], perl=TRUE)
					if (r != -1)
						b_par[t] <- as.numeric(substr(covar.options_nonpar[[t]], r[1]+2, attr(r, "match.length")+r[1]-1))
					else 
						b_par[t] <- 0.001								
				
					# Create ranks of K
					rank_K[t]		<- qr(K_par[[t]])$rank
				}
				t <- t + 1
			} # end while (t <= s)
		} # end if(s>0)
		
#		cat("\n\ncovar.names_nonpar\n")
#		print(covar.names_nonpar)
#		cat("\ncovar.types_nonpar\n")
#		print(covar.types_nonpar)
#		cat("\ncovar.options_nonpar\n")
#		print(covar.options_nonpar)		
						
    if (q>0) {
	    # Matrix X für fixen Effekte mit Anzahl q erstellen
 	   	covar.terms <- terms(as.formula(covar.form_fixed))
		  if (attr(covar.terms, "response") > 0) 
 		 	  stop("Response not allowed in \"covar.ind\" in MCMCm().\n")
 	  	nrCovar <- length(attr(covar.terms, "variables"))-1
	 	 	attributes(covar.terms)$intercept <- TRUE
 		 	X <- model.matrix(covar.terms, data, contrasts)
 			X <- as.matrix(X[,-1])						# Intercept der ersten Column rausschmeissen
      if(is.null(dimnames(X)[[2]]))
 		 		dimnames(X)[[2]] <- as.vector(attr(covar.terms, "term.labels")) 
 		 	
 		 	Xvars <- dimnames(X)[[2]] 
 		 	Xobs <- dimnames(X)[[1]]
 		 	Xrow <- nrow(X)	      # number of observations      
	 	 	q <- ncol(X)        # number of covariates q  
		} else {   ## q=0
 			X <- matrix(-999, nrow=Yrow, ncol=1)
 			Xvars <- Xobs <- NULL
			Xrow <- q <- nrCovar <- 0
 		}	
 		if ( q>0) X_exists 		 <- 1 else 		 X_exists <- 0
 		if ( s>0) X_par_exists <- 1 else X_par_exists <- 0
 		

 		##### Hier temporäre Effektkodierung eingebaut !!! #####
# 		for (xz in 1:dim(X)[1]) {
# 			if (sum(X[xz,])==0) {
# 				X[xz,] <- -1
 			#cat("Umgewandelt für i = ", i)
# 			}
# 		}
 		########################################################		
 		
 		
 		# Output wesentlicher Größen des Covariate-Prädiktors für latente Variablen
 		cat ("Information about covariates x:\n\tNumber of observations: \t\t",Xrow, "\n")
    cat ("\tTotal number of covariates q+s: \t\t", q+s,"\n")
    cat ("\tNumber of nonparametric covariates s: \t\t", s,"\n")
    cat ("\tNumber of fixed effects covariates q: \t\t", q,"\n")
    #return(X)
#    cat ("\tNumber of fixed effects covariate variables: \t\t", nrCovar, "\n" )  
#    cat ("\tNonparametric effects variable names: \t\t", Xvars,"\n\n")


#############################################################################
######  3.1	-	Prioris validieren und festlegen 												 ######
#############################################################################

		## Parametervektor gamma priori festlegen (enthält GAMMA)		
		ngamma <- factors*q
		gamma_prior <- form.mvn.prior(gamma_stern, Gamma_stern, q)
		# gamma_prior <- form.mvn.prior(gamma_stern, Gamma_stern, q+1) #Änderung gamma_intercept
		gamma_stern_call <- gamma_stern
		Gamma_stern_call <- Gamma_stern
		
		gamma_stern <- gamma_prior[[1]]
		Gamma_stern <- gamma_prior[[2]]

		## Parametervektor beta priori festlegen (enthält b_0, A, B)
		if (factors < 1) stop("Number of factors have to integer and >= 1!\n\n") 
	 	if (factors >= p) stop("Number of factors have to be smaller than number of indicators p !\n\n")         # hier auskommentieren falls nur ein Indikator 
    nbeta <- p*(factors+d+1)
		beta_prior <- form.mvn.prior(beta_stern, Beta_stern, nbeta)
		beta_stern_call <- beta_stern
		Beta_stern_call <- Beta_stern
		beta_stern <- beta_prior[[1]]
		Beta_stern <- beta_prior[[2]]
	
		
		## Parametervektor phi priori festlegen (enthält phi_j für metrische Var.)
		if (p3 > 0) {
			phi_prior <- form.mvn.prior(v_stern, 0, p3)
			v_stern <- phi_prior[[1]]
			phi_prior = form.mvn.prior(s_stern, 0, p3)
			s_stern <- phi_prior[[1]]
		} else {
			v_stern <- matrix(-999, nrow=1, ncol=1)
			s_stern <- matrix(-999, nrow=1, ncol=1)
		}
		
		## Print-Ausgabe Prioris
#		cat ("Number of latent factors m: \t", factors, "\n\n")
#		cat ("Length of parameter vector beta (p*(m+d+1)): \t", nbeta, "\n")
#		cat ("Length of parameter vector gamma (m*q): \t", ngamma, "\n")	
#		cat ("Length of parameter vector phi (2*p3): \t", 2*p3, "\n\n")

#############################################################################
######  3.2	-	Startwerte für Param.vektoren beta und gamma festlegen 	 ######
#############################################################################

##### Matrizen für Startwerte initialisieren (für gamma0 nicht notwendig)
		beta_start 										<- matrix(		0,	 nbeta,						1)
		gamma_start 									<- matrix(		0,  ngamma,				 		1)
		b0_start 											<- matrix(		0,			 p,           1)
    
#    if (Omega_exists) {
#      A_start 		<- matrix(0, p, d)
#      if(offset) {A_start[,1] <- 1}
#      }
#		else	A_start <- matrix( -999, 1, 1)
					
		if (Omega_exists)	A_start 		<- matrix(		0, 			 p, 					d)
			else						A_start 		<- matrix( -999,  		 1, 					1)
			
		B_start 											<- matrix(		0,   		 p,     factors)			
		phi_start 										<- matrix(		0,			 p,           p)
		if (p1==0)				nu_start		<- matrix( -999,			 p, 					1) 
			else 						nu_start		<- matrix( -999,			 p, max(ncat)+1)
		eta_start											<- matrix(		0, 	  Yrow,  		factors)

spa <- matrix(NA,2,p)
spa[1,] <- rep(1,p)
spa[2,] <- (1:p)
if(p2 >0){
 	for(j in c((p1+1):(p1 + p2))) {
		for(i in c(1:Yrow)){
  			if(spa[1,j] < max(Y[1:Yrow,j]+1)){
				spa[1,j] <- max(Y[1:Yrow,j]+1)
			}
		}
 	}
	for(k in (p1+1): p){
	 spa[2,k] <- sum(spa[1,1:k])
	}
}
spaanz <- spa[2,p]
z_start <- matrix(0, Yrow, spaanz)
##############################################################################
	
##### Wenn Startwerte gegeben sind, diese verwenden für gesamt beta und gamma
		beta.start.data <- form.mvn.prior(beta.start, 0, nbeta)
		beta_start <- beta.start.data[[1]]
			
		gamma.start.data <- form.mvn.prior(gamma.start, 0, ngamma)
		gamma_start <- gamma.start.data[[1]]			
		
##### Startwerte für b0, A und B aus beta.start ableiten
		for (j in 1:p) {
			b0_start[j,1] <- beta_start[(j-1)*(factors+d+1)+1,1]
			if (Omega_exists) {
				for (i in 1:d)
					A_start[j,i] <- beta_start[(j-1)*(factors+d+1)+1+i,1]
			}
			for (i in 1:factors)
				B_start[j,i] <- beta_start[(j-1)*(factors+d+1)+1+d+i,1]				
			## For B_start wird 0.7 gewählt für die erste freie Faktorladung, wenn beta.start=0 ist (default)
			if (beta.start == 0) 
				for (j in 1:min(p,factors))
					B_start[j,j] <- 0.7
			
		}	

####### Startwerte für gamma erzeugen
		if (is.na(gamma.start)[1]) {
			if (X_exists){
				for (i in 1:factors) {
					gamma_start[(i-1)*q+1,1] = 0
					for (j in 1:q) 
						gamma_start[(i-1)*q+1+j,1] <- 0
				}
			}
		}
				
####### Startwerte für eta (latente Faktoren) berechnen
		
		if (X_exists && q>0)    # Hier fehlen noch die Additionen für die nonparametr. Effekte
			eta_start <- X %*% t(matrix(gamma_start, nrow=factors, ncol=q, byrow=TRUE))

####### Startwerte für phi (Varianzen für metrische Var.) berechnen
		for (i in 1:p)
				phi_start[i,i] = 1  	# werden 1 gesetzt
		
####### Startwerte for cutpoints nu (cutpoints go from 0 to ncat[i]-2)
    if (p1>0) {
    	for (i in 1:p1){
     	  nu_start[i,1] <- -10000000
     	  nu_start[i,2] <- 0
     		nu_start[i,ncat[i]+1] <- 10000000
     	 	if(ncat[i] > 2)
     	  	nu_start[i,3:(ncat[i])] <- 1:(ncat[i]-2)
    	}
    }	

####### Startwerte für z (unterliegende Variablen) berechnen
		# Für ordinale und poissonverteilte Variablen nicht notwendig, da Durchführung im ersten Gibbs-Schritt
		# For continuous indicators the actual response is saved in z_start
		if (p3>0){
			for (j in (p1 + p2 +1):p){
				for (i in 1:Yrow)
					z_start[i,sum(spa[1,1:j])] = Y[i,j]
			}
		}

#######################################################################################
####### Startwerte für beta aus b0_start, A_start, B_start erzeugen
		if (is.na(beta.start)[1]) {
			for (i in 1:p) {
				beta_start[(i-1)*(factors+d+1)+1,1] <- b0_start[i,1]
				if (Omega_exists) 
					beta_start[(i-1)*(factors+d+1)+1+1:d,1] <- as.vector(A_start[i,])
				beta_start[(i-1)*(factors+d+1)+1+d+1:factors,1] <- as.vector(B_start[i,])
			}
		}

####### Startwerte ausdrucken
		#cat("b0_start:\n")
		#print(b0_start, quote=FALSE)
		#cat("B_start:\n")
		#print(B_start, quote=FALSE)
		#cat("A_start:\n")
		#print(A_start, quote=FALSE)
		#cat("phi_start:\n")
		#print(phi_start, quote=FALSE)
		#cat("nu_start:\n")
		#print(nu_start, quote=FALSE)
		#cat("beta_start:\n")
		#print(beta_start, quote=FALSE)
		#cat("gamma_start:\n")
		#print(gamma_start, quote=FALSE)
		#return(sem.object)
	
#############################################################################
######  4	-	Verwaltungsarbeiten(MCMC-Parameter festlegen, init rseed() ######
#############################################################################
		 
    # Check MCMC parameters
    check.mcmc.parameters(burnin, mcmc, thin)
    
    # seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]
    
    # Check and initialize MH tuning parameter
    if (p1>0) { 					# Nur für ordinale Variablen durchführen
   	 if (is.double(tune)){
   	 	 tune_original <- tune
   	   tune <- matrix(tune/ncat[1:p1], nrow=p1, ncol=1)
   	 }
   	 else stop("\"tune\"-parameter is not a double digit !\n") 
    }
          
    # define holder for posterior density sample
    sample_size <- 0
    if (store.scores == TRUE) sample_size <- sample_size + Yrow*factors    	
    sample_size <- sample_size + nbeta
################################################################################
#
#     Änderung um semiparametrische Schätzung der direkten Effekte zu ermöglichen
   if (Xdir_par_exists) {
      for (i in 1:sdir) sample_size <- sample_size + pdir_nr[i]*p + p # +p für tau-Variablen
    }
################################################################################
    if (X_exists) sample_size <- sample_size + ngamma
    if (X_par_exists) {
    	for (i in 1:s) sample_size <- sample_size + p_nr[i]*factors + 1*factors   # 1*factors für tau-Variablen 
    }
    sample_size <- sample_size + p		# theta is saved for ordinal and metric variables
    if ((p1 > 0) && (ncol(nu_start) > 3)){
  		for (i in 1:p)
  			if (ncat[i] > 2) sample_size <- sample_size + (ncat[i] - 2);
  	}

	if(store.z == TRUE){
		sample_size <- sample_size + Yrow*(spaanz-p3)
	}

##########################################################################################

  	# EIN UM EIN ERHÖHTES SAMPLESIZE FÜR DEN GAMMA-INTERCEPT
  	if (X_exists || X_par_exists) sample_size <- sample_size + 1
   	if(DIC == TRUE) sample_size <- sample_size + 1
   	 	
    sample <- matrix(data=0, mcmc/thin, sample_size)
    cat("R - Size of sample: nrows = ", mcmc/thin, "   ncols = ", sample_size, "\n")
    
#############################################################################
######  5	-	Aufruf C++ Code																						 ######
#############################################################################
        														 
    ### Matrices for nonparametric effects
		X_par_total <- NULL
		K_par_total <- NULL
		nonpar_par_total <- NULL
		if (s>0) {
			max_parameter <- max(p_nr)
			nonpar_par_total <- rbind(p_nr, a_par, b_par, rank_K, degr, var_cof) # contains prioris a/b, nr of parameters, ranks, degrees, var_cof
			
			for (j in 1:s) {
				X_par_total <- cbind(X_par_total, X_par[[j]])
				K_par_temp <- matrix(-876, nrow=max_parameter, ncol=p_nr[j])
				K_par_temp[1:p_nr[j], 1:p_nr[j]] <- K_par[[j]]
				K_par_total <- cbind(K_par_total, K_par_temp)
			}
		} else {
			  X_par_total <- matrix(-999, nrow=Yrow, ncol=1)
			  K_par_total <- matrix(-999, nrow=1, ncol=1)
			  nonpar_par_total <- matrix(-999, nrow=1, ncol=1)
		}
################################################################################
#
#     Änderung um semiparametrische Schätzung der direkten Effekte zu ermöglichen
#

    ### Matrices for nonparametric direct effects
		Xdir_par_total <- NULL
		Kdir_par_total <- NULL
		nonpardir_par_total <- NULL
		if (sdir>0) {
			maxdir_parameter <- max(pdir_nr)
			nonpardir_par_total <- rbind(pdir_nr, adir_par, bdir_par, rank_Kdir, degrdir, var_cofdir) # contains prioris adir/bdir, nr of parameters, ranks, degrees, var_cof
			
			for (j in 1:sdir) {
				Xdir_par_total <- cbind(Xdir_par_total, Xdir_par[[j]])
				Kdir_par_temp <- matrix(-876, nrow=maxdir_parameter, ncol=pdir_nr[j])
				Kdir_par_temp[1:pdir_nr[j], 1:pdir_nr[j]] <- Kdir_par[[j]]
				Kdir_par_total <- cbind(Kdir_par_total, Kdir_par_temp)
			}
		} else {
			  Xdir_par_total <- matrix(-999, nrow=Yrow, ncol=1)
			  Kdir_par_total <- matrix(-999, nrow=1, ncol=1)
			  nonpardir_par_total <- matrix(-999, nrow=1, ncol=1)
		}

################################################################################
	
	if (  !((nrow(Y) == nrow(X)) && (nrow(Y) == nrow(X_par_total)) && (nrow(Y) == nrow(Omega)) &&
			(nrow(Y) == nrow(z_start)) && (nrow(Y) == nrow(eta_start))) )
	stop ("Number of rows of Y, X, Omega, X_par_total, z_start and eta_start are not all equal !!!\n\n")
		if (  !((nrow(Y) == nrow(Xdir_par_total))) )
	stop ("Number of rows of Y and Xdir_par_total are not all equal !!!\n\n")
	if (  !((nrow(Y) == nrow(Offs))) )
	stop ("Number of rows of Y and Offs are not all equal !!!\n\n")

    if (is.null(lambda.constraints)){
      Lncol <- factors
      Lnrow <- nrow(B_start)
      lambda.constraints <- matrix(1,ncol=Lncol,nrow=Lnrow)
      if(Lncol != 1) {
        for(i in 2:ncol(lambda.constraints)){
          lambda.constraints[1:(i-1),i] <- 0
        }
      }
    } else {
      lambda.constraints <- do.call("cbind", lambda.constraints)
    }
		### MCMC parameter Übergabevektor definieren
		MCMCparameter <- as.vector(c(as.integer(burnin),
                    as.integer(mcmc),
                    as.integer(thin),
                    as.integer(lecuyer),
                    as.integer(lecuyer.stream),
                    as.integer(verbose),
                    as.integer(store.scores),
                    as.integer(store.z),
                    as.integer(gm),
                    as.integer(mh),
                    as.integer(p1),
                    as.integer(factors),
                    as.integer(d),
                    as.integer(q),
                    as.integer(s),
                    as.integer(Omega_exists),
                    as.integer(X_exists),
                    as.integer(X_par_exists),
                    as.integer(DIC),
                    as.integer(sim),
                    as.integer(nrow(ncat)),
                    as.integer(ncol(ncat)),
                    as.integer(p3),
                    as.integer(p2),
                    as.integer(ncol(z_start)),     # Spaltenzal von z_start als MCMC-Parameter übergeben
                    as.integer(Xdir_par_exists),
                    as.integer(nrow(beta_stern)),
                    as.integer(ncol(beta_stern)),  
                    as.integer(nrow(Beta_stern)),
                    as.integer(ncol(Beta_stern)),
                    as.integer(nrow(A_start)),
                    as.integer(ncol(A_start)),
                    as.integer(nrow(B_start)),
                    as.integer(ncol(B_start)),
                    as.integer(sdir),  
                    as.integer(offset),
                    as.integer(Offsrow),
                    as.integer(Offscol),
                    as.integer(ncol(lambda.constraints)),
                    as.integer(nrow(lambda.constraints)),
                    as.integer(nrow(sample)),
                    as.integer(ncol(sample))
					   ))
    cat("VERSION MIT CONSTRAINTS WIRD VERWENDET", fill=TRUE)
    posterior <- .C("MCMCm",
                    sampledata = as.double(sample),											### Sample-Parameter
										Y	=	as.double(Y),																		###	Datenmatrizen	Y, Omega,	X
										Yrow = as.integer(nrow(Y)),
										p = as.integer(ncol(Y)),	## Ycol
										Omega	=	as.double(Omega),
#										Omegarow = as.integer(nrow(Omega)),
										Omegacol = as.integer(ncol(Omega)),
										Offs = as.double(Offs),
										X	=	as.double(X),
#										Xrow = as.integer(nrow(X)),
										Xcol = as.integer(ncol(X)),
										MCMCparameter = as.integer(MCMCparameter), 					### MCMC-Parameter
										seedarray = as.integer(seed.array),
										beta_stern = as.double(beta_stern),									###	Prior-Information	der	3	Parameterblöcke
#										beta_sternrow	=	as.integer(nrow(beta_stern)),				###	beta,	gamma, phi
#										beta_sterncol	=	as.integer(ncol(beta_stern)),
										Beta_stern	=	as.double(Beta_stern),
#										Beta_sternrow = as.integer(nrow(Beta_stern)),
#										Beta_sterncol = as.integer(ncol(Beta_stern)),
										Xdir_par_total = as.double(Xdir_par_total),								### Alle Designmatrizen für nonparam. direkte Effekte
#										Xdir_par_totalrow	=	as.integer(nrow(Xdir_par_total)),
										Xdir_par_totalcol	=	as.integer(ncol(Xdir_par_total)),
										Kdir_par_total	=	as.double(Kdir_par_total),								### Alle Strafmatrizen für nonparam. direkte Effekte
										Kdir_par_totalrow = as.integer(nrow(Kdir_par_total)),
										Kdir_par_totalcol = as.integer(ncol(Kdir_par_total)),
										nonpardir_par_total	=	as.double(nonpardir_par_total),		### Alle Strafmatrizen für nonparam. direkte Effekte
										nonpardir_par_totalrow = as.integer(nrow(nonpardir_par_total)),
										nonpardir_par_totalcol = as.integer(ncol(nonpardir_par_total)),
										X_par_total = as.double(X_par_total),								### Alle Designmatrizen für nonparam. Effekte
#										X_par_totalrow	=	as.integer(nrow(X_par_total)),
										X_par_totalcol	=	as.integer(ncol(X_par_total)),
										K_par_total	=	as.double(K_par_total),								### Alle Strafmatrizen für nonparam. Effekte
										K_par_totalrow = as.integer(nrow(K_par_total)),
										K_par_totalcol = as.integer(ncol(K_par_total)),
										nonpar_par_total	=	as.double(nonpar_par_total),		### Alle Strafmatrizen für nonparam. Effekte
										nonpar_par_totalrow = as.integer(nrow(nonpar_par_total)),
										nonpar_par_totalcol = as.integer(ncol(nonpar_par_total)),	
								   	v_stern	= as.double(v_stern),												
										v_sternrow = as.integer(nrow(v_stern)),
										v_sterncol = as.integer(ncol(v_stern)),
										s_stern = as.double(s_stern),
										s_sternrow = as.integer(nrow(s_stern)),
										s_sterncol = as.integer(ncol(s_stern)),
										z_start = as.double(z_start),												### Startwerte für Parametervektoren 
#										z_startrow = as.integer(nrow(z_start)),							### z, eta, gamma, beta, phi, nu
                    z_spa = as.integer(spa),
										eta_start = as.double(eta_start),
#										eta_startrow = as.integer(nrow(eta_start)),
										eta_startcol = as.integer(ncol(eta_start)),
										gamma_start = as.double(gamma_start),
										gamma_startrow = as.integer(nrow(gamma_start)),
										gamma_startcol = as.integer(ncol(gamma_start)),
										beta_start = as.double(beta_start),
										beta_startrow = as.integer(nrow(beta_start)),
										beta_startcol = as.integer(ncol(beta_start)),
										phi_start = as.double(phi_start),
										phi_startrow = as.integer(nrow(phi_start)),
										phi_startcol = as.integer(ncol(phi_start)),
										nu_start = as.double(nu_start),
										nu_startrow = as.integer(nrow(nu_start)),
										nu_startcol = as.integer(ncol(nu_start)),
										b0_start = as.double(b0_start),											### Zusätzlich Startwerte der einzelnen
										b0_startrow = as.integer(nrow(b0_start)),						### Matrizen b0, A, B, Gamma 
										b0_startcol = as.integer(ncol(b0_start)),
										A_start = as.double(A_start),
#										A_startrow = as.integer(nrow(A_start)),
#										A_startcol = as.integer(ncol(A_start)),
										B_start = as.double(B_start),
#										B_startrow = as.integer(nrow(B_start)),
#										B_startcol = as.integer(ncol(B_start)),
										ncat = as.integer(ncat),														### Anzahl Kategorien je Responsevariable
										tune = as.double(tune),
								##  OtherConstants = as.integer(OtherConstants),
                    accepts = as.integer(0),
                    DICav = as.double(0),
                    DICestimate = as.double(0),
                    gamma_stern = as.double(gamma_stern),
                    Gamma_stern = as.double(Gamma_stern),
                    lambda_constraints = as.double(lambda.constraints),
										PACKAGE="MCMCm"
										)
		cat("C++ Aufruf beendet !\n")
		# return(0)
		
    # put together matrix and build MCMC object to return
    sample <- matrix(posterior$sampledata, MCMCparameter[41], MCMCparameter[42],
                     byrow=TRUE)
    output <- mcmc(data=sample,start=1, end=mcmc, thin=thin)

    ####### Create parameter names
    par.names <- NULL
    #eta names
    eta.names <- NULL
    if(store.scores == TRUE) {
    	eta.names <- paste("eta", 1:Yrow, rep(1:factors, each=Yrow), sep=".")
    }
    par.names <- c(par.names, eta.names)
    #beta names
    beta.names <- NULL
    beta.names <- c(beta.names, paste("b", 1:p, "0", sep="."))
 		#beta.names <- c(beta.names, paste("b", Yvars[1:p], "0", sep="."))
    for (j in 1:p) {
    	if (Omega_exists)
    		beta.names <- c(beta.names, paste("a", j, 1:d, sep="."))
    		#beta.names <- c(beta.names, paste("a", Yvars[j], rep(Omegavars,1), sep="."))
    	beta.names <- c(beta.names, paste("b", j, 1:factors, sep="."))
	   	#beta.names <- c(beta.names, paste("b", Yvars[j], 1:factors, sep="."))
    }
    par.names <- c(par.names, beta.names)
    # ALS TEST KOMMT DER GAMMA INTERCEPT HIER REIN !!!!!!!!
    if (X_exists || X_par_exists) par.names <- c(par.names, "gamma.intercept")
    #gamma names
    gamma.names <- NULL
    if (X_exists) {
    	for (lat in 1:factors) {
    		gamma.names <- c(gamma.names, paste("gamma_par", lat, 1:q, sep="."))
    		#gamma.names <- c(gamma.names, paste("gamma", lat, rep(Xvars,1), sep="."))
    	}
    }
    par.names <- c(par.names, gamma.names)
    #gamma_par names
    gamma_par.names <- NULL
    if (X_par_exists) {
    	for (lat in 1:factors) {    
    		for (t in 1:s) {
    			gamma_par.names <- c(gamma_par.names, paste("gamma_f-", lat, "-", t,"-", covar.names_nonpar[t], "|",
    																									x_sorted[[t]][1:p_nr[t]], sep=""))
    		}
    	}	
    	par.names <- c(par.names, gamma_par.names)
    }
    
    #tau_par names
    tau_par.names <- NULL
    if (X_par_exists) {
    	for (lat in 1:factors)
    		tau_par.names <- c(tau_par.names, paste("tau_par", lat, 1:s, sep="."))
    	par.names <- c(par.names, tau_par.names)
    }
#################################################################################
# gammadir_par names    
    gammadir_par.names <- NULL
    if (Xdir_par_exists) {
    	for (lat in 1:p) {    
    		for (t in 1:sdir) {
    			gammadir_par.names <- c(gammadir_par.names, paste("gammadir_f-", lat, "-", t,"-", fixed.names_nonpar[t], "|",
    																									xdir_sorted[[t]][1:pdir_nr[t]], sep=""))
    		}
    	}	
    	par.names <- c(par.names, gammadir_par.names)
    }
    
    #tau_par names
    taudir_par.names <- NULL
    if (Xdir_par_exists) {
    	for (lat in 1:p)
    		taudir_par.names <- c(taudir_par.names, paste("taudir_par", lat, 1:sdir, sep="."))
    	par.names <- c(par.names, taudir_par.names)
    }


#################################################################################

    ## theta.names
    theta.names <- NULL
    #if (p3>0) {
    	theta.names <- paste("theta", 1:p, sep=".")
    	par.names <- c(par.names, theta.names)
    #}
    ## nu.names
    nu.names <- NULL
    if ((p1 > 0) && (ncol(nu_start) > 3)){
  		for (j in 1:p)
  			if (ncat[j] > 2) nu.names <- c(nu.names, paste("nu", j, 2:(ncat[j]-1),sep=".")) 
  	}
  	par.names <- c(par.names, nu.names)

    z.names <-NULL
    if(store.z == TRUE){
	     z.names <- paste("z",1:Yrow, rep(1:(spaanz-p3),each=Yrow), sep=".")
    }
    par.names <- c(par.names, z.names)
    ## DIC name
    if(DIC == TRUE) {
    	par.names <- c(par.names, "DIC")
    }

    varnames(output) <- par.names
    
    # add attributes info
    #attr(output, "constraints") <- lambda.constraints
    attr(output, "bspline.basis") <- bspline
    attr(output, "n.factors") <- factors
    attr(output, "n.indicators") <- p
    attr(output, "n.ordinal.indicators") <- p1
    attr(output, "n.poisson.indicators") <- p2
    attr(output, "n.metric.indicators") <- p3
		attr(output, "vars.indicators") <- Yvars
    attr(output, "n.fixed.effects") <- d
    attr(output, "vars.fixed.effects") <- Omegavars
    attr(output, "n.par_covar") <- q
    attr(output, "vars.par_covar") <- Xvars
    attr(output, "n.nonpar_covar") <- s
    attr(output, "n.nonpar_direct_covar") <- sdir
    if (X_par_exists)
    	attr(output, "vars.nonpar_covar") <- covar.names_nonpar
    else 
    	attr(output, "vars.nonpar_covar") <- NULL
    if (X_par_exists)
    	attr(output, "nonpar_covar.types") <- covar.types_nonpar
    else 
    	attr(output, "nonpar_covar.types") <- NULL
    attr(output, "n.degree.psplines") <- degr
    ###########################################################################
    attr(output, "bsplinedir.basis") <- bsplinedir
    if (Xdir_par_exists)
    	attr(output, "vars.nonpar_direct_covar") <- fixed.names_nonpar
    else 
    	attr(output, "vars.nonpar_direct_covar") <- NULL
    if (Xdir_par_exists)
    	attr(output, "nonpar_direct_covar.types") <- fixed.types_nonpar
    else 
    	attr(output, "nonpar_direct_covar.types") <- NULL
    attr(output, "n.degreedir.psplines") <- degrdir
    ###########################################################################
    
    attr(output, "n.cat") <- ncat
    attr(output, "gamma_stern") <- gamma_stern_call    
    attr(output, "Gamma_stern") <- Gamma_stern_call
    attr(output, "beta_stern") <- beta_stern_call
    attr(output, "Beta_stern") <- Beta_stern_call
    attr(output, "beta.start") <- beta.start
    attr(output, "gamma.start") <- gamma.start
    if (mh) {
    	attr(output, "accept.rate") <- posterior$accepts  / ((burnin+mcmc)*p1)
    	attr(output, "mh.tune") <- tune_original
    	}
    else attr(output, "accept.rate") <- "MH-algorithm for cutpoints not used."
    if (gm) attr(output, "GM-step") <- "GM-step used."
    	else attr(output, "GM-step") <- "GM-step not used."
    attr(output,"title") <-
        "MCMCpack MIMIC Posterior Density Sample"
    attr(output, "function.call") <- call
    
    if (DIC == TRUE) {
    	DICF <- data.frame(c("DIC_averaged", "DIC_estimate", "free parameters", "DIC"), 
    										 c(posterior$DICav, posterior$DICestimate, 
    										 	 posterior$DICav-posterior$DICestimate, 2*posterior$DICav-posterior$DICestimate ))
    	colnames(DICF) <- c("Measure", "Value")
    	attr(output, "DIC") <- DICF
    }
    
    return(output)
  }

tpower <- function(x, t, pow) {
	(x - t) ^ pow * (x > t)
}

bbase <- function(x, x_min, x_max, intervals, degr) {
	# Construct a B-spline basis of degree 'degr'
	dx <- (x_max-x_min) / intervals
	knots <- seq(x_min - degr * dx, x_max + degr * dx, by = dx)
	P <- outer(x, knots, tpower, degr)
	n <- dim(P)[2]
	D <- diff(diag(n), diff = degr + 1) / (gamma(degr + 1) * dx ^ degr)
	B <- (-1) ^(degr+1) * P %*% t(D)
	return(B)
}

create_K_r1 <- function(knots) {
	NR <- length(knots)
	K_temp <- matrix(0, nrow=NR, ncol=NR)
	diff <- knots[2:NR] - knots[1:(NR-1)]
	for (i in 1:NR) {
		if (i==1) K_temp[i,i] <- diff[1]^(-1)
		else {
			if (i==NR) K_temp[i,i] <- diff[i-1]^(-1)
			 	else K_temp[i,i] <- diff[i-1]^(-1) + diff[i]^(-1)
			K_temp[i,i-1] <- K_temp[i-1,i] <- -diff[i-1]^(-1)
		} 
	}
	return(K_temp)
}

create_K_r2 <- function(knots) {
	NR <- length(knots)
	K_temp <- matrix(0, nrow=NR, ncol=NR)
	diff <- knots[2:NR] - knots[1:(NR-1)]			
	diff <- c(0, diff)
	weights <- NULL
	for (i in 3:NR) {   
		weights[i] <- diff[i]	# Alternativ: weights[i] <- diff[i]*(1+diff[i]/diff[i-1])  
	}
	K_temp[1,1] <- (diff[3]/diff[2])^2 / weights[3]
	K_temp[1,2] <- K_temp[2,1] <- -(1+diff[3]/diff[2])*diff[3]/(diff[2]*weights[3])
	K_temp[1,3] <- diff[3]/(diff[2]*weights[3])
	K_temp[2,2] <- (1+diff[3]/diff[2])^2/weights[3] + (diff[4]/diff[3])^2 / weights[4]
	K_temp[2,3] <- -(1+diff[3]/diff[2])/weights[3] - (1+diff[4]/diff[3])*diff[4]/(diff[3]*weights[4])
	K_temp[2,4] <- diff[4]/(diff[3]*weights[4])
	for (i in 3:NR) {
		K_temp[i,i-2] <- diff[i]/(diff[i-1]*weights[i])
		if (i < NR) 
			K_temp[i,i-1] <- -(1+diff[i]/diff[i-1])/weights[i] - (1+diff[i+1]/diff[i])*diff[i+1]/(diff[i]*weights[i+1])
		else																													
			K_temp[i,i-1] <- -(1+diff[i]/diff[i-1])/weights[i]
		if (i < NR-1)
			K_temp[i,i] 	<- 1/weights[i] + (1+diff[i+1]/diff[i])^2/weights[i+1] + (diff[i+2]/diff[i+1])^2/weights[i+2]
		else if (i < NR)
			K_temp[i,i] 	<- 1/weights[i] + (1+diff[i+1]/diff[i])^2/weights[i+1]
		else
			K_temp[i,i] 	<- 1/weights[i]
		if (i < NR-1)
			K_temp[i,i+1] <- -(1+diff[i+1]/diff[i])/weights[i+1] -(1+diff[i+2]/diff[i+1])*diff[i+2]/(diff[i+1]*weights[i+2])
		else if (i < NR)
			K_temp[i,i+1] <- -(1+diff[i+1]/diff[i])/weights[i+1]
		if (i < NR-1)
			K_temp[i,i+2] <- diff[i+2]/(diff[i+1]*weights[i+2])
	}
	return(K_temp)
}
