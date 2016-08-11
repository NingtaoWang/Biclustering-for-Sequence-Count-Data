REPEAT_LIMIT <<- 10001
RELATIVE_DIFF <<- 0.0000001
Method_H0 <<- 1

library(Rcpp)
library(RcppArmadillo)

cppFunction('
Rcpp::NumericVector MyOuterNB( 	Rcpp::NumericMatrix DataR,
					Rcpp::NumericMatrix rR,
					Rcpp::NumericMatrix probR) {
	int n = DataR.nrow(),
	m = DataR.ncol(),
	K = probR.nrow(),
	L = probR.ncol();

	int n_iter = n * m * K * L, 
	n_iter_ijk = n * m * K, 
	n_iter_ij = n * m, 
	n_iter_i = n;
	Rcpp::NumericVector Output(n_iter);
	for(int h = 0; h < n_iter; h++){
		int l = floor(h / n_iter_ijk);
		int ijk = h - l * n_iter_ijk;
		int k = floor(ijk / n_iter_ij);
		int ij = ijk - k * n_iter_ij;
		int j = floor(ij / n_iter_i);
		int i = ij - j * n_iter_i;

		double tmp = lgamma(DataR(i,j) + rR(k,l)) - lgamma(rR(k,l)) - lgamma(DataR(i,j) + 1) + rR(k,l) * log(probR(k,l)) + DataR(i,j) * log(1 - probR(k,l));
		Output[h] = exp(tmp);
	}
	return Output;      
}')

cppFunction('
Rcpp::NumericVector MyOuterPsiPsi( 	Rcpp::NumericMatrix Psiik,
						Rcpp::NumericVector Psiijkl) {
	Rcpp::IntegerVector DimsPsiijkl = Psiijkl.attr("dim");
	
	int n_iter = DimsPsiijkl[0] * DimsPsiijkl[1] * DimsPsiijkl[2] * DimsPsiijkl[3], 
	n_iter_ijk = DimsPsiijkl[0] * DimsPsiijkl[1] * DimsPsiijkl[2], 
	n_iter_ij = DimsPsiijkl[0] * DimsPsiijkl[1], 
	n_iter_i = DimsPsiijkl[0];
	Rcpp::NumericVector Output(n_iter);
	for(int h = 0; h < n_iter; h++){
		int l = floor(h / n_iter_ijk);
		int ijk = h - l * n_iter_ijk;
		int k = floor(ijk / n_iter_ij);
		int ij = ijk - k * n_iter_ij;
		int j = floor(ij / n_iter_i);
		int i = ij - j * n_iter_i;

		Output[h] = Psiik(i, k) * Psiijkl(h);
	}
	return Output;      
}')

cppFunction('
Rcpp::NumericVector MyOuterOmegaOmega( 	Rcpp::NumericMatrix Omegajl,
							Rcpp::NumericVector Omegaijkl) {
	Rcpp::IntegerVector DimsOmegaijkl = Omegaijkl.attr("dim");
	
	int n_iter = DimsOmegaijkl[0] * DimsOmegaijkl[1] * DimsOmegaijkl[2] * DimsOmegaijkl[3], 
	n_iter_ijk = DimsOmegaijkl[0] * DimsOmegaijkl[1] * DimsOmegaijkl[2], 
	n_iter_ij = DimsOmegaijkl[0] * DimsOmegaijkl[1], 
	n_iter_i = DimsOmegaijkl[0];
	Rcpp::NumericVector Output(n_iter);
	for(int h = 0; h < n_iter; h++){
		int l = floor(h / n_iter_ijk);
		int ijk = h - l * n_iter_ijk;
		int k = floor(ijk / n_iter_ij);
		int ij = ijk - k * n_iter_ij;
		int j = floor(ij / n_iter_i);
		int i = ij - j * n_iter_i;

		Output[h] = Omegajl(j, l) * Omegaijkl(h);
	}
	return Output;      
}')

cppFunction('
Rcpp::NumericVector MyOuterMu( 	Rcpp::NumericMatrix DataR,
						Rcpp::NumericVector PsiPsi,
						Rcpp::NumericVector OmegaOmega) {
	Rcpp::IntegerVector DimsPsiPsi = PsiPsi.attr("dim");
	
	int n_iter = DimsPsiPsi[0] * DimsPsiPsi[1] * DimsPsiPsi[2] * DimsPsiPsi[3], 
	n_iter_ijk = DimsPsiPsi[0] * DimsPsiPsi[1] * DimsPsiPsi[2], 
	n_iter_ij = DimsPsiPsi[0] * DimsPsiPsi[1], 
	n_iter_i = DimsPsiPsi[0];
	Rcpp::NumericVector Output(n_iter);
	for(int h = 0; h < n_iter; h++){
		int l = floor(h / n_iter_ijk);
		int ijk = h - l * n_iter_ijk;
		int k = floor(ijk / n_iter_ij);
		int ij = ijk - k * n_iter_ij;
		int j = floor(ij / n_iter_i);
		int i = ij - j * n_iter_i;

		Output(h) = DataR(i,j) * (PsiPsi(h) + OmegaOmega(h));
	}
	return Output;      
}')

cppFunction('
Rcpp::NumericVector MyOuterPartial1( 	Rcpp::NumericMatrix DataR,
					Rcpp::NumericVector PsiPsi,
					Rcpp::NumericVector OmegaOmega,					
					Rcpp::NumericMatrix muR,
					Rcpp::NumericMatrix rR) {
	Rcpp::IntegerVector DimsPsiPsi = PsiPsi.attr("dim");
	
	int n_iter = DimsPsiPsi[0] * DimsPsiPsi[1] * DimsPsiPsi[2] * DimsPsiPsi[3], 
	n_iter_ijk = DimsPsiPsi[0] * DimsPsiPsi[1] * DimsPsiPsi[2], 
	n_iter_ij = DimsPsiPsi[0] * DimsPsiPsi[1], 
	n_iter_i = DimsPsiPsi[0];
	Rcpp::NumericVector Output(n_iter);
	for(int h = 0; h < n_iter; h++){
		int l = floor(h / n_iter_ijk);
		int ijk = h - l * n_iter_ijk;
		int k = floor(ijk / n_iter_ij);
		int ij = ijk - k * n_iter_ij;
		int j = floor(ij / n_iter_i);
		int i = ij - j * n_iter_i;

		double tmp = R::digamma(DataR(i,j) + rR(k,l)) - R::digamma(rR(k,l)) - (DataR(i,j) + rR(k,l))/(muR(k,l) + rR(k,l)) + 1 + log(rR(k,l)) - log(muR(k,l) + rR(k,l));
		Output(h) = tmp * (PsiPsi(h) + OmegaOmega(h));
	}
	return Output;      
}')

cppFunction('
Rcpp::NumericVector MyOuterPartial2( 	Rcpp::NumericMatrix DataR,
					Rcpp::NumericVector PsiPsi,
					Rcpp::NumericVector OmegaOmega,					
					Rcpp::NumericMatrix muR,
					Rcpp::NumericMatrix rR) {
	Rcpp::IntegerVector DimsPsiPsi = PsiPsi.attr("dim");
	
	int n_iter = DimsPsiPsi[0] * DimsPsiPsi[1] * DimsPsiPsi[2] * DimsPsiPsi[3], 
	n_iter_ijk = DimsPsiPsi[0] * DimsPsiPsi[1] * DimsPsiPsi[2], 
	n_iter_ij = DimsPsiPsi[0] * DimsPsiPsi[1], 
	n_iter_i = DimsPsiPsi[0];
	Rcpp::NumericVector Output(n_iter);
	for(int h = 0; h < n_iter; h++){
		int l = floor(h / n_iter_ijk);
		int ijk = h - l * n_iter_ijk;
		int k = floor(ijk / n_iter_ij);
		int ij = ijk - k * n_iter_ij;
		int j = floor(ij / n_iter_i);
		int i = ij - j * n_iter_i;

		double MuRR = muR(k,l) + rR(k,l);
		double tmp = R::trigamma(DataR(i,j) + rR(k,l)) - R::trigamma(rR(k,l)) - 2/MuRR + (DataR(i,j) + rR(k,l))/MuRR/MuRR + 1/rR(k,l);
		Output(h) = tmp * (PsiPsi(h) + OmegaOmega(h));
	}
	return Output;      
}')

cppFunction('
Rcpp::NumericVector MyOuterDMu( 	Rcpp::NumericMatrix DataR,
						Rcpp::NumericMatrix MuR,
						Rcpp::NumericMatrix rR,
						Rcpp::NumericVector Fijkl) {
	Rcpp::IntegerVector DimsFijkl = Fijkl.attr("dim");
	
	int n_iter = DimsFijkl[0] * DimsFijkl[1] * DimsFijkl[2] * DimsFijkl[3], 
	n_iter_ijk = DimsFijkl[0] * DimsFijkl[1] * DimsFijkl[2], 
	n_iter_ij = DimsFijkl[0] * DimsFijkl[1], 
	n_iter_i = DimsFijkl[0];
	Rcpp::NumericVector Output(n_iter);
	for(int h = 0; h < n_iter; h++){
		int l = floor(h / n_iter_ijk);
		int ijk = h - l * n_iter_ijk;
		int k = floor(ijk / n_iter_ij);
		int ij = ijk - k * n_iter_ij;
		int j = floor(ij / n_iter_i);
		int i = ij - j * n_iter_i;

		double DataMu = DataR(i,j) - MuR(k,l);
		Output(h) = Fijkl(h) * rR(k,l) * DataMu / MuR(k,l) / (rR(k,l) + MuR(k,l));
	}
	return Output;      
}')

cppFunction('
Rcpp::NumericVector MyOuterDSigma( 	Rcpp::NumericMatrix DataR,
						Rcpp::NumericMatrix MuR,
						Rcpp::NumericMatrix rR,
						Rcpp::NumericVector Fijkl) {
	Rcpp::IntegerVector DimsFijkl = Fijkl.attr("dim");
	
	int n_iter = DimsFijkl[0] * DimsFijkl[1] * DimsFijkl[2] * DimsFijkl[3], 
	n_iter_ijk = DimsFijkl[0] * DimsFijkl[1] * DimsFijkl[2], 
	n_iter_ij = DimsFijkl[0] * DimsFijkl[1], 
	n_iter_i = DimsFijkl[0];
	Rcpp::NumericVector Output(n_iter);
	for(int h = 0; h < n_iter; h++){
		int l = floor(h / n_iter_ijk);
		int ijk = h - l * n_iter_ijk;
		int k = floor(ijk / n_iter_ij);
		int ij = ijk - k * n_iter_ij;
		int j = floor(ij / n_iter_i);
		int i = ij - j * n_iter_i;

		double tmp = R::digamma(DataR(i,j) + rR(k,l)) - R::digamma(rR(k,l)) - (DataR(i,j) + rR(k,l))/(MuR(k,l) + rR(k,l)) + 1 + log(rR(k,l)) - log(MuR(k,l) + rR(k,l));
		Output(h) = Fijkl(h) * tmp;
	}
	return Output;      
}')

BM.StepE <- function(data, mu, r, p, q, K, L, n, m)
{
	prob <- r / (r + mu)
	Fijkl <- MyOuterNB(data, r, prob)
	Fijkl <- array(Fijkl, dim = c(n, m, K, L))
	Fijkl[is.nan(Fijkl)] <- 1e-300
	Fijkl[Fijkl < 1e-300] <- 1e-300

	Psi.ijkl <- apply(Fijkl,1:3,FUN=function(x){x*q/(x%*%q)})
	Psi.ijkl <- aperm(Psi.ijkl, c(2,3,4,1))
	F1.ijk <- apply(Fijkl,1:3,FUN=function(x){x%*%q})
	F1.log.ijk <- log(F1.ijk)
	F1.log.ik <- apply(F1.log.ijk,c(1,3),sum)
	F1.log.ik.adjtmp <- apply(F1.log.ik,1,
		FUN=function(x){
			tmpx <- max(x)-1400
			x[which(x<tmpx)]<-tmpx
			return(x)})	
	F1.log.ik.adj <- apply(F1.log.ik.adjtmp,2,FUN=function(x){x-max(x)/2-min(x)/2})
	F1.log.ik.adj <- t(F1.log.ik.adj)
	F1.ik <- apply(F1.log.ik.adj,1,FUN=function(x){exp(x)*p/(exp(x)%*%p)})
	Psi.ik <- t(F1.ik)

	Omega.ijkl <- apply(Fijkl,c(1,2,4),FUN=function(x){x*p/(x%*%p)})
	Omega.ijkl <- aperm(Omega.ijkl, c(2,3,1,4))
	F2.ijl <- apply(Fijkl,c(1,2,4),FUN=function(x){x%*%p})
	F2.log.ijl <- log(F2.ijl)
	F2.log.jl <- apply(F2.log.ijl,c(2,3),sum)
	F2.log.jl.adjtmp <- apply(F2.log.jl,1,
		FUN=function(x){
			tmpx <- max(x)-1400
			x[which(x<tmpx)]<-tmpx
			return(x)})	
	F2.log.jl.adj <- apply(F2.log.jl.adjtmp,2,FUN=function(x){x-max(x)/2-min(x)/2})
	F2.log.jl.adj <- t(F2.log.jl.adj)
	F2.jl <- apply(F2.log.jl.adj,1,FUN=function(x){exp(x)*q/(exp(x)%*%q)})
	Omega.jl <- t(F2.jl)

	F1.log.i.mean <- apply(F1.log.ik.adjtmp,2,FUN=function(x){max(x)/2+min(x)/2})
	F1.ik.tmp <- apply(F1.log.ik.adj,1,FUN=function(x){exp(x)*p})
	F1.i.tmp <- apply(F1.ik.tmp,2,sum)
	F1.i.log <- log(F1.i.tmp)
	F2.log.j.mean <- apply(F2.log.jl.adjtmp,2,FUN=function(x){max(x)/2+min(x)/2})
	F2.jl.tmp <- apply(F2.log.jl.adj,1,FUN=function(x){exp(x)*q})
	F2.j.tmp <- apply(F2.jl.tmp,2,sum)
	F2.j.log <- log(F2.j.tmp)
	LogLikelihood <- sum(F1.i.log)+sum(F2.j.log)+sum(F1.log.i.mean)+sum(F2.log.j.mean)

	return(list(Psi.ik = Psi.ik, Psi.ijkl = Psi.ijkl, 
			Omega.jl = Omega.jl, Omega.ijkl = Omega.ijkl, LogLikelihood = LogLikelihood))
}

BM.StepM <- function(data, Psi.ik, Psi.ijkl, Omega.jl, Omega.ijkl, mu, r, K, L, n, m, method)
{
	PsiPsi <- MyOuterPsiPsi(Psi.ik, Psi.ijkl)
	PsiPsi <- array(PsiPsi, dim = c(n, m, K, L))
	PsiPsi.jl <- apply(PsiPsi,c(2,4),sum)

	OmegaOmega <- MyOuterOmegaOmega(Omega.jl, Omega.ijkl)
	OmegaOmega <- array(OmegaOmega, dim = c(n, m, K, L))
	OmegaOmega.ik <- apply(OmegaOmega,c(1,3),sum)

	PsiOmega.kl <- apply(OmegaOmega+PsiPsi, 3:4, sum)

	p <- apply(Psi.ik+OmegaOmega.ik,2,sum)/n/(m+1)
	q <- apply(Omega.jl+PsiPsi.jl,2,sum)/m/(n+1)

	PsiOmegaData <- MyOuterMu(data, PsiPsi, OmegaOmega)
	PsiOmegaData <- array(PsiOmegaData, dim=c(n, m, K, L))
	PsiOmegaData.kl <- apply(PsiOmegaData, 3:4, sum)
	new.mu <- PsiOmegaData.kl/PsiOmega.kl

	PsiOmegaPartial1 <- MyOuterPartial1(data, PsiPsi, OmegaOmega, mu, r)
	PsiOmegaPartial1 <- array(PsiOmegaPartial1, dim=c(n, m, K, L))
	PsiOmegaPartial1.kl <- apply(PsiOmegaPartial1, 3:4, sum)

	PsiOmegaPartial2 <- MyOuterPartial2(data, PsiPsi, OmegaOmega, mu, r)
	PsiOmegaPartial2 <- array(PsiOmegaPartial2, dim=c(n, m, K, L))
	PsiOmegaPartial2.kl <- apply(PsiOmegaPartial2, 3:4, sum)
	new.r <- r - PsiOmegaPartial1.kl / PsiOmegaPartial2.kl
	#new.r[new.r < 1] <- 1
	#new.r[new.r < 1] <- r[new.r < 1]
	new.r[is.na(new.r)] <- r[is.na(new.r)]
	new.r[new.r > 1e+10] <- 1e+10
	new.r[new.r <= 0] <- r[new.r <= 0]


	return(list(p=p, q=q, mu=new.mu, r=new.r) )
}

BM.InitTheta <- function(data, K, L)
{
	res1 <- kmeans(data,K)
	res2 <- kmeans(t(data),L)
	
	p <- res1$size/sum(res1$size)
	q <- res2$size/sum(res2$size)	
	
	mu <- array(0,dim=c(K,L))
	r <- array(0,dim=c(K,L))
	for(k in 1:K)
	{
		for(l in 1:L)
		{
			tmp <- data[which(res1$cluster==k),which(res2$cluster==l)]
			tmp <- as.vector(tmp)
			mu[k,l] <- mean(tmp)

			if(is.na(var(tmp))){
				r[k, l] <- 1
			}
			else{
				if(var(tmp) > mean(tmp)){
					r[k, l] <- mu[k,l]^2 / {var(tmp) - mu[k,l]}
				}
				else{
					r[k, l] <- max(tmp)
				}
			}
			
			r[r<=0] <- 1
		}
	}

	return(list(mu=mu, r = r, p=p, q=q, row.class=res1$cluster, col.class=res2$cluster) )
}

BM.RunEM <- function(data, K, L, method=Method_H0)
{
	n <- dim(data)[1]
	m <- dim(data)[2]

	#theta <- BM.InitTheta(data, K, L)
	theta <- result0
	mu <- theta$mu
	r <- theta$r
	p <- theta$p
	q <- theta$q

	rpt <- 1
	LogLikelihood <- -Inf
	while(TRUE){
		OldLogLikelihood <- LogLikelihood

		EResult <- BM.StepE(data, mu, r, p, q, K, L, n, m)

		Psi.ik <- EResult$Psi.ik
		Psi.ijkl <- EResult$Psi.ijkl
		Omega.jl <- EResult$Omega.jl
		Omega.ijkl <- EResult$Omega.ijkl
		LogLikelihood <- EResult$LogLikelihood
	
		cat("K:",K,"  L:",L,"	rpt:", rpt, "\n")
		cat("mu:", mu, "\n")
		cat("r:", r, "\n")
		cat("p:", p, "\n")
		cat("q:", q, "\n")
		cat("LogLikelihood:",LogLikelihood,"\n")
		cat("\n")
		
		if (is.infinite(LogLikelihood))
			break
		if( (abs(1 - OldLogLikelihood/LogLikelihood) < RELATIVE_DIFF) ){
			cat("quit due to likelihood\n")
			cat("RELATIVE_DIFF:",RELATIVE_DIFF,"\n")
			break
		}
		if( rpt >= REPEAT_LIMIT ){
			cat("quit due to rpt\n")
			break
		}
		if( OldLogLikelihood > LogLikelihood ){
			cat("quit due to likelihood decreasing\n")
			#break
		}
		oldtheta <- theta

		theta <- BM.StepM(data, Psi.ik, Psi.ijkl, Omega.jl, Omega.ijkl, mu, r, K, L, n, m, method)
		mu <- theta$mu
		r <- theta$r
		p <- theta$p
		q <- theta$q

		if( sum(is.na(mu)) > 0 ){
			cat("mu na\n")
			tmp.old.mu <- as.vector(oldtheta$mu)
			tmp.mu <- as.vector(mu)
			ind.na <- which(is.na(mu) == T)
			tmp.mu[ind.na] <- tmp.old.mu[ind.na]
			mu <- array(tmp.mu, dim = c(K, L))
			theta$mu <- mu
		}


		#identifiable
		#tmp.ind1 <- sort(colSums(mu),index.return=T)$ix
		#mu <- mu[,tmp.ind1]
		#r <- r[,tmp.ind1]
		#q <- q[tmp.ind1]

		#tmp.ind2 <- sort(rowSums(mu),index.return=T)$ix
		#mu <- mu[tmp.ind2,]
		#r <- r[tmp.ind2,]
		#p <- p[tmp.ind2]

		rpt <- rpt + 1
	}

	row.class <- max.col(Psi.ik)
	col.class <- max.col(Omega.jl)

	return(list(p=p, q=q, row.class = row.class, col.class = col.class, 
				LogLikelihood = LogLikelihood, mu = mu, r=r) )
}

BM.BICsimu <- function(data, mu, r, p, q, K, L)
{
	n <- dim(data)[1]
	m <- dim(data)[2]
	#K <- dim(mu)[1]
	#L <- dim(mu)[2]
	
	prob <- r / (r + mu)
	Fijkl <- MyOuterNB(data, r, prob)
	Fijkl <- array(Fijkl, dim = c(n, m, K, L))

	Fijkl.dmu <- MyOuterDMu(data, mu, r, Fijkl)
	Fijkl.dmu <- array(Fijkl.dmu, dim = c(n, m, K, L))

	Fijkl.dr <- MyOuterDSigma(data, mu, r, Fijkl)
	Fijkl.dr <- array(Fijkl.dr, dim = c(n, m, K, L))

	#Fijkl[is.nan(Fijkl)] <- 1e-300
	Fijkl[Fijkl < 1e-300] <- 1e-300
	#Fijkl.dmu[is.nan(Fijkl.dmu)] <- 1e-300
	#Fijkl.dr[is.nan(Fijkl.dr)] <- 1e-300

	F1.ijk <- apply(Fijkl,1:3,FUN=function(x){x%*%q})
	F1.log.ijk <- log(F1.ijk)
	F1.log.ik <- apply(F1.log.ijk,c(1,3),sum)
	F1.log.ik.adjtmp <- apply(F1.log.ik,1,
		FUN=function(x){
			tmpx <- max(x)-1400
			x[which(x<tmpx)]<-tmpx
			return(x)})	
	F1.log.ik.adj <- apply(F1.log.ik.adjtmp,2,FUN=function(x){x-max(x)/2-min(x)/2})
	F1.log.ik.adj <- t(F1.log.ik.adj)
	F1.ik <- apply(F1.log.ik.adj,1,FUN=function(x){exp(x)/(exp(x)%*%p)})
	F1.ik <- t(F1.ik)
 	drow.dk <- apply(as.matrix(F1.ik[,-K]),2,FUN=function(x){x-F1.ik[,K]})


	F1.i <- apply(F1.log.ik.adj,1,FUN=function(x){exp(x)%*%p})

	F1.ijkl <- apply(Fijkl,1:3,FUN=function(x){x/(x%*%q)})
	F1.ijkl <- aperm(F1.ijkl, c(2,3,4,1))
	F1.ikl <- apply(F1.ijkl,c(1,3,4),sum)
	F1.ikl <- apply(F1.ikl,3,FUN=function(x){x*exp(F1.log.ik.adj)})
	F1.ikl <- array(F1.ikl,dim=c(n,K,L))
	F1.il <- apply(F1.ikl,c(1,3),FUN=function(x){x%*%p})
	F1.il <- apply(F1.il,2,FUN=function(x){x/F1.i})
	drow.dl <- apply(as.matrix(F1.il[,-L]),2,FUN=function(x){x-F1.il[,L]})

	
	F1.ijkl.dmu <- apply(Fijkl.dmu,c(1,2,3),FUN=function(x){x*q})
	F1.ijkl.dmu <- apply(F1.ijkl.dmu,1,FUN=function(x){x/F1.ijk})
	F1.ijkl.dmu <- array(F1.ijkl.dmu,dim=c(n,m,K,L))
	F1.ikl.dmu <- apply(F1.ijkl.dmu,c(1,3,4),sum)
	F1.ikl.dmu <- apply(F1.ikl.dmu,3,FUN=function(x){x*exp(F1.log.ik.adj)})
	F1.ikl.dmu <- array(F1.ikl.dmu,dim=c(n,K,L))
	F1.ikl.dmu <- apply(F1.ikl.dmu,c(1,3),FUN=function(x){x*p})
	F1.ikl.dmu <- apply(F1.ikl.dmu,c(1,3),FUN=function(x){x/F1.i})
	drow.dmu <- apply(F1.ikl.dmu,1,FUN=function(x){as.vector(x)})
	drow.dmu <- t(drow.dmu)

	F1.ijkl.dr <- apply(Fijkl.dr,c(1,2,3),FUN=function(x){x*q})
	F1.ijkl.dr <- apply(F1.ijkl.dr,1,FUN=function(x){x/F1.ijk})
	F1.ijkl.dr <- array(F1.ijkl.dr,dim=c(n,m,K,L))
	F1.ikl.dr <- apply(F1.ijkl.dr,c(1,3,4),sum)
	F1.ikl.dr <- apply(F1.ikl.dr,3,FUN=function(x){x*exp(F1.log.ik.adj)})
	F1.ikl.dr <- array(F1.ikl.dr,dim=c(n,K,L))
	F1.ikl.dr <- apply(F1.ikl.dr,c(1,3),FUN=function(x){x*p})
	F1.ikl.dr <- apply(F1.ikl.dr,c(1,3),FUN=function(x){x/F1.i})
	drow.dr <- apply(F1.ikl.dr,1,FUN=function(x){as.vector(x)})
	drow.dr <- t(drow.dr)



	F2.ijl <- apply(Fijkl,c(1,2,4),FUN=function(x){x%*%p})
	F2.log.ijl <- log(F2.ijl)
	F2.log.jl <- apply(F2.log.ijl,c(2,3),sum)
	F2.log.jl.adjtmp <- apply(F2.log.jl,1,
		FUN=function(x){
			tmpx <- max(x)-1400
			x[which(x<tmpx)]<-tmpx
			return(x)})	
	F2.log.jl.adj <- apply(F2.log.jl.adjtmp,2,FUN=function(x){x-max(x)/2-min(x)/2})
	F2.log.jl.adj <- t(F2.log.jl.adj)
	F2.jl <- apply(F2.log.jl.adj,1,FUN=function(x){exp(x)/(exp(x)%*%q)})
	F2.jl <- t(F2.jl)
	dcol.dl <- apply(as.matrix(F2.jl[,-L]),2,FUN=function(x){x-F2.jl[,L]})


	F2.j <- apply(F2.log.jl.adj,1,FUN=function(x){exp(x)%*%q})

	F2.ijkl <- apply(Fijkl,c(1,2,4),FUN=function(x){x/(x%*%p)})
	F2.ijkl <- aperm(F2.ijkl, c(2,3,1,4))
	F2.jkl <- apply(F2.ijkl,c(2,3,4),sum)
	F2.jkl <- apply(F2.jkl,2,FUN=function(x){x*exp(F2.log.jl.adj)})
	F2.jkl <- array(F2.jkl,dim=c(m,L,K))
	F2.jk <- apply(F2.jkl,c(1,3),FUN=function(x){x%*%q})
	F2.jk <- apply(F2.jk,2,FUN=function(x){x/F2.j})
	dcol.dk <- apply(as.matrix(F2.jk[,-K]),2,FUN=function(x){x-F2.jk[,K]})


	F2.ijkl.dmu <- apply(Fijkl.dmu,c(1,2,4),FUN=function(x){x*p})
	F2.ijkl.dmu <- apply(F2.ijkl.dmu,1,FUN=function(x){x/F2.ijl})
	F2.ijkl.dmu <- array(F2.ijkl.dmu,dim=c(n,m,L,K))
	F2.jkl.dmu <- apply(F2.ijkl.dmu,c(2,3,4),sum)
	F2.jkl.dmu <- apply(F2.jkl.dmu,3,FUN=function(x){x*exp(F2.log.jl.adj)})
	F2.jkl.dmu <- array(F2.jkl.dmu,dim=c(m,L,K))
	F2.jkl.dmu <- apply(F2.jkl.dmu,c(1,3),FUN=function(x){x*q})
	F2.jkl.dmu <- apply(F2.jkl.dmu,c(1,3),FUN=function(x){x/F2.j})
	F2.jkl.dmu <- aperm(F2.jkl.dmu,c(1,3,2))
	dcol.dmu <- apply(F2.jkl.dmu, 1, FUN = function(x){as.vector(x)})
	dcol.dmu <- t(dcol.dmu)

	F2.ijkl.dr <- apply(Fijkl.dr,c(1,2,4),FUN=function(x){x*p})
	F2.ijkl.dr <- apply(F2.ijkl.dr,1,FUN=function(x){x/F2.ijl})
	F2.ijkl.dr <- array(F2.ijkl.dr,dim=c(n,m,L,K))
	F2.jkl.dr <- apply(F2.ijkl.dr,c(2,3,4),sum)
	F2.jkl.dr <- apply(F2.jkl.dr,3,FUN=function(x){x*exp(F2.log.jl.adj)})
	F2.jkl.dr <- array(F2.jkl.dr,dim=c(m,L,K))
	F2.jkl.dr <- apply(F2.jkl.dr,c(1,3),FUN=function(x){x*q})
	F2.jkl.dr <- apply(F2.jkl.dr,c(1,3),FUN=function(x){x/F2.j})
	F2.jkl.dr <- aperm(F2.jkl.dr,c(1,3,2))
	dcol.dr <- apply(F2.jkl.dr, 1, FUN = function(x){as.vector(x)})
	dcol.dr <- t(dcol.dr)


	drow <- cbind(drow.dk,drow.dl,drow.dmu,drow.dr)
	dcol <- cbind(dcol.dk,dcol.dl,dcol.dmu,dcol.dr)

	return(list(DROW=drow,DCOL=dcol))
}

BM.SimuData <- function(n, m, mu, r, p, q)
{
	K <- length(p)
	L <- length(q)


	rep.row <- rmultinom(1, size=n, prob=p)
	rep.col <- rmultinom(1, size=m, prob=q)



	Z <- NULL
	for (k in 1:K)
	{
		if(rep.row[k]!=0)
		{
			Y <- NULL
			for (l in 1:L)
			{
				tmp1 <- rnbinom(rep.row[k]*rep.col[l], size = r[k,l], mu = mu[k,l])
				tmp2 <- array(tmp1,dim=c(rep.row[k],rep.col[l]))
				Y <- cbind(Y,tmp2)
			}
			Z <- rbind(Z,Y)
		}

	}

	Z <- as.array(Z)
	return(list(data = Z, reprow = rep.row, repcol = rep.col))
}


BM.BIC <- function(data, K, L, rep=30)
{
	n <- dim(data)[1]
	m <- dim(data)[2]

	emResult <- BM.RunEM(data, K, L)

	mu <- emResult$mu
	r <- emResult$r
	p <- emResult$p
	q <- emResult$q

	H <- NULL
	J <- NULL
	for(i in 1:rep)
	{
		tmpdata <- BM.SimuData.alt(n, m, mu, r, p, q)
		#tmpdata <- BM.SimuData(n, m, mu, r, p, q)
		BICResult <- BM.BICsimu(tmpdata, mu, r, p, q, K, L)
		drow <- BICResult$DROW
		dcol <- BICResult$DCOL

		cat("Simulation:",i,"\n")
		#cat("dRow:",drow,"\n")
		#cat("dCol:",dcol,"\n")
		#cat("\n")

		tmp1 <- t(drow)%*%drow+t(dcol)%*%dcol
		H <- cbind(H,as.vector(tmp1))
		tmp2 <- (colSums(drow)+colSums(dcol))%*%t(colSums(drow)+colSums(dcol))
		J <- cbind(J,as.vector(tmp2))
	}
	tmp <- K+L-2+2*K*L
	H.simu <- apply(H,1,mean)
	J.simu <- apply(J,1,mean)
	H.simu <- array(H.simu,dim=c(tmp,tmp))
	J.simu <- array(J.simu,dim=c(tmp,tmp))

	lambda <- 5e-5
	H.simu.1 <- H.simu + diag(rep(lambda, tmp)) 
	tmpdiag <- diag(J.simu%*%solve(H.simu.1))

	BIC <- -2*emResult$LogLikelihood + (log(n)+log(m))*sum(tmpdiag)

	cat("trace:",sum(tmpdiag),"\n")
	cat("relative DF:",sum(tmpdiag)/tmp,"\n")
	cat("BIC:",BIC,"\n")
	cat("LogLikelihood:",emResult$LogLikelihood,"\n")
	cat("\n")
	return(list(H = H, J = J, Hsimu=H.simu, Jsimu=J.simu, BIC = BIC , est = emResult))
}

PlotBIC <- function(BIC, startN=2)
{
	par(mar=(c(4,4,2,2)+.1),cex=1.3)
	Nk <- length(BIC)
	plot(startN :  Nk, BIC[startN:Nk] / 1e+4 , type="b",pch = 5 , col= "blue",
		lty = "solid",main="",xlab = "# of cluster", ylab = "BIC", cex=0.9)

	mtext(expression("x"*10^4),line = 0.3,adj = 0,cex=1.2)
	
}


BM.SimuData.alt <- function(n, m, mu, r, p, q)
{
	K <- length(p)
	L <- length(q)


	rep.row <- rmultinom(1, size=n, prob=p)
	rep.col <- rmultinom(1, size=m, prob=q)



	Z <- NULL
	for (k in 1:K)
	{
		if(rep.row[k]!=0)
		{
			Y <- NULL
			for (l in 1:L)
			{
				tmp1 <- rnbinom(rep.row[k]*rep.col[l], size = r[k,l], mu = mu[k,l])
				tmp2 <- array(tmp1,dim=c(rep.row[k],rep.col[l]))
				Y <- cbind(Y,tmp2)
			}
			Z <- rbind(Z,Y)
		}

	}

	Z <- as.array(Z)
	return(Z)
}

