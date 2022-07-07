data {
	int<lower=0> N;										// number of observations
	int<lower=0> Npat;									// number of individuals	
	real<lower=0> y[N];									// outcome data
	real<lower=0> age[N];								// explanatory variable (age) data
	int<lower=1,upper=Npat> id[N];						// patient id for each observation
	real<lower=0> fixedkp;								// fixed knot point (if determined a priori)
}

parameters {
	vector[3] beta;										// fixed effects
	vector<lower=0>[3] u_sd;							// level 2 error sd (sds of the random effects u[j])
	real<lower=0> y_sd;									// level 1 error sd
	vector[3] uRaw[Npat];								// working parameters for random effects (level 2) errors
	cholesky_factor_corr[3] L_u_Corr;					// cholesky factor for the random effects correlation matrix
}

transformed parameters {	
	vector[3] u[Npat];									// random effects (level 2) errors
	vector[3] alpha[Npat];								// random effects
	real y_mu[N];										// regression equation

	//==========================
	// calculate random effects
	//==========================
	
	// matts trick
	for (i in 1:Npat) u[i] <- u_sd .* (L_u_Corr * uRaw[i]); 	// equivalent to: diag_matrix(u_sd) * L_u_Corr * uRaw[i];

	for (i in 1:Npat) for (k in 1:3) alpha[i,k] <- beta[k] + u[i,k];

	//=====================
	// regression equation
	//=====================

	for (j in 1:N) {
		if (age[j] < fixedkp) 
			y_mu[j] <- alpha[id[j],1] + alpha[id[j],2] * (age[j] - fixedkp);
		else 	
			y_mu[j] <- alpha[id[j],1] + alpha[id[j],3] * (age[j] - fixedkp);			
	}	
}

model {

	//========
	// priors
	//========
	
	beta[1] ~ normal(20, 20);							// prior: fixed effect, intercept
	beta[2] ~ normal(0, 4);								// prior: fixed effect, slope before knot
	beta[3] ~ normal(0, 4);								// prior: fixed effect, slope after knot

	u_sd[1] ~ cauchy(0,5);								// prior: random effect sd, intercept
	u_sd[2] ~ cauchy(0,5);								// prior: random effect sd, slope before knot
	u_sd[3] ~ cauchy(0,5);								// prior: random effect sd, slope after knot
	
	y_sd ~ cauchy(0,5);									// prior: level 1 error sd
	
	L_u_Corr ~ lkj_corr_cholesky(1);					// prior: cholesky factor for random effects correlation matrix
														// NB. this prior is the "lkj correlation distribution" with shape parameter 1 
														// which is equivalent to a uniform distribution over the possible correlation 
														// matrices (where a shape parameter > 1 would have resulted in an upside down
														// U-shaped distribution with the mode being located at the identity matrix)
	
	//=============================
	// random effects distribution
	//=============================
	
	// matts trick
	for (i in 1:Npat) uRaw[i] ~ normal(0,1); 			// implies: u ~ multi_normal(0, diag_matrix(u_sd) * L_u_Corr * L_u_Corr' * diag_matrix(u_sd))
  	
	//==================
	// model likelihood
	//==================
	
	y ~ normal(y_mu, y_sd);								// likelihood for the observed data
	
}

generated quantities {
	corr_matrix[3] u_Corr;								// random effects correlation matrix
	matrix[3,3] u_Sigma;								// random effects covariance matrix
	
	//=====================================================
	// recover the correlation and covariance matrices
	// using the cholesky factor of the correlation matrix
	//=====================================================
	
	u_Corr <- multiply_lower_tri_self_transpose(L_u_Corr);		// correlation matrix: u_Corr = L_u_Corr * L_u_Corr'
	u_Sigma <- quad_form_diag(u_Corr, u_sd);					// covariance matrix: u_Sigma = diag(u_sd) * u_Corr * diag(u_sd)
		
}

