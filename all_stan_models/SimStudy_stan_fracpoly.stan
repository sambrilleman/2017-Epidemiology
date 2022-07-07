data {
	int<lower=0> N;										// number of observations
	int<lower=0> Npat;									// number of individuals	
	real<lower=0> y[N];									// outcome data
	real<lower=0> age[N];								// explanatory variable (age)
	real<lower=0> var2[N];								// explanatory variable (age * ln(age))
	int<lower=1,upper=Npat> id[N];						// patient id for each observation
}

parameters {
	vector[3] beta;										// fixed effects, intercept and slopes
	vector<lower=0>[3] u_sd;							// level 2 error sd (sds of the random effects u[j])
	real<lower=0> y_sd;									// level 1 error sd
	vector[3] uRaw[Npat];								// working parameters for random effects (level 2) errors
	cholesky_factor_corr[3] L_u_Corr;					// cholesky factor for the random effects correlation matrix
}

transformed parameters {	
	vector[3] u[Npat];									// random effects (level 2) errors
	vector[3] alpha[Npat];								// random effects
	real y_mu[N];										// mean parameter based on regression equation

	//==========================
	// calculate random effects
	//==========================
	
	// matts trick
	for (i in 1:Npat) u[i] <- u_sd .* (L_u_Corr * uRaw[i]); 	// equivalent to: diag_matrix(u_sd) * L_u_Corr * uRaw[i];

	for (i in 1:Npat) {
		for (k in 1:3) alpha[i,k] <- beta[k] + u[i,k];
	}
	
	//=====================
	// regression equation
	//=====================
	
	for (j in 1:N)
		y_mu[j] <- alpha[id[j],1] + (alpha[id[j],2] * age[j]) + (alpha[id[j],3] * var2[j]);				
}

model {

	//========
	// priors
	//========
	
	beta[1] ~ normal(20, 20);							// prior: fixed effect, intercept
	beta[2] ~ normal(0, 4);								// prior: fixed effect, slope term 1
	beta[3] ~ normal(0, 4);								// prior: fixed effect, slope term 2

	u_sd[1] ~ cauchy(0, 5);								// prior: random effect sd, intercept
	u_sd[2] ~ cauchy(0, 5);								// prior: random effect sd, slope term 1
	u_sd[3] ~ cauchy(0, 5);								// prior: random effect sd, slope term 2
	
	y_sd ~ cauchy(0, 5);								// prior: level 1 error sd
	
	L_u_Corr ~ lkj_corr_cholesky(1);					// prior: cholesky factor for random effects correlation matrix
														// NB. this prior is the "lkj correlation distribution" with shape parameter 1 
														// which is equivalent to a uniform distribution over the possible correlation 
														// matrices (where a shape parameter > 1 would have resulted in an upside down
														// U-shaped distribution with the mode being located at the identity matrix)

	//=============================
	// random effects distribution
	//=============================

	// matts trick
	for (i in 1:Npat) uRaw[i] ~ normal(0,1); 		// implies: u ~ multi_normal(0, diag_matrix(u_sd) * L_u_Corr * L_u_Corr' * diag_matrix(u_sd))
  	
	//==================
	// model likelihood
	//==================
	
	y ~ normal(y_mu, y_sd);								// likelihood for the observed data
	
}

generated quantities {
	corr_matrix[3] u_Corr;								// random effects correlation matrix
	matrix[3,3] u_Sigma;								// random effects covariance matrix
	real meankp;										// knot point for curve based only on fixed effects
	int meankp_minimum;									// indicator for whether the 2nd derivative is positive
	real kp[Npat];										// knot point for each individual
	int kp_minimum[Npat];								// indicator for whether the 2nd derivative is positive
	real meankp2;										// knot point for curve based on averaging knot points across individuals
	
	//=====================================================
	// recover the correlation and covariance matrices
	// using the cholesky factor of the correlation matrix
	//=====================================================
	
	u_Corr <- multiply_lower_tri_self_transpose(L_u_Corr);		// correlation matrix: u_Corr = L_u_Corr * L_u_Corr'
	u_Sigma <- quad_form_diag(u_Corr, u_sd);					// covariance matrix: u_Sigma = diag(u_sd) * u_Corr * diag(u_sd)

	//==========================================================================
	// calculate knot point based on analytical derivative for the fitted curve
	//==========================================================================
	
	meankp <- exp(-(beta[2] + beta[3]) / beta[3]);
	if (beta[3] > 0) meankp_minimum <- 1;
	else meankp_minimum <- 0;
	
	for (i in 1:Npat) {
		kp[i] <- exp(-(alpha[i,2] + alpha[i,3]) / alpha[i,3]);
		if (alpha[i,3] > 0) kp_minimum[i] <- 1; 
		else kp_minimum[i] <- 0; 
	}
	meankp2 <- sum(kp) / Npat; 
	
}









