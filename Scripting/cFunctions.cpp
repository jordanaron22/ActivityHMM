//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp; 
using namespace arma;



// [[Rcpp::export]]
NumericVector vectorEqBool(NumericVector vec, double lod) {
	NumericVector vecBool;
	for(NumericVector::iterator i = vec.begin(); i != vec.end(); ++i) {
		if (*i == lod) {
			vecBool.push_back(1);
		} else {
			vecBool.push_back(0);
		}
  }
  return vecBool;
}


// [[Rcpp::export]]
NumericVector logClassificationC(int current_state, NumericVector act_obs, double mu, double sig, double act_binom, double lod) {
	NumericVector temp;
	if (current_state == 0) {
		temp = Rcpp::dnorm( act_obs, mu, sig, true );
	} else {
		NumericVector vec_eq = vectorEqBool(act_obs, lod);
		temp = (log((1-act_binom) * Rcpp::dnorm( act_obs, mu, sig, false )) * (1 - vec_eq)) + (log(act_binom)*vec_eq);
		
	}

	LogicalVector nan_vec = is_na(act_obs);
	int n = act_obs.length();

	for (int i = 0; i < n; i++) {
		if(nan_vec[i] == true){
			temp[i] = 0;
		}
	}


	return temp;
}

// [[Rcpp::export]]
vec logClassificationC2(int current_state, NumericVector act_obs, double mu, double sig, double act_binom, double lod) {
  vec temp;
  if (current_state == 0) {
    temp = Rcpp::dnorm( act_obs, mu, sig, true );
  } else {
    NumericVector vec_eq = vectorEqBool(act_obs, lod);
    temp = (log((1-act_binom) * Rcpp::dnorm( act_obs, mu, sig, false )) * (1 - vec_eq)) + (log(act_binom)*vec_eq);
  }

  LogicalVector nan_vec = is_na(act_obs);
	int n = act_obs.length();

	for (int i = 0; i < n; i++) {
		if(nan_vec[i] == true){
			temp[i] = 0;
		}
	}

  return temp;
}

// [[Rcpp::export]]
double logSumExpC2(NumericVector lx){

  // Obtain environment containing function
  Rcpp::Environment package_env("package:matrixStats"); 

  // Make function callable from C++
  Rcpp::Function rfunction = package_env["logSumExp"];    

  // Call the function and receive output (might not be list)
  Rcpp::NumericVector test_out = rfunction(Rcpp::_["lx"] = lx);
  double lse = test_out(0);
  return lse;
}

/* This is from the seqHMM github*/
#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#  define EXPL expl
#else
#  define LDOUBLE double
#  define EXPL exp
#endif

// [[Rcpp::export]]
double logSumExpC(const arma::vec& x) {
  unsigned int maxi = x.index_max();
  LDOUBLE maxv = x(maxi);
  if (!(maxv > -arma::datum::inf)) {
    return -arma::datum::inf;
  }
  LDOUBLE cumsum = 0.0;
  for (unsigned int i = 0; i < x.n_elem; i++) {
    if ((i != maxi) && (x(i) > -arma::datum::inf)) {
      cumsum += EXPL(x(i) - maxv);
    }
  }
  
  return maxv + log1p(cumsum);
}

// [[Rcpp::export]]
mat ForwardIndC(const NumericVector& act_ind, NumericVector init, Rcpp::List tran_list, cube emit_act,int tran_ind, double act_light_binom, int clust_i, double lepsilon, double log_sweight){

	mat alpha( act_ind.length(), 2 );

	NumericVector log_class_0 = logClassificationC( 0, act_ind, emit_act(0,0,clust_i), emit_act(0,1,clust_i), act_light_binom, lepsilon );

	NumericVector log_class_1 = logClassificationC( 1, act_ind, emit_act(1,0,clust_i), emit_act(1,1,clust_i), act_light_binom, lepsilon );

	alpha(0,0) = log(init(0)) + log_class_0[0];
	alpha(0,1) = log(init(1)) + log_class_1[0];
	
	List tran_time_list = tran_list[tran_ind - 1]; 
	
	for (int i = 1; i < act_ind.length(); i++) {
	  
	  NumericMatrix tran = tran_time_list[i%96];
	  double fp_00 = alpha(i-1,0) + log(tran(0,0)) + log_class_0[i];
	  
	  double fp_10 = alpha(i-1,1) + log(tran(1,0)) + log_class_0[i];
	  
	  double fp_01 = alpha(i-1,0) + log(tran(0,1)) + log_class_1[i];
	  
	  double fp_11 = alpha(i-1,1) + log(tran(1,1)) + log_class_1[i];
	  
	  NumericVector fp_0 = NumericVector::create(fp_00,fp_10);
	  NumericVector fp_1 = NumericVector::create(fp_01,fp_11);
	  
	  alpha(i,0) = logSumExpC(fp_0);
	  alpha(i,1) = logSumExpC(fp_1);
	}

	alpha = alpha + log_sweight;

	return(alpha);
}

// [[Rcpp::export]]
mat BackwardIndC(const NumericVector& act_ind, Rcpp::List tran_list, cube emit_act,int tran_ind, double act_light_binom, int clust_i, double lepsilon){
  
  int n = act_ind.length(); 
  mat beta( n, 2 );
  
  NumericVector log_class_0 = logClassificationC( 0, act_ind, emit_act(0,0,clust_i), emit_act(0,1,clust_i), act_light_binom, lepsilon );
  
  NumericVector log_class_1 = logClassificationC( 1, act_ind, emit_act(1,0,clust_i), emit_act(1,1,clust_i), act_light_binom, lepsilon );
  
  beta(n-1,0) = log(1);
  beta(n-1,1) = log(1);
  
  List tran_time_list = tran_list[tran_ind - 1]; 
  
  for (int i = n-2; i >= 0; i--) {
    
    NumericMatrix tran = tran_time_list[(i+1)%96];
    
    double bp_00 = log(tran(0,0)) + log_class_0[i+1] + beta(i+1,0);
    
    double bp_01 = log(tran(0,1)) + log_class_1[i+1] + beta(i+1,1);
    
    double bp_10 = log(tran(1,0)) + log_class_0[i+1] + beta(i+1,0);
    
    double bp_11 = log(tran(1,1)) + log_class_1[i+1] + beta(i+1,1);
    
    NumericVector bp_0 = NumericVector::create(bp_00,bp_01);
    NumericVector bp_1 = NumericVector::create(bp_10,bp_11);
    
    beta(i,0) = logSumExpC(bp_0);
    beta(i,1) = logSumExpC(bp_1);

  }
  
  return beta;

}

// [[Rcpp::export]]
List ForwardC(const NumericMatrix& act, NumericVector init, List tran_list, cube emit_act, NumericVector tran_ind_vec,double act_light_binom, double lepsilon, NumericVector log_sweights_vec){
	int num_people = act.ncol();
	int len = act.nrow();
	int num_re = emit_act.n_slices;
	List alpha_list(num_people);
	
	for (int ind = 0; ind < num_people; ind++) {
		arma::cube Cube1(len, 2, num_re);
		NumericVector act_ind = act.column(ind);
		int tran_ind = tran_ind_vec(ind);
		double log_sweight = log_sweights_vec(ind);

		
		for (int clust_i = 0; clust_i < num_re; clust_i++){
			Cube1.slice(clust_i) = ForwardIndC(act_ind, init, tran_list, emit_act, tran_ind, act_light_binom, clust_i, lepsilon, log_sweight);

		}

		alpha_list(ind) = Cube1;
	}
	return(alpha_list);
}

// [[Rcpp::export]]
List BackwardC(const NumericMatrix& act, List tran_list, cube emit_act, NumericVector tran_ind_vec,double act_light_binom, double lepsilon){
	int num_people = act.ncol();
	int len = act.nrow();
	int num_re = emit_act.n_slices;
	List beta_list(num_people);
	
	for (int ind = 0; ind < num_people; ind++) {
		arma::cube Cube1(len, 2, num_re);
		NumericVector act_ind = act.column(ind);
		int tran_ind = tran_ind_vec(ind);
		
		for (int clust_i = 0; clust_i < num_re; clust_i++){
			Cube1.slice(clust_i) = BackwardIndC(act_ind, tran_list, emit_act, tran_ind, act_light_binom, clust_i, lepsilon);

		}

		beta_list(ind) = Cube1;
	}
	return(beta_list);
}  

// [[Rcpp::export]]
cube CalcTranHelperC(int init_state, int new_state, NumericMatrix act, List tran_list_mat, NumericVector tran_ind_vec, cube emit_act, NumericVector ind_like_vec, List alpha, List beta,double act_light_binom, double lepsilon, vec pi_l){
  int num_people = act.ncol();
  int len = act.nrow();
  int num_re = emit_act.n_slices;

  arma::cube tran_vals_re_cube( len-1, num_people, num_re );

  for (int clust_i = 0; clust_i < num_re; clust_i++){

	mat tran_vals_re_mat( len-1, num_people );

    for (int ind = 0; ind < num_people; ind++) {
      
      int tran_ind = tran_ind_vec(ind);
	  mat tran_mat = tran_list_mat(tran_ind-1);
	  

	  //0,0->0 & 1,0->1 & 0,1->2 & 1,1->3
	  int tran_vec_ind = init_state + (new_state * 2);

	  arma::cube alpha_ind = alpha(ind); 
	  arma::cube beta_ind = beta(ind);
	  double likelihood = ind_like_vec(ind);
	  



	  NumericMatrix act_ind = act( Range(1,len-1) , Range(ind,ind) );
	  NumericVector act_ind_m1 = act_ind.column(0); 
	  
	  


	  vec class_vec = logClassificationC2( new_state, act_ind_m1, emit_act(new_state,0,clust_i), emit_act(new_state,1,clust_i), act_light_binom, lepsilon );
	  
	  vec alpha_ind_slice = alpha_ind(span(0,len-2),span(init_state,init_state),span(clust_i,clust_i));
	  vec beta_ind_slice = beta_ind(span(1,len-1),span(new_state,new_state),span(clust_i,clust_i));
	  vec tran_vec_slice = tran_mat.col(tran_vec_ind);

	  vec temp = alpha_ind_slice + beta_ind_slice + log(tran_vec_slice) + log(pi_l(clust_i)) + class_vec - likelihood;
	  vec tran_vals_re_ind = arma::exp(temp);

	  tran_vals_re_mat.col(ind) = tran_vals_re_ind;
	
    }

	tran_vals_re_cube.slice(clust_i) = tran_vals_re_mat;
    
  }

  return tran_vals_re_cube;
}

