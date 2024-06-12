#include <Rcpp.h>
#include <iostream>

#include <typeinfo>
using namespace Rcpp;


double ks_sample_gsva(IntegerVector geneset_idxs,
                      IntegerVector expr,
                      IntegerVector sort_idxs,
                      double tau,
                      int mx_diff,
                      int abs_rank){

  int n_genes = expr.size();
  int n_geneset = geneset_idxs.size();

  IntegerVector geneset_mask(n_genes);
  for(int i = 0; i < n_genes; i++){
    geneset_mask[i] = 0;
  }

  for(int i = 0; i < n_geneset; i++){
    geneset_mask[geneset_idxs[i]-1] = 1;
  }

  double dec = 1.0 / (n_genes - n_geneset);

  double sum_gset = 0.0;
  //Rprintf("%f\n", dec);
  for(int i = 0; i < n_geneset; i++){
    //Rprintf("%f\n",pow(x[geneset_idxs[i]-1], tau));
    sum_gset += pow(expr[geneset_idxs[i]-1], tau);
  }

  double mx_value_sign = 0.0;
  double cum_sum = 0.0;

  double mx_pos = 0.0;
  double mx_neg = 0.0;

  int idx;
  for(int i = 0; i < n_genes; i++){
    idx = sort_idxs[i]-1;

    if(geneset_mask[idx] == 1){
      //Rprintf("%d\t", i);
      cum_sum += pow(expr[idx], tau) / sum_gset;
    }else{
      cum_sum -= dec;
    }
    //walk[i] = cum_sum;

    if(cum_sum > mx_pos){ mx_pos = cum_sum; }
    if(cum_sum < mx_neg){ mx_neg = cum_sum; }
  }

  if (mx_diff != 0) {
    mx_value_sign = mx_pos + mx_neg;
    if (abs_rank != 0)
      mx_value_sign = mx_pos - mx_neg;
  } else {
    mx_value_sign = (mx_pos > fabs(mx_neg)) ? mx_pos : mx_neg;
  }

  return mx_value_sign;
}






// [[Rcpp::export]]
NumericVector ks_gset_gsva(List geneset_idxs_lis,
                           IntegerVector expr,
                           IntegerVector sort_idxs,
                           int n_gsets,
                           double tau,
                           int mx_diff,
                           int abs_rank){

  NumericVector pas_cell(n_gsets);
  for(int i =0; i< n_gsets; i++){
    IntegerVector gset = as<IntegerVector>(geneset_idxs_lis[i]);
    pas_cell[i] = ks_sample_gsva(gset,
                                 expr,
                                 sort_idxs,
                                 tau,
                                 mx_diff,
                                 abs_rank);
  }

  return pas_cell;

}
