#include <Rcpp.h>
#include <iostream>

#include <typeinfo>
using namespace Rcpp;





double ks_sample_ssgsea(IntegerVector geneset_idxs,
                        IntegerVector expr,
                        IntegerVector sort_idxs,
                        double tau){

  int n_genes = expr.size();
  //int n_samples = expr.ncol();
  int n_geneset = geneset_idxs.size();
  IntegerVector geneset_mask(n_genes);
  //NumericVector es_geneset(n_samples);
  for(int i = 0; i < n_genes; i++){
    geneset_mask[i] = 0;
  }

  for(int i = 0; i < n_geneset; i++){
    geneset_mask[geneset_idxs[i]-1] = 1;
  }

  double dec = 1.0 / (n_genes - n_geneset);


  double walk_sum=0;

  double sum_gset = 0.0;

  for(int i = 0; i < n_geneset; i++){
    //Rprintf("%f\n",pow(x[geneset_idxs[i]-1], tau));
    sum_gset += pow(expr[geneset_idxs[i]-1], tau);
  }

  double cum_sum = 0.0;


  int idx;
  for(int i = 0; i < n_genes; i++){
    idx = sort_idxs[i]-1;
    if(geneset_mask[idx] == 1){
      cum_sum += pow(expr[idx], tau) / sum_gset;
    }else{
      cum_sum -= dec;
    }
    walk_sum += cum_sum;
  }

  return walk_sum;

}



// [[Rcpp::export]]
NumericVector ks_gset_ssgsea(List geneset_idxs_lis,
                             IntegerVector expr,
                             IntegerVector sort_idxs,
                             int n_gsets,
                             double tau){
  //int n_gsets = geneset_idxs_lis.size();
  NumericVector pas_cell(n_gsets);
  for(int i =0; i<n_gsets; i++){
    //Rprintf(geneset_idxs_lis["Glycolysis_Gluconeogenesis"]);
    //char name=gsets_name[i];
    //Rprintf("%s\n", typeid(gsets_name[i]).name());
    IntegerVector gset = as<IntegerVector>(geneset_idxs_lis[i]);
    pas_cell[i] = ks_sample_ssgsea(gset,
                                   expr,
                                   sort_idxs,
                                   tau);
    //NumericVector vsc = as<NumericVector>(geneset_idxs_lis[i]);
    //pas_cell[i] = vsc[0];
  }

  return pas_cell;

}
