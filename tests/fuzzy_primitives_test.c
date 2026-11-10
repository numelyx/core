#include "../src/fuzzy/numelyx_core_fuzzy.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
static void close_to(double a,double b,double t){
 if (!(fabs(a-b)<=t)) fprintf(stderr,"fuzzy close_to failed actual=%.17g expected=%.17g tolerance=%.17g\\n",a,b,t);
 assert(fabs(a-b)<=t);
}
static void test_membership_functions(void){
 double mu;
 assert(numelyx_core_fuzzy_triangular_membership(1,0,1,2,&mu)==NUMELYX_FUZZY_SUCCESS);
 assert(mu==1);
 assert(numelyx_core_fuzzy_triangular_membership(.5,0,1,2,&mu)==NUMELYX_FUZZY_SUCCESS);
 close_to(mu,.5,1e-15);
 assert(numelyx_core_fuzzy_trapezoidal_membership(1.5,0,1,2,3,&mu)==NUMELYX_FUZZY_SUCCESS);
 assert(mu==1);
 assert(numelyx_core_fuzzy_triangular_membership(0,0,0,2,&mu)==NUMELYX_FUZZY_SUCCESS);
 assert(mu==1);
 assert(numelyx_core_fuzzy_trapezoidal_membership(0,0,0,1,2,&mu)==NUMELYX_FUZZY_SUCCESS);
 assert(mu==1);
 assert(numelyx_core_fuzzy_trapezoidal_membership(.5,0,1,2,3,&mu)==NUMELYX_FUZZY_SUCCESS);
 close_to(mu,.5,1e-15);
 assert(numelyx_core_fuzzy_gaussian_membership(1,0,1,&mu)==NUMELYX_FUZZY_SUCCESS);
 close_to(mu,exp(-.5),1e-15);
 assert(numelyx_core_fuzzy_bell_membership(2,2,2,0,&mu)==NUMELYX_FUZZY_SUCCESS);
 close_to(mu,.5,1e-15);
 assert(numelyx_core_fuzzy_sigmoid_membership(0,1,0,&mu)==NUMELYX_FUZZY_SUCCESS);
 assert(mu==.5);
 assert(numelyx_core_fuzzy_sigmoid_membership(1000,1,0,&mu)==NUMELYX_FUZZY_SUCCESS);
 assert(mu==1.0);
}
static void test_fuzzy_operators(void){
 double value;
 assert(numelyx_core_fuzzy_complement(.2,&value)==NUMELYX_FUZZY_SUCCESS);
 close_to(value,.8,1e-15);
 assert(numelyx_core_fuzzy_tnorm_min(.2,.7,&value)==NUMELYX_FUZZY_SUCCESS);
 close_to(value,.2,1e-15);
 assert(numelyx_core_fuzzy_tnorm_product(.2,.7,&value)==NUMELYX_FUZZY_SUCCESS);
 close_to(value,.14,1e-15);
 assert(numelyx_core_fuzzy_tnorm_lukasiewicz(.2,.7,&value)==NUMELYX_FUZZY_SUCCESS);
 assert(value==0.0);
 assert(numelyx_core_fuzzy_snorm_max(.2,.7,&value)==NUMELYX_FUZZY_SUCCESS);
 close_to(value,.7,1e-15);
 assert(numelyx_core_fuzzy_snorm_probabilistic_sum(.2,.7,&value)==NUMELYX_FUZZY_SUCCESS);
 close_to(value,.76,1e-15);
 assert(numelyx_core_fuzzy_snorm_bounded_sum(.6,.7,&value)==NUMELYX_FUZZY_SUCCESS);
 assert(value==1.0);
}
static void test_alpha_cuts_and_validation(void){
 double lower=9,upper=10;
 assert(numelyx_core_fuzzy_triangular_alpha_cut(0,1,2,.5,&lower,&upper)==NUMELYX_FUZZY_SUCCESS);
 close_to(lower,.5,1e-15); close_to(upper,1.5,1e-15);
 assert(numelyx_core_fuzzy_triangular_alpha_cut(-DBL_MAX,0,DBL_MAX,.5,&lower,&upper)==NUMELYX_FUZZY_SUCCESS);
 assert(isfinite(lower)&&isfinite(upper));
 assert(numelyx_core_fuzzy_trapezoidal_alpha_cut(0,1,2,3,1,&lower,&upper)==NUMELYX_FUZZY_SUCCESS);
 close_to(lower,1,1e-15); close_to(upper,2,1e-15);
 lower=9;upper=10;
 assert(numelyx_core_fuzzy_triangular_alpha_cut(0,1,2,1.5,&lower,&upper)==NUMELYX_FUZZY_DOMAIN_ERROR);
 assert(lower==9&&upper==10);
 assert(numelyx_core_fuzzy_gaussian_membership(0,0,0,&lower)==NUMELYX_FUZZY_DOMAIN_ERROR);
 assert(numelyx_core_fuzzy_tnorm_product(NAN,.5,&lower)==NUMELYX_FUZZY_NON_FINITE_INPUT);
}

static void test_fuzzy_interval_arithmetic(void){
 double lo=91,hi=92;
 assert(numelyx_core_fuzzy_interval_add(1,2,3,4,&lo,&hi)==NUMELYX_FUZZY_SUCCESS);
 close_to(lo,4,1e-15); close_to(hi,6,1e-15);
 assert(numelyx_core_fuzzy_interval_subtract(1,2,3,4,&lo,&hi)==NUMELYX_FUZZY_SUCCESS);
 close_to(lo,-3,1e-15); close_to(hi,-1,1e-15);
 assert(numelyx_core_fuzzy_interval_multiply(-2,3,-4,5,&lo,&hi)==NUMELYX_FUZZY_SUCCESS);
 close_to(lo,-12,1e-15); close_to(hi,15,1e-15);
 assert(numelyx_core_fuzzy_interval_divide(1,2,3,4,&lo,&hi)==NUMELYX_FUZZY_SUCCESS);
 close_to(lo,.25,1e-15); close_to(hi,2.0/3.0,1e-15);
 lo=91;hi=92;
 assert(numelyx_core_fuzzy_interval_divide(1,2,-1,1,&lo,&hi)==NUMELYX_FUZZY_DOMAIN_ERROR);
 assert(lo==91&&hi==92);
 assert(numelyx_core_fuzzy_interval_add(DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,&lo,&hi)==NUMELYX_FUZZY_NUMERICAL_FAILURE);
 assert(lo==91&&hi==92);
}


static void test_fuzzy_defuzzification(void){
 const double x[]={0,1,2,3,4};
 const double mu[]={0,.5,1,.5,0};
 const double plateau[]={0,.8,1,1,.8};
 double value=77;
 assert(numelyx_core_fuzzy_defuzzify_centroid(x,mu,5,&value)==NUMELYX_FUZZY_SUCCESS);
 close_to(value,2,1e-15);
 assert(numelyx_core_fuzzy_defuzzify_mean_of_maxima(x,plateau,5,&value)==NUMELYX_FUZZY_SUCCESS);
 close_to(value,2.5,1e-15);
 value=77;
 {
  const double zero[]={0,0,0,0,0};
  assert(numelyx_core_fuzzy_defuzzify_centroid(x,zero,5,&value)==NUMELYX_FUZZY_DOMAIN_ERROR);
  assert(value==77);
 }
}

static void test_mamdani_and_sugeno_inference_aggregation(void){
 const double strengths[]={.5,.8};
 const double consequents[]={0,.6,1,.4,.5,0};
 double output[]={91,92,93};
 double sugeno=94;
 assert(numelyx_core_fuzzy_mamdani_aggregate(strengths,consequents,2,3,
  NUMELYX_FUZZY_IMPLICATION_MIN,NUMELYX_FUZZY_AGGREGATION_MAX,output)==NUMELYX_FUZZY_SUCCESS);
 close_to(output[0],.4,1e-15); close_to(output[1],.5,1e-15); close_to(output[2],.5,1e-15);
 assert(numelyx_core_fuzzy_mamdani_aggregate(strengths,consequents,2,3,
  NUMELYX_FUZZY_IMPLICATION_PRODUCT,NUMELYX_FUZZY_AGGREGATION_PROBABILISTIC_SUM,output)==NUMELYX_FUZZY_SUCCESS);
 close_to(output[0],.32,1e-15); close_to(output[1],.58,1e-15); close_to(output[2],.5,1e-15);
 {
  const double crisp[]={10,20};
  assert(numelyx_core_fuzzy_sugeno_aggregate(strengths,crisp,2,&sugeno)==NUMELYX_FUZZY_SUCCESS);
  close_to(sugeno,21.0/1.3,1e-14);
 }
 {
  const double invalid[]={.5,1.1};
  output[0]=91;
  assert(numelyx_core_fuzzy_mamdani_aggregate(invalid,consequents,2,3,
   NUMELYX_FUZZY_IMPLICATION_MIN,NUMELYX_FUZZY_AGGREGATION_MAX,output)==NUMELYX_FUZZY_DOMAIN_ERROR);
  assert(output[0]==91);
 }
 {
  const double zero[]={0,0};
  sugeno=94;
  assert(numelyx_core_fuzzy_sugeno_aggregate(zero,strengths,2,&sugeno)==NUMELYX_FUZZY_DOMAIN_ERROR);
  assert(sugeno==94);
 }
}

int main(void){
 test_membership_functions(); test_fuzzy_operators(); test_alpha_cuts_and_validation(); test_fuzzy_interval_arithmetic(); test_fuzzy_defuzzification(); test_mamdani_and_sugeno_inference_aggregation();
 puts("Fuzzy primitive tests passed."); return 0;
}
