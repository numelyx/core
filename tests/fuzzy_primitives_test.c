#include "../src/fuzzy/numelyx_core_fuzzy.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
static void close_to(double a,double b,double t){assert(fabs(a-b)<=t);}
static void test_membership_functions(void){
 double mu;
 assert(numelyx_core_fuzzy_triangular_membership(1,0,1,2,&mu)==NUMELYX_FUZZY_SUCCESS);
 assert(mu==1);
 assert(numelyx_core_fuzzy_triangular_membership(.5,0,1,2,&mu)==NUMELYX_FUZZY_SUCCESS);
 close_to(mu,.5,1e-15);
 assert(numelyx_core_fuzzy_trapezoidal_membership(1.5,0,1,2,3,&mu)==NUMELYX_FUZZY_SUCCESS);
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
 assert(numelyx_core_fuzzy_trapezoidal_alpha_cut(0,1,2,3,1,&lower,&upper)==NUMELYX_FUZZY_SUCCESS);
 close_to(lower,1,1e-15); close_to(upper,2,1e-15);
 lower=9;upper=10;
 assert(numelyx_core_fuzzy_triangular_alpha_cut(0,1,2,1.5,&lower,&upper)==NUMELYX_FUZZY_DOMAIN_ERROR);
 assert(lower==9&&upper==10);
 assert(numelyx_core_fuzzy_gaussian_membership(0,0,0,&lower)==NUMELYX_FUZZY_DOMAIN_ERROR);
 assert(numelyx_core_fuzzy_tnorm_product(NAN,.5,&lower)==NUMELYX_FUZZY_NON_FINITE_INPUT);
}
int main(void){
 test_membership_functions(); test_fuzzy_operators(); test_alpha_cuts_and_validation();
 puts("Fuzzy primitive tests passed."); return 0;
}
