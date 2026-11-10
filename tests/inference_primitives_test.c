#include "../src/statistics/numelyx_core_inference.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
static void close_to(double a,double b,double t){assert(fabs(a-b)<=t);}
static void test_delta_method(void){
 const double gradient[]={1,2};
 const double covariance[]={1,.5,.5,4};
 double variance=99,se=99;
 assert(numelyx_core_delta_method_variance(gradient,covariance,2,&variance)==NUMELYX_INFERENCE_SUCCESS);
 close_to(variance,19,1e-14);
 assert(numelyx_core_delta_method_standard_error(gradient,covariance,2,&se)==NUMELYX_INFERENCE_SUCCESS);
 close_to(se,sqrt(19),1e-14);
}
static void test_wald_scalar_and_joint(void){
 const double estimate[]={1,2};
 const double null_scalar[]={0};
 const double restriction_scalar[]={1,1};
 const double covariance[]={1,0,0,4};
 const double null_joint[]={0,0};
 const double restriction_joint[]={1,0,0,1};
 double statistic=99,logp=99;
 assert(numelyx_core_wald_test(estimate,null_scalar,restriction_scalar,
  covariance,2,1,&statistic,&logp)==NUMELYX_INFERENCE_SUCCESS);
 close_to(statistic,9.0/5.0,1e-14);
 {
  double joint_statistic=99,joint_logp=99;
  assert(numelyx_core_wald_test(estimate,null_joint,restriction_joint,
   covariance,2,2,&joint_statistic,&joint_logp)==NUMELYX_INFERENCE_SUCCESS);
  close_to(joint_statistic,2.0,1e-14);
  close_to(exp(joint_logp),exp(-1.0),2e-14);
 }
}
static void test_inference_validation(void){
 const double estimate[]={1,2};
 const double null_value[]={0,0};
 const double restriction[]={1,0,2,0};
 const double singular_covariance[]={1,0,0,0};
 const double gradient[]={1,0};
 const double covariance[]={-1,0,0,1};
 double statistic=7,logp=8,variance=9;
 assert(numelyx_core_wald_test(estimate,null_value,restriction,
  singular_covariance,2,2,&statistic,&logp)==NUMELYX_INFERENCE_SINGULAR_COVARIANCE);
 assert(statistic==7&&logp==8);
 assert(numelyx_core_delta_method_variance(gradient,covariance,2,&variance)
  ==NUMELYX_INFERENCE_NUMERICAL_FAILURE);
 assert(variance==9);
}
int main(void){
 test_delta_method();test_wald_scalar_and_joint();test_inference_validation();
 puts("Inference primitive tests passed.");return 0;
}
