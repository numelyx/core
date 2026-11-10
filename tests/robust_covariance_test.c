#include "../src/statistics/numelyx_core_covariance.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
static void close_to(double a,double b,double t){assert(fabs(a-b)<=t);}
static void test_hc0_hc1_sandwich(void){
 const double bread[]={.5};
 const double scores[]={1,2,3};
 double covariance=99;
 assert(numelyx_core_sandwich_covariance(bread,scores,3,1,
  NUMELYX_COVARIANCE_HC0,&covariance)==NUMELYX_COVARIANCE_SUCCESS);
 close_to(covariance,3.5,1e-14);
 assert(numelyx_core_sandwich_covariance(bread,scores,3,1,
  NUMELYX_COVARIANCE_HC1,&covariance)==NUMELYX_COVARIANCE_SUCCESS);
 close_to(covariance,5.25,1e-14);
}
static void test_multivariate_sandwich_orientation(void){
 const double bread[]={1,.5,0,2};
 const double scores[]={1,2,3,4};
 double covariance[4]={0,0,0,0};
 assert(numelyx_core_sandwich_covariance(
  bread,scores,2,2,NUMELYX_COVARIANCE_HC0,covariance)
  ==NUMELYX_COVARIANCE_SUCCESS);
 close_to(covariance[0],29,1e-14);
 close_to(covariance[1],48,1e-14);
 close_to(covariance[2],48,1e-14);
 close_to(covariance[3],80,1e-14);
}
static void test_cluster_and_hac_sandwich(void){
 const double bread[]={.5};
 const double scores[]={1,2,3};
 const size_t clusters[]={0,0,1};
 double covariance=99;
 assert(numelyx_core_cluster_sandwich_covariance(
  bread,scores,clusters,3,1,2,NUMELYX_COVARIANCE_CLUSTER_CR0,&covariance)
  ==NUMELYX_COVARIANCE_SUCCESS);
 close_to(covariance,4.5,1e-14);
 assert(numelyx_core_cluster_sandwich_covariance(
  bread,scores,clusters,3,1,2,NUMELYX_COVARIANCE_CLUSTER_CR1,&covariance)
  ==NUMELYX_COVARIANCE_SUCCESS);
 close_to(covariance,9.0,1e-14);
 assert(numelyx_core_hac_sandwich_covariance(
  bread,scores,3,1,1,&covariance)==NUMELYX_COVARIANCE_SUCCESS);
 close_to(covariance,5.5,1e-14);
}
static void test_hc2_hc3_regression_covariance(void){
 const double bread[]={.5};
 const double design[]={1,1,1};
 const double residuals[]={1,2,3};
 const double leverage[]={0,.5,.75};
 double covariance=99;
 assert(numelyx_core_hc_sandwich_covariance(bread,design,residuals,leverage,
  3,1,NUMELYX_COVARIANCE_HC_TYPE_0,&covariance)==NUMELYX_COVARIANCE_SUCCESS);
 close_to(covariance,3.5,1e-14);
 assert(numelyx_core_hc_sandwich_covariance(bread,design,residuals,leverage,
  3,1,NUMELYX_COVARIANCE_HC_TYPE_1,&covariance)==NUMELYX_COVARIANCE_SUCCESS);
 close_to(covariance,5.25,1e-14);
 assert(numelyx_core_hc_sandwich_covariance(bread,design,residuals,leverage,
  3,1,NUMELYX_COVARIANCE_HC_TYPE_2,&covariance)==NUMELYX_COVARIANCE_SUCCESS);
 close_to(covariance,11.25,1e-14);
 assert(numelyx_core_hc_sandwich_covariance(bread,design,residuals,leverage,
  3,1,NUMELYX_COVARIANCE_HC_TYPE_3,&covariance)==NUMELYX_COVARIANCE_SUCCESS);
 close_to(covariance,40.25,1e-14);
}
static void test_covariance_validation_and_output_preservation(void){
 const double bread[]={1};
 const double scores[]={1,2};
 const size_t invalid_clusters[]={0,2};
 double covariance=17;
 assert(numelyx_core_sandwich_covariance(
  bread,scores,2,1,NUMELYX_COVARIANCE_HC0,&covariance)
  ==NUMELYX_COVARIANCE_SUCCESS);
 close_to(covariance,5.0,1e-14);
 covariance=17;
 assert(numelyx_core_cluster_sandwich_covariance(
  bread,scores,invalid_clusters,2,1,2,NUMELYX_COVARIANCE_CLUSTER_CR0,&covariance)
  ==NUMELYX_COVARIANCE_DOMAIN_ERROR);
 assert(covariance==17);
 {
  const double nonfinite[]={NAN,2};
  assert(numelyx_core_sandwich_covariance(
   bread,nonfinite,2,1,NUMELYX_COVARIANCE_HC0,&covariance)
   ==NUMELYX_COVARIANCE_NON_FINITE_INPUT);
  assert(covariance==17);
 }
 assert(numelyx_core_sandwich_covariance(
  bread,scores,1,1,NUMELYX_COVARIANCE_HC1,&covariance)
  ==NUMELYX_COVARIANCE_DOMAIN_ERROR);
 assert(covariance==17);
}
int main(void){
 test_hc0_hc1_sandwich(); test_multivariate_sandwich_orientation(); test_hc2_hc3_regression_covariance(); test_cluster_and_hac_sandwich();
 test_covariance_validation_and_output_preservation();
 puts("Robust covariance tests passed."); return 0;
}
