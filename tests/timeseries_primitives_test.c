#include "../src/statistics/numelyx_core_timeseries.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
static void close_to(double a,double b,double t){
 if (!(fabs(a-b)<=t)) fprintf(stderr,"timeseries close_to actual=%.17g expected=%.17g tolerance=%.17g\\n",a,b,t);
 assert(fabs(a-b)<=t);
}
static void test_autocovariance_and_acf(void){
 const double series[]={1,2,3,4};
 double value=99;
 double acf[]={-9,-9,-9};
 assert(numelyx_core_autocovariance(series,4,0,false,&value)==NUMELYX_TIMESERIES_SUCCESS);
 close_to(value,1.25,1e-14);
 assert(numelyx_core_autocovariance(series,4,1,false,&value)==NUMELYX_TIMESERIES_SUCCESS);
 close_to(value,.3125,1e-14);
 assert(numelyx_core_autocovariance(series,4,1,true,&value)==NUMELYX_TIMESERIES_SUCCESS);
 close_to(value,5.0/12.0,1e-14);
 assert(numelyx_core_autocorrelation(series,4,1,&value)==NUMELYX_TIMESERIES_SUCCESS);
 close_to(value,.25,1e-14);
 assert(numelyx_core_acf(series,4,2,acf)==NUMELYX_TIMESERIES_SUCCESS);
 close_to(acf[0],1,1e-15);close_to(acf[1],.25,1e-14);close_to(acf[2],-.3,1e-14);
}
static void test_differencing_and_validation(void){
 const double series[]={1,2,3,4};
 const double constant[]={5,5,5};
 double output[]={77,78,79};
 size_t output_count=99;
 double value=88;
 assert(numelyx_core_difference(series,4,1,output,&output_count)==NUMELYX_TIMESERIES_SUCCESS);
 assert(output_count==3);close_to(output[0],1,1e-15);close_to(output[2],1,1e-15);
 output[0]=77;output_count=99;
 assert(numelyx_core_difference(series,4,0,output,&output_count)==NUMELYX_TIMESERIES_DOMAIN_ERROR);
 assert(output[0]==77&&output_count==99);
 assert(numelyx_core_autocorrelation(constant,3,1,&value)==NUMELYX_TIMESERIES_DOMAIN_ERROR);
 assert(value==88);
}
int main(void){
 test_autocovariance_and_acf();test_differencing_and_validation();
 puts("Time-series primitive tests passed.");return 0;
}
