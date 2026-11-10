#include "../src/optimization/numelyx_core_box_qp.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
static void close_to(double a,double b,double t){assert(fabs(a-b)<=t);}
static void test_box_qp_interior_and_active_bounds(void){
 const double h[]={2,0,0,4},g[]={-4,-8},lo[]={0,0},hi[]={1,3},x0[]={0,0};
 numelyx_core_box_qp_options o; double x[]={-7,-7},f=-7; size_t it=0;
 assert(numelyx_core_box_qp_default_options(&o)==NUMELYX_QP_CONVERGED);
 assert(numelyx_core_box_qp_solve(h,g,lo,hi,x0,2,&o,x,&f,&it)==NUMELYX_QP_CONVERGED);
 close_to(x[0],1,2e-7); close_to(x[1],2,2e-7); close_to(f,-10,2e-6); assert(it>0);
}
static void test_box_qp_linear_and_unbounded_cases(void){
 const double h[]={0},g[]={-2},lo[]={0},hi[]={5},x0[]={1},uhi[]={INFINITY};
 numelyx_core_box_qp_options o; double x=0,f=0;
 assert(numelyx_core_box_qp_default_options(&o)==NUMELYX_QP_CONVERGED);
 assert(numelyx_core_box_qp_solve(h,g,lo,hi,x0,1,&o,&x,&f,NULL)==NUMELYX_QP_CONVERGED);
 assert(x==5); close_to(f,-10,1e-14);
 assert(numelyx_core_box_qp_solve(h,g,lo,uhi,x0,1,&o,&x,&f,NULL)==NUMELYX_QP_UNBOUNDED);
}
static void test_box_qp_semidefinite_and_validation(void){
 const double h[]={1,1,1,1},g[]={-2,-2},lo[]={0,0},hi[]={2,2},x0[]={0,0};
 const double bad[]={1,2,2,1};
 numelyx_core_box_qp_options o; double x[]={19,20},f=21;
 assert(numelyx_core_box_qp_default_options(&o)==NUMELYX_QP_CONVERGED);
 assert(numelyx_core_box_qp_solve(h,g,lo,hi,x0,2,&o,x,&f,NULL)==NUMELYX_QP_CONVERGED);
 close_to(x[0]+x[1],2,2e-6);
 x[0]=19;x[1]=20;f=21;
 assert(numelyx_core_box_qp_solve(bad,g,lo,hi,x0,2,&o,x,&f,NULL)==NUMELYX_QP_NOT_POSITIVE_SEMIDEFINITE);
 assert(x[0]==19&&x[1]==20&&f==21);
 assert(numelyx_core_box_qp_solve(h,g,lo,hi,x0,2,&o,NULL,&f,NULL)==NUMELYX_QP_INVALID_ARGUMENT);
}
int main(void){
 test_box_qp_interior_and_active_bounds();
 test_box_qp_linear_and_unbounded_cases();
 test_box_qp_semidefinite_and_validation();
 puts("Box-constrained QP tests passed."); return 0;
}
