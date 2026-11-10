#include "../src/statistics/numelyx_core_gmm.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
static void close_to(double a,double b,double t){assert(fabs(a-b)<=t);}
static void test_gmm_objective_gradient_and_curvature(void){
 const double moments[]={1,2};
 const double weighting[]={2,0,0,1};
 const double jacobian[]={1,2};
 double objective=99,gradient=99,hessian=99;
 assert(numelyx_core_gmm_evaluate(moments,weighting,jacobian,2,1,
  &objective,&gradient,&hessian)==NUMELYX_INFERENCE_SUCCESS);
 close_to(objective,6,1e-14);
 close_to(gradient,12,1e-14);
 close_to(hessian,12,1e-14);
}
static void test_gmm_multivariate_and_output_preservation(void){
 const double moments[]={1,2};
 const double weighting[]={1,.5,.5,2};
 const double jacobian[]={1,0,0,1};
 double objective=77,gradient[]={78,79},hessian[]={80,81,82,83};
 const double nonsymmetric[]={1,1,0,1};
 assert(numelyx_core_gmm_evaluate(moments,weighting,jacobian,2,2,
  &objective,gradient,hessian)==NUMELYX_INFERENCE_SUCCESS);
 close_to(objective,11,1e-14);
 close_to(gradient[0],4,1e-14);
 close_to(gradient[1],9,1e-14);
 close_to(hessian[0],2,1e-14);
 close_to(hessian[1],1,1e-14);
 close_to(hessian[2],1,1e-14);
 close_to(hessian[3],4,1e-14);
 objective=77;gradient[0]=78;gradient[1]=79;
 assert(numelyx_core_gmm_evaluate(moments,nonsymmetric,jacobian,2,2,
  &objective,gradient,hessian)==NUMELYX_INFERENCE_INVALID_ARGUMENT);
 assert(objective==77&&gradient[0]==78&&gradient[1]==79);
 {
  const double one_moment[]={1};
  const double negative_weight[]={-1};
  const double one_jacobian[]={1};
  double one_objective=17,one_gradient=18,one_hessian=19;
  assert(numelyx_core_gmm_evaluate(one_moment,negative_weight,one_jacobian,1,1,
   &one_objective,&one_gradient,&one_hessian)==NUMELYX_INFERENCE_NUMERICAL_FAILURE);
  assert(one_objective==17&&one_gradient==18&&one_hessian==19);
 }
}
int main(void){
 test_gmm_objective_gradient_and_curvature();
 test_gmm_multivariate_and_output_preservation();
 puts("GMM primitive tests passed.");return 0;
}
