//Echebarria B, Karma A.
//Mechanisms for initiation of cardiac discordant alternans.
//The European Physical Journal Special Topics 2007; 146: 217Å]31.
//http://link.springer.com/article/10.1140/epjst/e2007-00182-y
//

#include <iostream>
#include <fstream>
using namespace std;


__global__ void odeKernel( float* v, float* h, float* f, float* stim, int num)
{
  __shared__ float vsm[256];

  const unsigned int tid = threadIdx.x;
  const unsigned int bid = blockIdx.x;
  const unsigned int bdim = blockDim.x;
  const unsigned int gdim = gridDim.x;
  int step=bdim*gdim;

  for (int id=bid * bdim + tid;id<num;id+=step)
  {
    vsm[tid]=v[id];
    const float dt=0.1;

    //Euler method
    const float tauso=15;
    const float taufi=0.8;
    const float tauh1=4.8;
    const float tauh2=10.0;
    const float tausi=4.0;
    const float tauf1=100;
    const float tauf2=30;

    float minf=pow((vsm[tid]/0.2),6)/(1+pow((vsm[tid]/0.2),6));
    float hinf=1/(1+pow((vsm[tid]/0.1),6));
    float dinf=pow((vsm[tid]/0.4),4)/(1+pow((vsm[tid]/0.4),4));
    float finf=1/(1+pow((vsm[tid]/0.1),4));
    float tauh=tauh1+tauh2*exp(-20*pow((vsm[tid]-0.1),2));
    float tauf=tauf2+(tauf1-tauf2)*vsm[tid]*vsm[tid]*vsm[tid];

    float jfi=h[id]*minf*(vsm[tid]-1.3)/taufi;//Fast inward current (Na current)
    float jsi=f[id]*dinf*(vsm[tid]-1.4)/tausi;//Slow inward current (Ca current)
    float jso=(1-exp(-4*vsm[tid]))/tauso;//outward current (K current)
    float ion=-(jfi+jsi+jso-stim[id]);//total transmembrane current

    float dh=(hinf-h[id])/tauh;
    float df=(finf-f[id])/tauf;

    //update variables
    v[id]+=ion*dt;
    h[id]+=dh*dt;
    f[id]+=df*dt;
  }
}

__global__ void diffKernel( float* v, float* vold ,int num)
{
  const unsigned int tid = threadIdx.x;
  const unsigned int bid = blockIdx.x;
  const unsigned int bdim = blockDim.x;
  const unsigned int gdim = gridDim.x;
  int step=bdim*gdim;

  const float dt=0.1;
  const float dfu=0.0005;
  const float dx=0.015;
  for (int id=bid * bdim + tid;id<num;id+=step)
  {
    if (id==0)
      v[id]=vold[id]+(vold[id+1]+vold[id+1]-2*vold[id])*dfu*dt/(dx*dx)/2;
    else if (id==num-1)
      v[id]=vold[id]+(vold[id-1]+vold[id-1]-2*vold[id])*dfu*dt/(dx*dx)/2;
    else
      v[id]=vold[id]+(vold[id-1]+vold[id+1]-2*vold[id])*dfu*dt/(dx*dx)/2;
  }
}


int main(void) {
  const int num=400;
  
  //array
  float *h_v;//membrane voltage
  float *h_h;//Na channel inactivation gate
  float *h_f;//Ca channel inactivation gate
  float *h_stim;
  float *h_tmp;//temporary variable to solve diffusion eqn

  h_v = (float*) malloc(sizeof( float)*num);
  h_h = (float*) malloc(sizeof( float)*num);
  h_f = (float*) malloc(sizeof( float)*num);
  h_stim = (float*) malloc(sizeof( float)*num);
  h_tmp = (float*) malloc(sizeof( float)*num);

  float *d_v;//membrane voltage
  float *d_h;//Na channel inactivation gate
  float *d_f;//Ca channel inactivation gate
  float *d_stim;
  float *d_tmp;//temporary variable to solve diffusion eqn
  cudaMalloc( (void**) &d_v, sizeof( float) * num );
  cudaMalloc( (void**) &d_h, sizeof( float) * num );
  cudaMalloc( (void**) &d_f, sizeof( float) * num );
  cudaMalloc( (void**) &d_stim, sizeof( float) * num );
  cudaMalloc( (void**) &d_tmp, sizeof( float) * num );


  //initial values
  for (int i=0;i<num;i++) {
    h_v[i]=0;
    h_h[i]=0.8;
    h_f[i]=0.5;
    h_stim[i]=0;
  }

  cudaMemcpy( d_v, h_v, sizeof( float) * num , cudaMemcpyHostToDevice);
  cudaMemcpy( d_h, h_h, sizeof( float) * num , cudaMemcpyHostToDevice);
  cudaMemcpy( d_f, h_f, sizeof( float) * num , cudaMemcpyHostToDevice);
  cudaMemcpy( d_stim, h_stim, sizeof( float) * num , cudaMemcpyHostToDevice);



  //constants
  const float dt=0.1;// time step (0.1 ms)
  const float pcl=140;//pacing cycle length (ms)
  const int itr=10;//# of beats
  float tmax=pcl*itr;//time to simulate

  int tnmax=tmax/dt;//convert to int
  int pcln=pcl/dt;//convert to int
  int durn=1.0/dt;//duration of stimulation (1.0 ms)

  dim3 grid( 256, 1, 1);
  dim3 threads(256, 1, 1);

  
  //main loop
  for (int tn=0;tn<tnmax;tn++) {
    if (tn%100==0) {//write results every 10 ms
      //copy from GPU to CPU
      cudaMemcpy( h_v, d_v, sizeof( float) * num, cudaMemcpyDeviceToHost);
      for (int i=0;i<num;i++) {
        cout<<h_v[i]<<"\t";
      }
      cout<<endl;
    }

    //stimlate the end of the cable
    if (tn%pcln == 0) {
      h_stim[0]=h_stim[1]=h_stim[2]=h_stim[3]=h_stim[4]=0.3;
      cudaMemcpy( d_stim, h_stim, sizeof( float) * num , cudaMemcpyHostToDevice);
    }
    else if (tn%pcln == durn) {
      h_stim[0]=h_stim[1]=h_stim[2]=h_stim[3]=h_stim[4]=0;
      cudaMemcpy( d_stim, h_stim, sizeof( float) * num , cudaMemcpyHostToDevice);
    }
    //solve ODE
    odeKernel<<< grid, threads >>>( d_v, d_h, d_f, d_stim,num);
    cudaThreadSynchronize();
    //solve Diffusion
    diffKernel<<< grid, threads >>>( d_tmp, d_v, num);
    cudaThreadSynchronize();
    diffKernel<<< grid, threads >>>( d_v, d_tmp, num);
    cudaThreadSynchronize();
  }

  free(h_v);
  free(h_h);
  free(h_f);
  free(h_stim);
  free(h_tmp);
  cudaFree(d_v);
  cudaFree(d_h);
  cudaFree(d_f);
  cudaFree(d_stim);
  cudaFree(d_tmp);

  return 0;
}
