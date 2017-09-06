#include <iostream>
#include <fstream>
using namespace std;

__global__ void addition( float* x, float* y, float* z, int num)
{
  const unsigned int tid = threadIdx.x;
  const unsigned int bid = blockIdx.x;
  const unsigned int bdim = blockDim.x;
  const unsigned int gdim = gridDim.x;
  int step=bdim*gdim;

  for (int id=bid * bdim + tid;id<num;id+=step)
  {
    z[id]=x[id]+y[id];
  }
}

int main(void) {
  const int num=100;
  
  //array
  float *h_x;//variable x(host)
  float *h_y;//variable y(host)
  float *h_z;//variable z(host)

  h_x = (float*) malloc(sizeof( float)*num);
  h_y = (float*) malloc(sizeof( float)*num);
  h_z = (float*) malloc(sizeof( float)*num);

  float *d_x;//variable x(device)
  float *d_y;//variable y(device)
  float *d_z;//variable z(device)
  cudaMalloc( (void**) &d_x, sizeof( float) * num );
  cudaMalloc( (void**) &d_y, sizeof( float) * num );
  cudaMalloc( (void**) &d_z, sizeof( float) * num );


  //initial values
  for (int i=0;i<num;i++) {
    h_x[i]=i;
    h_y[i]=i*2;
  }

  //copy x and y from CPU to GPU
  cudaMemcpy( d_x, h_x, sizeof( float) * num , cudaMemcpyHostToDevice);
  cudaMemcpy( d_y, h_y, sizeof( float) * num , cudaMemcpyHostToDevice);

  dim3 grid(256, 1, 1);
  dim3 threads(256, 1, 1);
  //execute GPU function
  addition<<< grid, threads >>>( d_x, d_y, d_z, num);

  //copy z from GPU to CPU
  cudaMemcpy( h_z, d_z, sizeof( float) * num, cudaMemcpyDeviceToHost);
  for (int i=0;i<num;i++) {
    cout<<h_x[i]<<"\t"<<h_y[i]<<"\t"<<h_z[i]<<"\n";
  }

  free(h_x);
  free(h_y);
  free(h_z);
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_z);

  return 0;
}
