#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main(void) {
  const int xnum=600;
  const int ynum=600;
  const int num=xnum*ynum;
  
  //array
  float *v=new float [num];//membrane voltage
  float *h=new float [num];//Na channel inactivation gate
  float *f=new float [num];//Ca channel inactivation gate
  float *stim=new float [num];
  float *tmp=new float [num];//temporary variable to solve diffusion eqn

  //initial values
  for (int i=0;i<num;i++) {
    v[i]=0;
    h[i]=0.8;
    f[i]=0.5;
    stim[i]=0;
  }
  //stimlate the corner
  for(int i=0;i<10;i++) {
    for(int j=0;j<10;j++) {
      v[i*ynum+j]=1;
    }
  }

  //constants
  const float dt=0.1;// time step (0.1 ms)
  float tmax=200;//total simulation time

  //convert to int
  int tnmax=tmax/dt;
  
  //main loop
  for (int tn=0;tn<tnmax;tn++) {
    
    //write results every 10 ms
    if (tn%100==0) {
      char fname[255];
      sprintf(fname,"v%d.txt",tn);
      ofstream os(fname);
      for(int i=0;i<xnum;i+=6) {
        for(int j=0;j<ynum;j+=6) {
          os<<v[i*ynum+j]<<"\t";
        }
        os<<endl;
      }
    }

    //Euler method
    for (int i=0;i<num;i++) {
      const float tauso=15;
      const float taufi=0.8;
      const float tauh1=4.8;
      const float tauh2=10.0;
      const float tausi=4.0;
      const float tauf1=100;
      const float tauf2=30;
      float minf=pow((v[i]/0.2),6)/(1+pow((v[i]/0.2),6));
      float hinf=1/(1+pow((v[i]/0.1),6));
      float dinf=pow((v[i]/0.4),4)/(1+pow((v[i]/0.4),4));
      float finf=1/(1+pow((v[i]/0.1),4));
      float tauh=tauh1+tauh2*exp(-20*pow((v[i]-0.1),2));
      float tauf=tauf2+(tauf1-tauf2)*v[i]*v[i]*v[i];

      float jfi=h[i]*minf*(v[i]-1.3)/taufi;//Fast inward current (Na current)
      float jsi=f[i]*dinf*(v[i]-1.4)/tausi;//Slow inward current (Ca current)
      float jso=(1-exp(-4*v[i]))/tauso;//outward current (K current)
      float ion=-(jfi+jsi+jso-stim[i]);//total transmembrane current

      float dh=(hinf-h[i])/tauh;
      float df=(finf-f[i])/tauf;

      //update variables
      v[i]+=ion*dt;
      h[i]+=dh*dt;
      f[i]+=df*dt;

    }
    //solve Diffusion
    const float dfu=0.001;
    const float dx=0.015;//0.15 mm
    //non-flux boundary
    for (int i=0;i<xnum;i++) {
      v[i*ynum+0]=v[i*ynum+2];
      v[i*ynum+ynum-1]=v[i*ynum+ynum-3];
    }
    for (int j=0;j<ynum;j++) {
      v[0*ynum+j]=v[2*ynum+j];
      v[(xnum-1)*ynum+j]=v[(xnum-3)*ynum+j];
    }
    for (int i=1;i<xnum-1;i++) {
      for (int j=1;j<ynum-1;j++) {
        tmp[i*ynum+j]=v[i*ynum+j]+(v[(i-1)*ynum+j]+v[(i+1)*ynum+j]+v[i*ynum+(j-1)]+v[i*ynum+(j+1)]-4*v[i*ynum+j])*dfu*dt/(dx*dx)/2;
      }
    }

    for (int i=0;i<xnum;i++) {
      tmp[i*ynum+0]=tmp[i*ynum+2];
      tmp[i*ynum+ynum-1]=tmp[i*ynum+ynum-3];
    }
    for (int j=0;j<ynum;j++) {
      tmp[0*ynum+j]=tmp[2*ynum+j];
      tmp[(xnum-1)*ynum+j]=tmp[(xnum-3)*ynum+j];
    }
    for (int i=1;i<xnum-1;i++) {
      for (int j=1;j<ynum-1;j++) {
        v[i*ynum+j]=tmp[i*ynum+j]+(tmp[(i-1)*ynum+j]+tmp[(i+1)*ynum+j]+tmp[i*ynum+(j-1)]+tmp[i*ynum+(j+1)]-4*tmp[i*ynum+j])*dfu*dt/(dx*dx)/2;
      }
    }
  }
  
  delete[] v;
  delete[] h;
  delete[] f;
  delete[] stim;
  delete[] tmp;

  return 0;
}
