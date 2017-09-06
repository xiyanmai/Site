#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

int main(void) {
  const int num=100;
  
  //array
  float v[num];//membrane voltage
  float h[num];//Na channel inactivation gate
  float f[num];//Ca channel inactivation gate
  float stim[num];
  float tmp[num];//temporary variable to solve diffusion eqn

  //initial values
  for (int i=0;i<num;i++) {
    v[i]=0;
    h[i]=0.8;
    f[i]=0.5;
    stim[i]=0;
  }

  //constants
  const float dt=0.1;// time step (0.1 ms)
  const float pcl=260;//pacing cycle length (ms)
  const int itr=5;//# of beats
  float tmax=pcl*itr;//time to simulate

  int tnmax=tmax/dt;//convert to int
  int pcln=pcl/dt;//convert to int
  int durn=1.0/dt;//duration of stimulation (1.0 ms)
  
  //main loop
  for (int tn=0;tn<tnmax;tn++) {
    if (tn%10==0) {//write results every 1 ms
      float t=tn*dt;
      for (int i=0;i<num;i++) {
        cout<<v[i]<<"\t";
      }
      cout<<endl;
    }

    //stimlate the end of the cable
    if (tn%pcln < durn) {
      stim[0]=stim[1]=stim[2]=stim[3]=stim[4]=0.3;
    }
    else {
      stim[0]=stim[1]=stim[2]=stim[3]=stim[4]=0;
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
    const float dfu=0.001;//diffusion coefficient
    const float dx=0.15;//0.15 mm
    //non-flux boundary
    v[0]=v[2];
    v[num-1]=v[num-3];
    for (int i=1;i<num-1;i++) {
      tmp[i]=v[i]+(v[i-1]+v[i+1]-2*v[i])*dfu*dt/(dx*dx);
    }
    for (int i=1;i<num-1;i++) {
      v[i]=tmp[i];
    }
  }
  return 0;
}
