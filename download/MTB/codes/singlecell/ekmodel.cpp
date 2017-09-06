//Echebarria B, Karma A.
//Mechanisms for initiation of cardiac discordant alternans.
//The European Physical Journal Special Topics 2007; 146: 217‐31.
//http://link.springer.com/article/10.1140/epjst/e2007-00182-y

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

int main(void) {
  //initial values
  float v=0.0;//membrane voltage
  float h=0.8;//Na channel inactivation gate
  float f=0.5;//Ca channel inactivation gate

  //constants
  const float dt=0.1;// time step (0.1 ms)
  const float pcl=200;//pacing cycle length (ms)
  const int itr=20;//the number of beats
  float tmax=pcl*itr;//time to simulate

  int tnmax=tmax/dt;//convert to integer
  int pcln=pcl/dt;//convert to integer
  int durn=1.0/dt;//duration of stimulation (1.0 ms)
  
  //main loop
  for (int tn=0;tn<tnmax;tn++) {

    float stim=0;//stimulation current
    if (tn%pcln < durn) {
      stim=0.3;
    }

    //Euler method
    const float tauso=15;
    const float taufi=0.8;
    const float tauh1=4.8;
    const float tauh2=10.0;
    const float tausi=4.0;
    const float tauf1=100;
    const float tauf2=30;

    float minf=pow((v/0.2),6)/(1+pow((v/0.2),6));
    float hinf=1/(1+pow((v/0.1),6));
    float dinf=pow((v/0.4),4)/(1+pow((v/0.4),4));
    float finf=1/(1+pow((v/0.1),4));

    float tauh=tauh1+tauh2*exp(-20*pow((v-0.1),2));
    float tauf=tauf2+(tauf1-tauf2)*v*v*v;

    float jfi=h*minf*(v-1.3)/taufi;//Fast inward current (Na current)
    float jsi=f*dinf*(v-1.4)/tausi;//Slow inward current (Ca current)
    float jso=(1-exp(-4*v))/tauso;//outward current (K current)
    float ion=-(jfi+jsi+jso-stim);//total transmembrane current

    float dh=(hinf-h)/tauh;
    float df=(finf-f)/tauf;

    //update variables
    v+=ion*dt;
    h+=dh*dt;
    f+=df*dt;

    if (tn%10==0) {//write results every 1 ms
      float t=tn*dt;
      cout<<t<<"\t"<<v<<"\t"<<h<<"\t"<<f<<endl;
    }
  }
  return 0;
}
