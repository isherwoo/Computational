#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <new>
#include "planet.h"
#include "lib.h"
//#include "time.h"


using namespace std;

int main()
{
    double **R,**X,**Y,**Xhalf,**Yhalf,**Rhalf;
    double **Vxhalf,**Vyhalf,**Vx,**Vy;
   // double *vx,*vy,*x,*y,*r;
    //double *vxhalf,*vyhalf,*xhalf,*yhalf,*rhalf;
    double gravnowx,gravnowy,gravnextx,gravnexty,gravhalfx,gravhalfy,tmax,t0,h,h2,pi,gravnow,gravnext,test,gravhalf;
    double *kx1,*kx2,*kx3,*kx4,*kvx1,*kvx2,*kvx3,*kvx4;
    double *ky1,*ky2,*ky3,*ky4,*kvy1,*kvy2,*kvy3,*kvy4;
    int i,j,n,k,m;
    t0 = 0.0;
    tmax = 3.0;
    n = 20500;
    h = (tmax - t0)/(n+1);
    h2 = h*h;
    pi = acos(-1.0);
  /*  r = new double[n];
    kvx1 = new double[n];
    kvy1 = new double[n];
    kx1 = new double[n];
    ky1 = new double[n];

    rhalf = new double[n];
    vxhalf = new double[n];
    vyhalf = new double[n];
    xhalf = new double[n];
    yhalf = new double[n];*/


        kvx1 = new double[3];
        kvy1 = new double[3];
        kx1 = new double[3];
        ky1 = new double[3];
        kvx2 = new double[3];
        kvy2 = new double[3];
        kx2 = new double[3];
        ky2 = new double[3];
        kvx3 = new double[3];
        kvy3 = new double[3];
        kx3 = new double[3];
        ky3 = new double[3];
        kvx4 = new double[3];
        kvy4 = new double[3];
        kx4 = new double[3];
        ky4 = new double[3];

    R = (double **) matrix(n,3,sizeof(double));
    X = (double **) matrix(n,3,sizeof(double));
    Y = (double **) matrix(n,3,sizeof(double));
    Rhalf = (double **) matrix(n,3,sizeof(double));
    Xhalf = (double **) matrix(n,3,sizeof(double));
    Yhalf = (double **) matrix(n,3,sizeof(double));

    Vx = (double **) matrix(n,3,sizeof(double));
    Vy = (double **) matrix(n,3,sizeof(double));
    Vxhalf = (double **) matrix(n,3,sizeof(double));
    Vyhalf = (double **) matrix(n,3,sizeof(double));




   // planet Earth(0.000003,-9.687087657503362E-01, -2.304478999525212E-01,0.0,3.715293823075145E-03*365.0,-1.679555302972549E-02*365.0,0.0);
   // planet Sun(1.0,0.00510182404,0.0008987727027,0.0,-0.0004829298902,-0.002481983714,0.0); //(mass,x,y,z,vx,vy,vz)
   // planet Jupiter(.000954,-5.344777687381148,.9428344302147893,0.0,-1.398575965353750E-03*365.0,-7.075019624039646E-03*365.0,0.0);
    //planet all_planets[3] = {Sun,Earth,Jupiter};
planet Earth(0.000003,1.0000,0,0,0,2*pi,0);
planet Sun(1.0,0,0,0,0,0,0);
//planet CM(0,0,0,0,0,0,0);
planet Jupiter(0.000954*10,5.2,0,0,0,-2*pi/sqrt(5.2),0);
planet all_planets[3] = {Sun,Earth,Jupiter};
m = 3;


/*x[0] = 1.0;
vy[0] = 2*pi+2.5885754;
vx[0] = 0.0;
y[0] = 0.0;*/

for(i=0;i<m;i++){
X[0][i] = all_planets[i].position[0];
Y[0][i] = all_planets[i].position[1];
Vx[0][i] = all_planets[i].velocity[0];
Vy[0][i] = all_planets[i].velocity[1];

}
  //  cout << X[0][1] << endl;





    //Verlet-------------------------------------------------------------------------------------------------------------

   for(i = 0;i<n-1;i++){
       for(j = 0;j<m;j++){
         //  cout << R[0][1] << endl;
        R[i][j] = sqrt(X[i][j]*X[i][j] + Y[i][j]*Y[i][j]);
        gravnow = 4.0*pi*pi/pow(R[i][j],3);
        kx1[j] = 0.0;
        ky1[j] = 0.0;
           for(k=0;k<m;k++){

            if(k!=j){

               kx1[j] += 4*pi*pi*(X[i][j]-X[i][k])*all_planets[k].mass/pow(pow(X[i][j]-X[i][k],2.0)+pow(Y[i][j]-Y[i][k],2.0),1.5);
               ky1[j] += 4*pi*pi*(Y[i][j]-Y[i][k])*all_planets[k].mass/pow(pow(X[i][j]-X[i][k],2.0)+pow(Y[i][j]-Y[i][k],2.0),1.5);}
               else{
               kx1[j]+= 0.0;
               ky1[j]+= 0.0;
}           }

         //  cout << gravnow*X[i][1] << "  " << kx1[1] << endl;
    X[i+1][j] = X[i][j] + h*Vx[i][j] - h2/2.0*kx1[j];
    Y[i+1][j] = Y[i][j] + h*Vy[i][j] - h2/2.0*ky1[j];
       }
    //all_planets[j].position[0] = X[i+1][j];
    //all_planets[];
       for(j=0;j<m;j++){
           kx2[j] = 0.0;
           ky2[j] = 0.0;
    R[i+1][j] = sqrt(X[i+1][j]*X[i+1][j] + Y[i+1][j]*Y[i+1][j]);
    gravnext = 4.0*pi*pi/(R[i+1][j]*R[i+1][j]*R[i+1][j]);
    for(k=0;k<m;k++){
     if(k!=j){
         //cout<<X[i+1][k] << endl;
        kx2[j] += 4*pi*pi*(X[i+1][j]-X[i+1][k])*all_planets[k].mass/pow(sqrt(pow(X[i+1][j]-X[i+1][k],2.0)+pow(Y[i+1][j]-Y[i+1][k],2.0)),3);
        ky2[j] += 4*pi*pi*(Y[i+1][j]-X[i+1][k])*all_planets[k].mass/pow(sqrt(pow(X[i+1][j]-X[i+1][k],2.0)+pow(Y[i+1][j]-Y[i+1][k],2.0)),3);
        //cout << gravnextx << " " << j << " " << k << endl;
     }else{
        kx2[j]+= 0.0;
        ky2[j]+= 0.0;
}           }

    Vx[i+1][j] = Vx[i][j] -h/2.0*(kx2[j] + kx1[j]);
    test =Vy[i][j] - h/2.0*(gravnext*Y[i+1][1]+gravnow*Y[i][1]);
    Vy[i+1][j] = Vy[i][j] -h/2.0*(ky2[j] + ky1[j]);

   // cout << Vy[i+1][1] << "  "<< test << endl;
}
}
    ofstream myfile;
        myfile.open("proj3_jup_circle_n=6000_t=12_VV_3plans.txt");
    for(i = 0; i<n;i++){
     myfile << i << "  " << X[i][2] << "  " << Y[i][2] << "  "<< Vx[i][2] << "  " << Vy[i][2] << endl;
      //"  " << i << endl;
    }
//-------------------------------------------------------------------------------------------------------------------------

    //Runge-kutta---------------------------------------------------------------------------------------------------------
/*for(i=0;i<n-1;i++){

   //First
    for(j=1;j<m;j++){
    R[i][j] = sqrt(X[i][j]*X[i][j] + Y[i][j]*Y[i][j]);
    gravnowx = 0.0;
    gravnowy = 0.0;
       for(k=0;k<m;k++){

        if(k!=j){

           gravnowx += 4*pi*pi*(X[i][j]-X[i][k])*all_planets[k].mass/pow(pow(X[i][j]-X[i][k],2.0)+pow(Y[i][j]-Y[i][k],2.0),1.5);
           gravnowy += 4*pi*pi*(Y[i][j]-Y[i][k])*all_planets[k].mass/pow(pow(X[i][j]-X[i][k],2.0)+pow(Y[i][j]-Y[i][k],2.0),1.5);}
           else{
           gravnowx+= 0.0;
           gravnowy+= 0.0;
}           }



    gravnow = 4.0*pi*pi/(R[i][j]*R[i][j]*R[i][j]);


    kvx1[j]= gravnowx;
    kvy1[j] = gravnowy;
    Xhalf[i][j] = X[i][j] + h/2.0*Vx[i][j];
    Yhalf[i][j] = Y[i][j] + h/2.0*Vy[i][j];
    Vxhalf[i][j] = Vx[i][j] - kvx1[j]*h/2.0;
    Vyhalf[i][j] = Vy[i][j] - kvy1[j]*h/2.0;
    //test = gravnow*X[i][j];
    //cout << test << "  " << kvx1[1] << endl;
}
    //Second

    for(j=1;j<m;j++){
    R[i][j] = sqrt(X[i][j]*X[i][j] + Y[i][j]*Y[i][j]);
    gravhalfx = 0.0;
    gravhalfy = 0.0;
       for(k=0;k<m;k++){

        if(k!=j){

           gravhalfx += 4*pi*pi*(Xhalf[i][j]-Xhalf[i][k])*all_planets[k].mass/pow(pow(Xhalf[i][j]-Xhalf[i][k],2.0)+pow(Yhalf[i][j]-Yhalf[i][k],2.0),1.5);
           gravhalfy += 4*pi*pi*(Yhalf[i][j]-Yhalf[i][k])*all_planets[k].mass/pow(pow(Xhalf[i][j]-Xhalf[i][k],2.0)+pow(Yhalf[i][j]-Yhalf[i][k],2.0),1.5);}
           else{
           gravhalfx+= 0.0;
           gravhalfy+= 0.0;
}           }

        Rhalf[i][j] = sqrt(Xhalf[i][j]*Xhalf[i][j] +Yhalf[i][j]*Yhalf[i][j] );
        gravhalf = 4.0*pi*pi/(Rhalf[i][j]*Rhalf[i][j]*Rhalf[i][j]);
        kvx2[j] = gravhalfx;
        kvy2[j] = gravhalfy;
        kx2[j] = Vxhalf[i][j];
        ky2[j] = Vyhalf[i][j];
        //test = gravhalf*Yhalf[i][1];
        //cout << test << "   " << kvy2[1] << endl;
}
    //Third
    for(j=1;j<m;j++){
        Xhalf[i][j] = X[i][j] + h/2.0*kx2[j];
        Yhalf[i][j] = Y[i][j] + h/2.0*ky2[j];
        Vxhalf[i][j] = Vx[i][j] - h/2.0*kvx2[j];
        Vyhalf[i][j] = Vy[i][j] - h/2.0*kvy2[j];}

    //Fourth

    for(j=1;j<m;j++){
    Rhalf[i][j] = sqrt(Xhalf[i][j]*Xhalf[i][j] + Yhalf[i][j]*Yhalf[i][j]);
    gravhalfx = 0.0;
    gravhalfy = 0.0;
       for(k=0;k<m;k++){

        if(k!=j){

           gravhalfx += 4*pi*pi*(Xhalf[i][j]-Xhalf[i][k])*all_planets[k].mass/pow(pow(Xhalf[i][j]-Xhalf[i][k],2.0)+pow(Yhalf[i][j]-Yhalf[i][k],2.0),1.5);
           gravhalfy += 4*pi*pi*(Yhalf[i][j]-Yhalf[i][k])*all_planets[k].mass/pow(pow(Xhalf[i][j]-Xhalf[i][k],2.0)+pow(Yhalf[i][j]-Yhalf[i][k],2.0),1.5);}
           else{
           gravhalfx+= 0.0;
           gravhalfy+= 0.0;
}           }






        gravhalf = 4.0*pi*pi/(Rhalf[i][1]*Rhalf[i][1]*Rhalf[i][1]);
        kvx3[j] = gravhalfx;
        kvy3[j] = gravhalfy;
        kx3[j] = Vxhalf[i][j];
        ky3[j] = Vyhalf[i][j];
        //test = gravhalf*X[i][1];
        //cout << test << "  " << kvx3[1] << endl;


    }


    //Fifth

    for(j=1;j<m;j++){
    Xhalf[i][j] = X[i][j] + h*kx3[j];
        Yhalf[i][j] = Y[i][j] + h*ky3[j];
        Vxhalf[i][j] = Vx[i][j] - h*kvx3[j];
        Vyhalf[i][j] = Vy[i][j] - h*kvy3[j];
        Rhalf[i][j] = sqrt(Xhalf[i][j]*Xhalf[i][j] +Yhalf[i][j]*Yhalf[i][j] );
                gravhalf = 4.0*pi*pi/(Rhalf[i][j]*Rhalf[i][j]*Rhalf[i][j]);


                gravhalfx = 0.0;
                gravhalfy = 0.0;
                   for(k=0;k<m;k++){

                    if(k!=j){

                       gravhalfx += 4*pi*pi*(Xhalf[i][j]-Xhalf[i][k])*all_planets[k].mass/pow(pow(Xhalf[i][j]-Xhalf[i][k],2.0)+pow(Yhalf[i][j]-Yhalf[i][k],2.0),1.5);
                       gravhalfy += 4*pi*pi*(Yhalf[i][j]-Yhalf[i][k])*all_planets[k].mass/pow(pow(Xhalf[i][j]-Xhalf[i][k],2.0)+pow(Yhalf[i][j]-Yhalf[i][k],2.0),1.5);}
                       else{
                       gravhalfx+= 0.0;
                       gravhalfy+= 0.0;
            }           }




                kvx4[j] = gravhalfx;
                kvy4[j] = gravhalfy;
                kx4[j] = Vxhalf[i][j];
                ky4[j] = Vyhalf[i][j];
                //test = gravhalf*X[i][j];
                //cout <<test << "  " << kvx4[1] << endl;
}

    //Sixth
    for(j=1;j<m;j++){
        X[i+1][j] = X[i][j] + h/6.0*(Vx[i][j] +2.0*(kx2[j]+kx3[j]) + kx4[j]);
        Y[i+1][j] = Y[i][j] + h/6.0*(Vy[i][j] +2.0*(ky2[j]+ky3[j]) + ky4[j]);
        Vx[i+1][j] = Vx[i][j] - h/6.0*(kvx1[j] +2.0*(kvx2[j]+kvx3[j]) + kvx4[j]);
        Vy[i+1][j] = Vy[i][j] - h/6.0*(kvy1[j] +2.0*(kvy2[j]+kvy3[j]) + kvy4[j]);
        R[i+1][j] = sqrt(X[i+1][j]*X[i+1][j] + Y[i+1][j]*Y[i+1][j]);

}

    //cout << r[i] << endl;

//cout << R[i][1] << endl;

}

    ofstream myfile;
        myfile.open("proj3_earth1000_circle_n=6000_t=12_RK4_2plans.txt");
    for(i = 0; i<n;i++){
   cout << X[i][1] << "  " << Y[i][1]
                   << endl;
      //"  " << i << endl;
    }*/
//cout << r[n-1] << " " << x[n-1] << " " << y[n-1] << " " << vx[n-1] << " " << vy[n-1] - 2*pi << endl;

    //---------------------------------------------------------------------------------------------------------------------
    return 0;
}
