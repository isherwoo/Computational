#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <new>
#include "lib.h"
#include "time.h"

using namespace std;



int main()
{
    //double *p;
    double **A,**V, **R;
    int i,j,k,l,n,jt,it, count = 0;
    double akl,akk,all,aik,ail,rik,ril,tempu,tempb,max = 1.0, t,tau,c,s;
    double rmin = 0.0,rmax,h,h2,c2,s2,omeg,omeg2,on,off,r0,r1,r2;
    double *B,*e,*U;
    string answer;
    cout<< "Please give a stepsize" << endl;
    cin >> n ;
    cout << "Please give a max radial distance" << endl;
    cin << rmax;
    cout << "Would you like to turn on electrostatic forces?" << "YES | NO" << endl;
    cin >> answer;
    if(answer == "yes" || answer == "Yes" || answer == "YES"){
        cout << "Please enter an oscillation frequency" << endl;
        cin >> omeg;
    }

    omeg2 = omeg*omeg;




            U = new double[n];
            e = new double[n];
            B = new double[n];
            A = (double **) matrix(n,n,sizeof(double));
            V = (double **) matrix(n,n,sizeof(double));
            R = (double **) matrix(n,n,sizeof(double));



    h = (rmax-rmin)/(n+1);
    h2 = h*h;
    on = 2.0/h2;
    off = -1.0/h2;




    for (i=0; i<n; i++){
        for (k=0; k<n; k++){
            if(i == k){
                if(answer == "yes" || answer == "Yes" || answer == "YES"){
                A[i][k] = on+omeg2*(i+1.0)*(i+1.0)*h2 - 1.0/((i+1)*h);}
                else{
                    A[i][k] = on+(i+1.0)*(i+1.0)*h2;
                }
              //  V[i][k] = (i+1.0)*(i+1.0)*h2;
                R[i][k] = 1.0;
                B[i] = A[i][k];
                e[i] = off;}
            else if (k == i-1){
                A[i][k] = off;
                R[i][k] = 0.0;}
            else if (k == i+1){
                A[i][k] = off;
                R[i][k] = 0.0;}
            else{
                A[i][k] = 0.0;
                R[i][k] = 0.0;}


       //  A[i][k] = A[i][k] + V[i][k];

        }
       }
    free_matrix((void **) V);


    //tqli(B,e,n,R);


  //  for (i = 0; i<n; i++){
    //    cout<< A[i][0] << " "<<A[i][1]<<"  "<<A[i][2] <<endl;}

   while(fabs(max) >= 0.0000000001){
   // for(k = 0; k<2;k++){
       max = 0.0;

    for(i = 0; i<n;i++){
        for(j = i+1; j<n; j++)
            if(fabs(A[i][j]) >= fabs(max)){
                max = A[i][j]; it = i; jt = j; }}

akl = 2.0*A[it][jt];
    if (A[it][jt] != 0){
        tau = (A[jt][jt]-A[it][it])/akl;
    if(tau >= 0){
        t = -tau + sqrt(1+tau*tau);}
    else{
        t = -tau - sqrt(1+tau*tau);}


    c = 1/sqrt(1+t*t);
    s = c*t;
    //cout << max<< " " << fabs(max) << endl;
    }
    else{
        c = 1.0;
        s = 0.0;
    }


    c2 = c*c;
    s2 = s*s;
    akk = A[it][it];
    all = A[jt][jt];
    A[it][it] = c2*akk-c*s*akl+s2*all;
    A[jt][jt] = s2*akk + c*s*akl+c2*all;
    A[it][jt] = 0.0;
    A[jt][it] = 0.0;
    for(i = 0;i<n;i++){
        if(i != it && i != jt){
            aik = A[i][it];
            ail = A[i][jt];
            A[i][it] = c*aik -s*ail;
            A[it][i] = A[i][it];
            A[i][jt] = c*ail + s*aik;
            A[jt][i] = A[i][jt];
            //cout << aik << endl;

        }
 rik = R[i][it];
 ril = R[i][jt];

 R[i][it] = c*rik - s*ril;
 R[i][jt] = c*ril + s*rik;
}

   // cout<< fabs(max)<< " " << it << " " << jt <<endl;
    count += 1;
}








   // akl = off_diag(A,n); //returns pointer array whose elements are the largest off diagonal and its respective coordinates
    //it = p[1];
   // jt = p[2];

for (i = 0; i<n; i++){
    U[i] = A[i][i];
   // cout << R[i][39] << " " << i << endl;}
}





//Sort Eigenvalues and vectors
    for (i = 0; i<n; i++){
        for(j=0;j<n-1;j++){
            if(U[j]>U[j+1]){
                tempu = U[j+1];
                U[j+1] = U[j];
                U[j] = tempu;
                for(k=0;k<n;k++){
                 tempu = R[k][j];
                  R[k][j] = R[k][j+1];
                  R[k][j+1] = tempu;
               }
            }
            if(B[j]>B[j+1]){
                tempb = B[j+1];
                B[j+1] = B[j];
                B[j] = tempb;
            }

        }}
   // cout<< B[0] << " " << B[1] << " " << B[2] << endl;
   // cout << A[i][i]<< endl;


    ofstream myfile;
        myfile.open("proj2d.txt");

        for( i = 0; i<n; i++){
            r0 += R[i][0]*R[i][0];
            r1 += R[i][1]*R[i][1];
            r2 += R[i][2]*R[i][2];

        }
       /* for(i=0; i< n ; i++){

            myfile << R[i][0]/sqrt(r0) << " "<<R[i][1]/sqrt(r1) <<" " <<R[i][2]/sqrt(r2) << "\n";



    }*/

for (i = 0;i<n;i++){
    tempu += R[i][0]*R[i][n-1];
}



cout<< "n = " << n << "  Rmax = "<< rmax << endl;
cout << "Groundstate eigenvalue is " <<U[0]<< endl;
cout << "First excited state eigenvalue is " << U[1] << endl;
cout << "Second excited state eigenvalue is " << U[2] << endl;
cout << "This calculation required " << count - 1 << " Jacobi Rotations" << endl;

delete [] U;
delete [] B;
delete [] e;
free_matrix((void **) R);
free_matrix((void **) A);








    return 0;

}



