#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

/*int matrix()
{
    int i;
    int k;
    int j;
double A[5][5];
j = 5;



 for (i=0; i<j; i++){
     for (k=0; k<j; k++){
         if(i == k)
             A[i][k] = 2.0;
         else if (k == i-1)
             A[i][k] = -1.0;
         else if (k == i+1)
             A[i][k] = -1.0;
         else
             A[i][k] = 0.0;


     }

 }


for (i = 0; i<j; i++){

    cout << A[i][0] <<" " << A[i][1] <<" " << A[i][2] <<" " << A[i][3] << " " << A[i][4] << endl;
}



    return 0;



}*/

int vectorb(){
    int n;

    cout << "Give step number" << endl;
    cin >> n;
    int i;
    double h = 1.0/(n+1.0);
    double *bt,*temp, *v;
    double *a,*b,*c,*u,*err;
    a = new double[n];
    b = new double[n];
    c = new double[n];
    u = new double[n];
    err = new double[n];
    bt = new double[n];
    temp = new double[n];
    v = new double[n];
    double * A;
    double btemp,atemp = 0.0 ;



    for(i=0; i<n ; i++){
        bt[i] = h*h*100.0*exp(-10*i*h);
        b[i] = 2.0;
        a[i] = -1.0;
        c[i] = -1.0;
    //cout << bt[i] << endl;
}
    a[0] = 0.0;
    c[n] = 0.0;

    btemp = b[0];
    u[0] = bt[0]/btemp;
    for(i=1 ; i<n ; i++){
        temp[i] = c[i-1]/btemp;
        btemp = b[i]-a[i]*temp[i];
        v[i] = 1-(1-exp(-10.0))*(i+1)*h-exp(-10.0*(i+1)*h);
        u[i] = (bt[i] - a[i]*u[i-1])/btemp;
        // cout << u[i] << endl;
    }

    for(i = n-1 ; i>= 1; i--){
        u[i] -= temp[i+1]*u[i+1];
    }

    ofstream myfile;
    myfile.open("proj1.txt");

    for(i=0; i< n ; i++){
        err[i] = log10(abs((u[i] - v[i])/v[i] ));
        myfile << u[i] << "\n";
        //if(err[i]<atemp){
          //  atemp = err[i];
        cout << u[i] << endl;
//}
    }

    delete [] a;
    delete [] b;
    delete [] c;
    delete [] bt;
    delete [] err;
    delete [] temp;
    delete [] v;



    return 0;



}


int main(){
    vectorb();

    return 0;



}

