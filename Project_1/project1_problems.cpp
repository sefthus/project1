#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "lib.h"
#include "time.h"
#include <chrono>
#include "armadillo"
using namespace arma;


using namespace std;
ofstream ofile;
// Declaring two functions that will be used:

//The analytic solution
double analytic_func(double x){
    return (1-(1-exp(-10))*x - exp(-10*x)); //5 FLOPS
}
// the source term
double source_func(double x){
    return 100*exp(-10*x); //2 FLOPS
}


void problem_1b(int n, char* outfilename){

    double h = 1.0/(n+1); //step size 2 FLOPS

    //diagonals of the tridiagonalmatrix
    double *a = new double[n+2];
    double *b = new double[n+2];
    double *c = new double[n+2];

    // steps
    double *x = new double[n+2];

    // right hand side
    double *f = new double[n+2];
    double *f_tilde= new double[n+2];
    // f_tilde[0]=0

    double *u = new double[n+2]; // analytical
    double *v = new double[n+2]; // numerical solution
    u[0]=0;
    v[0]=0;

    // we dont use i=0,i=n+2 for a,b,c,f,u so it's fine having all in one loop
    for (int i=0;i<n+2;i++){
        x[i] = i*h; //1 FLOP

        a[i] = -1;
        b[i] = 2;
        c[i] = -1;

        f[i] = h*h*source_func(x[i]);   //2
        u[i] = analytic_func(x[i]);
    }
    //these are in the l.notes, but not with trygve
    //Dirch-something Boundary conditions
    c[n]=0;
    a[1]=0;

    clock_t start, finish;
    start = clock();

    //Forward substitution
    // a(i)*v(i-1) + b(i)*v(i) + c(i)*v(i+1) = f_tilde(i)
    f_tilde[0]=0; //unnecesary
    f_tilde[1]=f[1];
    for(int i=2;i<n+1;i++){
        //update RHS
        b[i] = b[i]-a[i]*c[i-1]/b[i-1]; // 3
        f_tilde[i] = f[i] -a[i]*f_tilde[i-1]/b[i-1]; // 3
        // (n-2)*6
    }
    //Backward substitution
    v[n]=f_tilde[n]/b[n]; // 1
    for(int i=n-1;i>=1;i--){
        v[i]=(f_tilde[i]-c[i]*v[i+1])/b[i]; // 3
        // (n-1+1)*3
    }

    finish=clock();
    double timeused = (double)(finish-start)/(CLOCKS_PER_SEC);
    cout <<"time used to run program: "<<timeused<<endl;

    // write to file
    ofile.open(outfilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << n << endl;
    ofile << "       x:             u(x):          v(x):  " << endl;
    for (int i=1;i<n+1;i++) {
        ofile << setw(15) << setprecision(8) << x[i];
        ofile << setw(15) << setprecision(8) << u[i];
        ofile << setw(15) << setprecision(8) << v[i] << endl;
    }

    ofile.close();

    delete [] a;
    delete [] b;
    delete [] c;
    delete [] u;
    delete [] v;
    delete [] x;
    delete [] f;
    delete [] f_tilde;
    //TOTAL FLOPS ( n-2)*6 + n*3 + 11 above
    double TFLOPS = ( n-2)*6 + n*3 + 11+5; // 5 on this line
    cout <<"Total number of FLOP: "<< TFLOPS << endl;
}

void problem_1cd(int n, char* outfilename){
    double h = 1.0/(n+1); //step size 2 FLOPS

    //diagonal of the tridiagonalmatrix
    double *b = new double[n+2];

    // steps
    double *x = new double[n+2];

    // right hand side
    double *f = new double[n+2];
    double *f_tilde= new double[n+2];
    // f_tilde[0]=0

    double *u = new double[n+2]; // analytical
    double *v = new double[n+2]; // numerical solution

    // we dont use i=0,i=n+2 for a,b,c,f,u so it's fine having all in one loop
    for (int i=0;i<n+2;i++){
        x[i] = i*h; //1 FLOP

        b[i] = 2;

        f[i] = h*h*source_func(x[i]);   //2
        u[i] = analytic_func(x[i]);
    }

    clock_t start, finish;
    start = clock();

    //Forward substitution
    f_tilde[1]=f[1];
    for(int i=2;i<n+1;i++){
        //update RHS
        b[i] = (double)(i+1)/i; // 2
        f_tilde[i] = (double)f[i] +(i-1)*f_tilde[i-1]/i; // 4
        //b[i] -= 1/b[i-1];
        //f_tilde[i] = f[i] +f_tilde[i-1]/b[i-1];
        // (n-2)*6
    }
    //Backward substitution
    v[n]=f_tilde[n]/b[n]; // 1
    for(int i=n;i>=1;i--){
        v[i-1]=(double)(i-1)/i*(f_tilde[i-1]+v[i]); // 3
        //v[i] = (f_tilde[i]+v[i+1])/b[i];
        // (n-1+1)*3
    }

    finish=clock();
    double timeused = (double)(finish-start)/(CLOCKS_PER_SEC);
    cout <<"time used to run program: "<<timeused<<endl;

    ofile.open(outfilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << n << endl;
    ofile << "       x:             u(x):          v(x):  " << endl;
    for (int i=1;i<n+1;i++) {
        ofile << setw(15) << setprecision(8) << x[i];
        ofile << setw(15) << setprecision(8) << u[i];
        ofile << setw(15) << setprecision(8) << v[i] << endl;
    }

    ofile.close();

    //TOTAL FLOPS ( n-2)*6 + n*3 + 11 above
    double TFLOPS = ( n-2)*6 + n*3 + 11+5; // 5 on this line
    cout <<"Total number of FLOP: "<< TFLOPS << endl;//at this point in the algorithm

    //problem 1d relative error
    double *eps = new double[n+2];

    for(int i=2;i<n+1;i++){
        eps[i]=log10(fabs((v[i]-u[i])/u[i]));
    }

    double max = eps[2];

    for (int i=2; i<n+1; i++){
        if(fabs(eps[i])>fabs(max))
            max=eps[i];
    }

    //cout <<max<<endl;
    ofile.open("project_1d_error.txt",ofstream::app);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    //ofile << n << endl;
    //ofile << "       n:             h:        max error:  " << endl;
    ofile << setw(25) << setprecision(10) << n;
    ofile << setw(25) << setprecision(10) << h;
    ofile << setw(25) << setprecision(10) << max << endl;

    ofile.close();


    delete [] b;
    delete [] u;
    delete [] v;
    delete [] x;
    delete [] f;
    delete [] f_tilde;

}
void problem_1e(int n){

    mat A = zeros<mat>(n,n);

    double h = 1.0/(n+1); //step size 2 FLOPS

    //diagonal of the tridiagonalmatrix
    vec a(n); a.fill(-1);
    vec b(n); b.fill(2);
    vec c(n); c.fill(-1);

    // steps
    vec x(n);

    // right hand side
    vec f(n);

    // making matrix A
    for (int i=0;i<n;i++){
        if (i>0) A(i,i-1)=a(i);
        A(i,i)=b(i);
        if (i < n-1) A(i,i+1)=c(i);
    }

    // fill in vectors x and f
    for(int i=0;i<n;i++){
        x[i] = i*h;
        f[i] = h*h*source_func(x[i]);
    }

    clock_t start, finish;
    start = clock();

    // solve equation using LU-decomposition Av=f
    vec v = solve(A,f);  // FLOPS goes as n^3

    finish=clock();
    double timeused = (double)(finish-start)/(CLOCKS_PER_SEC);
    cout <<"time used to run program: "<<timeused<<endl;

    cout<<"total number of FLOP:"<<" "<<pow(n,3)<<endl;

    // test LU-decomposition:
    mat L, U, P;
    lu(L,U,P,A);

    //Check that A = LU
    (A-P*L*U).print("test");

}
