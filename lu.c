#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define q 1.6e-19 //in coulomb
#define kt 0.026 //in eV
#define mun 1400 //in cm^2/V-s
#define ni 1.5e10 //in 1/cc
#define Tn 2e-6 //in 1/s
#define MS 3

void gaussian(int n, double x[], double A[][n], double b[]){
    double p[n];
    for (int i = 0; i < n; i++) {
        p[i] = 0;
    }

    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            double factor = A[j][i] / A[i][i];
            b[j] -= factor * b[i];
            for (int k = i; k < n; k++) {
                A[j][k] -= factor * A[i][k];
            }
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i] / A[i][i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= (A[i][j] / A[i][i]) * x[j];
        }
    }

    for (int i = 0; i < n; i++) {
        printf("%lf\n", x[i]);
    }
}

double n(double x, double t) {
    return 0;
}

double Fn(double x, double t) {
    double fn = kt * log(n(x, t) / ni);
    return fn;
}

int main() {
    double t = 10e-6; //in s
    double l = 4e-6; // in m

    double to = t / MS; //in s - incremental time
    double h = l / MS; //in m - incremental length

    double Jf[MS][MS], Jg[MS][MS], Jh[MS][MS], Uf[MS], Ug[MS], Uh[MS], Rf[MS], Rg[MS], Rh[MS] = {0};

    Jf[0][0] = Jg[0][0] = Jh[0][0] = 1;
    Jf[MS - 1][MS - 1] = Jg[MS - 1][MS - 1] = Jh[MS - 1][MS - 1] = 1;

    for (int i = 1; i < MS - 1; i++) {
        Jf[i][i - 1] = Jg[i][i - 1] = 1 / h;
        Jf[i][i] = Jg[i][i] = -1 / h;
        Jh[i][i - 1] = Jh[i][i + 1] = 1 / h / h;
        Jh[i][i] = -2 / h / h;
    }
    double A[][MS] = {{10, -7, 3},
                       {-6, 8, 4},
                       {2, 6, 9}};
    
    double b[MS] = {5, 7, -1};
    double x[MS] = {0, 0, 0};

    gaussian(MS, x, A, b);
    return 0;
}
