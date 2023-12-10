#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define q 1.6e-19 //in coulomb
#define kt 0.026 //in eV
#define mun 1400 //in cm^2/V-s
#define ni 1.5e10 //in 1/cc
#define Tn 2e-2 //in 1/s
#define MS 200
#define T 1e-8 //in s
#define l 6e-6 // in cm
#define no 1.7e10 //in 1/cc

double X[MS][2], Y[MS];

double to = T/MS; //in s - incremental time
double h = l/MS; //in cm - incremental length

void WriteToFile(int N, double y[], double Y[][2], FILE *fptr, int c){
    fptr = fopen("Continuity.txt", "a");
    if(c==0){
        for (int i = 0; i < N; i++) {
        fprintf(fptr, "%e ", y[i]);
        }
        fprintf(fptr, "\n");
    }
    
    for (int i = 0; i < N; i++) {
        fprintf(fptr, "%e ", Y[i][0]);
    }
    fprintf(fptr, "\n");

    fclose(fptr);
}

void gaussian(int n,  double x[], double A[][n], double b[]){
    int count = 0;
    double p = 0;
    double temp[n], a[n][n];
    for(int i=0; i<n; i++){
        for(int j=0; j<MS; j++){
            a[i][j] = A[i][j];
        }
        temp[i] = b[i];
    }
    for(int i=0; i<n-1; i++){
        for(int j=i+1; j<n; j++){
            count = 0;
            if(count == 0)
                p = a[j][i];
            temp[j] = (a[i][i]*temp[j] - p*temp[i])/a[i][i];
            count = count++;
            for(int k=0; k<n; k++){
                a[j][k] = (a[i][i]*a[j][k] - p*a[i][k])/a[i][i];
            }
        }
    }
    for(int i=0; i<n; i++){
        x[i] = 0;
    }
    for (int i = n - 1; i >= 0; i--) {
        x[i] = temp[i] / a[i][i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= (a[i][j] / a[i][i]) * x[j];
        }
    }
    // for(int i=0; i<MS; i++){
    //    for(int j=0; j<MS; j++){
    //         printf("%e ",a[i][j]);
    //     }
    //     printf("\n");
    // }
    // for(int i=0; i<MS; i++){
    //     printf("%e %e %e\n",b[i], X[i][0], x[i]);
    // }
}

double n(double x) {
    int i = (int)floor(x/h);
    double nxt = X[i][0];
    return nxt;
}

double Fn(double x) {
    double fn = kt * log(n(x) / ni);
    return fn;
}

int main() {

    double Jf[MS][MS], Jg[MS][MS], Jh[MS][MS], Uf[MS], Ug[MS], Uh[MS], Rf[MS], Rg[MS], Rh[MS], x[MS];
    double result[MS];
    double k = 1.5;
    double max, length;
    FILE *fptr; // Declare a FILE pointer

    for(int i=0; i<MS; i++){
        double constant = (k*rand()+1)/RAND_MAX;
        X[i][0] = (k+constant)*no;
        X[i][1] = 0;
        x[i] = length;
        length += h;
    }
    X[0][0] = X[MS-1][0] = no;

    for (int i=0; i<MS; i++) {
        for(int j=0; j<MS; j++){
            Jf[i][j] = 0.0;
            Jg[i][j] = 0.0;
            Jh[i][j] = 0.0;
        }if(i > 0 && i < MS-1){
            Jf[i][i - 1] = Jg[i][i - 1] = 1 / h;
            Jf[i][i] = Jg[i][i] = -1 / h;
            Jh[i][i - 1] = Jh[i][i + 1] = 1 / h / h;
            Jh[i][i] = -2 / h / h;
        }
        Uf[i] = Ug[i] = Uh[i] = 0.0;
    }

    Jf[0][0] = Jg[0][0] = Jh[0][0] = 1;
    Jf[MS - 1][MS - 1] = Jg[MS - 1][MS - 1] = Jh[MS - 1][MS - 1] = 1;

    int c,d = 0;
    double tol = 1e-4;
    fptr = fopen("Continuity.txt", "w"); // Open the file in write mode

    for(double j=0.0; j<l; j+=h){
            if(c==0 || c==MS-1){
                Rf[c] = 0.0;
                Rg[c] = 0.0;
                Rh[c] = 0.0;
            }else{
                Rf[c] = (n(j) - n(j-h))/h;
                Rg[c] = (Fn(j) - Fn(j-h))/h;
                Rh[c] = (Fn(j-h) - 2*Fn(j) + Fn(j+h))/h/h;
            }
            Y[c] = Fn(j);
            c++;
    }
    // for(int i=0; i<MS; i++){
    //     printf("%e %e %e\n",Rf[i], Rg[i], Rh[i]);
    // }
    // printf("\n");
    //gaussian(MS, Uf, Jf, Rf);
    for(double i=0.0; i<T; i+=to){
        WriteToFile(MS, x, X, fptr, d); // Call the modified WriteToFile function
        d++;
        
        for(int k=0; k<MS; k++){
            if(i==0)
                result[k] = mun*Rf[k]*Rg[k] + mun*X[k][0]*Rh[k];
            if(i>=0)
                result[k] = mun*Rf[k]*Rg[k] + mun*X[k][0]*Rh[k] - (1/to + 1/Tn)*(X[k][0] - X[k][1]);
        }

        // gaussian(MS, Uf, Jf, Rf);
        // gaussian(MS, Ug, Jg, Rg);
        // gaussian(MS, Uh, Jh, Rh);
        
        // printf("At t: %es\n", i);
        // for(int j=0; j<MS; j++){
        //     printf("%e %e %e\n",Uf[j], Ug[j], Uh[j]);
        // }

        gaussian(MS, Uf, Jf, Rf);

        for(int j=0; j<MS; j++){
                X[j][1] = X[j][0];
                X[j][0] = X[j][0] + Uf[j];
                Y[j] = Fn(j*h);
        }
    //     for(int i=0; i<MS; i++){
    //     printf("%e %e\n",Y[i], X[i][0]);
    // }
    // printf("\n");
        // for(int j=0; j<MS; j++){
        //     printf("%e %e\n",X[j][0], Y[j]);
        // }
        // printf("\n");
        for(int j=1; j<MS-1; j++){
                Rf[j] = (X[j] - X[j-1])/h;
                Rg[j] = (Y[j] - Y[j-1])/h;
                Rh[j] = (Y[j+1] - 2*Y[j] + Y[j-1])/h/h;
        }

        max = fabs(result[0]);
        for(int j=0; j<MS; j++){
            if (fabs(result[j]) > max) {
                max = fabs(result[j]);
            }    
        }
        if(max<tol){
            printf("%e\n%e\n", max, i);
            break;
        }
    }
    return 0;
}
