#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define q 1.6e-19 //in coulomb
#define kt 0.026 //in eV
#define mun 1400 //in cm^2/V-s - mobility of e-
#define ni 1.5e10 //in 1/cc - steady state e- concn
#define Tn 2e-2 //in 1/s - recombinational lifetime
#define MS 100 //no of nodes
#define T 1e-5 //in s - total time
#define l 6e-6 // in cm - total length
#define no 1.7e10 //in 1/cc - equilibrium conc of e-

double X[MS][2], Y[MS];
/*X[][0]-(Contains entries of e- conc at present time node: t_i), X[][1]-(Contains entries of e- conc at one time node behind the present time-node: t_i-1)
and Y[]-(Contains entries of the function Fn(): quasi-fermi levels)*/

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

    if(c<99)
        fprintf(fptr, "\n");
    fclose(fptr);
}

//function to solve the system of linear equations using gaussian elimination
void gaussian(int n,  double x[], double A[][n], double b[]){
    int count = 0;
    double p = 0;
    double temp[n], a[n][n]; 
    //creating a copy of the residual matrix-(b[]) and Jacobian matrix-(A[][n]) respectively.
    for(int i=0; i<n; i++){
        for(int j=0; j<MS; j++){
            a[i][j] = A[i][j];
        }
        temp[i] = b[i];
    }
    //converting jacobian into upper-triangular matrix-forward elimination
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
    //solving for x[]-(update matrix)-backward substitution
    for(int i=0; i<n; i++){
        x[i] = 0;
    }
    for (int i = n - 1; i >= 0; i--) {
        x[i] = temp[i] / a[i][i];
        for (int j = i + 1; j < n; j++) {
            x[i] = x[i] - (a[i][j] / a[i][i]) * x[j];
        }
    }
}

//function for the e- carrier concentration
double n(double x) {
    int i = (int)floor(x/h);
    double nxt = X[i][0]; //here X[] matrix has entries of the e- conc at every node...
    return nxt;
}

//function for quasi-fermi level
double Fn(double x) {
    double fn = kt * log(n(x) / ni);
    return fn;
}

int main() {
    //declaring the Jacobian, Update and Residual matrix for each function defined
    double Jf[MS][MS], Jg[MS][MS], Jh[MS][MS], Uf[MS], Ug[MS], Uh[MS], Rf[MS], Rg[MS], Rh[MS], x[MS];
    double result[MS];
    double k = 1.5;
    double max, length;
    FILE *fptr; // Declare a FILE pointer
    //generating the initial e- concn (non-uniform distribution using rand())
    for(int i=0; i<MS; i++){
        double constant = (k*rand()+1)/RAND_MAX;
        X[i][0] = (k+constant)*no; 
        X[i][1] = 0;
        x[i] = length;
        length = length + h;
    }
    //Defining the boundary conditions
    X[0][0] = X[MS-1][0] = no;
    //filling the jacobian as well as the Update matrices of all the functions.
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
        Uf[i] = Ug[i] = Uh[i] = 0.0; //Update matrix is filled with zeroes initially
    }
    Jf[0][0] = Jg[0][0] = Jh[0][0] = 1; 
    Jf[MS - 1][MS - 1] = Jg[MS - 1][MS - 1] = Jh[MS - 1][MS - 1] = 1;

    int c,d = 0;
    double tol = 1e-4; //defining tolerence for convergence criteria

    fptr = fopen("Continuity.txt", "w"); // Open the file in write mode
    //Filling the Residual Matrices of each function
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

    for(double i=0.0; i<T; i+=to){
        WriteToFile(MS, x, X, fptr, d); // Call the modified WriteToFile function
        d++;
        // calculating the eqn of continuity:-
        for(int k=0; k<MS; k++){
            if(i==0)
                result[k] = mun*Rf[k]*Rg[k] + mun*X[k][0]*Rh[k];
            if(i>0)
                result[k] = mun*Rf[k]*Rg[k] + mun*X[k][0]*Rh[k] - (1/to + 1/Tn)*(X[k][0] - X[k][1]);
        }

        gaussian(MS, Uf, Jf, Rf); //calculting the update matrix-Uf[]
        //Updating X[] and also the Y[]-(the Y[] contains the entries for the function Fn at each node...):
        for(int j=0; j<MS; j++){
                X[j][1] = X[j][0];
                X[j][0] = X[j][0] + Uf[j];
                Y[j] = Fn(j*h);
        }
        //Updating the Residual matrices accordingly:
        for(int j=1; j<MS-1; j++){
                Rf[j] = (X[j] - X[j-1])/h;
                Rg[j] = (Y[j] - Y[j-1])/h;
                Rh[j] = (Y[j+1] - 2*Y[j] + Y[j-1])/h/h;
        }

        max = fabs(result[0]);
        //having a convergence check
        for(int j=0; j<MS; j++){
            if (fabs(result[j]) > max) {
                max = fabs(result[j]);
            }    
        }
        //breaking the loop if convergence criteria is met...
        if(max<tol){
            printf("%e\n%e\n", max, i);
            break;
        }

    }
    return 0;
}
