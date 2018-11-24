#include <math.h>
#include <string.h>
#include "timer.h"
#include <stdio.h>

#define NN 128
#define NM 128

float A[NN][NM];
float Anew[NN][NM];

int main(int argc, char** argv)
{
    int i,j;
    const int n = NN;
    const int m = NM;
    const int iter_max = 1000000;
    
    const double tol = 1.0e-7;
    double error     = 1.0;
    
    memset(A, 0, n * m * sizeof(float));
    memset(Anew, 0, n * m * sizeof(float));
        
    for (j = 0; j < n; j++)
    {
        A[j][0]    = 1.0;
        Anew[j][0] = 1.0;
    }
    
    printf("Jacobi relaxation Calculation: %d x %d mesh\n", n, m);
    
    StartTimer();
    int iter = 0;
    
    while ( error > tol && iter < iter_max )
    {
        error = 0.0;

        for( j = 1; j < n-1; j++)
        {
            for( i = 1; i < m-1; i++ )
            {
                Anew[j][i] = 0.25 * ( A[j][i+1] + A[j][i-1]
                                    + A[j-1][i] + A[j+1][i]);
                //error= fmax( error, fabs(Anew[j][i] - A[j][i]));

                //if(error<(fabs(Anew[j][i] - A[j][i])))error=fabs(Anew[j][i] - A[j][i]);

                //if(error<((Anew[j][i] - A[j][i])>0?(Anew[j][i] - A[j][i]):-(Anew[j][i] - A[j][i])))error=((Anew[j][i] - A[j][i])>0?(Anew[j][i] - A[j][i]):-(Anew[j][i] - A[j][i]));

                error=error<(fabs(Anew[j][i] - A[j][i]))?fabs(Anew[j][i] - A[j][i]):error;
            }
        }
        
        for( j = 1; j < n-1; j++)
        {
            for( i = 1; i < m-1; i++ )
            {
                A[j][i] = Anew[j][i];    
            }
        }

        if(iter % 4000 == 0) printf("%5d, %0.10f\n", iter, error);
        
        iter++;
    }

    double runtime = GetTimer();
 
    printf(" total: %f s\n", runtime / 1000);

    return 0;
}