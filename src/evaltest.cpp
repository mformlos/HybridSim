#include <iostream> 
#include <cmath>


int eigenvals(double **matrix, double *evals) {
    double sum_diag, sum_offdiag;      // Sums of diagonal resp. off-diagonal elements
    double sin_phi, cos_phi, t;  // sin(phi), cos(phi), tan(phi) and temporary storage
    double g, h, z, theta;             // More temporary storage
    double thresh;
    
    //initialize evals to diag(matrix)
    for (int i = 0; i < 3; i++) {
        evals[i] = matrix[i][i]; 
    }
    
    //calculate sum of diagonal elements and square it 
    sum_diag = 0.0; 
    for (int i = 0; i < 3; i++) {
        sum_diag += fabs(evals[i]); 
    }
    sum_diag *= sum_diag; 
    
    //iteration
    for (int n = 0; n < 50; n++) {
        //test convergence
        sum_offdiag = 0.0; 
        for (int p = 0; p < 3; p++) {
            for (int q = p+1; q < 3; q++) {
                sum_offdiag += fabs(matrix[p][q]); 
            }
        }
        if (sum_offdiag == 0.0) {
            std::cout << "took " << n << " iterations" << std::endl; 
            return 0;
        } 
        
        if (n < 4) thresh = 0.2*sum_offdiag / 9; 
        else thresh = 0.0; 
        
        // sweep
        for (int p = 0; p < 3; p++) {
            for (int q = p+1; q < 3; q++) {
                g = 100.0*fabs(matrix[p][q]); 
                if (n > 4 && fabs(evals[p])+g == fabs(evals[p]) && fabs(evals[q])+g == fabs(evals[q])) {
                    matrix[p][q] = 0.0; 
                }
                else if (fabs(matrix[p][q]) > thresh) {
                    //calculate Jacobi iteration
                    h = evals[q] - evals[p]; 
                    if (fabs(h)+g == fabs(h)) {
                        t = matrix[p][q] / h; 
                    }
                    else {
                        theta = 0.5 * h / matrix[p][q]; 
                        if (theta < 0.0) t = -1.0 /(sqrt(1.0 + theta*theta)- theta); 
                        else t = -1.0 /(sqrt(1.0 + theta*theta) + theta); 
                    }
                    cos_phi = 1.0/sqrt(1.0+t*t); 
                    sin_phi = t*cos_phi; 
                    z = t*matrix[p][q]; 
                    
                    //apply trafo
                    matrix[p][q] = 0.0; 
                    evals[p] -= z; 
                    evals[q] += z; 
                    for (int r = 0; r < p; r++) {
                        t = matrix[r][p]; 
                        matrix[r][p] = cos_phi*t - sin_phi*matrix[r][q]; 
                        matrix[r][q] = sin_phi*t + cos_phi*matrix[r][q]; 
                    }
                    for (int r = p+1; r < q; r++) {
                        t = matrix[p][r]; 
                        matrix[p][r] = cos_phi*t - sin_phi*matrix[r][q]; 
                        matrix[r][q] = sin_phi*t + cos_phi*matrix[r][q];
                    }
                    for (int r = q+1; r < 3; r++) {
                        t = matrix[p][r]; 
                        matrix[p][r] = cos_phi*t - sin_phi*matrix[q][r]; 
                        matrix[q][r] = sin_phi*t + cos_phi*matrix[q][r];
                    }
                }
            }
        }
    }
    return -1; 
    
}

void straightsort(int n, double* arr) {
    double a; 
    
    for (int j = 1; j < n; j++) {
        a = arr[j]; 
        int i = j - 1; 
        while(i >= 0 && arr[i] > a) {
            arr[i+1] = arr[i]; 
            i--; 
        }
        arr[i+1] = a; 
    }
    
}

#define SQR(x)      ((x)*(x))                        // x^2 


// ----------------------------------------------------------------------------
int jacobi(double A[3][3], double Q[3][3], double w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
// matrix A using the Jacobi algorithm.
// The upper triangular part of A is destroyed during the calculation,
// the diagonal elements are read but not destroyed, and the lower
// triangular elements are not referenced at all.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error (no convergence)
// ----------------------------------------------------------------------------
{
  const int n = 3;
  double sd, so;                  // Sums of diagonal resp. off-diagonal elements
  double s, c, t;                 // sin(phi), cos(phi), tan(phi) and temporary storage
  double g, h, z, theta;          // More temporary storage
  double thresh;
  
  // Initialize Q to the identitity matrix
#ifndef EVALS_ONLY
  for (int i=0; i < n; i++)
  {
    Q[i][i] = 1.0;
    for (int j=0; j < i; j++)
      Q[i][j] = Q[j][i] = 0.0;
  }
#endif

  // Initialize w to diag(A)
  for (int i=0; i < n; i++)
    w[i] = A[i][i];

  // Calculate SQR(tr(A))  
  sd = 0.0;
  for (int i=0; i < n; i++)
    sd += fabs(w[i]);
  sd = SQR(sd);
 
  // Main iteration loop
  for (int nIter=0; nIter < 50; nIter++)
  {
    // Test for convergence 
    so = 0.0;
    for (int p=0; p < n; p++)
      for (int q=p+1; q < n; q++)
        so += fabs(A[p][q]);
    if (so == 0.0)
      return 0;

    if (nIter < 4)
      thresh = 0.2 * so / SQR(n);
    else
      thresh = 0.0;

    // Do sweep
    for (int p=0; p < n; p++)
      for (int q=p+1; q < n; q++)
      {
        g = 100.0 * fabs(A[p][q]);
        if (nIter > 4  &&  fabs(w[p]) + g == fabs(w[p])
                       &&  fabs(w[q]) + g == fabs(w[q]))
        {
          A[p][q] = 0.0;
        }
        else if (fabs(A[p][q]) > thresh)
        {
          // Calculate Jacobi transformation
          h = w[q] - w[p];
          if (fabs(h) + g == fabs(h))
          {
            t = A[p][q] / h;
          }
          else
          {
            theta = 0.5 * h / A[p][q];
            if (theta < 0.0)
              t = -1.0 / (sqrt(1.0 + SQR(theta)) - theta);
            else
              t = 1.0 / (sqrt(1.0 + SQR(theta)) + theta);
          }
          c = 1.0/sqrt(1.0 + SQR(t));
          s = t * c;
          z = t * A[p][q];

          // Apply Jacobi transformation
          A[p][q] = 0.0;
          w[p] -= z;
          w[q] += z;
          for (int r=0; r < p; r++)
          {
            t = A[r][p];
            A[r][p] = c*t - s*A[r][q];
            A[r][q] = s*t + c*A[r][q];
          }
          for (int r=p+1; r < q; r++)
          {
            t = A[p][r];
            A[p][r] = c*t - s*A[r][q];
            A[r][q] = s*t + c*A[r][q];
          }
          for (int r=q+1; r < n; r++)
          {
            t = A[p][r];
            A[p][r] = c*t - s*A[q][r];
            A[q][r] = s*t + c*A[q][r];
          }

          // Update eigenvectors
#ifndef EVALS_ONLY          
          for (int r=0; r < n; r++)
          {
            t = Q[r][p];
            Q[r][p] = c*t - s*Q[r][q];
            Q[r][q] = s*t + c*Q[r][q];
          }
#endif
        }
      }
  }

  return -1;
}


int main() {
    /*double **matrix; 
    double *evals; 
    
    matrix = (double**)calloc(3, sizeof(double*)); 
    for (int i = 0; i < 3; i++) {
        matrix[i] = (double*)calloc(3, sizeof(double)); 
    }
    evals = (double*)calloc(3, sizeof(double)); 
    */
    double matrix[3][3], Q[3][3]; 
    double evals[3]; 
    
    matrix[0][0] = 1.0; 
    matrix[0][1] = 3.0;
    matrix[0][2] = 0.0;
    matrix[1][0] = matrix[0][1];
    matrix[1][1] = 2.0;
    matrix[1][2] = 6.0; 
    matrix[2][0] = matrix[0][2];
    matrix[2][1] = matrix[1][2];
    matrix[2][2] = 5.0;  
    
    jacobi(matrix, Q, evals);
    //eigenvals(matrix,evals); 
    
    double a; 
    
    a =  evals[2]; 
    evals[2] = evals[0]; 
    evals[0] = a;   
    
    std::cout << "eigenvalues are: " << evals[0] << " " << evals[1] << " " << evals[2] << std::endl; 
    
    straightsort(3, evals); 
    
    
    std::cout << "sorted eigenvalues are: " << evals[0] << " " << evals[1] << " " << evals[2] << std::endl; 
    
    
    
    return 0; 
}
