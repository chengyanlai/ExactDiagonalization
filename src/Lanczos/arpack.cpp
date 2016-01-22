#include <iostream>
#include "src/Lanczos/arpack_wrapper.hpp"
#include "src/Hamiltonian/Hamiltonian.hpp"

template<>
void Hamiltonian<ComplexType>::arpackDiagonalize(int n,
  ComplexType* input_ptr, std::vector<RealType> &evals, int nev, RealType tol){
  /*
  // n        : dimension of the matrix
  // input_ptr: input trail vector pointer
  // nev      : number of eigenvalues to calculate
  // tol      : tolerance. 0 - calculate until machine precision
  */
  // ido: reverse communication parameter, must be zero on first iteration
  int ido = 0;
  // bmat: standard eigenvalue problem A*x=lambda*x
  char bmat = 'I';
  // which: calculate the smallest real part eigenvalue
  char which[] = {'S','R'};
  // resid: the residual vector
  ComplexType *resid = new ComplexType[n];
  memcpy(resid, input_ptr, n * sizeof(ComplexType));
  // the number of columns in v: the number of lanczos vector
  // generated at each iteration, ncv <= n
  // We use the answer to life, the universe and everything, if possible
  int ncv = 42;
  if( n < ncv )
    ncv = n;
  // v containts the lanczos basis vectors
  int ldv = n;
  ComplexType *v = new ComplexType[ldv*ncv];

  int *iparam = new int[11];
  iparam[0] = 1;   // Specifies the shift strategy (1->exact)
  iparam[2] = 3*n; // Maximum number of iterations
  iparam[6] = 1;   /* Sets the mode of dsaupd.
                      1 is exact shifting,
                      2 is user-supplied shifts,
                      3 is shift-invert mode,
                      4 is buckling mode,
                      5 is Cayley mode. */

  int *ipntr = new int[14];
  /* IPNTR   Integer array of length 14.  (OUTPUT)
             Pointer to mark the starting locations in the WORKD and WORKL
             arrays for matrices/vectors used by the Arnoldi iteration.
             -------------------------------------------------------------
             IPNTR(1): pointer to the current operand vector X in WORKD.
             IPNTR(2): pointer to the current result vector Y in WORKD.
             IPNTR(3): pointer to the vector B * X in WORKD when used in
                       the shift-and-invert mode.
             IPNTR(4): pointer to the next available location in WORKL
                       that is untouched by the program.
             IPNTR(5): pointer to the NCV by NCV upper Hessenberg
                       matrix H in WORKL.
             IPNTR(6): pointer to the  ritz value array  RITZ
             IPNTR(7): pointer to the (projected) ritz vector array Q
             IPNTR(8): pointer to the error BOUNDS array in WORKL.
             IPNTR(14): pointer to the NP shifts in WORKL. See Remark 5 below.
             Note: IPNTR(9:13) is only referenced by zneupd. See Remark 2 below.
             IPNTR(9): pointer to the NCV RITZ values of the
                       original system.
             IPNTR(10): Not Used
             IPNTR(11): pointer to the NCV corresponding error bounds.
             IPNTR(12): pointer to the NCV by NCV upper triangular
                        Schur matrix for H.
             IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
                        of the upper Hessenberg matrix H. Only referenced by
                        zneupd if RVEC = .TRUE. See Remark 2 below.
        -------------------------------------------------------------*/
  ComplexType *workd = new ComplexType[3*n];
  /* WORKD   Complex*16 work array of length 3*N.  (REVERSE COMMUNICATION)
             Distributed array to be used in the basic Arnoldi iteration
             for reverse communication.  The user should not use WORKD
             as temporary workspace during the iteration !!!!!!!!!!
             See Data Distribution Note below.  */
  int lworkl = 3*ncv*(ncv+2);
  /* LWORKL  Integer.  (INPUT)
             LWORKL must be at least 3*NCV**2 + 5*NCV.*/
  ComplexType *workl = new ComplexType[lworkl];
  /* WORKL   Complex*16 work array of length LWORKL.  (OUTPUT/WORKSPACE)
             Private (replicated) array on each PE or array allocated on
             the front end.  See Data Distribution Note below.*/
  RealType *rwork = new RealType[ncv];
  /* RWORK   Double precision  work array of length NCV (WORKSPACE)
             Private (replicated) array on each PE or array allocated on
             the front end. */
  int info = 1;
  /* INFO    Integer.  (INPUT/OUTPUT)
    If INFO .EQ. 0, a randomly initial residual vector is used.
    If INFO .NE. 0, RESID contains the initial residual vector,
                    possibly from a previous run.
    Error flag on output.
    =  0: Normal exit.
    =  1: Maximum number of iterations taken.
          All possible eigenvalues of OP has been found. IPARAM(5)
          returns the number of wanted converged Ritz values.
    =  2: No longer an informational error. Deprecated starting
          with release 2 of ARPACK.
    =  3: No shifts could be applied during a cycle of the
          Implicitly restarted Arnoldi iteration. One possibility
          is to increase the size of NCV relative to NEV.
          See remark 4 below.
    = -1: N must be positive.
    = -2: NEV must be positive.
    = -3: NCV-NEV >= 2 and less than or equal to N.
    = -4: The maximum number of Arnoldi update iteration
        must be greater than zero.
    = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
    = -6: BMAT must be one of 'I' or 'G'.
    = -7: Length of private work array is not sufficient.
    = -8: Error return from LAPACK eigenvalue calculation;
    = -9: Starting vector is zero.
    = -10: IPARAM(7) must be 1,2,3.
    = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
    = -12: IPARAM(1) must be equal to 0 or 1.
    = -9999: Could not build an Arnoldi factorization.
             User input error highly likely.  Please
           check actual array dimensions and layout.
             IPARAM(5) returns the size of the current Arnoldi
             factorization.
  */
  /* dneupd parameters
     rvec == 0 : calculate only eigenvalue
     rvec > 0 : calculate eigenvalue and eigenvector */
  int rvec = 1;

  // how many eigenvectors to calculate: 'A' => nev eigenvectors
  char howmny = 'A';

  int *select;
  // when howmny == 'A', this is used as workspace to reorder the eigenvectors
  select = new int[ncv];

  // This vector will return the eigenvalues from the second routine, dseupd.
  ComplexType *d = new ComplexType[nev+1];

  ComplexType *z = 0;
  z = new ComplexType[n*nev];
  /* On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of
       Z represent approximate eigenvectors (Ritz vectors) corresponding
       to the NCONV=IPARAM(5) Ritz values for eigensystem
       A*z = lambda*B*z. */

  // not used if iparam[6] == 1
  ComplexType sigma;
  ComplexType *workev = new ComplexType[3*ncv];
  /*WORKEV  Double precision work array of dimension 3*NCV.  (WORKSPACE)*/

  // first iteration
  znaupd_(&ido, &bmat, &n, &which[0], &nev, &tol, resid, &ncv, v, &ldv,
          iparam, ipntr, workd, workl, &lworkl, rwork, &info);

  while( ido != 99 ){
    mvprod(workd+ipntr[0]-1, workd+ipntr[1]-1, 0.0e0);
    znaupd_(&ido, &bmat, &n, &which[0], &nev, &tol, resid, &ncv, v, &ldv,
            iparam, ipntr, workd, workl, &lworkl, rwork, &info);
  }

  if( info < 0 )
    std::cerr << "Error with znaupd, info = " << info << std::endl;
  else if ( info == 1 )
    std::cerr << "Maximum number of Lanczos iterations reached." << std::endl;
  else if ( info == 3 )
    std::cerr << "No shifts could be applied during implicit Arnoldi update," <<
                 " try increasing NCV." << std::endl;
  zneupd_(&rvec, &howmny, select, d, z, &ldv, &sigma, workev,
          &bmat, &n, &which[0], &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr,
          workd, workl, &lworkl, rwork, &info);
  if ( info != 0 )
    std::cerr << "Error with dneupd, info = " << info << std::endl;
  evals.clear();
  evals.push_back(d[0].real());
  evals.push_back(d[1].real());
  std::cout << d[0] << std::endl;
  std::cout << d[1] << std::endl;
  /* TODO: make complex eigenvector */
  memcpy(input_ptr, z, n * sizeof(ComplexType));
  delete [] workev;
  delete [] z;
  delete [] d;
  delete [] select;
  delete [] rwork;
  delete [] workl;
  delete [] workd;
  delete [] ipntr;
  delete [] iparam;
  delete [] v;
  delete [] resid;
}

template<>
void Hamiltonian<RealType>::arpackDiagonalize(int n, RealType* input_ptr,
  std::vector<RealType> &evals, int nev, RealType tol){
  int ido = 0;
  char bmat = 'I';
  char which[] = {'S','A'};
  RealType *resid = new RealType[n];
  memcpy(resid, input_ptr, n * sizeof(RealType));
  int ncv = 42;
  if( n < ncv )
    ncv = n;
  int ldv = n;
  RealType *v = new RealType[ldv*ncv];
  int *iparam = new int[11];
  iparam[0] = 1;
  iparam[2] = 3*n;
  iparam[6] = 1;
  int *ipntr = new int[11];
  RealType *workd = new RealType[3*n];
  int lworkl = ncv*(ncv+8);
  RealType *workl = new RealType[lworkl];
  int info = 1;
  int rvec = 1;
  char howmny = 'A';
  int *select;
  select = new int[ncv];
  RealType *d = new RealType[nev];
  RealType *z = 0;
  z = new RealType[n*nev];
  RealType sigma;
  dsaupd_(&ido, &bmat, &n, &which[0], &nev, &tol, resid, &ncv, v, &ldv,
          iparam, ipntr, workd, workl, &lworkl, &info);
  while( ido != 99 ){
    mvprod(workd+ipntr[0]-1, workd+ipntr[1]-1, 0.0e0);
    dsaupd_(&ido, &bmat, &n, &which[0], &nev, &tol, resid, &ncv, v, &ldv,
            iparam, ipntr, workd, workl, &lworkl, &info);
  }
  if( info < 0 )
    std::cerr << "Error with dsaupd, info = " << info << std::endl;
  else if ( info == 1 )
    std::cerr << "Maximum number of Lanczos iterations reached." << std::endl;
  else if ( info == 3 )
    std::cerr << "No shifts could be applied during implicit Arnoldi update," <<
                 "try increasing NCV." << std::endl;
  dseupd_(&rvec, &howmny, select, d, z, &ldv, &sigma, &bmat, &n, which, &nev,
          &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  if ( info != 0 )
    std::cerr << "Error with dseupd, info = " << info << std::endl;
  evals.clear();
  evals.push_back(d[0]);
  evals.push_back(d[1]);
  memcpy(input_ptr, z, n * sizeof(RealType));
  delete [] resid;
  delete [] v;
  delete [] iparam;
  delete [] ipntr;
  delete [] workd;
  delete [] workl;
  delete [] d;
  delete [] z;
  delete [] select;
}
