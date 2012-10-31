//// HEADER GUARD ///////////////////////////
// If automatically generated, keep above
// comment as first line in file.
#ifndef __LAPACK_HPP__
#define __LAPACK_HPP__
//// HEADER GUARD ///////////////////////////

extern "C" {
  void dgetri_(const int* N, 
               double A[], 
               const int* LDA, 
               int IPIV[], 
               double WORK[], 
               int* LWORK, 
               int* INFO);

  void dgetrf_(const int* M, 
               const int* N, 
               double A[], 
               const int* LDA, 
               int IPIV[], 
               int* INFO);
  void dsytrd_(char* uplo, int* n, double* a, int* lda, 
               double * d, double* e, double* tau, 
               double* work, int* lwork, int* info);
  void dorgtr_(char* uplo, int* n, double* a, int* lda, double* tau,
               double* work, int* lwork, int* info);
  void dsteqr_(char* compz, int* n, double *d, double *e, double *z, int* ldz,
               double* work, int* info);
  void dgelsd_(const int* M,
               const int* N,
               const int* NRHS,
               const double A[],
               const int* LDA,
               double B[],
               const int* LDB,
               double S[],
               const double* RCOND,
               int* RANK,
               double WORK[],
               const int* LWORK,
               const int IWORK[],
               int* INFO);
};

//// HEADER GUARD ///////////////////////////
#endif
//// HEADER GUARD ///////////////////////////
