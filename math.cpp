
#include "lapack.hpp"
#include "math.hpp"
#include <cassert>
#include <vector>

using namespace std;
using namespace blitz;
using namespace blitz::tensor;

bool
invert
/////////////////////////////////////////////////////////
/** Inverts a matrix.
 * \note Changes the input!
 *
 * Note: lapack sees A^T, not A.  We don't have to transpose because
 * ((A^T)^-1)^T = A^-1
 *
 * \return true or false based on whether the 
 * inverse exists
 */
(
 FloatData2D& mat
 ///< \pre matrix to invert
 ///< \post inverted matrix
 ) {
  assert(mat.rows() == mat.columns());
  assert(mat.rows() != 0);

  //mat may not be continuous.  Force it to be for this routine.
  FloatData2D working = mat.copy();

  //do the LU decomposition.
  int N = working.rows();
  int LDA = working.rows();
  std::vector<int> index(mat.rows());
  int* IPIV = &index[0];
  Real* A = &working(0,0);
  int INFO = -1;
  assert(sizeof(Real) == sizeof(double));
  dgetrf_(&N,&N,A,&LDA,IPIV,&INFO);
  assert(INFO >= 0);
  if (INFO != 0) {
    return false;
  }
  
  //calculate workspace size.
  int LWORK = -1;
  int optimal_work;
  {
    double WORK;
    INFO = -1;
    dgetri_(&N,A, &LDA,IPIV, &WORK, &LWORK,&INFO);
    assert(INFO == 0);
    optimal_work = (int)WORK;
  }
  //allocate workspace
  LWORK = optimal_work;
  std::vector<Real> work(3*LWORK);
  Real* WORK = &work[0];

  //invert matrix.
  INFO = -1;
  dgetri_(&N,A, &LDA,IPIV, WORK, &LWORK,&INFO);
  assert(INFO >= 0);
  if (INFO == 0) {
    mat = working;
    return true;
  } else {
    return false;
  }
}

Real
det
/////////////////////////////////////////////////////////
(
 const FloatData2D& mat
 ///< \pre matrix to invert
 ///< \post inverted matrix
 ) {
  assert(mat.rows() == mat.columns());
  assert(mat.rows() != 0);

  //mat may not be continuous.  Force it to be for this routine.
  FloatData2D working = mat.copy();

  //do the LU decomposition.
  int N = working.rows();
  int LDA = working.rows();
  std::vector<int> index(mat.rows());
  int* IPIV = &index[0];
  Real* A = &working(0,0);
  int INFO = -1;
  assert(sizeof(Real) == sizeof(double));
  dgetrf_(&N,&N,A,&LDA,IPIV,&INFO);
  assert(INFO >= 0);
  if (INFO != 0) {
    return false;
  }
  
  //calculate the product of the U on the diagonal.
  Real product = 1;
  for (int ii=0; ii<working.rows(); ii++) {
    product *= working(ii,ii);
  }
  //Now figure out the sign.  
  int inversions = 0;
  for (uint ii=0; ii<index.size(); ii++) {
    for (uint jj=ii+1; jj<index.size(); jj++) {
      if (index[ii] > index[jj]) {
        inversions++;
      }
    }
  }
  inversions = (inversions % 2);
  if (inversions) {
    product *= -1;
  }

  return product;
}

void
solve_least_squares
/////////////////////////////////////////////////////////
/** Solves Ax = b, where A has more rows than b
 *
 */
(
 const FloatData2D& A, 
 const FloatData1D& b, 
 FloatData1D& x,
 const Real rcond
 ) {
  assert(A.columns() == x.rows());
  assert(A.rows() == b.rows());
  assert(A.rows() >= A.columns());
  assert(x.rows() > 0);
  assert(b.rows() > 0);


  FloatData2D copy_A(A.shape(), fortranArray);
  copy_A = A;
  FloatData1D copy_b(b.shape());
  copy_b = b;

  int M = A.rows();
  int N = A.columns();
  int NRHS = 1;
  int LDA = M;
  int LDB = M;
  double RCOND = rcond;
  
  //returned
  FloatData1D S(x.rows());
  int RANK = -1;
  int INFO = -1;
  
  //working memory
  int SMLSIZ = 30;
  int NLVL = 1;
  int LWORK = 12*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS +
    (SMLSIZ+1)*(SMLSIZ+1);
  double WORK[LWORK];
  int LIWORK = 3 * 4 * NLVL + 11 * 4;
  int IWORK[LIWORK];
  dgelsd_(&M, &N, &NRHS, &copy_A(1,1), &LDA, &copy_b(0), &LDB, &S(0), &RCOND,
          &RANK, WORK, &LWORK, IWORK, &INFO );
  if (INFO != 0) {
    cout << INFO << endl;
    assert(0 && "the algorithm for computing the SVD failed to converge");
  }
  assert(INFO == 0);
  //assert(RANK == N);

  x = copy_b(Range(0,N-1));
}

bool
find_interpolation_coeff
//---------------------------------------------------
/** Computes interpolation coefficients for two points in a
 *  neighborhood.
 */
(
 const FloatData2D& given_points,
 const FloatData1D& interpolated_point,
 FloatData1D& output_coeff,
 Real rcond
 ) {
  int dim = given_points.columns();
  int neighborhood_size = given_points.rows();

  assert(neighborhood_size == output_coeff.size());
  assert(dim == interpolated_point.size());
  
  
  FloatData2D A(dim, dim);
  A = matmult(transpose(given_points),given_points);
  FloatData1D b(dim);
  b = interpolated_point;

  FloatData1D coeff(neighborhood_size);
  FloatData1D x(dim);
  bool solution_sane = false;
  while (!solution_sane && rcond<=1) {
    solve_least_squares(A, b, x, rcond);
    
    coeff = matmult(x,transpose(given_points));
    
    solution_sane = true;
    for (int ii=0; ii<neighborhood_size; ii++) {
      if (fabs(coeff(ii)) > 1.1) {
        solution_sane = false;
        break;
      }
    }
    Real soln_sum=0;
    for (int ii=0; ii<neighborhood_size; ii++) {
      soln_sum += fabs(coeff(ii));
    }
    solution_sane &= (0.9 < soln_sum);
    solution_sane &= (soln_sum < 1.1);

    if (!solution_sane) {
      rcond *= 10;
    }
  }
  if (solution_sane) {
  } else {  
  }
  output_coeff = coeff;
  return solution_sane;
}

void adjoint(FloatData2D& mat) {
  assert(mat.rows() == mat.columns());
  Real D = det(mat);
  invert(mat);
  mat *= D;
}

void
cross
//---------------------------------------------------------
/** Computes the cross product of two vectors.
 */
(
 const FloatData1D& a,
 const FloatData1D& b,
 FloatData1D& c
 ) {
  assert(a.size() == 3);
  assert(b.size() == 3);
  assert(c.size() == 3);

  c(0) = a(1)*b(2)-a(2)*b(1);
  c(1) = a(2)*b(0)-a(0)*b(2);
  c(2) = a(0)*b(1)-a(1)*b(0);
}


bool disassemble_matrix(FloatData2D& A, FloatData1D& eig) {
  assert(A.rows() == A.columns());
  assert(eig.size() == A.rows());
  char compz='V', uplo='L';
  int n=A.rows(), info;
  double e[n-1], tau[n-1];
  
  FloatData2D working(A.shape());
  working = A.transpose(secondDim,firstDim);

  int lwork;
  {
    lwork = -1;
    double work[1];
    dsytrd_(&uplo,&n, working.data(), &n, eig.data(), e, tau, work, &lwork, &info);
    assert(info == 0);
    lwork = work[0];
  }
  {
    double work[lwork];
    dsytrd_(&uplo,&n, working.data(), &n, eig.data(), e, tau, work, &lwork, &info);
    assert(info == 0);
  }
  {
    lwork = -1;
    double work[1];
    dorgtr_(&uplo,&n, working.data(), &n, tau, work, &lwork, &info);
    assert(info == 0);
    lwork = work[0];
  }
  {
    double work[lwork];
    dorgtr_(&uplo,&n, working.data(), &n, tau, work, &lwork, &info);
    assert(info == 0);
  }    
  {
    double work[2*n-1];
    dsteqr_(&compz, &n, eig.data(), e, working.data(), &n, work, &info);
    assert(info == 0);
  }

  //transpose working, assign to A.
  A = working.transpose(secondDim,firstDim);
  return true;
}

void
normalize
//-----------------------------------------------------------
/** Normalizes a vector
 */
(
 FloatData1D& vec
 ) {
  vec /= magnitude(vec);
}

Real
magnitude
//------------------------------------------------------------
/** Computes the magnitude of a vector
 */
(
 const FloatData1D& vec
 ) {
  return sum(vec*vec);
}

void 
reassemble_matrix
//-----------------------------------------------------------
/** Reassembles a sym pos def matrix from eigenvalues and vectors. 
 * This is the inverse of disassemble matrix.
 */
(
 FloatData2D& Q,
 const FloatData1D& eig
 ) {
  assert(Q.rows() == Q.columns());
  assert(Q.rows() == eig.size());

  // compute Q * D * Q^T
  FloatData2D A(Q.shape());
  A = sum(Q(tensor::i,tensor::k)*eig(tensor::k)*Q(tensor::k,tensor::j),tensor::k);
  Q = A;
}

//void exp_matrix(float a[3][3]);
//void ln_matrix(float a[3][3]);


#if 0 //comment

#include <cblas.h>

inline
void
my_dgemm
//------------------------------------------------------
/** Computes C <= alpha *A*B + beta*C
 */
(
 Real alpha, 
 //! Row major matrix
 const FloatData2D& A, 
 //! Row major matrix
 const FloatData2D& B, 
 Real beta, 
 //! Row major matrix
 FloatData2D& C
 //---------------------------------------------------------
 ) {
  assert(C.rows() == A.rows());
  assert(A.columns() == B.rows());
  assert(C.columns() == B.columns());
  cblas_dgemm(CblasRowMajor,
              CblasNoTrans,
              CblasNoTrans,
              C.rows(),
              C.columns(),
              A.columns(),
              alpha,
              &A(0,0),
              A.columns(),
              &B(0,0),
              B.columns(),
              beta,
              &C(0,0),
              C.columns()
              );
}

inline
void
matmat_mult(const FloatData2D& A, const FloatData2D& B, FloatData2D& C) {
  my_dgemm(1,A,B,0,C);
}

void matmatmat_mult(const FloatData2D& A, const FloatData2D& B, const FloatData2D& C, FloatData2D& D) {
  assert(D.rows() == A.rows());
  assert(A.columns() == B.rows());
  assert(B.columns() == C.rows());
  assert(D.columns() == C.columns());
  
  //stupid blitz.  We have to do this one manually if we don't want to create temp storage.
  FloatData2D temp(B.rows(),C.columns());
  matmat_mult(B,C,temp);
  matmat_mult(A,temp,D);
}

#endif //comment
