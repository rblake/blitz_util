//// HEADER GUARD ///////////////////////////
// If automatically generated, keep above
// comment as first line in file.
#ifndef __BLITZARRAY_HPP__
#define __BLITZARRAY_HPP__
//// HEADER GUARD ///////////////////////////

#include <blitz/array.h>

typedef double Real;
typedef unsigned int uint;

typedef blitz::Array<Real, 1> FloatData1D;
typedef blitz::Array<Real, 2> FloatData2D;
typedef blitz::Array<Real, 3> FloatData3D;
typedef blitz::Array<Real, 4> FloatData4D;
typedef blitz::Array<Real, 5> FloatData5D;
 
typedef blitz::Array<int, 1> IntData1D;
typedef blitz::Array<int, 2> IntData2D;
typedef blitz::Array<int, 3> IntData3D;
typedef blitz::Array<int, 4> IntData4D;
typedef blitz::Array<int, 5> IntData5D;

#define ALL blitz::Range::all()

#define tensor_sum(A,B) sum(A,B)
#define total_sum(A) sum(A)

//template<typename TTT>
//typedef typename blitz::Array<TTT,2> Array2D<TTT>;

//template<typename TTT>
//typedef blitz::Array<TTT,1> Array1D<TTT>;


template<typename TTT>
inline const blitz::Array<TTT,2> transpose(blitz::Array<TTT,2> A) {
  return A.transpose(blitz::secondDim,blitz::firstDim);
}

template<typename TTT>
inline 
blitz::_bz_ArrayExpr<blitz::_bz_ArrayExprReduce<blitz::_bz_ArrayExpr<blitz::_bz_ArrayExprBinaryOp<blitz::_bz_ArrayExpr<blitz::ArrayIndexMapping<TTT, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0> >, blitz::_bz_ArrayExpr<blitz::ArrayIndexMapping<TTT, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0> >, blitz::Multiply<TTT, TTT> > >, 2, blitz::ReduceSum<TTT, TTT> > >
//const Array2D 
matmult(const blitz::Array<TTT,2>& A, const blitz::Array<TTT,2>& B) {
  return blitz::tensor_sum(A(blitz::tensor::i, blitz::tensor::k)*B(blitz::tensor::k, blitz::tensor::j),blitz::tensor::k);
}

template<typename TTT>
inline 
blitz::_bz_ArrayExpr<blitz::_bz_ArrayExprReduce<blitz::_bz_ArrayExpr<blitz::_bz_ArrayExprReduce<blitz::_bz_ArrayExpr<blitz::_bz_ArrayExprBinaryOp<blitz::_bz_ArrayExpr<blitz::_bz_ArrayExprBinaryOp<blitz::_bz_ArrayExpr<blitz::ArrayIndexMapping<TTT, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0> >, blitz::_bz_ArrayExpr<blitz::ArrayIndexMapping<TTT, 2, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0> >, blitz::Multiply<TTT, TTT> > >, blitz::_bz_ArrayExpr<blitz::ArrayIndexMapping<TTT, 2, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0> >, blitz::Multiply<TTT, TTT> > >, 3, blitz::ReduceSum<TTT, TTT> > >, 2, blitz::ReduceSum<TTT, TTT> > >
//const blitz::Array<TTT,2> 
matmult(const blitz::Array<TTT,2>& A, const blitz::Array<TTT,2>& B, const blitz::Array<TTT,2>& C) {
  return blitz::tensor_sum(blitz::tensor_sum(A(blitz::tensor::i, blitz::tensor::k)*B(blitz::tensor::k, blitz::tensor::l)*C(blitz::tensor::l,blitz::tensor::j),blitz::tensor::l),blitz::tensor::k);
}

template<typename TTT>
inline 
blitz::_bz_ArrayExpr<blitz::_bz_ArrayExprReduce<blitz::_bz_ArrayExpr<blitz::_bz_ArrayExprBinaryOp<blitz::_bz_ArrayExpr<blitz::ArrayIndexMapping<TTT, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0> >, blitz::_bz_ArrayExpr<blitz::ArrayIndexMapping<TTT, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0> >, blitz::Multiply<TTT, TTT> > >, 1, blitz::ReduceSum<TTT, TTT> > >
//const blitz::Array<TTT,1>
matmult(const blitz::Array<TTT, 1>& A, const blitz::Array<TTT, 2>& B) {
  return blitz::tensor_sum(A(blitz::tensor::j)*B(blitz::tensor::j, blitz::tensor::i),blitz::tensor::j);
}

template<typename TTT>
inline 
blitz::_bz_ArrayExpr<blitz::_bz_ArrayExprReduce<blitz::_bz_ArrayExpr<blitz::_bz_ArrayExprBinaryOp<blitz::_bz_ArrayExpr<blitz::ArrayIndexMapping<TTT, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0> >, blitz::_bz_ArrayExpr<blitz::ArrayIndexMapping<TTT, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0> >, blitz::Multiply<TTT, TTT> > >, 1, blitz::ReduceSum<TTT, TTT> > >
//const blitz::Array<TTT,1> 
matmult(const blitz::Array<TTT,2>& A, const blitz::Array<TTT,1>& B) {
  return blitz::tensor_sum(A(blitz::tensor::i, blitz::tensor::j)*B(blitz::tensor::j),blitz::tensor::j);
}

#define ARRAY_SIZE(x) (sizeof(x)/sizeof(x[0]))

//// HEADER GUARD ///////////////////////////
#endif
//// HEADER GUARD ///////////////////////////
