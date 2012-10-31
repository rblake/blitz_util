//// HEADER GUARD ///////////////////////////
// If automatically generated, keep above
// comment as first line in file.
#ifndef __MATH_HPP__
#define __MATH_HPP__
//// HEADER GUARD ///////////////////////////
#include "blitzarray.hpp"

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
 );

void
adjoint
//////////////////////////////////////////////////////////
/** Computes the adjoint of a matrix in place
 *  WARNING, should work for singular matrices, but doesn't.
 */
(
 FloatData2D& mat
 ///< \pre matrix to compute adjoint for
 ///< \post matrix adjoint
 );

Real
det
/////////////////////////////////////////////////////////
(
 const FloatData2D& mat
 );

void
cross
//---------------------------------------------------------
/** Computes the cross product of two vectors.
 */
(
 const FloatData1D& a,
 const FloatData1D& b,
 FloatData1D& c
 );

void
normalize
//-----------------------------------------------------------
/** Normalizes a vector
 */
(
 FloatData1D& vec
 );

Real
magnitude
//------------------------------------------------------------
/** Computes the magnitude of a vector
 */
(
 const FloatData1D& vec
 );

bool
disassemble_matrix
//-----------------------------------------------------------
/** Disassembles a symmetric positive definite matrix into it's
 * eigenvalues and eigenvectors.
 */
(
 FloatData2D& A, 
 FloatData1D& eig
 );

void
reassemble_matrix
//-----------------------------------------------------------
/** Reassembles a sym pos def matrix from eigenvalues and vectors. 
 * This is the inverse of disassemble matrix.
 */
(
 FloatData2D& A,
 const FloatData1D& eig
 );


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
 const Real rcond=1e-7
 );

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
 Real rcond=1e-4
 );

//// HEADER GUARD ///////////////////////////
#endif
//// HEADER GUARD ///////////////////////////
