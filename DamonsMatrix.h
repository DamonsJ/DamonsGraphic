#ifndef _DAMONS_MATRIX_H_
#define _DAMONS_MATRIX_H_

#include "DamonsVector.h"

/// @cond MATHFU_INTERNAL
/// This will unroll loops for matrices with <= 4 columns
#define DAMONSMATH_MAT_OPERATION(OP) DAMONSMATH_UNROLLED_LOOP(i, columns, OP)
/// @endcond

/// @cond MATHFU_INTERNAL
/// This will perform a given OP on each matrix column and return the result
#define DAMONSMATH_MAT_OPERATOR(OP)                  \
  {                                                  \
    Matrix<T, rows, columns> result;                 \
    DAMONSMATH_MAT_OPERATION(result.data_[i] = (OP)); \
    return result;                                    \
  }
/// @endcond

/// @cond MATHFU_INTERNAL
/// This will perform a given OP on each matrix column
#define DAMONSMATH_MAT_SELF_OPERATOR(OP) \
  {                                      \
    DAMONSMATH_MAT_OPERATION(OP);        \
    return *this;                        \
  }


#define DAMONSMATH_MATRIX_4X4_DOT(data1, data2, r)               \
  ((data1)[r] * (data2)[0] + (data1)[(r) + 4] * (data2)[1] + \
   (data1)[(r) + 8] * (data2)[2] + (data1)[(r) + 12] * (data2)[3])



#define DAMONSMATH_MATRIX_3X3_DOT(data1, data2, r, size)              \
  ((data1)[r] * (data2)[0] + (data1)[(r) + (size)] * (data2)[1] + \
   (data1)[(r) + 2 * (size)] * (data2)[2])

#define DAMONSMATH_VECTOR_STRIDE_FLOATS(vector) (sizeof(vector) / sizeof(float))

/// @endcond
namespace DMath {

	/// @{
	/// @class Matrix
	/// @brief Matrix stores a set of "rows" by "columns" elements of type T
	/// and provides functions that operate on the set of elements.
	///
	/// @tparam T type of each element in the matrix.
	/// @tparam rows Number of rows in the matrix.
	/// @tparam columns Number of columns in the matrix.
	template<class T,int rows,int columns = rows>
	class DMatrix {

	public:
		/// @brief Construct a Matrix of uninitialized values.
		inline DMatrix() {}

		/// @brief Construct a Matrix from another Matrix copying each element.
		////
		/// @param m Matrix that the data will be copied from.
		inline DMatrix(const DMatrix<T, rows, columns>& m) {
			DAMONSMATH_MAT_OPERATION(data_[i] = m.data_[i]);
		}

		/// @brief Construct a Matrix from a single float.
		///
		/// @param s Scalar value used to initialize each element of the matrix.
		explicit inline DMatrix(const T& s) {
			DAMONSMATH_MAT_OPERATION((data_[i] = DVector<T, rows>(s)));
		}

		/// @brief Construct a Matrix from four floats.
		///
		/// @note This method only works with a 2x2 Matrix.
		///
		/// @param s00 Value of the first row and column.
		/// @param s10 Value of the second row, first column.
		/// @param s01 Value of the first row, second column.
		/// @param s11 Value of the second row and column.
		inline DMatrix(const T& s00, const T& s10, const T& s01, const T& s11) {
			assert(rows == 2 && columns == 2);
			data_[0] = DVector<T, rows>(s00, s10);
			data_[1] = DVector<T, rows>(s01, s11);
		}

		/// @brief Create a Matrix from nine floats.
		///
		/// @note This method only works with a 3x3 Matrix.
		///
		/// @param s00 Value of the first row and column.
		/// @param s10 Value of the second row, first column.
		/// @param s20 Value of the third row, first column.
		/// @param s01 Value of the first row, second column.
		/// @param s11 Value of the second row and column.
		/// @param s21 Value of the third row, second column.
		/// @param s02 Value of the first row, third column.
		/// @param s12 Value of the second row, third column.
		/// @param s22 Value of the third row and column.
		inline DMatrix( const T& s00, const T& s10, const T& s20, 
						const T& s01, const T& s11, const T& s21, 
						const T& s02, const T& s12, const T& s22) {
			assert(rows == 3 && columns == 3);
			data_[0] = DVector<T, rows>(s00, s10, s20);
			data_[1] = DVector<T, rows>(s01, s11, s21);
			data_[2] = DVector<T, rows>(s02, s12, s22);
		}

		/// @brief Creates a Matrix from twelve floats.
		///
		/// @note This method only works with Matrix<float, 4, 3>.
		///
		///
		/// @param s00 Value of the first row and column.
		/// @param s10 Value of the second row, first column.
		/// @param s20 Value of the third row, first column.
		/// @param s30 Value of the fourth row, first column.
		/// @param s01 Value of the first row, second column.
		/// @param s11 Value of the second row and column.
		/// @param s21 Value of the third row, second column.
		/// @param s31 Value of the fourth row, second column.
		/// @param s02 Value of the first row, third column.
		/// @param s12 Value of the second row, third column.
		/// @param s22 Value of the third row and column.
		/// @param s32 Value of the fourth row, third column.
		inline DMatrix(const T& s00, const T& s10, const T& s20, const T& s30,
					  const T& s01, const T& s11, const T& s21, const T& s31,
					  const T& s02, const T& s12, const T& s22, const T& s32) {

			assert(rows == 4 && columns == 3);
			data_[0] = DVector<T, rows>(s00, s10, s20, s30);
			data_[1] = DVector<T, rows>(s01, s11, s21, s31);
			data_[2] = DVector<T, rows>(s02, s12, s22, s32);
		}

		/// @brief Create a Matrix from sixteen floats.
		///
		/// @note This method only works with a 4x4 Matrix.
		///
		/// @param s00 Value of the first row and column.
		/// @param s10 Value of the second row, first column.
		/// @param s20 Value of the third row, first column.
		/// @param s30 Value of the fourth row, first column.
		/// @param s01 Value of the first row, second column.
		/// @param s11 Value of the second row and column.
		/// @param s21 Value of the third row, second column.
		/// @param s31 Value of the fourth row, second column.
		/// @param s02 Value of the first row, third column.
		/// @param s12 Value of the second row, third column.
		/// @param s22 Value of the third row and column.
		/// @param s32 Value of the fourth row, third column.
		/// @param s03 Value of the first row, fourth column.
		/// @param s13 Value of the second row, fourth column.
		/// @param s23 Value of the third row, fourth column.
		/// @param s33 Value of the fourth row and column.
		inline DMatrix( const T& s00, const T& s10, const T& s20, const T& s30,
					    const T& s01, const T& s11, const T& s21, const T& s31,
						const T& s02, const T& s12, const T& s22, const T& s32,
						const T& s03, const T& s13, const T& s23, const T& s33) {
			assert(rows == 4 && columns == 4);
			data_[0] = DVector<T, rows>(s00, s10, s20, s30);
			data_[1] = DVector<T, rows>(s01, s11, s21, s31);
			data_[2] = DVector<T, rows>(s02, s12, s22, s32);
			data_[3] = DVector<T, rows>(s03, s13, s23, s33);
		}

		/// @brief Create 4x4 Matrix from 4, 4 element vectors.
		///
		/// @note This method only works with a 4x4 Matrix.
		///
		/// @param column0 Vector used for the first column.
		/// @param column1 Vector used for the second column.
		/// @param column2 Vector used for the third column.
		/// @param column3 Vector used for the fourth column.
		inline DMatrix( const DVector<T, 4>& column0, const DVector<T, 4>& column1,
						const DVector<T, 4>& column2, const DVector<T, 4>& column3) {
			assert(rows == 4 && columns == 4);
			data_[0] = column0;
			data_[1] = column1;
			data_[2] = column2;
			data_[3] = column3;
		}

		/// @brief Create a Matrix from the first row * column elements of an array.
		///
		/// @param a Array of values that the matrix will be iniitlized to.
		explicit inline DMatrix(const T* const a) {
			DAMONSMATH_MAT_OPERATION((data_[i] = DVector<T, rows>(&a[i * columns])));
		}

		/// @brief Access an element of the matrix.
		///
		/// @param row Index of the row to access.
		/// @param column Index of the column to access.
		/// @return Const reference to the element.
		inline const T& operator()(const int row, const int column) const {
			return data_[column][row];
		}

		/// @brief Access an element of the Matrix.
		///
		/// @param row Index of the row to access.
		/// @param column Index of the column to access.
		/// @return Reference to the data that can be modified by the caller.
		inline T& operator()(const int row, const int column) {
			return data_[column][row];
		}

		/// @brief Access an element of the Matrix.
		///
		/// @param i Index of the element to access in flattened memory.  Where
		/// the column accessed is i / rows and the row is i % rows.
		/// @return Reference to the data that can be modified by the caller.
		inline const T& operator()(const int i) const { return operator[](i); }

		/// @brief Access an element of the Matrix.
		///
		/// @param i Index of the element to access in flattened memory.  Where
		/// the column accessed is i / rows and the row is i % rows.
		/// @return Reference to the data that can be modified by the caller.
		inline T& operator()(const int i) { return operator[](i); }

		/// @brief Access an element of the Matrix.
		///
		/// @param i Index of the element to access in flattened memory.  Where
		/// the column accessed is i / rows and the row is i % rows.
		/// @return Const reference to the data.
		inline const T& operator[](const int i) const {
			return const_cast<DMatrix<T, rows, columns>*>(this)->operator[](i);
		}

		/// @brief Access an element of the Matrix.
		///
		/// @param i Index of the element to access in flattened memory.  Where
		/// the column accessed is i / rows and the row is i % rows.
		/// @return Reference to the data that can be modified by the caller.
		inline T& operator[](const int i) {
			return reinterpret_cast<T*>(data_)[i];
		}

		/// @cond MATHFU_INTERNAL
		/// @brief Access a column vector of the Matrix.
		///
		/// @param i Index of the column to access.
		/// @return Reference to the data that can be modified by the caller.
		inline DVector<T, rows>& GetColumn(const int i) { return data_[i]; }

		/// @brief Access a column vector of the Matrix.
		///
		/// @param i Index of the column to access.
		/// @return Const reference to the data.
		inline const DVector<T, rows>& GetColumn(const int i) const {
			return data_[i];
		}

		/// @brief Negate this Matrix.
		///
		/// @return Matrix containing the result.
		inline DMatrix<T, rows, columns> operator-() const {
			DAMONSMATH_MAT_OPERATOR(-data_[i]);
		}

		/// @brief Add a Matrix to this Matrix.
		///
		/// @param m Matrix to add to this Matrix.
		/// @return Matrix containing the result.
		inline DMatrix<T, rows, columns> operator+(const DMatrix<T, rows, columns>& m) const {
			DAMONSMATH_MAT_OPERATOR(data_[i] + m.data_[i]);
		}

		/// @brief Subtract a Matrix from this Matrix.
		///
		/// @param m Matrix to subtract from this Matrix.
		/// @return Matrix containing the result.
		inline DMatrix<T, rows, columns> operator-(const DMatrix<T, rows, columns>& m) const {
			DAMONSMATH_MAT_OPERATOR(data_[i] - m.data_[i]);
		}

		/// @brief Add a scalar to each element of this Matrix.
		///
		/// @param s Scalar to add to this Matrix.
		/// @return Matrix containing the result.
		inline DMatrix<T, rows, columns> operator+(const T& s) const {
			DAMONSMATH_MAT_OPERATOR(data_[i] + s);
		}

		/// @brief Subtract a scalar from each element of this Matrix.
		///
		/// @param s Scalar to subtract from this matrix.
		/// @return Matrix containing the result.
		inline DMatrix<T, rows, columns> operator-(const T& s) const {
			DAMONSMATH_MAT_OPERATOR(data_[i] - s);
		}

		/// @brief Multiply each element of this Matrix with a scalar.
		///
		/// @param s Scalar to multiply with this Matrix.
		/// @return Matrix containing the result.
		inline DMatrix<T, rows, columns> operator*(const T& s) const {
			DAMONSMATH_MAT_OPERATOR(data_[i] * s);
		}

		/// @brief Divide each element of this Matrix with a scalar.
		///
		/// @param s Scalar to divide this Matrix with.
		/// @return Matrix containing the result.
		inline DMatrix<T, rows, columns> operator/(const T& s) const {
			return (*this) * (T(1) / s);
		}

		/// @brief Multiply this Matrix with another Matrix.
		///
		/// @param m Matrix to multiply with this Matrix.
		/// @return Matrix containing the result.
		inline DMatrix<T, rows, columns> operator*(const DMatrix<T, rows, columns>& m) const {
			DMatrix<T, rows, columns> result;
			TimesHelper(*this, m, &result);
			return result;
		}

		/// @brief Add a Matrix to this Matrix (in-place).
		///
		/// @param m Matrix to add to this Matrix.
		/// @return Reference to this class.
		inline DMatrix<T, rows, columns>& operator+=(const DMatrix<T, rows, columns>& m) {
			DAMONSMATH_MAT_SELF_OPERATOR(data_[i] += m.data_[i]);
		}

		/// @brief Subtract a Matrix from this Matrix (in-place).
		///
		/// @param m Matrix to subtract from this Matrix.
		/// @return Reference to this class.
		inline DMatrix<T, rows, columns>& operator-=(const DMatrix<T, rows, columns>& m) {
			DAMONSMATH_MAT_SELF_OPERATOR(data_[i] -= m.data_[i]);
		}

		/// @brief Add a scalar to each element of this Matrix (in-place).
		///
		/// @param s Scalar to add to each element of this Matrix.
		/// @return Reference to this class.
		inline DMatrix<T, rows, columns>& operator+=(const T& s) {
			DAMONSMATH_MAT_SELF_OPERATOR(data_[i] += s);
		}

		/// @brief Subtract a scalar from each element of this Matrix (in-place).
		///
		/// @param s Scalar to subtract from each element of this Matrix.
		/// @return Reference to this class.
		inline DMatrix<T, rows, columns>& operator-=(const T& s) {
			DAMONSMATH_MAT_SELF_OPERATOR(data_[i] -= s);
		}

		/// @brief Multiply each element of this Matrix with a scalar (in-place).
		///
		/// @param s Scalar to multiply with each element of this Matrix.
		/// @return Reference to this class.
		inline DMatrix<T, rows, columns>& operator*=(const T& s) {
			DAMONSMATH_MAT_SELF_OPERATOR(data_[i] *= s);
		}

		/// @brief Divide each element of this Matrix by a scalar (in-place).
		///
		/// @param s Scalar to divide this Matrix by.
		/// @return Reference to this class.
		inline DMatrix<T, rows, columns>& operator/=(const T& s) {
			return (*this) *= (T(1) / s);
		}

		/// @brief Multiply this Matrix with another Matrix (in-place).
		///
		/// @param m Matrix to multiply with this Matrix.
		/// @return Reference to this class.
		inline DMatrix<T, rows, columns>& operator*=(const DMatrix<T, rows, columns>& m) {
			const DMatrix<T, rows, columns> copy_of_this(*this);
			TimesHelper(copy_of_this, m, this);
			return *this;
		}

		/// @brief Multiply a Vector by a Matrix.
		///
		/// @param v Vector to multiply.
		/// @param m Matrix to multiply.
		/// @return Matrix containing the result.
		friend inline DVector<T, columns> operator*(const DVector<T, rows>& v, const DMatrix<T, rows, columns>& m) {
			const int d = columns;
			DAMONSMATH_VECTOR_OPERATOR((DVector<T, rows>::DotProduct(m.data_[i], v)));
		}

		/// @brief Calculate the inverse of this Matrix.
		///
		/// This calculates the inverse Matrix such that
		/// <code>(m * m.Inverse())</code> is the identity.
		/// @return Matrix containing the result.
		inline DMatrix<T, rows, columns> Inverse() const {
			DMatrix<T, rows, columns> inverse;
			InverseHelper<false>(*this, &inverse);
			return inverse;
		}

		/// @brief Calculate the inverse of this Matrix.
		///
		/// This calculates the inverse Matrix such that
		/// <code>(m * m.Inverse())</code> is the identity.
		/// By contrast to Inverse() this returns whether the matrix is invertible.
		///
		/// The invertible check simply compares the calculated determinant with
		/// Constants<T>::GetDeterminantThreshold() to roughly determine whether the
		/// matrix is invertible.  This simple check works in common cases but will
		/// fail for corner cases where the matrix is a combination of huge and tiny
		/// values that can't be accurately represented by the floating point
		/// datatype T.  More extensive checks (relative to the input values) are
		/// possible but <b>far</b> more expensive, complicated and difficult to
		/// test.
		/// @return Whether the matrix is invertible.
		inline bool InverseWithDeterminantCheck(DMatrix<T, rows, columns>* const inverse)const {
			return InverseHelper<true>(*this, inverse);
		}

		/// @brief Calculate the transpose of this Matrix.
		///
		/// @return The transpose of the specified Matrix.
		inline DMatrix<T, columns, rows> Transpose() const {
			DMatrix<T, columns, rows> transpose;

			DAMONSMATH_UNROLLED_LOOP(i, columns, 
			DAMONSMATH_UNROLLED_LOOP(j, rows, transpose.GetColumn(j)[i] = GetColumn(i)[j])
			)

				return transpose;
		}

		/// @brief Get the 2-dimensional translation of a 2-dimensional affine
		/// transform.
		///
		/// @note 2-dimensional affine transforms are represented by 3x3 matrices.
		/// @return Vector with the first two components of column 2 of this Matrix.
		inline DVector<T, 2> TranslationVector2D() const {
			assert(rows == 3 && columns == 3);
			return DVector<T, 2>(data_[2][0], data_[2][1]);
		}

		/// @brief Get the 3-dimensional translation of a 3-dimensional affine
		/// transform.
		///
		/// @note 3-dimensional affine transforms are represented by 4x4 matrices.
		/// @return Vector with the first three components of column 3.
		inline DVector<T, 3> TranslationVector3D() const {
			assert(rows == 4 && columns == 4);
			return DVector<T, 3>(data_[3][0], data_[3][1], data_[3][2]);
		}

		/// @brief Calculate the outer product of two Vectors.
		///
		/// @return Matrix containing the result.
		static inline DMatrix<T, rows, columns> OuterProduct(const DVector<T, rows>& v1, 
															 const DVector<T, columns>& v2) {
			return OuterProductHelper(v1, v2);
		}

		/// @brief Calculate the hadamard / component-wise product of two matrices.
		///
		/// @param m1 First Matrix.
		/// @param m2 Second Matrix.
		/// @return Matrix containing the result.
		static inline DMatrix<T, rows, columns> HadamardProduct(const DMatrix<T, rows, columns>& m1, 
																const DMatrix<T, rows, columns>& m2) {
			DAMONSMATH_MAT_OPERATOR(m1[i] * m2[i]);
		}

		/// @brief Calculate the identity Matrix.
		///
		/// @return Matrix containing the result.
		static inline DMatrix<T, rows, columns> Identity() {
			return IdentityHelper<T, rows, columns>();
		}

		/// @brief Create a 3x3 translation Matrix from a 2-dimensional Vector.
		///
		/// This matrix will have an empty or zero rotation component.
		///
		/// @param v Vector of size 2.
		/// @return Matrix containing the result.
		static inline DMatrix<T, 3> FromTranslationVector(const DVector<T, 2>& v) {
			return DMatrix<T, 3>(1, 0, 0, 0, 1, 0, v[0], v[1], 1);
		}

		/// @brief Create a 4x4 translation Matrix from a 3-dimensional Vector.
		///
		/// This matrix will have an empty or zero rotation component.
		///
		/// @param v The vector of size 3.
		/// @return Matrix containing the result.
		static inline DMatrix<T, 4> FromTranslationVector(const DVector<T, 3>& v) {
			return Matrix<T, 4>(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, v[0], v[1], v[2],
				1);
		}

		/// @brief Create a square Matrix with the diagonal component set to v.
		///
		/// This is an affine transform matrix, so the dimension of the vector is
		/// one less than the dimension of the matrix.
		///
		/// @param v Vector containing components for scaling.
		/// @return Matrix with v along the diagonal, and 1 in the bottom right.
		static inline DMatrix<T, rows> FromScaleVector(const DVector<T, rows - 1>& v) {
			// TODO OPT: Use a helper function in a similar way to Identity to
			// construct the matrix for the specialized cases 2, 3, 4, and only run
			// this method in the general case. This will also allow you to use the
			// helper methods from specialized classes like Matrix<T, 4, 4>.
			Matrix<T, rows> return_matrix(Identity());
			for (int i = 0; i < rows - 1; ++i) 
				return_matrix(i, i) = v[i];

			return return_matrix;
		}

		/// @brief Create a 4x4 Matrix from a 3x3 rotation Matrix.
		///
		/// This Matrix will have an empty or zero translation component.
		///
		/// @param m 3x3 rotation Matrix.
		/// @return Matrix containing the result.
		static inline DMatrix<T, 4> FromRotationMatrix(const DMatrix<T, 3>& m) {
			return DMatrix<T, 4>(m[0], m[1], m[2], 0, m[3], m[4], m[5], 0, m[6], m[7],
				m[8], 0, 0, 0, 0, 1);
		}

		/// @brief Extracts the 3x3 rotation Matrix from a 4x4 Matrix.
		///
		/// This resulting Matrix will contain the upper-left 3x3 sub-matrix of the
		/// input Matrix.
		///
		/// @param m 4x4 Matrix.
		/// @return rotation Matrix containing the result.
		static inline DMatrix<T, 3> ToRotationMatrix(const DMatrix<T, 4>& m) {
			return DMatrix<T, 3>(m[0], m[1], m[2], m[4], m[5], m[6], m[8], m[9],
				m[10]);
		}

		/// @brief Constructs a Matrix<float, 4> from an AffineTransform.
		///
		/// @param affine An AffineTransform reference to be used to construct
		/// a Matrix<float, 4> by adding in the 'w' row of [0, 0, 0, 1].
		static inline DMatrix<T, 4> FromAffineTransform(const DMatrix<T, 4, 3>& affine) {
			return DMatrix<T, 4>(affine[0], affine[4], affine[8], static_cast<T>(0),
				affine[1], affine[5], affine[9], static_cast<T>(0),
				affine[2], affine[6], affine[10], static_cast<T>(0),
				affine[3], affine[7], affine[11], static_cast<T>(1));
		}

		/// @brief Converts a Matrix<float, 4> into an AffineTransform.
		///
		/// @param m A Matrix<float, 4> reference to be converted into an
		/// AffineTransform by dropping the fixed 'w' row.
		///
		/// @return Returns an AffineTransform that contains the essential
		/// transformation data from the Matrix<float, 4>.
		static inline DMatrix<T, 4, 3> ToAffineTransform(const DMatrix<T, 4>& m) {
			return DMatrix<T, 4, 3>(m[0], m[4], m[8], m[12], m[1], m[5], m[9], m[13],
				m[2], m[6], m[10], m[14]);
		}

		/// @brief Create a 3x3 rotation Matrix from a 2D normalized directional
		/// Vector around the X axis.
		///
		/// @param v 2D normalized directional Vector.
		/// @return Matrix containing the result.
		static inline DMatrix<T, 3> RotationX(const DVector<T, 2>& v) {
			return DMatrix<T, 3>(1, 0, 0, 0, v.x, v.y, 0, -v.y, v.x);
		}

		/// @brief Create a 3x3 rotation Matrix from a 2D normalized directional
		/// Vector around the Y axis.
		///
		/// @param v 2D normalized directional Vector.
		/// @return Matrix containing the result.
		static inline DMatrix<T, 3> RotationY(const DVector<T, 2>& v) {
			return DMatrix<T, 3>(v.x, 0, -v.y, 0, 1, 0, v.y, 0, v.x);
		}

		/// @brief Create a 3x3 rotation Matrix from a 2D normalized directional
		/// Vector around the Z axis.
		///
		/// @param v 2D normalized directional Vector.
		/// @return Matrix containing the result.
		static inline DMatrix<T, 3> RotationZ(const DVector<T, 2>& v) {
			return DMatrix<T, 3>(v.x, v.y, 0, -v.y, v.x, 0, 0, 0, 1);
		}

		/// @brief Create a 3x3 rotation Matrix from an angle (in radians) around
		/// the X axis.
		///
		/// @param angle Angle (in radians).
		/// @return Matrix containing the result.
		static inline DMatrix<T, 3> RotationX(T angle) {
			return RotationX(DVector<T, 2>(cosf(angle), sinf(angle)));
		}

		/// @brief Create a 3x3 rotation Matrix from an angle (in radians) around
		/// the Y axis.
		///
		/// @param angle Angle (in radians).
		/// @return Matrix containing the result.
		static inline DMatrix<T, 3> RotationY(T angle) {
			return RotationY(DVector<T, 2>(cosf(angle), sinf(angle)));
		}

		/// @brief Create a 3x3 rotation Matrix from an angle (in radians)
		/// around the Z axis.
		///
		/// @param angle Angle (in radians).
		/// @return Matrix containing the result.
		static inline DMatrix<T, 3> RotationZ(T angle) {
			return RotationZ(DVector<T, 2>(cosf(angle), sinf(angle)));
		}

		/// @brief Create a 4x4 perspective Matrix.
		///
		/// @param fovy Field of view.
		/// @param aspect Aspect ratio.
		/// @param znear Near plane location.
		/// @param zfar Far plane location.
		/// @param handedness 1.0f for RH, -1.0f for LH
		/// @return 4x4 perspective Matrix.
		static inline DMatrix<T, 4, 4> Perspective(T fovy, T aspect, T znear, T zfar,T handedness = 1) {
			return PerspectiveHelper(fovy, aspect, znear, zfar, handedness);
		}

		/// @brief Create a 4x4 orthographic Matrix.
		///
		/// @param left Left extent.
		/// @param right Right extent.
		/// @param bottom Bottom extent.
		/// @param top Top extent.
		/// @param znear Near plane location.
		/// @param zfar Far plane location.
		/// @param handedness 1.0f for RH, -1.0f for LH
		/// @return 4x4 orthographic Matrix.
		static inline DMatrix<T, 4, 4> Ortho(T left, T right, T bottom, T top,
											T znear,T zfar, T handedness = 1) {
			return OrthoHelper(left, right, bottom, top, znear, zfar, handedness);
		}

		/// @brief Create a 3-dimensional camera Matrix.
		///
		/// @param at The look-at target of the camera.
		/// @param eye The position of the camera.
		/// @param up The up vector in the world, for example (0, 1, 0) if the
		/// y-axis is up.
		/// @param handedness 1.0f for RH, -1.0f for LH.
		/// @return 3-dimensional camera Matrix.
		/// TODO: Change default handedness to +1 so that it matches Perspective().
		static inline DMatrix<T, 4, 4> LookAt(const DVector<T, 3>& at,
											  const DVector<T, 3>& eye,
											  const DVector<T, 3>& up,
													T handedness = -1) {
			return LookAtHelper(at, eye, up, handedness);
		}

		/// @brief Get the 3D position in object space from a window coordinate.
		///
		/// @param window_coord The window coordinate. The z value is for depth.
		/// A window coordinate on the near plane will have 0 as the z value.
		/// And a window coordinate on the far plane will have 1 as the z value.
		/// z value should be with in [0, 1] here.
		/// @param model_view The Model View matrix.
		/// @param projection The projection matrix.
		/// @param window_width Width of the window.
		/// @param window_height Height of the window.
		/// @return the mapped 3D position in object space.
		static inline DVector<T, 3> UnProject(const DVector<T, 3>& window_coord,
											  const DMatrix<T, 4, 4>& model_view,
											  const DMatrix<T, 4, 4>& projection,
											  const float window_width,
											  const float window_height) {
			Vector<T, 3> result;
			UnProjectHelper(window_coord, model_view, projection, window_width,
				window_height, result);
			return result;
		}

		// Dimensions of the matrix.
		/// Number of rows in the matrix.
		static const int kRows = rows;
		/// Number of columns in the matrix.
		static const int kColumns = columns;
		/// Total number of elements in the matrix.
		static const int kElements = rows * columns;

		/// @brief get string of vector data seperate by space.
		///  vector data must be basic type like :float int double
		///
		/// @return string of vector data.
		inline std::string ToString() const {

			std::string str = "mat"+std::to_string(rows)+ "X" + std::to_string(columns) + "\n";
			for (int r = 0; r < rows; ++r)
			{
				str += std::string("(");
				 for (int c = 0; c < columns; ++c)
				 {
					 str += (std::to_string((*this)(r, c)) + std::string(" "));
				 }
				 str += std::string(")\n");
			 }

			return str;
		}

		/// @brief get wstring of vector data seperate by space.
		///  vector data must be basic type like :float int double
		///
		/// @return wstring of vector data.
		inline std::wstring ToWString() const {
			std::wstring str = L"mat" + std::to_wstring(rows) + L"X" + std::to_wstring(columns) + L"\n";
			for (int r = 0; r < rows; ++r)
			{
				str += std::wstring(L"(");
					for (int c = 0; c < columns; ++c)
					{
						str += (std::to_wstring((*this)(r, c)) + std::wstring(L" "));
					}
				str += std::wstring(L")\n");
			}


			return str;
		}

		DAMONS_DEFINE_CLASS_SIMD_AWARE_NEW_DELETE

	private:
		DVector<T, rows> data_[columns];
	};


	/// @brief Multiply each element of a Matrix by a scalar.
	///
	/// @param s Scalar to multiply by.
	/// @param m Matrix to multiply.
	/// @return Matrix containing the result.
	/// @tparam T Type of each element in the Matrix and the scalar type.
	/// @tparam rows Number of rows in the matrix.
	/// @tparam columns Number of columns in the matrix.
	///
	/// @related mathfu::Matrix
	template <class T, int rows, int columns>
	inline DMatrix<T, rows, columns> operator*(const T& s,
		const DMatrix<T, columns, rows>& m) {
		return m * s;
	}

	/// @brief Multiply a Matrix by a Vector.
	///
	/// @note Template specialized versions are implemented for 2x2, 3x3, and 4x4
	/// matrices to increase performance.  The 3x3 float is also specialized
	/// to supported padded the 3-dimensional Vector in SIMD build configurations.
	///
	/// @param m Matrix to multiply.
	/// @param v Vector to multiply.
	/// @return Vector containing the result.
	///

	template <class T, int rows, int columns>
	inline DVector<T, rows> operator*(const DMatrix<T, rows, columns>& m,
									  const DVector<T, columns>& v) {
		const DVector<T, rows> result(0);
		int offset = 0;
		for (int column = 0; column < columns; column++) {
			for (int row = 0; row < rows; row++) {
				result[row] += m[offset + row] * v[column];
			}
			offset += rows;
		}
		return result;
	}

	/// @cond MATHFU_INTERNAL
	template <class T>
	inline DVector<T, 2> operator*(const DMatrix<T, 2, 2>& m, const DVector<T, 2>& v) {
		return DVector<T, 2>(m[0] * v[0] + m[2] * v[1], m[1] * v[0] + m[3] * v[1]);
	}
	/// @endcond

	/// @cond MATHFU_INTERNAL
	template <class T>
	inline DVector<T, 3> operator*(const DMatrix<T, 3, 3>& m, const DVector<T, 3>& v) {
		return DVector<T, 3>(DAMONSMATH_MATRIX_3X3_DOT(&m[0], v, 0, 3),
							 DAMONSMATH_MATRIX_3X3_DOT(&m[0], v, 1, 3),
							 DAMONSMATH_MATRIX_3X3_DOT(&m[0], v, 2, 3));
	}
	/// @endcond

	/// @cond MATHFU_INTERNAL
	template <>
	inline DVector<float, 3> operator*(const DMatrix<float, 3, 3>& m,
		const DVector<float, 3>& v) {
		return DVector<float, 3>(
			DAMONSMATH_MATRIX_3X3_DOT(&m[0], v, 0, DAMONSMATH_VECTOR_STRIDE_FLOATS(v)),
			DAMONSMATH_MATRIX_3X3_DOT(&m[0], v, 1, DAMONSMATH_VECTOR_STRIDE_FLOATS(v)),
			DAMONSMATH_MATRIX_3X3_DOT(&m[0], v, 2, DAMONSMATH_VECTOR_STRIDE_FLOATS(v)));
	}
	/// @endcond

	/// @cond MATHFU_INTERNAL
	template <class T>
	inline DVector<T, 4> operator*(const DMatrix<T, 4, 4>& m, const DVector<T, 4>& v) {
		return DVector<T, 4>(
			DAMONSMATH_MATRIX_4X4_DOT(&m[0], v, 0), DAMONSMATH_MATRIX_4X4_DOT(&m[0], v, 1),
			DAMONSMATH_MATRIX_4X4_DOT(&m[0], v, 2), DAMONSMATH_MATRIX_4X4_DOT(&m[0], v, 3));
	}
	/// @endcond

	/// @brief Multiply a 4x4 Matrix by a 3-dimensional Vector.
	///
	/// This is provided as a convenience and assumes the vector has a fourth
	/// component equal to 1.
	///
	/// @param m 4x4 Matrix.
	/// @param v 3-dimensional Vector.
	/// @return 3-dimensional Vector result.
	///
	/// @related mathfu::Matrix
	template <class T>
	inline DVector<T, 3> operator*(const DMatrix<T, 4, 4>& m, const DVector<T, 3>& v) {
		DVector<T, 4> v4(v[0], v[1], v[2], 1);
		v4 = m * v4;
		return DVector<T, 3>(v4[0] / v4[3], v4[1] / v4[3], v4[2] / v4[3]);
	}

	/// @brief Multiply a Matrix with another Matrix.
	///
	/// @note Template specialized versions are implemented for 2x2, 3x3, and 4x4
	/// matrices to improve performance. 3x3 float is also specialized because if
	/// SIMD is used the vectors of this type of length 4.
	///
	/// @param m1 Matrix to multiply.
	/// @param m2 Matrix to multiply.
	/// @param out_m Pointer to a Matrix which receives the result.
	///
	/// @tparam T Type of each element in the returned Matrix.
	/// @tparam size1 Number of rows in the returned Matrix and columns in m1.
	/// @tparam size2 Number of columns in the returned Matrix and rows in m2.
	/// @tparam size3 Number of columns in m3.
	template <class T, int size1, int size2, int size3>
	inline void TimesHelper(const DMatrix<T, size1, size2>& m1,
		const DMatrix<T, size2, size3>& m2,
		DMatrix<T, size1, size3>* out_m) {
		for (int i = 0; i < size1; i++) {
			for (int j = 0; j < size3; j++) {
				DVector<T, size2> row;
				for (int k = 0; k < size2; k++) {
					row[k] = m1(i, k);
				}
				(*out_m)(i, j) = DVector<T, size2>::DotProduct(m2.GetColumn(j), row);
			}
		}
	}

	/// @cond MATHFU_INTERNAL
	template <class T>
	inline void TimesHelper(const DMatrix<T, 2, 2>& m1, const DMatrix<T, 2, 2>& m2,
		DMatrix<T, 2, 2>* out_m) {
		DMatrix<T, 2, 2>& out = *out_m;
		out[0] = m1[0] * m2[0] + m1[2] * m2[1];
		out[1] = m1[1] * m2[0] + m1[3] * m2[1];
		out[2] = m1[0] * m2[2] + m1[2] * m2[3];
		out[3] = m1[1] * m2[2] + m1[3] * m2[3];
	}
	/// @endcond

	/// @cond MATHFU_INTERNAL
	template <typename T>
	inline void TimesHelper(const DMatrix<T, 3, 3>& m1, const DMatrix<T, 3, 3>& m2,
		DMatrix<T, 3, 3>* out_m) {
		DMatrix<T, 3, 3>& out = *out_m;
		{
			DVector<T, 3> row(m1[0], m1[3], m1[6]);
			out[0] = DVector<T, 3>::DotProduct(m2.GetColumn(0), row);
			out[3] = DVector<T, 3>::DotProduct(m2.GetColumn(1), row);
			out[6] = DVector<T, 3>::DotProduct(m2.GetColumn(2), row);
		}
		{
			DVector<T, 3> row(m1[1], m1[4], m1[7]);
			out[1] = DVector<T, 3>::DotProduct(m2.GetColumn(0), row);
			out[4] = DVector<T, 3>::DotProduct(m2.GetColumn(1), row);
			out[7] = DVector<T, 3>::DotProduct(m2.GetColumn(2), row);
		}
		{
			DVector<T, 3> row(m1[2], m1[5], m1[8]);
			out[2] = DVector<T, 3>::DotProduct(m2.GetColumn(0), row);
			out[5] = DVector<T, 3>::DotProduct(m2.GetColumn(1), row);
			out[8] = DVector<T, 3>::DotProduct(m2.GetColumn(2), row);
		}
	}
	/// @endcond

	/// @cond MATHFU_INTERNAL
	template <class T>
	inline void TimesHelper(const DMatrix<T, 4, 4>& m1, const DMatrix<T, 4, 4>& m2,
		DMatrix<T, 4, 4>* out_m) {
		DMatrix<T, 4, 4>& out = *out_m;
		{
			Vector<T, 4> row(m1[0], m1[4], m1[8], m1[12]);
			out[0]  = DVector<T, 4>::DotProduct(m2.GetColumn(0), row);
			out[4]  = DVector<T, 4>::DotProduct(m2.GetColumn(1), row);
			out[8]  = DVector<T, 4>::DotProduct(m2.GetColumn(2), row);
			out[12] = DVector<T, 4>::DotProduct(m2.GetColumn(3), row);
		}
		{
			Vector<T, 4> row(m1[1], m1[5], m1[9], m1[13]);
			out[1]  = DVector<T, 4>::DotProduct(m2.GetColumn(0), row);
			out[5]  = DVector<T, 4>::DotProduct(m2.GetColumn(1), row);
			out[9]  = DVector<T, 4>::DotProduct(m2.GetColumn(2), row);
			out[13] = DVector<T, 4>::DotProduct(m2.GetColumn(3), row);
		}
		{
			Vector<T, 4> row(m1[2], m1[6], m1[10], m1[14]);
			out[2]  = DVector<T, 4>::DotProduct(m2.GetColumn(0), row);
			out[6]  = DVector<T, 4>::DotProduct(m2.GetColumn(1), row);
			out[10] = DVector<T, 4>::DotProduct(m2.GetColumn(2), row);
			out[14] = DVector<T, 4>::DotProduct(m2.GetColumn(3), row);
		}
		{
			Vector<T, 4> row(m1[3], m1[7], m1[11], m1[15]);
			out[3]  = DVector<T, 4>::DotProduct(m2.GetColumn(0), row);
			out[7]  = DVector<T, 4>::DotProduct(m2.GetColumn(1), row);
			out[11] = DVector<T, 4>::DotProduct(m2.GetColumn(2), row);
			out[15] = DVector<T, 4>::DotProduct(m2.GetColumn(3), row);
		}
	}
	/// @endcond

	/// @brief Compute the identity matrix.
	///
	/// @note There are template specializations for 2x2, 3x3, and 4x4 matrices to
	/// increase performance.
	///
	/// @return Identity Matrix.
	/// @tparam T Type of each element in the returned Matrix.
	/// @tparam rows Number of rows in the returned Matrix.
	/// @tparam columns Number of columns in the returned Matrix.
	template <class T, int rows, int columns>
	inline DMatrix<T, rows, columns> IdentityHelper() {
		DMatrix<T, rows, columns> return_matrix(0.f);
		int min_d = rows < columns ? rows : columns;
		for (int i = 0; i < min_d; ++i) return_matrix(i, i) = 1;
		return return_matrix;
	}
	/// @endcond

	/// @cond MATHFU_INTERNAL
	template <class T>
	inline DMatrix<T, 2, 2> IdentityHelper() {
		return DMatrix<T, 2, 2>(1, 0, 0, 1);
	}
	/// @endcond

	/// @cond MATHFU_INTERNAL
	template <class T>
	inline DMatrix<T, 3, 3> IdentityHelper() {
		return DMatrix<T, 3, 3>(1, 0, 0, 0, 1, 0, 0, 0, 1);
	}
	/// @endcond

	/// @cond MATHFU_INTERNAL
	template <class T>
	inline DMatrix<T, 4, 4> IdentityHelper() {
		return DMatrix<T, 4, 4>(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
	}
	/// @endcond

	/// @brief Compute the outer product of two vectors.
	///
	/// @note There are template specialization for 2x2, 3x3, and 4x4 matrices to
	/// increase performance.
	template <class T, int rows, int columns>
	static inline DMatrix<T, rows, columns> OuterProductHelper(
		const DVector<T, rows>& v1, const DVector<T, columns>& v2) {
		DMatrix<T, rows, columns> result(0);
		int offset = 0;
		for (int column = 0; column < columns; column++) {
			for (int row = 0; row < rows; row++) {
				result[row + offset] = v1[row] * v2[column];
			}
			offset += rows;
		}
		return result;
	}
	/// @endcond

	/// @cond MATHFU_INTERNAL
	template <class T>
	static inline DMatrix<T, 2, 2> OuterProductHelper(const DVector<T, 2>& v1,
													  const DVector<T, 2>& v2) {
		return DMatrix<T, 2, 2>(v1[0] * v2[0], v1[1] * v2[0], v1[0] * v2[1],
			v1[1] * v2[1]);
	}
	/// @endcond

	/// @cond MATHFU_INTERNAL
	template <class T>
	static inline DMatrix<T, 3, 3> OuterProductHelper(const DVector<T, 3>& v1,
													  const DVector<T, 3>& v2) {
		return DMatrix<T, 3, 3>(v1[0] * v2[0], v1[1] * v2[0], v1[2] * v2[0],
			v1[0] * v2[1], v1[1] * v2[1], v1[2] * v2[1],
			v1[0] * v2[2], v1[1] * v2[2], v1[2] * v2[2]);
	}
	/// @endcond

	/// @cond MATHFU_INTERNAL
	template <class T>
	static inline DMatrix<T, 4, 4> OuterProductHelper(const DVector<T, 4>& v1,
		const DVector<T, 4>& v2) {
		return DMatrix<T, 4, 4>(
			v1[0] * v2[0], v1[1] * v2[0], v1[2] * v2[0], v1[3] * v2[0], v1[0] * v2[1],
			v1[1] * v2[1], v1[2] * v2[1], v1[3] * v2[1], v1[0] * v2[2], v1[1] * v2[2],
			v1[2] * v2[2], v1[3] * v2[2], v1[0] * v2[3], v1[1] * v2[3], v1[2] * v2[3],
			v1[3] * v2[3]);
	}
	/// @endcond

	/// Struct used for template specialization for functions that
	/// returns constants.
	template <class T>
	class Constants {
	public:
		/// Minimum absolute value of the determinant of an invertible matrix.
		static T GetDeterminantThreshold() {
			// No constant defined for the general case.
			assert(false);
			return 0;
		}
	};
	/// @endcond

	/// Functions that return constants for <code>float</code> values.
	template <>
	class Constants<float> {
	public:
		/// @brief Minimum absolute value of the determinant of an invertible
		/// <code>float</code> Matrix.
		///
		/// <code>float</code> values have 23 bits of precision which is roughly
		/// 1e7f, given that the final step of matrix inversion is multiplication
		/// with the inverse of the determinant, the minimum value of the
		/// determinant is 1e-7f before the precision too low to accurately
		/// calculate the inverse.
		/// @returns Minimum absolute value of the determinant of an invertible
		/// <code>float</code> Matrix.
		///
		/// @related mathfu::Matrix::InverseWithDeterminantCheck()
		static float GetDeterminantThreshold() { return 1e-7f; }
	};

	/// Functions that return constants for <code>double</code> values.
	template <>
	class Constants<double> {
	public:
		/// @brief Minimum absolute value of the determinant of an invertible
		/// <code>double</code> Matrix.
		///
		/// <code>double</code> values have 46 bits of precision which is roughly
		/// 1e15, given that the final step of matrix inversion is multiplication
		/// with the inverse of the determinant, the minimum value of the
		/// determinant is 1e-15 before the precision too low to accurately
		/// calculate the inverse.
		/// @returns Minimum absolute value of the determinant of an invertible
		/// <code>double</code> Matrix.
		///
		/// @related mathfu::Matrix::InverseWithDeterminantCheck()
		static double GetDeterminantThreshold() { return 1e-15; }
	};

	/// @brief Compute the inverse of a matrix.
	///
	/// There is template specialization  for 2x2, 3x3, and 4x4 matrices to
	/// increase performance. Inverse is not implemented for dense matrices that
	/// are not of size 2x2, 3x3, and 4x4.  If check_invertible is true the
	/// determine of the matrix is compared with
	/// Constants<T>::GetDeterminantThreshold() to roughly determine whether the
	/// Matrix is invertible.
	template <bool check_invertible, class T, int rows, int columns>
	inline bool InverseHelper(const DMatrix<T, rows, columns>& m,
							 DMatrix<T, rows, columns>* const inverse) {
		assert(false);
		(void)m;
		*inverse = T::Identity();
		return false;
	}
	/// @endcond

	/// @cond MATHFU_INTERNAL
	template <bool check_invertible, class T>
	inline bool InverseHelper(const DMatrix<T, 2, 2>& m,
							  DMatrix<T, 2, 2>* const inverse) {
		T determinant = m[0] * m[3] - m[1] * m[2];
		if (check_invertible &&
			fabs(determinant) < Constants<T>::GetDeterminantThreshold()) {
			return false;
		}
		T inverseDeterminant = 1 / determinant;
		(*inverse)[0] = inverseDeterminant * m[3];
		(*inverse)[1] = -inverseDeterminant * m[1];
		(*inverse)[2] = -inverseDeterminant * m[2];
		(*inverse)[3] = inverseDeterminant * m[0];
		return true;
	}
	/// @endcond

	/// @cond MATHFU_INTERNAL
	template <bool check_invertible, class T>
	inline bool InverseHelper(const DMatrix<T, 3, 3>& m,
		DMatrix<T, 3, 3>* const inverse) {
		// Find determinant of matrix.
		T sub11 =  m[4] * m[8] - m[5] * m[7];
		T sub12 = -m[1] * m[8] + m[2] * m[7];
		T sub13 =  m[1] * m[5] - m[2] * m[4];
		T determinant = m[0] * sub11 + m[3] * sub12 + m[6] * sub13;
		if (check_invertible &&
			fabs(determinant) < Constants<T>::GetDeterminantThreshold()) {
			return false;
		}
		// Find determinants of 2x2 submatrices for the elements of the inverse.
		*inverse = DMatrix<T, 3, 3>(
			sub11, sub12, sub13, m[6] * m[5] - m[3] * m[8], m[0] * m[8] - m[6] * m[2],
			m[3] * m[2] - m[0] * m[5], m[3] * m[7] - m[6] * m[4],
			m[6] * m[1] - m[0] * m[7], m[0] * m[4] - m[3] * m[1]);
		*(inverse) *= 1 / determinant;
		return true;
	}
	/// @endcond

	template <class T>
	inline int FindLargestPivotElem(const DMatrix<T, 4, 4>& m) {
		DVector<T, 4> fabs_column(fabs(m[0]), fabs(m[1]), fabs(m[2]), fabs(m[3]));
		if (fabs_column[0] > fabs_column[1]) {
			if (fabs_column[0] > fabs_column[2]) {
				if (fabs_column[0] > fabs_column[3]) {
					return 0;
				}
				else {
					return 3;
				}
			}
			else if (fabs_column[2] > fabs_column[3]) {
				return 2;
			}
			else {
				return 3;
			}
		}
		else if (fabs_column[1] > fabs_column[2]) {
			if (fabs_column[1] > fabs_column[3]) {
				return 1;
			}
			else {
				return 3;
			}
		}
		else if (fabs_column[2] > fabs_column[3]) {
			return 2;
		}
		else {
			return 3;
		}
	}
	/// @endcond

	/// @cond MATHFU_INTERNAL
	template <bool check_invertible, class T>
	bool InverseHelper(const DMatrix<T, 4, 4>& m, DMatrix<T, 4, 4>* const inverse) {
		// This will find the pivot element.
		int pivot_elem = FindLargestPivotElem(m);
		// This will perform the pivot and find the row, column, and 3x3 submatrix
		// for this pivot.
		DVector<T, 3> row, column;
		DMatrix<T, 3> matrix;
		if (pivot_elem == 0) {
			row    = DVector<T, 3>(m[4], m[8], m[12]);
			column = DVector<T, 3>(m[1], m[2], m[3]);
			matrix = DMatrix<T, 3>(m[5], m[6], m[7], m[9], m[10], m[11], m[13], m[14], m[15]);
		}
		else if (pivot_elem == 1) {
			row    = DVector<T, 3>(m[5], m[9], m[13]);
			column = DVector<T, 3>(m[0], m[2], m[3]);
			matrix = DMatrix<T, 3>(m[4], m[6], m[7], m[8], m[10], m[11], m[12], m[14], m[15]);
		}
		else if (pivot_elem == 2) {
			row    = DVector<T, 3>(m[6], m[10], m[14]);
			column = DVector<T, 3>(m[0], m[1], m[3]);
			matrix = DMatrix<T, 3>(m[4], m[5], m[7], m[8], m[9], m[11], m[12], m[13], m[15]);
		}
		else {
			row    = DVector<T, 3>(m[7], m[11], m[15]);
			column = DVector<T, 3>(m[0], m[1], m[2]);
			matrix = DMatrix<T, 3>(m[4], m[5], m[6], m[8], m[9], m[10], m[12], m[13], m[14]);
		}
		T pivot_value = m[pivot_elem];
		if (check_invertible &&
			fabs(pivot_value) < Constants<T>::GetDeterminantThreshold()) {
			return false;
		}
		// This will compute the inverse using the row, column, and 3x3 submatrix.
		T inv = -1 / pivot_value;
		row *= inv;
		matrix += DMatrix<T, 3>::OuterProduct(column, row);
		DMatrix<T, 3> mat_inverse;
		if (!InverseHelper<check_invertible>(matrix, &mat_inverse) &&
			check_invertible) {
			return false;
		}
		DVector<T, 3> col_inverse = mat_inverse * (column * inv);
		DVector<T, 3> row_inverse = row * mat_inverse;
		T pivot_inverse = DVector<T, 3>::DotProduct(row, col_inverse) - inv;
		if (pivot_elem == 0) {
			*inverse = DMatrix<T, 4, 4>(
				pivot_inverse, col_inverse[0], col_inverse[1], col_inverse[2],
				row_inverse[0], mat_inverse[0], mat_inverse[1], mat_inverse[2],
				row_inverse[1], mat_inverse[3], mat_inverse[4], mat_inverse[5],
				row_inverse[2], mat_inverse[6], mat_inverse[7], mat_inverse[8]);
		}
		else if (pivot_elem == 1) {
			*inverse = DMatrix<T, 4, 4>(
				row_inverse[0], mat_inverse[0], mat_inverse[1], mat_inverse[2],
				pivot_inverse, col_inverse[0], col_inverse[1], col_inverse[2],
				row_inverse[1], mat_inverse[3], mat_inverse[4], mat_inverse[5],
				row_inverse[2], mat_inverse[6], mat_inverse[7], mat_inverse[8]);
		}
		else if (pivot_elem == 2) {
			*inverse = DMatrix<T, 4, 4>(
				row_inverse[0], mat_inverse[0], mat_inverse[1], mat_inverse[2],
				row_inverse[1], mat_inverse[3], mat_inverse[4], mat_inverse[5],
				pivot_inverse, col_inverse[0], col_inverse[1], col_inverse[2],
				row_inverse[2], mat_inverse[6], mat_inverse[7], mat_inverse[8]);
		}
		else {
			*inverse = DMatrix<T, 4, 4>(
				row_inverse[0], mat_inverse[0], mat_inverse[1], mat_inverse[2],
				row_inverse[1], mat_inverse[3], mat_inverse[4], mat_inverse[5],
				row_inverse[2], mat_inverse[6], mat_inverse[7], mat_inverse[8],
				pivot_inverse, col_inverse[0], col_inverse[1], col_inverse[2]);
		}
		return true;
	}

	/// Create a 4x4 perpective matrix.
	template <class T>
	inline DMatrix<T, 4, 4> PerspectiveHelper(T fovy, T aspect, T znear, T zfar,T handedness) {
		const T y = 1 / std::tan(fovy * static_cast<T>(.5));
		const T x = y / aspect;
		const T zdist = (znear - zfar);
		const T zfar_per_zdist = zfar / zdist;
		return DMatrix<T, 4, 4>(x, 0, 0, 0, 0, y, 0, 0, 0, 0,
			zfar_per_zdist * handedness, -1 * handedness, 0, 0,
			2.0f * znear * zfar_per_zdist, 0);
	}
	/// @endcond

	/// @cond MATHFU_INTERNAL
	/// Create a 4x4 orthographic matrix.
	template <class T>
	static inline DMatrix<T, 4, 4> OrthoHelper( T left, T right, T bottom, T top,
												T znear, T zfar, T handedness) {
		return DMatrix<T, 4, 4>(static_cast<T>(2) / (right - left), 0, 0, 0, 0,
			static_cast<T>(2) / (top - bottom), 0, 0, 0, 0,
			-handedness * static_cast<T>(2) / (zfar - znear), 0,
			-(right + left) / (right - left),
			-(top + bottom) / (top - bottom),
			-(zfar + znear) / (zfar - znear), static_cast<T>(1));
	}
	/// @endcond

	/// @cond MATHFU_INTERNAL
	/// Calculate the axes required to construct a 3-dimensional camera matrix that
	/// looks at "at" from eye position "eye" with the up vector "up".  The axes
	/// are returned in a 4 element "axes" array.
	template <class T>
	static void LookAtHelperCalculateAxes(const DVector<T, 3>& at,const DVector<T, 3>& eye,
										 const DVector<T, 3>& up, T handedness,
										 DVector<T, 3>* const axes) {
		// Notice that y-axis is always the same regardless of handedness.
		axes[2] = (at - eye).Normalized();
		axes[0] = DVector<T, 3>::CrossProduct(up, axes[2]).Normalized();
		axes[1] = DVector<T, 3>::CrossProduct(axes[2], axes[0]);
		axes[3] = DVector<T, 3>(handedness * DVector<T, 3>::DotProduct(axes[0], eye),
			-DVector<T, 3>::DotProduct(axes[1], eye),
			handedness * DVector<T, 3>::DotProduct(axes[2], eye));

		// Default calculation is left-handed (i.e. handedness=-1).
		// Negate x and z axes for right-handed (i.e. handedness=+1) case.
		const T neg = -handedness;
		axes[0] *= neg;
		axes[2] *= neg;
	}
	/// @endcond

	/// @cond MATHFU_INTERNAL
	/// Create a 3-dimensional camera matrix.
	template <class T>
	static inline DMatrix<T, 4, 4> LookAtHelper(const DVector<T, 3>& at,
		const DVector<T, 3>& eye,
		const DVector<T, 3>& up,
		T handedness) {
		DVector<T, 3> axes[4];
		LookAtHelperCalculateAxes(at, eye, up, handedness, axes);
		const DVector<T, 4> column0(axes[0][0], axes[1][0], axes[2][0], 0);
		const DVector<T, 4> column1(axes[0][1], axes[1][1], axes[2][1], 0);
		const DVector<T, 4> column2(axes[0][2], axes[1][2], axes[2][2], 0);
		const DVector<T, 4> column3(axes[3], 1);
		return DMatrix<T, 4, 4>(column0, column1, column2, column3);
	}
	/// @endcond

	/// @cond MATHFU_INTERNAL
	/// Get the 3D position in object space from a window coordinate.
	template <class T>
	static inline bool UnProjectHelper( const DVector<T, 3>& window_coord,
										const DMatrix<T, 4, 4>& model_view,
										const DMatrix<T, 4, 4>& projection,
										const float window_width,
										const float window_height,
										DVector<T, 3>& result) {
		if (window_coord.z < static_cast<T>(0) ||
			window_coord.z > static_cast<T>(1)) {
			// window_coord.z should be with in [0, 1]
			// 0: near plane
			// 1: far plane
			return false;
		}
		DMatrix<T, 4, 4> matrix = (projection * model_view).Inverse();
		DVector<T, 4> standardized = DVector<T, 4>(
			static_cast<T>(2) * (window_coord.x - window_width) / window_width +
			static_cast<T>(1),
			static_cast<T>(2) * (window_coord.y - window_height) / window_height +
			static_cast<T>(1),
			static_cast<T>(2) * window_coord.z - static_cast<T>(1),
			static_cast<T>(1));

		DVector<T, 4> multiply = matrix * standardized;
		if (multiply.w == static_cast<T>(0)) {
			return false;
		}
		result = multiply.xyz() / multiply.w;
		return true;
	}
	/// @endcond


}; /// namespace DMath end

#endif // !_DAMONS_MATRIX_H_

