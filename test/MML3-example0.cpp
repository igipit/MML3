// MML3-example1.cpp : an example on usng MML3
//


#include"MML3-Matrix.h"
#include"MML3-Timer.h"
#include"MML3-Math.h"
#include<iostream>
#include<iomanip>
#include"MML3-MatrixAlgorithms.h"



using MML3::iSet;
using MML3::Mat;
typedef MML3::index_type index_t;

typedef Mat<>::type 	 Matrix;
// equivalent to the finer but longer
//  typedef Mat<double>::General::Rectangular::RowMajor::type 	 Matrix;
// or, directly , to
//  typedef MML3::Matrix<double, MML3::M_PROP::GE, MML3::M_SHAPE::RE, MML3::M_ORD::ROW> Matrix;


typedef Mat<>::sym_type   symMatrix;
// equivalent to
//typedef MML3::Mat<double>::Symmetric::Rectangular::RowMajor::type 	 Matrix;



int main()
{

	// demo: Matrix CTOR's
	{
		// Matrix of 7x6 doubles 
		Matrix B(7, 6);
		// Accessing and assigning some components !!! INDICES ARE 1-BASED !!!
		B(1, 1) = 1.0;
		B(3, 5) = 3.5;
		std::cout << "B:" << std::scientific << std::setw(20) << std::setprecision(12) << B;

		// matrix 5x6 initialized with unitary components
		Matrix K(5, 6);
		K.fill(3.14);
		std::cout << "K:" << std::scientific << std::setw(10) << std::setprecision(3) << K;

		// resizing a matrix and filling with a desired value
		K.resize(4, 4).fill(0.0);

		// matrix 8x3 initialized  with  values assigned by a C like constant array. Note that the unspecified components remain uninitialized
		// the matrix will have a number of columns equal to the maximum number of columns in the rows of the constant C array
		Matrix A = {
			{ 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0 },
			{ 21, 22, 23, 24, 25, 26, 27, 28 },
			{ 31, 32 } /* components {3,...,8 } are left unitialized*/
		};
		std::cout << "A:" << std::scientific << std::setw(10) << std::setprecision(1) << A;


		// matrix initialized  with serial data pointed by p; note that it isn't necessary that p points to the exact type of the matrix template
		// but if the types are equal, initialization is more efficient (via memcpy)
		double p[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
		Matrix C(p,10, 5, 2);
		std::cout << "C:" << std::setw(8) << std::setprecision(2) << C << std::endl;
	}

	{
		// sets the matrix to the  5x5 identity
		Matrix K=Matrix::get_eye(5);
		std::cout << "K:" << std::scientific << std::setw(10) << std::setprecision(1) << K;
		//  transpose the matrix
		Matrix A(3,5);
		Matrix B=transpose(A);
		//std::cout << "A:" << std::scientific << std::setw(10) << std::setprecision(1) << A;
		// A is transformed into a 4x5 random Matrix with components in the range [1.3 - 0.2,  1.3 +0.2]
		A = Matrix::get_random_matrix(4, 5, 1.3, 0.2);

		auto X = Matrix::get_random_matrix(4, 4, 1.0, 2.0);
		std::cout << std::setw(12) << std::setprecision(3) << X << std::endl;
		// creating a matrix by applying a function to another matrix
		auto Z = X.apply(sin);
		std::cout << std::setw(12) << std::setprecision(3) << Z << std::flush;
		// or applying a function to tranform the a matrix
		X.transform(cos);
		std::cout << std::setw(12) << std::setprecision(3) << X << std::flush;

	}


	// Demo: accessing matrices with SubMatrix
	// A SubMatrix is a view of the original matrix : modifying the SubMatrix  we modify the original matrix(only the Submatrix components)
	{


		// construct a 6x7 matrix filled with random numbers in the range [1.3 - 0.2,  1.3 +0.2]
		Matrix K = Matrix::get_random_matrix(6, 7, 1.0, 0.2);

		// a SubMatrix formed by the columns {2,3} of  K
		auto c23 = K.column({ 2, 3 });

		// setting all the SubMatrix components to a specified value
		c23 = 2.0;
		// changes the components of the referred matrix K
		std::cout << "K:" << std::scientific << std::setw(10) << std::setprecision(12) << K;

		// a SubMatrix formed by row 3 of K
		K.row({ 3 }).fill(3.0);
		std::cout << "K:" << std::scientific << std::setw(20) << std::setprecision(12) << K;

		// a SubMatrix formed by rows 2,3 and columns 4,5,6 of K
		auto ssK = K({ 2, 3 }, { 4, 5, 6 });

		auto KT = transpose(K);

		std::cout << "K':" << std::scientific << std::setw(20) << std::setprecision(12) << KT;

		// more generally sub matrices are created passing two set of row and column indexes  to K by means of the operator ()
		// a set of indexes is represented by objects of type iSet or convertible to an iSet
		//-------------------------------------------------
		// Let us explore index sets construction ...

		// creates a set of indices (not necessarily orderd)
		iSet   a({ 3, 1, 2 });
		// or, equivalently
		int q[] = { 3, 1, 2 }; iSet   a_bis(q, 3);
		// or, equivalently
		iSet a_tris; a.resize(3); a(1) = 3; a(2) = 1; a(3) = 2;
		// or equivalently with any std::container (vector,maps,lists, stacks, hash...)
		std::vector<MML3::index_type> a_v({ 3, 1, 2 });
		iSet a_quad(a_v);
		std::cout << a << std::endl;
		// creates an empty  index set
		iSet   b;
		// fills the index set one component at a time increasing dynamically the size
		b.add(4).add(5);
		std::cout << b << std::endl;

		// since iSet stores the indexes  consecutively, the use of add() causes reallocation. To prevent reallocation
		// you can reserve an amount of memory :
		b.reserve(2).add(4).add(5);
		// or, alternatively you can resize the iSet and then assign the components
		b.resize(2); b(1) = 4; b(2) = 5;

		// creates the index set {1,2,3,4,5}
		iSet   c(1, 5); // equivalent to iSet c=iSet::range(1,5);
		std::cout << c << std::endl;

		// creates the set of  indexes  {1,2,...,100}
		iSet d(1, 100);
		// modifies manually the components
		d(5) = 7; d(7) = 5; //...


		// a SubMatrix can be changed by assignement to another SubMatrix
		auto ssK1 = K(b, b);
		ssK1.fill(2.222);
		std::cout << K << std::endl;


		// full impementation  of const Matrix sub matrices (const_SubMatrix)
		const Matrix& cK = K;
		auto crK = cK(b, b);

	}


	// iterating over matrix components
	{
		
		/////////////////////////
		// standard iteration

		
		std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);

		// on a General Rectangular matrix
		// is the same as MML3::Mat<double>::RE
		typedef MML3::Mat<double>::type  Matrix;
		int nr = 4,  nc = 7;
		Matrix K(nr, nc);
		for (Matrix::index_t r = 1; r != nr + 1; ++r)
		for (Matrix::index_t c = 1; c != nc + 1; ++c)
			K(r, c) = r*c;
		std::cout <<"K: " << std::setw(6) << std::setprecision(1) << K << std::endl;


		// only for Row major rectangular matrices there is also the standard C 0-based direct access
		nr = 7;
		nc = 3;
		MML3::Mat<double>::General::type A(nr, nc);
		for (index_t i = 0; i != nr; ++i)
		for (index_t j = 0; j != nc; ++j)
			A[i][j] = (i + 1)*(j + 1);

		// on a General Lower triangular   matrix
		nr = nc = 6;
		typedef MML3::Mat<double>::General::LowerTria::type ltMatrix;
		ltMatrix K1(nr, nr);
		for (Matrix::index_t r = 1; r != nr + 1; ++r)
		for (Matrix::index_t c = 1; c != r + 1; ++c)
			K1(r, c) = r*c;
		std::cout << "K1: " << std::setw(6) << std::setprecision(1)  << K1 << std::endl;


		// on a General Upper triangular   matrix
		nr = nc = 3;
		typedef MML3::Mat<double>::General::UpperTria::type  utMatrix;
		utMatrix K2(nr, nr);
		for (Matrix::index_t r = 1; r != nr + 1; ++r)
		for (Matrix::index_t c = r; c != nr + 1; ++c)
			K2(r, c) = r*c;
		std::cout << "K2: " << std::setw(6) << std::setprecision(1)  << K2 << std::endl;


		// to iterate over the proper column range for every matrix type:
		for (Matrix::index_t r = 1; r != nr+1; ++r)
		for (Matrix::index_t c = K2.col_beg_idx(r); c != K2.col_end_idx(r); ++c)
			K2(r, c) = r*c*c;
		std::cout << "K: " << std::setw(6) << std::setprecision(1) << K << std::endl;
		
		// or ( best practice)
		for (auto r : K2.row_range())
		for (auto c : K2.col_range(r)) // dont forget the argument r in the inner loop since for same matrices it depends on the row
			K2(r, c) = r*c*c;
		std::cout << "K2: " << std::setw(6) << std::setprecision(1) << K2 << std::endl;

	}

	

return 0;
}

