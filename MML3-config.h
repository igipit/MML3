#pragma once

#include<cstdint>
#include<limits>
#include<complex>


// #undef MML3_LAPACK  if you dont want lapack support
// else define as
// 1: intel MKL LAPACKE interface
// 2: netlib    LAPACKE interface
#define MML3_LAPACK 1
#define MML3_LAPACK_INT_TYPE  int



#define MML3_INT_TYPE	 std::int32_t

// the type used for matrix indexing
#define MML3_INDEX_TYPE  std::int32_t

// the triangle in which the data of a rectangular symmetric matrix are stores 'U' for upper triangle
//'L' for Lower triangle
#define MML3_RE_SYM_MAT_PART 'U'

// the first index of a matrix
#define MML3_BASE_INDEX_OFFSET 1


#define MML3_TEST_INDEX_ON_ACCESS

// when this macro is defined access of matrices by indexes is tested
#if defined(_DEBUG) && !defined(MML3_TEST_INDEX_ON_ACCESS)
	#define MML3_TEST_INDEX_ON_ACCESS
#endif


namespace MML3
{
   //the type of dense matrix indexes, used to access the matric components
	typedef		MML3_INDEX_TYPE				index_type;
	typedef     MML3_LAPACK_INT_TYPE		lapack_int_t;
	typedef     MML3_INT_TYPE				int_t;


	// Cblas standardized options
	struct Option
	{
		//typedef unsigned char option_t;
		enum Order		{ RowMajor = 101, ColMajor = 102 };
		enum UpLo		{ Upper = 121, Lower = 122 };
		enum Transpose	{ NoTrans = 111, Trans = 112, ConjTrans = 113 };
		enum Side		{ Left = 141, Right = 142 };
	};

		

   // the type of the complex numbers
   template<typename T>
   using complex = std::complex<T>;



   struct M_PROP {
	   struct GE
	   {
		   enum{
			   ID = 0,
			   SHAPE_GE = 1,
			   PROP_HER = 0
		   };
		   static const char* name(){ return  "GEneral"; }

	   };
	   struct SYM
	   {
		   enum{
			   ID = 1,
			   SHAPE_GE = 0,
			   PROP_HER = 0
		   };
		   static const char* name(){ return  "SYMmetric"; }

	   };
	   struct HER
	   {
		   enum{
			   ID = 2,
			   SHAPE_GE = 0,
			   PROP_HER = 0
		   };
		   static const char* name(){ return  "HERmitian"; }

	   };


   };


   struct M_SHAPE {
	   struct RE
	   {

		   static const char* name(){ return  "REctangular"; }

#if MML3_RE_SYM_MAT_PART=='U'
		   enum { ID=0, SYMM_ID = Option::Upper};
		   enum :bool{ SymUpper = true };
		   static void sym_swap(index_type& r, index_type& c){ if (c < r) std::swap(r, c); }
#else
		   enum {ID=0, SYMM_ID = Option::Lower};
		   enum :bool{ SymUpper = false };
		   static void sym_swap(index_type& r, index_type& c){ if (r < c) std::swap(r, c); }
#endif
	   };


	   struct LT
	   {
		   enum { ID = Option::Lower, SYMM_ID = Option::Lower };
		   static const char* name(){ return  "packed LT"; }
		   static void sym_swap(index_type& r, index_type& c){ if (r < c) std::swap(r, c); }
	   };

	   struct UT
	   {
		   enum { ID = Option::Upper, SYMM_ID = Option::Upper };
		   static const char* name(){ return  "packed UT"; }
		   static void sym_swap(index_type& r, index_type& c){ if (c < r) std::swap(r, c); }
	   };
   };



   struct M_ORD {
	   struct ROW
	   {
		   enum  { ID = Option::RowMajor };
		   static const char* name(){ return  "ROW major"; }
	   };
	   struct COL
	   {
		   enum { ID = Option::ColMajor };
		   static const char* name(){ return  "COL major"; }
	   };
   };

   template< typename MP, typename MS>
   struct if_necessary;


   template<typename MS>
   struct if_necessary<M_PROP::GE, MS>
   {
	   template<typename I>
	   static void swap(I& r, I& c){}
   };

   template<typename MS>
   struct if_necessary<M_PROP::SYM, MS>
   {
	   template<typename I>
	   static void swap(I& r, I& c){ MS::sym_swap(r, c); }
   };








   struct BASE
   {
	   enum :index_type{ OFFSET = MML3_BASE_INDEX_OFFSET };
//
	   static const index_type		first(){ return OFFSET; }
	   static const index_type		last(index_type dim){ return dim + (OFFSET-1); }
	   static const index_type		end(index_type dim){ return dim + OFFSET; }
	   static bool					test_index(index_type i, index_type dim)	{ return ((i< OFFSET+dim) && i >= OFFSET) ? true : false; }
	   static index_type			op(index_type i){ return i- OFFSET; }




   };






}// end namespace
