#pragma once
#include<type_traits>
#include<algorithm>
#include<numeric>
#include<initializer_list>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<valarray>
#include<vector>
#include<stdexcept>
#include<cassert>
#include"MML3-memory_array_.h"
#include"MML3-config.h"
#include"MML3-iSet.h"
#include"MML3-SubMatrix.h"
#include"MML3-TypeSupport.h"

namespace MML3
{

	// the complete Matrix template
	template<typename T=double, typename MP=M_PROP::GE, typename MS=M_SHAPE::RE, typename MO=M_ORD::ROW>
	class Matrix;

	// the matrix selector
	template<typename T=double>
	struct Mat
	{

		typedef		Matrix<T, M_PROP::GE, M_SHAPE::RE, M_ORD::ROW>			type;
		typedef		Matrix<T, M_PROP::SYM, M_SHAPE::RE, M_ORD::ROW>			sym_type;

		// the whole variety of General matrices
		struct General
		{
			typedef M_PROP::GE  P;
			// default General
			typedef		Matrix<T, P, M_SHAPE::RE, M_ORD::ROW> type;

			struct Rectangular
			{
				typedef M_SHAPE::RE S;
				// default rectangular
				typedef		Matrix<T, P, S, M_ORD::ROW> type;

				struct RowMajor
				{
					typedef 	Matrix<T, P, S, M_ORD::ROW> type;
				};

				struct ColMajor
				{
					typedef 	Matrix<T, P, S, M_ORD::COL> type;
				};
			};

			struct LowerTria
			{
				typedef M_SHAPE::LT S;
				// default rectangular
				typedef		Matrix<T, P, S, M_ORD::ROW> type;

				struct RowMajor
				{
					typedef 	Matrix<T, P, S, M_ORD::ROW> type;
				};

				struct ColMajor
				{
					typedef 	Matrix<T, P, S, M_ORD::COL> type;
				};

			};

			struct UpperTria
			{
				typedef M_SHAPE::UT S;
				// default rectangular
				typedef		Matrix<T, P, S, M_ORD::ROW> type;

				struct RowMajor
				{
					typedef 	Matrix<T, P, S, M_ORD::ROW> type;
				};

				struct ColMajor
				{
					typedef 	Matrix<T, P, S, M_ORD::COL> type;
				};

			};
		};


		// the whole variety of symmetric matrices
		struct Symmetric
		{
			typedef M_PROP::SYM  P;
			// default General
			typedef		Matrix<T, P, M_SHAPE::RE, M_ORD::ROW> type;

			struct Rectangular
			{
				typedef M_SHAPE::RE S;
				// default rectangular
				typedef		Matrix<T, P, S, M_ORD::ROW> type;

				struct RowMajor
				{
					typedef 	Matrix<T, P, S, M_ORD::ROW> type;
				};

				struct ColMajor
				{
					typedef 	Matrix<T, P, S, M_ORD::COL> type;
				};
			};

			struct LowerTria
			{
				typedef M_SHAPE::LT S;
				// default rectangular
				typedef		Matrix<T, P, S, M_ORD::ROW> type;

				struct RowMajor
				{
					typedef 	Matrix<T, P, S, M_ORD::ROW> type;
				};

				struct ColMajor
				{
					typedef 	Matrix<T, P, S, M_ORD::COL> type;
				};

			};

			struct UpperTria
			{
				typedef M_SHAPE::UT S;
				// default rectangular
				typedef		Matrix<T, P, S, M_ORD::ROW> type;

				struct RowMajor
				{
					typedef 	Matrix<T, P, S, M_ORD::ROW> type;
				};

				struct ColMajor
				{
					typedef 	Matrix<T, P, S, M_ORD::COL> type;
				};

			};

		};

	};





	template<typename T, typename MP, typename MS, typename MO>
	class Matrix
	{
		typedef std::valarray<T>			array_;
		// the index swapper for symmetric matrices
		typedef MML3::if_necessary<MP, MS>	if_necessary;

	public:
		typedef T							value_t;
		typedef	index_type					index_t;
		typedef Matrix						matrix_t;
		typedef SubMatrix<T, Matrix>		sub_matrix_t;
		typedef const_SubMatrix<T, Matrix>	const_sub_mat_t;
		typedef MP							property_t;
		typedef MS							shape_t;
		typedef MO							order_t;
		typedef iSet						iset_t;


		enum IS :bool{	RowMajor= std::is_same<MO, M_ORD::ROW>::value,
						ColMajor= !RowMajor,
						RE		= std::is_same<M_SHAPE::RE, MS>::value,
						LT		= std::is_same<M_SHAPE::LT, MS>::value,
						UT		= std::is_same<M_SHAPE::UT, MS>::value,
						SYM		= std::is_same<M_PROP::SYM, MP>::value,
						GE		= std::is_same<M_PROP::GE, MP>::value,
						Packed  = LT||UT
		};
		enum :index_t {ONE_=BASE::OFFSET};


		//------------------------------------//
		//           CTORS                    //
		//------------------------------------//
		Matrix()					=default;
		Matrix(const Matrix& rval)	=default;
		Matrix(Matrix&& rval) :data_(std::move(rval.data_)), nr_(rval.nr_), nc_(rval.nc_){ rval.nr_ = 0; rval.nc_ = 0; }
		explicit Matrix(index_t nr, index_t nc = 1)							{resize(nr, nc);}
		Matrix(const std::initializer_list<T>& s)							{assign(s); }
		Matrix(const std::initializer_list<std::initializer_list<T>>& s)	{assign(s);}
		Matrix(const sub_matrix_t& o)										{assign(o);}
		Matrix(const const_sub_mat_t& o)									{assign(o); }
		Matrix(const T* src, index_t src_size, index_t nr, index_t nc)		{assign(src,src_size,nr,nc); }
		~Matrix() =default;

		Matrix&		resize(index_t nr, index_t nc = 1);
		void		swap(Matrix& o)											{data_.swap(o.data_); std::swap(nr_, o.nr_); std::swap(nc_, o.nc_); }
		void		free()													{data_.free(); nr_ = nc_ = 0; }

		//------------------------------------//
		//      INFO                          //
		//------------------------------------//
		index_t		nrows()const											{return nr_; }
		index_t		ncols()const											{return nc_; }
		size_t		size()const												{return data_.size(); }
		index_t		base_index()const										{return ONE_; }

		static const std::string	name();

		index_t		leading_dim()const										{return IS::RowMajor?nc_:nr_;}
		static Option::Order	CblasOrder()								{return Option::Order(MO::ID); }
		// tells where the components of a symmetric matrix are stored (Lower or Upper triangle)
		static Option::UpLo	CblasSymUpLo()									{return Option::UpLo(MS::SYMM_ID); }

		//  REctangualr GEneral matrices can have any dimension. Other combinations must be square
		bool		test_size(index_type nr, index_type nc)const			{ return  (IS::RE && IS::GE) ? true : (nr == nc); }

		bool		test_indexes(index_t r, index_t c)const;

		index_t		col_beg_idx(index_t r)const								{return IS::UT?r:ONE_; }
		index_t		col_end_idx(index_t r)const								{return (IS::LT?r:nc_)+ONE_; }
		index_t		row_beg_idx(index_t c)const								{return IS::LT?c:ONE_; }
		index_t		row_end_idx(index_t c)const								{return (IS::UT?c:nr_)+ONE_; }
		// the 0b position into the array data_ of the components [r][c] (0b)
		index_t     pos_0b(index_t r, index_t c)const						{ return pos_0b_(r, c); }

		// these are safe ranges to iterate over matrix components
		iRange		row_range()const										{return iRange(ONE_, nr_+ONE_); }
		iRange		col_range()const{ return iRange(ONE_, nc_+ONE_); }
		iRange		row_range(index_t col)const								{return iRange(row_beg_idx(col),row_end_idx(col)); }
		iRange		col_range(index_t row)const								{return iRange(col_beg_idx(row), col_end_idx(row)); }

		iSet		all_ridx()const											{ return iSet(ONE_, nr_ - (1 - ONE_)); }
		iSet		all_cidx()const											{ return iSet(ONE_, nc_ - (1 - ONE_)); }

		//------------------------------------//
		//      ASSIGNMENT                    //
		//------------------------------------//
		Matrix&		operator=(const Matrix&  rhs) = default;
		Matrix&		operator=(Matrix&&  o)									{ data_ = (std::move(o.data_)); nr_ = o.nr_; nc_ = o.nc_; return *this; }
		Matrix&		operator=(value_t  val)									{fill(val); return *this; }
		Matrix&		operator=(const std::initializer_list<T>& s)			{return assign(s); }
		Matrix&		operator=(const std::initializer_list<std::initializer_list<T>>& s){ return assign(s); }
		Matrix&		operator=(const sub_matrix_t&  o)						{return assign(o); }
		Matrix&		operator=(const const_sub_mat_t&  o)					{return assign(o); }

		Matrix&		assign(const std::initializer_list<std::initializer_list<T>> & s);
		Matrix&		assign(const std::initializer_list<T> & s);
		Matrix&		assign(const sub_matrix_t& o);
		Matrix&		assign(const const_sub_mat_t& o);
		template<typename G>
		Matrix&		assign(const G*src, index_t src_size, index_t nr, index_t nc);

		void		fill(T val)												{ data_ = val; }

		//------------------------------------//
		//      ACCESSORS                     //
		//------------------------------------//
		//  pointers to the first Matrix component
		T*			begin()													{ return data_.size() ? &(data_[0]) : nullptr; }
		const T*	begin()const											{ return data_.size() ? &(data_[0]) : nullptr; }
		//  pointers to the last+1 Matrix component
		T*			end()													{ return data_.size() ? (&(data_[0]) + data_.size()) : nullptr; }
		const T*	end()const												{ return data_.size() ? (&(data_[0]) + data_.size()) : nullptr; }
		// BASE::OFFSET-based indices (by default 1)
		T&			operator()(index_t r, index_t c);
		const T&	operator()(index_t r, index_t c)const;
		T& 			operator()(index_t r);
		const T& 	operator()(index_t r)const;

		// access with 0-based indexes
		T&			at_(index_t r, index_t c)								{ return data_[pos_0b_(r, c)]; }
		const T&	at_(index_t r, index_t c)const							{ return data_[pos_0b_(r, c)]; }
		// C access operator 0-b defined only for Row major Rectangular matrices
		T*			operator[](index_t r)									{ static_assert(IS::RowMajor && IS::RE, " C access operator [] is defined only for Row major rectangular matrices"); return begin() + r*nc_; }
		const T*	operator[](index_t r)const								{ static_assert(IS::RowMajor && IS::RE, " C access operator [] is defined only for Row major rectangular matrices"); return begin() + r*nc_; }


		//------------------------------------//
		//      SUB MATRIX ACCESS             //
		//------------------------------------//

		// Matrix 2D access
		sub_matrix_t	operator()(const iSet& 	r, const iSet& 	c={1})		{ static_assert(IS::RE, ""); return sub_matrix_t(*this, r, c); }
		sub_matrix_t	operator()(iSet&& 		r, const iSet&  c={1})		{ static_assert(IS::RE, ""); return sub_matrix_t(*this, std::move(r), c); }
		sub_matrix_t	operator()(const iSet& 	r, iSet&& 		c)			{ static_assert(IS::RE, ""); return sub_matrix_t(*this, r, std::move(c)); }
		sub_matrix_t	operator()(iSet&& 		r, iSet&& 		c)			{ static_assert(IS::RE, ""); return sub_matrix_t(*this, std::move(r), std::move(c)); }

		const_sub_mat_t	operator()(const iSet& r, const iSet& c={1})const	{ static_assert(IS::RE, ""); return const_sub_mat_t(*this, r, c); }
		const_sub_mat_t	operator()(iSet&& r, const iSet& c={1})const		{ static_assert(IS::RE, ""); return const_sub_mat_t(*this, std::move(r), c); }
		const_sub_mat_t	operator()(const iSet& r, iSet&& c)const			{ static_assert(IS::RE, ""); return const_sub_mat_t(*this, r, std::move(c)); }
		const_sub_mat_t	operator()(iSet&& r, iSet&& c)const					{ static_assert(IS::RE, ""); return const_sub_mat_t(*this, std::move(r), std::move(c)); }


		sub_matrix_t	row(const iSet&	r)              					{ static_assert(IS::RE, ""); return sub_matrix_t(*this, r, all_cidx()); }
		sub_matrix_t	row(iSet&&	r)              						{ static_assert(IS::RE, ""); return sub_matrix_t(*this, std::move(r), all_cidx()); }
		sub_matrix_t	column(const iSet& c)           					{ static_assert(IS::RE, ""); return sub_matrix_t(*this, all_ridx(), c); }
		sub_matrix_t	column(iSet&&	c)              					{ static_assert(IS::RE, ""); return sub_matrix_t(*this, all_ridx(), std::move(c)); }

		const_sub_mat_t	row(const iSet&	r)const     						{ static_assert(IS::RE, ""); return const_sub_mat_t(*this, r, all_cidx()); }
		const_sub_mat_t	row(iSet&& r)const          						{ static_assert(IS::RE, ""); return const_sub_mat_t(*this, std::move(r), all_cidx()); }
		const_sub_mat_t	column(const iSet& c)const							{ static_assert(IS::RE, ""); return const_sub_mat_t(*this, all_ridx(), c); }
		const_sub_mat_t	column(iSet&&	c)const								{ static_assert(IS::RE, ""); return const_sub_mat_t(*this, all_ridx(), std::move(c)); }


		//------------------------------------//
		//           IO                       //
		//------------------------------------//
		int				fwrite(const char* fname);
		int				fread(const char* fname);
		bool			print(std::ostream& os)const;

		//------------------------------------//
		//  MATH  OPS                   //
		//------------------------------------//
		Matrix&			operator+=(const Matrix& o)						{ data_ += o.data_;	return *this; }
		Matrix&			operator-=(const Matrix& o)						{ data_ -= o.data_;	return *this; }
		Matrix&			operator+()										{ return *this; }
		Matrix			operator-()const								{ Matrix tmp(*this); tmp.data_ = -data_;	return tmp; }
		Matrix&			operator*=(const T& val)						{ data_ *= val;	return *this; }
		Matrix&			operator/=(const T& val)						{ data_ /= val;	return *this; }
		//returns a copy of the Matrix transformed by func
		Matrix			apply(T func(T)) const							{ Matrix tmp(*this); return tmp.transform(func); }
		Matrix&			transform(T func(T))							{ for (auto& v : data_)v = func(v);return *this;}

		//------------------------------------//
		//  GENERATORS                        //
		//------------------------------------//
		static   Matrix get_random_matrix(index_type nr, index_type nc, T mean, T delta);
		static   Matrix get_eye(index_type n);

	private:
		//functions that need specialization

		// the 0b position of the component [r,c]0b into data_
		index_t		pos_0b_(index_t r, index_t c)const;
		// number of components for every  Shape
		index_type	num_elements_(index_type nr, index_type nc)const;

    private:
		// the class data
		array_ data_;
		index_t nr_ = 0, nc_ = 0;


	};


	//------------------------------------//
	//   Matrix NON MEMBER OPERATORS      //
	//------------------------------------//


	template<typename T, typename MP, typename MS, typename MO>
	std::ostream& operator << (std::ostream& os, const Matrix<T, MP, MS, MO>& obj)
	{
		obj.print(os);
		return os;
	}



	// tutte le simmetriche non cambiano con la trasposizione
	template<typename T, typename MP, typename MS, typename MO, typename MP1, typename MS1, typename MO1>
	Matrix<T, MP1, MS1, MO1>     transpose(const Matrix<T, MP, MS, MO>&);


	template<typename T, typename MP, typename MS, typename MO>
	inline Matrix<T, MP, MS, MO> operator +(const Matrix<T, MP, MS, MO>& x, const Matrix<T, MP, MS, MO>& y) { return (Matrix<T, MP, MS, MO>(x) += y); }

	template<typename T, typename MP, typename MS, typename MO>
	inline Matrix<T, MP, MS, MO> operator -(Matrix<T, MP, MS, MO> x, const Matrix<T, MP, MS, MO>& y) { return (x -= y); }

	template<typename T, typename MP, typename MS, typename MO>
	inline Matrix<T, MP, MS, MO> operator *(Matrix<T, MP, MS, MO> x, const T& val) { return (x *= val); }

	template<typename T, typename MP, typename MS, typename MO>
	inline Matrix<T, MP, MS, MO> operator *(const T& val, Matrix<T, MP, MS, MO> x) { return (x *= val); }

	template<typename T, typename MP, typename MS, typename MO>
	inline Matrix<T, MP, MS, MO> operator /(Matrix<T, MP, MS, MO> x, const T& val){ return (x /= val); }







	////////////////////////////////////////////////////////////////////////////////
	// IMPLEMENTATION




	template<typename T, typename MP, typename MS, typename MO>
	inline auto Matrix<T, MP, MS, MO>::pos_0b_(index_t r, index_t c)const->index_t
	{
		if (IS::RowMajor)
		{
			if (IS::RE)
				return nc_*r + c;
			else if (IS::LT)
				return (r*(r + 1)) / 2 + c;
			else if (IS::UT)
				return r*nc_ - (r*(r - 1)) / 2 + (c - r);

		}
		else if (IS::ColMajor)
		{
			if (IS::RE)
				return r + c*nr_;
			else if (IS::LT)
				return r + ((2 * nr_ - c - 1)*c) / 2;
			else if (IS::UT)
				return (c*(c + 1)) / 2 + r;
		}
		else
			assert(false);
		return 0;
	}






	// BASE::OFFSET-based indices (by default 1)
	template<typename T, typename MP, typename MS, typename MO>
	 inline T& Matrix<T, MP, MS, MO>::operator()(index_t r, index_t c)
	{
#ifdef MML3_TEST_INDEX_ON_ACCESS
		 test_indexes(r, c);
#endif
		if_necessary::swap(r, c);
		return at_(r - ONE_, c - ONE_);
	}
	template<typename T, typename MP, typename MS, typename MO>
	inline const T&	Matrix<T, MP, MS, MO>::operator()(index_t r, index_t c)const
	{
#ifdef MML3_TEST_INDEX_ON_ACCESS
		test_indexes(r, c);
#endif

		if_necessary::swap(r, c);
		return at_(r - ONE_, c - ONE_);
	}
	template<typename T, typename MP, typename MS, typename MO>
	inline T& 			Matrix<T, MP, MS, MO>::operator()(index_t r)
	{
#ifdef MML3_TEST_INDEX_ON_ACCESS
		test_indexes(r, 1);
#endif

		return at_(r - ONE_, 0);
	}
	template<typename T, typename MP, typename MS, typename MO>
	inline const T& 	Matrix<T, MP, MS, MO>::operator()(index_t r)const
	{
#ifdef MML3_TEST_INDEX_ON_ACCESS
		test_indexes(r, 1);
#endif
		return at_(r - ONE_, 0);
	}



	template<typename T, typename MP, typename MS, typename MO>
	inline bool 	Matrix<T, MP, MS, MO>::test_indexes(index_t r, index_t c)const
	{
		if (r < ONE_ || !(r <(nr_ + ONE_)))
			return false;
		if (IS::SYM || (IS::GE && IS::RE))
			return (c < (nc_ + ONE_) && c >= ONE_);
		else if (IS::LT)
			return  (c <= r && c >= ONE_);
		else if (IS::UT)
			return (c >= r && c <(nc_+ONE_));
	}

template<typename T, typename MP, typename MS, typename MO>
inline auto Matrix<T, MP, MS, MO>::num_elements_(index_type nr, index_type nc)const->index_t
{

    if(IS::RE)
        return nr*nc;
	else    if (IS::LT || IS::UT)
		return nr == nc ? (nr*(nr + 1)) / 2 : 0;
	else
		assert(false);
    return 0;
}


	template<typename T, typename MP, typename MS, typename MO>
	template<typename G>
	inline Matrix<T, MP, MS, MO>& Matrix<T, MP, MS, MO>::assign(const G*src, index_t src_size, index_t nr, index_t nc)
	{
		if (num_elements_(nr, nc) != src_size)
			throw std::runtime_error("Matrix::assign(*src, src_size, nr, nc)");
		resize(nr, nc);
		std::copy(src, src + src_size, begin());
		return *this;
	}

	// tutte le simmetriche non cambiano con la trasposizione
	template<typename T, typename MS,typename MO>
	inline Matrix<T, M_PROP::SYM, MS, MO>     transpose(const Matrix<T, M_PROP::SYM, MS, MO>& M)
	{
		return M;
	}




	// le rettangolari (non simmetriche ) si trattano in modo standard
	template<typename T, typename MO>
	Matrix<T, M_PROP::GE, M_SHAPE::RE, MO>     transpose(const Matrix<T, M_PROP::GE, M_SHAPE::RE, MO>& M)
	{
		typedef Matrix<T, M_PROP::GE, M_SHAPE::RE, MO> Mat_t;
		Mat_t R(M.ncols(), M.nrows());
		size_t ldin = M.leading_dim(), ldout = R.leading_dim();
		size_t  x = M.nrows(), y = M.ncols();
		if (Mat_t::IS::ColMajor)
			std::swap(x, y);
		size_t end_i = std::min(y, ldin),
			   end_j = std::min(x, ldout);
		const T* in = M.begin();
		T* out = R.begin();
		for (size_t i = 0; i < end_i; i++)
				for (size_t j = 0; j < end_j; j++)
					out[i*ldout + j] = in[j*ldin + i];
		return R;
	}

	// la trasposta di una triangolare (non simmetriche ) è sul triangolo opposto
	template<typename T>
	Matrix<T, M_PROP::GE, M_SHAPE::UT, M_ORD::COL>     transpose(const Matrix<T, M_PROP::GE, M_SHAPE::LT, M_ORD::ROW>& M)
	{
		Matrix<T, M_PROP::GE, M_SHAPE::UT, M_ORD::COL> R(M.ncols(), M.nrows());
		std::copy(M.begin(), M.end(), R.begin());
		return R;
	}
	template<typename T>
	Matrix<T, M_PROP::GE, M_SHAPE::UT, M_ORD::ROW>     transpose(const Matrix<T, M_PROP::GE, M_SHAPE::LT, M_ORD::COL>& M)
	{
		Matrix<T, M_PROP::GE, M_SHAPE::UT, M_ORD::ROW> R(M.ncols(), M.nrows());
		std::copy(M.begin(), M.end(), R.begin());
		return R;
	}

	template<typename T>
	Matrix<T, M_PROP::GE, M_SHAPE::LT, M_ORD::ROW>     transpose(const Matrix<T, M_PROP::GE, M_SHAPE::UT, M_ORD::COL>& M)
	{
		Matrix<T, M_PROP::GE, M_SHAPE::LT, M_ORD::ROW> R(M.ncols(), M.nrows());
		std::copy(M.begin(), M.end(), R.begin());
		return R;
	}

	template<typename T>
	Matrix<T, M_PROP::GE, M_SHAPE::LT, M_ORD::COL>     transpose(const Matrix<T, M_PROP::GE, M_SHAPE::UT, M_ORD::ROW>& M)
	{
		Matrix<T, M_PROP::GE, M_SHAPE::LT, M_ORD::COL> R(M.ncols(), M.nrows());
		std::copy(M.begin(), M.end(), R.begin());
		return R;
	}



	template<typename T, typename MP, typename MS, typename MO>
	auto Matrix<T, MP, MS, MO>::resize(index_t nr, index_t nc )->Matrix&
	{
		if (!test_size(nr, nc))
			throw std::runtime_error("MML3::Matrix::resize(nr,nc): dimensions incompatible with the matrix type");
		index_t new_sz = num_elements_(nr, nc);
		if (new_sz != data_.size())
			data_.resize(new_sz);
		nr_=nr;
		nc_=nc;
		return *this;
	}


	//print a matrix in text format (human readable),
	template<typename T, typename MP, typename MS, typename MO>
	bool  Matrix<T, MP, MS, MO>::print(std::ostream& os)const
	{
		std::streamsize ssize = os.width();
		std::streamsize sprecision = os.precision();
		os << nrows() << " * " << ncols() << std::endl;

		for (index_t r = ONE_; r != (nr_ + 1); ++r)
		{
			index_t first_cidx = col_beg_idx(r);
			index_t end_cidx = col_end_idx(r);
			for (index_t c = ONE_; c != (nc_ + ONE_); ++c)
			{
			if ( c >= first_cidx && c < end_cidx)
				os << std::setw(ssize) << std::setprecision(sprecision) << this->operator()(r, c) << ' ';
			else
				os << std::setw(ssize) << std::setprecision(sprecision) << '.' << ' ';
			}
			os << std::endl;
		}
		return os.good();
	}




	template<typename T, typename MP, typename MS, typename MO>
	auto Matrix<T, MP, MS, MO>::assign(const sub_matrix_t& o)->Matrix&
	{

		index_t nr = o.nrows(), nc = o.ncols();
		resize(nr, nc);
		for (index_t r = base_index(); r != nr + 1; ++r)
		{
			auto end = col_end_idx(r);
			for (auto c = col_beg_idx(r); c != end; ++c)
				operator()(r, c) = o(r, c);
		}

		return *this;
	}

	template<typename T, typename MP, typename MS, typename MO>
	auto Matrix<T, MP, MS, MO>::assign(const const_sub_mat_t& o)->Matrix&
	{

		index_t nr = o.nrows(), nc = o.ncols();
		resize(nr, nc);
		for (index_t r = base_index(); r != nr + 1; ++r)
		{
			auto end = col_end_idx(r);
			for (auto c = col_beg_idx(r); c != end; ++c)
				operator()(r, c) = o(r, c);
		}

		return *this;
	}




	/// fwrite/fread : Scrivono e Leggono la matrice su/da un file binario
	///@param fname nome del file
	///@return:  0 in caso di successo
	///@return: -1 nel caso il file non possa essere aperto
	///@return: -2 se si sono verificati errori in lettura
	///@return: -3 (solo fread) se i tipi  del'indice o dei valori su file non coincidono con quelli della Matrix
	template<typename T, typename MP, typename MS, typename MO>
	int  Matrix<T, MP, MS, MO>::fwrite(const char* fname)
	{
		std::ofstream f(fname, std::ios_base::binary);
		if (!f.is_open())
			return -1;
		// l'intestazione è costituita da 7 uint64 che indicano:
		// numero di  righe,
		// numero di colonne
		// ID del tipo degli indici
		// ID del tipo delle componenti
		// ID di proprieta'
		// ID di forma
		// ID di orientamento

		std::int64_t type_[7] = { nr_, nc_, type<index_t>::id(), type<T>::id(), MP::ID, MS::ID,MO::ID };
		f.write(reinterpret_cast<const char*>(type_), sizeof(type_));
		f.write(reinterpret_cast<const char*>(begin()), size()*sizeof(T));
		return f.good() ? 0 : -2;
	}


	template<typename T, typename MP, typename MS, typename MO>
	int  Matrix<T, MP, MS, MO>::fread(const char* fname)
	{
		std::ifstream f(fname, std::fstream::binary);
		if (!f.is_open())
			return -1;
		// leggo l'intestazione
		std::int64_t type_[7];
		f.read(reinterpret_cast<char*>(type_), sizeof(type_));
		if (!f.good())
			return -2;
		// controllo di consistenza dei tipi

		if (type_[2] != std::uint64_t(type<index_t>::id()) ||
			type_[3] != std::uint64_t(type<T>::id()) ||
			type_[4] != std::uint64_t(MP::ID) ||
			type_[5] != std::uint64_t(MS::ID) ||
			type_[6] != std::uint64_t(MO::ID)
			)
			return -3;
		resize(index_t(type_[0]), index_t(type_[1]));
		f.read(reinterpret_cast<char*>(begin()), size()*sizeof(T));
		return f.good() ? 0 : -2;
	}





	template<typename T, typename MP, typename MS, typename MO>
	auto  Matrix<T, MP, MS, MO>::assign(const std::initializer_list<std::initializer_list<T>> & s)->Matrix&

	{
		size_t nr = s.size();
		size_t nc = 0;
		// the number of columns in the initializer list matrix is determined as the maximum of the row lengts
		for (size_t r = 0; r != nr; ++r)
		{
			const std::initializer_list<T>* r_ptr = s.begin() + r;
			nc = std::max(nc, r_ptr->size());
		}

		// resize this
		resize(nr, nc);
		const std::initializer_list<T>* r_ptr = s.begin() ;

		for (size_t r = BASE::OFFSET; r != (nr + BASE::OFFSET); ++r)
		{
			// pointer to row r of the initializer list
			index_t end_idx = col_end_idx(r);
			const T*	c_init_ptr     = (*r_ptr).begin();
			const T*    c_init_end_ptr = (*r_ptr).end();

			for (index_t c = col_beg_idx(r); (c != end_idx) && (c_init_ptr != c_init_end_ptr); ++c)
				operator()(r,c) = *c_init_ptr++;

			++r_ptr;
		}
		return *this;
	}

	template<typename T, typename MP, typename MS, typename MO>
	auto  Matrix<T, MP, MS, MO>::assign(const std::initializer_list<T> & s)->Matrix&
	{
		size_t sz = s.size();
		resize(sz, 1);
		std::copy(s.begin(), s.end(), begin());
		return *this;
	}





	// genera una Matrix nr x nc con componenti random distribuite nell'intervallo (mean - delta/2, mean+delta/2)
	template<typename T, typename MP, typename MS, typename MO>
	auto  Matrix<T, MP, MS, MO>::get_random_matrix(index_type nr, index_type nc, T mean, T delta)->Matrix
	{

		Matrix A(nr, nc);
		// inizializzo il generatore
		srand((unsigned)time(NULL));
		T*        p = A.begin();
		const T*  end = A.end();

		while (p != end)
			*p++ = mean + (T(rand()) / (RAND_MAX ) - 0.5) * delta;
		return A;
	}

	template<typename T, typename MP, typename MS, typename MO>
	auto  Matrix<T, MP, MS, MO>::get_eye(index_type n)->Matrix
	{
		Matrix tmp(n, n);
		tmp.fill(0);
		for (index_t i = ONE_; i != n + 1; ++i)
			tmp(n, n) = T(1);
		return tmp;
	}

	template<typename T, typename MP, typename MS, typename MO>
	const std::string	  Matrix<T, MP, MS, MO>::name()
	{

		std::string tmp = std::string("Matrix<") + typeid(T).name() + ","
			+ MP::name() + ","
			+ MS::name() + ","
			+ MO::name() + ">";
		return tmp;

	}






}// end namespace
