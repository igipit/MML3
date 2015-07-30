#pragma once
#include<type_traits>
#include<algorithm>
#include<initializer_list>
#include"MML3-config.h"
#include"MML3-types.h"

#include "MML3-array_.h"



namespace MML3
{


	template<typename T>
	class mat_T
	{
	public:
		typedef T			value_t;
		typedef	index_type  index_t;


		index_t							nrows()const{ return nr_; }
		index_t							ncols()const{ return nc_; }
		index_t							size()const { return data_.size; }
		index_t							leading_dimension()const                { return nc_; }
		bool							is_row_major()const						{ return true; }

		index_t							start_index()const{ return BASE::OFFSET;}
	//	virtual index_t					first_col_index(index_t r)const = 0;
		//virtual index_t					end_col_index(index_t r)const = 0;


		void							fill(T val)	{ data_ = val; }
		void							free()		{ data_.destroy(); nr_ = nc_ = 0; }
		



		//  pointers to the first Matrix component
		T*					begin()                         	{ return data_.size()?&(data_[0]):nullptr; }
		const T*			begin()const                    	{ return data_.size()?&(data_[0]):nullptr; }
		//  pointers to the last Matrix component +1
		T*					end()                           	{ return data_.size() ? (&(data_[0]) + data_.size()) : nullptr; }
		const T*			end()const                      	{ return data_.size() ? (&(data_[0]) + data_.size()) : nullptr; }


	protected:
		// CTORs
		mat_T() = default;
		mat_T(const mat_T& rhs) = default;
		mat_T(mat_T&& rhs) :data_(std::move(rhs.data_)), nr_(o.nr_), nc_(o.nc_){}
		mat_T(index_t nr, index_t nc, index_t sz) :nr_(nr), nc_(nc), data_(sz){}
			
		mat_T& operator =(const mat_T&) = default;
		mat_T& operator =(mat_T&& o){data_ = std::move(o.data_);nr_ = o.nr_;nc_ = o.nc_;return *this;}

		void resize(index_t nr, index_t nc, index_t sz) 
		{
			if (sz != data_.size())
				data_.resize(sz);
			nr_ = nr;
			nc_ = nc;
		}
		void swap(mat_T& o){ data_.swap(o.data_); std::swap(nr_, o.nr_); std::swap(nc_, o.nc_); }
		
		index_t nr_ = 0, nc_ = 0;
		std::valarray<T> data_;
	};


	template<typename T, typename ORD, typename SHP>
	class mat_T_O_S: public mat_T<T>
	{
	
	};


	// ROW MAJOR ORDER SPECIALIZATIONS
	// RE
	template<typename T>
	class mat_T_O_S<T,MAT_ORD::ROW, MAT_SHAPE::RE> : public mat_T<T>
	{
		
	protected:
		using mat_T<T>::mat_T<T>;
		using mat_T<T>::operator=;

		T*			get_p_0b(index_t r, index_t c){ return begin() + nc_ * r + c; };
		const T*	get_p_0b(index_t r, index_t c)const{ return begin() + nc_ * r + c; };
		index_t		get_size_(index_t nr, index_t nc)const{ return nr_*nc_; }
	};

	// LT
	template<typename T>
	class mat_T_O_S<T, MAT_ORD::ROW, MAT_SHAPE::LT> : public mat_T<T>
	{
	protected:
		using mat_T<T>::mat_T<T>;
		using mat_T<T>::operator=;
		T*			get_p_0b(index_t r, index_t c){ return begin() + (r*(r + 1))/2; };
		const T*	get_p_0b(index_t r, index_t c)const { return begin() + (r*(r + 1)) / 2; };
		index_t		get_size_(index_t nr, index_t nc)const{ return (nr_*(nr_-1))/2; }
	};

	// UT
	template<typename T>
	class mat_T_O_S<T, MAT_ORD::ROW, MAT_SHAPE::UT> : public mat_T<T>
	{
	protected:
		using mat_T<T>::mat_T<T>;
		using mat_T<T>::operator=;
		T*			get_p_0b(index_t r, index_t c){ return begin()+ r*nc_ - (r*(r - 1)) / 2; };
		const T*	get_p_0b(index_t r, index_t c)const{ return begin() + r*nc_ - (r*(r - 1)) / 2; };
		index_t		get_size_(index_t nr, index_t nc)const{ return (nr_*(nr_ - 1)) / 2; }
	};

	// COL MAJOR ORDER SPECIALIZATIONS
	//RE
	template<typename T>
	class mat_T_O_S<T, MAT_ORD::COL, MAT_SHAPE::RE> : public mat_T<T>
	{
	protected:
		using mat_T<T>::mat_T<T>;
		using mat_T<T>::operator=;
		T*			get_p_0b(index_t r, index_t c){ return begin() + nr_ * c + r; };
		const T*	get_p_0b(index_t r, index_t c)const { return begin() + nr_ * c + r; };
		index_t		get_size_(index_t nr, index_t nc)const{ return nr_*nc_; }
	};
	//LT
	template<typename T>
	class mat_T_O_S<T, MAT_ORD::COL, MAT_SHAPE::LT> : public mat_T<T>
	{
	protected:
		using mat_T<T>::mat_T<T>;
		using mat_T<T>::operator=;
		T*			get_p_0b(index_t r, index_t c){ return begin() + c*nr_ - (c*(c - 1))/2; };
		const T*	get_p_0b(index_t r, index_t c)const { return begin() + c*nr_ - (c*(c - 1)) / 2; };
		index_t		get_size_(index_t nr, index_t nc)const{ return (nr_*(nr_ - 1)) / 2; }
	};

	//UT
	template<typename T>
	class mat_T_O_S<T, MAT_ORD::COL, MAT_SHAPE::UT> : public mat_T<T>
	{
	protected:
		using mat_T<T>::mat_T<T>;
		using mat_T<T>::operator=;
		T*			get_p_0b(index_t r, index_t c){ return + begin() + (c*(c + 1)) / 2; };
		const T*	get_p_0b(index_t r, index_t c)const{ return +begin() + (c*(c + 1)) / 2; };
		index_t		get_size_(index_t nr, index_t nc)const{ return (nr_*(nr_ - 1)) / 2; }
	};

	

	template<typename T, typename MSHP>
	class Mat:public mat_T<T>
	{
		
	};

	template<typename T>
	class Mat<T,MAT_SHAPE::RE>:public mat_T<T>
	{
	protected:
		typedef mat_T<T> mat_t;
		T*			p_row_0b(index_t r){ return &(data_[r*nc_]); }
		const T*	p_row_0b(index_t r)const { return &(data_[r*nc_]); }
	
		Mat() = default;
		Mat(const Mat& ) = default;
		Mat(Mat&& o) :mat_t(std::move(o)){}
		Mat(index_t nr, index_t nc) :mat_t(nr, nc, nr*nc){}
		Mat& operator=(const Mat& rhs) = default;
		Mat& operator=(Mat&& o){ mat_t::operator =(std::move(o)); return *this; }

		index_t		first_col_index(index_t r)const { return BASE::OFFSET; }
		index_t		end_col_index(index_t r)const	{ return nc_ + BASE::OFFSET; }


		T&			at_0b(index_t r)			{ return data_[r*nc_]; }
		const T&	at_0b(index_t r)const		{ return data_[r*nc_]; }
		T&			at(index_t r)				{ return data_[(r - BASE::OFFSET)*nc_]; }
		const T&    at(index_t r)const			{ return data_[(r - BASE::OFFSET)*nc_]; }
		T&			at(index_t r, index_t c)	{ return  data_[(r - BASE::OFFSET)*nc_ + c - BASE::OFFSET]; }
		const T&	at(index_t r, index_t c)const{ return data_[(r - BASE::OFFSET)*nc_ + c - BASE::OFFSET]; }

		
		

	public:
		void resize(index_t nr, index_t nc){ mat_t::resize(nr, nc, nr*nc);}
		void swap(Mat& o){ mat_t::swap(o); }
	};




	

	

	template<typename T>
	class Mat<T, MAT_SHAPE::LT>:public mat_T<T>
	{
	protected:
		typedef mat_T<T> mat_t;

		T*			p_row_0b(index_t r){ return &(data_[(r*(r + 1)) / 2]); }
		const T*	p_row_0b(index_t r)const { return &(data_[(r*(r + 1)) / 2]); }
		Mat() = default;
		Mat(const Mat& rhs) = default;
		
		Mat(Mat&& o) :mat_t(std::move(o)){}
		Mat(index_t nr, index_t nc) :mat_t(nr, nc, (nr*(nr+1))/2){ if (nr != nc) throw std::runtime_error("Mat<T, MAT_SHAPE::LT> must be square"); }
		Mat& operator=(const Mat& rhs) = default;
		Mat& operator=(Mat&& o){ mat_t::operator =(std::move(o)); return *this; }


		index_t		first_col_index(index_t r)const { return BASE::OFFSET; }
		index_t		end_col_index(index_t r)const	{ return r; }

		T&			at_0b(index_t r)		{ return *p_row_0b(r); }
		const T&	at_0b(index_t r)const	{ return *p_row_0b(r);; }
		T&			at(index_t r)			{ return *p_row_0b(r-BASE::OFFSET); }
		const T&    at(index_t r)const		{ return *p_row_0b(r - BASE::OFFSET); }
		T&			at(index_t r, index_t c){return  * (p_row_0b(r - BASE::OFFSET) + c-BASE::OFFSET); }
		const T&	at(index_t r, index_t c)const{return  *(p_row_0b(r - BASE::OFFSET) + c - BASE::OFFSET);}

	public:
		void resize(index_t nr, index_t nc){ mat_t::resize(nr, nc, (nr*(nr + 1)) / 2); }
		void swap(Mat& o){ mat_t::swap(o); }
	};



	template<typename T>
	class Mat<T, MAT_SHAPE::UT>:public mat_T<T>
	{
	protected:
		typedef mat_T<T> mat_t;
		T*			p_row_0b(index_t r)		 { return &(data_[r*nc_ - (r*(r - 1))/2]); }
		const T*	p_row_0b(index_t r)const { return &(data_[r*nc_ - (r*(r - 1))/2]); }
	
		Mat() = default;
		Mat(const Mat& rhs) = default;
		Mat(Mat&& o) :mat_t(std::move(o)){}
		Mat(index_t nr, index_t nc) :mat_t(nr, nc, (nr*(nr + 1)) / 2){ if (nr != nc) throw std::runtime_error("Mat<T, MAT_SHAPE::LT> must be square"); }
		Mat& operator=(const Mat& rhs) = default;
		Mat& operator=(Mat&& o){ mat_t::operator =(std::move(o)); return *this; }


		index_t		first_col_index(index_t r)const { return r; }
		index_t		end_col_index(index_t r)const	{ return nc_ + BASE::OFFSET; }


		T&			at_0b(index_t r)		{ return *p_row_0b(r); }
		const T&	at_0b(index_t r)const	{ return *p_row_0b(r); }
		T&			at(index_t r)			{ return *p_row_0b(r - BASE::OFFSET); }
		const T&    at(index_t r)const		{ return *p_row_0b(r - BASE::OFFSET); }
		T&			at(index_t r, index_t c){return  *(p_row_0b(r - BASE::OFFSET) + c-r);}
		const T&	at(index_t r, index_t c)const{return  *(p_row_0b(r - BASE::OFFSET) + c-r);}

	public:
		void resize(index_t nr, index_t nc){ mat_t::resize(nr, nc, (nr*(nr + 1)) / 2); }
		void swap(Mat& o){ mat_t::swap(o); }
	};


	
	template< typename MP,typename MS>
	struct MPS_;
	// these are 3 of the possible combinations of properties and shapes
	template<>
	struct MPS_< MAT_PROP::GE, MAT_SHAPE::RE>
	{
		static inline void swap(index_type& r, index_type& c){}
	};

	template<>
	struct MPS_< MAT_PROP::SYM, MAT_SHAPE::RE>
	{
		static inline void swap(index_type& r, index_type& c){ if (c > r) std::swap(r, c); }
	};


	template<>
	struct MPS_<MAT_PROP::SYM, MAT_SHAPE::LT>
	{
		static inline void swap(index_type& r, index_type& c){ if (c > r) std::swap(r, c); }
	};

	template<>
	struct MPS_<MAT_PROP::SYM, MAT_SHAPE::UT>
	{
		static inline void swap(index_type& r, index_type& c){ if (c < r) std::swap(r, c); }
	};



	// there are many others ... left
	
	


	template<typename T, typename MP, typename MSH>
	class Matrixe :public Mat<T,MSH>
	{
		typedef Mat<T, MSH>				Mty_;
		typedef MPS_<MP, MSH>			ICT_;
	public:
		typedef MP						property_t;
		typedef MSH					    shape_t;
		typedef Matrixe					Matrix;
		
		//------------------------------------//
		//           CTORS                    //
		//------------------------------------//
		Matrix() = default;
		Matrix(const Matrix& rval) = default;
		Matrix(Matrix&& rval) :Mty_(std::move(rval)){}
		explicit Matrix(index_t nr, index_t nc=1) :Mty_(nr,nc){}
		Matrix(const std::initializer_list<T>& s){ assign(s); }
		Matrix(const std::initializer_list<std::initializer_list<T>>& s){ assign(s); }


		Matrix& resize(index_t nr, index_t nc=1){Mty_::resize(nr, nc); return *this; }
		void	swap(Matrix& o){ Mty_::swap(o); }
				
		//------------------------------------//
		//      ASSIGNMENT                    //
		//------------------------------------//
		Matrix&		operator=(const Matrix&  rhs) = default;
		Matrix&		operator=(Matrix&&  o)  { Mty_::operator=(std::move(o)); return *this; }
		Matrix&		operator=(value_t  val) { fill(val); return *this;}
		Matrix&		operator=(const mat_T<T>& o);
		Matrix&		operator=(const std::initializer_list<T>& s){return assign(s);}
		Matrix&		operator=(const std::initializer_list<std::initializer_list<T>>& s){ return assign(s); }
		Matrix&		assign(const std::initializer_list<std::initializer_list<T>> & s);
		Matrix&		assign(const std::initializer_list<T> & s);

		//------------------------------------//
		//      ACCESSORS                     //
		//------------------------------------//

		// WARNING!!! direct access to the first element of the r-th row (0-base)
		//            for REctangular and Lower Triangular matrices the first element of row r has column index 0 (0-b)
		//            but for Upper Triangular matrices it hase  column index r (0-b)
		T&			operator[](index_t r){ return at_0b(r);}
		const T&	operator[](index_t r)const { return at_0b(r); }
		
		// 1-based indices
		T&			operator()(index_t r, index_t c);
		const T&	operator()(index_t r, index_t c)const;
		T& 			operator()(index_t r);
		const T& 	operator()(index_t r)const;


		//------------------------------------//
		//      SUBMATRIX ACCESSORS           //
		//------------------------------------//
		/*
		SubMatrix_t		    operator()(const iSet& 	r, const iSet& 	c)	{ return SubMatrix_t(*this, r, c); }
		SubMatrix_t		    operator()(iSet&& 		r, const iSet&  c)	{ return SubMatrix_t(*this, std::move(r), c); }
		SubMatrix_t		    operator()(const iSet& 	r, iSet&& 		c)	{ return SubMatrix_t(*this, r, std::move(c)); }
		SubMatrix_t		    operator()(iSet&& 		r, iSet&& 		c)	{ return SubMatrix_t(*this, std::move(r), std::move(c)); }

		const_SubMatrix_t	operator()(const iSet& r, const iSet& c)const	{ return const_SubMatrix_t(*this, r, c); }
		const_SubMatrix_t	operator()(iSet&& r, const iSet& c)const		{ return const_SubMatrix_t(*this, std::move(r), c); }
		const_SubMatrix_t	operator()(const iSet& r, iSet&& c)const		{ return const_SubMatrix_t(*this, r, std::move(c)); }
		const_SubMatrix_t	operator()(iSet&& r, iSet&& c)const				{ return const_SubMatrix_t(*this, std::move(r), std::move(c)); }

		SubMatrix_t			operator()(const iSet&  r)			{ return SubMatrix_t(*this, r, BASE::first()); }
		SubMatrix_t			operator()(const iSet&& r)			{ return SubMatrix_t(*this, std::move(r), BASE::first()); }
		SubMatrix_t			row(const iSet&	r)              	{ return SubMatrix_t(*this, r, iSet(BASE::first(), BASE::last(nc_))); }
		SubMatrix_t			row(iSet&&	r)              		{ return SubMatrix_t(*this, std::move(r), iSet(BASE::first(), BASE::last(nc_))); }
		SubMatrix_t			column(const iSet& c)           	{ return SubMatrix_t(*this, iSet(BASE::first(), BASE::last(nr_)), c); }
		SubMatrix_t			column(iSet&&	c)              	{ return SubMatrix_t(*this, iSet(BASE::first(), BASE::last(nr_)), std::move(c)); }
		const_SubMatrix_t	operator()(const iSet&  r)const		{ return const_SubMatrix_t(*this, r, BASE::first()); }
		const_SubMatrix_t	operator()(const iSet&& r)const		{ return const_SubMatrix_t(*this, std::move(r), BASE::first()); }
		// 1D array access
		const_SubMatrix_t	row(const iSet&	r)const     		{ return const_SubMatrix_t(*this, r, iSet(BASE::first(), BASE::last(nc_))); }
		const_SubMatrix_t	row(iSet&& r)const          		{ return const_SubMatrix_t(*this, std::move(r), iSet(BASE::first(), BASE::last(nc_))); }
		const_SubMatrix_t	column(const iSet& c)const			{ return const_SubMatrix_t(*this, iSet(BASE::first(), BASE::last(nr_)), c); }
		const_SubMatrix_t	column(iSet&&	c)const				{ return const_SubMatrix_t(*this, iSet(BASE::first(), BASE::last(nr_)), std::move(c)); }
		*/
		
		//------------------------------------//
		//           IO                       //
		//------------------------------------//
		int		fwrite(const char* fname);
		int		fread(const char* fname);
		//------------------------------------//
		//  MATH  OPS                   //
		//------------------------------------//
		Matrix& operator+=(const Matrix& o){ data_ += o.data_;	return *this; }
		Matrix& operator-=(const Matrix& o){ data_ -= o.data_;	return *this; }
		Matrix& operator+(){ return *this; }
		Matrix  operator-()const{ Matrix tmp(*this); tmp.data_=-data_;	return tmp; }
		Matrix& operator*=(const T& val){ data_ *= val;	return *this; }
		Matrix& operator/=(const T& val){ data_ /= val;	return *this; }
		//returns a copy of the Matrix transformed by func
		Matrix	apply(T func(T)) const{ Matrix tmp(*this); tmp.data_.apply(func);return tmp; }
		Matrix& transform(T func(T)) const{ data_.transform(func); return *this; }


		

	};






	template<typename T, typename MP, typename MS, typename MO>
	class mat_T_P_S_O :public mat_T_O_S<T,MO, MS>
	{
		typedef mat_T_O_S<T, MO, MS>	Mty_;
		typedef MPS_<MP, MS>			ICT_;
	public:
		typedef MP						property_t;
		typedef MS					    shape_t;
		typedef MO					    order_t;
		typedef mat_T_P_S_O				Matrix;

		//------------------------------------//
		//           CTORS                    //
		//------------------------------------//
		Matrix() = default;
		Matrix(const Matrix& rval) = default;
		Matrix(Matrix&& rval) :Mty_(std::move(rval)){}

		explicit Matrix(index_t nr, index_t nc = 1) :Mty_(nr, nc,get_size_(nr,nc)){}
		Matrix(const std::initializer_list<T>& s){ assign(s); }
		Matrix(const std::initializer_list<std::initializer_list<T>>& s){ assign(s); }


		Matrix& resize(index_t nr, index_t nc = 1){ Mty_::resize(nr, nc); return *this; }
		void	swap(Matrix& o){ Mty_::swap(o); }

		//------------------------------------//
		//      ASSIGNMENT                    //
		//------------------------------------//
		Matrix&		operator=(const Matrix&  rhs) = default;
		Matrix&		operator=(Matrix&&  o)  { Mty_::operator=(std::move(o)); return *this; }
		Matrix&		operator=(value_t  val) { fill(val); return *this; }
		

		Matrix&		operator=(const std::initializer_list<T>& s){ return assign(s); }
		Matrix&		operator=(const std::initializer_list<std::initializer_list<T>>& s){ return assign(s); }
		Matrix&		assign(const std::initializer_list<std::initializer_list<T>> & s);
		Matrix&		assign(const std::initializer_list<T> & s);

		//------------------------------------//
		//      ACCESSORS                     //
		//------------------------------------//

		// WARNING!!! direct access to the first element of the r-th row (0-base)
		//            for REctangular and Lower Triangular matrices the first element of row r has column index 0 (0-b)
		//            but for Upper Triangular matrices it hase  column index r (0-b)
		T&			operator[](index_t r){ return at_0b(r); }
		const T&	operator[](index_t r)const { return at_0b(r); }

		// 1-based indices
		T&			operator()(index_t r, index_t c);
		const T&	operator()(index_t r, index_t c)const;
		T& 			operator()(index_t r);
		const T& 	operator()(index_t r)const;


		//------------------------------------//
		//      SUBMATRIX ACCESSORS           //
		//------------------------------------//
		/*
		SubMatrix_t		    operator()(const iSet& 	r, const iSet& 	c)	{ return SubMatrix_t(*this, r, c); }
		SubMatrix_t		    operator()(iSet&& 		r, const iSet&  c)	{ return SubMatrix_t(*this, std::move(r), c); }
		SubMatrix_t		    operator()(const iSet& 	r, iSet&& 		c)	{ return SubMatrix_t(*this, r, std::move(c)); }
		SubMatrix_t		    operator()(iSet&& 		r, iSet&& 		c)	{ return SubMatrix_t(*this, std::move(r), std::move(c)); }

		const_SubMatrix_t	operator()(const iSet& r, const iSet& c)const	{ return const_SubMatrix_t(*this, r, c); }
		const_SubMatrix_t	operator()(iSet&& r, const iSet& c)const		{ return const_SubMatrix_t(*this, std::move(r), c); }
		const_SubMatrix_t	operator()(const iSet& r, iSet&& c)const		{ return const_SubMatrix_t(*this, r, std::move(c)); }
		const_SubMatrix_t	operator()(iSet&& r, iSet&& c)const				{ return const_SubMatrix_t(*this, std::move(r), std::move(c)); }

		SubMatrix_t			operator()(const iSet&  r)			{ return SubMatrix_t(*this, r, BASE::first()); }
		SubMatrix_t			operator()(const iSet&& r)			{ return SubMatrix_t(*this, std::move(r), BASE::first()); }
		SubMatrix_t			row(const iSet&	r)              	{ return SubMatrix_t(*this, r, iSet(BASE::first(), BASE::last(nc_))); }
		SubMatrix_t			row(iSet&&	r)              		{ return SubMatrix_t(*this, std::move(r), iSet(BASE::first(), BASE::last(nc_))); }
		SubMatrix_t			column(const iSet& c)           	{ return SubMatrix_t(*this, iSet(BASE::first(), BASE::last(nr_)), c); }
		SubMatrix_t			column(iSet&&	c)              	{ return SubMatrix_t(*this, iSet(BASE::first(), BASE::last(nr_)), std::move(c)); }
		const_SubMatrix_t	operator()(const iSet&  r)const		{ return const_SubMatrix_t(*this, r, BASE::first()); }
		const_SubMatrix_t	operator()(const iSet&& r)const		{ return const_SubMatrix_t(*this, std::move(r), BASE::first()); }
		// 1D array access
		const_SubMatrix_t	row(const iSet&	r)const     		{ return const_SubMatrix_t(*this, r, iSet(BASE::first(), BASE::last(nc_))); }
		const_SubMatrix_t	row(iSet&& r)const          		{ return const_SubMatrix_t(*this, std::move(r), iSet(BASE::first(), BASE::last(nc_))); }
		const_SubMatrix_t	column(const iSet& c)const			{ return const_SubMatrix_t(*this, iSet(BASE::first(), BASE::last(nr_)), c); }
		const_SubMatrix_t	column(iSet&&	c)const				{ return const_SubMatrix_t(*this, iSet(BASE::first(), BASE::last(nr_)), std::move(c)); }
		*/

		//------------------------------------//
		//           IO                       //
		//------------------------------------//
		int		fwrite(const char* fname);
		int		fread(const char* fname);
		//------------------------------------//
		//  MATH  OPS                   //
		//------------------------------------//
		Matrix& operator+=(const Matrix& o){ data_ += o.data_;	return *this; }
		Matrix& operator-=(const Matrix& o){ data_ -= o.data_;	return *this; }
		Matrix& operator+(){ return *this; }
		Matrix  operator-()const{ Matrix tmp(*this); tmp.data_ = -data_;	return tmp; }
		Matrix& operator*=(const T& val){ data_ *= val;	return *this; }
		Matrix& operator/=(const T& val){ data_ /= val;	return *this; }
		//returns a copy of the Matrix transformed by func
		Matrix	apply(T func(T)) const{ Matrix tmp(*this); tmp.data_.apply(func); return tmp; }
		Matrix& transform(T func(T)) const{ data_.transform(func); return *this; }




	};

























	////////////////////////////////////////////////////////////
	/// IMPLEMENTATION

	template<typename T, typename MP, typename MS>
	T&	Matrixe<T,MP,MS>::operator()(index_t r, index_t c)
	{
		ICT_::swap(r, c);
		return at(r, c);
	}

	template<typename T, typename MP, typename MS>
	T&	Matrixe<T, MP, MS>::operator()(index_t r)
	{
		static_assert(std::is_same<MSH, MAT_SHAPE::RE>::value, "");
		return at(r);
	}

	

	template<typename T, typename MP, typename MS>
	Matrixe<T, MP, MS>&  Matrixe<T, MP, MS>::assign(const std::initializer_list<std::initializer_list<T>> & s)
	{
		size_t nr = s.size();
		size_t nc = 0;
		// determino il numero di colonne della matrice in base al numero massimo di colonne 
		// delle varie righe
		for (size_t r = 0; r != nr; ++r)
		{
			const std::initializer_list<T>* r_ptr = s.begin() + r;
			nc = std::max(nc, r_ptr->size());
		}

		resize(nr, nc);

		for (size_t r = 0; r != nr; ++r)
		{
			const std::initializer_list<T>* r_ptr = s.begin() + r;
			size_t this_row_len = end_col_index(r) - first_col_index(r);

			size_t cpy_sz = std::min(r_ptr->size(),this_row_len);

			std::copy(r_ptr->begin(), r_ptr->begin() + cpy_sz, p_row_0b(r));
		}
		return *this;
	}


	template<typename T, typename MP, typename MS>
	Matrixe<T, MP, MS>&  Matrixe<T, MP, MS>::assign(const std::initializer_list<T> & s)
	{
		size_t nr = s.size();
		resize(nr, 1);
		std::copy(s.begin(), s.end(), p_row_0b(0));
		return *this;
	}


	/// fwrite/fread : Scrivono e Leggono la matrice su/da un file binario
	///@param fname nome del file
	///@return:  0 in caso di successo
	///@return: -1 nel caso il file non possa essere aperto
	///@return: -2 se si sono verificati errori in lettura
	///@return: -3 (solo fread) se i tipi  del'indice o dei valori su file non coincidono con quelli della Matrix
	template<typename T, typename MP, typename MS>
	int Matrixe<T, MP,MS>::fwrite(const char* fname)
	{
		std::ofstream f(fname, std::ios_base::binary);
		if (!f.is_open())
			return -1;
		// l'intestazione è costituita da 6 uint64 che indicano:
		// numero di  righe, 
		// numero di colonne
		// ID del tipo degli indici
		// ID del tipo delle componenti
		// ID di proprieta' 
		// ID di forma
		std::uint64_t type_[6] = { nr_, nc_, std::uint64_t(type<index_t>::id()), type<T>::id(), MP::ID, MS::ID };
		f.write(reinterpret_cast<const char*>(type_), sizeof(type_));
		f.write(reinterpret_cast<const char*>(begin()), size()*sizeof(T));
		return f.good() ? 0 : -2;
	}


	template<typename T, typename MP, typename MS>
	int Matrixe<T, MP, MS>::fread(const char* fname)
	{
		std::ifstream f(fname, std::fstream::binary);
		if (!f.is_open())
			return -1;
		// leggo l'intestazione
		std::uint64_t type_[6];
		f.read(reinterpret_cast<char*>(type_), sizeof(type_));
		if (!f.good())
			return -2;
		// controllo di consistenza dei tipi
		if (type_[2] != std::uint64_t(type<index_t>::id()) ||
			type_[3] != std::uint64_t(type<T>::id()) ||
			type_[4] != std::uint64_t(MP::ID)||
			type_[5] != std::uint64_t(MS::ID) ||
			)
			return -3;
		resize(index_t(type_[0]), index_t(type_[1]));
		f.read(reinterpret_cast<char*>(data_.begin()), size()*sizeof(T));
		return f.good() ? 0 : -2;
	}





}// end namespace