#ifndef _MML3_DIAGONAL_BLOCK_MATRIX_H_
#define _MML3_DIAGONAL_BLOCK_MATRIX_H_

#include"MML3-Matrix.h"
#include<vector>

namespace MML3
{

	template<typename BMAT>
	struct Elem_
	{
		Elem_() = default;

		~Elem_()
		{
			if (allocated_)
				delete mat_;
		}

		void copy(const BMAT& o)
		{
			if (allocated_)
				delete mat_;
			mat_ = new BMAT(o);
			allocated_ = true;
		}
		void reference(const BMAT& o)
		{
			if (allocated_)
			{
				allocated_ = false;
				delete mat_;
			}
			mat_ = &o;
		}

		const BMAT& get()const{ return *mat_;}


	private:
		bool  allocated_ = false;
		const BMAT* mat_ = nullptr;

		
	};
	

	template<typename BMAT>
	class DiagonalBlockMatrix
	{
	public:
		typedef typename BMAT::value_t		value_t;

		DiagonalBlockMatrix(size_t nb):block_(nb){}


		size_t	nblocks()const{ return block_.size(); }
		const BMAT& operator ()(size_t i)const{ return block_[i - 1].get(); }
		const BMAT& operator [](size_t i)const{ return block_[i].get(); }
		
		DiagonalBlockMatrix&   copy_block(size_t i, const BMAT& block)
		{
			block_[i - 1].copy(block);
			return *this;
		}

		DiagonalBlockMatrix&   reference_block(size_t i, const BMAT& block)
		{
			block_[i - 1].reference(block);
			return *this;
		}

		size_t		nrows()const
		{
			size_t count = 0;
			for (auto& b : block_)
				count += b.get().nrows();
			return count;
		}

		size_t		ncols()const
		{
			size_t count = 0;
			for (auto& b : block_)
				count += b.get().ncols();
			return count;
		}

		

		
	private:
		std::vector<Elem_<BMAT>> block_;
	};


	namespace Algorithm{

		//C= Diag * A  
		template<typename BMAT, typename MAT>
		MAT&  product_Diag_A(const DiagonalBlockMatrix<BMAT>& D, const MAT& A, MAT& C)
		{

			typedef typename BMAT::value_t value_t;
			size_t Dsz = D.ncols();
			size_t Anc = A.ncols();
			if (Dsz != A.nrows())
				throw std::range_error("product_Diag_A(...): incompatible dimensions");

			assert(&A != &C);
			try
			{
				C.resize(Dsz, Anc);
			}
			catch (...)
			{
				throw std::range_error("product_Diag_A(...): destination matrix C not resizeable");
			}

			size_t nblocks = D.nblocks();

			size_t I = 0;
			for (size_t i = 1; i <= nblocks; ++i)
			{
				const BMAT& Di = D(i);
				size_t bsize = Di.nrows();
				for (size_t a = 1; a <= bsize; ++a)
				{
					for (size_t b = 1; b <= Anc; ++b)
					{
						value_t acc = 0.;
						for (size_t c = 1; c <= bsize; ++c)
							acc += Di(a, c)*A(I + c, b);
						C(I + a, b) = acc;
					}
				}
				I += bsize;
			}
			return C;
		}



		//C= Diag' * A  
		template<typename BMAT, typename MAT>
		MAT&  product_DiagT_A(const DiagonalBlockMatrix<BMAT>& D, const MAT& A, MAT& C=MAT())
		{

			typedef typename BMAT::value_t value_t;
			size_t Dsz = D.ncols();
			if (Dsz != A.nrows())
				throw std::range_error("product_DiagT_A(...): incompatible dimensions");
			assert(&A != &C);
			try
			{
				C.resize(Dsz, A.ncols());
			}
			catch (...)
			{
				throw std::range_error("product_DiagT_A(...): destination matrix C not resizeable");
			}


			size_t nblocks = D.nblocks();

			size_t I = 1;
			for (size_t i = 0; i != nblocks; ++i)
			{
				const BMAT& Di = D(i + 1);
				size_t bsize = Di.nrows();
				size_t J = 1;
				for (size_t j = 0; j != nblocks; ++j)
				{
					for (size_t a = 0; a != bsize; ++a)
					{
						for (size_t b = 0; b != bsize; ++b)
						{
							value_t acc = 0.;
							for (size_t c = 0; c != bsize; ++c)
								acc += Di(c + 1, a + 1)*A(I + c, J + b);
							C(I + a, J + b) = acc;
						}
					}


					J += bsize;
				}
				I += bsize;
			}
			return C;

		}

		// C = Diag' * A * Diag
		template<typename BMAT, typename MAT>
		MAT& product_DiagT_A_Diag(const DiagonalBlockMatrix<BMAT>& D, const MAT& A, MAT& C=MAT())
		{

			typedef typename BMAT::value_t value_t;
			size_t Dsz = D.ncols();

			if (Dsz != A.nrows() || Dsz != A.ncols())
				throw std::range_error("product_DiagT_A_Diag(const DiagonalBlockMatrix<BMAT>& D, const MAT& B, MAT& C): incompatible dimensions");

			assert(&A != &C);
			try
			{
				C.resize(Dsz, Dsz);
			}
			catch (...)
			{
				throw std::range_error("product_DiagT_A_Diag(...): destination matrix C not resizeable");
			}


			size_t nblocks = D.nblocks();

			size_t I = 1;
			for (size_t i = 0; i != nblocks; ++i)
			{
				const BMAT& Di = D(i + 1);
				size_t bsize = Di.nrows();
				size_t J = 1;
				for (size_t j = 0; j != nblocks; ++j)
				{
					for (size_t a = 0; a != bsize; ++a)
					{
						for (size_t b = 0; b != bsize; ++b)
						{
							value_t acc = 0;
							for (size_t c = 0; c != bsize; ++c)
								for (size_t d = 0; d != bsize; ++d)
									acc += Di(c + 1, a + 1)*A(I + c, J + d) * Di(d + 1, b + 1);
							C(I + a, J + b) = acc;
						}
					}
					J += bsize;
				}
				I += bsize;
			}
			return C;
		}

	}// end namespace Algorithm
	


} // end namespace MML3

#endif