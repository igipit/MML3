#ifndef _MML3_CSR3_MATRIX_H_
#define _MML3_CSR3_MATRIX_H_
#include<cstdio>
#include<iostream>
#include<iomanip>
#include<vector>
#include<limits>
#include"MML3-Matrix.h"
#include"MML3-basic_sparse_matrix.h"


namespace MML3
{


	template<  typename	VALUE_TYPE,	/* il tipo delle componenti della matrice    */
               typename IDX_TYPE= int_t,    /* il tipo intero degli indici della matrice */
               class    MP=M_PROP::GE>
	class sparseCSR3Matrix;

	template<typename T, typename IDX_TYPE = int_t>
	using CSR3_matrix = sparseCSR3Matrix<T, IDX_TYPE, M_PROP::GE>;

	template<typename T, typename IDX_TYPE = int_t>
	using CSR3_sym_matrix = sparseCSR3Matrix<T, IDX_TYPE, M_PROP::SYM>;


/** 
Static    Compressed Sparse Row Matrix (3 vector format)
Static refers to the sparsity structure:  once the sparsity structure is formed it is not possible to add compoments 
*/
template<typename VALUE_TYPE,typename IDX_TYPE,class MP>
class sparseCSR3Matrix : public basic_sparse_matrix<VALUE_TYPE,IDX_TYPE>
{
public:
	typedef VALUE_TYPE									value_t;
	typedef IDX_TYPE									index_t;
	// Since row_position_array contains the offset position of the beginning element of each row, it should be
	// of type size_t to conform to the machine architecture (32 bit and 64 bit have consistent differences in addressable memory space)
	// But, since sparse solvers (pardiso) expect the same type for column indexes and row positions  here i forced this constraint
	// and thus the maximum nonzeros of a matrix equals std::numeric_limits<index_t>::max() that is 2'147'483'647 if an int32 is used
	typedef Vector<index_t>						row_position_array_t;

	typedef Vector<index_t>						column_idx_array_t;
	typedef Vector<value_t>						column_val_array_t;

	

	
	class col_iterator;
	class const_col_iterator;

	


	enum IS :bool {
		RowMajor = true,
		ColMajor = false,
		SYM = std::is_same<M_PROP::SYM, MP>::value,
		GE=!SYM,
		RE = false,
		LT = false,
		UT = SYM // if it is symmetric it is also Upper triangular
	};

	//-------------------------------------------------------------------------
	// COMMON INTERFACE
	//-------------------------------------------------------------------------
	typedef basic_sparse_matrix<VALUE_TYPE, IDX_TYPE>			base_matrix_t;
	typedef iSet										iset_t;
	typedef base_matrix_t::inserter_matrix_t			inserter_matrix_t;
	typedef base_matrix_t::inserter_sym_matrix_t		inserter_sym_matrix_t;

	
	sparseCSR3Matrix() = default;
	sparseCSR3Matrix(const sparseCSR3Matrix& lhs) = default;
	~sparseCSR3Matrix(){ destroy(); }

	void				destroy();
	index_t				nrows()const{ return nrows_; }
	index_t				ncols()const{ return ncols_; }
	bool				is_row_major()const													{ return IS::RowMajor; }
	bool				is_symmetric()const													{ return IS::SYM; }
	size_t				nonzeros()const{ return col_val_.size(); }
	bool				test_subscripts(index_t r, index_t c)const;
	void				fill(const value_t& val)											{ col_val_=val; }

	
	value_t&			put(index_t i, index_t j, const value_t& val);
	value_t&			sum(index_t i, index_t j, const value_t& val);
	/// Ritorna il puntatore alla componente in posizione (i,j) se questa esiste, nullptr altrimenti
	const value_t*		get_p(index_t i, index_t j)const;
	value_t*			get_p(index_t i, index_t j);
	
	void				put(const iSet& ir, const iSet& ic, const inserter_matrix_t& K)		{ put_(ir, ic, K); }

	void				sum(const iSet& ir, const iSet& ic, const inserter_matrix_t& K)		{ sum_(ir, ic, K);}
	void				sum(const iSet& irc, const inserter_matrix_t& K)					{ sum_(irc, K); }
	void				sum(const iSet& irc, const inserter_sym_matrix_t& K)				{ sum_(irc, K); }

	void				print(std::ostream& os)const;
	int					fwrite(const std::string& fname)const;
	int					fread(const std::string&);
	//-------------------------------------------------------------------------
	// EXTENDED INTERFACE
	//-------------------------------------------------------------------------

	/// costruttore che alloca le risorse necessarie
	sparseCSR3Matrix(index_t nr, index_t nc, size_t nnz);
	/// ridimensiona la matrice perdendo contenuto e struttura
	void resize(index_t nr, index_t nc, size_t nnz);
	void swap(sparseCSR3Matrix& lhs);
	
	///------------------------------------------------------ 
	/// per accedere dall'esterno direttamante ai 3 vettori
	///------------------------------------------------------
	/// ritorna il puntatore agli indici di colonna
	index_t*	column_index(){return &(col_idx_[0]);}
	/// ritorna il puntatore all'vettore di puntatori alle righe 
	index_t*	row_pos(){return & row_pos_[0];}
	/// ritorna il puntatore ai valori diversi da zero
	value_t*	column_value(){return & col_val_[0];}





	col_iterator		row_begin(index_t r){ return col_iterator(col_idx_, col_val_, row_pos_[r-1]-1); }
	col_iterator		row_end(index_t r){ return col_iterator(col_idx_, col_val_, row_pos_[r]-1); }
	const_col_iterator	row_begin(index_t r)const{ return const_col_iterator(col_idx_, col_val_, row_pos_[r - 1] - 1); }
	const_col_iterator	row_end(index_t r)const{ return const_col_iterator(col_idx_, col_val_, row_pos_[r] - 1); }

	/// determina il valore massimo e minimo sulla diagonale della matrice
	///@param row_idx_min indice di riga del valore minimo
	///@param row_idx_max indice di riga del valore massimo
	///@param val_min     valore minimo
	///@param val_max     valore massimo
	//void	get_min_max_diag(index_t& row_idx_min, value_type& val_min , index_t& row_idx_max, value_type& val_max)const;


	/// determina il valore massimo e minimo  della matrice e la loro posizione
	//void	get_min_max(index_t& row_idx_min, index_t& col_idx_min, value_type& val_min , 
	//									   index_t& row_idx_max, index_t& col_idx_max, value_type& val_max)const;

	
	//int		copy2(index_t rpos_sz, index_t nnz, index_t rpos[rpos_sz], index_t cidx[nnz], value_t  val[nnz])const;
	

	///@name OPERATORS
	//@{
	/// prodotto matrice vettore
	void product(const value_t* B, index_t szB, value_t* C, index_t szC)const;


	class col_iterator
	{
	public:
		col_iterator(const column_idx_array_t& idxp, column_val_array_t& valp, index_t pos) :col_idx_(idxp), col_val_(valp), pos_(pos){}
		col_iterator() = delete;
		col_iterator(const col_iterator&) = default;

		col_iterator&		operator++(){ ++pos_; return *this; }
		bool				operator!=(const col_iterator& o)const{ return pos_ != o.pos_; }
		bool				operator!=(const const_col_iterator& o)const{ return pos != o.pos_; }
		const	index_t&	idx()const{ return col_idx_[pos_]; }
		value_t&			value(){ return col_val_[pos_]; }
	private:
		const column_idx_array_t&	col_idx_ ;
		column_val_array_t&		col_val_;
		index_t pos_ = 0;

	};

	class const_col_iterator
	{
		friend class col_iterator;
	public:
		const_col_iterator(const column_idx_array_t idxp, const column_val_array_t& valp, index_t pos) :col_idx_(idxp), col_val_(valp), pos_(pos){}
		const_col_iterator() = delete;
		const_col_iterator(const const_col_iterator&) = default;

		const_col_iterator&		operator++(){ ++pos_; }
		bool					operator!=(const col_iterator& o)const{ return pos != o.pos_; }
		const	index_t&		idx()const{ return	col_idx_[pos_]; }
		const   value_t&		value(){ return		col_val_[pos]; }
	protected:
		const column_idx_array_t&  col_idx_;
		const column_val_array_t&	col_val_;
		index_t pos_ = 0;

	};


private:
	index_t			nrows_ = 0, ncols_ = 0;
	//size_t			nnz_ = 0;
	row_position_array_t	row_pos_;
	column_idx_array_t		col_idx_;
	column_val_array_t		col_val_;

	// row_pos_[] ha dimensione nrows+1 e contiene la posizione ( a partire da 1) in col_idx_ e col_val_ del primo elemento di ogni riga:
	// row_pos_[0]-1 � la posizione in  col_idx_ del primo indice di colonna della riga 0
	// row_pos_[1]-1 � la posizione in  col_idx_ del primo indice di colonna della riga 1
	// ...
	// row_pos_[nrows()] contiene il numero totale di elementi diversi da zero
	// il numero di elementi diversi da zero nella riga i-esima (0-base) � sempre dato da rsz_[i+1] - rsz_[i]
	
	
	/// ritorna la posizione 1-based della componente (r,c) all'interno di col_idx se l'elemento esiste, ritorna 0  altrimenti
	size_t				find_position_1b_(index_t r, index_t c)const;

	template<typename MAT>
	void				put_(const iSet& ir, const iSet& ic, const MAT& K);

	template<typename MAT>
	void				sum_(const iSet& ir, const iSet& ic, const MAT& K);
	
	template<typename MAT>
	void				sum_(const iSet& icr, const MAT& K);
	
	

};


//@}


////////////////////////////////////////////////////////////////////////////////
// Implementazione


template<typename	VAL, typename	IDX, typename MP >
inline bool sparseCSR3Matrix<VAL, IDX, MP>::test_subscripts(index_t r, index_t c)const
{
	if ((r > nrows_ || r<1) || (c > ncols_ || c<1))
		return false;
	return true;
}



/// costruttore che alloca le risorse necessarie
template<typename	VAL, typename	IDX, typename MP >
sparseCSR3Matrix<VAL, IDX, MP>::sparseCSR3Matrix(index_t nr, index_t nc, size_t nnz) :
		nrows_(nr),
		ncols_(nc),
		row_pos_(nr+1),
		col_idx_(nnz),
		col_val_(nnz)
	{
		
	}


	/// deallocates resources
template<typename	VAL, typename	IDX, typename MP >
void sparseCSR3Matrix<VAL, IDX, MP>::destroy()
{
	row_pos_.clear();
	col_idx_.clear();
	col_val_.clear();
	nrows_=0;
	ncols_=0;
}



/// ritorna la posizione 1-based della componente (r,c) all'interno di col_idx se l'elemento esiste, ritorna 0  altrimenti
/// essendo gli indici di colonna ordinati (per ogni riga) si pu� fare una ricerca piu' efficente, perlomeno quando gli 
/// elementi sono molti. Probabilmente con pochi elementi per riga la ricerca sequanziale e' la piu' veloce
template<typename	VAL, typename	IDX, typename MP >
size_t sparseCSR3Matrix<VAL, IDX, MP>::find_position_1b_(index_t r, index_t c)const
{
	// does not test indexes correctness
	size_t pos		= row_pos_[r - 1] - 1;
	size_t i_end	= row_pos_[r] - 1;
	for (; pos != i_end && col_idx_[pos]<c; ++pos);
	if (pos != i_end && col_idx_[pos] == c)
		return pos+1;
	else
		return 0;
}


/// ridimensiona la matrice perdendo contenuto e struttura
template<typename	VAL, typename	IDX, typename MP >
void sparseCSR3Matrix<VAL, IDX, MP>::resize(index_t nr, index_t nc, size_t nnz)
	{
		sparseCSR3Matrix	tmp(nr,nc,nnz);
		swap(tmp);
	}


template<typename	VAL, typename	IDX, typename MP >
void sparseCSR3Matrix<VAL, IDX, MP>::swap(sparseCSR3Matrix& rhs)
{
		std::swap(nrows_,rhs.nrows_);
		std::swap(ncols_,rhs.ncols_);
		row_pos_.swap(rhs.row_pos_);
		col_idx_.swap(rhs.col_idx_);
		col_val_.swap(rhs.col_val_);
}

	

/// Ritorna il puntatore alla componente in posizione (i,j) se questa esiste, nullptr altrimenti
template<typename	VAL, typename	IDX, typename MP >
const VAL*  sparseCSR3Matrix<VAL, IDX, MP>::get_p(index_t i, index_t j)const
{
	if (!test_subscripts(i, j))
		throw std::out_of_range("static_sparse_CR_Matrix<VAL,IDX>::get()const");
	size_t pos = find_position_1b_(i, j);
	if (pos)
		return & col_val_[pos - 1];
	else
		return nullptr;
}


/// Ritorna il puntatore alla componente in posizione (i,j) se questa esiste, nullptr altrimenti
template<typename	VAL, typename	IDX, typename MP >
VAL*  sparseCSR3Matrix<VAL, IDX, MP>::get_p(index_t i, index_t j)
{
	if (!test_subscripts(i, j))
		throw std::out_of_range("static_sparse_CR_Matrix<VAL,IDX>::get()const");
	size_t pos = find_position_1b_(i, j);
	if (pos)
		return &col_val_[pos - 1];
	else
		return nullptr;
}



/*template<int MT, typename VAL, typename	IDX>
void	static_sparse_CR_Matrix<MT,VAL,IDX>::get_min_max_diag(index_t& row_idx_min, value_type& val_min , index_t& row_idx_max, value_type& val_max)const
{
		row_idx_min=0;
		row_idx_max=0;
		val_min=std::numeric_limits<VAL>::max();
		val_max=std::numeric_limits<VAL>::min();
		VAL value;
		
		for (index_t i = 0; i != nrows_; ++i)
		{
			if( get(i+1,i+1,value)==end_pos())
				value=VAL(0);

			if(value > val_max)
			{
				val_max=value;
				row_idx_max=i+1;
			}
			if(value < val_min)
			{
				val_min=value;
				row_idx_min=i+1;
			}

		}

}

template<int MT, typename VAL, typename	IDX>
void	static_sparse_CR_Matrix<MT,VAL,IDX>::get_min_max(index_t& row_idx_min, index_t& col_idx_min, value_type& val_min , 
										   index_t& row_idx_max, index_t& col_idx_max, value_type& val_max)const
{
		row_idx_min=0;
		col_idx_min=0;
		row_idx_max=0;
		col_idx_max=0;
		val_min=std::numeric_limits<VAL>::max();
		val_max=std::numeric_limits<VAL>::min();
		//VAL value;

		
		int pos_min=MML::imin_(&col_val_[0],col_val_.size());
		val_min=col_val_[pos_min];
		col_idx_min=col_idx_[pos_min];

		for(int i=0; i!=nrows_; ++i)
			if(row_pos_[i+1] > pos_min)
			{
				row_idx_min=i+1;
				break;
			}


		int pos_max=MML::imax_(&col_val_[0], col_val_.size());
		val_max=col_val_[pos_max];
		col_idx_max=col_idx_[pos_max];
		for(int i=0; i!=nrows_; ++i)
			if(row_pos_[i+1] > pos_max)
			{
				row_idx_max=i+1;
				break;
			}

}
*/

	
template<typename	VAL, typename	IDX, typename MP>
auto sparseCSR3Matrix<VAL, IDX,MP>::put(index_t R, index_t C, const value_t& val)->value_t&
{
	value_t* p=get_p(R, C);
	if(p)
		*p=val;
	else
		throw std::out_of_range("sparseCSR3Matrix<VAL,IDX>::put():2");
	return *p;
}


template<typename	VAL, typename	IDX, typename MP >
auto	sparseCSR3Matrix<VAL, IDX,MP>::sum(index_t R, index_t C, const value_t& val)->value_t& 
{
	value_t* p = get_p(R, C);
	if (p)
		*p += val;
	else
		throw std::out_of_range("sparseCSR3Matrix<VAL,IDX>::sum():2");
	return *p;
}






	

template<typename	VAL, typename	IDX, typename MP >
template<typename MAT>
void sparseCSR3Matrix< VAL, IDX, MP>::put_(const iSet& ir, const iSet& ic, const MAT& K)
{
	if (index_t(ir.max()) > nrows() || index_t(ic.max())>ncols())
		throw std::out_of_range("sparseCSR3Matrix<VAL,IDX>::put_():6");
	index_t rsz = (index_t)ir.size(),
		csz = (index_t)ic.size();

	for (index_t r = 1; r <= rsz; ++r)
	{
		index_t R = ir(r);
		if (!R)
			continue;

		for (index_t c = 1; c <= csz; ++c)
		{
			index_t C = ic(c);
			if (!C || (IS::SYM && C<R))
				continue;
			value_t * v = get_p(R, C);
			if (v)
				*v = K(r, c);
			else
				throw std::out_of_range("sparseCSR3Matrix<VAL,IDX>::put_():6");
		}
	}
}


// Somma una intera sottomatrice. Le componenti K[i][j] vengono sommate alle componenti (RI[i],CI[j]) a meno che 
// RI[i] o CI[j] siano zero.
// Se symut==true, tratta sia K che la matrice sparsa come simmetriche triangolari superiori
// e aggiunge solo le componenti del triangolo superiore di K nel triangolo superiore. 
// Assume che K punti alle componenti di una matrice nr x nc orientata alle righe.  
// Tutte le  componenti devono gi� esistere.

template<typename	VAL, typename	IDX, typename MP >
template<typename MAT>
void sparseCSR3Matrix< VAL, IDX, MP>::sum_(const iSet& ir, const iSet& ic, const MAT& K)
{
	if (index_t(ir.max()) > nrows() || index_t(ic.max())>ncols())
		throw std::out_of_range("sparseCSR3Matrix<VAL,IDX>::sum():6");
	index_t rsz = (index_t)ir.size(),
		csz = (index_t)ic.size();
	for (index_t r = 1; r <= rsz; ++r)
	{
		index_t R = ir(r);
		if (!R)
			continue;

		for (index_t c = 1; c <= csz; ++c)
		{
			index_t C = ic(c);
			if (!C || (IS::SYM && C<R))
				continue;
			value_t * v = get_p(R, C);
			if (v)
				*v += K(r, c);
			else
				throw std::out_of_range("sparseCSR3Matrix<VAL,IDX>::sum():6");
		}
	}
}

	


template<typename	VAL, typename	IDX, typename MP >
template<typename MAT>
void	sparseCSR3Matrix<VAL, IDX, MP>::sum_(const iSet& irc, const MAT& K)
{
	if (index_t(irc.max())> nrows())
		throw std::out_of_range("sparseCSR3Matrix::gen_sum_");
	for (index_t r = 1, iend = irc.size() + 1; r != iend; ++r)
	{
		index_t R = irc(r);
		if (R == IDX(0))
			continue;
		for (index_t c = IS::SYM?r:1; c != iend; ++c)
		{
			index_t C = irc(c);
			if (C == IDX(0))
				continue;
			value_t * v = (!IS::SYM || C >= R  ) ? get_p(R, C) : get_p(C, R);
			if (v)
				*v += K(r, c);
			else
				throw std::out_of_range("sparseCSR3Matrix<VAL,IDX>::sum_():7");
		}
	}
}








template< typename VAL, typename	IDX, typename MP>
int	sparseCSR3Matrix<VAL,IDX, MP>::fwrite(const std::string& fname)const
{
	std::ofstream os(fname, std::ios_base::binary | std::ios_base::trunc);
	if (!os)
		return -1;
	// Scrive un header con le dimensioni ed il tipo dei dati e della matrice
	std::uint64_t header_[8] = 
	{ 
		(std::uint64_t)nrows_,
		(std::uint64_t)ncols_,
		(std::uint64_t)nonzeros(),
		(std::uint64_t)(is_symmetric() ? 1 : 0),
		(std::uint64_t)(is_row_major() ? 1 : 0),
		(std::uint64_t)(IS::UT ? 1 : 0),
		(std::uint64_t)type<index_t>::id(),
		(std::uint64_t)type<value_t>::id()
	};
	os.write(reinterpret_cast<const char*>(header_), sizeof(header_));
	if (!os)
		return -2;
	// scrive il vettore posizione delle righe
	os.write(reinterpret_cast<const char*>(row_pos_.data()), sizeof(index_t)* row_pos_.size());
	if (!os)
		return -2;
	// scrive il vettore degli indici di colonna 
	os.write(reinterpret_cast<const char*>(col_idx_.data()), sizeof(index_t)* col_idx_.size());
	if (!os)
		return -2;
	// scrive il vettore dei valori
	os.write(reinterpret_cast<const char*>(col_val_.data()), sizeof(value_t)* col_val_.size());
	if (!os)
		return -2;
	return 0;
}




template< typename VAL, typename	IDX, typename MP>
int	sparseCSR3Matrix<VAL,IDX, MP>::fread(const std::string& fname)
{
	std::ifstream is(fname, std::ios_base::binary);
	if (!is.is_open())
		return -1;

	std::uint64_t this_header_[8] =
	{
		(std::uint64_t)0,
		(std::uint64_t)0,
		(std::uint64_t)0,
		(std::uint64_t)(is_symmetric() ? 1 : 0),
		(std::uint64_t)(is_row_major() ? 1 : 0),
		(std::uint64_t)(IS::UT ? 1 : 0),
		(std::uint64_t)type<index_t>::id(),
		(std::uint64_t)type<value_t>::id()
	};

	std::uint64_t onfile_header_[8];

	is.read(reinterpret_cast< char*>(onfile_header_), sizeof(onfile_header_));
	if (!is)
		return -2;
	for (int i = 3; i != 8; ++i)
	if (this_header_[i] != onfile_header_[i])
		return -3;

	index_t nr = (index_t)onfile_header_[0];
	index_t nc = (index_t)onfile_header_[1];
	size_t nnz = (size_t)onfile_header_[2];
	if (nnz > size_t(std::numeric_limits<index_t>::max()))
		return -4;

	resize(nr, nc, nnz);

	is.read(reinterpret_cast< char*>(row_pos_.data()), sizeof(index_t)* row_pos_.size());
	if (!is)
		return -2;

	is.read(reinterpret_cast< char*>(col_idx_.data()), sizeof(index_t)* col_idx_.size());
	if (!is)
		return -2;

	is.read(reinterpret_cast< char*>(col_val_.data()), sizeof(value_t)* col_val_.size());
	if (!is)
		return -2;
	return 0;
}

	


	


template<typename VAL, typename	IDX, typename MP>
void	sparseCSR3Matrix<VAL,IDX,MP>::print(std::ostream& os)const
	{

		std::streamsize width=os.width(0);
		int count=1;
		int idx_width=std::max( nrows(),ncols())/10;
		for(;idx_width; idx_width/=10)
			++count;
		idx_width+=1;
		os	<< (unsigned int) nrows() << " * "  << (unsigned int) ncols() << "\tnnz="	<< (unsigned int) nonzeros() << std::endl;
		for(index_t i=1; i<=nrows();++i) 
		{
			
			index_t j1=row_pos_[i-1]-1;
			index_t j2=row_pos_[i]-1;
			for (index_t j = j1; j != j2; ++j)
				os << "(" << std::setw(idx_width) << i << "," << std::setw(idx_width) << col_idx_[j] << ")= " << std::setw(width) << col_val_[j] << std::endl;
			
		}
	}



template<typename VAL, typename	IDX, typename MP>
void	sparseCSR3Matrix< VAL, IDX, MP>::product(const value_t* vecB, index_t szB, value_t* vecC, index_t szC)const
{
	index_t NR = nrows(), NC = ncols();
	if (NC != szB || NR!= szC)
		throw std::runtime_error("sparse product");
	if (IS::GE)
	{
		for (index_t r = 0; r != NR; ++r)
		{
			index_t j1 = row_pos_[r] - 1;
			index_t j2 = row_pos_[r + 1] - 1;
			value_t tmp = 0.0;
			for (size_t j = j1; j != j2; ++j)
				tmp += col_val_[j] * vecB[col_idx_[j] - 1];
			vecC[r] = tmp;
		}
	}
	else if (IS::SYM)
	{
		value_t val;
		index_t J;
		std::fill(vecC, vecC + szC, value_t(0));
		for (index_t r = 0; r != NR; ++r)
		{
			index_t j1 = row_pos_[r] - 1;
			index_t j2 = row_pos_[r + 1] - 1;
			value_t tmp = 0.0;

			if (j1 == j2)
				continue;
			if (col_idx_[j1] == r + 1)
			{
				vecC[r] += col_val_[j1] * vecB[col_idx_[j1] - 1];
				j1++;
			}
			for (size_t j = j1; j != j2; ++j)
			{
				val = col_val_[j];
				J = col_idx_[j]-1;
				tmp += val * vecB[J];
				vecC[J] += val*vecB[r];
			}
			vecC[r]+= tmp;
		}

	}
	else
		throw std::runtime_error("matrix type product not supported");
}


	/*
template<int MT, typename VAL, typename	IDX>
void product(const static_sparse_CR_Matrix<MT, VAL,IDX>& MAT, const Vector<VAL>& U, Vector<VAL>& V, bool transpA=false)
{
	if(!transpA)
	{
		if(MAT.ncols() != U.size())
			throw std::runtime_error("product(const static_sparse_CR_Matrix<MT, VAL,IDX>& MAT, const Vector<VAL>& U, Vector<VAL>& V): size error");
		if(MAT.nrows() != V.size())
			V.resize(MAT.nrows());
	}
	else
	{
		if(MAT.nrows() != U.size())
			throw std::runtime_error("product(const static_sparse_CR_Matrix<MT, VAL,IDX>& MAT, const Vector<VAL>& U, Vector<VAL>& V): size error-2");
		if(MAT.ncols() != V.size())
			V.resize(MAT.ncols());

	}

	MAT.product(U.begin(),U.size(),V.begin(), transpA);

}
*/


} // end namespace MML
#endif