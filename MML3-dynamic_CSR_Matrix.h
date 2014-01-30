#pragma once
#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

#include<stdexcept> 
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<iomanip>
#include<vector>
#include<stdexcept>
#include<limits>
#include<algorithm>
#include<random>

#include"MML3-matrix.h"
#include"MML3-sparseMatrix.h"
#include"MML3-memory-block-allocator.h"
#include"MML3-ordered-single-linked-list.h"




namespace MML3{

/// \addtogroup MATRICES
// @{

// Dynamic Sparse Compressed Row  Matrix (index base 1)
// "Dynamic" refers to the sparsity structure : it is always possible to add components
	
	template<	typename	VALUE_TYPE,				// the matrix component type
				typename	IDX_TYPE = int_t,		// the matrix index type
				class       MP = M_PROP::GE,		// GE for general, SYM for symmetric matrices
				class		ALLOCATOR = MML3::MemoryBlockAllocator< MML3::sparse_component<VALUE_TYPE, IDX_TYPE> >>
	class dynamic_sparse_CSR_Matrix;

	template<typename T>
	using 	dynamic_sparse_matrix		= dynamic_sparse_CSR_Matrix<T,int_t, M_PROP::GE>;
			
	template<typename T>
	using 	dynamic_sparse_sym_matrix	= dynamic_sparse_CSR_Matrix<T,int_t, M_PROP::SYM>;
	
// WARNING  !!!
// be careful in changing the index_type (IDX_TYPE). Sparse solver algorithms generally require the type (int) that is int32 or int64 
// depending on the machine

	

/** Class dynamic_sparse_CSR_Matrix<VAL, IDX, ALLOCATOR> : Matrici sparse con allocazione dinamica ottimizzata organizzate per righe.
 Gestisce  matrici sparse con componenti di tipo VAL e indici di tipo IDX mediante un array
 di liste concatenate semplici ordinate (ogni riga è una lista). E' dotata di un allocatore di memoria di default ottimizzato che gestisce 
 la memoria in blocchi di allocazione. Il blocco di allocazione  è il numero di elementi (di tipo VAL) che vengono allocati  ogni qualvolta 
 è esaurita la memoria riservata ed è  richiesto l'inserimento di un nuovo elemento.
 Gli indici partono da 1.
 */
template<typename VALUE_TYPE,typename IDX_TYPE,class MP,class ALLOCATOR>
class dynamic_sparse_CSR_Matrix : public sparseMatrix<VALUE_TYPE,IDX_TYPE>
{
	typedef sparseMatrix<VALUE_TYPE, IDX_TYPE>										base_matrix_t;
	typedef  MML3::ordered_single_linked_list<VALUE_TYPE, IDX_TYPE, ALLOCATOR>		OSL_list;
	
public:
	

	typedef typename OSL_list::iterator					col_iterator;	///< il tipo dell'iteratore sugli elementi delle righe
	typedef typename OSL_list::const_iterator			const_col_iterator;
	
	
	enum IS :bool {
		RowMajor = true,
		ColMajor = false,
		SYM = std::is_same<M_PROP::SYM, MP>::value,
		RE = false,
		LT = false,
		UT =SYM // if it is symmetric it is also Upper triangular
	};

	//-------------------------------------------------------------------------
	// COMMON INTERFACE
	//-------------------------------------------------------------------------
	typedef VALUE_TYPE									value_t;
	typedef IDX_TYPE									index_t;
	typedef iSet										iset_t;
	typedef base_matrix_t::inserter_matrix_t			inserter_matrix_t;
	typedef base_matrix_t::inserter_sym_matrix_t		inserter_sym_matrix_t;
	

	dynamic_sparse_CSR_Matrix()=default;
	dynamic_sparse_CSR_Matrix(const dynamic_sparse_CSR_Matrix& );
	dynamic_sparse_CSR_Matrix(dynamic_sparse_CSR_Matrix&& );
	~dynamic_sparse_CSR_Matrix()																{ destroy(); }
	void				destroy();
	index_t				nrows()const														{ return row_.size(); }
	index_t				ncols()const														{ return ncols_; }
	bool				is_row_major()const													{ return IS::RowMajor; }
	bool				is_symmetric()const													{ return IS::SYM; }
	size_t				nonzeros()const;
	bool				test_subscripts(index_t r, index_t c)const;
	void				fill(const value_t&);

	// get the pointer to the component (i,j) if it exists, nullptr elsewhere
	const	value_t*	get_p(index_t i, index_t j)const;
			value_t*	get_p(index_t i, index_t j);
	
	value_t&			put(index_t i, index_t j, const value_t& val);
	value_t&			sum(index_t i, index_t j, const value_t& val);

	void				put(const iSet& ir, const iSet& ic, const inserter_matrix_t& K);

	void				sum(const iSet& ir, const iSet& ic, const inserter_matrix_t& K);
	void				sum(const iSet& irc, const inserter_matrix_t& K);
	void				sum(const iSet& irc, const inserter_sym_matrix_t& K);


	//int					copy2(index_t rpos_sz, index_t nnz, index_t rpos[], index_t cidx[], value_t  val[])const;
	int					fwrite(const std::string& fname)const;
	int					fread(const std::string&);
	void				print(std::ostream& os)const;

	
	//-------------------------------------------------------------------------
	// EXTENDED INTERFACE
	//-------------------------------------------------------------------------	
	
	///sparse matrix  nr X nc with allocation blocks of bsz elements (every time memory is exausted a memory block of bsz elements is created) 
	//@param nr number of rows
	//@param nc number of cols
	//@param bsz size of the allocation block in number of componentsstored in the block. if bsz is left unspecified than it is assumed bsz= max(nr,nc)+1  
						dynamic_sparse_CSR_Matrix(index_t nr, index_t nc, size_t bsz = 0);
	///Ridimensiona la matrice sparsa, deallocando prima tutto il suo contenuto precedente.
	///@param R numero di righe
	///@param C numero di colonne
	///@dimensione (in numero di componenti) del blocco di allocazione. Se bsz non viene specificato allora il blocco di allocazione è di max(C,R)+1 elementi 
	void				resize(index_t R, index_t C, size_t bsz=0);
	/// Scambio
	void				swap(dynamic_sparse_CSR_Matrix& rhs);
	/// Ritorna un iteratore che punta al primo elemento della lista i-esima
	col_iterator		row_begin(index_t i){return row_[i-1]->begin();}
	const_col_iterator	row_begin(index_t i)const{return row_[i-1]->begin();}
	/// Ritorna un iteratore che punta all'elemento dopo l'ultimo
	const_col_iterator	row_end(index_t i)const{ return row_[i - 1]->end(); }
	col_iterator		row_end(index_t i){ return row_[i - 1]->end(); }
	
	///Informazioni sulla struttura della matrice
	///@param NCR in uscita ha dimensione pari al numero di righe della matrice sparsa e contiene  il numero di elementi diversi da zero di ogni riga
	///@return il numero totale di elementi 
	size_t				nonzeros(std::vector<index_t>& nr_nnz)const;
	void				get_min_max_diag(index_t& row_idx_min, value_t& val_min , index_t& row_idx_max, value_t& val_max)const;


	static dynamic_sparse_CSR_Matrix  get_random_matrix(index_t nr, index_t nc, value_t sparsity, value_t mean, value_t range);
private:

	void				sym_put_(const iSet& ir, const iSet& ic, const inserter_matrix_t& K);
	void				gen_put_(const iSet& ir, const iSet& ic, const inserter_matrix_t& K);

	void				sym_sum_(const iSet& ir, const iSet& ic, const inserter_matrix_t& K);
	void				gen_sum_(const iSet& ir, const iSet& ic, const inserter_matrix_t& K);
	

	template<typename MAT>
	void				sym_sum_(const iSet& icr, const MAT& K);
	template<typename MAT>
	void				gen_sum_(const iSet& irc, const MAT& K);
	
	index_t					ncols_=0;
	std::vector<OSL_list*>	row_;
	ALLOCATOR*				allocator_=nullptr;
};






////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// IMPLEMENTAZIONE 



template<typename T>
inline T get_max_val(const T* p, size_t sz)
{
	if (!sz)
		return 0;
	T max = *p++;
	for (int i = 1; i != sz; ++i)
		if (*p > max)
			max = *p;
	return max;
}







template<typename	VAL,typename	IDX, class MP, class		ALLOCATOR>
dynamic_sparse_CSR_Matrix<VAL, IDX, MP,ALLOCATOR>::dynamic_sparse_CSR_Matrix(index_t R, index_t C, size_t bsz):
	ncols_(C), 
	row_(R), 
	allocator_( new ALLOCATOR(bsz == 0 ? (std::max(R, C) + 1): bsz))
{
	for (auto& r:row_)
		r = new OSL_list(allocator_);
}

template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::dynamic_sparse_CSR_Matrix(const dynamic_sparse_CSR_Matrix& rhs)
:ncols_(rhs.ncols_), row_(rhs.nrows()), allocator_(new ALLOCATOR(rhs.allocator_->block_size()))
	{
		size_t nr = nrows();
		for (size_t r = 0; r != nr; ++r)
			row_[r] = new OSL_list(allocator_, *rhs.row_[r]);
	}


template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::dynamic_sparse_CSR_Matrix(dynamic_sparse_CSR_Matrix&& rhs)
:ncols_(rhs.ncols_), allocator_(rhs.allocator_), row_(std::move(rhs.row_))
{
	rhs.allocator_=nullptr;
	rhs.ncols_ = 0;
}


template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
inline bool dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::test_subscripts(index_t r, index_t c)const
{
	if ((r > nrows() || r<1) || (c > ncols() || c<1))
		return false;
	return true;
}


		/// Ridimensiona il vettore di liste deallocando prima tutto il suo contenuto precedente
template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
void dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::resize(index_t R, index_t C, size_t bsz)
	{
		destroy();
		dynamic_sparse_CSR_Matrix	tmp(R,C,bsz);
		swap(tmp);
	}



	/// Distrugge la matrice deallocando le risorse
template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
void dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::destroy()
{
		for(auto& r:row_) 
			delete r;
		row_.clear();
		delete allocator_;
		allocator_=nullptr;
		ncols_ = 0;
}

	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	void dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::swap(dynamic_sparse_CSR_Matrix& rhs)
	{
		std::swap(ncols_,rhs.ncols_);
		std::swap(row_,rhs.row_);
		std::swap(allocator_,rhs.allocator_);
	}


	

	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	inline auto  dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::get_p(index_t i, index_t j)const ->const value_t*
	{
		if (i <1 || i> nrows())
			throw std::out_of_range("dynamic_sparse_CSR_Matrix<>::get_p");
		if (IS::SYM && j < i)
			std::swap(i, j);
		auto e = row_[i - 1]->get(j);
		return e ? &(e->value) : nullptr;
	}

	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	inline auto  dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::get_p(index_t i, index_t j)-> value_t*
	{
		if (i <1 || i> nrows())
		throw std::out_of_range("dynamic_sparse_CSR_Matrix<>::get_p");
		if (IS::SYM && j < i)
			std::swap(i, j);
		auto e = row_[i - 1]->get(j);
		return e ? &(e->value) : nullptr;
	}

	
	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	inline void dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::put(const iSet& ir, const iSet& ic, const inserter_matrix_t& K)
	{
		if (IS::SYM)
			sym_put_(ir, ic, K);
		else
			gen_put_(ir,  ic,  K);
	}


	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	inline void dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::gen_put_(const iSet& ir, const iSet& ic, const inserter_matrix_t& K)
	{
		size_t NR = nrows(), NC = ncols();
		if (index_t(ir.max()) > nrows() || index_t(ic.max()) > ncols())
			throw std::out_of_range("dynamic_sparse_CSR_Matrix::gen_put_");
		size_t R, C;
		for (size_t r = 1, rend=ir.size()+1;	r != rend; ++r)
		{
			R = ir(r);
			if (R == IDX(0))
				continue;
			for (size_t c = 1, cend=ic.size()+1;	c != cend; ++c)
			{
				C = ic(c);
				if (C == IDX(0))
					continue;
				row_[R - 1]->insert(C, K(r, c));
			}
		}
	}
	


	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	inline void dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::sym_put_(const iSet& ir, const iSet& ic, const inserter_matrix_t& K)
	{
		size_t NR = nrows(), NC = ncols(), R, C;
		if (index_t(ir.max())>nrows() || index_t(ic.max())> ncols())
			throw std::out_of_range("dynamic_sparse_CSR_Matrix::sym_put_");


		for (size_t r = 1, rend=ir.size()+1; r != rend; ++r)
		{
			R = ir(r);
			if (R == IDX(0))
				continue;
			for (size_t c = 1, cend=ic.size(); c != cend; ++c)
			{
				C = ic(c);
				if (C == IDX(0) || C < R)
					continue;
				row_[R - 1]->insert(C, K(r, c));
			}
		}
	}



	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	inline void dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::sum(const iSet& ir, const iSet& ic, const inserter_matrix_t& K)
	{
		if (IS::SYM)
			sym_sum_(ir, ic, K);
		else
			gen_sum_(ir, ic, K);
	}

	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	inline void dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::sum(const iSet& irc, const inserter_matrix_t& K)
	{
		if (IS::SYM)
			sym_sum_(irc, K);
		else
			gen_sum_(irc, K);
	}

	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	inline void dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::sum(const iSet& irc, const inserter_sym_matrix_t& K)
	{
		if (IS::SYM)
			sym_sum_(irc, K);
		else
			gen_sum_(irc, K);
	}

	
	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	inline void dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::gen_sum_(const iSet& ir, const iSet& ic, const inserter_matrix_t& K)
	{
		
		if (index_t(ir.max()) > nrows() || index_t(ic.max())>ncols())
			throw std::out_of_range("dynamic_sparse_CSR_Matrix::gen_sum_");
		size_t R, C;
		for (size_t r = 1, rend=ir.size()+1; r != rend; ++r)
		{
			R = ir(r);
			if (R == IDX(0))
				continue;
			for (size_t c = 1, cend=ic.size()+1; c != cend; ++c)
			{
				C = ic(c);
				if (C == IDX(0))
					continue;
				row_[R - 1]->add_at(C, K(r, c));
			}
		}

	}

	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	template<typename MAT>
	inline void dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::gen_sum_(const iSet& irc,  const MAT& K)
	{
		index_t max_idx = irc.max();
		if (max_idx> nrows() || max_idx> ncols())
			throw std::out_of_range("dynamic_sparse_CSR_Matrix::gen_sum_");
		
		for (index_t r = 1, end = irc.size() + 1; r != end; ++r)
		{
			index_t R = irc(r);
			if (R == IDX(0))
				continue;
			for (index_t c = 1; c != end; ++c)
			{
				index_t C = irc(c);
				if (C == IDX(0))
					continue;
				row_[R - 1]->add_at(C, K(r, c));
			}
		}
	}


	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	inline void dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::sym_sum_(const iSet& ir, const iSet& ic, const inserter_matrix_t& K)
	{
		
		index_t  R, C;
		if (index_t(ir.max()) > nrows() || index_t(ic.max()) > ncols())
			throw std::out_of_range("dynamic_sparse_CSR_Matrix::sym_sum_");

		for (index_t r = 1, rend = (index_t)ir.size()+1; r != rend; ++r)
		{
			R = ir(r);
			if (R == IDX(0))
				continue;
			for (index_t c = 1, cend = (index_t)ic.size()+1; c != cend; ++c)
			{
				C = ic(c);
				if (C == IDX(0) || C < R)
					continue;
				row_[R - 1]->add_at(C, K(r, c));
			}
		}
	}

	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	template<typename MAT>
	inline void dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::sym_sum_(const iSet& irc, const MAT& K)
	{
		if (index_t(irc.max())> nrows())
			throw std::out_of_range("dynamic_sparse_CSR_Matrix::sym_sum_");
		for (index_t r = 1, iend = irc.size() + 1; r != iend; ++r)
		{
			index_t R=irc(r);
			if (R == IDX(0))
				continue;
			for (index_t c = r; c !=iend; ++c)
			{
				index_t C = irc(c);
				if (C == IDX(0))
					continue;
				if ( C >= R)
					row_[R - 1]->add_at(C, K(r, c));
				else
					row_[C - 1]->add_at(R, K(r, c));
			}
		}
	}




	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	void	dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::fill(const value_t& val)
	{
		for (size_t i = 0; i != nrows(); ++i)
			row_[i]->fill(val);
	}




	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	void	dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::get_min_max_diag(index_t& idx_min, value_t& val_min, index_t& idx_max, value_t& val_max)const
	{
		idx_min=0;
		idx_max=0;
		val_min=std::numeric_limits<VAL>::max();
		val_max=std::numeric_limits<VAL>::min();

		VAL value;
		bool found;
		for(size_t i=0; i!=nrows(); ++i)
		{
			row_[i]->get(i+1, value, found);
			if(!found)
				value=VAL(0);
			if(value > val_max)
			{
				val_max=value;
				idx_max=i+1;
			}
			if(value < val_min)
			{
				val_min=value;
				idx_min=i+1;
			}

		}
		
	}

	/** Inserisce nella componente (i,j) il valore val; se la componente non c'è la crea.
	    ritorna il puntatore alla componente inserita  in caso di successo e 
		lancia un'eccezione  std::out_of_range nel caso gli indici siano fuori campo, .
	 */
	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	inline VAL& dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR > ::put(index_t i, index_t j, const value_t& val)
	{
		if(!  test_subscripts(i,j))
			throw std::out_of_range("DSC_Matrix<VAL,IDX,ALLOCATOR>::put");

		if (IS::SYM && j < i)
			std::swap(i, j);
		return *(row_[i-1]->insert(j,val));
	}

	
	template<typename	VAL, typename	IDX, class MP, class	ALLOCATOR>
	inline VAL& dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::sum(index_t i, index_t j, const value_t& val)
	{
		if(!  test_subscripts(i,j))
			throw std::out_of_range("");
		if (IS::SYM && j < i)
			std::swap(i, j);
		return	*(row_[i-1]->add_at(j,val));
	}
	
//}



	
	// Informazioni sulla struttura della matrice
	// In uscita NCR ha dimensione nrows e contiene  il numero di elementi diversi da zero di ogni riga 
	// Ritorna il numero totale di elementi 
	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	inline size_t	dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::nonzeros(std::vector<index_t>& NCR)const
	{
		NCR.clear();
		NCR.reserve(nrows());
		size_t s, counter=0;
		for (auto r : row_)
		{
			 NCR.push_back( s=r->size() );
			 counter += s;
		}
		return counter;
	}


	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	inline size_t	dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::nonzeros()const
	{
		size_t counter=0;
		for(auto cp:row_)
			 counter +=cp->size();
		return counter;
	}



	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	void	dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::print(std::ostream& os)const
	{
		std::streamsize width=os.width(0);
		int count=1;
		int idx_width=std::max( nrows(),ncols())/10;
		for(;idx_width; idx_width/=10)
			++count;
		idx_width=count +1;

		os	<< (unsigned int) nrows() << " * "  << (unsigned int) ncols() << "\tnnz="	<< (unsigned int) nonzeros() << std::endl;
		
		for(index_t r=0; r!=nrows();++r)
		{	
			for( auto it=row_begin(r+1); it.is_valid(); ++it)
				os << "(" << std::setw(idx_width) << r + 1 << "," << std::setw(idx_width) << it->idx << ")= " << std::setw(width) << it->value << std::endl;
		}
	}

	
			

		
		///@param fname nome del file
		///@return:  0 in caso di successo
		///@return: -1 nel caso il file non possa essere aperto
		///@return: -2 se si sono verificati errori in scrittura
	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	int	dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::fwrite(const std::string& fname)const
	{
		std::ofstream os(fname, std::ios_base::binary | std::ios_base::trunc);

		if(!os)
			return -1; 
		std::vector<index_t>	RC;
		size_t	nnz = nonzeros(RC);

		std::vector<int> row_pos(nrows() + 1);
		row_pos[0] = 1;
		for (size_t i = 0; i != RC.size(); ++i)
			row_pos[i + 1] = row_pos[i] + RC[i];

		if ((row_pos[nrows()]-1) != nnz)
			throw std::runtime_error("dynamic_sparse_CSR_Matrix::fwrite");

		// Scrive un header con le dimensioni ed il tipo dei dati e della matrice
		std::uint64_t header_[8] = {nrows(), 
									ncols(), 
									nnz, 
									is_symmetric()?1:0,
									is_row_major()?1:0,
									IS::UT?1:0,
									std::uint64_t(type<index_t>::id()), 
									type<value_t>::id()
		};

		os.write(reinterpret_cast<const char*>(header_), sizeof(header_));
		if (!os)
			return -2;
		// scrive il vettore posizione delle righe
		os.write(reinterpret_cast<const char*>(&row_pos[0]), sizeof(index_t)* row_pos.size());
		if (!os)
			return -2;
		// scrive tutti gli nnz indici di colonna
		std::vector<IDX>	i_buffer;
		for( index_t r=0; r!=nrows();++r)
		{
			row_[r]->get_indexes(i_buffer);
			if(i_buffer.size())
			{
				os.write( reinterpret_cast<const char*>(&i_buffer[0]), sizeof(index_t)* i_buffer.size());
				if (!os.good())
					return -1;
			}
		}

		// scrive tutte le nnz componenti
		std::vector<value_t>	v_buffer;
		// scrivo i valori 
		for( index_t r=0; r!=nrows();++r)
		{
			row_[r]->get_values(v_buffer);
			if(v_buffer.size())
			{
				os.write(reinterpret_cast<const char*>(&v_buffer[0]), sizeof(value_t)* v_buffer.size());
				if (!os)
					return -2;
			}
		}
		return 0;
	}
	
	/* int read(const char* fname, index_t bsize=0);
	Acquisisce da file una matrice sparsa 
	*/
		///@param fname nome del file
		///@param bsize
		///@return:  0 in caso di successo
		///@return: -1 nel caso il file non possa essere aperto
		///@return: -2 se si sono verificati errori in scrittura
		///@return: -3 se si sono verificati errori di incompatibilita' dei tipi
		
	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	int	dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::fread(const std::string& fname)
	{

		std::ifstream is(fname,std::ios_base::binary);
		if(!is.is_open())
			return -1; 



		std::uint64_t this_header_[8] = { 
			nrows(),
			ncols(),
			nonzeros(),
			is_symmetric() ? 1 : 0,
			is_row_major() ? 1 : 0,
			IS::UT ? 1 : 0,
			std::uint64_t(type<index_t>::id()),
			type<value_t>::id()
		};


		std::uint64_t onfile_header_[8];

		is.read(reinterpret_cast< char*>(onfile_header_), sizeof(onfile_header_));
		if(!is)
				return -2;
		for (int i = 3; i != 8;++i)
		if (this_header_[i] != onfile_header_[i])
			return -3;

		index_t nr = (index_t)onfile_header_[0];
		index_t nc = (index_t)onfile_header_[1];
		index_t nnz = (index_t)onfile_header_[2];

		std::vector<index_t>	row_pos(nr+1);
		is.read(reinterpret_cast< char*>(&row_pos[0]), sizeof(index_t)* row_pos.size());
		if(!is)
				return -2;
		index_t bsize = std::max(std::max(nr, nc), nnz / 4);

		resize(nr,nc, bsize);

		std::vector<index_t>	i_buff(nnz);
		std::vector<value_t>	v_buff(nr);
		is.read(reinterpret_cast< char*>(&i_buff[0]), sizeof(index_t)*nnz);
		if(!is)
				return -2;
		index_t* pi=&(i_buff[0]);
		for(index_t r=0; r!=nr; ++r)
		{
			index_t r_sz = row_pos[r + 1] - row_pos[r];
			if(index_t(v_buff.size()) < r_sz )
				v_buff.resize(r_sz);
			is.read(reinterpret_cast< char*>(&v_buff[0]), sizeof(value_t)*r_sz );
			if(!is )
				return -2;

			row_[r]->set(pi, &v_buff[0], r_sz);
			pi+=r_sz;
		}
		return 0;
	}



	
	template<typename	VAL, typename	IDX, class MP, class		ALLOCATOR>
	
	
	auto dynamic_sparse_CSR_Matrix<VAL, IDX, MP, ALLOCATOR>::get_random_matrix(index_t nr, index_t nc, value_t sparsity, value_t mean, value_t range)->dynamic_sparse_CSR_Matrix
	{
		dynamic_sparse_CSR_Matrix A(nr, nc);
		// Seed with a real random value, if available
		std::random_device rd;
		std::default_random_engine e1(rd());
		std::uniform_int_distribution<int> uniform_dist_cidx(1, nc);
		std::uniform_real_distribution<double> uniform_dist_val(mean-range/2, mean+range/2);
		index_t ncs = size_t(nc*sparsity) + 1;
		for (index_t R = 1; R <= nr; ++R) 
		{
			for (index_t j = 0; j != ncs; ++j)
			{
				index_t C = index_t(uniform_dist_cidx(e1));
				A.sum(R, C, uniform_dist_val(e1));
			}
			// inserts always the diagonal element
			if ( R <=nc)
				A.put(R, R, uniform_dist_val(e1));
		}
		return A;
	}


	
} // end namespace MML

