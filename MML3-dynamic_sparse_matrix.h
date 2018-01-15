#ifndef _MML3_DYNAMIC_SPARSE_MATRIX_H_
#define _MML3_DYNAMIC_SPARSE_MATRIX_H_


#include<stdexcept> 
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<iomanip>
#include<vector>
#include<map>
#include<stdexcept>
#include<limits>
#include<algorithm>
#include<random>

#include"MML3-matrix.h"
#include"MML3-basic_sparse_matrix.h"
#include"MML3-ordered_list_map.h"
#include"MML3-ordered_list_slist.h"



namespace MML3{

/// \addtogroup MATRICES
// @{
/// Dynamic Sparse (Row ordered)  Matrix (index base 1)
/// "Dynamic" refers to the sparsity structure : it is always possible to add components
/// Template parameter   VAL	: the type of components
/// Template parameter   IDX	: the type of the idexes
/// Template parameter   MP		: M_PROP::GE for GEneral matrices
///								  M_PROP::SYM for SYMmetric (only the upper triangular part is stored)
/// Template parameter   ROW_T	: the structure used to manage the elements of each row
///								  see the header "MML3-ordered_list_map.h" or "MML3-ordered_list_slist.h"
	
	template<	typename	VAL,					
				typename	IDX = int_t,			
				class       MP = M_PROP::GE,		
				class       ROW_T=ordered_list_map<IDX,VAL>	>
	class dynamic_sparse_matrix;

			
	template<typename T, typename IDX_TYPE = int_t>
	using 	sparse_matrix = dynamic_sparse_matrix<T, IDX_TYPE, M_PROP::SYM>;

	template<typename T, typename IDX_TYPE = int_t>
	using 	sparse_sym_matrix = dynamic_sparse_matrix<T, IDX_TYPE, M_PROP::SYM>;
	


template<typename VAL,typename IDX,class MP, class ROW_T>
class dynamic_sparse_matrix : public basic_sparse_matrix<VAL,IDX>
{
	typedef basic_sparse_matrix<VAL, IDX>		base_matrix_t;
	typedef  ROW_T						row_t;
		
public:
	
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
	typedef VAL										value_t;
	typedef IDX										index_t;
	typedef iSet									iset_t;
	typedef base_matrix_t::inserter_matrix_t		inserter_matrix_t;
	typedef base_matrix_t::inserter_sym_matrix_t	inserter_sym_matrix_t;
	

	dynamic_sparse_matrix()=default;
	dynamic_sparse_matrix(const dynamic_sparse_matrix& );
	dynamic_sparse_matrix(dynamic_sparse_matrix&& );
	~dynamic_sparse_matrix()																{ destroy(); }
	void				destroy();
	index_t				nrows()const														{ return row_.size(); }
	index_t				ncols()const														{ return ncols_; }
	bool				is_row_major()const													{ return IS::RowMajor; }
	bool				is_symmetric()const													{ return IS::SYM; }
	size_t				nonzeros()const;
	bool				test_subscripts(index_t r, index_t c)const;
	void				fill(const value_t&);

	/// get the pointer to the component (i,j) if it exists, nullptr elsewhere
	const	value_t*	get_p(index_t i, index_t j)const;
			value_t*	get_p(index_t i, index_t j);
	
	/// inserts the value val at position (i,j); if the component does not exists it is created
	value_t&			put(index_t i, index_t j, const value_t& val);
	/// sum value val at the compoent (i,j); if the component does not exists it is created
	value_t&			sum(index_t i, index_t j, const value_t& val);

	void				put(const iSet& ir, const iSet& ic, const inserter_matrix_t& K);
	void				sum(const iSet& ir, const iSet& ic, const inserter_matrix_t& K);
	void				sum(const iSet& irc, const inserter_matrix_t& K);
	void				sum(const iSet& irc, const inserter_sym_matrix_t& K);



	//int				copy2(index_t rpos_sz, index_t nnz, index_t rpos[], index_t cidx[], value_t  val[])const;
	int					fwrite(const std::string& fname)const;
	int					fread(const std::string&);
	void				print(std::ostream& os)const;

	
	//-------------------------------------------------------------------------
	// EXTENDED INTERFACE
	//-------------------------------------------------------------------------	
	
	typedef typename row_t::iterator		col_iterator;
	typedef typename row_t::const_iterator	const_col_iterator;

	///sparse matrix  nr X nc with allocation blocks of bsz elements (every time memory is exausted a memory block of bsz elements is created) 
	//@param nr number of rows
	//@param nc number of cols
	//@param bsz size of the allocation block in number of componentsstored in the block. if bsz is left unspecified than it is assumed bsz= max(nr,nc)+1  
						dynamic_sparse_matrix(index_t nr, index_t nc, size_t bsz = 0);
	///Ridimensiona la matrice sparsa, deallocando prima tutto il suo contenuto precedente.
	///@param R numero di righe
	///@param C numero di colonne
	///@dimensione (in numero di componenti) del blocco di allocazione. Se bsz non viene specificato allora il blocco di allocazione è di max(C,R)+1 elementi 
	void				resize(index_t R, index_t C, size_t bsz=0);
	/// Scambio
	void				swap(dynamic_sparse_matrix& rhs);
	/// Ritorna un iteratore che punta al primo elemento della lista i-esima
	col_iterator		row_begin(index_t i){return row_[i-1]->begin();}
	const_col_iterator	row_begin(index_t i)const{return row_[i-1]->begin();}
	/// Ritorna un iteratore che punta all'elemento dopo l'ultimo
	const_col_iterator	row_end(index_t i)const{ return row_[i - 1]->end(); }
	
	
	///Informazioni sulla struttura della matrice
	///@param NCR in uscita ha dimensione pari al numero di righe della matrice sparsa e contiene  il numero di elementi diversi da zero di ogni riga
	///@return il numero totale di elementi 
	size_t				nonzeros(std::vector<index_t>& nr_nnz)const;
	void				get_min_max_diag(index_t& row_idx_min, value_t& val_min , index_t& row_idx_max, value_t& val_max)const;



	/// Copia il contenuto della matrice nelle strutture di una matrice compressa (formato a tre vettori)
	/// torna 0 in caso di successo, un codice di errore !=0 altrimenti altrimenti
	/// sz1:  (IN) dimensione di rpos, deve essere = numero delle righe +1
	/// sz2:  (IN) dimensione di cidx e val, deve essere = nonzeros()
	/// rpos: (OUT) vettore dimensione =sz1 che conterrà le posizioni di inizio delle righe 
	/// cidx: (OUT) vettore dimensione =sz2 che conterrà gli indici di colonna 
	/// val : (OUT) vettore dimensione =sz2 che conterrà i valori 
	int copy2CSR3(index_t sz1, index_t sz2, index_t* rpos, index_t* cidx, value_t*  val)const;


	/// Inserisce il contenuto della matrice nelle strutture di una matrice compressa (formato a tre vettori)
	/// torna 1 in caso di successo, 0 altrimenti
	int put2CSR3(index_t sz1, index_t	sz2, index_t* rpos, index_t* cidx, value_t*  val)const;


	//---------------------------------------------------------------
	// Aggiunte per compatibilita' con MML2::basic_sparse_matrix

	/// Inserisce una intera sottomatrice. Le componenti K[i][j] vengono inserite	nelle componenti (ir[i],ic[j]) a meno che 
	/// ir[i] o ic[j] siano zero.
	/// Se symut==true, tratta sia K che la matrice sparsa come simmetriche triangolari superiori
	/// e aggiunge solo le componenti del triangolo superiore di K nel triangolo superiore. 
	/// Assume che K punti alle componenti di una matrice ir_sz x ic_sz orientata alle righe.  
	/// Se gli indici sono fuori limite lancia un 'eccezione std::out_of_range()
	///@param ir vettore di indici di riga
	///@param dimensione del vettore ir coincidente con il numero di righe di K
	///@param ic vettore di indici di colonna
	///@param dimensione del vettore ic coincidente con il numero di colonne di K
	///@param K puntatore alla matrice K orientata alle righe di dimensione ir*ic
	template<typename I_T, typename V_T>
	void				put_mat(const I_T* ir, size_t ir_sz, const I_T* ic, size_t ic_sz, const V_T* K, bool symut = false);

	

	// Somma una intera sottomatrice. Le componenti K[i][j] vengono sommate alle componenti (ir[i],ic[j]) a meno che 
	// ir[i] o ic[j] siano zero.
	/// Se symut==true, tratta sia K che la matrice sparsa come simmetriche triangolari superiori
	/// e aggiunge solo le componenti del triangolo superiore di K nel triangolo superiore. 
	/// Assume che K punti alle componenti di una matrice ir_sz x ic_sz orientata alle righe.  
	/// Se gli indici sono fuori limite lancia un 'eccezione std::out_of_range()
	///@param ir vettore di indici di riga
	///@param dimensione del vettore ir coincidente con il numero di righe di K
	///@param ic vettore di indici di colonna
	///@param dimensione del vettore ic coincidente con il numero di colonne di K
	///@param K puntatore alla matrice K orientata alle righe di dimensione ir*ic
	template<typename I_T, typename V_T>
	void				add_mat(const I_T* ir, size_t ir_sz, const I_T* ic, size_t ic_sz, const V_T* K, bool symut = false);







	static dynamic_sparse_matrix  get_random_matrix(index_t nr, index_t nc, value_t sparsity, value_t mean, value_t range);
private:

	void				sym_put_(const iSet& ir, const iSet& ic, const inserter_matrix_t& K);
	void				gen_put_(const iSet& ir, const iSet& ic, const inserter_matrix_t& K);

	void				sym_sum_(const iSet& ir, const iSet& ic, const inserter_matrix_t& K);
	void				gen_sum_(const iSet& ir, const iSet& ic, const inserter_matrix_t& K);
	

	template<typename MAT>
	void				sym_sum_(const iSet& icr, const MAT& K);
	template<typename MAT>
	void				gen_sum_(const iSet& irc, const MAT& K);
	
	index_t				ncols_=0;
	std::vector<row_t*>	row_;
	
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







template<typename VAL, typename IDX, class MP, class ROW_T>
dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::dynamic_sparse_matrix(index_t R, index_t C, size_t bsz) 
	
{
	if (std::is_same<MP, M_PROP::SYM>::value && C!=R)
	{
		std::cerr << "dynamic_sparse_matrix<VAL,IDX,SYM,ROW>(index_t R, index_t C, size_t bsz): symmetric matrices must be square \n";
		assert(false);
	}
	

	ncols_ = C;
	row_.resize(R);
	for (size_t i = 0; i != R; ++i)
	{
		row_[i] = new row_t;
		row_[i]->value_at(i + 1) = VAL();
	}
	
}

template<typename VAL, typename IDX, class MP, class ROW_T>
dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::dynamic_sparse_matrix(const dynamic_sparse_matrix& rhs)
:ncols_(rhs.ncols_),  row_(rhs.nrows())
{
	for (size_t i = 0; i != row_.size(); ++i)
	{
		row_[i] = new row_t(*rhs.row_[i]);
	}

		
}


template<typename VAL, typename IDX, class MP, class ROW_T>
dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::dynamic_sparse_matrix(dynamic_sparse_matrix&& rhs)
:ncols_(rhs.ncols_), row_(std::move(rhs.row_))
{
	rhs.ncols_ = 0;
}


template<typename VAL, typename IDX, class MP, class ROW_T>
inline bool dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::test_subscripts(index_t r, index_t c)const
{
	if ((r > nrows() || r<1) || (c > ncols() || c<1))
		return false;
	return true;
}


		/// Ridimensiona il vettore di liste deallocando prima tutto il suo contenuto precedente
template<typename VAL, typename IDX, class MP, class ROW_T>
void dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::resize(index_t R, index_t C, size_t bsz)
	{
		destroy();
		dynamic_sparse_matrix	tmp(R,C,bsz);
		swap(tmp);
	}



	/// Distrugge la matrice deallocando le risorse
template<typename VAL, typename IDX, class MP, class ROW_T>
void dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::destroy()
{
	for (size_t i = 0; i != row_.size(); ++i)
	{
		delete row_[i];
		row_[i] = nullptr;
	}
	row_.clear();
	ncols_ = 0;
}

template<typename VAL, typename IDX, class MP, class ROW_T>
	void dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::swap(dynamic_sparse_matrix& rhs)
	{
		std::swap(ncols_,rhs.ncols_);
		std::swap(row_,rhs.row_);
	}


	

	template<typename VAL, typename IDX, class MP, class ROW_T>
	inline auto  dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::get_p(index_t i, index_t j)const ->const value_t*
	{
		if (i <1 || i> nrows())
			throw std::out_of_range("dynamic_sparse_matrix<>::get_p");
		if (IS::SYM && j < i)
			std::swap(i, j);
		return row_[i - 1]->find(j);
	}

	template<typename VAL, typename IDX, class MP, class ROW_T>
	inline auto  dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::get_p(index_t i, index_t j)-> value_t*
	{
		if (i <1 || i> nrows())
			throw std::out_of_range("dynamic_sparse_matrix<>::get_p");
		if (IS::SYM && j < i)
			std::swap(i, j);
		return row_[i - 1]->find(j);
	}

	
	template<typename VAL, typename IDX, class MP, class ROW_T>
	inline void dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::put(const iSet& ir, const iSet& ic, const inserter_matrix_t& K)
	{
		if (IS::SYM)
			sym_put_(ir, ic, K);
		else
			gen_put_(ir,  ic,  K);
	}


	template<typename VAL, typename IDX, class MP, class ROW_T>
	inline void dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::gen_put_(const iSet& ir, const iSet& ic, const inserter_matrix_t& K)
	{
		size_t NR = nrows(), NC = ncols();
		if (index_t(ir.max()) > nrows() || index_t(ic.max()) > ncols())
			throw std::out_of_range("dynamic_sparse_matrix::gen_put_");
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
				row_[R - 1]->value_at(C) = K(r, c);
			}
		}
	}
	


	template<typename VAL, typename IDX, class MP, class ROW_T>
	inline void dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::sym_put_(const iSet& ir, const iSet& ic, const inserter_matrix_t& K)
	{
		size_t NR = nrows(), NC = ncols(), R, C;
		if (index_t(ir.max())>nrows() || index_t(ic.max())> ncols())
			throw std::out_of_range("dynamic_sparse_matrix::sym_put_");


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
				row_[R - 1]->value_at(C) = K(r, c);
			}
		}
	}



	template<typename VAL, typename IDX, class MP, class ROW_T>
	inline void dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::sum(const iSet& ir, const iSet& ic, const inserter_matrix_t& K)
	{
		if (IS::SYM)
			sym_sum_(ir, ic, K);
		else
			gen_sum_(ir, ic, K);
	}

	template<typename VAL, typename IDX, class MP, class ROW_T>
	inline void dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::sum(const iSet& irc, const inserter_matrix_t& K)
	{
		if (IS::SYM)
			sym_sum_(irc, K);
		else
			gen_sum_(irc, K);
	}

	template<typename VAL, typename IDX, class MP, class ROW_T>
	inline void dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::sum(const iSet& irc, const inserter_sym_matrix_t& K)
	{
		if (IS::SYM)
			sym_sum_(irc, K);
		else
			gen_sum_(irc, K);
	}

	
	template<typename VAL, typename IDX, class MP, class ROW_T>
	inline void dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::gen_sum_(const iSet& ir, const iSet& ic, const inserter_matrix_t& K)
	{
		
		if (index_t(ir.max()) > nrows() || index_t(ic.max())>ncols())
			throw std::out_of_range("dynamic_sparse_matrix::gen_sum_");
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
				row_[R - 1]->value_at(C) += K(r, c);
			}
		}

	}

	template<typename VAL, typename IDX, class MP, class ROW_T>
	template<typename MAT>
	inline void dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::gen_sum_(const iSet& irc,  const MAT& K)
	{
		typedef typename MAT::value_t K_value_t;
		static_assert(std::is_same<VAL, K_value_t>::value, "type mismatch");
		index_t max_idx = irc.max();
		if (max_idx> nrows() || max_idx> ncols())
			throw std::out_of_range("dynamic_sparse_matrix::gen_sum_");
		
		//col_iterator it;
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
				row_[R - 1]->value_at(C) += K(r, c);
			}
		}
	}


	template<typename	VAL, typename	IDX, class MP, class ROW_T>
	inline void dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::sym_sum_(const iSet& ir, const iSet& ic, const inserter_matrix_t& K)
	{
		
		index_t  R, C;
		if (index_t(ir.max()) > nrows() || index_t(ic.max()) > ncols())
			throw std::out_of_range("dynamic_sparse_matrix::sym_sum_");

		col_iterator it;
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
				row_[R - 1]->value_at(C) += K(r, c);
				
			}
		}
	}

	template<typename VAL, typename IDX, class MP, class ROW_T>
	template<typename MAT>
	inline void dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::sym_sum_(const iSet& irc, const MAT& K)
	{
		if (index_t(irc.max())> nrows())
			throw std::out_of_range("dynamic_sparse_matrix::sym_sum_");

		index_t RR=0, CC = 0;
		col_iterator it;
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
				if (C >= R)
				{
					CC = C; RR = R;
				}
				else
				{
					CC = R; RR = C;
				}
				row_[RR - 1]->value_at(CC) += K(r, c);
			}
		}
	}




	template<typename VAL, typename IDX, class MP, class ROW_T>
	void	dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::fill(const value_t& val)
	{
		for (size_t i = 0; i != nrows(); ++i)
		{
			auto end = row_[i]->end();
			for (col_iterator beg = row_[i]->begin(); beg != end; ++beg)
				row_t::value(beg) = val;
		}
		
	}




	template<typename VAL, typename IDX, class MP, class ROW_T>
	void	dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::get_min_max_diag(index_t& idx_min, value_t& val_min, index_t& idx_max, value_t& val_max)const
	{
		idx_min=0;
		idx_max=0;
		val_min=std::numeric_limits<VAL>::max();
		val_max=std::numeric_limits<VAL>::min();

		VAL value;
		bool found;
		for(size_t i=0; i!=nrows(); ++i)
		{
			const value_t * v=row_[i].find(i + 1);
			
			if (!v)
				value = VAL(0);
			else
				value = *v;

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
	template<typename VAL, typename IDX, class MP, class ROW_T>
	inline VAL& dynamic_sparse_matrix<VAL, IDX, MP, ROW_T> ::put(index_t i, index_t j, const value_t& val)
	{
		if(!  test_subscripts(i,j))
			throw std::out_of_range("DSC_Matrix<VAL,IDX,ALLOCATOR>::put");

		if (IS::SYM && j < i)
			std::swap(i, j);
		return (row_[i - 1]->value_at(j) = val);
		
	}

	
	template<typename VAL, typename IDX, class MP, class ROW_T>
	inline VAL& dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::sum(index_t i, index_t j, const value_t& val)
	{
		if(!  test_subscripts(i,j))
			throw std::out_of_range("");
		if (IS::SYM && j < i)
			std::swap(i, j);

		return row_[i - 1]->value_at(j) += val;
		
	}
	
//}



	
	// Informazioni sulla struttura della matrice
	// In uscita NCR ha dimensione nrows e contiene  il numero di elementi diversi da zero di ogni riga 
	// Ritorna il numero totale di elementi 
	template<typename VAL, typename IDX, class MP, class ROW_T>
	inline size_t	dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::nonzeros(std::vector<index_t>& NCR)const
	{
		NCR.clear();
		NCR.reserve(nrows());
		size_t s, counter=0;
		for (auto r : row_)
		{
			NCR.push_back(s = r->size());
			 counter += s;
		}
		return counter;
	}


	template<typename VAL, typename IDX, class MP, class ROW_T>
	inline size_t	dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::nonzeros()const
	{
		size_t counter=0;
		for(auto r:row_)
			counter += r->size();
		return counter;
	}



	template<typename VAL, typename IDX, class MP, class ROW_T>
	void	dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::print(std::ostream& os)const
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
			auto end = row_[r]->end();
			for (auto it = row_[r]->begin(); it != end; ++it)
				os << "(" << std::setw(idx_width) << r + 1 << "," << std::setw(idx_width) << row_t::key(it) << ")= " << std::setw(width) << row_t::value(it) << std::endl;
		}
	}

	
			

		
		///@param fname nome del file
		///@return:  0 in caso di successo
		///@return: -1 nel caso il file non possa essere aperto
		///@return: -2 se si sono verificati errori in scrittura
	template<typename VAL, typename IDX, class MP, class ROW_T>
	int	dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::fwrite(const std::string& fname)const
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
			throw std::runtime_error("dynamic_sparse_matrix::fwrite");

		// Scrive un header con le dimensioni ed il tipo dei dati e della matrice
		std::uint64_t header_[8] = {
			(std::uint64_t)nrows(),
			(std::uint64_t)ncols(),
			(std::uint64_t)nnz,
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
		os.write(reinterpret_cast<const char*>(&row_pos[0]), sizeof(index_t)* row_pos.size());
		if (!os)
			return -2;
		// scrive tutti gli nnz indici di colonna
		std::vector<IDX>	i_buffer;
		for( index_t r=0; r!=nrows();++r)
		{
			auto end = row_[r]->end();
			i_buffer.clear();
			i_buffer.reserve(row_[r]->size());
			for (auto it = row_[r]->begin(); it != end; ++it)
				i_buffer.push_back(row_t::key(it));
			if(i_buffer.size())
			{
				os.write( reinterpret_cast<const char*>(i_buffer.data()), sizeof(index_t)* i_buffer.size());
				if (!os)
					return -1;
			}
		}



		// scrive tutte le nnz componenti
		std::vector<value_t>	v_buffer;
		// scrivo i valori 
		for( index_t r=0; r!=nrows();++r)
		{
			auto end = row_[r]->end();
			v_buffer.clear();
			v_buffer.reserve(row_[r]->size());
			for (auto it = row_[r]->begin(); it != end; ++it)
				v_buffer.push_back(row_t::value(it));
			if(v_buffer.size())
			{
				os.write(reinterpret_cast<const char*>(v_buffer.data()), sizeof(value_t)* v_buffer.size());
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
		
	template<typename VAL, typename IDX, class MP, class ROW_T>
	int	dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::fread(const std::string& fname)
	{

		std::ifstream is(fname,std::ios_base::binary);
		if(!is.is_open())
			return -1; 
		destroy();

		std::uint64_t this_header_[8] = { 
			(std::uint64_t)0,
			(std::uint64_t)0,
			(std::uint64_t)0,
			(std::uint64_t)(is_symmetric() ? 1 : 0),
			(std::uint64_t)(is_row_major() ? 1 : 0),
			(std::uint64_t)(IS::UT ? 1 : 0),
			std::uint64_t(type<index_t>::id()),
			(std::uint64_t)type<value_t>::id()
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
			row_[r]->put(pi, v_buff.data(), r_sz);
			/*for (size_t j = 0; j != r_sz; ++j)
				row_[r].value_at(pi[j]) = v_buff[j];*/
			pi+=r_sz;
		}
		return 0;
	}



	
	template<typename VAL, typename IDX, class MP, class ROW_T>
	auto dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::get_random_matrix(index_t nr, index_t nc, value_t sparsity, value_t mean, value_t range)->dynamic_sparse_matrix
	{
		dynamic_sparse_matrix A(nr, nc);
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
				A.put(R, C, uniform_dist_val(e1));
			}
			// inserts always the diagonal element
			if ( R <=nc)
				A.put(R, R, uniform_dist_val(e1));
		}
		return A;
	}



	template<typename VAL, typename IDX, class MP, class ROW_T>
	int	dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::copy2CSR3(index_t sz1, index_t sz2, index_t* rpos, index_t* cidx, value_t*  val)const
	{
		size_t	nnz = 0;
		size_t nr = nrows();
		if (nr != (sz1 - 1))
			return -1;
		
		rpos[0] = 1;
		for (size_t i = 0; i != nr; ++i)
		{
			auto sz = row_[i]->size();
			rpos[i + 1] = rpos[i] + sz;
			nnz += sz;
		}
		if (nnz != sz2)
			return -2;

		
		for (index_t r = 0; r != nr; ++r)
		{
			for (row_t::const_iterator it = row_[r]->begin(), end = row_[r]->end(); it != end; ++it)
			{
				*cidx++ = it->first;
				*val++  = it->second;
			}
		}
		return 0;
	}



	template<typename VAL, typename IDX, class MP, class ROW_T>
	int	dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::put2CSR3(index_t NRP1, index_t NNZ, index_t* RPOS, index_t* CIDX, value_t*  VAL)const
	{
		size_t nr = nrows();
		if (NRP1 != (nr + 1))
			return -1;

		for (index_t r = 0; r != nr; ++r)
		{
			index_type csrPOS = RPOS[r] - 1;
			index_type csrEND = RPOS[r + 1] - 1;
			for (auto it = row_[r]->begin(), end = row_[r]->end(); it != end; ++it)
			{
				index_t cidx = it->first;
				bool found = false;
				for (; csrPOS < csrEND; ++csrPOS)
				if (CIDX[csrPOS] == cidx)
				{
					found = true;
					break;
				}
				if (found)
					VAL[csrPOS] = it->second;
				else
					return -3;
			}
		}
		return 0;
	}



	template<typename VAL, typename IDX, class MP, class ROW_T>
	template<typename I_T, typename V_T>
	void	dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::put_mat(const I_T* ir, size_t ir_sz, const I_T* ic, size_t ic_sz, const V_T* K, bool symut)
	{

		for (size_t i = 0; i != ir_sz; ++i)
		{
			index_t R = (index_t)ir[i];
			if (!R)
				continue;
			for (size_t j = 0; j != ic_sz; ++j)
			{
					index_t C = ic[j];
					if (!C)
						continue;
					if (C>=R)
						put(R, C, K[i*ic_sz + j]);
			}

		}
	}

	template<typename VAL, typename IDX, class MP, class ROW_T>
	template<typename I_T, typename V_T>
	void	dynamic_sparse_matrix<VAL, IDX, MP, ROW_T>::add_mat(const I_T* ir, size_t ir_sz, const I_T* ic, size_t ic_sz, const V_T* K, bool symut)
	{

		for (size_t i = 0; i != ir_sz; ++i)
		{
			index_t R = (index_t)ir[i];
			if (!R)
				continue;
			for (size_t j = 0; j != ic_sz; ++j)
			{
				index_t C = ic[j];
				if (!C)
					continue;
				if (C>=R)
					sum(R, C, K[i*ic_sz + j]);
			}

		}
	}



} // end namespace MML

#endif