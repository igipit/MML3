#pragma once
#include"MML3-matrix.h"

namespace MML3
{

	// Common base and interface of sparse matrices
	template<typename VAL_TYPE, typename IDX_TYPE>
	class sparseMatrix
	{
	public:

		typedef IDX_TYPE														index_t;	///< il tipo dell'indice
		typedef VAL_TYPE														value_t;	///< il tipo della'componente
		typedef iSet															iset_t;
		typedef Matrix<VAL_TYPE, M_PROP::GE, M_SHAPE::RE, M_ORD::ROW>			inserter_matrix_t;
		typedef Matrix<VAL_TYPE, M_PROP::SYM, M_SHAPE::RE, M_ORD::ROW>			inserter_sym_matrix_t;
		//-------------------------------------------
		// CTORs
		//-------------------------------------------
		sparseMatrix()=default;
		sparseMatrix(const sparseMatrix& )=default;
		virtual ~sparseMatrix() = default;

		
		
		// deallocates all the resources and resize ther matrix to 0x0
		virtual void				destroy()=0;

		///@return il numero totale di elementi  diversi da zero
		virtual size_t				nonzeros()const=0;
		///@return il numero di colonne
		virtual index_t				ncols()const=0;
		///@return il numero di righe
		virtual index_t				nrows()const=0;

		virtual bool				is_row_major()const = 0;

		virtual bool				is_symmetric()const = 0;
		///
		virtual bool				test_subscripts(index_t r, index_t c)const=0;

		virtual void				fill(const value_t&)=0;
		


		/// Ritorna il puntatore alla componente in posizione (i,j) se questa esiste, nullptr altrimenti
		virtual const	value_t*	get_p(index_t i, index_t j)const=0;
		virtual			value_t*	get_p(index_t i, index_t j)=0;



		/// Inserisce nella componente (i,j), creandola se non c'e' il valore val
		/// Se gli indici sono fuori limite lancia un 'eccezione std::out_of_range()
		///@param  i indice di riga
		///@param  j indice di colonna
		///@param value il valore da inserire
		///@return  componente inserita per riferimento. 
		virtual value_t&			put(index_t i, index_t j, const value_t& value)=0;

		/// Somma  alla componente (i,j) il valore val. Se la componente non c'è la crea.
		/// Se gli indici sono fuori limite lancia un 'eccezione std::out_of_range()
		///@param  i indice di riga
		///@param  j indice di colonna
		///@param value il valore da inserire
		///@return  componente modificata per riferimento. 
		virtual  value_t&			sum(index_t i, index_t j, const value_t& val)=0;

		/// Inserisce nella sparsa una intera sottomatrice densa K. Le componenti K[i][j] vengono inserite	nelle componenti (ir[i],ic[j]) a meno che 
		/// ir[i] o ic[j] siano zero.
		/// Se gli indici sono fuori limite lancia un 'eccezione std::out_of_range()
		///@param ir vettore<I> di indici di riga
		///@param ic vettore<I> di indici di colonna
		///@param K  matrice<T> densa  ir.size x ic.size
		virtual void				put(const iSet& ir, const iSet& ic, const inserter_matrix_t& K) = 0;


		/// sum(...,matrix k)
		///
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
		virtual void				sum(const iSet& ir, const iSet& ic, const inserter_matrix_t& K) = 0;
		// specialization for the frequent case  ir=ic=irc
		virtual void				sum(const iSet& irc, const inserter_matrix_t& K) = 0;
		// specialization for the frequent case where the inserted matrix is symmetric
		virtual void				sum(const iSet& irc, const inserter_sym_matrix_t& K) = 0;
		
		/// Copia il contenuto della matrice nelle strutture di una matrice compressa (formato a tre vettori)
		/// torna 1 in caso di successo, 0 altrimenti
		///@param rpos_sz:(IN)  dimensione di rpos, deve essere >= numero delle righe +1
		///@param nnz:    (IN)  dimensione di cidx e val, deve essere >= nonzeros()
		///@param rpos:   (OUT) vettore dimensione >=sz1 che conterrà le posizioni di inizio delle righe 
		///@param cidx:   (OUT) vettore dimensione >=sz2 che conterrà gli indici di colonna 
		///@param val :   (OUT) vettore dimensione >=sz2 che conterrà i valori 
		///@return   0 in case of success, <0 else
		//virtual  int					copy2(index_t rpos_sz, index_t nnz, index_t rpos[/*rpos_sz*/], index_t cidx[/*nnz*/], value_t  val[/*nnz*/])const=0;

		
		
		
		virtual void				print(std::ostream& os)const=0;

		
		///fwrite(const std::string& fname)
		///Scrive su file in formato binario universale la matrice sparsa .
		///@param fname nome del file
		///@return:  0 in caso di successo
		///@return: -1 nel caso il file non possa essere aperto
		///@return: -2 se si sono verificati errori in scrittura
		virtual int fwrite(const std::string& fname)const=0;

		///fread(const std::string& fname)
		///legge da  file in formato binario universale la matrice sparsa .
		///@param fname nome del file
		///@param bsize
		///@return:  0 in caso di successo
		///@return: -1 nel caso il file non possa essere aperto
		///@return: -2 se si sono verificati errori in lettura
		///@return: -3 se si sono verificati errori di incompatibilita' dei tipi
		///@return: -4 se la dimensione della matrice su file non e' gestibile 
		virtual int fread(const std::string& fname) = 0;
	};


	template<typename	VAL, typename	IDX>
	std::ostream& operator<<(std::ostream& os, const sparseMatrix<VAL, IDX>& V)
	{
		V.print(os);
		return os;
	}




} // end namespace
