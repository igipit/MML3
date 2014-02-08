#pragma once

#include"MML3-Matrix.h"
#include"MML3-static_CSR_Matrix.h"
#include"MML3-Vector.h"
#include<map>
#include<array>
#include<iostream>
#include<string>

namespace MML3
{

	/**
	* solutore per la determinazione dei primi modi propri generalizzati (l_i, Z_i) che soddisfano le eq. A Z_i = li B Z_i per matrici CSR3  simmetriche </summary>
	* @param A : matrice sparsa nxn
	* @param B : matrice sparsa nxn
	* @param Z : matrice densa ordinata per colonne  n x m0, dove m0 >=m e' la dimensione del sottospazio ed m e' il numero di autovettori contenuti nel range specificato
	* @param L : vettore m0 che in uscita, in caso di successo, conterra' m autovalori (la routine ridimensiona L ad m)
	* @eigen_error : vettore m0 che in uscita conterra' l'errore su cianscun autovalore
	* @param toll_exp  : esponente della tolleranza con cui cercare gli autovalori (0 per il valore di default)
	* @param max_iter  : numero massimo di iterazioni interne (0 per il valore di default)
	* @log	   : stream su cui mandare  messaggi a runtime
	* @return  : 0 in caso di successo
	*			 1 se non sono stati trovati autovalori nell'intervallo prescritto
	*			 2 se non e' stata raggiunta laq convergenza nel numero massimo di cicli prescritto
	*            3 se la dimensione del sottospazio e' troppo piccola per contenere gli autovettori
	*           -3 se la matrice B non e' definita positiva
	*           altri valori in caso di errori severi
	*/
	int SparseEigenSolver(static_sparse_CSR_Matrix<double, std::int32_t, M_PROP::SYM>& A,
		static_sparse_CSR_Matrix<double, std::int32_t, M_PROP::SYM>& B,
		Matrix<double, M_PROP::GE, M_SHAPE::RE, M_ORD::COL>& eigen_vector,
		Vector<double>& eigen_value,
		Vector<double>& eigen_error,
		double eigenv_max,
		int    toll_exp,
		int    max_iter,
		std::ostream& log);


} // end namespace MML3