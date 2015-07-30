#pragma once
#include"MML3-CSR3_Matrix.h"
namespace MML3
{

	// Risolve il sistema lineare A X=B e sovrascrive B con la soluzione X
	// INPUT
	//		A:			matrice simmetrica sparsa statica
	//		B:			matrice dei termini noti, ogni colonna è un termine noto
	//      log:		std::ostream su cui verrà diretto il log
	//		threads:	numero massimo di threads che l'algoritmo dovrà usare
	// OUTPUT
	//      B:			soluzione del sistema lineare
	//		log:		il log
	//      RET_VAL		torna 0 in caso di successo, un codice di errore in caso contrario
	int LinearSparseSolver(CSR3_sym_matrix<double, std::int32_t>& A,
		Matrix<double, M_PROP::GE, M_SHAPE::RE, M_ORD::COL>& B,
		std::ostream& log,
		int threads);


	int LinearSparseSolver(CSR3_matrix<double, std::int32_t>& A,
		Matrix<double, M_PROP::GE, M_SHAPE::RE, M_ORD::COL>& B,
		std::ostream& log,
		int threads);

}// end namespace