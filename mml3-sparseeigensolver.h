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

	// solutore per la determinazione dei primi modi propri generalizzati (l_i, Z_i) per matrici CSR3  simmetriche 
	// A Z_i = li B Z_i
	int SparseEigenSolver(static_sparse_CSR_Matrix<double, std::int32_t, M_PROP::SYM>& A,
		static_sparse_CSR_Matrix<double, std::int32_t, M_PROP::SYM>& B,
		Matrix<double, M_PROP::GE, M_SHAPE::RE, M_ORD::COL>& eigen_vector,
		Vector<double>& eigen_value,
		Vector<double>& eigen_error,
		double eigenv_max,
		std::array<int, 32>& param,
		std::ostream& log);


} // end namespace MML3