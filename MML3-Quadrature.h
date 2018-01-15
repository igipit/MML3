#pragma once

namespace MML3
{


/** Formule di quadratura di Gauss nel dominio 1D [-1,+1]*/
class gauss_quadrature
{
public:

	// copia nella matrice A i punti (prima riga)  ed i pesi (seconda riga) della formula di quadratura ad ngp punti
	// ritorna il numero di punti in caso di successo e <=0 se la formula non e' implementata
	// ARRAY deve essere ridimensionabile con ARRAY::resize(NR,NC) e deve avere l'operatore di accesso ARRAY::operator()(i,j)
	template< class ARRAY>
	static int get(int ngp, ARRAY& A)
	{
		if(ngp >size_ || ngp<1)
			return -1;
		A.resize(2,ngp);
		for(int i=0; i!=ngp;++i)
		{
			A(1,i+1)=X_[ngp-1][i];
			A(2,i+1)=W_[ngp-1][i];
		}
		return ngp;
	}
	// torna un puntatore agli ngp punti della formula di quadratura o zero se la formula non e' implementata 
	static const double*	getX(int ngp);
	// torna un puntatore agli ngp pesi della formula di quadratura o zero se la formula non e' implementata 
	static const double*	getW(int ngp);

	// numero massimo di punti di quadratura per cui la formula è implementata
	static int	max_points(){return size_;}
private:
	enum {size_=10};
	static const double X_[size_][size_];
	static const double W_[size_][size_];
};




/** Formule di quadratura di Gauss nel dominio 2D [-1,+1] X [-1,+1] */
class gauss_quadrature_quad
{
public:
	// copia nella matrice A le coordinate X (prima riga), Y(seconda riga)  ed i pesi (terza riga) della formula di quadratura ad ngpx * ngpy punti
	// ritorna il numero di punti totali in caso di successo e <=0 se la formula non e' implementata
	// ARRAY deve essere ridimensionabile con ARRAY::resize(NR,NC) e deve avere l'operatore di accesso ARRAY::operator()(i,j)
	template< class ARRAY>	
	static int	get(int ngpx, int ngpy, ARRAY& A)
	{
	
		if((ngpx> gauss_quadrature::max_points() || ngpx<1) || (ngpy> gauss_quadrature::max_points() || ngpy<1))
			return -1;
		A.resize(3,ngpx*ngpy);
		const double * x=gauss_quadrature::getX(ngpx);
		const double * y=gauss_quadrature::getX(ngpy);
		const double * wx=gauss_quadrature::getW(ngpx);
		const double * wy=gauss_quadrature::getW(ngpy);
		int pos=1;
		for(int i=0; i!=ngpy; ++i)
			for(int j=0; j!=ngpx; ++j)
			{
				A(1,pos)=x[j];
				A(2,pos)=y[i];
				A(3,pos)=wy[i]*wx[j];
				++pos;
			}
		return ngpx*ngpy;
	}
};



/** Formule di quadratura  nel dominio 2D triangolare (triangolo rettangolo con lati unitari) */
// copia nella matrice A le coordinate X (prima riga), Y(seconda riga)  ed i pesi (terza riga) della formula di quadratura ad ngp  punti
// ritorna il numero di punti totali in caso di successo e <= 0se la formula non e' implementata
// ARRAY deve essere ridimensionabile con ARRAY::resize(NR,NC) e deve avere l'operatore di accesso ARRAY::operator()(i,j)
class gauss_quadrature_tria
{
public:

	template< class ARRAY>
	static int get(int ngp, ARRAY& A)
	{

		if((ngp >size_ || ngp<1) || triangleGaussianQuadratureWeight[ngp-1][0]==0.0 )
			return -1;
		A.resize(3,ngp);
		for(int i=0; i!=ngp;++i)
		{
			A(1,i+1)=triangleGaussianQuadraturePoint[ngp-1][i][0];
			A(2,i+1)=triangleGaussianQuadraturePoint[ngp-1][i][1];
			A(3,i+1)=triangleGaussianQuadratureWeight[ngp-1][i];
		}
		return ngp;

	}

private:
	enum{size_=7};
	static double triangleGaussianQuadraturePoint [size_][size_][2];
	static double triangleGaussianQuadratureWeight[size_][size_];


}; // end class 





} // end namespaces MML3

