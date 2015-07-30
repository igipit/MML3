#pragma once
#include"MML3-math.h"
#include"MML3-array-algorithms.h"
#include"MML3-geometry.h"
#include<vector>



namespace MML3
{
	

	
	//===================================================================================================
	// Ramo cubico di una curva nello spazio N dimensionale, parametrizzata con t \in [t_a, t_b]
	template<int N>
	class CubicSplineCurve
	{
	public:
		CubicSplineCurve() :x1_(0), h_(1.){}
		CubicSplineCurve(double ta, double tb, const double z[/*N*4*/]){ configure(ta,tb,z); }
		CubicSplineCurve(const double z[/*N*4*/]){ configure(z); }


		//-------------------------------------------------------------------------------------------------------
		// Calcola  coefficienti della curva r(t) :[0, 1] -> R^N, a partire dai valori di r=(r_x, r_y, ...)  e delle sue derivate negli estremi
		// Siccome r è un vettore di dimensione N, devo passare un array N*4 di valori organizzati come segue
		//        [ r_x(0),  r'_x(0), r_x(1), r'_x(1) ]
		//     z= [ r_y(t0), r'_y(0), r_y(1), r'_y(1) ]
		//        [ ...                                      ]
		//  Se i contributi di grado 2 e\o 3 sono  trascurabili nel senso che la massima variazione Dr della curva rispetto all'andamento rettilineo
		//  è tale che || Dp || < tolerance *||  r(1) - r(0)||
		//  allora i relativi coefficienti vengono annullati
		void	configure(const double z[/*N*4*/], double tolerance=std::numeric_limits<double>::epsilon());

		// Calcola  coefficienti della curva r(t) :[t_a, t_b] -> R^N, a partire dai valori di r ed r' negli estremi
		// Siccome r è un vettore di dimensione N, devo passare un array N*4 di valori organizzati come segue
		//        [ r_x(t_a), r'_x(t_a), r_x(t_b), r'_x(t_b) ]
		//     z= [ r_y(t_a), r'_y(t_a), r_y(t_b), r'_y(t_b) ]
		//        [ ...                                      ]
		// 
		void	configure(double ta, double tb, const double z[/*N*4*/]);


		//-------------------------------------------------------------------------------------------------
		// una volta chiamato configure, questa funzione calcola efficientemente il valore, la derivata prima 
		// e seconda della curva  nel punto t specificato
		void	eval(double t, double value[N], double deriv[/*N*/], double deriv2[/*N*/])const;


		// ritorna i coefficienti 
		const double* cx()const{ return c_; }
		const double* cy()const{ static_assert(N > 1,""); return c_ + 4; }
		const double* cz()const{ static_assert(N > 2,""); return c_ + 8; }

	private:
		double	x1_,  /* estremo sinistro dell'intervallo*/
				h_;  /* ampiezza dell'intervallo */

		// array N*4 dei coefficienti tali che (convenzione d Einstein)
		// r_i(t) = c_[i][j] * t^j,      i \in [0,N), j\in [0,4)
		double c_[N * 4];
	};


	



	//===========================================================================================================================
	// Interpolatore di tipo Nagata C0 sul triangolo
	// Nagata T.,  Simple local interpolation of surfaces using normal vectors, Computer Aided Geometric Design 22 (2005) 327–347
	// Esegue l'interpolazione di una superficie in R^3 parametrizzata con due parametri t={t_0,t_}  sul triangolo  unitario (u,v) \in [0,1]x[0,1] : v <=1-u
	// a partire dalla conoscenza dei punti e delle normali nei vertici della patch triangolare 

	class NagataC0Triangle
	{

		typedef MML3::StaticArrayOp<3, double>	V3Op;
		typedef MML3::GEO::Point<3, double>		Point3;
		typedef MML3::GEO::Vector<3, double>	Vector3;
	public:
		//-----------------------------------------------------------------------------------------------
		// Configura una patch  alla Nagata a partire dalle coordinate e dalle normali dei 3  vertici 
		// ritorna >=0 in caso di successo
		// Deve essere chiamata prima di ogni altra funzione
		int config(const double* P[3], const double* n[3]);
		

		//-----------------------------------------------------------------------------------------------
		// Determina le coordinate e la normale della superficie  nel punto immagine di (u,v)=(t[0],t[1]) 
		// se n==NULL non determina le normali
		void get_at(const double t[2], double r[3], double n[3])const;

		// //Calcola le derivate prime rispetto alle cordinate naturali nel punto (u,v)=(t[0],t[1])
		// OUT
		//	dp_u=[dx/du, dy/du, dz/du]
		//  dp_v=[dx/dv, dy/dv, dz/dv]
		void get_derivatives_at(const double t[2], double dp_u[3], double dp_v[3])const;

		//Calcola le derivate seconde rispetto alle cordinate naturali nel punto t[]
		// IN
		//		t[]=(u,v)
		// OUT
		//	dd=[	d^2x/du^2, d^2x/dv^2, d^2x/dudv, 
		//			d^2y/du^2, d^2y/dv^2, d^2y/dudv,
		//			d^2z/du^2, d^2z/dv^2, d^2z/dudv ]
		void get_second_derivatives(const double t[2], double dd[9])const;


		//==============================================================================================
		// Determina l'intersezione del raggio  [0,0,0] + alpha r[]  con la superficie d'interpolazione 
		// INPUT
		//		r[] la direzione  in Re^3 del raggio, non necessariamente unitaria
		// OUTPUT
		//		alpha:	moltiplicatore di r[]
		//		t[] :	le coordinate naturali  del punto di intersezione
		//		RetVal: >=0 in caso di successo
		//
		// Questa funzione non è altro che un wrapper che chiama una delle sottostanti
		int get_ray_intersection(const double r[3], double& alpha, double t[2], double parameter_bounds_tol, double  tolerance )const;


		// Versione 0
		int ray_intersection_Newton_(const double r[3], double& alpha, double t[2], double  tolerance)const;
		// Versione 1
		int ray_intersection_first_order_(const double r[3], double& alpha, double t[2], double  tolerance)const;



	private:


		// coefficienti del polinomio quadratico completo nei termini delle coordinate naturali  (x,y) \in [0,1] x [0,1]:
		// r_i(x,y) =	 a_[0][0][i] +  a_[1][0][i]* x + a_[0][1][i]* y + 
		//				+a_[2][0][i]* x^2 + a_[1][1][i]* x * y + a_[0][2][i]* y^2
		double a_[3][3][3];

		// coordinate dei tre vertici. Queste sono ridondanti in quanto possono essere ricavate da a_[][][]
		// mediante operazioni algebriche, ma ne tengo copia per velocizzae le operazioni
		double P_[3][3];


		// work horse per il calccolo delle curvature sui lati del triangolo
		void get_c_(const double* d, const double* n0, const double* n1, double* c)const;
	};



	

		//===========================================================================================================================
		// Interpolatore di tipo Nagata C0 
		// Nagata T.,  Simple local interpolation of surfaces using normal vectors, Computer Aided Geometric Design 22 (2005) 327–347
		// Esegue l'interpolazione di una superficie in R^3 parametrizzata con due parametri t={t_0,t_}  sul quadrato unitario [0,1]x[0,1]
		// a partire dalla conoscenza dei punti e delle normali nei vertici della patch triangolare o quadrilatera

		class NagataC0_Patch
		{

			typedef MML3::StaticArrayOp<3, double>	V3Op;
			typedef MML3::GEO::Point<3,double>		Point3;
			typedef MML3::GEO::Vector<3, double>		Vector3;
		public:
			//-----------------------------------------------------------------------------------------------
			// Configura una patch  alla Nagata a partire dalle coordinate e dalle normali dei 3 o 4 vertici 
			// chiama una delle due versioni specializzate:
			//      se P[3]==nullptr  chiama config_tria
			//		altrimenti chiama config_quad
			// ritorna >=0 in caso di successo
			// Deve essere chiamata prima di ogni altra funzione
			int config(const double* P[4], const double* n[4]);
			int config_quad(const double* P[4], const double* n[4]);
			int config_tria(const double* P[3], const double* n[3]);


			//-----------------------------------------------------------------------------------------------
			// Determina le coordinate e la normale della superficie  nel punto immagine di (t[0],t[1]) 
			// le coordinate in input t[] sono quelle naturali in [0,1]x[0,1]
			// se n==NULL non determina le normali
			void get_at(const double t[2], double r[3], double n[3])const;

			// //Calcola le derivate prime rispetto alle cordinate naturali nel punto t[]
			// OUT
			//	d=[dr_1/dt_1, dr_2/dt_1, dr_3/dt_1, dr_1/dt_2, dr_2/dt_2, dr_3/dt_2]
			void get_derivatives_at(const double t[2], double d[6])const;

			//Calcola le derivate seconde rispetto alle cordinate naturali nel punto t[]
			// OUT
			//	dd
			void get_second_derivatives(const double t[2], double dd[9])const;


			//==============================================================================================
			// Determina l'intersezione del raggio  [0,0,0] + alpha r[]  con la superficie d'interpolazione 
			// INPUT
			//      P[] le coordinate dei 3 o 4  vertici
			//		r[] la direzione  in Re^3 del raggio, non necessariamente unitaria
			// OUTPUT
			//		alpha:	moltiplicatore di r[]
			//		t[] :	le coordinate naturali  del punto di intersezione
			//		RetVal: >=0 in caso di successo
			//
			// Questa funzione non è altro che un wrapper che chiama una delle sottostanti
			int get_ray_intersection(const double* P[4], const double r[3], double& alpha, double t[2], double  tolerance=1e-7)const;


			// Versione 0
			int ray_intersection_Newton_(const double* P[4], const double r[3], double& alpha, double t[2], double  tolerance)const;
			// Versione 1
			int ray_intersection_first_order_(const double* P[4], const double r[3], double& alpha, double t[2], double  tolerance)const;



		private:

			
			// coefficienti del polinomio quadratico completo nei termini delle coordinate naturali  (x,y) \in [0,1] x [0,1]:
			// r_i(x,y) =	 a_[0][0][i] +  a_[1][0][i]* x + a_[0][1][i]* y + 
			//				+a_[2][0][i]* x^2 + a_[1][1][i]* x * y + a_[0][2][i]* y^2
			double a_[3][3][3];

			bool  is_triangle_;


			void get_c_(const double* d, const double* n0, const double* n1, double* c)const;
		};




		//  Cubiche di Hermite sull'intervallo [0,1]
		struct CubicHermite
		{

			static double base(int i, double t)
			{
				double tt = t*t;
				switch (i)
				{
				case 0:
					return 1.0 + tt*(-3. + 2. * t);
				case 1:
					return  t * (1. + t*(-2. + t));
				case 2:
					return  tt*(3. - 2. * t);
				case 3:
					return  tt*(t - 1.);
				default:
					throw std::out_of_range("MML3::CubicHermite::base");
				}
			}


			static void eval_all_basis(double t, double y[], double dy[], double ddy[])
			{

				if (y != nullptr)
				{
					double tt = t*t;

					y[0] = 1 + tt*(2 * t - 3);
					y[1] = t * (1 + t*(t - 2));
					y[2] = 1 - y[0];
					y[3] = tt*(t - 1);
				}
				if (dy != nullptr)
				{

					// derivate
					dy[0] = (6 * t) * (t - 1);
					dy[1] = 1 + t*(3 * t - 4);
					dy[2] = -dy[0];
					dy[3] = t*(3 * t - 2);
				}
				if (ddy != nullptr)
				{
					// derivate seconde
					double six_t = 6.*t;

					ddy[0] = six_t + six_t - 6;
					ddy[1] = six_t - 4;
					ddy[2] = -ddy[0];
					ddy[3] = six_t - 2;
				}
			}
		};








		//////////////////////////////////////////////////////////////////////////////////////////////
		// NAGATAC0_TRIA IMPL

		inline void  NagataC0Triangle::get_c_(const double* d, const double* n0, const double* n1, double* c)const
		{
			// coseno dell'angolo tra le due normali 
			double co = n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2]; 
			// sin^2 dell'angolo
			double si2 = 1. - co*co; 
			// il  test è positivo se l'angolo tra le due normali è inferiore ad 1/100 di grado sessagesimale
			// in tal caso le due normali vengono considerate parallele e la curva trea i due punti un segmento 
			// rettilineo (quindi il vettore curvatura = 0)
			if (si2 < 3E-8)
			{
				V3Op::zero(c);
				return;
			}
			double dn0 = V3Op::dot(d, n0);
			double dn1 = V3Op::dot(d, n1);
			double alpha_1 = -1 / si2 * (dn1 + co*dn0);
			double alpha_0 = dn0 - alpha_1 * co;
			V3Op::lin_comb(alpha_0, n0, alpha_1, n1, c);
		}



		// versione modificata con i parametri (u,v) : 0<= u <=1,   v <=1-u
		inline int NagataC0Triangle::config(const double* P[3], const double* n[3])
		{
			double d[3][3]; //vettori dei lati
			double c[3][3]; //vettori curvatura dei lati

			// calcolo dei vettori curvatura per ogni lato
			get_c_( V3Op::sub(P[1], P[0], d[0]), n[0], n[1], c[0]);
			get_c_( V3Op::sub(P[2], P[1], d[1]), n[1], n[2], c[1]);
			get_c_( V3Op::sub(P[2], P[0], d[2]), n[0], n[2], c[2]);

			// copia interna delle coordinate dei vertici
			for (auto i = 0; i != 3; ++i)
				V3Op::copy(P[i], P_[i]);


			V3Op::copy(P[0], a_[0][0]);			//a_00=P_0
			V3Op::sub(d[0], c[0], a_[1][0]);	//a_10=d_0-c_0
			V3Op::sub(d[2], c[2], a_[0][1]);	//a_01=d_2-c_2
			V3Op::copy(c[0], a_[2][0]);			//a_20=c_0
			V3Op::copy(c[2], a_[0][2]);			//a_02=c_2
						
			// a_11=c_0 -c_1 +c_2    in two steps
			V3Op::sum(c[2], c[0], a_[1][1]);		
			V3Op::sub(a_[1][1], c[1], a_[1][1]); 

			// NB: a_[2][2],a_[1][2],a_[2][1] not used

			return 0;
		}


		// versione originale di Nagata con i parametri u,v:   0 <= v <= u <=1
		//inline int NagataC0Triangle::config(const double* P[3], const double* n[3])
		//{
		//	double d[3][3]; //vettori dei lati
		//	double c[3][3]; //vettori curvatura dei lati

		//	// calcolo dei vettori curvatura per ogni lato
		//	get_c_(V3Op::sub(P[1], P[0], d[0]), n[0], n[1], c[0]);
		//	get_c_(V3Op::sub(P[2], P[1], d[1]), n[1], n[2], c[1]);
		//	get_c_(V3Op::sub(P[2], P[0], d[2]), n[0], n[2], c[2]);

		//	// copia interna delle coordinate dei vertici
		//	for (auto i = 0; i != 3; ++i)
		//		V3Op::copy(P[i], P_[i]);


		//	V3Op::copy(P[0], a_[0][0]);	//a_00=P_0
		//	V3Op::copy(c[0], a_[2][0]);	//a_20=c_0
		//	V3Op::copy(c[2], a_[0][2]);	//a_02=c_2

		//	V3Op::sub(d[0], c[0], a_[1][0]);	//a_10=d_0-c_0

		//	// a_01= d_1 + c_0 - c_2  in two steps
		//	V3Op::sum(d[1], c[0], a_[0][1]);	 //a_01 =d_1 +c_0
		//	V3Op::sub(a_[0][1], c[2], a_[0][1]);// a_01 =a_01 -c_2

		//	// a_11=c_2 -c_0-c_1    in two steps
		//	V3Op::sub(c[2], c[0], a_[1][1]);
		//	V3Op::sub(a_[1][1], c[1], a_[1][1]);


		//	// NB: a_[2][2],a_[1][2],a_[2][1] not used

		//	return 0;
		//}


		inline void NagataC0Triangle::get_at(const double t[2], double r[3], double n[3])const
		{
			const double u = t[0], v = t[1];
			const double	uu = u*u,
							vv = v*v,
							uv = u*v;
			
			if (r != nullptr)
			{
				for (auto i = 0; i != 3; ++i)
				{
					r[i] = a_[0][0][i] + a_[1][0][i] * u + a_[0][1][i] * v +
						a_[2][0][i] * uu + a_[1][1][i] * uv + a_[0][2][i] * vv;
				}
			}

			if (n != nullptr)
			{
				double t0[3], t1[3]; // vettori coordinati nel punto
				double _2u = 2 * u, _2v=2*v;
				for (auto i = 0; i != 3; ++i)
				{
					t0[i] = a_[1][0][i] + _2u*a_[2][0][i] + v*a_[1][1][i];
					t1[i] = a_[0][1][i] + _2v*a_[0][2][i] + u*a_[1][1][i];
				}
				V3Op::cross_product(t0, t1, n);
				V3Op::unitize(n);
			}

		}


		inline void NagataC0Triangle::get_derivatives_at(const double t[2], double dp_u[3], double dp_v[3])const
		{

			double u = t[0], v = t[1];
			double _2u = 2 * u, _2v = 2 * v;

			
			// p[i](u,v)= a_[0][0][i] + a_[1][0][i] * u + a_[0][1][i] * v
			//			+ a_[2][0][i] * u^2 + a_[1][1][i] * u*v  + a_[0][2][i] * v^2
			//			+ a_[2][1][i] * u^2*v + a_[1][2][i] * u*v^2;
			for (auto i = 0; i != 3; ++i)
			{
				dp_u[i] = a_[1][0][i] + _2u * a_[2][0][i] + a_[1][1][i] * v;
				dp_v[i] = a_[0][1][i] + _2v * a_[0][2][i] + a_[1][1][i] * u;
			}

		}




		inline void NagataC0Triangle::get_second_derivatives(const double t[2], double dd2[9])const
		{
			for (auto i = 0; i != 3; ++i)
			{
				double* p = dd2 + 3 * i;
				p[0] = 2 * a_[2][0][i] ;
				p[1] = 2 * a_[0][2][i] ;
				p[2] = a_[1][1][i] ;
			}

		}

		// versione con dominio dei parametri modificato
		inline int NagataC0Triangle::get_ray_intersection(const double r[3], double& alpha, double t[2], double parameter_bounds_tol, double  tolerance)const
		{

			int it=ray_intersection_Newton_( r, alpha, t, tolerance);
			// in caso di convergenza eseguo il test sui parametri
			
			if (it >= 0)
			{
				if ((t[0] > 1 + parameter_bounds_tol) ||
					(t[0] < -parameter_bounds_tol) ||
					(t[1] > 1- t[0] + parameter_bounds_tol) ||
					(t[1] < -parameter_bounds_tol))

					return -it;
			}
			return it;
		}


		// versione originale di Nagata
		//inline int NagataC0Triangle::get_ray_intersection(const double r[3], double& alpha, double t[2], double parameter_bounds_tol, double  tolerance)const
		//{

		//	int it = ray_intersection_Newton_(r, alpha, t, tolerance);
		//	// in caso di convergenza eseguo il test sui parametri

		//	if (it >= 0)
		//	{
		//		if ((t[0] > 1 + parameter_bounds_tol) ||
		//			(t[0] < -parameter_bounds_tol) ||
		//			(t[1] > t[0] + parameter_bounds_tol) ||
		//			(t[1] < -parameter_bounds_tol))

		//			return -it;
		//	}
		//	return it;
		//}


		inline int NagataC0Triangle::ray_intersection_Newton_( const double r[3], double& alpha, double t[2], double  tolerance)const
		{
			

			//DEBUG
			//  verifico che tipo di curva e' l'intersezione tra il piano span{ P0 ,r} (che contiene l'intersezione) e la superficie
			// determino la normale al piano


			Vector3 nu;
			double P0_norm = V3Op::norm2(P_[0]);
			V3Op::cross_product(r, P_[0], (double*)nu);
			double detB,  gamma[3][3];
			if (nu.norm()/P0_norm < tolerance)
			{
			// il punto intersezione coincide con il vertice p0
				t[0] = t[1] = 0;
				alpha = P0_norm / V3Op::norm2(r);
				return 0;
			}
			else
			{
				nu /= P0_norm;
				gamma[2][0] = V3Op::dot((double*)nu, a_[2][0]);
				gamma[0][2] = V3Op::dot((double*)nu, a_[0][2]);
				gamma[1][1] = V3Op::dot((double*)nu, a_[1][1]);
				detB = gamma[2][0] * gamma[0][2] - 0.25*gamma[1][1] * gamma[1][1];
				// se detB >0   ->  ellisse
				// se detB ==0	-> parabola
				// se detB <0	-> iperbole
			}



			int MAXIT = 50;
			// temps
			double a[3], b[3], c[3], d[3];

			// distanza di riferimento utile per vari test successivi
			V3Op::sub(P_[2], P_[0], d);
			double r_ref = V3Op::norm2(d);


			//=======================================================================================
			// comincio determinando una stima iniziale dell'intersezione basata su un'aprossimazione
			// affine della supercie del tipo
			//		r[](u,v)= a[] + b[] u + c[] v
			// in cui  i vettori a[],b[],c[] sono ottenuti per interpolazione dei vertici 

			V3Op::sub(P_[1], P_[0], b);
			V3Op::sub(P_[2], P_[1], c);
			V3Op::copy(P_[0], a);
			

			double k[3][3] = {
				{ -b[0], -c[0], r[0] },
				{ -b[1], -c[1], r[1] },
				{ -b[2], -c[2], r[2] } };
			double x[3];
			double det = MML3::Math::Matrix<3, double>::solve(k[0], a, x);


			// adesso ho la stima x[]=[x,y,alpha] dell'intersezione
			// verifico l'entita' dell'errore

			//Siccome  adesso a,b,c non servono piu'le riciclo
			// stima del punto calcolato con la stima delle coordinate
			double* r_est = c;
			// residuo
			double* res = a;
			// incremento delle variabili
			double* dx = b;

			// calcolo la stima del punto
			get_at(x, r_est, nullptr);
			// calcolo il residuo
			V3Op::lin_comb(1, r_est, -x[2], r, res);
			// errore relativo
			double err = V3Op::norm2(res) / r_ref;

			//=======================================================================================
			// adesso,se necessario, rifinisco la stima iterando con Newton sulla superficie quadratica reale
			int it = 0;
			for (; it != MAXIT && (err > tolerance); ++it)
			{
				// calcolo il gradiente  del residuo
				get_derivatives_at(x, k[0],k[1]);
				k[2][0] = -r[0];
				k[2][1] = -r[1];
				k[2][2] = -r[2];
				MML3::Math::Matrix<3, double>::transpose(k[0]);

				// determino l'incremento delle variabili
				det = MML3::Math::Matrix<3, double>::solve(k[0], res, dx);
				if (det == 0)
				{
					++it;
					break;
				}
				// aggiorno la stima delle variabili
				V3Op::axpy(-1.0, dx, x);
				// aggiorno la stima del punto sulla superficie
				get_at(x, r_est, nullptr);
				// aggiorno il residuo
				V3Op::lin_comb(1, r_est, -x[2], r, res);
				// aggiorno l'errore
				err = V3Op::norm2(res) / r_ref;
			}


			alpha = x[2];
			t[0] = x[0];
			t[1] = x[1];
			if (err > tolerance)
				return -it;
			else
				return it;

		}





		inline int NagataC0Triangle::ray_intersection_first_order_(const double r[3], double& alpha, double t[2], double  tolerance)const
		{


			int MAXIT = 50;

			typedef MML3::StaticArrayOp<2, double> V2Op;

			
			// determino il vertice più vicino al raggio
			double cos_max = 0;
			int idx_nearest_vertex = 0;

			for (int i = 0; i != 3; ++i)
			{
				double co = V3Op::dot(r, P_[i]) / V3Op::norm2(P_[i]);
				if (co > cos_max)
				{
					cos_max = co;
					idx_nearest_vertex = i;
				}
			}

			// coordinate parametriche  superficiali della stima di partenza del punto 
			t[0] = t[1] = 0;
			switch (idx_nearest_vertex)
			{
			case 0: break;
			case 1: t[0] = 1.0; break;
			case 2:	t[1] = 1; 	break;
			case 3: t[1] = 1;
			default:
				throw std::exception("");
			};

			// iterazioni
			//  stima del punto
			Point3 P;
			// e della normale 
			Vector3 n, residual;

			double der[2][3], dt[2];
			double k[2][2], b[2];


			// stima del punto e relativa normale sulla superficie

			get_at(t, (double*)P, (double*)n);
			alpha = V3Op::dot((double *)P, (double *)n) / V3Op::dot(r, (double *)n);

			int it = 0;
			for (; it != MAXIT; ++it)
			{
				// derivate sulla superficie
				get_derivatives_at(t, der[0], der[1]);
				V3Op::lin_comb(-1., (double*)P, alpha, r, (double*)residual);
				k[0][0] = V3Op::dot(der[0], der[0]);
				k[1][1] = V3Op::dot(der[1], der[1]);
				k[0][1] = k[1][0] = V3Op::dot(der[0], der[1]);
				b[0] = V3Op::dot((double*)residual, der[0]);
				b[1] = V3Op::dot((double*)residual, der[1]);

				double det = MML3::Math::Matrix<2, double>::solve(k[0], b, dt);

				t[0] += dt[0];
				t[1] += dt[1];
				get_at(t, (double*)P, (double*)n);
				alpha = V3Op::dot((double *)P, (double *)n) / V3Op::dot(r, (double *)n);

				if (V2Op::norm2(dt) < tolerance)
					break;
			}

			if (V2Op::norm2(dt) > tolerance)
				return -it;
			else
				return it;

		}

		







		//////////////////////////////////////////////////////////////////////////////////////////////
		// NAGATAC0_PATCH IMPL

		inline void  NagataC0_Patch::get_c_(const double* d, const double* n0, const double* n1, double* c)const
		{
			double co = n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2]; // coseno dell'angolo tra n0 ed n1
			
			double si2 = 1. - co*co; // sin^2 dell'angolo

			// il  test è positivo se l'angolo tra le due normali è inferiore ad 1/100 di grado sessagesimale
			if (si2 < 3E-8)
			{
				V3Op::zero(c);
				return;
			}


			double dn0 = V3Op::dot(d, n0);
			double dn1 = V3Op::dot(d, n1);
			double alpha_1 = -1 / si2 * ( dn1 + co*dn0 );
			double alpha_0 = dn0 - alpha_1 * co;


			V3Op::lin_comb(alpha_0,n0, alpha_1, n1, c);
			
		}

		inline int NagataC0_Patch::config(const double* P[4], const double* n[4])
		{
			if (P[3] != nullptr)
				return config_quad(P, n);
			else
				return config_tria(P, n);
	
		}

		inline int NagataC0_Patch::config_quad(const double* P[4], const double* n[4])
		{
			double d[4][3]; //vettori dei lati
			double c[4][3]; //vettori curvatura dei lati

			// i lati sono numerati da 0 a 3 in senso antiorario. Lid contiene per ogni lato
			// gli id del  nodo di inizio e quello di fine 
			static int Lid[4][2] = {
				{ 0, 1 },
				{ 1, 2 },
				{ 3, 2 },
				{ 0, 3 } };


			is_triangle_ = false;

			for (auto i = 0; i!=4; ++i)
			{
				int i0 =Lid[i][0], i1 =Lid[i][1];
				V3Op::sub(P[i1], P[i0], d[i]); 
				get_c_(d[i], n[i0], n[i1], c[i]);
			}

			
			

			V3Op::copy(P[0], a_[0][0]); //a_00=P_0
			V3Op::copy(c[0], a_[2][0]);	//a_20=c_0
			V3Op::copy(c[3], a_[0][2]);	//a_02=c_3

			V3Op::sub(d[0], c[0], a_[1][0]);	//a_10=d_0-c_0
			V3Op::sub(d[3], c[3], a_[0][1]);	//a_01=d_3-c_3
			V3Op::sub(c[1], c[3], a_[1][2]);	//a_12=c_1-c_3
			V3Op::sub(c[2], c[0], a_[2][1]);	//a_21=c_2-c_0

			// infine il piu' complicato a_11=(d_1-d_3) - (c_1-c_3) - (c_2-c_0)
			V3Op::sub(d[1], d[3], a_[1][1]);
			V3Op::sub(a_[1][1], a_[1][2], a_[1][1]);
			V3Op::sub(a_[1][1], a_[2][1], a_[1][1]);

			// NB: a_[2][2] not used

		
			return 0;
		}


		inline int NagataC0_Patch::config_tria(const double* P[3], const double* n[3])
		{
			double d[3][3]; //vettori dei lati
			double c[3][3]; //vettori curvatura dei lati

			// i lati sono numerati da 0 a 2 in senso antiorario. Lid contiene per ogni lato
			// gli id del  nodo di inizio e quello di fine 
			//
			static int Lid[3][2] = {
				{ 0, 1 },
				{ 1, 2 },
				{ 0, 2 } };


			is_triangle_ = true;

			for (auto i = 0; i != 3; ++i)
			{
				int i0 = Lid[i][0], i1 = Lid[i][1];
				V3Op::sub(P[i1], P[i0], d[i]);
				get_c_(d[i], n[i0], n[i1], c[i]);
			}


			V3Op::zero(a_[1][2]);
			V3Op::zero(a_[2][1]);


			V3Op::copy(P[0], a_[0][0]);	//a_00=P_0
			V3Op::copy(c[0], a_[2][0]);	//a_20=c_0
			V3Op::copy(c[2], a_[0][2]);	//a_02=c_2

			V3Op::sub(d[0], c[0], a_[1][0]);	//a_10=d_0-c_0

			V3Op::sub(d[2], c[2], a_[0][1]);	//a_01=d_2-c_2
			
				
			V3Op::sum(c[2], c[0], a_[1][1]);		// a_11=c_0+c_2-c_1
			V3Op::sub(a_[1][1], c[1], a_[1][1]); // ...
			

			// NB: a_[2][2] not used


			return 0;
		}


		inline void NagataC0_Patch::get_at(const double t[2], double r[3], double n[3])const
		{
			
			
			
			const double u =t[0], v = t[1];
			const double	uu = u*u,
							vv = v*v, 
							uv = u*v;
			const double	uuv = uu*v, 
							uvv = u*vv;

			if( r!=nullptr)
			{
				
				for (auto i = 0; i != 3; ++i)
				{
					r[i] =		a_[0][0][i] + a_[1][0][i] * u + a_[0][1][i] * v 
							+	a_[2][0][i] * uu + a_[1][1][i] * uv + a_[0][2][i] * vv 
							+	a_[2][1][i] * uuv + a_[1][2][i] * uvv;
				}
				
			}


			if (n!=nullptr)
			{
				
				double t0[3],t1[3]; // vettori coordinati nel punto

				for (auto i = 0; i != 3; ++i)
				{


					t0[i] = a_[1][0][i] + 2 * u*a_[2][0][i] + v*a_[1][1][i] + 2 * uv *a_[2][1][i] + vv* a_[1][2][i];
					t1[i] = a_[0][1][i] + 2 * v*a_[0][2][i] + u*a_[1][1][i] +  uu *a_[2][1][i] + 2*uv*a_[1][2][i] ;
				}
				V3Op::cross_product(t0, t1, n);
				V3Op::unitize(n);
				
			}
			

		}


		inline void NagataC0_Patch::get_derivatives_at(const double t[2], double d[6])const
		{

			double u = t[0], v =t[1];
			double _2u = 2 * u, _2v = 2 * v;

			double* dp_u = d;			// dp_u[i] := dp[i]/du
			double* dp_v = d + 3;		// dp_v[i] := dp[i]/dv
			

			// p[i](u,v)= a_[0][0][i] + a_[1][0][i] * u + a_[0][1][i] * v
			//			+ a_[2][0][i] * u^2 + a_[1][1][i] * u*v  + a_[0][2][i] * v^2
			//			+ a_[2][1][i] * u^2*v + a_[1][2][i] * u*v^2;
			for (auto i = 0; i != 3; ++i)
			{
				dp_u[i]	= a_[1][0][i] + _2u * a_[2][0][i] + (a_[1][1][i] + a_[2][1][i] * _2u	+ a_[1][2][i] * v	) * v;
				dp_v[i]	= a_[0][1][i] + _2v * a_[0][2][i] + (a_[1][1][i] + a_[2][1][i] * u		+ a_[1][2][i] * _2v	) * u;
			}
	
		}


		

		inline void NagataC0_Patch::get_second_derivatives(const double t[2], double dd2[9])const
		{
			for (auto i = 0; i != 3; ++i)
			{
				double* p = dd2 + 3 * i;
				p[0] = 2 * (a_[2][0][i] + a_[2][1][i]*t[1] );
				p[1] = 2 * (a_[0][2][i] + a_[1][2][i]*t[0]);
				p[2] = a_[1][1][i] + 2 * (a_[2][1][i] * t[0] + a_[1][2][i] * t[1]);
			}
			 
		}


		inline int NagataC0_Patch::get_ray_intersection(const double* p[4], const double r[3], double& alpha, double t[2], double  tolerance)const
		{

			return ray_intersection_Newton_(p, r, alpha, t, tolerance);
		}
		 

		inline int NagataC0_Patch::ray_intersection_Newton_(const double* p[4], const double r[3], double& alpha, double t[2], double  tolerance)const
		{

			
			int MAXIT = 50;
			// temps
			double a[3], b[3], c[3], d[3];
			
			// distanza di riferimento utile per vari test successivi
			V3Op::sub(p[2], p[0], d);
			double r_ref = V3Op::norm2(d);


			//=======================================================================================
			// comincio determinando una stima iniziale dell'intersezione basata su un'aprossimazione
			// affine della supercie del tipo
			//		r[](x,y)= a[] + b[] x + c[] y
			// in cui  i vettori a[],b[],c[] sono ottenuti per interpolazione dei vertici nel caso del
			// triangolo e con la tecnica dei minimi quadrati,  minimizzando cioe' l'errore
			//		err= sum_{i=1:4} ( r(x_i,y_i)- P_i)^2, 
			// nel caso del quadrilatero

			

			if (is_triangle_)
			{
				V3Op::sub(p[1], p[0], b);
				V3Op::sub(p[2], p[0], c);
				V3Op::copy(p[0], a);
			}
			else
			{
				const double* P[] = { p[0], p[1], p[2], p[3] };
				double mult_a[] = { 3. / 4., 0.25, -0.25, 0.25 };
				V3Op::lin_comb(4,  mult_a, P, a);
				double mult_b[] = { -0.5, 0.5, 0.5, -0.5 };
				V3Op::lin_comb(4, mult_b,P,  b);
				double mult_c[] = { -0.5, -0.5, 0.5, 0.5 };
				V3Op::lin_comb(4,  mult_c,P, c);
			}

			double k[3][3] = {
				{-b[0],- c[0], r[0] },
				{-b[1],- c[1], r[1] },
				{-b[2],- c[2], r[2] } };
			double x[3];
			double det=MML3::Math::Matrix<3, double>::solve(k[0], a, x);


						
			

			// adesso ho la stima x[]=[x,y,alpha] dell'intersezione
			// verifico l'entita' dell'errore
			
			//Siccome  adesso a,b,c non servono piu'le riciclo
			// stima del punto calcolato con la stima delle coordinate
			double* r_est=c;
			// residuo
			double* res = a;
			// incremento delle variabili
			double* dx = b;

			// calcolo la stima del punto
			get_at(x, r_est, nullptr);
			// calcolo il residuo
			V3Op::lin_comb( 1, r_est,  -x[2],r, res);
			// errore relativo
			double err = V3Op::norm2(res) / r_ref;
			
			//=======================================================================================
			// adesso,se necessario, rifinisco la stima iterando con Newton sulla superficie quadratica reale
			
			int it = 0;
			for( ; it != MAXIT && (err > tolerance); ++it)
			{
				
				
				// calcolo il gradiente  del residuo
				get_derivatives_at(x, k[0]);
				k[2][0] = -r[0];
				k[2][1] = -r[1];
				k[2][2] = -r[2];
				MML3::Math::Matrix<3, double>::transpose(k[0]);
				
				// determino l'incremento delle variabili
				det = MML3::Math::Matrix<3, double>::solve(k[0],res, dx);
				if (det == 0)
				{
					++it;
					break;
				}
				// aggiorno la stima delle variabili
				V3Op::axpy(-1.0, dx, x);
				// aggiorno la stima del punto sulla superficie
				get_at(x, r_est, nullptr);
				// aggiorno il residuo
				V3Op::lin_comb(1,r_est,-x[2], r, res);
				// aggiorno l'errore
				err = V3Op::norm2(res) / r_ref;
			}


			alpha = x[2];
			t[0] = x[0];
			t[1] = x[1];
			if (err > tolerance)
				return -it;
			else
				return it;

		}





		inline int NagataC0_Patch::ray_intersection_first_order_(const double* p[4], const double r[3], double& alpha, double t[2], double  tolerance)const
		{

			
			int MAXIT = 50;
		
			typedef MML3::StaticArrayOp<2, double> V2Op;
			
			int num_vertices = is_triangle_ ? 3 : 4;
			
			// determino il vertice più vicino al raggio
			double cos_max = 0;
			int idx_nearest_vertex = 0;
			
			for (int i = 0; i != num_vertices; ++i)
			{
				double co=V3Op::dot(r, p[i])/V3Op::norm2(p[i]);
				if (co > cos_max)
				{
					cos_max = co;
					idx_nearest_vertex = i;
				}
			}
			
			// coordinate parametriche  superficiali della stima di partenza del punto 
			t[0] = t[1] = 0;
			switch (idx_nearest_vertex)
				{
				case 0: break;
				case 1: t[0] = 1.0; break;
				case 2: 
					if (is_triangle_)
						t[1] = 1;
					else
						t[0] = t[1] = 1; 
					break;
				case 3: t[1] = 1;
				default:
					throw std::exception("");
				};

			// iterazioni
			//  stima del punto
			Point3 P;
			// e della normale 
			Vector3 n, residual;

			double der[2][3], dt[2];
			double k[2][2],b[2];
			

			// stima del punto e relativa normale sulla superficie
			
			get_at(t, (double*)P, (double*)n);
			alpha = V3Op::dot((double *)P, (double *)n) / V3Op::dot(r, (double *)n);
			
			int it = 0;
			for (; it != MAXIT; ++it)
			{
				// derivate sulla superficie
				get_derivatives_at(t, der[0]);
				V3Op::lin_comb(-1., (double*)P, alpha, r,(double*)residual);
				k[0][0] = V3Op::dot(der[0], der[0]);
				k[1][1] = V3Op::dot(der[1], der[1]);
				k[0][1] = k[1][0]=V3Op::dot(der[0], der[1]);
				b[0] = V3Op::dot((double*)residual, der[0]);
				b[1] = V3Op::dot((double*)residual, der[1]);

				double det=MML3::Math::Matrix<2, double>::solve(k[0], b, dt);

				t[0] += dt[0];
				t[1] += dt[1];
				get_at(t, (double*)P, (double*)n);
				alpha = V3Op::dot((double *)P, (double *)n) / V3Op::dot(r, (double *)n);

				if (V2Op::norm2(dt) < tolerance)
					break;
			}

			if (V2Op::norm2(dt) > tolerance)
				return -it;
			else
				return it;

		}



		template<int N>
		void CubicSplineCurve<N>::configure( const double p_[/*N*4*/], double tolerance)
		{
			

			// calcolo un parametro r con cui confrontare i coefficienti della cubica
			// per stabilire se devono essere nulli
			const double *p = p_;
			// calcolo la norma di P2-P1
			double r = 0.0;
			for (int i = 0; i != N; ++i)
			{
				r+=pow(p[2] - p[0], 2);
				p += 4;
			}
			r = sqrt(r);



			
			p = p_;
			double *C = c_;
			for (int a = 0; a != N; ++a)
			{
				// parte essenziale 
				p = p_ + a * 4;
				C = c_ + a * 4;
				C[0] = p[0];
				C[1] = p[1];
				C[2] = 3 * (p[2] - p[0]) - (2 * p[1] + p[3]);
				C[3] = 2 * (p[0] - p[2]) + (p[1] + p[3]);
				
				// cancellazione dei termini irrilevanti  in modo che vengano rispettati i vincoli
				// P(0)=P_0  e P(1)=P_1
				// tutto dipende dall' equivalenza delle due rappresentazioni 
				// P(t) = P_0 + (d-c) t  + (c+g) t^2   - g   t^3
				//      = C_0  + C_1  t   + C_2  t^2   + C_3 t^3
				// essendo d=P_1 - P_0


				// se il coefficiente cubico g=C_3 è irrilevante, oltre che annullare C_3
				// devo correggere C_2
				if (fabs(C[3]) / r < tolerance)
				{
					C[2] += C[3];
					C[3] = 0.;
				}
				// se il coefficiente quadratico C_2=c+g è trascurabile, porlo a zero significa assumere che c=-g , 
				// ovvero c= C_3 e questo implica che C_1 = d-C_3
				if (fabs(C[2]) / r < tolerance)
				{
					C[2] = 0.;
					C[1] = p[2] - p[0] - C[3];
				}

			}
		}



		template<int N>
		void CubicSplineCurve<N>::configure(double ta, double tb, const double p_[/*N*4*/])
		{
			x1_ = ta;
			h_ = tb - ta;

			const double *p = p_;
			double *c = c_;
			for (int a = 0; a != N; ++a)
			{
				p = p_ + a * 4;
				c = c_ + a * 4;
				c[0] = p[0];
				c[1] = h_* p[1];
				c[2] = 3 * (p[2] - p[0]) - h_*(2 * p[1] + p[3]);
				c[3] = 2 * (p[0] - p[2]) + h_*(p[1] + p[3]);
			}
		}


		template<int N>
		void CubicSplineCurve<N>::eval(double t, double f[N], double df[/*N*/], double ddf[/*N*/])const
		{
			t = (t - x1_) / h_;
			const double *c = c_;
			for (int a = 0; a != N; ++a)
			{
				c = c_ + a * 4;
				if (f)
					f[a] = c[0] + t*(c[1] + t*(c[2] + c[3] * t));
				if (df)
					df[a] = (c[1] + t*(2 * c[2] + 3 * t* c[3])) / h_;
				if (ddf)
					ddf[a] = (2 * c[2] + 6 * t*c[3]) / (h_*h_);
			}


		}

		/*inline void CubicSplines::set(double x1, double x2, const double p[4])
		{
			x1_ = x1;
			h_ = x2 - x1;
			double pp1 = p[1] * h_;
			double pp3 = p[3] * h_;

			c_[0] = p[0];
			c_[1] = h_* p[1];
			c_[2] = 3 * (p[2] - p[0]) - (2 * pp1 + pp3);
			c_[3] = 2 * (p[0] - p[2]) + (pp1 + pp3);
		}

		inline double	CubicSplines::value(double x)const
		{
			double t = (x - x1_) / h_;
			return c_[0] + t*(c_[1] + t*(c_[2] + c_[3] * t));
		}

		inline double	CubicSplines::derivative(double x)const
		{
			double t = (x - x1_) / h_;
			return (c_[1] + t*(2 * c_[2] + 3 * t* c_[3])) / h_;
		}

		inline double	CubicSplines::second_derivative(double x) const
		{
			double t = (x - x1_) / h_;
			return (2 * c_[2] + 6 * t*c_[3]) / (h_*h_);

		}

		inline void	CubicSplines::eval(double x, double& value, double& deriv, double& deriv2)const
		{
			double t = (x - x1_) / h_;
			value = c_[0] + t*(c_[1] + t*(c_[2] + c_[3] * t));
			deriv = (c_[1] + t*(2 * c_[2] + 3 * t* c_[3])) / h_;
			deriv2 = (2 * c_[2] + 6 * t*c_[3]) / (h_*h_);


		}*/



	
}//END MML


