#pragma once
#include<cmath>

namespace MML3
{ 

namespace  Math
{

	///////////////////////////////////////////////////////////////////////////////////////
	// costanti caratteristiche 
	namespace Const
	{
		static double Pi = 3.141592653589793238462;
		static double E = 2.718281828459045235360;
	}

	using namespace Const;
	
	///////////////////////////////////////////////////////////////////////////////////////
	// funzioni caratteristiche
		//% sign of a real
		inline double sign(double x ){return x < 0.? -1. : 1.;}
		//% positive part of a real (Macauley brackets)
		inline double positive_part(double x ){return x > 0.? x : 0.;}
		//% heavisyde function 
		inline double heavyside(double x, double x0=0.0){return x >= x0 ? 1.0 : 0.;}

		// conversione di angoli da gradi a radianti e viceversa
		inline  double deg2rad(double deg){	return deg/180. * Pi;	}
		inline  double rad2deg(double rad){ return rad/Pi * 180.;}

		/// delta di Kronecker
		inline  double delta(int a, int b){return (a==b?1.:0.);}	
		/// delta di Ricci
		inline  double delta(int a, int b, int c)
		{
			static signed char D[3][3][3]=
			{
				{{0,0,0} , {0,0,1},{0,-1,0}},
				{{0,0,-1}, {0,0,0},{1,0,0}},
				{{0,1,0} , {-1,0,0},{0,0,0}},
			};
			return double(D[a-1][b-1][c-1]);
		}	

		// prodotto scalare tra vettori
		inline double scalar_product(const double a[], const double b[], size_t sz)
		{
			double sum = 0;
			for (size_t i = 0; i != sz; ++i)
				sum += a[i] * b[i];
			return sum;
		}
		// norma euclidea di un vettore
		inline double norm2(const double a[], size_t sz)
		{
			return sqrt(scalar_product(a,a,sz));
		}
                
                
               
                
                
		// prodotto vettore 3D
                template<typename T>
		inline void vector_product_3D(const T a[], const T b[], T c[])
		{
			c[0] = a[1] * b[2] - a[2] * b[1];
			c[1] = a[2] * b[0] - a[0] * b[2];
			c[2] = a[0] * b[1] - a[1] * b[0];
		}

		///////////////////////////////////////////////////////////////////////////////////////	
		// algoritmi 
		
		
		/**
        Trova le radici di un equazione di secondo grado nella forma
        a x^2 + b x + c=0
        Le radici sono $\frac{ - b +- \sqrt{b^2 - 4 a c} }{2 a}
        Se le radici sono reali (distine o coincidenti) x1 ed x2 in uscita
        contengono le radici.
        Se le radici sono complesse la funzione ritorna la parte immaginaria
        e le radici sono
        x1 + j retval,    x2 - j retval
        dove retval indica il valore di ritorno della funzione.
		L'utente pu� discriminare tra i due casi in base al valore di ritorno della funzione,
		se questop � 0.0 allora le radici sono reali
        */
		inline double poly2_root(double a, double b, double c, double& x1, double& x2)
		{

			double Delta=b*b - 4. * a * c;
            
            if(Delta < 0.)
            {
                x1=x2=- b / (2.0 * a);
            	return sqrt(-Delta)/(2.0 * a);
            }
            
			Delta=sqrt(Delta);
			x1 = ( -b - Delta)/( 2. * a);
			x2 = ( -b + Delta)/( 2. * a);
			return 0.0;
		}

		/** 
		Trova le radici reali dell'equazione 
		x^3  + a_2 x^2 + a_1 x + a_0 = 0

		Siccome una cubica o ha 3 radici reali oppure ne ha 1 reale e 2
		complesse coniugate, la soluzione � sempre esprimibile con 4 
		numeri reali z_1, z_2, z_3, Im:
		z_1
		z_2 + j Im
		z_3 - j Im

		Le tre parti reali sono parametri della funzione passati 
		per riferimento, la eventuale parte immaginaria viene tornata
		dalla funzione

		Torna la parte immaginaria delle sue radici complesse. 
		Se torna 0.0 tutte tre le radici sono reali.
		*/

		inline  double poly3_root(double a2, double  a1, double  a0,   double&  z1, double& z2, double&  z3 ) 
		{    
			double Q, R, D, S, T, Q3;
			double im, th;
		
			Q = (3.*a1 - a2*a2)/9.0;
			Q3=Q*Q*Q;
			R = (9.*a1*a2 - 27.*a0 - 2. * a2*a2*a2)/54.0;
			D = Q3 + R*R;									//polynomial discriminant

			if (D >= 0) //  complex or duplicate roots
			{ 
				
				S = sign(R + sqrt(D)) *  pow( fabs(R + sqrt(D) ), 1./3.);
				T = sign(R - sqrt(D)) *  pow( fabs(R - sqrt(D) ), 1./3.);
		        
				z1 = -a2/3 + (S + T);         //       % real root
				z2 = -a2/3 - (S + T)/2.;      //       % real part of complex root
				z3 = -a2/3 - (S + T)/2.;      //       % real part of complex root
				im = fabs(sqrt(3.)*(S - T)/2.); //     % complex part of root pair

			}
			else  // distinct real roots
			{
				th = acos(R/sqrt(- Q3)) ;
				z1 = 2.*sqrt(-Q) * cos(th/3.) - a2/3. ;
				z2 = 2.*sqrt(-Q) * cos((th + 2.* Pi )/3.) - a2/3. ;
				z3 = 2.*sqrt(-Q) * cos((th + 4.* Pi )/3.) - a2/3. ;
				im = 0.;

			}
			return im;// imaginary part
		}




		/// costruisce la matrice ortogonale di Rodriguez (IPMEC-NLFEM ) che ruota un generico vettore attorno a w
		/// dell'angolo |w| 
		/// G= 1 + sin(a)/a <w> + (1-cos(a))/a^2  <w>^2, con w=mult * v , a=|w| e <w> matrice emissimetrica tale 
		/// che <w> x = w ^ x per ogni vettore x 
		///@param v vettore di dimensione 3 
		///@param mult moltiplicatore di v
		///@param G matrice 3x3 output
		template <typename T>
		void rodriguez_rotation_matrix( T* G, const T* v, T mult=1.)
		{
				

				T w[3],tmp=0.;
				for(int i=0;i!=3;++i)
				{
					w[i]=v[i]*mult;
					tmp+=w[i]*w[i];
				}
				tmp=sqrt(tmp);

				T a=sin(tmp)/tmp ;
				T b=(1 - cos(tmp))/( tmp * tmp) ;
				
				
				G[0]=1.0          - b* ( w[1]*w[1] + w[2]*w[2] );
				G[1]=    -a*w[2]  + b* w[0]*w[1];
				G[2]=    +a*w[1]  + b* w[0]*w[2];
				G[3]=    +a*w[2]  + b* w[0]*w[1];
				G[4]=1.0          - b* ( w[0]*w[0] + w[2]*w[2] );
				G[5]=    -a*w[0]  + b* w[1]*w[2];
				G[6]=    +a*w[1]  + b* w[0]*w[2];
				G[7]=    +a*w[0]  + b* w[1]*w[2];
				G[8]=1.0          - b* ( w[1]*w[1] + w[0]*w[0] );
		}



		///////////////////////////////////////////////////////////////////////////////////////
		// Algoritmi per matrici di dimensione 2,3,4 
		template<int N, typename T>
		class Matrix;

		template<typename T>
		class Matrix<2,T>
		{
		public:
			/// calcola il determinante
			static  T   det(const T* pA)
			{
				return pA[0] * pA[3] - pA[1] * pA[2];
			}
			
			// Inverte la matrice e ne ritorna il determinante
			static  T   inv( T* pA)
			{
				const T de=det(pA);
				if(de==T(0))
					return T(0);
				T A11=pA[3]/de;
				pA[3]=pA[0]/de;
				pA[0]=A11;
				pA[1]/=(-de);
				pA[2]/=(-de);
				return de;
			}

			// risolve il sistema lineare A x=b e torna il determinante di A
			static T	solve(const T* A, const T* b, T* x)
			{
				const T de=det(A);
				if(de==T(0))
					return T(0);
				x[0]= ( A[3] * b[0] - A[1] * b[1])/de;
				x[1]= (-A[2] * b[0] + A[0] * b[1])/de;
				return det;
			}


		};


		template<typename T>
		class Matrix<3,T>
		{
		public:
			static T   det(const double* A)
			{

				return	A[4] * ( A[0]*A[8] - A[2]*A[6] )+
					A[1] * ( A[5]*A[6] - A[3]*A[8] )+
					A[7] * ( A[2]*A[3] - A[0]*A[5] );
			}


			static T   inv( double* A)
			{
					const T A11=A[0]; 
					const T A12=A[1];
					const T A13=A[2];
					const T A21=A[3]; 
					const T A22=A[4];
					const T A23=A[5];
					const T A31=A[6]; 
					const T A32=A[7];
					const T A33=A[8];

					const T de= det(A);
					if(de==T(0))
						return T(0);
					
					A[0]= (A22*A33 - A23*A32)/de;
					A[1]=-(A12*A33 - A13*A32)/de;
					A[2]= (A12*A23 - A13*A22)/de;
					A[3]=-(A21*A33 - A23*A31)/de;
					A[4]= (A11*A33 - A13*A31)/de;
					A[5]=-(A11*A23 - A13*A21)/de;
					A[6]= (A21*A32 - A22*A31)/de;
					A[7]=-(A11*A32 - A12*A31)/de;
					A[8]= (A11*A22 - A12*A21)/de;
					return de;

				}


				static T solve(const T* A, const T* b, T* x)
				{

					const T A11=A[0]; 
					const T A12=A[1];
					const T A13=A[2];
					const T A21=A[3]; 
					const T A22=A[4];
					const T A23=A[5];
					const T A31=A[6]; 
					const T A32=A[7];
					const T A33=A[8];

					const T de= det(A);
					if(de==T(0))
						return T(0);
					
					x[0]=(   (A22*A33 - A23*A32)*b[0] - (A12*A33 - A13*A32)*b[1] +(A12*A23 - A13*A22)*b[2] )/de;
					x[1]=(  -(A21*A33 - A23*A31)*b[0] + (A11*A33 - A13*A31)*b[1] -(A11*A23 - A13*A21)*b[2] )/de;
					x[2]=(   (A21*A32 - A22*A31)*b[0] - (A11*A32 - A12*A31)*b[1] +(A11*A22 - A12*A21)*b[2] )/de;

					return de;
				}

		}; // end class



		template<typename T>
		class Matrix<4,T>
		{
		public:
			static T det(const T* A)
			{
				const T A11=A[0]; 
				const T A12=A[1];
				const T A13=A[2];
				const T A14=A[3];
				const T A21=A[4]; 
				const T A22=A[5];
				const T A23=A[6];
				const T A24=A[7];
				const T A31=A[8]; 
				const T A32=A[9];
				const T A33=A[10];
				const T A34=A[11];
				const T A41=A[12]; 
				const T A42=A[13];
				const T A43=A[14];
				const T A44=A[15];
				

				return 	A14*(A23*A32*A41-A22*A33*A41-A23*A31*A42+A21*A33*A42+A22*A31*A43-A21*A32*A43)+
					A13*(-A24*A32*A41+A22*A34*A41+A24*A31*A42-A21*A34*A42-A22*A31*A44+A21*A32*A44)+
					A12*(A24*A33*A41-A23*A34*A41-A24*A31*A43+A21*A34*A43+A23*A31*A44-A21*A33*A44)+
					A11*(-A24*A33*A42+A23*A34*A42+A24*A32*A43-A22*A34*A43-A23*A32*A44+A22*A33*A44);


			};


			static T   inv( T* A)
			{
				const T A11=A[0]; 
				const T A12=A[1];
				const T A13=A[2];
				const T A14=A[3];
				const T A21=A[4]; 
				const T A22=A[5];
				const T A23=A[6];
				const T A24=A[7];
				const T A31=A[8]; 
				const T A32=A[9];
				const T A33=A[10];
				const T A34=A[11];
				const T A41=A[12]; 
				const T A42=A[13];
				const T A43=A[14];
				const T A44=A[15];

				const T de=det(A);
				if(de!=T(0))
				{
					A[0]=(-A24*A33*A42 +A23*A34*A42 +A24*A32*A43 -A22*A34*A43 -A23*A32*A44 +A22*A33*A44)/de;
					A[1]=( A14*A33*A42 -A13*A34*A42 -A14*A32*A43 +A12*A34*A43 +A13*A32*A44 -A12*A33*A44)/de;
					A[2]=(-A14*A23*A42 +A13*A24*A42 +A14*A22*A43 -A12*A24*A43 -A13*A22*A44 +A12*A23*A44)/de;
					A[3]=( A14*A23*A32 -A13*A24*A32 -A14*A22*A33 +A12*A24*A33 +A13*A22*A34 -A12*A23*A34)/de;

					A[4]=( A24*A33*A41 -A23*A34*A41 -A24*A31*A43 +A21*A34*A43 +A23*A31*A44 -A21*A33*A44)/de;
					A[5]=(-A14*A33*A41 +A13*A34*A41 +A14*A31*A43 -A11*A34*A43 -A13*A31*A44 +A11*A33*A44)/de;
					A[6]=( A14*A23*A41 -A13*A24*A41 -A14*A21*A43 +A11*A24*A43 +A13*A21*A44 -A11*A23*A44)/de;
					A[7]=(-A14*A23*A31 +A13*A24*A31 +A14*A21*A33 -A11*A24*A33 -A13*A21*A34 +A11*A23*A34)/de;

					A[8]=(-A24*A32*A41 +A22*A34*A41 +A24*A31*A42 -A21*A34*A42 -A22*A31*A44 +A21*A32*A44)/de;
					A[9]=( A14*A32*A41 -A12*A34*A41 -A14*A31*A42 +A11*A34*A42 +A12*A31*A44 -A11*A32*A44)/de;
					A[10]=(-A14*A22*A41 +A12*A24*A41 +A14*A21*A42 -A11*A24*A42 -A12*A21*A44 +A11*A22*A44)/de;
					A[11]=( A14*A22*A31 -A12*A24*A31 -A14*A21*A32 +A11*A24*A32 +A12*A21*A34 -A11*A22*A34)/de;

					A[12]=( A23*A32*A41 -A22*A33*A41 -A23*A31*A42 +A21*A33*A42 +A22*A31*A43 -A21*A32*A43)/de;
					A[13]=(-A13*A32*A41 +A12*A33*A41 +A13*A31*A42 -A11*A33*A42 -A12*A31*A43 +A11*A32*A43)/de;
					A[14]=( A13*A22*A41 -A12*A23*A41 -A13*A21*A42 +A11*A23*A42 +A12*A21*A43 -A11*A22*A43)/de;
					A[15]=(-A13*A22*A31 +A12*A23*A31 +A13*A21*A32 -A11*A23*A32 -A12*A21*A33 +A11*A22*A33)/de;
				}
				return de;
			}

			static T   solve(const T* A, const T* b, T* x)
			{
				const T A11=A[0]; 
				const T A12=A[1];
				const T A13=A[2];
				const T A14=A[3];
				const T A21=A[4]; 
				const T A22=A[5];
				const T A23=A[6];
				const T A24=A[7];
				const T A31=A[8]; 
				const T A32=A[9];
				const T A33=A[10];
				const T A34=A[11];
				const T A41=A[12]; 
				const T A42=A[13];
				const T A43=A[14];
				const T A44=A[15];

				const T de=det(A);
				if(de!=T(0))
				{
					x[0]=(	(-A24*A33*A42 +A23*A34*A42 +A24*A32*A43 -A22*A34*A43 -A23*A32*A44 +A22*A33*A44)*b[0] +
							( A14*A33*A42 -A13*A34*A42 -A14*A32*A43 +A12*A34*A43 +A13*A32*A44 -A12*A33*A44)*b[1] +
							(-A14*A23*A42 +A13*A24*A42 +A14*A22*A43 -A12*A24*A43 -A13*A22*A44 +A12*A23*A44)*b[2] +
							( A14*A23*A32 -A13*A24*A32 -A14*A22*A33 +A12*A24*A33 +A13*A22*A34 -A12*A23*A34)*b[3] )/de;

					x[1]=(	( A24*A33*A41 -A23*A34*A41 -A24*A31*A43 +A21*A34*A43 +A23*A31*A44 -A21*A33*A44)*b[0]+
							(-A14*A33*A41 +A13*A34*A41 +A14*A31*A43 -A11*A34*A43 -A13*A31*A44 +A11*A33*A44)*b[1]+
							( A14*A23*A41 -A13*A24*A41 -A14*A21*A43 +A11*A24*A43 +A13*A21*A44 -A11*A23*A44)*b[2]+
							(-A14*A23*A31 +A13*A24*A31 +A14*A21*A33 -A11*A24*A33 -A13*A21*A34 +A11*A23*A34)*b[3] )/de;

					x[2]=(	(-A24*A32*A41 +A22*A34*A41 +A24*A31*A42 -A21*A34*A42 -A22*A31*A44 +A21*A32*A44)*b[0]+
							( A14*A32*A41 -A12*A34*A41 -A14*A31*A42 +A11*A34*A42 +A12*A31*A44 -A11*A32*A44)*b[1]+
							(-A14*A22*A41 +A12*A24*A41 +A14*A21*A42 -A11*A24*A42 -A12*A21*A44 +A11*A22*A44)*b[2]+
							( A14*A22*A31 -A12*A24*A31 -A14*A21*A32 +A11*A24*A32 +A12*A21*A34 -A11*A22*A34)*b[3] )/de;

					x[3]=(	( A23*A32*A41 -A22*A33*A41 -A23*A31*A42 +A21*A33*A42 +A22*A31*A43 -A21*A32*A43)*b[0]+
							(-A13*A32*A41 +A12*A33*A41 +A13*A31*A42 -A11*A33*A42 -A12*A31*A43 +A11*A32*A43)*b[1]+
							( A13*A22*A41 -A12*A23*A41 -A13*A21*A42 +A11*A23*A42 +A12*A21*A43 -A11*A22*A43)*b[2]+
							(-A13*A22*A31 +A12*A23*A31 +A13*A21*A32 -A11*A23*A32 -A12*A21*A33 +A11*A22*A33)*b[3] )/de;
				}

				return de;
			}

		}; // end class


		




}  // end namespace Math

}  // end namespace MML
    


