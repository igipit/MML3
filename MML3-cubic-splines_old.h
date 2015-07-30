#ifndef __MML_C_SPLINES_H__
#define __MML_C_SPLINES_H__

#include<cmath>

namespace MML3
{ 

	class CubicSplines
	{
	public:
		CubicSplines() :x1_(0), h_(1.){}
		void			set(double x1, double x2, const double p[4]);
		inline double	value(double x)const;
		inline double	derivative(double x)const;
		inline double	second_derivative(double x) const;
		inline	void	eval(double x, double& value, double& deriv, double& deriv2)const;
	private:

		
		double x1_, h_;
		double c_[4];
	};


	template<int N>
	class CubicSplineCurve
	{
	public:
		CubicSplineCurve() :x1_(0), h_(1.){}
		// qui p è un array Nx4 che nella riga a contiene i 4 valori della componente a-esima [p_a(t_1),p'_a(t_1),p_a(t_2),p'_a(t_2)]
		void	set(double t1, double t2, const double p[/*N*4*/]);

		void	eval(double t, double value[N], double deriv[/*N*/], double deriv2[/*N*/])const;

		const double* cx()const{ return c_; }
		const double* cy()const{ assert(N > 1); return c_ + 4; }
		const double* cz()const{ assert(N > 2); return c_ + 8; }

	private:
		double x1_, h_;
		double c_[N*4];
	};


	template<int N>
	void CubicSplineCurve<N>::set(double x1, double x2, const double p_[/*N*4*/])
	{
		x1_ = x1;
		h_ = x2 - x1;

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
	void CubicSplineCurve<N>::eval(double x, double value[N], double deriv[/*N*/], double deriv2[/*N*/])const
	{
		double t = (x - x1_) / h_;
		const double *c = c_;
		for (int a = 0; a != N; ++a)
		{
			c = c_ + a * 4;
			if (value)
				value[a] = c[0] + t*(c[1] + t*(c[2] + c[3] * t));
			if (deriv)
				deriv[a] = (c[1] + t*(2 * c[2] + 3 * t* c[3])) / h_;
			if (deriv2)
				deriv2[a] = (2 * c[2] + 6 * t*c[3]) / (h_*h_);
		}


	}

	inline void CubicSplines::set(double x1, double x2, const double p[4])
	{
		x1_ = x1; 
		h_= x2-x1;
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
		return (c_[1] + t*( 2* c_[2] + 3*t* c_[3]))/h_;
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


	}



	class Splines3Basis
	{
	public:
		/// Costruttore splines cubiche sul segmento [x1,x2]
		Splines3Basis(double x1, double x2) :x1_(x1), x2_(x2){}
		Splines3Basis() :x1_(), x2_(){}

		void set(double x1, double x2)
		{
			x1_ = x1;
			x2_ = x2;
		}

				

		
		// calcola il valore delle 4 splines di base sul dominio [x1_,x2_], delle loro derivate prime e seconde nel punto di coordinata x
		void get(double x, double* y, double* dy, double* ddy)const
		{
			const double h = x2_ - x1_;
			
			const double eta = (x - x1_) / h;


			get_standard(eta, y, dy, ddy);
			

			// valori
			if (y)
			{
				y[1] *= h;
				y[3] *= h;
			}

			if (dy)
			{
				// derivate
				dy[0] /= h;
				dy[2] = -dy[0];
			}

			if (ddy)
			{
				ddy[0] /=(h*h);
				ddy[1] /=h;
				ddy[2] = -ddy[0];
				ddy[3] /=h;
			}

		}

		// calcola il valore delle 4 splines di base definite sul dominio unitario [0,1], delle loro derivate prime e seconde nel punto di coordinata t
		static void  get_standard(double eta, double* y, double* dy, double* ddy)
		{
			double tmp;
			// valori
			if (y)
			{
				tmp = eta*eta;
				y[0] = 1.0 + tmp*(-3. + 2. * eta);
				y[1] = eta * (1. + eta*(-2. + eta));
				y[2] = 1. - y[0];
				y[3] = tmp*(eta - 1.);
			}

			if (dy)
			{

				// derivate
				dy[0] = (6.* eta) * (-1. + eta);
				dy[1] = 1.0 + eta*(-4. + 3. * eta);
				dy[2] = -dy[0];
				dy[3] =eta*(-2. + 3. * eta);
			}

			if (ddy)
			{
				tmp = 6.*eta;
				ddy[0] = 12.* eta - 6.0;
				ddy[1] =tmp- 4.0;
				ddy[2] = -ddy[0];
				ddy[3] = tmp - 2.0;
			}

		}

	private:
		double x1_, x2_;

	};


	


}  // end namespace MML
    





#endif