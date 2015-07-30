#pragma once

namespace MML3{



	namespace Math
	{

		//===================================================================================================
		// Algoritmi  per determinare lo zero di una funzione F:R ->R nell'intervallo [x1, x2]
		// La funzione deve avere segno diverso nei due estremi, altrimenti viene tornato un codice di errore

		// RETVAL:
		// N>=0	:   Successo, N=numero di iterazioni eseguite per trovare la soluzione
		// -1   :   Errore,  il segno della funzione  nei due estremi è lo stesso
		// -iter:   Errore, esaurite le iterazioni consentite senza aver trovato la soluzione con la precisione richiesta



		//--------------------------------------------------------------
		// Algoritmo di bisezione con la possibilità di scegliere se la bisezione
		// avviene con la regola del medio o per interplazione lneare

		enum class BisectionType { MEAN, LINEAR };

		template<  typename T, typename FUNCTOR, BisectionType BT= BisectionType::LINEAR>
		int bisection_root(T x1, T x2, T& sol, const FUNCTOR& F,  T precision = 1e-6, int max_iter = 1000)
		{


			T y1 = F(x1);
			T y2 = F(x2);

			// la soluzione potrebbe essere uno dei due estremi
			if (fabs(y1) < precision)
			{
				sol = x1;
				return 1;
			}
			else if (fabs(y2) < precision)
			{
				sol = x2;
				return 1;
			}

			// l'algoritmo è efficente solo se parte da due punti x1 e x2 per i quali il segno della funzione è diverso
			// per cui ritorno un codice d'errore se non è così
			if (y1*y2 > 0)
				return -1;




			// algoritmo di bisezione con interpolazione lineare
			int j;
			for (j = 0; j != max_iter; ++j)
			{

				if (BT==BisectionType::LINEAR)
				{
					// interpolazione lineare
					sol = x1 - y1*(x2 - x1) / (y2 - y1);
				}
				else
				{
					// regola del medio
					sol = 0.5*(x1 + x2);
				}
				
				double y = F(sol);
				if (fabs(y) < precision)
					return j + 1;

				// se y ha lo stesso segno di y1, allora 
				// sol diventa il nuovo x1
				if (y * y1 > 0)
				{
					x1 = sol;
					y1 = y;
				}
				// altrimenti sol diventa il nuovo x2
				else
				{
					x2 = sol;
					y2 = y;
				}
			}

			return -j;

		}

		

		



		template<  typename FUNCTOR>
		int brent_root(double a, double b, double& root, const FUNCTOR& F, double t = 1e-6 /*toll*/, int max_it = 1000)

			//****************************************************************************80
			//
			//  Purpose:
			//
			//    brent_root seeks the root of a function F(X) in an interval [a,b].
			//
			//  Discussion:
			//
			//    The interval [a,b] must be a change of sign interval for F.
			//    That is, F(a) and F(b) must be of opposite signs.  Then
			//    assuming that F is continuous implies the existence of at least
			//    one value c between a and b for which F(c) = 0.
			//
			//    The location of the zero is determined to within an accuracy
			//    of 6 * MACHEPS * abs ( c ) + 2 * t.
			//
			//    Thanks to Thomas Secretin for pointing out a transcription error in the
			//    setting of the value of P, 11 February 2013.
			//
			//  Licensing:
			//
			//    This code is distributed under the GNU LGPL license.
			//
			//  Modified:
			//
			//    11 February 2013
			//
			//  Author:
			//
			//    Original FORTRAN77 version by Richard Brent.
			//    C++ version by John Burkardt.
			//
			//  Reference:
			//
			//    Richard Brent,
			//    Algorithms for Minimization Without Derivatives,
			//    Dover, 2002,
			//    ISBN: 0-486-41998-3,
			//    LC: QA402.5.B74.
			//
			//  Parameters:
			//
			//    Input, double a, b, the endpoints of the change of sign interval.
			//
			//    Input, double t, a positive error tolerance.
			//
			//    Input, FUNCTOR& F, the name of a user-supplied c++ functor
			//    whose zero is being sought.  The input and output
			//    of F() are of type double.
			//
			//    Input, max_it
			//
			//    Output, double root, the estimated value of a zero of
			//    the function F.
			//    Output return value: 
			//				the numer of iteration in case of success, 
			//				a negative value in case of failure:
			//					-1 if the function values at initial guesses a,b have the same sign
			//					-iter if convergence was not reached in the maximum prescribed iterations
			//
		{
			double c;
			double d;
			double e;
			double fa;
			double fb;
			double fc;
			double m;
			double macheps;
			double p;
			double q;
			double r;
			double s;
			double sa;
			double& sb = root;
			double tol;
			//
			//  Make local copies of A and B.
			//
			sa = a;
			sb = b;
			fa = F(sa);
			fb = F(sb);



			// la soluzione potrebbe essere uno dei due estremi
			if (fabs(fa) < t)
			{
				sb = a;
				return 1;
			}
			else if (fabs(fb) < t)
			{
				sb = b;
				return 1;
			}

			// l'algoritmo è efficente solo se parte da due punti a e b per i quali il segno della funzione è diverso
			// per cui ritorno un codice d'errore se non è così
			if (sign(fa) == sign(fb))
				return -1;

			c = sa;
			fc = fa;
			e = sb - sa;
			d = e;

			macheps = std::numeric_limits<double>::epsilon();
			int iter = 0;
			for (; iter != max_it; ++iter)
			{
				if (fabs(fc) < fabs(fb))
				{
					sa = sb;
					sb = c;
					c = sa;
					fa = fb;
					fb = fc;
					fc = fa;
				}

				tol = 2.0 * macheps * fabs(sb) + t;
				m = 0.5 * (c - sb);

				if (fabs(m) <= tol || fb == 0.0)
				{
					break;
				}

				if (fabs(e) < tol || fabs(fa) <= fabs(fb))
				{
					e = m;
					d = e;
				}
				else
				{
					s = fb / fa;

					if (sa == c)
					{
						p = 2.0 * m * s;
						q = 1.0 - s;
					}
					else
					{
						q = fa / fc;
						r = fb / fc;
						p = s * (2.0 * m * q * (q - r) - (sb - sa) * (r - 1.0));
						q = (q - 1.0) * (r - 1.0) * (s - 1.0);
					}

					if (0.0 < p)
					{
						q = -q;
					}
					else
					{
						p = -p;
					}

					s = e;
					e = d;

					if (2.0 * p < 3.0 * m * q - fabs(tol * q) &&
						p < fabs(0.5 * s * q))
					{
						d = p / q;
					}
					else
					{
						e = m;
						d = e;
					}
				}
				sa = sb;
				fa = fb;

				if (tol < fabs(d))
				{
					sb = sb + d;
				}
				else if (0.0 < m)
				{
					sb = sb + tol;
				}
				else
				{
					sb = sb - tol;
				}

				fb = F(sb);

				if ((0.0 < fb && 0.0 < fc) || (fb <= 0.0 && fc <= 0.0))
				{
					c = sa;
					fc = fa;
					e = sb - sa;
					d = e;
				}
			}
			if (iter == max_it)
				return -2;
			return iter;

		}

	} // end namespace Math

}// end namespace MML3