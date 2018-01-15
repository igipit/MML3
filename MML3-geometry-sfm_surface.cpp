#pragma once
#include<algorithm>
#include<vector>
#include<map>
#include<assert.h>
#include<fstream>
#include"MML3-math.h"
#include"MML3-geometry.h"
#include"MML3-interpolation.h"
#include"MML3-geometry-sfm_surface.h"
#include<iostream>


namespace MML3
{

	namespace GEO{


		


		// Calcola le coordinate geografiche di un punto r[], con la latitudine nel range [-pi/2, + pi/2] 
		// e la longitudine nel range [0, 2 pi)
		// ritorna 0 se il punto r[] e' identicamente nullo
		int SphericalMappedSurface::get_geo_coord(const double * r, double& latitude, double& longitude)const
		{
			double R=MML3::StaticArrayOp<3, double>::norm2(r);
			if (R == 0.)
				return 0;
			latitude = asin(r[2] / R);
			double cos_lat = cos(latitude);

			// se il punto specificato sta' su un polo la longitudine e' indetermnata
			// e la pongo uguale a zero
			if (cos_lat < 1e-12)
				longitude = 0;
			else
			{
				longitude = atan2(r[1], r[0]);

				// l'immagine di atan2 e' [-pi, +pi] e pertanto, se necessario traslo
				if (longitude < 0)
					longitude +=  MML3::Math::two_Pi;

				// se la longitudine e' 2*pi la pongo a zero
				if (fabs(longitude - MML3::Math::two_Pi) < 1e-12)
					longitude = 0.;
			}
			return 1;
		}



		// ritorna l'identificatore del  parallelo con latitudne <= a quella specificata
		int SphericalMappedSurface::get_base_parallel(double latitude) const
		{

			// ricerca sequenziale; se i punti fossero tanti, essendo ordinati converrebbe 
			// un algoritmo diverso
			if (latitude < -MML3::Math::half_Pi || latitude> MML3::Math::half_Pi)
				throw std::out_of_range("SphericalMappedSurface::get_base_parallel");
			// se siamo al polo nord il parallelo di base e' il penultimo
			if (MML3::Math::half_Pi - latitude < 1e-6)
				return GLatitude_.size() - 2;

			for (auto i = 0; i != GLatitude_.size(); ++i)
				if (GLatitude_[i + 1] > latitude)
					return i;
			return -1;
		}






		
		int SphericalMappedSurface::get_index(int meridian, int parallel)const
		{
			// se necessario normalizzo il meridiano
			if (meridian < 0)
				meridian += meridians();
			else if (meridian >= int(meridians()))
				meridian = meridian % meridians();

			// la normalizzazione del parallelo è  più complessa ed implica anche una traslazione di pi equilvalenti del meridiano

			if (parallel < 0)
			{
				parallel = -parallel;
				meridian = (meridian + meridians() / 2) % meridians();
			}
			else if (parallel >= int(parallels()))
			{
				parallel = (parallels() - 2) - parallel%parallels();
				meridian = (meridian + meridians() / 2) % meridians();
			}
		
			return 		meridian*parallels() + parallel;

		}

		 void SphericalMappedSurface::get_meridian_and_parallel_id(size_t i, int& meridian, int& parallel)const
		{
			meridian = i / parallels();
			parallel = i - meridian*parallels();

		}

		void  SphericalMappedSurface::get_quad_based_at_grid(int meridian, int parallel, Quad& quad)const
		{
			// non c'è una quad con base al polo nord
			assert(parallel >= 0 && parallel < int(parallels()) - 1);

			quad.set_pole(Quad::INTERNAL);
			
			if (parallel == 0)
				quad.set_pole(Quad::SOUTH_POLE);
			else if (parallel == parallels() - 2)
				quad.set_pole(Quad::NORTH_POLE);
			


			quad[0] = get_index(meridian, parallel);
			quad[1] = get_index(meridian + 1, parallel);
			quad[2] = quad[1] + 1;
			quad[3] = quad[0] + 1;
		}


		void  SphericalMappedSurface::get_enclosing_quad(double latitude, double longitude, Quad& quad)const
		{
			 get_quad_based_at_grid(get_base_meridian(longitude), get_base_parallel(latitude), quad);
		}

		




		double SphericalMappedSurface::get_distance(const double* r, double tolerance)const
		{
	
			double alpha = evaluate_at_ray_intersection(r, nullptr,tolerance);
			return 1 / alpha - 1;


		}





		 double SphericalMappedSurface::evaluate_at_ray_intersection(const double* r,  double * n, double tolerance)const
		{

			const int MAX_PATCH_ITER = 9;
			int near_quad_strides[9][2] = {
				 { 0, 0 },		{ -1, 0 }, { +1, 0 } ,
				 { -1, -1 },	{ 0, -1 }, { 1, -1 } ,
				 { -1, +1 },	{ 0, +1 }, { 1, +1 }  
			};


			

			double upper_t_lim = 1.5;
			double lower_t_lim = -0.5;

			double latitude, longitude;

			// la direzione del raggio non può essere identicamente nulla 
			if (Vec3Op::norm2(r) < tolerance)
				return 0;

			// calcola le coordinate geografiche sferiche del raggio  r[]
			get_geo_coord(r, latitude, longitude);
			

			// se ci troviamo ai poli ..
			if (cos(latitude) < 1E-9)
			{
				int idx = (latitude > 0) ? parallels() - 1 : 0;
				if (n != nullptr)
				{
					n[0] = n[1] = 0;
					n[2] = (idx == 0 )? -1. : +1.;
				}
				return pnt_[idx][2]/r[2];
			}
			
		
			Quad quad;
			int base_meridian = get_base_meridian(longitude);
			int base_parallel = get_base_parallel(latitude);
			bool intersection_found = false;

			/*int trial_meridian = base_meridian;
			int trial_parallel = base_parallel;*/

			double t[2], alpha;
			int RV;


			NagataC0Triangle NagataInterpolator;
			Triangle  T3[2];

			
			for (int patch_iter = 0; patch_iter != MAX_PATCH_ITER; ++patch_iter)
			{


				// i meridiani si mettono a osto da soli?
				int trial_meridian = base_meridian + near_quad_strides[patch_iter][0];
				

				int trial_parallel = base_parallel + near_quad_strides[patch_iter][1];

				if (trial_parallel < 0 || trial_parallel >= int(parallels()) - 1)
					continue;


				get_quad_based_at_grid(trial_meridian, trial_parallel, quad);

				
				quad.get_triangles(T3[0], T3[1]);

				int n_tria = quad.at_pole() ? 1 : 2;
				RV = -1;
				for (int tria = 0; tria != n_tria && RV < 0; ++tria)
				{
					int I = T3[tria][0],
						J = T3[tria][1],
						K = T3[tria][2];
					const double* P[] = { pnt_[I].data(), pnt_[J].data(), pnt_[K].data() };
					const double* N[] = { normal_[I].data(), normal_[J].data(), normal_[K].data() };
					NagataInterpolator.config(P, N);
					// questa può fallire se l'intersezione cade in realtà su una patch adiacente
					RV = NagataInterpolator.get_ray_intersection(r, alpha, t, 0.01, tolerance);

					if (RV >= 0)
					{
						search_iterations_count_ += RV;
						if (n != nullptr)
							NagataInterpolator.get_at(t, nullptr, n);
						double nr = Vec3Op::dot(n, r);
						if (nr <= 0.)
						{
							std::cout << "opsss... " << std::endl;
						}
						return alpha;
					}
					else
						search_iterations_count_ += -RV;
				}
				// devo cercare sugli elementi vicini


				
			}

			if ( RV <0)
				throw std::exception("SphericalMappedSurface::evaluate_at_ray_intersection: failure");
			
			return 0;
		}




		


		int SphericalMappedSurface::create(int n_meridians, const std::vector<double>& latitude, const std::vector<double>& point, const std::vector<double>& normal )
		{

			if (n_meridians*latitude.size()*3 != point.size() )
				return -1;

			if (normal.size() != 0 && point.size() != normal.size())
				return -2;

			meridians_ = n_meridians;
			meridian_step_ = MML3::Math::two_Pi / meridians_;
			GLatitude_=latitude;
			pnt_.resize(meridians_ * GLatitude_.size());
			normal_.resize(meridians_ * GLatitude_.size());
			const double* p = point.data();
			for (auto& el : pnt_)
			{
				el.copy(p);
				p += 3;
			}
			if (normal.size())
			{
				p = normal.data();
				for (auto& el : normal_)
				{
					el.copy(p);
					p += 3;
				}
			}
			else
				evaluate_normals_();
				

			// un test di cnvessità non sarebbe male
			return 0;
		}








		

		 void  SphericalMappedSurface::evaluate_normals_()
		{
			Point zero{ 0, 0, 0 };
			normal_.resize(pnt_.size());

			Vector a, b, c, d, tmp;
			
			for (int mer = 0; mer != meridians_; ++mer)
			{
				// polo sud
				normal_[get_index(mer, 0)] = { 0, 0, -1. };
				for (int par = 1; par != GLatitude_.size() - 1; ++par)
				{

					Vector& N = normal_[get_index(mer, par)];

					const Point& P = grid_point(mer , par);
	
					a = grid_point(mer + 1, par) - P;
					b = grid_point(mer , par+1)	- P;
					c = grid_point(mer - 1, par) - P;
					d = grid_point(mer , par-1)	- P;

					MML3::GEO::cross_product(a, b, N);

					MML3::GEO::cross_product(b, c, tmp);
					N += tmp;
					MML3::GEO::cross_product(c, d, tmp);
					N += tmp;
					MML3::GEO::cross_product(d, a, tmp);
					N += tmp;
					N.unitize();
				}
				// polo nord
				normal_[get_index(mer, GLatitude_.size() - 1)] = { 0, 0, +1. };
			}
			

			

		}



		


		void SphericalMappedSurface::write_vtk_polydata(std::ofstream& os)const
		{

			os << "# vtk DataFile Version 2.0\n"
				<< "created by SphericalMappedSurface::write_vtk_unstructured\n"
				<< "ASCII\n"
				<< "DATASET POLYDATA\n";

			os << "POINTS " << pnt_.size() << " DOUBLE\n";
			for (auto el : pnt_)
				os << el[0] << " " << el[1] << " " << el[2] << std::endl;
			
			int cells = meridians_ * (parallels() - 1);
			os << "POLYGONS " << cells << " " << cells *5 << std::endl;
			int idx[4];
			for (int mer = 0; mer != meridians(); ++mer)
			{

				idx[0] = get_index(mer, 0);
				idx[1] = get_index(mer + 1,0);
				
				for (int par = 0; par != parallels()-1; ++par)
				{
					idx[2] = idx[1] + 1;
					idx[3] = idx[0] + 1;
					os << 4 << " " << idx[0] << " " << idx[1] << " " << idx[2] << " " << idx[3] << std::endl;
					idx[0]++;
					idx[1]++;
				}

				idx[0]++;
				idx[1]++;
			}
	
			os << std::endl
				<< "POINT_DATA " << pnt_.size() << "\n"
				<< "NORMALS normal_unit_vectors double\n";
			for (auto el : normal_)
				os << el[0] << " " << el[1] << " " << el[2] << std::endl;
			os << std::endl;

		}






		 void write_vtk_polydata(std::ofstream& os, const SphericalMappedSurface& surf, size_t latitude_intervals, size_t longitude_intervals)
		{

			// genero l'array di punti e normali su una griglia (latitude_intervals+1)*(longitude_intervals +1)

			
			size_t num_points = (latitude_intervals + 1)*(longitude_intervals + 1);
			std::vector< SphericalMappedSurface::Point> point(num_points), normal(num_points);

			double Dangle[2] = { MML3::Math::two_Pi / longitude_intervals, MML3::Math::Pi / latitude_intervals };
			point.reserve(num_points);
			normal.reserve(num_points);
			double latitude, longitude = 0;


			SphericalMappedSurface::Point  P;
			size_t count=0;
			for (size_t mer = 0; mer < (longitude_intervals+1); ++mer)
			{
				longitude = mer*Dangle[0];
				latitude = -MML3::Math::half_Pi;
				for (size_t par = 0; par < latitude_intervals+1; ++par)
				{
					set_spherical(P, 1., latitude, longitude);
					double	alpha = surf.evaluate_at_ray_intersection((double*)P, (double*)normal[count], 1e-6);
					point[count++] = (P*=alpha);
					latitude += Dangle[1];
				}
			}



			os << "# vtk DataFile Version 2.0\n"
				<< "created by SphericalMappedSurface::write_vtk_unstructured\n"
				<< "ASCII\n"
				<< "DATASET POLYDATA\n";

			os << "POINTS " << point.size() << " DOUBLE\n";
			for (auto el :point)
				os << el[0] << " " << el[1] << " " << el[2] << std::endl;

			int cells = latitude_intervals * longitude_intervals;
			os << "POLYGONS " << cells << " " << cells * 5 << std::endl;
			// inizializzo solo le prime due componenti in quanto le ultime due vengono
			// calcolate aggiungendo 1 alle precedenti
			size_t idx[4] = { 0, latitude_intervals + 1 };
			for (int mer = 0; mer != longitude_intervals; ++mer)
			{
				for (int par = 0; par != latitude_intervals; ++par)
				{
					
					idx[2] = idx[1] + 1;
					idx[3] = idx[0] + 1;
					os << 4 << " " << idx[0] << " " << idx[1] << " " << idx[2] << " " << idx[3] << std::endl;
					idx[0]++;
					idx[1]++;
				}
				idx[0]++;
				idx[1]++;
			}

			os << std::endl
				<< "POINT_DATA " << normal.size() << "\n"
				<< "NORMALS normal_unit_vectors double\n";
			for (auto el : normal)
				os << el[0] << " " << el[1] << " " << el[2] << std::endl;
			os << std::endl;

		}

		



 bool  CircleMapped2DCurve::set_by_polygonal(const std::vector<REAL>& x, const std::vector<REAL>& y)
{

	if (x.size() < 3 || x.size() != y.size())
		return false;

	pnt_.clear();
	// memorizzo i punti in coordinate polari in una mappa ordinandoli per angoli crescenti
	std::map<REAL, REAL> vertex;
	for (size_t i = 0; i != x.size(); ++i)
	{
		REAL	X = x[i],
				Y = y[i];
		REAL    theta = angle(X, Y);
		REAL    r = sqrt(X*X + Y*Y);
		vertex[theta] = r;
	}


	
	// verifico la convessità
	REAL x1[2], x2[2], t[2]/*, theta0*/;
	
	if (verify_convexity(vertex) != 1)
	{
		return false;
	}
	

	// per ogni segmento definito da una coppia di vertici adiacenti calcolo 
	// il punto medio,
	// la tangente al lato
	// ds/dt che assumo uguale alla lunghezza del segmento diviso 1
	cpoint pt;

	auto it = vertex.begin();
	double r = it->second;
	double theta = it->first;
	x1[0] = r * cos(theta);
	x1[1] = r * sin(theta);
	bool end = false;
	for (; !end;)
	{
		++it;
		if (it == vertex.end())
		{
			it = vertex.begin();
			end = true;
		}

		r = it->second;
		theta = it->first;
		x2[0] = r * cos(theta);
		x2[1] = r * sin(theta);

		pt.set_cartesian_coord((x1[0] + x2[0]) / 2., (x1[1] + x2[1]) / 2.);
		t[0] = x2[0] - x1[0];
		t[1] = x2[1] - x1[1];
		pt.set_tangent(t[0], t[1]);


		pt.set_ds_dt(sqrt(t[0] * t[0] + t[1] * t[1]));

		pnt_[pt.theta()] = pt;
		x1[0] = x2[0];
		x1[1] = x2[1];
	}


	PointLocator P_y_max, P_y_min;
	search_vertical_extremes(&P_y_min, &P_y_max);

	y_lim_[0] = P_y_min.x[1];
	y_lim_[1] = P_y_max.x[1];
	
	return true;
}




 int				CircleMapped2DCurve::set(const std::vector<REAL>& x, const std::vector<REAL>& y)
{

	if (x.size() < 6 || x.size() != y.size())
		return 0;


	// memorizzo i punti in coordinate polari in una mappa ordinata per angoli crescenti nel range [0,2*pi)
	std::map<REAL, REAL> vertex;
	double max_y = -1e20, min_y = 1E20;
	for (size_t i = 0; i != x.size(); ++i)
	{
		REAL	X = x[i],
				Y = y[i];

		if (Y > max_y)
			max_y = Y;
		if (Y < min_y)
			min_y = Y;
		REAL    theta = angle(X, Y);
		REAL    r = sqrt(X*X + Y*Y);
		vertex[theta] = r;
	}

	// verifico la convessità
	if (verify_convexity(vertex) != 1)
		return -1;

	double y_range = max_y - min_y;
	REAL x1[2], x2[2], t[2], theta, r;

	// verifico l'esistenza e la consistenza del taglio superiore
	auto it = vertex.lower_bound(MML3::Math::half_Pi);
	auto it_ = it;
	if (it == vertex.begin())
		it = vertex.end();
	--it;

	r = it->second;
	theta = it->first;
	x1[0] = r * cos(theta);
	x1[1] = r * sin(theta);

	r = it_->second;
	theta = it_->first;
	x2[0] = r * cos(theta);
	x2[1] = r * sin(theta);


	if (fabs(x1[1] - x2[1]) > 1e-12*y_range)
		return -2;
	if (x2[0] >= 0 || x1[0] <= 0)
		return -3;


	// verifico l'esistenza e la consistenza del taglio inferiore
	it = vertex.lower_bound(3. / 2. * MML3::Math::Pi);
	it_ = it;
	if (it == vertex.begin())
		it = vertex.end();
	--it;

	r = it->second;
	theta = it->first;
	x1[0] = r * cos(theta);
	x1[1] = r * sin(theta);

	r = it_->second;
	theta = it_->first;
	x2[0] = r * cos(theta);
	x2[1] = r * sin(theta);


	if (fabs(x1[1] - x2[1]) > 1e-12*y_range)
		return -2;
	if (x2[0] <= 0 || x1[0] >= 0)
		return -3;



	// per ogni segmento della poligonale definito da una coppia di vertici adiacenti calcolo 
	// a) il punto medio,
	// b) la tangente al segmento
	// c) ds/dt che assumo uguale alla lunghezza del segmento, come se al variare del parametro t tra 0 ed 1
	//    l'ascissa curviliea vari da 0 alla lunghezza del segmento



	cpoint pt;
	it = vertex.begin();
	r = it->second;
	theta = it->first;
	x1[0] = r * cos(theta);
	x1[1] = r * sin(theta);
	bool end = false;
	for (; !end;)
	{
		++it;
		if (it == vertex.end())
		{
			it = vertex.begin();
			end = true;
		}

		r = it->second;
		theta = it->first;
		x2[0] = r * cos(theta);
		x2[1] = r * sin(theta);

		// se x1[1]==x2[1] allora siamo su un taglio orizzontale che attraversa l'asse verticale
		if (fabs(x2[1] - x1[1]) < 1e-12*y_range)
			pt.set_cartesian_coord(0., (x1[1] + x2[1]) / 2.);
		else
			pt.set_cartesian_coord((x1[0] + x2[0]) / 2., (x1[1] + x2[1]) / 2.);

		t[0] = x2[0] - x1[0];
		t[1] = x2[1] - x1[1];
		pt.set_tangent(t[0], t[1]);
		pt.set_ds_dt(sqrt(t[0] * t[0] + t[1] * t[1]));
		pnt_[theta] = pt;
		x1[0] = x2[0];
		x1[1] = x2[1];
	}


	PointLocator P_y_max, P_y_min;
	search_vertical_extremes(&P_y_min, &P_y_max);

	y_lim_[0] = P_y_min.x[1];
	y_lim_[1] = P_y_max.x[1];

	return 1;
}


 int CircleMapped2DCurve::verify_convexity(const std::map<REAL, REAL>&  poly_vertex)
{
	// verifico la convessità
	REAL x1[2], x2[2], theta0;
	auto it = poly_vertex.begin();
	x1[0] = it->second * cos(it->first);
	x1[1] = it->second * sin(it->first);
	++it;
	x2[0] = it->second * cos(it->first);
	x2[1] = it->second * sin(it->first);
	theta0 = angle(x2[0] - x1[0], x2[1] - x1[1]);

	Vector2 V1 = { x2[0] - x1[0], x2[1] - x1[1] }, V2;

	V1.unitize();
	bool end(false);
	for (; !end;)
	{
		++it;
		if (it == poly_vertex.end())
		{
			it = poly_vertex.begin();
			end = true;
		}

		x2[0] = it->second * cos(it->first);
		x2[1] = it->second * sin(it->first);
		V2[0] = x2[0] - x1[0];
		V2[1] = x2[1] - x1[1];
		V2.unitize();

		// calcolo il seno dell'angolo tra due segmenti successivi
		double cr = cross(V1, V2);
		// la tolleranza sulla non convessità di due segmenti consecutivi è di 1  grado
		// corrispondente a 0.017 radianti
		if (cr < -0.017)
			return -1;
		V1 = V2;
		x1[0] = x2[0];
		x1[1] = x2[1];
	}
	return 1;
}


 void CircleMapped2DCurve::evaluate_at(const PointLocator* lc, double v[2], double dv[2], double ddv[2] )
{
	static  MML3::CubicSplineCurve<2> CS2;
	
	if (v)
	{
		v[0] = lc->x[0];
		v[1] = lc->x[1];
	}
	// se è richiesto solo il punto c'è poco da fare
	if (!dv && !ddv)
		return;


	configure_spline_(*(lc->P1), *(lc->P2), CS2);
	CS2.eval(lc->t, nullptr, dv, ddv);
}



 void CircleMapped2DCurve::configure_spline_(const cpoint& P1, const cpoint& P2, MML3::CubicSplineCurve<2>& CSpline)
{
	
	double D1[2], D2[2];
	P1.get_dp_dt(D1);
	P2.get_dp_dt(D2);

	double pp[] = {
		P1.x(), D1[0], P2.x(), D2[0],
		P1.y(), D1[1], P2.y(), D2[1] };

	CSpline.configure(pp,1e-7);
}



 double CircleMapped2DCurve::intersect(REAL x, REAL y, PointLocator* lc,  double tolerance)const
{

	// numero massimo di iterazioni per la rifinitura dell'intersezione
	int MAXIT = 50;
	


	REAL theta = angle(x, y); 
	REAL r = sqrt(x*x + y*y);

	point_map_t::const_iterator it_P1, it_P2;
	find_bounds_(theta, it_P1, it_P2);



	double DT = it_P2->first - it_P1->first;
	if (DT < 0)
		DT += MML3::Math::two_Pi;

	MML3::CubicSplineCurve<2> cubic_curve;
	configure_spline_(it_P1->second, it_P2->second, cubic_curve);

	// stima della pseudo coordinata angolare local_t in base all'angolo angolo del punto (x,y) misurato a partire dall'angolo di P1;
	// NB è solo una stima grezza in quanto la coordinata della spline non dipende linearmente dalla coordinata angolare
	double  local_t = (theta - it_P1->first) / DT;
	 	
	double X[2],  res[2], k[2][2],dt[2];
	
	// stima iniziale del punto intersezione
	cubic_curve.eval(local_t, X, nullptr, nullptr);
	// stima iniziale di alpha
	double alpha = radius(X[0], X[1])/r;
	

	// iterazioni alla Newton per rifinire la stima 
	double err ;
	// parte costante della matrice delle derivate el residuo
	k[0][1] = -x;
	k[1][1] = -y;
	for (int n = 0; n != MAXIT; ++n)
	{
		// calcolo il residuo
		res[0] = X[0] - alpha *x;
		res[1] = X[1] - alpha *y;
		// e l'errore relativo
		err = sqrt(res[0] * res[0] + res[1] * res[1])/(alpha*r);
		// termino le iterazioni se l'errore è piccolo
		if (err < tolerance)
			break;
		// calcolo le derivate del residuo rispetto a local_t
		cubic_curve.eval(local_t, nullptr, dt, nullptr);
		k[0][0] = dt[0];
		k[1][0] = dt[1];
		
		// calcolo la variazione delle incognite
		MML3::Math::Matrix<2, double>::solve(k[0], res,  dt );
		// aggiorno le stime 
		local_t -= dt[0];
		alpha	-= dt[1];
		// aggiorno X[]
		cubic_curve.eval(local_t, X, nullptr, nullptr);
	}

	if (err > tolerance)
		throw std::exception("MML3::GEO::CircleMapped2DCurve::position_: convergence failure");
	if (lc)
	{
		lc->P1= &it_P1->second;
		lc->P2= &it_P2->second;
		lc->t = local_t;
		lc->x[0] = X[0];
		lc->x[1] = X[1];
	}
	
	return alpha;

	


}






//inline bool CircleMapped2DCurve::horizontal_intersect(REAL Y, REAL X, PointLocator* lc)const
//{
//
//
//	MML3::CubicSplineCurve<2> cubic_curve;
//
//	
//	// ciclo sulle splines cubiche che compongono la frontiera
//	for (auto it1 = pnt_.begin(); it1 != pnt_.end(); ++it1)
//	{
//		//-------------------------------------------------
//		// individuo i due punti estremi della spline
//		auto it2 = it1;
//		++it2;
//		if (it2 == pnt_.end())
//			it2 = pnt_.begin();
//
//		const cpoint& P1 = it1->second;
//		const cpoint& P2 = it2->second;
//
//		//  ----------------------------------------------------------------------------
//		//  test rapidi per escludere che l'intersezione stia nell'intervallo corrente
//
//		// 1) se le coordinate x dei due punti hanno entrambe segno opposto ad X, oppure 
//		// una ha segno opposto e l'altra e' nulla, escludo la spline corrente
//		if (P1.x() *X <= 0 && P2.x() *X <= 0)
//			continue;
//
//		// 2) se Y>=0 ed entrambi i punti P1 e P2 sono sopra Y , oppure
//		//    se Y<0  ed entrambi i punti P1 e P2 sono sotto Y escludo
//		if (Y >= 0)
//		{
//			if (P1.y() > Y && P2.y() > Y)
//				continue;
//		}
//		else
//		{
//			if (P1.y()< Y && P2.y() < Y)
//				continue;
//
//
//		}
//		
//		// angolo compreso tra P1 e P2
//		REAL DT = P2.theta() - P1.theta();
//		// l'ultima spline potrebbe essere a cavallo dell'asse orizzontale
//		if (DT < 0)
//			DT += MML3::Math::two_Pi;
//
//		
//
//
//		//---------------------------------------------------------------------------------
//		// cerco  le intersezioni della retta orizzontale y=Y con la spline corrente
//	
//		configure_spline_(P1, P2, cubic_curve);
//		// coefficienti polinomiali della componente y  della spline  nei termini del pseudo-angolo 
//		// tali cioè che p_y(t) = c[0] + c[1] * t + c[2]*t^2 + c[3]*t^3
//		const double* c = cubic_curve.cy();
//		
//		double z[3];
//		// radici reali dell'equazione p_y(z) - Y=0
//		// attenzione:: i coefficienti sono relativi al parametro della spline
//		int nr = MML3::Math::poly_real_root(c[3], c[2], c[1], c[0] - Y, z[0], z[1], z[2]);
//
//
//
//		// ciclo sulle radici reali
//		for (int i = 0; i != nr; ++i)
//		{
//			double root = z[i];
//			// se le radici della cubica cadono fuori dall'intervallo  [0,DT] non sono intersezioni
//			if (root > 1 || root < 0)
//				continue;
//
//			
//			double r[2];
//
//			cubic_curve.eval(root, r, NULL, NULL);
//			if (fabs(r[1] - Y) > 1e-6)
//				throw
//				std::exception("MML3::GEO::CircleMapped2DCurve::horizontal_intersect");
//
//			if (r[0] * X >= 0) /* radice trovata */
//			{
//				lc->P1 = &it1->second;
//				lc->P2 = &it2->second;
//				lc->x[0] = r[0];
//				lc->x[1] = r[1];
//				lc->t = root;
//				return true;
//			}
//		}
//	}
//
//
//	return false;
//
//}










 bool CircleMapped2DCurve::horizontal_intersect(REAL Y, PointLocator* left_lc, PointLocator* right_lc)const
{


	MML3::CubicSplineCurve<2> cubic_curve;


	if(left_lc)
		left_lc->x[0] = std::numeric_limits<double>::max();
	if (right_lc)
		right_lc->x[0] = std::numeric_limits<double>::lowest();

	double Y_ref = y_lim_[1] - y_lim_[0];

	

	// ciclo sulle splines cubiche che compongono la frontiera
	for (auto it1 = pnt_.begin(); it1 != pnt_.end(); ++it1)
	{
		//-------------------------------------------------
		// individuo i due punti estremi della spline
		auto it2 = it1;
		++it2;
		if (it2 == pnt_.end())
			it2 = pnt_.begin();

		const cpoint& P1 = it1->second;
		const cpoint& P2 = it2->second;

		//  ----------------------------------------------------------------------------
		//  test rapidi per escludere che l'intersezione stia nell'intervallo corrente
		double y1 = P1.y();
		double y2 = P2.y();

		double x1 = P1.x();
		double x2 = P2.x();
		

		//   se Y>=0 ed entrambi i punti P1 e P2 sono sopra Y , oppure
		//   se Y<0  ed entrambi i punti P1 e P2 sono sotto Y escludo
		if (Y >= 0)
		{
			if (P1.y() > Y && P2.y() > Y)
				continue;
		}
		else
		{
			if (P1.y()< Y && P2.y() < Y)
				continue;


		}


		double theta2 = P2.theta();
		// angolo compreso tra P1 e P2
		REAL DT = P2.theta() - P1.theta();
		// l'ultima spline potrebbe essere a cavallo dell'asse orizzontale
		if (DT < 0)
			DT +=  MML3::Math::two_Pi;




		//---------------------------------------------------------------------------------
		// cerco  le intersezioni della retta orizzontale y=Y con la spline corrente

		configure_spline_(P1, P2, cubic_curve);
		// coefficienti polinomiali della componente y  della spline  nei termini del pseudo-angolo 
		// tali cioè che p_y(t) = c[0] + c[1] * t + c[2]*t^2 + c[3]*t^3
		const double* c = cubic_curve.cy();

		double z[3];
		// radici reali dell'equazione p_y(z) - Y=0
		// attenzione:: i coefficienti sono relativi al parametro della spline
		int nr = MML3::Math::poly_real_root(c[3], c[2], c[1], c[0] - Y, z[0], z[1], z[2]);



		// ciclo sulle radici reali
		for (int i = 0; i != nr; ++i)
		{
			double root = z[i];
			// se le radici della cubica cadono fuori dall'intervallo  [0,DT] non sono intersezioni
			if (root > (1+1e-6) || root < -1E-6)
				continue;
			
			double r[2];
			cubic_curve.eval(root, r, NULL, NULL);
			if (fabs(r[1] - Y)/Y_ref > 1e-6)
				throw
				std::exception("MML3::GEO::CircleMapped2DCurve::horizontal_intersect");

			if (right_lc && r[0] > right_lc->x[0])
				right_lc->set(&it1->second, &it2->second, r, root);
				
			if (left_lc && r[0] < left_lc->x[0])
				left_lc->set(&it1->second, &it2->second, r, root);
		}

		// se tutte le intersezioni richieste sono state trovate  posso terminare precocemente

		bool left_ok	= left_lc  ? (left_lc->x[0]  <=0):true;
		bool right_ok	= right_lc ? (right_lc->x[0] >=0):true;
		bool ok = left_ok && right_ok;
		
		if (ok && (left_lc &&  right_lc))
			ok = left_lc->x[0] < right_lc->x[0];
		if (ok)
			return true;
	}

		bool left_ok	= left_lc	? (left_lc->x[0]	<=0) : true;
		bool right_ok	= right_lc	? (right_lc->x[0]	>=0) : true;
	if (left_ok && right_ok)
		return true;
	else
		return false;

}





 void CircleMapped2DCurve::search_vertical_extremes(PointLocator* pl_bottom, PointLocator* pl_top)const
{


	MML3::CubicSplineCurve<2> cubic_curve;


	pl_bottom->x[1] = std::numeric_limits<double>::max();
	pl_top->x[1] = std::numeric_limits<double>::lowest();



	// inizializzo entrambi gli estremi con il punto iniziale
	auto it = pnt_.begin();
	auto ito = it;
	++ito;
	pl_top->set(&(it->second), &(ito->second), it->second.x(), it->second.y(), 0);
	pl_bottom->set(&(it->second), &(ito->second), it->second.x(), it->second.y(), 0);


	// ciclo sulle splines cubiche che compongono la frontiera
	for (auto it1 = pnt_.begin(); it1 != pnt_.end(); ++it1)
	{
		//-------------------------------------------------
		// individuo i due punti estremi della spline
		auto it2 = it1;
		++it2;
		if (it2 == pnt_.end())
			it2 = pnt_.begin();

		const cpoint& P1 = it1->second;
		const cpoint& P2 = it2->second;

		//  ----------------------------------------------------------------------------
		//  test rapidi di esclusione dell'intervallo corrente
		// .... TODO! 


		

		// angolo compreso tra P1 e P2
		REAL DT = P2.theta() - P1.theta();
		// l'ultima spline potrebbe essere a cavallo dell'asse orizzontale
		if (DT < 0)
			DT += MML3::Math::two_Pi;




		//---------------------------------------------------------------------------------
		// cerco  i punti estremi in y sulla  spline corrente

		configure_spline_(P1, P2, cubic_curve);
		// coefficienti polinomiali della componente y  della spline  nei termini del pseudo-angolo 
		// tali cioè che p_y(t) = c[0] + c[1] * t + c[2]*t^2 + c[3]*t^3
		const double* c = cubic_curve.cy();

		double z[3];
		// radici reali dell'equazione p'_y(z)=0
		// attenzione:: i coefficienti sono relativi al parametrodella spline
		int nr = MML3::Math::poly_real_root(0, 3*c[3], 2*c[2], c[1], z[0], z[1], z[2]);

		// ciclo sulle radici reali
		for (int i = 0; i != nr; ++i)
		{
			double root = z[i];
			// se le radici  cadono fuori dall'intervallo  (0,1) non interessano
			// ricordo che il test negli estremi degli ntervalli è svolto fuori
			if (root >= 1 || root <= 0)
				continue;
			
			double r[2];
			cubic_curve.eval(root, r, NULL, NULL);

			if (r[1] > pl_top->x[1])
				pl_top->set(&it1->second, &it2->second, r, root);

			else if (r[1] < pl_bottom->x[1])
				pl_bottom->set(&it1->second, &it2->second, r, root);
		}

		// controllo se l'estremo locale è sul punto finale . ricordo che il punto iniziale è già stato testato
		// come punto finale del ramo precedente

		if (P2.y() > pl_top->x[1])
			pl_top->set(&it1->second, &it2->second, P2.x(),P2.y(), 1);
		else if (P2.y() < pl_bottom->x[1])
			pl_bottom->set(&it1->second, &it2->second, P2.x(), P2.y(), 1);
	
	}
}




 void		CircleMapped2DCurve::find_bounds_(const double theta, point_map_t::const_iterator& It1, point_map_t::const_iterator& It2)const
{

	It2 = pnt_.lower_bound(theta);
	// se it2 è il primo vertice
	if (It2 == pnt_.begin())
		It1 = --pnt_.end();
	else if (It2 == pnt_.end())
	{
		It1 = --pnt_.end();
		It2 = pnt_.begin();

	}
	else
		(It1 = It2)--;

}




 void write_vtk_polydata(std::ofstream& os, const  CircleMapped2DCurve& curve, size_t vertical_steps)
{

	double y_min = curve.get_y_min();
	double y_max = curve.get_y_max();
	double Dy = (curve.get_y_max() - curve.get_y_min()) / vertical_steps;
	MML3::GEO::CircleMapped2DCurve::PointLocator  PntL, PntR;

	std::vector<double> x_left(vertical_steps + 1), x_right(vertical_steps + 1), Y(vertical_steps + 1);

	// il numero di punti della griglia verticale è vertical_steps +1


	
	for (size_t i = 0; i <= vertical_steps;++i)
	{
		double y = curve.get_y_min() + i*Dy;
		if (i == 0)
			y = curve.get_y_min() + Dy / 10;
		else if (i == vertical_steps)
			y = curve.get_y_max() - Dy/10;
		
		bool rvb_ = curve.horizontal_intersect(y, &PntL, &PntR);

		if (!rvb_)
			throw std::exception("MML3::GEO:: write_vtk_polydata(std::ofstream& , const  CircleMapped2DCurve& ...): failure");

		Y[i] = y;
		x_left[i]  = PntL.x[0];
		x_right[i] = PntR.x[0];
	}

	// adesso determino i due punti estremi
	curve.search_vertical_extremes(&PntL, &PntR);



	// costruisco i punti della poligonale partendo dall'estremo inferiore e girando in senso antiorario
	typedef MML3::GEO::Point<2, double> APoint;
	std::vector<APoint> Point;
	Point.reserve(2 * vertical_steps + 4);

	Point.push_back(APoint(PntL.x));
	for (auto i = 0; i != Y.size(); ++i)
			Point.push_back({ x_right[i], Y[i] });
	Point.push_back(APoint(PntR.x));
	for (int i = Y.size()-1; i >=0; --i)
		Point.push_back({ x_left[i], Y[i] });


	// VTK output


	os << "# vtk DataFile Version 2.0\n"
		<< "created by write_vtk_polydata(... CircleMapped2DCurve ...)\n"
		<< "ASCII\n"
		<< "DATASET POLYDATA\n";

	os << "POINTS " <<Point.size() << " DOUBLE\n";
	for (auto el :Point)
		os << el[0] << " " << el[1]  << " 0" << std::endl;

	
	os << "\nLINES " << Point.size()  << " " << (Point.size()  *3) << std::endl;
	for (int line = 0; line != Point.size()-1 ; ++line)
		os << "2 " << line << " " << line + 1 << std::endl;

	// linea di chiusura della poligonale
	os << "2 " << Point.size()-1 << " " << 0 << std::endl;
	


}


 void write_vtk_polygonal(std::ofstream& os, const char* ztstring, size_t size, const  double* x, const  double* y, const double* z, bool close)
{

	// VTK output
	os << "# vtk DataFile Version 2.0\n"
		<< ztstring
		<< "\nASCII\n"
		<< "DATASET POLYDATA\n";

	os << "POINTS " << size << " DOUBLE\n";
	for (size_t i = 0; i != size; ++i)
		os <<x[i] << " " << y[i] << " " << z[i] << std::endl;

	int num_lines = close ? size : size - 1;
	os << "\nLINES " << num_lines << " " << (num_lines * 3) << std::endl;
	for (int line = 0; line != size - 1; ++line)
		os << "2 " << line << " " << line + 1 << std::endl;
	
	// linea di chiusura della poligonale
	if (close)
		os << "2 " << size-1 << " " << 0 << std::endl;

}
 



	}// end namespace GEO

}// end namespace MML