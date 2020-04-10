#include <math.h>
#include <iostream>
#include <fstream>
#include "vfield.hpp"
#include <ctime>

#define lambda 1.4e-7 /* to do: parameter */
#define B 1. /* ---""--- */

#define XMIN (-5)
#define XMAX   5
#define YMIN (-3)
#define YMAX   3

using namespace std;

// def. Konstruktor :
vortex_field_direct::
vortex_field_direct(int numpar  // Zahl der FLn
, double xl, double yl
, double a0, double h0   // h0==0 bedeutet random feld!
, double xfrac, double yfrac)
: num(numpar), xl(xl), yl(yl), a0(a0), h0(h0), xfrac(xfrac), yfrac(yfrac)  //TH Fixed SVN-17
{

	// int i;
	// E.Zero(); //hier noch keine Länge....
	// dxy = 0;
	// fxy = 0;
	// prop = 0;
	// //FMSrem: subnn=0;   //ONLY FOR TEST in relax_step
	allout = true;

  time_t timx;
	time(&timx);
  std::cout << "Random seed: " << timx << std::endl;
  // gsl_rng_default_seed = timx; // Für echtes random, sonst immer dieselbe init...
  gsl_rng_default_seed =1533633677;
    rndm = gsl_rng_alloc (gsl_rng_taus2);  //

	if (xl == 0 || yl == 0) {
		std::cerr << "vortex_field (a0): xl==0 or yl==0" << std::endl;
		exit(89);
	}

	int cnt = 0;

	if (num == 0) {
		if (a0 != 0) {
			num = (int)((xl*yl) / (a0*a0*r_h00) + 1e-5);
		}
		else {
			std::cerr << "vortex_field (num): num==0 and a0==0" << std::endl;
			exit(88);
		}
	}
	else { // num>0 oder num<0
		if (a0 == 0) {
			a0 = sqrt((xl*yl) / (num*r_h00));
		}
	}

	xh = 0.5*xl; yh = 0.5*yl;
	ww_max_dist = (yh < 9) ? yh : 9; // realer Wert in Einh. von lambda

	dx = dy = a0*FACdx;

	nummatch = (int)((xl*yl) / (a0*a0*r_h00) + 0.5);


	if (h0 == 0.0) { // random array
		resize(nummatch);
		int i = 0;
		double minborder = a0*0.6;
		while (i < nummatch) {
			int ix = 2 * i;
			int iy = ix + 1;
			xy(ix) = gsl_rng_uniform(rndm) *xl*xfrac - xh;
			xy(iy) = gsl_rng_uniform(rndm) *yl*yfrac - yh; // -xh und -yh aufheben
			for (int i1 = 0; i1<i; i1++) {
				double distx = xy(ix) - xy(2 * i1);
				if (distx>xh) distx -= xl;
				else if (distx<-xh) distx += xl;
				double disty = xy(iy) - xy(2 * i1 + 1);
				if (disty>yh) disty -= yl;
				else if (disty < -yh) disty += yl;
				if (hypot(distx, disty) < minborder) { // schließe zu nahe F-n aus
					std::cerr << '*';
					i--;
					break;
				}
			}
			i++;
		}

	}
	else if (num < 0) {
		num = 0;
		resize(nummatch);
		fill(nummatch, xfrac, yfrac);

	}
	else { // regular lattice

		if (nummatch != num) {
			std::cerr << "Achtung, num( " << num
				<< " -> " << nummatch << ") verändert!" << std::endl;
			num = nummatch;
		}
		resize(nummatch);
		int ix = 0, iy = 0
			, even = 1; /* Verschiebung für hex Gitter aus Dreieck bei jeder 2. Reihe */
		double xstep = h0*0.5, ystep = a0*0.25;
		double xstart = -xh + xstep, ystart = -yh + 2 * ystep;
		double x, y;
		double *poi = xy.data();
		x = xstart; // int i00=0;

		cnt = 0;
		do  { // Reihen...
			iy = 0; y = ystart; // y-Koordinate vorbereiten
			do { // weitere Punkte in der Spalte
				cnt++;
				*poi++ = x; *poi++ = y + even*ystep; // Punkt schreiben
				iy++; y = ystart + a0*iy; // nächste y-Koordinate
				//cout << i00++ << " (" << xy(i00*2)<< ", "<<xy(i00*2+1)<<")"<< endl; //ix, iy
			} while (y < yh); // Ende y, wnn nächste y über dem Rand ist
			even *= -1; // nächste Spalte vorbereiten
			ix++; x = xstart + h0*ix;
		} while (x < xh);
		std::cout << "\nreg. lattice: " << num << " reserviert, und "
			<< cnt << " belegt. (" << ix << " x " << iy << ")\n" << std::endl;
	} // end else .... regular lattice


	if (h0 == 0) h0 = r_h00*a0;

	cout << "done init " << num * 2 << " doubles" << endl;
	cout << "  xl: " << xl << endl;
	cout << "  yl: " << yl << endl;

	min_dis = MIN_DIS*xl;

  ffac0 = 0.3 * 0.01 / (gsl_sf_bessel_K1(a0)-gsl_sf_bessel_K1(a0+0.01));  // Startwert für Kraft->Weg
  cout << "ffac0 (*0.3) berechnet: " << ffac0 << endl; // weniger bringt auch nichts bei fixem ffac, getestet...

	init_energy_file();

	return;
}


// // Konstruktor mit eingefügter Versetzung
// vortex_field_direct::
// vortex_field_direct(dismode_t disloc_mode, int numpar
// 	     , double xl, double yl
// 	     , double a0, double h0
// 	     , double xfrac, double yfrac)
//   : num(numpar), xl(xl), yl(yl), a0(a0), h0(h0), xfrac(xfrac), yfrac(yfrac)  //TH Fixed SVN-17
// {
//   int i, ix=0, iy=0, even=1; /* Verschiebung für hex Gitter aus Dreieck bei jeder 2. Reihe */
//   double xstep = h0*0.5, ystep=a0*0.25;
//   double xstart, ystart;
//   double x, y;
//   double *poi;
//   double y_temp;
//   vector<int>::iterator iter1, iter2;
//   int xint, xh_2;
//   double delta[-XMIN+XMAX+1] =
//   {0.1, 0.25, 0.35, 0.4, 0.5, -0.5, -0.4, -0.35, -0.25, -0.1, 0.0};
//   double verschiebung;
//   bool versetzt = false;
//   wwf={0};
//   dxy={0};
//   fxy={0};


//   if (xl==0 || yl==0) {
//     std::cerr << "vortex_field (a0): xl==0 or yl==0"<< std::endl;
//     exit(89);
//   }

//   if (num==0) {
//     if (a0!=0) {
//       num = (int)((xl*yl)/(a0*a0*r_h00)+1e-5);
//     } else {
//       std::cerr << "vortex_field (num): num==0 and a0==0"<< std::endl;
//       exit(88);
//     }
//   } else {
//     if (a0==0) {
//       a0 = sqrt((xl*yl)/(num*r_h00));
//     }
//   }

//   xh = 0.5*xl; yh = 0.5*yl;
//   ww_max_dist = (yh<9)?yh:9; // realer Wert in Einh. von lambda

//   dx = dy =a0*FACdx;

//   nummatch = (int)((xl*yl)/(a0*a0*r_h00)+1e-5);

//   switch(disloc_mode)
//   {
//   case vertical:
//    // vertikale Versetzung erzeugen
//    resize(nummatch);
//    poi = xy.data();
//    xstart = -xh+xstep; ystart = -yh+2*ystep;
//    x=xstart;
//    xh_2 = lrint(xh/(h0*2.));
//    i = 0;
//    versetzt = false;

//    do
//    {
//      xint = lrint(x / h0);
//      iy=0; y=ystart;
//      do
//      {
//       y_temp = y+even*ystep;
//       verschiebung = 0.;
//       if((-yh/2. <= y_temp) && (y_temp <= yh/2.))
//       {
//        if((XMIN-xh_2 <= xint) && (xint <= XMAX-xh_2))
//        {
//         verschiebung = -delta[xint-XMIN+xh_2] * a0;
//        }
//        else if((XMIN+xh_2 <= xint) && (xint <= XMAX+xh_2))
//        {
//         verschiebung = delta[xint-XMIN-xh_2] * a0;
//        }
//        if(!versetzt && ((xint == -1-xh_2) || (xint == -xh_2) ||
//                         (xint == -1+xh_2) || (xint ==  xh_2)))
//        {
//         prop(i-2) = -2;
//         prop(i-1) = -2;
//         prop(i  ) = -2; /* Pinnen der Versetzung */
//         prop(i+1) = -2;
//        }
//        versetzt = true;
//       }
//       else
//       {
//        if(versetzt && ((xint == -1-xh_2) || (xint == -xh_2) ||
//                        (xint == -1+xh_2) || (xint ==  xh_2)))
//        {
//         prop(i-2) = -1;
//         prop(i-1) = -1;
//         prop(i  ) = -1; /* Pinnen der Versetzung */
//         prop(i+1) = -1;
//        }
//        versetzt = false;
//       }
//       y_temp += verschiebung;
//       *poi++ = x; *poi++ = y_temp;
//       iy++; i++; y=ystart+a0*iy;
//      } while(y<yh);
//      even *= -1;
//      ix++; x=xstart+h0*ix;
//    } while(x<xh);
//    ext_force = 10 * (sqrt(B / Fi0) * lambda);
//    force_delta = 0;
//    pinnedtostart();

//    break;

//   case none:

//    /* keine Versetzung */

//    resize(nummatch);
//    poi = xy.data();
//    xstart = -xh+xstep; ystart = -yh+2*ystep;
//    x=xstart;
//    xh_2 = lrint(xh/(h0*2.));
//    i = 0;
//    versetzt = false;

//    do
//    {
//      xint = lrint(x / h0);
//      iy=0; y=ystart;
//      do
//      {
//       y_temp = y+even*ystep;
//       *poi++ = x; *poi++ = y_temp;
//       iy++; i++; y=ystart+a0*iy;
//      } while(y<yh);
//      even *= -1;
//      ix++; x=xstart+h0*ix;
//    } while(x<xh);
//    ext_force = 0.;
//    force_delta = 0.;

//    break;

//   default:

//    cerr << "Modus " << disloc_mode << " nicht bekannt!" << endl;
//    abort();
//   }
//   if(h0==0) h0 = r_h00*a0;

//   cout<<"done init "<< num*2 <<" doubles" <<endl;
//   cout<< "  xl: " << xl <<endl;
//   cout<< "  yl: " << yl <<endl;

//   min_dis=MIN_DIS*xl;
//   return;
// }



vortex_field_direct::vortex_field_direct()
{
	//resize(0);
	return;
}

// n n+1 n Konstruktor
// vortex_field_direct::vortex_field_direct(int nh_num_c, int nv_num_c, int n_1_h_num_c, double h0_c, double a0_c, double gsmult_c)
// 	: nh_num(nh_num_c), nv_num(nv_num_c), n_1_h_num(n_1_h_num_c), h0(h0_c), a0(a0_c),gsmult(gsmult_c)
// {
// 	apply_force = false;
//         allout = false;
// 	sum0 = 0.0;
// 	sum_last=0.0;
// 	sumfabs=0.;
// 	n_count = 0;
// 	grenze_rechts = 0.0;
// 	grenze_links = 0.0;
// 	outfil1.str("");
// 	force_param = 2.0*M_PI*lambda;
// 	force_param = force_param*force_param*2.0e11*2.0e-7*1.4e-7/2.07e-15; /* force auf K0 normiert */

// 	n_1_v_num = nv_num + 1;						//Zeilenzahl des (n+1)-Gitters [auch immer GERADE]


// 	n_1_v_dist = (nv_num * a0) / n_1_v_num;				//Vertikaler Abstand im (n+1)-Gitter
// 	n_1_h_dist = sqrt(3.0) * 0.5 * n_1_v_dist;			//Horizontaler Abstand im (n+1)-Gitter


// 	n_1_h_width = n_1_h_dist * (n_1_h_num - 1);			//Breite des (n+1) Gitters

// 	gapsize = (h0 + n_1_h_dist) * 0.5 * gsmult;

// 	xl = ((nh_num - 1) * h0) + (gapsize * 2.0) + n_1_h_width;	//Breite des Gesamtsystems
// 	yl = (nv_num) * a0;						//Höhe des Gesamtsystems

// 	num = ((nv_num) * nh_num) + (n_1_v_num * n_1_h_num);

// 	first_half = (nh_num ) / 2;

// 	xfrac = 1.0;
// 	yfrac = 1.0;
// 	point_cnt = 0;

// 	xh = 0.5*xl; yh = 0.5*yl;
// 	ww_max_dist = (9>yh)?yh:9; // realer Wert in Einh. von lambda

// 	dx = dy = a0*FACdx;

// 	resize(num);


// 	for(int i = 0; i < first_half; i++)				//füllt 1.Hälfte des n-Gitters auf
// 	{
// 		if(i % 2 == 0) //gerade Spalten
// 		{
// 			for(int j = 0; j < nv_num; j++)
// 			{
// 				xy(point_cnt)		= (i * h0) - xh;
// 				xy(point_cnt + 1)	= yh - (j * a0);
// 				point_cnt += 2;
// 			}
// 		}
// 		if(i % 2 != 0) //ungerade Spalten
// 		{
// 			for(int j = 0; j < nv_num; j++)
// 			{
// 				xy(point_cnt)		= (i * h0) - xh;
// 				xy(point_cnt + 1)	= yh - (a0 * 0.5) - (j * a0);
// 				point_cnt += 2;
// 			}
// 		}
// 	}


// 	for(int i = 0; i < n_1_h_num; i++) 				//füllt n+1-Gitter auf
// 	{
// 		if(i % 2 == 0) //gerade Spalten
// 		{
// 			for(int j = 0; j < n_1_v_num; j++)
// 			{
// 				xy(point_cnt)		= ((i * n_1_h_dist) + ((first_half-1) * h0) + gapsize) - xh;
// 				xy(point_cnt + 1)	= yh - (j * n_1_v_dist);
// 				point_cnt += 2;
// 			}
// 		}
// 		if(i % 2 != 0) //ungerade Spalten
// 		{
// 			for(int j = 0; j < n_1_v_num; j++)
// 			{
// 				xy(point_cnt)		= ((i * n_1_h_dist) + ((first_half-1) * h0) + gapsize) - xh;
// 				xy(point_cnt + 1)	= yh - (n_1_v_dist * 0.5) - (j * n_1_v_dist);
// 				point_cnt += 2;
// 			}
// 		}
// 	}

// 	for(int i = 0; i < (first_half); i++)				//füllt 2.Hälfte des n-Gitters auf
// 	{
// 		if(i % 2 == 0) //gerade Spalten
// 		{
// 			for(int j = 0; j < nv_num; j++)
// 			{
// 				xy(point_cnt)		= ((i * h0) + ((n_1_h_num-1) * n_1_h_dist) + ((first_half-1) * h0) + (2.0 * gapsize)) - xh;
// 				xy(point_cnt + 1)	= yh - (j * a0);
// 				point_cnt += 2;
// 			}
// 		}
// 		if(i % 2 != 0) //ungerade Spalten
// 		{
// 			for(int j = 0; j < nv_num; j++)
// 			{
// 				xy(point_cnt)		= ((i * h0) + ((n_1_h_num-1) * n_1_h_dist) + ((first_half-1) * h0) + (2.0 * gapsize)) - xh;
// 				xy(point_cnt + 1)	= yh - (a0 * 0.5) - (j * a0);
// 				point_cnt += 2;
// 			}
// 		}
// 	}


// 	for(int i = 0; i < (point_cnt - 1); i = i + 2)
// 		{
// 			cout << xy(i) << " " << xy(i + 1) << endl;
// 		}


// 	cout << "Breite: " << xl << " Hoehe: " << yl << " (n+1)-Abstand: " << n_1_h_dist << " Teilchen (gezaehlt): " << (double) point_cnt / 2 << "Teilchen (berechnet): " << num << endl;
// 	num = point_cnt / 2;


// }


/* Alles, was in diesen zwei Funktionen gemacht wird, wird nirgends wirklich verwendet.
void vortex_field_direct::toggle_force(bool set)
{
	apply_force = set;
	return;
}

void vortex_field_direct::set_force(double nf)
{
	if (nf == ext_force) return;
	fscale = 0;
	fprefac_reset(); // set fprefac to 1
	force_delta = nf - ext_force; // welcher Kraft-Schritt, für neues ffac
	ext_force = nf;
	return;
}
*/

void vortex_field_direct::resize(int n)
{
	int xymax = n * 2;
	// int alt_num = prop.numElements();
	nummax = n;
	wwmax = (n*(n + 1)) / 2;

  E.conservativeResize(n); // Energie(max_num)
  oE.conservativeResize(n); // Energie(max_num)
  // Ecurr.conservativeResize(n); // Energie(max_num)
  cosa.conservativeResize(n); // Hilfsgröße für adaptive Schritte (array nur für debug)
  prop.conservativeResize(n); // Eigenschaften, ob geändert o.ä.
	xy.conservativeResize(xymax); // Vektor für FL-Koordinaten (max_xy):
	oxy.conservativeResize(xymax); // Vektor für FL-Koordinaten (max_xy):
	// x0, y0, x1, y1, ...
	fxy.conservativeResize(xymax);  // Vektor für FL-Kraft oder Bewegung(max_xy)
	ofxy.conservativeResize(xymax);  // Vektor für FL-Kraft oder Bewegung(max_xy)
	fext_xy.conservativeResize(xymax);

	dxy.conservativeResize(xymax);
	odxy.conservativeResize(xymax);
	fval.conservativeResize(n);
	ofval.conservativeResize(n);
  cout << "wwmax: " << wwmax << endl;

  ffac.conservativeResize(n);

  min_dis = MIN_DIS*xl;

	// diffi.conservativeResize(xymax);

  // int ix, iy;

  // for (int i = alt_num; i < n; i++)	// todo: Müssen die alle initialisiert werden?
	// 	{
	// 		ix = 2 * i;
	// 		iy = ix + 1;

	// 		oE(i) = 0;
	// 		E(i) = 0;
	// 		prop(i) = 0;

	// 		xy(ix) = 0;
	// 		xy(iy) = 0;

	// 		fxy(ix) = 0;
	// 		fxy(iy) = 0;
	// 		fext_xy(ix) = 0;
	// 		fext_xy(iy) = 0;

	// 		dxy(ix) = 0;
	// 		dxy(iy) = 0;

	// 		fval(i) = 0;
	// 		diffi(i) = 0;
	// 	}

	return;
}


vortex_field_direct::~vortex_field_direct(void)
{}

int vortex_field_direct::getnum(void)
{
	return num;
}

double vortex_field_direct::getXsize(void)
{
	return xl;
}

double vortex_field_direct::getYsize(void)
{
	return yl;
}

// int vortex_field_direct::verzerr (double meas)
// {
//   double fac = 1.0 +meas;

//   xl = xl * fac; xh = xl*0.5;
//   yl = yl / fac; yh = yl*0.5;

//   return 0;
// }


int vortex_field_direct::ww_save(string filnam, int prec)
{
	ostringstream outfilx;

	int i = 0;
	double Esum = 0;

	outfilx << setprecision(prec);

	outfilx << "x y fx fy dx dy ffac oE dE cosa fext_x fext_y prop" /* << " dx dy oE -fval*ffac0*a0 E?? prop force_x force_y" */ << endl;

	int ix, iy;

	for (i = 0; i < num; i++) {
		ix = 2 * i;
		iy = ix + 1;

		outfilx << xy(ix);
		outfilx << " " << xy(iy) << " "
            << fxy(ix) << " "
            << fxy(iy) << " "
            << dxy(ix) << " "
            << dxy(iy) << " "
            << ffac(i) << " "
			// << max_force << " "
            << E(i) << ' '
			//<< -fval(i)*ffac0*a0 << ' '
            << (E(i) - oE(i) + odxy(i)*fval(i)) << " "
            << cosa(i) << " "
            << fext_xy(ix) << " " << fext_xy(iy)	// David
            << " " << prop(i)
            << endl;
		Esum += E(i);
	}

	outfilx << "# " << num << " " << xl << " " << yl << " " << getNpinned() << endl;
	ofstream filestream;
	filestream.open(filnam.c_str(), ios::out);
	filestream << outfilx.str();
	filestream.close();
	std::cout << "E saved (" << filnam << ")  = " << Esum << std::endl;
	return 1;
}



int vortex_field_direct::xy_save(string filnam, int prec)
{
	double Ex;
	ofstream outfilx;
	int i = 0;
	E.fill(0); // reset all ww !!!!
	fxy.fill(0);
	double Esum = 0;

	outfilx.open(filnam.c_str(), ios::out);
	outfilx << setprecision(prec);

	for (i = 0; i < num; i++) {
		Ex = f_analyze(i, true, false);
    int ix=2*i;
    int iy=ix+1;
		outfilx << xy(ix);
		outfilx << " " << xy(iy) << " "
            << fxy(ix) << " "
            << fxy(iy) << " "
            << Ex << endl;
		Esum += Ex;
	}

	outfilx << "# " << num << " " << xl << " " << yl << " " << getNpinned() << endl;
	outfilx.close();
	std::cout << "E saved (" << filnam << ")  = " << Esum << std::endl;

	return 1;
}


/*int vortex_field_direct::fxy_save( string filnam, int prec )
{
ofstream outfil;
outfil.open(filnam.c_str(), ios::out);
outfil<<setprecision(prec);

int i=0;
//  cout<<filnam<<": "<<endl;
for (int i0=0; i0<num; i0++) {
// outfil <<
outfil << fxy(i++)/a0;
outfil << " " << fxy(i++)/a0 << endl;
}
outfil << "# " << num << " " << xl << " " << yl << " " << getNpinned() << endl;
outfil.close();
return 1;
}*/

bool  vortex_field_direct::read_xy(string filnam)
{
	ifstream infil;
	char istr[60];
	cout << "Versuche relaxiertes disloc-file einzulesen: ";
	infil.open(filnam.c_str(), ios::in);
	if (!infil) { cout << "no file!" << endl; return false; }

	infil.seekg(-60, std::ios::end);
	infil.get(istr, 60, '/');
	if (infil.eof()) {
		cout << "unerwartetes file-Ende: " << istr << endl;
		return false;
	}
	infil.getline(istr, 60);
	if (istr[0] == '#') {
		string iistr(string(istr + 1));
		ostringstream ostr;
		ostr << " " << num << " " << xl << " " << yl << endl;
		if (!iistr.compare(ostr.str())) {
			cout << "falsche Dimensionen, muss neu rechnen" << endl;
			return false;
		}
		infil.seekg(0, std::ios::beg);
		int i = 0;
		for (int i0 = 0; i0 < num; i0++) {
			infil >> xy(i++);
			infil >> xy(i++);
		}
		cout << "... ok" << endl << endl;
		return true;
	}
	else {
		cout << "kein #: " << istr << endl;
		return false;
	}

}

int vortex_field_direct::save(string filnam)
{
	ofstream ofi;
	long sz;

	cout << "status speichern " << endl;

	ofi.open(filnam.c_str(), ios::out | ios::binary);
	sz = sizeof(*this);
	ofi.write((char*)&sz, sizeof(sz));
	ofi.write((char*)this, sizeof(*this));

	sz = prop.size();
	ofi.write((char*)&sz, sizeof(sz));
	ofi.write((char*)&prop(0), sizeof(prop(0))*sz);

	sz = oE.size();
	ofi.write((char*)&sz, sizeof(sz));
	ofi.write((char*)&oE(0), sizeof(oE(0))*sz);

	sz = xy.size();
	ofi.write((char*)&sz, sizeof(sz));
	ofi.write((char*)&xy(0), sizeof(xy(0))*sz);

	sz = fxy.size();
	ofi.write((char*)&sz, sizeof(sz));
	ofi.write((char*)&fxy(0), sizeof(fxy(0))*sz);

	sz = E.size();
	ofi.write((char*)&sz, sizeof(sz));
	ofi.write((char*)&E(0), sizeof(E(0))*sz);

	sz = fval.size();
	ofi.write((char*)&sz, sizeof(sz));
	ofi.write((char*)&fval(0), sizeof(fval(0))*sz);

	sz = dxy.size();
	ofi.write((char*)&sz, sizeof(sz));
	ofi.write((char*)&dxy(0), sizeof(dxy(0))*sz);

	//sz=dxy_sum.size();    // Oleg experiment with dxy_sums
	//ofi.write((char*)&sz,sizeof(sz));
	//ofi.write((char*)&dxy_sum(0),sizeof(dxy_sum(0))*sz);

	ofi.close();
	return 0;
}

int vortex_field_direct::read(string filnam)
{
	ifstream ifi;
	long sz = 0;

	cout << "Status aus Datei einlesen: " << filnam.c_str() << endl;

	ifi.open(filnam.c_str(), ios::in | ios::binary);

	if (ifi.good()){

		//ifi.read((char*)&temp,sizeof(*this));
		//Geht nicht da die Array Pointer nicht überschrieben werden dürfen
		ifi.read((char*)&sz, sizeof(sz));
		ifi.seekg(sz + sizeof(sz), ios::beg);

		ifi.read((char*)&sz, sizeof(sz));
		prop.resize(sz);
		ifi.read((char*)&prop(0), sizeof(prop(0))*sz);

		ifi.read((char*)&sz, sizeof(sz));
		oE.resize(sz);
		ifi.read((char*)&oE(0), sizeof(oE(0))*sz);

		ifi.read((char*)&sz, sizeof(sz));
		xy.resize(sz);
		ifi.read((char*)&xy(0), sizeof(xy(0))*sz);

		ifi.read((char*)&sz, sizeof(sz));
		fxy.resize(sz);
		ifi.read((char*)&fxy(0), sizeof(fxy(0))*sz);

		ifi.read((char*)&sz, sizeof(sz));
		E.resize(sz);
		ifi.read((char*)&E(0), sizeof(E(0))*sz);

		ifi.read((char*)&sz, sizeof(sz));
		fval.resize(sz);
		ifi.read((char*)&fval(0), sizeof(fval(0))*sz);

		ifi.read((char*)&sz, sizeof(sz));
		dxy.resize(sz);
		ifi.read((char*)&dxy(0), sizeof(dxy(0))*sz);

		ifi.read((char*)&sz, sizeof(sz));   // Oleg experiment with dxy_sums
		// dxy_sum.resize(sz);
		// ifi.read((char*)&dxy_sum(0),sizeof(dxy_sum(0))*sz);

		//ifi.read((char*)&sz,sizeof(sz));   // Oleg experiment with dxy_sums
		//dxy_temp.resize(sz);
		//ifi.read((char*)&dxy_temp(0),sizeof(dxy_temp(0))*sz);

	}
	else
		cout << "Fehler beim Einlesen der Datei!" << endl;
	ifi.close();

	return 0;
}

bool  vortex_field_direct::read_new_xy(string filnam)
{
	ifstream infil;
	char istr[60];
	cout << "Versuche neues Positions-file einzulesen: ";
	infil.open(filnam.c_str(), ios::in);
	if (!infil) { cout << "no file!" << endl; return false; }

	infil.seekg(-60, std::ios::end);
	infil.get(istr, 60, '/');
	if (infil.eof()) {
		cout << "unerwartetes file-Ende: " << istr << endl;
		return false;
	}
	infil.getline(istr, 60);
	if (istr[0] == '/' && istr[1] == '/') {
		string iistr(string(istr + 2));
		istringstream iiistr(iistr);
		iiistr >> num >> xl >> yl;
		cout << num << ": " << xl << ", " << yl << endl;
		nummax = num;

		resize(num);

		infil.seekg(0, std::ios::beg);
		int i = 0;
		for (int i0 = 0; i0 < num; i0++) {
			infil >> xy(i++);
			infil >> xy(i++);
		}
		cout << "... ok" << endl << endl;
		return true;
	}
	else {
		cout << "kein //: " << istr << endl;
		return false;
	}

}


int vortex_field_direct::v_add(double xadd = 0, double yadd = 0, int pinned = 0)
// add new vortex to field ... Neue FLL dazugeben
{
	int resizei = -1;

	//Feld anpassen falls zu klein*****
	if (num + 2 > nummax) {
		nummax = nummax < nummatch ? nummatch : nummax * 2;
		resize(nummax);
		resizei = 0;
	}

	double closeval = a0*0.5; //Minimaler Abstand zu anderen FL
	bool close = false;

	//Wähle zufällige Position**********
	if (xadd == 0 && yadd == 0) {
		do {

			xadd = gsl_rng_uniform(rndm) * xl; // sehr grob...,
			// hoffentlich trifft es keine existierende FL
			yadd = gsl_rng_uniform(rndm) * yl;

			close = false;
			if (pinned == 1) {
				for (int i = 0; i<num; i++) { // zu nahe aneinander?
					double xx = fabs(xadd - xy(2 * i));
					if (xx>xh) xx = xl - xx;
					double yy = fabs(yadd - xy(2 * i + 1));
					if (yy > yh) yy = yl - yy;

					if (hypot(xx, yy) < closeval) {
						close = true;
						break;
					}
				}
			}
		} while (close); //Solange zu nahe
	}

	//Koordinaten eintragen***********
	xy(2 * num) = xadd;
	xy(1 + 2 * num) = yadd;
	if (pinned == 1) prop(num) = -1;
	num++;
	cout << "xy+1 (" << num << ") ";
	return resizei;
}

int vortex_field_direct::fill(int addn, double xfrac, double yfrac)
// füllt in Bereich 0..xfrac*xl und 0..yfrac*yl
// addn FL-n ein, oder bis num==nummatch
{
	int end = addn == 0 ? nummatch : num + addn;

	//Feld anpassen falls addn zu groß*****
	if (end > nummax) {
		nummax = end;
		resize(nummax);
	}

	//Fügt FL bis Ende auf******************
	while (num < end) {
		xy(2 * num) = (drand48() - 0.5)*xl*xfrac - xh;     //Gerade mit mit zufälliger x Koordinaten rund um (-xh) Wird mit checkbounds in richtigen Bereich gebracht
		xy(1 + 2 * num++) = (drand48() - 0.5)*yl*yfrac - yh; //Ungerade mit y Koordinaten
	}
	check_bounds();
	cout << "filled to " << end << endl;
	return addn;
}

void vortex_field_direct::swapfl(const int& fl1, const int& fl2)
{
	swap(oE(fl1), oE(fl2));
	swap(xy(fl1 * 2), xy(fl2 * 2));
	swap(xy(fl1 * 2 + 1), xy(fl2 * 2 + 1));
	swap(fxy(fl1 * 2), fxy(fl2 * 2));
	swap(prop(fl1), prop(fl2));
	return;
}

int vortex_field_direct::pinn_rand(const int pinn_num)
{
	int i;
	int pinnfl = 0;
	int cntdoublepinn = 0; //Sicherheit um endlosschleife im fall zu vieler pinned FL zu vermeiden
	for (i = 0; i < pinn_num; i++){
		pinnfl = lrint(drand48()*(num - 1));
		if (prop(pinnfl) < 0){
			cntdoublepinn++;
			if (cntdoublepinn > 10000){
				cout << "zu viele pinned FL, Anzahl neugepinnt:" << i << endl;
				return -1;
			}
		}
		else {
			prop(pinnfl) = -1;
		}
	}
	cout << "FL Pinned: " << i << endl;
	pinnedtostart();
	return 0;
}

int vortex_field_direct::pinnedtostart()
{
	unsigned int ipinned = 0;
  while (prop(ipinned)<0) ipinned++;
	for (int i = ipinned+1; i < num; i++){

		if (prop(i) < 0 ) {
			swapfl(ipinned, i);
			ipinned++;
		}
	}
	return 0;
}

void vortex_field_direct::setNpinned(int p)
{
	for (int i = 0; i < p; i++) prop(i) = -1;
	return;
}

int vortex_field_direct::getNpinned()
{
	int n = 0;
	for (int i = 0; i < num; i++) if (prop(i) < 0) n++;
	return n;
}

void vortex_field_direct::rm_pin()
{
	int i;

	for (i = 0; i < num; i++)
		if (prop(i) == -2)
		{
		prop(i) = 0;
    // FMSrem:		pin2.push_back(i); // fehlt: Sollte aus den gepinnten (am Anfang) entfernt werden
		}

	return;
}

void vortex_field_direct::out_debug(int subn)
{
	double Ex;
	double sum0, delta;
	#ifdef ALLOW_EXT_FORCES
		double sumf = 0;
	#endif
	int i = 0, ix, iy;;

	fxy.fill(0); // fxy muss auf 0 gesetzt werden, bevor f_analyze mit calcforce=true aufgerufen wird, sonst steht ein Blödsinn in fxy

	ostringstream outfilename("");
	outfilename << base_out << "_detail_" << setfill('0') << setw(5) << subn << "_.dat";
	std::ofstream testout(outfilename.str().c_str());

	testout << setprecision(11);

	for (sum0 = 0, i = 0; i < num; i++) { // lokale dE nach gesamter Verschiebung
		ix = 2 * i; iy = ix + 1;

		Ex = f_analyze(i, true, false); //Energie nach Verschiebung
		sum0 += Ex - oE(i);
		#ifdef ALLOW_EXT_FORCES
			sumf += abs(dxy(ix)*fxy(ix)) + abs(dxy(iy)*fxy(iy)); //Energie korrigieren
		#endif
			delta = Ex - oE(i) - dxy(ix)*fext_xy(ix) - dxy(iy)*fext_xy(iy);
		testout << xy(ix) << ' ' << xy(iy) << ' '
			<< fxy(ix) << ' ' << fxy(iy) << ' '
			<< Ex - oE(i) << ' '
			#ifdef ALLOW_EXT_FORCES
			<< -dxy(ix)*fext_xy(ix) - dxy(iy)*fext_xy(iy)
			#endif
			<< dxy(ix) << ' ' << dxy(iy) << ' ' << delta << ' ' << Ex << ' ' << prop(i) << std::endl;
	} // Ende Berechnung lokale dE mit Ausgabe Details
	testout.close();

	return;
}



/// Initialisiert das energy_file.
void vortex_field_direct::init_energy_file()
{
	energy_file.open("energy_graph.dat", ios::out);
	energy_file << setprecision(11);

	energy_file << "stepnr delta_E delta_e_ww delta_e_force e_ges e_ges_vorher ffacmed" << endl;
}

/// Speichert die Gesamtenergie und delta_E für jeden Schritt.
void vortex_field_direct::save_energy_graph(double stepnr, double delta_e, double delta_e_ww, double delta_e_force)
{

	if (!plot_energy) return;

	double eges = 0;
	double eges_vor = 0;

	for (int i = 0; i < num; i++)
	{
		eges += E(i);
		eges_vor += oE(i);
	}

	energy_file << stepnr << " " << delta_e << " " << delta_e_ww << " " << delta_e_force << " " << eges << " "
              << eges_vor << " " << ffac.mean() << endl;
}


bool  vortex_field_direct::pickup(string filnam, int n0)
// reads positions from old data-file, x=col(1), y=col(2) discards the rest
// n0 is the step number of the starting step

{
	ifstream infil;
	string istr;
	cout << "Versuche neues Positions-file einzulesen: "<< endl;
	infil.open(filnam.c_str(), ios::in);
	if (!infil) { cout << "no file " << filnam << "!\n" << endl; return false; }

	if (getline(infil, istr)) {
    cout << "First line: " << istr << endl;

    int i = 0;
    while(getline(infil, istr) && i<2*num) {
      istringstream doneline(istr);
      doneline >> xy(i++);
      doneline >> xy(i++);
    }

    if(infil.good()) {
      cout << "Read " << num << "data. Last line: " << istr << endl;
    } else {
      cout << "Nicht genug Datenpunkte in der Datei oder keine Schlusszeile"<< endl;
      cout << "Read " << i/2 << "data. Last line: " << istr << endl << endl;
    }

    if (num==i/2) {
      cout << "... ok" << endl << endl;
      istringstream doneline(istr);
      char c;
      double d0;
      int pinned0, num0;
      doneline >> c >> num0 >> d0 >> d0 >> pinned0;
      if (pinned0>num) {
        std::cerr << "schwerer Fehler! Mehr gepinnte als Flusslinien!"<< std::endl;
        return true; // ???? ok oder nicht?
        }
      std::cout << "number of pinned vortices: " << pinned0 << std::endl;
      setNpinned(pinned0);
      return true;
    }	else {
      cout << "gelesene Zahl: " << i/2 << " vs. die erwartete: " << num << endl;
      return false;
    }
  } else {
    cout << "konnte von Datei nicht lesen :( " << endl;
    return false;
  }

}

