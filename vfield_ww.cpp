#include "vfield.hpp"



void vortex_field::bess_generate(double bessmin)
{
	// jeweils ein zweidimensionales Feld mit Stützwerten für die Bessel-K0, K1-x und K1-y, die nicht nur den Abstand,
	// sondern auch die Wechselwirkung über die Periode hinweg beinhaltet

	ddx = 1e-2; // Bruchteil von lambda für die Schrittweite der Stützpunkte Soll: 1e-3, aber langsam...
	ddxh = ddx / 2;
	dd1x = 1. / ddx;

	double x = 0, y = 0;

	int ixnum = xh / ddx + 3;
	int iynum = yh / ddx + 2;

	// Tabellen:
	bessk0tbl.resize(ixnum, iynum);
	bessk1xtbl.resize(ixnum, iynum);
	bessk1ytbl.resize(ixnum, iynum);
	std::cerr << "bessk?tbl: " << ixnum << " x " << iynum << std::endl;
	std::cerr << "xh: " << xh << "   yh: " << yh << "   ddx: " << ddx << std::endl << std::endl;
	// std::cerr << "xh+ddx  = " << xh+ddx << " < " << ixnum*ddx << "?     yh+ddx =  "
	// << yh+ddx << " < " << iynum*ddx  << std::endl<< std::endl;

	double xr, r = 0;
	double corrr0 = 0, corr0 = 0;
	double corrr1, corrr1x = 0, corr1x = 0;
	double corrr1y = 0, corr1y = 0;
	int ir;

	for (ir = 1; ir < 100; ir++) { // Zahl der Ringe: 100 ist die absolute Schranke, sollte nie erreicht werden
		// (ev. bei Felder < 1T und zwei FL-n im Bereich
		corrr0 = 0;
		std::cerr << "0,0" << std::endl;
#pragma omp parallel for private(xr) reduction(+: corr0)
		for (int is = -ir; is < ir; is++) { // jeder Ring!
			xr = hypot(x - ir*xl, y + is*yl);
			corrr0 += gsl_sf_bessel_K0(xr);
			xr = hypot(x + is*xl, y + ir*yl);
			corrr0 += gsl_sf_bessel_K0(xr);
			xr = hypot(x + ir*xl, y - is*yl);
			corrr0 += gsl_sf_bessel_K0(xr);
			xr = hypot(x - is*xl, y - ir*yl);
			corrr0 += gsl_sf_bessel_K0(xr);
		}
		corr0 += corrr0; // Ringe zusammenzählen
		if (abs(corrr0) < bessmin) break;
	}

	bessk0tbl(0, 0) = corr0;
	bessk1xtbl(0, 0) = 0;
	bessk1ytbl(0, 0) = 0;
	cout << x << ' ' << y << ' ' << 0 << ' ' << corr0 << ' ' << ir
		<< std::endl;

	int ix, iy;
	double x0;
	std::cerr << "y: ";
	for (iy = 0; iy < iynum; iy++) {
		std::cerr << ' ' << iy;
		y = ddx*iy;
		if (1 == iy % 2) x0 = -ddxh; else x0 = 0; // Versetzen jeder zweiten Reihe um -0.5 ddx
		for (ix = 0; ix < ixnum; ix++) {
			x = x0 + ddx*ix;
			r = hypot(x, y);
			if (r < 0.05*a0) {
				if (ix != 0 || iy != 0) {
					bessk0tbl(ix, iy) = bessk0tbl(0, 0);
					bessk1xtbl(ix, iy) = 0;
					bessk1ytbl(ix, iy) = 0;
				}

     }
			else {
				corr0 = gsl_sf_bessel_K0(r);  // unkorrigierte Werte
				corr1x = gsl_sf_bessel_K1(r) / r; // unkorrigierte Werte
				corr1y = corr1x*y;  // unkorrigierte Werte
				corr1x *= x;  // unkorrigierte Werte


				double xx, yy;
				for (int ir = 1; ir < 100; ir++) {
					corrr0 = 0; corrr1x = 0; corrr1y = 0;
#pragma omp parallel for private(xx, yy, xr, corrr1) reduction(+: corrr0, corrr1x, corrr1y)
					for (int is = -ir; is < ir; is++) {

						xx = x - ir*xl, yy = y + is*yl;
						xr = hypot(xx, yy);
						corrr0 += gsl_sf_bessel_K0(xr);
						corrr1 = gsl_sf_bessel_K1(xr) / xr;
						corrr1x += corrr1*xx;
						corrr1y += corrr1*yy;
						// if(i==0) std::cerr << "corr0: " << ir << ' ' << is << ' ' << xx << ' ' << yy << ' ' << corrr0 << ' ' << corrr1x << ' ' << corrr1y << std::endl;

						yy = y + ir*yl; xx = x + is*xl;
						xr = hypot(xx, yy);
						corrr0 += gsl_sf_bessel_K0(xr);
						corrr1 = gsl_sf_bessel_K1(xr) / xr;
						corrr1x += corrr1*xx;
						corrr1y += corrr1*yy;
						// if(i==0) std::cerr << "corr1: " << ir << ' ' << is << ' '<< xx << ' ' << yy << ' ' << corrr0 << ' ' << corrr1x << ' ' << corrr1y << std::endl;

						xx = x + ir*xl; yy = y - is*yl;
						xr = hypot(xx, yy);
						corrr0 += gsl_sf_bessel_K0(xr);
						corrr1 = gsl_sf_bessel_K1(xr) / xr;
						corrr1x += corrr1*xx;
						corrr1y += corrr1*yy;
						// if(i==0) std::cerr << "corr2: " << ir << ' ' << is << ' '<< xx << ' ' << yy << ' ' << corrr0 << ' ' << corrr1x << ' ' << corrr1y << std::endl;

						yy = y - ir*yl; xx = x - is*xl;
						xr = hypot(xx, yy);
						corrr0 += gsl_sf_bessel_K0(xr);
						corrr1 = gsl_sf_bessel_K1(xr) / xr;
						corrr1x += corrr1*xx;
						corrr1y += corrr1*yy;
						// if(i==0) std::cerr << "corr3: " << ir << ' ' << is << ' '<< xx << ' ' << yy << ' ' << corrr0 << ' ' << corrr1x << ' ' << corrr1y << std::endl;
					}
					corr0 += corrr0;
					corr1x += corrr1x;
					corr1y += corrr1y;
					if (abs(corrr0) < bessmin) break;
				}

				bessk0tbl(ix, iy) = corr0;
				bessk1xtbl(ix, iy) = corr1x;
				bessk1ytbl(ix, iy) = corr1y;
				// out_any << x << ' ' << y << ' ' <<  gsl_sf_bessel_K0(r)
				// 	<< ' ' << corr0 << ' ' << ir
				// 	<< std::endl;
        }
		}
	}

	std::cerr << "\nWechselwirkungsenergie (0,0): " << bessk0tbl(0, 0)
		<< "  ir= " << ir - 1 << std::endl;
	std::cerr << "Wechselwirkungsenergie (a0): " << gsl_sf_bessel_K0(a0) << std::endl;
	std::cerr << "\nTabelle (ixnum/2,iynum/2): " << bessk0tbl(ixnum / 2, iynum / 2)
		<< " statt " << gsl_sf_bessel_K0(ddx*hypot(ixnum / 2, iynum / 2)) << std::endl;
	// Testausgabe für die ergänzte Wechselwirkung
	// out_any.open("bessk0idx.dat");

	// for(ix=0; ix<ixnum; ix++) {
	//   x=ix*ddx;
	//   if (x>0.5*a0) iy=0;
	//   else iy=sqrt(a0*a0*0.25-x*x)/ddx+3;
	//   for(; iy<iynum; iy++) {
	//     y = iy*ddx;
	//     r = hypot(x,y);
	//     // out_any << x << ' ' << y << ' '
	//     // 	      << 0.25*(bessk0tbl(ix-1,iy-1)+bessk0tbl(ix+1,iy-1)
	//     // 		       +bessk0tbl(ix+1,iy+1)+bessk0tbl(ix-1,iy+1))
	//     // 	- bessk0tbl(ix,iy)
	//     // 	      <<  ' ' << 0.25*(bessk1xtbl(ix-1,iy-1)+bessk1xtbl(ix+1,iy-1)
	//     // 		       +bessk1xtbl(ix+1,iy+1)+bessk1xtbl(ix-1,iy+1))
	//     // 	-  bessk1xtbl(ix,iy)
	//     // 	      <<  ' ' << 0.25*(bessk1ytbl(ix-1,iy-1)+bessk1ytbl(ix+1,iy-1)
	//     // 		       +bessk1ytbl(ix+1,iy+1)+bessk1ytbl(ix-1,iy+1))
	//     // 	-  bessk1ytbl(ix,iy)
	//           out_any << ix << ' ' << iy << ' ' << x << ' ' << y << ' '
	// 		    <<  bessk0tbl(ix,iy)
	// 		    <<  ' '  << bessk1xtbl(ix,iy)
	// 		    <<  ' ' << bessk1ytbl(ix,iy)
	// 		    << std::endl;
	//   }
	// }

	// out_any.close();
	// exit(99);

	return;
}

double vortex_field_direct::calcWW(int pnr, int i, bool force_calc)
// Oleg: this funktion replaces bessk0_tab() and both calcWW_k#().
//
// ACHTUNG: diese Funktion ist speziell für die openMP in f_analyze() angepasst
// berechnet K0 und K1 Summen direkt, keine Tabellen nötig!

{
	int ix, iy, ipx, ipy;
	double xip, yip;


	ix = 2 * i; iy = ix + 1;
	ipx = 2 * pnr; ipy = ipx + 1;

	xip = xy(ix) - xy(ipx); // xip und yip steht für den Vektor von p nach i
	yip = xy(iy) - xy(ipy);

	if (xip > xh) xip -= xl;    // muss translatierte FL (Klon) nehmen
	else if (xip<-xh) xip += xl;
	if (yip>yh) yip -= yl;    // muss translatierte FL (Klon) nehmen
	else if (yip < -yh) yip += yl;

	int fx = 1, fy = 1;
	if (xip < 0) { fx = -1; xip = -xip; }
	if (yip < 0) { fy = -1; yip = -yip; }
	// Differenzvektor im ersten Quadranten, damit die Interpolation funktioniert
	// , Vorzeichen gemerkt in fx, fy

	double bk0, bkx, bky; // Energie, Kraft in x-, Kraft in y-Richtung

	double r = hypot(xip, yip);
	double corr0 = gsl_sf_bessel_K0(r);  // unkorrigierte Werte
	double corr1x = gsl_sf_bessel_K1(r) / r; // unkorrigierte Werte
	double corr1y = corr1x*yip;  // unkorrigierte Werte
	corr1x *= xip;  // unkorrigierte Werte

	// if(i==0) std::cerr << pnr << ' ' << xip << ' ' << yip << ' ' << corr0 << ' ' << corr1x << ' ' << corr1y << ' ';



	for (int ir = 1; ir < 100; ir++) {
		double corrr0 = 0, corrr1x = 0, corrr1y = 0, corrr1;
		for (int is = -ir; is < ir; is++) {

			double xx = xip - ir*xl, yy = yip + is*yl;
			double xr = hypot(xx, yy);
			corrr0 += gsl_sf_bessel_K0(xr);
			corrr1 = gsl_sf_bessel_K1(xr) / xr;
			corrr1x += corrr1*xx;
			corrr1y += corrr1*yy;
			// if(i==0) std::cerr << "corr0: " << ir << ' ' << is << ' ' << xx << ' ' << yy << ' ' << corrr0 << ' ' << corrr1x << ' ' << corrr1y << std::endl;

			yy = yip + ir*yl; xx = xip + is*xl;
			xr = hypot(xx, yy);
			corrr0 += gsl_sf_bessel_K0(xr);
			corrr1 = gsl_sf_bessel_K1(xr) / xr;
			corrr1x += corrr1*xx;
			corrr1y += corrr1*yy;
			// if(i==0) std::cerr << "corr1: " << ir << ' ' << is << ' '<< xx << ' ' << yy << ' ' << corrr0 << ' ' << corrr1x << ' ' << corrr1y << std::endl;

			xx = xip + ir*xl; yy = yip - is*yl;
			xr = hypot(xx, yy);
			corrr0 += gsl_sf_bessel_K0(xr);
			corrr1 = gsl_sf_bessel_K1(xr) / xr;
			corrr1x += corrr1*xx;
			corrr1y += corrr1*yy;
			// if(i==0) std::cerr << "corr2: " << ir << ' ' << is << ' '<< xx << ' ' << yy << ' ' << corrr0 << ' ' << corrr1x << ' ' << corrr1y << std::endl;

			yy = yip - ir*yl; xx = xip - is*xl;
			xr = hypot(xx, yy);
			corrr0 += gsl_sf_bessel_K0(xr);
			corrr1 = gsl_sf_bessel_K1(xr) / xr;
			corrr1x += corrr1*xx;
			corrr1y += corrr1*yy;
			// if(i==0) std::cerr << "corr3: " << ir << ' ' << is << ' '<< xx << ' ' << yy << ' ' << corrr0 << ' ' << corrr1x << ' ' << corrr1y << std::endl;
		}
		corr0 += corrr0;
		corr1x += corrr1x;
		corr1y += corrr1y;
		if (abs(corrr0) < 1e-16) break;
	}
	//  if(i==0) std::cerr  << corr0 << ' ' << corr1x << ' ' << corr1y << std::endl;


	bkx = fx*corr1x;  // Zuweisungen für direkte Berechnung
	bky = fy*corr1y;
	bk0 = corr0;

	if (force_calc) {
		//Kräfte speichern. Gleich für beide Punkte i und p
#pragma omp atomic
		fxy(ix) += bkx;
#pragma omp atomic
		fxy(iy) += bky;
#pragma omp atomic
		fxy(ipx) -= bkx;

#pragma omp atomic
		fxy(ipy) -= bky;


	}

	return(bk0);
}


double vortex_field::calcWW(int pnr, int i, bool force_calc)
// Oleg: this funktion replaces bessk0_tab() and both calcWW_k#().
// ACHTUNG: diese Funktion ist speziel für die openMP in f_analyze() angepasst

{
	int ix, iy, ipx, ipy;
	double xip, yip;
	double g1anst = 0.5;

	// if (i==pnr) return(0); // bei berechnung der WW einer FL mit sich selber soll E=0 und die Kraft +=0 rauskommen
	//:: entfernt, weil das nicht vorkommen darf

	ix = 2 * i; iy = ix + 1;
	ipx = 2 * pnr; ipy = ipx + 1;

	xip = xy(ix) - xy(ipx); // xip und yip steht für den Vektor von p nach i
	yip = xy(iy) - xy(ipy);

	if (xip > xh) xip -= xl;    // muss translatierte FL (Klon) nehmen
	else if (xip<-xh) xip += xl;
	if (yip>yh) yip -= yl;    // muss translatierte FL (Klon) nehmen
	else if (yip < -yh) yip += yl;

	int fx = 1, fy = 1;
	if (xip < 0) { fx = -1; xip = -xip; }
	if (yip < 0) { fy = -1; yip = -yip; }
	// Differenzvektor im ersten Quadranten, damit die Interpolation funktioniert
	// , Vorzeichen gemerkt in fx, fy

	double bk0, bkx, bky; // Energie, Kraft in x-, Kraft in y-Richtung

	// Indizes im 2D Lookup-table:

	double rti0x, rti0y = yip*dd1x;  // Division für Index
	int ti0x, ti0y = rti0y;  // Index unterhalb des Wertes
	double cx, cy = rti0y - ti0y;  // Rest/ddx, gut für weitere Verarbeitung
	int ti1x, ti2x, ti2y; // (ti2x,ti2y) , (ti1x,ti1y) vervollständigen das Dreieck

	//bool down = false;
	if (cy > 0.5) { // FL näher zur oberen Reihe!
		cy = 1. - cy;
		ti2y = ti0y++;  // obere Reihe und nächste Nachbarn in der Reihe darunter
		// down = true;
	}
	else {
		ti2y = ti0y + 1;   // nächste Nachbarn in der Reihe darüber
	}


	if (0 == ti0y % 2) { // gerade Reihen, ohne x-Versetzung
		ti0x = rti0x = xip*dd1x;
		ti2x = ti0x + 1;
	}
	else {
		ti0x = rti0x = xip*dd1x + 0.5;  // ungerade Reihen, um -0.5 ddx versetzt
		ti2x = ti0x;
	}

	bool left = false;
	cx = rti0x - ti0x;
	if (cx > 0.5) { // näher zum nächsten Wert!
		cx = 1 - cx;
		ti1x = ti0x++; // der Wert vorher in der selben Reihe
		left = true;
	}
	else {
		ti1x = ti0x + 1; // der Wert danach in der selben Reihe
	}

	// jetzt sind (ti0x,ti0y) die Koordinaten des nächsten Punktes
	// und (tinx, tiny) die Koordinaten des passenden Punktes in der nächsten Reihe
	// und (ti1x,ti0y) der zugehörige Punkt in derselben Reihe

	// int wtch=0;

	double cx1 = cx - g1anst*cy;  // Restbewegung in x-Richtung nach vorhergehender in Richtung 2, um y zu erreichen
	//double x0, y0, x1, x2, y2;
	//double x=0, y=0;

	if (cx1 >= 0) {
		// ti1y wäre =ti0y ....
		// Interpolation:
		// Basis des Dreiecks auf der Reihe des nächsten Punktes:
		bk0 = bessk0tbl(ti0x, ti0y) + cy*(bessk0tbl(ti2x, ti2y) - bessk0tbl(ti0x, ti0y)) + cx1* (bessk0tbl(ti1x, ti0y) - bessk0tbl(ti0x, ti0y));
    if (force_calc) {
      bkx = bessk1xtbl(ti0x, ti0y) + cy*(bessk1xtbl(ti2x, ti2y) - bessk1xtbl(ti0x, ti0y)) + cx1* (bessk1xtbl(ti1x, ti0y) - bessk1xtbl(ti0x, ti0y));
      bky = bessk1ytbl(ti0x, ti0y) + cy*(bessk1ytbl(ti2x, ti2y) - bessk1ytbl(ti0x, ti0y)) + cx1* (bessk1ytbl(ti1x, ti0y) - bessk1ytbl(ti0x, ti0y));
		}
	}	else {

		// der nächste Punkt ist die Spitze des Dreieckes
		ti1x = ti2x + (left ? 1 : -1); // ti1y wäre =ti2y....

		bk0 = bessk0tbl(ti0x, ti0y) + cy*(bessk0tbl(ti2x, ti2y) - bessk0tbl(ti0x, ti0y)) + cx1* (bessk0tbl(ti2x, ti2y) - bessk0tbl(ti1x, ti2y));
    if (force_calc) {
      bkx = bessk1xtbl(ti0x, ti0y) + cy*(bessk1xtbl(ti2x, ti2y) - bessk1xtbl(ti0x, ti0y)) + cx1* (bessk1xtbl(ti2x, ti2y) - bessk1xtbl(ti1x, ti2y));
      bky = bessk1ytbl(ti0x, ti0y) + cy*(bessk1ytbl(ti2x, ti2y) - bessk1ytbl(ti0x, ti0y)) + cx1* (bessk1ytbl(ti2x, ti2y) - bessk1ytbl(ti1x, ti2y));
    }
		// Vielleicht numerisch besser mit den beiden Vektoren vom Punkt?
	}

  // cout << "calcWW " << pnr << ':' << i << ' '
  //      << bk0 << ' ' << bkx << ' ' << bky << std::endl;


	if (force_calc) {
    if (fx < 0) bkx = -bkx;
    if (fy < 0) bky = -bky;
		//Kräfte speichern. Gleich für beide Punkte i und p
#pragma omp atomic
		fxy(ix) += bkx;
#pragma omp atomic
		fxy(iy) += bky;

#pragma omp atomic
		fxy(ipx) -= bkx;
#pragma omp atomic
		fxy(ipy) -= bky;
	}

	return(bk0);
}

