#include <cmath>
#include <iostream>
#include <fstream>
#include "vfield.hpp"
#include <sstream>
// #include <algorithm>
// algorithm für max()

#define FACSTEP 5
#define lambda __LAMBDA // 1.4e-7
#define signum(x) ((x) >= 0 ? 1 : (-1))

// Andreas Burger: ad_analyze() für zwischenf aufgerufen


double vortex_field_direct::point_analyze(double x, double y)
{
	// Oleg: Es wird nur die Energie der WW berechnet für einen beliebigen Punkt (x,y)
	// Man kann nicht von (x,y) auf (i,pnr) übergehen, also wiederholen wir das f_analyze mit abhängigkeit von (x,y)

	// BROKEN! wait for corrected calcWW()!!!!!!!!!!!!!!!!!!!!!!! FMS

	double ix, iy;
	double Ex=0;
	// int k=bessk0_num, i;
	int i;


	for (i = 0; i < num; i++) {
		ix = 2 * i; iy = ix + 1;
		ix = xy(ix) - x;
		iy = xy(iy) - y;
		if ((ix = fabs(ix)) > xh) ix -= xl; // muss translatierte FL (Klon) nehmen
		if ((iy = fabs(iy)) > yh) iy -= yl; // muss translatierte FL (Klon) nehmen

		// look up in the bessk0 table (replacement for  bessk0_tab())
		// while(bessk0start(k)>r2) //TH: Bereich für r2 Suchen
		//   k--;

		// c=((r2-bessk0start(k))/bessk0step(k));
		// n=lrint(c);        //TH: Zahl der Schritte vom Start der Tabellarisierungs Bereiches gezählt.
		// c-=n;              //TH:c=Rundungsfehler
		// n+=bessk0n(k);     //Erstes Element im Bereich
		// bn0 = bessk0tbl(n); //Damit man nicht doppelt in der Tabelle nachschlagt

		// if ((c)<0)         //TH: Wenn aufgerundet linearisieren zwischen vorhergehendem Stützpunkt
		//   bk0 = (bn0 + (bn0 - bessk0tbl(n-1)*c));
		// else               //TH: Sonst mit nächstem
		//   bk0 = (bn0 + (bn0 - bessk0tbl(n+1)*c));
		// // end of bessk0_tab()

		// wwsum.E += bk0;
	}
	return Ex;
}


void vortex_field_direct::f_analyze_setup(double force, bool force_calc)
{
  ofstream testout;
  unsigned int curr=0;
  ostringstream outfil; 
  icurr=0; iold=1;
  oxy=xy;  // bei step-back passiert sonst ein Fehler
  
  // std::cout << "f_analyze_setup Anfang: ";
  // std::cout << num << " last: " << xy(2*num-2) << ", " << xy(2*num-1)
  //           << " auf: " << &(xy(2*num-2)) << ", " << &(xy(2*num-1)) 
  //           <<  " von: " << &(xy(0)) << std::endl;

	E.fill(0); fxy.fill(0); fext_xy.fill(0); // Alle Energien und Kräfte zurücksetzen
  ofxy.fill(0); ofval.fill(1);
  ffac.fill(ffac0);
  std::cout << ffac(0) << ' ' << ffac(1) << ' ' << ffac(num-1) << std::endl;
  
  minffac = ffac0*.05;
  maxffac = ffac0*20.;
  maxdx = 0.05;

  std::cout << "f_analyze_setup() maxdx =  " << maxdx <<
    " ffac=" << ffac.mean() << " ("<< ffac0  << ") " << std::endl;

  calc_set_forces(); //  bei jedem Schritt in relax() inkludierten Kraftberechnung
		// (->David Bader im subdir analyze)

  for (int i = 0; i < num; i++) { // berechnen von Energie und WW-Kräften inklusive externe Kräfte
		f_analyze(i, force_calc, true);
    int ix = i * 2; 
		fval(i) = hypot(fxy(ix), fxy(ix+1));
	}
  // std::cout << "nach f_analyze: ";
  // std::cout << num << " last: " << xy(2*num-2) << ", " << xy(2*num-1) << std::endl;

  do {
    std::cerr << "fxy: " << fxy.size() << "   dxy: " << dxy.size() << "      num: " << num << std::endl;
    dxy = fxy*ffac0*0.1; // damit der erste step funktioniert
    while(dxy.head(2*num).maxCoeff()>maxdx || dxy.head(2*num).minCoeff()<-maxdx) {
      std::cerr << "fd "<< dxy.maxCoeff() << ' ' << dxy.minCoeff() << std::endl;
      dxy*=.3;
    }
    // dxy_shift(testout);
    outfil1.str("");
    outfil1 << "ffset_" << curr++ << ".dat";
    // std::cerr << "writing "<< outfil1.str() << std::endl;
    ww_save(outfil1.str());
  } while (relax_step(force)==0); // muss neu, falls dn+rev in relax_step() 
  return;
}


double vortex_field_direct::f_analyze(int pnr, bool force_calc, bool recalc)
	// Oleg: f_analyze für Berechnung der WW-Energie und des WW-Kraftvektors auf den Punkt pnr
{
  // !!!!!!!
  //Achtung: Setzt voraus, dass funktion mit Parameter pnr nacheineinder von 0 bis num aufgerufen wird!
  // !!!!!!!
  
  // Ei ist die WW-Energie zw Punkten i und p
	// xi und yi definieren die WW-Kraft von i auf p

	//bool moved=(prop(pnr)==0);

	//TH: Nur Variablen für Parallelisierung (Struct nicht möglich!)
	double E00;
#pragma omp parallel for  private(E00)
	for (int i = pnr + 1; i < num; i++) {  //  über alle Punkte >pnr ... , die darunter sind schon gerechnet
		
    E00 = calcWW(pnr, i, force_calc);
    // hier werden auch die Kräfte ins fxy(2*i) fxy(2*i+1, und fxy(2*pnr...) gespeichert, wenn force_calc
    // if class vortex_field lookup in 2d-table, in class vortex_field_direct WW is directly calculated
				
    // if (E00 > 0) {	// todo: für was ist die Abfrage gut??
#pragma omp atomic
    E(pnr) += E00; // addiert für Punkte von 0 bis pnr
#pragma omp atomic
    E(i) += E00;
	}

	//	Jetzt ist pnr wirklich fertig und die externe Kraft kann eingerechnet werden!
	#ifdef ALLOW_EXT_FORCES
		if (force_calc)
		{
			int ipx = 2 * pnr;
			int ipy = ipx + 1;
			double xpos = xy(ipx);
			double ypos = xy(ipy);

			double fx = fxy(ipx);
			double fy = fxy(ipy);
      

			fext_xy(ipx) = force_x(xpos, ypos, fx, fy, prop(pnr));
			fext_xy(ipy) = force_y(xpos, ypos, fx, fy, prop(pnr));

			fxy(ipx) += fext_xy(ipx) - fpinx; // Für Rest-gepinntes korrigieren...  Test FMS
			fxy(ipy) += fext_xy(ipy) - fpiny;

      // std::cerr << ' ' << fext_xy(ipx) << '_'; // check ob erreicht

		}
#endif

	return E(pnr);
}


int vortex_field_direct::check_bounds(void)
{
	int ipt = 0;
	int icnt = 0;

	#ifdef RAND_FEST
		
		double alt;

		for (int i = 0; i < num; i++) {
			ipt = 2 * i;
			if (xy(ipt) < -xh) 
			{
				alt = xy(ipt) - dxy(ipt);
				xy(ipt) = -xh;
				dxy(ipt) = xy(ipt) - alt;
				// alt = dxy(ipt);	// Zum Testen
			}
			if (xy(ipt) > xh)
			{
				alt = xy(ipt) - dxy(ipt);
				xy(ipt) = xh;
				dxy(ipt) = xy(ipt) - alt;
				// alt = dxy(ipt);	// Zum Testen
			}

			ipt++;

			if (xy(ipt) < -yh)
			{
				alt = xy(ipt) - dxy(ipt);
				xy(ipt) = -yh;
				dxy(ipt) = xy(ipt) - alt;
				// alt = dxy(ipt);	// Zun Testen
			}
			if (xy(ipt) > yh)
			{
				alt = xy(ipt) - dxy(ipt);
				xy(ipt) = yh;
				dxy(ipt) = xy(ipt) - alt;
				// alt = dxy(ipt);	// Zun Testen
			}
		}
	#else
		#pragma omp parallel for private(ipt) reduction(+: icnt)
		for (int i=0; i<num; i++) {
			ipt=2*i;
			while (xy(ipt) < -xh) {xy(ipt)+=xl;icnt++;cout<<"x";}
			while (xy(ipt) > xh) {xy(ipt)-=xl;icnt++;cout<<"X";}
			ipt++;
			while (xy(ipt) < -yh) {xy(ipt)+=yl;icnt++;cout<<"y";}
			while (xy(ipt) > yh) {xy(ipt)-=yl;icnt++;cout<<"Y";}
		}
	#endif


	return icnt;
}



double vortex_field_direct::getE(void)
// berechnet die Gesamtenergie des Feldes
{
	int i;
	double sumE = 0.0;
	double get_E;

	E.fill(0);  //die E soll auf 0 gesetzt werden, sonst wird E akkumuliert

	for (i = 0; i < num; i++) {
		get_E = f_analyze(i, false, false);
		sumE += get_E;
	}
	return sumE;
}


double vortex_field_direct::dxy_shift(std::ofstream &testout)
// um dxy() relativ verschoben
// Gesamtverschiebung seit oE()wird mit dem Schritt-count sv gezählt
// und testout nur für die Ausgabe...
{
	// int i, ix, iy;

	// //  #pragma omp parallel for 
	// for (int i = 0; i < num; i++) { // Verschieben um dxy in Kraft-Richtung
	// 	int ix = 2 * i; int iy = ix + 1;
	// 	if (prop(i) >= 0) {      // falls nicht verschiebbar, nicht verschieben....
	// 		xy(ix) += dxy(ix);
	// 		xy(iy) += dxy(iy);
	// 	}
	// }

  // tausche alt gegen neu
  int irem= iold;
  iold = icurr; // alt-neu vertauschen nach der Verschiebung 
  icurr = irem; // icurr und iold wählen aus fxy_base und fval_base den richtigen Index aus

  // neues xy aus altem+ altem dxy 
  xy = oxy + odxy;  // prop schon in relax_step() bei dxy-Berechnung berücksichtigt.
  
	check_bounds();

  double dE_sum = 0.0;
	double dE_force_sum = 0;
	double dE_WW_sum = 0;

  double dE_WW, dE_force;
  
	E.fill(0); fxy.fill(0); fext_xy.fill(0); // Alle Energien und Kräfte zurücksetzen

	for (int i = 0; i < num; i++) {      // berechne Energie-gewinn nach gesamter Verschiebung
		// int ix = 2 * i; int iy = ix + 1;
		f_analyze(i, true, false); // neue E und f nach Verschiebung,
                               // wird im f_analyze() parallelisiert....

    dE_WW = E(i) - oE(i);
		dE_WW_sum += dE_WW;	// stimmt jetzt! oE und E sind ohne ext. Kräfte

#ifdef ALLOW_EXT_FORCES		//Energie korrigieren mit zusätzl. Kraft
    int ix = 2 * i; int iy = ix + 1;
    dE_force =  - dxy(ix)*fext_xy(ix) - dxy(iy)*fext_xy(iy);
    dE_force_sum -= dE_force;
    dE_sum -= dE_force;
#endif

		dE_sum += dE_WW;  // Achtung: negatives deltaE bedeutet das E kleiner geworden ist (weiter gleich so)
	}

	if (testout)  testout << ffac0 << ' ' << (dE_sum / num) << std::endl;

	save_energy_graph(active_step_nr, dE_sum, dE_WW_sum, dE_force_sum);
	active_step_nr += 0.2;

	return dE_sum;
}



int vortex_field_direct::relax_step(double force)
// definiert einen Relaxationsschritt für das Feld
{

	ostringstream outfilename("");
	outfilename << base_out << "_ffacdElarge_" << setfill('0') << setw(5) << stepn << "_.dat";
	std::ofstream testout;

	int i, ix, iy;          // temporary variables for cycles
	double  dE_sum = 0.0;

	
  // ****************
  dE_sum = dxy_shift(testout); // verschiebt um dxy und gibt die Änderung bez. oE() zurück.
                               // wechselt auch alte nach aktuelle: xy, fxy, dxy, fval
  // ****************

  // Jetzt wird der Step ausgewertet und weiter verarbeitet:
  
  double ffac_mean=ffac.mean();
  std::cerr << "dE_sum: " << dE_sum << " (" << ffac_mean << ") " << std::endl;
 
  max_force = 0;
  // double cosa; // removed, because vector for monitoring in testing
	//#pragma omp parallel for private (ix,iy)
	for (i = 0; i < num; i++) {     // max Kraft und Betrag der neuen Kraft bestimmen
                                  // , ffac für Umrechnung in dxy

		ix = i * 2; iy = ix + 1;
		fval(i) = hypot(fxy(ix), fxy(iy)); // Länge der Kraft berechnen

		if (prop(i) >= 0 && max_force < fval(i)) max_force = fval(i); // max_force beim Parallelisieren?
    // maxforce nur für die Flusslinen berechnen, die verschiebbar sind.
	}

  double mean_force=fval.mean(); // not really necessary, Ausgabe und f_irrelev, siehe unten
  
  std::cerr << "f: " << "<" << mean_force << "> "<< max_force << "<=? " << force << std::endl;
  double f_irrelev=force*.001;  // minimum relevant force, ?0.1? or tie to 0.01*max_force oder mean_force?
  // FMS: more aggressive set, old one in comments, brackets
  static const double upbnd= 0.8; // boundary for increase of ffac, ALL above! (.99)
  static const double upfac= 1.05;  // increase grow step when in the same direction
  static const double downbnd= 0.5;  // boundary for decrease of ffac  (.8)
  static const double downfac= 0.95;  // decrease grow step ffac
  static const double revbnd =0.1;  // boundary for revert step (.4) 
  static const double revfac =0.5;  // defrease ffac in case of revert step

  double mincosa=1.;  // Find min of cosa, necessary for upstep and revert 
  double ffac_lowest=ffac0*0.001; // no decrease of ffac below this value
  double ffac_highest=ffac0*1000; // no decrease of ffac below this value
  
  int dnfrac = 0; // zählt, ob retour notwendig
  int dncnt = 0;  // Zähler für down
  static int revcnt=0;  // zählt wie oft hintereinander rev war 
  static int urevcnt=0;  // zählt wie oft hintereinander rev mit unify war 
  
  int movcnt =0;
	if (max_force > force) { // Kräfte behandeln, falls signifikant
                           // Setup für neue Bewegung, oder retour
    dnfrac=0;

    unsigned int ipin=0;
    fpinx=0, fpiny=0, ipin=0;  // Summen über gepinnte FL-n initialisieren
    for (i = 0; i < num; i++) {     // FL nur verschieben
      // , wenn die Kraft größer als die Pinning-Kraft ist
      ix = i * 2; iy = ix + 1;
      cosa(i) = (ofxy(ix)*fxy(ix) + ofxy(iy)*fxy(iy)) / (ofval(i)*fval(i)); // detect nonlinearity
      if (fval(i) < force || prop(i)<0 ) {  // force too small or pinned
        dxy(ix) = 0;
        dxy(iy) = 0;
        // if (prop(i) == -1) {  // remove fpinx....
        //   ipin++;  // Zähler für gepinnte im Aussenbereich
        //   fpinx += fxy(ix); // Kraftsumme über die gepinnten, _nur_ wenn nicht im Innenbereich
        //   fpiny += fxy(iy);
        // }
      } else {
        if (fval(i)>f_irrelev) { // Herumspielen mit relevanten Kräften ...&& fval(i)<fbreak
                                      // Kernbegrenzung f_irrelev verhindert Steckenbleiben(, _aber_ ...???)  

          if (cosa(i)<mincosa) mincosa=cosa(i);
          if (cosa(i) < downbnd) {  // Achtung, kurvig, slow down
            if (ffac(i)>ffac_lowest) {
              ffac(i) *=downfac;  // if() to avoid getting unrecoverably low
              std::cerr << 'd';
            } else {
              std::cerr << 'D';
              dnfrac=1;  //single slow down not good, too far out of others, general slowdown
            }
            dncnt++;
            // if(cosa(i) < 0.3) { // change too large, slow down all  
            //   std::cerr << 'R';
            //   dncnt+=5;  // Kurve zu groß
            // }
          }
          dxy(ix) = fxy(ix)* ffac(i);
          dxy(iy) = fxy(iy)* ffac(i);
          movcnt++;
        } else { // fval<irrelev, ohne Zählung der +/- cosa
          std::cout << 'i';
          dxy(ix) = 0; //fxy(ix)* ffac(i);
          dxy(iy) = 0; //fxy(iy)* ffac(i);
          movcnt++;
        } 
      } // Ende dxy berechnen
      while (dxy.head(2*num).minCoeff() < -maxdx || dxy.head(2*num).maxCoeff() >maxdx) {
        std::cerr << "fd "<< dxy.maxCoeff() << ' ' << dxy.minCoeff() << std::endl;
        ffac*=.3; dxy*=.3;
      }
        
    } // Ende checken ob Bewegung

    if (ipin!=0) {
      fpinx /= ipin;  // mittlere gepinnte Kraft 
      fpiny /= ipin;
      std::cout << "fpin (" << fpinx << ", " << fpiny <<", "<< ipin << ")"<<std::endl;
    } // else sind fpinx und fpiny sowieso 0...
    
    if (innen_force!=0) std::cout << "max. angepasste Kraft/innen_force(abs): "
      << max_force << "/" << innen_force << " = " << max_force/innen_force << std::endl;
    if (movcnt>0 && max_force>innen_force*0.3) { // Test Franz innen_force*0.3 !! max_force meist am Rand, wo's wurscht ist
      dnfrac=dncnt/movcnt;
      std::cerr << "  moved: " << movcnt << "(dn: "<< dncnt << ") ";
      
      if (mincosa>upbnd && ffac_mean<ffac_highest) {
        ffac *= upfac;
        std::cerr << " up"<< endl;
        urevcnt=revcnt=0;
        return 3;
      } else if(dnfrac<0.5) { // all is well, continue
        // nix tun, nächster step beginnt mit dxy_shift()
        std::cerr << " fwd"<< endl;
        second_try = false;
        urevcnt=revcnt=0;
        return 2;
      } else {  // dnfrac>0.5  , etwas muss passieren 
        ffac *= downfac;
        std::cerr << " dn ";
        
        if(mincosa < revbnd) {
          if(urevcnt<4 ) { // retour und neu ansetzen mit verringerter Verstärkung
          int irem= iold;
          iold = icurr; // alt-neu vertauschen nach der Verschiebung 
          icurr = irem; // icurr und iold wählen aus fxy_base und fval_base den richtigen Index aus
        
          if (revcnt<3) {
            ffac *= revfac;  // decrease ohne unify global ffac
            // urevcnt = 0; // besser nicht zurücksetzen, sonst 
          } else {
            ffac.fill(ffac_mean*revfac);
            // zu oft ohne unify... jetzt devcrease mit unify
            std::cerr<< " U";
            revcnt=0; 
            urevcnt++;
          }
          for (i = 0; i < num; i++) { // adaptiertes dxy
            ix = i * 2; iy = ix + 1;
            if (prop(i)>=0 && fval(i)>force) {
              dxy(ix) = fxy(ix)*ffac(i);
              dxy(iy) = fxy(iy)*ffac(i);
            }
          }
        std::cerr << "rev "  << endl;
        revcnt++;
        return 1; // end revert
        } else {  // should revert, overridden
          urevcnt = revcnt = 0;
          std::cerr << "rev override";
          second_try = false;
          return 7;  // flag für override, zwecks Ausgabe

          }// end revert? or override
        }
        return 0; // dn ohne revert
      } // end up, fwd, dn, rev
      
      second_try = false;
      return 0;
    }
  }

  // nix gegangen, keine signifikante Kraft mehr da  
  std::cerr << " no significant forces left! " << ffac0 << "  dE_sum:  " << dE_sum << std::endl;
  if (second_try) {
    second_try = false;
    return -1;
  }	else {
    second_try = true;
    return 8;
  }
} // Ende relax_step


int vortex_field_direct::relax(double force, int actstep)
{

	//ostringstream outfil1;
	int ret;
	int currstep = actstep;
	unsigned int tries = 0;

	outfil1.str("");

	cout << "FL_num: " << num 
		// << "    fprefac: " << fprefac 
       << setprecision(9) << "    av_E: " << getE() / num << endl;

	//int i; for (i=0;i<num;i++) {if (prop(i)<0) cout<<"prop("<<i<<"): "<<prop(i)<<endl;}  //pinned check

	double forcex = force;  // bis hierher wird relaxiert...
	// clock_t start;
	// time_t startsum;
	// startsum = time(NULL);

	// Anfangswerte: (fac wird immer auf den aktuellen Wert gesetzt, wenn
	// die Minimumsuche in relax_step() durchgeführt wird. In regelmäßigen Abständen
	// FACSTEP wird die Minimumsuche erzwungen, gezählt durch faccnt)

	//fac =1; // in relax_step() ist die Auslenkung genau ins lokale Minimum
	faccnt = 0;
	second_try = false;
	recalc = false;

	// std::cout << "relax() vor analyze_setup: ";
  // std::cout << num << " last: " << xy(2*num-2) << ", " << xy(2*num-1) 
  //           << " auf: " << &(xy(2*num-2)) << ", " << &(xy(2*num-1))
  //           <<  " von: " << &(xy(0)) << std::endl << xy << std::endl;
  stepn = 0;
  f_analyze_setup(force, true); // -> innen_force_x, -y
  unsigned int secondcnt=0;
	do {
		// start = clock();

		cout << endl << "___stepnr " << stepn << "___" << endl;

    breaknumber =0;
		calc_set_forces(); // Plugin für bei jedem Schritt in relax() inkludierten Kraftberechnung
		// (->David Bader, Fabio)
    if(breaknumber == 1) break;
    
		active_step_nr = stepn + 1;	// damit ich im energy_file in dxy_shift weiß, bei welchem step ich bin.
		if (allout) {
		  //if (tries < 220 || (tries<5000 && ((tries%10==0)||tries%10==1)) || tries%100==0 || tries%100==1) {
      if(tries%25==0) {
			outfil1.str("");
						outfil1 << "zwischenf_" << currstep << "_" << tries << ".dat";
			// std::cerr << "writing "<< outfil1.str() << std::endl;
            ww_save(outfil1.str());            
            }
           
            if (tries%100==0) { ad_analyze(-1., -1., -1., -1., "zwischenf_"+to_string(currstep)+"_"+to_string(tries), force); } //call ad_analyze
            
            
      tries++;
      }
    if (innen_force_med!=0 && forcex!=innen_force_med) {
      std::cout << " ( diff = " << forcex-innen_force_med << ") ";
      forcex = innen_force_med;
      std::cout << "new forcex: " << forcex  <<  std::endl;
    }
		ret = relax_step(forcex);    // re-calc mit forcex als 
    std::cerr << " ** " << ret << " ** ";

		if (allout) {
      tries--;
		  if (ret==7) {  // rev override signal
        outfil1.str("");
        outfil1 << "overrideb_" << currstep << "_" << tries << ".dat";
        // std::cerr << "writing "<< outfil1.str() << std::endl;
        ww_save(outfil1.str());
      }
      tries++;
    }
    
		stepn++;

		// cout << "Schritt gestartet.\nDIAG: " << "n_vorher: " << sum_last << " n_jetzt: " << n_count << endl;
		// alles ausgegeben
    if(ret==8) {
      secondcnt++;
      if (secondcnt>6) {
        ret= -1;
        std::cerr << "ret 8 wiederholt. no movement..." << std::endl;
      }
    }

	} while (ret < 9 && ret >= 0); {
		/* TH:Solange failcnt< 9 bzw. bis keine Bewegung mehr (-1) und nicht abgebrochen */

		if (ww_max_dist < 10) ww_max_dist += 1.0;
		cout << "(" << stepn << ", " << ret << ")" << setprecision(7) << endl;  //<<"ww_max_dist = "<<ww_max_dist<<endl;

		//while ((ww_max_dist< (10))); // xh<yh?xh:yh*
		//cout<< "Fertig: (" << stepn << ", " << ret << ")" <<endl;
	}
	// fprefac *= 2;

  outfil1.str("");
  outfil1 << "final_0T6a_" << currstep << ".dat";
  // std::cerr << "writing "<< outfil1.str() << std::endl;
  ww_save(outfil1.str());

	return stepn;
}


int vortex_field_direct::pullout(const int mode)

{
	//Select FL
	double Ex;
	double max_E0 = 0;
	int xy_max = 0;

	double min_E0 = 0, y;
	double y_min = 0;

	const int ysteps = 30;
	double y_position = gsl_rng_uniform(rndm) * yl; //Zufällige Y-Position auswählen für mode=1.
	//FL mit max Feld im Bereich +/-(a0)
	// Test mit min...
	switch (mode){
		//Feldlinie mit maximaler Energie herausziehen:
	case 0:
		for (int i = 0; i < num; i++){
			if ((fabs(xy(i * 2)) < a0) && (prop(i) > -1)) {
				if (xy_max == 0)//Erste Feldlinie die in der Mitte gefunden wird zum vergleichen nehmen
					xy_max = i * 2;
				Ex = f_analyze(i, false, false);
				if (max_E0 == 0) max_E0 = Ex;
				if (Ex < max_E0){
					max_E0 = Ex;
					xy_max = i * 2;
				}
			}
		}

		min_E0 = 0;
		//FLPunkt mit min Feld bei xh-(3*a0/4)
		for (int i = -ysteps; i < ysteps; i++){
			y = (yh - a0)*i / ysteps;
			Ex = point_analyze(xh - 3 * a0 / 4, y);
			if (Ex < min_E0 || min_E0 == 0){
				min_E0 = Ex;
				y_min = y;
			}
		}
		break;

		//Zufällige Feldline aus Mitte herausziehen:
	case 1:

		for (int i = 0; i < num; i++){
			if ((fabs(xy(i * 2)) < a0) && (prop(i) > -1)) {
				if (xy_max == 0)//Erste Feldlinie die in der Mitte gefunden wird nehmen, damit sicher eine erwischt wird
					xy_max = i * 2;
				//Wenn unterhalb von y_position und näher bei y_position als letztes gefundenes
				if ((xy(i * 2 + 1) > y_position) && ((xy(xy_max) < xy(i * 2 + 1)) || (xy(xy_max) < y_position))){
					xy_max = i * 2;
				}
			}
		}


		min_E0 = 0;
		//FLPunkt mit min Feld bei xh-(3*a0/4)
		for (int i = -ysteps; i < ysteps; i++){
			y = (yh - a0)*i / ysteps;
			Ex = point_analyze(xh - 3 * a0 / 4, y);
			if (Ex < min_E0 || min_E0 == 0){
				min_E0 = Ex;
				y_min = y;
			}
		}
		break;

	default:
		cout << "Pullout(): Mode nicht definiert!";
	}

	cout << "max_E0:" << max_E0 << " minE0:" << min_E0 << endl;
	cout << "Pullout from (" << xy_max << "): x" << xy(xy_max) << " y" << xy(xy_max + 1);
	xy(xy_max) = xh - 3 * a0 / 4;
	xy(xy_max + 1) = y_min;
	check_bounds();
	cout << " to: x" << xy(xy_max) << " y" << xy(xy_max + 1) << endl;
	return 0;
}


int vortex_field_direct::verschiebe(const double x1, const double y1, const double x2, const double y2)
{
	/* naechste Flussline zu (x1, y1) suchen und die nach (x2, y2) verschieben */
	int xy_min;
	double d_min, d_neu;

	xy_min = 0;
	d_min = sqrt((xy(0) - x1)*(xy(0) - x1) + (xy(1) - y1)*(xy(1) - y1));
	for (int i = 1; i < num; i++)
	{
		d_neu = sqrt((xy(2 * i) - x1)*(xy(2 * i) - x1) + (xy(2 * i + 1) - y1)*(xy(2 * i + 1) - y1));
		if (d_neu < d_min)
		{
			d_min = d_neu;
			xy_min = i * 2;
		}
	}
	xy(xy_min) = x2;
	xy(xy_min + 1) = y2;

	return 0;
}


void vortex_field_direct::toggle_allout()
{
	if (allout == false)
	{
		allout = true;
	}
	else
	{
		allout = false;
	}
}

