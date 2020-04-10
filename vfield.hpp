#ifndef vfield_h
#define vfield_h

/* definiert die Klasse vortex_field und vortex_field_direct */
// vortex_field_direct: arbeitet mit direkter Berechnung, langsam...zum Testen der Tabelle
// vortex_field: davon abgeleitet, WW zwische Flusslinien ist 2-dimensional tabelliert bess_generate()

// Andreas Burger: ad_analyze() in vortex_field_direct, public deklariert 

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_rng.h>

#define r_h00 (sqrt(3.0)*0.5)
#define DISLOC_LEN 20
#define Fi0 2.07e-15
#define MIN_DIS 3.6e-15*16

#define FACdx 0.01

#define __LAMBDA 1.4e-7

using namespace std;

typedef struct { double xf, yf; } force_t;

// Um externe Kräfte zu erlauben, müssen diese mit
#define ALLOW_EXT_FORCES
// erlaubt werden. Sonst werden sie nicht einkompiliert.

// #define START_RANDOM	// Wenn aktiviert, wird mit einem zufälligem FL-Gitter gestartet, sonst mit einem regelmäßigen.

// #define RAND_GEGENKRAFT	// gibt an, dass im Randbereich eine Kraft proportional zur Kraft auf die Flusslinie und dem Abstand zum Rand wirkt.
// #define RAND_FEST	// gibt an, dass Flusslinien den Rand nicht durchdringen können.

//#define FORCE_INWARDS // Gibt an, dass die externe Kraft in der Mitte das Vorzeichen wechselt und somit nach innen zeigt.

//#define TEST_FORCE_X 1 // Kräfte zum Testen, werden zu force_x bzw. force_y dazugerechnet.
//#define TEST_FORCE_Y 0.5


class vortex_field_direct  // Flusslinien-Feld in periodischen Randbedingungen
{

private:
	// eventuell kommt noch was von protected hier herein...

protected:

	int num;             // Anzahl der FL-n
	long wwmax;
	double xl, yl;       // Feldgrenzen
	double xh, yh;       // halbe Länge x und y
	double ww_max_dist;  // max. WW-Länge (min(xh, yh)

	double dx, dy;       // Auslenkung zur E-Berechnung

	int nummax;          // max. Anzahl der FLL (progn. maximales num zum Reservieren)
	int nummatch;        // Genau passende Zahl von FL-n

	double a0, h0;       // FL-Abstand und Höhe des Dreiecks
	double xfrac, yfrac; // Bruchteil von xhl und yl, der mit FL-n gefüllt wird !? Hier notwendig??
	// FMSrem:   double clos_dx;      // kürzester FL-Abstand im Feld

	// FMSrem:   //Zusätzliche Variablen für n - n+1 - n - Geometrie:
	// FMSrem:
	// FMSrem:   int nh_num;		//Spaltenzahl des n-Gitters [immer UNGERADE => RB!!]
	// FMSrem:   int nv_num;		//Zeilenzahl des n-Gitters [immer GERADE => RB!!]
	// FMSrem:   int n_1_h_num;	//Spaltenzahl des (n+1)-Gitters.
	// FMSrem:   int n_1_v_num;	//Zeilenzahl des (n+1)-Gitters [auch immer GERADE]
	// FMSrem:   double n_1_v_dist; 	//Vertikaler Abstand im (n+1)-Gitter
	// FMSrem:   double n_1_h_dist;	//Horizontaler Abstand im (n+1)-Gitter (immer gleich n_1_v_dist)
	// FMSrem:   double n_1_h_width;	//Breite des (n+1) Gitters
	// FMSrem:   double gsmult;
	// FMSrem:   double gapsize;
	// FMSrem:   int first_half;	//Erste Hälfte
	// FMSrem:   int point_cnt;	//FL-Counter


	// bool apply_force;     //Kraft ein-/ausschalten // gelöscht, wird über #ifdef gemacht
	// FMSrem:
	// FMSrem:   double ext_force2;    //Externer Kraftparameter

  const gsl_rng_type *Trndm;
  gsl_rng *rndm;

	//Für relax und WechselWirkung funktionen:
	//double fac;         // Faktor für Weg aus Kraft ... wird in relax_step adaptiert
	int faccnt;         // TH: Zählt die Zahl der Relaxationsschritte
	bool second_try;    // Sicherheirsschritt für Ende Relaxation ... Achtung, Langsam....
	bool recalc;        // Für gezwungene Rekalkulation

	// double fprefac;     // Krafterhöhungsfaktor für multiplikation im relax
	// double E_sum;       // [wird nicht verwendet] Enegrie Summe nach dxy_shift. Eingeführt um nicht mehrmals in einem relax_step zu rechnen

	bool debugfail;     // Für debugging von Kräften und Auslenkungen (relax_step(), relax())
	bool allout;		//Zwischenschrittausgabe
	const double resE = 2e-7;  // sinnvolle Energieauflösung -> in relax_step() abgeschätzt aus Tests...


  Eigen::VectorXd prop;     // FL Eigenschaft: pinned=-1, wird bewegt=0, sonst=1
                                 // -2: eingelesen

  int icurr=0; // zeigt auf das jeweils aktuelle Feld von fxy_base[icurr], das andere ist das alte
  int iold=1;  // zeigt auf das jeweils alte Feld von fxy_base[iold]

  Eigen::VectorXd E_base[2];    // Energie(max_num)
#define E E_base[icurr]
#define oE E_base[iold]

	Eigen::VectorXd xy_base[2];    // Vektor für FL-Koordinaten (max_xy): x0, y0, x1, y1, ..q
                                         // zweifach, für alte und neue Werte
#define xy xy_base[icurr]
#define oxy xy_base[iold]

	Eigen::VectorXd fxy_base[2];   // Vektor für WW-Kraft + externe Kraft: x0, y0, x1, y1, ..
                                         // zweifach, für alte und neue Werte
#define fxy fxy_base[icurr]
#define ofxy fxy_base[iold]

  Eigen::VectorXd fext_xy;   // Vektor für externe Kraft: x0, y0, x1, y1, ..

  Eigen::VectorXd fval_base[2];  //  fval ist ein Vektor(max_num) für die Absolutwerte der Kraft fxy()
#define fval fval_base[icurr]
#define ofval fval_base[iold]

  	Eigen::VectorXd ffac;        // Krafterhöhungsfaktor für multiplikation im relax_step
  double ffac0; // Startverschiebung je Kraft: ffac0 (gesetzt im Konstruktor!)
  double minffac, maxffac; // begrenzt die Dynamik der ffac-s

	// Lookup-tables für die Besselfunktionen K0 und K1
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>  bessk0tbl;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> bessk1xtbl;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> bessk1ytbl;
	double ddx, ddxh, dd1x; // Schrittweite für 2-dimensionales mesh für die bessk?-Tabellen, und die Hälfte, und 1/ddx
	virtual void bess_generate(double bessmin) {bessmin++; return;};

	void resize(int n);  // Ändere Zahl der Flusslinien (alle arrays. Achtung bei zusätzlichen blitz-Feldern! Einfügen!)

	int check_bounds(void);  // ist die FL innerhalb der Simulationsbox?

	virtual double calcWW(int i, int pnr, bool force_calc);
  // berechnet die Energie durch bessk0 und Kraft durch bessk1 mittels 2D-Tabelle


	double dsum0 = 0, dsum48 = 0;

  Eigen::VectorXd cosa;   // lokaler Faktor: Kraft zu Verschiebung
	//externe Kraft
	double ext_force;  // FMS Test
	// double force_delta;
	// double force_param;

  double upfac, dnfac; // für FMS Test in force_x, gesetzt in calc_set_forces()
  double maxdx; // max movement Wert, gesetzt im analyze_setup()
  double ftwidth=0.01;

	virtual inline double force_x(const double &x, const double &y, const double &fx, const double &fy, const int &flprop){
    // double fxx=0;
    // if(x<-ftwidth)
    //   if (x>-xh+ftwidth) fxx=ext_force*dnfac;
    //   else fxx=ext_force*dnfac*(xh+x)/ftwidth;
    // else if(x>ftwidth)
    //   if (x<xh-ftwidth) fxx=-ext_force*upfac;
    //   else fxx=ext_force*upfac*(x-xh)/ftwidth;
    // else fxx = -ext_force/ftwidth*x* (x<0?dnfac:upfac);

    // return fxx;}	// Kraft an Position (x, y)
    return ext_force;}

  
	virtual inline double force_y(const double &x, const double &y, const double &fx, const double &fy, const int &flprop)
	{
		// if(apply_force == false) return 0.;

		// if((x < grenze_links) || (x > grenze_rechts))
		//   {
		// 	return  0; // 2*ext_force;
		//   }
		// else
		//   {
		return  0;
		// }
	}

  // Folgendes gehört eigentlich zu analyze()
  double innen_force=0; // Referenzwert für Abbruch (innen absolut) bei großen inneren Kräften (calc_set_forces() und vfield_dyn() in asim)
  double innen_force_x, innen_force_y; // mittlere Innenkraft für Kompensationsrechnung
  double delta_force_x, delta_force_y; // mittlere Differenzkraft zwischen innnen und Aussen
	double innen_force_med=0; // Referenzwert für Abbruch (innen l des MW) bei großen inneren Kräften (calc_set_forces() und vfield_dyn() in asim)
	virtual inline void calc_set_forces() {
    // FMS Test: Ausbalancieren der Gesamtkräfte
    // double up=0, dn=0;
    // // std::cerr << " [";
    // for (int i=0; i< 2*num; i+=2){
    //   if (xy(i)<-ftwidth)
    //     if (xy(i)>-xh+ftwidth) dn++;
    //     else {dn+=(xy(i)+xh)/ftwidth;}// std::cerr<<'.';}
    //   else if (xy(i)>ftwidth)
    //     if (xy(i)<xh-ftwidth) up++;
    //     else {up += (xh-xy(i))/ftwidth; } // std::cerr<<':';}
    //   else if (xy(i)<0) {dn-=xy(i)/ftwidth;} // std::cerr<<'-';}
    //   else {up+=xy(i)/ftwidth; } // std::cerr<<'+';}
    // }
    // dnfac=num*0.5;
    // upfac=dnfac/up;
    // dnfac/=dn;
    // std::cerr << " | " << dn<< ", " << up<< "] ";
    return;
  } // stub für Plugin für bei jedem Schritt in relax() inkludierten Kraftberechnung
	  // (->David Bader)
  // Ende des Blocks für analyze;

	// double fscale; // damit bei kleinen force_delta nicht immer so viel gerechnet werden muss
	// wird beim ersten Mal gemerkt, und geht dann immer von dort weg....
	// reset in set_force()

	unsigned int retcount; // zählt wie oft nicht mit der externen Kraft global verschoben wurde
	//  wird in relax_step() benötigt, um unnötige Verschiebungen zu vermeiden
	unsigned int waitretcount; // wie oft die externe Feldberechnung ausgelassen wird
	// (aus retcount berechnet....)
	// Mechanismus zum Aussetzen der Verschiebung(ext Kraft) bei Sinnlosigkeit:
	// Bei einmaligem "retour" wird eine Runde ausgesetzt, falls danach ok ist wieder ok.
	// Falls danach wieder "retour", wird 2mal ausgesetzt
	// mit jedem signfikanten Verschieben wird retcount wieder gesenkt, mit jedem
	// insignifikanten erhöht, und dann retcount lange gewartet bis wieder mit
	// der externen Kraft verschoben wird


	// Energie berechnen***********************************************************
	struct eVector { float Exy; float x; float y; };

	double f_analyze(int pnr, bool force_calc, bool recalc);
	double point_analyze(double x, double y); //Energie eines Punktes mit Testflusslinie
    

  Eigen::VectorXd dxy_base[2]; // der dr zum berechnen der Arbeit gegenüber der äusseren Kraft
#define dxy dxy_base[icurr]
#define odxy dxy_base[iold]
  


	void resetww(void);
	string base_out; // signifikanter Anteil des Ausgabe Dateinames

	double dxy_shift(std::ofstream &testout);



	// Pinning:
	void swapfl(const int& fl1, const int& fl2); //Vertauscht elemente fl1,fl2 derobigen Arrays
	// TH: Achtung: ww wird von swapfl Funktion NICHT berücksichtigt!!
	int pinnedtostart();
	// FMSrem vector<int> pin2;

	ofstream energy_file;
	void init_energy_file();
	void save_energy_graph(double stepnr, double delta_e, double delta_e_ww, double delta_e_force); // David, zum Testen
	double active_step_nr;
	bool plot_energy = true;

	double max_force;

	void f_analyze_setup(double force, bool force_calc);

public:
	vortex_field_direct(int number // = 200
		, double xlen = 10, double ylen = 10
		, double a0 = 0.3071, double h0 = 2660
		, double xfrac = 1.0, double yfrac = 1.0); /* def. Konstruktor */
	// TH: h0=0: Füllt Feld mit zufälligen Feldlinien
	// TH: number<0: Füllt Feld mit zufälligen Feldlinien mit fill() Methode.
	// FMSrem:   vortex_field(dismode_t disloc_mode, int number=200
	// FMSrem: 	       , double xlen= 10, double ylen= 10
	// FMSrem: 	       , double a0=0.3071, double h0=2660
	// FMSrem: 	       , double xfrac= 1.0, double yfrac= 1.0); /* Konstruktor mit Versetzung */

	vortex_field_direct(); //erstellt dummy
	stringstream outfil1;

	// FMSrem:     vortex_field(int nh_num, int nv_num, int n_1_h_num, double h0, double a0, double gsmult);
	// FMSrem:       //Konstruktor für n - (n+1) - n - Geometrie.


	~vortex_field_direct(void);


	double sum0, sum_last, sumfabs, grenze_links, grenze_rechts; // ?????
	int n_count;

	double min_dis;  // ergibt minimale signifikante Verschiebung: MIN_DIS*xl (in Konstr.)

	// Manipulation****************************************************************
	// FMSrem:     int compress (double clen);   // selbe Fläche, x verkürzt
	// FMSrem:     int verzerr (double measure); // ändert Länge- zu Breite um (1+measure)


	// I/O Functions***************************************************************
	int xy_save(string filnam = "default_xy.dat", int prec = 5);
	bool read_xy(string filnam);
	bool read_new_xy(string filnam);

    bool pickup(string filnam, int n0=0);
// reads positions from old data-file, x=col(1), y=col(2) discards the rest
// n0 is the step number of the starting step

	int ww_save(string filnam, int prec = 11);
	int fxysave(string filnam, int prec = 11);
	int save(string filnam = "default.bin");
	int read(string filnam = "default.bin");

	// void toggle_force(bool);
	// void set_force(double nf);
	void toggle_allout();

	double getE(void);     // Gesamtenergie des Feldes
	int getnum(void);      // Zahl der FL im Feld
	double getXsize(void); // Länge xl des Sim-bereiches
	double getYsize(void); // Höhe yl des Sim-bereiches
	// void show_connected(int id, const char *path = "/proc/self/fd/1");

	// Flusslinien hinzufügen******************************************************
	int v_add(double xadd, double yadd, int pinned);
	// Neue FLL dazugeben; gibt (-1) zurück wenn Feld NICHT vergrößert wurde, sonst (0)
	//  für xadd=yadd=0 zufälliger Platz im ganzen Feld

	int fill(int addn, double xfrac, double yfrac); // füllt von num bis num+addn oder nummatch (bei addn==0) auf
	// TH: (bzw. fügt addn FL hinzu)

	int pinn_rand(const int pinn_num); //TH: Macht zufällige pinn_num FL pinned
	int pullout(const int mode);
	int verschiebe(const double x1, const double y1, const double x2, const double y2);
	void rm_pin();

	int getNpinned();  // return number of pinned vortices
    void setNpinned(int p); // pin the first p vortices
    
  
    void ad_analyze( double rxu=(-1.), double rxl=(-1.), double ryu=(-1.), double ryl=(-1.), string filename="", double extForce=(0.) ); 
    //Andreas Mach-Alles Funktion deklariert; r(x/y)(upper/lower) ist Punkte-Bereich für den E berechnet wird
    
  
	int relax(double force, int currstep); // relaxiert das Feld ins Gleichgewicht
	int relax_step(double force); // definiert einen Relaxationsschritt für das Feld
	// TH: (if return(-1): Keine Bewegung mehr)
	int stepn; // Zähler der Aufrufe von relax_step() je relax;

	void set_base_out(string bo) { base_out = bo; return; }
	string Base_out() { return base_out; }

	void out_debug(int subn);  // gibt das ganze Feld aus....

	// void fprefac_reset() { fprefac = 1; return; }; // reset Vorfaktor für ext. Feld-> Weg auf dx

	int unpin() {
		int i = 0, unp = 0;
		for (; i < num; i++) {
			if (prop(i) < 0) {
				prop(i) = 0;
				unp++;
			}
		}
		return unp;
	}

	// Debugging output
	std::ofstream out_any;  // Damit Routinen das kennen...

  void set_ext_force(double force) {ext_force=force; return;}  // FMS Test

  double fpinx, fpiny; // mittlere Kraft aller gepinnten /num (!)
  
  int breaknumber = 0; //um bei Bedarf relax zum weitermachen zu zwingen (mit eingelesenem txt-file)

};



class vortex_field : public vortex_field_direct  // Flusslinien-Feld in periodischen Randbedingungen,
{
    protected:
        virtual void bess_generate(double  = 1e-14);  // 1e-16...

        virtual double calcWW(int i, int pnr, bool force_calc); // berechnet die Energie und Kraft aus Tabelle (-> bess_generate)

    public:
        vortex_field() : vortex_field_direct() {}

        vortex_field(int number // = 200
                , double xlen = 10, double ylen = 10
                , double a0 = 0.3071, double h0 = 2660
                , double xfrac = 1.0, double yfrac = 1.0) : 
                vortex_field_direct(number, xlen, ylen, a0, h0, xfrac, yfrac)
        {
            double bessmin = 1e-16*gsl_sf_bessel_K0(a0); // wegen Auflösung geht nicht weniger...

            // generiere Lookup-Tables für die gegebene Geometrie:
            bess_generate(bessmin);  // only stub for derived class vortex_field_direct, works in class vortex_field
        }
};

#endif
