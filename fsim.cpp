
// Andreas Burger: ad_analyze() für relax steps aufgerufen

// #define _GNU_SOURCE

#include "vfield.hpp"

ofstream outfil_E;
ofstream res("R_result.dat",ios::app); //global progress output
 


int main(int argc, char*argv[])
{
  //IO Streams&Names********************************************************
  ostringstream outfilename;
  //string filnam; /* zur Ausgabe */ ?!
  //ifstream infile; /* zum Einlesen gespeicherter Info */
    
  int min_lambda=9; /* minimale Länge in Vielfachen von Lambda */

  /* reale Größen */
  double B = 1.0;// 1.0; /* Feld... in T */
  double lambda = 1.4e-7; /* in m */
  
  double r_Jc0 = 2.0e11;  /* maximales Jc in A/m2 */
                        /* fp x Flussliniendicht = Pv = Jc*B */
  //  double force =r_Jc0*Fi0; /* elem. Pinning-Kraft je Längeneinheit N/m */
                        /* r_force = force, Längen kürzen sich weg. */

  double r_a0 = sqrt(Fi0/B/r_h00); /* a0 */
  double r_h0 = r_a0*r_h00; /* Höhenfaktor r_h00 defined in vfield.h*/ 

  /* reduzierte Längen... in Einheiten von lambda */
  /* force ist jetzt die Kraft je Länge lambda einer FL */
  double Kforce = 2.0*M_PI*lambda;
  Kforce = Kforce*Kforce*r_Jc0*2.0e-7*lambda/Fi0; /* force auf K0 normiert */
  double len = 3*min_lambda;
  double hght = 2*min_lambda;
  double a0 = r_a0/lambda;
  double h0 = r_h0/lambda;
  int n_len = 2*(int)((len/h0+1)/2);
  int n_hght = (int)(hght/a0)+1;
  len = n_len*h0;
  hght = n_hght*a0;
  // double Jc0 = r_Jc0*lambda*lambda; /* in A/lambda^2 */
  
  double sepfac=1.8;
  if (argc==2) {
    istringstream insepf(argv[1]);
    insepf >> sepfac;
  }

  double sep=h0*sepfac;

  if (3*sep>hght) hght = (int)((3*sep+a0)/a0)*a0;
  
  int xcnt=26; // muss gerade sein 40x34
  int ycnt=20; // 26x20
  int flc = xcnt*ycnt;
  
  omp_set_num_threads(3);
  
  //vortex_field vf(flc,70,72,h0,a0,1.0); // Versetzg
  vortex_field vf(flc, xcnt*h0, ycnt*a0, a0, h0, 1., 1.); // clean, dep on xcnt, ycnt, flc
  //vortex_field vf(20*(18+20), 38*h0, 20*a0, a0, 0.0, 1., 1.); // pinned small random
  //vortex_field vf(flc, xcnt*h0, ycnt*a0, a0, 0.0, 1., 1.); // random, dep on xcnt, ycnt, flc
  //vortex_field vf(72*(70+72), 142*h0, 72*a0, a0, h0, 1., 1.); // clean
  //vortex_field vf(4680, 27.2185, 18.1591, 0, h0, 1., 1.); // clean David
  

  double external_force = 0; 
  vf.set_ext_force( external_force );  
  cout << "argc: " << argc << endl << endl;
  // if (argc<2) vf.relax(1e-6,0); // 
 
  outfilename.str("");
  outfilename << "clnstart_" << flc;
  // std::cerr << outfilename.str() << std::endl;
  vf.set_base_out(outfilename.str());
  outfilename << "_re-3.dat";  
  vf.xy_save(outfilename.str(),11);
  // Bis hierher ok, auch parallel    
  outfilename.str("");
 
  //vf.read("relaxed_force_exp0_72_70_72_1000mT_binary_00000.bin");
  //vf.read("relaxed_klein.bin");
  outfilename << "force_start_" << flc;
  vf.set_base_out(outfilename.str());
  outfilename << ".dat";

  // clean
  //vf.grenze_links = -10.58; //-11.65;-11.04 //???
  //vf.grenze_rechts = 10.58; // 10.73; //11.65; //???
  // Kforce = 1.e-8; //0.16107; //0.0095206;
  
  ofstream out;
  out.open("force_list_1000mT.txt", ios::app); 
  
  // vf.pinn_rand(6); // 6 zufällige pins...
  

  
  // vf.set_force(external_force); 
  // vf.toggle_force(false);
  // vf.fprefac_reset();  // reset force prefac to 1
  
  // vf.verschiebe(0.151, -0.262, 0.151, -0.1);
  // vf.xy_save("verschoben.dat",11);

  cout << "argc: " << argc << endl;
  if (argc>1) {
    cout << "Picking up: " << argv[1] << endl;
    if (argc==3){
      if (!vf.pickup(argv[1], atoi(argv[2]))) return 66;
      cout << "Starting at " << atoi(argv[2]) << endl << endl;
    } else {
      if (!vf.pickup(argv[1])) return 67;
      cout << "Starting at " << 0 << endl << endl;
    }
    vf.relax(Kforce,0);
    vf.pinn_rand(36); // 36 zufällige pins...    
  } else {
    cout << "No pickup, continue new..."<< endl;
    cout <<  "     ext-force: " <<external_force << endl<<endl;
    vf.set_ext_force(0.0);
    vf.relax(1e-6,0);
    vf.pinn_rand(36); // 36 zufällige pins...
    vf.ww_save("../largerelax.dat",11);
  }
  
  external_force = 0.1; // 5.12e-3; // 1.e-5; // Andreas 0.1 oder 0.01
  vf.set_ext_force(external_force);  
  Kforce = 2*external_force; // Andreas 2* oder 1*

  
  for(int i = 1; i < 12; Kforce*=0.5, i++)  {  //  Kforce Änderung!!!!! // Andreas *=0.5 oder *=0.9
    // external_force*=2;
    // vf.set_ext_force(external_force);
    
    cout<<endl<< "========";
    cout<<endl<< "STEP: " <<i<< "     ext-force: " <<external_force;
    cout<<endl<< "--------" <<endl;
    
    out << i << " " << external_force << endl;
    outfilename.str("");
    outfilename << "cln_" << flc << "_" << setfill('0') << setw(4) << i ;
    vf.set_base_out(outfilename.str());
    outfilename << ".dat";
    
    //vf.fprefac_reset();  // reset force prefac to 1
    vf.relax(Kforce,i);
    
    
    //Andreas: ad_analyze aufrufen.
    //vf.ad_analyze(-1., -1., -1., -1., "step_"+to_string(i), Kforce); //Default=(-1) will calc all points. 
    
    
    vf.xy_save(outfilename.str(),11);
    outfilename.str("");
    outfilename << vf.Base_out()<< ".bin";
    vf.save(outfilename.str());
    outfilename.str("");
    
  }
  
  
  outfilename << "new_" << flc << "last.dat";
  vf.xy_save(outfilename.str(),11);
  outfilename.str("");
  outfilename << "new_" << flc << "last.bin";
  vf.save(outfilename.str());
  outfilename.str("");
  out.close();
  
  return 0;
}
