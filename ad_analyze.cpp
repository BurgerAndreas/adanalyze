/*  Writen by ANDREAS BURGER in March 2020.
 * 
 * This program primarly finds adjacencies from a set of initial points
 * and calculates energies with that information.
 * It utalizes external functions for the triangulation and the 
 * calculation of energies.
 * 
 * It is seperated into 6 segments: 
 * SETUP, TRANSLATION, TRIANGULATION,
 * ADJACENCY, SHIFTS AND ENERGIES, OUTPUT
*/

#include "delabella.h" //Andreas: must be bevore vfield.hpp, becouse of (define xy).
#include "vfield.hpp"
#include <iostream>
#include <fstream>
#include <math.h>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <string>
#include <chrono>

using namespace std;
using namespace std::chrono;

void vortex_field_direct::ad_analyze( double rxu, double rxl, double ryu, double ryl, string filename, double Kforce ) // Default values and declaration in vfield.hpp.
{
    cout <<"\n\n\t=== ad_analyze started ==="<<endl; // Startup message.
    cout <<"\tpassed " <<rxu <<" " <<rxl <<" " <<ryu <<" " <<ryl <<" " <<Kforce <<endl; // Passed parameters.
    if( rxu!=-1 || rxl!=-1 || ryu!=-1 || ryl !=-1) filename+="_roi"; // If region of interest spezified, mark output files.
    
    // ====================== SETUP ======================
    struct MyPoint // Needed for delabella.
    {
        double x; // Coordinates.
        double y;
        int ind; // Index. For translated Points: equal to index of initial Point it came from.
    };
    
    // Parameters.
    unsigned int nbrange = 8; // Max number of adjacent points. Defines matrices sizes. Better solution?
    const double edge=1.5*a0; // Definition of near-boundary for translation. 
    const double lx = xl, ly = yl; // Boundary containing initial points.
    const int flc = num; // Number of initial points.
    const double norm=a0/100; // Defines length of shift towards adjacent points.
    int ix , iy; // Temp variable used for iteration in various loops.
    string destination = "./Output/"; // Location of output files.  
    
    // Containers.
    Eigen::VectorXd adNum = Eigen::VectorXd::Zero( flc ); // Number of adjacencies for a point.
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> adIndex( flc, nbrange ); // Indices of adjacent points.
    Eigen::MatrixXd shiftE = Eigen::MatrixXd::Zero( flc, nbrange ); // Energies after shift towards adjacent points.           
    Eigen::MatrixXd adVec = Eigen::MatrixXd::Zero( flc*2, nbrange ); // Vectors towards adjacent points. 
    Eigen::VectorXd regE = Eigen::VectorXd::Zero( flc ); // Energies un-shifted, in case not already calculated.
    
    
    // ====================== TRANSLATION ======================  
    auto t1 = high_resolution_clock::now(); // Timekeeping.
    MyPoint* extendCloud = new MyPoint[9*flc]; // Temp Array of Points for triangulation.
    
    unsigned int transPointCount=0; // Temp variables for translation loop.
    unsigned int copyCount=0;    
    
    for (int i=0; i<flc; i++ ) // Will translate near-boundary-points.
    {
        ix= 2*i;
        iy= ix+1;
        
        extendCloud[i].x = xy(ix); // Copying of initial points.
        extendCloud[i].y = xy(iy); 
        extendCloud[i].ind = i; // Index of copied point = index initial point.
        copyCount++;
        
        if ( xy(iy) > ((ly/2)-edge) ) // Near top.
        {
            extendCloud[ flc + transPointCount ].x = ( xy(ix) ); // Translate x.
            extendCloud[ flc + transPointCount ].y = ( xy(iy) - ly ); // Translate y.
            extendCloud[ flc + transPointCount ].ind = ( i ); // Same Index as initial point.
            transPointCount++;
            
            if ( xy(ix) < -((lx/2)-edge) ){ // Top left corner.
                extendCloud[ flc + transPointCount ].x = ( xy(ix) + lx );
                extendCloud[ flc + transPointCount ].y = ( xy(iy) - ly );
                extendCloud[ flc + transPointCount ].ind = ( i );
                transPointCount++;                
            }else if ( xy(ix) > ((lx/2)-edge) ){ // Top right corner.
                extendCloud[ flc + transPointCount ].x = ( xy(ix) - lx );
                extendCloud[ flc + transPointCount ].y = ( xy(iy) - ly );
                extendCloud[ flc + transPointCount ].ind = ( i );
                transPointCount++;                
            }              
        }else if ( xy(iy) < -((ly/2)-edge) ){ // Near bottom.
            extendCloud[ flc + transPointCount ].x = ( xy(ix) );
            extendCloud[ flc + transPointCount ].y = ( xy(iy) + ly );
            extendCloud[ flc + transPointCount ].ind = ( i );
            transPointCount++;            
            if ( xy(ix) < -((lx/2)-edge) ){ // Bottom left corner.
                extendCloud[ flc + transPointCount ].x = ( xy(ix) + lx );
                extendCloud[ flc + transPointCount ].y = ( xy(iy) + ly );
                extendCloud[ flc + transPointCount ].ind = ( i );
                transPointCount++;                
            }else if ( xy(ix) > ((lx/2)-edge) ){ // Bottom right corner.
                extendCloud[ flc + transPointCount ].x = ( xy(ix) - lx );
                extendCloud[ flc + transPointCount ].y = ( xy(iy) + ly );
                extendCloud[ flc + transPointCount ].ind = ( i );
                transPointCount++;                
            }            
        }if ( xy(ix) < -((lx/2)-edge) ){ // Near left.
            extendCloud[ flc + transPointCount ].x = ( xy(ix) + lx );
            extendCloud[ flc + transPointCount ].y = ( xy(iy) );
            extendCloud[ flc + transPointCount ].ind = ( i );
            transPointCount++;            
        }if ( xy(ix) > ((lx/2)-edge) ){ // Near right.
            extendCloud[ flc + transPointCount ].x = ( xy(ix) - lx );
            extendCloud[ flc + transPointCount ].y = ( xy(iy) );
            extendCloud[ flc + transPointCount ].ind = ( i );
            transPointCount++;            
        }
    } // End for() translate near-boundary-points.     
    const unsigned int triaNum = transPointCount+copyCount; // Number of points used for Triangulation.    
    
    size_t newSize = triaNum; // Resize "extendCloud" into new array "cloud". 
    MyPoint* cloud = new MyPoint[newSize];                
    memcpy( cloud, extendCloud, newSize * sizeof(MyPoint) );            
    delete [] extendCloud;    
    
    
    // ====================== TRIANGULATION ======================        
    auto t2 = high_resolution_clock::now(); // Timekeeping.
    IDelaBella* idb = IDelaBella::Create(); // Call for triangulation.
    unsigned int verts = idb->Triangulate(triaNum, &cloud->x, &cloud->y, sizeof(MyPoint));        
    if (verts<=0) {cout <<"Error Triangulation: verts !> 0" <<endl;} // No points given or all points are colinear. 
    
    const DelaBella_Triangle* dela = idb->GetFirstDelaunayTriangle();    
    #define d0i (dela->v[0]->i) // v[] are verteces of a given triangle.
    #define d1i (dela->v[1]->i) // i are the indices of those points.
    #define d2i (dela->v[2]->i)
        
    // ====================== ADJACENCY ====================== 
    auto t3 = high_resolution_clock::now(); // Timekeeping.
    for (unsigned int j = 0; j<(verts/3); j++) // All triangles. verts are number of verteces from delabella.
    {   
        // Explanation: verteces in triangles are ordered clockwise.
        // Every edge is contained in two triangles: once as P1->P2 and once as P2->P1 (because they are ordered clockwise).
        // Every edge needs to be noted twice: once as "P1 is Nb of P2", once as "P2 is Nb of P1".
        // So we note the edges "in a circle".
        
        if( (d0i < flc) && (adNum(d0i) < nbrange) ){ // Only information for initial points noted.
            adIndex( d0i, adNum(d0i) ) = cloud[ d1i ].ind; // Save Indices of adjacent points.    
            
            adVec( 2*d0i, adNum(d0i) ) = (cloud[d1i].x - cloud[d0i].x); // Shift towards adjacent point.    
            adVec( (2*d0i)+1, adNum(d0i) ) = (cloud[d1i].y - cloud[d0i].y); 
            
            adNum( d0i )++; // Count number of adjacencies. 
        }  
        if( (d1i < flc) && (adNum(d1i) < nbrange) ){
            adIndex( d1i, adNum(d1i) ) = cloud[ d2i ].ind;   
            
            adVec( 2*d1i, adNum(d1i) ) = (cloud[d2i].x - cloud[d1i].x);
            adVec( (2*d1i)+1, adNum(d1i) ) = (cloud[d2i].y - cloud[d1i].y);
            
            adNum( d1i )++; 
        }
        if( (d2i < flc) && (adNum(d2i) < nbrange) ){
            adIndex( d2i, adNum(d2i) ) = cloud[ d0i ].ind;
            
            adVec( 2*d2i, adNum(d2i) ) = (cloud[d0i].x - cloud[d2i].x);
            adVec( (2*d2i)+1, adNum(d2i) ) = (cloud[d0i].y - cloud[d2i].y);
            
            adNum( d2i )++; 
        }
        dela = dela->next; // Next Triangle.
    }      
    
               
    // ====================== SHIFTS AND ENERGIES ====================== 
    auto t4 = high_resolution_clock::now(); // Timekeeping.
    double xsav, ysav; // Temp save before shift.
    double shiftlen; // Temp variable for length of shift / length of vector.
    if( rxu==(-1.) ) rxu=  lx/2; // Region of interest for shifted energy calculation.
    if( rxl==(-1.) ) rxl= -lx/2; // If passed as default value, calc whole region.
    if( ryu==(-1.) ) ryu=  ly/2;
    if( ryl==(-1.) ) ryl= -ly/2;
    
    for (int pnr=0; pnr<flc; pnr++) // All initial points.
    {
        ix= 2*pnr;
        iy= ix+1;
        
        for (int j=0; j<pnr; j++) // Calc E un-shifted for further calculation.
        {
            regE(pnr) += calcWW(pnr, j, false); // WW between pnr and j.
        }
        for (int j=pnr+1; j<flc; j++) 
        {
            regE(pnr) +=calcWW(pnr, j, false);
        } 
                
        if( (xy(ix)>rxu) || (xy(ix)<rxl) || (xy(iy)>ryu) || (xy(iy)<ryl) ) continue; // Is pnr outside region of interest?
        
        xsav=xy(ix); // Save coordinates before shift.
        ysav=xy(iy);  
        
        for (int k=0; k<adNum(pnr); k++) // Will calc energies after shift towards adjacent points. 
        {
            shiftlen = norm / hypot(adVec(ix,k),adVec(iy,k));
            xy(ix)= xsav + (adVec(ix,k) * shiftlen); // Execute shift with length norm.
            xy(iy)= ysav + (adVec(iy,k) * shiftlen);              
            
            for (int j=0; j<pnr; j++) // WW between pnr and all other points.
            {
                shiftE(pnr, k) += calcWW(pnr, j, false); // WW between pnr and k.
            }
            for (int j=pnr+1; j<flc; j++) 
            {
                shiftE(pnr, k) += calcWW(pnr, j, false);
            }             
        }
        xy(ix)=xsav; // Return coordinates to regular position.
        xy(iy)=ysav; // Regular means before shift.           
    }
    
    
    // ====================== OUTPUT ======================
    auto t5 = high_resolution_clock::now(); // Timekeeping.
   
    
    
    // Output difference in regularE between adjacent points
    ofstream odregEnb(destination+filename+"_dEad.dat");
    odregEnb <<"# x y xvec_nb yvec_nb regE(nb)-regE(xy)" <<endl;
    for (int j=0; j<flc; j++) // All points.
    {
        for (int k=0; k<adNum(j); k++) // All adjacencies of j.
        {
            ix = j * 2;
            iy = ix + 1;
            odregEnb <<xy(ix) <<" " <<xy(iy) <<" "
            << adVec(ix,k)/2.3 <<" " <<adVec(iy,k)/2.3 <<" " // Length of vector in plot.
            <<( regE(adIndex(j,k)) - regE(j) ) / hypot(adVec(ix,k),adVec(iy,k))<<endl;
            // Note that with the logic used in segment ADJACENCY, we have every adjacency twice.
        }
    } odregEnb <<"# Kforce= " <<Kforce <<endl;
    odregEnb.close();
    
    
    
    // Output (shiftedE-regularE)/(shiftlength) with un-normed direction vector.
    ofstream odiffE(destination+filename+"_dEshift.dat");
    odiffE << "# x y shiftx shifty (regularE-shiftedE)/(shiftlength)" <<endl;
    for(int j=0; j<flc; j++) //Points
    {
        ix = 2 * j;
        iy = ix + 1;
        
        //Is pnr outside region of interest?
        if( (xy(ix)>rxu) || (xy(ix)<rxl) || (xy(iy)>ryu) || (xy(iy)<ryl) ) continue;
        
        for(int k=0; k<adNum(j); k++)
        {           
            odiffE <<xy(ix) <<" " <<xy(iy) 
            <<" " <<adVec(ix,k)/2.3 <<" " <<adVec(iy,k)/2.3 // Length of vector in plot.
            <<" " <<( shiftE(j,k)-regE(j) )/norm <<endl;  
            // Calcs in the Kforce in x direction, over length of shift in x direction.
        }
    } odiffE <<"# rxu/l ryu/l norm Kforce: " <<rxu <<" " <<rxl <<" " <<ryu <<" " <<ryl <<" " <<norm <<" " <<Kforce <<endl; // Boundary of points.
    odiffE.close();
    

    
    // Output triangulation with number of adjacencies.
    dela = idb->GetFirstDelaunayTriangle();    
    ofstream otrian(destination+filename+"_triangulation.dat"); 
    otrian <<"# x y adNum" <<endl;
    for (unsigned int j = 0; j<(verts/3); j++)
    {
        if( d0i < flc ){ // Translated points are set to 0.
            otrian << dela->v[0]->x <<" " << dela->v[0]->y <<" " <<adNum( d0i ) <<endl;
        }else{otrian << dela->v[0]->x <<" " << dela->v[0]->y <<" " <<0<<endl;}
        
        if( d1i < flc ){
            otrian << dela->v[1]->x <<" " << dela->v[1]->y <<" " <<adNum( d1i ) <<endl;       
        }else{otrian << dela->v[1]->x <<" " << dela->v[1]->y <<" " <<0<<endl;}
        
        if( d2i < flc ){
            otrian << dela->v[2]->x <<" " << dela->v[2]->y <<" " <<adNum( d2i ) <<"\n" <<endl;
        }else{otrian << dela->v[2]->x <<" " << dela->v[2]->y <<" " <<0 <<"\n"<<endl;}    
        
        // Points of the same triangle will form block of threes. Enables gnuplot linespoints command.
        dela = dela->next;
    } otrian <<"# " <<lx/2 <<" " <<-lx/2 <<" " <<ly/2 <<" " <<-ly/2 <<endl; // Boundary of points.
    otrian.close(); 
    
    
    // Output pinned points.
    ofstream opinn(destination+filename+"_pinned.dat");
    for(int j = 0; j<flc; j++)
    {
        ix = 2 * j;
        iy = ix + 1;
        if (prop(j) < 0) opinn << xy(ix) <<" " << xy(iy) <<" " << prop(j) <<endl;
    } opinn.close();
    
    
    
    // Finish.
    delete[] cloud; 
    idb->Destroy(); // DelaBella triangulation.
    
    auto t6 = high_resolution_clock::now(); 
    auto d1 = duration_cast<milliseconds>(t2 - t1); // Time the segments took.
    auto d2 = duration_cast<milliseconds>(t3 - t2);
    auto d3 = duration_cast<milliseconds>(t4 - t3);
    auto d4 = duration_cast<milliseconds>(t5 - t4);
    auto d5 = duration_cast<milliseconds>(t6 - t5);
    
    cout<<"\tOutput saved to: " <<destination+filename+"..." <<endl;
    cout<<"\tTime taken for Translation, Triangulation, Adjacency, Energies, Output: " 
    <<"\n\t" <<d1.count() <<"ms " <<d2.count() <<"ms " <<d3.count() <<"ms " <<d4.count() <<"ms " <<d5.count() <<"ms " <<endl;
    cout<<"\t=== ad_analyze done ===\n\n" <<endl;
    return;
} 
