#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

//
//  benchmarking program
//
//
//
//



void HalfKick(int n, particle_t *p) {
/*------------------------------------------------------------------------------
    Accelerates atomic velocities, rv, by half the time step.
------------------------------------------------------------------------------*/
    double DeltaT = 0.0005;
    double DeltaTH = .5*DeltaT;
    for (int i=0; i<n; i++){
        p[i].vx = p[i].vx + DeltaTH*p[i].ax;
        p[i].vy = p[i].vy + DeltaTH*p[i].ay;
                            }
                }

double SignR(double v,double x) {if (x > 0) return v; else return -v;}
void ApplyBoundaryCond(int n, particle_t *p) {
/*------------------------------------------------------------------------------
    Applies periodic boundary conditions to atomic coordinates.
------------------------------------------------------------------------------*/
    double Density = 0.0005;
    int InitUcell[2];  //number of unit cells
    InitUcell[0] = 10;     
    InitUcell[1] = 10;
    double Region[2];   //MD box lengths
    double RegionH[2];  //Half the box lengths
////// Compute basic parameters
    for (int k=0;k<2; k++){
        Region[k] = InitUcell[k]/pow(Density, 1.0/2.0);
        RegionH[k] = 0.5*Region[k];
              } 
    for (int i=0; i<n; i++){
        p[i].x = p[i].x - SignR(RegionH[0], p[i].x) - SignR(RegionH[0], p[i].x-Region[0]);
        p[i].y = p[i].y - SignR(RegionH[1], p[i].y) - SignR(RegionH[1], p[i].y-Region[1]);
                         } 
                        }
/*----------------------------------------------------------------------------*/
void ApplyForce(int n, particle_t *p) {
/*------------------------------------------------------------------------------
    Acceleration, are computed as a function of atomic coordinates.
------------------------------------------------------------------------------*/
    int i,j,a,lcy,lcxy,mc[2],c,mc1[2],c1;
    double dr[2],rr,ri2,ri6,r1,rrCut,fcVal,f,rshift[2];

    #define RCUT 0.01   //cutoff length
    #define NCLMAX n/10    //maximum number of linked-list cells
    #define EMPTY -1
    #define DeltaT 0.0005
    #define Density 0.0005



    /////////Variables
    int head[NCLMAX];   //Headers for the linked cell lists
    int lscl[n];     //linked cell lists
    int lc[2];      //Number of cells in the x|y direction
    double rc[2];   //Length of a cell in the x|y direction

    int navg = 0;
    double davg = 0.0;
    double dmin = 1.0;
    int InitUcell[2];     //number of unit cells
    InitUcell[0] = 10;     
    InitUcell[1] = 10;
    double Region[2];   //MD box lengths
    double RegionH[2];  //Half the box lengths
////// Compute basic parameters
    for (int k=0;k<2; k++){
        Region[k] = InitUcell[k]/pow(Density, 1.0/2.0);
        RegionH[k] = 0.5*Region[k];
              } 
    ///// Compute the # of cells for linked cell lists
    for (int k=0; k<2; k++){
    lc[k] = Region[k]/RCUT;
    rc[k] = Region[k]/lc[k];
               }

    /* Reset the potential & forces */
    for (i=0; i<n; i++) p[i].ax = p[i].ay = 0;


  /* Make a linked-cell list, lscl--------------------------------------------*/

    lcy = lc[1];
    lcxy = lc[0]*lcy;

    /* Reset the headers, head */
    for (c=0; c<lcxy; c++) head[c] = EMPTY;

    /* Scan atoms to construct headers, head, & linked lists, lscl */

    for (i=0; i<n; i++) {
        mc[0] = p[i].x/rc[0];
        mc[1] = p[i].y/rc[1];

        /* Translate the vector cell index, mc, to a scalar cell index */
        c = mc[0]*lcy+mc[1];

        /* Link to the previous occupant (or EMPTY if you're the 1st) */
        lscl[i] = head[c];

        /* The last one goes to the header */
        head[c] = i;
    } /* Endfor atom i */

  /* Calculate pair interaction-----------------------------------------------*/

    rrCut = RCUT*RCUT;

    /* Scan inner cells */
    for (mc[0]=0; mc[0]<lc[0]; (mc[0])++)
    for (mc[1]=0; mc[1]<lc[1]; (mc[1])++){

        /* Calculate a scalar cell index */
        c = mc[0]*lcy+mc[1];
        /* Skip this cell if empty */
        if (head[c] == EMPTY) continue;

        /* Scan the neighbor cells (including itself) of cell c */
        for (mc1[0]=mc[0]-1; mc1[0]<=mc[0]+1; (mc1[0])++)
        for (mc1[1]=mc[1]-1; mc1[1]<=mc[1]+1; (mc1[1])++){
            /* Periodic boundary condition by shifting coordinates */
            for (a=0; a<2; a++) {
                if (mc1[a] < 0)
                    rshift[a] = -Region[a];
                else if (mc1[a]>=lc[a])
                    rshift[a] = Region[a];
                else
                    rshift[a] = 0.0;
            }
            /* Calculate the scalar cell index of the neighbor cell */
            c1 = ((mc1[0]+lc[0])%lc[0])*lcy
                +((mc1[1]+lc[1])%lc[1]);
            /* Skip this neighbor cell if empty */
            if (head[c1] == EMPTY) continue;

            /* Scan atom i in cell c */
            i = head[c];
            while (i != EMPTY) {

                /* Scan atom j in cell c1 */
                j = head[c1];
                while (j != EMPTY) {

                    /* Avoid double counting of pairs */
                    if (i < j) {
                        /* Pair vector dr = r[i]-r[j] */
                        rr = 0.0;
                        dr[0] = p[i].x - (p[j].x+rshift[0]);
                        dr[1] = p[i].y - (p[j].y+rshift[1]); 
                        rr += dr[0]*dr[0];
                        rr += dr[1]*dr[1];

                        /* Calculate potential & forces if rij<RCUT */
                        if (rr < rrCut) {
                           apply_force( p[i], p[j],&dmin,&davg,&navg);
                        }
                    } /* Endif i<j */

                    j = lscl[j];
                } /* Endwhile j not empty */

                i = lscl[i];
            } /* Endwhile i not empty */

        } /* Endfor neighbor cells, c1 */
    } /* Endfor central cell, c */
}

void SingleStep(int n, particle_t *p) {
/*------------------------------------------------------------------------------
    r & rv are propagated by DeltaT in time using the velocity-Verlet method.
------------------------------------------------------------------------------*/
    HalfKick(n, p); /* First half kick to obtain v(t+Dt/2) */
    for (int i=0; i<n; i++){ /* Update atomic coordinates to r(t+Dt) */
        p[i].x = p[i].x + DeltaT*p[i].vx;
        p[i].y = p[i].y + DeltaT*p[i].vy;
                            }
    ApplyBoundaryCond(n, p);
    ApplyForce(n, p); /* Computes new force, a(t+Dt) */
    HalfKick(n, p); /* Second half kick to obtain v(t+Dt) */
}


int main( int argc, char **argv )
{    
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );
    /////// add constants

    
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );  //initialize particles
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );



    for( int step = 0; step < NSTEPS; step++ )
    { 
        navg = 0;
        davg = 0.0;
        dmin = 1.0;

        SingleStep(n, particles);

        //  move particles
        //
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );		

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;
		
          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -the minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    //
    // Printing summary data
    //
    if( fsum) 
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
