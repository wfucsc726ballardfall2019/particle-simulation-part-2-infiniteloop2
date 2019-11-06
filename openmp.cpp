#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"

//
//  benchmarking program
//
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
    printf( "-p to specify the number of threads for OpenMP\n");
    return 0;
  }

  int n = read_int( argc, argv, "-n", 1000 );

  char *savename = read_string( argc, argv, "-o", NULL );
  char *sumname = read_string( argc, argv, "-s", NULL );

  FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
  FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

  int num_threads = read_int ( argc, argv, "-p", 1 );
  omp_set_num_threads(num_threads);

  particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
  set_size( n );
  init_particles( n, particles );
  int num_bins = init_grid();

  int index;
  omp_lock_t * bins = (omp_lock_t*) malloc(num_bins * sizeof(omp_lock_t));
  for (int i = 0; i < num_bins; i++) {
    omp_init_lock(&bins[i]);
  }
  
  //
  //  simulate a number of time steps
  //
  double simulation_time = read_timer( );

  for( int step = 0; step < NSTEPS; step++ )
  {
    navg = 0;
    davg = 0.0;
    dmin = 1.0;
    
    //
    //  bin particles
    //
#pragma omp parallel for private(index)
    for (int i = 0; i < n; i++) 
    {
      particles[i].ax = particles[i].ay = 0;
      index = grid_index(&particles[i]);
      omp_set_lock(&bins[index]);
      bin_particle(&particles[i]);
      omp_unset_lock(&bins[index]);
    }


    //
    //  compute forces
    //
#pragma omp parallel for reduction(+:davg,navg) reduction(min:dmin) schedule(static)
    for (int i = 0; i < num_bins; i++) {
      bin_forces(i, &dmin, &davg, &navg);
    }

    //
    //  unbin particles
    //
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
      unbin_particle(&particles[i]);
    }

    //
    //  reset bins
    //
#pragma omp parallel for
    for (int i = 0; i < num_bins; i++)
    {
      reset_bin(i);
    }

    //
    //  move particles
    //
#pragma omp parallel for
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

  printf( "n = %d, threads = %d, simulation time = %g seconds", n, num_threads, simulation_time);

  if( find_option( argc, argv, "-no" ) == -1 )
  {
    if (nabsavg) absavg /= nabsavg;
    // 
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
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
    fprintf(fsum,"%d %d %g\n",n,num_threads,simulation_time);

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
