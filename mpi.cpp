#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include "common.h"

using namespace std;

int n_proc, rank;
MPI_Datatype PARTICLE;

void dist_particles(vector<particle_t> & main, vector<particle_t> & boundary, int * ybounds) {
  vector< vector<int> > counts(n_proc, vector<int>(2));
  vector<particle_t> smain[n_proc];
  vector<particle_t> sboundary[n_proc];

  int x, y;

  for (int i = 0; i < main.size(); i++) {
    x = grid_x(&main[i]);
    y = grid_y(&main[i]);

    for (int j = 0; j < n_proc; j++) {
      if (y < ybounds[j]) {
        smain[j].push_back(main[i]);
        counts[j][0]++;
        if (j > 0 && y == ybounds[j-1]) {
          sboundary[j-1].push_back(main[i]);
          counts[j-1][1]++;
        }
        if (j < n_proc-1 && y == ybounds[j]-1) {
          sboundary[j+1].push_back(main[i]);
          counts[j+1][1]++;
        } 
        break;
      }
    }
  }

  vector< vector<int> > rcounts = counts;
  MPI_Request sen[3][n_proc];
  MPI_Request rec[3][n_proc];

  for (int i = 0; i < n_proc; i++) {
    MPI_Isend(&counts[i][0], 2, MPI_INT, i, 0, MPI_COMM_WORLD, &sen[0][i]);
    MPI_Isend(&smain[i][0], counts[i][0], PARTICLE, i, 0, MPI_COMM_WORLD, &sen[1][i]);
    MPI_Isend(&sboundary[i][0], counts[i][1], PARTICLE, i, 0, MPI_COMM_WORLD, &sen[2][i]);
  }

  main.resize(0);
  boundary.resize(0);
  for (int i = 0; i < n_proc; i++) {
    MPI_Recv(&rcounts[i][0], 2, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    int msize = main.size();
    int bsize = boundary.size();
    main.resize(msize + rcounts[i][0]);
    boundary.resize(bsize + rcounts[i][1]);
    MPI_Recv(&main[msize], rcounts[i][0], PARTICLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&boundary[bsize], rcounts[i][1], PARTICLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }


  MPI_Waitall(3*n_proc, sen[0], MPI_STATUSES_IGNORE);
  return;
}









/*
   void dist_particles( vector<particle_t> & main, vector<particle_t & boundary, 
   int counts[n_proc];
   for (int i = 0; i < n_proc; i++) {
   counts[i] = 0;
   }
   vector<particle_t> sparticles[n_proc];


   for (int i = 0; i < n; i++) {
   if (particles[i].valid) {
   int x = grid_x(&particles[i]);
   int y = grid_y(&particles[i]);

   for (int j = 0; j < n_proc; j++) {
   if (y < ybounds[j]) {
   counts[j]++;
   sparticles[j].push_back(particles[i]);
   particles[i].valid = false;
   break;
   }
   }
   }
   }

   for (int i = 0; i < n_proc; i++) {
   printf("%d: ", i);
   for (int j = 0; j < counts[i]; j++) {
   printf("%d ", indices[i][j]);
   }
   printf("\n");
   }

   MPI_Request sents[n_proc];
   MPI_Status status[n_proc];
   for (int i = 0; i < n_proc; i++) {
   MPI_Isend(&counts[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &sents[i]);
   if (counts[i] > 0) {
   MPI_Isend(&sparticles[i][0], counts[i], PARTICLE, i, 0, MPI_COMM_WORLD, &sents[i]);
   }
   }

   vector<particle_t> nlocal;

   for (int i = 0; i < n_proc; i++) {
   int count;
   int size = nlocal.size();
   MPI_Recv(&count, 1, MPI_INT, i, 0, MPI_COMM_WORLD, NULL);
   if (count > 0) {
   nlocal.resize(size + count);
   MPI_Recv(&nlocal[size], count, PARTICLE, i, 0, MPI_COMM_WORLD, NULL);
   }
   }

   return nlocal;
   }
   */

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
  int navg, nabsavg=0;
  double dmin, absmin=1.0,davg,absavg=0.0;
  double rdavg,rdmin;
  int rnavg; 

  //
  //  process command line parameters
  //
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
  char *savename = read_string( argc, argv, "-o", NULL );
  char *sumname = read_string( argc, argv, "-s", NULL );

  //
  //  set up MPI
  //
  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  //
  //  allocate generic resources
  //
  FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
  FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;


  particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
  for (int i = 0; i < n; i++) {
    particles[i].valid = false;
    particles[i].next = NULL;
  }

  MPI_Type_contiguous( sizeof(particle_t), MPI_BYTE, &PARTICLE );
  MPI_Type_commit( &PARTICLE );

  //
  //  initialize and distribute the particles (that's fine to leave it unoptimized)
  //
  set_size( n );
  int num_bins = init_grid();
  int num_cols = sqrt(num_bins);

  int ystarts[n_proc];
  int ybounds[n_proc];
  for (int i = 0; i < n_proc; i++) {
    ystarts[i] = i*ceil((double)num_cols/n_proc);
    ybounds[i] = (i+1)*ceil((double)num_cols/n_proc);
  }
  ybounds[n_proc-1] = num_cols;

  int ystart = ystarts[rank];
  int ybound = ybounds[rank];

  vector<particle_t> local;
  vector<particle_t> boundary;

  if (rank == 0) {
    init_particles(n, particles);
    local.insert(local.end(), &particles[0], &particles[n]);
  }




  //
  //  simulate a number of time steps
  //
  double simulation_time = read_timer( );
  double dist_time = 0;
  double bin_time = 0;
  double force_time = 0;
  double move_time = 0;
  double start;
  for( int step = 0; step < NSTEPS; step++ )
  {
    navg = 0;
    dmin = 1.0;
    davg = 0.0;

    start = read_timer();
    dist_particles(local, boundary, ybounds);
    dist_time += read_timer() - start;
    
    int size = local.size();
    start = read_timer();
    for (int i = 0; i < local.size(); i++) {
      local[i].next = NULL;
      local[i].ax = local[i].ay = 0;
      bin_particle(&local[i]);
    }

    for (int i = 0; i < boundary.size(); i++) {
      boundary[i].next = NULL;
      bin_particle(&boundary[i]);
    }
    bin_time += read_timer() - start;
    
    start = read_timer();
    for (int i = 0; i < num_bins; i++) {
      if (particles_per_bin(i) > 0) {
        bin_forces(i, &dmin, &davg, &navg);
      }
    }
    force_time += read_timer() - start;

    for (int i = 0; i < local.size(); i++) {  
      unbin_particle(&local[i]);
    }

    for (int i = 0; i < boundary.size(); i++) {
      unbin_particle(&boundary[i]);
    }

    for (int i = 0; i < num_bins; i++) {
      reset_bin(i);
    }
    
    start = read_timer();
    for (int i = 0; i < local.size(); i++) {
      move(local[i]);
    }
    move_time += read_timer() - start;

    if (find_option(argc, argv, "-no") == -1) {
      printf("entered\n");
      MPI_Request sen[2];
      MPI_Isend(&size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &sen[0]);
      MPI_Isend(&local[0], size, PARTICLE, 0, 0, MPI_COMM_WORLD, &sen[1]);


      MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);

      if (rank == 0) {
        int index = 0;
        for (int i = 0; i < n_proc; i++) {
          int count;
          MPI_Recv(&count, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Recv(&particles[index], count, PARTICLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          index += count;
        }
        sort(particles, particles+n, pcompare);

        if ( fsave && (step%SAVEFREQ) == 0 ) {
          save( fsave, n, particles );
        }

        if (rnavg) {
          absavg += rdavg/rnavg;
          nabsavg++;
        }
        if (rdmin < absmin) absmin = rdmin;
      }

      MPI_Waitall(2, sen, MPI_STATUSES_IGNORE);
    }

  }
  simulation_time = read_timer( ) - simulation_time;

  if (rank == 0) {  
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    printf("\n%g\n%g\n%g\n%g\n", dist_time, bin_time, force_time, move_time);

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
      fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
  }

  //
  //  release resources
  //
  if ( fsum )
    fclose( fsum );
  free( particles );
  if( fsave )
    fclose( fsave );

  MPI_Finalize( );

  return 0;
}
