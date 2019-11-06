#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;

//
// particle data structure
//
typedef struct particle_t
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
  particle_t * next;
}particle_t;

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
void init_particles( int n, particle_t *p );
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void move( particle_t &p );


// 
// bin routines
//
int init_grid();
int grid_index(particle_t *p);
int grid_index(int x, int y);
int particles_per_bin(int bin);
void bin_particle(particle_t *p);
void unbin_particle(particle_t *p);
void reset_bin(int bin);
void bin_forces(int bin, double *dmin, double *davg, int *navg);

//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
