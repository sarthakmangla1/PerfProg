/*
 * This file defines static arrays that contains the primary coordinates
 * of the particles,
 *
 *  Nbody	  Number of particles
 *  Npair	  Number of particle pairs
 *  pos		  Position of the particles
 *  r             Distance of partice from central mass 
 *  vel		  Velocity of the particles
 *  f		  Forces acting on each particle
 *  vis           Viscosity coefficient for each particle
 *  mass	  Mass of each particle
 */

#ifdef DECL
#define DEF
#else
#define DEF extern
#endif
#define Nbody 4*1024
#define  Npair ((Nbody*(Nbody-1))/2)

enum{ Xcoord=0, Ycoord, Zcoord, Ndim };
      
DEF double *pos[Ndim], *velo[Ndim];
DEF double *f[Ndim], *vis, *mass, *radius, *inv_mass, *Gm;
DEF double *r;
DEF double wind[Ndim];
DEF int collisions;

#define G 2.0
#define M_central 1000.0

void evolve(int Nstep, double dt);
void vis_forces(int N, double *f, double *vis, double *velo);
void wind_forces(int N, double *f, double *vis, double velo);
void compute_forces(double dt);
