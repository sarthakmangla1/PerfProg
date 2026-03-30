#include <stdio.h>
#include <math.h>
#include "coord.h"

#ifndef DEBUG_PRINT
#define DEBUG_PRINT 1
#endif

void evolve(int count, double dt){
    int step;
    int i, j;

    for(step = 1; step <= count; step++){
#if DEBUG_PRINT
        printf("timestep %d\n", step);
        printf("collisions %d\n", collisions);
#endif

/* set the viscosity term in the force calculation */
        for(j=0;j<Ndim;j++){
          vis_forces(Nbody,f[j],vis,velo[j]);
        }
/* add the wind term in the force calculation */
        for(j=0;j<Ndim;j++){
          wind_forces(Nbody,f[j],vis,wind[j]);
        }

        compute_forces(dt);

        for (i = 0; i < Nbody; i++) {
            for (j = 0; j < Ndim; j++) {
                pos[j][i] = pos[j][i] + dt * velo[j][i];
            }
        }

        for (i = 0; i < Nbody; i++) {
            for (j = 0; j < Ndim; j++) {
                velo[j][i] = velo[j][i] + dt * (f[j][i] * inv_mass[i]);
            }
        }
    }
}
