#include <math.h>
#include "coord.h"

void vis_forces(int N,double *f, double *vis, double *velo)
{
  int i;
          for(i=0;i<N;i++){
            f[i] = -vis[i] * velo[i];
          }
}
void wind_forces(int N,double *f, double *vis, double velo)
{
  int i;
          for(i=0;i<N;i++){
            f[i] = f[i] -vis[i] * velo;
          }
}
void add_norms(int N,double *r, double *delta)
{
  int k;
        for(k=0;k<N;k++){
          r[k] += (delta[k] * delta[k]);
        }
}

double forces(double Wv, double deltav, double rv){
  return Wv*deltav/(rv * rv * rv);
}

void compute_forces(double dt){
    int i, j;
    double size;
    double dx, dy, dz;
    double r_sq, r_val;
    double force_scalar;

#define BLOCK 128

    /* central force — use precomputed Gm[] */
    for (i = 0; i < Nbody; i++) {
        double px = pos[0][i];
        double py = pos[1][i];
        double pz = pos[2][i];
        r_sq = px * px + py * py + pz * pz;
        r_val = sqrt(r_sq);
        double r3 = r_sq * r_val;
        double gm = Gm[i] * M_central / r3;
        f[0][i] -= gm * px;
        f[1][i] -= gm * py;
        f[2][i] -= gm * pz;
    }

    /* particle-particle forces — cache blocked, local accumulator, local collision counter */
    for (i = 0; i < Nbody; i++) {
        double gmi = Gm[i];
        double ri  = radius[i];
        double fix = 0.0, fiy = 0.0, fiz = 0.0;
        int local_col = 0;

        for (int jb = i + 1; jb < Nbody; jb += BLOCK) {
            int jend = jb + BLOCK < Nbody ? jb + BLOCK : Nbody;
            for (j = jb; j < jend; j++) {
                dx = pos[0][i] - pos[0][j];
                dy = pos[1][i] - pos[1][j];
                dz = pos[2][i] - pos[2][j];

                r_sq = dx * dx + dy * dy + dz * dz;
                r_val = sqrt(r_sq);

                force_scalar = gmi * mass[j] / (r_sq * r_val);
                size = ri + radius[j];

                double sign = (r_val < size) ? 1.0 : -1.0;
                local_col += (r_val < size);

                fix += sign * force_scalar * dx;
                fiy += sign * force_scalar * dy;
                fiz += sign * force_scalar * dz;

                f[0][j] -= sign * force_scalar * dx;
                f[1][j] -= sign * force_scalar * dy;
                f[2][j] -= sign * force_scalar * dz;
            }
        }

        f[0][i] += fix;
        f[1][i] += fiy;
        f[2][i] += fiz;
        collisions += local_col;
    }

}



