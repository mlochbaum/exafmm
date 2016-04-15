#ifndef modify_tree_h
#define modify_tree_h
#include "morton_ops.h"

namespace exafmm {
  namespace modify {
    // Initialize each body's KEY value with its Morton key
    void init_ibody(Bounds bounds, Bodies bodies) {
      B_iter Bi0 = bodies.begin();
      real_t diam = max(bounds.Xmax - bounds.Xmin);
      vec3 Xmin = ((bounds.Xmin + bounds.Xmax) - diam) / 2;
      for (int ibody = 0; ibody < bodies.size(); ibody++) {
        B_iter Bi = Bi0 + ibody;
        vec3 X = Bi->X;
        int iX[3] = {0, 0, 0};
        for (int d=0; d<3; d++) iX[d] = (int)((1<<21) * ((X[d] - Xmin[d]) / diam));
        Bi->KEY = morton::interleave(iX[0], iX[1], iX[2]);
      }
    }

    // Assume KEY has been initialized for each body.
    // Perform a radix sort on the bodies.
    void sort_bodies(Bodies bodies) {
      // Radix sort
      int numBodies = bodies.size();
      int *perm  = new int[numBodies]; // Final permutation
      int *p     = new int[numBodies]; // Permutation buffer
      int *perm1 = new int[numBodies]; // Alternate permutation
      for (int i=0; i<numBodies; i++) perm[i] = i;
      // Radix is 2^8 = 256, requiring 64/8 = 8 sorts.
      int nbody[256];
      for (int n=0; n<8; n++) {
// Get the bin of body i
#define BIN(i) (bodies[perm[i]].KEY >> (n*8)) & 255
        // After both loops, nbody[s] is the number of bodies in bin s,
        // and p[i] is the position of body i in its bin.
        for (int i=0; i<256; i++) { nbody[i]=0; }
        for (int i=0; i<numBodies; i++) {
          int s = BIN(i); p[i] = nbody[s]; nbody[s]++;
        }
        // Turn nbody into a list of cumulative sums
        for (int i=0,s=0; i<256; i++) { int t=nbody[i]; nbody[i]=s; s+=t; }
        // Permute perm to give perm1
        for (int i=0; i<numBodies; i++) {
          perm1[p[i] + nbody[BIN(i)]] = perm[i];
        }
        // Swap perm and perm1
        { int *t = perm; perm = perm1; perm1 = t; }
#undef BIN
      }
      delete[] p; delete[] perm1;
      Body *buffer = new Body[numBodies];
      for (int i=0; i<numBodies; i++) { buffer[i] = bodies[perm[i]]; }
      for (int i=0; i<numBodies; i++) { bodies[i] = buffer[i]; }
      delete[] perm;
      delete[] buffer;
    }

    // Subdivide a cell
    // Assumes bodies are all sorted by Morton key
    void subdivide(Bodies bodies, Cells cells, int icell) {
      B_iter B = bodies.begin();
      // Allocate enough space for 8 cells
      // We will shrink to the actual number later
      cells.resize(cells.size() + 8);
      C_iter C = cells.begin();
      C_iter Ci = C + icell;
      if (Ci->NCHILD) return; // Ci is already subdivided

      int jcell = icell; C_iter Cj; // Child cell
      int iX; // Index of child relative to parent
      real_t Rc = Ci->R/2; // Child cell radius
      int shift = 3 * (20 - morton::getLevel(Ci->ICELL));
      // Loop over bodies to place each in a cell
      for (int ibody = 0; ibody < Ci->NBODY; ibody++) {
        Body Bi = bodies[Ci->IBODY + ibody];
        vec3 Xb = Bi.X;
        // If ibody is zero or Bi is not in Cj create a new cell
        int jX = (Bi.KEY >> shift) & 7;
        if (ibody==0 || iX!=jX) {
          jcell++; Cj = C + jcell;
          iX = jX;
          Cj->ICELL = morton::append(Ci->ICELL, iX);
          Cj->IPARENT = icell;
          Cj->NCHILD = 0;
          Cj->IBODY = Ci->IBODY + ibody;
          Cj->NBODY = 0;
          Cj->BODY = B + Cj->IBODY;
          for (int i=0; i<3; i++) {
            Cj->X[i] = iX&(1<<i) ? Ci->X[i]+Rc : Ci->X[i]-Rc;
          }
          Cj->R = Rc;
          if (ibody==0) { Ci->ICHILD = jcell; }
          Ci->NCHILD++;
        }
        Cj->NBODY++;
      }
      cells.resize(icell+1+Ci->NCHILD);
    }
  }
}
#endif
