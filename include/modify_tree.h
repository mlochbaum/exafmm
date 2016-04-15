#ifndef modify_tree_h
#define modify_tree_h
#include "morton_ops.h"

namespace exafmm {
  namespace modify {
    // Initialize each body's IBODY value with its Morton key
    // Assume Ci0[0] is the root cell
    void init_ibody(C_iter Ci0, B_iter Bi0) {
      vec3 Xmin = Ci0->X;
      real_t diameter = Ci0->R;
      for (int ibody = 0; ibody < Ci0->NBODY; ibody++) {
        B_iter Bi = Bi0 + ibody;
        vec3 X = Bi->X;
        int iX[3] = {0, 0, 0};
        for (int d=0; d<3; d++) iX[d] = int((1<<20) * (X[d] - Xmin[d]) / diameter);
        Bi->IBODY = morton::interleave(iX[0], iX[1], iX[2]);
      }
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
      int shift = 3 * (20 - morton::getLevel(Cj->ICELL));
      // Loop over bodies to place each in a cell
      for (int ibody = 0; ibody < Ci->NBODY; ibody++) {
        Body Bi = bodies[Ci->IBODY + ibody];
        vec3 Xb = Bi.X;
        // If ibody is zero or Bi is not in Cj create a new cell
        int jX = (Bi.IBODY >> shift) & 7;
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

    /*
    void sortbodies(Body * bodies, int len) {
      int n[8]; for (int i=0; i<8; i++) n[i]=0;
      Body Bt[Ci->NBODY]; int sub[Ci->NBODY]; int ind[Ci->NBODY];
      for (int ibody = 0; ibody < Ci->NBODY; ibody++) {
        Bt[ibody] = bodies[Ci->IBODY + ibody];
        vec3 Xb = Bt[ibody].X;
        int s=0; for(int i=0;i<3;i++) s += (Xb[i]>=Ci->X[i]) << i;
        sub[ibody] = s;
        ind[ibody] = n[s];
        n[s]++;
      }
      int sum=0; for (int i=0; i<8; i++) { int t=n[i]; n[i]=sum; sum+=t; }
      for (int ibody = 0; ibody < Ci->NBODY; ibody++) {
        bodies[Ci->IBODY + n[sub[ibody]] + ind[ibody]] = Bt[ibody];
      }
    }
    */
  }
}
#endif
