#ifndef morton_ops_h
#define morton_ops_h
#include "types.h"
#include "morton_key.h"

#define MAX_DEPTH 21
#define SHRINK(B, p, d) ((Morton){B.h+d, p>>(21-B.h-d)*3})
#define FLS(a) (64-__builtin_clzll(a))

namespace exafmm {
  namespace morton {
    typedef struct { uint8_t h; uint64_t k; } Morton;

    // Convert back and forth between struct representation
    Morton expandMorton(uint64_t m) {
      Morton M;
      /*
      int n = FLS(m);
      M.h = n/3;
      M.k = m & ((1<<n-1) - 1);
      */
      int n = getLevel(m);
      M.h = n;
      M.k = m - getOffset(n);
      return M;
    }
    Morton expandMorton(C_iter C) {
      return expandMorton(C->ICELL);
    }

    uint64_t compressMorton(Morton M) {
      // return M.k + (1<<(3*M.h));
      return M.k + getOffset(3*M.h);
    }


    // Parent of box B
    Morton parent(Morton B) {
      return (Morton){ B.h-1, B.k>>3 };
    }

    // Determine whether two boxes are colleagues
    char collk(uint64_t b, uint64_t c) {
      if (b == c) return 0;
      for (int i=0; i<3; i++) {
        int diff = unint(b>>i) - unint(c>>i);
        if (diff > 1 || diff < -1) return 0;
      }
      return 1;
    }
    char coll(Morton B, Morton C) {
      return B.h == C.h && collk(B.k, C.k);
    }

    // Determine whether two boxes are adjacent
    char adjk(uint64_t b, uint64_t c, int dh) {
      if (dh == 0) return collk(b, c);
      if (b>>(3*dh) == c) return 0;
      for (int i=0; i<3; i++) {
        int uB = unint(b>>i), uC = unint(c>>i);
        int diff = (uB >> dh) - uC;
        if (diff == 0) continue;
        if (diff != 1 && diff != -1) return 0;
        if (uB + (diff==-1)  &  (1<<dh) - 1) return 0;
      }
      return 1;
    }
    char adj(Morton B, Morton C) {
      int dh = B.h - C.h;
      return (dh >= 0) ? adjk(B.k, C.k, dh)
                       : adjk(C.k, B.k, -dh);
    }

    typedef enum {
      U_list,
      V_list,
      W_list,
      X_list,
      N_list
    } List_Name;

    // Which list of B the box C is in
    char wlist(uint64_t b, uint64_t c, int dh) {
      return collk(b, c >> dh*3) && adjk(c>>3, b, dh-1);
    }
    List_Name whichlist(Morton B, Morton C) {
      if (B.h == C.h) {
        return (B.k==C.k || collk(B.k,C.k)) ? U_list
          : collk(B.k>>3, C.k>>3) ? V_list : N_list;
      } else if (B.h < C.h) {
        return adj(B,C) ? U_list
          : wlist(B.k, C.k, C.h-B.h) ? W_list : N_list;
      } else {
        return adj(C,B) ? U_list
          : wlist(C.k, B.k, B.h-C.h) ? X_list : N_list;
      }
    }

    // Assume b is in A and A is adjacent to C.
    uint8_t nonadjno(uint64_t b, Morton A, Morton C) {
      uint64_t c = C.k;
      uint8_t dh = 21 - C.h;
      uint8_t r = 21;
      for (int i=0; i<3; i++) {
        int uB = unint(b>>i), uC = unint(c>>i);
        int diff = (uB >> dh) - uC;
        uint8_t hh = (diff == 0) ? 21 :
          21 - FLS((uint64_t)((diff==1 ? 0 : -1) ^ uB) & ((1<<dh)-1));
        r = (hh < r) ? hh : r;
      }
      return r - A.h + 1;
    }

    // Subdivision number of p in A with respect to B
    uint8_t subdivno(uint64_t p, Morton A, Morton B) {
      switch(whichlist(A,B)) {
        case U_list: { uint8_t n = nonadjno(p,A,B);
                       return n + subdivno(p, SHRINK(B,p,n), B); }
        case V_list:
        case X_list: return 1;
        case W_list: { uint8_t n = nonadjno(p,A,parent(B)), dh = B.h-A.h+1;
                       return (n < dh) ? n : dh; }
        case N_list: return 0;
      }
    }
  }
}
#endif
