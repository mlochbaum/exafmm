#ifndef traversal_uvwx_h
#define traversal_uvwx_h
#include "kernel.h"
#include "logger.h"
#include "thread.h"
#if EXAFMM_OPTIMIZE_TREE
#include "morton_ops.h"
#else
#include "morton_key.h"
#endif

#define U_LIST 0
#define V_LIST 1
#define W_LIST 2
#define X_LIST 3

#if EXAFMM_COUNT_KERNEL
#define countKernel(N) N++
#else
#define countKernel(N)
#endif

namespace exafmm {
  using namespace morton;
  class Traversal {
  private:
    const int nspawn;                                           //!< Threshold of NBODY for spawning new threads
    const int images;                                           //!< Number of periodic image sublevels
    // Interaction lists are stored as linked lists:
    // listOffset[icell][itype] gives the index (in lists) of the last
    // element of icell's list itype.
    // lists[ilast][0] gives the index of the element before ilast in
    // the same list (determined by cell and type).
    // lists[ilast][1] gives the cell index (in Cj0) of the list cell.
    // lists[ilast][2] gives that cell's periodic key.
    int (* listOffset)[4];
    int (* lists)[3];
#if EXAFMM_COUNT_KERNEL
    // Counters for kernel calls of each type
    real_t numP2P; // U_LIST
    real_t numM2L; // V_LIST
    real_t numM2P; // W_LIST
    real_t numP2L; // X_LIST
#endif
    // Iterators of target and source cells
    C_iter Ci0;
    C_iter Cj0;

  private:
#if EXAFMM_COUNT_LIST
    //! Accumulate interaction list size of cells
    //! Combines V, W, and X into M2L; this should be fixed.
    void countList(C_iter Ci, C_iter Cj, bool mutual, int type) {
      bool isP2P = type==U_LIST;
      if (isP2P) Ci->numP2P++;                                  // If P2P, increment P2P counter of target cell
      else Ci->numM2L++;                                        // Else, increment M2L counter of target cell
      if (mutual) {                                             // If mutual interaction in on
	if (isP2P) Cj->numP2P++;                                //  If P2P, increment P2P counter of source cell
	else Cj->numM2L++;                                      //  Else, increment M2L counter of source cell
      }                                                         // End if for mutual interaction
    }
#else
    void countList(C_iter, C_iter, bool, bool) {}
#endif

#if EXAFMM_USE_WEIGHT
    //! Accumulate interaction weights of cells
    void countWeight(C_iter Ci, C_iter Cj, bool mutual, real_t weight) {
      Ci->WEIGHT += weight;                                     // Increment weight of target cell
      if (mutual) Cj->WEIGHT += weight;                         // Increment weight of source cell
    }
#else
    void countWeight(C_iter, C_iter, bool, real_t) {}
#endif

    // Add jcell to icell's list itype
    void setList(int itype, int icell, int jcell, int periodicKey, int & numLists) {
      lists[numLists][0] = listOffset[icell][itype]; // Index of previous
      lists[numLists][1] = jcell;                    // Cell value
      lists[numLists][2] = periodicKey;              // Periodic key
      listOffset[icell][itype] = numLists;           // Index of current
      numLists++;
    }

    // Find whether two cells are adjacent
    int isAdj(C_iter Ci, C_iter Cj) {
      uint64_t ki = Ci->ICELL, kj = Cj->ICELL;
      return morton::adj(ki,kj);
    }

    // Get the periodic key of Cj with respect to Ci
    int getPeriodicKey(C_iter Ci, C_iter Cj) {
      uint64_t ki = Ci->ICELL, kj = Cj->ICELL;
      return morton::getPeriodicKey(ki,kj);
    }

    // Get the level (height) of cell Ci
    int getLevel(C_iter Ci) {
      return morton::getLevel(Ci->ICELL);
    }

    // Get 3-D index from periodic key
    ivec3 getPeriodicIndex(int key) {
      ivec3 iX;
      iX[0] = key % 3;
      iX[1] = (key / 3) % 3;
      iX[2] = key / 9;
      iX -= 1;   // {0,1,2} -> {-1,0,1}
      return iX;
    }

    // Generate all interaction lists
    void setLists(Cells & icells) {
      int numCells = icells.size();
      // colleagues will store all adjacent boxes of equal or greater height
      int (* colleagues)[27] = new int[numCells][27]();
      int *parentColls, *childColls;
      for (int i=0; i<numCells; i++) {
        for (int j=0; j<4; j++) {
          listOffset[i][j] = -1;
        }
      }
      int numLists = 0;
      // Initialize root with no neighbors
      for (int j=0; j<27; j++) { colleagues[0][j] = -1; }
      // Traverse in order: parents are visited before children
      for (int icell=1; icell<numCells; icell++) {
        C_iter Ci = Ci0 + icell;
        int iparent = Ci->IPARENT;
        C_iter Cp = Ci0 + iparent;
        parentColls = colleagues[iparent];
        childColls  = colleagues[icell];
        // loop over parent colleagues
        for (int j=0; j<27; j++) {
          int jcell = parentColls[j];
          childColls[j] = -1;
          if (jcell == -1) { continue; }
          C_iter Cj = Cj0 + jcell;
          if (Cj->NCHILD == 0) { // If j is a leaf:
            if (isAdj(Ci, Cj)) { // If it is still adjacent, add to colleagues
              childColls[j] = parentColls[j];
            } else { // If it is no longer adjacent to i, add to i's X list
              setList(X_LIST, icell, jcell, j, numLists);
              setList(W_LIST, jcell, icell, 26-j, numLists);
            }
          } else { // If j has children:
            // Add all children which are not colleagues to i's V list
            // Add colleague children to i's colleague list
            for (int k=0; k<Cj->NCHILD; k++) {
              int kcell = Cj->ICHILD+k;
              if (!isAdj(Ci, Cj0+kcell)) {
                setList(V_LIST, icell, kcell, j, numLists);
              } else {
                childColls[j] = kcell;
              }
            }
          }
        }
        // Add parent's children to colleague list
        for (int k=0; k<Cp->NCHILD; k++) {
          int kcell = Cp->ICHILD+k;
          C_iter Ck = Ci0 + kcell;
          if (icell != kcell) {
            childColls[getPeriodicKey(Ci, Ck)] = kcell;
          }
        }
        // If i is a leaf, add its colleagues to appropriate lists
        if (Ci->NCHILD==0) {
          for (int j=0; j<27; j++) {
            int jcell = childColls[j];
            if (jcell == -1) { continue; }
            C_iter Cj = Cj0 + jcell;
            if (Cj->NCHILD == 0) { // If j is a leaf, add to U
              int periodicKey = getPeriodicKey(Ci, Cj);
              setList(U_LIST, icell, jcell, periodicKey, numLists);
              // Since larger cells don't inspect smaller ones, the
              // smaller cell must add itself to the larger cell's U list.
              if (getLevel(Cj) < getLevel(Ci)) {
                setList(U_LIST, jcell, icell, 26-periodicKey, numLists);
              }
            }
          }
        }
      }
      delete[] colleagues;
    }

    // After lists have been constructed, execute the corresponding kernels
    void listBasedTraversal(int numCells, vec3 cycle, bool mutual, real_t remote) {
      for (int itype=0; itype<4; itype++) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (int icell=0; icell<numCells; icell++) {
          C_iter Ci = Ci0 + icell;
          // Last cell in list
          int ilast = listOffset[icell][itype];
          while (ilast >= 0) {
            int jcell = lists[ilast][1];
            int periodicKey = lists[ilast][2];
            C_iter Cj = Cj0 + jcell;
            ivec3 pX = getPeriodicIndex(periodicKey);
            for (int d=0; d<3; d++) {
              kernel::Xperiodic[d] = pX[d] * cycle[d];
            }

            switch (itype) {
              case U_LIST: kernel::P2P(Ci, Cj, mutual); countKernel(numP2P); break;
              case V_LIST: kernel::M2L(Ci, Cj, mutual); countKernel(numM2L); break;
              case W_LIST: kernel::M2P(Ci, Cj, mutual); countKernel(numM2P); break;
              case X_LIST: kernel::P2L(Ci, Cj, mutual); countKernel(numP2L); break;
            }
            countList(Ci, Cj, mutual, itype);
            countWeight(Ci, Cj, mutual, remote);
            // Move to previous cell in list
            ilast = lists[ilast][0];
          }
        }
      }
    }

  public:
    //! Constructor
    Traversal(int _nspawn, int _images) :                       // Constructor
      nspawn(_nspawn), images(_images)                          // Initialize variables
#if EXAFMM_COUNT_KERNEL
      , numP2P(0), numM2L(0), numM2P(0), numP2L(0)
#endif
    {}

#if EXAFMM_COUNT_LIST
    //! Initialize size of P2P and M2L interaction lists per cell
    void initListCount(Cells & cells) {
      for (C_iter C=cells.begin(); C!=cells.end(); C++) {       // Loop over cells
	C->numP2P = C->numM2L = 0;                              //  Initialize size of interaction list
      }                                                         // End loop over cells
    }
#else
    void initListCount(Cells) {}
#endif

#if EXAFMM_USE_WEIGHT
    //! Initialize interaction weights of bodies and cells
    void initWeight(Cells & cells) {
      for (C_iter C=cells.begin(); C!=cells.end(); C++) {       // Loop over cells
	C->WEIGHT = 0;                                          //  Initialize cell weights
	if (C->NCHILD==0) {                                     //  If leaf cell
	  for (B_iter B=C->BODY; B!=C->BODY+C->NBODY; B++) {    //   Loop over bodies in cell
	    B->WEIGHT = 0;                                      //    Initialize body weights
	  }                                                     //   End loop over bodies in cell
	}                                                       //  End if for leaf cell
      }                                                         // End loop over cells
    }
#else
    void initWeight(Cells) {}
#endif

    //! Evaluate P2P and M2L using list based traversal
    void traverse(Cells & icells, Cells & jcells, vec3 cycle, bool dual, bool mutual, real_t remote=1) {
      if (icells.empty() || jcells.empty()) return;
      logger::startTimer("Traverse");
      logger::initTracer();
      Ci0 = icells.begin();
      Cj0 = jcells.begin();
      kernel::Xperiodic = 0;

      int numCells = icells.size();
      // Do not allow mutual interaction
      mutual = false;
      // Initialize lists
      listOffset = new int [numCells][4]();
      lists = new int [(216+27)*numCells][3]();
      // Construct interaction lists
      setLists(icells);

      // Perform interactions based on lists
      listBasedTraversal(numCells, cycle, mutual, remote);

      delete[] listOffset;
      delete[] lists;

      logger::stopTimer("Traverse");
      logger::writeTracer();
    }

    // Print number of each kernel type used
    void printTraversalData() {
#if EXAFMM_COUNT_KERNEL
      if (logger::verbose) {
#define PRINT(X2X) \
     std::setw(logger::stringLength) << std::left << #X2X " calls" << " : " \
  << std::setprecision(0) << std::fixed << num##X2X << std::endl
        std::cout << "--- Traversal stats --------------" << std::endl
                  << PRINT(P2P)
                  << PRINT(M2L)
                  << PRINT(M2P)
                  << PRINT(P2L);
#undef PRINT
      }
#endif // EXAFMM_COUNT_KERNEL
    }

#if EXAFMM_COUNT_LIST
    void writeList(Cells cells, int mpirank) {
      std::stringstream name;                                   // File name
      name << "list" << std::setfill('0') << std::setw(6)       // Set format
	   << mpirank << ".dat";                                // Create file name for list
      std::ofstream listFile(name.str().c_str());               // Open list log file
      for (C_iter C=cells.begin(); C!=cells.end(); C++) {       // Loop over all lists
	listFile << std::setw(logger::stringLength) << std::left//  Set format
		 << C->ICELL << " " << C->numP2P << " " << C->numM2L << std::endl; // Print list size
      }                                                         // End loop over all lists
    }
#else
    void writeList(Cells, int) {}
#endif


#if EXAFMM_OPTIMIZE_TREE
  private:
    // Constants of proportionality for kernel functions
    double K_P2P, K_L2L, K_M2P, K_P2L, K_M2M;
    // Morton keys of bodies
    uint64_t* keys;
    // Cost of subdivision at each depth
    double (* costs)[MAX_DEPTH];

  public:
    void initialize_costs(double P2P, double L2L, double M2P, double P2L, double M2M) {
      K_P2P = P2P; K_L2L = L2L; K_M2P = M2P; K_P2L = P2L; K_M2M = M2M;
    }

    void optimize_tree(Cells cells) {
      int numCells = cells.size();
      int *depth = new int [cells[0].NBODY];
      int *subdivcost = new int [numCells]; // Cost of optimal subdivision
      for (int icell=0; icell<numCells; icell++) {
        C_iter Ci = Ci0 + icell; Morton Mi = expandMorton(Ci);
        int ilast; int numBodies = Ci->NBODY;
        if (numBodies == 0) { subdivcost[icell] = 0; continue; }
        keys = new uint64_t [numBodies];
        costs = new double [numBodies][MAX_DEPTH];
        // Initialize body keys
        for (int i=0; i<numBodies; i++) {
          keys[i] = (Ci0->BODY + i)->KEY;
          for (int j=0; j<MAX_DEPTH; j++) costs[i][j]=0;
        }
        // Add all costs
        // Traverse U and W lists
#define ADD_LIST(LIST, L) \
        ilast = listOffset[icell][LIST];                 \
        while (ilast >= 0) {                             \
          C_iter Cj = Cj0+lists[ilast][1];               \
          add_cost_##L(Mi, expandMorton(Cj), Cj->NBODY, numBodies, 0); \
          ilast = lists[ilast][0];                       \
        }
        ADD_LIST(U_LIST, U); ADD_LIST(W_LIST, W);
#undef ADD_LIST
        // Find the optimal subdivision depth of each body in Ci and
        // add to depth
        int p = 0;
        subdivcost[icell] = find_subdivision_cost(depth+Ci->ICELL, Mi.h, &p, numBodies);
        delete[] costs;
        delete[] keys;
      }
      for (int icell=0; icell<numCells; icell++) {
        if (subdivcost[icell] < 0) {
          // TODO subdivide cell
        }
      }
      delete[] depth;
      delete[] subdivcost;
    }

  private:
    /* Cost chart
        UV      L2L - c*d*P2P
        UW    c*M2P - c*d*P2P
        UX    d*P2L - c*d*P2P
        UN          - c*d*P2P
        WV      L2L - c * M2P
        WN          - c * M2P
    */

    void add_cost_W(Morton B, Morton D, int size_D, int end, int p) {
      add_cost_W(B, D, size_D, end, p, 0);
    }
    void add_cost_W(Morton B, Morton D, int size_D, int end, int p, int o) {
      int sh = (MAX_DEPTH-B.h)*3; // Amount to shift to test that p is still in B
      do {
        int p0=p; uint64_t pk = keys[p]; double cost = 0;

        // Find subdivision depth n and move p past corresponding box
        uint8_t n = nonadjno(pk,B,parent(D)), dh = D.h-B.h+1;
        n = (n < dh) ? n : dh;
        uint8_t k = (MAX_DEPTH-n)*3; while (p < end && keys[p]>>k == pk>>k) p++;

        // Compute costs and add to cost array
        cost -= K_M2P*(p-p0); // For both W->V and W->N
        if (n > dh) cost += K_L2L; // For W->V only
        costs[p0][n+o] += cost;

      } while (p < end && keys[p]>>sh == B.k);
    }

    void add_cost_U(Morton B, Morton D, int size_D, int end, int p) {
      int sh = (MAX_DEPTH-B.h)*3; // Amount to shift to test that p is still in B
      do {
        int p0=p; uint64_t pk = keys[p]; double cost = 0;

        // Find subdivision depth n and move p past corresponding box
        uint8_t n = nonadjno(pk,B,D);
        uint8_t k = (MAX_DEPTH-n)*3; while (p < end && keys[p]>>k == pk>>k) p++;

        // Compute costs and add to cost array
        cost -= K_P2P*(p-p0)*size_D; // All transitions (GPU)
        switch(whichlist(B,D)) {
          case V_list: { cost += K_L2L; break; }
          case X_list: { cost += K_M2P * (p-p0); break; }
          case W_list: { cost += K_P2L * size_D;
                         add_cost_W(SHRINK(B,pk,n),D,end,p,n); }
          case N_list: break;
          default: break; // not possible
        }
        costs[p0][n] += cost;

      } while (p < end && keys[p]>>sh == B.k);
    }

    double find_subdivision_cost(int* depth, int d, int* p, int end) {
      uint64_t p0 = keys[*p];
      int n = (MAX_DEPTH-d)*3;
      // q is the first body with a different Morton key than p
      int q = *p; while (q < end && keys[q] == p0) q++;
      if (q < end && keys[q]>>n == p0>>n) {
        double cost = K_M2M + costs[*p][d];
        do {
          cost += find_subdivision_cost(depth, d+1, p, end);
        } while (keys[*p]>>n == p0>>n);
        if (cost >= 0) {
          depth[*p] = d;
          return 0;
        } else {
          return cost;
        }
      } else {
        double cost = 0, c = 0;
        for (; d < MAX_DEPTH; d++) {
          c += K_M2M + costs[*p][d];
          if (c < cost) { depth[*p] = d; cost = c; }
        }
        *p = q;
        return cost;
      }
    }
#endif // EXAFMM_OPTIMIZE_TREE
  };
}
#endif
