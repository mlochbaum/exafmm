/*

This file contain all the functions used for data rearrangement based on the octree indexing

author: Nikos Sismanis
date: Jul 2014  

*/

#include "stdio.h"
#include "stdlib.h"
#include "cilk/cilk.h"
#include "cilk/cilk_api.h"
#include "math.h"
#include "sys/time.h"
#include <string.h>
#include <stdint.h>

#define NSPLITS 64
#define DIM 3
#ifndef LDIM
#define LDIM 12
#endif
#define ADIM 12
#define NP 512

void relocateTL(float (*restrict Y), float (*restrict X), uint32_t (*restrict index), uint32_t N){
  for(uint i=0; i<N; i++){
    Y[i*LDIM:LDIM] = X[index[i]*LDIM:LDIM];
  }
}

void rearrange_dataTL(float (*restrict Y), float (*restrict X), uint32_t (*restrict Ids), int N){

  if(N>100000){

    uint64_t M = (uint32_t)ceil((float)N / (float)NP);
    for(uint64_t i=0; i<NP; i++){
      uint64_t size = ((i+1)*M < N) ? M : N - i*M;
      cilk_spawn relocateTL(&Y[i*M*LDIM], X, &Ids[i*M], size);
    }
    cilk_sync;
  }
  else{

    cilk_for(uint64_t i=0; i<N; i++){
      Y[i*LDIM:LDIM] = X[Ids[i]*LDIM:LDIM];
    }    

  }
}

void rearrangeDirect(float (*restrict Y),
		     float (*restrict X),
		     uint32_t (*restrict Ids), int N){

 cilk_for(uint i=0; i<N; i++){
    // Y[i*LDIM:LDIM] = X[Ids[i]*LDIM:LDIM];
    memcpy(&Y[i*LDIM], &X[Ids[i]*LDIM], LDIM*sizeof(float));
  }    
}

