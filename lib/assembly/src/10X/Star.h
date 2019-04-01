// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

#ifndef TENX_STARX_H
#define TENX_STARX_H

#include "CoreTools.h"
#include "graph/Digraph.h"
#include "math/HoInterval.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/paths/ReadPathVecX.h"

void Star( const HyperBasevectorX& hb, const vec<int>& inv, const vec<Bool>& dup,
     const vec<int32_t>& bc, const ReadPathVecX& paths, 
     digraphE<vec<int>>& D, 
     vec<int>& dinv, ReadPathVec& dpaths,
     const VecIntVec& ebcx, const vec< triple<int,int,int> >& qept,
     const double MIN_ADVANTAGE, const Bool OO, 
     const Bool DJANGO, const Bool verbose = False );

void Star2( const HyperBasevectorX& hb, const vec<int>& inv, const vec<Bool>& dup,
     const vec<int32_t>& bc, const ReadPathVecX& paths, 
     digraphE<vec<int>>& D, 
     vec<int>& dinv, ReadPathVec& dpaths,
     const VecIntVec& ebcx, const vec< triple<int,int,int> >& qept,
     const double MIN_ADVANTAGE, const Bool OO, 
     const Bool DJANGO, const Bool verbose = False );

void Star3( const HyperBasevectorX& hb, const vec<int>& inv, const vec<Bool>& dup,
     const vec<int32_t>& bc, const ReadPathVecX& paths, 
     digraphE<vec<int>>& D, 
     vec<int>& dinv, ReadPathVec& dpaths,
     const VecIntVec& ebcx, const vec< triple<int,int,int> >& qept,
     const double MIN_ADVANTAGE, const Bool OO, 
     const Bool DJANGO, const Bool verbose = False );

void Star4( const HyperBasevectorX& hb, const vec<int>& inv, const vec<Bool>& dup,
     const vec<int32_t>& bc, const ReadPathVecX& paths, 
     digraphE<vec<int>>& D, 
     vec<int>& dinv, ReadPathVec& dpaths,
     const VecIntVec& ebcx, const vec< triple<int,int,int> >& qept,
     const double MIN_ADVANTAGE, const Bool OO, 
     const Bool DJANGO, const Bool verbose = False );

void ScaffoldMe( const HyperBasevectorX& hb, const vec<int>& inv, const vec<Bool>& dup,
     const vec<int32_t>& bc, const ReadPathVecX& paths, 
     digraphE<vec<int> >& D, vec<int>& dinv, 
     ReadPathVec& dpaths, const VecIntVec& ebcx, 
     const vec< triple<int,int,int> >& qept, 
     const double MIN_ADVANTAGE, const Bool OO, 
     const Bool DJANGO, 
     const int BC_REQUIRE, int BC_FLANK, int BC_IGNORE,
     const Bool require_symmetric, const Bool verbose = False);
#endif
