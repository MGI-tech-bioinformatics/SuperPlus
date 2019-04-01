// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/Lines.h"
#include "10X/Heuristics.h"
#include "10X/IntIndex.h"
#include "10X/LineOO.h"
#include "10X/PlaceReads.h"
#include "10X/Star.h"
#include "10X/Super.h"
#include "10X/Gap.h"

#define USE_PIVOT 1

Bool TooMuchDupContent(const digraphE<vec<int> >& D,
        const vec<vec<vec<int> > >& line_L,
        const vec<vec<vec<int> > >& line_R,
        const vec<int>& inv)
{
    set<int> edge_L, edge_R;
    for ( int j = 0; j < line_L.isize( ); j++ )
    {   
        for ( int k = 0; k < line_L[j].isize( ); k++ )
        {   
            for ( int l = 0; l < line_L[j][k].isize( ); l++ )
            {    int d = line_L[j][k][l];
                if ( D.O(d)[0] < 0 ) continue;
                for ( int m = 0; m < D.O(d).isize( ); m++ ) {
                    edge_L.insert(D.O(d)[m]);
                }
            }
        }
    }
    for ( int j = 0; j < line_R.isize( ); j++ )
    {   
        for ( int k = 0; k < line_R[j].isize( ); k++ )
        {   
            for ( int l = 0; l < line_R[j][k].isize( ); l++ )
            {    int d = line_R[j][k][l];
                if ( D.O(d)[0] < 0 ) continue;
                for ( int m = 0; m < D.O(d).isize( ); m++ ) {
                    edge_R.insert(D.O(d)[m]);
                }
            }
        }
    }
    int cnt = 0;
    for (auto & edge : edge_L) {
        if (edge_R.find(edge) != edge_R.end() ||
            edge_R.find(inv[edge]) != edge_R.end()) {
            ++cnt;
        }
    }
    size_t total_cnt = Min(edge_L.size(), edge_R.size());
    double ratio = double(cnt) / total_cnt;
    if (ratio >= 0.5) {
        // #pragma omp critical
        // cout << "WTF!: " << ratio << endl;
        return True;
    }
    return False;
}

Bool TooMuchDupContentSingle(const digraphE<vec<int> >& D,
        const vec<vec<vec<int> > >& line_L,
        const vec<int>& dinv) {
    set<int> edge_L;
    for ( int j = 0; j < line_L.isize( ); j++ )
    {   
        for ( int k = 0; k < line_L[j].isize( ); k++ )
        {   
            for ( int l = 0; l < line_L[j][k].isize( ); l++ )
            {    
                int d = line_L[j][k][l];
                ForceAssertGe(d, 0);
                ForceAssertGe(D.O(d).size(), 0);
                if ( D.O(d)[0] < 0 ) continue;
                edge_L.insert(d);
            }
        }
    }
    int cnt = 0;
    for (auto &edge : edge_L) {
        if (edge  == dinv[edge]) {
            ++cnt;
            continue;
        }
        if (edge_L.find(dinv[edge]) != edge_L.end()) {
            ++cnt;
        }
    }
    double ratio = double(cnt) / edge_L.size();
    if (cnt > 0) {
#pragma omp critical
        cout << "nimei: " << cnt << ':' << edge_L.size() << endl;
    }
    if (ratio >= 0.25) {
        return True;
    }
    return False;
}

void Star( const HyperBasevectorX& hb, const vec<int>& inv, const vec<Bool>& dup,
     const vec<int32_t>& bc, const ReadPathVecX& paths, 
     digraphE<vec<int>>& D, vec<int>& dinv, 
     ReadPathVec& dpaths, const VecIntVec& ebcx, 
     const vec< triple<int,int,int> >& qept, 
     const double MIN_ADVANTAGE, const Bool OO, 
     const Bool DJANGO, const Bool verbose )
{
     cout << Date( ) << ": begin Star" << endl;

     // Heuristics.

     const int MIN_STAR = 5000;
     const double MAX_CN_DIFF = 0.5;
     const int MAX_VIEW = 10;
     int MAX_RIGHTS = 6;
     if (OO) MAX_RIGHTS = 4;
     // const int MIN_BAR_TO = 2000;
     const int MIN_BAR_TO = 2000;
     const int BC_VIEW = 50000;

     // Place reads.

     cout << Date( ) << ": placing reads" << endl;
     PlaceReads( hb, paths, dup, D, dpaths, True, False );
     cout << Date( ) << ": making index" << endl;
     IntIndex dpaths_index( dpaths, D.E( ) );

     // Find lines.

     cout << Date( ) << ": finding lines" << endl;
     vec<vec<vec<vec<int>>>> dlines;
     FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH, True );

     // Compute kmers.

     vec<int> kmers( hb.E( ) );
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
          kmers[e] = hb.Kmers(e);

     // Compute coverage.

     vec<int> llens;
     GetLineLengths( hb, D, dlines, llens );
     vec<double> cov;
     {    vec<vec<pair<int,int>>> lbp;
          BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, 0 );
          MasterVec<SerfVec<pair<int,int>>> lbpx;
          for ( auto x : lbp )
          {    SerfVec<pair<int,int>> y( x.size( ) );
               for ( int j = 0; j < x.isize( ); j++ )
                    y[j] = x[j];
               lbpx.push_back(y);    }
          LineCN( kmers, lbpx, D, dlines, llens, cov );    }

     // Start passes.

     int sp = 0;
     while(1)
     {    cout << Date( ) << ": start star pass " << ++sp << endl;
          cout << Date( ) << ": N50 line = " 
             << ToStringAddCommas( LineN50( hb, D, dlines ) ) << endl;

          // Introduce barcode-only gaps.  First step: find line neighbors.

          cout << Date( ) << ": get line ancillary data" << endl;
          vec<int> llens, to_left, to_right, linv;
          vec< vec< pair<int,int> > > lhood;
          GetLineLengths( hb, D, dlines, llens );
          LineProx( hb, inv, ebcx, D, dinv, dlines, qept, lhood );
          D.ToLeft(to_left), D.ToRight(to_right);
          LineInv( dlines, dinv, linv );

          // Get barcode positions.

          cout << Date( ) << ": getting barcode positions" << endl;
          vec<vec<pair<int,int>>> lbp;
          BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, BC_VIEW );
     
          // Go through long lines and find their nearest right neighbors.

          cout << Date( ) << ": begin star stage" << endl;
          vec< pair<int,int> > stars;
          vec<int> l1s;
          for ( int L1 = 0; L1 < dlines.isize( ); L1++ )
          {    if ( llens[L1] < MIN_STAR ) continue;
               if ( linv[L1] == L1) continue;
               int v = to_right[ dlines[L1].back( )[0][0] ];
               if ( v < 0 ) { continue; }
               if ( D.From(v).nonempty( ) || !D.To(v).solo( ) ) continue;
               l1s.push_back(L1);    }
          vec<int> l1slen( l1s.isize( ) );
          for ( int i = 0; i < l1s.isize( ); i++ )
               l1slen[i] = llens[ l1s[i] ];
          ReverseSortSync( l1slen, l1s );
          cout << Date( ) << ": start star loop on " << l1s.size( )
               << " L1 values" << endl;
          map< vec<int>, double > memory;
          double clock = WallClockTime( );
          const int batch = Max( 1, l1s.isize( )/500 );
          #pragma omp parallel for schedule(dynamic, batch)
          for ( int li = 0; li < l1s.isize( ); li++ )
          {    int L1 = l1s[li];
               vec< triple<int,int,int> > M;

               // Set up logging.

               ostringstream* outp;
               if (verbose) outp = new ostringstream;
               auto DumpOut = [&]
               {    if (verbose)
                    {
                         #pragma omp critical
                         {    cout << outp->str( );    }    
                         delete outp;    }    };

               // Find right neighbors of L1.

               vec<int> rights;
               for ( int i = 0; i < Min( MAX_VIEW, lhood[L1].isize( ) ); i++ )
               {    
                    int L2 = lhood[L1][i].second;
                    ForceAssert(linv[L2] != L1);
                    ForceAssert(L2 != L1);
                    if (linv[L2] == L2) continue;
                    if ( AbsDiff( cov[L1], cov[L2] ) > MAX_CN_DIFF ) continue;
                    if ( llens[L2] < MIN_BAR_TO ) continue;
                    vec<double> scores;
                    scores.push_back( MemoryScoreOrder( 
                         { L2, L1 }, lbp, llens, M, memory ) );
                    scores.push_back( MemoryScoreOrder( 
                         { linv[L2], L1 }, lbp, llens, M, memory ) );
                    scores.push_back( MemoryScoreOrder( 
                         { L1, L2 }, lbp, llens, M, memory ) );
                    scores.push_back( MemoryScoreOrder( 
                         { L1, linv[L2] }, lbp, llens, M, memory ) );
                    vec<int> ids( 4, vec<int>::IDENTITY );
                    SortSync( scores, ids );
                    if (verbose) {
                        (*outp) << "tmp score: " << printSeq(scores) 
                            << " tmp index: " << printSeq(ids) << endl;
                    }
                    double ad = scores[1] - scores[0];
                    if ( ad < MIN_ADVANTAGE ) continue;
                    if ( ids[0] <= 1 ) continue;
                    if ( ids[0] == 2 ) rights.push_back(L2);
                    else rights.push_back( linv[L2] );    }

               // Cap number of rights.  Go with the biggest lines.

               if ( rights.isize( ) > MAX_RIGHTS )
               {    vec<int> lens( rights.size( ) );
                    for ( int i = 0; i < rights.isize( ); i++ )
                         lens[i] = llens[ rights[i] ];
                    ReverseSortSync( lens, rights );
                    rights.resize(MAX_RIGHTS);    }
               if (verbose)
               {    (*outp) << "\nlooking right from L" << L1 << endl;
                    (*outp) << "see " << printSeq(rights) << endl;    }
               if ( rights.empty( ) )
               {    DumpOut( );
                    continue;    }
     
               // Find leftmost right.

               int R = -1;
               if ( rights.size( ) == 1 ) R = rights[0];
               else if ( !DJANGO )
               {    vec<int> brights;
                    double ad;
                    if (OO)
                    {    
#if !USE_PIVOT
                        MemoryOrderAndOrientN(rights, lbp, llens, 
                                linv, M, memory, ad, brights );    
#else                        
                        MemoryOrderAndOrientN( L1, rights, lbp, llens, 
                               linv, M, memory, ad, brights );    
#endif
                              
                    }
                    else {
#if !USE_PIVOT
                        MemoryOrderN(rights, lbp, llens, M, memory, ad, brights );
#else
                        MemoryOrderN( L1, rights, lbp, llens, M, memory, ad, brights );
#endif
                    }
                    if (verbose)
                    {    (*outp) << "\nlooking right from L" << L1 << endl
                              << "see " << printSeq(rights) << endl
                              << "winner = L" << brights[0]
                              << " [" << ad << "]" << endl;    }
                    if ( ad < MIN_ADVANTAGE )
                    {    DumpOut( );
                         continue;    }
                    R = brights[0];    }
               else
               {    vec<int> lens( rights.size( ) );
                    for ( int i = 0; i < rights.isize( ); i++ )
                         lens[i] = llens[ rights[i] ];
                    ReverseSortSync( lens, rights );
                    while( rights.size( ) > 1 )
                    {    vec<int> brights;
                         double ad;
                         if (OO)
                         {    
#if !USE_PIVOT
                             MemoryOrderAndOrientN(rights, lbp, llens, 
                                     linv, M, memory, ad, brights );
#else
                             MemoryOrderAndOrientN( L1, rights, lbp, llens, 
                                     linv, M, memory, ad, brights );
#endif
                                  
                         }
                         else 
                         {    
#if !USE_PIVOT
                             MemoryOrderN( rights, lbp, llens, M, memory, ad, brights );
#else
                             MemoryOrderN( L1, rights, lbp, llens, M, memory, ad, brights );
#endif
                         }
                         if (verbose)
                         {    (*outp) << "\nlooking right from L" << L1 << endl
                                   << "see " << printSeq(rights) << endl
                                   << "winner = L" << brights[0]
                                   << " [" << ad << "]" << endl;    }
                         if ( ad >= MIN_ADVANTAGE )
                         {    R = brights[0];
                              break;    }
                         rights.pop_back( );    }    }

               // Check for done and otherwise save.

               if ( R == -1 )
               {    DumpOut( );
                    continue;    }
               int w = to_left[ dlines[R].front( )[0][0] ]; 
               if ( w < 0 )
               {    DumpOut( );
                    continue;    }
               if ( D.To(w).nonempty( ) || !D.From(w).solo( ) )
               {    DumpOut( );
                    continue;    }
               if (verbose) (*outp) << "saving " << L1 << " --> " << R << endl;
               DumpOut( );
               #pragma omp critical
               {    stars.push( L1, R );    }    }
     
          // Introduce joins.
     
          cout << Date( ) << ": done, time used = " << TimeSince(clock) << endl;
          Sort(stars);
          vec<pair<int,int>> stars2;
          for ( int i = 0; i < stars.isize( ); i++ )
          {    int L = stars[i].first, R = stars[i].second;
               if (verbose) cout << "\nexamining " << L << " --> " << R << endl;
               int ri = BinPosition( stars, make_pair( linv[R], linv[L] ) );
               if ( !DJANGO && ri < i ) continue;
               int d1 = dlines[L].back( )[0][0], d2 = dlines[R].front( )[0][0];
               int v = to_right[d1], w = to_left[d2];
               int rv = to_right[ dinv[d2] ], rw = to_left[ dinv[d1] ];
               if (DJANGO)
               {    if ( D.From(v).nonempty( ) || D.From(rv).nonempty( )
                         || D.To(w).nonempty( ) || D.To(rw).nonempty( ) )
                    {    continue;    }    }
               dinv.push_back( D.E( ) + 1, D.E( ) );
               D.AddEdge( v, w, {-2} );
               D.AddEdge( rv, rw, {-2} );
               stars2.push( L, R );
               stars2.push( linv[R], linv[L] );
               if (verbose)
               {    cout << "linking from L" << L << " to L" << R << endl;
                    cout << "linking from L" << linv[R] << " to L" << linv[L]
                         << endl;    }    }
          cout << Date( ) << ": made " << stars2.size( ) << " star joins" << endl;
          if ( stars2.empty( ) ) break;
          
          // Find lines of lines.

          vec<vec<int>> linelines;
          vec<int> ll_inv;
          vec<Bool> touched( dlines.size( ), False );
          {    Sort(stars2);
               vec<pair<int,int>> stars2r( stars2.size( ) );
               for ( int i = 0; i < stars2.isize( ); i++ )
                    stars2r[i] = make_pair( stars2[i].second, stars2[i].first );
               Sort(stars2r);
               for ( int i = 0; i < stars2.isize( ); i++ )
               {    int L = stars2[i].first, R = stars2[i].second;
                    if ( touched[L] || touched[R] ) continue;
                    vec<int> x = {L,R};
                    touched[L] = touched[R] = True;
                    Bool circle = False;
                    while(1)
                    {    int p = BinPosition1( stars2, x.back( ) );
                         if ( p < 0 ) break;
                         int R = stars2[p].second;
                         if ( R == x.front( ) ) 
                         {    circle = True;
                              break;    }
                         x.push_back(R);
                         touched[R] = True;    }
                    if ( !circle )
                    {    while(1)
                         {    int p = BinPosition1( stars2r, x.front( ) );
                              if ( p < 0 ) break;
                              int L = stars2r[p].second;
                              // if ( L == x.back( ) ) 
                              // {    circle = True;
                              //      break;    }
                              x.push_front(L);
                              touched[L] = True;    }    }
                    vec<int> y;
                    for ( int j = x.isize( ) - 1; j >= 0; j-- )
                    {    y.push_back( linv[ x[j] ] );
                         touched[ linv[ x[j] ] ] = True;    }
                    vec<int> xs(x), ys(y);
                    Sort(xs), Sort(ys);
                    if (circle)
                    {    x.push_back( x[0] ), y.push_back( y[0] );    }
                    int nll = linelines.size( );
                    linelines.push_back(x);
                    if ( ys != xs ) 
                    {    linelines.push_back(y);
                         ll_inv.push_back( nll+1, nll );    }
                    else ll_inv.push_back(nll);    }    }

          // Update lines and parallel data structures.

          cout << Date( ) << ": updating lines" << endl;
          vec<vec<int>> gap = { { } };
          for ( int i = 0; i < linelines.isize( ); i++ )
          {    const vec<int>& x = linelines[i];
               vec<vec<vec<int>>> y;
               if ( ll_inv[i] == i - 1 )
               {    y = dlines.back( );
                    y.ReverseMe( );
                    for ( int j = 0; j < y.isize( ); j++ )
                    for ( int k = 0; k < y[j].isize( ); k++ )
                    for ( int l = 0; l < y[j][k].isize( ); l++ )
                         y[j][k][l] = dinv[ y[j][k][l] ];    }
               else {
                   y = dlines[ x[0] ];
                   for ( int j = 1; j < x.isize( ); j++ )
                   {    y.push_back(gap);
                       if ( x[j] == x[0] ) y.push_back( dlines[ x[j] ][0] );
                       else y.append( dlines[ x[j] ] );    }
               }
               
               dlines.push_back(y);
               int64_t lensum = 0;
               double cnsum = 0.0;
               for ( int j = 0; j < x.isize( ); j++ )
               {    lensum += llens[ x[j] ];
                    cnsum += llens[ x[j] ] * cov[ x[j] ];    }
               cov.push_back( cnsum/lensum );
               touched.push_back(False);    }
          EraseIf( dlines, touched );
          EraseIf( cov, touched );    }

     // Clean up.

     cout << Date( ) << ": cleaning up at end of Star" << endl;
     // added later.
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );    }



// Star2
// 1. given L, find its lefts and rights
// 2. take its largest right R
// 3. get lefts and rights of R
// 4. must have L in lefts(R)
// 5. any member of rights(L) >= 10 kb must be present in R+lefts(R)+rights(R)
// 6. any member of lefts(R) >= 10 kb must be present in L+lefts(L)+rights(L)
// 7. pool all, drop any lefts of L and rights of R
// 8. drop tiny stuff until <= 5
// 9. O&O, dropping tiny stuff until confident
// 10. save
// 11. traverse largest by Min(L,R) and don’t join if touched

void Star2( const HyperBasevectorX& hb, const vec<int>& inv, const vec<Bool>& dup,
     const vec<int32_t>& bc, const ReadPathVecX& paths, 
     digraphE<vec<int>>& D, vec<int>& dinv, 
     ReadPathVec& dpaths, const VecIntVec& ebcx, 
     const vec< triple<int,int,int> >& qept, 
     const double MIN_ADVANTAGE, const Bool OO, 
     const Bool DJANGO, const Bool verbose )
{
     cout << Date( ) << ": begin Star2" << endl;

     // Heuristics.

     const int MIN_STAR = 5000;
     const double MAX_CN_DIFF = 0.5;
     const int MAX_VIEW = 10;
     const int MIN_BAR_TO = 2000;
     const int BC_VIEW = 50000;

     // Place reads.

     cout << Date( ) << ": placing reads" << endl;
     PlaceReads( hb, paths, dup, D, dpaths, True, False );
     cout << Date( ) << ": making index" << endl;
     IntIndex dpaths_index( dpaths, D.E( ) );

     // Find lines.

     cout << Date( ) << ": finding lines" << endl;
     vec<vec<vec<vec<int>>>> dlines;
     FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH, True );

     // Compute kmers.

     vec<int> kmers( hb.E( ) );
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
          kmers[e] = hb.Kmers(e);

     // Compute coverage.

     cout << Date( ) << ": mult before CN" << endl;
     vec<int> mult;
     vec<int> llens;
     // ComputeMult( hb, D, mult );
     GetLineLengths( hb, D, dlines, llens );
     vec<double> cov;
     {    vec<vec<pair<int,int>>> lbp;
          BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, 0 );
          MasterVec<SerfVec<pair<int,int>>> lbpx;
          for ( auto x : lbp )
          {    SerfVec<pair<int,int>> y( x.size( ) );
               for ( int j = 0; j < x.isize( ); j++ )
                    y[j] = x[j];
               lbpx.push_back(y);    }
          LineCN( kmers, lbpx, D, dlines, llens, cov );    }

     // Start passes.

     int sp = 0;
     while(1)
     {    cout << Date( ) << ": start star pass " << ++sp << endl;

          // Introduce barcode-only gaps.  First step: find line neighbors.

          cout << Date( ) << ": get line ancillary data" << endl;
          vec<int> llens, to_left, to_right, linv;
          vec< vec< pair<int,int> > > lhood;
          GetLineLengths( hb, D, dlines, llens );
          LineProx( hb, inv, ebcx, D, dinv, dlines, qept, lhood );
          D.ToLeft(to_left), D.ToRight(to_right);
          LineInv( dlines, dinv, linv );

          // Get barcode positions.

          cout << Date( ) << ": getting barcode positions" << endl;
          vec<vec<pair<int,int>>> lbp;
          BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, BC_VIEW );
     
          // Find long lines that terminate at their right end.

          cout << Date( ) << ": begin star stage" << endl;
          vec<vec<int> > stars;
          vec<int> stars_lens;
          vec<int> l1s;
          for ( int L1 = 0; L1 < dlines.isize( ); L1++ )
          {    if ( llens[L1] < MIN_STAR ) continue;
               int v = to_right[ dlines[L1].back( )[0][0] ];
               if ( v < 0 ) continue;
               if ( D.From(v).nonempty( ) || !D.To(v).solo( ) ) continue;
               l1s.push_back(L1);    }
          vec<int> l1slen( l1s.isize( ) );
          for ( int i = 0; i < l1s.isize( ); i++ )
               l1slen[i] = llens[ l1s[i] ];
     
          // Go through long lines and find their nearest left and right neighbors.

          ReverseSortSync( l1slen, l1s );
          cout << Date( ) << ": start star loop on " << l1s.size( )
               << " L1 values" << endl;
          map< vec<int>, double > memory;
          double clock = WallClockTime( );
          const int batch = Max( 1, l1s.isize( )/500 );
          #pragma omp parallel for schedule(dynamic, batch)
          for ( int li = 0; li < l1s.isize( ); li++ )
          {    int L1 = l1s[li];
               vec< triple<int,int,int> > M;

               // Set up logging.

               ostringstream* outp;
               if (verbose) outp = new ostringstream;
               auto DumpOut = [&]
               {    if (verbose)
                    {
                         #pragma omp critical
                         {    cout << outp->str( );    }    
                         delete outp;    }    };

               // Find left and right neighbors of L1.

               vec<int> lefts, rights;
               for ( int i = 0; i < Min( MAX_VIEW, lhood[L1].isize( ) ); i++ )
               {    int L2 = lhood[L1][i].second;
                    if ( AbsDiff( cov[L1], cov[L2] ) > MAX_CN_DIFF ) continue;
                    if ( llens[L2] < MIN_BAR_TO ) continue;
                    vec<double> scores;
                    scores.push_back( MemoryScoreOrder( 
                         { L2, L1 }, lbp, llens, M, memory ) );
                    scores.push_back( MemoryScoreOrder( 
                         { linv[L2], L1 }, lbp, llens, M, memory ) );
                    scores.push_back( MemoryScoreOrder( 
                         { L1, L2 }, lbp, llens, M, memory ) );
                    scores.push_back( MemoryScoreOrder( 
                         { L1, linv[L2] }, lbp, llens, M, memory ) );
                    vec<int> ids( 4, vec<int>::IDENTITY );
                    SortSync( scores, ids );
                    double ad = scores[1] - scores[0];
                    // if ( L2 != linv[L2] && ad < MIN_ADVANTAGE ) continue;
                    if ( ad < MIN_ADVANTAGE ) continue;
                    if (verbose) (*outp) << "advantage for best right = " << ad << endl;
                    if ( ids[0] == 0 ) lefts.push_back(L2);
                    else if ( ids[0] == 1 ) lefts.push_back( linv[L2] );
                    else if ( ids[0] == 2 ) rights.push_back(L2);
                    else rights.push_back( linv[L2] );    }

               // Find the largest right neighbor L2.

               
               if (verbose)
               {    (*outp) << "\nlooking right from L " << L1 << endl;
                    (*outp) << "see " << printSeq(rights) << endl;    }
               
               if ( rights.empty( ) )
               {    DumpOut( );
                    continue;    }
               int maxlen = -1, L2 = -1;
               for ( auto l : rights )
               {    
                   int w = to_left[ dlines[l].front( )[0][0] ]; 
                   if ( D.To(w).nonempty( ) || !D.From(w).solo( ) )
                   {  
                       continue;    }
                   if ( llens[l] > maxlen )
                   {    maxlen = llens[l];
                        L2 = l;    }    }
               // if (L2 == -1 || maxlen < MIN_STAR) {
               if (L2 == -1) {
                   DumpOut( );
                   continue;
               }
               if (verbose) (*outp) << "largest right = " << L2 << endl;

//====================================================================================
               
               // find lefts and rights for R(L2)

               vec<int> lefts_L2, rights_L2;
               for ( int i = 0; i < Min( MAX_VIEW, lhood[L2].isize( ) ); i++ )
               {    int L1 = lhood[L2][i].second;
                    if ( AbsDiff( cov[L2], cov[L1] ) > MAX_CN_DIFF ) continue;
                    if ( llens[L1] < MIN_BAR_TO ) continue;
                    vec<double> scores;
                    scores.push_back( MemoryScoreOrder( 
                         { L1, L2 }, lbp, llens, M, memory ) );
                    scores.push_back( MemoryScoreOrder( 
                         { linv[L1], L2 }, lbp, llens, M, memory ) );
                    scores.push_back( MemoryScoreOrder( 
                         { L2, L1 }, lbp, llens, M, memory ) );
                    scores.push_back( MemoryScoreOrder( 
                         { L2, linv[L1] }, lbp, llens, M, memory ) );
                    vec<int> ids( 4, vec<int>::IDENTITY );
                    SortSync( scores, ids );
                    double ad = scores[1] - scores[0];
                    if (verbose) (*outp) << "advantage for best left = " << ad << endl;
                    // if ( ad != -1.0 && ad < MIN_ADVANTAGE ) continue;
                    if ( ad < MIN_ADVANTAGE ) continue;
                    if ( ids[0] == 0 ) lefts_L2.push_back(L1);
                    else if ( ids[0] == 1 ) lefts_L2.push_back( linv[L1] );
                    else if ( ids[0] == 2 ) rights_L2.push_back(L1);
                    else rights_L2.push_back( linv[L1] );    }

               if (lefts_L2.empty()) {
                   DumpOut( );
                   continue;
               }
               int L = L1, R = L2;
               // 4. must have L in lefts(R)
               if (Position(lefts_L2, L) == -1) {
                   DumpOut( );
                   continue;
               }
               Bool give_up = False;
               // 5. any member of rights(L) >= 10 kb must be present in R+lefts(R)+rights(R)
               for (int i = rights.isize() - 1; i >= 0; --i) {
                   int l = rights[i];
                   if (llens[l] < 10000) {
                       continue;
                   }
                   if ((l != R) &&
                       (Position(lefts_L2, l) == -1) &&
                       (Position(rights_L2, l) == -1)) {
                       // rights.erase(rights.begin() + i);
                       give_up = True;
                       break;
                   }
               }
               if (give_up) {
                   DumpOut( );
                   continue;
               }
               // 6. any member of lefts(R) >= 10 kb must be present in L+lefts(L)+rights(L)
               for (int i = lefts_L2.isize() - 1; i >= 0; --i) {
                   int l = lefts_L2[i];
                   if (llens[l] < 10000) {
                       continue;
                   }
                   if ((l != L) &&
                       (Position(lefts, l) == -1) &&
                       (Position(rights, l) == -1)) {
                       // lefts_L2.erase(rights.begin() + i);
                       give_up = True;
                       break;
                   }
               }
               if (give_up) {
                   DumpOut( );
                   continue;
               }
               // 7. pool all, drop any lefts of L and rights of R
               vec<int> line_pool;
               for (auto l : rights) {
                   if (l == R || Position(rights_L2, l) != -1) {
                       continue;
                   }
                   line_pool.push_back(l);
               }
               for (auto l : lefts_L2) {
                   if (l == L || Position(lefts, l) != -1) {
                       continue;
                   }
                   line_pool.push_back(l);
               }
               if ( line_pool.empty() ) {
                   DumpOut();
                   continue;
               }
               UniqueSort(line_pool);
               for (int i = line_pool.isize() - 1; i >= 0; --i) {
                   int l = line_pool[i];
                   if ( l > linv[l] && Position(line_pool, linv[l]) != -1) {
                       line_pool.erase(line_pool.begin() + i);
                   }
               }

               // 8. drop tiny stuff until <= 5
               int MAX_BETWEEN = 5;
               if ( line_pool.isize( ) > MAX_BETWEEN )
               {    vec<int> lens( line_pool.size( ) );
                    for ( int i = 0; i < line_pool.isize( ); i++ )
                         lens[i] = llens[ line_pool[i] ];
                    ReverseSortSync( lens, line_pool );
                    line_pool.resize(MAX_BETWEEN);    }
               // 9. O&O, dropping tiny stuff until confident
               vec<int> best_oo;
               /*
               if (line_pool.size() != 1) {
                   DumpOut( );
                   continue;
               } else {
                   best_oo = line_pool;
               }
               */
               
               if (line_pool.size() == 1) {
                   best_oo = line_pool;
               }
               else if ( !DJANGO )
               {    vec<int> best_between;
                    double ad;
                    vec<int> pivot({L, R});
                    if (OO)
                    {    MemoryOrderAndOrientN( pivot, line_pool, lbp, llens, 
                              linv, M, memory, ad, best_between );    }
                    else MemoryOrderN( pivot, line_pool, lbp, llens, M, memory, ad, best_between );
                    if (verbose)
                    {    (*outp) << "\nlooking right from L" << L << endl
                              << "see " << printSeq(line_pool) << endl
                              << "winner = L" << printSeq(best_between)
                              << " [" << ad << "]" << endl;    }
                    if (verbose) (*outp) << "advantage for best bewteen = " << ad << endl;
                    // if ( ad != -1.0 && ad < MIN_ADVANTAGE )
                    if ( ad < MIN_ADVANTAGE )
                    {    DumpOut( );
                         continue;    }
                    best_oo = best_between;    }
               else
               {    vec<int> lens( line_pool.size( ) );
                    for ( int i = 0; i < line_pool.isize( ); i++ )
                         lens[i] = llens[ line_pool[i] ];
                    ReverseSortSync( lens, line_pool );
                    while( line_pool.size( ) > 1 )
                    {    vec<int> best_between;
                         double ad;
                         vec<int> pivot({L, R});
                         if (OO)
                         {    MemoryOrderAndOrientN( pivot, line_pool, lbp, llens, 
                                   linv, M, memory, ad, best_between );    }
                         else 
                         {    MemoryOrderN( pivot, line_pool, lbp, llens, 
                                   M, memory, ad, best_between );    }
                         if (verbose)
                         {    (*outp) << "\nlooking right from L" << L << endl
                                   << "see " << printSeq(line_pool) << endl
                                   << "winner = L" << printSeq(best_between)
                                   << " [" << ad << "]" << endl;    }
                         if (verbose) (*outp) << "advantage for best bewteen = " << ad << endl;
                         // if ( ad == -1.0 || ad >= MIN_ADVANTAGE )
                         if ( ad >= MIN_ADVANTAGE )
                         {    best_oo = best_between;
                              break;    }
                         line_pool.pop_back( );    }    }
               
               // 10. save
               // Check for done and otherwise save.

               if ( best_oo.empty() )
               {    DumpOut( );
                    continue;    }
               give_up = False;
               for (auto l : best_oo) {
                   int v = to_left[ dlines[l].front( )[0][0] ]; 
                   int w = to_right[ dlines[l].back( )[0][0] ]; 
                   if (v < 0 || w < 0) {
                       give_up = True;
                       break;
                   }
                   if ( D.To(v).nonempty( ) || !D.From(v).solo( ) ) {
                       give_up = True;
                       break;
                   }
                   if ( D.From(w).nonempty( ) || !D.To(w).solo( ) ) {
                       give_up = True;
                       break;
                   }
               }
               if (give_up)
               {    DumpOut( );
                    continue;    }
               if (verbose) (*outp) << "saving " << L << " --> " << R << endl;
               DumpOut( );
               vec<int> star({L});
               star.append(best_oo);
               star.push_back(R);
               int len = llens[L];
               for (auto l : best_oo) {
                   len += llens[l];
               }
               len += llens[R];
               #pragma omp critical
               {    
                   stars.push_back(star);    
                   stars_lens.push_back(len);
               }
          }

          // 11. traverse largest by Min(L,R) and don’t join if touched
          // Introduce joins.
          cout << Date( ) << ": done, time used = " << TimeSince(clock) << endl;
          ReverseSortSync(stars_lens, stars);
          vec<Bool> touched_L( dlines.size( ), False );
          vec<Bool> touched_R( dlines.size( ), False );
          vec<pair<int,int>> stars2;
          for ( int i = 0; i < stars.isize( ); i++ ) {
              vec<int> &star = stars[i];
              for (int j = 0; j < star.isize() - 1; ++j) {
                  int L = star[j], R = star[j + 1];
                  if (touched_R[L] || touched_L[R]) {
                      continue;
                  }
                  int d1 = dlines[L].back( )[0][0], d2 = dlines[R].front( )[0][0];
                  int v = to_right[d1], w = to_left[d2];
                  int rv = to_right[ dinv[d2] ], rw = to_left[ dinv[d1] ];
                  if ( D.From(v).nonempty( ) || D.From(rv).nonempty( )
                      || D.To(w).nonempty( ) || D.To(rw).nonempty( ) )
                  {    continue;    }
                  dinv.push_back( D.E( ) + 1, D.E( ) );
                  D.AddEdge( v, w, {-2} );
                  D.AddEdge( rv, rw, {-2} );
                  touched_R[L] = touched_L[R] = True;
                  touched_R[linv[L]] = touched_L[linv[R]] = True;
                  stars2.push( L, R );
                  stars2.push( linv[R], linv[L] );
                  if (verbose)
                  {    cout << "linking from L" << L << " to L" << R << endl;
                      cout << "linking from L" << linv[R] << " to L" << linv[L]
                          << endl;    }
              }
          }
          cout << Date( ) << ": made " << stars2.size( ) << " star joins" << endl;
          if ( stars2.empty( ) ) break;
     
          // Find lines of lines.

          vec<vec<int>> linelines;
          vec<int> ll_inv;
          vec<Bool> touched( dlines.size( ), False );
          {    Sort(stars2);
               vec<pair<int,int>> stars2r( stars2.size( ) );
               for ( int i = 0; i < stars2.isize( ); i++ )
                    stars2r[i] = make_pair( stars2[i].second, stars2[i].first );
               Sort(stars2r);
               for ( int i = 0; i < stars2.isize( ); i++ )
               {    int L = stars2[i].first, R = stars2[i].second;
                    if ( touched[L] || touched[R] ) continue;
                    vec<int> x = {L,R};
                    touched[L] = touched[R] = True;
                    Bool circle = False;
                    while(1)
                    {    int p = BinPosition1( stars2, x.back( ) );
                         if ( p < 0 ) break;
                         int R = stars2[p].second;
                         if ( R == x.front( ) ) 
                         {    circle = True;
                              break;    }
                         x.push_back(R);
                         touched[R] = True;    }
                    if ( !circle )
                    {    while(1)
                         {    int p = BinPosition1( stars2r, x.front( ) );
                              if ( p < 0 ) break;
                              int L = stars2r[p].second;
                              // if ( L == x.back( ) ) 
                              // {    circle = True;
                              //      break;    }
                              x.push_front(L);
                              touched[L] = True;    }    }
                    else {
                        if (verbose) cout << "Fuck, circle detected!!!" << endl;
                    }
                    vec<int> y;
                    for ( int j = x.isize( ) - 1; j >= 0; j-- )
                    {    y.push_back( linv[ x[j] ] );
                         touched[ linv[ x[j] ] ] = True;    }
                    vec<int> xs(x), ys(y);
                    Sort(xs), Sort(ys);
                    if (circle)
                    {    x.push_back( x[0] ), y.push_back( y[0] );    }
                    int nll = linelines.size( );
                    linelines.push_back(x);
                    if ( ys != xs ) 
                    {    linelines.push_back(y);
                         ll_inv.push_back( nll+1, nll );    }
                    else ll_inv.push_back(nll);    }    }

          // Update lines and parallel data structures.

          cout << Date( ) << ": updating lines" << endl;
          vec<vec<int>> gap = { { } };
          for ( int i = 0; i < linelines.isize( ); i++ )
          {    const vec<int>& x = linelines[i];
               vec<vec<vec<int>>> y = dlines[ x[0] ];
               for ( int j = 1; j < x.isize( ); j++ )
               {    y.push_back(gap);
                    if ( x[j] == x[0] ) y.push_back( dlines[ x[j] ][0] );
                    else y.append( dlines[ x[j] ] );    }
               if ( ll_inv[i] == i - 1 )
               {    y = dlines.back( );
                    y.ReverseMe( );
                    for ( int j = 0; j < y.isize( ); j++ )
                    for ( int k = 0; k < y[j].isize( ); k++ )
                    for ( int l = 0; l < y[j][k].isize( ); l++ )
                         y[j][k][l] = dinv[ y[j][k][l] ];    }
               dlines.push_back(y);
               int64_t lensum = 0;
               double cnsum = 0.0;
               for ( int j = 0; j < x.isize( ); j++ )
               {    lensum += llens[ x[j] ];
                    cnsum += llens[ x[j] ] * cov[ x[j] ];    }
               cov.push_back( cnsum/lensum );
               touched.push_back(False);    }
          EraseIf( dlines, touched );
          EraseIf( cov, touched );    }

     // Clean up.

     cout << Date( ) << ": cleaning up at end of Star2" << endl;
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );    }

//==============================================================================================     
void reverse_line_content(const vec<vec<vec<int>>> &line, const vec<int> &dinv, vec<vec<vec<int>>> &reversed_line)
{   
    reversed_line = line;
    reversed_line.ReverseMe( );
    for ( int j = 0; j < reversed_line.isize( ); j++ )
        for ( int k = 0; k < reversed_line[j].isize( ); k++ )
            for ( int l = 0; l < reversed_line[j][k].isize( ); l++ )
                reversed_line[j][k][l] = dinv[ reversed_line[j][k][l] ];    }

void Star3_bk( const HyperBasevectorX& hb, const vec<int>& inv, const vec<Bool>& dup,
     const vec<int32_t>& bc, const ReadPathVecX& paths, 
     digraphE<vec<int>>& D, vec<int>& dinv, 
     ReadPathVec& dpaths, const VecIntVec& ebcx, 
     const vec< triple<int,int,int> >& qept, 
     const double MIN_ADVANTAGE, const Bool OO, 
     const Bool DJANGO, const Bool verbose ) {

     cout << Date( ) << ": begin Star3" << endl;

     // Heuristics.

     const int MIN_STAR = 5000;
     const double MAX_CN_DIFF = 0.5;
     const int MAX_VIEW = 10;
     const int MIN_BAR_TO = 2000;
     const int BC_VIEW = 50000;

     // Place reads.

     cout << Date( ) << ": placing reads" << endl;
     PlaceReads( hb, paths, dup, D, dpaths, True, False );
     cout << Date( ) << ": making index" << endl;
     IntIndex dpaths_index( dpaths, D.E( ) );


     // Compute kmers.

     vec<int> kmers( hb.E( ) );
     #pragma omp parallel for schedule(dynamic)
     for ( int e = 0; e < hb.E( ); e++ )
         kmers[e] = hb.Kmers(e);

     vec<vec<vec<vec<int>>>> dlines;
     int sp = 0;
     while(1)
     {    
         cout << Date( ) << ": start star pass " << ++sp << endl;
         // Introduce barcode-only gaps.  First step: find line neighbors.
         // Find lines.
         cout << Date( ) << ": finding lines" << endl;
         FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH, True );


         vec<int> llens;
         GetLineLengths( hb, D, dlines, llens );

         cout << Date( ) << ": get line ancillary data" << endl;
         vec<int> linv;
         LineInv( dlines, dinv, linv );
         vec<int> to_left, to_right;
         D.ToLeft(to_left), D.ToRight(to_right);
         vec<double> cov;
         {    
             vec<vec<pair<int,int>>> lbp;
             BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, 0 );
             MasterVec<SerfVec<pair<int,int>>> lbpx;
             for ( auto x : lbp )
             {    
                 SerfVec<pair<int,int>> y( x.size( ) );
                 for ( int j = 0; j < x.isize( ); j++ )
                     y[j] = x[j];
                 lbpx.push_back(y);
             }
             LineCN( kmers, lbpx, D, dlines, llens, cov );    
         }
         // Get barcode positions.

         vec<vec<pair<int,int>>> lbp;
         BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, BC_VIEW );

         // Find (long) lines with bc-only gaps.

         cout << Date( ) << ": begin star stage" << endl;
         vec<vec<int> > stars;
         vec<int> stars_lens;
         vec<int> l1s;
         vec<vec<int> > break_cell_ids(dlines.isize());
         for ( int L1 = 0; L1 < dlines.isize( ); L1++ )
         {    
             // if ( llens[L1] < MIN_STAR ) continue;
             const vec<vec<vec<int>>>& L = dlines[L1];
             for (int j = 0; j < L.isize(); ++j) {
                 const vec<vec<int>> &cell = L[j];
                 if (cell.solo() && cell[0].empty()) {
                     int d = L[j - 1][0][0];
                     int e = D.IFrom(to_right[d], 0);
                     if (IsBarcodeOnlyGap( D.O(e))) {
                         break_cell_ids[L1].push_back(j);
                     }
                 }
             }
             if (!break_cell_ids[L1].empty()) {
                 l1s.push_back(L1);
             }
         }
         vec<int> l1slen( l1s.isize( ) );
         for ( int i = 0; i < l1s.isize( ); i++ )
             l1slen[i] = llens[ l1s[i] ];

         // Go through long lines and find their nearest left and right neighbors.
         vec<vec<triple<double, int, int> > > line_scores(dlines.size());
         ReverseSortSync( l1slen, l1s );
         cout << Date( ) << ": start star loop on " << l1s.size( )
              << " L1 values" << endl;
         map< vec<int>, double > memory;


         for ( int li = 0; li < l1s.isize( ); li++ ) {

             int lid = l1s[li];
             vec<int> &break_cell_id = break_cell_ids[lid];
             vec<vec<vec<int>>> line = dlines[lid];
             vec<vec<vec<vec<int>>>> broken_lines;
             for (int i = 0; i < break_cell_id.isize(); ++i) {
                 int cell_id = break_cell_id[i];
                 if (i == 0) {
                     const vec<vec<vec<int>>> broken_line(line.begin(), line.begin() + cell_id);
                     broken_lines.push_back(broken_line);
                 } else {
                     int cell_id_prev = break_cell_id[i - 1];
                     ForceAssertGt(cell_id, cell_id_prev);
                     const vec<vec<vec<int>>> broken_line(line.begin() + cell_id_prev + 1, line.begin() + cell_id);
                     broken_lines.push_back(broken_line);
                 }
             }
             // add last broken line.
             const vec<vec<vec<int> > > broken_line(line.begin() + break_cell_id.back() + 1, line.end());
             broken_lines.push_back(broken_line);

             // create broken lines.
             vec<vec<vec<vec<int>>>> dlines_extra;
             for (int broken_lid = 0; broken_lid < broken_lines.isize(); ++ broken_lid) {
                 vec<vec<vec<int>>> &L = broken_lines[broken_lid];
                 vec<vec<vec<int>>> reversed_line;
                 dlines_extra.push_back(L);
                 reverse_line_content(L, dinv, reversed_line);
                 dlines_extra.push_back(reversed_line);
             }
             vec<int> llens_extra;
             GetLineLengths( hb, D, dlines_extra, llens_extra );

             vec<double> cov_extra;
             {    
                 vec<vec<pair<int,int>>> lbp;
                 BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines_extra, lbp, 0 );
                 MasterVec<SerfVec<pair<int,int>>> lbpx;
                 for ( auto x : lbp )
                 {    
                     SerfVec<pair<int,int>> y( x.size( ) );
                     for ( int j = 0; j < x.isize( ); j++ )
                         y[j] = x[j];
                     lbpx.push_back(y);
                 }
                 LineCN( kmers, lbpx, D, dlines_extra, llens_extra, cov_extra );    
             }
             vec<vec<pair<int,int>>> lbp_extra;
             BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines_extra, lbp_extra, BC_VIEW );

             #pragma omp parallel for schedule(dynamic)
             for (int broken_cell_id = 0; broken_cell_id < break_cell_id.isize(); ++broken_cell_id) {
                 vec< triple<int,int,int> > M;
                 // Set up logging.

                 ostringstream* outp;
                 if (verbose) outp = new ostringstream;
                 auto DumpOut = [&]
                 {    if (verbose)
                     {
                         #pragma omp critical
                         {    cout << outp->str( );    }    
                         delete outp;    }    };
                 // make a copy of lines
                 // remember the difference between lines before current line
                 // and after it.
                 vec<int> new_to_old_lid(dlines.size() + 4, -1);
                 vec<vec<vec<vec<int>>>> dlines_copy;
                 vec<int> llens_copy;
                 vec<double> cov_copy;
                 vec<vec<pair<int,int>>> lbp_copy;
                 for (int j = 0; j < dlines.isize(); ++j) {
                     if (j == lid || j == linv[lid]) {
                         continue;
                     }
                     new_to_old_lid[dlines_copy.size()] = j;
                     dlines_copy.push_back(dlines[j]);
                     llens_copy.push_back(llens[j]);
                     cov_copy.push_back(cov[j]);
                     lbp_copy.push_back(lbp[j]);
                 }
                 // assign ids for lhs and rhs
                 int L = dlines_copy.isize();
                 dlines_copy.push_back(dlines_extra[2*broken_cell_id]);
                 dlines_copy.push_back(dlines_extra[2*broken_cell_id + 1]);
                 int R = dlines_copy.isize();
                 dlines_copy.push_back(dlines_extra[2*broken_cell_id + 2]);
                 dlines_copy.push_back(dlines_extra[2*broken_cell_id + 3]);

                 llens_copy.push_back(llens_extra[2*broken_cell_id],
                                      llens_extra[2*broken_cell_id + 1],
                                      llens_extra[2*broken_cell_id + 2],
                                      llens_extra[2*broken_cell_id + 3]);

                 cov_copy.push_back(cov_extra[2*broken_cell_id],
                                      cov_extra[2*broken_cell_id + 1],
                                      cov_extra[2*broken_cell_id + 2],
                                      cov_extra[2*broken_cell_id + 3]);
                 // Get barcode positions.
                 lbp_copy.push_back(lbp_extra[2*broken_cell_id],
                                      lbp_extra[2*broken_cell_id + 1],
                                      lbp_extra[2*broken_cell_id + 2],
                                      lbp_extra[2*broken_cell_id + 3]);


                 vec<int> linv_copy;
                 LineInv( dlines_copy, dinv, linv_copy );

                 vec< vec< pair<int,int> > > lhood;
                 LineProx( hb, inv, ebcx, D, dinv, dlines_copy, qept, lhood );


                 // find intersection of L's and R's lhood.
                 vec<int> line_pool, lefts, rights;

                 for (int i = 0; i < lhood[L].isize(); ++i) {
                     rights.push_back(lhood[L][i].second);
                 }
                 for (int i = 0; i < lhood[R].isize(); ++i) {
                     lefts.push_back(lhood[R][i].second);
                 }
                 for (int i = 0; i < lefts.isize(); ++i) {
                     int L1 = lefts[i];
                     if (L1 == R || Position(rights, L1) == -1) {
                         continue;
                     }
                     line_pool.push_back(L1);
                 }

                 if ( line_pool.empty() )
                 {    DumpOut( );
                     continue;    }

                 vec<double> scores;
                 for (int i = 0; i < line_pool.isize(); ++i) {
                     // LB means Line Between
                     int LB = line_pool[i];
                     if ( llens_copy[LB] < MIN_BAR_TO ) continue;
                     if ( AbsDiff( cov_copy[L], cov_copy[LB] ) > MAX_CN_DIFF ) continue;
                     if ( AbsDiff( cov_copy[R], cov_copy[LB] ) > MAX_CN_DIFF ) continue;
                     scores.push_back( MemoryScoreOrder( 
                                 { L, LB , R }, lbp_copy, llens_copy, M, memory ) );
                     scores.push_back( MemoryScoreOrder( 
                                 { L, linv_copy[LB] , R }, lbp_copy, llens_copy, M, memory ) );
                 }
                 if (scores.empty()) {
                     DumpOut( );
                     continue;    
                 }
                 vec<int> ids( scores.size(), vec<int>::IDENTITY );
                 SortSync( scores, ids );
                 double ad = scores[1] - scores[0];
                 if (ad == 0 && (line_pool[ids[0]/2] == linv_copy[line_pool[ids[1]/2]])) {
                     if (verbose) (*outp) << "REALLY?" << endl;
                 } else if ( ad < MIN_ADVANTAGE ) {
                     DumpOut();
                     continue;
                 }
                 int best_oo = line_pool[ids[0] / 2];
                 if (ids[0] % 2 != 0) {
                     best_oo = linv_copy[best_oo];
                 }
                 int v = to_left[ dlines_copy[best_oo].front( )[0][0] ]; 
                 int w = to_right[ dlines_copy[best_oo].back( )[0][0] ]; 
                 if (v < 0 || w < 0) {
                     DumpOut();
                     continue;
                 }
                 if ( D.To(v).nonempty( ) || !D.From(v).solo( ) ) {
                     DumpOut();
                     continue;
                 }
                 if ( D.From(w).nonempty( ) || !D.To(w).solo( ) ) {
                     DumpOut();
                     continue;
                 }
                 // transfer new lid to old one
                 best_oo = new_to_old_lid[best_oo];
                 // shouldn't happen
                 if (best_oo == -1) {
                     if (verbose) (*outp) << "can not find proper old lid"  << endl;
                     DumpOut( );
                     continue;
                 }

                 if (verbose) {
                     (*outp) << "line " << lid << " gap cell "
                             << break_cell_id[broken_cell_id]
                             << " can fill line " << best_oo << endl;;
                 }
                 DumpOut( );
                 #pragma omp critical
                 {
                     // <score, lid, cell_id>
                     line_scores[best_oo].push(scores[0], lid, break_cell_id[broken_cell_id]);
                 }
             }
         }
         // introduce joins.
         vec<triple<double, int, int> > best_line_scores(dlines.size(), triple<double, int, int> (-1, -1, -1));
         for (int i = 0; i < line_scores.isize(); ++i) {
             vec<triple<double, int, int> > &line_score = line_scores[i];
             if (verbose && line_score.nonempty()) {
                 cout << "scores for line " << i << " is (are): "
                      << printSeq(line_score) << endl;
             }
             if (line_score.size() > 1) {
                 Sort(line_score);
                 best_line_scores[i] = line_score[0];
             } else if (line_score.solo()) {
                 best_line_scores[i] = line_score[0];
             }
         }
         // find non-circle (dead lock) line relationship and insert bc-only
         // gaps, also break down the old ones.
         vec<int> lines_to_fill;
         for (int i = 0; i < best_line_scores.isize(); ++i) {
             // lines with no place to fill at all.
             if (best_line_scores[i].second == -1 ||
                 best_line_scores[i].third == -1) {
                 continue;
             }
             // detect if there is a circle.
             int start = best_line_scores[i].second;
             int next = start;
             Bool is_circle = False;
             vec<int> path({start});
             while (1) {
                 next = best_line_scores[next].second;
                 path.push_back(next);
                 if (next == -1) {
                     break;
                 } else if (next == start) {
                     is_circle = True;
                     if (verbose) {
                         cout << "circle detected: " 
                             << printSeq(path) << endl;
                     }
                     break;
                 }
             }
             if (is_circle) {
                 continue;
             }
             // if reach here, there is no dead lock.
             lines_to_fill.push_back(i);
         }
         
         if (verbose) {
             cout << "there are " << lines_to_fill.size()
                  << " lines to fill (candidates)" << endl;
         }
         if (lines_to_fill.empty()) {
             break;
         }

         // check if there is any collision with lines to fill.
         // first, create a matrix of line-cell-line_to_fill
         vec<map<int, int> > line_cell(dlines.size());
         // vec<vec<int> > break_cell_ids(dlines.isize());
         for (int i = 0; i < best_line_scores.isize(); ++i) {
             triple<double, int, int> &best_line_score = best_line_scores[i];
             int lid = best_line_score.second;
             int cell_id = best_line_score.third;
             if (lid == -1 || cell_id == -1) {
                 continue;
             }
             line_cell[lid][cell_id] = i;
         }

         Bool give_up = True;
         vec<int> dels;
         for (int i = 0; i < lines_to_fill.isize(); ++i) {
             int lid_to_fill = lines_to_fill[i];
             triple<double, int, int> &best_line_score = best_line_scores[lid_to_fill];
             int lid = best_line_score.second;
             int cell_id = best_line_score.third;
             // gap cell id on reverse strand line
             int cell_id_re = dlines[lid].isize() - cell_id - 1;
             if (line_cell[linv[lid]].find(cell_id_re) == line_cell[linv[lid]].end()) {
                 if (verbose) {
                     cout << "warning type 1! line " << lid_to_fill
                         << " only fill into one strand of line "
                         << lid << endl;
                 }
             } else if (line_cell[linv[lid]][cell_id_re] != linv[lid_to_fill]) {
                 if (verbose) {
                     cout << "warning type 2! line " << lid
                         << " have different lines to fill in on two strands: "
                         << lid_to_fill << " and "
                         << line_cell[linv[lid]][cell_id_re] << endl;
                 }
                 if (llens[lid_to_fill] <= llens[line_cell[linv[lid]][cell_id_re]]) {
                     continue;
                 }
             }
             // another check
             int lid_to_fill_re = linv[lid_to_fill];
             triple<double, int, int> &best_line_score_re = best_line_scores[lid_to_fill_re];
             if (lid_to_fill_re != lid_to_fill &&
                 best_line_score_re.second != -1 &&
                 best_line_score_re.third != -1) {

                 // collision can't be solved.
                 if (best_line_score_re.second != linv[lid] ||
                     best_line_score_re.third != cell_id_re) {
                     cout << "warning type 3! line" << lid_to_fill << " and "
                         <<  lid_to_fill_re << " have different gap cell to fill"
                         << endl;
                     continue;
                 }
             }

             give_up = False;
             // break original bc-only gaps.
             int d1 = dlines[lid][cell_id - 1][0][0];
             int d2 = dlines[lid][cell_id + 1][0][0];
             int v = to_right[d1], w = to_left[d2];
             int rv = to_right[ dinv[d2] ], rw = to_left[ dinv[d1] ];
             ForceAssert((D.IFrom(v).solo() &&
                          D.ITo(w).solo() &&
                          D.IFrom(rv).solo() &&
                          D.ITo(rw).solo()));
             int bc_gap_d = D.IFrom(v, 0);
             ForceAssert(IsBarcodeOnlyGap( D.O(bc_gap_d)));
             dels.push_back(bc_gap_d);
             dels.push_back(dinv[bc_gap_d]);

             // insert new bc-only gaps.
             // get v1, w1 for line_to_fill
             int d_left = dlines[lid_to_fill].front()[0][0];
             int d_right = dlines[lid_to_fill].back()[0][0];
             int v1 = to_left[d_left], w1 = to_right[d_right];
             int rv1 = to_left[ dinv[d_right] ], rw1 = to_right[ dinv[d_left] ];
             ForceAssert((D.ITo(v1).empty() &&
                          D.IFrom(w1).empty() &&
                          D.ITo(rv1).empty() &&
                          D.IFrom(rw1).empty()));

             // add 4 new superedges of bc-only gap.
             dinv.push_back( D.E( ) + 1, D.E( ), D.E() + 3, D.E() + 2 );
             D.AddEdge( v, v1, {-2} );
             D.AddEdge( rw1, rw, {-2} );
             D.AddEdge( w1, w, {-2} );
             D.AddEdge( rv, rv1, {-2} );
         }
         if (give_up) {
             ForceAssert(dels.empty());
             if (verbose) {
                 cout << "no lines can  really be filled. Terminate." << endl;
             }
             break;
         }
         // clean up current supergraph D.
         D.DeleteEdges(dels);
         RemoveUnneededVertices( D, dinv );
         CleanupCore( D, dinv );
     }
}

void LineProxSP( const HyperBasevectorX& hb, const vec<int>& inv,
     const VecIntVec& ebcx, const digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines, const vec<vec<vec<vec<int>>>>& dlines_broken,
     const vec<int> &broken_to_whole,
     const vec< triple<int,int,int> >& qept, vec< vec< pair<int,int> > >& lhood )
{
     // Heuristics.

     const int LOOK_IN = 10000;

     // Stuff in parallel.

     vec<int> mult, tod( hb.E( ), -1 ), tol( D.E( ), -1 ), linv;
     #pragma omp parallel sections
     {
          // Index edges.

          #pragma omp section
          {    ComputeMult( hb, D, mult );
               for ( int d = 0; d < D.E( ); d++ )
               {    if ( D.O(d)[0] < 0 ) continue;
                    for ( int i = 0; i < D.O(d).isize( ); i++ )
                         tod[ D.O(d)[i] ] = d;    }    }

          // Index lines.

          #pragma omp section
          {    for ( int i = 0; i < dlines.isize( ); i++ )
               for ( int j = 0; j < dlines[i].isize( ); j++ )
               for ( int k = 0; k < dlines[i][j].isize( ); k++ )
               for ( int l = 0; l < dlines[i][j][k].isize( ); l++ )
               {    int e = dlines[i][j][k][l];
                    if ( e >= 0 ) tol[e] = i;    }    }

          // Get line inversion.

          #pragma omp section
          {    LineInv( dlines, dinv, linv );    }    }

     // Compute barcode sets for each line, looking only at the ends.

     cout << Date( ) << ": compute barcode sets for lines" << endl;
     vec<vec<int>> lbc( dlines.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    const vec<vec<vec<int>>>& L = dlines[i];
          int pos = 0;
          for ( int j = 0; j < L.isize( ); j++ )
          {    vec<int> lensj;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    int len = 0;
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                    {    int d = L[j][k][l];
                         if ( D.O(d)[0] < 0 ) continue;
                         for ( int m = 0; m < D.O(d).isize( ); m++ )
                         {    int e = D.O(d)[m];
                              if ( pos + len <= LOOK_IN && mult[e] == 1 )
                              {    for ( auto b : ebcx[e] ) lbc[i].push_back(b);    }
                              len += hb.Kmers(e);    }    }
                    lensj.push_back(len);    }
               Sort(lensj);
               if ( lensj.nonempty( ) ) pos += Median(lensj);    }
          pos = 0;
          for ( int j = L.isize( ) - 1; j >= 0; j-- )
          {    vec<int> lensj;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    int len = 0;
                    for ( int l = L[j][k].isize( ) - 1; l >= 0; l-- )
                    {    int d = L[j][k][l];
                         if ( D.O(d)[0] < 0 ) continue;
                         for ( int m = D.O(d).isize( ) - 1; m >= 0; m-- )
                         {    int e = D.O(d)[m];
                              if ( pos + len <= LOOK_IN && mult[e] == 1 )
                              {    for ( auto b : ebcx[e] ) lbc[i].push_back(b);    }
                              len += hb.Kmers(e);    }    }
                    lensj.push_back(len);    }
               Sort(lensj);
               if ( lensj.nonempty( ) ) pos += Median(lensj);    }
          UniqueSort( lbc[i] );    }

     cout << Date( ) << ": compute barcode sets for broken lines" << endl;
     vec<vec<int>> lbc_broken( dlines_broken.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < dlines_broken.isize( ); i++ )
     {    const vec<vec<vec<int>>>& L = dlines_broken[i];
          int pos = 0;
          for ( int j = 0; j < L.isize( ); j++ )
          {    vec<int> lensj;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    int len = 0;
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                    {    int d = L[j][k][l];
                         if ( D.O(d)[0] < 0 ) continue;
                         for ( int m = 0; m < D.O(d).isize( ); m++ )
                         {    int e = D.O(d)[m];
                              if ( pos + len <= LOOK_IN && mult[e] == 1 )
                              {    for ( auto b : ebcx[e] ) lbc_broken[i].push_back(b);    }
                              len += hb.Kmers(e);    }    }
                    lensj.push_back(len);    }
               Sort(lensj);
               if ( lensj.nonempty( ) ) pos += Median(lensj);    }
          pos = 0;
          for ( int j = L.isize( ) - 1; j >= 0; j-- )
          {    vec<int> lensj;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    int len = 0;
                    for ( int l = L[j][k].isize( ) - 1; l >= 0; l-- )
                    {    int d = L[j][k][l];
                         if ( D.O(d)[0] < 0 ) continue;
                         for ( int m = D.O(d).isize( ) - 1; m >= 0; m-- )
                         {    int e = D.O(d)[m];
                              if ( pos + len <= LOOK_IN && mult[e] == 1 )
                              {    for ( auto b : ebcx[e] ) lbc_broken[i].push_back(b);    }
                              len += hb.Kmers(e);    }    }
                    lensj.push_back(len);    }
               Sort(lensj);
               if ( lensj.nonempty( ) ) pos += Median(lensj);    }
          UniqueSort( lbc_broken[i] );    }
     // Index barcode links.

     cout << Date( ) << ": indexing barcode links" << endl;
     vec<int> qep_index( hb.E( ) + 1, -1 );
     qep_index.back( ) = qept.size( );
     for ( int64_t i = qept.jsize( ) - 1; i >= 0; i-- )
          qep_index[ qept[i].first ] = i;
     for ( int e = hb.E( ) - 1; e >= 0; e-- )
          if ( qep_index[e] < 0 ) qep_index[e] = qep_index[e+1];

     // Find barcode neighbors.

     cout << Date( ) << ": finding barcode neighbors" << endl;
     lhood.clear_and_resize( dlines_broken.size( ) );
     const int MIN_LEN = 100;
     const int MIN_LINKS = 6;
     const double MIN_NHOOD_FRAC = 0.1;
     #pragma omp parallel for schedule (dynamic, 1)
     for ( int i1 = 0; i1 < dlines_broken.isize( ); i1++ )
     {    // if ( linv[i1] < i1 ) continue;
          vec<int> n;
          const vec<vec<vec<int>>>& L = dlines_broken[i1];
          for ( int j = 0; j < L.isize( ); j++ )
          for ( int k = 0; k < L[j].isize( ); k++ )
          for ( int l = 0; l < L[j][k].isize( ); l++ )
          {    int d = L[j][k][l];
               if ( D.O(d)[0] < 0 ) continue;
               for ( int m = 0; m < D.O(d).isize( ); m++ )
               {    int e = D.O(d)[m];
                    if ( mult[e] != 1 ) continue;
                    e = Min( e, inv[e] );
                    for ( int64_t l = qep_index[e]; l < qep_index[e+1]; l++ )
                    {    int f = qept[l].second;
                         if ( hb.Kmers(f) < MIN_LEN || mult[f] != 1 ) continue;
                         int i2a = tol[tod[f]];
                         if ( i2a < 0 ) continue;
                         int i2 = Min( i2a, linv[i2a] );
                         if (i2 != broken_to_whole[i1] &&
                             i2 != linv[broken_to_whole[i1]] ) {
                             n.push_back(i2);
                         }
                    }
               }
          }
          UniqueSort(n);
          for ( int j = 0; j < n.isize( ); j++ )
          {    int c = MeetSize( lbc_broken[i1], lbc[ n[j] ] );
               if ( c >= MIN_LINKS ) lhood[i1].push( c, n[j] );    }
          ReverseSort( lhood[i1] );
          for ( int j = 1; j < lhood[i1].isize( ); j++ )
          {    if ( double( lhood[i1][j].first ) / lhood[i1][0].first
                    < MIN_NHOOD_FRAC )
               {    lhood[i1].resize(j);
                    break;    }    }
          /*if ( linv[i1] > i1 ) lhood[ linv[i1] ] = lhood[i1];*/    }    }

void Star3( const HyperBasevectorX& hb, const vec<int>& inv, const vec<Bool>& dup,
     const vec<int32_t>& bc, const ReadPathVecX& paths, 
     digraphE<vec<int>>& D, vec<int>& dinv, 
     ReadPathVec& dpaths, const VecIntVec& ebcx, 
     const vec< triple<int,int,int> >& qept, 
     const double MIN_ADVANTAGE, const Bool OO, 
     const Bool DJANGO, const Bool verbose ) {

     cout << Date( ) << ": begin Star3" << endl;

     // Heuristics.

     const int MIN_STAR = 5000;
     const double MAX_CN_DIFF = 0.5;
     const int MAX_VIEW = 10;
     const int MIN_BAR_TO = 2000;
     const int BC_VIEW = 50000;

     // Place reads.

     cout << Date( ) << ": placing reads" << endl;
     PlaceReads( hb, paths, dup, D, dpaths, True, False );
     cout << Date( ) << ": making index" << endl;
     IntIndex dpaths_index( dpaths, D.E( ) );


     // Compute kmers.

     vec<int> kmers( hb.E( ) );
     #pragma omp parallel for schedule(dynamic)
     for ( int e = 0; e < hb.E( ); e++ )
         kmers[e] = hb.Kmers(e);

     vec<vec<vec<vec<int>>>> dlines;
     int sp = 0;
     while(1)
     {    
         cout << Date( ) << ": start star pass " << ++sp << endl;
         // Introduce barcode-only gaps.  First step: find line neighbors.
         // Find lines.
         cout << Date( ) << ": finding lines" << endl;
         FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH, True );


         vec<int> llens;
         GetLineLengths( hb, D, dlines, llens );

         cout << Date( ) << ": get line ancillary data" << endl;
         vec<int> linv;
         LineInv( dlines, dinv, linv );
         vec<int> to_left, to_right;
         D.ToLeft(to_left), D.ToRight(to_right);

         // Find (long) lines with bc-only gaps.

         cout << Date( ) << ": begin star stage" << endl;
         vec<vec<int> > stars;
         vec<int> stars_lens;
         vec<int> l1s;
         vec<vec<int> > break_cell_ids(dlines.size());
         for ( int L1 = 0; L1 < dlines.isize( ); L1++ )
         {    
             // if ( llens[L1] < MIN_STAR ) continue;
             const vec<vec<vec<int>>>& L = dlines[L1];
             for (int j = 0; j < L.isize(); ++j) {
                 const vec<vec<int>> &cell = L[j];
                 if (cell.solo() && cell[0].empty()) {
                     int d = L[j - 1][0][0];
                     int e = D.IFrom(to_right[d], 0);
                     if (IsBarcodeOnlyGap( D.O(e))) {
                         break_cell_ids[L1].push_back(j);
                     }
                 }
             }
             if (!break_cell_ids[L1].empty()) {
                 l1s.push_back(L1);
             }
         }
         if (l1s.empty()) {
             if (verbose) {
                 cout << "there is no line with bc-only barcode. Terminate."
                      << endl;
             }
             break;
         }
         vec<int> l1slen( l1s.isize( ) );
         for ( int i = 0; i < l1s.isize( ); i++ )
             l1slen[i] = llens[ l1s[i] ];

         // Go through long lines and find their nearest left and right neighbors.
         vec<vec<triple<double, int, int> > > line_scores(dlines.size());
         ReverseSortSync( l1slen, l1s );
         cout << Date( ) << ": start star loop on " << l1s.size( )
              << " L1 values" << endl;

         vec<int> dlines_broken_start_pos;
         vec<vec<vec<vec<int>>>> dlines_broken;
         vec<int> broken_to_whole;
         for ( int li = 0; li < l1s.isize( ); li++ ) {
             dlines_broken_start_pos.push_back(dlines_broken.size());
             int lid = l1s[li];
             vec<int> &break_cell_id = break_cell_ids[lid];
             vec<vec<vec<int>>> line = dlines[lid];
             for (int i = 0; i < break_cell_id.isize(); ++i) {
                 int cell_id = break_cell_id[i];
                 if (i == 0) {
                     const vec<vec<vec<int>>> broken_line(line.begin(), line.begin() + cell_id);
                     dlines_broken.push_back(broken_line);
                     broken_to_whole.push_back(lid);
                 } else {
                     int cell_id_prev = break_cell_id[i - 1];
                     ForceAssertGt(cell_id, cell_id_prev);
                     const vec<vec<vec<int>>> broken_line(line.begin() + cell_id_prev + 1, line.begin() + cell_id);
                     dlines_broken.push_back(broken_line);
                     broken_to_whole.push_back(lid);
                 }
             }
             // add last broken line.
             const vec<vec<vec<int> > > broken_line(line.begin() + break_cell_id.back() + 1, line.end());
             dlines_broken.push_back(broken_line);
             broken_to_whole.push_back(lid);
         }
         // add broken lines size to star_pos
         dlines_broken_start_pos.push_back(dlines_broken.size());

         // compute cov for original dlines
         vec<double> cov;
         {    
             vec<vec<pair<int,int>>> lbp;
             BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, 0 );
             MasterVec<SerfVec<pair<int,int>>> lbpx;
             for ( auto x : lbp )
             {    
                 SerfVec<pair<int,int>> y( x.size( ) );
                 for ( int j = 0; j < x.isize( ); j++ )
                     y[j] = x[j];
                 lbpx.push_back(y);
             }
             LineCN( kmers, lbpx, D, dlines, llens, cov );    
         }
         // Get barcode positions for original lines
         vec<vec<pair<int,int>>> lbp;
         BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, BC_VIEW );

         // compute cov for broken dlines
         vec<int> llens_broken;
         GetLineLengths( hb, D, dlines_broken, llens_broken );
         vec<double> cov_broken;
         {    
             vec<vec<pair<int,int>>> lbp_broken;
             BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines_broken, lbp_broken, 0 );
             MasterVec<SerfVec<pair<int,int>>> lbpx;
             for ( auto x : lbp_broken )
             {    
                 SerfVec<pair<int,int>> y( x.size( ) );
                 for ( int j = 0; j < x.isize( ); j++ )
                     y[j] = x[j];
                 lbpx.push_back(y);
             }
             LineCN( kmers, lbpx, D, dlines_broken, llens_broken, cov_broken );    
         }
         // Get barcode positions for broken lines
         vec<vec<pair<int,int>>> lbp_broken;
         BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines_broken, lbp_broken, BC_VIEW );

         // compute barcode neighbor
         vec< vec< pair<int,int> > > lhood;
         LineProxSP( hb, inv, ebcx, D, dinv, dlines, dlines_broken, broken_to_whole, qept, lhood );

         for ( int li = 0; li < l1s.isize( ); li++ ) {
             int start = dlines_broken_start_pos[li];
             int end = dlines_broken_start_pos[li + 1];
             int lid = l1s[li];
             vec<int> &break_cell_id = break_cell_ids[lid];


             ForceAssertEq(break_cell_id.isize(), end - start - 1);
             #pragma omp parallel for schedule(dynamic)
             for (int broken_cell_id = 0; broken_cell_id < break_cell_id.isize(); ++broken_cell_id) {
                 vec< triple<int,int,int> > M;
                 map< vec<int>, double > memory;
                 // Set up logging.

                 ostringstream* outp;
                 if (verbose) outp = new ostringstream;
                 auto DumpOut = [&]
                 {    if (verbose)
                     {
                         #pragma omp critical
                         {    cout << outp->str( );    }    
                         delete outp;    }    };
                 // make a copy of lines
                 // remember the difference between lines before current line
                 // and after it.
                 // vec<vec<vec<vec<int>>>> dlines_copy(dlines);
                 vec<int> llens_copy(llens);
                 // vec<double> cov_copy(cov);
                 vec<vec<pair<int,int>>> lbp_copy(lbp);

                 // assign ids for lhs and rhs
                 int llid = start + broken_cell_id;
                 int rlid = llid + 1;
                 int L = dlines.isize();
                 // dlines_copy.push_back(dlines_broken[llid]);
                 int R = L + 1;
                 // dlines_copy.push_back(dlines_broken[rlid]);

                 llens_copy.push_back(llens_broken[llid],
                                      llens_broken[rlid]);

                 /*
                 cov_copy.push_back(cov_broken[llid],
                                    cov_broken[rlid]);
                                    */
                 // Get barcode positions.
                 lbp_copy.push_back(lbp_broken[llid],
                                    lbp_broken[rlid]);


                 // vec<int> linv_copy;
                 // LineInv( dlines_copy, dinv, linv_copy );



                 // find intersection of L's and R's lhood.
                 vec<int> line_pool, lefts, rights;

                 for (int i = 0; i < lhood[llid].isize(); ++i) {
                     rights.push_back(lhood[llid][i].second);
                 }
                 for (int i = 0; i < lhood[rlid].isize(); ++i) {
                     lefts.push_back(lhood[rlid][i].second);
                 }
                 for (int i = 0; i < lefts.isize(); ++i) {
                     int L1 = lefts[i];
                     if (Position(rights, L1) == -1 || L1 == linv[L1]) {
                         continue;
                     }
                     line_pool.push_back(L1);
                 }

                 if ( line_pool.empty() )
                 {    DumpOut();
                     continue;    }

                 vec<double> scores;
                 for (int i = 0; i < line_pool.isize(); ++i) {
                     // LB means Line Between
                     int LB = line_pool[i];
                     if ( llens_copy[LB] < MIN_BAR_TO ) continue;
                     if ( AbsDiff( cov_broken[llid], cov[LB] ) > MAX_CN_DIFF ) continue;
                     if ( AbsDiff( cov_broken[rlid], cov[LB] ) > MAX_CN_DIFF ) continue;
                     scores.push_back( MemoryScoreOrder( 
                                 { L, LB , R }, lbp_copy, llens_copy, M, memory ) );
                     scores.push_back( MemoryScoreOrder( 
                                 { L, linv[LB] , R }, lbp_copy, llens_copy, M, memory ) );
                 }
                 if (scores.empty()) {
                     DumpOut();
                     continue;    
                 }
                 vec<int> ids( scores.size(), vec<int>::IDENTITY );
                 SortSync( scores, ids );
                 double ad = scores[1] - scores[0];
                 if (ad == 0 && (line_pool[ids[0]/2] == linv[line_pool[ids[1]/2]])) {
                     if (verbose) (*outp) << "REALLY?" << endl;
                     DumpOut();
                     continue;
                 } else if ( ad < MIN_ADVANTAGE ) {
                     DumpOut();
                     continue;
                 }
                 int best_oo = line_pool[ids[0] / 2];
                 if (ids[0] % 2 != 0) {
                     best_oo = linv[best_oo];
                 }
                 ForceAssertLt(best_oo, dlines.isize());
                 int v = to_left[ dlines[best_oo].front( )[0][0] ]; 
                 int w = to_right[ dlines[best_oo].back( )[0][0] ]; 
                 if (v < 0 || w < 0) {
                     DumpOut();
                     continue;
                 }
                 if ( D.To(v).nonempty( ) || !D.From(v).solo( ) ) {
                     DumpOut();
                     continue;
                 }
                 if ( D.From(w).nonempty( ) || !D.To(w).solo( ) ) {
                     DumpOut();
                     continue;
                 }

                 // shouldn't happen
                 /*
                 if (best_oo == -1) {
                     if (verbose) (*outp) << "can not find proper old lid"  << endl;
                     DumpOut( );
                     continue;
                 }
                 */

                 if (verbose) {
                     (*outp) << "line " << lid << " gap cell "
                             << break_cell_id[broken_cell_id]
                             << " can be filled with line " << best_oo << endl;;
                 }
                 DumpOut();
                 #pragma omp critical
                 {
                     // <score, lid, cell_id>
                     line_scores[best_oo].push(scores[0], lid, break_cell_id[broken_cell_id]);
                 }
             }
         }
         // introduce joins.
         vec<triple<double, int, int> > best_line_scores(dlines.size(), triple<double, int, int> (-1, -1, -1));
         for (int i = 0; i < line_scores.isize(); ++i) {
             vec<triple<double, int, int> > &line_score = line_scores[i];
             if (verbose && line_score.nonempty()) {
                 cout << "scores for line " << i << " is (are): "
                      << printSeq(line_score) << endl;
             }
             if (line_score.size() > 1) {
                 Sort(line_score);
                 best_line_scores[i] = line_score[0];
             } else if (line_score.solo()) {
                 best_line_scores[i] = line_score[0];
             }
         }
         // find non-circle (dead lock) line relationship and insert bc-only
         // gaps, also break down the old ones.
         vec<int> lines_to_fill;
         for (int i = 0; i < best_line_scores.isize(); ++i) {
             // lines with no place to fill at all.
             if (best_line_scores[i].second == -1 ||
                 best_line_scores[i].third == -1) {
                 continue;
             }
             // detect if there is a circle.
             int start = best_line_scores[i].second;
             int next = start;
             Bool is_circle = False;
             vec<int> path({start});
             while (1) {
                 next = best_line_scores[next].second;
                 if (next == -1) {
                     break;
                 } else if (next == start) {
                     is_circle = True;
                     if (verbose) {
                         cout << "circle detected: " 
                             << printSeq(path) << endl;
                     }
                     break;
                 }
                 path.push_back(next);
             }
             if (is_circle) {
                 continue;
             }
             // if reach here, there is no dead lock.
             lines_to_fill.push_back(i);
         }
         
         if (verbose) {
             cout << "there are " << lines_to_fill.size()
                  << " lines to fill (candidates)" << endl;
         }
         if (lines_to_fill.empty()) {
             break;
         }

         // check if there is any collision with lines to fill.
         // first, create a matrix of line-cell-line_to_fill
         vec<map<int, int> > line_cell(dlines.size());
         // vec<vec<int> > break_cell_ids(dlines.isize());
         for (int i = 0; i < best_line_scores.isize(); ++i) {
             triple<double, int, int> &best_line_score = best_line_scores[i];
             int lid = best_line_score.second;
             int cell_id = best_line_score.third;
             if (lid == -1 || cell_id == -1) {
                 continue;
             }
             line_cell[lid][cell_id] = i;
         }

         Bool give_up = True;
         vec<int> dels;
         vec<Bool> touched( D.E(), False );
         for (int i = 0; i < lines_to_fill.isize(); ++i) {
             int lid_to_fill = lines_to_fill[i];
             triple<double, int, int> &best_line_score = best_line_scores[lid_to_fill];
             int lid = best_line_score.second;
             int cell_id = best_line_score.third;
             // gap cell id on reverse strand line
             int cell_id_re = dlines[lid].isize() - cell_id - 1;
             if (line_cell[linv[lid]].find(cell_id_re) == line_cell[linv[lid]].end()) {
                 if (verbose) {
                     cout << "warning type 1! line " << lid_to_fill
                         << " only fill into one strand of line "
                         << lid << endl;
                 }
             } else if (line_cell[linv[lid]][cell_id_re] != linv[lid_to_fill]) {
                 if (verbose) {
                     cout << "warning type 2! line " << lid
                         << " have different lines to fill in on two strands: "
                         << lid_to_fill << " and "
                         << line_cell[linv[lid]][cell_id_re] << endl;
                 }
                 if (llens[lid_to_fill] <= llens[line_cell[linv[lid]][cell_id_re]]) {
                     continue;
                 }
             }
             // another check
             int lid_to_fill_re = linv[lid_to_fill];
             triple<double, int, int> &best_line_score_re = best_line_scores[lid_to_fill_re];
             if (lid_to_fill_re != lid_to_fill &&
                 best_line_score_re.second != -1 &&
                 best_line_score_re.third != -1) {

                 // collision can't be solved.
                 if (best_line_score_re.second != linv[lid] ||
                     best_line_score_re.third != cell_id_re) {
                     if (verbose) {
                         cout << "warning type 3! line " << lid_to_fill << " and "
                             <<  lid_to_fill_re << " have different gap cell to fill"
                             << endl;
                     }
                     continue;
                 }
             }

             give_up = False;
             // break original bc-only gaps.
             int d1 = dlines[lid][cell_id - 1][0][0];
             int d2 = dlines[lid][cell_id + 1][0][0];
             int d3 = dinv[d1], d4 = dinv[d2];
             if (touched[d1] || touched[d2] || touched[d3] || touched[d4]) {
                 continue;
             }
             touched[d1] = touched[d2] = touched[d3] = touched[d4] = True;
             int v = to_right[d1], w = to_left[d2];
             int rv = to_right[ dinv[d2] ], rw = to_left[ dinv[d1] ];
             
             ForceAssert((D.IFrom(v).solo() &&
                          D.ITo(w).solo() &&
                          D.IFrom(rv).solo() &&
                          D.ITo(rw).solo()));
                          
             /*
             if (!(D.IFrom(v).solo() &&
                          D.ITo(w).solo() &&
                          D.IFrom(rv).solo() &&
                          D.ITo(rw).solo()))
             {
                 cout << "nimei: " << D.IFrom(v).size() << ' '
                      << D.ITo(w).size() << ' ' 
                      << D.IFrom(rv).size() << ' '
                      << D.ITo(rw).size() << endl;
                 Scram(1);
             }
             */
             int bc_gap_d = D.IFrom(v, 0);
             ForceAssert(IsBarcodeOnlyGap( D.O(bc_gap_d)));
             dels.push_back(bc_gap_d);
             dels.push_back(dinv[bc_gap_d]);

             // insert new bc-only gaps.
             // get v1, w1 for line_to_fill
             int d_left = dlines[lid_to_fill].front()[0][0];
             int d_right = dlines[lid_to_fill].back()[0][0];
             int v1 = to_left[d_left], w1 = to_right[d_right];
             int rv1 = to_left[ dinv[d_right] ], rw1 = to_right[ dinv[d_left] ];
             ForceAssert((D.ITo(v1).empty() &&
                          D.IFrom(w1).empty() &&
                          D.ITo(rv1).empty() &&
                          D.IFrom(rw1).empty()));

             // add 4 new superedges of bc-only gap.
             dinv.push_back( D.E( ) + 1, D.E( ), D.E() + 3, D.E() + 2 );
             D.AddEdge( v, v1, {-2} );
             D.AddEdge( rw1, rw, {-2} );
             D.AddEdge( w1, w, {-2} );
             D.AddEdge( rv, rv1, {-2} );
         }
         if (give_up) {
             ForceAssert(dels.empty());
             if (verbose) {
                 cout << "no lines can be filled. Terminate." << endl;
             }
             break;
         }
         // clean up current supergraph D.
         D.DeleteEdges(dels);
         RemoveUnneededVertices( D, dinv );
         CleanupCore( D, dinv );
     }
}

void Star4( const HyperBasevectorX& hb, const vec<int>& inv, const vec<Bool>& dup,
     const vec<int32_t>& bc, const ReadPathVecX& paths, 
     digraphE<vec<int>>& D, vec<int>& dinv, 
     ReadPathVec& dpaths, const VecIntVec& ebcx, 
     const vec< triple<int,int,int> >& qept, 
     const double MIN_ADVANTAGE, const Bool OO, 
     const Bool DJANGO, const Bool verbose ) {

     cout << Date( ) << ": begin Star4" << endl;

     // Heuristics.

     const int MIN_STAR = 5000;
     const double MAX_CN_DIFF = 0.5;
     const int MAX_VIEW = 10;
     int MAX_RIGHTS = 6;
     if (OO) MAX_RIGHTS = 4;
     int MAX_LEFTS = 6;
     if (OO) MAX_LEFTS = 4;
     const int MIN_BAR_TO = 2000;
     const int MAX_BAR_TO = 1000000;
     const int BC_VIEW = 50000;

     // Place reads.

     cout << Date( ) << ": placing reads" << endl;
     PlaceReads( hb, paths, dup, D, dpaths, True, False );
     cout << Date( ) << ": making index" << endl;
     IntIndex dpaths_index( dpaths, D.E( ) );


     // Compute kmers.

     vec<int> kmers( hb.E( ) );
     #pragma omp parallel for schedule(dynamic)
     for ( int e = 0; e < hb.E( ); e++ )
         kmers[e] = hb.Kmers(e);

     vec<vec<vec<vec<int>>>> dlines;
     int sp = 0;
     while(1)
     {    
         cout << Date( ) << ": start star pass " << ++sp << endl;
         // Introduce barcode-only gaps.  First step: find line neighbors.
         // Find lines.
         cout << Date( ) << ": finding lines" << endl;
         FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH, True );


         vec<int> llens;
         GetLineLengths( hb, D, dlines, llens );

         cout << Date( ) << ": get line ancillary data" << endl;
         vec<int> linv;
         LineInv( dlines, dinv, linv );
         vec<int> to_left, to_right;
         D.ToLeft(to_left), D.ToRight(to_right);

         // Find (long) lines with bc-only gaps.

         cout << Date( ) << ": begin star stage" << endl;
         vec<int> l1s;
         vec<vec<int> > break_cell_ids(dlines.size());
         for ( int L1 = 0; L1 < dlines.isize( ); L1++ )
         {    
             if ( llens[L1] < MIN_STAR ) continue;
             const vec<vec<vec<int>>>& L = dlines[L1];
             break_cell_ids[L1].push_back(0);
             for (int j = 0; j < L.isize(); ++j) {
                 const vec<vec<int>> &cell = L[j];
                 if (cell.solo() && cell[0].empty()) {
                     int d = L[j - 1][0][0];
                     int e = D.IFrom(to_right[d], 0);
                     if (IsBarcodeOnlyGap( D.O(e))) {
                         break_cell_ids[L1].push_back(j);
                     }
                 }
             }
             break_cell_ids[L1].push_back(L.isize());
             
             if (!break_cell_ids[L1].empty()) {
                 l1s.push_back(L1);
             }
             
             // l1s.push_back(L1);
         }
         if (l1s.empty()) {
             if (verbose) {
                 cout << "there is no line with bc-only barcode. Terminate."
                      << endl;
             }
             break;
         }
         vec<int> l1slen( l1s.isize( ) );
         for ( int i = 0; i < l1s.isize( ); i++ )
             l1slen[i] = llens[ l1s[i] ];

         // Go through long lines and find their nearest left and right neighbors.
         vec<vec<triple<double, int, int> > > line_scores(dlines.size());
         ReverseSortSync( l1slen, l1s );
         cout << Date( ) << ": start star loop on " << l1s.size( )
              << " L1 values" << endl;

         vec<int> dlines_broken_start_pos;
         vec<vec<vec<vec<int>>>> dlines_broken;
         vec<int> broken_to_whole;
         for ( int li = 0; li < l1s.isize( ); li++ ) {
             dlines_broken_start_pos.push_back(dlines_broken.size());
             int lid = l1s[li];
             vec<int> &break_cell_id = break_cell_ids[lid];
             ForceAssertGe(break_cell_id.isize(), 2);
             vec<vec<vec<int>>> &line = dlines[lid];
             for (int i = 0; i < break_cell_id.isize() - 1; ++i) {
                 int cell_id = break_cell_id[i];
                 int cell_id_next = break_cell_id[i + 1];
                 const vec<vec<vec<int>>> broken_line(line.begin() + cell_id,
                                                      line.begin() + cell_id_next);
                 ForceAssert(broken_line.nonempty());
                 dlines_broken.push_back(broken_line);
                 broken_to_whole.push_back(lid);
             }
         }
         // add broken lines size to start_pos
         dlines_broken_start_pos.push_back(dlines_broken.size());
         if (verbose) {
             cout << Date() << " there are " << dlines_broken.size()
                  << " broken lines" << endl;
         }

         // compute cov for original dlines
         vec<double> cov;
         {    
             vec<vec<pair<int,int>>> lbp;
             BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, 0 );
             MasterVec<SerfVec<pair<int,int>>> lbpx;
             for ( auto x : lbp )
             {    
                 SerfVec<pair<int,int>> y( x.size( ) );
                 for ( int j = 0; j < x.isize( ); j++ )
                     y[j] = x[j];
                 lbpx.push_back(y);
             }
             LineCN( kmers, lbpx, D, dlines, llens, cov );    
         }
         // Get barcode positions for original lines
         vec<vec<pair<int,int>>> lbp;
         BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, BC_VIEW );

         // compute cov for broken dlines
         vec<int> llens_broken;
         GetLineLengths( hb, D, dlines_broken, llens_broken );
         vec<double> cov_broken;
         {    
             vec<vec<pair<int,int>>> lbp_broken;
             BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines_broken, lbp_broken, 0 );
             MasterVec<SerfVec<pair<int,int>>> lbpx;
             for ( auto x : lbp_broken )
             {    
                 SerfVec<pair<int,int>> y( x.size( ) );
                 for ( int j = 0; j < x.isize( ); j++ )
                     y[j] = x[j];
                 lbpx.push_back(y);
             }
             LineCN( kmers, lbpx, D, dlines_broken, llens_broken, cov_broken );    
         }
         // Get barcode positions for broken lines
         vec<vec<pair<int,int>>> lbp_broken;
         BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines_broken, lbp_broken, BC_VIEW );

         // compute barcode neighbor
         vec< vec< pair<int,int> > > lhood;
         LineProxSP( hb, inv, ebcx, D, dinv, dlines, dlines_broken, broken_to_whole, qept, lhood );
         // make a copy of lines
         // remember the difference between lines before current line
         // and after it.
         // vec<vec<vec<vec<int>>>> dlines_copy(dlines);
         vec<int> llens_copy(llens);
         llens_copy.append(llens_broken);
         // vec<double> cov_copy(cov);
         vec<vec<pair<int,int>>> lbp_copy(lbp);
         lbp_copy.append(lbp_broken);

         for ( int li = 0; li < l1s.isize( ); li++ ) {
             int start = dlines_broken_start_pos[li];
             int end = dlines_broken_start_pos[li + 1];
             int lid = l1s[li];
             vec<int> &break_cell_id = break_cell_ids[lid];


             ForceAssertEq(break_cell_id.isize(), end - start + 1);
             #pragma omp parallel for schedule(dynamic, 1)
             for (int broken_cell_id = 0; broken_cell_id < break_cell_id.isize(); ++broken_cell_id) {
                 vec< triple<int,int,int> > M;
                 // CAUTION: memory has to be here !!!!!
                 map< vec<int>, double > memory;
                 // Set up logging.

                 ostringstream* outp;
                 if (verbose) outp = new ostringstream;
                 auto DumpOut = [&]
                 {    if (verbose)
                     {
                         #pragma omp critical
                         {    cout << outp->str( );    }    
                         delete outp;    }    };


                 if (broken_cell_id == 0) {
                     int rlid = start;
                     int v = to_left[ dlines_broken[rlid].front( )[0][0] ]; 
                     if ( v < 0 )
                     {    DumpOut( );
                         continue;    }
                     if ( D.To(v).nonempty( ) || !D.From(v).solo( ) )
                     {    DumpOut( );
                         continue;    }
                     // only consider left
                     // assign ids for rhs
                     int R = dlines.isize() + rlid;
                     // llens_copy.push_back(llens_broken[rlid]);
                     // lbp_copy.push_back(lbp_broken[rlid]);
                     vec<int> lefts;
                     for (int i = 0; i < lhood[rlid].isize(); ++i) {
                         int lid = lhood[rlid][i].second;
                         if ( lid == linv[lid]) continue;
                         if ( AbsDiff( cov_broken[rlid], cov[lid] ) > MAX_CN_DIFF ) continue;
                         if ( llens_copy[lid] < MIN_BAR_TO ) continue;
                         vec<double> scores;
                         scores.push_back( MemoryScoreOrder( 
                                     { R, lid }, lbp_copy, llens_copy, M, memory ) );
                         scores.push_back( MemoryScoreOrder( 
                                     { R, linv[lid] }, lbp_copy, llens_copy, M, memory ) );
                         scores.push_back( MemoryScoreOrder( 
                                     { lid, R }, lbp_copy, llens_copy, M, memory ) );
                         scores.push_back( MemoryScoreOrder( 
                                     { linv[lid], R }, lbp_copy, llens_copy, M, memory ) );
                         vec<int> ids( 4, vec<int>::IDENTITY );
                         SortSync( scores, ids );
                         double ad = scores[1] - scores[0];
                         if ( ad < MIN_ADVANTAGE ) continue;
                         if ( ids[0] <= 1 ) continue;
                         if ( ids[0] == 2 ) lefts.push_back(lid);
                         else lefts.push_back( linv[lid] );
                     }

                     if ( lefts.isize( ) > MAX_LEFTS )
                     {    vec<int> lens( lefts.size( ) );
                         for ( int i = 0; i < lefts.isize( ); i++ )
                             lens[i] = llens_copy[ lefts[i] ];
                         ReverseSortSync( lens, lefts );
                         lefts.resize(MAX_LEFTS);    }

                     if ( lefts.empty() )
                     {   DumpOut();
                         continue;    }
                     int L = -1;
                     if ( lefts.size( ) == 1 ) L = lefts[0];
                     else if ( !DJANGO )
                     {    
                         vec<int> blefts;
                         double ad;
                         if (OO)
                         {    
#if !USE_PIVOT
                             MemoryOrderAndOrientN(lefts, lbp_copy, llens_copy, 
                                     linv, M, memory, ad, blefts );    
#else                        
                             MemoryOrderAndOrientN( R, lefts, lbp_copy, llens_copy, 
                                     linv, M, memory, ad, blefts, True );    
#endif

                         }
                         else {
#if !USE_PIVOT
                             MemoryOrderN(lefts, lbp_copy, llens_copy, M, memory, ad, blefts );
#else
                             MemoryOrderN( R, lefts, lbp_copy, llens_copy, M, memory, ad, blefts, True );
#endif
                         }
                         if ( ad < MIN_ADVANTAGE )
                         {    DumpOut( );
                             continue;    }
                         L = blefts[0];    }
                     else
                     {    
                         vec<int> lens( lefts.size( ) );
                         for ( int i = 0; i < lefts.isize( ); i++ )
                             lens[i] = llens[ lefts[i] ];
                         ReverseSortSync( lens, lefts );
                         while( lefts.size( ) > 1 )
                         {    vec<int> blefts;
                             double ad;
                             if (OO)
                             {    
#if !USE_PIVOT
                                 MemoryOrderAndOrientN(lefts, lbp_copy, llens_copy, 
                                         linv, M, memory, ad, blefts );
#else
                                 MemoryOrderAndOrientN( R, lefts, lbp_copy, llens_copy, 
                                         linv, M, memory, ad, blefts, True );
#endif

                             }
                             else 
                             {    
#if !USE_PIVOT
                                 MemoryOrderN( lefts, lbp_copy, llens_copy, M, memory, ad, blefts );
#else
                                 MemoryOrderN( R, lefts, lbp_copy, llens_copy, M, memory, ad, blefts, True );
#endif
                             }
                             if ( ad >= MIN_ADVANTAGE )
                             {    L = blefts[0];
                                 break;    }
                             lefts.pop_back( );    }    }

                     // Check for done and otherwise save.

                     if ( L == -1 )
                     {    DumpOut( );
                         continue;    }
                     int w = to_right[ dlines[L].back( )[0][0] ]; 
                     if ( w < 0 )
                     {    DumpOut( );
                         continue;    }
                     if ( !D.To(w).solo( ) || D.From(w).nonempty( ) )
                     {    DumpOut( );
                         continue;    }
                     if (verbose) {
                         (*outp) << "line " << lid << " gap cell "
                             << break_cell_id[broken_cell_id]
                             << " can be left inserted with line " << L << endl;
                     }
                     DumpOut( );
                     #pragma omp critical
                     {   
                         // dummy score
                         line_scores[L].push(0.0, lid, break_cell_id[broken_cell_id]);
                     }

                 } else if (broken_cell_id == break_cell_id.isize() - 1) {
                     // check if right side of left is sink
                     int llid = start + break_cell_id.isize() - 1;
                     int v = to_right[ dlines_broken[llid].back( )[0][0] ]; 
                     if ( v < 0 )
                     {    DumpOut( );
                         continue;    }
                     if ( !D.To(v).solo( ) || D.From(v).nonempty( ) )
                     {    DumpOut( );
                         continue;    }

                     int L = dlines.isize() + llid;
                     // llens_copy.push_back(llens_broken[llid]);
                     // lbp_copy.push_back(lbp_broken[llid]);
                     vec<int> rights;
                     for (int i = 0; i < lhood[llid].isize(); ++i) {
                         int lid = lhood[llid][i].second;
                         if ( lid == linv[lid]) continue;
                         if ( AbsDiff( cov_broken[llid], cov[lid] ) > MAX_CN_DIFF ) continue;
                         if ( llens_copy[lid] < MIN_BAR_TO ) {
                             continue;
                         }
                         vec<double> scores;
                         scores.push_back( MemoryScoreOrder( 
                                     { lid, L }, lbp_copy, llens_copy, M, memory ) );
                         scores.push_back( MemoryScoreOrder( 
                                     { linv[lid], L }, lbp_copy, llens_copy, M, memory ) );
                         scores.push_back( MemoryScoreOrder( 
                                     { L, lid }, lbp_copy, llens_copy, M, memory ) );
                         scores.push_back( MemoryScoreOrder( 
                                     { L, linv[lid] }, lbp_copy, llens_copy, M, memory ) );
                         vec<int> ids( 4, vec<int>::IDENTITY );
                         SortSync( scores, ids );
                         double ad = scores[1] - scores[0];
                         if ( ad < MIN_ADVANTAGE ) continue;
                         if ( ids[0] <= 1 ) continue;
                         if ( ids[0] == 2 ) rights.push_back(lid);
                         else rights.push_back( linv[lid] );
                     }
                     if ( rights.isize( ) > MAX_RIGHTS )
                     {    vec<int> lens( rights.size( ) );
                         for ( int i = 0; i < rights.isize( ); i++ )
                             lens[i] = llens_copy[ rights[i] ];
                         ReverseSortSync( lens, rights );
                         rights.resize(MAX_RIGHTS);    }

                     if ( rights.empty() )
                     {   DumpOut();
                         continue;    }
                     int R = -1;
                     if ( rights.size( ) == 1 ) R = rights[0];
                     else if ( !DJANGO )
                     {    vec<int> brights;
                         double ad;
                         if (OO)
                         {    
#if !USE_PIVOT
                             MemoryOrderAndOrientN(rights, lbp_copy, llens_copy, 
                                     linv, M, memory, ad, brights );    
#else                        
                             MemoryOrderAndOrientN( L, rights, lbp_copy, llens_copy, 
                                     linv, M, memory, ad, brights );    
#endif

                         }
                         else {
#if !USE_PIVOT
                             MemoryOrderN(rights, lbp_copy, llens_copy, M, memory, ad, brights );
#else
                             MemoryOrderN( L, rights, lbp_copy, llens_copy, M, memory, ad, brights );
#endif
                         }
                         if ( ad < MIN_ADVANTAGE )
                         {    DumpOut( );
                             continue;    }
                         R = brights[0];    }
                     else
                     {    vec<int> lens( rights.size( ) );
                         for ( int i = 0; i < rights.isize( ); i++ )
                             lens[i] = llens[ rights[i] ];
                         ReverseSortSync( lens, rights );
                         while( rights.size( ) > 1 )
                         {    vec<int> brights;
                             double ad;
                             if (OO)
                             {    
#if !USE_PIVOT
                                 MemoryOrderAndOrientN(rights, lbp_copy, llens_copy, 
                                         linv, M, memory, ad, brights );
#else
                                 MemoryOrderAndOrientN( L, rights, lbp_copy, llens_copy, 
                                         linv, M, memory, ad, brights );
#endif

                             }
                             else 
                             {    
#if !USE_PIVOT
                                 MemoryOrderN( rights, lbp_copy, llens_copy, M, memory, ad, brights );
#else
                                 MemoryOrderN( L, rights, lbp_copy, llens_copy, M, memory, ad, brights );
#endif
                             }
                             if ( ad >= MIN_ADVANTAGE )
                             {    R = brights[0];
                                 break;    }
                             rights.pop_back( );    }    }

                     // Check for done and otherwise save.

                     if ( R == -1 )
                     {    DumpOut( );
                         continue;    }
                     int w = to_left[ dlines[R].front( )[0][0] ]; 
                     if ( w < 0 )
                     {    DumpOut( );
                         continue;    }
                     if ( D.To(w).nonempty( ) || !D.From(w).solo( ) )
                     {    DumpOut( );
                         continue;    }
                     if (verbose) {
                         (*outp) << "line " << lid << " gap cell "
                             << break_cell_id[broken_cell_id]
                             << " can be right inserted with line " << R << endl;
                     }
                     DumpOut( );
                     #pragma omp critical
                     {   
                         // dummy score
                         line_scores[R].push(0.0, lid, break_cell_id[broken_cell_id]);
                     }
                 } else {
                     // assign ids for lhs and rhs
                     int llid = start + broken_cell_id - 1;
                     int rlid = llid + 1;
                     int L = dlines.isize() + llid;
                     // dlines_copy.push_back(dlines_broken[llid]);
                     int R = L + 1;
                     // dlines_copy.push_back(dlines_broken[rlid]);

                     // llens_copy.push_back(llens_broken[llid],
                     //         llens_broken[rlid]);

                     /*
                        cov_copy.push_back(cov_broken[llid],
                        cov_broken[rlid]);
                        */
                     // Get barcode positions.
                     // lbp_copy.push_back(lbp_broken[llid],
                     //         lbp_broken[rlid]);

                     // vec<int> linv_copy;
                     // LineInv( dlines_copy, dinv, linv_copy );

                     // find intersection of L's and R's lhood.
                     vec<int> line_pool, lefts, rights;

                     for (int i = 0; i < lhood[llid].isize(); ++i) {
                         int lid = lhood[llid][i].second;
                         /*
                            int w = to_left[ dlines[lid].front( )[0][0] ]; 
                            if ( w < 0 )
                            {    
                            continue;    }
                            if ( D.To(w).nonempty( ) || !D.From(w).solo( ) )
                            {   
                            continue;    }
                            */
                         if ( lid == linv[lid]) continue;
                         if ( AbsDiff( cov_broken[llid], cov[lid] ) > MAX_CN_DIFF ) continue;
                         if ( llens_copy[lid] < MIN_BAR_TO ||
                                 llens_copy[lid] > MAX_BAR_TO ) {
                             continue;
                         }
                         vec<double> scores;
                         scores.push_back( MemoryScoreOrder( 
                                     { lid, L }, lbp_copy, llens_copy, M, memory ) );
                         scores.push_back( MemoryScoreOrder( 
                                     { linv[lid], L }, lbp_copy, llens_copy, M, memory ) );
                         scores.push_back( MemoryScoreOrder( 
                                     { L, lid }, lbp_copy, llens_copy, M, memory ) );
                         scores.push_back( MemoryScoreOrder( 
                                     { L, linv[lid] }, lbp_copy, llens_copy, M, memory ) );
                         vec<int> ids( 4, vec<int>::IDENTITY );
                         SortSync( scores, ids );
                         /*
                            double ad = scores[1] - scores[0];
                            if ( ad < MIN_ADVANTAGE ) continue;
                            */
                         if ( ids[0] <= 1 ) continue;
                         rights.push_back(lid);
                         /*
                            if ( ids[0] == 2 ) rights.push_back(lid);
                            else rights.push_back( linv[lid] );
                            */
                     }

                     if ( rights.empty() )
                     {   DumpOut();
                         continue;    }

                     for (int i = 0; i < lhood[rlid].isize(); ++i) {
                         int lid = lhood[rlid][i].second;
                         if ( lid == linv[lid]) continue;
                         if ( AbsDiff( cov_broken[rlid], cov[lid] ) > MAX_CN_DIFF ) continue;
                         if ( llens_copy[lid] < MIN_BAR_TO || 
                                 llens_copy[lid] > MAX_BAR_TO ) {
                             continue;
                         }
                         vec<double> scores;
                         scores.push_back( MemoryScoreOrder( 
                                     { R, lid }, lbp_copy, llens_copy, M, memory ) );
                         scores.push_back( MemoryScoreOrder( 
                                     { R, linv[lid] }, lbp_copy, llens_copy, M, memory ) );
                         scores.push_back( MemoryScoreOrder( 
                                     { lid, R }, lbp_copy, llens_copy, M, memory ) );
                         scores.push_back( MemoryScoreOrder( 
                                     { linv[lid], R }, lbp_copy, llens_copy, M, memory ) );
                         vec<int> ids( 4, vec<int>::IDENTITY );
                         SortSync( scores, ids );
                         /*
                            double ad = scores[1] - scores[0];
                            if ( ad < MIN_ADVANTAGE ) continue;
                            */
                         if ( ids[0] <= 1 ) continue;
                         lefts.push_back(lid);
                         /*
                            if ( ids[0] == 2 ) lefts.push_back(lid);
                            else lefts.push_back( linv[lid] );
                            */
                     }

                     if ( lefts.empty() )
                     {   DumpOut();
                         continue;    }

                     for (int i = 0; i < lefts.isize(); ++i) {
                         int lid = lefts[i];
                         if (Position(rights, lid) == -1 || lid == linv[lid]) {
                             continue;
                         }
                         line_pool.push_back(lid);
                     }

                     if ( line_pool.empty() )
                     {   DumpOut();
                         continue;    }

                     vec<double> scores;
                     for (int i = 0; i < line_pool.isize(); ++i) {
                         // LB means Line Between
                         int LB = line_pool[i];
                         // if ( llens_copy[LB] < MIN_BAR_TO ) continue;
                         // if ( AbsDiff( cov_broken[llid], cov[LB] ) > MAX_CN_DIFF ) continue;
                         // if ( AbsDiff( cov_broken[rlid], cov[LB] ) > MAX_CN_DIFF ) continue;
                         scores.push_back( MemoryScoreOrder( 
                                     { L, LB , R }, lbp_copy, llens_copy, M, memory ) );
                         scores.push_back( MemoryScoreOrder( 
                                     { L, linv[LB] , R }, lbp_copy, llens_copy, M, memory ) );
                     }
                     if (scores.empty()) {
                         DumpOut();
                         continue;    
                     }
                     vec<int> ids( scores.size(), vec<int>::IDENTITY );
                     SortSync( scores, ids );
                     double ad = scores[1] - scores[0];
                     if (ad == 0 && (line_pool[ids[0]/2] == linv[line_pool[ids[1]/2]])) {
                         if (verbose) (*outp) << "REALLY?" << endl;
                         DumpOut();
                         continue;
                     } else if ( ad < MIN_ADVANTAGE ) {
                         DumpOut();
                         continue;
                     }
                     int best_oo = line_pool[ids[0] / 2];
                     if (ids[0] % 2 != 0) {
                         best_oo = linv[best_oo];
                     }
                     ForceAssertLt(best_oo, dlines.isize());
                     int v = to_left[ dlines[best_oo].front( )[0][0] ]; 
                     int w = to_right[ dlines[best_oo].back( )[0][0] ]; 
                     if (v < 0 || w < 0) {
                         DumpOut();
                         continue;
                     }
                     if ( D.To(v).nonempty( ) || !D.From(v).solo( ) ) {
                         DumpOut();
                         continue;
                     }
                     if ( D.From(w).nonempty( ) || !D.To(w).solo( ) ) {
                         DumpOut();
                         continue;
                     }

                     if (verbose) {
                         (*outp) << "line " << lid << " gap cell "
                             << break_cell_id[broken_cell_id]
                             << " can be filled with line " << best_oo << endl;
                     }
                     DumpOut();
                     #pragma omp critical
                     {
                         // <score, lid, cell_id>
                         line_scores[best_oo].push(scores[0], lid, break_cell_id[broken_cell_id]);
                     }
                 }
             }
         }
         // introduce joins.
         if (verbose) {
             cout << Date() << " start to introduce joins" << endl;
         }
         vec<triple<double, int, int> > best_line_scores(dlines.size(), triple<double, int, int> (-1, -1, -1));
         for (int i = 0; i < line_scores.isize(); ++i) {
             vec<triple<double, int, int> > &line_score = line_scores[i];
             if (verbose && line_score.nonempty()) {
                 cout << "scores for line " << i << " is (are): "
                     << printSeq(line_score) << endl;
             }
             /*
                if (line_score.size() > 1) {
                Sort(line_score);
                best_line_scores[i] = line_score[0];
                } else if (line_score.solo()) {
                best_line_scores[i] = line_score[0];
                }
                */
             if (line_score.solo()) {
                 best_line_scores[i] = line_score[0];
             }
         }
         // find non-circle (dead lock) line relationship and insert bc-only
         // gaps, also break down the old ones.
         vec<int> lines_to_fill;
         for (int i = 0; i < best_line_scores.isize(); ++i) {
             vec<Bool> touched_left( D.E(), False );
             vec<Bool> touched_right( D.E(), False );
             // lines with no place to fill at all.
             if (best_line_scores[i].second == -1 ||
                     best_line_scores[i].third == -1) {
                 continue;
             }
             // detect if there is a circle.
             int start = i;
             triple<double, int, int> &line_score = best_line_scores[start];
             int next = line_score.second;
             int cell_id = line_score.third;
             if (cell_id == 0) {
                 touched_right[i] = True;
             } else if (cell_id == dlines[next].isize()) {
                 touched_left[i] = True;
             } else {
                 touched_left[i] = True;
                 touched_right[i] = True;
             }
             Bool is_circle = False;
             vec<int> path({start});
             while (1) {
                 triple<double, int, int>  &line_score = best_line_scores[next];
                 next = best_line_scores[next].second;
                 if (next == -1) {
                     break;
                 } else if (next == start) {
                     is_circle = True;
                     if (verbose) {
                         cout << "circle detected: " 
                             << printSeq(path) << endl;
                     }
                     break;
                 }
                 path.push_back(next);
             }
             if (is_circle) {
                 continue;
             }
             // if reach here, there is no dead lock.
             lines_to_fill.push_back(i);
         }

         if (verbose) {
             cout << "there are " << lines_to_fill.size()
                 << " lines to fill (candidates)" << endl;
         }
         if (lines_to_fill.empty()) {
             break;
         }

         // check if there is any collision with lines to fill.
         // first, create a matrix of line-cell-line_to_fill
         vec<map<int, int> > line_cell(dlines.size());
         // vec<vec<int> > break_cell_ids(dlines.isize());
         for (int i = 0; i < best_line_scores.isize(); ++i) {
             triple<double, int, int> &best_line_score = best_line_scores[i];
             int lid = best_line_score.second;
             int cell_id = best_line_score.third;
             if (lid == -1 || cell_id == -1) {
                 continue;
             }
             line_cell[lid][cell_id] = i;
         }

         Bool give_up = True;
         vec<int> dels;
         vec<Bool> touched( D.E(), False );
         for (int i = 0; i < lines_to_fill.isize(); ++i) {
             int lid_to_fill = lines_to_fill[i];
             triple<double, int, int> &best_line_score = best_line_scores[lid_to_fill];
             int lid = best_line_score.second;
             int cell_id = best_line_score.third;
             // gap cell id on reverse strand line
             int cell_id_re = dlines[lid].isize() - cell_id - 1;
             if (line_cell[linv[lid]].find(cell_id_re) == line_cell[linv[lid]].end()) {
                 if (verbose) {
                     cout << "warning type 1! line " << lid_to_fill
                         << " only fill into one strand of line "
                         << lid << endl;
                 }
             } else if (line_cell[linv[lid]][cell_id_re] != linv[lid_to_fill]) {
                 if (verbose) {
                     cout << "warning type 2! line " << lid
                         << " have different lines to fill in on two strands: "
                         << lid_to_fill << " and "
                         << line_cell[linv[lid]][cell_id_re] << endl;
                 }
                 if (llens[lid_to_fill] <= llens[line_cell[linv[lid]][cell_id_re]]) {
                     continue;
                 }
             }
             // another check
             int lid_to_fill_re = linv[lid_to_fill];
             triple<double, int, int> &best_line_score_re = best_line_scores[lid_to_fill_re];
             if (lid_to_fill_re != lid_to_fill &&
                     best_line_score_re.second != -1 &&
                     best_line_score_re.third != -1) {

                 // collision can't be solved.
                 if (best_line_score_re.second != linv[lid] ||
                         best_line_score_re.third != cell_id_re) {
                     if (verbose) {
                         cout << "warning type 3! line " << lid_to_fill << " and "
                             <<  lid_to_fill_re << " have different gap cell to fill"
                             << endl;
                     }
                     continue;
                 }
             }

             give_up = False;
             if (cell_id == 0) {
                 // insert line_to_fill on the left
                 int d_left = dlines[lid_to_fill].back()[0][0]; 
                 int d_left_r = dinv[d_left];
                 int d_right = dlines[lid].front()[0][0];
                 int d_right_r = dinv[d_right];
                 if (touched[d_left] || touched[d_left_r] || touched[d_right] || touched[d_right_r]) {
                     continue;
                 }
                 touched[d_left] = touched[d_left_r] = touched[d_right] = touched[d_right_r] = True;
                 int v = to_right[d_left], w = to_left[d_right];
                 int rv = to_right[d_right_r], rw = to_left[d_left_r];
                 ForceAssert(D.IFrom(v).empty() && D.ITo(v).solo());
                 ForceAssert(D.ITo(w).empty() && D.IFrom(w).solo());
                 ForceAssert(D.IFrom(rv).empty() && D.ITo(rv).solo());
                 ForceAssert(D.ITo(rw).empty() && D.IFrom(rw).solo());
                 /*
                 if (D.IFrom(rv).nonempty() || !D.ITo(rv).solo()) {
                     cout << "left nimei rv" << D.IFrom(rv).size() << ' ' << D.ITo(rv).size() << endl;
                 }
                 if (D.ITo(rw).nonempty() || !D.IFrom(rw).solo()) {
                     cout << "left nimei rw" << D.IFrom(rv).size() << ' ' << D.ITo(rv).size() << endl;
                 }
                 */
                 dinv.push_back( D.E( ) + 1, D.E( ) );
                 D.AddEdge( v, w, {-2} );
                 D.AddEdge( rv, rw, {-2} );
             } else if (cell_id == break_cell_ids[lid].back()) {
                 // insert line_to fill on the right
                 int d_left = dlines[lid].back()[0][0]; 
                 int d_left_r = dinv[d_left];
                 int d_right = dlines[lid_to_fill].front()[0][0];
                 int d_right_r = dinv[d_right];
                 if (touched[d_left] || touched[d_left_r] || touched[d_right] || touched[d_right_r]) {
                     continue;
                 }
                 touched[d_left] = touched[d_left_r] = touched[d_right] = touched[d_right_r] = True;
                 int v = to_right[d_left], w = to_left[d_right];
                 int rv = to_right[d_right_r], rw = to_left[d_left_r];
                 ForceAssert(D.IFrom(v).empty() && D.ITo(v).solo());
                 ForceAssert(D.ITo(w).empty() && D.IFrom(w).solo());
                 // ForceAssert(D.IFrom(rv).empty() && D.ITo(rv).solo());
                 // ForceAssert(D.ITo(rw).empty() && D.IFrom(rw).solo());
                 if (D.IFrom(rv).nonempty() || !D.ITo(rv).solo()) {
                     cout << "right nimei rv" << D.IFrom(rv).size() << ' ' << D.ITo(rv).size() << endl;
                 }
                 if (D.ITo(rw).nonempty() || !D.IFrom(rw).solo()) {
                     cout << "right nimei rw" << D.IFrom(rv).size() << ' ' << D.ITo(rv).size() << endl;
                 }
                 dinv.push_back( D.E( ) + 1, D.E( ) );
                 D.AddEdge( v, w, {-2} );
                 D.AddEdge( rv, rw, {-2} );
             } else {
                 int d1 = dlines[lid][cell_id - 1][0][0];
                 int d2 = dlines[lid][cell_id + 1][0][0];
                 int d3 = dinv[d1], d4 = dinv[d2];
                 if (touched[d1] || touched[d2] || touched[d3] || touched[d4]) {
                     continue;
                 }
                 touched[d1] = touched[d2] = touched[d3] = touched[d4] = True;
                 int v = to_right[d1], w = to_left[d2];
                 int rv = to_right[ dinv[d2] ], rw = to_left[ dinv[d1] ];
                 ForceAssert((D.IFrom(v).solo() &&
                             D.ITo(w).solo() &&
                             D.IFrom(rv).solo() &&
                             D.ITo(rw).solo()));

                 int d_left = dlines[lid_to_fill].front()[0][0];
                 int d_left_r = dinv[d_left];
                 int d_right = dlines[lid_to_fill].back()[0][0];
                 int d_right_r = dinv[d_right];
                 if (touched[d_left] || touched[d_left_r] || touched[d_right] || touched[d_right_r]) {
                     continue;
                 }
                 touched[d_left] = touched[d_left_r] = touched[d_right] = touched[d_right_r] = True;
                 // get v1, w1 for line_to_fill
                 int v1 = to_left[d_left], w1 = to_right[d_right];
                 int rv1 = to_left[ dinv[d_right] ], rw1 = to_right[ dinv[d_left] ];
                 ForceAssert((D.ITo(v1).empty() &&
                             D.IFrom(w1).empty() &&
                             D.ITo(rv1).empty() &&
                             D.IFrom(rw1).empty()));

                 // break original bc-only gaps.
                 int bc_gap_d = D.IFrom(v, 0);
                 ForceAssert(IsBarcodeOnlyGap( D.O(bc_gap_d)));
                 dels.push_back(bc_gap_d);
                 dels.push_back(dinv[bc_gap_d]);

                 // insert new bc-only gaps.
                 // add 4 new superedges of bc-only gap.
                 dinv.push_back( D.E( ) + 1, D.E( ), D.E() + 3, D.E() + 2 );
                 D.AddEdge( v, v1, {-2} );
                 D.AddEdge( rw1, rw, {-2} );
                 D.AddEdge( w1, w, {-2} );
                 D.AddEdge( rv, rv1, {-2} );
             }
         }
         if (give_up) {
             ForceAssert(dels.empty());
             if (verbose) {
                 cout << "no lines can be filled. Terminate." << endl;
             }
             break;
         }
         // clean up current supergraph D.
         D.DeleteEdges(dels);
         RemoveUnneededVertices( D, dinv );
         CleanupCore( D, dinv );
     }
}

Bool IsKillable(vec< pair<int,int> >  &orc_L,
                vec< pair<int,int> >  &orc_R,
                int len_L,
                int len_R,
                double winpos,
                const int BC_REQUIRE,
                const int BC_FLANK,
                const int BC_IGNORE,
                const int BC_MIN) {
    vec<int> bleft, bright;
    if (len_L >= BC_REQUIRE && len_R >= BC_REQUIRE)
    {
        // Compute number of bridging barcodes.

        int64_t low = LowerBound1( orc_L, len_L - BC_FLANK );
        int64_t high = UpperBound1( orc_L, len_L - BC_IGNORE );
        for ( int64_t j = low; j < high; j++ )
        {    
            int start = orc_L[j].first, b = orc_L[j].second;
            bleft.push_back(b);
        }
        low = LowerBound1( orc_R, BC_IGNORE );
        high = UpperBound1( orc_R, BC_FLANK );
        for ( int64_t j = low; j < high; j++ )
        {    
            int start = orc_R[j].first, b = orc_R[j].second;
            bright.push_back(b);
        }
        int nleft = bleft.size( ), nright = bright.size( );
        #pragma omp critical
        cout << "nleft vs. nright " << nleft << ':' << nright << endl;
        int n = Min( nleft, nright );
        UniqueSort(bleft), UniqueSort(bright);
        int bridge = MeetSize( bleft, bright );
        int expect = Min( 1.0, n/winpos ) * BC_MIN;
        #pragma omp critical
        cout << Date() << ": expect vs. bridge " << expect << ':' << bridge << endl;
        if ( bridge < expect ) {
            #pragma omp critical
            cout << "KILL BEFOREHAND!: " << bridge << ':' << expect << endl;
            return True;
        }
    }
    return False;
}

void ScaffoldMe( const HyperBasevectorX& hb, const vec<int>& inv, const vec<Bool>& dup,
     const vec<int32_t>& bc, const ReadPathVecX& paths, 
     digraphE<vec<int> >& D, vec<int>& dinv, 
     ReadPathVec& dpaths, const VecIntVec& ebcx, 
     const vec< triple<int,int,int> >& qept, 
     const double MIN_ADVANTAGE, const Bool OO, 
     const Bool DJANGO, 
     const int BC_REQUIRE, int BC_FLANK, int BC_IGNORE,
     const Bool require_symmetric, const Bool verbose ) {

     cout << Date( ) << ": begin ScaffoldMe" << endl;

     // Heuristics.

     const int MIN_STAR = 4000;
     const double MAX_CN_DIFF = 0.5;
     const int MAX_VIEW = 10;
     int MAX_RIGHTS = 6;
     if (OO) MAX_RIGHTS = 4;
     int MAX_LEFTS = 6;
     if (OO) MAX_LEFTS = 4;
     const int MIN_BAR_TO = 2000;
     const int BC_VIEW = 50000;
     int BC_MIN = 10;


     // Compute kmers.

     vec<int> kmers( hb.E( ) );
     #pragma omp parallel for schedule(dynamic)
     for ( int e = 0; e < hb.E( ); e++ )
         kmers[e] = hb.Kmers(e);

     vec<vec<vec<vec<int>>>> dlines;
     int sp = 0;
     while(1)
     {    
         cout << Date( ) << ": start ScaffoldMe pass " << ++sp << endl;
         // Introduce barcode-only gaps.  First step: find line neighbors.

         // Place reads.

         cout << Date( ) << ": placing reads" << endl;
         PlaceReads( hb, paths, dup, D, dpaths, True, False );
         cout << Date( ) << ": making index" << endl;
         IntIndex dpaths_index( dpaths, D.E( ) );
         // Find lines.
         cout << Date( ) << ": finding lines" << endl;
         FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH, True );
         cout << Date( ) << ": N50 line = " 
             << ToStringAddCommas( LineN50( hb, D, dlines ) ) << endl;


         vec<int> llens;
         GetLineLengths( hb, D, dlines, llens );

         cout << Date( ) << ": get line ancillary data" << endl;
         vec<int> linv;
         LineInv( dlines, dinv, linv );
         vec<int> to_left, to_right;
         D.ToLeft(to_left), D.ToRight(to_right);
         // compute cov for original dlines
         vec<double> cov;
         vec< vec< pair<int,int> > > orcs(dlines.size());
         vec<vec<pair<int,int>>> lbp_all;
         {    
             BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp_all, 0 );
             MasterVec<SerfVec<pair<int,int>>> lbpx;
             for ( auto x : lbp_all )
             {    
                 SerfVec<pair<int,int>> y( x.size( ) );
                 for ( int j = 0; j < x.isize( ); j++ )
                     y[j] = x[j];
                 lbpx.push_back(y);
             }
             LineCN( kmers, lbpx, D, dlines, llens, cov );    

             for (int i = 0; i < orcs.isize(); ++i) {
                 vec< pair<int,int> > &orc = orcs[i];
                 orc.resize( lbp_all[i].size( ) );
                 for ( int j = 0; j < orc.isize( ); j++ )
                     orc[j] = make_pair( lbp_all[i][j].second, lbp_all[i][j].first );
                 Sort(orc);
             }
         }
         // Get barcode positions for original lines
         vec<vec<pair<int,int>>> lbp;
         BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, BC_VIEW );
         vec< vec< pair<int,int> > > lhood;
         LineProx( hb, inv, ebcx, D, dinv, dlines, qept, lhood );

         // Estimate fragment length distribution.

         const int MAX_GAP_FRAG = 60000;
         const int MIN_READS_FRAG = 4;
         const int MIN_FRAG_LEN = 1000;
         int mol_lwml=0; 
         // Define local block so lrpb can be destroyed
         // by going out of scope
         {
             vec <double> fhist;
             vec <pair<float, int> > lr;
             // compute read, pos on lines
             vec< vec<triple<int32_t, int, int64_t> > > lrpb;
             ReadPosLine( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lrpb, 0 );

             FindMoleculesOnLines( dlines, llens, lrpb, fhist, lr, mol_lwml,
                     MAX_GAP_FRAG, MIN_READS_FRAG, MIN_FRAG_LEN);
         }
         cout << Date() << ": mol_lwml = " << mol_lwml << endl;
         const int MAX_BAR_TO = mol_lwml / 3;

         // Scale down if parameters are not supported by the fragment size.

         BC_IGNORE = Min( BC_IGNORE, mol_lwml/4 );

         // Get mean number of positions in a window.

         int64_t total_bases = 0, total_pos = 0;
         for ( int i = 0; i < dlines.isize( ); i++ )
         {    if ( llens[i] < BC_FLANK ) continue;
             total_bases += llens[i];
             total_pos += lbp_all[i].size( );    }
         double winpos
             = double(BC_FLANK-BC_IGNORE) * double(total_pos) / double(total_bases);
         cout << Date( ) << ": " << winpos << " positions per window" << endl;

         // Find (long) lines with bc-only gaps.

         vec<int> l1s;
         for ( int L1 = 0; L1 < dlines.isize( ); L1++ )
         {    if ( llens[L1] < MIN_STAR ) continue;
             if ( linv[L1] == L1) continue;
             int v = to_right[ dlines[L1].back( )[0][0] ];
             if ( v < 0 ) { continue; }
             if ( D.From(v).nonempty( ) || !D.To(v).solo( ) ) continue;
             /*
             if ( TooMuchDupContentSingle(D, dlines[L1], dinv)) {
                 continue;
             }
             */
             l1s.push_back(L1);
         }
         vec<int> l1slen( l1s.isize( ) );
         for ( int i = 0; i < l1s.isize( ); i++ )
             l1slen[i] = llens[ l1s[i] ];
         ReverseSortSync( l1slen, l1s );
         cout << Date( ) << ": start star loop on " << l1s.size( )
             << " L1 values" << endl;
         map< vec<int>, double > memory;
         const int batch = Max( 1, l1s.isize( )/500 );
         vec<vec<pair<double, int> > > confident_lefts(dlines.size());
         vec<int> dels;
         #pragma omp parallel for schedule(dynamic, batch)
         for ( int li = 0; li < l1s.isize( ); li++ )
         {    
             // Set up logging.

             ostringstream* outp;
             if (verbose) outp = new ostringstream;
             auto DumpOut = [&]
             {    if (verbose)
                 {
                     #pragma omp critical
                     {    cout << outp->str( );    }    
                     delete outp;    }    };
             // check if right side of left is sink
             int L = l1s[li];
             vec< triple<int,int,int> > M;
             vec<int> rights;
             for (int i = 0; i < lhood[L].isize(); ++i) {
                 int lid = lhood[L][i].second;
                 ForceAssert(linv[lid] != L);
                 ForceAssert(lid != L);
                 if ( lid == linv[lid]) continue;
                 if ( AbsDiff( cov[L], cov[lid] ) > MAX_CN_DIFF ) continue;
                 if ( llens[lid] < MIN_BAR_TO ) {
                     continue;
                 }
                 /*
                 if ( TooMuchDupContentSingle(D, dlines[lid], dinv)) {
                     continue;
                 }
                 */
                 
                 if (TooMuchDupContent(D, dlines[L], dlines[lid], inv)) {
                     continue;
                 }
                 
                 vec<double> scores;
                 scores.push_back( MemoryScoreOrder( 
                             { lid, L }, lbp, llens, M, memory ) );
                 scores.push_back( MemoryScoreOrder( 
                             { linv[lid], L }, lbp, llens, M, memory ) );
                 scores.push_back( MemoryScoreOrder( 
                             { L, lid }, lbp, llens, M, memory ) );
                 scores.push_back( MemoryScoreOrder( 
                             { L, linv[lid] }, lbp, llens, M, memory ) );
                 vec<int> ids( 4, vec<int>::IDENTITY );
                 SortSync( scores, ids );
                 double ad = scores[1] - scores[0];
                 if ( ad < MIN_ADVANTAGE ) continue;
                 if ( ids[0] <= 1 ) continue;
                 if ( ids[0] == 2 ) rights.push_back(lid);
                 else rights.push_back( linv[lid] );
             }
             if ( rights.isize( ) > MAX_RIGHTS )
             {    vec<int> lens( rights.size( ) );
                 for ( int i = 0; i < rights.isize( ); i++ )
                     lens[i] = llens[ rights[i] ];
                 ReverseSortSync( lens, rights );
                 rights.resize(MAX_RIGHTS);    }

             if ( rights.empty() )
             {   DumpOut();
                 continue;    }
             int R = -1;
             if ( rights.size( ) == 1 ) R = rights[0];
             else if ( !DJANGO )
             {    vec<int> brights;
                 double ad;
                 if (OO)
                 {    
#if !USE_PIVOT
                     MemoryOrderAndOrientN(rights, lbp, llens, 
                             linv, M, memory, ad, brights );    
#else                        
                     MemoryOrderAndOrientN( L, rights, lbp, llens, 
                             linv, M, memory, ad, brights );    
#endif

                 }
                 else {
#if !USE_PIVOT
                     MemoryOrderN(rights, lbp, llens, M, memory, ad, brights );
#else
                     MemoryOrderN( L, rights, lbp, llens, M, memory, ad, brights );
#endif
                 }
                 if ( ad < MIN_ADVANTAGE )
                 {    DumpOut( );
                     continue;    }
                 R = brights[0];    }
             else
             {    vec<int> lens( rights.size( ) );
                 for ( int i = 0; i < rights.isize( ); i++ )
                     lens[i] = llens[ rights[i] ];
                 ReverseSortSync( lens, rights );
                 while( rights.size( ) > 1 )
                 {    vec<int> brights;
                     double ad;
                     if (OO)
                     {    
#if !USE_PIVOT
                         MemoryOrderAndOrientN(rights, lbp, llens, 
                                 linv, M, memory, ad, brights );
#else
                         MemoryOrderAndOrientN( L, rights, lbp, llens, 
                                 linv, M, memory, ad, brights );
#endif

                     }
                     else 
                     {    
#if !USE_PIVOT
                         MemoryOrderN( rights, lbp, llens, M, memory, ad, brights );
#else
                         MemoryOrderN( L, rights, lbp, llens, M, memory, ad, brights );
#endif
                     }
                     if ( ad >= MIN_ADVANTAGE )
                     {    R = brights[0];
                         break;    }
                     rights.pop_back( );    }    }

             // Check for done and otherwise save.

             if ( R == -1 )
             {    DumpOut( );
                 continue;    }
             int w = to_left[ dlines[R].front( )[0][0] ]; 
             if ( w < 0 )
             {    DumpOut( );
                 continue;    }
             if ( D.To(w).nonempty( ) || !D.From(w).solo( ) )
             {    DumpOut( );
                 continue;    }
             if ( IsKillable (orcs[L], orcs[R], llens[L], llens[R], winpos, BC_REQUIRE, BC_FLANK,
                              BC_IGNORE, BC_MIN))
             {    DumpOut( );
                 continue;    }
             if (verbose) {
                 (*outp) << "line " << L
                     << " can be right connected with line " << R << endl;
             }
             DumpOut( );
             #pragma omp critical
             {   
                 // dummy score
                 confident_lefts[R].push(0.0, L);
             }
         }


         // put global rights aside, take care of broken lines.
         set<int> line_between;
         for (int pass = 0; pass < 1; ++pass) {
             vec<int> l1s;
             vec<vec<int> > break_cell_ids(dlines.size());
             vec<Bool> is_line_broken(dlines.size(), False);
             for ( int L1 = 0; L1 < dlines.isize( ); L1++ )
             {    
                 // if ( llens[L1] < MIN_STAR ) continue;
                 const vec<vec<vec<int>>>& L = dlines[L1];
                 break_cell_ids[L1].push_back(0);
                 for (int j = 0; j < L.isize(); ++j) {
                     const vec<vec<int>> &cell = L[j];
                     if (cell.solo() && cell[0].empty()) {
                         int d = L[j - 1][0][0];
                         int e = D.IFrom(to_right[d], 0);
                         if (IsBarcodeOnlyGap( D.O(e))) {
                             break_cell_ids[L1].push_back(j);
                         }
                     }
                 }
                 break_cell_ids[L1].push_back(L.isize());

                 // if it has real gaps
                 if (break_cell_ids[L1].size() > 2) {
                     l1s.push_back(L1);
                     is_line_broken[L1] = True;
                 }

             }
             if (l1s.empty()) {
                 if (verbose) {
                     cout << "there is no line with bc-only barcode. Terminate."
                         << endl;
                 }
                 break;
             }
             vec<int> l1slen( l1s.isize( ) );
             for ( int i = 0; i < l1s.isize( ); i++ )
                 l1slen[i] = llens[ l1s[i] ];

             // Go through long lines and find their nearest left and right neighbors.
             vec<vec<triple<double, int, int> > > line_scores(dlines.size());
             ReverseSortSync( l1slen, l1s );
             cout << Date( ) << ": start ScaffoldMe loop on " << l1s.size( )
                 << " broken lines" << endl;

             vec<int> dlines_broken_start_pos;
             vec<vec<vec<vec<int>>>> dlines_broken;
             vec<int> broken_to_whole;
             for ( int li = 0; li < l1s.isize( ); li++ ) {
                 dlines_broken_start_pos.push_back(dlines_broken.size());
                 int lid = l1s[li];
                 vec<int> &break_cell_id = break_cell_ids[lid];
                 ForceAssertGe(break_cell_id.isize(), 2);
                 vec<vec<vec<int>>> &line = dlines[lid];
                 for (int i = 0; i < break_cell_id.isize() - 1; ++i) {
                     int cell_id = break_cell_id[i];
                     int cell_id_next = break_cell_id[i + 1];
                     
                     if (i != 0) {
                         cell_id += 1;
                     }
                     
                     const vec<vec<vec<int> > > broken_line(line.begin() + cell_id,
                             line.begin() + cell_id_next);
                     /*
                     if (broken_line.solo()) {
                         const vec<vec<int> > &cell = broken_line[0];
                         if (cell.solo()) {
                             ForceAssert(cell[0].nonempty());
                         }
                     }
                     */
                     ForceAssert(broken_line.nonempty());
                     dlines_broken.push_back(broken_line);
                     broken_to_whole.push_back(lid);
                 }
             }
             // add broken lines size to start_pos
             dlines_broken_start_pos.push_back(dlines_broken.size());
             if (verbose) {
                 cout << Date() << ": there are " << dlines_broken.size()
                     << " broken lines" << endl;
             }


             // compute cov for broken dlines
             vec<int> llens_broken;
             GetLineLengths( hb, D, dlines_broken, llens_broken );
             vec<double> cov_broken;
             vec< vec< pair<int,int> > > orcs_broken(dlines_broken.size());
             {    
                 vec<vec<pair<int,int>>> lbp_broken;
                 BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines_broken, lbp_broken, 0 );
                 MasterVec<SerfVec<pair<int,int>>> lbpx;
                 for ( auto x : lbp_broken )
                 {    
                     SerfVec<pair<int,int>> y( x.size( ) );
                     for ( int j = 0; j < x.isize( ); j++ )
                         y[j] = x[j];
                     lbpx.push_back(y);
                 }
                 LineCN( kmers, lbpx, D, dlines_broken, llens_broken, cov_broken );    

                 for (int i = 0; i < orcs_broken.isize(); ++i) {
                     vec< pair<int,int> > &orc = orcs_broken[i];
                     orc.resize( lbp_broken[i].size( ) );
                     for ( int j = 0; j < orc.isize( ); j++ )
                         orc[j] = make_pair( lbp_broken[i][j].second, lbp_broken[i][j].first );
                     Sort(orc);
                 }
             }
             // Get barcode positions for broken lines
             vec<vec<pair<int,int>>> lbp_broken;
             BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines_broken, lbp_broken, BC_VIEW );

             // compute barcode neighbor
             vec< vec< pair<int,int> > > lhood_broken;
             LineProxSP( hb, inv, ebcx, D, dinv, dlines, dlines_broken, broken_to_whole, qept, lhood_broken );
             // make a copy of lines
             // remember the difference between lines before current line
             // and after it.
             // vec<vec<vec<vec<int>>>> dlines_copy(dlines);
             vec<int> llens_copy(llens);
             llens_copy.append(llens_broken);
             // vec<double> cov_copy(cov);
             vec<vec<pair<int,int>>> lbp_copy(lbp);
             lbp_copy.append(lbp_broken);

             for ( int li = 0; li < l1s.isize( ); li++ ) {
                 int start = dlines_broken_start_pos[li];
                 int end = dlines_broken_start_pos[li + 1];
                 int lid = l1s[li];
                 vec<int> &break_cell_id = break_cell_ids[lid];


                 ForceAssertEq(break_cell_id.isize(), end - start + 1);
#pragma omp parallel for schedule(dynamic, 1)
                 for (int broken_cell_id = 0; broken_cell_id < break_cell_id.isize(); ++broken_cell_id) {
                     vec< triple<int,int,int> > M;
                     // CAUTION: memory has to be here !!!!!
                     map< vec<int>, double > memory;
                     // Set up logging.

                     ostringstream* outp;
                     if (verbose) outp = new ostringstream;
                     auto DumpOut = [&]
                     {    if (verbose)
                         {
#pragma omp critical
                             {    cout << outp->str( );    }    
                             delete outp;    }    };


                     if (broken_cell_id == 0 || broken_cell_id == break_cell_id.isize() - 1) {
                         continue;
                     } else {
                         // assign ids for lhs and rhs
                         int llid = start + broken_cell_id - 1;
                         int rlid = llid + 1;
                         ForceAssertLt(llid, end);
                         ForceAssertLt(rlid, end);
                         int L = dlines.isize() + llid;
                         int R = L + 1;
                         // find intersection of L's and R's lhood.
                         vec<int> line_pool, lefts, rights;
                         
                         /*
                         if ( TooMuchDupContentSingle(D, dlines_broken[llid], dinv)) {
                             DumpOut();
                             continue;
                         }
                         */
                         /*
                         if ( TooMuchDupContentSingle(D, dlines_broken[rlid], dinv)) {
                             DumpOut();
                             continue;
                         }
                         */
                         

                         for (int i = 0; i < lhood_broken[llid].isize(); ++i) {
                             int lid = lhood_broken[llid][i].second;
                             // if (is_line_broken[lid]) continue; 
                             if (lid == linv[lid]) continue;
                             if ( AbsDiff( cov_broken[llid], cov[lid] ) > MAX_CN_DIFF ) continue;
                             if ( llens_copy[lid] < MIN_BAR_TO ||
                                     llens_copy[lid] > MAX_BAR_TO ) {
                                 continue;
                             }
                             /*
                             if ( TooMuchDupContentSingle(D, dlines[lid], dinv)) {
                                 continue;
                             }
                             */
                             
                             if (TooMuchDupContent(D, dlines_broken[llid], dlines[lid], inv)) {
                                 continue;
                             }
                             
                             vec<double> scores;
                             scores.push_back( MemoryScoreOrder( 
                                         { lid, L }, lbp_copy, llens_copy, M, memory ) );
                             scores.push_back( MemoryScoreOrder( 
                                         { linv[lid], L }, lbp_copy, llens_copy, M, memory ) );
                             scores.push_back( MemoryScoreOrder( 
                                         { L, lid }, lbp_copy, llens_copy, M, memory ) );
                             scores.push_back( MemoryScoreOrder( 
                                         { L, linv[lid] }, lbp_copy, llens_copy, M, memory ) );
                             vec<int> ids( 4, vec<int>::IDENTITY );
                             SortSync( scores, ids );
                             /*
                                double ad = scores[1] - scores[0];
                                if ( ad < MIN_ADVANTAGE ) continue;
                                */
                             if ( ids[0] <= 1 ) continue;
                             rights.push_back(lid);
                             /*
                                if ( ids[0] == 2 ) rights.push_back(lid);
                                else rights.push_back( linv[lid] );
                                */
                         }

                         if ( rights.empty() )
                         {   DumpOut();
                             continue;    }

                         for (int i = 0; i < lhood_broken[rlid].isize(); ++i) {
                             int lid = lhood_broken[rlid][i].second;
                             // if (is_line_broken[lid]) continue; 
                             if ( lid == linv[lid]) continue;
                             if ( AbsDiff( cov_broken[rlid], cov[lid] ) > MAX_CN_DIFF ) continue;
                             if ( llens_copy[lid] < MIN_BAR_TO || 
                                     llens_copy[lid] > MAX_BAR_TO ) {
                                 continue;
                             }
                             /*
                             if ( TooMuchDupContentSingle(D, dlines[lid], dinv)) {
                                 continue;
                             }
                             */
                             
                             if (TooMuchDupContent(D, dlines[lid], dlines_broken[rlid], inv)) {
                                 continue;
                             }
                             
                             vec<double> scores;
                             scores.push_back( MemoryScoreOrder( 
                                         { R, lid }, lbp_copy, llens_copy, M, memory ) );
                             scores.push_back( MemoryScoreOrder( 
                                         { R, linv[lid] }, lbp_copy, llens_copy, M, memory ) );
                             scores.push_back( MemoryScoreOrder( 
                                         { lid, R }, lbp_copy, llens_copy, M, memory ) );
                             scores.push_back( MemoryScoreOrder( 
                                         { linv[lid], R }, lbp_copy, llens_copy, M, memory ) );
                             vec<int> ids( 4, vec<int>::IDENTITY );
                             SortSync( scores, ids );
                             /*
                                double ad = scores[1] - scores[0];
                                if ( ad < MIN_ADVANTAGE ) continue;
                                */
                             if ( ids[0] <= 1 ) continue;
                             lefts.push_back(lid);
                             /*
                                if ( ids[0] == 2 ) lefts.push_back(lid);
                                else lefts.push_back( linv[lid] );
                                */
                         }

                         if ( lefts.empty() )
                         {   DumpOut();
                             continue;    }

                         for (int i = 0; i < lefts.isize(); ++i) {
                             int lid = lefts[i];
                             if (Position(rights, lid) == -1 || lid == linv[lid]) {
                                 continue;
                             }
                             line_pool.push_back(lid);
                         }

                         if ( line_pool.empty() )
                         {   DumpOut();
                             continue;    }

                         vec<double> scores;
                         for (int i = 0; i < line_pool.isize(); ++i) {
                             // LB means Line Between
                             int LB = line_pool[i];
                             // if ( llens_copy[LB] < MIN_BAR_TO ) continue;
                             // if ( AbsDiff( cov_broken[llid], cov[LB] ) > MAX_CN_DIFF ) continue;
                             // if ( AbsDiff( cov_broken[rlid], cov[LB] ) > MAX_CN_DIFF ) continue;
                             scores.push_back( MemoryScoreOrder( 
                                         { L, LB , R }, lbp_copy, llens_copy, M, memory ) );
                             scores.push_back( MemoryScoreOrder( 
                                         { L, linv[LB] , R }, lbp_copy, llens_copy, M, memory ) );
                         }
                         if (scores.empty()) {
                             DumpOut();
                             continue;    
                         }
                         vec<int> ids( scores.size(), vec<int>::IDENTITY );
                         SortSync( scores, ids );
                         double ad = scores[1] - scores[0];
                         if (ad == 0 && (line_pool[ids[0]/2] == linv[line_pool[ids[1]/2]])) {
                             if (verbose) (*outp) << "REALLY?" << endl;
                             DumpOut();
                             continue;
                         } else if ( ad < MIN_ADVANTAGE ) {
                             DumpOut();
                             continue;
                         }
                         if (verbose) {
                             (*outp) << "tmp score: " << printSeq(scores) 
                                 << " tmp index: " << printSeq(ids) << endl;
                         }
                         int best_oo = line_pool[ids[0] / 2];
                         if (ids[0] % 2 != 0) {
                             best_oo = linv[best_oo];
                         }
                         ForceAssertLt(best_oo, dlines.isize());
                         int v = to_left[ dlines[best_oo].front( )[0][0] ]; 
                         int w = to_right[ dlines[best_oo].back( )[0][0] ]; 
                         if (v < 0 || w < 0) {
                             DumpOut();
                             continue;
                         }
                         if ( D.To(v).nonempty( ) || !D.From(v).solo( ) ) {
                             DumpOut();
                             continue;
                         }
                         if ( D.From(w).nonempty( ) || !D.To(w).solo( ) ) {
                             DumpOut();
                             continue;
                         }
                         if ( IsKillable (orcs_broken[llid], orcs[best_oo], llens_copy[L], llens[best_oo], winpos, BC_REQUIRE, BC_FLANK,
                                     BC_IGNORE, BC_MIN))
                         {    DumpOut( );
                             continue;    }
                         if ( IsKillable (orcs[best_oo], orcs_broken[rlid], llens[best_oo], llens_copy[R], winpos, BC_REQUIRE, BC_FLANK,
                                     BC_IGNORE, BC_MIN))
                         {    DumpOut( );
                             continue;    }

                         if (verbose) {
                             (*outp) << "line " << lid << " gap cell "
                                 << break_cell_id[broken_cell_id]
                                 << " can be filled with line " << best_oo << endl;
                         }
                         DumpOut();
#pragma omp critical
                         {
                             // <score, lid, cell_id>
                             line_scores[best_oo].push(scores[0], lid, break_cell_id[broken_cell_id]);
                         }
                     }
                 }
             }
             // introduce joins.
             if (verbose) {
                 cout << Date() << " start to introduce joins" << endl;
             }
             vec<triple<double, int, int> > best_line_scores(dlines.size(), triple<double, int, int> (-1, -1, -1));
             for (int i = 0; i < line_scores.isize(); ++i) {
                 vec<triple<double, int, int> > &line_score = line_scores[i];
                 if (verbose && line_score.nonempty()) {
                     cout << "scores for line " << i << " is (are): "
                         << printSeq(line_score) << endl;
                 }
                 if (line_score.solo()) {
                     best_line_scores[i] = line_score[0];
                 }
             }
             // find non-circle (dead lock) line relationship and insert bc-only
             // gaps, also break down the old ones.
             vec<int> lines_to_fill;
             for (int i = 0; i < best_line_scores.isize(); ++i) {
                 // lines with no place to fill at all.
                 if (best_line_scores[i].second == -1 ||
                         best_line_scores[i].third == -1) {
                     continue;
                 }
                 // detect if there is a circle.
                 int start = best_line_scores[i].second;
                 int next = start;
                 Bool is_circle = False;
                 vec<int> path({start});
                 while (1) {
                     next = best_line_scores[next].second;
                     if (next == -1) {
                         break;
                     } else if (next == start) {
                         is_circle = True;
                         if (verbose) {
                             cout << "circle detected: " 
                                 << printSeq(path) << endl;
                         }
                         break;
                     }
                     path.push_back(next);
                 }
                 if (is_circle) {
                     continue;
                 }
                 // if reach here, there is no dead lock.
                 lines_to_fill.push_back(i);
             }
             if (verbose) {
                 cout << "there are " << lines_to_fill.size()
                     << " lines to fill (candidates)" << endl;
             }
             if (lines_to_fill.empty()) {
                 break;
             }

             // check if there is any collision with lines to fill.
             // first, create a matrix of line-cell-line_to_fill
             vec<map<int, int> > line_cell(dlines.size());
             // vec<vec<int> > break_cell_ids(dlines.isize());
             for (int i = 0; i < best_line_scores.isize(); ++i) {
                 triple<double, int, int> &best_line_score = best_line_scores[i];
                 int lid = best_line_score.second;
                 int cell_id = best_line_score.third;
                 if (lid == -1 || cell_id == -1) {
                     continue;
                 }
                 line_cell[lid][cell_id] = i;
             }

             vec<Bool> touched( D.E(), False );
             int cnt = 0;
             for (int i = 0; i < lines_to_fill.isize(); ++i) {
                 int lid_to_fill = lines_to_fill[i];
                 triple<double, int, int> &best_line_score = best_line_scores[lid_to_fill];
                 int lid = best_line_score.second;
                 int cell_id = best_line_score.third;
                 // gap cell id on reverse strand line
                 int cell_id_re = dlines[lid].isize() - cell_id - 1;
                 if (line_cell[linv[lid]].find(cell_id_re) == line_cell[linv[lid]].end()) {
                     if (verbose) {
                         cout << "warning type 1! line " << lid_to_fill
                             << " only fill into one strand of line "
                             << lid << endl;
                     }
                     if (require_symmetric) {
                         continue;
                     }
                 } else if (line_cell[linv[lid]][cell_id_re] != linv[lid_to_fill]) {
                     if (verbose) {
                         cout << "warning type 2! line " << lid
                             << " have different lines to fill in on two strands: "
                             << lid_to_fill << " and "
                             << line_cell[linv[lid]][cell_id_re] << endl;
                     }
                     /*
                        if (llens[lid_to_fill] <= llens[line_cell[linv[lid]][cell_id_re]]) {
                        continue;
                        }
                        */
                     continue;
                 }
                 // another check
                 int lid_to_fill_re = linv[lid_to_fill];
                 triple<double, int, int> &best_line_score_re = best_line_scores[lid_to_fill_re];
                 if (lid_to_fill_re != lid_to_fill &&
                         best_line_score_re.second != -1 &&
                         best_line_score_re.third != -1) {

                     // collision can't be solved.
                     if (best_line_score_re.second != linv[lid] ||
                             best_line_score_re.third != cell_id_re) {
                         if (verbose) {
                             cout << "warning type 3! line " << lid_to_fill << " and "
                                 <<  lid_to_fill_re << " have different gap cell to fill"
                                 << endl;
                         }
                         continue;
                     }
                 }

                 // two dummies for now
                 if (cell_id == 0 || cell_id == break_cell_ids[lid].back()) {
                     continue;
                 } else {
                     int d1 = dlines[lid][cell_id - 1][0][0];
                     int d2 = dlines[lid][cell_id + 1][0][0];
                     int d3 = dinv[d1], d4 = dinv[d2];
                     if (touched[d1] || touched[d2] || touched[d3] || touched[d4]) {
                         continue;
                     }
                     touched[d1] = touched[d2] = touched[d3] = touched[d4] = True;
                     int v = to_right[d1], w = to_left[d2];
                     int rv = to_right[ dinv[d2] ], rw = to_left[ dinv[d1] ];
                     ForceAssert((D.IFrom(v).solo() &&
                                 D.ITo(w).solo() &&
                                 D.IFrom(rv).solo() &&
                                 D.ITo(rw).solo()));

                     int d_left = dlines[lid_to_fill].front()[0][0];
                     int d_left_r = dinv[d_left];
                     int d_right = dlines[lid_to_fill].back()[0][0];
                     int d_right_r = dinv[d_right];
                     if (touched[d_left] || touched[d_left_r] || touched[d_right] || touched[d_right_r]) {
                         continue;
                     }
                     touched[d_left] = touched[d_left_r] = touched[d_right] = touched[d_right_r] = True;
                     // get v1, w1 for line_to_fill
                     int v1 = to_left[d_left], w1 = to_right[d_right];
                     int rv1 = to_left[ dinv[d_right] ], rw1 = to_right[ dinv[d_left] ];
                     ForceAssert((D.ITo(v1).empty() &&
                                 D.IFrom(w1).empty() &&
                                 D.ITo(rv1).empty() &&
                                 D.IFrom(rw1).empty()));

                     // break original bc-only gaps.
                     int bc_gap_d = D.IFrom(v, 0);
                     ForceAssert(IsBarcodeOnlyGap( D.O(bc_gap_d)));
                     dels.push_back(bc_gap_d);
                     dels.push_back(dinv[bc_gap_d]);

                     // insert new bc-only gaps.
                     // add 4 new superedges of bc-only gap.
                     dinv.push_back( D.E( ) + 1, D.E( ), D.E() + 3, D.E() + 2 );
                     D.AddEdge( v, v1, {-2} );
                     D.AddEdge( rw1, rw, {-2} );
                     D.AddEdge( w1, w, {-2} );
                     D.AddEdge( rv, rv1, {-2} );
                     cnt += 2;
                     // mark line_to_fill
                     line_between.insert(lid_to_fill);
                     line_between.insert(linv[lid_to_fill]);
                 }
             }
             if (verbose) {
                 cout << cnt << " lines can be filled." << endl;
             }
         }

         // come back and take care of global rights
         vec<int> most_confident_rights(dlines.size(), -1);
         for (int i = 0; i < confident_lefts.isize(); ++i) {
             vec<pair<double, int> > &lefts = confident_lefts[i];
             if (lefts.solo() &&
                 line_between.find(lefts[0].second) == line_between.end() &&
                 line_between.find(i) == line_between.end()) {
                 most_confident_rights[lefts[0].second] = i;
             }
         }
         // check for assymmetric rights
         for (int i = 0; i < most_confident_rights.isize(); ++i) {
             int &right = most_confident_rights[i];
             if (right != -1 && most_confident_rights[linv[right]] != linv[i]) {
                 if (require_symmetric) {
                     if (verbose) {
                         cout << "line " << i << " --> " << right
                             << " not symmetric" << endl;
                     }
                     most_confident_rights[linv[right]] = -1;
                     right = -1;
                 } else {
                     
                     if (most_confident_rights[linv[right]] == -1) {
                         most_confident_rights[linv[right]] == linv[i];
                     } else {
                         if (verbose) {
                             cout << "line " << i << " --> " << right
                                 << " not symmetric" << endl;
                         }
                         most_confident_rights[linv[right]] = -1;
                         right = -1;
                     }
                        
                 }
             }
         }
         
         // make joins
         vec<Bool> touched( D.E(), False );
         int cnt = 0;
         for (int i = 0; i < most_confident_rights.isize(); ++i) {
             int &right = most_confident_rights[i];
             if (right == -1) {
                 continue;
             }
             // insert line_to fill on the right
             int d_left = dlines[i].back()[0][0]; 
             int d_left_r = dinv[d_left];
             int d_right = dlines[right].front()[0][0];
             int d_right_r = dinv[d_right];
             if (touched[d_left] || touched[d_left_r] || touched[d_right] || touched[d_right_r]) {
                 continue;
             }
             touched[d_left] = touched[d_left_r] = touched[d_right] = touched[d_right_r] = True;
             int v = to_right[d_left], w = to_left[d_right];
             int rv = to_right[d_right_r], rw = to_left[d_left_r];
             ForceAssert(D.IFrom(v).empty() && D.ITo(v).solo());
             ForceAssert(D.ITo(w).empty() && D.IFrom(w).solo());
             ForceAssert(D.IFrom(rv).empty() && D.ITo(rv).solo());
             ForceAssert(D.ITo(rw).empty() && D.IFrom(rw).solo());
             /*
             if (D.IFrom(rv).nonempty() || !D.ITo(rv).solo()) {
                 cout << "right nimei rv" << D.IFrom(rv).size() << ' ' << D.ITo(rv).size() << endl;
             }
             if (D.ITo(rw).nonempty() || !D.IFrom(rw).solo()) {
                 cout << "right nimei rw" << D.IFrom(rv).size() << ' ' << D.ITo(rv).size() << endl;
             }
             */
             dinv.push_back( D.E( ) + 1, D.E( ) );
             D.AddEdge( v, w, {-2} );
             D.AddEdge( rv, rw, {-2} );
             cnt += 2;
         }
         if (verbose) {
             cout << Date( ) << ": made " << cnt << " joins" << endl;
         }


         // clean up current supergraph D.
         if (dels.nonempty()) {
             D.DeleteEdges(dels);
         }
         RemoveUnneededVertices( D, dinv );
         CleanupCore( D, dinv );
         if (dels.empty() && cnt == 0) {
             break;
         }
     }
}
