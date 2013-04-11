// $Id: ssm_csia.h,v 1.1.1.1 2004/11/23 16:24:37 keb Exp $
// =================================================================
//
//    05.04.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// -----------------------------------------------------------------
//
//  **** Module  :  ssm_csia       <interface>
//       ~~~~~~~~~
//  **** Classes :  ssm::GraphMatch ( matching SS graphs )
//       ~~~~~~~~~
//
//  E. Krissinel 2001-2013
//
//  When used, please cite:
//
//   Krissinel, E. and Henrick, K. (2004)
//   Common subgraph isomorphism detection by backtracking search.
//   Software - Practice and Experience, 34, 591-607.
//
// =================================================================
//

#ifndef  __SSM_CSIA__
#define  __SSM_CSIA__

#include "ssm_graph.h"

namespace ssm  {

  //  =========================  Match  ===========================

  DefineClass(Match);

  class Match : public CStream  {

    friend class GraphMatch;

    public :

      Match ();
      Match ( RPCStream Object );
      Match ( ivector FV1, ivector FV2, int nv, int n, int m );
      ~Match();

      void SetMatch ( ivector FV1, ivector FV2,
                      int nv, int n, int m );   // FV1[], FV2[] are copied

      void    Swap();

      Boolean isMatch ( ivector FV1, ivector FV2, int nv );

      int  isSubMatch ( ivector FV1, ivector FV2, int nv );
      // return 0 <=> no submatch relations
      //        1 <=> "this" is submatch of (FV1,FV2)
      //       -1 <=> (FV1,FV2) is submatch of "this"

      void GetMatch ( ivector  & FV1,  // do not allocate or
                      ivector  & FV2,  // dispose FV1 and FV2 in
                      int      & nv ); // application!

      void GetMatch ( ivector  & FV1, // do not allocate or
                      ivector  & FV2, // dispose FV1 and FV2 in
                      int      & nv,  // application!
                      realtype & p1,
                      realtype & p2 );

      void read  ( RCFile f );
      void write ( RCFile f );

    protected :
      int     mlength,n1,n2;
      ivector F1,F2;

      void InitMatch();

    private :
      int nAlloc;

  };

  DefineStreamFunctions(Match)

  //  =========================  GraphMatch  ===========================


  DefineClass(GraphMatch);

  class GraphMatch : public CStream  {

    public :

      GraphMatch ();
      GraphMatch ( RPCStream Object );
      ~GraphMatch();

      void  SetUniqueMatch ( Boolean unique_match );
      void  SetBestMatch   ( Boolean best_match   );
      void  SetMatchBufferLength ( int matchBufLen );
      void  SetFlags       ( word Flags );
      void  RemoveFlags    ( word Flags );

      void  MatchGraphs    ( PGraph Gh1, PGraph Gh2, int minMatch );

      PGraph  GetGraph1 ();
      PGraph  GetGraph2 ();
      void  GetMatches     ( PPMatch & SSMatch, int & nOfMatches );
      void  GetMatch       ( int   matchNo, int & matchLen,
                             ivector  & F1, ivector  & F2,
                             realtype & p1, realtype & p2 );
      inline int GetMaxRecursionLevel()  { return maxRecursionLevel; }
      inline int GetNofMatches       ()  { return nMatches;          }
      int   GetNofMatches     ( realtype p1, realtype p2 );

      int   CheckConnectivity ( int matchNo );

      void  read  ( RCFile f );
      void  write ( RCFile f );

    protected :

      PGraph    G1,G2;
      PPVertex  V1;
      PPVertex  V2;
      PPEdge    E1;
      PPEdge    E2;
      imatrix   c1,c2;
      Boolean   swap;
      word      flags;
      int       n,m;

      imatrix3  P;
      imatrix   iF1;
      ivector   F1,F2,ix;

      int       nMatches,maxNofMatches;
      PPMatch   match;
      Boolean   UniqueMatch,BestMatch,wasFullMatch,Stop;
      int       maxMatch,maxCollectedMatch,maxRecursionLevel;

      void  InitGraphMatch   ();
      void  FreeMemory       ();
      void  FreeRecHeap      ();
      void  GetMemory        ();
      void  GetRecHeap       ();
      int   Initialize       ();
      void  DoMatch          ( int minMatch );
      void  MatchSingleVertex();
      void  Backtrack        ( int i );
      void  Backtrack1       ( int i, int k0 );
      void  CollectMatch     ( int nm );

    private :
      int nAlloc,mAlloc,nMAlloc;

  };

  DefineStreamFunctions(GraphMatch)

}

#endif
