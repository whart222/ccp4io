// $Id: ssm_graph.h,v 1.2 2005/12/20 12:09:14 keb Exp $
// =================================================================
//
//    10.04.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// -----------------------------------------------------------------
//
//  **** Module  :  ssm_graph  <interface>
//       ~~~~~~~~~
//  **** Classes :  ssm::Graph  ( secondary structure graph )
//       ~~~~~~~~~
//
//  E. Krissinel 2002-2013
//
// =================================================================
//


#ifndef  __SSM_Graph__
#define  __SSM_Graph__

#include "mmdb/mmdb_manager.h"

#include "ssm_vxedge.h"

namespace ssm  {

  //  ==========================  Graph  ===========================

  DefineClass(Graph);
  DefineStreamFunctions(Graph)

  class Graph : public CStream  {

    friend class GraphMatch;

    public :

      Graph ();
      Graph ( RPCStream Object );
      ~Graph();

      void  Reset();  // must be called before building a graph
                      // The sequence of calls is:
                      //    Graph.Reset();
                      //    for (....)  {
                      //      V = new CSSVertex();
                      //      .....
                      //      SSGraph.AddVertex ( V );
                      //    }
                      //    SSGraph.Build();

      void  SetGraphName  ( cpstr gname );

      void  SelectCalphas ( PCMMDBManager MMDB, int & selHnd,
                            cpstr selstring );

      //   AddVertex(..) do not copy the objects, but take them over.
      // This means that application should forget about pointers to
      // V once they were given to CSSGraph. All vertices must be
      // allocated newly prior each call to AddVertex(..).
      void  AddVertex  ( PVertex Vx );

      int   MakeGraph ( PCMMDBManager MMDB );

      void  CalcVertexOrder();
      void  RepairSS ( PCMMDBManager MMDB );

      //   BuildGraph() calculates all edges and builds the graph.
      void    BuildGraph();
      Boolean isBuild   ();

      void  calcVTypes();  // calculates nHelices and nStrands only

      //   ReleaseEdges() deallocates all graph edges and
      //  the connectivity matrix
      void  ReleaseEdges();

      void  RemoveShortVertices   ( int nmin_hx, int nmin_sd );

      //   LeaveVertices(..) removes all vertices from the graph
      // except those having numbers listed in vector vlist. Thus,
      // if vlist[i]=j, 1<=i<=vllen,  1<=j, then jth vertex will
      // not be removed.
      void  LeaveVertices         ( ivector vlist, int vllen );

      //   LeaveVertices(..) removes all vertices from the graph
      // except those found in the specified range. 'select' is of
      // the following format:
      //    "*", "(all)"            - take all file
      //    "-"                     - take chain without chain ID
      //    "a:Ni-Mj,b:Kp-Lq,..."   - take chain a residue number N
      //                              insertion code i to residue number M
      //                              insertion code j plus chain b
      //                              residue number K insertion code p to
      //                              residue number L insertion code q and
      //                              so on.
      //    "a:,b:..."              - take whole chains a and b and so on
      //    "a:,b:Kp-Lq,..."        - any combination of the above.
      void  LeaveVertices ( cpstr select, PCMMDBManager M );

      //    LeaveVertices ( selHnd,MMDB ) leaves only vertices that are
      // covered by the given selection. selHnd may refer to the
      // selection of atoms, residues or chains.
      void  LeaveVertices ( int selHnd, PCMMDBManager M );

      void  RemoveVertex  ( int vertex_no );  // 1..nVertices

      Boolean    inRange  ( cpstr chainID, int initPos, int endPos );
      inline cpstr    GetGraphName   () { return name;        }
      inline cpstr    GetDevChain    () { return devChain;    }
      pstr            GetChainList   ( pstr S );
      inline int      GetNofVertices () { return nVertices;   }
      inline PPVertex GetVertices    () { return V;           }
      inline int      GetNofEdges    () { return nEdges;      }
      inline int      GetNofHelices  () { return nHelices;    }
      inline int      GetNofStrands  () { return nStrands;    }
      void     GetAllChains     ( PChainID & chain, int & nchains );
      int      GetNofChains     ();
      Boolean  GetEdgeDirection ( int v1, int v2, vect3 & v );
      VERTEX_TYPE GetVertexType ( int vertex_no  ); // 1..nVertices
      int      GetVertexClass   ( int vertex_no  ); // 1..nVertices
      Boolean  GetVertexDirection ( int vertex_no, vect3 & v );
      int      GetSeqLength     ( int vertex_no  ); // 1..nVertices
      realtype GetMass          ( int vertex_no  ); // 1..nVertices
      PVertex  GetGraphVertex   ( int vertex_no  ); // 1..nVertices
      pstr     GetVertexChainID ( int vertex_no  ); // 1..nVertices
      pstr     GetVertexInitRes ( int vertex_no  ); // 1..nVertices
      pstr     GetVertexEndRes  ( int vertex_no  ); // 1..nVertices
      void     GetVertexRange   ( int     vertex_no,  // 1..nVertices
                                  ChainID chID,
                                  int &   initSeqNum,
                                  InsCode initICode,
                                  int &   endSeqNum,
                                  InsCode endICode  );
      void     GetVertexRange   ( int     vertex_no,  // 1..nVertices
                                  ChainID chID,
                                  int &   initPos,
                                  int &   endPos );
      VERTEX_TYPE GetSSEType    ( pstr chainID, int atomPos );
      VERTEX_TYPE GetSSEType    ( PCAtom A );

      PEdge   GetGraphEdge   ( int edge_no );     // 1..nEdges
      PEdge   GetGraphEdge   ( int v1, int v2 );  // 1..nVertices

      realtype CalcCombinations ( ivector F, int nm );

      void  DevelopChainGraphs ( PPGraph & G, int & nGraphs );

      //  Superpose(..) returns TMatrix - a transformation matrix for
      // G's coordinates, such that TMatrix*{G} ~= {this}
      //  F1 is for 'this' graph, F2 = for G.
      void  Superpose   ( PGraph G, ivector F1, ivector F2,
                          int nMatch, mat44 & TMatrix );

      void  Copy  ( PGraph G );

      void  read  ( RCFile f );
      void  write ( RCFile f );

    protected :
      pstr     name;        // graph name
      ChainID  devChain;    // chain of a developed graph
      int      nVertices,nEdges;
      int      nHelices,nStrands;

      PPVertex V;
      PPEdge   E;
      imatrix  graph;

      void  InitGraph      ();
      void  FreeMemory     ();
      void  _leaveVertices ( PCMMDBManager M, int selHnd1 );

      //   CompareEdges(..) compares edge (ij) of the graph with
      // edge (kl) of graph G. i may be either less or greater
      // than j, same about k and l. If edges compare, the function
      // returns 0. Edges with equal indices (i.e. (ii) and (kk))
      // are considered as comparable (returns 0).
      //   The function may be used only after both graphs have
      // been built.
      int CompareEdges ( int i, int j, PGraph G, int k, int l );

      int CheckEdgeConnectivity ( int i, int j, PGraph G, int k, int l );

    private :
      int  nVAlloc,nEAlloc,nGAlloc;

  };


  //  ==================================================================

  //   In SelectDomain(..) and CutOutDomain(..), select is of the
  // following format:
  //    "*", "(all)"            - take all file
  //    "-"                     - take chain without chain ID
  //    "a:Ni-Mj,b:Kp-Lq,..."   - take chain a residue number N
  //                              insertion code i to residue number M
  //                              insertion code j plus chain b
  //                              residue number K insertion code p to
  //                              residue number L insertion code q and
  //                              so on.
  //    "a:,b:..."              - take whole chains a and b and so on
  //    "a:,b:Kp-Lq,..."        - any combination of the above.
  extern int SelectDomain ( PCMMDBManager MMDB, int & selHnd,
                            cpstr select, int selType );
  extern int CutOutDomain ( PCMMDBManager MMDB, cpstr select );

  extern PGraph GetSSGraph ( PCMMDBManager M, int selHnd, int & rc );

  extern void DisposeGraphs ( PPGraph & G, int & nGraphs );

  extern int  SuperposeGraphs ( PGraph G1, ivector F1,
                                PGraph G2, ivector F2,
                                int     matchlen,
                                mat44 & TMatrix );

}

/*

extern realtype  GetTorsion ( rvector U, rvector W, rvector V );
//      U     W      V
//   o<----o----->o----->o
//

extern realtype  GetAngle   ( rvector v1, rvector v2 );
//  returns angle between v1 and v2



extern void  CalcCombinations ( rvector & combs, int & vlen,
                                PCSSGraph G1, PCSSGraph G2 );

*/

#endif
