// $Id: ssm_vxedge.h,v 1.1.1.1 2004/11/23 16:24:37 keb Exp $
// =================================================================
//
//    05.04.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// -----------------------------------------------------------------
//
//  **** Module  :  ssm_vxedge  <interface>
//       ~~~~~~~~~
//  **** Classes :  ssm::Vertex  ( secondary structure graph vertex )
//       ~~~~~~~~~  ssm::Edge    ( secondary structure graph edge   )
//
//  E. Krissinel 2002-2013
//
// =================================================================
//


#ifndef  __SSM_VxEdge__
#define  __SSM_VxEdge__

#include "mmdb/mmdb_manager.h"
#include "ssm_defs.h"

//  ==========================  Tune-up  ============================

namespace ssm  {

  extern int hx_min_len;
  extern int sd_min_len;

  extern void InitGraph();  // should be called on top of application

  extern void SetMatchPrecision    ( PRECISION    prec      );
  extern void writeMatchParameters ( cpstr        FileName  );
  extern int  readMatchParameters  ( cpstr        FileName  );
  extern void SetConnectivityCheck ( CONNECTIVITY checkMode );


  //  =========================  CSSVertex  ===========================

  DefineClass(Vertex);
  DefineStreamFunctions(Vertex)

  class Vertex : public CStream  {

    friend class Edge;
    friend class Graph;

    public :

      Vertex ();
      Vertex ( RPCStream Object );
      ~Vertex();

      int  SetVertex ( PCMMDBManager MMDB, PCHelix  Helix  );
      int  SetVertex ( PCMMDBManager MMDB, PCStrand Strand );
      int  SetVertex ( PCMMDBManager MMDB, VERTEX_TYPE v_type,
                       int sNum, int  iclass, ChainID chID,
                       int seqNum1, InsCode iCode1,
                       int seqNum2, InsCode iCode2 );

      inline void SetID ( int vid ) { id = vid; }

      realtype GetAngle  ( PVertex v );
      realtype GetCosine ( PVertex v );
      realtype GetAngle  ( realtype vx, realtype vy, realtype vz );

      pstr     GetShortVertexDesc ( pstr S );
      pstr     GetFullVertexDesc  ( pstr S );

      Boolean  Compare ( PVertex v ); // True if vertices compare

      realtype GetLengthDeviation ( PVertex v );

      void     GetDirection ( vect3 & v );
      void     GetPosition  ( vect3 & p );
      void     GetPosition  ( realtype & vx0, realtype & vy0,
                              realtype & vz0 );

      inline realtype GetLength    ()  { return length; }
      inline int      GetSeqLength ()  { return nres;   }
      inline realtype GetMass      ()  { return mass;   }

      inline realtype GetX1        ()  { return x1;     };
      inline realtype GetX2        ()  { return x2;     };
      inline realtype GetY1        ()  { return y1;     };
      inline realtype GetY2        ()  { return y2;     };
      inline realtype GetZ1        ()  { return z1;     };
      inline realtype GetZ2        ()  { return z2;     };

      Boolean  inRange ( cpstr chID, int Pos1, int Pos2 );

      inline int   GetVertexType   () { return type;    }
      inline int   GetVertexChainNo() { return VNo;     }
      inline cpstr GetChainID      () { return chainID; }
      void  GetVertexRange  ( ChainID chID,
                              ResName name1,
                              int &   seqNum1,
                              InsCode insCode1,
                              ResName name2,
                              int &   seqNum2,
                              InsCode insCode2 );

      void  Copy  ( PVertex v );

      void  read  ( RCFile f );
      void  write ( RCFile f );

    protected :

      //  matching info
      int       id;          //!< unique identifier that MUST be the vertex
                             /// number starting from 1 on
      VERTEX_TYPE  type;     //!< a V_XXXXX constant
      int       classID;     //!< class ID for helices
      int       nres;        //!< number of residues
      realtype  x0,y0,z0;    //!< center of mass
      realtype  mass;        //!< the mass
      realtype  ex,ey,ez;    //!< direction vector
      realtype  dalpha;      //!< uncertainty angle
      realtype  length;      //!< vertex length

      //  identification info
      pstr     name;        //!< composed name for short identification
      int      serNum;      //!< helix serial number
      int      strandNo;    //!< strand number
      maxMMDBName vertexID; //!< helix ID or sheet ID
      ChainID  chainID;     //!< chain ID (only for identification)
      ResName  initResName; //!< name of the strand's initial residue
      int      initSeqNum;  //!< sequence number of the initial residue
      int      initPos;     //!< sequence position of the initial residue
      InsCode  initICode;   //!< insertion code of the initial residue
      ResName  endResName;  //!< name of the strand's terminal residue
      int      endSeqNum;   //!< sequence number of the terminal residue
      int      endPos;      //!< sequence position of the terminal residue
      InsCode  endICode;    //!< insertion code of the terminal residue
      int      VNo;         //!< number of vertex in the chain

      realtype x1,x2;       //!< coordinates
      realtype y1,y2;       ///   SSE
      realtype z1,z2;       ///     ends

      void  InitVertex ();
      void  FreeMemory   ();
      void  CalcGeometry ( PPCAtom CA );
      int   GetPositions ( PCMMDBManager MMDB, int minlen );
      realtype  GetCoor1 ( PPCAtom CA, int coor_key );
      realtype  GetCoor2 ( PPCAtom CA, int coor_key );

  };


  //  ==========================  CSSEdge  ============================

  DefineClass(Edge);
  DefineStreamFunctions(Edge)

  class Edge : public CStream  {

    friend class Graph;
    friend class GraphMatch;

    public :

      Edge ();
      Edge ( RPCStream Object );
      ~Edge();

      void     SetEdge  ( PVertex v1, PVertex v2 );

      realtype GetAngle ( PVertex v );  // returns angle between
                                        // the edge and vertex
      realtype GetCosine ( PEdge E );   // returns cosine angle between
                                        // the edges
      realtype GetAngle ( rvector V1, rvector V2 );

      // Compare(..) returns 0 if edges compare, that is:
      //   1. edge lengths compare within relative precision
      //      edge_len_tol
      //   2. angles alpha1, alpha2 and alpha3 compare within
      //      absolute deviations edge_alphaX_tol .
      int   Compare ( Boolean swap_this, PEdge edge,
                      Boolean swap_edge );

      int   CheckConnectivity ( Boolean swap_this, PEdge edge,
                                Boolean swap_edge );

      void  GetDirection ( vect3 & v );
      inline realtype GetLength () { return length; }

      void  read  ( RCFile f );
      void  write ( RCFile f );

    protected :
      int      id1,id2;  //!< linked vertices
      int      vtype1;   //!< type of 1st linked vertex
      int      vtype2;   //!< type of 2nd linked vertex
      int      bdir;     //!< bond direction along the chain
      realtype length;   //!< length of edge (between v1 and v2 mass centers)
      realtype ex,ey,ez; //!< direction vector from v1 to v2
      realtype alpha1;   //!< angle V1E between v1 and the edge
      realtype alpha2;   //!< angle V2E between v2 and the edge
      realtype alpha3;   //!< angle V1V2 between v1 and v2
      realtype alpha4;   //!< torsion angle V1EV2 of v1, edge and v2
      realtype dalpha1;  //!< uncertainty in alpha1
      realtype dalpha2;  //!< uncertainty in alpha2
      realtype dalpha3;  //!< uncertainty in alpha3
      realtype dalpha4;  //!< uncertainty in alpha4
      realtype dr12;
      Boolean  GoodTorsion; //!< True if the VEV torsion angle is well defined

      void  InitEdge();

  };

}

#endif
