// $Id$
// =================================================================
//
//    05.04.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  SSM_Align <interface>
//       ~~~~~~~~~
//  **** Project :  Structure alignment in 3D
//       ~~~~~~~~~
//  **** Classes :  ssm::Align   ( Secondary Structure Matching )
//       ~~~~~~~~~  ssm::XAlign  ( Output alignment             )
//                  ssm::XTAlign ( Text output alignment        )
//
//  E. Krissinel, 2002-2013
//
// =================================================================
//

#ifndef  __SSM_Align__
#define  __SSM_Align__

#include "mmdb/mmdb_manager.h"

#include "ssm_superpose.h"
#include "ssm_csia.h"


//  ---------------------------  CSSMAlign  ------------------------

namespace ssm {

  DefineClass(Align);
  DefineStreamFunctions(Align)

  class Align : public CStream  {

    public :
      mat44    TMatrix; //!< superposition matrix to be applied to 1st structure
      realtype rmsd;         //!< core rmsd achieved
      realtype Qscore;       //!< core Q achieved
      int      cnCheck;      //!< connectivity option used
      int      nres1,nres2;  //!< number of residues in structures
      int      nsel1,nsel2;  //!< number of residues in aligned selections
      int      nalgn;        //!< number of aligned residues
      int      ngaps;        //!< number of gaps
      int      nmd;          //!< number of misdirections
      realtype ncombs;       //!< number of SSE combinations
      realtype seqIdentity;  //!< sequence identity
      int      selHndCa1,selHndCa2; //!< selection handles to used C-alphas
      ivector  Ca1,Ca2;      //!< C-alpha correspondence vectors
                             /// Ca1[i] corresponds to a[i], where a is
                             /// selection identified by selHndCa1
      rvector  dist1;        //!< optimizedd distances between the query
                             /// and target C-alphas
      PGraph   G1,G2;        //!< retained SSE graphs

      Align ();
      Align ( RPCStream Object );
      ~Align();

      int align ( PCMMDBManager M1, PCMMDBManager M2,
                  PRECISION     precision,
                  CONNECTIVITY  connectivity,
                  int selHnd1=0, int selHnd2=0 );

      int AlignSelectedMatch ( PCMMDBManager M1, PCMMDBManager M2,
                               PRECISION     precision,
                               CONNECTIVITY  connectivity,
                               int selHnd1=0, int selHnd2=0,
                               int nselect=0 );

      rvector GetQvalues () const { return pqvalues; }
      int     GetNMatches() const { return nMatches; }

      PSuperpose GetSuperpose() { return &superpose; }

      void  read  ( RCFile f );
      void  write ( RCFile f );

    protected :
      GraphMatch U;
      Superpose  superpose;
      rvector    pqvalues;
      int        nMatches;

      void  InitAlign ();
      void  FreeMemory();
      void  MapSelections  ( int & selHndCa, PCMMDBManager M,
                             PGraph G, int selHnd, ivector & newID );
      void  MakeSelections ( PCMMDBManager M1, int selHnd1,
                             PCMMDBManager M2, int selHnd2 );

  };


  //  -----------------------------  CXAlign --------------------------

  DefineStructure(XBlock);

  struct XBlock  {
    int      i1,i2;    //!< the outer block boundaries
    int      ip1,ip2;  //!< the alignment boundaries (ip1>=i1, ip2<=i2)
    int      icol;     //!< the block "column" number
    realtype mc;       //!< center of "index mass"
  };


  DefineClass(XAlign);

  class XAlign  {

    public :
      XAlign();
      virtual ~XAlign();

      void align ( PGraph g1, PPCAtom Calpha1, ivector Ca1, int nat1,
                   PGraph g2, PPCAtom Calpha2, ivector Ca2, int nat2,
                   rvector dist1, int & nr );

      int  GetNCols2() { return nCols2; }

    protected :
      PXBlock  XBlock1,XBlock2;
      int       nBlock1,nBlock2;
      int       na1,na2,nCols1,nCols2,nRows,algnLen;

      ivector   a1,a2;
      PPCAtom   alpha1,alpha2;
      PGraph sg1,sg2;
      rvector   d1;
      realtype  maxdist;

      virtual void FreeMemory();
      virtual void customInit();
      int   makeXBlocks  ( ivector Ca, int nat, RPXBlock xBlock,
                           int & nBlocks );
      void  alignXBlocks ( RXBlock B1, RXBlock B2, int & nr );

      virtual void makeRow ( PCAtom A1, int sseType1,
                             PCAtom A2, int sseType2,
                             realtype dist, int rowNo, int icol,
                             Boolean aligned );
  };


  //  ----------------------------  CXTAlign --------------------------

  DefineStructure(XTAlign);

  struct XTAlign  {
    realtype hydropathy1,hydropathy2,dist;
    ChainID  chID1,chID2;
    ResName  resName1,resName2;
    InsCode  insCode1,insCode2;
    int      alignKey; //!< 0: aligned, 1: not aligned, 2: NULL 1, 3: NULL 2
    int      loopNo;
    int      sseType1,sseType2;
    int      seqNum1,seqNum2;
    int      simindex;
    void  Print ( RCFile f );
  };

  DefineClass(XAlignText);

  class XAlignText : public XAlign  {

    public :
      XAlignText ();
      ~XAlignText();

      PXTAlign GetTextRows   () { return R; }
      void     GetAlignments ( pstr & algn1, pstr & algn2 );
      void     WipeTextRows  ();

    protected :
      PXTAlign R;

      void customFree();
      void customInit();
      void makeRow   ( PCAtom A1, int sseType1,
                       PCAtom A2, int sseType2,
                       realtype dist, int rowNo, int icol,
                       Boolean aligned );
  };

  extern void PrintAlignTable ( RCFile f,
                                PCMMDBManager M1, PCMMDBManager M2,
                                PAlign SSMAlign );
}

#endif
