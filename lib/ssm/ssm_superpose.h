// $Id$
// =================================================================
//
//    05.04.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// -----------------------------------------------------------------
//
//  **** Module  :  ssm_spose  <interface>
//       ~~~~~~~~~
//  **** Functions : SuperposeCalphas ( superposing protein structures )
//       ~~~~~~~~~~~
//
//  E. Krissinel 2002-2013
//
// =================================================================
//


#ifndef  __SSM_Superpose__
#define  __SSM_Superpose__

#include "ssm_graph.h"


//  =================================================================

namespace ssm  {

  DefineStructure(SpAtom);

  struct SpAtom  {
    ChainID  chID;
    int      c,sse,c0;
    realtype dist,dist0;
    int      unmap1,unmap2;
    Boolean  excluded;
    Boolean  CompatibleSSE ( RSpAtom a );
  };


  DefineStructure(SectionDist);

  struct SectionDist  {
    realtype dist,rmsd,cosine;
    int      core_pos1,core_pos2,core_e1,core_e2;
    int      na,pos1,pos2,e1,e2;
    int      sse1,sse2;
    void Copy ( RSectionDist D );
  };


  DefineStructure(SSEDesc);

  struct SSEDesc  {
    realtype x1,y1,z1,x2,y2,z2;      // transformed start/end coordinates
    realtype xs1,ys1,zs1,xs2,ys2,zs2;   // original start/end coordinates
    realtype score,Qscore,Rscore,Xscore; // overlaping scores
    int      pos,len,pend, type,classID;
    int      m,match;
    void  Transform ( mat44  & T );
    void  CalcScore ( RSSEDesc D );
    realtype Cosine ( RSSEDesc D );
    void       Copy ( RSSEDesc D );
  };


  DefineStructure(SortDistData);

  struct SortDistData  {
    realtype dist;
    int      index,unmap1,unmap2;
  };

  DefineClass(SortDist);

  class SortDist : public CQuickSort  {
    public :
      SortDist() : CQuickSort() {}
      int  Compare ( int i, int j );
      void Swap    ( int i, int j );
      void Sort    ( PSortDistData sdata, int len );
    protected :
      PSortDistData sd;
  };


  DefineStructure(SuperposeData);

  struct SuperposeData  {
    PGraph        G;  // SSE graph
    PCMMDBManager M;  // the structure
    PSpAtom       a;  // atom superposition vector
    PPCAtom  Calpha;  // selected C-alphas
    PSSEDesc   SSED;  // SSE description vector
    pstr  selstring;  // C-alpha selection string
    int      selHnd;  // C-alpha selection handle
    int  selHndIncl;  // selection handle of inculded C-alphas
    int        nres;  // number of residues (C-alphas)
    int       nSSEs;  // number of SSEs
    void  Init   ();
    void  Dispose();
    void  DeselectCalphas();
    void  SelectCalphas  ();
  };


  DefineClass(Superpose);

  class Superpose  {

    public :
      Superpose();
      ~Superpose();

      void SetAllowMC         ( Boolean allowMisconnections );
      void SetIterationLimits ( int iter_max, int iter_min,
                                int max_hollow );
      void SetCaSelections    ( cpstr selection1, cpstr selection2 );

      int  SuperposeSSGraphs  ( PGraph G1, ivector F1,
                                PGraph G2, ivector F2,
                                int matchlen );

      //  driver #1
      int  SuperposeCalphas  (
              PGraph        G1,   // SSE graph of 1st structure
              PGraph        G2,   // SSE graph of 2nd structure
              ivector       F1,   // matched vertices of G1 [1..mlen]
              ivector       F2,   // matched vertices of G2 [1..mlen]
              int         mlen,   // length of match (F1,F2)
              PCMMDBManager M1,   // 1st structure
              PCMMDBManager M2,   // 2nd structure
              int  selHndIncl1=0, // sel handle to include atoms from M1
              int  selHndIncl2=0  // sel handle to include atoms from M2
                             );

      //  driver #2
      int  SuperposeCalphas  (
              PSuperposeData SD1, // superposition data of 1st structure
              PSuperposeData SD2, // superposition data of 2nd structure
              ivector        F1, // matched vertices of SD1.G [1..mlen]
              ivector        F2, // matched vertices of SD2.G [1..mlen]
              int            mlen  // length of match (F1,F2)
                             );

      void  GetTMatrix ( mat44 & TMat ); // to be applied to 1st struct.
      mat44 *  GetTMatrix    ();    // to be applied to 1st structure
      realtype GetRMSD       ();
      int      GetNAlign     ();
      void  GetSuperposition ( ivector  & Ca1  ,
                               rvector  & dist1, int & nCa1,
                               ivector  & Ca2  , int & nCa2,
                               mat44    & TMat ,
                               realtype & rmsdAchieved,
                               int & nAligned,   int & nGaps,
                               realtype & seqIdentity,
                               int & nMisD, realtype & nCombs );

      void GetCalphas1 ( PPCAtom & Calpha, int & numRes );
      void GetCalphas2 ( PPCAtom & Calpha, int & numRes );

      void GetSSEDesc1 ( RPSSEDesc sseDesc, int & numSSEs );
      void GetSSEDesc2 ( RPSSEDesc sseDesc, int & numSSEs );
      PSSEDesc GetSSEDesc1();
      PSSEDesc GetSSEDesc2();

      void GetSuperposedSSEs ( ivector v1, ivector v2, int & nSupSSEs );

      realtype GetCalphaQ   ()  { return Q_achieved; }
      realtype MatchQuality ( int Nalign, realtype Rmsd );

    protected :
      mat44     TMatrix,TMx;
      PSpAtom   a1,a2;
      realtype  Rmsd0;       // optimization parameter
      realtype  minContact;  // minimal Calpha-pair contact parameter
      realtype  maxContact;  // maximal Calpha-pair contact parameter
      realtype  maxRMSD;     // maximal RMSD allowed
      realtype  minQStep;    // minimal quality improvement that counts
      realtype  minCosine;   // min cosine between co-directional SSEs
      realtype  SSEweight;   // additional weight for SSE atoms
      int       sseGray;     // gray zone on the ends of SSEs allowed for
                             // matching to non-SSE atoms
      int       selInclHnd1; // selection handle for included Calpha1
      int       selInclHnd2; // selection handle for included Calpha2
      int       driverID;    // ID of the used Superpose driver
      pstr      selString1;  // optional sel-n string for 1st structure
      pstr      selString2;  // optional sel-n string for 2nd structure


      realtype  rmsd_achieved,Q_achieved,ncombs,seqIdent;
      int       shortSect1,shortSect2;
      int       iterMax,iterMin;
      int       maxHollowIt;  // maximal allowed number of consequtive
                              // iterations without quality improvement

      int       nres1,nres2,nalgn,ngaps,nmd,nmisdr;
      Boolean   allowMC;     // allowing for misconnection

      rmatrix   A,U,V, AD;
      rvector   W,RV1;

      ivector   copyF1,copyF2;   // copy pointers to input F1,F2
      int       copyFlen;        // length of FF1,FF2
      rvector   cax0,cay0,caz0;  // working arrays
      PSortDistData sdata;

      PCMMDBManager MMDB1,MMDB2; // copies of 1st and 2nd structure MMDBs

      PPCAtom   Calpha1,Calpha2;
      PSSEDesc  SSED1,SSED2;
      ivector   FH1,FS1,FH2,FS2;
      int       nSSEs1,nSSEs2;
      int       nFH1,nFS1,nFH2,nFS2;
      PPSectionDist SDist;
      int       SDistAlloc;

      SortDist sortDist;


      void  InitSuperpose       ();
      void  FreeMemory          ();
      void  SelectCalphas       ( PCMMDBManager MMDB, PGraph G,
                                  PPCAtom & Calpha, PSpAtom & a,
                                  int & nres, int & selHnd,
                                  int selInclHnd, cpstr selString );
      void  MapSSEs             ( PPCAtom Calpha, PSpAtom a, int nres,
                                  PGraph G, RPSSEDesc SSED,
                                  int & nSSEs );
      void  IdentifyUnmatchedSSEs ( ivector & FH, int & nFH,
                                  ivector & FS, int & nFS,
                                  ivector F, int mlen,
                                  PGraph G );
      void  GetSSESpseCenters   ( RSSEDesc Q1, RSSEDesc Q2,
                                  RSSEDesc T1, RSSEDesc T2,
                                  realtype & qc1, realtype & qc2,
                                  realtype & tc1, realtype & tc2 );
      int   FirstGuess          ( ivector F1, ivector F2, int mlen );
      void  ChooseFirstRotation ( int rotSSE1, int rotSSE2 );
      void  CalcDistance        ( int SSE1, int SSSE2, RSectionDist D );
      void  AlignSSEs           ( RSectionDist D, int unmap );
      Boolean isMC              ( int pos1, int pos2 );
      void  CorrespondSSEs      ( ivector F1, int nF1, ivector F2,
                                  int nF2, realtype rmsd_est );
      void  CorrespondContacts  ( PCMMDBManager M1, realtype rmsd_est );
      void  ExpandContact       ( RSContact c, int & ip, int & im,
                                  realtype maxDist2 );
      void  RecoverGaps         ( PPCAtom Ca1, PSpAtom at1, int nat1,
                                  PPCAtom Ca2, PSpAtom at2, int nat2,
                                  realtype thresh );
      void  CleanShortSections  ( PSpAtom at1, int nat1, PSpAtom at2 );

      int   CalculateTMatrix    ();
      void     CalcNGaps        ( PSpAtom a, int nres, int & Ng, int & Nm );
      realtype CalcNCombs       ( PGraph G, PSSEDesc SSED, int nSSEs,
                                  PSpAtom  a, int nres );
      realtype MatchQuality2    ( int Nalign, realtype dist2 );
      void     CalcQScore       ( RSSEDesc SSE1 );
      int   OptimizeNalign      ();
      void  UnmapExcluded       ( PSpAtom a1, PSpAtom a2, int nres1 );

      void  _superpose ( PGraph G1, PGraph G2, int & rc );

  };

}

#endif
