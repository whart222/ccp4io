// $Id$
// =================================================================
//
//    26.06.06   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  SSM_Align <interface>
//       ~~~~~~~~~
//  **** Project :  Structure alignment in 3D
//       ~~~~~~~~~
//  **** Classes :  ssm::Align ( Secondary Structure Matching )
//       ~~~~~~~~~  ssm::XAlign   ( Output alignment             )
//                  CXTAlign  ( Text output alignment        )
//
//  E. Krissinel, 2002-2006
//
// =================================================================
//

#include <string.h>
#include <ctype.h>

#include "mmdb/mmdb_tables.h"
#include "ssm_align.h"


//  ---------------------------  ssm::Align ------------------------

ssm::Align::Align() : CStream()  {
  InitAlign();
}

ssm::Align::Align ( RPCStream Object ) : CStream ( Object )  {
  InitAlign();
}

ssm::Align::~Align()  {
  FreeMemory();
}

void ssm::Align::FreeMemory()  {
  FreeVectorMemory ( Ca1     ,0 );
  FreeVectorMemory ( dist1   ,0 );
  FreeVectorMemory ( Ca2     ,0 );
  FreeVectorMemory ( pqvalues,0 );
  if (G1)  delete G1;
  if (G2)  delete G2;
  G1 = NULL;
  G2 = NULL;
  nMatches = 0;
}

void ssm::Align::InitAlign()  {

  Mat4Init ( TMatrix ); // transformation matrix

  cnCheck     = CONNECT_Flexible;
  rmsd        = 0.0;  // core rmsd achieved
  Qscore      = 0.0;  // core Q achieved
  nres1       = 0;    // number of residues in query  structure
  nres2       = 0;    // number of residues in target structure
  nsel1       = 0;    // number of residues in query  selection
  nsel2       = 0;    // number of residues in target selection
  nalgn       = 0;    // number of aligned residues
  ngaps       = 0;    // number of gaps
  nmd         = 0;    // number of misdirections
  seqIdentity = 0;    // sequence identity
  ncombs      = 1.0;  // number of SSE combinations

  selHndCa1   = 0;    // selection handle to used C-alphas in
                      //   query structure
  Ca1         = NULL; // C-alpha correspondence vector for
                      //   query structure
  dist1       = NULL; // optimized distances between the query
                      //   and target C-alphas
  selHndCa2   = 0;    // selection handle to used C-alphas in
                      //   target structure
  Ca2         = NULL; // C-alpha correspondence vector for
                      //   target structure

  G1          = NULL;
  G2          = NULL; // retained SSE graphs

  pqvalues    = NULL;
  nMatches    = 0;

}


void ssm::Align::MapSelections ( int & selHndCa, PCMMDBManager M,
                                 PGraph G, int selHnd,
                                 ivector & newID )  {
PPCAtom a;
int     nr,i,k;
  G->SelectCalphas ( M,selHndCa,NULL );
  if (selHnd>0)  {
    M->GetSelIndex   ( selHndCa,a,nr   );
    GetVectorMemory  ( newID,nr,0      );
    k = 0;
    for (i=0;i<nr;i++)
      if (a[i]->isInSelection(selHnd)) newID[i] = k++;
                                  else newID[i] = -1;
    M->Select ( selHndCa,STYPE_ATOM,selHnd,SKEY_AND );
  } else
    newID = NULL;
}


void ssm::Align::MakeSelections ( PCMMDBManager M1, int selHnd1,
                                 PCMMDBManager M2, int selHnd2 )  {
ivector newID1,newID2;
int     i,k;

  MapSelections ( selHndCa1,M1,G1,selHnd1,newID1 );
  MapSelections ( selHndCa2,M2,G2,selHnd2,newID2 );

  if (newID2)  {
    k = 0;
    if (newID1)  {
      for (i=0;i<nres2;i++)
        if (newID2[i]>=0)  {
          if (Ca2[i]>=0)  Ca2[k++] = newID1[Ca2[i]];
                    else  Ca2[k++] = -1;
        }
    } else  {
      for (i=0;i<nres2;i++)
        if (newID2[i]>=0)  Ca2[k++] = Ca2[i];
    }
  }

  if (newID1)  {
    k = 0;
    if (newID2)  {
      for (i=0;i<nres1;i++)
        if (newID1[i]>=0)  {
          if (Ca1[i]>=0)  Ca1[k] = newID2[Ca1[i]];
                    else  Ca1[k] = -1;
          dist1[k] = dist1[i];
          k++;
        }
    } else  {
      for (i=0;i<nres1;i++)
        if (newID1[i]>=0)  {
          Ca1[k]   = Ca1[i];
          dist1[k] = dist1[i];
          k++;
        }
    }
  }

  nsel1 = M1->GetSelLength ( selHndCa1 );
  nsel2 = M2->GetSelLength ( selHndCa2 );

  FreeVectorMemory ( newID1,0 );
  FreeVectorMemory ( newID2,0 );

}


int ssm::Align::align ( PCMMDBManager M1, PCMMDBManager M2,
                        PRECISION     precision,
                        CONNECTIVITY  connectivity,
                        int selHnd1,  int selHnd2 )  {
PPMatch  Match;
ivector  F1,F2;
realtype Q1;
int      i,nm;

  FreeMemory();

  SetMatchPrecision    ( precision    );
  SetConnectivityCheck ( connectivity );
  cnCheck = connectivity;

  U.SetUniqueMatch ( True );
  U.SetBestMatch   ( True );

  G1 = GetSSGraph ( M1,selHnd1,i );
  if (!G1)  return i;

  G2 = GetSSGraph ( M2,selHnd2,i );
  if (!G2)  return i+2;

  U.MatchGraphs ( G1,G2,1 );

  U.GetMatches ( Match,nMatches );
  if (nMatches<=0)  return RC_NoHits;

  GetVectorMemory ( pqvalues,nMatches,0 );
  for (i=0;i<nMatches;i++)
    pqvalues[i] = -1.0;

  Qscore = -0.5;
  for (i=0;i<nMatches;i++)
    if (Match[i])  {
      Match[i]->GetMatch ( F1,F2,nm );
      superpose.SuperposeCalphas ( G1,G2,F1,F2,nm,M1,M2,
                                   selHnd1,selHnd2 );
      Q1 = superpose.GetCalphaQ();
      if ((Q1>0.0) && (Q1>Qscore))  {
        Qscore = Q1;
        superpose.GetSuperposition ( Ca1,dist1,nres1,Ca2,nres2,
                                     TMatrix,rmsd,nalgn,ngaps,
                                     seqIdentity, nmd,ncombs );
      }
      pqvalues[i] = Q1;
    }

  if (Qscore>0.0)  {
    MakeSelections ( M1,selHnd1, M2,selHnd2 );
    return RC_Ok;
  }

  return RC_NoSuperposition;

}


int ssm::Align::AlignSelectedMatch ( PCMMDBManager M1,
                                     PCMMDBManager M2,
                                     PRECISION    precision,
                                     CONNECTIVITY connectivity,
                                     int          selHnd1,
                                     int          selHnd2,
                                     int          nselect )  {
PPMatch   Match;
ivector   F1,F2;
realtype  Q1;
int       i,nGMatches,nm;

  FreeMemory();

  SetMatchPrecision    ( precision    );
  SetConnectivityCheck ( connectivity );
  cnCheck = connectivity;

  U.SetUniqueMatch ( True );
  U.SetBestMatch   ( True );

  G1 = GetSSGraph ( M1,selHnd1,i );
  if (!G1)  return i;

  G2 = GetSSGraph ( M2,selHnd2,i );
  if (!G2)  return i+2;

  U.MatchGraphs ( G1,G2,1 );

  U.GetMatches ( Match,nGMatches );
  if (nGMatches<=0)  return RC_NoHits;

  if (nselect>=nGMatches)  {
//  printf(" There are only %d matches for this alignment", nGMatches);
    return RC_TooFewMatches;
  }

  Qscore = -0.5;

  if (Match[nselect])  {

    Match[nselect]->GetMatch   ( F1,F2,nm );
    superpose.SuperposeCalphas ( G1,G2,F1,F2,nm,M1,M2,
                                 selHnd1,selHnd2 );

    Q1 = superpose.GetCalphaQ();
    if (Q1>0.0)  {
      superpose.GetSuperposition ( Ca1,dist1,nres1,Ca2,nres2,TMatrix,
                                   rmsd,nalgn,ngaps,seqIdentity,
                                   nmd,ncombs );
      MakeSelections ( M1,selHnd1, M2,selHnd2 );
      return RC_Ok;
    }

  }

  return RC_NoSuperposition;

}


void ssm::Align::write ( RCFile f )  {
int i,j;

  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      f.WriteReal ( &(TMatrix[i][j]) );

  f.WriteInt  ( &cnCheck     );
  f.WriteReal ( &rmsd        );
  f.WriteInt  ( &nres1       );
  f.WriteInt  ( &nres2       );
  f.WriteInt  ( &nsel1       );
  f.WriteInt  ( &nsel2       );
  f.WriteInt  ( &nalgn       );
  f.WriteInt  ( &ngaps       );
  f.WriteInt  ( &nmd         );
  f.WriteReal ( &seqIdentity );
  f.WriteReal ( &ncombs      );

  if (Ca1 && (nsel1>0))
    for (i=0;i<nsel1;i++)  {
      f.WriteInt  ( &(Ca1  [i]) );
      f.WriteReal ( &(dist1[i]) );
    }
  if (Ca2 && (nsel2>0))
    for (i=0;i<nsel2;i++)
      f.WriteInt ( &(Ca2[i]) );

  StreamWrite ( f,G1 );
  StreamWrite ( f,G2 );

}

void ssm::Align::read ( RCFile f )  {
int i,j;

  FreeMemory();

  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      f.ReadReal ( &(TMatrix[i][j]) );

  f.ReadInt  ( &cnCheck     );
  f.ReadReal ( &rmsd        );
  f.ReadInt  ( &nres1       );
  f.ReadInt  ( &nres2       );
  f.ReadInt  ( &nsel1       );
  f.ReadInt  ( &nsel2       );
  f.ReadInt  ( &nalgn       );
  f.ReadInt  ( &ngaps       );
  f.ReadInt  ( &nmd         );
  f.ReadReal ( &seqIdentity );
  f.ReadReal ( &ncombs      );

  if (nsel1>0)  {
    GetVectorMemory ( Ca1  ,nsel1,0 );
    GetVectorMemory ( dist1,nsel1,0 );
    for (i=0;i<nsel1;i++)  {
      f.ReadInt  ( &(Ca1  [i]) );
      f.ReadReal ( &(dist1[i]) );
    }
  }
  if (nsel2>0)  {
    GetVectorMemory ( Ca2,nsel2,0 );
    for (i=0;i<nsel2;i++)
      f.ReadInt ( &(Ca2[i]) );
  }

  StreamRead ( f,G1 );
  StreamRead ( f,G2 );

}

namespace ssm  {
  MakeStreamFunctions(Align)
}


//  -----------------------------  ssm::XAlign --------------------------

ssm::XAlign::XAlign()  {
  XBlock1 = NULL;
  nBlock1 = 0;
  XBlock2 = NULL;
  nBlock2 = 0;
  algnLen = 0;
}

ssm::XAlign::~XAlign()  {
  FreeMemory();
}

void ssm::XAlign::FreeMemory()  {
  if (XBlock1)  delete[] XBlock1;
  if (XBlock2)  delete[] XBlock2;
  XBlock1 = NULL;
  nBlock1 = 0;
  XBlock2 = NULL;
  nBlock2 = 0;
  algnLen = 0;
}


void ssm::XAlign::customInit() {}

void ssm::XAlign::align (
                PGraph g1, PPCAtom Calpha1, ivector Ca1, int nat1,
                PGraph g2, PPCAtom Calpha2, ivector Ca2, int nat2,
                rvector dist1, int & nr )  {
int i,j;

  FreeMemory();

  a1     = Ca1;
  a2     = Ca2;
  alpha1 = Calpha1;
  alpha2 = Calpha2;
  sg1    = g1;
  sg2    = g2;
  d1     = dist1;
  na1    = nat1;
  na2    = nat2;

  nCols1 = makeXBlocks ( Ca1,nat1,XBlock1,nBlock1 );
  nCols2 = makeXBlocks ( Ca2,nat2,XBlock2,nBlock2 );
  nRows  = nat1 + nat2 + 2;

  maxdist = 0.0;
  for (i=0;i<nat1;i++)
    if (Ca1[i]>=0)  {
      if (dist1[i]>maxdist)  maxdist = dist1[i];
    }
  if (maxdist<=1.0e-2)  maxdist = 1.0;

  customInit();
  nr = 0;
  for (i=0;i<nBlock1;i++)
    for (j=0;j<nBlock2;j++)
      alignXBlocks ( XBlock1[i],XBlock2[j],nr );

  algnLen = nr;

}


int  ssm::XAlign::makeXBlocks ( ivector Ca, int nat, RPXBlock xBlock,
                                int & nBlocks )  {
//    Ca is considered as blocks of non-negative,
//  increasing-by-one numbers Ca[i]>=0, and negative
//  Ca[i]<0 surrounding them. Block boundaries are drawn
//  at the middle of negative-Ca[i] gaps.
//    nBlocks returns the number of such blocks, each block is
//  identified by the initial and final indices i1 and i2, and
//  by "index mass center" mc used for sorting.
//    Returns the number of fold-columns.
PXBlock  XB1;
int      nAlloc,i,j,i1,i2,ip1,ip2,iv,k,icol;
realtype mc;

  if (xBlock)  delete[] xBlock;
  xBlock  = NULL;
  nBlocks = 0;
  nAlloc  = 0;

  i  = 0;
  i1 = 0;  // begining of a block
  // begining of first block, check for leading negatives
  while (i<nat)
    if (Ca[i]<0)  i++;
            else  break;
  do  {
    if (i<nat)  {
      // check for increasing-by-one positives
      ip1 = i;
      iv  = Ca[i++];
      mc  = iv;
      k   = 1;
      while (i<nat)
        if (Ca[i]==iv+1)  {
          iv  = Ca[i++];
          mc += iv;
          k++;
        } else
          break;
      mc /= k;
      ip2 = i-1;  // increasing-by-one has stopped, check for negatives
      while (i<nat)
        if (Ca[i]<0)  i++;
                else  break;
      //  get i2 as end of a block
      if (i>=nat)       i2 = nat-1;     // the last block
      else if (i-ip2>1) i2 = (ip2+i)/2; // take the medium
                   else i2 = ip2;
    } else  {
      i2  = nat-1;  // the only empty (all-negatives) block
      ip1 = -1;
      ip2 = -1;
      mc  = 0.0;
    }
    // create new block
    if (nBlocks>=nAlloc)  {
      nAlloc += 20;
      XB1 = new XBlock[nAlloc];
      for (j=0;j<nBlocks;j++)  {
        XB1[j].i1   = xBlock[j].i1;
        XB1[j].i2   = xBlock[j].i2;
        XB1[j].ip1  = xBlock[j].ip1;
        XB1[j].ip2  = xBlock[j].ip2;
        XB1[j].mc   = xBlock[j].mc;
        XB1[j].icol = xBlock[j].icol;
      }
      delete[] xBlock;
      xBlock = XB1;
    }
    xBlock[nBlocks].i1   = i1;
    xBlock[nBlocks].i2   = i2;
    xBlock[nBlocks].ip1  = ip1;
    xBlock[nBlocks].ip2  = ip2;
    xBlock[nBlocks].mc   = mc;
    xBlock[nBlocks].icol = 0;
    nBlocks++;
    i1 = i2+1;
  } while (i<nat);

  // assign fold-columns to the blocks
  icol = 0;
  do  {
    icol++;
    i  = 0;
    iv = 0;
    while (i<nBlocks)  {
      mc = MaxReal;
      k  = -1;
      for (j=i;j<nBlocks;j++)
        if ((xBlock[j].icol==0) && (xBlock[j].mc<mc))  {
          mc = xBlock[j].mc;
          k  = j;
        }
      if (k>=0)  {
        xBlock[k].icol = icol;
        i  = k+1;
        iv = 1;
      } else
        i = nBlocks;
    }
  } while (iv);

  return icol-1;

}


void ssm::XAlign::alignXBlocks ( RXBlock B1, RXBlock B2, int & nr )  {
int  l1,l2, i1,i2, sseType1,sseType2, icol;

  if (((a1[B1.ip1]>=B2.ip1) && (a1[B1.ip1]<=B2.ip2)) ||
      ((a1[B1.ip2]>=B2.ip1) && (a1[B1.ip2]<=B2.ip2)) ||
      ((a2[B2.ip1]>=B1.ip1) && (a2[B2.ip1]<=B1.ip2)) ||
      ((a2[B2.ip2]>=B1.ip1) && (a2[B2.ip2]<=B1.ip2)))  {

    if (a1[B1.ip1]<B2.ip1)  {
      l1 = 0;
      l2 = B2.ip1 - B2.i1;
      i1 = a2[B2.ip1];
      i2 = B2.i1;
    } else if (a1[B1.ip1]==B2.ip1)  {
      l1 = B1.ip1 - B1.i1;
      l2 = B2.ip1 - B2.i1;
      i1 = B1.i1;
      i2 = B2.i1;
    } else  {
      l1 = B1.ip1 - B1.i1;  // number of leading unmappings
      l2 = 0;
      i1 = B1.i1;
      i2 = a1[B1.ip1];
    }

    icol = B1.icol;

    while (l1>l2)  {
      if (alpha1[i1])
           sseType1 = sg1->GetSSEType ( alpha1[i1]->GetChainID(),i1 );
      else sseType1 = V_UNKNOWN;
      makeRow ( alpha1[i1],sseType1,NULL,V_UNKNOWN,
                d1[i1],nr++,icol,False );
      i1++;
      l1--;
    }
    while (l2>l1)  {
      if (alpha2[i2])
           sseType2 = sg2->GetSSEType ( alpha2[i2]->GetChainID(),i2 );
      else sseType2 = V_UNKNOWN;
      makeRow ( NULL,V_UNKNOWN,alpha2[i2++],sseType2,
                -1.0,nr++,icol,False );
      l2--;
    }
    while (l2>0)  {
      if (alpha1[i1])
           sseType1 = sg1->GetSSEType ( alpha1[i1]->GetChainID(),i1 );
      else sseType1 = V_UNKNOWN;
      if (alpha2[i2])
           sseType2 = sg2->GetSSEType ( alpha2[i2]->GetChainID(),i2 );
      else sseType2 = V_UNKNOWN;
      makeRow ( alpha1[i1],sseType1,alpha2[i2++],sseType2,
                d1[i1],nr++,icol,False );
      i1++;
      l2--;
    }

    l1 = IMin ( B1.ip2-i1, B2.ip2-i2 ) + 1;
    while (l1>0)  {
      if (alpha1[i1])
           sseType1 = sg1->GetSSEType ( alpha1[i1]->GetChainID(),i1 );
      else sseType1 = V_UNKNOWN;
      if (alpha2[i2])
           sseType2 = sg2->GetSSEType ( alpha2[i2]->GetChainID(),i2 );
      else sseType2 = V_UNKNOWN;
      makeRow ( alpha1[i1],sseType1,alpha2[i2++],sseType2,
                d1[i1],nr++,icol,True );
      i1++;
      l1--;
    }

    if (i1<=B1.ip2)  {
      l1 = 0;
      l2 = B2.i2 - i2 + 1;
    } else if (i2<=B2.ip2)  {
      l1 = B1.i2 - i1 + 1;
      l2 = 0;
    } else  {
      l1 = B1.i2 - i1 + 1;
      l2 = B2.i2 - i2 + 1;
    }
    while ((l1>0) && (l2>0))  {
      if (alpha1[i1])
           sseType1 = sg1->GetSSEType ( alpha1[i1]->GetChainID(),i1 );
      else sseType1 = V_UNKNOWN;
      if (alpha2[i2])
           sseType2 = sg2->GetSSEType ( alpha2[i2]->GetChainID(),i2 );
      else sseType2 = V_UNKNOWN;
      makeRow ( alpha1[i1],sseType1,alpha2[i2++],sseType2,
                d1[i1],nr++,icol,False );
      i1++;
      l1--;
      l2--;
    }
    while (l1>0)  {
      if (alpha1[i1])
           sseType1 = sg1->GetSSEType ( alpha1[i1]->GetChainID(),i1 );
      else sseType1 = V_UNKNOWN;
      makeRow ( alpha1[i1],sseType1,NULL,V_UNKNOWN,
                d1[i1],nr++,icol,False );
      i1++;
      l1--;
    }
    while (l2>0)  {
      if (alpha2[i2])
           sseType2 = sg2->GetSSEType ( alpha2[i2]->GetChainID(),i2 );
      else sseType2 = V_UNKNOWN;
      makeRow ( NULL,V_UNKNOWN,alpha2[i2++],sseType2,
                -1.0,nr++,icol,False );
      l2--;
    }

  }

}


void ssm::XAlign::makeRow ( PCAtom A1, int sseType1,
                        PCAtom A2, int sseType2,
                        realtype dist, int rowNo,
                        int icol, Boolean aligned )  {
UNUSED_ARGUMENT(A1);
UNUSED_ARGUMENT(sseType1);
UNUSED_ARGUMENT(A2);
UNUSED_ARGUMENT(sseType2);
UNUSED_ARGUMENT(dist);
UNUSED_ARGUMENT(rowNo);
UNUSED_ARGUMENT(icol);
UNUSED_ARGUMENT(aligned);
}



//  ----------------------------  CXTAlign --------------------------

namespace ssm  {

  void  PrintAtom ( RCFile f, int sseType, realtype hydropathy,
                    ChainID chID, ResName resName, int seqNum,
                    InsCode insCode )  {
  char sse[2],hp[2],ch[3],S[200];

    if (sseType==V_HELIX)       sse[0] = 'H';
    else if (sseType==V_STRAND) sse[0] = 'S';
                           else sse[0] = ' ';
    sse[1] = char(0);

    if ((-5.0<hydropathy) && (hydropathy<5.0))  {
      if (hydropathy>=-0.5)      hp[0] = '-';
      else if (hydropathy<=-1.5) hp[0] = '+';
                            else hp[0] = '.';
    } else
      hp[0] = ' ';
    hp[1] = char(0);

    if ((!chID[0]) || (chID[0]==' '))  {
      ch[0] = ' ';      ch[1] = ' ';
    } else  {
      ch[0] = chID[0];  ch[1] = ':';
    }
    ch[2] = char(0);

    sprintf ( S," |%1s%1s %2s%3s%4i%1s|",
                sse,hp,ch,resName,seqNum,insCode );
    f.Write ( S );

  }

}

void  ssm::XTAlign::Print ( RCFile f )  {
char S[100],SI[10];
int  i;

  if (alignKey<4)  {

    if (alignKey!=2)
      PrintAtom ( f,sseType1,hydropathy1,chID1,
                    resName1,seqNum1,insCode1 );
    else
      f.Write ( " |             |" );

    if (alignKey==0)  {
      switch (simindex)  {
        case 5 :  strcpy ( SI,"*****" );  break;
        case 4 :  strcpy ( SI,"+++++" );  break;
        case 3 :  strcpy ( SI,"=====" );  break;
        case 2 :  strcpy ( SI,"-----" );  break;
        case 1 :  strcpy ( SI,":::::" );  break;
        default:
        case 0 :  strcpy ( SI,"....." );  break;
      }
      SI[1] = char(0);
      sprintf ( S," %1s%5.2f%1s%1s",SI,dist,SI,SI );
      if (S[3]==' ')  S[3] = SI[0];
      f.Write ( S );
      SI[1] = SI[0];
      for (i=1;i<loopNo;i++)  f.Write ( SI );
//      f.Write ( " " );
    } else  {
      f.Write ( "        " );
      for (i=1;i<loopNo;i++)  f.Write ( "     " );
      f.Write ( " " );
    }

    if (alignKey!=3)
      PrintAtom ( f,sseType2,hydropathy2,chID2,
                    resName2,seqNum2,insCode2 );
    else
      f.Write ( " |             |" );

    f.LF();

  }

}


ssm::XAlignText::XAlignText() : ssm::XAlign() {
  R = NULL;
}

ssm::XAlignText::~XAlignText() {
  customFree();
}

void ssm::XAlignText::customFree()  {
  if (R)  delete[] R;
  R = NULL;
}

void ssm::XAlignText::customInit()  {
int i;
  customFree();
  R = new ssm::XTAlign[nRows];
  for (i=0;i<nRows;i++)
    R[i].alignKey = 5;
}

void  ssm::XAlignText::WipeTextRows()  {
  R = NULL;
}

void ssm::XAlignText::makeRow ( PCAtom A1, int sseType1,
                            PCAtom A2, int sseType2,
                            realtype dist, int rowNo, int icol,
                            Boolean aligned )  {

  if (aligned)  R[rowNo].alignKey = 0;
          else  R[rowNo].alignKey = 1;

  if (A1)  {
    R[rowNo].sseType1    = sseType1;
    R[rowNo].hydropathy1 = A1->GetAAHydropathy();
    R[rowNo].seqNum1     = A1->GetSeqNum      ();
    strcpy ( R[rowNo].chID1   ,A1->GetChainID() );
    strcpy ( R[rowNo].resName1,A1->GetResName() );
    strcpy ( R[rowNo].insCode1,A1->GetInsCode() );
  } else
    R[rowNo].alignKey = 2;

  if (A2)  {
    R[rowNo].sseType2    = sseType2;
    R[rowNo].hydropathy2 = A2->GetAAHydropathy();
    R[rowNo].seqNum2     = A2->GetSeqNum      ();
    strcpy ( R[rowNo].chID2   ,A2->GetChainID() );
    strcpy ( R[rowNo].resName2,A2->GetResName() );
    strcpy ( R[rowNo].insCode2,A2->GetInsCode() );
  } else
    R[rowNo].alignKey = 3;

  if ((!A1) && (!A2))  R[rowNo].alignKey = 4;

  R[rowNo].simindex = -5;
  R[rowNo].dist     = -1.0;
  if (aligned)  {
    if (A1 && A2)  R[rowNo].simindex = A1->GetAASimilarity ( A2 );
             else  R[rowNo].simindex = -5;
    R[rowNo].dist = dist;
  }

  R[rowNo].loopNo = icol;

}

void  ssm::XAlignText::GetAlignments ( pstr & algn1, pstr & algn2 )  {
char rn1[10];
char rn2[10];
int i;
  if (algn1)  delete[] algn1;
  if (algn2)  delete[] algn2;
  if (algnLen>0)  {
    algn1 = new char[algnLen+1];
    algn2 = new char[algnLen+1];
    for (i=0;i<algnLen;i++)  {
      if (R[i].alignKey<=3)  {
        if (R[i].alignKey!=2)
              Get1LetterCode ( R[i].resName1,rn1 );
        else  strcpy ( rn1,"-" );
        if (R[i].alignKey!=3)
              Get1LetterCode ( R[i].resName2,rn2 );
        else  strcpy ( rn2,"-" );
        if (R[i].alignKey==0)  {
          rn1[0] = char(toupper(int(rn1[0])));
          rn2[0] = char(toupper(int(rn2[0])));
        } else  {
          rn1[0] = char(tolower(int(rn1[0])));
          rn2[0] = char(tolower(int(rn2[0])));
        }
      } else  {
        strcpy ( rn1,"-" );
        strcpy ( rn2,"-" );
      }
      algn1[i] = rn1[0];
      algn2[i] = rn2[0];
    }
    algn1[algnLen] = char(0);
    algn2[algnLen] = char(0);
  } else  {
    algn1 = NULL;
    algn2 = NULL;
  }
}


namespace ssm  {

  void PrintAlignTable ( RCFile f, PCMMDBManager M1, PCMMDBManager M2,
                                   PAlign SSMAlign )  {
  XAlignText CXA;
  PXTAlign   XTA;
  PPCAtom     Calpha1,Calpha2;
  int         nat1,nat2,nr,j;

    M1->GetSelIndex ( SSMAlign->selHndCa1,Calpha1,nat1 );
    M2->GetSelIndex ( SSMAlign->selHndCa2,Calpha2,nat2 );

    CXA.align ( SSMAlign->G1,Calpha1,SSMAlign->Ca1,nat1,
                SSMAlign->G2,Calpha2,SSMAlign->Ca2,nat2,
                SSMAlign->dist1,nr );
    f.LF();

    if (SSMAlign->cnCheck!=CONNECT_None)  {
      f.WriteLine ( " .-------------.----------.-------------." );
      f.WriteLine ( " |    Query    | Dist.(A) |   Target    |" );
      f.WriteLine ( " |-------------+----------+-------------|" );
    } else  {
      f.WriteLine (
      " .-------------.----------.-----------------------------------");
      f.WriteLine (
      " |    Query    | Dist.(A) |   Target"                          );
      f.WriteLine (
      " |-------------+----------+-----------------------------------");
    }

    XTA = CXA.GetTextRows();
    for (j=0;j<nr;j++)
      XTA[j].Print ( f );

    if (SSMAlign->cnCheck!=CONNECT_None)
      f.WriteLine ( " `-------------'----------'-------------'" );
    else
      f.WriteLine (
      " `-------------'----------'-----------------------------------");
    f.LF();

    f.WriteLine ( " Notations:" );
    f.WriteLine ( " S/H   residue belongs to a strand/helix" );
    f.WriteLine ( " +/-/. hydrophylic/hydrophobic/neutral residue" );
    f.WriteLine ( " **    identical residues matched: similarity 5" );
    f.WriteLine ( " ++    similarity 4" );
    f.WriteLine ( " ==    similarity 3" );
    f.WriteLine ( " --    similarity 2" );
    f.WriteLine ( " ::    similarity 1" );
    f.WriteLine ( " ..    dissimilar residues: similarity 0" );

  }

}
