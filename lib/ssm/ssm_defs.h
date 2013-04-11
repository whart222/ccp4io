// $Id: ssm_defs.h,v 1.2 2005/12/20 12:09:14 keb Exp $
// =================================================================
//
//    05.04.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// -----------------------------------------------------------------
//
//  **** Module  :  ssm_defs  <interface>
//       ~~~~~~~~~
//
//  E. Krissinel 2002-2013
//
// =================================================================
//


#ifndef  __SSM_Defs__
#define  __SSM_Defs__


namespace ssm  {

  #define SSGP_Distance   0
  #define SSGP_Alpha1     1
  #define SSGP_Alpha2     2
  #define SSGP_Alpha3     3
  #define SSGP_Alpha4     4
  #define SSGP_dAlpha1    5
  #define SSGP_dAlpha2    6
  #define SSGP_dAlpha3    7
  #define SSGP_dAlpha4    8

  #define SSGE_Ok                     0
  #define SSGE_NoVertices             70
  #define SSGE_UnmatchedConnectivity  5001
  #define SSGE_AlignError             5002
  #define SSGE_WrongSelLine1          5003
  #define SSGE_WrongSelLine2          5004
  #define SSGE_WrongSelLine3          5005

  #define SSGT_None       0
  #define SSGT_PDB        1
  #define SSGT_SCOP       2
  #define SSGT_PDBDOMAIN  3
  #define SSGT_PDBRANGE   4
  #define SSGT_CFDOMAIN   5
  #define SSGT_CFRANGE    6

  #define SSMF_UniqueMatch       0x00000001
  #define SSMF_BestMatch         0x00000002
  #define SSMF_WrongConnectOnly  0x00000004

  #define MALIGN_Ok             0
  #define MALIGN_BadInput       1
  #define MALIGN_NoStructure    2
  #define MALIGN_NoAlignment    3
  #define MALIGN_NoGraph     1000

  #define UNMAP_YES          (-2)
  #define UNMAP_NO           (-1)

  enum SUPERPOSITION_RESULT {
    SPOSE_Ok,SPOSE_BadData,SPOSE_NoCalphas1,SPOSE_NoCalphas2,
    SPOSE_RemoteStruct,SPOSE_SVDFail
  };

  enum RETURN_CODE {
    RC_Ok,RC_NoHits,RC_NoSuperposition,RC_NoGraph,RC_NoVertices,
    RC_NoGraph2,RC_NoVertices2,RC_TooFewMatches
  };

  //  precision level conatsnts
  enum PRECISION { PREC_Highest,PREC_High,PREC_Normal,
                   PREC_Low,PREC_Lowest };

  //  regimes of checking the SS connectivity
  enum CONNECTIVITY { CONNECT_None,CONNECT_Flexible,CONNECT_Strict };

  enum VERTEX_TYPE { V_UNKNOWN=-1,V_HELIX,V_STRAND };

}


#endif
