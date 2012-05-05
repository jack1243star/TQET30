/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.  
 *
 * Copyright (c) 2010-2012, ITU/ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the ITU/ISO/IEC nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/** \file     TComSlice.h
    \brief    slice header and SPS class (header)
*/

#ifndef __TCOMSLICE__
#define __TCOMSLICE__

#include <cstring>
#include <map>
#include "CommonDef.h"
#include "TComRom.h"
#include "TComList.h"

//! \ingroup TLibCommon
//! \{

class TComPic;
class TComTrQuant;

// ====================================================================================================================
// Constants
// ====================================================================================================================

/// max number of supported APS in software
#define MAX_NUM_SUPPORTED_APS 1

// ====================================================================================================================
// Class definition
// ====================================================================================================================

#if RPS_IN_SPS
/// Reference Picture Set class
class TComReferencePictureSet
{
private:
  Int  m_numberOfPictures;
  Int  m_numberOfNegativePictures;
  Int  m_numberOfPositivePictures;
  Int  m_numberOfLongtermPictures;
  Int  m_deltaPOC[MAX_NUM_REF_PICS];
  Int  m_POC[MAX_NUM_REF_PICS];
  Bool m_used[MAX_NUM_REF_PICS];
  Bool m_interRPSPrediction;
  Int  m_deltaRIdxMinus1;   
  Int  m_deltaRPS; 
  Int  m_numRefIdc; 
  Int  m_refIdc[MAX_NUM_REF_PICS+1];
#if LTRP_MULT
  Bool m_bCheckLTMSB[MAX_NUM_REF_PICS];
#endif

public:
  TComReferencePictureSet();
  virtual ~TComReferencePictureSet();

  Void setUsed(Int bufferNum, Bool used);
  Void setDeltaPOC(Int bufferNum, Int deltaPOC);
  Void setPOC(Int bufferNum, Int deltaPOC);
  Void setNumberOfPictures(Int numberOfPictures);
#if LTRP_MULT
  Void      setCheckLTMSBPresent     (Int bufferNum, Bool b );
  Bool      getCheckLTMSBPresent     (Int bufferNum);
#endif

  Int  getUsed(Int bufferNum);
  Int  getDeltaPOC(Int bufferNum);
  Int  getPOC(Int bufferNum);
  Int  getNumberOfPictures();

  Void setNumberOfNegativePictures(Int number)  { m_numberOfNegativePictures = number; }
  Int  getNumberOfNegativePictures()            { return m_numberOfNegativePictures; }
  Void setNumberOfPositivePictures(Int number)  { m_numberOfPositivePictures = number; }
  Int  getNumberOfPositivePictures()            { return m_numberOfPositivePictures; }
  Void setNumberOfLongtermPictures(Int number)  { m_numberOfLongtermPictures = number; }
  Int  getNumberOfLongtermPictures()            { return m_numberOfLongtermPictures; }

  Void setInterRPSPrediction(Bool flag)         { m_interRPSPrediction = flag; }
  Bool getInterRPSPrediction()                  { return m_interRPSPrediction; }
  Void setDeltaRIdxMinus1(Int x)                { m_deltaRIdxMinus1 = x; }
  Int  getDeltaRIdxMinus1()                     { return m_deltaRIdxMinus1; }
  Void setDeltaRPS(Int x)                       { m_deltaRPS = x; }
  Int  getDeltaRPS()                            { return m_deltaRPS; }
  Void setNumRefIdc(Int x)                      { m_numRefIdc = x; }
  Int  getNumRefIdc()                           { return m_numRefIdc; }

  Void setRefIdc(Int bufferNum, Int refIdc);
  Int  getRefIdc(Int bufferNum);

  Void sortDeltaPOC();
  Void printDeltaPOC();
};

/// Reference Picture Set set class
class TComRPSList
{
private:
  Int  m_numberOfReferencePictureSets;
  TComReferencePictureSet* m_referencePictureSets;
  
public:
  TComRPSList();
  virtual ~TComRPSList();
  
  Void  create  (Int numberOfEntries);
  Void  destroy ();


  TComReferencePictureSet* getReferencePictureSet(Int referencePictureSetNum);
  Int getNumberOfReferencePictureSets();
  Void setNumberOfReferencePictureSets(Int numberOfReferencePictureSets);
};
#endif
/// SPS class
class TComSPS
{
private:
  Int         m_SPSId;
  Int         m_ProfileIdc;
  Int         m_LevelIdc;
  Int         m_chromaFormatIdc;

  UInt        m_uiMaxTLayers;           // maximum number of temporal layers

  // Structure
  UInt        m_picWidthInLumaSamples;
  UInt        m_picHeightInLumaSamples;
  Bool        m_picCroppingFlag;
  Int         m_picCropLeftOffset;
  Int         m_picCropRightOffset;
  Int         m_picCropTopOffset;
  Int         m_picCropBottomOffset;
  UInt        m_uiMaxCUWidth;
  UInt        m_uiMaxCUHeight;
  UInt        m_uiMaxCUDepth;
  UInt        m_uiMinTrDepth;
  UInt        m_uiMaxTrDepth;
#if RPS_IN_SPS
  TComRPSList* m_RPSList;
  Bool        m_bLongTermRefsPresent;
#endif
#if H0567_DPB_PARAMETERS_PER_TEMPORAL_LAYER
  Int         m_numReorderPics[MAX_TLAYER];
#else
  Int         m_maxNumberOfReferencePictures;
  Int         m_numReorderFrames;
#endif
  
  // Tool list
  UInt        m_uiQuadtreeTULog2MaxSize;
  UInt        m_uiQuadtreeTULog2MinSize;
  UInt        m_uiQuadtreeTUMaxDepthInter;
  UInt        m_uiQuadtreeTUMaxDepthIntra;
  Bool        m_usePCM;
  UInt        m_pcmLog2MaxSize;
  UInt        m_uiPCMLog2MinSize;
  Bool        m_bDisInter4x4;
  Bool        m_useAMP;
  Bool        m_bUseALF;
#if LCU_SYNTAX_ALF
  Bool        m_bALFCoefInSlice;
#endif
  Bool        m_bUseLMChroma; // JL:

  Bool        m_bUseLComb;
  Bool        m_bLCMod;
  Bool        m_useNSQT;
  
#if H0412_REF_PIC_LIST_RESTRICTION
  Bool        m_restrictedRefPicListsFlag;
  Bool        m_listsModificationPresentFlag;
#endif

  // Parameter
  AMVP_MODE   m_aeAMVPMode[MAX_CU_DEPTH];
  UInt        m_uiBitDepth;
  UInt        m_uiBitIncrement;
#if H0736_AVC_STYLE_QP_RANGE
  Int         m_qpBDOffsetY;
  Int         m_qpBDOffsetC;
#endif

#if LOSSLESS_CODING
  Bool        m_useLossless;
#endif

  UInt        m_uiPCMBitDepthLuma;
  UInt        m_uiPCMBitDepthChroma;
  Bool        m_bPCMFilterDisableFlag;

  UInt        m_uiBitsForPOC;
  // Max physical transform size
  UInt        m_uiMaxTrSize;
  
  Int m_iAMPAcc[MAX_CU_DEPTH];

  Bool        m_bLFCrossSliceBoundaryFlag;
  Bool        m_bUseSAO; 

  Bool     m_bLFCrossTileBoundaryFlag;
  Int      m_iUniformSpacingIdr;
  Int      m_iTileBoundaryIndependenceIdr;
  Int      m_iNumColumnsMinus1;
  UInt*    m_puiColumnWidth;
  Int      m_iNumRowsMinus1;
  UInt*    m_puiRowHeight;
  
  Bool        m_bTemporalIdNestingFlag; // temporal_id_nesting_flag

  Bool        m_scalingListEnabledFlag;
#if H0567_DPB_PARAMETERS_PER_TEMPORAL_LAYER
  UInt        m_uiMaxDecPicBuffering[MAX_TLAYER]; 
  UInt        m_uiMaxLatencyIncrease[MAX_TLAYER];
#else
  UInt        m_uiMaxDecFrameBuffering; 
  UInt        m_uiMaxLatencyIncrease;
#endif

  Bool        m_useDF;

#if TILES_WPP_ENTRY_POINT_SIGNALLING
  UInt        m_tilesOrEntropyCodingSyncIdc;
  Int         m_numSubstreams;
#endif

public:
  TComSPS();
  virtual ~TComSPS();

  Int  getSPSId       ()         { return m_SPSId;          }
  Void setSPSId       (Int i)    { m_SPSId = i;             }
  Int  getProfileIdc  ()         { return m_ProfileIdc;     }
  Void setProfileIdc  (Int i)    { m_ProfileIdc = i;        }
  Int  getLevelIdc    ()         { return m_LevelIdc;       }
  Void setLevelIdc    (Int i)    { m_LevelIdc = i;          }

  Int  getChromaFormatIdc ()         { return m_chromaFormatIdc;       }
  Void setChromaFormatIdc (Int i)    { m_chromaFormatIdc = i;          }
  
  // structure
  Void setPicWidthInLumaSamples       ( UInt u ) { m_picWidthInLumaSamples = u;        }
  UInt getPicWidthInLumaSamples       ()         { return  m_picWidthInLumaSamples;    }
  Void setPicHeightInLumaSamples      ( UInt u ) { m_picHeightInLumaSamples = u;       }
  UInt getPicHeightInLumaSamples      ()         { return  m_picHeightInLumaSamples;   }

  Bool getPicCroppingFlag() const          { return m_picCroppingFlag; }
  Void setPicCroppingFlag(Bool val)        { m_picCroppingFlag = val; }
  Int  getPicCropLeftOffset() const        { return m_picCropLeftOffset; }
  Void setPicCropLeftOffset(Int val)       { m_picCropLeftOffset = val; }
  Int  getPicCropRightOffset() const       { return m_picCropRightOffset; }
  Void setPicCropRightOffset(Int val)      { m_picCropRightOffset = val; }
  Int  getPicCropTopOffset() const         { return m_picCropTopOffset; }
  Void setPicCropTopOffset(Int val)        { m_picCropTopOffset = val; }
  Int  getPicCropBottomOffset() const      { return m_picCropBottomOffset; }
  Void setPicCropBottomOffset(Int val)     { m_picCropBottomOffset = val; }

  Void setMaxCUWidth  ( UInt u ) { m_uiMaxCUWidth = u;      }
  UInt getMaxCUWidth  ()         { return  m_uiMaxCUWidth;  }
  Void setMaxCUHeight ( UInt u ) { m_uiMaxCUHeight = u;     }
  UInt getMaxCUHeight ()         { return  m_uiMaxCUHeight; }
  Void setMaxCUDepth  ( UInt u ) { m_uiMaxCUDepth = u;      }
  UInt getMaxCUDepth  ()         { return  m_uiMaxCUDepth;  }
  Void setUsePCM      ( Bool b ) { m_usePCM = b;           }
  Bool getUsePCM      ()         { return m_usePCM;        }
  Void setPCMLog2MaxSize  ( UInt u ) { m_pcmLog2MaxSize = u;      }
  UInt getPCMLog2MaxSize  ()         { return  m_pcmLog2MaxSize;  }
  Void setPCMLog2MinSize  ( UInt u ) { m_uiPCMLog2MinSize = u;      }
  UInt getPCMLog2MinSize  ()         { return  m_uiPCMLog2MinSize;  }
  Void setBitsForPOC  ( UInt u ) { m_uiBitsForPOC = u;      }
  UInt getBitsForPOC  ()         { return m_uiBitsForPOC;   }
  Bool getDisInter4x4()         { return m_bDisInter4x4;        }
  Void setDisInter4x4      ( Bool b ) { m_bDisInter4x4  = b;          }
  Bool getUseAMP() { return m_useAMP; }
  Void setUseAMP( Bool b ) { m_useAMP = b; }
  Void setMinTrDepth  ( UInt u ) { m_uiMinTrDepth = u;      }
  UInt getMinTrDepth  ()         { return  m_uiMinTrDepth;  }
  Void setMaxTrDepth  ( UInt u ) { m_uiMaxTrDepth = u;      }
  UInt getMaxTrDepth  ()         { return  m_uiMaxTrDepth;  }
  Void setQuadtreeTULog2MaxSize( UInt u ) { m_uiQuadtreeTULog2MaxSize = u;    }
  UInt getQuadtreeTULog2MaxSize()         { return m_uiQuadtreeTULog2MaxSize; }
  Void setQuadtreeTULog2MinSize( UInt u ) { m_uiQuadtreeTULog2MinSize = u;    }
  UInt getQuadtreeTULog2MinSize()         { return m_uiQuadtreeTULog2MinSize; }
  Void setQuadtreeTUMaxDepthInter( UInt u ) { m_uiQuadtreeTUMaxDepthInter = u;    }
  Void setQuadtreeTUMaxDepthIntra( UInt u ) { m_uiQuadtreeTUMaxDepthIntra = u;    }
  UInt getQuadtreeTUMaxDepthInter()         { return m_uiQuadtreeTUMaxDepthInter; }
  UInt getQuadtreeTUMaxDepthIntra()         { return m_uiQuadtreeTUMaxDepthIntra; }
#if H0567_DPB_PARAMETERS_PER_TEMPORAL_LAYER
  Void setNumReorderPics(Int i, UInt tlayer)              { m_numReorderPics[tlayer] = i;    }
  Int  getNumReorderPics(UInt tlayer)                     { return m_numReorderPics[tlayer]; }
#else
  Void setMaxNumberOfReferencePictures( Int u )  { m_maxNumberOfReferencePictures = u;    }
  Int  getMaxNumberOfReferencePictures()         { return m_maxNumberOfReferencePictures; }
  Void setNumReorderFrames( Int i )              { m_numReorderFrames = i;    }
  Int  getNumReorderFrames()                     { return m_numReorderFrames; }
#endif
#if RPS_IN_SPS
  Void      setRPSList( TComRPSList* RPSList )   { m_RPSList = RPSList;       }
  TComRPSList* getRPSList()                      { return m_RPSList;          }
  Bool      getLongTermRefsPresent()         { return m_bLongTermRefsPresent; }
  Void      setLongTermRefsPresent(Bool b)   { m_bLongTermRefsPresent=b;      }
#endif
  
  // physical transform
  Void setMaxTrSize   ( UInt u ) { m_uiMaxTrSize = u;       }
  UInt getMaxTrSize   ()         { return  m_uiMaxTrSize;   }
  
  // Tool list
  Bool getUseALF      ()         { return m_bUseALF;        }
#if LCU_SYNTAX_ALF
  Void setUseALFCoefInSlice(Bool b) {m_bALFCoefInSlice = b;}
  Bool getUseALFCoefInSlice()    {return m_bALFCoefInSlice;}
#endif

  Void setUseALF      ( Bool b ) { m_bUseALF  = b;          }
  Void setUseLComb    (Bool b)   { m_bUseLComb = b;         }
  Bool getUseLComb    ()         { return m_bUseLComb;      }
  Void setLCMod       (Bool b)   { m_bLCMod = b;     }
  Bool getLCMod       ()         { return m_bLCMod;  }

  Bool getUseLMChroma ()         { return m_bUseLMChroma;        }
  Void setUseLMChroma ( Bool b ) { m_bUseLMChroma  = b;          }

#if LOSSLESS_CODING
  Bool getUseLossless ()         { return m_useLossless; }
  Void setUseLossless ( Bool b ) { m_useLossless  = b; }
#endif
  Bool getUseNSQT() { return m_useNSQT; }
  Void setUseNSQT( Bool b ) { m_useNSQT = b; }
  
#if H0412_REF_PIC_LIST_RESTRICTION
  Bool getRestrictedRefPicListsFlag    ()          { return m_restrictedRefPicListsFlag;   }
  Void setRestrictedRefPicListsFlag    ( Bool b )  { m_restrictedRefPicListsFlag = b;      }
  Bool getListsModificationPresentFlag ()          { return m_listsModificationPresentFlag; }
  Void setListsModificationPresentFlag ( Bool b )  { m_listsModificationPresentFlag = b;    }
#endif

  // AMVP mode (for each depth)
  AMVP_MODE getAMVPMode ( UInt uiDepth ) { assert(uiDepth < g_uiMaxCUDepth);  return m_aeAMVPMode[uiDepth]; }
  Void      setAMVPMode ( UInt uiDepth, AMVP_MODE eMode) { assert(uiDepth < g_uiMaxCUDepth);  m_aeAMVPMode[uiDepth] = eMode; }
  
  // AMP accuracy
  Int       getAMPAcc   ( UInt uiDepth ) { return m_iAMPAcc[uiDepth]; }
  Void      setAMPAcc   ( UInt uiDepth, Int iAccu ) { assert( uiDepth < g_uiMaxCUDepth);  m_iAMPAcc[uiDepth] = iAccu; }

  // Bit-depth
  UInt      getBitDepth     ()         { return m_uiBitDepth;     }
  Void      setBitDepth     ( UInt u ) { m_uiBitDepth = u;        }
  UInt      getBitIncrement ()         { return m_uiBitIncrement; }
  Void      setBitIncrement ( UInt u ) { m_uiBitIncrement = u;    }
#if H0736_AVC_STYLE_QP_RANGE
  Int       getQpBDOffsetY  ()             { return m_qpBDOffsetY;   }
  Void      setQpBDOffsetY  ( Int value  ) { m_qpBDOffsetY = value;  }
  Int       getQpBDOffsetC  ()             { return m_qpBDOffsetC;   }
  Void      setQpBDOffsetC  ( Int value  ) { m_qpBDOffsetC = value;  }
#endif

  Void      setLFCrossSliceBoundaryFlag     ( Bool   bValue  )    { m_bLFCrossSliceBoundaryFlag = bValue; }
  Bool      getLFCrossSliceBoundaryFlag     ()                    { return m_bLFCrossSliceBoundaryFlag;   } 

  Void setUseDF                   ( Bool b ) { m_useDF = b; }
  Bool getUseDF                   ()         { return m_useDF; }

  Void setUseSAO                  (Bool bVal)  {m_bUseSAO = bVal;}
  Bool getUseSAO                  ()           {return m_bUseSAO;}

  UInt      getMaxTLayers()                           { return m_uiMaxTLayers; }
  Void      setMaxTLayers( UInt uiMaxTLayers )        { assert( uiMaxTLayers <= MAX_TLAYER ); m_uiMaxTLayers = uiMaxTLayers; }

  Bool      getTemporalIdNestingFlag()                { return m_bTemporalIdNestingFlag; }
  Void      setTemporalIdNestingFlag( Bool bValue )   { m_bTemporalIdNestingFlag = bValue; }
  UInt      getPCMBitDepthLuma     ()         { return m_uiPCMBitDepthLuma;     }
  Void      setPCMBitDepthLuma     ( UInt u ) { m_uiPCMBitDepthLuma = u;        }
  UInt      getPCMBitDepthChroma   ()         { return m_uiPCMBitDepthChroma;   }
  Void      setPCMBitDepthChroma   ( UInt u ) { m_uiPCMBitDepthChroma = u;      }
  Void      setPCMFilterDisableFlag     ( Bool   bValue  )    { m_bPCMFilterDisableFlag = bValue; }
  Bool      getPCMFilterDisableFlag     ()                    { return m_bPCMFilterDisableFlag;   } 

  Void    setLFCrossTileBoundaryFlag               ( Bool   bValue  )    { m_bLFCrossTileBoundaryFlag = bValue; }
  Bool    getLFCrossTileBoundaryFlag               ()                    { return m_bLFCrossTileBoundaryFlag;   }
  Void     setUniformSpacingIdr             ( Int i )           { m_iUniformSpacingIdr = i; }
  Int      getUniformSpacingIdr             ()                  { return m_iUniformSpacingIdr; }
#if !REMOVE_TILE_DEPENDENCE
  Void     setTileBoundaryIndependenceIdr   ( Int i )           { m_iTileBoundaryIndependenceIdr = i; }
  Int      getTileBoundaryIndependenceIdr   ()                  { return m_iTileBoundaryIndependenceIdr; }
#endif
  Void     setNumColumnsMinus1              ( Int i )           { m_iNumColumnsMinus1 = i; }
  Int      getNumColumnsMinus1              ()                  { return m_iNumColumnsMinus1; }
  Void     setColumnWidth ( UInt* columnWidth )
  {
    if( m_iUniformSpacingIdr == 0 && m_iNumColumnsMinus1 > 0 )
    {
      m_puiColumnWidth = new UInt[ m_iNumColumnsMinus1 ];

      for(Int i=0; i<m_iNumColumnsMinus1; i++)
      {
        m_puiColumnWidth[i] = columnWidth[i];
     }
    }
  }
  UInt     getColumnWidth  (UInt columnIdx) { return *( m_puiColumnWidth + columnIdx ); }
  Void     setNumRowsMinus1( Int i )        { m_iNumRowsMinus1 = i; }
  Int      getNumRowsMinus1()               { return m_iNumRowsMinus1; }
  Void     setRowHeight    ( UInt* rowHeight )
  {
    if( m_iUniformSpacingIdr == 0 && m_iNumRowsMinus1 > 0 )
    {
      m_puiRowHeight = new UInt[ m_iNumRowsMinus1 ];

      for(Int i=0; i<m_iNumRowsMinus1; i++)
      {
        m_puiRowHeight[i] = rowHeight[i];
      }
    }
  }
  UInt     getRowHeight           (UInt rowIdx)    { return *( m_puiRowHeight + rowIdx ); }
  Bool getScalingListFlag       ()         { return m_scalingListEnabledFlag;     }
  Void setScalingListFlag       ( Bool b ) { m_scalingListEnabledFlag  = b;       }
#if H0567_DPB_PARAMETERS_PER_TEMPORAL_LAYER
  UInt getMaxDecPicBuffering  (UInt tlayer)            { return m_uiMaxDecPicBuffering[tlayer]; }
  Void setMaxDecPicBuffering  ( UInt ui, UInt tlayer ) { m_uiMaxDecPicBuffering[tlayer] = ui;   }
  UInt getMaxLatencyIncrease  (UInt tlayer)            { return m_uiMaxLatencyIncrease[tlayer];   }
  Void setMaxLatencyIncrease  ( UInt ui , UInt tlayer) { m_uiMaxLatencyIncrease[tlayer] = ui;      }
#else
  UInt getMaxDecFrameBuffering  ()            { return m_uiMaxDecFrameBuffering; }
  Void setMaxDecFrameBuffering  ( UInt ui )   { m_uiMaxDecFrameBuffering = ui;   }
  UInt getMaxLatencyIncrease    ()            { return m_uiMaxLatencyIncrease;   }
  Void setMaxLatencyIncrease    ( UInt ui )   { m_uiMaxLatencyIncrease= ui;      }
#endif
#if TILES_WPP_ENTRY_POINT_SIGNALLING
  UInt getTilesOrEntropyCodingSyncIdc ()                    { return m_tilesOrEntropyCodingSyncIdc;   }
  Void setTilesOrEntropyCodingSyncIdc ( UInt val )          { m_tilesOrEntropyCodingSyncIdc = val;    }
  Int  getNumSubstreams               ()                    { return m_numSubstreams;                 }
  Void setNumSubstreams               ( Int numSubstreams ) { m_numSubstreams = numSubstreams;        }
#endif
};

#if !RPS_IN_SPS
/// Reference Picture Set class
class TComReferencePictureSet
{
private:
  Int m_numberOfPictures;
  Int m_numberOfNegativePictures;
  Int m_numberOfPositivePictures;
  Int m_numberOfLongtermPictures;
  Int  m_deltaPOC[MAX_NUM_REF_PICS];
  Int  m_POC[MAX_NUM_REF_PICS];
  Bool m_used[MAX_NUM_REF_PICS];
  Bool m_interRPSPrediction;
  Int  m_deltaRIdxMinus1;   
  Int  m_deltaRPS; 
  Int  m_numRefIdc; 
  Int  m_refIdc[MAX_NUM_REF_PICS+1];

public:
  TComReferencePictureSet();
  virtual ~TComReferencePictureSet();

  Void setUsed(Int bufferNum, Bool used);
  Void setDeltaPOC(Int bufferNum, Int deltaPOC);
  Void setPOC(Int bufferNum, Int deltaPOC);
  Void setNumberOfPictures(Int numberOfPictures);

  Int  getUsed(Int bufferNum);
  Int  getDeltaPOC(Int bufferNum);
  Int  getPOC(Int bufferNum);
  Int  getNumberOfPictures();

  Void setNumberOfNegativePictures(Int number)  { m_numberOfNegativePictures = number; }
  Int  getNumberOfNegativePictures()            { return m_numberOfNegativePictures; }
  Void setNumberOfPositivePictures(Int number)  { m_numberOfPositivePictures = number; }
  Int  getNumberOfPositivePictures()            { return m_numberOfPositivePictures; }
  Void setNumberOfLongtermPictures(Int number)  { m_numberOfLongtermPictures = number; }
  Int  getNumberOfLongtermPictures()            { return m_numberOfLongtermPictures; }

  Void setInterRPSPrediction(Bool flag)         { m_interRPSPrediction = flag; }
  Bool getInterRPSPrediction()                  { return m_interRPSPrediction; }
  Void setDeltaRIdxMinus1(Int x)                { m_deltaRIdxMinus1 = x; }
  Int  getDeltaRIdxMinus1()                     { return m_deltaRIdxMinus1; }
  Void setDeltaRPS(Int x)                       { m_deltaRPS = x; }
  Int  getDeltaRPS()                            { return m_deltaRPS; }
  Void setNumRefIdc(Int x)                      { m_numRefIdc = x; }
  Int  getNumRefIdc()                           { return m_numRefIdc; }

  Void setRefIdc(Int bufferNum, Int refIdc);
  Int  getRefIdc(Int bufferNum);

  Void sortDeltaPOC();
  Void printDeltaPOC();
};

/// Reference Picture Set set class
class TComRPSList
{
private:
  Int  m_numberOfReferencePictureSets;
  TComReferencePictureSet* m_referencePictureSets;
  
public:
  TComRPSList();
  virtual ~TComRPSList();
  
  Void  create  (Int numberOfEntries);
  Void  destroy ();


  TComReferencePictureSet* getReferencePictureSet(Int referencePictureSetNum);
  Int getNumberOfReferencePictureSets();
  Void setNumberOfReferencePictureSets(Int numberOfReferencePictureSets);
};
#endif

/// Reference Picture Lists class
class TComRefPicListModification
{
private:
  UInt      m_bRefPicListModificationFlagL0;  
  UInt      m_bRefPicListModificationFlagL1;  
#if !H0137_0138_LIST_MODIFICATION
  UInt      m_uiNumberOfRefPicListModificationsL0;
  UInt      m_uiNumberOfRefPicListModificationsL1;
  UInt      m_ListIdcL0[32];
#endif
  UInt      m_RefPicSetIdxL0[32];
#if !H0137_0138_LIST_MODIFICATION
  UInt      m_ListIdcL1[32];
#endif
  UInt      m_RefPicSetIdxL1[32];
    
public:
  TComRefPicListModification();
  virtual ~TComRefPicListModification();
  
  Void  create                    ();
  Void  destroy                   ();

  Bool       getRefPicListModificationFlagL0() { return m_bRefPicListModificationFlagL0; }
  Void       setRefPicListModificationFlagL0(Bool flag) { m_bRefPicListModificationFlagL0 = flag; }
  Bool       getRefPicListModificationFlagL1() { return m_bRefPicListModificationFlagL1; }
  Void       setRefPicListModificationFlagL1(Bool flag) { m_bRefPicListModificationFlagL1 = flag; }
#if !H0137_0138_LIST_MODIFICATION
  UInt       getNumberOfRefPicListModificationsL0() { return m_uiNumberOfRefPicListModificationsL0; }
  Void       setNumberOfRefPicListModificationsL0(UInt nr) { m_uiNumberOfRefPicListModificationsL0 = nr; }
  UInt       getNumberOfRefPicListModificationsL1() { return m_uiNumberOfRefPicListModificationsL1; }
  Void       setNumberOfRefPicListModificationsL1(UInt nr) { m_uiNumberOfRefPicListModificationsL1 = nr; }
  Void       setListIdcL0(UInt idx, UInt idc) { m_ListIdcL0[idx] = idc; }
  UInt       getListIdcL0(UInt idx) { return m_ListIdcL0[idx]; }
#endif
  Void       setRefPicSetIdxL0(UInt idx, UInt refPicSetIdx) { m_RefPicSetIdxL0[idx] = refPicSetIdx; }
  UInt       getRefPicSetIdxL0(UInt idx) { return m_RefPicSetIdxL0[idx]; }
#if !H0137_0138_LIST_MODIFICATION
  Void       setListIdcL1(UInt idx, UInt idc) { m_ListIdcL1[idx] = idc; }
  UInt       getListIdcL1(UInt idx) { return m_ListIdcL1[idx]; }
#endif
  Void       setRefPicSetIdxL1(UInt idx, UInt refPicSetIdx) { m_RefPicSetIdxL1[idx] = refPicSetIdx; }
  UInt       getRefPicSetIdxL1(UInt idx) { return m_RefPicSetIdxL1[idx]; }
};

/// PPS class
class TComPPS
{
private:
  Int         m_PPSId;                    // pic_parameter_set_id
  Int         m_SPSId;                    // seq_parameter_set_id
  Int         m_picInitQPMinus26;
  Bool        m_useDQP;
  Bool        m_bConstrainedIntraPred;    // constrained_intra_pred_flag
 
  // access channel
  TComSPS*    m_pcSPS;
#if !RPS_IN_SPS
  TComRPSList* m_RPSList;
#endif
  UInt        m_uiMaxCuDQPDepth;
  UInt        m_uiMinCuDQPSize;

  Int        m_iChromaQpOffset;
  Int        m_iChromaQpOffset2nd;

  UInt        m_numRefIdxL0DefaultActive;
  UInt        m_numRefIdxL1DefaultActive;

#if !RPS_IN_SPS
  Bool        m_bLongTermRefsPresent;
#endif

#if !H0566_TLA
  UInt        m_uiNumTlayerSwitchingFlags;            // num_temporal_layer_switching_point_flags
  Bool        m_abTLayerSwitchingFlag[ MAX_TLAYER ];  // temporal_layer_switching_point_flag
#endif

  Int         m_iSliceGranularity;

  Bool        m_bUseWeightPred;           // Use of Weighting Prediction (P_SLICE)
  UInt        m_uiBiPredIdc;              // Use of Weighting Bi-Prediction (B_SLICE)

#if H0388
  Bool        m_OutputFlagPresentFlag;   // Indicates the presence of output_flag in slice header
#endif

  Int      m_iTileBehaviorControlPresentFlag;
  Bool     m_bLFCrossTileBoundaryFlag;
  Int      m_iColumnRowInfoPresent;
  Int      m_iUniformSpacingIdr;
#if !REMOVE_TILE_DEPENDENCE
  Int      m_iTileBoundaryIndependenceIdr;
#endif
  Int      m_iNumColumnsMinus1;
  UInt*    m_puiColumnWidth;
  Int      m_iNumRowsMinus1;
  UInt*    m_puiRowHeight;
  
  Int      m_iEntropyCodingMode; // !!! in PPS now, but also remains in slice header!
#if !WPP_SIMPLIFICATION  
  Int      m_iEntropyCodingSynchro;
  Bool     m_bCabacIstateReset;
#endif
  Int      m_iNumSubstreams;

  Bool     m_enableTMVPFlag;

#if MULTIBITS_DATA_HIDING
  Int      m_signHideFlag;
  Int      m_signHidingThreshold;
#endif

#if CABAC_INIT_FLAG
  Bool     m_cabacInitPresentFlag;
  UInt     m_encCABACTableIdx;           // Used to transmit table selection across slices
#endif
#if DBL_CONTROL
  Bool     m_DeblockingFilterControlPresent;
#endif
#if PARALLEL_MERGE
  UInt m_log2ParallelMergeLevelMinus2;
#endif
public:
  TComPPS();
  virtual ~TComPPS();
  
  Int       getPPSId ()      { return m_PPSId; }
  Void      setPPSId (Int i) { m_PPSId = i; }
  Int       getSPSId ()      { return m_SPSId; }
  Void      setSPSId (Int i) { m_SPSId = i; }
  
  Int       getSliceGranularity()        { return m_iSliceGranularity; }
  Void      setSliceGranularity( Int i ) { m_iSliceGranularity = i;    }
  Int       getPicInitQPMinus26 ()         { return  m_picInitQPMinus26; }
  Void      setPicInitQPMinus26 ( Int i )  { m_picInitQPMinus26 = i;     }
  Bool      getUseDQP ()                   { return m_useDQP;        }
  Void      setUseDQP ( Bool b )           { m_useDQP   = b;         }
  Bool      getConstrainedIntraPred ()         { return  m_bConstrainedIntraPred; }
  Void      setConstrainedIntraPred ( Bool b ) { m_bConstrainedIntraPred = b;     }

#if !H0566_TLA
  UInt      getNumTLayerSwitchingFlags()                                  { return m_uiNumTlayerSwitchingFlags; }
  Void      setNumTLayerSwitchingFlags( UInt uiNumTlayerSwitchingFlags )  { assert( uiNumTlayerSwitchingFlags < MAX_TLAYER ); m_uiNumTlayerSwitchingFlags = uiNumTlayerSwitchingFlags; }

  Bool      getTLayerSwitchingFlag( UInt uiTLayer )                       { assert( uiTLayer < MAX_TLAYER ); return m_abTLayerSwitchingFlag[ uiTLayer ]; }
  Void      setTLayerSwitchingFlag( UInt uiTLayer, Bool bValue )          { m_abTLayerSwitchingFlag[ uiTLayer ] = bValue; }
#endif

#if !RPS_IN_SPS
  Bool      getLongTermRefsPresent()         { return m_bLongTermRefsPresent; }
  Void      setLongTermRefsPresent(Bool b)   { m_bLongTermRefsPresent=b;      }
#endif
  Void      setSPS              ( TComSPS* pcSPS ) { m_pcSPS = pcSPS; }
  TComSPS*  getSPS              ()         { return m_pcSPS;          }
#if !RPS_IN_SPS
  Void      setRPSList          ( TComRPSList* RPSList ) { m_RPSList = RPSList; }
  TComRPSList* getRPSList       ()         { return m_RPSList;        }
#endif
  Void      setMaxCuDQPDepth    ( UInt u ) { m_uiMaxCuDQPDepth = u;   }
  UInt      getMaxCuDQPDepth    ()         { return m_uiMaxCuDQPDepth;}
  Void      setMinCuDQPSize     ( UInt u ) { m_uiMinCuDQPSize = u;    }
  UInt      getMinCuDQPSize     ()         { return m_uiMinCuDQPSize; }

  Void      setChromaQpOffset   ( Int i ) { m_iChromaQpOffset = i; }
  Int       getChromaQpOffset   () { return m_iChromaQpOffset;}
  Void      setChromaQpOffset2nd( Int i ) { m_iChromaQpOffset2nd = i; }
  Int       getChromaQpOffset2nd() { return m_iChromaQpOffset2nd;}

  Void      setNumRefIdxL0DefaultActive(UInt ui)    { m_numRefIdxL0DefaultActive=ui;     }
  UInt      getNumRefIdxL0DefaultActive()           { return m_numRefIdxL0DefaultActive; }
  Void      setNumRefIdxL1DefaultActive(UInt ui)    { m_numRefIdxL1DefaultActive=ui;     }
  UInt      getNumRefIdxL1DefaultActive()           { return m_numRefIdxL1DefaultActive; }

  Bool getUseWP                     ()          { return m_bUseWeightPred;  }
  UInt getWPBiPredIdc               ()          { return m_uiBiPredIdc;     }

  Void setUseWP                     ( Bool b )  { m_bUseWeightPred = b;     }
  Void setWPBiPredIdc               ( UInt u )  { m_uiBiPredIdc = u;        }

#if H0388
  Void      setOutputFlagPresentFlag( Bool b )  { m_OutputFlagPresentFlag = b;    }
  Bool      getOutputFlagPresentFlag()          { return m_OutputFlagPresentFlag; }
#endif 

  Void    setTileBehaviorControlPresentFlag        ( Int i )             { m_iTileBehaviorControlPresentFlag = i;    }
  Int     getTileBehaviorControlPresentFlag        ()                    { return m_iTileBehaviorControlPresentFlag; }
  Void    setLFCrossTileBoundaryFlag               ( Bool   bValue  )    { m_bLFCrossTileBoundaryFlag = bValue; }
  Bool    getLFCrossTileBoundaryFlag               ()                    { return m_bLFCrossTileBoundaryFlag;   }
  Void     setColumnRowInfoPresent          ( Int i )           { m_iColumnRowInfoPresent = i; }
  Int      getColumnRowInfoPresent          ()                  { return m_iColumnRowInfoPresent; }
  Void     setUniformSpacingIdr             ( Int i )           { m_iUniformSpacingIdr = i; }
  Int      getUniformSpacingIdr             ()                  { return m_iUniformSpacingIdr; }
#if !REMOVE_TILE_DEPENDENCE
  Void     setTileBoundaryIndependenceIdr   ( Int i )           { m_iTileBoundaryIndependenceIdr = i; }
  Int      getTileBoundaryIndependenceIdr   ()                  { return m_iTileBoundaryIndependenceIdr; }
#endif
  Void     setNumColumnsMinus1              ( Int i )           { m_iNumColumnsMinus1 = i; }
  Int      getNumColumnsMinus1              ()                  { return m_iNumColumnsMinus1; }
  Void     setColumnWidth ( UInt* columnWidth )
  {
    if( m_iUniformSpacingIdr == 0 && m_iNumColumnsMinus1 > 0 )
    {
      m_puiColumnWidth = new UInt[ m_iNumColumnsMinus1 ];

      for(Int i=0; i<m_iNumColumnsMinus1; i++)
      {
        m_puiColumnWidth[i] = columnWidth[i];
      }
    }
  }
  UInt     getColumnWidth  (UInt columnIdx) { return *( m_puiColumnWidth + columnIdx ); }
  Void     setNumRowsMinus1( Int i )        { m_iNumRowsMinus1 = i; }
  Int      getNumRowsMinus1()               { return m_iNumRowsMinus1; }
  Void     setRowHeight    ( UInt* rowHeight )
  {
    if( m_iUniformSpacingIdr == 0 && m_iNumRowsMinus1 > 0 )
    {
      m_puiRowHeight = new UInt[ m_iNumRowsMinus1 ];

      for(Int i=0; i<m_iNumRowsMinus1; i++)
      {
        m_puiRowHeight[i] = rowHeight[i];
      }
    }
  }
  UInt     getRowHeight           (UInt rowIdx)    { return *( m_puiRowHeight + rowIdx ); }
  Void     setEntropyCodingMode(Int iEntropyCodingMode)       { m_iEntropyCodingMode = iEntropyCodingMode; }
  Int      getEntropyCodingMode()                             { return m_iEntropyCodingMode; }
#if !WPP_SIMPLIFICATION
  Void     setEntropyCodingSynchro(Int iEntropyCodingSynchro) { m_iEntropyCodingSynchro = iEntropyCodingSynchro; }
  Int      getEntropyCodingSynchro()                          { return m_iEntropyCodingSynchro; }
  Void     setCabacIstateReset(Bool bCabacIstateReset)        { m_bCabacIstateReset = bCabacIstateReset; }
  Bool     getCabacIstateReset()                              { return m_bCabacIstateReset; }
#endif
  Void     setNumSubstreams(Int iNumSubstreams)               { m_iNumSubstreams = iNumSubstreams; }
  Int      getNumSubstreams()                                 { return m_iNumSubstreams; }

#if MULTIBITS_DATA_HIDING
  Void      setSignHideFlag( Int signHideFlag ) { m_signHideFlag = signHideFlag; }
  Void      setTSIG( Int tsig )                 { m_signHidingThreshold = tsig; }
  Int       getSignHideFlag()                    { return m_signHideFlag; }
  Int       getTSIG()                            { return m_signHidingThreshold; }
#endif

  Void     setEnableTMVPFlag( Bool b )  { m_enableTMVPFlag = b;    }
  Bool     getEnableTMVPFlag()          { return m_enableTMVPFlag; }

#if CABAC_INIT_FLAG
  Void     setCabacInitPresentFlag( Bool flag )     { m_cabacInitPresentFlag = flag;    }
  Void     setEncCABACTableIdx( Int idx )           { m_encCABACTableIdx = idx;         }
  Bool     getCabacInitPresentFlag()                { return m_cabacInitPresentFlag;    }
  UInt     getEncCABACTableIdx()                    { return m_encCABACTableIdx;        }
#endif
#if DBL_CONTROL
  Void setDeblockingFilterControlPresent    ( Bool bValue )       { m_DeblockingFilterControlPresent = bValue; }
  Bool getDeblockingFilterControlPresent    ()                    { return m_DeblockingFilterControlPresent; }
#endif
#if PARALLEL_MERGE
  UInt getLog2ParallelMergeLevelMinus2      ()                    { return m_log2ParallelMergeLevelMinus2; }
  Void setLog2ParallelMergeLevelMinus2      (UInt mrgLevel)       { m_log2ParallelMergeLevelMinus2 = mrgLevel; }
#endif
};

/// SCALING_LIST class
class TComScalingList
{
public:
  TComScalingList();
  virtual ~TComScalingList();
  Void     setScalingListPresentFlag    (Bool b)                               { m_scalingListPresentFlag = b;    }
  Bool     getScalingListPresentFlag    ()                                     { return m_scalingListPresentFlag; }
  Int*     getScalingListAddress          (UInt sizeId, UInt listId)           { return m_scalingListCoef[sizeId][listId]; } //!< get matrix coefficient
  Bool     checkPredMode                  (UInt sizeId, UInt listId);
  Void     setRefMatrixId                 (UInt sizeId, UInt listId, UInt u)   { m_refMatrixId[sizeId][listId] = u;    }     //!< set reference matrix ID
  UInt     getRefMatrixId                 (UInt sizeId, UInt listId)           { return m_refMatrixId[sizeId][listId]; }     //!< get reference matrix ID
  Int*     getScalingListDefaultAddress   (UInt sizeId, UInt listId);                                                        //!< get default matrix coefficient
  Void     processDefaultMarix            (UInt sizeId, UInt listId);
#if SCALING_LIST
  Void     setScalingListDC               (UInt sizeId, UInt listId, UInt u)   { m_scalingListDC[sizeId][listId] = u; }      //!< set DC value
  Int      getScalingListDC               (UInt sizeId, UInt listId)           { return m_scalingListDC[sizeId][listId]; }   //!< get DC value
  Void     checkDcOfMatrix                ();
  Void     setUseDefaultScalingMatrixFlag (UInt sizeId, UInt listId, Bool b)   { m_useDefaultScalingMatrixFlag[sizeId][listId] = b;    } //!< set default matrix enabled/disabled in each matrix
  Bool     getUseDefaultScalingMatrixFlag (UInt sizeId, UInt listId)           { return m_useDefaultScalingMatrixFlag[sizeId][listId]; } //!< get default matrix enabled/disabled in each matrix
#endif
  Void     processRefMatrix               (UInt sizeId, UInt listId , UInt refListId );
  Bool     xParseScalingList              (char* pchFile);

private:
  Void     init                    ();
  Void     destroy                 ();
#if SCALING_LIST
  Int      m_scalingListDC               [SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM]; //!< the DC value of the matrix coefficient for 16x16
  Bool     m_useDefaultScalingMatrixFlag [SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM]; //!< UseDefaultScalingMatrixFlag
#endif
  UInt     m_refMatrixId                 [SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM]; //!< RefMatrixID
  Bool     m_scalingListPresentFlag;                                                //!< flag for using default matrix
  UInt     m_predMatrixId                [SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM]; //!< reference list index
  Int      *m_scalingListCoef            [SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM]; //!< quantization matrix
};

/// APS class
class TComAPS
{
public:
  TComAPS();
  virtual ~TComAPS();

  Void      setAPSID      (Int iID)   {m_apsID = iID;            }  //!< set APS ID 
  Int       getAPSID      ()          {return m_apsID;           }  //!< get APS ID
  Void      setSaoEnabled (Bool bVal) {m_bSaoEnabled = bVal;     }  //!< set SAO enabled/disabled in APS
  Bool      getSaoEnabled ()          {return m_bSaoEnabled;     }  //!< get SAO enabled/disabled in APS
  Void      setAlfEnabled (Bool bVal) {m_bAlfEnabled = bVal;     }  //!< set ALF enabled/disabled in APS
  Bool      getAlfEnabled ()          {return m_bAlfEnabled;     }  //!< get ALF enabled/disabled in APS

#if LCU_SYNTAX_ALF
  AlfParamSet* getAlfParam   ()          {return m_alfParamSet;}
#else
  ALFParam* getAlfParam   ()          {return m_pAlfParam;       }  //!< get ALF parameters in APS
#endif
  SAOParam* getSaoParam   ()          {return m_pSaoParam;       }  //!< get SAO parameters in APS

  Void      createSaoParam();   //!< create SAO parameter object
  Void      destroySaoParam();  //!< destroy SAO parameter object

  Void      createAlfParam();   //!< create ALF parameter object
  Void      destroyAlfParam();  //!< destroy ALF parameter object

  Void      setLoopFilterOffsetInAPS(Bool val)  {m_loopFilterOffsetInAPS = val; }      //!< set offset for deblocking filter enabled/disabled in APS
  Bool      getLoopFilterOffsetInAPS()          {return m_loopFilterOffsetInAPS; }     //!< get offset for deblocking filter enabled/disabled in APS
  Void      setLoopFilterDisable(Bool val)      {m_loopFilterDisable = val; }          //!< set offset for deblocking filter disabled
  Bool      getLoopFilterDisable()              {return m_loopFilterDisable; }         //!< get offset for deblocking filter disabled
  Void      setLoopFilterBetaOffset(Int val)    {m_loopFilterBetaOffsetDiv2 = val; }    //!< set beta offset for deblocking filter
  Int       getLoopFilterBetaOffset()           {return m_loopFilterBetaOffsetDiv2; }   //!< get beta offset for deblocking filter
  Void      setLoopFilterTcOffset(Int val)      {m_loopFilterTcOffsetDiv2 = val; }      //!< set tc offset for deblocking filter
  Int       getLoopFilterTcOffset()             {return m_loopFilterTcOffsetDiv2; }     //!< get tc offset for deblocking filter

  Void      createScalingList();
  Void      destroyScalingList();
  Void      setScalingListEnabled (Bool bVal) { m_scalingListEnabled = bVal; }  //!< set ScalingList enabled/disabled in APS
  Bool      getScalingListEnabled ()          { return m_scalingListEnabled; }  //!< get ScalingList enabled/disabled in APS
  TComScalingList* getScalingList ()          { return m_scalingList; }         //!< get ScalingList class pointer in APS
#if SAO_UNIT_INTERLEAVING
  Bool     getSaoInterleavingFlag() {return m_saoInterleavingFlag;}             //!< get SAO interleaving flag in APS
  Void     setSaoInterleavingFlag(Bool bVal) {m_saoInterleavingFlag = bVal;}    //!< set SAO interleaving flag in APS
#endif

private:
  Int         m_apsID;        //!< APS ID
  Bool        m_bSaoEnabled;  //!< SAO enabled/disabled in APS (true for enabled)
  Bool        m_bAlfEnabled;  //!< ALF enabled/disabled in APS (true for enabled)
  SAOParam*   m_pSaoParam;    //!< SAO parameter object pointer 
#if LCU_SYNTAX_ALF
  AlfParamSet*   m_alfParamSet;
#else
  ALFParam*   m_pAlfParam;    //!< ALF parameter object pointer
#endif
  Bool        m_loopFilterOffsetInAPS;       //< offset for deblocking filter in 0 = slice header, 1 = APS
  Bool        m_loopFilterDisable;           //< Deblocking filter enabled/disabled in APS
  Int         m_loopFilterBetaOffsetDiv2;    //< beta offset for deblocking filter
  Int         m_loopFilterTcOffsetDiv2;      //< tc offset for deblocking filter
  Bool        m_scalingListEnabled;     //!< ScalingList enabled/disabled in APS (true for enabled)
  TComScalingList*     m_scalingList;   //!< ScalingList class pointer
#if SAO_UNIT_INTERLEAVING
  Bool        m_saoInterleavingFlag;    //!< SAO interleaving flag
#endif

public:
  TComAPS& operator= (const TComAPS& src);  //!< "=" operator for APS object
};

typedef struct {
  // Explicit weighted prediction parameters parsed in slice header,
  // or Implicit weighted prediction parameters (8 bits depth values).
  Bool        bPresentFlag;
  UInt        uiLog2WeightDenom;
  Int         iWeight;
  Int         iOffset;

  // Weighted prediction scaling values built from above parameters (bitdepth scaled):
  Int         w, o, offset, shift, round;
} wpScalingParam;

typedef struct {
  Int64 iAC;
  Int64 iDC;
} wpACDCParam;

/// slice header class
class TComSlice
{
  
private:
  //  Bitstream writing
  Int         m_iAPSId; //!< APS ID in slice header
  bool       m_alfEnabledFlag;
  bool       m_saoEnabledFlag;
#if SAO_UNIT_INTERLEAVING
  bool       m_saoInterleavingFlag;   ///< SAO interleaving flag
  bool       m_saoEnabledFlagCb;      ///< SAO Cb enabled flag
  bool       m_saoEnabledFlagCr;      ///< SAO Cr enabled flag
#endif
  Int         m_iPPSId;               ///< picture parameter set ID
#if H0388
  Bool        m_PicOutputFlag;        ///< pic_output_flag 
#endif
  Int         m_iPOC;
  Int         m_iLastIDR;
  static Int  m_prevPOC;
  TComReferencePictureSet *m_pcRPS;
  TComReferencePictureSet m_LocalRPS;
  Int         m_iBDidx; 
  Int         m_iCombinationBDidx;
  Bool        m_bCombineWithReferenceFlag;
  TComRefPicListModification m_RefPicListModification;
  NalUnitType m_eNalUnitType;         ///< Nal unit type for the slice
  SliceType   m_eSliceType;
  Int         m_iSliceQp;
#if ADAPTIVE_QP_SELECTION
  Int         m_iSliceQpBase;
#endif
  Bool        m_bLoopFilterDisable;
  Bool        m_loopFilterOffsetInAPS;
  Bool        m_inheritDblParamFromAPS;      //< offsets for deblocking filter inherit from APS
  Int         m_loopFilterBetaOffsetDiv2;    //< beta offset for deblocking filter
  Int         m_loopFilterTcOffsetDiv2;      //< tc offset for deblocking filter
  
  Int         m_aiNumRefIdx   [3];    //  for multiple reference of current slice

  Int         m_iRefIdxOfLC[2][MAX_NUM_REF_LC];
  Int         m_eListIdFromIdxOfLC[MAX_NUM_REF_LC];
  Int         m_iRefIdxFromIdxOfLC[MAX_NUM_REF_LC];
  Int         m_iRefIdxOfL1FromRefIdxOfL0[MAX_NUM_REF_LC];
  Int         m_iRefIdxOfL0FromRefIdxOfL1[MAX_NUM_REF_LC];
  Bool        m_bRefPicListModificationFlagLC;
  Bool        m_bRefPicListCombinationFlag;

  Bool        m_bCheckLDC;

  //  Data
  Int         m_iSliceQpDelta;
  TComPic*    m_apcRefPicList [2][MAX_NUM_REF+1];
  Int         m_aiRefPOCList  [2][MAX_NUM_REF+1];
  Int         m_iDepth;
  
  // referenced slice?
  Bool        m_bRefenced;
  
  // access channel
  TComSPS*    m_pcSPS;
  TComPPS*    m_pcPPS;
  TComPic*    m_pcPic;
#if ADAPTIVE_QP_SELECTION
  TComTrQuant* m_pcTrQuant;
#endif  
  TComAPS*    m_pcAPS;  //!< pointer to APS parameter object

  UInt        m_uiColDir;  // direction to get colocated CUs
  
  UInt        m_colRefIdx;

#if ALF_CHROMA_LAMBDA || SAO_CHROMA_LAMBDA
  Double      m_dLambdaLuma;
  Double      m_dLambdaChroma;
#else
  Double      m_dLambda;
#endif

  Bool        m_abEqualRef  [2][MAX_NUM_REF][MAX_NUM_REF];
  
  Bool        m_bNoBackPredFlag;
  Bool        m_bRefIdxCombineCoding;

  UInt        m_uiTLayer;
  Bool        m_bTLayerSwitchingFlag;

  UInt        m_uiSliceMode;
  UInt        m_uiSliceArgument;
  UInt        m_uiSliceCurStartCUAddr;
  UInt        m_uiSliceCurEndCUAddr;
  UInt        m_uiSliceIdx;
  UInt        m_uiEntropySliceMode;
  UInt        m_uiEntropySliceArgument;
  UInt        m_uiEntropySliceCurStartCUAddr;
  UInt        m_uiEntropySliceCurEndCUAddr;
  Bool        m_bNextSlice;
  Bool        m_bNextEntropySlice;
  UInt        m_uiSliceBits;
  UInt        m_uiEntropySliceCounter;
  Bool        m_bFinalized;

  wpScalingParam  m_weightPredTable[2][MAX_NUM_REF][3]; // [REF_PIC_LIST_0 or REF_PIC_LIST_1][refIdx][0:Y, 1:U, 2:V]
  wpACDCParam    m_weightACDCParam[3];                 // [0:Y, 1:U, 2:V]
  wpScalingParam  m_weightPredTableLC[2*MAX_NUM_REF][3]; // [refIdxLC][0:Y, 1:U, 2:V]

  UInt        *m_uiTileByteLocation;
  UInt        m_uiTileCount;
  Int         m_iTileMarkerFlag;
  UInt        m_uiTileOffstForMultES;

  UInt*       m_puiSubstreamSizes;
  TComScalingList*     m_scalingList;                 //!< pointer of quantization matrix
#if CABAC_INIT_FLAG
  Bool        m_cabacInitFlag; 
#else
  Int         m_cabacInitIdc; 
#endif

#if H0111_MVD_L1_ZERO
  Bool       m_bLMvdL1Zero;
#endif
#if TILES_WPP_ENTRY_POINT_SIGNALLING
  Int         m_numEntryPointOffsets;
#endif

public:
  TComSlice();
  virtual ~TComSlice();
  
  Void      initSlice       ();
  Void      initTiles();

  
  Void      setSPS          ( TComSPS* pcSPS ) { m_pcSPS = pcSPS; }
  TComSPS*  getSPS          () { return m_pcSPS; }
  
  Void      setPPS          ( TComPPS* pcPPS )         { assert(pcPPS!=NULL); m_pcPPS = pcPPS; m_iPPSId = pcPPS->getPPSId(); }
  TComPPS*  getPPS          () { return m_pcPPS; }

#if ADAPTIVE_QP_SELECTION
  Void          setTrQuant          ( TComTrQuant* pcTrQuant ) { m_pcTrQuant = pcTrQuant; }
  TComTrQuant*  getTrQuant          () { return m_pcTrQuant; }
#endif

  Void      setPPSId        ( Int PPSId )         { m_iPPSId = PPSId; }
  Int       getPPSId        () { return m_iPPSId; }
  Void      setAPS          ( TComAPS* pcAPS ) { m_pcAPS = pcAPS; } //!< set APS pointer
  TComAPS*  getAPS          ()                 { return m_pcAPS;  } //!< get APS pointer
  Void      setAPSId        ( Int Id)          { m_iAPSId =Id;    } //!< set APS ID
  Int       getAPSId        ()                 { return m_iAPSId; } //!< get APS ID
#if H0388
  Void      setPicOutputFlag( Bool b )         { m_PicOutputFlag = b;    }
  Bool      getPicOutputFlag()                 { return m_PicOutputFlag; }
#endif
  Void      setAlfEnabledFlag(Bool s) {m_alfEnabledFlag =s; }
  Bool      getAlfEnabledFlag() { return m_alfEnabledFlag; }
  Void      setSaoEnabledFlag(Bool s) {m_saoEnabledFlag =s; }
  Bool      getSaoEnabledFlag() { return m_saoEnabledFlag; }
#if SAO_UNIT_INTERLEAVING
  Void      setSaoInterleavingFlag(Bool s) {m_saoInterleavingFlag =s; } //!< set SAO interleaving flag
  Bool      getSaoInterleavingFlag() { return m_saoInterleavingFlag;  } //!< get SAO interleaving flag
  Void      setSaoEnabledFlagCb(Bool s) {m_saoEnabledFlagCb =s; }       //!< set SAO Cb enabled flag
  Bool      getSaoEnabledFlagCb() { return m_saoEnabledFlagCb; }        //!< get SAO Cb enabled flag
  Void      setSaoEnabledFlagCr(Bool s) {m_saoEnabledFlagCr =s; }       //!< set SAO Cr enabled flag
  Bool      getSaoEnabledFlagCr() { return m_saoEnabledFlagCr; }        //!< get SAO Cr enabled flag
#endif
  Void      setRPS          ( TComReferencePictureSet *pcRPS ) { m_pcRPS = pcRPS; }
  TComReferencePictureSet*  getRPS          () { return m_pcRPS; }
  TComReferencePictureSet*  getLocalRPS     () { return &m_LocalRPS; }

  Void      setRPSidx          ( Int iBDidx ) { m_iBDidx = iBDidx; }
  Int       getRPSidx          () { return m_iBDidx; }
  Void      setCombinationBDidx          ( Int iCombinationBDidx ) { m_iCombinationBDidx = iCombinationBDidx; }
  Int       getCombinationBDidx          () { return m_iCombinationBDidx; }
  Void      setCombineWithReferenceFlag          ( Bool bCombineWithReferenceFlag ) { m_bCombineWithReferenceFlag = bCombineWithReferenceFlag; }
  Bool      getCombineWithReferenceFlag          () { return m_bCombineWithReferenceFlag; }
  Int       getPrevPOC      ()                          { return  m_prevPOC;       }

  TComRefPicListModification* getRefPicListModification() { return &m_RefPicListModification; }
  Void      setLastIDR(Int iIDRPOC)                       { m_iLastIDR = iIDRPOC; }
  Int       getLastIDR()                                  { return m_iLastIDR; }
  SliceType getSliceType    ()                          { return  m_eSliceType;         }
  Int       getPOC          ()                          { return  m_iPOC;           }
  Int       getSliceQp      ()                          { return  m_iSliceQp;           }
#if ADAPTIVE_QP_SELECTION
  Int       getSliceQpBase  ()                          { return  m_iSliceQpBase;       }
#endif
  Int       getSliceQpDelta ()                          { return  m_iSliceQpDelta;      }
  Bool      getLoopFilterDisable()                      { return  m_bLoopFilterDisable; }
  Bool      getLoopFilterOffsetInAPS()                  { return  m_loopFilterOffsetInAPS;}
  Bool      getInheritDblParamFromAPS()                 { return  m_inheritDblParamFromAPS; }
  Int       getLoopFilterBetaOffset()                   { return  m_loopFilterBetaOffsetDiv2; }
  Int       getLoopFilterTcOffset()                     { return  m_loopFilterTcOffsetDiv2; }

  Int       getNumRefIdx        ( RefPicList e )                { return  m_aiNumRefIdx[e];             }
  TComPic*  getPic              ()                              { return  m_pcPic;                      }
  TComPic*  getRefPic           ( RefPicList e, Int iRefIdx)    { return  m_apcRefPicList[e][iRefIdx];  }
  Int       getRefPOC           ( RefPicList e, Int iRefIdx)    { return  m_aiRefPOCList[e][iRefIdx];   }
  Int       getDepth            ()                              { return  m_iDepth;                     }
  UInt      getColDir           ()                              { return  m_uiColDir;                   }
  Bool      getColRefIdx        ()                              { return  m_colRefIdx;                  }
  Void      checkColRefIdx      (UInt curSliceIdx, TComPic* pic);
  Bool      getCheckLDC     ()                                  { return m_bCheckLDC; }
#if H0111_MVD_L1_ZERO
  Bool      getMvdL1ZeroFlag ()                                  { return m_bLMvdL1Zero;    }
#endif
#if H0137_0138_LIST_MODIFICATION
  Int       getNumRpsCurrTempList();
#endif
  Int       getRefIdxOfLC       (RefPicList e, Int iRefIdx)     { return m_iRefIdxOfLC[e][iRefIdx];           }
  Int       getListIdFromIdxOfLC(Int iRefIdx)                   { return m_eListIdFromIdxOfLC[iRefIdx];       }
  Int       getRefIdxFromIdxOfLC(Int iRefIdx)                   { return m_iRefIdxFromIdxOfLC[iRefIdx];       }
  Int       getRefIdxOfL0FromRefIdxOfL1(Int iRefIdx)            { return m_iRefIdxOfL0FromRefIdxOfL1[iRefIdx];}
  Int       getRefIdxOfL1FromRefIdxOfL0(Int iRefIdx)            { return m_iRefIdxOfL1FromRefIdxOfL0[iRefIdx];}
  Bool      getRefPicListModificationFlagLC()                   {return m_bRefPicListModificationFlagLC;}
  Void      setRefPicListModificationFlagLC(Bool bflag)         {m_bRefPicListModificationFlagLC=bflag;}     
  Bool      getRefPicListCombinationFlag()                      {return m_bRefPicListCombinationFlag;}
  Void      setRefPicListCombinationFlag(Bool bflag)            {m_bRefPicListCombinationFlag=bflag;}     
  Void      setListIdFromIdxOfLC(Int  iRefIdx, UInt uiVal)      { m_eListIdFromIdxOfLC[iRefIdx]=uiVal; }
  Void      setRefIdxFromIdxOfLC(Int  iRefIdx, UInt uiVal)      { m_iRefIdxFromIdxOfLC[iRefIdx]=uiVal; }
  Void      setRefIdxOfLC       (RefPicList e, Int iRefIdx, Int RefIdxLC)     { m_iRefIdxOfLC[e][iRefIdx]=RefIdxLC;}

  Void      setReferenced(Bool b)                               { m_bRefenced = b; }
  Bool      isReferenced()                                      { return m_bRefenced; }
  
  Void      setPOC              ( Int i )                       { m_iPOC              = i; if(getTLayer()==0) m_prevPOC=i; }
  Void      setNalUnitType      ( NalUnitType e )               { m_eNalUnitType      = e;      }
  NalUnitType getNalUnitType    ()                              { return m_eNalUnitType;        }
  Void      checkCRA(TComReferencePictureSet *pReferencePictureSet, Int& pocCRA, TComList<TComPic*>& rcListPic);
  Void      decodingRefreshMarking(Int& pocCRA, Bool& bRefreshPending, TComList<TComPic*>& rcListPic);
  Void      setSliceType        ( SliceType e )                 { m_eSliceType        = e;      }
  Void      setSliceQp          ( Int i )                       { m_iSliceQp          = i;      }
#if ADAPTIVE_QP_SELECTION
  Void      setSliceQpBase      ( Int i )                       { m_iSliceQpBase      = i;      }
#endif
  Void      setSliceQpDelta     ( Int i )                       { m_iSliceQpDelta     = i;      }
  Void      setLoopFilterDisable( Bool b )                      { m_bLoopFilterDisable= b;      }
  Void      setLoopFilterOffsetInAPS( Bool b )                  { m_loopFilterOffsetInAPS = b;}
  Void      setInheritDblParamFromAPS( Bool b )                 { m_inheritDblParamFromAPS = b; }
  Void      setLoopFilterBetaOffset( Int i )                    { m_loopFilterBetaOffsetDiv2 = i; }
  Void      setLoopFilterTcOffset( Int i )                      { m_loopFilterTcOffsetDiv2 = i; }
  
  Void      setRefPic           ( TComPic* p, RefPicList e, Int iRefIdx ) { m_apcRefPicList[e][iRefIdx] = p; }
  Void      setRefPOC           ( Int i, RefPicList e, Int iRefIdx ) { m_aiRefPOCList[e][iRefIdx] = i; }
  Void      setNumRefIdx        ( RefPicList e, Int i )         { m_aiNumRefIdx[e]    = i;      }
  Void      setPic              ( TComPic* p )                  { m_pcPic             = p;      }
  Void      setDepth            ( Int iDepth )                  { m_iDepth            = iDepth; }
  
  Void      setRefPicList       ( TComList<TComPic*>& rcListPic );
  Void      setRefPOCList       ();
  Void      setColDir           ( UInt uiDir ) { m_uiColDir = uiDir; }
  Void      setColRefIdx        ( UInt refIdx) { m_colRefIdx = refIdx; }
  Void      setCheckLDC         ( Bool b )                      { m_bCheckLDC = b; }
#if H0111_MVD_L1_ZERO
  Void      setMvdL1ZeroFlag     ( Bool b)                       { m_bLMvdL1Zero = b; }
#endif  

  Bool      isIntra         ()                          { return  m_eSliceType == I_SLICE;  }
  Bool      isInterB        ()                          { return  m_eSliceType == B_SLICE;  }
  Bool      isInterP        ()                          { return  m_eSliceType == P_SLICE;  }
  
#if ALF_CHROMA_LAMBDA || SAO_CHROMA_LAMBDA  
  Void      setLambda( Double d, Double e ) { m_dLambdaLuma = d; m_dLambdaChroma = e;}
  Double    getLambdaLuma() { return m_dLambdaLuma;        }
  Double    getLambdaChroma() { return m_dLambdaChroma;        }
#else
  Void      setLambda( Double d ) { m_dLambda = d; }
  Double    getLambda() { return m_dLambda;        }
#endif
  
  Void      initEqualRef();
  Bool      isEqualRef  ( RefPicList e, Int iRefIdx1, Int iRefIdx2 )
  {
    if (iRefIdx1 < 0 || iRefIdx2 < 0) return false;
    return m_abEqualRef[e][iRefIdx1][iRefIdx2];
  }
  
  Void setEqualRef( RefPicList e, Int iRefIdx1, Int iRefIdx2, Bool b)
  {
    m_abEqualRef[e][iRefIdx1][iRefIdx2] = m_abEqualRef[e][iRefIdx2][iRefIdx1] = b;
  }
  
  static Void      sortPicList         ( TComList<TComPic*>& rcListPic );
  
  Bool getNoBackPredFlag() { return m_bNoBackPredFlag; }
  Void setNoBackPredFlag( Bool b ) { m_bNoBackPredFlag = b; }
  Bool getRefIdxCombineCoding() { return m_bRefIdxCombineCoding; }
  Void setRefIdxCombineCoding( Bool b ) { m_bRefIdxCombineCoding = b; }
  Void generateCombinedList       ();

  UInt getTLayer             ()                            { return m_uiTLayer;                      }
  Void setTLayer             ( UInt uiTLayer )             { m_uiTLayer = uiTLayer;                  }

#if !H0566_TLA
  Bool getTLayerSwitchingFlag()                            { return m_bTLayerSwitchingFlag;          }
  Void setTLayerSwitchingFlag( Bool bValue )               { m_bTLayerSwitchingFlag = bValue;        }
#endif

  Void setTLayerInfo( UInt uiTLayer );
  Void decodingMarking( TComList<TComPic*>& rcListPic, Int iGOPSIze, Int& iMaxRefPicNum ); 
  Void applyReferencePictureSet( TComList<TComPic*>& rcListPic, TComReferencePictureSet *RPSList);
#if H0566_TLA && H0566_TLA_SET_FOR_SWITCHING_POINTS
  Bool isTemporalLayerSwitchingPoint( TComList<TComPic*>& rcListPic, TComReferencePictureSet *RPSList);
#endif
#if START_DECODING_AT_CRA
  Int       checkThatAllRefPicsAreAvailable( TComList<TComPic*>& rcListPic, TComReferencePictureSet *pReferencePictureSet, Bool printErrors, Int pocRandomAccess = 0);
#else
  Int       checkThatAllRefPicsAreAvailable( TComList<TComPic*>& rcListPic, TComReferencePictureSet *pReferencePictureSet, Bool printErrors);
#endif
  Void      createExplicitReferencePictureSetFromReference( TComList<TComPic*>& rcListPic, TComReferencePictureSet *pReferencePictureSet);

  Void decodingMarkingForNoTMVP( TComList<TComPic*>& rcListPic, Int currentPOC );

  UInt m_uiMaxNumMergeCand;
  Void setMaxNumMergeCand               (UInt maxNumMergeCand ) { m_uiMaxNumMergeCand = maxNumMergeCand;  }
  UInt getMaxNumMergeCand               ()                  {return m_uiMaxNumMergeCand;                  }

  Void setSliceMode                     ( UInt uiMode )     { m_uiSliceMode = uiMode;                     }
  UInt getSliceMode                     ()                  { return m_uiSliceMode;                       }
  Void setSliceArgument                 ( UInt uiArgument ) { m_uiSliceArgument = uiArgument;             }
  UInt getSliceArgument                 ()                  { return m_uiSliceArgument;                   }
  Void setSliceCurStartCUAddr           ( UInt uiAddr )     { m_uiSliceCurStartCUAddr = uiAddr;           }
  UInt getSliceCurStartCUAddr           ()                  { return m_uiSliceCurStartCUAddr;             }
  Void setSliceCurEndCUAddr             ( UInt uiAddr )     { m_uiSliceCurEndCUAddr = uiAddr;             }
  UInt getSliceCurEndCUAddr             ()                  { return m_uiSliceCurEndCUAddr;               }
  Void setSliceIdx                      ( UInt i)           { m_uiSliceIdx = i;                           }
  UInt getSliceIdx                      ()                  { return  m_uiSliceIdx;                       }
  Void copySliceInfo                    (TComSlice *pcSliceSrc);
  Void setEntropySliceMode              ( UInt uiMode )     { m_uiEntropySliceMode = uiMode;              }
  UInt getEntropySliceMode              ()                  { return m_uiEntropySliceMode;                }
  Void setEntropySliceArgument          ( UInt uiArgument ) { m_uiEntropySliceArgument = uiArgument;      }
  UInt getEntropySliceArgument          ()                  { return m_uiEntropySliceArgument;            }
  Void setEntropySliceCurStartCUAddr    ( UInt uiAddr )     { m_uiEntropySliceCurStartCUAddr = uiAddr;    }
  UInt getEntropySliceCurStartCUAddr    ()                  { return m_uiEntropySliceCurStartCUAddr;      }
  Void setEntropySliceCurEndCUAddr      ( UInt uiAddr )     { m_uiEntropySliceCurEndCUAddr = uiAddr;      }
  UInt getEntropySliceCurEndCUAddr      ()                  { return m_uiEntropySliceCurEndCUAddr;        }
  Void setNextSlice                     ( Bool b )          { m_bNextSlice = b;                           }
  Bool isNextSlice                      ()                  { return m_bNextSlice;                        }
  Void setNextEntropySlice              ( Bool b )          { m_bNextEntropySlice = b;                    }
  Bool isNextEntropySlice               ()                  { return m_bNextEntropySlice;                 }
  Void setSliceBits                     ( UInt uiVal )      { m_uiSliceBits = uiVal;                      }
  UInt getSliceBits                     ()                  { return m_uiSliceBits;                       }  
  Void setEntropySliceCounter           ( UInt uiVal )      { m_uiEntropySliceCounter = uiVal;            }
  UInt getEntropySliceCounter           ()                  { return m_uiEntropySliceCounter;             }
  Void setFinalized                     ( Bool uiVal )      { m_bFinalized = uiVal;                       }
  Bool getFinalized                     ()                  { return m_bFinalized;                        }
  Void  setWpScaling    ( wpScalingParam  wp[2][MAX_NUM_REF][3] ) { memcpy(m_weightPredTable, wp, sizeof(wpScalingParam)*2*MAX_NUM_REF*3); }
  Void  getWpScaling    ( RefPicList e, Int iRefIdx, wpScalingParam *&wp);

  Void  resetWpScaling  (wpScalingParam  wp[2][MAX_NUM_REF][3]);
  Void  initWpScaling    (wpScalingParam  wp[2][MAX_NUM_REF][3]);
  Void  initWpScaling   ();
  inline Bool applyWP   () { return( (m_eSliceType==P_SLICE && m_pcPPS->getUseWP()) || (m_eSliceType==B_SLICE && m_pcPPS->getWPBiPredIdc()) ); }
  
  Void  setWpAcDcParam  ( wpACDCParam wp[3] ) { memcpy(m_weightACDCParam, wp, sizeof(wpACDCParam)*3); }
  Void  getWpAcDcParam  ( wpACDCParam *&wp );
  Void  initWpAcDcParam ();
  Void  copyWPtable     (wpScalingParam *&wp_src, wpScalingParam *&wp_dst);
  Void  getWpScalingLC  ( Int iRefIdx, wpScalingParam *&wp);
  Void  resetWpScalingLC(wpScalingParam  wp[2*MAX_NUM_REF][3]);
  Void  setWpParamforLC();
  Void setTileLocationCount             ( UInt uiCount )      { m_uiTileCount = uiCount;                  }
  UInt getTileLocationCount             ()                    { return m_uiTileCount;                     }
  Void setTileLocation                  ( Int i, UInt uiLOC ) { m_uiTileByteLocation[i] = uiLOC;          }
  UInt getTileLocation                  ( Int i )             { return m_uiTileByteLocation[i];           }
  Void setTileMarkerFlag                ( Int iFlag )         { m_iTileMarkerFlag = iFlag;                }
  Int  getTileMarkerFlag                ()                    { return m_iTileMarkerFlag;                 }
  Void setTileOffstForMultES            (UInt uiOffset )      { m_uiTileOffstForMultES = uiOffset;        }
  UInt getTileOffstForMultES            ()                    { return m_uiTileOffstForMultES;            }
  Void allocSubstreamSizes              ( UInt uiNumSubstreams );
  UInt* getSubstreamSizes               ()                  { return m_puiSubstreamSizes; }
  Void  setScalingList              ( TComScalingList* scalingList ) { m_scalingList = scalingList; }
  TComScalingList*   getScalingList ()                               { return m_scalingList; }
  Void  setDefaultScalingList       ();
  Bool  checkDefaultScalingList     ();
#if CABAC_INIT_FLAG
  Void      setCabacInitFlag  ( Bool val ) { m_cabacInitFlag = val;      }  //!< set CABAC initial flag 
  Bool      getCabacInitFlag  ()           { return m_cabacInitFlag;     }  //!< get CABAC initial flag 
#else
  Void      setCABACinitIDC(Int iVal) {m_cabacInitIdc = iVal;    }  //!< set CABAC initial IDC number 
  Int       getCABACinitIDC()         {return m_cabacInitIdc;    }  //!< get CABAC initial IDC number 
#endif
#if TILES_WPP_ENTRY_POINT_SIGNALLING
  Void      setNumEntryPointOffsets(Int val)  { m_numEntryPointOffsets = val;     }
  Int       getNumEntryPointOffsets()         { return m_numEntryPointOffsets;    }
#endif

protected:
  TComPic*  xGetRefPic  (TComList<TComPic*>& rcListPic,
                         UInt                uiPOC);
  TComPic*  xGetLongTermRefPic  (TComList<TComPic*>& rcListPic,
                         UInt                uiPOC);
};// END CLASS DEFINITION TComSlice


template <class T> class ParameterSetMap
{
public:
  ParameterSetMap(Int maxId)
  :m_maxId (maxId)
  {}

  ~ParameterSetMap()
  {
    for (typename std::map<Int,T *>::iterator i = m_paramsetMap.begin(); i!= m_paramsetMap.end(); i++)
    {
      delete (*i).second;
    }
  }

  Void storePS(Int psId, T *ps)
  {
    assert ( psId < m_maxId );
    if ( m_paramsetMap.find(psId) != m_paramsetMap.end() )
    {
      delete m_paramsetMap[psId];
    }
    m_paramsetMap[psId] = ps; 
  }

  Void mergePSList(ParameterSetMap<T> &rPsList)
  {
    for (typename std::map<Int,T *>::iterator i = rPsList.m_paramsetMap.begin(); i!= rPsList.m_paramsetMap.end(); i++)
    {
      storePS(i->first, i->second);
    }
    rPsList.m_paramsetMap.clear();
  }


  T* getPS(Int psId)
  {
    return ( m_paramsetMap.find(psId) == m_paramsetMap.end() ) ? NULL : m_paramsetMap[psId];
  }

  T* getFirstPS()
  {
    return (m_paramsetMap.begin() == m_paramsetMap.end() ) ? NULL : m_paramsetMap.begin()->second;
  }

private:
  std::map<Int,T *> m_paramsetMap;
  Int               m_maxId;
};

class ParameterSetManager
{
public:
  ParameterSetManager();
  virtual ~ParameterSetManager();

  //! store sequence parameter set and take ownership of it 
  Void storeSPS(TComSPS *sps) { m_spsMap.storePS( sps->getSPSId(), sps); };
  //! get pointer to existing sequence parameter set  
  TComSPS* getSPS(Int spsId)  { return m_spsMap.getPS(spsId); };
  TComSPS* getFirstSPS()      { return m_spsMap.getFirstPS(); };

  //! store picture parameter set and take ownership of it 
  Void storePPS(TComPPS *pps) { m_ppsMap.storePS( pps->getPPSId(), pps); };
  //! get pointer to existing picture parameter set  
  TComPPS* getPPS(Int ppsId)  { return m_ppsMap.getPS(ppsId); };
  TComPPS* getFirstPPS()      { return m_ppsMap.getFirstPS(); };

  //! store adaptation parameter set and take ownership of it 
  Void storeAPS(TComAPS *aps) { m_apsMap.storePS( aps->getAPSID(), aps); };
  //! getPointer to existing adaptation parameter set  
  TComAPS* getAPS(Int apsId)  { return m_apsMap.getPS(apsId); };

protected:
  ParameterSetMap<TComSPS> m_spsMap; 
  ParameterSetMap<TComPPS> m_ppsMap; 
  ParameterSetMap<TComAPS> m_apsMap; 
};

//! \}

#endif // __TCOMSLICE__
