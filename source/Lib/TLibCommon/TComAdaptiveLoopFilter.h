/* ====================================================================================================================

  The copyright in this software is being made available under the License included below.
  This software may be subject to other third party and   contributor rights, including patent rights, and no such
  rights are granted under this license.

  Copyright (c) 2010, SAMSUNG ELECTRONICS CO., LTD. and BRITISH BROADCASTING CORPORATION
  All rights reserved.

  Redistribution and use in source and binary forms, with or without modification, are permitted only for
  the purpose of developing standards within the Joint Collaborative Team on Video Coding and for testing and
  promoting such standards. The following conditions are required to be met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions and
      the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and
      the following disclaimer in the documentation and/or other materials provided with the distribution.
    * Neither the name of SAMSUNG ELECTRONICS CO., LTD. nor the name of the BRITISH BROADCASTING CORPORATION
      may be used to endorse or promote products derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
  INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 * ====================================================================================================================
*/

/** \file     TComAdaptiveLoopFilter.h
    \brief    adaptive loop filter class (header)
*/

#ifndef __TCOMADAPTIVELOOPFILTER__
#define __TCOMADAPTIVELOOPFILTER__

#include "TComPic.h"

// ====================================================================================================================
// Constants
// ====================================================================================================================

#define ALF_MAX_NUM_TAP       9                                       ///< maximum number of filter taps (9x9)
#define ALF_MIN_NUM_TAP       5                                       ///< minimum number of filter taps (5x5)
#define ALF_MAX_NUM_TAP_C     5                                       ///< number of filter taps for chroma (5x5)
#define ALF_MAX_NUM_COEF      42                                      ///< maximum number of filter coefficients
#define ALF_MIN_NUM_COEF      14                                      ///< minimum number of filter coefficients
#define ALF_MAX_NUM_COEF_C    14                                      ///< number of filter taps for chroma
#define ALF_NUM_BIT_SHIFT     8                                       ///< bit shift parameter for quantization of ALF param.
#define ALF_ROUND_OFFSET      ( 1 << ( ALF_NUM_BIT_SHIFT - 1 ) )      ///< rounding offset for ALF quantization

#include "../TLibCommon/CommonDef.h"

#define NUM_BITS               9
#define NO_TEST_FILT           3       // Filter supports (5/7/9)
#define NO_VAR_BINS           16 
#define NO_FILTERS            16
#define VAR_SIZE               3
#define FILTER_LENGTH          9

#define MAX_SQR_FILT_LENGTH   ((FILTER_LENGTH*FILTER_LENGTH) / 2 + 2)
#if TI_ALF_MAX_VSIZE_7
#define SQR_FILT_LENGTH_9SYM  ((9*9) / 4 + 2 - 1) 
#else
#define SQR_FILT_LENGTH_9SYM  ((9*9) / 4 + 2) 
#endif
#define SQR_FILT_LENGTH_7SYM  ((7*7) / 4 + 2) 
#define SQR_FILT_LENGTH_5SYM  ((5*5) / 4 + 2) 
#define MAX_SCAN_VAL    11
#define MAX_EXP_GOLOMB  16

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))
#define imgpel  unsigned short

#if TI_ALF_MAX_VSIZE_7
extern Int depthInt9x9Sym[21];
#else
extern Int depthInt9x9Sym[22];
#endif
extern Int depthInt7x7Sym[14];
extern Int depthInt5x5Sym[8];
extern Int *pDepthIntTab[NO_TEST_FILT];
void destroyMatrix_int(int **m2D);
void initMatrix_int(int ***m2D, int d1, int d2);

#if MTK_NONCROSS_INLOOP_FILTER
#define EXTEND_NUM_PEL    (UInt)(ALF_MAX_NUM_TAP/2)
#define EXTEND_NUM_PEL_C  (UInt)(ALF_MAX_NUM_TAP_C/2)
enum PaddingRegionPosition
{
  PRP_L = 0,
  PRP_R,
  PRP_T,
  PRP_B,
  NUM_PADDING_TILE,
  PRP_LT = NUM_PADDING_TILE,
  PRP_RT,
  PRP_LB,
  PRP_RB,
  NUM_PADDING_REGION
};
enum AlfChromaID
{
  ALF_Cb = 0,
  ALF_Cr = 1
};
#endif

// ====================================================================================================================
// Class definition
// ====================================================================================================================

#if MTK_NONCROSS_INLOOP_FILTER

class CAlfCU
{
public:
  CAlfCU() {}
  ~CAlfCU() {}
public:
  //  Void        init(TComDataCU* pcCU, UInt uiCUAddr, UInt uiStartCU, UInt uiEndCU, UInt uiNumCUWidth, UInt uiNumCUHeight);
  Void        init(TComPic* pcPic, UInt uiCUAddr, UInt uiStartCU, UInt uiEndCU, UInt uiNumCUWidth, UInt uiNumCUHeight);
  TComDataCU* getCU()           {return m_pcCU;}
  UInt        getWidth()        {return m_uiWidth;}
  UInt        getHeight()       {return m_uiHeight;}
  Int*        getCUBorderFlag() {return m_aiCUBorderFlag;}
  Void        extendCUBorder(Pel* pCUPel, UInt uiCUWidth, UInt uiCUHeight, Int iStride, UInt uiExtSize);
private:
  Void assignBorderStatus(UInt uiStartCU, UInt uiEndCU, UInt uiNumCUWidth, UInt uiNumCUHeight);

private:
  TComDataCU* m_pcCU;

  UInt        m_uiCUAddr;
  Int         m_aiCUBorderFlag[NUM_PADDING_REGION];
  UInt        m_uiWidth;
  UInt        m_uiHeight;
};

class CAlfSlice
{
public:
  CAlfSlice()
  {
    m_pcAlfCU= NULL;
  }
  ~CAlfSlice()
  {
    destroy(); 
  }
public:  //operator to access CAlfCU
  CAlfCU& operator[] (Int idx)
  {
    assert(idx < m_uiNumLCUs);
    return m_pcAlfCU[idx];
  }

public:
  Void init(UInt uiNumLCUsInPicWidth, UInt uiNumLCUsInPicHeight);
  Void create(TComPic* pcPic, Int iSliceID, UInt uiStartLCU, UInt uiEndLCU);
  Void destroy();

  UInt getNumLCUs     ()         {return m_uiNumLCUs;}

  Void extendSliceBorderLuma(Pel* pPelSrc, Int iStride, UInt uiExtSize);
  Void extendSliceBorderChroma(Pel* pPelSrc, Int iStride, UInt uiExtSize);
  Void copySliceLuma(Pel* pPicDst, Pel* pPicSrc, Int iStride);
  Void copySliceChroma(Pel* pPicDst, Pel* pPicSrc, Int iStride );

private: 

  Int    m_iSliceID;
  UInt   m_uiStartLCU;
  UInt   m_uiEndLCU;
  UInt   m_uiNumLCUs;

  CAlfCU* m_pcAlfCU;

  //----------------------------------------//
  UInt   m_uiNumLCUsInPicWidth;
  UInt   m_uiNumLCUsInPicHeight;
};

#endif 


/// adaptive loop filter class
class TComAdaptiveLoopFilter
{
protected:
  // quantized filter coefficients
  static const  Int m_aiSymmetricMag9x9[41];             ///< quantization scaling factor for 9x9 filter
  static const  Int m_aiSymmetricMag7x7[25];             ///< quantization scaling factor for 7x7 filter
  static const  Int m_aiSymmetricMag5x5[13];             ///< quantization scaling factor for 5x5 filter
#if TI_ALF_MAX_VSIZE_7
  static const  Int m_aiSymmetricMag9x7[32];             ///< quantization scaling factor for 9x7 filter
#endif
  
  // temporary picture buffer
  TComPicYuv*   m_pcTempPicYuv;                          ///< temporary picture buffer for ALF processing
  
  // ------------------------------------------------------------------------------------------------------------------
  // For luma component
  // ------------------------------------------------------------------------------------------------------------------
#if TI_ALF_MAX_VSIZE_7
  static Int m_pattern9x9Sym[39];
  static Int m_weights9x9Sym[21];
#else
  static Int m_pattern9x9Sym[41];
  static Int m_weights9x9Sym[22];
#endif
  static Int m_pattern9x9Sym_Quart[42];
  static Int m_pattern7x7Sym[25];
  static Int m_weights7x7Sym[14];
  static Int m_pattern7x7Sym_Quart[42];
  static Int m_pattern5x5Sym[13];
  static Int m_weights5x5Sym[8];
  static Int m_pattern5x5Sym_Quart[45];
#if TI_ALF_MAX_VSIZE_7
  static Int m_pattern9x9Sym_9[39];
#else
  static Int m_pattern9x9Sym_9[41];
#endif
  static Int m_pattern9x9Sym_7[25];
  static Int m_pattern9x9Sym_5[13];
  
  static Int *m_patternTab_filt[NO_TEST_FILT];
  static Int m_flTab[NO_TEST_FILT];
  static Int *m_patternTab[NO_TEST_FILT]; 
  static Int *m_patternMapTab[NO_TEST_FILT];
  static Int *m_weightsTab[NO_TEST_FILT];
  static Int m_sqrFiltLengthTab[NO_TEST_FILT];
  
  Int m_img_height,m_img_width;
  
  imgpel **m_imgY_pad;
  imgpel **m_imgY_var;
  Int    **m_imgY_temp;
  
  Int **m_filterCoeffSym;
  Int **m_filterCoeffPrevSelected;
  Int **m_filterCoeffTmp;
  Int **m_filterCoeffSymTmp;
  

#if MTK_NONCROSS_INLOOP_FILTER
  Bool        m_bUseNonCrossALF;
  UInt        m_uiNumLCUsInWidth;
  UInt        m_uiNumLCUsInHeight;
  UInt        m_uiNumSlicesInPic;
  CAlfSlice*  m_pSlice;
  Bool        m_bIsFirstDecodedSlice;

  Void xFilterOneSlice            (CAlfSlice* pSlice, imgpel* pDec, imgpel* pRest, Int iStride, ALFParam* pcAlfParam);
  Void calcVarforOneSlice         (CAlfSlice* pSlice, imgpel **imgY_var, imgpel *imgY_pad, Int pad_size, Int fl, Int img_stride);
  Void xFrameChromaforOneSlice    (CAlfSlice* pSlice, Int ComponentID, TComPicYuv* pcPicDec, TComPicYuv* pcPicRest, Int *qh, Int iTap);
  //for decoder CU on/off control
  Void setAlfCtrlFlagsforSlices   (ALFParam *pcAlfParam, UInt &idx);
  Void setAlfCtrlFlagsforOneSlice (CAlfSlice* pSlice, ALFParam *pcAlfParam, UInt &idx);
#endif


  /// ALF for luma component
  Void xALFLuma_qc( TComPic* pcPic, ALFParam* pcAlfParam, TComPicYuv* pcPicDec, TComPicYuv* pcPicRest );
  
  Void reconstructFilterCoeffs(ALFParam* pcAlfParam,int **pfilterCoeffSym, int bit_depth);
  Void getCurrentFilter(int **filterCoeffSym,ALFParam* pcAlfParam);
  // memory allocation
  Void destroyMatrix_imgpel(imgpel **m2D);
  Void destroyMatrix_int(int **m2D);
  Void initMatrix_int(int ***m2D, int d1, int d2);
  Void initMatrix_imgpel(imgpel ***m2D, int d1, int d2);
  Void destroyMatrix4D_double(double ****m4D, int d1, int d2);
  Void destroyMatrix3D_double(double ***m3D, int d1);
  Void destroyMatrix_double(double **m2D);
  Void initMatrix4D_double(double *****m4D, int d1, int d2, int d3, int d4);
  Void initMatrix3D_double(double ****m3D, int d1, int d2, int d3);
  Void initMatrix_double(double ***m2D, int d1, int d2);
  Void free_mem2Dpel(imgpel **array2D);
  Void get_mem2Dpel(imgpel ***array2D, int rows, int columns);
  Void no_mem_exit(const char *where);
  Void xError(const char *text, int code);
#if MTK_NONCROSS_INLOOP_FILTER
  Void calcVar(int ypos, int xpos, imgpel **imgY_var, imgpel *imgY_pad, int pad_size, int fl, int img_height, int img_width, int img_stride);
#else
  Void calcVar(imgpel **imgY_var, imgpel *imgY_pad, int pad_size, int fl, int img_height, int img_width, int img_stride);
#endif
  Void DecFilter_qc(imgpel* imgY_rec,ALFParam* pcAlfParam, int Stride);
  Void xSubCUAdaptive_qc(TComDataCU* pcCU, ALFParam* pcAlfParam, imgpel *imgY_rec_post, imgpel *imgY_rec, UInt uiAbsPartIdx, UInt uiDepth, Int Stride);
  Void xCUAdaptive_qc(TComPic* pcPic, ALFParam* pcAlfParam, imgpel *imgY_rec_post, imgpel *imgY_rec, Int Stride);
  Void subfilterFrame(imgpel *imgY_rec_post, imgpel *imgY_rec, int filtNo, int start_height, int end_height, int start_width, int end_width, int Stride);
  Void filterFrame(imgpel *imgY_rec_post, imgpel *imgY_rec, int filtNo, int Stride);
#if TSB_ALF_HEADER
  UInt  m_uiNumCUsInFrame;
  Void  setAlfCtrlFlags (ALFParam *pAlfParam, TComDataCU *pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt &idx);
#endif
  
  // ------------------------------------------------------------------------------------------------------------------
  // For chroma component
  // ------------------------------------------------------------------------------------------------------------------
  
  /// ALF for chroma component
  Void xALFChroma   ( ALFParam* pcAlfParam, TComPicYuv* pcPicDec, TComPicYuv* pcPicRest );
  
  /// sub function: non-adaptive ALF process for chroma
#if MTK_NONCROSS_INLOOP_FILTER
  Void xFrameChroma ( Int ypos, Int xpos, Int iHeight, Int iWidth, TComPicYuv* pcPicDec, TComPicYuv* pcPicRest, Int *qh, Int iTap, Int iColor );
#else
  Void xFrameChroma ( TComPicYuv* pcPicDec, TComPicYuv* pcPicRest, Int *qh, Int iTap, Int iColor );
#endif

public:
  TComAdaptiveLoopFilter();
  virtual ~TComAdaptiveLoopFilter() {}
  
  // initialize & destory temporary buffer
  Void create  ( Int iPicWidth, Int iPicHeight, UInt uiMaxCUWidth, UInt uiMaxCUHeight, UInt uiMaxCUDepth );
  Void destroy ();
  
  // alloc & free & set functions
  Void allocALFParam  ( ALFParam* pAlfParam );
  Void freeALFParam   ( ALFParam* pAlfParam );
  Void copyALFParam   ( ALFParam* pDesAlfParam, ALFParam* pSrcAlfParam );
#if TSB_ALF_HEADER
  Void  setNumCUsInFrame        (TComPic *pcPic);
#endif
  
  // predict filter coefficients
  Void predictALFCoeff        ( ALFParam* pAlfParam );                  ///< prediction of luma ALF coefficients
  Void predictALFCoeffChroma  ( ALFParam* pAlfParam );                  ///< prediction of chroma ALF coefficients
  
  // interface function
  Void ALFProcess             ( TComPic* pcPic, ALFParam* pcAlfParam ); ///< interface function for ALF process
  
#if TI_ALF_MAX_VSIZE_7
  static Int ALFTapHToTapV(Int tapH);
  static Int ALFTapHToNumCoeff(Int tapH);
  static Int ALFFlHToFlV(Int flH);
#endif


#if MTK_NONCROSS_INLOOP_FILTER
public:
  Void setNumSlicesInPic(UInt uiNum) {m_uiNumSlicesInPic = uiNum;}
  UInt getNumSlicesInPic()           {return m_uiNumSlicesInPic;}
  Void setUseNonCrossAlf(Bool bVal)  {m_bUseNonCrossALF = bVal;}
  Bool getUseNonCrossAlf()           {return m_bUseNonCrossALF;}
  Void createSlice      ();
  Void destroySlice     ();

public: //operator to access Alf slice
  CAlfSlice& operator[] (UInt i)
  {
    assert(i < m_uiNumSlicesInPic);
    return m_pSlice[i];
  }
#endif


};
#endif
