/**
 * @file gcbayes.h
 * Header file for gcbayes
 *
 * @author Jayanth Chennamangalam
 * @date 2009.12.09
 */

#ifndef __GCBAYES_H__
#define __GCBAYES_H__

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <signal.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#include <cpgplot.h>
#include "colourmap.h"

/**
 * @defgroup SPSBool Standard SPS booleans.
 */
/* @{ */
#define FALSE        0   /**< @brief Boolean false */
#define TRUE         1   /**< @brief Boolean true */
/* @} */

#define LEN_GENSTRING       256 /**< @brief Length of a generic string */

#define DEF_MIN_S_MIN       0.001
#define DEF_STEP_S_MIN      0.001

#define DEF_MIN_N           0
#define DEF_STEP_N          1
#define DEF_MAX_N           1000
#define DEF_N               1000

#define DEF_MIN_LMEAN       -1.19   /* Boyles et al. 2011 */
#define DEF_MAX_LMEAN       -1.04   /* Boyles et al. 2011 */
#define DEF_MIN_SD          0.91    /* Boyles et al. 2011 */
#define DEF_MAX_SD          0.98    /* Boyles et al. 2011 */

#define DEF_STEP_MEAN       0.01
#define DEF_STEP_SD         0.01

#define DEF_LEN_DIST        50

#define DEF_CONF_INT        95.45   /* 2-sigma */

#define SQRT_2              ((double) 1.4142135623730950488)

#define PG_DEV              "1/XS"  /**< @brief Device for screen plotting */
#define PG_DEV_PS_GRID      "gcbayes.ps/CPS"/**< @brief Device suffix for
                                                 plotting to a
                                                 PostScript file */
#define PG_VP_ML            0.17    /**< @brief Left margin */
#define PG_VP_MR            0.83    /**< @brief Right margin */
#define PG_VP_MB            0.17    /**< @brief Bottom margin */
#define PG_VP_MT            0.83    /**< @brief Top margin */
#define PG_SYMBOL           2
#define PG_CI_DEF           1
#define PG_CI_PLOT          11
#define PG_HBIN_LIM         200 
#define PG_TICK_STEPS_X     10      /**< @brief tick marks on the x-axis */
#define PG_TICK_STEPS_Y     5       /**< @brief tick marks on the y-axis */

typedef struct tagConf
{
    char acGCName[LEN_GENSTRING];
    char acFileFlux[LEN_GENSTRING];
    int in;
    float fSMin;    /* read from the above file, for now */
    float fDist;
    float fDistSD;
    float fSObs;
    float fMinSMean;
    float fMaxSMean;
} GC_CONF;

typedef struct tagStats
{
    int   iErrCode;
    float fMean;
    float fMode;
    float fMedian;
    float fMedianP;
    float fMinCI;
    float fMaxCI;
} STATS;

int RegisterSignalHandlers(void);
void HandleStopSignals(int iSigNo);
unsigned long int GetFreeMem(void);
STATS GetStats(float* pfX, float* pfP, float fStepX, int iLen,
               int iPlotLumAll);
int DoGrid(int iStartN, int iStepN, int iMaxN,
           float fStartLMean, float fStepLMean, float fMaxLMean,
           float fStartSD, float fStepSD, float fMaxSD,
           float fStartSMin, float fStepSMin,
           char acFileConf[],
           int iColourMap, int iNeedPS, int iPlotLumAll);
GC_CONF ReadGCConf(char* pcFileConf);
void CleanUp(void);
void PrintUsage(const char *pcProgName);
float CalcPObs(float fSMean, float fSD, float fSMin);
float CalcS1Lhd(int in, float *pfFlux,
                float fp_obs, float fSMean, float fSD);
float gammln(float xx);
float factln(int n);
float CalcS2Lhd(int in, float fp_obs, int iN);
float CalcS3Lhd(float fSObs, int iN, float fSMean, float fSD);
int GetIndex(float* pfArray, int iLen, float fEpsilon, float fVal);

#endif  /* __GCBAYES_H__ */

