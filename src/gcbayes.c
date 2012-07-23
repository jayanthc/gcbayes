/*
 * @file gcbayes.c
 * Chennamangalam, Lorimer, Mandel, and Bagchi (2012)
 *
 * @author Jayanth Chennamangalam
 * @date 2009.12.09
 */

#include "gcbayes.h"

/* NOTE: all runtime-allocated memory pointer variables are global so that the
         SIGINT-handler can free them */
float* g_pfSMean = NULL;
float* g_pfLMean = NULL;
float* g_pfSD = NULL;
float* g_pfN = NULL;
float* g_pfLogN = NULL;
float* g_pfSMin = NULL;
float* g_pfS_i = NULL;
float* g_pfDist = NULL;
float* g_pfPDist = NULL;
float* g_pfSTotLhd = NULL;
float* g_pfLTotLhd = NULL;
float* g_pfLPost = NULL;
float* g_pfLPostMarg_MeanSD = NULL;
float* g_pfLPostMarg_N = NULL;
float* g_pfLPostMarg_Mean = NULL;
float* g_pfLPostMarg_SD = NULL;
float* g_pfLPostMarg_SMin = NULL;
float* g_pfLPostMarg_d = NULL;

int main(int argc, char *argv[])
{
    int iStartN = DEF_MIN_N;
    int iMaxN = DEF_MAX_N;
    int iStepN = DEF_STEP_N;
    /* mean in log10 */
    float fStartLMean = DEF_MIN_LMEAN;
    float fMaxLMean = DEF_MAX_LMEAN;
    float fStepLMean = DEF_STEP_MEAN;
    /* standard deviation in log10 */
    float fStartSD = DEF_MIN_SD;
    float fMaxSD = DEF_MAX_SD;
    float fStepSD = DEF_STEP_SD;
    float fStartSMin = DEF_MIN_S_MIN;
    float fStepSMin = DEF_STEP_S_MIN;
    int iColourMap = DEF_CMAP;
    int iNeedPS = FALSE;
    int iPlotAll = FALSE;
    int iRet = EXIT_SUCCESS;
    char acFileConf[LEN_GENSTRING] = {0};
    const char *pcProgName = NULL;
    int iNextOpt = 0;
    /* valid short options */
    const char* const pcOptsShort = "hN:S:X:m:p:u:s:e:a:y:c:fw";
    /* valid long options */
    const struct option stOptsLong[] = {
        { "help",       0, NULL, 'h' },
        { "n-min",      1, NULL, 'N' },
        { "n-step",     1, NULL, 'S' },
        { "n-max",      1, NULL, 'X' },
        { "mu-min",     1, NULL, 'm' },
        { "mu-step",    1, NULL, 'p' },
        { "mu-max",     1, NULL, 'u' },
        { "sd-min",     1, NULL, 's' },
        { "sd-step",    1, NULL, 'e' },
        { "sd-max",     1, NULL, 'a' },
        { "conf-file",  1, NULL, 'y' },
        { "colour-map", 1, NULL, 'c' },
        { "plot-ps",    0, NULL, 'f' },
        { "plot-all",   0, NULL, 'w' },
        { NULL,         0, NULL, 0   }
    };

    /* get the filename of the program from the argument list */
    pcProgName = argv[0];

    /* parse the input */
    do
    {
        iNextOpt = getopt_long(argc, argv, pcOptsShort, stOptsLong, NULL);
        switch (iNextOpt)
        {
            case 'h':   /* -h or --help */
                /* print usage info and terminate */
                PrintUsage(pcProgName);
                return EXIT_SUCCESS;

            case 'N':   /* -N or --n-min */
                /* set option */
                iStartN = (int) atoi(optarg);
                break;

            case 'S':   /* -S or --n-step */
                /* set option */
                iStepN = (int) atoi(optarg);
                break;

            case 'X':   /* -X or --n-max */
                /* set option */
                iMaxN = (int) atoi(optarg);
                break;

            case 'm':   /* -m or --mu-min */
                /* set option */
                fStartLMean = (float) atof(optarg);
                break;

            case 'p':   /* -p or --mu-step */
                /* set option */
                fStepLMean = (float) atof(optarg);
                break;

            case 'u':   /* -u or --mu-max */
                /* set option */
                fMaxLMean = (float) atof(optarg);
                break;

            case 's':   /* -s or --sd-min */
                /* set option */
                fStartSD = (float) atof(optarg);
                break;

            case 'e':   /* -e or --sd-step */
                /* set option */
                fStepSD = (float) atof(optarg);
                break;

            case 'a':   /* -a or --sd-max */
                /* set option */
                fMaxSD = (float) atof(optarg);
                break;

            case 'y':  /* -y or --conf-file */
                /* set option */
                (void) strcpy(acFileConf, optarg);
                break;

            case 'c':   /* -c or --colour-map */
                /* set option */
                iColourMap = GetColourMapFromName(optarg);
                break;

            case 'f':  /* -f or --plot-ps */
                /* set option */
                iNeedPS = TRUE;
                break;

            case 'w':  /* -w or --plot-all */
                /* set option */
                iPlotAll = TRUE;
                break;

            case '?':   /* user specified an invalid option */
                /* print usage info and terminate with error */
                (void) fprintf(stderr, "ERROR: Invalid option!\n");
                PrintUsage(pcProgName);
                return EXIT_FAILURE;

            case -1:    /* done with options */
                break;

            default:    /* unexpected */
                assert(0);
        }
    } while (iNextOpt != -1);

    /* register SIGTERM and CTRL+C signal handler */
    iRet = RegisterSignalHandlers();
    if (iRet != EXIT_SUCCESS)
    {
        (void) fprintf(stderr, "ERROR: Registering signal handler failed!\n");
        return EXIT_FAILURE;
    }

    iRet = DoGrid(iStartN, iStepN, iMaxN,
                  fStartLMean, fStepLMean, fMaxLMean,
                  fStartSD, fStepSD, fMaxSD,
                  fStartSMin, fStepSMin,
                  acFileConf,
                  iColourMap, iNeedPS, iPlotAll);
    if (iRet != EXIT_SUCCESS)
    {
        (void) fprintf(stderr,
                       "ERROR: Grid failed!\n");
        return EXIT_FAILURE;
    }

    /* free resources */
    CleanUp();

    return EXIT_SUCCESS;
}

int DoGrid(int iStartN, int iStepN, int iMaxN,
           float fStartLMean, float fStepLMean, float fMaxLMean,
           float fStartSD, float fStepSD, float fMaxSD,
           float fStartSMin, float fStepSMin,
           char acFileConf[],
           int iColourMap, int iNeedPS, int iPlotAll)
{
    int iRet = EXIT_SUCCESS;
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int m = 0;
    int n = 0;
    int in = 0;
    int iLenN = (int) ((iMaxN - iStartN) / iStepN);
    int iLenLMean = (int) ((fMaxLMean - fStartLMean) / fStepLMean);
    int iLenSD = (int) ((fMaxSD - fStartSD) / fStepSD);
    int iLenSMin = 0;
    int iRangeN = iMaxN - iStartN;
    float fRangeLMean = fMaxLMean - fStartLMean;
    float fRangeSD = fMaxSD - fStartSD;
    float fRangeSMean = 0.0;
    float fRangeSMin = 0.0;
    float fDist = 0.0;
    float fDistSD = 0.0;
    float fSObs = 0.0;
    float fStartSMean = 0.0;
    float fStepSMean = fStepLMean;
    float fMaxSMean = 0.0;
    int iLenSMean = 0;
    float fStartDist = 0.0;
    float fMaxDist = 0.0;
    float fStepDist = 0.0;
    int iLenDist = 0;
    float fp_obs = 0.0;
    float fDataMin = FLT_MAX;
    float fDataMax = -(FLT_MAX);
    float fMaxSMin = 0.0;
    float afTM[6] = {0.0};
    FILE *pFFlux = NULL;
    char acFileFlux[LEN_GENSTRING] = {0};
    GC_CONF stConf = {{0}};
    float fOpt = 0.0;
    unsigned long iFreeBytes = 0;
    unsigned long iBytesNeeded = 0;

    stConf = ReadGCConf(acFileConf);
    if ('\0' == stConf.acGCName[0])
    {
        (void) fprintf(stderr,
                       "ERROR: Reading config file failed!\n");
        return EXIT_FAILURE;
    }
    (void) strncpy(acFileFlux, stConf.acFileFlux, LEN_GENSTRING);
    in = stConf.in;
    fDist = stConf.fDist;
    fDistSD = stConf.fDistSD;
    /* user input validation */
    if (0.0 == fDistSD)
    {
        (void) fprintf(stderr, "ERROR: Distance uncertainty unspecified!\n");
        return EXIT_FAILURE;
    }
    fSObs = stConf.fSObs;
    fStartSMean = stConf.fMinSMean;
    fMaxSMean = stConf.fMaxSMean;

    g_pfN = (float *) malloc((size_t) iLenN * sizeof(float));
    if (NULL == g_pfN)
    {
        (void) fprintf(stderr,
                       "ERROR: %s %d: Memory allocation failed! %s\n",
                       __FILE__,
                       __LINE__,
                       strerror(errno));
        return EXIT_FAILURE;
    }
    for (i = 0; i < iLenN; ++i)
    {
        g_pfN[i] = (float) (iStartN + (i * iStepN));
    }

    iLenSMean = (int) ((fMaxSMean - fStartSMean) / fStepSMean);
    g_pfSMean = (float *) malloc((size_t) iLenSMean * sizeof(float));
    if (NULL == g_pfSMean)
    {
        (void) fprintf(stderr,
                       "ERROR: %s %d: Memory allocation failed! %s\n",
                       __FILE__,
                       __LINE__,
                       strerror(errno));
        return EXIT_FAILURE;
    }
    for (i = 0; i < iLenSMean; ++i)
    {
        g_pfSMean[i] = fStartSMean + (i * fStepSMean);
    }
    fRangeSMean = fMaxSMean - fStartSMean;

    g_pfSD = (float*) malloc((size_t) iLenSD * sizeof(float));
    if (NULL == g_pfSD)
    {
        (void) fprintf(stderr,
                       "ERROR: %s %d: Memory allocation failed! %s\n",
                       __FILE__,
                       __LINE__,
                       strerror(errno));
        return EXIT_FAILURE;
    }
    for (i = 0; i < iLenSD; ++i)
    {
        g_pfSD[i] = fStartSD + (i * fStepSD);
    }

    /* open the flux file */
    pFFlux = fopen(acFileFlux, "r");
    if (NULL == pFFlux)
    {
        (void) fprintf(stderr,
                       "ERROR: Opening flux file failed! %s\n",
                       strerror(errno));
        return EXIT_FAILURE;
    }

    /* allocate memory for the observed flux array */
    g_pfS_i = (float *) malloc((size_t) in * sizeof(float));
    if (NULL == g_pfS_i)
    {
        (void) fprintf(stderr,
                       "ERROR: %s %d: Memory allocation failed! %s\n",
                       __FILE__,
                       __LINE__,
                       strerror(errno));
        return EXIT_FAILURE;
    }

    i = 0;
    while (feof(pFFlux) == 0)
    {
        if (i >= in)
        {
            (void) fprintf(stderr,
                           "ERROR: Pulsar count mismatch!\n");
            /* close the flux file */
            (void) fclose(pFFlux);
            return EXIT_FAILURE;
        }
        (void) fscanf(pFFlux, "%f\n", &g_pfS_i[i]);
        ++i;
    }

    /* close the flux file */
    (void) fclose(pFFlux);

    /* NOTE: assumes file with fluxes sorted in ascending order */
    stConf.fSMin = g_pfS_i[0];
    fMaxSMin = stConf.fSMin;
    iLenSMin = (int) ((fMaxSMin - fStartSMin) / fStepSMin);
    g_pfSMin = (float *) malloc((size_t) iLenSMin * sizeof(float));
    if (NULL == g_pfSMin)
    {
        (void) fprintf(stderr,
                       "ERROR: %s %d: Memory allocation failed! %s\n",
                       __FILE__,
                       __LINE__,
                       strerror(errno));
        return EXIT_FAILURE;
    }
    for (i = 0; i < iLenSMin; ++i)
    {
        g_pfSMin[i] = fStartSMin + (i * fStepSMin);
    }
    fRangeSMin = fMaxSMin - fStartSMin;

    /* check if we have enough memory - estimate the memory that will be needed
       and compare it to the free memory (not including reclaimable memory) */
    /* NOTE: this is just the memory required for the luminosity domain
       likelihood, but the actual memory needed will be greater than this. */
    iBytesNeeded = (unsigned long) DEF_LEN_DIST * iLenLMean * iLenSD * iLenSMin
                   * iLenN * sizeof(float);
    iFreeBytes = GetFreeMem();
    if (iBytesNeeded > iFreeBytes)
    {
        (void) fprintf(stderr,
                       "ERROR: Insufficient memory. "
                       "Try reducing ranges and/or step sizes.\n");
        return EXIT_FAILURE;
    }

    /* total likelihood (flux domain) */
    g_pfSTotLhd = (float*) malloc((size_t) iLenSMean * iLenSD * iLenSMin * iLenN
                                  * sizeof(float));
    if (NULL == g_pfSTotLhd)
    {
        (void) fprintf(stderr,
                       "ERROR: %s %d: Memory allocation failed! %s\n",
                       __FILE__,
                       __LINE__,
                       strerror(errno));
        return EXIT_FAILURE;
    }

    /* compute the total likelihood (flux domain) */
    for (i = 0; i < iLenSD; ++i)
    {
        for (j = 0; j < iLenSMean; ++j)
        {
            for (k = 0; k < iLenSMin; ++k)
            {
                fp_obs = CalcPObs(g_pfSMean[j], g_pfSD[i], g_pfSMin[k]);
                float fS1Lhd = CalcS1Lhd(in, g_pfS_i,
                                          fp_obs,
                                          g_pfSMean[j],
                                          g_pfSD[i]);
                for (l = 0; l < iLenN; ++l)
                {
                    (g_pfSTotLhd
                     + (i * iLenSMean * iLenSMin * iLenN)
                     + (j * iLenSMin * iLenN)
                     + (k * iLenN))[l]
                        = fS1Lhd
                          * CalcS2Lhd(in, fp_obs, (int) g_pfN[l])
                          * CalcS3Lhd(fSObs,
                                      (int) g_pfN[l],
                                      g_pfSMean[j],
                                      g_pfSD[i]);
                }
            }
        }
    }

    /* reclaim memory */
    free(g_pfS_i);
    g_pfS_i = NULL;

    /* create the prior on r - a gaussian */
    fStartDist = fDist - (3 * fDistSD);
    fMaxDist = fDist + (3 * fDistSD);
    iLenDist = DEF_LEN_DIST;
    fStepDist = (fMaxDist - fStartDist) / (iLenDist - 1);
    g_pfDist = (float*) malloc((size_t) iLenDist * sizeof(float));
    if (NULL == g_pfDist)
    {
        (void) fprintf(stderr,
                       "ERROR: %s %d: Memory allocation failed! %s\n",
                       __FILE__,
                       __LINE__,
                       strerror(errno));
        return EXIT_FAILURE;
    }
    for (i = 0; i < iLenDist; ++i)
    {
        g_pfDist[i] = fStartDist + (i * fStepDist);
    }
    g_pfPDist = (float *) malloc((size_t) iLenDist * sizeof(float));
    if (NULL == g_pfPDist)
    {
        (void) fprintf(stderr,
                       "ERROR: %s %d: Memory allocation failed! %s\n",
                       __FILE__,
                       __LINE__,
                       strerror(errno));
        return EXIT_FAILURE;
    }
    for (i = 0; i < iLenDist; ++i)
    {
        g_pfPDist[i] = expf(-powf((g_pfDist[i] - fDist), 2)
                            / (2 * fDistSD * fDistSD))
                       / (sqrtf(2 * M_PI) * fDistSD);
    }

    /* total likelihood (luminosity domain) */
    g_pfLTotLhd = (float*) calloc((size_t) iLenLMean * iLenSD * iLenSMin * iLenN
                                  * iLenDist,
                                  sizeof(float));
    if (NULL == g_pfLTotLhd)
    {
        (void) fprintf(stderr,
                       "ERROR: %s %d: Memory allocation failed! %s\n",
                       __FILE__,
                       __LINE__,
                       strerror(errno));
        return EXIT_FAILURE;
    }

    /* create mu and sigma arrays, using priors */
    g_pfLMean = (float *) malloc((size_t) iLenLMean * sizeof(float));
    if (NULL == g_pfLMean)
    {
        (void) fprintf(stderr,
                       "ERROR: %s %d: Memory allocation failed! %s\n",
                       __FILE__,
                       __LINE__,
                       strerror(errno));
        return EXIT_FAILURE;
    }
    for (i = 0; i < iLenLMean; ++i)
    {
        g_pfLMean[i] = fStartLMean + (i * fStepLMean);
    }

    /* convert likelihood to luminosity domain */
    for (i = 0; i < iLenDist; ++i)
    {
        for (j = 0; j < iLenSD; ++j)
        {
            for (k = 0; k < iLenSMean; ++k)
            {
                float fLMean = g_pfSMean[k] + (2 * log10f(g_pfDist[i]));
                /* get index in g_pfLMean corresponding to fLMean */
                n = GetIndex(g_pfLMean, iLenLMean, fStepLMean, fLMean);
                if (-1 == n)  /* not in prior, skip to next */
                {
                    continue;
                }
                assert(n < iLenLMean);
                assert(n > -1);
                for (l = 0; l < iLenSMin; ++l)
                {
                    for (m = 0; m < iLenN; ++m)
                    {
                        (g_pfLTotLhd
                         + (i * iLenSD * iLenLMean * iLenSMin * iLenN)
                         + (j * iLenLMean * iLenSMin * iLenN)
                         + (n * iLenSMin * iLenN)   /* n instead of k */
                         + (l * iLenN))[m]
                            = (g_pfSTotLhd
                                + (j * iLenSMean * iLenSMin * iLenN)
                                + (k * iLenSMin * iLenN)
                                + (l * iLenN))[m];
                    }
                }
            }
        }
    }

    /* no need for flux domain likelihood anymore - reclaim memory */
    free(g_pfSMean);
    g_pfSMean = NULL;
    free(g_pfSTotLhd);
    g_pfSTotLhd = NULL;

    g_pfLPost = (float*) calloc((size_t) iLenLMean * iLenSD * iLenN * iLenSMin
                                * iLenDist,
                                sizeof(float));
    if (NULL == g_pfLPost)
    {
        (void) fprintf(stderr,
                       "ERROR: %s %d: Memory allocation failed! %s\n",
                       __FILE__,
                       __LINE__,
                       strerror(errno));
        return EXIT_FAILURE;
    }

    g_pfLPostMarg_MeanSD = (float*) calloc((size_t) iLenLMean * iLenSD, sizeof(float));
    if (NULL == g_pfLPostMarg_MeanSD)
    {
        (void) fprintf(stderr,
                       "ERROR: %s %d: Memory allocation failed! %s\n",
                       __FILE__,
                       __LINE__,
                       strerror(errno));
        return EXIT_FAILURE;
    }
    g_pfLPostMarg_N = (float*) calloc((size_t) iLenN, sizeof(float));
    if (NULL == g_pfLPostMarg_N)
    {
        (void) fprintf(stderr,
                       "ERROR: %s %d: Memory allocation failed! %s\n",
                       __FILE__,
                       __LINE__,
                       strerror(errno));
        return EXIT_FAILURE;
    }
    g_pfLPostMarg_Mean = (float*) calloc((size_t) iLenLMean, sizeof(float));
    if (NULL == g_pfLPostMarg_Mean)
    {
        (void) fprintf(stderr,
                       "ERROR: %s %d: Memory allocation failed! %s\n",
                       __FILE__,
                       __LINE__,
                       strerror(errno));
        return EXIT_FAILURE;
    }
    g_pfLPostMarg_SD = (float*) calloc((size_t) iLenSD, sizeof(float));
    if (NULL == g_pfLPostMarg_SD)
    {
        (void) fprintf(stderr,
                       "ERROR: %s %d: Memory allocation failed! %s\n",
                       __FILE__,
                       __LINE__,
                       strerror(errno));
        return EXIT_FAILURE;
    }
    g_pfLPostMarg_SMin = (float*) calloc((size_t) iLenSMin, sizeof(float));
    if (NULL == g_pfLPostMarg_SMin)
    {
        (void) fprintf(stderr,
                       "ERROR: %s %d: Memory allocation failed! %s\n",
                       __FILE__,
                       __LINE__,
                       strerror(errno));
        return EXIT_FAILURE;
    }
    g_pfLPostMarg_d = (float*) calloc((size_t) iLenDist, sizeof(float));
    if (NULL == g_pfLPostMarg_d)
    {
        (void) fprintf(stderr,
                       "ERROR: %s %d: Memory allocation failed! %s\n",
                       __FILE__,
                       __LINE__,
                       strerror(errno));
        return EXIT_FAILURE;
    }

    /* constant, but for completion */
    fOpt = (1.0 / iRangeN) * (1.0 / fRangeLMean)
           * (1.0 / fRangeSD) * (1.0 / fRangeSMin);
    /* compute posterior (luminosity domain) */
    for (i = 0; i < iLenDist; ++i)
    {
        for (j = 0; j < iLenSD; ++j)
        {
            for (k = 0; k < iLenLMean; ++k)
            {
                for (l = 0; l < iLenSMin; ++l)
                {
                    for (m = 0; m < iLenN; ++m)
                    {
                        (g_pfLPost
                         + (i * iLenSD * iLenLMean * iLenSMin * iLenN)
                         + (j * iLenLMean * iLenSMin * iLenN)
                         + (k * iLenSMin * iLenN)
                         + (l * iLenN))[m]
                            = ((g_pfLTotLhd
                                + (i * iLenSD * iLenLMean * iLenSMin * iLenN)
                                + (j * iLenLMean * iLenSMin * iLenN)
                                + (k * iLenSMin * iLenN)
                                + (l * iLenN))[m] * g_pfPDist[i] * fOpt);
                        /* marginalise the posterior */
                        (g_pfLPostMarg_MeanSD + j * iLenLMean)[k]
                            += (g_pfLPost
                                + (i * iLenSD * iLenLMean * iLenSMin * iLenN)
                                + (j * iLenLMean * iLenSMin * iLenN)
                                + (k * iLenSMin * iLenN)
                                + (l * iLenN))[m];
                        g_pfLPostMarg_N[m]
                            += (g_pfLPost
                                + (i * iLenSD * iLenLMean * iLenSMin * iLenN)
                                + (j * iLenLMean * iLenSMin * iLenN)
                                + (k * iLenSMin * iLenN)
                                + (l * iLenN))[m];
                        g_pfLPostMarg_Mean[k]
                            += (g_pfLPost
                                + (i * iLenSD * iLenLMean * iLenSMin * iLenN)
                                + (j * iLenLMean * iLenSMin * iLenN)
                                + (k * iLenSMin * iLenN)
                                + (l * iLenN))[m];
                        g_pfLPostMarg_SD[j]
                            += (g_pfLPost
                                + (i * iLenSD * iLenLMean * iLenSMin * iLenN)
                                + (j * iLenLMean * iLenSMin * iLenN)
                                + (k * iLenSMin * iLenN)
                                + (l * iLenN))[m];
                        g_pfLPostMarg_SMin[l]
                            += (g_pfLPost
                                + (i * iLenSD * iLenLMean * iLenSMin * iLenN)
                                + (j * iLenLMean * iLenSMin * iLenN)
                                + (k * iLenSMin * iLenN)
                                + (l * iLenN))[m];
                        g_pfLPostMarg_d[i]
                            += (g_pfLPost
                                + (i * iLenSD * iLenLMean * iLenSMin * iLenN)
                                + (j * iLenLMean * iLenSMin * iLenN)
                                + (k * iLenSMin * iLenN)
                                + (l * iLenN))[m];
                    }
                }
            }
        }
        g_pfLPostMarg_d[i] *= (fStepLMean * fStepSD * fStepSMin * iStepN);
    }

    /* for completion */
    for (j = 0; j < iLenSD; ++j)
    {
        for (k = 0; k < iLenLMean; ++k)
        {
            (g_pfLPostMarg_MeanSD + j * iLenLMean)[k]
                *= (fStepDist * fStepSMin * iStepN);
        }
    }
    for (i = 0; i < iLenSD; ++i)
    {
        g_pfLPostMarg_SD[i] *= (fStepDist * fStepLMean * fStepSMin * iStepN);
    }
    for (i = 0; i < iLenLMean; ++i)
    {
        g_pfLPostMarg_Mean[i] *= (fStepDist * fStepSD * fStepSMin * iStepN);
    }
    for (i = 0; i < iLenN; ++i)
    {
        g_pfLPostMarg_N[i] *= (fStepDist * fStepLMean * fStepSMin * fStepSD);
    }
    for (i = 0; i < iLenSMin; ++i)
    {
        g_pfLPostMarg_SMin[i] *= (fStepDist * fStepLMean * fStepSD * iStepN);
    }

    /* reclaim memory */
    free(g_pfLTotLhd);
    g_pfLTotLhd = NULL;
    free(g_pfLPost);
    g_pfLPost = NULL;
    free(g_pfPDist);
    g_pfPDist = NULL;

    /* open the PGPLOT graphics device */
    if (FALSE == iNeedPS)    /* plot to screen */
    {
        iRet = cpgopen(PG_DEV);
        if (iRet <= 0)
        {
            (void) fprintf(stderr,
                           "ERROR: Opening graphics device %s failed!\n",
                           PG_DEV);
            return EXIT_FAILURE;
        }
    }
    else                            /* plot to PS file */
    {
        iRet = cpgopen(PG_DEV_PS_GRID);
        if (iRet <= 0)
        {
            (void) fprintf(stderr,
                           "ERROR: Opening graphics device %s failed!\n",
                           PG_DEV_PS_GRID);
            return EXIT_FAILURE;
        }
    }

    if (iPlotAll)
    {
        cpgpap(0.0, 0.667);
    }
    else
    {
        cpgpap(0.0, 1.0);
    }

    cpgsch(1.8);
    if (iPlotAll)
    {
        cpgsubp(3, 2);
    }
    else
    {
        cpgsubp(1, 1);
    }

    if (iPlotAll)
    {
        /* plot marginalised posterior for mean and SD */
        {
            /* find the volume for normalisation */
            float fVol = 0.0;
            for (i = 0; i < iLenSD; ++i)
            {
                for (j = 0; j < iLenLMean; ++j)
                {
                    fVol += (g_pfLPostMarg_MeanSD + i * iLenLMean)[j];
                }
            }
            fVol *= (fStepSD * fStepLMean);
            for (i = 0; i < iLenSD; ++i)
            {
                for (j = 0; j < iLenLMean; ++j)
                {
                    (g_pfLPostMarg_MeanSD + i * iLenLMean)[j] /= fVol;
                }
            }

            /* find min. and max. for plotting */
            fDataMin = FLT_MAX;
            fDataMax = -(FLT_MAX);
            for (i = 0; i < iLenSD; ++i)
            {
                for (j = 0; j < iLenLMean; ++j)
                {
                    if ((g_pfLPostMarg_MeanSD + i * iLenLMean)[j] > fDataMax)
                    {
                        fDataMax = (g_pfLPostMarg_MeanSD + i * iLenLMean)[j];
                    }
                    if ((g_pfLPostMarg_MeanSD + i * iLenLMean)[j] < fDataMin)
                    {
                        fDataMin = (g_pfLPostMarg_MeanSD + i * iLenLMean)[j];
                    }
                }
            }
            SetColourMap(iColourMap, 0, fDataMin, fDataMax);

            cpgpanl(1, 1);
            cpgsvp(PG_VP_ML, PG_VP_MR, PG_VP_MB, PG_VP_MT);
            cpgswin(g_pfLMean[0], g_pfLMean[iLenLMean-1],
                    g_pfSD[0], g_pfSD[iLenSD-1]);
            afTM[1] = (g_pfLMean[iLenLMean-1] - g_pfLMean[0] + fStepLMean)
                      / iLenLMean;
            afTM[0] = g_pfLMean[0] - afTM[1];
            afTM[2] = 0;
            afTM[5] = (g_pfSD[iLenSD-1] - g_pfSD[0] + fStepSD)
                      / iLenSD;
            afTM[3] = g_pfSD[0] - afTM[5];
            afTM[4] = 0;
            cpgimag((const float *) g_pfLPostMarg_MeanSD,
                    iLenLMean,
                    iLenSD,
                    1,
                    iLenLMean,
                    1,
                    iLenSD,
                    fDataMin,
                    fDataMax,
                    afTM);
            cpgwedg("RI", 1.0, 5.0, fDataMin, fDataMax, "");
            cpgbox("CST", 0.0, 0, "CST", 0.0, 0);
            cpgsci(PG_CI_DEF);
            cpgmtxt("T", -1.5, 0.85, 0.0, "(a)");
            cpgaxis("N",
                    g_pfLMean[0], g_pfSD[0],
                    g_pfLMean[iLenLMean-1], g_pfSD[0],
                    g_pfLMean[0], g_pfLMean[iLenLMean-1],
                    0.0,
                    0,
                    0.0,
                    -0.4,
                    0.5,
                    0.8,
                    0);
            cpgaxis("N",
                    g_pfLMean[0], g_pfSD[0],
                    g_pfLMean[0], g_pfSD[iLenSD-1],
                    g_pfSD[0], g_pfSD[iLenSD-1],
                    0.0,
                    0,
                    -0.4,
                    0.0,
                    0.5,
                    -0.8,
                    0);
            cpglab("\\gm", "\\gs", "");

            /* plot FK06 point */
            //cpgsci(PG_CI_DEF);
            /* make the marker bold */
            //cpgslw(4);
            //cpgpt1(-1.1, 0.9, PG_SYMBOL);
            /* reset the line width */
            //cpgslw(1);
            cpgsci(PG_CI_DEF);

            /* reclaim memory */
            free(g_pfLPostMarg_MeanSD);
            g_pfLPostMarg_MeanSD = NULL;
            free(g_pfLPostMarg_MeanSD);
            g_pfLPostMarg_MeanSD = NULL;
        }
    }

    /* plot marginalised posterior for N */
    {
        /* find the area for normalisation (to calculate the mean, etc.) */
        float fArea = 0.0;
        for (i = 0; i < iLenN; ++i)
        {
            fArea += g_pfLPostMarg_N[i];
        }
        fArea *= iStepN;
        for (i = 0; i < iLenN; ++i)
        {
            g_pfLPostMarg_N[i] /= fArea;
        }

        /* find min. and max. for plotting */
        fDataMin = FLT_MAX;
        fDataMax = -(FLT_MAX);
        for (i = 0; i < iLenN; ++i)
        {
            if (g_pfLPostMarg_N[i] > fDataMax)
            {
                fDataMax = g_pfLPostMarg_N[i];
            }
            if (g_pfLPostMarg_N[i] < fDataMin)
            {
                fDataMin = g_pfLPostMarg_N[i];
            }
        }
        /* keep 10% space between the top of the box and the plot */
        fDataMax = (float) (fDataMax + ((fDataMax - fDataMin) * 0.1));

        if (iPlotAll)
        {
            cpgpanl(2, 1);
        }
        else
        {
            cpgpanl(1, 1);
        }
        cpgsvp(PG_VP_ML, PG_VP_MR, PG_VP_MB, PG_VP_MT);
        g_pfLogN = (float *) malloc((size_t) iLenN * sizeof(float));
        if (NULL == g_pfLogN)
        {
            (void) fprintf(stderr,
                           "ERROR: %s %d: Memory allocation failed! %s\n",
                           __FILE__,
                           __LINE__,
                           strerror(errno));
            return EXIT_FAILURE;
        }
        for (i = 0; i < iLenN; ++i)
        {
            g_pfLogN[i] = log10f(g_pfN[i]);
        }
        cpgswin(g_pfLogN[0], g_pfLogN[iLenN-1], fDataMin, fDataMax);
        if (iPlotAll)
        {
            cpglab("N",
                   "",
                   "");
        }
        else
        {
            cpglab("N",
                   "",
                   "");
        }
        cpgsci(PG_CI_DEF);
        cpgline(iLenN, g_pfLogN, g_pfLPostMarg_N);
        if (iPlotAll)
        {
            cpgmtxt("T", -1.5, 0.85, 0.0, "(b)");
        }
        else
        {
            cpgmtxt("T", -1.5, 0.85, 0.0, "(a)"); /* Ter 5 */
            //cpgmtxt("T", -1.5, 0.85, 0.0, "(b)"); /* 47 Tuc */
            //cpgmtxt("T", -1.5, 0.85, 0.0, "(c)");   /* M 28 */
        }

        STATS stStats = {0};
        stStats = GetStats(g_pfLogN, g_pfLPostMarg_N, (float) iStepN, iLenN,
                           iPlotAll);
        if (stStats.iErrCode != EXIT_SUCCESS)
        {
            (void) fprintf(stderr, "ERROR: Getting statistics failed!\n");
            return EXIT_FAILURE;
        }
        /* convert values */
        stStats.fMean = powf(10, stStats.fMean);
        stStats.fMode = powf(10, stStats.fMode);
        stStats.fMedian = powf(10, stStats.fMedian);
        stStats.fMinCI = powf(10, stStats.fMinCI);
        stStats.fMaxCI = powf(10, stStats.fMaxCI);
        (void) printf("Mean of N = %g\n", stStats.fMean);
        (void) printf("Mode of N = %g\n", stStats.fMode);
        (void) printf("Median of N = %g, P = %g\n",
                      stStats.fMedian,
                      stStats.fMedianP);
        (void) printf("CI Min = %g\n", stStats.fMinCI);
        (void) printf("CI Max = %g\n", stStats.fMaxCI);
        (void) printf("CI Min - Median = %g\n",
                      stStats.fMinCI - stStats.fMedian);
        (void) printf("CI Max - Median = %g\n",
                      stStats.fMaxCI - stStats.fMedian);

        /* draw box */
        cpgbox("BCLNST", 0.0, 0, "BCNST", 0.0, 0);

        /* reclaim memory */
        free(g_pfLogN);
        g_pfLogN = NULL;
        free(g_pfN);
        g_pfN = NULL;
        free(g_pfLPostMarg_N);
        g_pfLPostMarg_N = NULL;
        free(g_pfLPostMarg_N);
        g_pfLPostMarg_N = NULL;
    }

    if (iPlotAll)
    {
        /* plot marginalised posterior for mean */
        {
            /* find the area for normalisation (to calculate the mean, etc.) */
            float fArea = 0.0;
            for (i = 0; i < iLenLMean; ++i)
            {
                fArea += g_pfLPostMarg_Mean[i];
            }
            fArea *= fStepLMean;
            for (i = 0; i < iLenLMean; ++i)
            {
                g_pfLPostMarg_Mean[i] /= fArea;
            }

            /* find min. and max. for plotting */
            fDataMin = FLT_MAX;
            fDataMax = -(FLT_MAX);
            for (i = 0; i < iLenLMean; ++i)
            {
                if (g_pfLPostMarg_Mean[i] > fDataMax)
                {
                    fDataMax = g_pfLPostMarg_Mean[i];
                }
                if (g_pfLPostMarg_Mean[i] < fDataMin)
                {
                    fDataMin = g_pfLPostMarg_Mean[i];
                }
            }
            /* keep 10% space between the top of the box and the plot */
            fDataMax = (float) (fDataMax + ((fDataMax - fDataMin) * 0.1));

            cpgpanl(3, 1);
            cpgsvp(PG_VP_ML, PG_VP_MR, PG_VP_MB, PG_VP_MT);
            cpgswin(g_pfLMean[0], g_pfLMean[iLenLMean-1], fDataMin, fDataMax);
            cpglab("\\gm",
                   "",
                   "");
            cpgsci(PG_CI_DEF);
            cpgline(iLenLMean, g_pfLMean, g_pfLPostMarg_Mean);
            cpgmtxt("T", -1.5, 0.85, 0.0, "(c)");

            STATS stStats = {0};
            stStats = GetStats(g_pfLMean, g_pfLPostMarg_Mean, fStepLMean,
                               iLenLMean, iPlotAll);
            if (stStats.iErrCode != EXIT_SUCCESS)
            {
                (void) fprintf(stderr, "ERROR: Getting statistics failed!\n");
                return EXIT_FAILURE;
            }
            (void) printf("Mean of mu = %g\n", stStats.fMean);
            (void) printf("Mode of mu = %g\n", stStats.fMode);
            (void) printf("Median of mu = %g, P = %g\n",
                          stStats.fMedian,
                          stStats.fMedianP);
            (void) printf("CI Min = %g\n", stStats.fMinCI);
            (void) printf("CI Max = %g\n", stStats.fMaxCI);
            (void) printf("CI Min - Median = %g\n",
                          stStats.fMinCI - stStats.fMedian);
            (void) printf("CI Max - Median = %g\n",
                          stStats.fMaxCI - stStats.fMedian);

            /* draw box */
            cpgbox("BCNST", 0.0, 0, "BCNST", 0.0, 0);

            /* reclaim memory */
            free(g_pfLMean);
            g_pfLMean = NULL;
            free(g_pfLPostMarg_Mean);
            g_pfLPostMarg_Mean = NULL;
            free(g_pfLPostMarg_Mean);
            g_pfLPostMarg_Mean = NULL;
        }

        /* plot marginalised posterior for SD */
        {
            /* find the area for normalisation (to calculate the mean, etc.) */
            float fArea = 0.0;
            for (i = 0; i < iLenSD; ++i)
            {
                fArea += g_pfLPostMarg_SD[i];
            }
            fArea *= fStepSD;
            for (i = 0; i < iLenSD; ++i)
            {
                g_pfLPostMarg_SD[i] /= fArea;
            }

            /* find min. and max. for plotting */
            fDataMin = FLT_MAX;
            fDataMax = -(FLT_MAX);
            for (i = 0; i < iLenSD; ++i)
            {
                if (g_pfLPostMarg_SD[i] > fDataMax)
                {
                    fDataMax = g_pfLPostMarg_SD[i];
                }
                if (g_pfLPostMarg_SD[i] < fDataMin)
                {
                    fDataMin = g_pfLPostMarg_SD[i];
                }
            }
            /* keep 10% space between the top of the box and the plot */
            fDataMax = (float) (fDataMax + ((fDataMax - fDataMin) * 0.1));

            cpgpanl(1, 2);
            cpgsvp(PG_VP_ML, PG_VP_MR, PG_VP_MB, PG_VP_MT);
            cpgswin(g_pfSD[0], g_pfSD[iLenSD-1], fDataMin, fDataMax);
            cpglab("\\gs",
                   "",
                   "");
            cpgsci(PG_CI_DEF);
            cpgline(iLenSD, g_pfSD, g_pfLPostMarg_SD);
            cpgmtxt("T", -1.5, 0.85, 0.0, "(d)");

            STATS stStats = {0};
            stStats = GetStats(g_pfSD, g_pfLPostMarg_SD, fStepSD, iLenSD,
                               iPlotAll);
            if (stStats.iErrCode != EXIT_SUCCESS)
            {
                (void) fprintf(stderr, "ERROR: Getting statistics failed!\n");
                return EXIT_FAILURE;
            }
            (void) printf("Mean of sigma = %g\n", stStats.fMean);
            (void) printf("Mode of sigma = %g\n", stStats.fMode);
            (void) printf("Median of sigma = %g, P = %g\n",
                          stStats.fMedian,
                          stStats.fMedianP);
            (void) printf("CI Min = %g\n", stStats.fMinCI);
            (void) printf("CI Max = %g\n", stStats.fMaxCI);
            (void) printf("CI Min - Median = %g\n",
                          stStats.fMinCI - stStats.fMedian);
            (void) printf("CI Max - Median = %g\n",
                          stStats.fMaxCI - stStats.fMedian);

            /* draw box */
            cpgbox("BCNST", 0.0, 0, "BCNST", 0.0, 0);

            /* reclaim memory */
            free(g_pfSD);
            g_pfSD = NULL;
            free(g_pfLPostMarg_SD);
            g_pfLPostMarg_SD = NULL;
            free(g_pfLPostMarg_SD);
            g_pfLPostMarg_SD = NULL;
        }

        /* plot marginalised posterior for SMin */
        {
            /* find the area for normalisation (to calculate the mean, etc.) */
            float fArea = 0.0;
            for (i = 0; i < iLenSMin; ++i)
            {
                fArea += g_pfLPostMarg_SMin[i];
            }
            fArea *= fStepSMin;
            for (i = 0; i < iLenSMin; ++i)
            {
                g_pfLPostMarg_SMin[i] /= fArea;
            }

            /* find min. and max. for plotting */
            fDataMin = FLT_MAX;
            fDataMax = -(FLT_MAX);
            for (i = 0; i < iLenSMin; ++i)
            {
                if (g_pfLPostMarg_SMin[i] > fDataMax)
                {
                    fDataMax = g_pfLPostMarg_SMin[i];
                }
                if (g_pfLPostMarg_SMin[i] < fDataMin)
                {
                    fDataMin = g_pfLPostMarg_SMin[i];
                }
            }
            /* keep 10% space between the top of the box and the plot */
            fDataMax = (float) (fDataMax + ((fDataMax - fDataMin) * 0.1));

            cpgpanl(2, 2);
            cpgsvp(PG_VP_ML, PG_VP_MR, PG_VP_MB, PG_VP_MT);
            cpgswin(g_pfSMin[0], g_pfSMin[iLenSMin-1], fDataMin, fDataMax);
            cpglab("S\\dmin\\u (mJy)",
                   "",
                   "");
            cpgsci(PG_CI_DEF);
            cpgline(iLenSMin, g_pfSMin, g_pfLPostMarg_SMin);
            cpgmtxt("T", -1.5, 0.85, 0.0, "(e)");

            STATS stStats = {0};
            stStats = GetStats(g_pfSMin, g_pfLPostMarg_SMin, fStepSMin,
                               iLenSMin, iPlotAll);
            if (stStats.iErrCode != EXIT_SUCCESS)
            {
                (void) fprintf(stderr, "ERROR: Getting statistics failed!\n");
                return EXIT_FAILURE;
            }
            (void) printf("Mean of Smin = %g\n", stStats.fMean);
            (void) printf("Mode of Smin = %g\n", stStats.fMode);
            (void) printf("Median of Smin = %g, P = %g\n",
                          stStats.fMedian,
                          stStats.fMedianP);
            (void) printf("CI Min = %g\n", stStats.fMinCI);
            (void) printf("CI Max = %g\n", stStats.fMaxCI);
            (void) printf("CI Min - Median = %g\n",
                          stStats.fMinCI - stStats.fMedian);
            (void) printf("CI Max - Median = %g\n",
                          stStats.fMaxCI - stStats.fMedian);

            /* draw box */
            cpgbox("BCNST", 0.0, 0, "BCNST", 0.0, 0);

            /* reclaim memory */
            free(g_pfSMin);
            g_pfSMin = NULL;
            free(g_pfLPostMarg_SMin);
            g_pfLPostMarg_SMin = NULL;
            free(g_pfLPostMarg_SMin);
            g_pfLPostMarg_SMin = NULL;
        }

        /* plot marginalised posterior for r */
        {
            /* find the area for normalisation (to calculate the mean, etc.) */
            float fArea = 0.0;
            for (i = 0; i < iLenDist; ++i)
            {
                fArea += g_pfLPostMarg_d[i];
            }
            fArea *= fStepDist;
            for (i = 0; i < iLenDist; ++i)
            {
                g_pfLPostMarg_d[i] /= fArea;
            }

            /* find min. and max. for plotting */
            fDataMin = FLT_MAX;
            fDataMax = -(FLT_MAX);
            for (i = 0; i < iLenDist; ++i)
            {
                if (g_pfLPostMarg_d[i] > fDataMax)
                {
                    fDataMax = g_pfLPostMarg_d[i];
                }
                if (g_pfLPostMarg_d[i] < fDataMin)
                {
                    fDataMin = g_pfLPostMarg_d[i];
                }
            }
            /* keep 10% space between the top of the box and the plot */
            fDataMax = (float) (fDataMax + ((fDataMax - fDataMin) * 0.1));

            cpgpanl(3, 2);
            cpgsvp(PG_VP_ML, PG_VP_MR, PG_VP_MB, PG_VP_MT);
            cpgswin(g_pfDist[0], g_pfDist[iLenDist-1], fDataMin, fDataMax);
            cpglab("r (kpc)",
                   "",
                   "");
            cpgsci(PG_CI_DEF);
            cpgline(iLenDist, g_pfDist, g_pfLPostMarg_d);
            cpgmtxt("T", -1.5, 0.85, 0.0, "(f)");

            STATS stStats = {0};
            stStats = GetStats(g_pfDist, g_pfLPostMarg_d, fStepDist, iLenDist,
                               iPlotAll);
            if (stStats.iErrCode != EXIT_SUCCESS)
            {
                (void) fprintf(stderr, "ERROR: Getting statistics failed!\n");
                return EXIT_FAILURE;
            }
            (void) printf("Mean of d = %g\n", stStats.fMean);
            (void) printf("Mode of d = %g\n", stStats.fMode);
            (void) printf("Median of d = %g, P = %g\n",
                          stStats.fMedian,
                          stStats.fMedianP);
            (void) printf("CI Min = %g\n", stStats.fMinCI);
            (void) printf("CI Max = %g\n", stStats.fMaxCI);
            (void) printf("CI Min - Median = %g\n",
                          stStats.fMinCI - stStats.fMedian);
            (void) printf("CI Max - Median = %g\n",
                          stStats.fMaxCI - stStats.fMedian);

            /* draw box */
            cpgbox("BCNST", 0.0, 0, "BCNST", 0.0, 0);

            /* reclaim memory */
            free(g_pfDist);
            g_pfDist = NULL;
            free(g_pfLPostMarg_d);
            g_pfLPostMarg_d = NULL;
            free(g_pfLPostMarg_d);
            g_pfLPostMarg_d = NULL;
        }
    }

    /* close PGPLOT device */
    cpgclos();

    return EXIT_SUCCESS;
}

/* calculate p_obs */
float CalcPObs(float fSMean, float fSD, float fSMin)
{
    return (erfc((log10(fSMin) - fSMean) / (fSD * SQRT_2)) / 2);
}


/* calculate step 1 likelihood */
float CalcS1Lhd(int in, float *pfFlux,
                float fp_obs, float fSMean, float fSD)
{
    float fOpt = 2 * fSD * fSD;
    float fS1Lhd = 1.0;

    for (int i = 0; i < in; ++i)
    {
        fS1Lhd *= (exp(-pow((log10(pfFlux[i]) - fSMean), 2) / fOpt)
                   / (fp_obs * fSD * sqrt(2 * M_PI)));
    }

    return fS1Lhd;
}

/* from Numerical Recipes */
float gammln(float xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}
float factln(int n)
{
    static float a[101];

    if (n < 0) {
      puts("Negative factorial in routine factln");
      exit(0);
    }
    if (n <= 1) return 0.0;
    if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
    else return gammln(n+1.0);
}
/* calculate step 2 likelihood */
float CalcS2Lhd(int in, float fp_obs, int iN)
{
    return (exp((float) (factln(iN) - factln(in) - factln(iN - in)))
            * pow(fp_obs, in)
            /* NOTE: pow(x, 0) is 1 in the math library, even when x is 0 or
               NaN! */
            * pow((1.0 - fp_obs), (iN - in)));
}

/* calculate step 3 likelihood */
float CalcS3Lhd(float fSObs, int iN, float fSMean, float fSD)
{
    float fSExp = pow(10, fSMean + ((log(10) * fSD * fSD) / 2));
    float fSStd = fSExp * sqrt(pow(10, log(10) * fSD * fSD) - 1);

    /* expectation of total flux */
    float fS_diff = iN * fSExp;
    /* standard deviation of total flux */
    float fsigma_diff = sqrt(iN) * fSStd;

    /* compute and return the step 3 likelihood */
    return (exp(-pow((fSObs - fS_diff), 2) / (2 * fsigma_diff * fsigma_diff))
            / (sqrt(2 * M_PI) * fsigma_diff));
}

int GetIndex(float* pfArray, int iLen, float fEpsilon, float fVal)
{
    int iIdx = -1;  /* to catch misses */
    int i = 0;

    for (i = 0; i < iLen; ++i)
    {
        if (fabsf(pfArray[i] - fVal) < fEpsilon)
        {
            iIdx = i;
            break;
        }
    }

    if (iIdx < 0)
    {
        return -1; 
    }

    return iIdx;
}

/*
 * Compute the mean, median, mode, and confidence intervals of a given
 * probability distribution
 */
STATS GetStats(float* pfX, float* pfP, float fStepX, int iLen, int iPlotAll)
{
    int iIdxCIMax = 0;
    int iIdxCIMin = 0;
    float* pfCDF = NULL;
    float fConfInt = DEF_CONF_INT;
    float fDiff = 0.0;
    float fDiffOld = FLT_MAX;
    float* pfCIPoly = NULL;
    float* pfCIPolyInt = NULL;
    int iIdxDiffMin = 0;
    int i = 0;
    int j = 0;
    STATS stStats = {0};
    float fMaxP = -(FLT_MAX);
    float fMinCI = 0.0;
    float fMaxCI = 0.0;

    stStats.iErrCode = EXIT_SUCCESS;

    /* find the mean and mode */
    for (i = 0; i < iLen; ++i)
    {
        stStats.fMean += (pfX[i] * pfP[i]);
        if (pfP[i] > fMaxP)
        {
            fMaxP = pfP[i];
            stStats.fMode = pfX[i];
        }
    }
    stStats.fMean *= fStepX;

    pfCDF = (float*) malloc((size_t) iLen * sizeof(float));
    if (NULL == pfCDF)
    {
        (void) fprintf(stderr,
                       "ERROR: Memory allocation failed in %s at %d! "
                       "%s\n",
                       __FILE__,
                       __LINE__,
                       strerror(errno));
        stStats.iErrCode = EXIT_FAILURE;
        return stStats;
    }

    /* calculate the CDF and apply the confidence interval */
    for (i = 0; i < iLen; ++i)
    {
        pfCDF[i] = 0.0;
        for (j = 0; j <= i; ++j)
        {
            pfCDF[i] += pfP[j];
        }
        pfCDF[i] *= fStepX;
    }

    /* calculate the confidence limits from the confidence interval */
    fMinCI = (100.0 - fConfInt) / 2;
    fMaxCI = (100.0 - fMinCI);
    fMinCI /= 100;
    fMaxCI /= 100;

    /* get the index for which pfCDF[i] = fMinCI */
    for (i = 0; i < iLen; ++i)
    {
        if (0 == pfCDF[i])
        {
            continue;
        }
        fDiff = fabs(pfCDF[i] - fMinCI);
        if (fDiff < fDiffOld)
        {
            iIdxCIMin = i;
        }
        fDiffOld = fDiff;
    }
    stStats.fMinCI = pfX[iIdxCIMin];
    fDiffOld = FLT_MAX;
    /* get the index for which pfCDF[i] = fMaxCI */
    for (i = 0; i < iLen; ++i)
    {
        fDiff = fabs(pfCDF[i] - fMaxCI);
        if (fDiff < fDiffOld)
        {
            iIdxCIMax = i;
        }
        fDiffOld = fDiff;
    }
    stStats.fMaxCI = pfX[iIdxCIMax];

    if (iIdxCIMin != 0) /* no need to draw CI otherwise */
    {
        /* allocate memory for the shaded polygon arrays */
        pfCIPoly = (float *) malloc((size_t) (iIdxCIMin + 2) * sizeof(float));
        if (NULL == pfCIPoly)
        {
            (void) fprintf(stderr,
                           "ERROR: Memory allocation failed in %s at %d! "
                           "%s\n",
                           __FILE__,
                           __LINE__,
                           strerror(errno));
            free(pfCDF);
            stStats.iErrCode = EXIT_FAILURE;
            return stStats;
        }
        pfCIPolyInt = (float *) malloc((size_t) (iIdxCIMin + 2) * sizeof(float));
        if (NULL == pfCIPolyInt)
        {
            (void) fprintf(stderr,
                           "ERROR: Memory allocation failed in %s at %d! "
                           "%s\n",
                           __FILE__,
                           __LINE__,
                           strerror(errno));
            free(pfCDF);
            free(pfCIPoly);
            stStats.iErrCode = EXIT_FAILURE;
            return stStats;
        }

        for (i = 0; i < iIdxCIMin; ++i)
        {
            pfCIPoly[i] = pfX[i];
            pfCIPolyInt[i] = pfP[i];
        }
        pfCIPoly[i] = pfCIPoly[i-1];
        pfCIPoly[i+1] = pfCIPoly[0];
        pfCIPolyInt[i] = 0;
        pfCIPolyInt[i+1] = 0;
        cpgsfs(3);
        cpgsci(PG_CI_DEF);
        cpgpoly(i+2, pfCIPoly, pfCIPolyInt);

        free(pfCIPoly);
        pfCIPoly = NULL;
        free(pfCIPolyInt);
        pfCIPolyInt = NULL;
    }

    if (iIdxCIMax != iLen)  /* no need to draw CI otherwise */
    {
        /* reallocate memory for the polygon arrays */
        /* backup the original pointer before realloc() */
        pfCIPoly = (float *) malloc((size_t) (iLen - iIdxCIMax + 2) * sizeof(float));
        if (NULL == pfCIPoly)
        {
            (void) fprintf(stderr,
                           "ERROR: Memory allocation failed! %s\n",
                           strerror(errno));
            free(pfCDF);
            stStats.iErrCode = EXIT_FAILURE;
            return stStats;
        }
        /* backup the original pointer before realloc() */
        pfCIPolyInt = (float *) malloc((size_t) (iLen - iIdxCIMax + 2) * sizeof(float));
        if (NULL == pfCIPolyInt)
        {
            (void) fprintf(stderr,
                           "ERROR: Memory allocation failed! %s\n",
                           strerror(errno));
            free(pfCDF);
            free(pfCIPoly);
            stStats.iErrCode = EXIT_FAILURE;
            return stStats;
        }

        for (i = 0, j = iIdxCIMax + 1; j < iLen; ++i, ++j)
        {
            pfCIPoly[i] = pfX[j];
            pfCIPolyInt[i] = pfP[j];
        }
        /* TODO: fix the problem when iIdxCIMax +1 = iLen */
        pfCIPoly[i] = pfCIPoly[i-1];
        pfCIPoly[i+1] = pfCIPoly[0];
        pfCIPolyInt[i] = 0;
        pfCIPolyInt[i+1] = 0;
        cpgsfs(3);
        cpgpoly(i+2, pfCIPoly, pfCIPolyInt);
    }

    /* find the median */
    /* get the value of N for which pfCDF[i] = 0.5 */
    fDiffOld = FLT_MAX;
    for (i = 0; i < iLen; ++i)
    {
        fDiff = fabs(pfCDF[i] - (float) 0.5);
        if (fDiff < fDiffOld)
        {
            iIdxDiffMin = i;
        }
        fDiffOld = fDiff;
    }

    stStats.fMedian = pfX[iIdxDiffMin];
    stStats.fMedianP = pfCDF[iIdxDiffMin];

    free(pfCDF);
    free(pfCIPoly);
    free(pfCIPolyInt);

    return stStats;
}

