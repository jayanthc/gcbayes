/**
 * @file gcbayes_func.c
 *
 * @author Jayanth Chennamangalam
 * @date 2010.05.26
 */

#include "gcbayes.h"

extern float* g_pfSMean;
extern float* g_pfLMean;
extern float* g_pfSD;
extern float* g_pfN;
extern float* g_pfLogN;
extern float* g_pfSMin;
extern float* g_pfS_i;
extern float* g_pfDist;
extern float* g_pfPDist;
extern float* g_pfSTotLhd;
extern double* g_pdSTotLhd;
extern float* g_pfLTotLhd;
extern float* g_pfLPost;
extern float* g_pfLPostMarg_MeanSD;
extern float* g_pfLPostMarg_N;
extern float* g_pfLPostMarg_Mean;
extern float* g_pfLPostMarg_SD;
extern float* g_pfLPostMarg_SMin;
extern float* g_pfLPostMarg_d;

double (*g_padErfTable)[][2] = NULL;
int g_iErfEntries = 0;

/*
 * Registers handlers for SIGTERM and CTRL+C
 */
int RegisterSignalHandlers(void)
{
    struct sigaction stSigHandler = {{0}};
    int iRet = EXIT_SUCCESS;

    /* register the CTRL+C-handling function */
    stSigHandler.sa_handler = HandleStopSignals;
    iRet = sigaction(SIGINT, &stSigHandler, NULL);
    if (iRet != EXIT_SUCCESS)
    {
        (void) fprintf(stderr,
                       "ERROR: Signal handler registration failed for "
                       "signal %d!\n",
                       SIGINT);
        return EXIT_FAILURE;
    }

    /* register the SIGTERM-handling function */
    stSigHandler.sa_handler = HandleStopSignals;
    iRet = sigaction(SIGTERM, &stSigHandler, NULL);
    if (iRet != EXIT_SUCCESS)
    {
        (void) fprintf(stderr,
                       "ERROR: Signal handler registration failed for "
                       "signal %d!\n",
                       SIGTERM);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

/*
 * Catches SIGTERM and CTRL+C and cleans up before exiting
 */
void HandleStopSignals(int iSigNo)
{
    /* clean up */
    CleanUp();

    /* exit */
    exit(EXIT_SUCCESS);

    /* never reached */
    return;
}

/*
 * Computes the amount of free memory
 */
unsigned long int GetFreeMem(void)
{
    return (sysconf(_SC_PAGESIZE) * sysconf(_SC_AVPHYS_PAGES));
}

/*
 * Reads the per-GC configuration file
 */
GC_CONF ReadGCConf(char* pcFileConf)
{
    GC_CONF stConf = {{0}};
    FILE* pFConf = NULL;
    char* pcLine = NULL;
    size_t iLen = 0;

    pFConf = fopen(pcFileConf, "r");
    if (NULL == pFConf)
    {
        (void) fprintf(stderr,
                       "ERROR: Opening config file failed! %s\n",
                       strerror(errno));
        return stConf;
    }

    (void) getline(&pcLine, &iLen, pFConf); 
    (void) sscanf(pcLine, "%s", stConf.acGCName);

    (void) getline(&pcLine, &iLen, pFConf); 
    (void) sscanf(pcLine, "%s", stConf.acFileFlux);

    (void) getline(&pcLine, &iLen, pFConf); 
    (void) sscanf(pcLine, "%d", &stConf.in);

    (void) getline(&pcLine, &iLen, pFConf); 
    (void) sscanf(pcLine, "%f", &stConf.fDist);

    (void) getline(&pcLine, &iLen, pFConf); 
    (void) sscanf(pcLine, "%f", &stConf.fDistSD);

    (void) getline(&pcLine, &iLen, pFConf); 
    (void) sscanf(pcLine, "%f", &stConf.fSObs);

    (void) getline(&pcLine, &iLen, pFConf); 
    (void) sscanf(pcLine, "%f", &stConf.fMinSMean);

    (void) getline(&pcLine, &iLen, pFConf); 
    (void) sscanf(pcLine, "%f", &stConf.fMaxSMean);

    free(pcLine);

    /* close the flux file */
    (void) fclose(pFConf);

    /* open the PGPLOT graphics device */
    return stConf;
}

/*
 * Cleans up all allocated memory
 */
void CleanUp()
{
    if (g_pfN != NULL)
    {
        free(g_pfN);
        g_pfN = NULL;
    }
    if (g_pfLogN != NULL)
    {
        free(g_pfLogN);
        g_pfLogN = NULL;
    }
    if (g_pfSMean != NULL)
    {
        free(g_pfSMean);
        g_pfSMean = NULL;
    }
    if (g_pfSD != NULL)
    {
        free(g_pfSD);
        g_pfSD = NULL;
    }
    if (g_pfSMin != NULL)
    {
        free(g_pfSMin);
        g_pfSMin = NULL;
    }
    if (g_pfS_i != NULL)
    {
        free(g_pfS_i);
        g_pfS_i = NULL;
    }
    if (g_pfDist != NULL)
    {
        free(g_pfDist);
        g_pfDist = NULL;
    }
    if (g_pfPDist != NULL)
    {
        free(g_pfPDist);
        g_pfPDist = NULL;
    }
    if (g_pfLMean != NULL)
    {
        free(g_pfLMean);
        g_pfLMean = NULL;
    }
    if (g_pfSTotLhd != NULL)
    {
        free(g_pfSTotLhd);
        g_pfSTotLhd = NULL;
    }
    if (g_pdSTotLhd != NULL)
    {
        free(g_pdSTotLhd);
        g_pdSTotLhd = NULL;
    }
    if (g_pfLTotLhd != NULL)
    {
        free(g_pfLTotLhd);
        g_pfLTotLhd = NULL;
    }
    if (g_pfLPost != NULL)
    {
        free(g_pfLPost);
        g_pfLPost = NULL;
    }
    if (g_pfLPostMarg_MeanSD != NULL)
    {
        free(g_pfLPostMarg_MeanSD);
        g_pfLPostMarg_MeanSD = NULL;
    }
    if (g_pfLPostMarg_N != NULL)
    {
        free(g_pfLPostMarg_N);
        g_pfLPostMarg_N = NULL;
    }
    if (g_pfLPostMarg_Mean != NULL)
    {
        free(g_pfLPostMarg_Mean);
        g_pfLPostMarg_Mean = NULL;
    }
    if (g_pfLPostMarg_SD != NULL)
    {
        free(g_pfLPostMarg_SD);
        g_pfLPostMarg_SD = NULL;
    }
    if (g_pfLPostMarg_SMin != NULL)
    {
        free(g_pfLPostMarg_SMin);
        g_pfLPostMarg_SMin = NULL;
    }
    if (g_pfLPostMarg_d != NULL)
    {
        free(g_pfLPostMarg_d);
        g_pfLPostMarg_d = NULL;
    }

    return;
}

/*
 * Prints usage information
 */
void PrintUsage(const char *pcProgName)
{
    (void) printf("Usage: %s [options]\n", pcProgName);
    (void) printf("    -h  --help                           ");
    (void) printf("Display this usage information\n");
    (void) printf("    -N  --n-min <value>                  ");
    (void) printf("Minimum N\n");
    (void) printf("                                         ");
    (void) printf("(default is %d)\n", DEF_MIN_N);
    (void) printf("    -S  --n-step <value>                 ");
    (void) printf("N step size\n");
    (void) printf("                                         ");
    (void) printf("(default is %d)\n", DEF_STEP_N);
    (void) printf("    -X  --n-max <value>                  ");
    (void) printf("Maximum N\n");
    (void) printf("                                         ");
    (void) printf("(default is %d)\n", DEF_MAX_N);
    (void) printf("    -m  --mean-min <value>               ");
    (void) printf("The minimum mean (<log S>)\n");
    (void) printf("                                         ");
    (void) printf("(default is %g)\n", DEF_MIN_LMEAN);
    (void) printf("    -p  --mean-step <value>              ");
    (void) printf("The mean (<log L>) step size\n");
    (void) printf("                                         ");
    (void) printf("(default is %g)\n", DEF_STEP_MEAN);
    (void) printf("    -u  --mean-max <value>               ");
    (void) printf("The maximum mean (<log L>)\n");
    (void) printf("                                         ");
    (void) printf("(default is %g)\n", DEF_MAX_LMEAN);
    (void) printf("    -s  --sd-min <value>                 ");
    (void) printf("The minimum standard deviation (of\n");
    (void) printf("                                         ");
    (void) printf("log L)\n");
    (void) printf("                                         ");
    (void) printf("(default is %g)\n", DEF_MIN_SD);
    (void) printf("    -e  --sd-step <value>                ");
    (void) printf("The standard deviation (of log L)\n");
    (void) printf("                                         ");
    (void) printf("step size\n");
    (void) printf("                                         ");
    (void) printf("(default is %g)\n", DEF_STEP_SD);
    (void) printf("    -a  --sd-max <value>                 ");
    (void) printf("The maximum standard deviation (of\n");
    (void) printf("                                         ");
    (void) printf("log L)\n");
    (void) printf("                                         ");
    (void) printf("(default is %g)\n", DEF_MAX_SD);
    (void) printf("    -l  --s-min-min <value>              ");
    (void) printf("The minimum flux in mJy\n");
    (void) printf("                                         ");
    (void) printf("(default is %g)\n", DEF_MIN_S_MIN);
    //(void) printf("    -t  --s-min-step <value>             ");
    //(void) printf("The flux step size in mJy\n");
    //(void) printf("                                         ");
    //(void) printf("(default is %g)\n", DEF_STEP_S_MIN);
    (void) printf("    -y  --conf-file <filename>           ");
    (void) printf("Config file\n");
    (void) printf("    -c  --colour-map <name>              ");
    (void) printf("Name of colour map\n");
    (void) printf("    -f  --plot-ps                        ");
    (void) printf("Print to PS file, instead of screen\n");
    (void) printf("    -z  --plot-all                       ");
    (void) printf("Plot both S and L posteriors\n");
    (void) printf("    -w  --plot-lall                      ");
    (void) printf("Plot all marginalised posteriors for L\n");

    return;
}

