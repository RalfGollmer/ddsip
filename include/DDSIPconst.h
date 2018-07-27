#ifndef DDSIPCONST_H
#define DDSIPCONST_H
#ifdef _WIN32
#pragma pack(push, 8)
#endif

#ifdef __cplusplus
extern "C" {
#endif

// Identifier should be unique within this number of letters
    extern const int    DDSIP_unique ;

// Filenames are not longer than
#define    DDSIP_ln_fname 1024

// Variable names in model file are not longer than
#define    DDSIP_ln_varname 65

// longest string occuring in file
#define    DDSIP_max_str_ln 512

// OUTLEV for first-stage solutions in LB
#define    DDSIP_first_stage_outlev 25

// OUTLEV for second-stage solutions in LB + UB
#define    DDSIP_second_stage_outlev 49

// OUTLEV for suggested first stage in UB
#define    DDSIP_suggested_first_stage_outlev 23

// OUTLEV for current lambda
#define    DDSIP_current_lambda_outlev 26

// Indicates 'Control C'
    extern int     DDSIP_killsignal;

// Maximal number of parameters
    extern const int DDSIP_maxparam;

// Number of implemented risk measures
    extern const int DDSIP_maxrisk;

// Large values
    extern const int    DDSIP_bigint ;
    extern const double DDSIP_bigvalue ;
    extern const double DDSIP_infty ;

// Small values
    extern const double DDSIP_brancheps ;
    extern const double DDSIP_epsilon ;
    extern const double DDSIP_accuracy ;

// Output directory
    extern const char   DDSIP_outdir[8];
    extern const char   DDSIP_outfname[16];

// Statistik of several runs
    extern const char   DDSIP_solfname[32];

// Debugging output
    extern const char   DDSIP_moreoutfname[16];

// Recourse function to file
    extern const char   DDSIP_recfunfname[32];

#ifdef __cplusplus
} // extern "C"
#endif

#ifdef _WIN32
#pragma pack(pop)
#endif

#endif
