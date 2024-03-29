/*  Authors:           Andreas M"arkert, Ralf Gollmer
    Copyright to:      University of Duisburg-Essen
    Language:          C

    Description:
    This file contains the main procedure of DDSIP --
    a decomposition method for two-stage stochastic mixed-integer programs
    based on the Lagrangian relaxation of nonanticipativity and on
    branch-and-bound.

    License:
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

#include <DDSIP.h>
#include <DDSIPconst.h>
#include <DDSIPversion.h>
#include <sys/resource.h>

#define CBFORALL 0
#define DELNODE1MULT

CPXENVptr    DDSIP_env = NULL;
CPXLPptr     DDSIP_lp  = NULL;

FILE      * DDSIP_outfile = NULL;

// These pointers have to be global if we want to use ConicBundle
para_t    * DDSIP_param= NULL;
data_t    * DDSIP_data = NULL;
bb_t      * DDSIP_bb   = NULL;
node_t   ** DDSIP_node = NULL;
// Identifiers are DDSIP_unique within ... letters
const int DDSIP_unique = 6;

// Maximal number of parameters (in each CPLEX section)
const int DDSIP_maxparam = 64;

// Number of risk measures implemented
const int DDSIP_maxrisk = 7;

// Large values
const int DDSIP_bigint      = 10000000;       // Upper bound on integer parameters
const double DDSIP_bigvalue = 1.0e9;       // Just to detect the print format
const double DDSIP_infty    = CPX_INFBOUND; // is 1.0e20; -- Infinity

// Output directory
const char DDSIP_outdir[8] = "sipout";

// Output file
const char DDSIP_outfname[16] = "sipout/sip.out";

// More output is printed if desired
const char DDSIP_moreoutfname[16] = "sipout/more.out";

// Solution output file
const char DDSIP_solfname[32] = "sipout/solution.out";

// Recourse function to file
//const char DDSIP_recfunfname[32] = "sipout/recfun.out";

// Indicates user termination ('Control-C')
int DDSIP_killsignal = 0;


//==========================================================================
int
main (int argc, char * argv[])
{
    // Temporary variables
    struct stat filestat;

    char astring[DDSIP_ln_fname];

    int status = 0, cont, boundstat, comb, i;

    // guarantee minimal stacksize limit
    const rlim_t kStackSize = 128 * 1024 * 1024;   // min stack size = 128 MB
    struct rlimit rl;
    int result;
    unsigned long cur_limit;
    double tmp_objval;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {
        cur_limit = rl.rlim_cur;
        printf ("original stacksize limit %lu\n", cur_limit);
        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
            result = getrlimit(RLIMIT_STACK, &rl);
            if (result == 0)
            {
                cur_limit = rl.rlim_cur;
                printf ("changed  stacksize limit %lu\n", cur_limit);
            }
        }
    }

    i = argc;
    //

    // Welcome
    printf ("#######################################################################################\n");
    printf ("#########  D D S I P -- Dual Decomposition In Stochastic Integer Programming  #########\n");
    printf ("#########  %-66s #########\n", DDSIP_version);
    printf ("For copyright, license agreements, help, comments, requests, ... ");
    printf ("see\n\thttp://www.uni-duisburg-essen.de/~hn215go/ddsip.shtml\n");
    printf ("\thttp://www.github.com/RalfGollmer/ddsip\n");
    /*    printf ("  Copyright (C) to University of Duisburg-Essen \n\n"); */
    /*    printf("  This program is free software; you can redistribute it and/or\n"); */
    /*    printf ("  modify it under the terms of the GNU General Public License\n"); */
    /*    printf("  as published by the Free Software Foundation; either version 2\n"); */
    /*    printf ("  of the License, or (at your option) any later version.\n\n"); */
    /*    printf ("  This program is distributed in the hope that it will be useful,\n"); */
    /*    printf("  but WITHOUT ANY WARRANTY; without even the implied warranty of\n"); */
    /*    printf ("  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"); */
    /*    printf ("  GNU General Public License for more details.\n\n"); */
    /*    printf ("  You should have received a copy of the GNU General Public License\n"); */
    /*    printf("  along with this program; if not, write to the Free Software\n"); */
    /*    printf ("  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.\n"); */
    /*    printf ("==============================================================================\n\n"); */

    // These four structures shall carry the whole problem data
    DDSIP_param = NULL;
    DDSIP_bb = NULL;
    DDSIP_data = NULL;
    DDSIP_node = NULL;

    // In DDSIP some signals are handled seperatly
    DDSIP_RegisterSignalHandlers ();

    // Create output directory if it doesn't already exist
    if (stat (DDSIP_outdir, &filestat) == -1)
    {
        sprintf (astring, "mkdir %s\n", DDSIP_outdir);
        i = system (astring);
        printf ("Creating subdirectory %s ...\n", DDSIP_outdir);
    }
    // remove possibly existing output file
    sprintf (astring, "rm -f %s\n", DDSIP_outfname);
    i = system (astring);
    // Open output file
    if ((DDSIP_outfile = fopen (DDSIP_outfname, "a")) == NULL)
    {
        fprintf (stderr, "ERROR: Cannot open '%s'. \n", DDSIP_outfname);
        status = 107;
        goto TERMINATE;
    }

    setbuf (DDSIP_outfile, 0);
    fprintf (DDSIP_outfile, "%s\n", argv[0]);
    fprintf (DDSIP_outfile, "DDSIP  %s\n", DDSIP_version);
    fprintf (DDSIP_outfile, "-----------------------------------------------------------\n");

    fprintf (DDSIP_outfile, "current system time: ");
    fflush (DDSIP_outfile);
#ifndef _WIN32
    i = system ("date");
    // Print time to output file
    sprintf (astring, "date >> %s\n", DDSIP_outfname);
    i = system (astring);
    fprintf (DDSIP_outfile, "host            : ");
    sprintf (astring, "hostname >> %s; cat /proc/cpuinfo | sed '/processor.*: 0/,/^$/!d' |grep -E 'vendor|cpu |model|stepping|MHz|cache |cores' >> %s\n", DDSIP_outfname, DDSIP_outfname);
    for (i = 0; i < 100; i++)
        exp(1.11*(-i-2));
    i = system (astring);
#else
    sprintf (astring, "date /T >> %s & time /T >> %s\n", DDSIP_outfname,DDSIP_outfname);
    i = system (astring);
#endif
    fprintf (DDSIP_outfile, "-----------------------------------------------------------\n");

    // Open cplex environment
    DDSIP_env = CPXopenCPLEX (&status);
    if (DDSIP_env == NULL)
    {
        fprintf (stderr, "ERROR: Failed to open cplex environment, CPLEX error code %d.\n",status);
        fprintf (DDSIP_outfile, "ERROR: Failed to open cplex environment, CPLEX error code %d.\n",status);
        return status;
    }
    sprintf (astring, "%s", CPXversion (DDSIP_env));
    printf ("CPLEX version is %s\n", astring);
    fprintf (DDSIP_outfile, "CPLEX version is %s\n", astring);
    fprintf (DDSIP_outfile, "CB    version is %s\n", CBV);

    // Allocate the structures to hold the information on the problem
    DDSIP_param = (para_t *) DDSIP_Alloc (sizeof (para_t), 1, "param(Main)");
    DDSIP_bb    = (bb_t *)   DDSIP_Alloc (sizeof (bb_t),   1, "bb(Main)");
    DDSIP_data  = (data_t *) DDSIP_Alloc (sizeof (data_t), 1, "data(Main)");

    // Read model (lp) file
    if ((status = DDSIP_ReadModel ()))
        goto TERMINATE;

    // Read specification file
    if ((status = DDSIP_ReadSpec ()))
        goto TERMINATE;

    // Set specified CPLEX parameters
    if ((status = DDSIP_InitCpxPara ()))
        goto TERMINATE;

    // if a Benders annotation file is given, read it
    if (DDSIP_param->annotationFile)
    {
        status = CPXreadcopyannotations (DDSIP_env, DDSIP_lp, DDSIP_param->annotationFile);
        if ( status )
        {
            fprintf (stderr, "ERROR: Failed to read and copy the annotation data from file %s.\n", DDSIP_param->annotationFile);
            fprintf (DDSIP_outfile, "ERROR: Failed to read and copy the annotation data from file %s.\n", DDSIP_param->annotationFile);
            goto TERMINATE;
        }
    }

    // data->cost contains: stoch. costs for all scenarios, followed by the original costs
    DDSIP_data->cost =
        (double *) DDSIP_Alloc (sizeof (double),
                                DDSIP_param->scenarios * DDSIP_param->stoccost + DDSIP_data->novar, "cost (ReadData)");
    // Store original objective coefficients
    status = CPXgetobj (DDSIP_env, DDSIP_lp, DDSIP_data->cost + DDSIP_param->scenarios * DDSIP_param->stoccost, 0, DDSIP_data->novar - 1);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to get objective coefficients\n");
        fprintf (DDSIP_outfile, "ERROR: Failed to get objective coefficients\n");
        goto TERMINATE;
    }

    DDSIP_bb->moreoutfile = NULL;
    if (DDSIP_param->outlev)
    {
        // Open debug output file
        if ((DDSIP_bb->moreoutfile = fopen (DDSIP_moreoutfname, "a")) == NULL)
        {
            fprintf (stderr, "ERROR: Cannot open '%s'. \n", DDSIP_moreoutfname);
            fprintf (DDSIP_outfile, "ERROR: Cannot open '%s'. \n", DDSIP_moreoutfname);
            status = 109;
            goto TERMINATE;
        }
        if (DDSIP_param->outlev > 10)
        {
            // Buffer size = 0
            setbuf (stdout, 0);
            if (DDSIP_param->outlev > 20)
                setbuf (DDSIP_bb->moreoutfile, 0);
        }
#ifndef _WIN32
        // Print time to output file
        sprintf (astring, "date > %s\n", DDSIP_moreoutfname);
        i = system (astring);
#else
        sprintf (astring, "date /T > %s & time /T >> %s\n", DDSIP_moreoutfname,DDSIP_moreoutfname);
        i = system (astring);
#endif
        fprintf (DDSIP_bb->moreoutfile, "--------------------------------------------------------------");
        fprintf (DDSIP_bb->moreoutfile, "---------\nThis is an additional output file of DDSIP. The ");
        fprintf (DDSIP_bb->moreoutfile, "actual amount of output\nis controlled by the parameter ");
        fprintf (DDSIP_bb->moreoutfile, "OUTLEV in the specification file.\n---------");
        fprintf (DDSIP_bb->moreoutfile, "--------------------------------------------------------------\n\n");

        if (DDSIP_param->cb)
            fprintf (DDSIP_bb->moreoutfile, "front node re-eval at nodes %d and %d.\n",
                  2*((DDSIP_param->cbContinuous + DDSIP_param->cbBreakIters)/2) + DDSIP_Imax (2, 2*(-DDSIP_param->cb/2-1)),
                  2*(DDSIP_param->cbContinuous + 2*(DDSIP_param->cbBreakIters/2) + 4)
                );
    }

    DDSIP_node = (node_t **) DDSIP_Alloc (sizeof (node_t *), DDSIP_param->nodelim + 3, "node(Main)");
    DDSIP_node[0] = (node_t *) DDSIP_Alloc (sizeof (node_t), 1, "node[0](Main)");
    DDSIP_node[0]->first_sol = (double **) DDSIP_Alloc (sizeof (double *), DDSIP_param->scenarios, "node[0]->first_sol(Main)");
    DDSIP_node[0]->cursubsol = (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "node[0]->cursubsol(Main)");
    DDSIP_node[0]->subbound = (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "node[0]->subbound(Main)");
    if (DDSIP_param->cb)
    {
        DDSIP_node[0]->scenBoundsNoLag = (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "node[0]->scenBoundsNoLag(Main)");
        DDSIP_node[0]->BoundNoLag = -DDSIP_infty;

    }
    DDSIP_node[0]->mipstatus = (int *) DDSIP_Alloc (sizeof (int), DDSIP_param->scenarios, "node[0]->mipstatus(Main)");
    DDSIP_node[0]->ref_scenobj = (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "node[0]->ref_scenobj(Main)");

    // Prepare first and second stage variables
    if ((status = DDSIP_InitStages ()))
    {
        fprintf (stderr, "ERROR: Failed to initialize stages\n");
        fprintf (DDSIP_outfile, "ERROR: Failed to initialize stages\n");
        goto TERMINATE;
    }

    // Read data file(s)
    if ((status = DDSIP_ReadData ()))
        goto TERMINATE;

    // Read order file if specified
    if (DDSIP_param->cpxorder)
    {
        status = CPXsetintparam (DDSIP_env, CPX_PARAM_MIPORDIND, 1);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to set cplex parameter CPX_PARAM_MIPORDIND.\n");
            return status;
        }
        if (DDSIP_ReadCPLEXOrder ())
            goto TERMINATE;
    }

    // Detect first and second stage constraints
    if ((status = DDSIP_DetectStageRows ()))
    {
        fprintf (stderr, "ERROR: Failed detection of row stages (BbInit)\n");
        fprintf (DDSIP_outfile, "ERROR: Failed detection of row stages (BbInit)\n");
        goto TERMINATE;
    }

    DDSIP_bb->DDSIP_step =  solve;

#ifdef CONIC_BUNDLE
    if (DDSIP_param->cb)
    {
        status = DDSIP_NonAnt ();
        if (status)
            goto TERMINATE;
    }
#endif

    if (DDSIP_param->outlev)
    {
        printf ("\t Total initialization time: %4.2f seconds.\n", DDSIP_GetCpuTime ());
        fprintf (DDSIP_bb->moreoutfile, " Total initialization time: %4.2f seconds.\n", DDSIP_GetCpuTime ());
    }
    fprintf (DDSIP_outfile, " Total initialization time: %4.2f seconds.\n", DDSIP_GetCpuTime ());
#ifndef _WIN32
    sprintf (astring, "grep 'MHz' /proc/cpuinfo|sort -k 4 -rn|head -1 >> %s\n", DDSIP_outfname);
    i = system (astring);
#endif

    // at the start there is no solution
    status = CPXsetintparam (DDSIP_env, CPX_PARAM_ADVIND, 0);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to set cplex parameter CPX_PARAM_ADVIND.\n");
        fprintf (DDSIP_outfile, "ERROR: Failed to set cplex parameter CPX_PARAM_ADVIND.\n");
        return status;
    }

#ifndef NEOS
    // Write deterministic equivalent (only expectation-based cases)
    if (DDSIP_param->write_detequ)
        DDSIP_DetEqu ();
#endif
    if (DDSIP_param->riskmod < 0)
    {
        fprintf (stderr, " pure risk models are disabled for now, exiting.\n");
        fprintf (DDSIP_outfile, " pure risk models are disabled for now, exiting.\n");
        goto TERMINATE;
    }
    if (DDSIP_param->stoccost && DDSIP_param->riskmod)
    {
        fprintf (DDSIP_outfile, "XXX Error: Risk optimization not implemented for stochastic cost coefficients.\n");
        fprintf (stderr, "XXX Error: Risk optimization not implemented for stochastic cost coefficients.\n");
        goto TERMINATE;
    }

    // Read advanced starting info
    if (DDSIP_param->advstart)
    {
        DDSIP_node[DDSIP_bb->curnode]->step = DDSIP_bb->DDSIP_step = adv;

        if ((status = DDSIP_AdvStart ()))
            goto TERMINATE;

        // Solution provided?
        if (DDSIP_param->advstart == 2)
        {
            DDSIP_bb->curnode = 0;
            DDSIP_bb->sug[DDSIP_param->nodelim + 2] = (struct sug_l *)DDSIP_Alloc (sizeof (sug_t), 1, "sug[0](Advanced start solution)");
            (DDSIP_bb->sug[DDSIP_param->nodelim + 2])->firstval =
                (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar, "sug[0]->firstval(adv start)");
            for (i = 0; i < (DDSIP_bb->firstvar); i++)
                (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[i] = DDSIP_bb->adv_sol[i];
            (DDSIP_bb->sug[DDSIP_param->nodelim + 2])->next = NULL;
            if ((status = DDSIP_UpperBound (DDSIP_param->scenarios, 0)) && status < 100000)
                goto TERMINATE;
        }
        // Only weak consistency check:
        if ((DDSIP_node[0]->bound) > DDSIP_bb->bestvalue)
        {
            status = 121;
            goto TERMINATE;
        }
        fprintf (DDSIP_outfile, "-----------------------------------------------------------\n\n");
    }
    // Solve EV and EEV if specified
    if (DDSIP_param->expected)
    {
        DDSIP_node[DDSIP_bb->curnode]->step = DDSIP_bb->DDSIP_step = eev;

        // status = 1 -> ExpValProb infeasible
        if ((status = DDSIP_ExpValProb ()) == 1)
        {
            if (DDSIP_param->outlev)
                fprintf (DDSIP_bb->moreoutfile, "    exp. val. prob: INFEASIBLE\n");
            fprintf (DDSIP_outfile, "    exp. val. prob: INFEASIBLE\n");
        }

        if (!status)
        {
            DDSIP_bb->curnode = 0;
            DDSIP_bb->sug[DDSIP_param->nodelim + 2] = (struct sug_l *)DDSIP_Alloc (sizeof (sug_t), 1, "sug[0](Advanced start solution)");
            DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval =
                (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar, "sug[i]->firstval(Heuristic)");
            for (i = 0; i < DDSIP_bb->firstvar; i++)
                (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[i] = DDSIP_bb->adv_sol[i];
            DDSIP_bb->sug[DDSIP_param->nodelim + 2]->next = NULL;
            if ((status = DDSIP_UpperBound (DDSIP_param->scenarios, 0)) && status < 100000)
                goto TERMINATE;
        }

        if (fabs (DDSIP_bb->expbest) < DDSIP_infty)
        {
            if (DDSIP_param->outlev)
                printf ("\t\t EEV:  %13.7f\n", DDSIP_bb->expbest);
            fprintf (DDSIP_outfile, "-EEV:     %13.7f\n", DDSIP_bb->expbest);
        }
        else
        {
            if (DDSIP_param->outlev)
                printf ("\t\t EEV:  No solution found.\n");
            fprintf (DDSIP_outfile, "-EEV:     No solution found.\n");
        }
    }                // END if (EV)
    // stop here if NODELIM is zero
    if (!(DDSIP_param->nodelim))
        goto TERMINATE;
    // Print cplex log to debugfile
    if (DDSIP_param->outlev > 51)
    {
#ifdef CPLEX_12_8
        if ((status = CPXsetlogfilename (DDSIP_env, DDSIP_moreoutfname, "a")))
            goto TERMINATE;
#else
        if ((status = CPXsetlogfile (DDSIP_env, DDSIP_bb->moreoutfile)))
            goto TERMINATE;
#endif
    }

    // No of DDSIP iterations
    DDSIP_bb->noiter = 0;
    // cont = 1 if no stop criteria is fullfilled
    cont = 1;
    // comb tells in case of a combined heuristic, which one to apply (3 = RoundNear)
    comb = 3;

    if (!DDSIP_param->nodelim)
        return 0;
    if (DDSIP_param->outlev)
        printf ("Starting branch-and-bound algorithm.\n");
    fprintf (DDSIP_outfile, "----------------------------------------------------------------------------------------\n");

    while (cont)
    {
        // the cuts from the root node are contained in every following node model, there is no need to check their violation
        // for the scenario solutions. But the rounding heuristics could violate a cut, so keep them.
        result = DDSIP_bb->cutNumber;
#ifdef CONIC_BUNDLE
#ifdef DEBUG
////////////////////////////////////////////
        if (DDSIP_bb->curnode && DDSIP_param->outlev)
        {
            if((DDSIP_node[DDSIP_bb->curnode-1])->step == dual)
                fprintf(DDSIP_bb->moreoutfile, "######## last node %d step=dual, leaf= %d,  dualdescitcnt = %d, DDSIP_bb->cutoff= %d\n",DDSIP_bb->curnode-1,(DDSIP_node[DDSIP_bb->curnode-1])->leaf,DDSIP_bb->dualdescitcnt,DDSIP_bb->cutoff);
        }
////////////////////////////////////////////
#endif

        if (DDSIP_param->outlev > 20)
        {
            if (DDSIP_node[DDSIP_bb->curnode]->cbReturn32)
                fprintf (DDSIP_bb->moreoutfile, "########## DDSIP_node[%d]->cbReturn32 = %d -> no ConicBundle\n", DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->cbReturn32);
        }
        // Dual method
        if ((!DDSIP_node[DDSIP_bb->curnode]->cbReturn32) &&
                ((DDSIP_param->cb > 0 && (!(DDSIP_bb->curnode % abs(DDSIP_param->cb))) && (abs(DDSIP_param->riskmod) != 4 || DDSIP_bb->curnode)) ||
                 (DDSIP_param->cb < 0 && (DDSIP_bb->curnode || (!DDSIP_bb->curnode && DDSIP_param->cbrootitlim >= 0)) &&
                                         (((DDSIP_node[DDSIP_bb->curnode]->depth <= DDSIP_param->cb_depth) && abs(DDSIP_param->riskmod) != 4) ||
                                          (abs(DDSIP_param->riskmod) == 5 && DDSIP_node[DDSIP_bb->curnode]->depth == 8) ||
                                          (DDSIP_bb->curnode > DDSIP_param->cbBreakIters &&
                                           ((!(DDSIP_bb->curnode % abs(DDSIP_param->cb))) ||
                                            (!((DDSIP_bb->curnode+1) % -DDSIP_param->cb)) ||
                                            (DDSIP_bb->curnode%200 > 199 - DDSIP_param->cbContinuous) ||
                                            (DDSIP_bb->curnode <= DDSIP_param->cbContinuous + DDSIP_param->cbBreakIters) ||
                                            ((DDSIP_bb->curnode  > 2*DDSIP_param->cbBreakIters) && (DDSIP_bb->curnode <= DDSIP_param->cbContinuous + 2*DDSIP_param->cbBreakIters)) ||
                                            ((DDSIP_bb->cutoff > 6) &&
                                             (((DDSIP_bb->no_reduced_front < 51) && (DDSIP_bb->curnode % -DDSIP_param->cb) < DDSIP_param->cbContinuous)
                                              || (((DDSIP_node[DDSIP_bb->curnode-1])->step == dual) && (DDSIP_node[DDSIP_bb->curnode-1])->leaf /*&& (DDSIP_bb->dualdescitcnt < 11)*/)
                                             )
                                            )
                                           )
                                          ) ||
                                          (abs(DDSIP_param->riskmod) != 5 && DDSIP_bb->curnode <= DDSIP_param->cbBreakIters && DDSIP_bb->curnode > DDSIP_param->cbBreakIters*.48 &&
                                           (CBFORALL || (DDSIP_node[DDSIP_bb->curnode]->numInheritedSols > (DDSIP_Imin(DDSIP_param->scenarios/20,2)+(DDSIP_param->scenarios+1)/2))))
                                         )
                 )
                )
           )
        {
            DDSIP_bb->DDSIP_step = DDSIP_node[DDSIP_bb->curnode]->step = dual;
            if ((status = DDSIP_DualOpt (0)))
            {
                if (!DDSIP_bb->curnode)
                    goto TERMINATE;
                else
                {
                    DDSIP_bb->skip = 1;
                    fprintf (DDSIP_outfile, "### return code of DualOpt is %d.\n", status);
                    if (status == 2 || DDSIP_bb->skip == 1)
                        DDSIP_node[DDSIP_bb->curnode]->bound = DDSIP_infty;
                }
            }
        }
        else
        {
            DDSIP_node[DDSIP_bb->curnode]->step = DDSIP_bb->DDSIP_step = solve;
            // status=1 means there was no solution found to a scenario problem
            if ((status = DDSIP_LowerBound (&tmp_objval, 0)))
                goto TERMINATE;
        }
#else
        DDSIP_node[DDSIP_bb->curnode]->step = DDSIP_bb->DDSIP_step = solve;
        // status=1 means there was no solution found to a subproblem
        if ((status = DDSIP_LowerBound (&tmp_objval, 0)))
            goto TERMINATE;
#endif

        if (!DDSIP_bb->skip || DDSIP_bb->skip == -1 || DDSIP_bb->skip == -11)
        {
            double old_bound;
            int cntr, maxCntr, nowAdded;

            DDSIP_bb->cutAdded = DDSIP_bb->cutNumber - result;
            nowAdded = DDSIP_bb->cutNumber;
            DDSIP_EvaluateScenarioSolutions (&comb);
            nowAdded = DDSIP_bb->cutNumber - nowAdded;
            cntr = 0;
            if (DDSIP_node[DDSIP_bb->curnode]->step != dual)
                maxCntr = DDSIP_param->numberReinits;
            else
                maxCntr = 0;
            old_bound = DDSIP_node[DDSIP_bb->curnode]->bound;
            boundstat = DDSIP_Bound ();
            if (!DDSIP_bb->curnode)
            {
                int cnt, j, noIncreaseCntr = 0;
                double lhs;
                cutpool_t *currentCut;
                DDSIP_PrintState (DDSIP_bb->noiter);
                if (nowAdded && DDSIP_param->outlev)
                {
                    fprintf (DDSIP_outfile, " %6d%101d cuts (UB)\n", DDSIP_bb->curnode, nowAdded);
                }
                while ((DDSIP_bb->cutAdded || DDSIP_node[0]->step == dual) && cntr < maxCntr)
                {
                    cntr++;
                    old_bound = DDSIP_node[0]->bound;
                    // Free the solutions from former LowerBound
                    for (i = 0; i < DDSIP_param->scenarios; i++)
                    {
                        if (((DDSIP_node[0])->first_sol)[i])
                        {
                            if (DDSIP_node[0]->step == solve)
                            {
                                currentCut = DDSIP_bb->cutpool;
                                while (currentCut)
                                {
                                    lhs = 0.;
                                    for (j = 0; j < DDSIP_bb->firstvar; j++)
                                    {
                                        lhs += (DDSIP_node[0])->first_sol[i][j] * currentCut->matval[j];
                                    }
                                    if (lhs < currentCut->rhs - 1.e-7)
                                    {
#ifdef DEBUG
                                        if (DDSIP_param->outlev > 50)
                                            fprintf (DDSIP_bb->moreoutfile, "scen %d solution violates cut %d.\n", i+1, currentCut->number);
#endif
                                        if ((cnt = (int) ((((DDSIP_node[0])->first_sol)[i])[DDSIP_bb->firstvar] - 0.9)))
                                            for (j = i + 1; cnt && j < DDSIP_param->scenarios; j++)
                                            {
                                                {
                                                    if (((DDSIP_node[0])->first_sol)[j]
                                                            && ((DDSIP_node[0])->first_sol)[i] == ((DDSIP_node[0])->first_sol)[j])
                                                    {
                                                        ((DDSIP_node[0])->first_sol)[j] = NULL;
                                                        cnt--;
                                                    }
                                                }
                                            }
                                        DDSIP_Free ((void **) &(((DDSIP_node[0])->first_sol)[i]));
                                        break;
                                    }
                                    currentCut = currentCut->prev;
                                }
                            }
                            else
                            {
                                if ((cnt = (int) ((((DDSIP_node[0])->first_sol)[i])[DDSIP_bb->firstvar] - 0.9)))
                                    for (j = i + 1; cnt && j < DDSIP_param->scenarios; j++)
                                    {
                                        {
                                            if (((DDSIP_node[0])->first_sol)[j]
                                                    && ((DDSIP_node[0])->first_sol)[i] == ((DDSIP_node[0])->first_sol)[j])
                                            {
                                                ((DDSIP_node[0])->first_sol)[j] = NULL;
                                                cnt--;
                                            }
                                        }
                                    }
                                DDSIP_Free ((void **) &(((DDSIP_node[0])->first_sol)[i]));
                            }
                        }
                    }
                    result = DDSIP_bb->cutNumber;
                    DDSIP_node[0]->step = DDSIP_bb->DDSIP_step = solve;
                    // status=1 means there was no solution found to a scenario problem
#ifdef GREATERMIPGAP
                    if ((status = DDSIP_LowerBound (&tmp_objval, 1)))
                        goto TERMINATE;
#else
                    if ((status = DDSIP_LowerBound (&tmp_objval, 0)))
                        goto TERMINATE;
#endif
                    nowAdded = DDSIP_bb->cutAdded = DDSIP_bb->cutNumber - result;
                    DDSIP_EvaluateScenarioSolutions (&comb);
                    nowAdded = DDSIP_bb->cutAdded - nowAdded;
                    DDSIP_bb->bestbound = DDSIP_node[0]->bound;
                    DDSIP_bb->noiter++;
                    DDSIP_PrintState (1);
                    if (/*nowAdded &&*/ DDSIP_param->outlev)
                    {
                        fprintf (DDSIP_outfile, " %6d %82d. reinit: %8d cuts (UB)\n", 0, cntr, nowAdded);
                    }
                    if ((DDSIP_node[0]->bound - old_bound)/(fabs(old_bound) + 1e-16) < 1.e-10)
                        noIncreaseCntr++;
                    if (noIncreaseCntr > 2 ||
                            (DDSIP_bb->bestvalue - DDSIP_node[0]->bound)/(fabs(DDSIP_bb->bestvalue) + 1e-16) < 0.5*DDSIP_param->relgap)
                        break;
                }
                if (DDSIP_param->deleteRedundantCuts)
                    DDSIP_CheckRedundancy(1);
// #ifdef GREATERMIPGAP
//                 // in case the computations werde done with reduced accuracy, repeat
//                 DDSIP_bb->cutAdded = 1;
//                 while (DDSIP_bb->cutAdded && cntr < maxCntr)
//                 {
//                     cntr++;
//                     old_bound = DDSIP_node[0]->bound;
//                     // Free the solutions from former LowerBound
//                     for (i = 0; i < DDSIP_param->scenarios; i++)
//                     {
//                         if (((DDSIP_node[0])->first_sol)[i])
//                         {
//                             double gap;
//                             gap = 100.0*((DDSIP_node[DDSIP_bb->curnode]->cursubsol)[i]-(DDSIP_node[DDSIP_bb->curnode]->subbound)[i])/
//                                   (fabs((DDSIP_node[DDSIP_bb->curnode]->cursubsol)[i])+1e-4);
//                             if (fabs(gap) < 1.e-10)
//                             {
//                                 currentCut = DDSIP_bb->cutpool;
//                                 while (currentCut)
//                                 {
//                                     lhs = 0.;
//                                     for (j = 0; j < DDSIP_bb->firstvar; j++)
//                                     {
//                                         lhs += (DDSIP_node[0])->first_sol[i][j] * currentCut->matval[j];
//                                     }
//                                     if (lhs < currentCut->rhs - 1.e-7)
//                                     {
// #ifdef DEBUG
//                                         if (DDSIP_param->outlev > 50)
//                                             fprintf (DDSIP_bb->moreoutfile, "scen %d solution violates cut %d.\n", i+1, currentCut->number);
// #endif
//                                         if ((cnt = (int) ((((DDSIP_node[0])->first_sol)[i])[DDSIP_bb->firstvar] - 0.9)))
//                                             for (j = i + 1; cnt && j < DDSIP_param->scenarios; j++)
//                                             {
//                                                 {
//                                                     if (((DDSIP_node[0])->first_sol)[j]
//                                                             && ((DDSIP_node[0])->first_sol)[i] == ((DDSIP_node[0])->first_sol)[j])
//                                                     {
//                                                         ((DDSIP_node[0])->first_sol)[j] = NULL;
//                                                         cnt--;
//                                                     }
//                                                 }
//                                             }
//                                         DDSIP_Free ((void **) &(((DDSIP_node[0])->first_sol)[i]));
//                                         break;
//                                     }
//                                     currentCut = currentCut->prev;
//                                 }
//                             }
//                             else
//                             {
//                                 if ((cnt = (int) ((((DDSIP_node[0])->first_sol)[i])[DDSIP_bb->firstvar] - 0.9)))
//                                     for (j = i + 1; cnt && j < DDSIP_param->scenarios; j++)
//                                     {
//                                         {
//                                             if (((DDSIP_node[0])->first_sol)[j]
//                                                     && ((DDSIP_node[0])->first_sol)[i] == ((DDSIP_node[0])->first_sol)[j])
//                                             {
//                                                 ((DDSIP_node[0])->first_sol)[j] = NULL;
//                                                 cnt--;
//                                             }
//                                         }
//                                     }
//                                 DDSIP_Free ((void **) &(((DDSIP_node[0])->first_sol)[i]));
//                             }
//                         }
//                     }
//                     result = DDSIP_bb->cutNumber;
//                     DDSIP_node[0]->step = DDSIP_bb->DDSIP_step = solve;
//                     // status=1 means there was no solution found to a scenario problem
//                     if ((status = DDSIP_LowerBound (&tmp_objval, 0)))
//                         goto TERMINATE;
//                     nowAdded = DDSIP_bb->cutAdded = DDSIP_bb->cutNumber - result;
//                     DDSIP_EvaluateScenarioSolutions (&comb);
//                     nowAdded = DDSIP_bb->cutAdded - nowAdded;
//                     DDSIP_bb->bestbound = DDSIP_node[0]->bound;
//                     DDSIP_bb->noiter++;
//                     DDSIP_PrintState (1);
//                     if (/*nowAdded &&*/ DDSIP_param->outlev)
//                     {
//                         fprintf (DDSIP_outfile, " %6d %82d. reinit: %8d cuts (UB)\n", 0, cntr, nowAdded);
//                     }
//                     if ((DDSIP_node[0]->bound - old_bound)/(fabs(old_bound) + 1e-16) < 1.e-10)
//                         noIncreaseCntr++;
//                     if (noIncreaseCntr > 2 ||
//                             (DDSIP_bb->bestvalue - DDSIP_node[0]->bound)/(fabs(DDSIP_bb->bestvalue) + 1e-16) < 0.5*DDSIP_param->relgap)
//                         break;
//                 }
// #endif
                for (i = 0; i < DDSIP_param->scenarios; i++)
                {
                    tmp_objval = 100.0*((DDSIP_node[DDSIP_bb->curnode]->cursubsol)[i]-(DDSIP_node[DDSIP_bb->curnode]->subbound)[i])/
                          (fabs((DDSIP_node[DDSIP_bb->curnode]->cursubsol)[i])+1e-4);
                    if (tmp_objval > DDSIP_param->relgap)
                    {
                        (DDSIP_node[DDSIP_bb->curnode]->mipstatus)[i] = CPXMIP_NODE_LIM_FEAS;
                    }
                }
            }
            else
            {
                // Print a line of output at the first, the last and each `ith' node
                if (!DDSIP_bb->noiter || (!((DDSIP_bb->noiter + 1) % DDSIP_param->logfreq)) || (DDSIP_param->outlev && (DDSIP_node[DDSIP_bb->curnode])->step == dual))
                    DDSIP_PrintState (DDSIP_bb->noiter);
                if (nowAdded && DDSIP_param->outlev > 1)
                {
                    fprintf (DDSIP_outfile, " %6d%101d cuts (UB)\n", DDSIP_bb->curnode, nowAdded);
                }
            }
        }
        else
        {
            boundstat = DDSIP_Bound ();
            // Print a line of output at the first, the last and each `ith' node
            if (!DDSIP_bb->noiter || (!((DDSIP_bb->noiter + 1) % DDSIP_param->logfreq)) || (DDSIP_param->outlev && (DDSIP_node[DDSIP_bb->curnode])->step == dual))
                DDSIP_PrintState (DDSIP_bb->noiter);
        }
        if (!DDSIP_bb->curnode)
        {
            DDSIP_bb->cutCntr0 = DDSIP_bb->cutNumber;
            DDSIP_bb->correct_bounding = fabs(DDSIP_node[0]->bound)*5.e-10;
        }
        else if (DDSIP_param->alwaysBendersCuts)
        {
            if (DDSIP_bb->curnode == 11)
            {
                if (DDSIP_bb->cutCntr0 == DDSIP_bb->cutNumber)
                {
                    if (DDSIP_param->outlev > 20)
                        fprintf (DDSIP_bb->moreoutfile, "### setting alwaysBendersCuts to 0\n");
                    DDSIP_param->alwaysBendersCuts = 0;
                }
            }
            if (DDSIP_bb->curnode == 24)
            {
                if (DDSIP_bb->cutCntr0 == DDSIP_bb->cutNumber)
                {
                    if (DDSIP_param->outlev > 20)
                        fprintf (DDSIP_bb->moreoutfile, "### setting alwaysBendersCuts to 0\n");
                    DDSIP_param->alwaysBendersCuts = 0;
                }
                else
                    DDSIP_bb->cutCntr0 = DDSIP_bb->cutNumber;
            }
            else if (DDSIP_bb->curnode == 49)
            {
                if (DDSIP_bb->cutCntr0 == DDSIP_bb->cutNumber)
                {
                    if (DDSIP_param->outlev > 20)
                        fprintf (DDSIP_bb->moreoutfile, "### setting alwaysBendersCuts to 0\n");
                    DDSIP_param->alwaysBendersCuts = 0;
                }
                else
                    DDSIP_bb->cutCntr0 = DDSIP_bb->cutNumber;
            }
            else if (DDSIP_bb->curnode == 99)
            {
                if (DDSIP_bb->cutCntr0 == DDSIP_bb->cutNumber)
                {
                    if (DDSIP_param->outlev > 20)
                        fprintf (DDSIP_bb->moreoutfile, "### setting alwaysBendersCuts to 0\n");
                    DDSIP_param->alwaysBendersCuts = 0;
                }
            }
        }
        DDSIP_bb->correct_bounding = fabs(DDSIP_bb->bestbound)*5.e-10;

        if (DDSIP_param->cb && DDSIP_param->cb_depth>0 && DDSIP_param->cb_depth_iters && DDSIP_bb->curnode == 2 && DDSIP_node[1]->dual && DDSIP_node[2]->dual)
        {
            double old_objval, *save_dual;
            int cnt = 0, k1;
            for (k1 = 0; k1 < DDSIP_bb->dimdual; k1++)
            {
                if (DDSIP_node[2]->dual[k1] != DDSIP_node[1]->dual[k1])
                {
                    cnt = 1;
                    break;
                }
            }
            if (cnt)
            {
                save_dual = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->dimdual + 3, "save_dual(DDSIP_main)");
                memcpy (save_dual, DDSIP_node[1]->dual, sizeof (double) * DDSIP_bb->dimdual + 3);
                memcpy (DDSIP_node[1]->dual, DDSIP_node[2]->dual, sizeof (double) * DDSIP_bb->dimdual + 3);
                DDSIP_bb->curnode = 1;
                status = DDSIP_SetBounds ();
                if (status)
                {
                    fprintf (stderr, "ERROR: Failed to add new bounds (Branch).\n");
                    return status;
                }
                DDSIP_bb->DDSIP_step = solve;
                for (i = 0; i < DDSIP_param->scenarios; i++)
                {
                    if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[i] != NULL)
                    {
                        if ((cnt = (int) ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i])[DDSIP_bb->firstvar] - 0.9)))
                            for (k1 = i + 1; cnt && k1 < DDSIP_param->scenarios; k1++)
                            {
                                if ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i]) == (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[k1]))
                                {
                                    ((DDSIP_node[DDSIP_bb->curnode])->first_sol)[k1] = NULL;
                                    cnt--;
                                }
                            }
                        DDSIP_Free ((void **) &(((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i]));
                    }
                }
                old_objval = DDSIP_node[1]->bound;
                DDSIP_bb->dualObjVal = DDSIP_node[1]->bound;
                DDSIP_LowerBound (&tmp_objval, 0);
                if (DDSIP_param->outlev > 10)
                    fprintf (DDSIP_bb->moreoutfile, " ## re-eval node %d : %20.14g (new) - %20.14g (old) = %g\n", DDSIP_bb->curnode, tmp_objval, old_objval, tmp_objval - old_objval);
                if (tmp_objval <= old_objval)
                {
                    memcpy (DDSIP_node[1]->dual, save_dual, sizeof (double) * DDSIP_bb->dimdual + 3);
                    DDSIP_node[1]->bound = old_objval;
                    DDSIP_bb->DDSIP_step = DDSIP_node[1]->step = dual;
                }
                else
                {
                    DDSIP_node[1]->bound = tmp_objval;
                    DDSIP_bb->DDSIP_step = DDSIP_node[1]->step = solve;
#ifdef DE1MULT
                    // delete node 1 multiplier from list
                    bbest_t * tmp_previous, * tmp_bestdual;
                    tmp_bestdual = DDSIP_bb->bestdual;
                    if (tmp_bestdual)
                    {
                        if (tmp_bestdual->node_nr == 1)
                        {
                            DDSIP_bb->bestdual = DDSIP_bb->bestdual->next;
                            DDSIP_bb->bestdual_cnt--;
                            if (DDSIP_param->outlev > 10)
                                fprintf (DDSIP_bb->moreoutfile, " ## delete entry from node %d from bestdual list, #entries: %d\n", tmp_bestdual->node_nr, DDSIP_bb->bestdual_cnt);
                            DDSIP_Free ((void **) &(tmp_bestdual->dual));
                            DDSIP_Free ((void **) &(tmp_bestdual));
                        }
                        else
                        {
                            while (tmp_bestdual)
                            {
                                tmp_previous = tmp_bestdual;
                                tmp_bestdual = tmp_bestdual->next;
                                if (tmp_bestdual && tmp_bestdual->node_nr == 1)
                                {
                                    tmp_previous->next = tmp_bestdual->next;
                                    DDSIP_bb->bestdual_cnt--;
                                    if (DDSIP_param->outlev > 10)
                                        fprintf (DDSIP_bb->moreoutfile, " ## delete entry from node %d from bestdual list, #entries: %d\n", tmp_bestdual->node_nr, DDSIP_bb->bestdual_cnt);
                                    DDSIP_Free ((void **) &(tmp_bestdual->dual));
                                    DDSIP_Free ((void **) &(tmp_bestdual));
                                    break;
                                }
                            }
                        }
                    }
#endif
                }
                DDSIP_Free ((void **) &(save_dual));
                boundstat = DDSIP_Bound ();
                if (!DDSIP_bb->skip)
                    DDSIP_EvaluateScenarioSolutions (&comb);
                DDSIP_PrintState (DDSIP_bb->noiter);
                DDSIP_bb->curnode = 2;
            }
        }
        // check if the bestdual multipliers improve bounds of the nodes in the front
        if (DDSIP_param->cb < 0 && DDSIP_param->cb_checkBestdual && DDSIP_bb->bestdual_cnt &&
                  (DDSIP_bb->curnode == DDSIP_Imin (2*((DDSIP_param->cbContinuous + DDSIP_param->cbBreakIters)/2) + DDSIP_Imax (2, 2*(-DDSIP_param->cb/2-1)), DDSIP_param->dive_start - 7) ||
                   DDSIP_bb->curnode == DDSIP_Imax (2*(DDSIP_param->cbContinuous + 2*(DDSIP_param->cbBreakIters/2) + 4), 62) ||
                   DDSIP_bb->curnode == DDSIP_param->nodelim - 2
                  )
                )
        {
            if (DDSIP_bb->bestdual_improvement > 0)
            {
                int ih, inode, nowAdded;

                double * front_node_bound, least_bound = DDSIP_infty;
                if (DDSIP_param->outlev > 20)
                    fprintf (DDSIP_bb->moreoutfile, " ## node %d, bestdual_improvement= %d, re-evaluation of front nodes with bestdual entries\n", DDSIP_bb->curnode, DDSIP_bb->bestdual_improvement); 
                // sort front nodes ascending wrt. bound
                DDSIP_Free ((void **) &DDSIP_bb->front_nodes_sorted);
                DDSIP_bb->front_nodes_sorted = (int *) DDSIP_Alloc(sizeof(int), DDSIP_bb->nofront, "DDSIP_bb->front_nodes_sorted(main)");
                front_node_bound = (double *) DDSIP_Alloc(sizeof(double), DDSIP_bb->nonode, "front_node_bound(main)");
                for  (ih = 0; ih < DDSIP_bb->nofront; ih++)
                {
                    DDSIP_bb->front_nodes_sorted[ih] = DDSIP_bb->front[ih];
                    front_node_bound[DDSIP_bb->front_nodes_sorted[ih]] = DDSIP_node[DDSIP_bb->front_nodes_sorted[ih]]->bound;
                }
                DDSIP_qsort_ins_A (front_node_bound, DDSIP_bb->front_nodes_sorted, 0, DDSIP_bb->nofront-1);
                tmp_objval = DDSIP_node[DDSIP_bb->front_nodes_sorted[DDSIP_bb->nofront - 1]]->bound;
                ih = DDSIP_bb->curnode;
                for (inode = 0; inode < DDSIP_bb->nofront; inode++)
                {
                    if (DDSIP_node[DDSIP_bb->front_nodes_sorted[inode]]->leaf)
                        continue;
                    DDSIP_bb->curnode = DDSIP_bb->front_nodes_sorted[inode];
                    if (DDSIP_node[DDSIP_bb->curnode]->bound > 0.30 * least_bound + 0.70 * tmp_objval)
                    {
                        break;
                    }
                    if (!DDSIP_node[DDSIP_bb->curnode]->leaf && DDSIP_node[DDSIP_bb->curnode]->solved)
                    {
                        double old_objval;
                        int cnt, k1;
                        nowAdded = DDSIP_bb->cutNumber;
                        status = DDSIP_SetBounds ();
                        if (status)
                        {
                            fprintf (stderr, "ERROR: Failed to add new bounds (main).\n");
                            return status;
                        }
                        for (i = 0; i < DDSIP_param->scenarios; i++)
                        {
                            if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[i] != NULL)
                            {
                                if ((cnt = (int) ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i])[DDSIP_bb->firstvar] - 0.9)))
                                    for (k1 = i + 1; cnt && k1 < DDSIP_param->scenarios; k1++)
                                    {
                                        if ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i]) == (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[k1]))
                                        {
                                            ((DDSIP_node[DDSIP_bb->curnode])->first_sol)[k1] = NULL;
                                            cnt--;
                                        }
                                    }
                                DDSIP_Free ((void **) &(((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i]));
                            }
                        }
                        old_objval = DDSIP_node[DDSIP_bb->curnode]->bound;
                        DDSIP_bb->DDSIP_step = DDSIP_node[DDSIP_bb->curnode]->step = dual;
                        if ((status = DDSIP_DualOpt (1)))
                        {
                            if (!DDSIP_bb->curnode)
                                goto TERMINATE;
                            else
                            {
                                DDSIP_bb->skip = 1;
                                fprintf (DDSIP_outfile, "### return code of DualOpt is %d.\n", status);
                                if (status == 2 || DDSIP_bb->skip == 1)
                                    DDSIP_node[DDSIP_bb->curnode]->bound = DDSIP_infty;
                            }
                        }
                        if (DDSIP_param->outlev > 10)
                            fprintf (DDSIP_bb->moreoutfile, " ## re-eval %2d node %d : %20.14g (new) - %20.14g (old) = %g\n", inode+1, DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound, old_objval, DDSIP_node[DDSIP_bb->curnode]->bound - old_objval);
                        if (DDSIP_node[DDSIP_bb->curnode]->bound < old_objval)
                        {
                            DDSIP_node[DDSIP_bb->curnode]->bound = old_objval;
                        }
                        least_bound = DDSIP_Dmin (least_bound, DDSIP_node[DDSIP_bb->curnode]->bound);
                        if (!DDSIP_bb->skip)
                            DDSIP_EvaluateScenarioSolutions (&comb);
                        nowAdded = DDSIP_bb->cutNumber - nowAdded;
                        DDSIP_PrintState (DDSIP_bb->noiter);
                        if (nowAdded && DDSIP_param->outlev > 1)
                        {
                            fprintf (DDSIP_outfile, " %6d%101d cuts (UB)\n", DDSIP_bb->curnode, nowAdded);
                        }
                    }
                }
                DDSIP_bb->curnode = ih;
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile, "##### node %d: set bestBound to -6\n", DDSIP_bb->curnode);
                DDSIP_bb->bestBound = -6;
                for  (ih = 0; ih < DDSIP_bb->nofront; ih++)
                {
                    DDSIP_bb->front_nodes_sorted[ih] = DDSIP_bb->front[ih];
                    front_node_bound[DDSIP_bb->front_nodes_sorted[ih]] = DDSIP_node[DDSIP_bb->front_nodes_sorted[ih]]->bound;
                }
                DDSIP_qsort_ins_A (front_node_bound, DDSIP_bb->front_nodes_sorted, 0, DDSIP_bb->nofront-1);
                boundstat = DDSIP_Bound ();
                DDSIP_Free ((void **) &front_node_bound);
                cont = DDSIP_Continue (&DDSIP_bb->noiter, &boundstat);
                if (!cont)
                {
                    status = boundstat;
                    goto TERMINATE;
                }

                if ((status = DDSIP_Branch ()))
                    goto TERMINATE;
            }
            else if (DDSIP_param->outlev > 20)
                fprintf (DDSIP_bb->moreoutfile, " ## node %d, bestdual_improvement= %d, no re-evaluation of front nodes with bestdual entries\n", DDSIP_bb->curnode, DDSIP_bb->bestdual_improvement); 
        }
        else if ( (DDSIP_bb->curnode == DDSIP_param->dive_start - 7 ||
                   DDSIP_bb->curnode == 62 ||
                   DDSIP_bb->curnode == DDSIP_param->nodelim - 6
                  )
                )
        {
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile, "##### node %d: set bestBound to -6\n", DDSIP_bb->curnode);
                DDSIP_bb->bestBound = -6;
        }


        // DDSIP_bb->skip > 0 indicates that UpperBound calls have been skipped
        // DDSIP_bb->heurval contains the objective value of the heuristic solution
        // Reset to initial values
        DDSIP_bb->heurval = DDSIP_infty;
        DDSIP_bb->skip = 0;

        cont = DDSIP_Continue (&DDSIP_bb->noiter, &boundstat);
        if (!cont)
        {
            status = boundstat;
            goto TERMINATE;
        }

        if ((status = DDSIP_Branch ()))
            goto TERMINATE;

        if (DDSIP_param->deleteRedundantCuts && !(DDSIP_bb->curnode % 10))
            DDSIP_CheckRedundancy(1);
    }

    // Termination
TERMINATE:

    DDSIP_PrintErrorMsg (status);

    printf ("\nOutput files in directory `%s'.\n", DDSIP_outdir);

    // Free up problem as allocated above
    //
    if (DDSIP_node != NULL)
    {
        DDSIP_FreeFrontNodes ();
        for (i = 0; i < DDSIP_bb->nonode; i++)
            DDSIP_Free ((void **) &(DDSIP_node[i]));
        DDSIP_Free ((void **) &(DDSIP_node));
    }

    if (DDSIP_data != NULL)
    {
        DDSIP_FreeData ();
        DDSIP_Free ((void **) &(DDSIP_data));
    }

    // Free up the problem as allocated by CPXcreateprob, if necessary
    if (DDSIP_lp != NULL)
    {
        status = CPXfreeprob (DDSIP_env, &DDSIP_lp);
        if (status)
            fprintf (stderr, "ERROR: CPXfreeprob failed, error code %d\n", status);
    }
    // Free up the CPLEX environment
    if (DDSIP_env != NULL)
    {
        status = CPXcloseCPLEX (&DDSIP_env);
        if (status)
        {
            char errmsg[1024];
            fprintf (stderr, "ERROR: Failed to close CPLEX environment.\n");
            CPXgeterrorstring (DDSIP_env, status, errmsg);
            fprintf (stderr, "%s\n", errmsg);
        }
    }

    if (DDSIP_bb != NULL)
    {
        DDSIP_FreeBb ();
        DDSIP_Free ((void **) &(DDSIP_bb));
    }

    if (DDSIP_param->outlev)
        printf ("Terminating DDSIP.\n");

    if (DDSIP_param != NULL)
    {
        DDSIP_FreeParam ();
        DDSIP_Free ((void **) &(DDSIP_param));
    }

    fprintf (DDSIP_outfile, "current system time: ");
#ifndef _WIN32
    i = system ("date");
    // Print time to output file
    sprintf (astring, "date >> %s\n", DDSIP_outfname);
    i = system (astring);
#else
    sprintf (astring, "date /T >> %s & time /T >> %s\n", DDSIP_outfname,DDSIP_outfname);
    i = system (astring);
#endif

    if (DDSIP_outfile != NULL)
        fclose (DDSIP_outfile);

    return 0;
}
