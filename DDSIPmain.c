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
#include <sys/resource.h>


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
const double DDSIP_bigint   = 1.0e7;	   // Upper bound on integer parameters
const double DDSIP_bigvalue = 1.0e9;	   // Just to detect the print format
const double DDSIP_infty    = CPX_INFBOUND; // is 1.0e20; -- Infinity

// Version
const char DDSIP_version[] = "2017-08-12 (with Benders feasibility cuts) ";

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
main (void)
{
    // Temporary variables
    struct stat filestat;

    char astring[DDSIP_ln_fname];

    int status = 0, cont, boundstat, comb, i;

    // guarantee minimal stacksize limit
    const rlim_t kStackSize = 128 * 1024 * 1024;   // min stack size = 128 MB
    struct rlimit rl;
    int result;
    unsigned cur_limit;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {
        cur_limit = rl.rlim_cur;
        printf ("original stacksize limit %u\n", cur_limit);
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
                printf ("changed  stacksize limit %u\n", cur_limit);
            }
        }
    }

    //

    // Welcome
    printf ("########################################################");
    printf ("################################\n##########  ");
    printf ("D D S I P -- Dual Decomposition In Stochastic Integer Programming");
    printf ("  #########\n");
    printf ("##########  Version: %-57.57s #########\n", DDSIP_version);
    printf ("For copyright, license agreements, help, comments, requests, ... ");
    printf ("see\n\thttp://www.uni-duisburg-essen.de/~hn215go/ddsip.shtml\n");
    printf ("\thttp://www.www.github.com/RalfGollmer/ddsip\n");
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
    fprintf (DDSIP_outfile, "-----------------------------------------------------------\n");
    fprintf (DDSIP_outfile, "Current system time: ");
    fflush (DDSIP_outfile);
#ifndef _WIN32
    i = system ("date");
    // Print time to output file
    sprintf (astring, "date >> %s\n", DDSIP_outfname);
    i = system (astring);
#else
    sprintf (astring, "date /T >> %s & time /T >> %s\n", DDSIP_outfname,DDSIP_outfname);
    i = system (astring);
#endif

    fprintf (DDSIP_outfile, "\nDDSIP VERSION %s\n", DDSIP_version);
    fprintf (DDSIP_outfile, "-----------------------------------------------------------\n");

    // Open cplex environment
    DDSIP_env = CPXopenCPLEX (&status);
    if (DDSIP_env == NULL)
    {
        fprintf (stderr, "ERROR: Failed to open cplex environment, CPLEX error code %d.\n",status);
        fprintf (DDSIP_outfile, "ERROR: Failed to open cplex environment, CPLEX error code %d.\n",status);
        return status;
    }
    else
        printf ("CPLEX version is %s\n", CPXversion (DDSIP_env));
    fprintf (DDSIP_outfile, "CPLEX version is %s\n\n", CPXversion (DDSIP_env));

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
        if ( status ) {
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
        // Buffer size = 0
        setbuf (DDSIP_bb->moreoutfile, 0);
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
    }

    DDSIP_node = (node_t **) DDSIP_Alloc (sizeof (node_t *), DDSIP_param->nodelim + 3, "node(Main)");
    DDSIP_node[0] = (node_t *) DDSIP_Alloc (sizeof (node_t), 1, "node[0](Main)");
    DDSIP_node[0]->first_sol = (double **) DDSIP_Alloc (sizeof (double *), DDSIP_param->scenarios, "node[0]->first_sol(Main)");
    DDSIP_node[0]->cursubsol = (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "node[0]->cursubsol(Main)");
    DDSIP_node[0]->subbound = (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "node[0]->subbound(Main)");
    if (DDSIP_param->cb)
        DDSIP_node[0]->subboundNoLag = (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "node[0]->subbound_NoLag(Main)");
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
        if ((status = DDSIP_ExpValProb ()) > 1)
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
    }				// END if (EV)

    if (DDSIP_bb->cutAdded && DDSIP_param->outlev > 1)
    {
        fprintf (DDSIP_outfile, " %6d%101d cuts\n", DDSIP_bb->curnode, DDSIP_bb->cutAdded);
    }

    // Print cplex log to debugfile
    if (DDSIP_param->outlev > 51)
        if ((status = CPXsetlogfile (DDSIP_env, DDSIP_bb->moreoutfile)))
            goto TERMINATE;

    // No of DDSIP iterations
    DDSIP_bb->noiter = 0;
    // cont = 1 if no stop criteria is fullfilled
    cont = 1;
    // comb tells in case of a combined heuristic, which one to apply (3 = RoundNear)
    comb = 3;

    if (DDSIP_param->outlev)
        printf ("Starting branch-and-bound algorithm.\n");
    fprintf (DDSIP_outfile, "----------------------------------------------------------------------------------------\n");

    while (cont)
    {
        // the cuts from the root node are contained in every following node model, there is no need to check their violation
        // for the scenario solutions. But the rounding heuristics could violate a cut, so keep them.
#ifdef CONIC_BUNDLE
        // Dual method
        if ((DDSIP_param->cb > 0 && (!(DDSIP_bb->noiter % abs(DDSIP_param->cb))) && (abs(DDSIP_param->riskmod) != 4 || DDSIP_bb->noiter)) ||
                (DDSIP_param->cb < 0 && ((!(DDSIP_bb->noiter) && abs(DDSIP_param->riskmod) != 4 && DDSIP_param->cbrootitlim) ||
                                         (abs(DDSIP_param->riskmod) == 5 && DDSIP_node[DDSIP_bb->curnode]->depth == 8) ||
                                         (DDSIP_bb->noiter > DDSIP_param->cbBreakIters &&
                                           ((!(DDSIP_bb->noiter % abs(DDSIP_param->cb))) || (!((DDSIP_bb->noiter+1) % -DDSIP_param->cb)) ||
                                           (DDSIP_bb->noiter < DDSIP_param->cbContinuous + DDSIP_param->cbBreakIters) ||
                                           ((DDSIP_bb->cutoff > 5) && (DDSIP_bb->no_reduced_front < 51) && (DDSIP_bb->noiter % -DDSIP_param->cb) < DDSIP_param->cbContinuous))) ||
                                         (DDSIP_bb->noiter%200 > 199 - DDSIP_param->cbContinuous) ||
                                         (abs(DDSIP_param->riskmod) != 5 && (DDSIP_bb->noiter <= DDSIP_param->cbBreakIters && DDSIP_bb->noiter > DDSIP_param->cbBreakIters*.6 && (DDSIP_node[DDSIP_bb->curnode]->numInheritedSols > (DDSIP_Imin(DDSIP_param->scenarios/20,2)+(DDSIP_param->scenarios+1)/2))))
                                        )
                )
           )
        {
            DDSIP_bb->DDSIP_step = DDSIP_node[DDSIP_bb->curnode]->step = dual;
            if ((status = DDSIP_DualOpt ()))
            {
                if (!DDSIP_bb->curnode)
                    goto TERMINATE;
                else
                {
                    DDSIP_bb->skip = 1;
                    DDSIP_node[DDSIP_bb->curnode]->bound = DDSIP_infty;
                }
            }
        }
        else
        {
            DDSIP_node[DDSIP_bb->curnode]->step = DDSIP_bb->DDSIP_step = solve;
            // status=1 means there was no solution found to a scenario problem
            if ((status = DDSIP_LowerBound ()))
                goto TERMINATE;
        }
#else
        DDSIP_node[DDSIP_bb->curnode]->step = DDSIP_bb->DDSIP_step = solve;
        // status=1 means there was no solution found to a subproblem
        if ((status = LowerBound ()))
            goto TERMINATE;
#endif

        if (!DDSIP_bb->skip)
        {
            double old_bound;
            int cntr, maxCntr;

            DDSIP_EvaluateScenarioSolutions (&comb);
            cntr = 0;
            if (DDSIP_node[DDSIP_bb->curnode]->step != dual)
                maxCntr = DDSIP_param->numberReinits;
            else
                maxCntr = 0;
            old_bound = DDSIP_node[DDSIP_bb->curnode]->bound;
            boundstat = DDSIP_Bound ();
            if (!DDSIP_bb->curnode)
            {
                int cnt, j;
                double lhs;
                cutpool_t *currentCut;
                DDSIP_PrintState (DDSIP_bb->noiter);
                if (DDSIP_bb->cutAdded && DDSIP_param->outlev > 1)
                {
                    fprintf (DDSIP_outfile, " %6d%101d cuts\n", DDSIP_bb->curnode, DDSIP_bb->cutAdded);
                }
                while ((DDSIP_bb->cutAdded || DDSIP_node[0]->step == dual) && cntr < maxCntr)
                {
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
                                        if (DDSIP_param->outlev > 50)
                                            fprintf (DDSIP_bb->moreoutfile, "scen %d solution violates cut %d.\n", i+1, currentCut->number);
                                        if ((cnt = (((DDSIP_node[0])->first_sol)[i])[DDSIP_bb->firstvar] - 1))
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
                               if ((cnt = (((DDSIP_node[0])->first_sol)[i])[DDSIP_bb->firstvar] - 1))
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
                    DDSIP_node[0]->step = DDSIP_bb->DDSIP_step = solve;
                    // status=1 means there was no solution found to a scenario problem
                    if ((status = DDSIP_LowerBound ()))
                        goto TERMINATE;
                    DDSIP_EvaluateScenarioSolutions (&comb);
                    cntr++;
                    DDSIP_bb->bestbound = DDSIP_node[0]->bound;
                    DDSIP_PrintState (1);
                    if (DDSIP_bb->cutAdded && DDSIP_param->outlev > 1)
                    {
                        fprintf (DDSIP_outfile, " %6d %82d. reinit: %8d cuts\n", 0, cntr, DDSIP_bb->cutAdded);
                    }
                    DDSIP_bb->noiter++;
                    if ((DDSIP_node[0]->bound - old_bound)/(fabs(old_bound) + 1e-16) < 1.e-7)
                        break;
                }
            }
            else
            {
                // Print a line of output at the first, the last and each `ith' node
                if (!DDSIP_bb->noiter || !((DDSIP_bb->noiter + 1) % DDSIP_param->logfreq))
                        DDSIP_PrintState (DDSIP_bb->noiter);
                if (DDSIP_bb->cutAdded && DDSIP_param->outlev > 1)
                {
                    fprintf (DDSIP_outfile, " %6d%101d cuts\n", DDSIP_bb->curnode, DDSIP_bb->cutAdded);
                }
            }
        }
        else
        {
            boundstat = DDSIP_Bound ();
            // Print a line of output at the first, the last and each `ith' node
            if (!DDSIP_bb->noiter || !((DDSIP_bb->noiter + 1) % DDSIP_param->logfreq))
                DDSIP_PrintState (DDSIP_bb->noiter);
        }

        // DDSIP_bb->skip indicates that UpperBound calls have been skipped
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
    }

    // Termination
TERMINATE:

    // check if some of the cuts are redundant
    if (DDSIP_param->redundancyCheck)
    {
        DDSIP_CheckRedundancy();
    }


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
            printf ("%s\n", errmsg);
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

    fprintf (DDSIP_outfile, "Current system time: ");
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
