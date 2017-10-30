/*  Authors:           Andreas M"arkert, Ralf Gollmer
	Copyright to:      University of Duisburg-Essen
    Language:          C
	Description:
	Here we check all termination criteria and perform a postprocessing
	if one criterion is fulfilled.

	License:
	This file is part of DDSIP.

    DDSIP is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    DDSIP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with DDSIP; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

#include <DDSIP.h>
#include <DDSIPconst.h>

//==========================================================================
// Function checks limits and tolerances
// solution status
// -1  user termination
//  0  continue
//  1  node limit reached
//  2  gap reached
//  3  time limit reached
//  4  maximal dispersion = 0
//  5  branching tree exhausted
//  6  no valid lower bound
//  7  problem infeasible
//  8  problem unbounded
// CPLEX errors are returned in boundstat
int
DDSIP_Continue (int *noiter, int *boundstat)
{
    int solstat;
    char message[DDSIP_max_str_ln];
    double waitandsee = 0.l;
    double wall_secs, cpu_secs;
    int    wall_hrs, wall_mins, cpu_hrs, cpu_mins;
    double sum;

    (*noiter)++;
    // User termination
    // Control C
    if (DDSIP_killsignal)
    {
        sprintf (message, "\nTermination signal received.\n");
        solstat = -1;
    }
    // No valid lower bound
    // This means a scenario problem is infeasible or insufficiently DDSIP_evaluated (e.g. low time limit)
    // We do not proceed in this case
    else if (fabs (DDSIP_bb->bestbound - DDSIP_infty) < DDSIP_param->accuracy)
    {
        if (DDSIP_bb->bestvalue < DDSIP_infty)
        {
            sprintf (message, "\nNo valid lower bound.\n");
            solstat = 6;
        }
        else
        {
            sprintf (message, "\nProblem infeasible.\n");
            solstat = 7;
        }
    }
    else if (fabs (DDSIP_bb->bestbound + DDSIP_infty) < DDSIP_param->accuracy)
    {
        sprintf (message, "\nProblem unbounded.\n");
        solstat = 8;
    }
    // Regular termination
    // Gap
    else if ((DDSIP_bb->bestvalue - DDSIP_bb->bestbound <= DDSIP_param->absgap + DDSIP_param->accuracy)
             || ((DDSIP_bb->bestvalue - DDSIP_bb->bestbound) / (fabs(DDSIP_bb->bestvalue)+1e-16) <= DDSIP_param->relgap))
    {
        sprintf (message, "\nGap reached.\n");
        solstat = 2;
    }
    // Nodelimit
    else if (DDSIP_bb->nonode > DDSIP_param->nodelim + 1)
    {
        sprintf (message, "\nNode limit reached (total number of nodes).\n");
        solstat = 1;
    }
    // Timelimit
    else if (DDSIP_GetCpuTime () > DDSIP_param->timelim)
    {
        sprintf (message, "\nTime limit reached.\n");
        solstat = 3;
    }
    // Front is empty
    else if (*boundstat == 2)
    {
        sprintf (message, "\nThe whole branching tree was backtracked.\nProbably due to MIP gaps (see below) or branching on continuous variables the specified gap tolerance could not be reached.\n");
        solstat = 5;
    }
    // Maximal dispersion norm of front nodes is zero
    else if (*boundstat == 1)
    {
        sprintf (message, "\nMaximal dispersion equals zero.\n");
        solstat = 4;
    }
    // Continue
    else
        solstat = 0;

    // If the end is near
    if (solstat)
    {
        int i, j, feas = 0, bound = 0;
        FILE *outputfile;

        // Infeasible ?
        if (DDSIP_bb->bestvalue < DDSIP_infty)
            feas = 1;

        // Valid bound ?
        if (fabs (DDSIP_bb->bestbound) < DDSIP_infty)
            bound = 1;

        // Print a line of output here if this hasn't be done before
        if ((*noiter + 1) % DDSIP_param->logfreq)
            DDSIP_PrintState (*noiter);

        // Print message according to solstat
        time (&DDSIP_bb->cur_time);
        printf ("%s", message);
        if (DDSIP_param->outlev)
            fprintf (DDSIP_bb->moreoutfile, "\n%s\n", message);
        DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
        DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
        printf ("Total CPU time: %dh %02d:%05.2f, wall time: %dh %02d:%02.0f\n",cpu_hrs,cpu_mins,cpu_secs,wall_hrs,wall_mins,wall_secs);
        fprintf (DDSIP_outfile,"%s", message);
        fprintf (DDSIP_outfile,"Total CPU time: %dh %02d:%05.2f, wall time: %dh %02d:%02.0f\n",cpu_hrs,cpu_mins,cpu_secs,wall_hrs,wall_mins,wall_secs);

        // If a feasible solution exists
        if (feas)
        {
            int *index = (int *) DDSIP_Alloc (sizeof (int), DDSIP_param->scenarios, "quantil (Continue)");
            // Display lower bound of best node
            if ((solstat == 4 || solstat == 5) && DDSIP_bb->feasbound > DDSIP_bb->bestbound)
                printf ("Bound in best node %12.7f\n", DDSIP_bb->feasbound);

            // Some more info on the solution
            if (DDSIP_param->files)
            {
                if ((outputfile = fopen (DDSIP_solfname, "w")) == NULL)
                    fprintf (stderr, "*Warning: Failed to write output to %s.\n", DDSIP_solfname);
                else
                {
                    char **colname = (char **) DDSIP_Alloc (sizeof (char *), DDSIP_bb->novar,
                                                            "colname (continue)");
                    char *colstore = (char *) DDSIP_Alloc (sizeof (char),
                                                           DDSIP_bb->novar * DDSIP_ln_varname,
                                                           "colstore (continue)");

                    fprintf (outputfile, "This file contains additional information on first- and second-stage ");
                    fprintf (outputfile, "scenario solutions, namely\n1. the best first-stage solution,\n");
                    fprintf (outputfile, "2. the lower bounds in the root node ");
                    fprintf (outputfile, "vs. the upper bound in the best node,\n3. the quantiles of the ");
                    fprintf (outputfile, "probability distribution w.r.t. the best first-stage solution, ");
                    fprintf (outputfile, "and \n4. the best second-stage solutions associated with ");
                    fprintf (outputfile, "the best first-stage solution.\n");
                    fprintf (outputfile, "The item 4. is only printed if the parameter ");
                    fprintf (outputfile, "OUTLEV is greater than 2.\n\n1. Best Solution\n");
                    fprintf (outputfile, "Variable name                Value\n");
                    fprintf (outputfile, "-----------------------------------\n");

                    // Get all variable names
                    (*boundstat) = CPXgetcolname (DDSIP_env, DDSIP_lp, colname, colstore, DDSIP_bb->novar * DDSIP_ln_varname, &i, 0, DDSIP_bb->novar - 1);
                    if (*boundstat)
                    {
                        fprintf (stderr, "ERROR: Failed to get column names (Continue)\n");
                        for (i = 0; i < DDSIP_bb->firstvar; i++)
                            if (DDSIP_bb->firsttype[i] == 'C')
                                fprintf (outputfile, "%18.16g\n", DDSIP_bb->bestsol[i]);
                            else
                                fprintf (outputfile, "%18.16g\n", DDSIP_bb->bestsol[i]);
                    }
                    else
                    {
                        for (i = 0; i < DDSIP_bb->firstvar; i++)
                            if (DDSIP_bb->firsttype[i] == 'C')
                                fprintf (outputfile, "%-25s %-18.16g\n", colname[DDSIP_bb->firstindex[i]], DDSIP_bb->bestsol[i]);
                            else
                                fprintf (outputfile, "%-25s %-18.16g\n", colname[DDSIP_bb->firstindex[i]], DDSIP_bb->bestsol[i]);
                    }

                    fprintf (outputfile, "\n2. Bounds\n");
                    fprintf (outputfile, "Scenario    Dual Bound (root)     Best Value\n");
                    fprintf (outputfile, "------------------------------------------------\n");
                    // 2. Bounds
                    // Some people may find this data file helpful
                    // It contains the lower bounds in the root node and the upper bounds
                    // in the best node for each scenario
                    for (i = 0; i < DDSIP_param->scenarios; i++)
                        fprintf (outputfile, "%5d    %13.7g         %13.7g\n", i + 1, DDSIP_bb->btlb[i], DDSIP_bb->subsol[i]);

                    fprintf (outputfile, "\nBest Value:        %20.12g (Risk model %d)\n", DDSIP_bb->bestvalue,DDSIP_param->riskmod);
                    if (DDSIP_bb->maximization)
                        fprintf (outputfile, "Upper Bound:       %20.12g\n", DDSIP_bb->bestbound);
                    else
                        fprintf (outputfile, "Lower Bound:       %20.12g\n", DDSIP_bb->bestbound);
                    // 3. Quantiles
                    fprintf (outputfile, "\n3. Quantiles\n");
                    // Sort index
                    for (i = 0; i < DDSIP_param->scenarios; i++)
                        index[i] = i;
                    //   DDSIP_QuickSort (DDSIP_bb->subsol, index, DDSIP_param->scenarios);
                    DDSIP_qsort_ins_A (DDSIP_bb->subsol, index, 0, DDSIP_param->scenarios-1);

                    fprintf (outputfile, "Quantile      Quantile value\n");
                    fprintf (outputfile, "-----------------------------\n");

                    // Print quantiles
                    i = 1;
                    while (i <= DDSIP_param->noquant)
                    {
                        j = 0;
                        sum = 0;
                        while (sum < (double) i / DDSIP_param->noquant && j < DDSIP_param->scenarios)
                            sum += DDSIP_data->prob[index[j++]];

                        fprintf (outputfile, " %4.2f         %12.7f \n", (double) i / DDSIP_param->noquant, DDSIP_bb->subsol[index[j - 1]]);
                        i++;
                    }

                    if (DDSIP_param->outlev > 2)
                    {
                        // 4. Second Stage
                        // Print optimal second-stage solutions of each scenario
                        fprintf (outputfile, "\n4. Second-stage solutions \n");
                        for (i = 0; i < DDSIP_param->scenarios; i++)
                        {
                            fprintf (outputfile, "Scenario %d:\n", i + 1);

                            for (j = 0; j < DDSIP_bb->secvar; j++)
                                fprintf (outputfile, "%20s  %18.12g\n", colname[DDSIP_bb->secondindex[j]], DDSIP_bb->secstage[j][i]);

                        }
                    }
                    DDSIP_Free ((void **) &(colstore));
                    DDSIP_Free ((void **) &(colname));
                    DDSIP_Free ((void **) &(index));
                    fclose (outputfile);
                }
            }

        }

        if (DDSIP_param->expected)
        {
            fprintf (DDSIP_outfile, "\n----------------------------------------------");
            fprintf (DDSIP_outfile, "------------------------------------------\n");
            if (fabs (DDSIP_bb->expbest) < DDSIP_infty)
            {
                fprintf (DDSIP_outfile, "EEV       %12.7f\n", DDSIP_bb->expbest);
                if (feas && DDSIP_param->riskmod == 0)
                    fprintf (DDSIP_outfile, "VSS       %12.7f    (%12.7f %%)\n",
                             DDSIP_bb->expbest - DDSIP_bb->bestvalue, 100. * (DDSIP_bb->expbest - DDSIP_bb->bestvalue) / (fabs (DDSIP_bb->bestvalue) + 1e-6));
            }
            else
                fprintf (DDSIP_outfile, "EEV       No solution found.\n");
        }

        if (DDSIP_param->riskmod == 0)
        {
            for (i = 0; i < DDSIP_param->scenarios; i++)
                waitandsee += DDSIP_data->prob[i] * DDSIP_bb->btlb[i];
            fprintf (DDSIP_outfile, "EVPI      %12.7f    (%12.7f %%)\n",
                     DDSIP_bb->bestvalue - waitandsee, (DDSIP_bb->bestvalue - waitandsee) / (fabs (DDSIP_bb->bestvalue) + 1e-6) * 100.0);
        }
        fprintf (DDSIP_outfile, "\n----------------------------------------------");
        fprintf (DDSIP_outfile, "------------------------------------------\n");
        fprintf (DDSIP_outfile, "Status           %6d\t ", solstat);
        fprintf (DDSIP_outfile, "Time             %6.2f\n", DDSIP_GetCpuTime ());
        fprintf (DDSIP_outfile, "Tree depth       %6d\n", DDSIP_bb->depth);
        fprintf (DDSIP_outfile, "Nodes            %6d\n", DDSIP_bb->nonode);
        fprintf (DDSIP_outfile, "Cuts             %6d\n", DDSIP_bb->cutCntr);
        fprintf (DDSIP_outfile, "max. mean MIP gap LB  %g%%\n", DDSIP_bb->meanGapLB);
        fprintf (DDSIP_outfile, "max. mean MIP gap UB  %g%%\n", DDSIP_bb->meanGapUB);
        if (DDSIP_param->cb)
        {
            fprintf (DDSIP_outfile, "max. mean MIP gap CB  %g%%\n", DDSIP_bb->meanGapCBLB);
            fprintf (DDSIP_outfile, "CB function eval %6d\n", DDSIP_bb->CBIters);
        }
        fprintf (DDSIP_outfile, "Lower Bound eval %6d\n", DDSIP_bb->LBIters);
        fprintf (DDSIP_outfile, "Upper Bound eval %6d\t ", DDSIP_bb->neobjcnt);

        // Risk Model
        if (DDSIP_param->riskmod && feas)
        {
            printf ("\nExpected value     %25.10f\n", DDSIP_bb->bestexp);
            printf ("Risk measure       %25.10f\n", DDSIP_bb->bestrisk);
            fprintf (DDSIP_outfile, "Expected value     %25.10f\n", DDSIP_bb->bestexp);
            fprintf (DDSIP_outfile, "Risk measure       %25.10f\n", DDSIP_bb->bestrisk);

            fprintf (DDSIP_outfile, "\nRisk measures (for target: %g, probability level: %g):\n", DDSIP_param->risktarget, DDSIP_param->risklevel);
            printf ("\nRisk measures (for target: %g, probability level: %g):\n", DDSIP_param->risktarget, DDSIP_param->risklevel);

            fprintf (DDSIP_outfile,
                     "         ExpExc         ExcProb          SemDev              WC            TVaR             VaR          StdDev \n");
            printf
            ("         ExpExc         ExcProb          SemDev              WC            TVaR             VaR          StdDev \n");
            for (i = 0; i < DDSIP_maxrisk; i++)
            {
                fprintf (DDSIP_outfile, " %15.13g", DDSIP_bb->bestriskval[i]);
                printf (" %15.13g", DDSIP_bb->bestriskval[i]);
            }
            fprintf (DDSIP_outfile, "\n");
            printf ("\n");
        }

        fprintf (DDSIP_outfile, "\n----------------------------------------------");
        fprintf (DDSIP_outfile, "------------------------------------------\n");

        if (feas)
        {
            if (DDSIP_bb->bestvalue < DDSIP_bigvalue)
            {
                printf ("\nBest Value  : %20.12f\n", DDSIP_bb->bestvalue);
                fprintf (DDSIP_outfile, "Best Value        %20.12f\n", DDSIP_bb->bestvalue);
            }
            else
            {
                printf ("\nBest Value  : %20.12g\n", DDSIP_bb->bestvalue);
                fprintf (DDSIP_outfile, "Best Value        %20.12g\n", DDSIP_bb->bestvalue);
            }
        }
        else
            printf ("\n");

        if (bound)
        {
            if (DDSIP_bb->bestbound < DDSIP_bigvalue)
            {
                printf ("Lower Bound : %20.12f\n", DDSIP_bb->bestbound);
                fprintf (DDSIP_outfile, "Lower Bound       %20.12f\n\n", DDSIP_bb->bestbound);
            }
            else
            {
                printf ("Lower Bound : %20.12g\n", DDSIP_bb->bestbound);
                fprintf (DDSIP_outfile, "Lower Bound       %20.12g\n\n", DDSIP_bb->bestbound);
            }
            // Gaps
            if (feas)
            {
                double rgap;

                if (fabs (DDSIP_bb->bestvalue) > DDSIP_param->accuracy)
                    rgap = ((DDSIP_bb->bestvalue - DDSIP_bb->bestbound) / fabs(DDSIP_bb->bestvalue));
                else
                    rgap = ((DDSIP_bb->bestvalue - DDSIP_bb->bestbound) / (fabs(DDSIP_bb->bestvalue) + DDSIP_param->accuracy));

                rgap = DDSIP_Dmin (1, rgap);

                if (fabs (DDSIP_bb->bestvalue - DDSIP_bb->bestbound) < DDSIP_bigvalue)
                {
                    printf ("Absolute Gap: %15.8g\n", (DDSIP_bb->bestvalue - DDSIP_bb->bestbound));
                    fprintf (DDSIP_outfile, "Absolute Gap      %15.8g\n", (DDSIP_bb->bestvalue - DDSIP_bb->bestbound));
                    printf ("Relative Gap: %15.8g\n", rgap);
                    fprintf (DDSIP_outfile, "Relative Gap      %15.8g\n", rgap);
                }
                else
                {
                    printf ("Absolute Gap: %15.8g\n", (DDSIP_bb->bestvalue - DDSIP_bb->bestbound));
                    fprintf (DDSIP_outfile, "Absolute Gap      %15.8g\n", (DDSIP_bb->bestvalue - DDSIP_bb->bestbound));
                    printf ("Relative Gap:  %15.8g\n", rgap);
                    fprintf (DDSIP_outfile, "Relative Gap      %15.8g\n", rgap);
                }
            }
        }

        fprintf (DDSIP_outfile, "\n----------------------------------------------");
        fprintf (DDSIP_outfile, "------------------------------------------\n");

        // Print recourse function values to separate file
        /*  	  if (DDSIP_param->outlev) */
        /*  		{ */
        /*  		  fprintf (DDSIP_param->recfunfile,"Recourse function evaluations:\n\n"); */
        /*  		  for (i=0; i<DDSIP_bb->neobjcnt; i++) */
        /*  			{ */
        /*  			  fprintf (DDSIP_param->recfunfile,"("); */
        /*  			  for (j=0; j<DDSIP_bb->firstcon+DDSIP_bb->seccon; j++) */
        /*  				fprintf (DDSIP_param->recfunfile,"%f, ",DDSIP_bb->Txph[i*(DDSIP_bb->firstcon+DDSIP_bb->seccon)+j]); */
        /*  			  fprintf (DDSIP_param->recfunfile,"%f)",DDSIP_bb->phiofTxph[i]); */
        /*  			  fprintf (DDSIP_param->recfunfile,"\n"); */
        /*  			} */
        /*  		} */

        return 0;
    }
    else
        return 1;
}
