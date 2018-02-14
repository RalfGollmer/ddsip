/*  Authors:           Andreas M"arkert, Ralf Gollmer
	Copyright to:      University of Duisburg-Essen
    Language:          C
	Description:
	Initializations: Set cplex parameters, bbtype, advanced solution, sort scenarios

	License:
	This file is part of DDSIP.

    DDSIP is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the DDSIP_Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    DDSIP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with DDSIP; if not, write to the DDSIP_Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

#include <DDSIP.h>
#include <DDSIPconst.h>
#include <math.h>
#include <limits.h>

int DDSIP_SortScen (void);

//==========================================================================
// Function sets cplex parameters, different parameters could be used in LB, UB etc.
int
DDSIP_SetCpxPara (const int cnt, const int * isdbl, const int * which, const double * what)
{
    int i;
    int status = 0;

    status = CPXsetdefaults (DDSIP_env);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to set default parameters\n");
        return status;
    }
    // if not changed by user, we require minimal feasibilty, optimality and integrality tolerances
    if ((status = CPXsetdblparam (DDSIP_env, CPX_PARAM_EPRHS, 1e-8)))
    {
        fprintf (stderr, "ERROR: Failed to set minimal feasibility tolerance\n");
        return status;
    }
    if ((status = CPXsetdblparam (DDSIP_env, CPX_PARAM_EPOPT, 1e-8)))
    {
        fprintf (stderr, "ERROR: Failed to set minimal optimality tolerance\n");
        return status;
    }
    if ((status = CPXsetdblparam (DDSIP_env, CPX_PARAM_EPINT, 0.0)))
    {
        fprintf (stderr, "ERROR: Failed to set minimal integrality tolerance\n");
        return status;
    }

    if (which != DDSIP_param->cpxwhich)
    {
        for (i = 0; i < DDSIP_param->cpxno; i++)
        {
            if (DDSIP_param->cpxisdbl[i] == 1)
                status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxwhich[i], DDSIP_param->cpxwhat[i]);
            else if (DDSIP_param->cpxisdbl[i] == 0)
                status = CPXsetintparam (DDSIP_env, DDSIP_param->cpxwhich[i], floor (DDSIP_param->cpxwhat[i] + 0.1));
            else if (DDSIP_param->cpxisdbl[i] == 3)
                status = CPXsetlongparam (DDSIP_env, DDSIP_param->cpxwhich[i], floor (DDSIP_param->cpxwhat[i] + 0.1));
            else
            {
                fprintf (stderr,"Error: unexpected parameter type for parameter (general) %d\n", DDSIP_param->cpxwhich[i]);
                status = 1;
            }

            if (status)
            {
                fprintf (stderr, "ERROR: Failed to set cplex parameter %d (general).\n", DDSIP_param->cpxwhich[i]);
                return status;
            }
        }
    }
    for (i = 0; i < cnt; i++)
    {
        if (isdbl[i] == 1)
            status = CPXsetdblparam (DDSIP_env, which[i], what[i]);
        else if (isdbl[i] == 0)
            status = CPXsetintparam (DDSIP_env, which[i], floor (what[i] + 0.1));
        else if (isdbl[i] == 3)
            status = CPXsetlongparam (DDSIP_env, which[i], floor (what[i] + 0.1));
        else
        {
            fprintf (stderr,"Error: unexpected parameter type for parameter (section) %d\n", which[i]);
            status = 1;
        }

        if (status)
        {
            fprintf (stderr, "ERROR: Failed to set cplex parameter %d (section).\n", which[i]);
            return status;
        }
    }
    if (DDSIP_param->watchkappa)
    {
        if ((status = CPXsetintparam (DDSIP_env, CPX_PARAM_MIPKAPPASTATS, DDSIP_param->watchkappa)))
        {
            fprintf (stderr, "ERROR: Failed to set kappastats\n");
            return status;
        }
    }
    return 0;
} // DDSIP_SetCpxPara

//==========================================================================
// Produce some output
int
DDSIP_CpxParaPrint (void)
{
    int status = 0, i;

    fprintf (DDSIP_outfile, "-----------------------------------------------------------\n");
    fprintf (DDSIP_outfile, "-CPLEX PARAMETERS (see CPLEX manual)\n\n");
    fprintf (DDSIP_outfile, "  Number                Value   \n");

    for (i = 0; i < DDSIP_param->cpxno; i++)
        if (DDSIP_param->cpxisdbl[i] > 0)
            fprintf (DDSIP_outfile, "  %4d           %12.4g\n", DDSIP_param->cpxwhich[i], DDSIP_param->cpxwhat[i]);
        else if (DDSIP_param->cpxisdbl[i] == 0)
            fprintf (DDSIP_outfile, "  %4d           %12.0f\n", DDSIP_param->cpxwhich[i], floor (DDSIP_param->cpxwhat[i] + 0.1));
        else
            fprintf (DDSIP_outfile, "  %4d   CPLEX complained!\n", DDSIP_param->cpxwhich[i]);

    if (DDSIP_param->cpxnolb)
        fprintf (DDSIP_outfile, "CLEXLB -   Parameters for evaluation of lower bounds \n");
    for (i = 0; i < DDSIP_param->cpxnolb; i++)
        if (DDSIP_param->cpxlbisdbl[i] > 0)
            fprintf (DDSIP_outfile, "  %4d           %12.4g\n", DDSIP_param->cpxlbwhich[i], DDSIP_param->cpxlbwhat[i]);
        else if (DDSIP_param->cpxlbisdbl[i] == 0)
            fprintf (DDSIP_outfile, "  %4d           %12.0f\n", DDSIP_param->cpxlbwhich[i], floor (DDSIP_param->cpxlbwhat[i] + 0.1));
        else
            fprintf (DDSIP_outfile, "  %4d   CPLEX complained!\n", DDSIP_param->cpxlbwhich[i]);
    if (DDSIP_param->cpxnolb2)
        fprintf (DDSIP_outfile, "CPLEX2LB - Parameters for second opt. in lower bounds \n");
    for (i = 0; i < DDSIP_param->cpxnolb2; i++)
        if (DDSIP_param->cpxlbisdbl2[i] > 0)
            fprintf (DDSIP_outfile, "  %4d           %12.4g\n", DDSIP_param->cpxlbwhich2[i], DDSIP_param->cpxlbwhat2[i]);
        else if (DDSIP_param->cpxlbisdbl2[i] == 0)
            fprintf (DDSIP_outfile, "  %4d           %12.0f\n", DDSIP_param->cpxlbwhich2[i], floor (DDSIP_param->cpxlbwhat2[i] + 0.1));
        else
            fprintf (DDSIP_outfile, "  %4d   CPLEX complained!\n", DDSIP_param->cpxlbwhich2[i]);

    if (DDSIP_param->cpxnoub)
        fprintf (DDSIP_outfile, "CPLEXUB -  Parameters for evaluation of upper bounds \n");
    for (i = 0; i < DDSIP_param->cpxnoub; i++)
        if (DDSIP_param->cpxubisdbl[i] > 0)
            fprintf (DDSIP_outfile, "  %4d           %12.4g\n", DDSIP_param->cpxubwhich[i], DDSIP_param->cpxubwhat[i]);
        else if (DDSIP_param->cpxubisdbl[i] == 0)
            fprintf (DDSIP_outfile, "  %4d           %12.0f\n", DDSIP_param->cpxubwhich[i], floor (DDSIP_param->cpxubwhat[i] + 0.1));
        else
            fprintf (DDSIP_outfile, "  %4d   CPLEX complained!\n", DDSIP_param->cpxubwhich[i]);
    if (DDSIP_param->cpxnoub2)
        fprintf (DDSIP_outfile, "CLPEX2UB - Parameters for second opt. in upper bounds \n");
    for (i = 0; i < DDSIP_param->cpxnoub2; i++)
        if (DDSIP_param->cpxubisdbl2[i] > 0)
            fprintf (DDSIP_outfile, "  %4d           %12.4g\n", DDSIP_param->cpxubwhich2[i], DDSIP_param->cpxubwhat2[i]);
        else if (DDSIP_param->cpxubisdbl2[i] == 0)
            fprintf (DDSIP_outfile, "  %4d           %12.0f\n", DDSIP_param->cpxubwhich2[i], floor (DDSIP_param->cpxubwhat2[i] + 0.1));
        else
            fprintf (DDSIP_outfile, "  %4d   CPLEX complained!\n", DDSIP_param->cpxubwhich2[i]);

    if (DDSIP_param->cpxnoeev)
        fprintf (DDSIP_outfile, "CPLEXEEV - Parameters for EV problem\n");
    for (i = 0; i < DDSIP_param->cpxnoeev; i++)
        if (DDSIP_param->cpxeevisdbl[i] > 0)
            fprintf (DDSIP_outfile, "  %4d           %12.4e\n", DDSIP_param->cpxeevwhich[i], DDSIP_param->cpxeevwhat[i]);
        else if (DDSIP_param->cpxeevisdbl[i] == 0)
            fprintf (DDSIP_outfile, "  %4d           %12.0f\n", DDSIP_param->cpxeevwhich[i], floor (DDSIP_param->cpxeevwhat[i]) + 0.1);
        else
            fprintf (DDSIP_outfile, "  %4d   CPLEX complained!\n", DDSIP_param->cpxeevwhich[i]);

    if (DDSIP_param->cpxnodual)
        fprintf (DDSIP_outfile, "CPLEXDUAL  Parameters for Lagrangian dual \n");
    for (i = 0; i < DDSIP_param->cpxnodual; i++)
        if (DDSIP_param->cpxdualisdbl[i] > 0)
            fprintf (DDSIP_outfile, "  %4d           %12.4e\n", DDSIP_param->cpxdualwhich[i], DDSIP_param->cpxdualwhat[i]);
        else if (DDSIP_param->cpxdualisdbl[i] == 0)
            fprintf (DDSIP_outfile, "  %4d           %12.0f\n", DDSIP_param->cpxdualwhich[i], floor (DDSIP_param->cpxdualwhat[i] + 0.1));
        else
            fprintf (DDSIP_outfile, "  %4d   CPLEX complained!\n", DDSIP_param->cpxdualwhich[i]);

    if (DDSIP_param->cpxnodual2)
        fprintf (DDSIP_outfile, "CPLEX2DUAL Parameters for second call in Lagrangian dual \n");
    for (i = 0; i < DDSIP_param->cpxnodual2; i++)
        if (DDSIP_param->cpxdualisdbl2[i] > 0)
            fprintf (DDSIP_outfile, "  %4d           %12.4e\n", DDSIP_param->cpxdualwhich2[i], DDSIP_param->cpxdualwhat2[i]);
        else if (DDSIP_param->cpxdualisdbl2[i] == 0)
            fprintf (DDSIP_outfile, "  %4d           %12.0f\n", DDSIP_param->cpxdualwhich2[i], floor (DDSIP_param->cpxdualwhat2[i] + 0.1));
        else
            fprintf (DDSIP_outfile, "  %4d   CPLEX complained!\n", DDSIP_param->cpxdualwhich2[i]);

    fprintf (DDSIP_outfile, "END of PARAMETER SECTION\n");
    fprintf (DDSIP_outfile, "-----------------------------------------------------------\n\n");

    return status;
} // DDSIP_CpxParaPrint




//==========================================================================
// Initialization of some cplex parameters, detecting integer/double parameters
// Initial setting of parameters
int
DDSIP_InitCpxPara (void)
{
    int i;
    int status = 0;
    double tmp;

    // Some parameters need an initialization
    // Screen indicator
    status = CPXinfointparam (DDSIP_env, CPX_PARAM_SCRIND, &i, NULL, NULL);
    DDSIP_param->cpxscr = i;
    DDSIP_param->cpxubscr = i;
    // Time limit
    status = CPXinfodblparam (DDSIP_env, CPX_PARAM_TILIM, &tmp, NULL, NULL);
    DDSIP_param->cpxtime = tmp;
    DDSIP_param->cpxubtime = tmp;
    // Relative gap
    status = CPXinfodblparam (DDSIP_env, CPX_PARAM_EPGAP, &tmp, NULL, NULL);
    DDSIP_param->cpxgap = tmp;
    DDSIP_param->cpxubgap = tmp;
    // Node limit
    status = CPXinfointparam (DDSIP_env, CPX_PARAM_NODELIM, &i, NULL, NULL);
    DDSIP_param->cpxnodelim = i;
    DDSIP_param->cpxubnodelim = i;

    // save some parameters seperately
    for (i = 0; i < DDSIP_param->cpxno; i++)
    {
        // Some parameters are saved seperately
        if (DDSIP_param->cpxwhich[i] == CPX_PARAM_TILIM)
        {
            DDSIP_param->cpxtime = DDSIP_param->cpxwhat[i];
            DDSIP_param->cpxubtime = DDSIP_param->cpxwhat[i];
        }
        else if (DDSIP_param->cpxwhich[i] == CPX_PARAM_SCRIND)
        {
            DDSIP_param->cpxscr = DDSIP_param->cpxwhat[i];
            DDSIP_param->cpxubscr = DDSIP_param->cpxwhat[i];
        }
        else if (DDSIP_param->cpxwhich[i] == CPX_PARAM_EPGAP)
        {
            DDSIP_param->cpxgap = DDSIP_param->cpxwhat[i];
            DDSIP_param->cpxubgap = DDSIP_param->cpxwhat[i];
            if (DDSIP_param->cpxgap > DDSIP_param->relgap)
            {
                fprintf (stderr, "*Warning: CPLEX GAP (%2.2e) greater than DUAL GAP (%2.2e).\n", DDSIP_param->cpxgap, DDSIP_param->relgap);
                fprintf (stderr, "*         This may lead to redundant branching.\n");
            }
        }
        else if (DDSIP_param->cpxwhich[i] == CPX_PARAM_NODELIM)
        {
            DDSIP_param->cpxnodelim = DDSIP_param->cpxwhat[i];
            DDSIP_param->cpxubnodelim = DDSIP_param->cpxwhat[i];
        }
    }

    // Print output to file
    DDSIP_CpxParaPrint ();

    status = DDSIP_SetCpxPara (DDSIP_param->cpxno, DDSIP_param->cpxisdbl, DDSIP_param->cpxwhich, DDSIP_param->cpxwhat);
    if (status)
        fprintf (stderr, "ERROR: Failed to set initial cplex parameters (SetCpxInitPara).\n");

    return status;
} // DDSIP_InitCpxPara

//==========================================================================
// The function sorts the scenario according to a certain criteria
int
DDSIP_SortScen (void)
{
    int i, j, chg;

    double tmp;

    double *sum = (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios,
                                          "sum(SortScen)");

    if (DDSIP_param->outlev > 3)
        fprintf (DDSIP_bb->moreoutfile, "Sorting scenarios according to sum of entries.\n");

    // Sum of entries
    if (DDSIP_param->outlev > 7)
        fprintf (DDSIP_bb->moreoutfile, "\tOriginal Sums\n");

    if (DDSIP_param->stocrhs)
    {
        for (i = 0; i < DDSIP_param->scenarios; i++)
        {
            sum[i] = DDSIP_data->rhs[i * DDSIP_param->stocrhs];
            for (j = 1; j < DDSIP_param->stocrhs; j++)
                sum[i] += DDSIP_data->rhs[i * DDSIP_param->stocrhs + j];
            if (DDSIP_param->outlev > 7)
                fprintf (DDSIP_bb->moreoutfile, "\tSum(%d)=%2.2f\n", i, sum[i]);
        }

        // Bubble sort
        chg = 1;
        while (chg)
        {
            chg = 0;
            for (i = 0; i < DDSIP_param->scenarios - 1; i++)
                if (sum[i] > sum[i + 1])
                {
                    tmp = sum[i];
                    sum[i] = sum[i + 1];
                    sum[i + 1] = tmp;

                    for (j = 0; j < DDSIP_param->stocrhs; j++)
                    {
                        tmp = DDSIP_data->rhs[i * DDSIP_param->stocrhs + j];
                        DDSIP_data->rhs[i * DDSIP_param->stocrhs + j] = DDSIP_data->rhs[ (i + 1) * DDSIP_param->stocrhs + j];
                        DDSIP_data->rhs[ (i + 1) * DDSIP_param->stocrhs + j] = tmp;
                    }
                    chg = 1;
                }
        }
    }

    DDSIP_Free ((void **) &(sum));
    return 0;
} // DDSIP_SortScen

//==========================================================================
// Initialization of bbtype members
int
DDSIP_BbTypeInit (void)
{
    int i, j, status = 0;

    time (&DDSIP_bb->start_time);
    // New = Original
    DDSIP_bb->firstvar = DDSIP_data->firstvar;
    DDSIP_bb->secvar = DDSIP_data->secvar;
    DDSIP_bb->firstcon = DDSIP_data->firstcon;
    DDSIP_bb->seccon = DDSIP_data->seccon;
    DDSIP_bb->correct_bounding = 0.;
    DDSIP_bb->LBIters = 0;
    DDSIP_bb->CBIters = 0;
    DDSIP_bb->from_scenario = -1;

    // Change according to risk model
    if (!DDSIP_param->riskalg && !DDSIP_param->scalarization)
        switch (abs (DDSIP_param->riskmod))
        {
            // ExpExc
        case 1:
            DDSIP_bb->secvar++;
            DDSIP_bb->seccon++;
            break;
            // ExcProb
        case 2:
            DDSIP_bb->secvar++;
            DDSIP_bb->seccon++;
            break;
            // SemDev
        case 3:
            DDSIP_bb->secvar++;
            DDSIP_bb->seccon++;
            break;
            // WorstCase
        case 4:
            DDSIP_bb->firstvar++;
            DDSIP_bb->seccon++;
            break;
            // TVaR
        case 5:
            DDSIP_bb->firstvar++;
            DDSIP_bb->secvar++;
            DDSIP_bb->seccon++;
            break;
        }

    DDSIP_bb->novar = DDSIP_bb->firstvar + DDSIP_bb->secvar;
    DDSIP_bb->nocon = DDSIP_bb->firstcon + DDSIP_bb->seccon;
    /////
    if(DDSIP_param->riskmod)
    {
        printf ("with risk model\n");
        printf ("\t\t No. of              variables:    %6d\n", DDSIP_bb->novar);
        printf ("\t\t No. of  first-stage variables:    %6d  (%d generals, %d binary, %d continuous)\n", DDSIP_bb->firstvar, DDSIP_bb->first_int - DDSIP_bb->first_bin, DDSIP_bb->first_bin, DDSIP_bb->firstvar - DDSIP_bb->first_int);
        printf ("\t\t No. of second-stage variables:    %6d  (%d integers)\n", DDSIP_bb->secvar, DDSIP_bb->total_int - DDSIP_bb->first_int);
        fprintf (DDSIP_outfile, "-----------------------------------------------------------\n");
        fprintf (DDSIP_outfile, "with risk model\n");
        fprintf (DDSIP_outfile, "\t\t No. of              variables:    %6d\n", DDSIP_bb->novar);
        fprintf (DDSIP_outfile, "\t\t No. of  first-stage variables:    %6d  (%d generals, %d binary, %d continuous)\n", DDSIP_bb->firstvar, DDSIP_bb->first_int - DDSIP_bb->first_bin, DDSIP_bb->first_bin, DDSIP_bb->firstvar - DDSIP_bb->first_int);
        fprintf (DDSIP_outfile, "\t\t No. of second-stage variables:    %6d  (%d integers)\n", DDSIP_bb->secvar, DDSIP_bb->total_int - DDSIP_bb->first_int);
        fprintf (DDSIP_outfile, "-----------------------------------------------------------\n");
    }

#ifdef ADDINTEGERCUTS
    if (DDSIP_bb->firstvar > DDSIP_bb->first_bin)
        DDSIP_param->addIntegerCuts = 0;
#endif

    // Memory allocation for bb-type-members
    DDSIP_bb->solstat = (int *) DDSIP_Alloc (sizeof (int), DDSIP_param->scenarios, "solstat(BbTypeInit)");
    DDSIP_bb->firstindex = (int *) DDSIP_Alloc (sizeof (int), DDSIP_bb->firstvar, "firstindex(BbTypeInit)");
    DDSIP_bb->firstindex_reverse = (int *) DDSIP_Alloc (sizeof (int), DDSIP_bb->novar, "firstindex(BbTypeInit)");
    DDSIP_bb->secondindex = (int *) DDSIP_Alloc (sizeof (int), DDSIP_bb->secvar, "secondindex(BbTypeInit)");
    DDSIP_bb->secondindex_reverse = (int *) DDSIP_Alloc (sizeof (int), DDSIP_bb->novar, "secondindex(BbTypeInit)");
    if (DDSIP_param->order)
        DDSIP_bb->order = (int *) DDSIP_Alloc (sizeof (int), DDSIP_bb->firstvar, "order(BbTypeInit)");
    DDSIP_bb->firsttype = (char *) DDSIP_Alloc (sizeof (char), DDSIP_bb->firstvar, "firsttype(BbTypeInit)");
    DDSIP_bb->sectype = (char *) DDSIP_Alloc (sizeof (char), DDSIP_bb->secvar, "sectype(BbTypeInit)");

    DDSIP_bb->lbident = (char *) DDSIP_Alloc (sizeof (char), DDSIP_bb->firstvar, "lbident(BbTypeInit)");
    DDSIP_bb->ubident = (char *) DDSIP_Alloc (sizeof (char), DDSIP_bb->firstvar, "ubident(BbTypeInit)");
    j = ceil (0.5 * DDSIP_param->nodelim) + 3;
    DDSIP_bb->front = (int *) DDSIP_Alloc (sizeof (int), j, "front(BbTypeInit)");

    DDSIP_bb->bestsol = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar, "bestsol(BbTypeInit)");
    DDSIP_bb->lborg = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar, "lborg(BbTypeInit)");
    DDSIP_bb->uborg = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar, "uborg(BbTypeInit)");

    DDSIP_bb->objcontrib = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->novar, "objcontrib(BbTypeInit)");
    DDSIP_bb->secstage = (double **) DDSIP_Alloc (sizeof (double *), DDSIP_bb->secvar, "secstage(BbTypeInit)");
    for (i = 0; i < DDSIP_bb->secvar; i++)
        DDSIP_bb->secstage[i] = (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "secstage[i](BbTypeInit)");

    DDSIP_bb->subsol = (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "subsol(BbTypeInit)");

    /*    DDSIP_bb->Txph = (double *) DDSIP_Alloc(sizeof(double), DDSIP_param->nodelim * (DDSIP_bb->firstcon + DDSIP_bb->seccon),"Txph(BbTypeInit)"); */
    /*    DDSIP_bb->phiofTxph = (double *) DDSIP_Alloc(sizeof(double), DDSIP_param->nodelim,"phiofTxph(BbTypeInit)"); */

    // Lower bounds in root node and upper bounds in best node are needed
    // to calculate best target
    DDSIP_bb->btlb = (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "btlb(BbTypeInit)");

    DDSIP_bb->curind = (int *) DDSIP_Alloc (sizeof (int), DDSIP_bb->firstvar, "curind(BbTypeInit)");
    DDSIP_bb->curlb = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar, "curlb(BbTypeInit)");
    DDSIP_bb->curub = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar, "curub(BbTypeInit)");

    DDSIP_bb->sug = (sug_t **) DDSIP_Alloc (sizeof (sug_t *), DDSIP_param->nodelim + 3, "sug(BbTypeInit)");
    for (i = 0; i < DDSIP_param->nodelim + 2; i++)
        DDSIP_bb->sug[i] = NULL;
    //For Conic Bundle: space to save best scenario solutions
    if (DDSIP_param->cb)
        DDSIP_bb->bestfirst = (bestfirst_t *) DDSIP_Alloc (sizeof (bestfirst_t), DDSIP_param->scenarios, "bestfirst(BbTypeInit)");

    if (DDSIP_param->riskmod)
    {
        DDSIP_bb->bestriskval = (double *) malloc (sizeof (double) * DDSIP_maxrisk);
        DDSIP_bb->curriskval = (double *) malloc (sizeof (double) * DDSIP_maxrisk);
    }
    if ((DDSIP_param->scalarization || DDSIP_param->cb))
        DDSIP_bb->ref_risk = (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "ref_risk(BbTypeInit)");
    DDSIP_bb->ref_scenobj = (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "ref_scenobj(BbTypeInit)");
    DDSIP_bb->n_buffer_len = 256;
    DDSIP_bb->n_buffer = (char *) DDSIP_Alloc (sizeof (char), DDSIP_bb->n_buffer_len, "n_buffer(Main)");
    // Initial settings
    DDSIP_bb->heurval = DDSIP_infty;
    DDSIP_bb->bestvalue = DDSIP_infty;
    DDSIP_bb->expbest = DDSIP_infty;
    DDSIP_bb->bestbound = -DDSIP_infty;
    DDSIP_bb->feasbound = -DDSIP_infty;
    DDSIP_bb->skip = 0;
    DDSIP_bb->nonode = 1;
    DDSIP_bb->nofront = 1;
    DDSIP_bb->no_reduced_front = 1;
    DDSIP_bb->curnode = 0;
    DDSIP_bb->curbdcnt = 0;

    DDSIP_bb->lboutcnt = 1;
    DDSIP_bb->uboutcnt = 1;
    // The  depth of the b&b tree
    DDSIP_bb->depth = 1;
    DDSIP_bb->front[0] = 0;

    DDSIP_bb->cutoff = 0;

    // Count the number of evaluated upper bounds
    DDSIP_bb->found_optimal_node = 0;
    DDSIP_bb->bound_optimal_node = DDSIP_infty;

    DDSIP_bb->neobjcnt = 0;

#ifdef CONIC_BUNDLE
    if (DDSIP_param->cb)
    {
        if (DDSIP_param->nonant == 3)
            DDSIP_bb->dimdual = DDSIP_bb->firstvar * DDSIP_param->scenarios;
        else
            DDSIP_bb->dimdual = DDSIP_bb->firstvar * (DDSIP_param->scenarios - 1);
    }
#endif

    DDSIP_node[0]->father = -1;
    DDSIP_node[0]->depth  = 0;
    DDSIP_node[0]->neoind = 0;
    DDSIP_node[0]->solved = 0;

    DDSIP_node[0]->bound = -DDSIP_infty;
    DDSIP_node[0]->dispnorm = DDSIP_infty;
    DDSIP_node[0]->leaf = 0;

    if (DDSIP_param->riskmod == 3)
        DDSIP_node[0]->target = -DDSIP_infty;
    //node[0]->target = DDSIP_param->risktarget;

    if (DDSIP_param->cb)
    {
        DDSIP_node[0]->dual = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->dimdual + 3, "dual(BbTypeInit)");
        DDSIP_node[0]->dual[DDSIP_bb->dimdual] = DDSIP_param->cbweight;
        DDSIP_node[0]->dual[DDSIP_bb->dimdual + 2] = -DDSIP_infty;
    }

    for (i = 0; i < DDSIP_bb->firstvar; DDSIP_bb->bestsol[i++] = 0.0);

    for (i = 0; i < DDSIP_param->scenarios; DDSIP_bb->btlb[i++] = -DDSIP_infty);

    // Initializing lower and upper bound identifiers used in cplex-routines
    //ASCII L=76, U=85
    for (i = 0; i < DDSIP_bb->firstvar; i++)
    {
        DDSIP_bb->lbident[i] = 76;
        DDSIP_bb->ubident[i] = 85;
    }
    // DDSIP_Allocate and initialize lb_scen_order for lower bound, scen_order for upper bound
    DDSIP_bb->lb_scen_order = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_param->scenarios), "DDSIP_bb->lb_scen_order(BbInit)");
    DDSIP_bb->ub_scen_order = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_param->scenarios), "DDSIP_bb->ub_scen_order(BbInit)");
    for(i=0; i<DDSIP_param->scenarios; i++)
    {
        DDSIP_bb->lb_scen_order[i] = DDSIP_bb->ub_scen_order[i]    = i;
    }
    DDSIP_bb->lb_sorted = 0;
    DDSIP_bb->ub_sorted = 0;
    DDSIP_bb->front_nodes_sorted = (int *) DDSIP_Alloc (sizeof (int), 1, "DDSIP_bb->front_nodes_sorted(BbInit)");
    DDSIP_bb->front_nodes_sorted[0] = 0;
    DDSIP_bb->meanGapLB = DDSIP_bb->meanGapCBLB = DDSIP_bb->meanGapUB = 0.;
    DDSIP_bb->bestBound = DDSIP_bb->newTry = DDSIP_bb->cutCntr = DDSIP_bb->cutAdded = 0;
    DDSIP_bb->bestsol_in_curnode = 1;
    if (DDSIP_param->cb)
    {
        DDSIP_bb->bestdual = NULL;
        DDSIP_bb->bestdual_cnt = 0;
        DDSIP_bb->bestdual_max = -DDSIP_infty;
        DDSIP_bb->local_bestdual = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->dimdual + 3, "bestdual(BbTypeInit)");
    }
    else
    {
        DDSIP_bb->bestdual = NULL;
        DDSIP_bb->local_bestdual = NULL;
    }
    DDSIP_bb->shifts = 0;

    return status;
} // DDSIP_BbTypeInit

//==========================================================================
// Arrange branching order for master bb (may include additional first-stage variables)
int
DDSIP_BranchOrder (void)
{
    int i;

    // Copy provided data, decrease by 1
    for (i = 0; i < DDSIP_data->firstvar; i++)
        DDSIP_bb->order[i] = DDSIP_data->order[i] - 1;
    // The original order information is not needed any more
    DDSIP_Free ((void **)&(DDSIP_data->order));
    // Some risk models introduce additional first-stage variables. Their branching order corresponds
    // to the parameter brancheta
    if (DDSIP_bb->firstvar - DDSIP_data->firstvar)
        DDSIP_bb->order[DDSIP_data->firstvar] = DDSIP_param->brancheta;

    return 0;
} // DDSIP_BranchOrder

//==========================================================================
// Initialize variables related to stages
int
DDSIP_InitStages (void)
{
    // Temporary variables
    int status = 0, cnt, i, j, length, ind;
    int seccnt;

    int *firindex = (int *) DDSIP_Alloc (sizeof (int), DDSIP_data->novar+1, "firindex(InitStages)");
    int *secondindex = (int *) DDSIP_Alloc (sizeof (int), DDSIP_data->novar+1, "secondindex(InitStages)");

    char *ctype = (char *) DDSIP_Alloc (sizeof (char), DDSIP_data->novar+2, "ctype(InitStages)");
    char **colname = (char **) DDSIP_Alloc (sizeof (char *), DDSIP_data->novar, "colname(InitStages)");
    char *colstore = (char *) DDSIP_Alloc (sizeof (char), DDSIP_data->novar * DDSIP_ln_varname,
                                           "colstore(InitStages)");
    double *lb, *ub;

    //Prepare first- and second-stage variables
    // Get all variable names
    status = CPXgetcolname (DDSIP_env, DDSIP_lp, colname, colstore, DDSIP_bb->novar * DDSIP_ln_varname, &i, 0, DDSIP_bb->novar - 1);
    if(i < 0)
    {
        fprintf (stderr, "ERROR: variable names require more storage than expected by DDSIP, missing extra storage: %d\n",-i);
        fprintf (stderr, "       reserved space: %d * %d = %d\n", DDSIP_bb->novar, DDSIP_ln_varname, DDSIP_bb->novar * DDSIP_ln_varname);
        cnt = 0;
        for (seccnt=0; seccnt<DDSIP_bb->novar; seccnt++)
        {
            status = CPXgetcolname (DDSIP_env, DDSIP_lp, colname, colstore, DDSIP_bb->novar * DDSIP_ln_varname, &i, seccnt, seccnt);
            cnt = DDSIP_Dmax(cnt, strlen(colname[0]));
        }
        fprintf (stderr, "       maximal variable name length: %d,  reserved for each variable name: %d chars.\n", cnt, DDSIP_ln_varname-1);
        return(cnt);
    }
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to get column names\n");
        return status;
    }
    // Get all types
    status = CPXgetctype (DDSIP_env, DDSIP_lp, ctype, 0, DDSIP_bb->novar - 1);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to get variables \n");
        return status;
    }
    // Which variables are in the first stage ?
    // --> defined by PREFIX or POSTFIX in the specs file
    cnt = 0;
    seccnt = 0;
    i = 0;
    DDSIP_bb->total_int = 0;
    DDSIP_bb->first_int  = 0;
    DDSIP_bb->first_bin  = 0;

    if (DDSIP_param->prefix)
    {
        if (!(length = strlen(DDSIP_param->prefix)))
        {
            fprintf (stderr," *** ERROR: The prefix for the first stage variables has to have a positive length.\n");
            exit (1);
        }
        while (i < DDSIP_bb->novar)
        {
            if (!strncmp(DDSIP_param->prefix,colname[i],length))
            {
                firindex[cnt] = i;
                if ((ctype[i] == 'B') || (ctype[i] == 'I'))
                {
                    DDSIP_bb->total_int++;
                    DDSIP_bb->first_int++;
                    if (ctype[i] == 'B')
                        DDSIP_bb->first_bin++;
                }
                cnt++;
                // Debug output of names of first-stage variables
                if (DDSIP_param->outlev > 30)
                {
                    printf("first-stage variable %3d:  %c  %s\n", cnt, ctype[i], colname[i]);
                }
            }
            else
            {
                if ((ctype[i] == 'B') || (ctype[i] == 'I'))
                {
                    DDSIP_bb->total_int++;
                }
                secondindex[seccnt++] = i;
            }
            i++;
        }
    }
    else
    {
        if (!(length = strlen(DDSIP_param->postfix)))
        {
            fprintf (stderr," *** ERROR: The postfix for the first stage variables has to have a positive length.\n");
            exit (1);
        }
        while (i < DDSIP_bb->novar)
        {
            j = strlen (colname[i])-length;
            if (j > 0 && !strcmp(DDSIP_param->postfix,colname[i]+j))
            {
                firindex[cnt] = i;
                if ((ctype[i] == 'B') || (ctype[i] == 'I'))
                {
                    DDSIP_bb->total_int++;
                    DDSIP_bb->first_int++;
                    if (ctype[i] == 'B')
                        DDSIP_bb->first_bin++;
                }
                cnt++;
                // Debug output of names of first-stage variables
                if (DDSIP_param->outlev > 30)
                {
                    printf("first-stage variable %3d:  %c  %s\n", cnt, ctype[i], colname[i]);
                }
            }
            else
            {
                if ((ctype[i] == 'B') || (ctype[i] == 'I'))
                {
                    DDSIP_bb->total_int++;
                }
                secondindex[seccnt++] = i;
            }
            i++;
        }
    }
    DDSIP_data->firstvar = cnt;
    DDSIP_data->secvar   = DDSIP_data->novar - cnt;
    if (!cnt)
    {
        printf ("XXX ERROR: no first-stage variables found. Wrong post- or prefix given? Problem cannot be handled by DDSIP.\n    Exiting.\n");
        fprintf (DDSIP_outfile, "XXX ERROR: no first-stage variables found. Wrong post- or prefix given? Problem cannot be handled by DDSIP.\n    Exiting.\n");
        return -1;
    }
    if (abs(DDSIP_param->riskmod) == 4)
    {
        firindex[cnt]             = DDSIP_data->novar;
        ctype[DDSIP_data->novar]  = 'C';
    }
    else if (abs(DDSIP_param->riskmod) == 5)
    {
        firindex[cnt]             = DDSIP_data->novar;
        secondindex[seccnt]       = DDSIP_data->novar+1;
        ctype[DDSIP_data->novar]  = 'C';
        ctype[DDSIP_data->novar+1]= 'C';
    }
    else if (DDSIP_param->riskmod)
    {
        if (abs(DDSIP_param->riskmod) == 1)
        {
            secondindex[seccnt]      = DDSIP_data->novar;
            ctype[DDSIP_data->novar] = 'B';
        }
        else
        {
            secondindex[seccnt]      = DDSIP_data->novar;
            ctype[DDSIP_data->novar] = 'C';
        }
    }

    // Initialize b&b variables
    if ((status = DDSIP_BbInit ()))
    {
        fprintf (stderr, "ERROR: Failed to initialize Bb data (InitStages), return code = %d\n", status);
        return status;
    }
    printf ("initialized Bb data (InitStages)\n");

    for (i = 0; i < DDSIP_bb->firstvar; i++)
    {
        DDSIP_bb->firstindex[i] = firindex[i];
        DDSIP_bb->firstindex_reverse[firindex[i]]  =  i;
        DDSIP_bb->secondindex_reverse[firindex[i]] = -1;
        DDSIP_bb->firsttype[i] = ctype[DDSIP_bb->firstindex[i]];
    }

    for (i = 0; i < DDSIP_bb->secvar; i++)
    {
        DDSIP_bb->secondindex[i] = secondindex[i];
        DDSIP_bb->secondindex_reverse[secondindex[i]] =  i;
        DDSIP_bb->firstindex_reverse[secondindex[i]]  = -1;
        DDSIP_bb->sectype[i] = ctype[DDSIP_bb->secondindex[i]];
    }

    DDSIP_Free ((void **) &(secondindex));
    DDSIP_Free ((void **) &(firindex));
    DDSIP_Free ((void **) &(ctype));


    DDSIP_bb->cost = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar, "DDSIP_bb->cost,InitStages");

    for (i = 0; i < DDSIP_data->firstvar; i++)
    {
        DDSIP_bb->cost[i] = DDSIP_data->cost[DDSIP_bb->firstindex[i]+ DDSIP_param->stoccost * DDSIP_param->scenarios];
    }
    if (abs(DDSIP_param->riskmod) == 4 || abs(DDSIP_param->riskmod) == 5)
    {
        if (DDSIP_param->riskmod > 0)
            DDSIP_bb->cost[DDSIP_data->firstvar] = DDSIP_param->riskweight;
        else
            DDSIP_bb->cost[DDSIP_data->firstvar] = 1.;
    }

    // Preserve lower and upper bounds on first stage variables
    lb = (double *) DDSIP_Alloc (sizeof (double), (DDSIP_bb->novar), "lb(InitStages)");
    ub = (double *) DDSIP_Alloc (sizeof (double), (DDSIP_bb->novar), "ub(InitStages)");
    status = CPXgetlb (DDSIP_env, DDSIP_lp, lb, 0, DDSIP_bb->novar - 1);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to get lower bounds \n");
        return status;
    }

    status = CPXgetub (DDSIP_env, DDSIP_lp, ub, 0, DDSIP_bb->novar - 1);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to get upper bounds \n");
        return status;
    }

    ind = 0;
    for (i = 0; i < DDSIP_bb->firstvar; i++)
    {
        DDSIP_bb->lborg[i] = lb[DDSIP_bb->firstindex[i]];
        DDSIP_bb->uborg[i] = ub[DDSIP_bb->firstindex[i]];

        if (DDSIP_param->cb)
        {
            if (!(DDSIP_bb->lborg[i] > -DDSIP_infty) || !(DDSIP_bb->uborg[i] < DDSIP_infty))
            {
                if (DDSIP_param->outlev > 2)
                    fprintf (DDSIP_bb->moreoutfile, "!!!! unbounded variable %25s [%g,%g]\n", colname[DDSIP_bb->firstindex[i]], DDSIP_bb->lborg[i], DDSIP_bb->uborg[i]);
                ind++;
            }
        }
    }

    DDSIP_Free ((void **) &(colstore));
    DDSIP_Free ((void **) &(colname));

    // Cplex does not always tell us whether a problem is infeasible or unbounded.
    // We treat both as infeasible. Unboundedness should be removed by imposing
    // upper bounds on all variables
    if (ind)
    {
        printf ("*Warning: %d first-stage variable(s) unbounded. This may cause problems with dual method (unboundedness of scenario problems due to Lagrangean term).\n", ind);
        fprintf (DDSIP_outfile, "*Warning: %d first-stage variable(s) unbounded. This may cause problems with dual method (unboundedness of scenario problems due to Lagrangean term).\n", ind);
        if (DDSIP_param->outlev)
            fprintf (DDSIP_bb->moreoutfile, "*Warning: %d first-stage variable(s) unbounded. This may cause problems with dual method (unboundedness of scenario problems due to Lagrangean term).\n", ind);
    }

    DDSIP_Free ((void **) &(lb));
    DDSIP_Free ((void **) &(ub));

    return status;
} // DDSIP_InitStages

//==========================================================================
// Identify constraints related to stages
int
DDSIP_DetectStageRows (void)
{
    // Now find out, which constraints are in the first stage
    int i, j, k, status = 0, nonzeros, surplus;
    int *RowSecondStage = (int *) DDSIP_Alloc (sizeof (int), DDSIP_data->nocon+1, "RowSecondStage(InitStages)");
    int *rmatbeg = (int *) DDSIP_Alloc (sizeof (int), 1, "rmatbeg(InitStages)");
    int *rmatind = (int *) DDSIP_Alloc (sizeof (int), DDSIP_data->novar, "rmatind(InitStages)");
    double *rmatval = (double *) DDSIP_Alloc (sizeof (double), DDSIP_data->novar+2, "rmatval(InitStages)");
    double time_start, time_end;
    time_start = DDSIP_GetCpuTime ();
    // first, all rows containing stochastic rhs or matrix entries are second stage
    for (i = 0; i < DDSIP_param->stocrhs; i++)
    {
        RowSecondStage[DDSIP_data->rhsind[i]] = 1;
    }
    for (i = 0; i < DDSIP_param->stocmat; i++)
    {
        RowSecondStage[DDSIP_data->matrow[i]] = 1;
    }
    // for the rest of the rows check whether there is a second-stage variable involved
    for (i = 0; i < DDSIP_data->nocon; i++)
    {
        if (!(RowSecondStage[i]))
        {
            // get row data from model
            if ((status = CPXgetrows (DDSIP_env, DDSIP_lp, &nonzeros, rmatbeg, rmatind, rmatval, DDSIP_data->novar, &surplus, i, i)))
                return status;
            for (j = 0; j < nonzeros; j++)
            {
                for (k = 0; k < DDSIP_data->firstvar; k++)
                {
                    if (rmatind[j] == DDSIP_bb->firstindex[k])
                        break;
                }
                if (k == DDSIP_data->firstvar)
                {
                    RowSecondStage[i] = 1;
                    break;
                }
            }
        }
    }
    for (i = 0; i < DDSIP_data->nocon; i++)
    {
        if (RowSecondStage[i])
        {
            DDSIP_data->seccon++;
        }
        else
        {
            DDSIP_data->firstcon++;
        }
    }

    DDSIP_bb->firstcon = DDSIP_data->firstcon;
    DDSIP_bb->seccon = DDSIP_data->seccon;

    // Change according to risk model
    if (DDSIP_param->riskmod)
    {
        DDSIP_bb->seccon++;
        RowSecondStage[DDSIP_data->nocon] = 1;
    }
    DDSIP_bb->nocon = DDSIP_bb->firstcon + DDSIP_bb->seccon;

    DDSIP_bb->firstrowind = (int *) DDSIP_Alloc(sizeof (int), DDSIP_bb->firstcon, "firstrowind(DDSIP_DetectStageRows)");
    DDSIP_bb->secondrowind = (int *) DDSIP_Alloc(sizeof (int), DDSIP_bb->seccon, "secondrowind(DDSIP_DetectStageRows)");
    DDSIP_bb->firstrowind_reverse = (int *) DDSIP_Alloc(sizeof (int), DDSIP_bb->nocon, "firstrowind_reverse(DDSIP_DetectStageRows)");
    DDSIP_bb->secondrowind_reverse = (int *) DDSIP_Alloc(sizeof (int), DDSIP_bb->nocon, "secondrowind_reverse(DDSIP_DetectStageRows)");
    j = k = 0;
    for (i = 0; i < DDSIP_bb->nocon; i++)
    {
        if (RowSecondStage[i])
        {
            DDSIP_bb->secondrowind[j] = i;
            DDSIP_bb->secondrowind_reverse[i] =  j;
            DDSIP_bb->firstrowind_reverse[i]  = -1;
            j++;
        }
        else
        {
            DDSIP_bb->firstrowind[k] = i;
            DDSIP_bb->firstrowind_reverse[i]  =  k;
            DDSIP_bb->secondrowind_reverse[i] = -1;
            k++;
        }
    }
    if (DDSIP_param->outlev > 30)
    {
        printf("First-stage constraint indices:");
        for (i = 0; i < DDSIP_data->firstcon; i++)
        {
            if (!(i%10))
                printf ("\n");
            printf(" %6d,",DDSIP_bb->firstrowind[i]);
        }
        printf ("\n");
    }
    time_end = DDSIP_GetCpuTime ();

    fprintf (DDSIP_outfile, "%6.2f sec  for detection of stages for constraints\n------------------------------------------------------------\n",time_end-time_start);
    fprintf (DDSIP_outfile, "in input model file:\n");
    fprintf (DDSIP_outfile, "\t\tNo. of              constraints:  %10d\n", DDSIP_data->nocon);
    fprintf (DDSIP_outfile, "\t\tNo. of  first-stage constraints:  %10d\n", DDSIP_data->firstcon);
    fprintf (DDSIP_outfile, "\t\tNo. of second-stage constraints:  %10d\n", DDSIP_data->seccon);
    if (DDSIP_param->riskmod)
    {
        fprintf (DDSIP_outfile, "with risk model: \n");
        fprintf (DDSIP_outfile, "\t\tNo. of              constraints:  %10d\n", DDSIP_bb->nocon);
        fprintf (DDSIP_outfile, "\t\tNo. of second-stage constraints:  %10d\n", DDSIP_bb->seccon);
    }
    fprintf (DDSIP_outfile, "------------------------------------------------------------\n");
    printf ("in input model file:\n");
    printf ("\t\tNo. of              constraints:  %10d\n", DDSIP_data->nocon);
    printf ("\t\tNo. of  first-stage constraints:  %10d\n", DDSIP_data->firstcon);
    printf ("\t\tNo. of second-stage constraints:  %10d\n", DDSIP_data->seccon);
    if (DDSIP_param->riskmod)
    {
        printf ("with risk model: \n");
        printf ("\t\tNo. of              constraints:  %10d\n", DDSIP_bb->nocon);
        printf ("\t\tNo. of second-stage constraints:  %10d\n", DDSIP_bb->seccon);
    }
    DDSIP_Free ((void **) &(RowSecondStage));
    DDSIP_Free ((void **) &(rmatbeg));
    DDSIP_Free ((void **) &(rmatind));
    DDSIP_Free ((void **) &(rmatval));
    return 0;   
} // DDSIP_DetectStageRows 

//==========================================================================
// For use of advanced solutions we need some initializations and an additional constraint
int
DDSIP_AdvSolInit (void)
{
    int i, j, status = 0;

    char *ctype = (char *) DDSIP_Alloc (sizeof (char), DDSIP_bb->novar, "ctype(AdvSolInit)");

    DDSIP_bb->objbndind = INT_MIN;
    // Get all types
    status = CPXgetctype (DDSIP_env, DDSIP_lp, ctype, 0, DDSIP_bb->novar - 1);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to get variable types \n");
        return status;
    }
    // Fill DDSIP_bb->intind with indices of all integer variables
    DDSIP_bb->intind = (int *) DDSIP_Alloc (sizeof (int), DDSIP_bb->total_int + 1, "intind(AdvSolInit)");
    j = 0;
    for (i = 0; i < DDSIP_bb->novar; i++)
        if ((ctype[i] == 'B') || (ctype[i] == 'I') || (ctype[i] == 'N'))
            DDSIP_bb->intind[j++] = i;
    if (DDSIP_param->cb)
    {
        DDSIP_bb->intsolvals = (double **) DDSIP_Alloc (sizeof (double *), DDSIP_param->scenarios, "intsolvals(AdvSolInit)");
        for (j = 0; j < DDSIP_param->scenarios; j++)
            DDSIP_bb->intsolvals[j] = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->total_int, "intsolvals[j](CBLowerBound)");
        DDSIP_bb->boundIncrease_intsolvals = (double **) DDSIP_Alloc (sizeof (double *), DDSIP_param->scenarios, "boundIncrease_intsolvals(AdvSolInit)");
        for (j = 0; j < DDSIP_param->scenarios; j++)
            DDSIP_bb->boundIncrease_intsolvals[j] = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->total_int, "boundIncrease_intsolvals[j](CBLowerBound)");
    }

    DDSIP_Free ((void **) &(ctype));

    // DDSIP_Allocate memory for solutions of node 0
    if (DDSIP_param->hot)
    {
        DDSIP_node[0]->solut = (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios * DDSIP_bb->total_int, "solut(AdvSolInit)");
        // DDSIP_Allocate memory for MIP start values
        if (DDSIP_param->hot != 2)
        {
            DDSIP_bb->values  = (double *)DDSIP_Alloc(sizeof(double), 2*DDSIP_param->scenarios*DDSIP_bb->total_int, "values(AdvSolInit)");
            DDSIP_bb->Names   = (char **) DDSIP_Alloc(sizeof(char *), 2*DDSIP_param->scenarios, "Names(AdvSolInit)");
            //DDSIP_bb->beg     = (int *)   DDSIP_Alloc(sizeof(int *),  2*DDSIP_param->scenarios, "beg(AdvSolInit)");
            DDSIP_bb->beg     = (int *)   DDSIP_Alloc(sizeof(int),    1, "beg(AdvSolInit)");
            DDSIP_bb->effort  = (int *)   DDSIP_Alloc(sizeof(int *),  1, "effort(AdvSolInit)");
            DDSIP_bb->effort[0] = 2;
            DDSIP_bb->beg[0] = 0;
            for(j=0; j<2*DDSIP_param->scenarios; j++)
            {
                DDSIP_bb->Names[j] = (char *) DDSIP_Alloc (sizeof (char), 24, "Names(LowerBound)");
            }
        }
    }
    if (DDSIP_param->hot == 2 || DDSIP_param->hot > 4)
    {
        int *rowlist = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_bb->firstvar + DDSIP_bb->secvar),
                                            "rowlist (AdvSolInit)");
        int *collist = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_bb->firstvar + DDSIP_bb->secvar),
                                            "collist (AdvSolInit)");
        double *vallist = (double *) DDSIP_Alloc (sizeof (double), (DDSIP_bb->firstvar + DDSIP_bb->secvar),
                          "vallist(AdvSolInit)");
        double *rhs = (double *) DDSIP_Alloc (sizeof (double), 1, "rhs(AdvSolInit)");
        char *sense = (char *) DDSIP_Alloc (sizeof (char), 1, "sense(AdvSolInit)");


        // Handle lower bound on objective function in scenario problems
        sense[0] = 'G';
        rhs[0] = -DDSIP_infty;

        // Add constraint to problem
        status = CPXnewrows (DDSIP_env, DDSIP_lp, 1, rhs, sense, NULL, NULL);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to add new constraint (BbInit) \n");
            return status;
        }
        // This results in one additional second stage constraint
        DDSIP_bb->seccon++;
        DDSIP_bb->nocon++;

        // Remember index of constraint
        DDSIP_bb->objbndind = DDSIP_bb->nocon - 1;

        // Change coefficients of new constraint
        for (i = 0; i < DDSIP_bb->novar; i++)
        {
            rowlist[i] = DDSIP_bb->nocon - 1;
            collist[i] = i;
        }

        status = CPXgetobj (DDSIP_env, DDSIP_lp, vallist, 0, DDSIP_bb->novar - 1);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to get objective function coefficients (AdvSolInit) \n");
            return status;
        }

        status = CPXchgcoeflist (DDSIP_env, DDSIP_lp, DDSIP_bb->novar, rowlist, collist, vallist);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change coefficients (AdvSolInit) \n");
            return status;
        }

        DDSIP_Free ((void **) &(sense));
        DDSIP_Free ((void **) &(rhs));
        DDSIP_Free ((void **) &(vallist));
        DDSIP_Free ((void **) &(collist));
        DDSIP_Free ((void **) &(rowlist));
    }

    return status;
} // DDSIP_AdvSolInit

//==========================================================================
// Function initializes b&b tree and detects indicies, lower bounds,
// upper bounds and types of first stage variables
// Scenario sorting, branching order, advanced solution, risk modelling
int
DDSIP_BbInit (void)
{
    // Temporary variables
    int status = 0, i;
    double *obj_coef = (double *) DDSIP_Alloc (sizeof (double), DDSIP_data->novar,"DDSIP_data->obj_coef(BbInit)");
    if (!obj_coef)
    {
        fprintf (stderr, "ERROR: Failed to allocate memory for objective coefficients\n");
        return status;
    }
    DDSIP_data->obj_coef = obj_coef;

    printf ("Initializing branch-and-bound tree.\n");

    if (DDSIP_param->outlev)
        fprintf (DDSIP_bb->moreoutfile, "Initializing branch-and-bound tree.\n");

    // Print infos
    printf ("\t\t Data from input model file:\n");
    printf ("\t\t No. of              variables:     %6d\n", DDSIP_data->novar);
    printf ("\t\t No. of  first-stage variables:     %6d  (%d generals, %d binary, %d continuous)\n", DDSIP_data->firstvar, DDSIP_bb->first_int - DDSIP_bb->first_bin, DDSIP_bb->first_bin, DDSIP_data->firstvar - DDSIP_bb->first_int);
    printf ("\t\t No. of second-stage variables:     %6d  (%d integers)\n", DDSIP_data->secvar, DDSIP_bb->total_int - DDSIP_bb->first_int);
    fprintf (DDSIP_outfile, "------------------------------------------------------------\n");
    fprintf (DDSIP_outfile, "in input model file:\n");
    fprintf (DDSIP_outfile, "\t\tNo. of              variables:    %10d\n", DDSIP_data->novar);
    fprintf (DDSIP_outfile, "\t\tNo. of  first-stage variables:    %10d  (%d generals, %d binary, %d continuous)\n", DDSIP_data->firstvar, DDSIP_bb->first_int - DDSIP_bb->first_bin, DDSIP_bb->first_bin, DDSIP_data->firstvar - DDSIP_bb->first_int);
    fprintf (DDSIP_outfile, "\t\tNo. of second-stage variables:    %10d  (%d integers)\n", DDSIP_data->secvar, DDSIP_bb->total_int - DDSIP_bb->first_int);
    fprintf (DDSIP_outfile, "------------------------------------------------------------\n");

#ifdef ADDINTEGERCUTS
    if (DDSIP_bb->first_int > DDSIP_bb->first_bin)
        DDSIP_param->addIntegerCuts = 0;
#endif

    status = CPXgetobj (DDSIP_env, DDSIP_lp, DDSIP_data->obj_coef, 0, DDSIP_data->novar - 1);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to get objective coefficients\n");
        return status;
    }

    DDSIP_bb->adv_sol = NULL;
    if (CPXgetobjsen (DDSIP_env, DDSIP_lp) == CPX_MAX)
    {
        DDSIP_bb->novar = DDSIP_data->novar;
        int *index = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_bb->novar), "index(BbInit)");

        DDSIP_bb->maximization = 1;
        printf ("Maximization becomes minimization (objective=-1*objective).\n");

        for (i = 0; i < DDSIP_bb->novar; i++)
        {
            if (DDSIP_data->obj_coef[i] != 0.0)
            {
                DDSIP_data->obj_coef[i] = -DDSIP_data->obj_coef[i];
            }
            index[i] = i;
        }

        status = CPXchgobj (DDSIP_env, DDSIP_lp, DDSIP_bb->novar, index, DDSIP_data->obj_coef);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change obj coefficients\n");
            return status;
        }
        // Change problem to minimization
        CPXchgobjsen (DDSIP_env, DDSIP_lp, CPX_MIN);
        for (i = 0; i < DDSIP_data->novar + DDSIP_param->stoccost * DDSIP_param->scenarios; i++)
            if (DDSIP_data->cost[i] != 0.0)
                DDSIP_data->cost[i] = -DDSIP_data->cost[i];

        DDSIP_Free ((void **) &(index));
    }
    else
        DDSIP_bb->maximization= 0;

    // Add a risk model
    if (DDSIP_param->riskmod)
    {
        status = DDSIP_RiskModel ();
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to add risk model\n");
            return status;
        }
    }
    // Initialize bbtype members
    status = DDSIP_BbTypeInit ();
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to initialize bb-type\n");
        return status;
    }
    // Sort scenarios
    if (DDSIP_param->prepro)
    {
        status = DDSIP_SortScen ();
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to sort scenarios (BbInit)\n");
            return status;
        }
    }

    // Initialisation, no meaning
    (DDSIP_node[0])->neolb = DDSIP_bb->lborg[0];
    (DDSIP_node[0])->neoub = DDSIP_bb->uborg[0];

    // Modification for advanced solution usage
    if (DDSIP_param->hot)
    {
        status = DDSIP_AdvSolInit ();
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to initialize advanced solution constraints\n");
            return status;
        }
    }

    return status;
} // DDSIP_BbInit
