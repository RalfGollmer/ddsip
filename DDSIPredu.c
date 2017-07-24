/*  Authors:           Andreas M"arkert, Ralf Gollmer
	Copyright to:      University of Duisburg-Essen
    Language:          C
	Description:
	The procedure in this file builds a redundancy check problem
        and checks the added cuts for redundancy

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

void
DDSIP_CheckRedundancy ()
{
    CPXLPptr redundancy;

    int status, i, j, k, nzcnt_row, ranged = 0, beg, ind;
    static int callNr = 1;
    double *scaled_obj_coef = NULL;
    char *sense = NULL, *sense_sorted = NULL;
    char **scen_spec_rowname = NULL;
    char **scen_spec_colname = NULL;
    char **rowname = NULL, *rownamestore = NULL;
    char **colname = NULL, *colnamestore = NULL;
    int rowstorespace, rowsurplus;
    int colstorespace, colsurplus;
    char *string2 = NULL;
    double *lb = NULL, *lb_sorted = NULL;
    double *ub = NULL, *ub_sorted = NULL;
    double *rng = NULL, *rng_sorted = NULL;
    char *vartype_sorted = NULL;
    int *colindex_sorted = NULL, *colindex_revers = NULL;
    double *redundancy_rhs = NULL;
    double *base_rhs = NULL;
    int nzcnt=0, *rmatbeg=NULL, *rmatind=NULL, *rmatbeg_stage=NULL, *rmatind_stage=NULL, *rowindex=NULL;
    double *rmatval=NULL, *rmatval_stage=NULL;
    double time_start, time_end, objval;
    char probname[50];
    cutpool_t *currentCut;

    if (DDSIP_bb->cutCntr < 50*callNr)
        return;
    if (!DDSIP_bb->cutpool)
    {
        fprintf(DDSIP_outfile, "*ERROR: cut counter is %d, but cutpool is empty\n", DDSIP_bb->cutCntr);
        return;
    }

    sprintf (probname, "sipout/redundancy%d.lp.gz", callNr);

    time_start = DDSIP_GetCpuTime ();

    colstorespace = (DDSIP_data->novar + DDSIP_bb->cutCntr) * 255;
    rowstorespace = (DDSIP_data->nocon + DDSIP_bb->cutCntr) * 255;
    if (!(sense = (char *) DDSIP_Alloc (sizeof (char), DDSIP_data->nocon, "sense(CheckRedundancy)")) ||
        !(sense_sorted = (char *) DDSIP_Alloc (sizeof (char), DDSIP_Imax(DDSIP_data->firstcon,DDSIP_bb->cutCntr), "sense_sorted(CheckRedundancy)")) ||
        !(base_rhs = (double *) DDSIP_Alloc(sizeof(double),DDSIP_data->nocon,"base_rhs(CheckRedundancy)")) ||
        !(redundancy_rhs = (double *) DDSIP_Alloc(sizeof(double),DDSIP_Imax(DDSIP_data->firstcon, DDSIP_bb->cutCntr),"redundancy_rhs(CheckRedundancy)")) ||
        !(scaled_obj_coef = (double *) DDSIP_Alloc (sizeof (double), DDSIP_Imax(DDSIP_data->firstvar, DDSIP_bb->cutCntr), "scaled_obj_coef(CheckRedundancy)")) ||
        !(colname = (char **) DDSIP_Alloc (sizeof (char *), DDSIP_Imax(DDSIP_data->novar, DDSIP_bb->cutCntr),"colname(CheckRedundancy)")) ||
        !(scen_spec_colname = (char **) DDSIP_Alloc (sizeof (char *), DDSIP_Imax(DDSIP_data->firstvar,DDSIP_bb->cutCntr), "scen_spec_colname(CheckRedundancy)")) ||
        !(colnamestore = (char *) DDSIP_Alloc (sizeof (char), colstorespace, "colnamestore(CheckRedundancy)")) ||
        !(rowname = (char **) DDSIP_Alloc (sizeof (char *), DDSIP_data->nocon + DDSIP_bb->cutCntr, "rowname(DetrEqu)")) ||
        !(scen_spec_rowname = (char **) DDSIP_Alloc (sizeof (char *), DDSIP_Imax(DDSIP_data->firstcon,DDSIP_bb->cutCntr), "scen_spec_rowname(CheckRedundancy)")) ||
        !(rownamestore = (char *) DDSIP_Alloc (sizeof (char), rowstorespace, "rownamestore(CheckRedundancy)")) ||
        !(lb = (double *) DDSIP_Alloc (sizeof (double), DDSIP_data->novar, "lb(CheckRedundancy)")) ||
        !(lb_sorted = (double *) DDSIP_Alloc (sizeof (double), DDSIP_Imax(DDSIP_data->firstvar,DDSIP_bb->cutCntr), "lb_sorted(CheckRedundancy)")) ||
        !(ub = (double *) DDSIP_Alloc (sizeof (double), DDSIP_data->novar, "ub(CheckRedundancy)")) ||
        !(ub_sorted = (double *) DDSIP_Alloc (sizeof (double), DDSIP_Imax(DDSIP_data->firstvar,DDSIP_bb->cutCntr), "ub_sorted(CheckRedundancy)")) ||
        !(vartype_sorted = (char *) DDSIP_Alloc (sizeof (double), DDSIP_Imax(DDSIP_data->firstvar,DDSIP_bb->cutCntr), "vartype_sorted(CheckRedundancy)")) ||
        !(colindex_sorted = (int *) DDSIP_Alloc (sizeof (int), DDSIP_data->novar, "colindex_sorted(CheckRedundancy)")) ||
        !(rowindex = (int *) DDSIP_Alloc (sizeof (int), DDSIP_data->firstcon, "rowindex(CheckRedundancy)")))
    {
        fprintf (stderr, "Not enough memory for building redundancy check problem\n");
        goto FREE;
    }

    // get problem data
    /*____________________________________________________________________________________*/
    if((status = CPXgetcolname (DDSIP_env, DDSIP_lp, colname, colnamestore,
                            colstorespace, &colsurplus, 0, DDSIP_data->novar - 1)) ||
       (status = CPXgetrowname (DDSIP_env, DDSIP_lp, rowname, rownamestore,
                            rowstorespace, &rowsurplus, 0, DDSIP_data->nocon - 1)) ||
       (status = CPXgetsense (DDSIP_env, DDSIP_lp, sense, 0, DDSIP_data->nocon - 1)) ||
       (status = CPXgetrhs (DDSIP_env, DDSIP_lp, base_rhs, 0, DDSIP_data->nocon - 1)) ||
       (status = CPXgetlb (DDSIP_env, DDSIP_lp, lb, 0, DDSIP_data->novar - 1)) ||
       (status = CPXgetub (DDSIP_env, DDSIP_lp, ub, 0, DDSIP_data->novar - 1)))
    {
        fprintf (stderr, "Could not get problem data, returned %d\n", status);
        goto FREE;
    }
    // check whether there are ranged
    for (j=0; j<DDSIP_data->nocon; j++)
    {
        if (sense[j] == 'R')
        {
            ranged = 1;
	    break;
        }
    } 
    if (ranged)
    {
        if (!(rng = (double *) DDSIP_Alloc (sizeof (double), DDSIP_data->nocon, "rng(CheckRedundancy)")) ||
            !(rng_sorted = (double *) DDSIP_Alloc (sizeof (double), DDSIP_Imax(DDSIP_data->firstcon,DDSIP_data->seccon), "rng_sorted(CheckRedundancy)")))
        {
            fprintf (stderr, "Not enough memory for building redundancy check problem\n");
            goto FREE;
        }
        if ((status = CPXgetrngval (DDSIP_env, DDSIP_lp, rng, 0, DDSIP_data->nocon-1)))
        {
            fprintf (stderr, "Coud not get problem ranges, returned %d\n", status);
            goto FREE;
        }
    }
    /*____________________________________________________________________________________*/

    // create empty problem
    redundancy = CPXcreateprob (DDSIP_env, &status, "redundancyCheck");
    if (status)
    {
        fprintf (stderr, "CPXcreateprob returned %d\n", status);
        goto FREE;
    }

    // add (original) first-stage variables
    for (j = 0; j < DDSIP_data->firstvar; j++)
    {
        vartype_sorted[j]   = 'C';
        lb_sorted[j]        = lb[DDSIP_bb->firstindex[j]];
        ub_sorted[j]        = ub[DDSIP_bb->firstindex[j]];
        scaled_obj_coef[j]  = 0.;
        scen_spec_colname[j]= colname[DDSIP_bb->firstindex[j]];
    }
    if ((status = CPXnewcols (DDSIP_env, redundancy, DDSIP_data->firstvar, scaled_obj_coef,
                        lb_sorted, ub_sorted, vartype_sorted, scen_spec_colname)))
    {
        fprintf (stderr, "CPXnewcols returned %d for first-stage variables\n", status);
        goto FREE;
    }
    DDSIP_Free ((void **) &(colnamestore));
    // add slack variables for all cuts
    for (j = 0;  j < DDSIP_bb->cutCntr; j++)
    {
        vartype_sorted[j] = 'C';
        lb_sorted[j]      = 0.;
        ub_sorted[j]      = DDSIP_infty;
        if (!(string2 = (char *) calloc (1, 255 * sizeof (char))))
        {
            fprintf (stderr, "Not enough memory for building redundancy check model\n");
            goto FREE;
        }
        sprintf (string2, "Cut%.03d", DDSIP_bb->cutCntr-j);
        scen_spec_colname[j] = string2;
        scaled_obj_coef[j] = 0.;
    }
    if ((status = CPXnewcols (DDSIP_env, redundancy, DDSIP_bb->cutCntr, scaled_obj_coef,
                    lb_sorted, ub_sorted, vartype_sorted, scen_spec_colname)))
    {
        fprintf (stderr, "CPXnewcols returned %d for slack variables\n", status);
        goto FREE;
    }
    for (j = 0; j < DDSIP_bb->cutCntr; j++)
        DDSIP_Free ((void **) &(scen_spec_colname[j]));
        
    // free arrays needeed only for columns
    DDSIP_Free ((void **) &(colname));
    DDSIP_Free ((void **) &(lb));
    DDSIP_Free ((void **) &(ub));
    //
    // get problem matrix coefficients
    // query the length needed for storage of coefficients
    CPXgetrows(DDSIP_env, DDSIP_lp, &nzcnt, rmatbeg, rmatind, rmatval, 0, &rowsurplus, 0, DDSIP_data->nocon - 1);
    nzcnt = -rowsurplus;
    if (!(rmatbeg = (int *) DDSIP_Alloc (sizeof (int), DDSIP_data->nocon, "rmatbeg(CheckRedundancy)")) ||
        !(rmatind = (int *) DDSIP_Alloc (sizeof (int), nzcnt, "rmatind(CheckRedundancy)")) ||
        !(rmatval = (double *) DDSIP_Alloc (sizeof (double), nzcnt, "rmatval(CheckRedundancy)")))
    {
        fprintf (stderr, "Not enough memory for building redundancy check problem\n");
        goto FREE;
    }
    CPXgetrows(DDSIP_env, DDSIP_lp, &nzcnt, rmatbeg, rmatind, rmatval, nzcnt, &rowsurplus, 0, DDSIP_data->nocon - 1);

    fprintf(stderr," got %d elements of the matrix\n", nzcnt);

    beg = DDSIP_Imax(DDSIP_data->firstcon, DDSIP_bb->cutCntr);
    ind = nzcnt + DDSIP_bb->cutCntr * (DDSIP_data->firstvar + 1);
    if (!(rmatbeg_stage = (int *) DDSIP_Alloc (sizeof (int), beg, "rmatbeg_stage(CheckRedundancy)")) ||
        !(rmatind_stage = (int *) DDSIP_Alloc (sizeof (int), ind, "rmatind_stage(CheckRedundancy)")) ||
        !(rmatval_stage = (double *) DDSIP_Alloc (sizeof (double), ind, "rmatval_stage(CheckRedundancy)")))
    {
        fprintf (stderr, "Not enough memory for building redundancy check problem\n");
        goto FREE;
    }


    // add first-stage constraints
    k = 0;
    for (j = 0; j < DDSIP_data->firstcon; j++)
    {
        for (i=0; i < beg; i++)
             rmatbeg_stage[j] = 0;
        for (i=0; i < ind; i++)
        {
             rmatind_stage[j] = 0;
             rmatval[j] = 0.;
        }
        sense_sorted[j] = sense[DDSIP_bb->firstrowind[j]];
        redundancy_rhs[j] = base_rhs[DDSIP_bb->firstrowind[j]];
        scen_spec_rowname[j] = rowname[DDSIP_bb->firstrowind[j]];
	rmatbeg_stage[j] = k;
	if (DDSIP_bb->firstrowind[j] == DDSIP_data->nocon -1)
            nzcnt_row = nzcnt - rmatbeg[DDSIP_data->nocon -1];
        else
            nzcnt_row = rmatbeg[DDSIP_bb->firstrowind[j]+1] - rmatbeg[DDSIP_bb->firstrowind[j]];
	for (i = 0; i < nzcnt_row; i++)
        {
            rmatind_stage[k + i] = DDSIP_bb->firstindex_reverse[rmatind[rmatbeg[DDSIP_bb->firstrowind[j]] + i]];
            rmatval_stage[k + i] = rmatval[rmatbeg[DDSIP_bb->firstrowind[j]] + i];
        }
	k += nzcnt_row;
    }
    if ((status = CPXaddrows(DDSIP_env, redundancy, 0, DDSIP_data->firstcon, k, redundancy_rhs, sense_sorted, rmatbeg_stage, rmatind_stage, rmatval_stage, NULL, scen_spec_rowname)))
    {
        fprintf (stderr, "CPXaddrows returned %d for first-stage constraints\n", status);
        goto FREE;
    }
    if (ranged)
    {
        for (j = 0; j < DDSIP_data->firstcon; j++)
        {
            rng_sorted[j] = rng[DDSIP_bb->firstrowind[j]];
            rowindex[j]   = j;
        }
        if((status = CPXchgrngval(DDSIP_env, redundancy, DDSIP_data->firstcon, rowindex, rng_sorted)))
        {
            fprintf (stderr, "CPXchgrngval returned %d for first-stage constraints\n", status);
            goto FREE;
        }
    }

    // add the cuts
    k = 0;
    j = 0;
    currentCut = DDSIP_bb->cutpool;
    while (currentCut)
    {
        for (i=0; i < beg; i++)
             rmatbeg_stage[j] = 0;
        for (i=0; i < ind; i++)
        {
             rmatind_stage[j] = 0;
             rmatval[j] = 0.;
        }
        sense_sorted[j] = 'E';
        redundancy_rhs[j] = currentCut->rhs;
        if (!(string2 = (char *) calloc (1, 255 * sizeof (char))))
        {
            fprintf (stderr, "Not enough memory for building redundancy check problem\n");
            goto FREE;
        }
        sprintf (string2, "DDSIPCut%.03d", currentCut->number);

        scen_spec_rowname[j] = string2;
        rmatbeg_stage[j] = k;
        nzcnt_row = DDSIP_data->novar + 1;
        for (i = 0; i < DDSIP_data->firstvar; i++)
        {
            rmatind_stage[k + i] = i;
            rmatval_stage[k + i] = currentCut->matval[i];
        }
        rmatind_stage[k + DDSIP_data->firstvar] = DDSIP_data->firstvar + j;
        rmatval_stage[k + DDSIP_data->firstvar] = -1.;
        k += DDSIP_bb->firstvar + 1;
        currentCut = currentCut->prev;
        j++;
    }
    if ((status = CPXaddrows(DDSIP_env, redundancy, 0, j, k, redundancy_rhs, sense_sorted, rmatbeg_stage, rmatind_stage, rmatval_stage, NULL, scen_spec_rowname)))
    {
        fprintf (stderr, "CPXaddrows returned %d for cut constraints\n", status);
        goto FREE;
    }
    if (j < DDSIP_bb->cutCntr)
    {
        fprintf (stderr, "  ------- while loop found only %d < %d cuts\n", j, DDSIP_bb->cutCntr);
    }
    for (i = 0; i < j; i++)
        DDSIP_Free ((void **) &(scen_spec_rowname[i]));

    time_end = DDSIP_GetCpuTime ();
    fprintf (DDSIP_outfile, " %6.2f sec  for building redundancy check problem\n",time_end-time_start);
    // now test for redundant cuts
    sense[0] = 'L';
    lb_sorted[0]       = -2.;
    ub_sorted[0]       = 0.;
    ind = 0;
    for (k = 0; k < j; k++)
    {
        // change lower bound of one slack variable and put it in the objective
        rmatind_stage[0] = DDSIP_data->firstvar + k;
        if ((status = CPXchgbds(DDSIP_env, redundancy, 1, rmatind_stage, sense, lb_sorted)))
        {
            fprintf (stderr, "CPXchgbds returned %d\n", status);
            goto FREE;
        }
        if ((status = CPXchgcoef(DDSIP_env, redundancy, -1, rmatind_stage[0], 1.)))
        {
            fprintf (stderr, "CPXchcoef returned %d\n", status);
            goto FREE;
        }
        // optimize
        status = CPXchgprobtype (DDSIP_env, redundancy, CPXPROB_LP);
        status = CPXdualopt (DDSIP_env, redundancy);
        if (status)
        {
             fprintf (stderr, "CPXdualopt retuned %d for k= %d\n", status, k);
             fprintf (DDSIP_outfile, "ERROR*: CPXdualopt returned %d for k= %d\n", status, k);
             if (DDSIP_param->outlev)
                 fprintf (DDSIP_bb->moreoutfile, "ERROR*: CPXdualopt status %d for k= %d\n", status, k);
             sprintf (probname, "sipout/redundancy%d_%d.lp.gz", callNr,k);
             status = CPXwriteprob (DDSIP_env, redundancy, probname, NULL);
             if (status)
             {
                 fprintf (DDSIP_outfile, " *** redundancy check problem not written successfully, status = %d\n", status);
                 printf  (" *** redundancy check problem not written successfully, status = %d\n", status);
             }
             else
             {
                 fprintf (DDSIP_outfile, " *** redundancy check problem %s written successfully\n", probname);
                 printf  (" *** redundancy check problem %s written successfully\n", probname);
             }
            continue;
        }
        status = CPXgetstat (DDSIP_env, redundancy);
        if (status != CPX_STAT_OPTIMAL)
        {
             fprintf (stderr, "CPXdualopt status %d not optimal\n", status);
             fprintf (DDSIP_outfile, "ERROR*: CPXdualopt status %d not optimal\n", status);
             if (DDSIP_param->outlev)
                 fprintf (DDSIP_bb->moreoutfile, "ERROR*: CPXdualopt status %d not optimal\n", status);
            continue;
        }
        status = CPXgetobjval (DDSIP_env, redundancy, &objval);
        if (status)
        {
             fprintf (stderr, "ERROR*: Failed to get objective value \n");
             fprintf (DDSIP_outfile, "ERROR*: Failed to get objective value \n");
             if (DDSIP_param->outlev)
                 fprintf (DDSIP_bb->moreoutfile, "ERROR*: Failed to get objective value \n");
             continue;
        }
        // if the objval is nonnegative, the cut is redundant
        if (objval >= 0.)
        {
             fprintf (stderr, "######## Cut %d is redundant\n", DDSIP_bb->cutCntr-k);
             fprintf (DDSIP_outfile, "######## Cut %d is redundant\n", DDSIP_bb->cutCntr-k);
             if (DDSIP_param->outlev)
                 fprintf (DDSIP_bb->moreoutfile, "######## Cut %d is redundant\n", DDSIP_bb->cutCntr-k);
             ind++;
        }
        
        if ((status = CPXchgbds(DDSIP_env, redundancy, 1, rmatind_stage, sense, ub_sorted)))
        {
            fprintf (stderr, "CPXchgbds returned %d\n", status);
            goto FREE;
        }
        if ((status = CPXchgcoef(DDSIP_env, redundancy, -1, rmatind_stage[0], 0.)))
        {
            fprintf (stderr, "CPXchcoef returned %d\n", status);
            goto FREE;
        }
    }
    if (!ind)
    {
         fprintf (stderr, "++++++++ none of the %d cuts is redundant\n", DDSIP_bb->cutCntr);
         fprintf (DDSIP_outfile, "++++++++ none of the %d cuts is redundant\n", DDSIP_bb->cutCntr);
         if (DDSIP_param->outlev)
             fprintf (DDSIP_bb->moreoutfile, "++++++++ none of the %d cuts is redundant\n", DDSIP_bb->cutCntr);
    }
    else
    {
         fprintf (stderr, "++++++++ %d of the %d cuts are redundant\n", ind, DDSIP_bb->cutCntr);
         fprintf (DDSIP_outfile, "++++++++ %d of the %d cuts are redundant\n", ind, DDSIP_bb->cutCntr);
         if (DDSIP_param->outlev)
             fprintf (DDSIP_bb->moreoutfile, "++++++++ %d of the %d cuts are redundant\n", ind, DDSIP_bb->cutCntr);
    }

    time_start = DDSIP_GetCpuTime ();
    fprintf (DDSIP_outfile, " %6.2f sec  for checking redundancy of cuts\n",time_start-time_end);

    status = CPXfreeprob (DDSIP_env, &redundancy);

    callNr++;

FREE:
    DDSIP_Free ((void **) &(vartype_sorted));
    DDSIP_Free ((void **) &(ub_sorted));
    DDSIP_Free ((void **) &(lb_sorted));
    DDSIP_Free ((void **) &(scaled_obj_coef));
    DDSIP_Free ((void **) &(sense));
    DDSIP_Free ((void **) &(sense_sorted));
    DDSIP_Free ((void **) &(rowname));
    DDSIP_Free ((void **) &(rownamestore));
    DDSIP_Free ((void **) &(colname));
    DDSIP_Free ((void **) &(colnamestore));
    DDSIP_Free ((void **) &(redundancy_rhs));
    DDSIP_Free ((void **) &(base_rhs));
    DDSIP_Free ((void **) &(lb));
    DDSIP_Free ((void **) &(ub));
    DDSIP_Free ((void **) &(vartype_sorted));
    DDSIP_Free ((void **) &(lb_sorted));
    DDSIP_Free ((void **) &(ub_sorted));
    DDSIP_Free ((void **) &(scaled_obj_coef));
    DDSIP_Free ((void **) &(colindex_sorted));
    DDSIP_Free ((void **) &(colindex_revers));
    DDSIP_Free ((void **) &(scen_spec_rowname));
    DDSIP_Free ((void **) &(scen_spec_colname));
    DDSIP_Free ((void **) &(rmatbeg));
    DDSIP_Free ((void **) &(rmatind));
    DDSIP_Free ((void **) &(rmatval));
    DDSIP_Free ((void **) &(rmatbeg_stage));
    DDSIP_Free ((void **) &(rmatind_stage));
    DDSIP_Free ((void **) &(rmatval_stage));
    DDSIP_Free ((void **) &(rowindex));
    if (ranged)
    {
        DDSIP_Free ((void **) &(rng));
        DDSIP_Free ((void **) &(rng_sorted));
    }
    return;
}
