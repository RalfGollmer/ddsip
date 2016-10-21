/*  Authors:           Andreas M"arkert, Ralf Gollmer
	Copyright to:      University of Duisburg-Essen
    Language:          C
    Last modification: 06.02.2016
	Description:
	The procedure in this file builds a deterministic equivalent and writes
        the lp.gz file via CPLEX calls -- This is implemented only for
        expectation-based models.

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

void
DDSIP_DetEqu ()
{
    CPXLPptr det_equ;

    int status, scen, i, j, k, nzcnt_row, ranged = 0;
    char probname[] = "sipout/det_equ.lp.gz";
    double *obj_coef = NULL;
    double *scaled_obj_coef = NULL;
    char *sense = NULL, *sense_sorted = NULL;
    char **scen_spec_rowname = NULL;
    char **scen_spec_colname = NULL;
    char **rowname = NULL, *rownamestore = NULL;
    char **colname = NULL, *colnamestore = NULL;
    int rowstorespace, rowsurplus;
    int colstorespace, colsurplus;
    char *string1 = NULL, *string2 = NULL;
    double *lb = NULL, *lb_sorted = NULL;
    double *ub = NULL, *ub_sorted = NULL;
    double *rng = NULL, *rng_sorted = NULL;
    char *vartype = NULL, *vartype_sorted = NULL;
    int *colindex_sorted = NULL, *colindex_revers = NULL;
    double *value = NULL;
    double *det_equ_rhs = NULL;
    double *base_rhs = NULL;
    int nzcnt=0, *rmatbeg=NULL, *rmatind=NULL, *rmatbeg_stage=NULL, *rmatind_stage=NULL, *rowindex=NULL;
    double *rmatval=NULL, *rmatval_stage=NULL;
    if (DDSIP_data->seccon)
        det_equ_rhs = (double *) DDSIP_Alloc(sizeof(double),DDSIP_data->seccon,"det_equ_rhs(DetEqu)");
    else
    {
        fprintf (stderr,"XXX ERROR: no second stage contraints, got DDSIP_data->seccon=%d.\n",DDSIP_data->seccon);
        return;
    }

    fprintf (stderr,
             "\nBuilding deterministic equivalent. This may take some time.\nWorks only for expectation model so far (no risk measure).\n");
             //"\nBuilding deterministic equivalent. This may take some time.\nWorks only for expectation-based model so far.\n");

    colstorespace = DDSIP_data->novar * 255;
    rowstorespace = DDSIP_data->nocon * 255;
    if (!(sense = (char *) DDSIP_Alloc (sizeof (char), DDSIP_data->nocon, "sense(DetEqu)")) ||
        !(sense_sorted = (char *) DDSIP_Alloc (sizeof (char), DDSIP_Imax(DDSIP_data->firstcon,DDSIP_data->seccon), "sense_sorted(DetEqu)")) ||
        !(base_rhs = (double *) DDSIP_Alloc(sizeof(double),DDSIP_data->nocon,"base_rhs(DetEqu)")) ||
        !(obj_coef = (double *) DDSIP_Alloc (sizeof (double), DDSIP_data->novar,"base_rhs(DetEqu)")) ||
        !(scaled_obj_coef = (double *) DDSIP_Alloc (sizeof (double), DDSIP_Imax(DDSIP_data->firstvar, DDSIP_data->secvar), "base_rhs(DetEqu)")) ||
        !(colname = (char **) DDSIP_Alloc (sizeof (char *), DDSIP_data->novar,"base_rhs(DetEqu)")) ||
        !(scen_spec_colname = (char **) DDSIP_Alloc (sizeof (char *), DDSIP_Imax(DDSIP_data->firstvar,DDSIP_data->secvar), "scen_spec_colname(DetEqu)")) ||
        !(colnamestore = (char *) DDSIP_Alloc (sizeof (char), colstorespace, "colnamestore(DetEqu)")) ||
        !(rowname = (char **) DDSIP_Alloc (sizeof (char *), DDSIP_data->nocon, "rowname(DetrEqu)")) ||
        !(scen_spec_rowname = (char **) DDSIP_Alloc (sizeof (char *), DDSIP_Imax(DDSIP_data->firstcon,DDSIP_data->seccon), "scen_spec_rowname(DetEqu)")) ||
        !(rownamestore = (char *) DDSIP_Alloc (sizeof (char), rowstorespace, "rownamestore(DetEqu)")) ||
        !(lb = (double *) DDSIP_Alloc (sizeof (double), DDSIP_data->novar, "lb(DetEqu)")) ||
        !(lb_sorted = (double *) DDSIP_Alloc (sizeof (double), DDSIP_Imax(DDSIP_data->firstvar,DDSIP_data->secvar), "lb_sorted(DetEqu)")) ||
        !(ub = (double *) DDSIP_Alloc (sizeof (double), DDSIP_data->novar, "ub(DetEqu)")) ||
        !(ub_sorted = (double *) DDSIP_Alloc (sizeof (double), DDSIP_Imax(DDSIP_data->firstvar,DDSIP_data->secvar), "ub_sorted(DetEqu)")) ||
        !(vartype = (char *) DDSIP_Alloc (sizeof (char), DDSIP_data->novar, "vartype(DetEqu)")) ||
        !(vartype_sorted = (char *) DDSIP_Alloc (sizeof (double), DDSIP_Imax(DDSIP_data->firstvar,DDSIP_data->secvar), "vartype_sorted(DetEqu)")) ||
        !(colindex_sorted = (int *) DDSIP_Alloc (sizeof (int), DDSIP_data->novar, "colindex_sorted(DetEqu)")) ||
        !(rowindex = (int *) DDSIP_Alloc (sizeof (int), DDSIP_Imax(DDSIP_data->firstcon, DDSIP_data->seccon), "rowindex(DetEqu)")))
    {
        fprintf (stderr, "Not enough memory for building deterministic equivalent\n");
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
       (status = CPXgetobj (DDSIP_env, DDSIP_lp, obj_coef, 0, DDSIP_data->novar - 1)) ||
       (status = CPXgetlb (DDSIP_env, DDSIP_lp, lb, 0, DDSIP_data->novar - 1)) ||
       (status = CPXgetub (DDSIP_env, DDSIP_lp, ub, 0, DDSIP_data->novar - 1)) ||
       (status = CPXgetctype (DDSIP_env, DDSIP_lp, vartype, 0, DDSIP_data->novar - 1)))
    {
        fprintf (stderr, "Coud not get problem data, returned %d\n", status);
        goto FREE;
    }
    // check whether there are ranged rows
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
        if (!(rng = (double *) DDSIP_Alloc (sizeof (double), DDSIP_data->nocon, "rng(DetEqu)")) ||
            !(rng_sorted = (double *) DDSIP_Alloc (sizeof (double), DDSIP_Imax(DDSIP_data->firstcon,DDSIP_data->seccon), "rng_sorted(DetEqu)")))
        {
            fprintf (stderr, "Not enough memory for building deterministic equivalent\n");
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
    det_equ = CPXcreateprob (DDSIP_env, &status, probname);
    if (status)
    {
        fprintf (stderr, "CPXcreateprob returned %d\n", status);
        goto FREE;
    }
    // add (original) first-stage variables
    for (j = 0; j < DDSIP_data->firstvar; j++)
    {
        vartype_sorted[j]   = vartype[DDSIP_bb->firstindex[j]];
        lb_sorted[j]        = lb[DDSIP_bb->firstindex[j]];
        ub_sorted[j]        = ub[DDSIP_bb->firstindex[j]];
        scaled_obj_coef[j]  = obj_coef[DDSIP_bb->firstindex[j]];
        scen_spec_colname[j]= colname[DDSIP_bb->firstindex[j]];
    }
    if ((status = CPXnewcols (DDSIP_env, det_equ, DDSIP_data->firstvar, scaled_obj_coef,
                        lb_sorted, ub_sorted, vartype_sorted, scen_spec_colname)))
    {
        fprintf (stderr, "CPXnewcols returned %d for first-stage variables\n", status);
        goto FREE;
    }
    // add (original) second-stage variables for all scenarios
    for (j = 0; j < DDSIP_data->secvar; j++)
    {
        vartype_sorted[j] = vartype[DDSIP_bb->secondindex[j]];
        lb_sorted[j]      = lb[DDSIP_bb->secondindex[j]];
        ub_sorted[j]      = ub[DDSIP_bb->secondindex[j]];
    }
    for (scen = 0; scen < DDSIP_param->scenarios; scen++)
    {
        for (j = 0; j < DDSIP_data->secvar; j++)
        {
            if (!(string2 = (char *) calloc (1, 255 * sizeof (char))))
            {
                fprintf (stderr, "Not enough memory for building deterministic equivalent\n");
                goto FREE;
            }
            // append scenario index to colname
            string1 = colname[DDSIP_bb->secondindex[j]];
            sprintf (string2, "%sSC%.3d", string1, scen+1);
            scen_spec_colname[j] = string2;
            scaled_obj_coef[j] = DDSIP_data->prob[scen] * obj_coef[DDSIP_bb->secondindex[j]];
        }
        if ((status = CPXnewcols (DDSIP_env, det_equ, DDSIP_data->secvar, scaled_obj_coef,
                        lb_sorted, ub_sorted, vartype_sorted, scen_spec_colname)))
        {
            fprintf (stderr, "CPXnewcols returned %d for second-stage variables of scenario %d\n", status, scen+1);
            goto FREE;
        }
        for (j = 0; j < DDSIP_data->secvar; j++)
            DDSIP_Free ((void **) &(scen_spec_colname[j]));
    }
    ///////enter stochastic cost coefficients //////////////////////////////////////////////////
    if (DDSIP_param->stoccost)
    {
        for (j = 0; j < DDSIP_param->stoccost; j++)
        {
            obj_coef[j] = 0.0;
            if ((colindex_sorted[j] = DDSIP_bb->firstindex_reverse[DDSIP_data->costind[j]]))
                 colindex_sorted[j] = DDSIP_data->firstvar + DDSIP_bb->secondindex_reverse[DDSIP_data->costind[j]];
        }
        for (scen = 0; scen < DDSIP_param->scenarios; scen++)
        {
            for (j = 0; j < DDSIP_param->stoccost; j++)
            {
                if (colindex_sorted[j] >= DDSIP_data->firstvar)
                    obj_coef[j] = DDSIP_data->prob[scen] * DDSIP_data->cost[scen * DDSIP_param->stoccost + j];
                else
                    obj_coef[j] += DDSIP_data->prob[scen] * DDSIP_data->cost[scen * DDSIP_param->stoccost + j];
            }
            status = CPXchgobj (DDSIP_env, det_equ, DDSIP_param->stoccost, colindex_sorted, obj_coef);
            if (status)
            {
                char errmsg[1024];
                CPXgeterrorstring (DDSIP_env, status, errmsg);
                fprintf (stderr, "in DetEqu: %s\n", errmsg);
            }
            for (j = 0; j < DDSIP_param->stoccost; j++)
            {
                if (colindex_sorted[j] >= DDSIP_data->firstvar)
                    colindex_sorted[j] += DDSIP_data->secvar;
            }
        }
    }
    //
    // free arrays needeed only for columns
    DDSIP_Free ((void **) &(vartype));
    DDSIP_Free ((void **) &(colname));
    DDSIP_Free ((void **) &(colnamestore));
    DDSIP_Free ((void **) &(lb));
    DDSIP_Free ((void **) &(ub));
    DDSIP_Free ((void **) &(vartype_sorted));
    DDSIP_Free ((void **) &(lb_sorted));
    DDSIP_Free ((void **) &(ub_sorted));
    DDSIP_Free ((void **) &(obj_coef));
    DDSIP_Free ((void **) &(scaled_obj_coef));
    //
    // get problem matrix coefficients
    // query the length needed for storage of coefficients
    CPXgetrows(DDSIP_env, DDSIP_lp, &nzcnt, rmatbeg, rmatind, rmatval, 0, &rowsurplus, 0, DDSIP_data->nocon-1);
    nzcnt = -rowsurplus;
    if (!(rmatbeg = (int *) DDSIP_Alloc (sizeof (int), DDSIP_data->nocon, "rmatbeg(DetEqu)")) ||
        !(rmatind = (int *) DDSIP_Alloc (sizeof (int), DDSIP_Imax(nzcnt, DDSIP_param->stocmat), "rmatind(DetEqu)")) ||
        !(rmatval = (double *) DDSIP_Alloc (sizeof (double), nzcnt, "rmatval(DetEqu)")))
    {
        fprintf (stderr, "Not enough memory for building deterministic equivalent\n");
        goto FREE;
    }
    CPXgetrows(DDSIP_env, DDSIP_lp, &nzcnt, rmatbeg, rmatind, rmatval, nzcnt, &rowsurplus, 0, DDSIP_data->nocon-1);

    printf(" got %d elements of the matrix\n", nzcnt);

    if (!(rmatbeg_stage = (int *) DDSIP_Alloc (sizeof (int), DDSIP_Imax(DDSIP_data->firstcon, DDSIP_data->seccon), "rmatbeg_stage(DetEqu)")) ||
        !(rmatind_stage = (int *) DDSIP_Alloc (sizeof (int), nzcnt + DDSIP_param->stocmat, "rmatind_stage(DetEqu)")) ||
        !(rmatval_stage = (double *) DDSIP_Alloc (sizeof (double), nzcnt + DDSIP_param->stocmat, "rmatval_stage(DetEqu)")))
    {
        fprintf (stderr, "Not enough memory for building deterministic equivalent\n");
        goto FREE;
    }

    // add first-stage constraints
    k = 0;
    for (j = 0; j < DDSIP_data->firstcon; j++)
    {
        sense_sorted[j] = sense[DDSIP_bb->firstrowind[j]];
        det_equ_rhs[j] = base_rhs[DDSIP_bb->firstrowind[j]];
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
//fprintf(stderr,"first-stage row %d: koeff %d in col %d -> rmatind[%d]= %d, value %g (an Stelle %d+%d)\n", j, i, rmatind_stage[k + i],rmatbeg[DDSIP_bb->firstrowind[j]] + i, rmatind[rmatbeg[DDSIP_bb->firstrowind[j]] + i],   rmatval_stage[k + i], k, i);
        }
	k += nzcnt_row;
    }
    if ((status = CPXaddrows(DDSIP_env, det_equ, 0, DDSIP_data->firstcon, k, det_equ_rhs, sense_sorted, rmatbeg_stage, rmatind_stage, rmatval_stage, NULL, scen_spec_rowname)))
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
        if((status = CPXchgrngval(DDSIP_env, det_equ, DDSIP_data->firstcon, rowindex, rng_sorted)))
        {
            fprintf (stderr, "CPXchgrngval returned %d for first-stage constraints\n", status);
            goto FREE;
        }
    }

    // add second-stage constraints
    for (scen = 0; scen < DDSIP_param->scenarios; scen++)
    {
        k = 0;
        for (j = 0; j < DDSIP_data->seccon; j++)
        {
            sense_sorted[j] = sense[DDSIP_bb->secondrowind[j]];
            det_equ_rhs[j] = base_rhs[DDSIP_bb->secondrowind[j]];
            if (!(string2 = (char *) calloc (1, 255 * sizeof (char))))
            {
                fprintf (stderr, "Not enough memory for building deterministic equivalent\n");
                goto FREE;
            }
            // append scenario index to colname
            string1 = rowname[DDSIP_bb->secondrowind[j]];
            sprintf (string2, "%sSC%.3d", string1, scen+1);
            scen_spec_rowname[j] = string2;
            rmatbeg_stage[j] = k;
            if (DDSIP_bb->secondrowind[j] == DDSIP_data->nocon -1)
                nzcnt_row = nzcnt - rmatbeg[DDSIP_data->nocon -1];
            else
            {
                nzcnt_row = rmatbeg[DDSIP_bb->secondrowind[j]+1] - rmatbeg[DDSIP_bb->secondrowind[j]];

            }
            for (i = 0; i < nzcnt_row; i++)
            {
                if (DDSIP_bb->firstindex_reverse[rmatind[rmatbeg[DDSIP_bb->secondrowind[j]] + i]] < 0)
                    rmatind_stage[k + i] = DDSIP_data->firstvar + scen*DDSIP_data->secvar + DDSIP_bb->secondindex_reverse[rmatind[rmatbeg[DDSIP_bb->secondrowind[j]] + i]];
                else
                    rmatind_stage[k + i] = DDSIP_bb->firstindex_reverse[rmatind[rmatbeg[DDSIP_bb->secondrowind[j]] + i]];
                rmatval_stage[k + i] = rmatval[rmatbeg[DDSIP_bb->secondrowind[j]] + i];
            }
            k += nzcnt_row;
        }
        ///////enter stochastic rhs entries//////////////////////////////////////////////////////
        for (j=0; j< DDSIP_param->stocrhs; j++)
        {
            det_equ_rhs[DDSIP_bb->secondrowind_reverse[DDSIP_data->rhsind[j]]] = DDSIP_data->rhs[scen * DDSIP_param->stocrhs + j];
        }

        if ((status = CPXaddrows(DDSIP_env, det_equ, 0, DDSIP_data->seccon, k, det_equ_rhs, sense_sorted, rmatbeg_stage, rmatind_stage, rmatval_stage, NULL, scen_spec_rowname)))
        {
            fprintf (stderr, "CPXaddrows returned %d for second-stage constraints scenario %d\n", status, scen+1);
            goto FREE;
        }
        for (j = 0; j < DDSIP_data->seccon; j++)
            DDSIP_Free ((void **) &(scen_spec_rowname[j]));
        if (ranged)
        {
            for (j = 0; j < DDSIP_data->seccon; j++)
            {
                rng_sorted[j] = rng[DDSIP_bb->secondrowind[j]];
                rowindex[j]   = DDSIP_data->firstcon + scen * DDSIP_data->seccon + j;
            }
            if ((status = CPXchgrngval(DDSIP_env, det_equ, DDSIP_data->seccon, rowindex, rng_sorted)))
            {
                fprintf (stderr, "CPXchgrngval returned %d for first-stage constraints\n", status);
                goto FREE;
            }
        }
    }

    ///////enter stochastic matrix entries//////////////////////////////////////////////////////
    if (DDSIP_param->stocmat)
    {

        if (!(value = (double *) calloc (DDSIP_param->stocmat, sizeof (double))))
        {
            fprintf (stderr, "Not enough memory for building deterministic equivalent\n");
            goto FREE;
        }
        for (j = 0; j < DDSIP_param->stocmat; j++)
        {
            if ((colindex_sorted[j] = DDSIP_bb->firstindex_reverse[DDSIP_data->matcol[j]]))
                 colindex_sorted[j] = DDSIP_data->firstvar + DDSIP_bb->secondindex_reverse[DDSIP_data->matcol[j]];
            rmatind[j] = DDSIP_data->firstcon + DDSIP_bb->secondrowind_reverse[DDSIP_data->matrow[j]];
        }
        for (scen = 0; scen < DDSIP_param->scenarios; scen++)
        {
            for (j = 0; j < DDSIP_param->stocmat; j++)
            {
                value[j] = DDSIP_data->matval[scen * DDSIP_param->stocmat + j];
            }
///////////////////////////
//fprintf(stderr,"matrix scen %d:\n",scen+1);
//for (j = 0; j < DDSIP_param->stocmat; j++)
//fprintf(stderr,"  row %d col %d : %g\n",rmatind[j],colindex_sorted[j],value[j]);
///////////////////////////
            status = CPXchgcoeflist (DDSIP_env, det_equ, DDSIP_param->stocmat, rmatind, colindex_sorted, value);
            if (status)
            {
                char errmsg[1024];
                CPXgeterrorstring (DDSIP_env, status, errmsg);
                fprintf (stderr, "in DetEqu chgcoeflist returned %d: %s\n", status, errmsg);
            }
            for (j = 0; j < DDSIP_param->stocmat; j++)
            {
                rmatind[j] += DDSIP_data->seccon;
                if (colindex_sorted[j] >= DDSIP_data->firstvar)
                    colindex_sorted[j] += DDSIP_data->secvar;
            }
        }
        DDSIP_Free ((void **) &(value));
    }

    status = CPXwriteprob (DDSIP_env, det_equ, probname, NULL);
    if (status)
        fprintf (DDSIP_outfile, " *** Deterministic equivalent not written successfully, status = %d\n", status);
    else
        fprintf (DDSIP_outfile, " *** Deterministic equivalent written successfully\n");
    status = CPXfreeprob (DDSIP_env, &det_equ);

FREE:
    DDSIP_Free ((void **) &(sense));
    DDSIP_Free ((void **) &(sense_sorted));
    DDSIP_Free ((void **) &(vartype));
    DDSIP_Free ((void **) &(rowname));
    DDSIP_Free ((void **) &(rownamestore));
    DDSIP_Free ((void **) &(colname));
    DDSIP_Free ((void **) &(colnamestore));
    DDSIP_Free ((void **) &(det_equ_rhs));
    DDSIP_Free ((void **) &(base_rhs));
    DDSIP_Free ((void **) &(lb));
    DDSIP_Free ((void **) &(ub));
    DDSIP_Free ((void **) &(vartype_sorted));
    DDSIP_Free ((void **) &(lb_sorted));
    DDSIP_Free ((void **) &(ub_sorted));
    DDSIP_Free ((void **) &(obj_coef));
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
