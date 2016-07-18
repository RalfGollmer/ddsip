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

    int status, scen, i, j, k;
    char probname[] = "sipout/det_equ.lp.gz";
    double *obj_coef;
    double *scaled_obj_coef;
    char *sense;
    char **scen_spec_rowname;
    char **scen_spec_colname;
    char **rowname, *rownamestore;
    char **colname, *colnamestore;
    int rowstorespace, rowsurplus_p;
    int colstorespace, colsurplus_p;
    char *string1, *string2;
    double coef;
    double *lb, *lb_sorted;
    double *ub, *ub_sorted;
    char *vartype, *vartype_sorted;
    int *colindex_sorted, *colindex_revers, *matcol_sorted;
    double *value;
    double *det_equ_rhs = NULL;
    double *non_stoc_rhs = NULL;
    if (DDSIP_data->seccon)
        det_equ_rhs = (double *) DDSIP_Alloc(sizeof(double),DDSIP_data->seccon,"det_equ_rhs(DetEqu)");
    else
    {
        fprintf (stderr,"XXX ERROR: no second stage contraints, got DDSIP_data->seccon=%d.\n",DDSIP_data->seccon);
        exit (1);
    }
    if (DDSIP_data->seccon - DDSIP_param->stocrhs>0)
        non_stoc_rhs = (double *) DDSIP_Alloc(sizeof(double),DDSIP_data->seccon - DDSIP_param->stocrhs,"non_stoc_rhs(DetEqu)");

    fprintf (stderr,
             "\nBuilding deterministic equivalent. This may take some time.\nWorks only for expectation-based model so far.\n");

    if (!(sense = (char *) calloc (DDSIP_data->seccon, sizeof (char))))
    {
        fprintf (stderr, "Not enough memory for building deterministic equivalent\n");
        return;
    }

    if (!(obj_coef = (double *) calloc (DDSIP_data->novar, sizeof (double))) ||
            !(scaled_obj_coef = (double *) calloc (DDSIP_bb->secvar, sizeof (double))))
    {
        fprintf (stderr, "Not enough memory for building deterministic equivalent\n");
        return;
    }


    det_equ = CPXcloneprob (DDSIP_env, DDSIP_lp, &status);
    CPXchgprobname (DDSIP_env, det_equ, probname);

    if (!(rowname = (char **) calloc (DDSIP_data->seccon, sizeof (char *)))
            || !(scen_spec_rowname = (char **) calloc (DDSIP_data->seccon, sizeof (char *)))
            || !(rownamestore = (char *) calloc (DDSIP_data->seccon * 255, sizeof (char))))
    {
        fprintf (stderr, "Not enough memory for building deterministic equivalent\n");
        return;
    }
    rowstorespace = DDSIP_data->seccon * 255;
    status = CPXgetrowname (DDSIP_env, DDSIP_lp, rowname, rownamestore,
                            rowstorespace, &rowsurplus_p, DDSIP_data->firstcon, DDSIP_data->nocon + DDSIP_data->seccon - 1);


    if (!(colname = (char **) calloc (DDSIP_data->novar, sizeof (char *)))
            || !(scen_spec_colname = (char **) calloc (DDSIP_data->secvar, sizeof (char *)))
            || !(colnamestore = (char *) calloc (DDSIP_data->secvar * 255, sizeof (char))))
    {
        fprintf (stderr, "Not enough memory for building deterministic equivalent\n");
        return;
    }
    colstorespace = (DDSIP_data->novar) * 255;
    status = CPXgetcolname (DDSIP_env, DDSIP_lp, colname, colnamestore,
                            colstorespace, &colsurplus_p, 0, DDSIP_data->novar - 1);

    /*____________________________________________________________________________________*/
    status = CPXgetsense (DDSIP_env, DDSIP_lp, sense, DDSIP_data->firstcon, DDSIP_data->nocon - 1);
    /*____________________________________________________________________________________*/
    status = CPXgetrhs (DDSIP_env, DDSIP_lp, non_stoc_rhs, DDSIP_data->firstcon + DDSIP_param->stocrhs, DDSIP_data->nocon - 1);
    /*____________________________________________________________________________________*/
    status = CPXgetobj (DDSIP_env, DDSIP_lp, obj_coef, 0, DDSIP_data->novar - 1);
    /*____________________________________________________________________________________*/
    //copy rownames scenario many times, append scenario index
    //and enter sense and rhs
    for (scen = 0; scen < DDSIP_param->scenarios; scen++)
    {
        for (j = 0; j < DDSIP_data->seccon; j++)
        {
            if (!(string2 = (char *) calloc (1, 255 * sizeof (char))))
            {
                fprintf (stderr, "Not enough memory for building deterministic equivalent\n");
                return;
            }
            string1 = rowname[j];
            sprintf (string2, "%sSC%.3d", string1, scen);
            scen_spec_rowname[j] = string2;

            if (j < DDSIP_param->stocrhs)
                det_equ_rhs[j] = DDSIP_data->rhs[DDSIP_param->stocrhs * scen + j];
            else
                det_equ_rhs[j] = non_stoc_rhs[j - DDSIP_param->stocrhs];
        }
        status = CPXnewrows (DDSIP_env, det_equ, DDSIP_data->seccon, det_equ_rhs, sense, NULL, scen_spec_rowname);
        for (j = 0; j < DDSIP_data->seccon; j++)
            DDSIP_Free ((void **) &(scen_spec_rowname[j]));
    }

    //copy colnames scenario many times, append scenario index
    //and enter into constraint matrix
    if (!(lb = (double *) calloc (DDSIP_data->novar, sizeof (double)))
            || !(lb_sorted = (double *) calloc (DDSIP_data->secvar, sizeof (double)))
            || !(ub = (double *) calloc (DDSIP_data->novar, sizeof (double)))
            || !(ub_sorted = (double *) calloc (DDSIP_data->secvar, sizeof (double)))
            || !(vartype = (char *) calloc (DDSIP_data->novar, sizeof (char)))
            || !(vartype_sorted = (char *) calloc (DDSIP_data->secvar, sizeof (double)))
            || !(colindex_revers = (int *) calloc (DDSIP_data->novar, sizeof (int)))
            || !(colindex_sorted = (int *) calloc (DDSIP_data->novar, sizeof (int))))
    {
        fprintf (stderr, "Not enough memory for building deterministic equivalent\n");
        return;
    }

    status = CPXgetlb (DDSIP_env, det_equ, lb, 0, DDSIP_data->novar - 1);
    status = CPXgetub (DDSIP_env, det_equ, ub, 0, DDSIP_data->novar - 1);
    status = CPXgetctype (DDSIP_env, det_equ, vartype, 0, DDSIP_data->novar - 1);
    for (j = 0; j < DDSIP_data->secvar; j++)
    {
        vartype_sorted[j] = vartype[DDSIP_bb->secondindex[j]];
        lb_sorted[j] = lb[DDSIP_bb->secondindex[j]];
        ub_sorted[j] = ub[DDSIP_bb->secondindex[j]];
    }


    for (scen = 0; scen < DDSIP_param->scenarios; scen++)
    {
        for (j = 0; j < DDSIP_data->secvar; j++)
        {
            if (!(string2 = (char *) calloc (1, 255 * sizeof (char))))
            {
                fprintf (stderr, "Not enough memory for building deterministic equivalent\n");
                return;
            }
            string1 = colname[DDSIP_bb->secondindex[j]];
            sprintf (string2, "%sSC%.3d", string1, scen);
            scen_spec_colname[j] = string2;
            scaled_obj_coef[j] = DDSIP_data->prob[scen] * obj_coef[DDSIP_bb->secondindex[j]];
        }

        status =
            CPXnewcols (DDSIP_env, det_equ, DDSIP_data->secvar, scaled_obj_coef,
                        lb_sorted, ub_sorted, vartype_sorted, scen_spec_colname);
        for (j = 0; j < DDSIP_data->secvar; j++)
            DDSIP_Free ((void **) &(scen_spec_colname[j]));

    }

    /////////////////////////////////////////////////
    for (j = 0; j < DDSIP_data->firstvar; j++)
    {
        colindex_sorted[j] = DDSIP_bb->firstindex[j];
    }
    for (j = 0; j < DDSIP_data->secvar; j++)
    {
        colindex_sorted[DDSIP_data->firstvar + j] = DDSIP_bb->secondindex[j];
    }
    for (j = 0; j < DDSIP_data->firstvar + DDSIP_data->secvar; j++)
    {
        colindex_revers[colindex_sorted[j]] = j;
    }

    k = DDSIP_data->seccon / 60;
    printf ("\n0%%                                                         100%%\n");
    for (i = 0; i < DDSIP_data->seccon; i++)
    {
        for (j = 0; j < DDSIP_data->firstvar; j++)
        {
            if ((status = CPXgetcoef (DDSIP_env, det_equ, DDSIP_data->firstcon + i, colindex_sorted[j], &coef)))
            {
                fprintf (stderr, " Build det. equivalent: Error retrieving coefficient of first-stage Variable %d.\n", j);
                exit (1);
            }
            if (coef)
            {
                for (scen = 0; scen < DDSIP_param->scenarios; scen++)
                {
                    status =
                        CPXchgcoef (DDSIP_env, det_equ, DDSIP_data->firstcon + DDSIP_bb->seccon + scen * DDSIP_data->seccon + i, colindex_sorted[j], coef);
                    if (status)
                    {
                        fprintf (stderr, " Build det. equivalent: Error setting coefficient of first-stage Variable %d.\n", j);
                        exit (1);
                    }
                }
            }
        }
        for (j = DDSIP_data->firstvar; j < DDSIP_data->firstvar + DDSIP_data->secvar; j++)
        {
            if ((status = CPXgetcoef (DDSIP_env, det_equ, DDSIP_data->firstcon + i, colindex_sorted[j], &coef)))
            {
                fprintf (stderr,
                         " Build det. equivalent: Error retrieving coefficient of second-stage Variable %d.\n",
                         j - DDSIP_data->firstvar);
                exit (1);
            }
            if (coef)
            {
                for (scen = 0; scen < DDSIP_param->scenarios; scen++)
                {
                    status =
                        CPXchgcoef (DDSIP_env, det_equ,
                                    DDSIP_data->firstcon + DDSIP_bb->seccon + scen * DDSIP_data->seccon + i, (scen + 1) * DDSIP_data->secvar + j, coef);
                }
                if (status)
                {
                    fprintf (stderr,
                             " Build det. equivalent: Error setting coefficient of second-stage Variable %d.\n",
                             j - DDSIP_data->firstvar);
                    exit (1);
                }
            }
        }
        if (!k)
        {
            for (j = 0; j <= 60 / DDSIP_data->seccon; j++)
                printf ("#");
        }
        else if (i % k == k - 1)
            printf ("#");
    }
    printf ("\n\n");

    ///////delete original second stage rows & cols ////////////////////////////////////////////

    status = CPXdelrows (DDSIP_env, det_equ, DDSIP_data->firstcon, DDSIP_data->firstcon + DDSIP_bb->seccon - 1);
    j = 0;
    for (i = 0; i < DDSIP_data->secvar; i++)
    {
        status = CPXdelcols (DDSIP_env, det_equ, DDSIP_bb->secondindex[i] - j, DDSIP_bb->secondindex[i] - j);
        j++;
    }

    ///////enter stochastic matrix entries//////////////////////////////////////////////////////
    if (DDSIP_param->stocmat)
    {

        if (!(value = (double *) calloc (DDSIP_param->stocmat, sizeof (double)))
                || !(matcol_sorted = (int *) calloc (DDSIP_param->stocmat, sizeof (int))))
        {
            fprintf (stderr, "Not enough memory for building deterministic equivalent\n");
            return;
        }
        for (j = 0; j < DDSIP_param->stocmat; j++)
        {
            matcol_sorted[j] = colindex_revers[DDSIP_data->matcol[j]];
        }
        for (scen = 0; scen < DDSIP_param->scenarios; scen++)
        {
            for (j = 0; j < DDSIP_param->stocmat; j++)
            {
                value[j] = DDSIP_data->matval[scen * DDSIP_param->stocmat + j];
            }
            status = CPXchgcoeflist (DDSIP_env, det_equ, DDSIP_param->stocmat, DDSIP_data->matrow, matcol_sorted, value);
            if (status)
            {
                char errmsg[1024];
                CPXgeterrorstring (DDSIP_env, status, errmsg);
                fprintf (stderr, "in DetEqu: %s\n", errmsg);
            }
            for (j = 0; j < DDSIP_param->stocmat; j++)
            {
                DDSIP_data->matrow[j] += DDSIP_data->seccon;
                if (matcol_sorted[j] >= DDSIP_data->firstvar)
                    matcol_sorted[j] += DDSIP_data->secvar;
            }
        }
        DDSIP_Free ((void **) &(value));
        DDSIP_Free ((void **) &(matcol_sorted));
        //set matrow to the old values
        for (j = 0; j < DDSIP_param->stocmat; j++)
        {
            DDSIP_data->matrow[j] -= DDSIP_param->scenarios * DDSIP_data->seccon;
        }

    }
    ///////enter stochastic cost coefficients //////////////////////////////////////////////////
    if (DDSIP_param->stoccost)
    {

        if (!(value = (double *) calloc (DDSIP_param->stoccost, sizeof (double)))
                || !(matcol_sorted = (int *) calloc (DDSIP_param->stoccost, sizeof (int))))
        {
            fprintf (stderr, "Not enough memory for building deterministic equivalent\n");
            return;
        }
        for (j = 0; j < DDSIP_param->stoccost; j++)
        {
            value[j] = 0.0;
            matcol_sorted[j] = colindex_revers[DDSIP_data->costind[j]];
        }
        for (scen = 0; scen < DDSIP_param->scenarios; scen++)
        {
            for (j = 0; j < DDSIP_param->stoccost; j++)
            {
                if (matcol_sorted[j] >= DDSIP_data->firstvar)
                    value[j] = DDSIP_data->prob[scen] * DDSIP_data->cost[scen * DDSIP_param->stoccost + j];
                else
                    value[j] += DDSIP_data->prob[scen] * DDSIP_data->cost[scen * DDSIP_param->stoccost + j];
            }
            status = CPXchgobj (DDSIP_env, det_equ, DDSIP_param->stoccost, matcol_sorted, value);
            if (status)
            {
                char errmsg[1024];
                CPXgeterrorstring (DDSIP_env, status, errmsg);
                fprintf (stderr, "in DetEqu: %s\n", errmsg);
            }
            for (j = 0; j < DDSIP_param->stoccost; j++)
            {
                if (matcol_sorted[j] >= DDSIP_data->firstvar)
                    matcol_sorted[j] += DDSIP_data->secvar;
            }
        }
        DDSIP_Free ((void **) &(value));
        DDSIP_Free ((void **) &(matcol_sorted));

    }
    ////////////////////////////////////////////////////////////////////////////////////////////

    status = CPXwriteprob (DDSIP_env, det_equ, probname, NULL);
    if (status)
        fprintf (DDSIP_outfile, " *** Deterministic equivalent not written successfully, status = %d\n", status);
    else
        fprintf (DDSIP_outfile, " *** Deterministic equivalent written successfully\n");
    status = CPXfreeprob (DDSIP_env, &det_equ);

    DDSIP_Free ((void **) &(sense));
    DDSIP_Free ((void **) &(vartype));
    DDSIP_Free ((void **) &(rowname));
    DDSIP_Free ((void **) &(rownamestore));
    DDSIP_Free ((void **) &(colname));
    DDSIP_Free ((void **) &(colnamestore));
    DDSIP_Free ((void **) &(det_equ_rhs));
    DDSIP_Free ((void **) &(non_stoc_rhs));
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
    return;
}
