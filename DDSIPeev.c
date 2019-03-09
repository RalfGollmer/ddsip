/*  Authors:           Andreas M"arkert, Ralf Gollmer
	Copyright to:      University of Duisburg-Essen
    Language:          C
	Description:
	LowerBound expected value problem (EV), i.e. the deterministic problem
	with the expected value of the scenarios.

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
// This function solves the expected value problem
int
DDSIP_ExpValProb (void)
{
    int status, j, mipstatus, nodes_1st = -1;
    int wall_hrs, wall_mins,cpu_hrs, cpu_mins;
    double objval, bobjval, time_start, time_end, wall_secs, cpu_secs, gap;

    char fname[DDSIP_ln_fname];

    double *mipx = (double *) DDSIP_Alloc (sizeof (double), (DDSIP_bb->firstvar + DDSIP_bb->secvar),
                                           "mipx(ExpValProb)");

    printf ("Solving expected value problem\n");
    if (DDSIP_param->outlev)
        fprintf (DDSIP_bb->moreoutfile, "Solving expected value problem...\n");

    status = DDSIP_ChgProb (-1, 0);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to change problem \n");
        goto TERMINATE;
    }

    if (DDSIP_param->files > 1)
    {
        sprintf (fname, "%s/ev%s", DDSIP_outdir, DDSIP_param->coretype);
        status = CPXwriteprob (DDSIP_env, DDSIP_lp, fname, NULL);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to write problem\n");
            goto TERMINATE;
        }
    }
    // New cplex parameters
    if (DDSIP_param->cpxnoeev)
    {
        status = DDSIP_SetCpxPara (DDSIP_param->cpxnoeev, DDSIP_param->cpxeevisdbl, DDSIP_param->cpxeevwhich, DDSIP_param->cpxeevwhat);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to set CPLEX parameters (ExpValProb) \n");
            goto TERMINATE;
        }
    }

    time_start = DDSIP_GetCpuTime ();
    mipstatus = CPXmipopt (DDSIP_env, DDSIP_lp);

    // Reset cplex parameters
    if (DDSIP_param->cpxnoeev)
    {
        status = DDSIP_SetCpxPara (DDSIP_param->cpxno, DDSIP_param->cpxisdbl, DDSIP_param->cpxwhich, DDSIP_param->cpxwhat);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to reset CPLEX parameters (ExpValProb) \n");
            goto TERMINATE;
        }
    }

    if (DDSIP_Error (mipstatus))
    {
        fprintf (stderr, "ERROR: Failed to optimize EXP\n");
        status = mipstatus;
        goto TERMINATE;
    }
    //  Error ? (blatant infeasible, scenario problem limit)
    if (DDSIP_Infeasible (mipstatus))
    {
        status = 1;
        goto TERMINATE;
    }
    // No solution found ? (integer infeasible, some limit (node, time))
    mipstatus = CPXgetstat (DDSIP_env, DDSIP_lp);
    if (DDSIP_NoSolution (mipstatus))
    {
        status = 1;
        goto TERMINATE;
    }

    status = CPXgetx (DDSIP_env, DDSIP_lp, mipx, 0, DDSIP_bb->firstvar + DDSIP_bb->secvar - 1);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to get solution \n");
        goto TERMINATE;
    }
// output of result
    if (DDSIP_param->outlev)
    {
        status = CPXgetobjval (DDSIP_env, DDSIP_lp, &objval);
        if (status)
        {
            fprintf (stderr, "ERROR*: Failed to get best objective value \n");
            fprintf (DDSIP_outfile, "ERROR*: Failed to get best objective value \n");
            if (DDSIP_param->outlev)
                fprintf (DDSIP_bb->moreoutfile, "ERROR*: Failed to get best objective value \n");
            goto TERMINATE;
        }
        if (mipstatus == CPXMIP_OPTIMAL)
        {
            bobjval = objval;
        }
        else
        {
            status = CPXgetbestobjval (DDSIP_env, DDSIP_lp, &bobjval);
            if (status)
            {
                fprintf (stderr, "ERROR: Failed to get value of best remaining node\n");
                fprintf (DDSIP_outfile, "ERROR: Failed to get value of best remaining node\n");
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to get value of best remaining node\n");
                goto TERMINATE;
            }
        }
        gap = 100.0*(objval-bobjval)/(fabs(objval)+1e-4);
        nodes_1st = CPXgetnodecnt (DDSIP_env,DDSIP_lp);
        time_end = DDSIP_GetCpuTime ();
        time_start = time_end-time_start;
        time (&DDSIP_bb->cur_time);
        DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
        DDSIP_translate_time (time_end,&cpu_hrs,&cpu_mins,&cpu_secs);
        if (mipstatus == CPXMIP_OPTIMAL)
            fprintf (DDSIP_bb->moreoutfile,
                 "    exp. val. prob:  Best=%-20.14g\tBound=%-20.14g (%9.4g%%)     \t %3dh %02d:%02.0f / %3dh %02d:%05.2f (%6.2fs n: %4d)",
                 objval, bobjval, gap,
                 wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs, time_start, nodes_1st);
        else if (mipstatus == CPXMIP_OPTIMAL_TOL)
            fprintf (DDSIP_bb->moreoutfile,
                 "    exp. val. prob:  Best=%-20.14g\tBound=%-20.14g (%9.4g%%) tol.\t %3dh %02d:%02.0f / %3dh %02d:%05.2f (%6.2fs n: %4d)",
                 objval, bobjval, gap,
                 wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs, time_start, nodes_1st);
        else if (mipstatus == CPXMIP_TIME_LIM_FEAS)
            fprintf (DDSIP_bb->moreoutfile,
                 "    exp. val. prob:  Best=%-20.14g\tBound=%-20.14g (%9.4g%%) TIME\t %3dh %02d:%02.0f / %3dh %02d:%05.2f (%6.2fs n: %4d)",
                 objval, bobjval, gap,
                 wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs, time_start, nodes_1st);
        else
            fprintf (DDSIP_bb->moreoutfile,
                 "    exp. val. prob:  Best=%-20.14g\tBound=%-20.14g (%9.4g%%) %-4d\t %3dh %02d:%02.0f / %3dh %02d:%05.2f (%6.2fs n: %4d)",
                 objval, bobjval, gap, mipstatus,
                 wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs, time_start, nodes_1st);
    }
    // Returns sometimes rubbish, don't know why..
    if (!DDSIP_bb->adv_sol)
        DDSIP_bb->adv_sol = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar, "DDSIP_bb->adv_sol(sipread)");
    for (j = 0; j < DDSIP_bb->firstvar; j++)
    {
        // Numerical errors ?
        if (DDSIP_Equal (mipx[DDSIP_bb->firstindex[j]], 0.0))
            mipx[DDSIP_bb->firstindex[j]] = 0.0;
        DDSIP_bb->adv_sol[j] = mipx[DDSIP_bb->firstindex[j]];
    }

TERMINATE:

    DDSIP_Free ((void **) &(mipx));
    return status;
}
