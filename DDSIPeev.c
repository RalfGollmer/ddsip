/*  Authors:           Andreas M"arkert, Ralf Gollmer
	Copyright to:      University of Duisburg-Essen
    Language:          C
    Last modification: 06.02.2016
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
    int status, j, mipstatus;

    char fname[DDSIP_ln_fname];

    double *mipx = (double *) DDSIP_Alloc (sizeof (double), (DDSIP_bb->firstvar + DDSIP_bb->secvar),
                                           "mipx(ExpValProb)");

    printf ("Solving expected value problem\n");
    if (DDSIP_param->outlev)
        fprintf (DDSIP_bb->moreoutfile, "Solving expected value problem...\n");

    status = DDSIP_ChgProb (-1);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to change problem \n");
        goto TERMINATE;
    }

    if (DDSIP_param->files > 1)
    {
        sprintf (fname, "%s/ev%s", DDSIP_outdir, DDSIP_param->coretype);
        status = CPXwriteprob (env, lp, fname, NULL);
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

    mipstatus = CPXmipopt (env, lp);

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
    mipstatus = CPXgetstat (env, lp);
    if (DDSIP_NoSolution (mipstatus))
    {
        status = 1;
        goto TERMINATE;
    }

    status = CPXgetx (env, lp, mipx, 0, DDSIP_bb->firstvar + DDSIP_bb->secvar - 1);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to get solution \n");
        goto TERMINATE;
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
