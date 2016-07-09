/*  Authors:           Andreas M"arkert, Ralf Gollmer
	Copyright to:      University of Duisburg-Essen
   Language:          C
   Last modification: 06.02.2016
	Description:
	This file contains procedures to restore the original problem during the
	branch-and-bound algorithm

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
// Restore bounds and type of the first stage variables
int
DDSIP_RestoreBoundAndType (void)
{
    int status = 0;

    if ((DDSIP_bb->DDSIP_step == neobj || DDSIP_bb->DDSIP_step == eev) && DDSIP_param->riskvar)
        DDSIP_UndeleteRiskObj ();

    /*   if (DDSIP_bb->DDSIP_step==solve || DDSIP_bb->DDSIP_step==dual) */
    /* 	{ */
    /* 	  status = CPXchgobj (DDSIP_env, DDSIP_lp, DDSIP_bb->firstvar, DDSIP_bb->firstindex, DDSIP_bb->cost); */
    /* 	  if ( status ) { */
    /* 		fprintf (stderr,"ERROR: Failed to update objective coefficients\n"); */
    /* 		return status; */
    /* 	  } */
    /* 	} */

    // Restore original bounds
    status = CPXchgbds (DDSIP_env, DDSIP_lp, DDSIP_bb->firstvar, DDSIP_bb->firstindex, DDSIP_bb->lbident, DDSIP_bb->lborg);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to change lower bounds \n");
        return status;
    }

    status = CPXchgbds (DDSIP_env, DDSIP_lp, DDSIP_bb->firstvar, DDSIP_bb->firstindex, DDSIP_bb->ubident, DDSIP_bb->uborg);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to change upper bounds \n");
        return status;
    }
    // probtype=0 (LP)
    if (!CPXgetprobtype (DDSIP_env, DDSIP_lp))
    {
        status = CPXchgprobtype (DDSIP_env, DDSIP_lp, CPXPROB_MILP);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change problem type (Restore) \n");
            return status;
        }
        status = CPXchgctype (DDSIP_env, DDSIP_lp, DDSIP_bb->secvar, DDSIP_bb->secondindex, DDSIP_bb->sectype);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change types of second-stage variables (Restore) \n");
            return status;
        }
    }
    //Restore ctypes if solving was processed with relaxed first stage
    if (DDSIP_bb->DDSIP_step == solve && DDSIP_param->relax == 1)
    {
        status = CPXchgctype (DDSIP_env, DDSIP_lp, DDSIP_bb->firstvar, DDSIP_bb->firstindex, DDSIP_bb->firsttype);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change types of first stage variables (Restore) \n");
            return status;
        }
    }
    return status;
}
