/*  Authors:           Andreas M"arkert, Ralf Gollmer
	Copyright to:      University of Duisburg-Essen
   Language:          C
   Last modification: 06.02.2016

	Description:
	This procedure prints a line of output when invoked.

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
#include <math.h>

// Print error message
void
DDSIP_PrintErrorMsg (int status)
{
    if (status == 113)
        fprintf (stderr, "ERROR: No node selected.\n");
    else if (status == 115)
        fprintf (stderr, "ERROR: Lower bound greater than upper bound (new nodes).\n");
    else if (status == 117)
        fprintf (stderr, "ERROR: Number of bounds exceeds number of first-stage variables.\n");
    else if (status == 119)
        fprintf (stderr, "ERROR: Lower bound greater than upper bound (current node).\n");
    else if (status == 121)
        fprintf (stderr, "ERROR: Inconsistent data.\n");
    else if (status == 123)
        fprintf (stderr, "ERROR: Wrong format.\n");
    else if (status == 127)
        fprintf (stderr, "ERROR: Inappropriate model.\n");
    else if (status == 129)
        fprintf (stderr, "ERROR: Inappropriate probability level.\n");
    else if (status == 131)
        fprintf (stderr, "ERROR: Failed to set branching index.\n");
    else if (status == 133)
        fprintf (stderr, "ERROR: FSD-algorithm and dual method cannot be combined.\n");
    /*    else if (status==135) */
    /*  DDSIP_	fprintf (stderr,"ERROR: Unsupported file extension of core file.\n"); */
    else if (status == 137)
        fprintf (stderr, "ERROR: Inconsistent start value(s).\n");
    else if (status > 1000)
    {
        char errmsg[1024];
        CPXgeterrorstring (env, status, errmsg);
        printf ("%s", errmsg);
    }

}

// Print message on stdout and to output file
void
DDSIP_Print2 (char b[], char e[], double d, int what)
{
    int ih;
    if (0 == what)
    {
        printf ("%s%s", b, e);
        fprintf (DDSIP_outfile, "%s%s", b, e);
    }
    else if (1 == what)
    {
        ih = floor (d + 0.01);
        printf ("%s%d%s", b, ih, e);
        fprintf (DDSIP_outfile, "%s%d%s", b, ih, e);
    }
    else
    {
        printf ("%s%6.2f%s", b, d, e);
        fprintf (DDSIP_outfile, "%s%6.2e%s", b, d, e);
    }
}

//==========================================================================
// Function prints status of optimization
// A headline is printed at every 10th printout

void
DDSIP_PrintState (int noiter)
{
    //  char    state[DDSIP_max_str_ln];
    double wall_secs, cpu_secs;
    int    wall_hrs, wall_mins, cpu_hrs, cpu_mins, print_violations = 1;
    double rgap, factor;
    //best<>0 ?

    if (!DDSIP_Equal (fabs (DDSIP_bb->bestvalue), 0.0))
        rgap = 100. * (DDSIP_bb->bestvalue - DDSIP_bb->bestbound) / fabs (DDSIP_bb->bestvalue);
    else
        rgap = 100. * (DDSIP_bb->bestvalue - DDSIP_bb->bestbound) / (fabs (DDSIP_bb->bestvalue) + DDSIP_param->accuracy);

    rgap = DDSIP_Dmin (rgap, 100.0);
    factor = (DDSIP_bb->bestvalue < 0.)? 1.-1.e-12 :  1.+1.e-12;

    // A headline is printed every 20th call
    if (!((noiter) % (DDSIP_param->logfreq * 20)) || (DDSIP_param->cb && !(noiter % DDSIP_Imax(DDSIP_param->cb,10)) && DDSIP_param->outlev > 1))
    {
        printf ("\n   Node   Nodes    Left  Objective         Heuristic");
        fprintf (DDSIP_outfile, "\n   Node   Nodes    Left  Objective         Heuristic");
        printf ("         Best Value        Bound             Viol./Dispersion        Gap   Wall Time    CPU Time   Father Depth\n");
        fprintf (DDSIP_outfile, "         Best Value        Bound             Viol./Dispersion        Gap   Wall Time    CPU Time   Father Depth\n");
    }

    if (!DDSIP_bb->violations && fabs(DDSIP_bb->bestvalue - DDSIP_node[DDSIP_bb->curnode]->bound)/(fabs(DDSIP_bb->bestvalue) + 3.e-16) < 5.e-16)
    {
        printf ("*%6d  %6d  %6d", DDSIP_bb->curnode, DDSIP_bb->nonode, DDSIP_bb->nofront);
        fprintf (DDSIP_outfile, "*%6d  %6d  %6d", DDSIP_bb->curnode, DDSIP_bb->nonode, DDSIP_bb->nofront);
    }
    else
    {
        printf (" %6d  %6d  %6d", DDSIP_bb->curnode, DDSIP_bb->nonode, DDSIP_bb->nofront);
        fprintf (DDSIP_outfile, " %6d  %6d  %6d", DDSIP_bb->curnode, DDSIP_bb->nonode, DDSIP_bb->nofront);
    }

    //bound=inf
    if (DDSIP_node[DDSIP_bb->curnode]->bound >= DDSIP_infty)
    {
        printf ("       infeasible ");
        fprintf (DDSIP_outfile, "       infeasible ");
        print_violations  = 0;
    }
    else if (DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bestvalue*factor)
    {
        printf ("       cutoff     ");
        fprintf (DDSIP_outfile, "       cutoff     ");
        print_violations  = 0;
    }
    else if (DDSIP_node[DDSIP_bb->curnode]->bound < DDSIP_infty)
    {
        printf ("  %-16.12g", DDSIP_node[DDSIP_bb->curnode]->bound);
        fprintf (DDSIP_outfile, "  %-16.12g", DDSIP_node[DDSIP_bb->curnode]->bound);
    }
    //bound=-inf
    else
    {
        printf ("        unbounded ");
        fprintf (DDSIP_outfile, "        unbounded ");
    }

    // DDSIP_bb->heurval contains the objective value of the heuristic solution
    // DDSIP_bb->skip indicates whether the evaluation of an heuristic solution was skipped for some reason
    if (fabs (DDSIP_bb->heurval) < DDSIP_infty)
    {
        printf ("  %-16.12g", DDSIP_bb->heurval);
        fprintf (DDSIP_outfile, "  %-16.12g", DDSIP_bb->heurval);
    }
    else
    {
        if (DDSIP_bb->skip >= 100)
        {
            // Print the number of the scenario which caused the stop
            printf ("         %4d-stop", DDSIP_bb->skip - 100 + 1);
            fprintf (DDSIP_outfile, "         %4d-stop", DDSIP_bb->skip - 100 + 1);
        }
        else if (DDSIP_bb->skip == 2 || DDSIP_bb->skip == 1 ||  DDSIP_bb->skip == 3)
        {
            printf ("                  ");
            fprintf (DDSIP_outfile, "                  ");
        }
        else if (DDSIP_bb->skip == 4)
        {
            printf ("     multiple     ");
            fprintf (DDSIP_outfile, "     multiple     ");
        }
        else
        {
            printf ("   infeasible     ");
            fprintf (DDSIP_outfile, "   infeasible     ");
        }
    }

    if (fabs (DDSIP_bb->bestvalue) < DDSIP_infty)
    {
        printf ("  %-16.12g", DDSIP_bb->bestvalue);
        fprintf (DDSIP_outfile, "  %-16.12g", DDSIP_bb->bestvalue);
    }
    else
    {
        printf ("                  ");
        fprintf (DDSIP_outfile, "                  ");
    }

    if (fabs (DDSIP_bb->bestbound - DDSIP_infty) < DDSIP_param->accuracy)
    {
        printf ("        infeasible");
        fprintf (DDSIP_outfile, "        infeasible");
    }
    else
    {
        printf ("  %-16.12g", DDSIP_bb->bestbound);
        fprintf (DDSIP_outfile, "  %-16.12g", DDSIP_bb->bestbound);
    }

    if (!print_violations || fabs (DDSIP_bb->bestbound - DDSIP_infty) < DDSIP_param->accuracy || fabs(DDSIP_node[DDSIP_bb->curnode]->dispnorm) >= DDSIP_infty)
    {
        printf ("                  ");
        fprintf (DDSIP_outfile, "                  ");
    }
    else
    {
        //number of violations of nonanticipativity
        printf ("  %4d",DDSIP_bb->violations);
        fprintf (DDSIP_outfile, "  %4d",DDSIP_bb->violations);
        //max dispersion of the variables
        printf ("  %-10.5g",DDSIP_node[DDSIP_bb->curnode]->dispnorm);
        fprintf (DDSIP_outfile, "  %-10.5g",DDSIP_node[DDSIP_bb->curnode]->dispnorm);
    }

    if (DDSIP_bb->bestvalue < DDSIP_infty && fabs (DDSIP_bb->bestbound) < DDSIP_infty)
    {
        printf (" %9.4g%%", rgap);
        fprintf (DDSIP_outfile, " %9.4g%%", rgap);
    }
    else
    {
        printf ("           ");
        fprintf (DDSIP_outfile, "           ");
    }

    DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
    time (&DDSIP_bb->cur_time);
    DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
    printf ("  %3dh %02d:%02.0f  %3dh %02d:%02.0f %8d %5d\n", wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs,DDSIP_node[DDSIP_bb->curnode]->father,DDSIP_node[DDSIP_bb->curnode]->depth);
    fprintf (DDSIP_outfile,"  %3dh %02d:%02.0f  %3dh %02d:%02.0f %8d %5d\n" ,wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs,DDSIP_node[DDSIP_bb->curnode]->father,DDSIP_node[DDSIP_bb->curnode]->depth);
}

//==========================================================================
// Function prints status of optimization
// in case of a new upper bound

void
DDSIP_PrintStateUB (int heur)
{
    double wall_secs, cpu_secs;
    int    wall_hrs, wall_mins, cpu_hrs, cpu_mins;

    printf ("*%6d  %6d        ", DDSIP_bb->curnode, DDSIP_bb->nonode);
    fprintf (DDSIP_outfile, "*%6d  %6d        ", DDSIP_bb->curnode, DDSIP_bb->nonode);

    if (DDSIP_bb->DDSIP_step == eev)
    {
        printf ("   exp. value sol ");
        fprintf (DDSIP_outfile, "   exp. value sol ");
    }
    else if (DDSIP_bb->DDSIP_step == adv)
    {
        fprintf (DDSIP_outfile, "   start point    ");
    }
    else
    {
        printf ("     Heuristic %2d ", heur);
        fprintf (DDSIP_outfile, "     Heuristic %2d ", heur);
    }

    // DDSIP_bb->heurval contains the objective value of the heuristic solution
    // DDSIP_bb->skip indicates whether the evaluation of an heuristic solution was skip for some reason
    printf ("  %-16.12g", DDSIP_bb->heurval);
    fprintf (DDSIP_outfile, "  %-16.12g", DDSIP_bb->heurval);

    if (fabs (DDSIP_bb->bestvalue) < DDSIP_infty)
    {
        printf ("  %-16.12g", DDSIP_bb->bestvalue);
        fprintf (DDSIP_outfile, "  %-16.12g", DDSIP_bb->bestvalue);
    }
    else
    {
        printf ("                  ");
        fprintf (DDSIP_outfile, "                  ");
    }

    printf ("                                               ");
    fprintf (DDSIP_outfile,"                                               ");

    DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
    time (&DDSIP_bb->cur_time);
    DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
    printf ("  %3dh %02d:%02.0f  %3dh %02d:%02.0f\n", wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs);
    fprintf (DDSIP_outfile,"  %3dh %02d:%02.0f  %3dh %02d:%02.0f\n", wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs);
}
