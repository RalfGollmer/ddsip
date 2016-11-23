/*  Authors:            Andreas M"arkert, Ralf Gollmer
	Copyright to:      University of Duisburg-Essen
    Language:          C
    Last modification: 12.03.2016
	Description: ,
	This file contains the implementation of the different heuristics to
	get a suggestion for an upper bound.

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

void DDSIP_RoundDown (double *);
void DDSIP_RoundUp (double *);
void DDSIP_RoundNear (double *);
void DDSIP_Frequent (void);
void DDSIP_Probable (void);
void DDSIP_CloseToAverage (double *);
void DDSIP_SmallValue (void);
void DDSIP_LargeValue (void);
void DDSIP_MinSum (void);
void DDSIP_MaxSum (void);
int  DDSIP_All (void);


//==========================================================================
// Rounding the average of integer components down
void
DDSIP_RoundDown (double *average)
{
    int j;

    for (j = 0; j < DDSIP_bb->firstvar; j++)
        if (DDSIP_bb->firsttype[j] == 'B' || DDSIP_bb->firsttype[j] == 'I' || DDSIP_bb->firsttype[j] == 'N') {
            ((DDSIP_bb->sug[DDSIP_param->nodelim + 2])->firstval)[j] = floor (average[j] + 0.00001);
        }
        else
            ((DDSIP_bb->sug[DDSIP_param->nodelim + 2])->firstval)[j] = average[j];
}

//==========================================================================
// Rounding the average of integer components up
void
DDSIP_RoundUp (double *average)
{
    int j;

    for (j = 0; j < DDSIP_bb->firstvar; j++)
        if (DDSIP_bb->firsttype[j] == 'B' || DDSIP_bb->firsttype[j] == 'I' || DDSIP_bb->firsttype[j] == 'N') {
            (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[j] = ceil (average[j]-0.00001);
        }
        else
            (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[j] = average[j];
}

//==========================================================================
// Rounding the average of integer components to nearest integer
void
DDSIP_RoundNear (double *average)
{
    int j;

    for (j = 0; j < DDSIP_bb->firstvar; j++)
    {
        if (DDSIP_bb->firsttype[j] == 'B' || DDSIP_bb->firsttype[j] == 'I' || DDSIP_bb->firsttype[j] == 'N') {
            (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[j] = floor (average[j] + 0.50001);
        }
        else
            (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[j] = average[j];
    }
}

//==========================================================================
// Most frequent scenario solution
void
DDSIP_Frequent (void)
{
    int i, j, k;
    k = 0;
    for (i = 1; i < DDSIP_param->scenarios; i++)
        if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][DDSIP_bb->firstvar] > (DDSIP_node[DDSIP_bb->curnode]->first_sol)[k][DDSIP_bb->firstvar])
            k = i;

    for (j = 0; j < DDSIP_bb->firstvar; j++)
        (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[j] = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[k][j];
}


//==========================================================================
// Most probable scenario solution
void
DDSIP_Probable (void)
{
    int i, j, k;
    double * prob;
    prob = (double *) DDSIP_Alloc(sizeof (double), DDSIP_param->scenarios, "prob(All)");

    //sum up probabilities of (equal) scenario solutions
    prob[0] += DDSIP_data->prob[0];
    for (i = 1; i < DDSIP_param->scenarios; i++)
    {
        k = i;
        if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][DDSIP_bb->firstvar] > 1)
        {
            for (j = 0; j < i; j++)
            {
                if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[i] == (DDSIP_node[DDSIP_bb->curnode]->first_sol)[j])
                {
                    k = j;
                    break;
                }
            }
        }
        prob[k] += DDSIP_data->prob[i];
    }
    // determine the solution with maximal probability
    k = 0;
    for (i = 1; i < DDSIP_param->scenarios; i++)
    {
        if (prob[i] > prob[k])
            k = i;
    }

    for (j = 0; j < DDSIP_bb->firstvar; j++)
        (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[j] = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[k][j];
    DDSIP_Free((void **) &prob);
}


//==========================================================================
// Scenario solution closest to average
void
DDSIP_CloseToAverage (double *average)
{
    int i, j, minind = 0;
    double tmp, mindev = DDSIP_infty, deviation;

    for (i = 0; i < DDSIP_param->scenarios; i++)
    {
        deviation = 0.0;
        for (j = 0; j < DDSIP_bb->firstvar; j++) {
            tmp = average[j] - (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][j];
            deviation = deviation + tmp*tmp;
        }
        if (deviation < mindev)
        {
            mindev = deviation;
            minind = i;
        }
    }
    for (j = 0; j < DDSIP_bb->firstvar; j++)
        (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[j] = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[minind][j];
}

//==========================================================================
// Scenario solution with smallest objective value
void
DDSIP_SmallValue (void)
{
    int j, i = 0;
    double tmp = (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[0];

    for (j = 1; j < DDSIP_param->scenarios; j++)
        if (tmp > (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j])
        {
            tmp = (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j];
            i = j;
        }
    for (j = 0; j < DDSIP_bb->firstvar; j++)
        (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[j] = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][j];
}

//==========================================================================
// Senario solution with largest objective value
void
DDSIP_LargeValue (void)
{
    int j, i = 0;
    double tmp = (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[0];

    for (j = 1; j < DDSIP_param->scenarios; j++)
        if (tmp < (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j])
        {
            tmp = (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j];
            i = j;
        }
    for (j = 0; j < DDSIP_bb->firstvar; j++)
        (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[j] = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][j];
}

//==========================================================================
// Senario solution with minimal sum of first-stage variables (originally to suit for EPS problem)
void
DDSIP_MinSum (void)
{
    int i, j, minind = 0;
    double min, tmp;

    tmp = 0.;
    if (DDSIP_bb->solstat[0])
        for (i = 0; i < DDSIP_bb->firstvar; i++)
            tmp += (DDSIP_node[DDSIP_bb->curnode]->first_sol)[0][i];

    min = tmp;

    for (j = 1; j < DDSIP_param->scenarios; j++)
    {
        if (DDSIP_bb->solstat[j])
        {
            tmp = 0;
            for (i = 0; i < DDSIP_bb->firstvar; i++)
                tmp += (DDSIP_node[DDSIP_bb->curnode]->first_sol)[j][i];

            if (tmp < min)
            {
                minind = j;
                min = tmp;
            }
        }
    }
    for (i = 0; i < DDSIP_bb->firstvar; i++)
        (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[i] = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[minind][i];
}

//==========================================================================
// Senario solution with maximal sum of first-stage variables (originally to suit for EPS problem)
void
DDSIP_MaxSum (void)
{
    int i, j, maxind = 0;
    double max, tmp;

    tmp = 0.;
    if (DDSIP_bb->solstat[0])
        for (i = 0; i < DDSIP_bb->firstvar; i++)
            tmp += (DDSIP_node[DDSIP_bb->curnode]->first_sol)[0][i];

    max = tmp;

    for (j = 1; j < DDSIP_param->scenarios; j++)
    {
        if (DDSIP_bb->solstat[j])
        {
            tmp = 0.;
            for (i = 0; i < DDSIP_bb->firstvar; i++)
                tmp += (DDSIP_node[DDSIP_bb->curnode]->first_sol)[j][i];

            if (tmp > max)
            {
                maxind = j;
                max = tmp;
            }
        }
    }
    for (i = 0; i < DDSIP_bb->firstvar; i++)
        (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[i] = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[maxind][i];
}

//==========================================================================
// All suggests successively all scenario problem solutions
int
DDSIP_All (void)
{
    int i, ii, j, status;
    int cnt;
    int *unind = (int *) DDSIP_Alloc (sizeof (int), DDSIP_param->scenarios, "unind(All)");
    sug_t *tmp;

    j = cnt = 0;
    for (i = DDSIP_param->scenarios-1; i > -1; i--)
    {
        if (!DDSIP_param->heuristic_auto)
        {
            if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][DDSIP_bb->firstvar + 2] < DDSIP_bb->curnode)
            {
                if (DDSIP_param->outlev > 49)
                    fprintf (DDSIP_bb->moreoutfile, "  Heuristic 12: solution of scenario %d is from father node\n", i);
                continue;
            }
        }
        for (j = i - 1; j > -1; j--)
        {
            if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[i] == (DDSIP_node[DDSIP_bb->curnode]->first_sol)[j])
            {
                if (DDSIP_param->outlev > 49)
                    fprintf (DDSIP_bb->moreoutfile, "  Heuristic 12: solution of scenario %d is the same as that of scenario %d\n", i, j);
                break;
            }
        }
        if (j==-1)
            unind[cnt++] = i;
    }

    if (!cnt)
    {
        if (DDSIP_param->outlev)
            fprintf(DDSIP_bb->moreoutfile, "  all scenario solutions did occur before.\n");
        DDSIP_bb->skip = -4;
        DDSIP_Free ((void **) &(unind));
        return 1;
    }

    for (i = cnt-1; i > 0; i--)
    {
        if (DDSIP_killsignal)
        {
            DDSIP_Free ((void **) &(unind));
            return 1;
        }
        if (!(tmp = DDSIP_bb->sug[DDSIP_param->nodelim + 2]))
        {
            // DDSIP_Allocate memory for heuristics suggestion if there is none yet
            tmp = (sug_t *)DDSIP_Alloc (sizeof (sug_t), 1, "sug[i](Heuristic)");
            tmp->firstval = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar, "sug[i]->firstval(Heuristic)");
            tmp->next = NULL;
            DDSIP_bb->sug[DDSIP_param->nodelim + 2] = tmp;
        }
        else
        {
            if (tmp->next != NULL)
                printf("\n\n FEHLER: !!!!!!!!!!!!!!!!!!!!! DDSIP_bb->sug[DDSIP_param->nodelim + 2]->next = %p\n\n",DDSIP_bb->sug[DDSIP_param->nodelim + 2]->next);
        }

        for (j = 0; j < DDSIP_bb->firstvar; j++)
            (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[j] = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[unind[i]][j];

        if (DDSIP_param->cpxubscr||DDSIP_param->outlev > 7)
            printf ("Heuristic 12: suggested first stage solution of scen %d\n", unind[i]+1);
        if(DDSIP_param->outlev)
            fprintf (DDSIP_bb->moreoutfile, "\nHeuristic 12: suggested first stage solution of scen %d", unind[i]+1);

        if (DDSIP_param->relax)
            for (ii = 0; ii < DDSIP_bb->firstvar; ii++)
                if (DDSIP_bb->firsttype[ii] == 'B' || DDSIP_bb->firsttype[ii] == 'I' || DDSIP_bb->firsttype[ii] == 'N')
                    (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[ii] = floor (((DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[ii] + 0.5));


        // In case of worst-case risk measure change value of aux var to worst_case_lb for heuristics
        if (abs(DDSIP_param->riskmod) == 4 && !DDSIP_bb->skip && !DDSIP_param->riskalg && !DDSIP_param->scalarization)
        {
            double worst_case_lb = -DDSIP_infty;
            for (ii=0; ii<DDSIP_param->scenarios; ii++)
            {
                worst_case_lb = DDSIP_Dmax((DDSIP_node[DDSIP_bb->curnode]->first_sol)[ii][DDSIP_bb->firstvar-1], worst_case_lb);
            }
            (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[DDSIP_bb->firstvar - 1] = worst_case_lb + DDSIP_param->accuracy;
        }

        if ((status = DDSIP_UpperBound ()))
        {
            if (status < 100000)
            {
                fprintf (stderr, "ERROR: Failed to perform UpperBound (All)\n");
                return status;
            }
            else if (DDSIP_param->interrupt_heur)
            {
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile, "\nGap reached, no further heuristics necessary.\n");
                break;
            }
        }
    }

// the last scenario
    if (cnt)
    {
        if (!(tmp = DDSIP_bb->sug[DDSIP_param->nodelim + 2]))
        {
            // DDSIP_Allocate memory for heuristics suggestion if there is none yet
            tmp = (sug_t *)DDSIP_Alloc (sizeof (sug_t), 1, "sug[i](Heuristic)");
            tmp->firstval = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar, "sug[i]->firstval(Heuristic)");
            tmp->next = NULL;
            DDSIP_bb->sug[DDSIP_param->nodelim + 2] = tmp;
        }
        else
        {
            if (tmp->next != NULL)
                printf("\n\n FEHLER: !!!!!!!!!!!!!!!!!!!!! DDSIP_bb->sug[DDSIP_param->nodelim + 2]->next = %p\n\n",DDSIP_bb->sug[DDSIP_param->nodelim + 2]->next);
        }
        for (j = 0; j < DDSIP_bb->firstvar; j++)
            (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[j] = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[0][j];

        if(DDSIP_param->outlev)
            fprintf (DDSIP_bb->moreoutfile, "\nHeuristic 12: suggested first stage solution of scen %d", unind[i]+1);
        if (DDSIP_param->cpxubscr||DDSIP_param->outlev > 7)
            printf ("Heuristic 12: suggested first stage solution of scen %d\n", unind[i]+1);
    }
    else
    {
        if (DDSIP_bb->sug[DDSIP_param->nodelim + 2])
            DDSIP_Free ((void **) &(DDSIP_bb->sug[DDSIP_param->nodelim + 2]));
        DDSIP_bb->sug[DDSIP_param->nodelim + 2] = NULL;
    }
    DDSIP_Free ((void **) &(unind));

    return 0;
}

//==========================================================================
// Function returns a suggestion in the first components of sug[DDSIP_param->nodelim + 2]
int
DDSIP_Heuristics (int *comb)
{
    int i, j, status = 0;
    sug_t *tmp;

    double *average = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar,
                      "average(Heuristic)");

    if (!(tmp = DDSIP_bb->sug[DDSIP_param->nodelim + 2]))
    {
        // DDSIP_Allocate memory for heuristics suggestion if there is none yet
        tmp = (sug_t *)DDSIP_Alloc (sizeof (sug_t), 1, "sug[i](Heuristic)");
        tmp->firstval = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar, "sug[i]->firstval(Heuristic)");
        for (i = 0; i < DDSIP_bb->firstvar; (tmp->firstval)[i++] = DDSIP_infty);
        tmp->next = NULL;
        DDSIP_bb->sug[DDSIP_param->nodelim + 2] = tmp;
    }
    else
    {
        for (i = 0; i < DDSIP_bb->firstvar; (tmp->firstval)[i++] = DDSIP_infty);
        if (tmp->next != NULL)
            printf("\n\n FEHLER: !!!!!!!!!!!!!!!!!!!!! DDSIP_bb->sug[DDSIP_param->nodelim + 2]->next = %p\n\n",DDSIP_bb->sug[DDSIP_param->nodelim + 2]->next);
        //tmp->next = NULL;
    }
    if (DDSIP_param->outlev)
    {
        fprintf (DDSIP_bb->moreoutfile, "\n----------------------\n");
        fprintf (DDSIP_bb->moreoutfile, "Invoking heuristic %d...", DDSIP_param->heuristic);
        if (DDSIP_param->cpxubscr||DDSIP_param->outlev > 10)
            printf ("Invoking heuristic %d...\n", DDSIP_param->heuristic);
    }
    // Calculate average of solutions of scenario problems
//if (DDSIP_param->outlev > 30)
//  fprintf (DDSIP_bb->moreoutfile, "\n__ 1st stage scen %d:", 0);
    for (j = 0; j < DDSIP_bb->firstvar; j++)
    {
//  if (DDSIP_param->outlev > 30)
//    fprintf (DDSIP_bb->moreoutfile, " %g,", (DDSIP_node[DDSIP_bb->curnode]->first_sol)[0][j]);
        average[j] = DDSIP_data->prob[0] * (DDSIP_node[DDSIP_bb->curnode]->first_sol)[0][j];
    }
    for (i = 1; i < DDSIP_param->scenarios; i++)
    {
//  if (DDSIP_param->outlev > 30)
//    fprintf (DDSIP_bb->moreoutfile, "\n__ 1st stage scen %d:", i);
        for (j = 0; j < DDSIP_bb->firstvar; j++)
        {
//    if (DDSIP_param->outlev > 30)
//	fprintf (DDSIP_bb->moreoutfile, " %g,", (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][j]);
            average[j] += DDSIP_data->prob[i] * (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][j];
        }
    }
//if (DDSIP_param->outlev > 30)
//  fprintf (DDSIP_bb->moreoutfile, "\n");

    switch (DDSIP_param->heuristic)
    {
    case 1:
        DDSIP_RoundDown (average);
        break;
    case 2:
        DDSIP_RoundUp (average);
        break;
    case 3:
        DDSIP_RoundNear (average);
        break;
    case 4:
        DDSIP_Frequent ();
        break;
    case 5:
        DDSIP_CloseToAverage (average);
        break;
    case 6:
        DDSIP_Probable();
        break;
    case 7:
        DDSIP_SmallValue ();
        break;
    case 8:
        DDSIP_LargeValue ();
        break;
    case 9:
        DDSIP_MinSum ();
        break;
    case 10:
        DDSIP_MaxSum ();
        break;
    case 11:			// Combined heuristic 3 and 5
        if (*comb == 3)
        {
            DDSIP_RoundNear (average);
            *comb = 5;
        }
        else
        {
            DDSIP_CloseToAverage (average);
            *comb = 3;
        }
        break;
    case 12:
        status = DDSIP_All ();
        break;
    case 13:			// Combined heuristic 4 and 5
        if (*comb == 5)
        {
            DDSIP_CloseToAverage (average);
            *comb = 4;
        }
        else
        {
            DDSIP_Frequent ();
            *comb = 5;
        }
        break;
    default:
        DDSIP_RoundNear (average);
        break;
    }

    DDSIP_Free ((void **) &(average));

    if (!status && DDSIP_bb->sug[DDSIP_param->nodelim + 2])
    {
        // round resulting values for integer first stage vars
        if (DDSIP_param->relax)
            for (i = 0; i < DDSIP_bb->firstvar; i++)
                if (DDSIP_bb->firsttype[i] == 'I' || DDSIP_bb->firsttype[i] == 'B' || DDSIP_bb->firsttype[i] == 'N')
                {
                    (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[i] = floor (((DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[i] + .5));
                }
        // Consistent ?
        for (i = 0; i < DDSIP_bb->firstvar; i++)
        {
            if ((DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval[i] - DDSIP_bb->uborg[i])/(fabs(DDSIP_bb->uborg[i])+ 1.) > DDSIP_param->accuracy)
            {
                printf ("   high suggestion for variable %d: sug=%20.18f ub=%20.18f, difference=%lg\n",i,
                        (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[i],DDSIP_bb->uborg[i],(DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[i]-DDSIP_bb->uborg[i]);
                if (DDSIP_param->outlev)
                {
                    fprintf (DDSIP_bb->moreoutfile,"   high suggestion for variable %d: sug=%20.18f ub=%20.18f, difference=%lg\n",i,
                             (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[i],DDSIP_bb->uborg[i],(DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[i]-DDSIP_bb->uborg[i]);
                }
                DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval[i] = DDSIP_bb->uborg[i];
            }
            if ((DDSIP_bb->lborg[i] - DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval[i])/(fabs(DDSIP_bb->lborg[i])+ 1.) > DDSIP_param->accuracy)
            {
                // in the root node the lower bound for the additional variable for worst case costs was updated
                if (DDSIP_bb->curnode || !(i == DDSIP_data->firstvar && (abs(DDSIP_param->riskmod) != 4 || abs(DDSIP_param->riskmod) != 5)))
                {
                    printf ("   low suggestion for variable %d: sug=%14.8f lb=%14.8f, difference=%lg\n",i,
                            (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[i],DDSIP_bb->lborg[i],(DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[i]-DDSIP_bb->lborg[i]);
                    if (DDSIP_param->outlev)
                    {
                        fprintf (DDSIP_bb->moreoutfile,"   low suggestion for variable %d: sug=%14.8f lb=%14.8f, difference=%lg\n",i,
                                 (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[i],DDSIP_bb->lborg[i],(DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[i]-DDSIP_bb->lborg[i]);
                    }
                }
                (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[i] = DDSIP_bb->lborg[i];
            }
        }

    }
    // In case of worst-case risk measure change value of aux var to worst_case_lb for heuristics
    if (abs(DDSIP_param->riskmod) == 4 && !DDSIP_bb->skip && !DDSIP_param->riskalg && !DDSIP_param->scalarization)
    {
        double worst_case_lb = -DDSIP_infty;
        for (i=0; i<DDSIP_param->scenarios; i++)
        {
            worst_case_lb = DDSIP_Dmax((DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][DDSIP_bb->firstvar-1], worst_case_lb);
        }
        (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[DDSIP_bb->firstvar - 1] = worst_case_lb + DDSIP_param->accuracy;
    }

    return 0;
}
