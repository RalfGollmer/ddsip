/*  Authors:           Andreas M"arkert, Ralf Gollmer
    Copyright to:      University of Duisburg-Essen
    Language:          C

    Description:
    This file contains the evaluation of the nodes in the
    branch-and-bound tree for the lower bound computation.

    License:
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

//#define DEBUG

#define SHIFT
#define CBHOTSTART
//#define CHECK_BOUNDS

#include <DDSIP.h>
#include <DDSIPconst.h>

static int DDSIP_SkipUB (void);
static int DDSIP_Warm (int);
static int DDSIP_LowerBoundWriteLp (int);


//==========================================================================
int
DDSIP_SkipUB (void)
{
    int i, status = 1;

    // There are some conditions which make a further proceeding unneccessary
    // If solstat=0 for all scenarios skip UB()
    for (i = 0; i < DDSIP_param->scenarios; i++)
        if (DDSIP_bb->solstat[i])
            status = 0;

    // Or if the node will be fathomed due to inferiority (bound > bestvalue)
    if (!DDSIP_bb->bestsol_in_curnode && DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bestvalue + DDSIP_bb->correct_bounding)
    {
        if (DDSIP_param->outlev > 20)
            fprintf (DDSIP_bb->moreoutfile, "########## skip=2 for DDSIP_bb->bestsol_in_curnode=%d, DDSIP_node[%d]->bound (%20.15g) - bestvalue = %.8g\n", DDSIP_bb->bestsol_in_curnode, DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bestvalue);
        return 2;
    }

    // Or if the node will be fathomed due to optimality (bound ~ bestvalue)
    //if (fabs(DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bestvalue) < DDSIP_param->accuracy)
    if (fabs(DDSIP_bb->bestvalue - DDSIP_node[DDSIP_bb->curnode]->bound)/(fabs(DDSIP_bb->bestvalue) + 3.e-16) < 2.e-14)
    {
        return 3;
    }

    // Or if bound=DDSIP_infty (this is different from the previous criteria if bestvalue=DDSIP_infty
    //  if (!(DDSIP_infty>DDSIP_node[DDSIP_bb->curnode]->bound))
    //    status = 1;

    return status;
} // DDSIP_SkipUB

//==========================================================================
// The component to branch on and the corresponding branch value are determined
// Maximal dispersion becomes dispnorm of current node
// Maximal dispersion = max_j (max_s x(js) - min_s x(js)), j in Firstvar, s in Scenarios, and
// x_s is the first-stage part of the solution of scenario problem s
int
DDSIP_GetBranchIndex (double *dispnorm)
{
    int i, j, k, cnt, below, above, diff, mindiff, maxdiff, cur_diff = 99999, maxord = -1;
    double hlb, hub, h, hbranch = 0., dist, dist_below = 0., dist_above = 0., maxdist = 0., factor = 1., cur_dist = DDSIP_infty;
    int *index = (int *) DDSIP_Alloc (sizeof (int), DDSIP_bb->firstvar, "index GetBI");
    // no branching on leaf
    if (DDSIP_node[DDSIP_bb->curnode]->leaf)
    {
        DDSIP_node[DDSIP_bb->curnode]->branchind = 0;
        DDSIP_node[DDSIP_bb->curnode]->branchval = 0.;
        return 0;
    }
    for (i=0; i<DDSIP_bb->firstvar; index[i++]=-100000);

    // Initialize branching index and dispersion norm of current node
    DDSIP_node[DDSIP_bb->curnode]->branchind = -1;
    h = -1.;
    mindiff = 1000000;
    maxdiff = -1;

    if (!DDSIP_param->riskalg && !DDSIP_param->scalarization)
    {
        // Some risk models (worst case, tvar) introduce an additional first-stage variable.
        // The parameter brancheta determines whether this variable (eta) will be branched.
        if ((abs (DDSIP_param->riskmod) == 4 || abs (DDSIP_param->riskmod) == 5) && !DDSIP_param->brancheta)
            dispnorm[DDSIP_bb->firstvar - 1] = -DDSIP_infty;

        if (DDSIP_param->cb && abs (DDSIP_param->riskmod) == 5)
        {
            if (DDSIP_param->cb < 0 && DDSIP_node[DDSIP_bb->curnode]->depth == 5)
                DDSIP_bb->local_bestdual[DDSIP_bb->dimdual] = DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_bb->dimdual] = DDSIP_Dmax(DDSIP_param->cbweight, 0.8*dispnorm[DDSIP_bb->firstvar - 1]);
            else if (DDSIP_param->cb > 0 && DDSIP_node[DDSIP_bb->curnode]->depth == 0)
                DDSIP_bb->local_bestdual[DDSIP_bb->dimdual] = DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_bb->dimdual] = DDSIP_Dmax(DDSIP_param->cbweight, 0.05*dispnorm[DDSIP_bb->firstvar - 1]);
        }

        // due to the changes in the bounds: do not branch in root node on additional variable for WC
        if (!(DDSIP_bb->curnode) && abs(DDSIP_param->riskmod) == 4)
        {
            dispnorm[DDSIP_bb->firstvar - 1] = 0.;
        }
        else if (DDSIP_param->brancheta && abs(DDSIP_param->riskmod) == 5)
        {
            // do branch in first nodes on additional variable for TVaR
            if (((DDSIP_bb->curnode<32) && (dispnorm[DDSIP_bb->firstvar - 1] > DDSIP_Dmax(0.99,0.01*(DDSIP_bb->uborg[DDSIP_bb->firstvar - 1] - DDSIP_bb->lborg[DDSIP_bb->firstvar - 1])))) ||
                    // and again at levels 7, 8, 10,11,20,21 etc.
                    ((((DDSIP_node[DDSIP_bb->curnode]->depth > 6 && DDSIP_node[DDSIP_bb->curnode]->depth < 9) || (DDSIP_node[DDSIP_bb->curnode]->depth%10) < 2)) && (dispnorm[DDSIP_bb->firstvar - 1] > DDSIP_Dmax(0.05,0.0001*(DDSIP_bb->uborg[DDSIP_bb->firstvar - 1] - DDSIP_bb->lborg[DDSIP_bb->firstvar - 1])))))
            {
                DDSIP_node[DDSIP_bb->curnode]->branchind = DDSIP_bb->firstvar - 1;
                hub = DDSIP_bb->uborg[DDSIP_node[DDSIP_bb->curnode]->branchind];
                hlb = DDSIP_bb->lborg[DDSIP_node[DDSIP_bb->curnode]->branchind];
                i = DDSIP_bb->curnode;
                while (i > 0)
                {
                    if (DDSIP_node[i]->neoind == DDSIP_node[DDSIP_bb->curnode]->branchind)
                    {
                        hlb = DDSIP_Dmax (DDSIP_node[i]->neolb, hlb);
                        hub = DDSIP_Dmin (DDSIP_node[i]->neoub, hub);
                    }
                    i = DDSIP_node[i]->father;
                }
                h = hub;
                hbranch = hlb;
                for (i = 0; i < DDSIP_param->scenarios; i++)
                {
                    if (hbranch  < (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][DDSIP_node[DDSIP_bb->curnode]->branchind])
                        hbranch  = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][DDSIP_node[DDSIP_bb->curnode]->branchind];
                    if (h        > (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][DDSIP_node[DDSIP_bb->curnode]->branchind])
                        h        = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][DDSIP_node[DDSIP_bb->curnode]->branchind];
                }
                h = DDSIP_Dmax (h, hlb);
                hbranch = DDSIP_Dmin (hbranch + 0.1, hub);
                //DDSIP_node[DDSIP_bb->curnode]->branchval = 0.5 * (h + hbranch);
                DDSIP_node[DDSIP_bb->curnode]->branchval = 0.48 * h + 0.52 * hbranch;
                maxdiff = mindiff = 0;
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile,"-- forcing branch on aux. variable: node %d, depth %d, branchval=%18.16g, (min= %.16g, max= %.16g)\n",DDSIP_bb->curnode,DDSIP_node[DDSIP_bb->curnode]->depth,DDSIP_node[DDSIP_bb->curnode]->branchval,h,hbranch);
                goto READY;
            }
            else
            {
                if ((DDSIP_bb->bestvalue == DDSIP_infty && DDSIP_bb->curnode<DDSIP_data->firstvar*10))
                {
                    h = -DDSIP_infty;
                    for(i=0; i<DDSIP_data->firstvar; i++)
                    {
                        if (h < dispnorm[i])
                            h = dispnorm[i];
                    }
                    if (h > 2.*DDSIP_param->nulldisp)
                        dispnorm[DDSIP_bb->firstvar - 1] = 0.;
                }
            }
        }
        // do branch in depth 3 on additional variable for WC plus depth equal multiples of 5
        else if (DDSIP_param->brancheta && (abs(DDSIP_param->riskmod) == 4)  && (DDSIP_node[DDSIP_bb->curnode]->depth == 3 || !(DDSIP_node[DDSIP_bb->curnode]->depth%5)))
        {
            DDSIP_node[DDSIP_bb->curnode]->branchind = DDSIP_bb->firstvar - 1;
            hub = - DDSIP_infty;
            hlb = DDSIP_infty;
            for (k = 0; k < DDSIP_param->scenarios; k++)
            {
                hlb = DDSIP_Dmin (hlb, (DDSIP_node[DDSIP_bb->curnode]->first_sol)[k][DDSIP_node[DDSIP_bb->curnode]->branchind]);
                hub = DDSIP_Dmax (hub, (DDSIP_node[DDSIP_bb->curnode]->first_sol)[k][DDSIP_node[DDSIP_bb->curnode]->branchind]);
            }
            if (hub - hlb > 1e-7*fabs(hub))
            {
                hbranch = hub - 1.0e-12*fabs(hub);
                if (hbranch < hlb + DDSIP_param->brancheps)
                    hbranch = 0.5* (hlb+hub);
                DDSIP_node[DDSIP_bb->curnode]->branchval = hbranch;
                maxdiff = mindiff = 0;
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile,"-- forcing branch on aux. variable: node %d, depth %d, branchval=%18.16g, (min= %.16g, max= %.16g)\n",DDSIP_bb->curnode,DDSIP_node[DDSIP_bb->curnode]->depth,DDSIP_node[DDSIP_bb->curnode]->branchval,hlb,hub);
                goto READY;
            }
        }
    }

    // We do not branch variables which have been fixed in the upper tree.
    // Continuous variables are considered fix if ub-lb <= DDSIP_param->brancheps
    j = DDSIP_bb->curnode;
    while (j > 0)
    {
        if ((DDSIP_bb->firsttype[DDSIP_node[j]->neoind] == 'B' && DDSIP_node[j]->neoub == DDSIP_node[j]->neolb)
                || (DDSIP_bb->firsttype[DDSIP_node[j]->neoind] == 'I' && DDSIP_node[j]->neoub == DDSIP_node[j]->neolb)
                || (DDSIP_bb->firsttype[DDSIP_node[j]->neoind] == 'N' && DDSIP_node[j]->neoub == DDSIP_node[j]->neolb)
                || (DDSIP_bb->firsttype[DDSIP_node[j]->neoind] == 'C' && !(DDSIP_node[j]->neoub - DDSIP_node[j]->neolb > DDSIP_param->brancheps))
                || (DDSIP_bb->firsttype[DDSIP_node[j]->neoind] == 'S' && !(DDSIP_node[j]->neoub - DDSIP_node[j]->neolb > DDSIP_param->brancheps)))
            dispnorm[DDSIP_node[j]->neoind] = -DDSIP_infty;
        j = DDSIP_node[j]->father;
    }

    h = -1.;
    mindiff = 1000000;
    maxdiff = -1;
    // Only some variables are candidates
    cnt = 0;
    if (DDSIP_param->boundstrat)
    {
        if ((DDSIP_bb->curnode%20)<9)
            factor = 0.999;
        else
            factor = 0.95;
    }
    else
    {
        factor = 0.9999;
    }
    maxord = -1;
    // With branching order
    if (DDSIP_param->order)
    {
        j = 0;
        // Get maximal branching order among variables with sufficient dispersion norm
        while (j < DDSIP_bb->firstvar)
        {
            if (dispnorm[j] > DDSIP_Dmax(DDSIP_param->nulldisp, (exp(-0.5*DDSIP_bb->order[j])+0.5)*DDSIP_node[DDSIP_bb->curnode]->dispnorm))
                maxord = DDSIP_Imax (maxord, DDSIP_bb->order[j]);
            j++;
        }
        // Index holds those variables with order equal to maxord and with sufficient dispersion norm
        if (maxord > -1)
        {
            h = 0.;
            for (j = 0; j < DDSIP_bb->firstvar; j++)
            {
                if (DDSIP_bb->order[j] == maxord)
                {
                    h = DDSIP_Dmax (h, dispnorm[j]);
                }
            }
            if ((h*=factor) > 0.)
            {
                for (j = 0; j < DDSIP_bb->firstvar; j++)
                {
                    if (DDSIP_bb->order[j] == maxord && dispnorm[j] > h)
                    {
                        index[cnt++] = j;
                    }
                }
            }
        }
        if (DDSIP_param->outlev > 4 && maxord > 1)
            fprintf (DDSIP_bb->moreoutfile,"\tGetBranchIndex: max.order = %d, count: %d\n", maxord, cnt);
    }
    // Without branching order
    if (!cnt || maxord == 1)
    {
        cnt = 0;
        h = 0.;
        if (DDSIP_param->intfirst)
        {
            for (j = 0; j < DDSIP_bb->firstvar; j++)
            {
                if (dispnorm[j] > DDSIP_param->nulldisp && (DDSIP_bb->firsttype[j] == 'B' || DDSIP_bb->firsttype[j] == 'I' || DDSIP_bb->firsttype[j] == 'N'))
                {
                    h = DDSIP_Dmax (h,dispnorm[j]);
                }
            }
            if (h >= DDSIP_param->nulldisp)
            {
                h*=factor;
                for (j = 0; j < DDSIP_bb->firstvar; j++)
                {
                    if (dispnorm[j] > h && (DDSIP_bb->firsttype[j] == 'B' || DDSIP_bb->firsttype[j] == 'I' || DDSIP_bb->firsttype[j] == 'N'))
                    {
                        index[cnt++] = j;
                    }
                }
            }
            if (cnt == 0)          /* no integer variables to branch on - include all vars */
            {
                h = 0.;
                for (j = 0; j < DDSIP_bb->firstvar; j++)
                {
                    if (dispnorm[j] > DDSIP_param->nulldisp)
                    {
                        h = DDSIP_Dmax (h,dispnorm[j]);
                    }
                }
                if (h >= DDSIP_param->nulldisp)
                {
                    h*=factor;
                    for (j = 0; j < DDSIP_bb->firstvar; j++)
                    {
                        if (dispnorm[j] > h)
                        {
                            index[cnt++] = j;
                        }
                    }
                }
            }
        }
        else
        {
            // Index holds those variables with sufficient dispersion norm
            h = 0.;
            for (j = 0; j < DDSIP_bb->firstvar; j++)
            {
                if (dispnorm[j] > DDSIP_param->nulldisp)
                {
                    h = DDSIP_Dmax (h, dispnorm[j]);
                }
            }
            if (h >= DDSIP_param->nulldisp)
            {
                h*=factor;
                for (j = 0; j < DDSIP_bb->firstvar; j++)
                {
                    if (dispnorm[j] > h)
                    {
                        index[cnt++] = j;
                    }
                }
            }
        }
    }

    //
    if (DDSIP_param->outlev > 4)
        fprintf (DDSIP_bb->moreoutfile,"\tGetBranchIndex: threshold= %g, number of candidates = %d\n",h,cnt);
    //

    // Among the variables in 'index': Get variable with maximal dispersion norm most equally dividing the scenario sol. values
    // (which means the minimal difference in the numbers of scenarios above or below the branching value)
    if (cnt>1)
    {
        mindiff =  DDSIP_bigint;
        maxdiff = -DDSIP_bigint;
        maxdist = -1.;
        DDSIP_node[DDSIP_bb->curnode]->branchind = index[0];
        for (j = 0; j < cnt; j++)
        {
            // Set tentative branching value
            // determine minimal and maximal value of the variable
            hlb =  DDSIP_infty;
            hub = -DDSIP_infty;
            for (k = 0; k < DDSIP_param->scenarios; k++)
            {
                hlb = DDSIP_Dmin (hlb, (DDSIP_node[DDSIP_bb->curnode]->first_sol)[k][index[j]]);
                hub = DDSIP_Dmax (hub, (DDSIP_node[DDSIP_bb->curnode]->first_sol)[k][index[j]]);
            }
            if (DDSIP_param->branchstrat == 0)
            {
                // branchstrat = 0: take the middle between maximal and minimal value
                hbranch = 0.5* (hlb+hub);
            }
            else
            {
                if ((abs(DDSIP_param->riskmod) == 4) && (index[j] == DDSIP_bb->firstvar - 1))
                {
                    hbranch = hub - DDSIP_Dmax(DDSIP_param->brancheps,1.0e-7)*fabs(hub);
                    if (hbranch < hlb + DDSIP_param->brancheps)
                        hbranch = 0.5* (hlb+hub);
                }
                else
                {
                    // branchstrat = 1: the expected value
                    hbranch = DDSIP_data->prob[0] * (DDSIP_node[DDSIP_bb->curnode]->first_sol)[0][index[j]];
                    for (i = 1; i < DDSIP_param->scenarios; i++)
                        hbranch += DDSIP_data->prob[i] * (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][index[j]];
                    if (DDSIP_param->branchstrat == 2)
                    {
                        // branchstrat = 2: take a point between the expected value and the middle
                        hbranch = 0.25*hbranch + 0.375* (hlb+hub);
                    }
                }
            }
            // now count the scenario sol. values below and above the branch value and measure the distances
            dist_below = dist_above = dist = 0.;
            below = above = 0;
            for (i = 0; i < DDSIP_param->scenarios; i++)
            {
                if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][index[j]] <= hbranch)
                {
                    below++;
                    hlb =  hbranch - (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][index[j]];
                    dist_below += DDSIP_data->prob[i] * hlb*hlb;
                }
                else if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][index[j]] > hbranch)
                {
                    above++;
                    hlb = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][index[j]] - hbranch;
                    dist_above += DDSIP_data->prob[i] * hlb*hlb;
                }
            }
            diff = abs (below-above);
            dist = DDSIP_Dmin (dist_below,dist_above) + 5e-1*(dist_below+dist_above);
            if (DDSIP_bb->firsttype[index[j]] == 'B' || DDSIP_bb->firsttype[index[j]] == 'I' || DDSIP_bb->firsttype[index[j]] == 'N')
            {
                // for integers: make hbranch end in .5
                hbranch=floor (hbranch + 1.e-6) + 0.5;
            }
            if (!((abs(DDSIP_param->riskmod) == 4) && (index[j] == DDSIP_bb->firstvar - 1)) &&
                    ((DDSIP_param->equalbranch == 1) ||
                     (!DDSIP_param->equalbranch && ((DDSIP_bb->curnode < 11) || ((DDSIP_bb->curnode > 25) && (DDSIP_bb->curnode%10 >= 6))))
                    )
               )
            {
                if ((diff <= mindiff && dist > maxdist) || (diff < mindiff && dist > 0.8*maxdist) || (diff <= mindiff && dispnorm[index[j]] > dispnorm[DDSIP_node[DDSIP_bb->curnode]->branchind] && dist > 0.9*maxdist) || (diff < mindiff+3 && dispnorm[index[j]] > factor*dispnorm[DDSIP_node[DDSIP_bb->curnode]->branchind] && dist > maxdist + 1e-1))
                {
                    if (diff < mindiff)
                        maxdist = dist;
                    else
                        maxdist = DDSIP_Dmax (maxdist,dist);
                    mindiff = DDSIP_Imin (mindiff,diff);
                    DDSIP_node[DDSIP_bb->curnode]->branchind = index[j];
                    DDSIP_node[DDSIP_bb->curnode]->branchval = hbranch;
                    cur_diff = diff;
                    cur_dist = dist;
                }
                //
                if (DDSIP_param->outlev > 21)
                    fprintf (DDSIP_bb->moreoutfile,"\tGetBranchIndex 1: index  %5d disp. %15.12g,  below: %3d above: %3d diff %3d  mindiff %3d dist %12.7g branchindex %d\n",DDSIP_bb->firstindex[index[j]],dispnorm[index[j]],below,above,diff,mindiff,dist,DDSIP_bb->firstindex[DDSIP_node[DDSIP_bb->curnode]->branchind]);
                //
            }
            else
            {
                if ((above && below) && (diff > maxdiff || (diff >= maxdiff && dist > maxdist) || (diff == maxdiff && dispnorm[index[j]] > dispnorm[DDSIP_node[DDSIP_bb->curnode]->branchind]) || (diff > maxdiff-4 && dist > 1.2*maxdist)))
                {
                    if (diff > maxdiff)
                        maxdist = dist;
                    else
                        maxdist = DDSIP_Dmax (maxdist,dist);
                    maxdiff = DDSIP_Imax (maxdiff,diff);
                    DDSIP_node[DDSIP_bb->curnode]->branchind = index[j];
                    DDSIP_node[DDSIP_bb->curnode]->branchval = hbranch;
                    cur_diff = diff;
                    cur_dist = dist;
                }
                //
                if (DDSIP_param->outlev > 21)
                    fprintf (DDSIP_bb->moreoutfile,"\tGetBranchIndex 2: index  %5d disp. %15.12g  below %3d above %3d diff %3d  maxdiff %3d dist %12.7g branchindex %d\n",DDSIP_bb->firstindex[index[j]],dispnorm[index[j]],below,above,diff,maxdiff,dist,DDSIP_bb->firstindex[DDSIP_node[DDSIP_bb->curnode]->branchind]);
                //
            }
        }
        if (DDSIP_param->outlev > 19)
            fprintf (DDSIP_bb->moreoutfile,"\tGetBranchIndex: index %d, diff %3d  dist %8.5g, maxdist %8.5g mindiff %d maxdiff %3d\n",DDSIP_bb->firstindex[DDSIP_node[DDSIP_bb->curnode]->branchind],cur_diff,cur_dist,maxdist,mindiff,maxdiff);
    }
    else if (cnt) // i.e. cnt == 1
        DDSIP_node[DDSIP_bb->curnode]->branchind = index[cnt-1];
    // If the branching index could not be set according to the dispnorm criterion
    // we branch the variable with the largest deviation among lower and upper bound
    // This may happen when no scenario problem yields a solution or when the maximal dispersion is 0
    if (DDSIP_node[DDSIP_bb->curnode]->branchind == -1)
    {
        double lb, ub;
        double *dev = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar, "dev (GetBI)");
        for (j = 0; j < DDSIP_bb->firstvar; j++)
        {
            ub = DDSIP_bb->uborg[j];
            lb = DDSIP_bb->lborg[j];
            cnt = DDSIP_bb->curnode;
            while (cnt > 0)
            {
                if (DDSIP_node[cnt]->neoind == j)
                {
                    ub = DDSIP_Dmin (DDSIP_node[cnt]->neoub, ub);
                    lb = DDSIP_Dmax (DDSIP_node[cnt]->neolb, lb);
                }
                cnt = DDSIP_node[cnt]->father;
            }
            dev[j] = ub - lb;
        }
        ub = dev[0];
        DDSIP_node[DDSIP_bb->curnode]->branchind = 0;
        for (j = 1; j < DDSIP_bb->firstvar; j++)
            if (dev[j] > ub)
            {
                ub = dev[j];
                DDSIP_node[DDSIP_bb->curnode]->branchind = j;
            }
        DDSIP_Free ((void **) &(dev));
    }

    // Set branching value
    if (!DDSIP_param->branchstrat|| ((abs(DDSIP_param->riskmod) == 4) && (DDSIP_node[DDSIP_bb->curnode]->branchind == DDSIP_bb->firstvar - 1)))
    {
        // the center of the values ((min+max)/2)
        hlb =  DDSIP_infty;
        hub = -DDSIP_infty;
        for (j = 0; j < DDSIP_param->scenarios; j++)
        {
            hlb = DDSIP_Dmin (hlb, (DDSIP_node[DDSIP_bb->curnode]->first_sol)[j][DDSIP_node[DDSIP_bb->curnode]->branchind]);
            hub = DDSIP_Dmax (hub, (DDSIP_node[DDSIP_bb->curnode]->first_sol)[j][DDSIP_node[DDSIP_bb->curnode]->branchind]);
        }
        if ((abs(DDSIP_param->riskmod) == 4) && (DDSIP_node[DDSIP_bb->curnode]->branchind == DDSIP_bb->firstvar - 1))
        {
            // DDSIP_node[DDSIP_bb->curnode]->branchval = 0.999*hub + 0.001*hlb;
            DDSIP_node[DDSIP_bb->curnode]->branchval = DDSIP_Dmax(hub - DDSIP_param->accuracy*fabs(hub), 0.5 * (hlb + hub));
        }
        else
            DDSIP_node[DDSIP_bb->curnode]->branchval = (0.5 * (hlb + hub));
    }
    else
    {
        // the expected value
        DDSIP_node[DDSIP_bb->curnode]->branchval = DDSIP_data->prob[0] * (DDSIP_node[DDSIP_bb->curnode]->first_sol)[0][DDSIP_node[DDSIP_bb->curnode]->branchind];
        for (j = 1; j < DDSIP_param->scenarios; j++)
            DDSIP_node[DDSIP_bb->curnode]->branchval += DDSIP_data->prob[j] * (DDSIP_node[DDSIP_bb->curnode]->first_sol)[j][DDSIP_node[DDSIP_bb->curnode]->branchind];
        if (DDSIP_param->branchstrat == 2)
        {
            // shift the expected value in direction of the center of the values ((min+max)/2)
            hlb =  DDSIP_infty;
            hub = -DDSIP_infty;
            for (j = 0; j < DDSIP_param->scenarios; j++)
            {
                hlb = DDSIP_Dmin (hlb, (DDSIP_node[DDSIP_bb->curnode]->first_sol)[j][DDSIP_node[DDSIP_bb->curnode]->branchind]);
                hub = DDSIP_Dmax (hub, (DDSIP_node[DDSIP_bb->curnode]->first_sol)[j][DDSIP_node[DDSIP_bb->curnode]->branchind]);
            }
            if (DDSIP_node[DDSIP_bb->curnode]->dispnorm < 1e-4)
            {
                // for continuous variables with difference small enough: take the middle of the values
                DDSIP_node[DDSIP_bb->curnode]->branchval = 0.5 * (hlb + hub);
            }
            else
                DDSIP_node[DDSIP_bb->curnode]->branchval = 0.25 * DDSIP_node[DDSIP_bb->curnode]->branchval + 0.375 * (hlb+hub);
        }
    }
    if (DDSIP_bb->firsttype[DDSIP_node[DDSIP_bb->curnode]->branchind] == 'B' || DDSIP_bb->firsttype[DDSIP_node[DDSIP_bb->curnode]->branchind] == 'I' || DDSIP_bb->firsttype[DDSIP_node[DDSIP_bb->curnode]->branchind] == 'N')
    {
        // for integers: make it end in .5
        DDSIP_node[DDSIP_bb->curnode]->branchval = floor (DDSIP_node[DDSIP_bb->curnode]->branchval + 1.e-6) + 0.5;
    }
    // Don't cut out the individual optimal solutions
    hlb = -DDSIP_infty;
    hub =  DDSIP_infty;
    for (i = 1; i < DDSIP_param->scenarios; i++)
    {
        if ((h = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][DDSIP_node[DDSIP_bb->curnode]->branchind]) <= DDSIP_node[DDSIP_bb->curnode]->branchval)
            hlb = DDSIP_Dmax (hlb,h);
        else
            hub = DDSIP_Dmin (hub,h);
    }
    if (hub < DDSIP_node[DDSIP_bb->curnode]->branchval + DDSIP_param->brancheps && hlb < DDSIP_node[DDSIP_bb->curnode]->branchval - DDSIP_param->brancheps - 1.e-14)
    {
#ifdef DEBUG
        if (DDSIP_param->outlev)
        {
            fprintf (DDSIP_bb->moreoutfile,"-------->shifting branching value from %22.16g (hlb=%22.16g, hub=%22.16g)",DDSIP_node[DDSIP_bb->curnode]->branchval,hlb,hub);
        }
#endif
        DDSIP_node[DDSIP_bb->curnode]->branchval = hub - DDSIP_param->brancheps - 1.e-14;
#ifdef DEBUG
        if (DDSIP_param->outlev)
        {
            fprintf (DDSIP_bb->moreoutfile,"  to %22.16g\n",DDSIP_node[DDSIP_bb->curnode]->branchval);
        }
#endif
    }
    // special case: branching variable is worst-case cost
    if (!DDSIP_param->riskalg && !DDSIP_param->scalarization && abs(DDSIP_param->riskmod) == 4 && DDSIP_node[DDSIP_bb->curnode]->branchind == DDSIP_bb->firstvar-1)
    {
        hub = DDSIP_bb->uborg[DDSIP_node[DDSIP_bb->curnode]->branchind];
        hlb = DDSIP_bb->lborg[DDSIP_node[DDSIP_bb->curnode]->branchind];
        i = DDSIP_bb->curnode;
        while (i > 0)
        {
            if (DDSIP_node[i]->neoind == DDSIP_node[DDSIP_bb->curnode]->branchind)
            {
                hlb = DDSIP_Dmax (DDSIP_node[i]->neolb, hlb);
                hub = DDSIP_Dmin (DDSIP_node[i]->neoub, hub);
            }
            i = DDSIP_node[i]->father;
        }
        for (j = 1; j < DDSIP_param->scenarios; j++)
        {
            if ((DDSIP_node[DDSIP_bb->curnode]->branchval < ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[j][DDSIP_bb->firstvar-1])) && ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[j][DDSIP_bb->firstvar-1] < hub - DDSIP_param->brancheps*0.2))
                DDSIP_node[DDSIP_bb->curnode]->branchval = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[j][DDSIP_node[DDSIP_bb->curnode]->branchind];
        }
        if (DDSIP_node[DDSIP_bb->curnode]->neoind == DDSIP_bb->firstvar-1)
        {
            if (DDSIP_node[DDSIP_bb->curnode]->neoub - DDSIP_node[DDSIP_bb->curnode]->neolb > 2.*DDSIP_param->brancheps)
            {
                DDSIP_node[DDSIP_bb->curnode]->branchval = DDSIP_Dmax (DDSIP_node[DDSIP_bb->curnode]->branchval,DDSIP_node[DDSIP_bb->curnode]->neolb + DDSIP_param->brancheps);
                DDSIP_node[DDSIP_bb->curnode]->branchval = DDSIP_Dmin (DDSIP_node[DDSIP_bb->curnode]->branchval,DDSIP_node[DDSIP_bb->curnode]->neoub - DDSIP_param->brancheps);
            }
            else
                DDSIP_node[DDSIP_bb->curnode]->branchval = 0.5* (DDSIP_node[DDSIP_bb->curnode]->neoub + DDSIP_node[DDSIP_bb->curnode]->neolb);
        }
    }

READY:
    if (DDSIP_param->outlev > 2)
    {
        hlb =  DDSIP_infty;
        hub = -DDSIP_infty;
        for (i = 0; i < DDSIP_param->scenarios; i++)
        {
            h = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][DDSIP_node[DDSIP_bb->curnode]->branchind];
            hlb = DDSIP_Dmin (hlb,h);
            hub = DDSIP_Dmax (hub,h);
        }
        fprintf (DDSIP_bb->moreoutfile,"\tchosen branch index: %d, branchval= %18.16g  (min= %.16g, max= %.16g) dispersion= %g\n",DDSIP_bb->firstindex[DDSIP_node[DDSIP_bb->curnode]->branchind],DDSIP_node[DDSIP_bb->curnode]->branchval,hlb,hub,dispnorm[DDSIP_node[DDSIP_bb->curnode]->branchind]);
    }

    DDSIP_Free ((void **) &(index));

    return 0;
} // DDSIP_GetBranchIndex

//==========================================================================
// Warm starts
int
DDSIP_Warm (int iscen)
{
    int scen, i, j, k, kf, status, jt, already_there;
    static int added, inserted;

    scen = DDSIP_bb->lb_scen_order[iscen];
    added = inserted = 0;

    if ((k = CPXgetnummipstarts(DDSIP_env, DDSIP_lp)) > 2)
    {
        status    = CPXdelmipstarts (DDSIP_env, DDSIP_lp, 1, k-2);
    }
    else if (k == 2)
    {
        status    = CPXdelmipstarts (DDSIP_env, DDSIP_lp, 1, 1);
    }

    if (DDSIP_param->hot && DDSIP_bb->DDSIP_step == dual && DDSIP_bb->dualitcnt)
    {
        DDSIP_bb->beg[0]=0;
        DDSIP_bb->effort[0]=3;

        if (memcmp(DDSIP_bb->boundIncrease_intsolvals[scen], DDSIP_bb->intsolvals[scen], DDSIP_bb->total_int*sizeof(double)))
        {
            // add solution of previous bound increase
            sprintf (DDSIP_bb->Names[0],"BoundIncr_%d",scen+1);
#ifdef DEBUG
            if(DDSIP_param->outlev > 90)
            {
                fprintf(DDSIP_bb->moreoutfile,"     ### MIP starts (in CB): %s  \n",DDSIP_bb->Names[0]);
                for (j = 0; j < DDSIP_bb->total_int; j++)
                {
                    fprintf(DDSIP_bb->moreoutfile,"%g, ", DDSIP_bb->boundIncrease_intsolvals[scen][j]);
                    if (!((j+1)%20))
                        fprintf(DDSIP_bb->moreoutfile,"\n");
                }
                fprintf(DDSIP_bb->moreoutfile,"\n");
            }
#endif
            // Copy starting informations to problem
            status = CPXaddmipstarts (DDSIP_env, DDSIP_lp, 1, DDSIP_bb->total_int, DDSIP_bb->beg, DDSIP_bb->intind, DDSIP_bb->boundIncrease_intsolvals[scen], DDSIP_bb->effort, DDSIP_bb->Names);
            if (status)
            {
                fprintf (stderr, "ERROR: Failed to copy mip start infos (Warm)\n");
                fprintf (DDSIP_outfile, "ERROR: Failed to copy mip start infos (Warm)\n");
                return status;
            }
#ifdef DEBUG
            else if (DDSIP_param->outlev > 25)
                fprintf (DDSIP_bb->moreoutfile, "##  Copied mip start info %s of previous bound increase (Warm)\n", DDSIP_bb->Names[0]);
#endif
            added++;
        }
        // add solution of previous iteration
        sprintf (DDSIP_bb->Names[0],"Prev_Iter_%d",scen+1);
        DDSIP_bb->effort[0]=3;
        // print debugging info
#ifdef DEBUG
        if(DDSIP_param->outlev > 90)
        {
            fprintf(DDSIP_bb->moreoutfile,"     ### MIP starts (in CB): %s  \n",DDSIP_bb->Names[0]);
            for (j = 0; j < DDSIP_bb->total_int; j++)
            {
                fprintf(DDSIP_bb->moreoutfile,"%g, ", DDSIP_bb->intsolvals[scen][j]);
                if (!((j+1)%20))
                    fprintf(DDSIP_bb->moreoutfile,"\n");
            }
            fprintf(DDSIP_bb->moreoutfile,"\n");
        }
#endif
        // Copy starting informations to problem
        status = CPXaddmipstarts (DDSIP_env, DDSIP_lp, 1, DDSIP_bb->total_int, DDSIP_bb->beg, DDSIP_bb->intind, DDSIP_bb->intsolvals[scen], DDSIP_bb->effort, DDSIP_bb->Names);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to copy mip start infos (Warm)\n");
            fprintf (DDSIP_outfile, "ERROR: Failed to copy mip start infos (Warm)\n");
            return status;
        }
#ifdef DEBUG
        else if (DDSIP_param->outlev > 30)
            fprintf (DDSIP_bb->moreoutfile, "##  Copied mip start info %s of previous iteration (Warm)\n", DDSIP_bb->Names[0]);
#endif
        added++;
//        if (memcmp(DDSIP_bb->boundIncrease_intsolvals[scen], DDSIP_bb->intsolvals[scen], DDSIP_bb->total_int*sizeof(double)))
//        {
//            // add solution of previous bound increase
//            sprintf (DDSIP_bb->Names[0],"BoundIncr_%d",scen+1);
//#ifdef DEBUG
//            if(DDSIP_param->outlev > 90)
//            {
//                fprintf(DDSIP_bb->moreoutfile,"     ### MIP starts (in CB): %s  \n",DDSIP_bb->Names[0]);
//                for (j = 0; j < DDSIP_bb->total_int; j++)
//                {
//                    fprintf(DDSIP_bb->moreoutfile,"%g, ", DDSIP_bb->boundIncrease_intsolvals[scen][j]);
//                    if (!((j+1)%20))
//                        fprintf(DDSIP_bb->moreoutfile,"\n");
//                }
//                fprintf(DDSIP_bb->moreoutfile,"\n");
//            }
//#endif
//            // Copy starting informations to problem
//            status = CPXaddmipstarts (DDSIP_env, DDSIP_lp, 1, DDSIP_bb->total_int, DDSIP_bb->beg, DDSIP_bb->intind, DDSIP_bb->boundIncrease_intsolvals[scen], DDSIP_bb->effort, DDSIP_bb->Names);
//            if (status)
//            {
//                fprintf (stderr, "ERROR: Failed to copy mip start infos (Warm)\n");
//                fprintf (DDSIP_outfile, "ERROR: Failed to copy mip start infos (Warm)\n");
//                return status;
//            }
//#ifdef DEBUG
//            else if (DDSIP_param->outlev > 21)
//                fprintf (DDSIP_bb->moreoutfile, "##  Copied mip start info %s of previous bound increase (Warm)\n", DDSIP_bb->Names[0]);
//#endif
//            added++;
//        }

        // ADVIND must be 2 when using start values
        status = CPXsetintparam (DDSIP_env, CPX_PARAM_ADVIND, 2);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to set cplex parameter CPX_PARAM_ADVIND.\n");
            fprintf (DDSIP_outfile, "ERROR: Failed to set cplex parameter CPX_PARAM_ADVIND.\n");
            return status;
        }
    }
    // Hotstart = 0 --> ADVIND  0 - start values should not be used
    if (DDSIP_param->hot == 0 && DDSIP_bb->DDSIP_step != dual)
    {
        status = CPXsetintparam (DDSIP_env, CPX_PARAM_ADVIND, 0);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to set cplex parameter CPX_PARAM_ADVIND.\n");
            fprintf (DDSIP_outfile, "ERROR: Failed to set cplex parameter CPX_PARAM_ADVIND.\n");
            return status;
        }
        return 0;
    }
    // Start Value
    // Hotstart = 4 or 6 --> use solutions from father node (from node 1 on), and from previous scenarios
    else if ((DDSIP_param->hot == 4 || DDSIP_param->hot == 6)
#ifdef CBHOTSTART
             || (DDSIP_node[DDSIP_bb->curnode]->step == dual &&
                 (DDSIP_param->cbhot > 2  ||
                  (DDSIP_param->cbhot == 2 && !DDSIP_bb->curnode) ||
                  (DDSIP_param->cbhot == 1 && !DDSIP_bb->curnode && DDSIP_bb->dualdescitcnt < 3))
                )
#endif
            )
    {
        int addMax = 14;
        if (DDSIP_param->cbhot > 2)
            addMax = 24;
        else if (DDSIP_param->cbhot == 2)
            addMax = 19;
#ifdef DEBUG
        // print debugging info
        if (DDSIP_param->outlev > 21)
        {
            fprintf (DDSIP_bb->moreoutfile,"      have %3d MIP starts (from previous iters)\n",added);
        }
#endif
        DDSIP_bb->effort[0]=2;
        added = 0;
        inserted = 0;
        // If not in the root node
        if (DDSIP_bb->curnode && DDSIP_node[DDSIP_bb->curnode]->step != dual)
        {
            // Copy solutions from father node to the problem (see initialization of subbound and solut)
            for(j=0; j<DDSIP_param->scenarios; j++)
            {
                sprintf (DDSIP_bb->Names[added],"Father_s_%d",j+1);
                for(k=0; k<DDSIP_bb->total_int; k++)
                {
                    DDSIP_bb->values[added*DDSIP_bb->total_int+k] = DDSIP_node[DDSIP_bb->curnode]->solut[j*DDSIP_bb->total_int + k];
                }
                already_there = 0;
                for(jt=0; jt<added; jt++)
                {
                    for(k=0; k<DDSIP_bb->total_int; k++)
                        if (!DDSIP_Equal (DDSIP_bb->values[added*DDSIP_bb->total_int+k], DDSIP_bb->values[jt*DDSIP_bb->total_int + k]))
                            break;
                    if (k >= DDSIP_bb->total_int)
                    {
                        already_there = 1;
                        break;
                    }
                }
                if (!already_there)
                    added++;
            }
            // A priori feasibility check (branching may lead to infeasibility)
            kf= DDSIP_node[DDSIP_bb->curnode]->neoind;
            k = DDSIP_bb->firstindex[kf];
            if (DDSIP_bb->firsttype[kf] == 'B' || DDSIP_bb->firsttype[kf] == 'I' || DDSIP_bb->firsttype[kf] == 'N')
            {
                for(i=0; i<DDSIP_bb->total_int; i++)
                    if (k == DDSIP_bb->intind[i])
                        break;
                if (i>DDSIP_bb->total_int-1)
                {
                    /* if this happens it is an error!!! */
                    fprintf (stderr,"XXX ERROR: variable %d branched on is of type %c, but not in DDSIP_bb->intind! DDSIP_bb->intind is:\n", k,DDSIP_bb->firsttype[kf]);
                    fprintf (DDSIP_outfile,"XXX ERROR: variable %d branched on is of type %c, but not in DDSIP_bb->intind! DDSIP_bb->intind is:\n", k,DDSIP_bb->firsttype[kf]);
                    exit (99);
                }
                for(j=0; j<added; j++)
                {
                    if (DDSIP_bb->values[j*DDSIP_bb->total_int + i] < DDSIP_node[DDSIP_bb->curnode]->neolb)
                    {
                        DDSIP_bb->values[j*DDSIP_bb->total_int + i] = DDSIP_node[DDSIP_bb->curnode]->neolb;
                    }
                    else if (DDSIP_bb->values[j*DDSIP_bb->total_int + i] > DDSIP_node[DDSIP_bb->curnode]->neoub)
                    {
                        DDSIP_bb->values[j*DDSIP_bb->total_int + i] = DDSIP_node[DDSIP_bb->curnode]->neolb;
                    }
                }
            }
        }
        // Add solutions of previous scenarios (the directly preceding scenario's are in m1 etc.)
        if (iscen > 1)
        {
//            for (kf=0; kf<iscen - 2; kf++)
            for (kf=iscen - 2; kf>-1; kf--)
            {
                if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[DDSIP_bb->lb_scen_order[kf]])
                {
                    if (!((DDSIP_node[DDSIP_bb->curnode]->first_sol)[DDSIP_bb->lb_scen_order[kf]][DDSIP_bb->firstvar + 1]) || DDSIP_bb->dualitcnt)
                    {
                        // if inherited, the solution is already there from the father
                        already_there = 0;
                        for(jt=0; jt<added; jt++)
                        {
                            for(k=0; k<DDSIP_bb->total_int; k++)
                                if (!DDSIP_Equal (DDSIP_node[DDSIP_bb->curnode]->solut[(DDSIP_bb->lb_scen_order[kf])*DDSIP_bb->total_int+k], DDSIP_bb->values[jt*DDSIP_bb->total_int + k]))
                                    break;
                            if (k >= DDSIP_bb->total_int)
                            {
                                already_there = 1;
                                break;
                            }
                        }
                        if (!already_there)
                        {
                            if (added < 2*DDSIP_param->scenarios)
                            {
                                sprintf (DDSIP_bb->Names[added],"Scen_%d",DDSIP_bb->lb_scen_order[kf]+1);
                                for(k=0; k<DDSIP_bb->total_int; k++)
                                    DDSIP_bb->values [added*DDSIP_bb->total_int+k] = DDSIP_node[DDSIP_bb->curnode]->solut[ (DDSIP_bb->lb_scen_order[kf])*DDSIP_bb->total_int + k];
                                added++;
                                if (DDSIP_node[DDSIP_bb->curnode]->step == dual && added > addMax)
                                    break;
                            }
                            else if (inserted < DDSIP_param->scenarios)
                            {
                                sprintf (DDSIP_bb->Names[inserted],"Scen_%d",DDSIP_bb->lb_scen_order[kf]+1);
                                for(k=0; k<DDSIP_bb->total_int; k++)
                                    DDSIP_bb->values [inserted*DDSIP_bb->total_int+k] = DDSIP_node[DDSIP_bb->curnode]->solut[ (DDSIP_bb->lb_scen_order[kf])*DDSIP_bb->total_int + k];
                                inserted++;
                            }
                            else
                            {
                                fprintf (stderr," XXX WARNING: Scenario %d - No space for MIP start info, added = %d, inserted = %d.\n",DDSIP_bb->lb_scen_order[kf], added, inserted);
                                fprintf (DDSIP_bb->moreoutfile," XXX WARNING: Scenario %d - No space for MIP start info, added = %d, inserted = %d.\n",DDSIP_bb->lb_scen_order[kf], added, inserted);
                            }
                        }
                    }
                }
            }
        }
        // Add start informations to problem
        if (added)
        {
            for (j = 0; j < added; j++)
            {
#ifdef DEBUG
                if (DDSIP_param->outlev > 99)
                {
                    fprintf (DDSIP_bb->moreoutfile,"Scenario %d  - %s:\n",j+1, DDSIP_bb->Names[j]);
                    for (k = 0; k < DDSIP_bb->total_int; k++)
                        fprintf (DDSIP_bb->moreoutfile,"%d:%g, ",DDSIP_bb->intind[k],DDSIP_bb->values[j*DDSIP_bb->total_int + k]);
                    fprintf (DDSIP_bb->moreoutfile,"\n");
                }
#endif
                status = CPXaddmipstarts (DDSIP_env, DDSIP_lp, 1, DDSIP_bb->total_int, DDSIP_bb->beg, DDSIP_bb->intind, DDSIP_bb->values + j*DDSIP_bb->total_int, DDSIP_bb->effort, DDSIP_bb->Names + j);
                if (status)
                {
                    fprintf (stderr, "ERROR: Failed to add mip start infos (Warm)\n");
                    fprintf (DDSIP_outfile, "ERROR: Failed to add mip start infos (Warm)\n");
                    return status;
                }
#ifdef DEBUG
                else if (DDSIP_param->outlev > 80)
                    fprintf (DDSIP_bb->moreoutfile, "  Copied mip start info %s (Warm)\n", DDSIP_bb->Names[j]);
#endif
            }
        }
#ifdef DEBUG
        // print debugging info
        if (DDSIP_param->outlev > 21)
        {
            fprintf (DDSIP_bb->moreoutfile,"      added %2d MIP starts (from other scens)\n",added);
        }
#endif
        // ADVIND must be 2 when using supplied start values
        status = CPXsetintparam (DDSIP_env, CPX_PARAM_ADVIND, 2);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to set cplex parameter CPX_PARAM_ADVIND.\n");
            fprintf (DDSIP_outfile, "ERROR: Failed to set cplex parameter CPX_PARAM_ADVIND.\n");
            return status;
        }
    }
    // Hotstart = 1 --> ADVIND  1 - start values of previous optimization should be used
    else if (DDSIP_param->hot == 1)
    {
        DDSIP_bb->effort[0]=2;
        // If not in the root node add solution of father
        if (DDSIP_bb->curnode && (DDSIP_node[DDSIP_bb->curnode]->step != dual || !DDSIP_bb->dualitcnt))
        {
            // Copy solution from father node to the problem (see initialization of subbound and solut)
            sprintf (DDSIP_bb->Names[0],"Father_%d",scen+1);
            for(k=0; k<DDSIP_bb->total_int; k++)
                DDSIP_bb->values[k] = DDSIP_node[DDSIP_bb->curnode]->solut[scen*DDSIP_bb->total_int + k];
            // A priori feasibility check (branching may lead to infeasibility)
            kf= DDSIP_node[DDSIP_bb->curnode]->neoind;
            k = DDSIP_bb->firstindex[kf];
            if (DDSIP_bb->firsttype[kf] == 'B' || DDSIP_bb->firsttype[kf] == 'I' || DDSIP_bb->firsttype[kf] == 'N')
            {
                for(i=0; i<DDSIP_bb->total_int; i++)
                    if (k == DDSIP_bb->intind[i])
                        break;
                //this should not happen, but for safety ...
                if (i>DDSIP_bb->total_int-1)
                {
                    fprintf (stderr,"XXX ERROR: variable %d branched on is of type %c, but not in DDSIP_bb->intind! DDSIP_bb->intind is:\n", k,DDSIP_bb->firsttype[k]);
                    fprintf (DDSIP_outfile,"XXX ERROR: variable %d branched on is of type %c, but not in DDSIP_bb->intind! DDSIP_bb->intind is:\n", k,DDSIP_bb->firsttype[k]);
                    exit (99);
                }
                if (DDSIP_bb->values[i] < DDSIP_node[DDSIP_bb->curnode]->neolb)
                {
                    DDSIP_bb->values[i] = DDSIP_node[DDSIP_bb->curnode]->neolb;
                }
                else if (DDSIP_bb->values[i] > DDSIP_node[DDSIP_bb->curnode]->neoub)
                {
                    DDSIP_bb->values[i] = DDSIP_node[DDSIP_bb->curnode]->neolb;
                }
            }
            // Add start informations to problem
            status = CPXaddmipstarts (DDSIP_env, DDSIP_lp, 1, DDSIP_bb->total_int, DDSIP_bb->beg, DDSIP_bb->intind, DDSIP_bb->values, DDSIP_bb->effort, DDSIP_bb->Names);
            status = CPXsetintparam (DDSIP_env, CPX_PARAM_ADVIND, 2);
            if (status)
            {
                fprintf (stderr, "ERROR: Failed to set cplex parameter CPX_PARAM_ADVIND.\n");
                fprintf (DDSIP_outfile, "ERROR: Failed to set cplex parameter CPX_PARAM_ADVIND.\n");
                return status;
            }
#ifdef DEBUG
            else if (DDSIP_param->outlev > 90)
                fprintf (DDSIP_bb->moreoutfile, "  Copied mip start info %s of father (Warm)\n", DDSIP_bb->Names[0]);
#endif
        }
        return 0;
    }
    // Hotstart = 3 or 5 --> solutions from previous scenarios in this node
    else if (DDSIP_param->hot == 3 || DDSIP_param->hot == 5)
    {
        status = 0;
        DDSIP_bb->effort[0]=2;
        // Add solution of previous scenarios
        added = 0;
        if (iscen>1)
        {
#ifdef DEBUG
            // print debugging info
            if (DDSIP_param->outlev > 99)
            {
                fprintf (DDSIP_bb->moreoutfile,"     ### MIP starts (from previous scenarios):\n");
            }
#endif
            j = DDSIP_bb->lb_scen_order[scen-2];
            sprintf (DDSIP_bb->Names[added],"Scen_%d",j+1);
            for(k=0; k<DDSIP_bb->total_int; k++)
                DDSIP_bb->values [added*DDSIP_bb->total_int+k] = DDSIP_node[DDSIP_bb->curnode]->solut[j*DDSIP_bb->total_int + k];
            already_there = 0;
            for(jt=0; jt<added; jt++)
            {
                for(k=0; k<DDSIP_bb->total_int; k++)
                    if (!DDSIP_Equal (DDSIP_bb->values [added*DDSIP_bb->total_int+k], DDSIP_bb->values[jt*DDSIP_bb->total_int + k]))
                        break;
                if (k >= DDSIP_bb->total_int)
                {
                    already_there = 1;
                    break;
                }
            }
            if (!already_there)
                added++;
            // Copy starting informations to problem
            for (j = 0; j < added; j++)
            {
#ifdef DEBUG
                // print debugging info
                if (DDSIP_param->outlev > 99)
                {
                    fprintf (DDSIP_bb->moreoutfile,"Scenario %d  - %s:\n",j+1, DDSIP_bb->Names[j]);
                    for (k = 0; k < DDSIP_bb->total_int; k++)
                        fprintf (DDSIP_bb->moreoutfile,"%d:%g, ",DDSIP_bb->intind[k],DDSIP_bb->values[j*DDSIP_bb->total_int + k]);
                    fprintf (DDSIP_bb->moreoutfile,"\n");
                }
#endif
                status = CPXaddmipstarts (DDSIP_env, DDSIP_lp, 1, DDSIP_bb->total_int, DDSIP_bb->beg, DDSIP_bb->intind, DDSIP_bb->values + j*DDSIP_bb->total_int, DDSIP_bb->effort, DDSIP_bb->Names + j);
                if (status)
                {
                    fprintf (stderr, "ERROR: Failed to add mip start infos (Warm)\n");
                    fprintf (DDSIP_outfile, "ERROR: Failed to add mip start infos (Warm)\n");
                    return status;
                }
#ifdef DEBUG
                else if (DDSIP_param->outlev > 90)
                    fprintf (DDSIP_bb->moreoutfile, "  Copied mip start info %s (Warm)\n", DDSIP_bb->Names[j]);
#endif
            }
            if (status)
            {
                fprintf (stderr, "ERROR: Failed to add mip start infos (Warm)\n");
                fprintf (DDSIP_outfile, "ERROR: Failed to add mip start infos (Warm)\n");
                return status;
            }
            // ADVIND must be 2 when using supplied start values
            status = CPXsetintparam (DDSIP_env, CPX_PARAM_ADVIND, 2);
            if (status)
            {
                fprintf (stderr, "ERROR: Failed to set cplex parameter CPX_PARAM_ADVIND.\n");
                fprintf (DDSIP_outfile, "ERROR: Failed to set cplex parameter CPX_PARAM_ADVIND.\n");
                return status;
            }
        }
        else
            added = 0;
    }

    // Hotstart = 2 or 5 or 6
    // Lower bounds are passed to sons and daughters
    // Doesn't work with risk models (needs additional implementations)
    // Doesn't work with conic bundle (when Lagrange multipliers change)
    if (DDSIP_bb->curnode && (DDSIP_param->hot == 2 || DDSIP_param->hot > 4) && !DDSIP_param->riskmod && !DDSIP_param->cb)
    {
        double *values = (double *) DDSIP_Alloc (sizeof (double), 1, "values(Warm)");
        int *rhsind = (int *) DDSIP_Alloc (sizeof (int), 1, "rhsind(Warm)");

        // If solution was feasible, modify lower-bound constraint
        values[0] = DDSIP_node[DDSIP_bb->curnode]->subbound[scen];

        // Debugging information
        if (DDSIP_param->outlev > 5)
            fprintf (DDSIP_bb->moreoutfile, "\tChanging lower bound on objective to %16.14g...\n", values[0]);

        rhsind[0] = DDSIP_bb->objbndind;

        status = CPXchgrhs (DDSIP_env, DDSIP_lp, 1, rhsind, values);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change rhs.(Warm)\n");
            fprintf (DDSIP_outfile, "ERROR: Failed to change rhs.(Warm)\n");
            return status;
        }
        DDSIP_Free ((void **) &(values));
        DDSIP_Free ((void **) &(rhsind));
    }

    return 0;
} // DDSIP_Warm


//==========================================================================
//
int
DDSIP_LowerBoundWriteLp (int scen)
{
    int status = 0;
    // LP-files for debugging are written under sophisticated conditions ;-)
    if ((DDSIP_param->files > 2 && DDSIP_param->files < 5 && ((DDSIP_param->cb && (DDSIP_bb->dualdescitcnt <= 1) && DDSIP_bb->nonode < 2)
            || (DDSIP_param->cb == 0 && DDSIP_bb->nonode == 1)))
            || (DDSIP_bb->DDSIP_step == solve && DDSIP_param->files == 4))
    {
        char fname[DDSIP_ln_fname];

        sprintf (fname, "%s/lb_sc%d_n%d_gc%d%s", DDSIP_outdir, scen + 1, DDSIP_bb->curnode, DDSIP_bb->lboutcnt++, DDSIP_param->coretype);

        status = CPXwriteprob (DDSIP_env, DDSIP_lp, fname, NULL);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to write problem\n");
            return status;
        }
        else if (DDSIP_param->outlev > 20)
        {
#ifdef DEBUG
            printf ("  LP file %s written.\n",fname);
#endif
            fprintf (DDSIP_bb->moreoutfile, "        LP file %s written.\n",fname);
        }
    }
    //
    else if(DDSIP_param->files == 5)
    {
        // write sav, mst and prm files - overwriting each other in cb - only the last ones for each scenario and node remain
        char fname[DDSIP_ln_fname];
        sprintf (fname, "%s/lb_sc%d_n%d%s", DDSIP_outdir, scen + 1, DDSIP_bb->curnode, ".sav.gz");
        status = CPXwriteprob (DDSIP_env, DDSIP_lp, fname, NULL);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to write problem\n");
            return status;
        }
        else if (DDSIP_param->outlev)
        {
            fprintf (DDSIP_bb->moreoutfile, "      problem file   %s written.\n",fname);
        }
        if (CPXgetnummipstarts(DDSIP_env,DDSIP_lp))
        {
            sprintf (fname, "%s/lb_sc%d_n%d%s", DDSIP_outdir, scen + 1, DDSIP_bb->curnode, ".mst.gz");
            status = CPXwritemipstarts(DDSIP_env,DDSIP_lp,fname,0, CPXgetnummipstarts(DDSIP_env,DDSIP_lp)-1);
            if (status)
            {
                fprintf (stderr, "ERROR: Failed to write startinfo\n");
                return status;
            }
            else if (DDSIP_param->outlev)
            {
                fprintf (DDSIP_bb->moreoutfile, "      mip start file %s written.\n",fname);
            }
        }
        sprintf (fname, "%s/lb_sc%d_n%d%s", DDSIP_outdir, scen + 1, DDSIP_bb->curnode, "_1.prm");
        status = CPXwriteparam(DDSIP_env,fname);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to write parameters\n");
            return status;
        }
        else if (DDSIP_param->outlev)
        {
            fprintf (DDSIP_bb->moreoutfile, "      parameter file %s written.\n",fname);
        }
    }
    //
    else if(DDSIP_param->files == 6)
    {
        // write sav, mst and prm files for each scenario, node and cb iteration
        char fname[DDSIP_ln_fname];
        sprintf (fname, "%s/lb_sc%d_n%d_gc%d%s", DDSIP_outdir, scen + 1, DDSIP_bb->curnode, DDSIP_bb->lboutcnt++, ".sav.gz");
        status = CPXwriteprob (DDSIP_env, DDSIP_lp, fname, NULL);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to write problem\n");
            return status;
        }
        else if (DDSIP_param->outlev)
        {
#ifdef DEBUG
            printf ("      problem file   %s written.\n",fname);
#endif
            fprintf (DDSIP_bb->moreoutfile, "      problem file   %s written.\n",fname);
        }
        if (CPXgetnummipstarts(DDSIP_env,DDSIP_lp))
        {
            sprintf (fname, "%s/lb_sc%d_n%d_gc%d%s", DDSIP_outdir, scen + 1, DDSIP_bb->curnode, DDSIP_bb->lboutcnt-1, ".mst.gz");
            status = CPXwritemipstarts(DDSIP_env,DDSIP_lp,fname,0, CPXgetnummipstarts(DDSIP_env,DDSIP_lp)-1);
            if (status)
            {
                fprintf (stderr, "ERROR: Failed to write startinfo\n");
                return status;
            }
            else if (DDSIP_param->outlev)
            {
#ifdef DEBUG
                printf ("      mip start file %s written.\n",fname);
#endif
                fprintf (DDSIP_bb->moreoutfile, "      mip start file %s written.\n",fname);
            }
        }
        sprintf (fname, "%s/lb_sc%d_n%d%s", DDSIP_outdir, scen + 1, DDSIP_bb->curnode, "_1.prm");
        status = CPXwriteparam(DDSIP_env,fname);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to write parameters\n");
            return status;
        }
        else if (DDSIP_param->outlev)
        {
#ifdef DEBUG
            printf ("      parameter file %s written.\n",fname);
#endif
            fprintf (DDSIP_bb->moreoutfile, "      parameter file %s written.\n",fname);
        }
    }

    return 0;
} // DDSIP_LowerBoundWriteLp

//==========================================================================
// Function solves scenario problems in b&b tree and calculates new bounds
int
DDSIP_LowerBound (void)
{
    double objval, bobjval, tmpbestbound = 0.0, tmpbestvalue = 0.0, tmpfeasbound = 0.0, worst_case_lb = -1.e+10;
    double we, wr, mipgap, d, time_start, time_end, time_help, wall_secs, cpu_secs, gap, meanGap, maxGap;
#ifdef DEBUG
    double time_lap;
#endif

    int iscen, i_scen, probs_solved=0, i, j, k, status = 0, optstatus, mipstatus, scen, allopt = 1, relax = 0;
    int wall_hrs, wall_mins,cpu_hrs, cpu_mins;
    int nodes_1st = -1, nodes_2nd = -1;
    int nrCuts;

    int *indices = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_bb->firstvar + DDSIP_bb->secvar),
                                        "indices (LowerBound)");
    double *mipx = (double *) DDSIP_Alloc (sizeof (double), (DDSIP_bb->firstvar + DDSIP_bb->secvar),
                                           "mipx (LowerBound)");
    double *minfirst = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar,
                       "minfirst (LowerBound)");
    double *maxfirst = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar,
                       "maxfirst (LowerBound)");
    char *type = (char *) DDSIP_Alloc (sizeof (char), DDSIP_bb->firstvar, "type(LowerBound)");
    int    *ordind  = NULL;
    double *scensol = NULL;
    double sumprob, maxdispersion, rest_bound, factor;
    char **colname;
    char *colstore;
    sug_t *tmp;

    DDSIP_bb->skip = 0;
    DDSIP_bb->LBIters++;
    // Output
    if (DDSIP_param->outlev)
    {
        fprintf (DDSIP_bb->moreoutfile, "\n----------------------\n");
        fprintf (DDSIP_bb->moreoutfile, "Solving scenario problems (Lower bounds) in node %d:\n", DDSIP_bb->curnode);
        printf ("Solving scenario problems (Lower bounds) in node %d:\n", DDSIP_bb->curnode);
    }
    // Initialization of indices, minfirst, and maxfirst
    for (j = 0; j < DDSIP_bb->firstvar + DDSIP_bb->secvar; j++)
        indices[j] = j;
    for (j = 0; j < DDSIP_bb->firstvar; j++)
    {
        minfirst[j] = DDSIP_infty;
        maxfirst[j] = -DDSIP_infty;
    }

    // CPLEX parameters
    if (DDSIP_bb->DDSIP_step == solve || DDSIP_bb->DDSIP_step == dual)
    {
        status = DDSIP_SetCpxPara (DDSIP_param->cpxnolb, DDSIP_param->cpxlbisdbl, DDSIP_param->cpxlbwhich, DDSIP_param->cpxlbwhat);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
            fprintf (DDSIP_outfile, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
            if (DDSIP_param->outlev)
                fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
            goto TERMINATE;
        }
    }
    // Absolute semideviation: change the lower bound on the risk variable
    // Details: See description of mean-asd-algorithm
    if (DDSIP_param->riskmod == 3 && !DDSIP_param->riskalg && DDSIP_bb->DDSIP_step != dual)
    {
        DDSIP_SemDevChgBd ();

        // Node can be pruned ?
        if (DDSIP_node[DDSIP_bb->curnode]->target > DDSIP_bb->bestvalue)
        {
            DDSIP_node[DDSIP_bb->curnode]->bound = DDSIP_node[DDSIP_bb->curnode]->target;
            if (DDSIP_param->outlev)
                fprintf (DDSIP_bb->moreoutfile, "Pruning node. \n");
            DDSIP_bb->skip=2;
            goto TERMINATE;
        }
    }
    // Relax
    relax = 0;
    if (DDSIP_param->relax > 2 && DDSIP_bb->curnode % DDSIP_param->relax)
        relax = 2;
    else if (DDSIP_param->relax == 2)
        relax = 2;
    else if (DDSIP_param->relax == 1)
        relax = 1;

    // Relax first stage
    if (relax == 1)
    {
        for (j = 0; j < DDSIP_bb->firstvar; j++)
        {
            if (DDSIP_bb->firsttype[j] == 'B' || DDSIP_bb->firsttype[j] == 'I')
                type[j] = 'C';
        }
        status = CPXchgctype (DDSIP_env, DDSIP_lp, DDSIP_bb->firstvar, DDSIP_bb->firstindex, type);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change types of first stage variables (LowerBound) \n");
            fprintf (DDSIP_outfile, "ERROR: Failed to change types of first stage variables (LowerBound) \n");
            if (DDSIP_param->outlev)
                fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to change types of first stage variables (LowerBound) \n");
            goto TERMINATE;
        }
    }
    // if specified, Relax integrality of all variables
    if (relax == 2)
        status = CPXchgprobtype (DDSIP_env, DDSIP_lp, CPXPROB_LP);

    // Change bounds according to node in bb tree
    if (DDSIP_bb->curbdcnt)
    {
        status = DDSIP_ChgBounds (1);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change bounds \n");
            fprintf (DDSIP_outfile, "ERROR: Failed to change bounds \n");
            return status;
        }
    }
    // This node has been solved (as soon as we enter the loop below)
    DDSIP_node[DDSIP_bb->curnode]->solved = 1;
    DDSIP_bb->violations = DDSIP_bb->firstvar;
    tmpbestbound = 0.0;
    tmpbestvalue = 0.0;
    // prepare for decision about stopping: rest_bound is the expectation of all the lower bounds in the father node
    if (DDSIP_bb->curnode /* && DDSIP_node[DDSIP_node[DDSIP_bb->curnode]->father]->step != dual */)
    {
        if(!DDSIP_param->cb || DDSIP_node[DDSIP_node[DDSIP_bb->curnode]->father]->step != dual)
        {
            rest_bound = DDSIP_node[DDSIP_node[DDSIP_bb->curnode]->father]->bound;
        }
        else
        {
            rest_bound = (DDSIP_node[DDSIP_bb->curnode]->subbound)[0] * DDSIP_data->prob[0] - DDSIP_param->accuracy;
            for (scen = 1; scen < DDSIP_param->scenarios; scen++)
            {
                rest_bound += (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen] * DDSIP_data->prob[scen];
            }
            rest_bound = DDSIP_Dmin(rest_bound, DDSIP_node[DDSIP_node[DDSIP_bb->curnode]->father]->bound);
        }
    }
    else
        rest_bound = -DDSIP_param->scenarios*DDSIP_infty;
    // LowerBound problem for each scenario
    //****************************************************************************
    maxGap = meanGap = 0.;
    nrCuts = DDSIP_bb->cutCntr;
    for (iscen = 0; iscen < DDSIP_param->scenarios; iscen++)
    {
        int nrFeasCheck = 0;
        scen=DDSIP_bb->lb_scen_order[iscen];
        // If even with the lower bounds of the father we would reach a greater value we may stop and save some time

        if (rest_bound > -DDSIP_infty)
        {
            rest_bound -= (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen] * DDSIP_data->prob[scen];
        }

        // User termination
        if (DDSIP_killsignal)
        {
            DDSIP_bb->skip = 1;
            goto TERMINATE;
        }
        //
        // Solve problem only if the solution is not inherited from father node
        //
        if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen])
        {
            if (DDSIP_param->cpxscr)
            {
                printf ("*****************************************************************************\n");
                printf ("Inherited solution for scenario problem %d  ....\n", scen + 1);
            }
            if (!scen)
                // Initialize Warm starts in case the solution of scen 0 is inherited
                status = DDSIP_Warm (iscen);
            gap = 100.0*((DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen]-(DDSIP_node[DDSIP_bb->curnode]->subbound)[scen])/
                  (fabs((DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen])+1e-4);
            if (fabs(gap) < 2.e-13)
                gap = 0.;
            meanGap += DDSIP_data->prob[scen] * gap;
            time (&DDSIP_bb->cur_time);
            if (DDSIP_param->outlev)
            {
                DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                if ((DDSIP_node[DDSIP_bb->curnode]->mipstatus)[scen] == CPXMIP_OPTIMAL)
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%)     ",
                             iscen + 1,scen + 1, (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen], (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen], gap);
                else if ((DDSIP_node[DDSIP_bb->curnode]->mipstatus)[scen] == CPXMIP_OPTIMAL_TOL)
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) tol.",
                             iscen + 1,scen + 1, (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen], (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen], gap);
                else if ((DDSIP_node[DDSIP_bb->curnode]->mipstatus)[scen] == CPXMIP_NODE_LIM_FEAS)
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) NODE",
                             iscen + 1,scen + 1, (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen], (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen], gap);
                else if ((DDSIP_node[DDSIP_bb->curnode]->mipstatus)[scen] == CPXMIP_TIME_LIM_FEAS)
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) TIME",
                             iscen + 1,scen + 1, (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen], (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen], gap);
                else
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) %-4.0d ",
                             iscen + 1,scen + 1, (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen], (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen], gap,
                             (DDSIP_node[DDSIP_bb->curnode]->mipstatus)[scen]);
                fprintf (DDSIP_bb->moreoutfile,
                         "\t %3dh %02d:%02.0f / %3dh %02d:%05.2f  > node %3g",
                         wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs, (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][DDSIP_bb->firstvar + 2]);
                /*
                                if (DDSIP_param->outlev > 8)
                                {
                                    printf ("%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%)\tStat.%3.0d\t%3dh %02d:%02.0f / %3dh %02d:%05.2f  > node %3g",
                                            iscen + 1 ,scen + 1, (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen], (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen],
                                            gap, (DDSIP_node[DDSIP_bb->curnode]->mipstatus)[scen], wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs,
                                            (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][DDSIP_bb->firstvar + 2]);
                                    if (DDSIP_param->cpxscr)
                                    {
                                        printf ("\n");
                                    }
                                }
                */
            }
            //**************************************************************************** output of first stage
            if (DDSIP_param->outlev)
            {
                for (i_scen = 0; i_scen < scen; i_scen++)
                {
                    // if the solutions have the same address, they were equal
                    if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[i_scen] == (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen])
                        break;
                }
                if (i_scen < scen)
                {
                    //first_sol[DDSIP_bb->firstvar]   equals the number of identical first stage scenario solutions
                    fprintf (DDSIP_bb->moreoutfile,"    \t\t->scen %4d (%3g ident.)\n", i_scen + 1, DDSIP_node[DDSIP_bb->curnode]->first_sol[scen][DDSIP_bb->firstvar]);
                }
                else if (DDSIP_param->outlev >= DDSIP_first_stage_outlev)
                {
                    fprintf (DDSIP_bb->moreoutfile, "\n    First-stage solution:\n");
                    for (j = 0; j < DDSIP_bb->firstvar; j++)
                    {
                        fprintf (DDSIP_bb->moreoutfile," %21.14g", (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j]);
                        if (!((j + 1) % 5))
                            fprintf (DDSIP_bb->moreoutfile, "\n");
                    }
                    fprintf (DDSIP_bb->moreoutfile, "\n");
                }
                else
                {
                    fprintf (DDSIP_bb->moreoutfile, "\n");
                }
            }
            //****************************************************************************
            DDSIP_bb->solstat[scen] = 1;
            if (!((DDSIP_node[DDSIP_bb->curnode]->mipstatus)[scen] == CPXMIP_OPTIMAL))
                allopt = 0;
            // Temporary bound is infinity if a scenario problem is infeasible
            if (DDSIP_Equal ((DDSIP_node[DDSIP_bb->curnode]->subbound)[scen], DDSIP_infty))
                tmpbestbound = DDSIP_infty;
            else
            {
                tmpbestbound += DDSIP_data->prob[scen] * (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen];
            }
            tmpbestvalue += DDSIP_data->prob[scen] * (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen];
            for (j = 0; j < DDSIP_bb->firstvar; j++)
            {
                // Calculate minimum and maximum of each component
                minfirst[j] = DDSIP_Dmin ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j], minfirst[j]);
                maxfirst[j] = DDSIP_Dmax ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j], maxfirst[j]);
            }
            // we don't know the second stage solution after branching
            for (j = 0; j < DDSIP_bb->secvar; j++)
            {
                DDSIP_bb->cur_secstage[scen][j] = DDSIP_infty;
            }
            if ((DDSIP_param->scalarization || DDSIP_param->cb) && (DDSIP_param->riskmod))
            {
                if ((DDSIP_node[DDSIP_bb->curnode]->ref_scenobj)[scen] > DDSIP_param->risktarget)
                    DDSIP_bb->ref_risk[scen] = 1;
                else
                    DDSIP_bb->ref_risk[scen] = 0;
            }
            factor = (DDSIP_bb->bestvalue < 0.)? 1.-2.e-15 :  1.+2.e-15;
            // Node can be fathomed ?
            if (rest_bound > -DDSIP_infty /*&& DDSIP_param->prematureStop*/ && !DDSIP_param->riskalg && iscen < DDSIP_param->scenarios-1 && ((tmpbestbound + rest_bound) > (DDSIP_bb->bestvalue*factor)) && !DDSIP_param->scalarization)
            {
                DDSIP_node[DDSIP_bb->curnode]->bound = tmpbestbound + rest_bound;
                //
#ifdef DEBUG
                if (DDSIP_param->outlev > 21)
                {
                    fprintf(DDSIP_bb->moreoutfile," --premature stop with tmpbestbound + rest_bound: %22.14g + %22.14g = %22.14g\n", tmpbestbound, rest_bound, DDSIP_node[DDSIP_bb->curnode]->bound);
                    fprintf(DDSIP_bb->moreoutfile,"                       tmpbestbound + rest_bound: %22.14g > %22.14g = %22.14g*%g, difference: %g\n", DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->bestvalue*factor, DDSIP_bb->bestvalue, factor, DDSIP_node[DDSIP_bb->curnode]->bound-(DDSIP_bb->bestvalue*factor));
                }
#endif
                //
                DDSIP_bb->skip = 2;
                //DDSIP_bb->cutoff++;
                if (DDSIP_param->outlev)
                {
                    fprintf (DDSIP_bb->moreoutfile, "\tLower bound of node %d >= %18.16g (-bestvalue = %g) cutoff after evaluation of %d scenarios, skip=%d\n", DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bestvalue,iscen+1,DDSIP_bb->skip);
                }
                DDSIP_node[DDSIP_bb->curnode]->leaf = 1;
                goto TERMINATE;
            }
        }
        else
        {
            nodes_1st = nodes_2nd = -1;
            // Output
            if (DDSIP_param->outlev || DDSIP_param->cpxscr)
            {
                if (DDSIP_param->cpxscr)
                    printf ("*****************************************************************************\n");
                printf ("Solving scenario problem %d (for lower bound) in node %d....\n", scen + 1, DDSIP_bb->curnode);
            }
NEXT_TRY:
            probs_solved++;
            // Change problem according to scenario
            status = DDSIP_ChgProb (scen, DDSIP_bb->multipliers);
            if (status)
            {
                fprintf (stderr, "ERROR: Failed to change problem \n");
                fprintf (DDSIP_outfile, "ERROR: Failed to change problem \n");
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to change problem \n");
                goto TERMINATE;
            }

            // Warm starts
            status = DDSIP_Warm (iscen);
            // Error?
            if (status)
                goto TERMINATE;

            status = DDSIP_LowerBoundWriteLp (scen);
            // Error?
            if (status)
                goto TERMINATE;
            time_start = DDSIP_GetCpuTime ();
            // Optimize
            if (relax < 2)
            {
#ifdef DEBUG
                // query time limit amd mip rel. gap parameters
                if (DDSIP_param->cpxscr || DDSIP_param->outlev > 21)
                {
                    status = CPXgetdblparam (DDSIP_env,CPX_PARAM_TILIM,&we);
                    status = CPXgetdblparam (DDSIP_env,CPX_PARAM_EPGAP,&wr);
                    printf ("   -- 1st optimization time limit: %gs, rel. gap: %g%% --\n",we,wr*100.0);
                }
#endif
                // Optimize MIP
                DDSIP_bb->scenLBIters++;
                optstatus = CPXmipopt (DDSIP_env, DDSIP_lp);
                mipstatus = CPXgetstat (DDSIP_env, DDSIP_lp);
                nodes_1st = CPXgetnodecnt (DDSIP_env,DDSIP_lp);
                if (!optstatus && !DDSIP_Error(optstatus) && !DDSIP_Infeasible (mipstatus))
                {
                    if (CPXgetmiprelgap(DDSIP_env, DDSIP_lp, &mipgap))
                    {
                        fprintf (stderr, "ERROR: Query of CPLEX mip gap from 1st optimization failed (LowerBound) \n");
                        fprintf (stderr, "       CPXgetstat returned: %d\n",mipstatus);
                        mipgap = 1.e+30;
                    }
#ifdef DEBUG
                    time_lap = DDSIP_GetCpuTime ();
                    if (DDSIP_param->cpxscr || DDSIP_param->outlev > 11)
                    {
                        printf ("      LB: after 1st optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)\n",mipgap*100.0,nodes_1st,time_lap-time_start);
                        if (DDSIP_param->outlev)
                            fprintf (DDSIP_bb->moreoutfile,"      LB: after 1st optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)\n",mipgap*100.0,nodes_1st,time_lap-time_start);
                    }
#endif
                    if (DDSIP_param->watchkappa)
                    {
                        double maxkappaval, stablekappaval, suspiciouskappaval, unstablekappaval, illposedkappaval;
                        status = (CPXgetdblquality(DDSIP_env,DDSIP_lp,&maxkappaval,CPX_KAPPA_MAX) ||
                                  CPXgetdblquality(DDSIP_env,DDSIP_lp,&illposedkappaval,CPX_KAPPA_ILLPOSED) ||
                                  CPXgetdblquality(DDSIP_env,DDSIP_lp,&stablekappaval,CPX_KAPPA_STABLE) ||
                                  CPXgetdblquality(DDSIP_env,DDSIP_lp,&suspiciouskappaval,CPX_KAPPA_SUSPICIOUS) ||
                                  CPXgetdblquality(DDSIP_env,DDSIP_lp,&unstablekappaval,CPX_KAPPA_UNSTABLE));
                        if (status)
                        {
                            fprintf (stderr, "Failed to obtain Kappa information. \n\n");
                        }
                        else
                        {
                            printf("maximal Kappa value        = %g\n",maxkappaval);
                            if (stablekappaval < 1.)
                            {
                                if (stablekappaval)
                                    printf("numerically stable bases     %7.3f%%\n",stablekappaval*1e2);
                                if (suspiciouskappaval)
                                    printf("numerically suspicious bases %7.3f%%\n",suspiciouskappaval*1e2);
                                if (unstablekappaval)
                                    printf("numerically unstable bases   %7.3f%%\n",unstablekappaval*1e2);
                                if (illposedkappaval)
                                    printf("numerically ill-posed bases  %7.3f%%\n",illposedkappaval*1e2);
                            }
                            if (DDSIP_param->outlev)
                            {
                                fprintf(DDSIP_bb->moreoutfile, "             maximal Kappa value        = %g\n",maxkappaval);
                                if (stablekappaval < 1.)
                                {
                                    if (stablekappaval)
                                        fprintf(DDSIP_bb->moreoutfile, "             numerically stable bases     %7.3f%%\n",stablekappaval*1e2);
                                    if (suspiciouskappaval)
                                        fprintf(DDSIP_bb->moreoutfile, "             numerically suspicious bases %7.3f%%\n",suspiciouskappaval*1e2);
                                    if (unstablekappaval)
                                        fprintf(DDSIP_bb->moreoutfile, "             numerically unstable bases   %7.3f%%\n",unstablekappaval*1e2);
                                    if (illposedkappaval)
                                        fprintf(DDSIP_bb->moreoutfile, "             numerically ill-posed bases  %7.3f%%\n",illposedkappaval*1e2);
                                }
                            }
                        }
                    }

                    if (DDSIP_param->cpxnolb2 && mipstatus != CPXMIP_OPTIMAL)
                    {
                        // more iterations with different settings
                        status = DDSIP_SetCpxPara (DDSIP_param->cpxnolb2, DDSIP_param->cpxlbisdbl2, DDSIP_param->cpxlbwhich2, DDSIP_param->cpxlbwhat2);
                        if (status)
                        {
                            fprintf (stderr, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
                            fprintf (DDSIP_outfile, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
                            if (DDSIP_param->outlev)
                                fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
                            goto TERMINATE;
                        }
                        // query time limit amd mip rel. gap parameters
                        status = CPXgetdblparam (DDSIP_env,CPX_PARAM_EPGAP,&wr);
#ifdef DEBUG
                        if (DDSIP_param->cpxscr || DDSIP_param->outlev > 21)
                        {
                            status = CPXgetdblparam (DDSIP_env,CPX_PARAM_TILIM,&we);
                            printf ("   -- 2nd optimization time limit: %gs, rel. gap: %g%% --\n",we,wr*100.0);
                        }
#endif
                        // continue if desired gap is not reached yet
                        if (mipgap > wr)
                        {
                            optstatus = CPXmipopt (DDSIP_env, DDSIP_lp);
                            mipstatus = CPXgetstat (DDSIP_env, DDSIP_lp);
                            nodes_2nd = CPXgetnodecnt (DDSIP_env,DDSIP_lp);
#ifdef DEBUG
                            if (DDSIP_param->cpxscr || DDSIP_param->outlev > 11)
                            {
                                if (CPXgetmiprelgap(DDSIP_env, DDSIP_lp, &mipgap))
                                {
                                    fprintf (stderr, "ERROR: Query of CPLEX mip gap from 2nd optimization failed (LowerBound) \n");
                                    fprintf (stderr, "       CPXgetstat returned: %d\n",mipstatus);
                                    mipgap = 1.e+30;
                                }
                                else
                                {
                                    time_end = DDSIP_GetCpuTime ();
                                    printf ("      LB: after 2nd optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)\n",mipgap*100.0,nodes_2nd,time_end-time_lap);
                                    if (DDSIP_param->outlev)
                                        fprintf (DDSIP_bb->moreoutfile,"      LB: after 2nd optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)\n",mipgap*100.0,nodes_2nd,time_end-time_lap);
                                }
                            }
#endif
                            if (DDSIP_param->watchkappa)
                            {
                                double maxkappaval, stablekappaval, suspiciouskappaval, unstablekappaval, illposedkappaval;
                                status = (CPXgetdblquality(DDSIP_env,DDSIP_lp,&maxkappaval,CPX_KAPPA_MAX) ||
                                          CPXgetdblquality(DDSIP_env,DDSIP_lp,&illposedkappaval,CPX_KAPPA_ILLPOSED) ||
                                          CPXgetdblquality(DDSIP_env,DDSIP_lp,&stablekappaval,CPX_KAPPA_STABLE) ||
                                          CPXgetdblquality(DDSIP_env,DDSIP_lp,&suspiciouskappaval,CPX_KAPPA_SUSPICIOUS) ||
                                          CPXgetdblquality(DDSIP_env,DDSIP_lp,&unstablekappaval,CPX_KAPPA_UNSTABLE));
                                if (status)
                                {
                                    fprintf (stderr, "Failed to obtain Kappa information. \n\n");
                                }
                                else
                                {
                                    printf("maximal Kappa value        = %g\n",maxkappaval);
                                    if (stablekappaval < 1.)
                                    {
                                        if (stablekappaval)
                                            printf("numerically stable bases     %7.3f%%\n",stablekappaval*1e2);
                                        if (suspiciouskappaval)
                                            printf("numerically suspicious bases %7.3f%%\n",suspiciouskappaval*1e2);
                                        if (unstablekappaval)
                                            printf("numerically unstable bases   %7.3f%%\n",unstablekappaval*1e2);
                                        if (illposedkappaval)
                                            printf("numerically ill-posed bases  %7.3f%%\n",illposedkappaval*1e2);
                                    }
                                    if (DDSIP_param->outlev)
                                    {
                                        fprintf(DDSIP_bb->moreoutfile, "             maximal Kappa value        = %g\n",maxkappaval);
                                        if (stablekappaval < 1.)
                                        {
                                            if (stablekappaval)
                                                fprintf(DDSIP_bb->moreoutfile, "             numerically stable bases     %7.3f%%\n",stablekappaval*1e2);
                                            if (suspiciouskappaval)
                                                fprintf(DDSIP_bb->moreoutfile, "             numerically suspicious bases %7.3f%%\n",suspiciouskappaval*1e2);
                                            if (unstablekappaval)
                                                fprintf(DDSIP_bb->moreoutfile, "             numerically unstable bases   %7.3f%%\n",unstablekappaval*1e2);
                                            if (illposedkappaval)
                                                fprintf(DDSIP_bb->moreoutfile, "             numerically ill-posed bases  %7.3f%%\n",illposedkappaval*1e2);
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            if (mipstatus == CPXMIP_TIME_LIM_FEAS || mipstatus == CPXMIP_NODE_LIM_FEAS)
                                mipstatus = CPXMIP_OPTIMAL_TOL;
                        }
                        // reset CPLEX parameters to lb
                        status = DDSIP_SetCpxPara (DDSIP_param->cpxnolb, DDSIP_param->cpxlbisdbl, DDSIP_param->cpxlbwhich, DDSIP_param->cpxlbwhat);
                        if (status)
                        {
                            fprintf (stderr, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
                            fprintf (DDSIP_outfile, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
                            if (DDSIP_param->outlev)
                                fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
                            goto TERMINATE;
                        }
                    }
                }

                if ((k = CPXgetnummipstarts(DDSIP_env, DDSIP_lp)) > 2)
                {
                    status    = CPXdelmipstarts (DDSIP_env, DDSIP_lp, 1, k-1);
                }
            }
            else
            {
                // Optimize relaxed MIP
                optstatus = CPXdualopt (DDSIP_env, DDSIP_lp);
                mipstatus = CPXgetstat (DDSIP_env, DDSIP_lp);
            }
            // Infeasible, unbounded ... ?
            // We handle some errors separately (blatant infeasible, error in scenario problem)
            if (DDSIP_Error (optstatus))
            {
                fprintf (stderr, "ERROR: Failed to optimize problem for scenario %d in node %d.(LowerBound), status=%d\n", scen + 1, DDSIP_bb->curnode, optstatus);
                fprintf (DDSIP_outfile, "ERROR: Failed to optimize problem for scenario %d in node %d.(LowerBound), status=%d\n", scen + 1, DDSIP_bb->curnode, optstatus);
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to optimize problem for scenario %d in node %d.(LowerBound) status=%d\n",
                             scen + 1, optstatus, DDSIP_bb->curnode);
                if (optstatus == CPXERR_SUBPROB_SOLVE)
                {
                    optstatus = CPXgetsubstat(DDSIP_env, DDSIP_lp);
                    fprintf (stderr, "ERROR:                                scenario %d.(LowerBound) subproblem status=%d\n", scen + 1, optstatus);
                    fprintf (DDSIP_outfile, "ERROR:                                scenario %d.(LowerBound) subproblem status=%d\n", scen + 1, optstatus);
                    if (DDSIP_param->outlev)
                        fprintf (DDSIP_bb->moreoutfile, "ERROR:                                scenario %d.(LowerBound) subproblem status=%d\n", scen + 1, optstatus);
                }
                status = optstatus;
                goto TERMINATE;
            }

            time (&DDSIP_bb->cur_time);
            // Infeasible, unbounded ... ?
            if (mipstatus==CPXMIP_TIME_LIM_INFEAS)
            {
                if (DDSIP_param->outlev)
                {
                    DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                    DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:                          \tinfeasible               \t             \tStat.%3.0d\t%3dh %02d:%02.0f / %3dh %02d:%05.2f\n",
                             iscen + 1, scen + 1, mipstatus, wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs);
                }
                DDSIP_node[DDSIP_bb->curnode]->bound = DDSIP_infty;
                CPXgetdblparam (DDSIP_env,CPX_PARAM_TILIM,&cpu_secs);
                printf ("ERROR: Problem infeasible within time limit of %g secs for scenario %d (LowerBound)\n", cpu_secs, scen + 1);
                fprintf (DDSIP_outfile, "ERROR: Problem infeasible within time limit of %g secs for scenario %d (LowerBound)\n", cpu_secs, scen + 1);
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile, "ERROR: Problem infeasible within time limit of %g secs for scenario %d (LowerBound)\n", cpu_secs, scen + 1);
                DDSIP_bb->skip = -2;
                for (k = scen; k < DDSIP_param->scenarios; k++)
                    DDSIP_bb->solstat[k] = 0;
                bobjval = DDSIP_infty;
                status = mipstatus;
                goto TERMINATE;
            }
            else if (DDSIP_Infeasible (mipstatus))
            {
                if (DDSIP_param->outlev)
                {
                    DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                    DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:                          \tinfeasible               \t             \tStat.%3.0d\t%3dh %02d:%02.0f / %3dh %02d:%05.2f\n",
                             iscen + 1, scen + 1, mipstatus, wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs);
                }
                // shift the infeasible scenario to be solved first
                if (iscen)
                {
                    i_scen = DDSIP_bb->lb_scen_order[iscen];
                    for(j = iscen; j>0; j--)
                        DDSIP_bb->lb_scen_order[j] = DDSIP_bb->lb_scen_order[j-1];
                    DDSIP_bb->lb_scen_order[0] = i_scen;
                }
                DDSIP_node[DDSIP_bb->curnode]->bound = DDSIP_infty;
                DDSIP_bb->skip = 1;
                for (k = scen; k < DDSIP_param->scenarios; k++)
                    DDSIP_bb->solstat[k] = 0;
                goto TERMINATE;
            }
            // We did not detect infeasibility
            // If we couldn't find a feasible solution we can at least obtain a lower bound
            if (DDSIP_NoSolution (mipstatus))
            {
                objval = DDSIP_infty;
                DDSIP_bb->solstat[scen] = 0;
                if (relax != 2)
                    status = CPXgetbestobjval (DDSIP_env, DDSIP_lp, &bobjval);
                else
                    bobjval = DDSIP_infty;
                if (status)
                {
                    fprintf (stderr, "ERROR: Failed to get value of best remaining node (LowerBound)\n");
                    fprintf (DDSIP_outfile, "ERROR: Failed to get value of best remaining node (LowerBound)\n");
                    if (DDSIP_param->outlev)
                        fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to get value of best remaining node (LowerBound)\n");
                    goto TERMINATE;
                }
            }
            else              // Solution exists
            {
                DDSIP_bb->solstat[scen] = 1;
                status = CPXgetobjval (DDSIP_env, DDSIP_lp, &objval);
                if (status)
                {
                    fprintf (stderr, "ERROR: Failed to get best objective value \n");
                    fprintf (DDSIP_outfile, "ERROR: Failed to get best objective value \n");
                    if (DDSIP_param->outlev)
                        fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to get best objective value \n");
                    goto TERMINATE;
                }

                status = CPXgetx (DDSIP_env, DDSIP_lp, mipx, 0, DDSIP_bb->firstvar + DDSIP_bb->secvar - 1);
                if (status)
                {
                    fprintf (stderr, "ERROR: Failed to get solution \n");
                    fprintf (DDSIP_outfile, "ERROR: Failed to get solution \n");
                    if (DDSIP_param->outlev)
                        fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to get solution \n");
                    goto TERMINATE;
                }
                // Numerical errors?
                for (j = 0; j < DDSIP_bb->firstvar; j++)
                {
#ifdef CHECK_BOUNDS
                    // outside of bounds?
                    // loose - restrict the first-stage variables inside of original bounds (disregarding branching)
                    if (DDSIP_bb->lborg[j] - mipx[DDSIP_bb->firstindex[j]] > (DDSIP_bb->lborg[j] + 1.) * 5.e-10)
                    {
//#ifdef DEBUG
                        if (DDSIP_param->outlev > 21)
                            fprintf (DDSIP_bb->moreoutfile, "   ## First-stage var. %6d (%c):  lb= %.16g > %.16g =mipx, diff= %g\n", DDSIP_bb->firstindex[j], DDSIP_bb->firsttype[j], DDSIP_bb->lborg[j], mipx[DDSIP_bb->firstindex[j]], DDSIP_bb->lborg[j] - mipx[DDSIP_bb->firstindex[j]]);
//#endif
                        mipx[DDSIP_bb->firstindex[j]] = DDSIP_bb->lborg[j];
                    }
                    else
                    {
                        if (mipx[DDSIP_bb->firstindex[j]] - DDSIP_bb->uborg[j] > (DDSIP_bb->uborg[j] + 1.) * 5.e-10)
                        {
//#ifdef DEBUG
                            if (DDSIP_param->outlev > 21)
                                fprintf (DDSIP_bb->moreoutfile, "   ## First-stage var. %6d (%c):  ub= %.16g < %.16g =mipx, diff= %g\n", DDSIP_bb->firstindex[j], DDSIP_bb->firsttype[j], DDSIP_bb->uborg[j], mipx[DDSIP_bb->firstindex[j]], mipx[DDSIP_bb->firstindex[j]] - DDSIP_bb->uborg[j]);
//#endif
                            mipx[DDSIP_bb->firstindex[j]] = DDSIP_bb->uborg[j];
                        }
                    }
#endif
                    // tiny value?
                    if (fabs (mipx[DDSIP_bb->firstindex[j]]) < 1e-18)
                    {
//#ifdef DEBUG
                        if (DDSIP_param->outlev > 21 && mipx[DDSIP_bb->firstindex[j]])
                            fprintf (DDSIP_bb->moreoutfile, "   ## First-stage var. %6d (%c):  tiny mipx %.6g -> 0.\n", DDSIP_bb->firstindex[j], DDSIP_bb->firsttype[j], mipx[DDSIP_bb->firstindex[j]]);
//#endif
                        mipx[DDSIP_bb->firstindex[j]] = 0.0;
                    }
                }
                if (relax != 2)
                {
                    //if (mipstatus == CPXMIP_OPTIMAL)
                    //{
                    //    bobjval = objval;
                    //}
                    //else
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
                }
                else
                    bobjval = DDSIP_infty;

            }                // end for else

            // If all scenario problems were solved to optimality ....
            if (!(mipstatus == CPXMIP_OPTIMAL))
            {
                allopt = 0;
                // If the tree was exhausted bobjval is huge
                if (fabs (bobjval) > DDSIP_infty)
                    bobjval = objval;
            }
            (DDSIP_node[DDSIP_bb->curnode]->mipstatus)[scen] = mipstatus;
            // Debugging information
            gap = 100.0*(objval-bobjval)/(fabs(objval)+1e-4);
            if (fabs(gap) < 2.e-13)
                gap = 0.;
            meanGap += DDSIP_data->prob[scen] * gap;
            maxGap   = DDSIP_Dmax (maxGap, gap);
            time_end = DDSIP_GetCpuTime ();
            time_help = time_end-time_start;
#ifdef SHIFT
            if (DDSIP_param->hot == 4)
            {
                // prepare for shifting difficult scenarios to the end such that they can make use of more mipstarts
                DDSIP_bb->aggregate_time[iscen] = time_help + 1e3*gap;
            }
#endif
            time (&DDSIP_bb->cur_time);
            if (DDSIP_param->outlev)
            {
                DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                DDSIP_translate_time (time_end,&cpu_hrs,&cpu_mins,&cpu_secs);
                if (mipstatus == CPXMIP_OPTIMAL)
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%)     ",
                             iscen + 1, scen + 1, objval, bobjval, gap);
                else if (mipstatus == CPXMIP_OPTIMAL_TOL)
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) tol.",
                             iscen + 1, scen + 1, objval, bobjval, gap);
                else if (mipstatus == CPXMIP_NODE_LIM_FEAS)
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) NODE",
                             iscen + 1, scen + 1, objval, bobjval, gap);
                else if (mipstatus == CPXMIP_TIME_LIM_FEAS)
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) TIME",
                             iscen + 1, scen + 1, objval, bobjval, gap);
                else
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) %-4.0d ",
                             iscen + 1, scen + 1, objval, bobjval, gap, mipstatus);
                fprintf (DDSIP_bb->moreoutfile,
                         "\t %3dh %02d:%02.0f / %3dh %02d:%05.2f (%7.2fs n:%5d",
                         wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs, time_help, nodes_1st);
                if (nodes_2nd >= 0)
                    fprintf (DDSIP_bb->moreoutfile, " +%4d)", nodes_2nd - nodes_1st);
                else
                    fprintf (DDSIP_bb->moreoutfile, ")");
            }
            // Remember values of integer variables for warm starts
            if (DDSIP_bb->solstat[scen] && DDSIP_param->hot)
            {
                for (j = 0; j < DDSIP_bb->total_int; j++)
                    DDSIP_node[DDSIP_bb->curnode]->solut[scen * DDSIP_bb->total_int + j] = mipx[DDSIP_bb->intind[j]];
            }

            wr = DDSIP_Dmin (bobjval, objval);
            (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen] = wr;
            (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen] = objval;
            // Remember objective function contribution
            DDSIP_Contrib_LB (mipx, scen);
            if (DDSIP_Equal ((DDSIP_node[DDSIP_bb->curnode]->subbound)[scen], DDSIP_infty))
                tmpbestbound = DDSIP_infty;


            // Wait-and-see lower bounds, should be set only once
            // Test if btlb contain still the initial value -DDSIP_infty
            if (!(DDSIP_bb->btlb[scen] > -DDSIP_infty))
                DDSIP_bb->btlb[scen] = wr;
            i_scen = 999999999;
            // If a feasible solution exists
            if (DDSIP_bb->solstat[scen])
            {
                for (i_scen = 0; i_scen < DDSIP_param->scenarios; i_scen++)
                {
                    if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[i_scen])
                    {
                        for (j = 0; j < DDSIP_bb->firstvar; j++)
                        {
                            if (!DDSIP_Equal (mipx[DDSIP_bb->firstindex[j]], (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i_scen][j]))
                                break;
                        }
                        if (j >= DDSIP_bb->firstvar)
                        {
                            break;
                        }
                    }
                }
                if (i_scen < DDSIP_param->scenarios)
                {
                    (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen] = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i_scen];
                    //first_sol[DDSIP_bb->firstvar]   equals the number of identical first stage scenario solutions
                    (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][DDSIP_bb->firstvar] += 1.0;
                    if (DDSIP_param->outlev)
                    {
                        if (nodes_2nd < 0)
                            fprintf (DDSIP_bb->moreoutfile," \t\t->scen %4d (%3g ident.)\n", i_scen + 1, DDSIP_node[DDSIP_bb->curnode]->first_sol[scen][DDSIP_bb->firstvar]);
                        else
                            fprintf (DDSIP_bb->moreoutfile," \t->scen %4d (%3g ident.)\n", i_scen + 1, DDSIP_node[DDSIP_bb->curnode]->first_sol[scen][DDSIP_bb->firstvar]);
                        if (DDSIP_param->outlev > 5)
                        {
                            //Print each solution once (also inherited ones)
                            if (DDSIP_param->outlev >= DDSIP_first_stage_outlev && i_scen > scen)
                            {
                                for (j = 0; j < DDSIP_bb->firstvar; j++)
                                {
                                    fprintf (DDSIP_bb->moreoutfile," %21.14g", (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j]);
                                    if (!((j + 1) % 5))
                                        fprintf (DDSIP_bb->moreoutfile, "\n");
                                }
                                fprintf (DDSIP_bb->moreoutfile, "\n");
                            }
                        }
                    }
                }
                else
                {
                    if (DDSIP_param->outlev)
                    {
                        fprintf (DDSIP_bb->moreoutfile, "\n");
                    }
                    //DDSIP_Allocate space for the first stage variables
                    //  first_sol[DDSIP_bb->firstvar]   equals the number of identical first stage scenario solutions
                    //  first_sol[DDSIP_bb->firstvar+1] counts the levels of inheriting
                    //  first_sol[DDSIP_bb->firstvar+2] is the number of the node where the solution was computed
                    if (!(DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen])
                        (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen] =
                            (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar + 3, "DDSIP_node[curode]->first_sol[scen](LowerBound)");
                    //first_sol[DDSIP_bb->firstvar] equals the number of identical first stage scenario solutions
                    (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][DDSIP_bb->firstvar] = 1.0;
                    //first_sol[DDSIP_bb->firstvar+1] equals the levels of inheritance
                    (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][DDSIP_bb->firstvar + 1] = 0.0;
                    //first_sol[DDSIP_bb->firstvar+2] is the number of the node where the solution was computed
                    (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][DDSIP_bb->firstvar + 2] = DDSIP_bb->curnode;
                    if (DDSIP_param->outlev >= DDSIP_first_stage_outlev && i_scen >= scen)
                        fprintf (DDSIP_bb->moreoutfile, "    First-stage solution:\n");
                    for (j = 0; j < DDSIP_bb->firstvar; j++)
                    {
                        // Collect solution vector
                        if (DDSIP_bb->firsttype[j] == 'B' || DDSIP_bb->firsttype[j] == 'I' || DDSIP_bb->firsttype[j] == 'N')
                            (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j] = floor (mipx[DDSIP_bb->firstindex[j]] + 0.1);
                        else
                            (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j] = mipx[DDSIP_bb->firstindex[j]];
                        // Print something
                        if (DDSIP_param->outlev >= DDSIP_first_stage_outlev && i_scen >= scen)
                        {
                            fprintf (DDSIP_bb->moreoutfile," %21.14g", (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j]);
                            if (!((j + 1) % 5))
                                fprintf (DDSIP_bb->moreoutfile, "\n");
                        }
                        // Calculate minimum and maximum of each component
                        minfirst[j] = DDSIP_Dmin ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j], minfirst[j]);
                        maxfirst[j] = DDSIP_Dmax ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j], maxfirst[j]);
                        ///////////////////
                        // Check for violations of bounds by the CPLEX solution
                        if (DDSIP_param->outlev > 5)
                        {
                            for (k = 0; k< DDSIP_bb->curbdcnt; k++)
                            {
                                if (DDSIP_bb->curind[k] == j)
                                {
                                    if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j] < DDSIP_bb->curlb[k] && (DDSIP_bb->curlb[k]-(DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j])/(fabs(DDSIP_bb->curlb[k])+1e-9)> 1.e-10)
                                    {
                                        if (DDSIP_param->outlev)
                                            fprintf(DDSIP_bb->moreoutfile, "\nXXX WARNING: first-stage variable %d of solution to scenario %d with value %18.16g violates its lower bound %18.16g by %g\n\n",j+1,scen+1,(DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j], DDSIP_bb->curlb[k], DDSIP_bb->curlb[k]-(DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j]);
                                        fprintf(DDSIP_outfile, "XXX WARNING: first-stage variable %d of solution to scenario %d with value %18.16g violates its lower bound %18.16g by %g\n",j+1,scen+1,(DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j], DDSIP_bb->curlb[k], DDSIP_bb->curlb[k]-(DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j]);
                                    }
                                    if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j] > DDSIP_bb->curub[k] && ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j]-DDSIP_bb->curub[k])/(fabs(DDSIP_bb->curub[k])+1e-9) > 1e-10)
                                    {
                                        if (DDSIP_param->outlev)
                                            fprintf(DDSIP_bb->moreoutfile, "\nXXX WARNING: first-stage variable %d of solution to scenario %d with value %18.16g violates its upper bound %18.16g by %g\n\n",j+1,scen+1,(DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j], DDSIP_bb->curub[k], (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j]-DDSIP_bb->curub[k]);
                                        fprintf(DDSIP_outfile, "XXX WARNING: first-stage variable %d of solution to scenario %d with value %18.16g violates its upper bound %18.16g by %g\n",j+1,scen+1,(DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j], DDSIP_bb->curub[k], (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j]-DDSIP_bb->curub[k]);
                                    }
                                }
                            }
                        }
                        ///////////////////
                    }
                    if (DDSIP_param->outlev>= DDSIP_first_stage_outlev)
                        fprintf (DDSIP_bb->moreoutfile, "\n");
                    //
                    if ((DDSIP_param->addBendersCuts || DDSIP_param->addIntegerCuts) && (DDSIP_param->numberScenReeval > -1) && (DDSIP_bb->curnode < 9) && (j = DDSIP_bb->cutCntr))
                    {
                        // check feasibility for other scenarios and possibly add cuts
                        tmp = DDSIP_bb->sug[DDSIP_param->nodelim + 2];
                        if (!tmp)
                        {
                            // DDSIP_Allocate memory for heuristics suggestion if there is none yet
                            tmp = (sug_t *)DDSIP_Alloc (sizeof (sug_t), 1, "sug(LowerBound)");
                            tmp->firstval = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar, "sug->firstval(LowerBound)");
                            tmp->next = NULL;
                            DDSIP_bb->sug[DDSIP_param->nodelim + 2] = tmp;
                        }
                        else
                        {
                            if (tmp->next != NULL)
                            {
                                printf("\n\n FEHLER: !!!!!!!!!!!!!!!!!!!!! DDSIP_bb->sug[DDSIP_param->nodelim + 2]->next = %p\n\n",DDSIP_bb->sug[DDSIP_param->nodelim + 2]->next);
                                exit (1);
                            }
                        }

                        memcpy (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval, (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen], sizeof (double) * DDSIP_bb->firstvar);

                        if ((status = DDSIP_UpperBound (0, 2)))
                        {
                            fprintf (stderr, "ERROR: Failed to perform feasibility check\n");
                            return status;
                        }
                        nrFeasCheck++;
                        // Change bounds according to node in bb tree
                        if (DDSIP_bb->curbdcnt)
                        {
                            status = DDSIP_ChgBounds (0);
                            if (status)
                            {
                                fprintf (stderr, "ERROR: Failed to change problem \n");

                                if (DDSIP_param->outlev > 11)
                                    fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to change problem \n");
                                return status;
                            }
                        }
                        if (j < DDSIP_bb->cutCntr && nrFeasCheck <= DDSIP_param->numberScenReeval)
                        {
                            if(DDSIP_param->outlev)
                                fprintf (DDSIP_bb->moreoutfile, " --> scenario reevaluation %d\n", nrFeasCheck);
                            goto NEXT_TRY;
                        }
                    }
                    //
                    for (j = 0; j < DDSIP_bb->firstvar; j++)
                    {
                        minfirst[j] = DDSIP_Dmin ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j], minfirst[j]);
                        maxfirst[j] = DDSIP_Dmax ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j], maxfirst[j]);
                    }
                }
                if (DDSIP_param->outlev >= DDSIP_second_stage_outlev)
                    fprintf (DDSIP_bb->moreoutfile, "    Second-stage solution:\n");
                for (j = 0; j < DDSIP_bb->secvar; j++)
                {
                    DDSIP_bb->cur_secstage[scen][j] = mipx[DDSIP_bb->secondindex[j]];
                    if (DDSIP_param->outlev >= DDSIP_second_stage_outlev)
                    {
                        fprintf (DDSIP_bb->moreoutfile, " %21.14g",  mipx[DDSIP_bb->secondindex[j]]);
                        if (!((j + 1) % 5))
                            fprintf (DDSIP_bb->moreoutfile, "\n");
                    }
                }
                if (DDSIP_param->outlev >= DDSIP_second_stage_outlev)
                    fprintf (DDSIP_bb->moreoutfile, "\n");
            }
            // Temporary bound is infinity if a scenario problem is infeasible
            if (DDSIP_Equal ((DDSIP_node[DDSIP_bb->curnode]->subbound)[scen], DDSIP_infty))
                tmpbestbound = DDSIP_infty;
            else
            {
                tmpbestbound += DDSIP_data->prob[scen] * (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen];
            }
            tmpbestvalue += DDSIP_data->prob[scen] * (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen];

            if ((DDSIP_param->scalarization || DDSIP_param->cb) && (DDSIP_param->riskmod))
            {
                if (DDSIP_bb->ref_scenobj[scen] > DDSIP_param->risktarget)
                    DDSIP_bb->ref_risk[scen] = 1;
                else
                    DDSIP_bb->ref_risk[scen] = 0;
            }
            factor = (DDSIP_bb->bestvalue < 0.)? 1.-2.e-15 :  1.+2.e-15;
            // Node can be fathomed ?
            if (rest_bound > -DDSIP_infty /*&& DDSIP_param->prematureStops*/ && !DDSIP_param->riskalg && iscen < DDSIP_param->scenarios-1 && ((tmpbestbound + rest_bound) > (DDSIP_bb->bestvalue*factor)) && !DDSIP_param->scalarization)
            {
                DDSIP_node[DDSIP_bb->curnode]->bound = tmpbestbound + rest_bound;
                //
#ifdef DEBUG
                if (DDSIP_param->outlev > 21)
                {
                    fprintf(DDSIP_bb->moreoutfile," --premature stop with tmpbestbound + rest_bound: %22.14g + %22.14g = %22.14g\n", tmpbestbound, rest_bound, DDSIP_node[DDSIP_bb->curnode]->bound);
                    fprintf(DDSIP_bb->moreoutfile,"                       tmpbestbound + rest_bound: %22.14g > %22.14g = %22.14g*%g, difference: %g\n", DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->bestvalue*factor, DDSIP_bb->bestvalue, factor, DDSIP_node[DDSIP_bb->curnode]->bound-(DDSIP_bb->bestvalue*factor));
                }
#endif
                //
                DDSIP_bb->skip = 2;
                //DDSIP_bb->cutoff++;
                if (DDSIP_param->outlev)
                {
                    fprintf (DDSIP_bb->moreoutfile, "\tLower bound of node %d >= %18.16g (-bestvalue = %g) cutoff after evaluation of %d scenarios, skip=%d\n", DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bestvalue,iscen+1,DDSIP_bb->skip);
                }
                DDSIP_node[DDSIP_bb->curnode]->leaf = 1;
                goto TERMINATE;
            }
        }                // end else
    }                // end for iscen
    if ((j = DDSIP_bb->cutCntr - nrCuts))
        fprintf (DDSIP_outfile, " %6d%101d cuts\n", DDSIP_bb->curnode, j);

    if (DDSIP_bb->bestvalue < 0.)
    {
        factor  = 1.+DDSIP_Dmax(0.5*DDSIP_param->relgap, 2.e-10);
    }
    else
    {
        factor  = 1.-DDSIP_Dmax(0.5*DDSIP_param->relgap, 2.e-10);
    }

    if (!(DDSIP_bb->lb_sorted))
    {
        // in order to allow for premature cutoff: sort scenarios according to lower bound in root node in descending order
        double * sort_array;
        sort_array = (double *) DDSIP_Alloc(sizeof(double), DDSIP_param->scenarios, "sort_array(LowerBound)");
        for (iscen = 0; iscen < DDSIP_param->scenarios; iscen++)
        {
            sort_array[DDSIP_bb->lb_scen_order[iscen]] = DDSIP_data->prob[DDSIP_bb->lb_scen_order[iscen]] * (DDSIP_node[DDSIP_bb->curnode]->subbound)[DDSIP_bb->lb_scen_order[iscen]];
        }
        DDSIP_qsort_ins_D (sort_array, DDSIP_bb->lb_scen_order, 0, DDSIP_param->scenarios-1);

        if (DDSIP_param->outlev > 21)
        {
            // debug output
            fprintf (DDSIP_bb->moreoutfile,"order of scenarios after sorting lb order\n");
            for (iscen = 0; iscen < DDSIP_param->scenarios; iscen++)
            {
                fprintf(DDSIP_bb->moreoutfile," %3d: Scen %3d  %21.14g =  %10.07g * %21.14g\n",
                        iscen+1, DDSIP_bb->lb_scen_order[iscen]+1, sort_array[DDSIP_bb->lb_scen_order[iscen]],
                        DDSIP_data->prob[DDSIP_bb->lb_scen_order[iscen]],
                        (DDSIP_node[DDSIP_bb->curnode]->subbound)[DDSIP_bb->lb_scen_order[iscen]]);
            }
            // debug output
        }

        DDSIP_Free ((void**) &sort_array);

        // for the case of worst case cost risk:
        // change lower bound of the risk variable in root node to be the maximal lower bound
        if (!DDSIP_param->riskalg && !DDSIP_param->scalarization)
        {
            if (abs(DDSIP_param->riskmod) == 4)
            {
                worst_case_lb = -DDSIP_infty;
                if (DDSIP_param->riskmod > 0)
                {
                    for (i_scen=0; i_scen<DDSIP_param->scenarios; i_scen++)
                    {
                        bobjval = DDSIP_bb->btlb[i_scen] - DDSIP_param->riskweight*(DDSIP_node[DDSIP_bb->curnode]->first_sol)[i_scen][DDSIP_bb->firstvar-1];
                        if (bobjval > worst_case_lb)
                            worst_case_lb = bobjval;
                    }
                }
                else
                {
                    for (i_scen=0; i_scen<DDSIP_param->scenarios; i_scen++)
                    {
                        bobjval = DDSIP_bb->btlb[i_scen];
                        if (bobjval > worst_case_lb)
                            worst_case_lb = bobjval;
                    }
                }
                DDSIP_bb->lborg[DDSIP_bb->firstvar-1] = worst_case_lb;
            }
            else if (abs(DDSIP_param->riskmod) == 5)
            {
                double tmp;
                ordind = (int *) DDSIP_Alloc(sizeof (int), DDSIP_param->scenarios, "ordind(RiskObjective)");
                scensol = (double *) DDSIP_Alloc(sizeof (double), DDSIP_param->scenarios, "scensol(RiskObjective)");
                if (DDSIP_param->riskmod > 0)
                {
                    for (i_scen=0; i_scen<DDSIP_param->scenarios; i_scen++)
                    {
                        scensol[i_scen] = DDSIP_bb->btlb[i_scen] - DDSIP_param->riskweight* (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i_scen][DDSIP_bb->firstvar-1];
                    }
                }
                else
                {
                    for (i_scen=0; i_scen<DDSIP_param->scenarios; i_scen++)
                    {
                        scensol[i_scen] = DDSIP_bb->btlb[i_scen];
                    }
                }
                for (i_scen = 0; i_scen < DDSIP_param->scenarios; i_scen++)
                {
                    ordind[i_scen] = i_scen;
                }

                DDSIP_qsort_ins_A (scensol, ordind, 0, DDSIP_param->scenarios-1);

                sumprob = DDSIP_data->prob[ordind[0]];
                i_scen = 1;
                while (sumprob < DDSIP_param->risklevel && i_scen < DDSIP_param->scenarios)
                {
                    sumprob += DDSIP_data->prob[ordind[i_scen]];
                    i_scen++;
                }
                i_scen--;

                if (!DDSIP_param->riskalg && !DDSIP_param->scalarization)
                {
                    DDSIP_bb->lborg[DDSIP_bb->firstvar-1] = scensol[ordind[i_scen]] - 1.;
                    //    if (i_scen < DDSIP_param->scenarios -1)
                    //      i_scen++;
                    DDSIP_bb->uborg[DDSIP_bb->firstvar-1] = scensol[ordind[DDSIP_param->scenarios -1]] + DDSIP_Dmax(1,0.1*fabs(scensol[ordind[DDSIP_param->scenarios -1]]));
                    if (DDSIP_param->outlev)
                        fprintf (DDSIP_bb->moreoutfile,"---- changed bounds for auxiliary variable to [%g , %g]\n",DDSIP_bb->lborg[DDSIP_bb->firstvar-1],DDSIP_bb->uborg[DDSIP_bb->firstvar-1]);
                }
                if (DDSIP_param->cb && DDSIP_param->cbweight < (tmp = fabs(scensol[ordind[i_scen]])*0.025))
                {
                    fprintf (DDSIP_outfile, "\n     CBWEIGHT less than 0.025*(lb of TVaR), resetting to %g, CBFACTOR to %g.\n", tmp, 3e-1/tmp);
                    if (DDSIP_param->outlev)
                    {
                        printf ("\n     CBWEIGHT less than 0.025*(lb of TVaR), resetting to %g, CBFACTOR to %g.\n\n", tmp, 3e-1/tmp);
                        fprintf (DDSIP_bb->moreoutfile, "\n     CBWEIGHT less than 0.025*(lb of TVaR), resetting to %g, CBFACTOR to %g.\n", tmp, 3e-1/tmp);
                    }
                    DDSIP_param->cbfactor = 3e-1/tmp;
                    DDSIP_param->cbweight = DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_bb->dimdual] = DDSIP_bb->local_bestdual[DDSIP_bb->dimdual] = tmp;
                }
                DDSIP_Free ((void **) &ordind);
                DDSIP_Free ((void **) &scensol);
            }
        }
        DDSIP_bb->lb_sorted = 1;
    }
#ifdef SHIFT
    if (DDSIP_param->hot == 4)
    {
        // shifting difficult scenarios to the end such that they can make use of more mipstarts
        wall_secs = cpu_secs = DDSIP_bb->aggregate_time[0];
        for (wall_hrs = 1; wall_hrs < DDSIP_param->scenarios; wall_hrs++)
        {
            wall_secs += DDSIP_bb->aggregate_time[wall_hrs];
            if (wall_hrs < DDSIP_param->scenarios - 2)
                cpu_secs = DDSIP_Dmax (cpu_secs, DDSIP_bb->aggregate_time[wall_hrs]);
        }
        wall_secs /= (0.01 + DDSIP_param->scenarios);
#ifdef DEBUG
        if (DDSIP_param->outlev > 20)
            fprintf (DDSIP_bb->moreoutfile, " ### max time %g, mean %g, max>4*mean: %d ###\n", cpu_secs, wall_secs, cpu_secs > 4.*wall_secs);
#endif
        if (cpu_secs > 4.*wall_secs)
        {
            cpu_hrs = 0;
            for (wall_hrs = 0; wall_hrs < DDSIP_param->scenarios-cpu_hrs; wall_hrs++)
            {
                if (DDSIP_bb->aggregate_time[wall_hrs] > 4.*wall_secs)
                {
                    cpu_mins = DDSIP_bb->lb_scen_order[wall_hrs];
                    cpu_secs = DDSIP_bb->aggregate_time[wall_hrs];
#ifdef DEBUG
                    if (DDSIP_param->outlev > 20)
                        fprintf (DDSIP_bb->moreoutfile, " ### shifting scenario %d with time %g to the end of lb_scen_order ###\n", cpu_mins+1, DDSIP_bb->aggregate_time[wall_hrs]);
#endif
                    for (wall_mins = wall_hrs+1; wall_mins < DDSIP_param->scenarios; wall_mins++)
                    {
                        DDSIP_bb->lb_scen_order[wall_mins-1] = DDSIP_bb->lb_scen_order[wall_mins];
                        DDSIP_bb->aggregate_time[wall_mins-1] = DDSIP_bb->aggregate_time[wall_mins];
                    }
                    DDSIP_bb->lb_scen_order[DDSIP_param->scenarios-1] = cpu_mins;
                    DDSIP_bb->aggregate_time[DDSIP_param->scenarios-1] = cpu_secs;
                    cpu_hrs++;
                    wall_hrs--;
                }
            }
        }
    }
#endif

    if (DDSIP_param->riskmod && (DDSIP_param->riskalg || DDSIP_param->scalarization))
    {
        if (DDSIP_param->scalarization)
        {
            double tmp1 = DDSIP_RiskLb (DDSIP_node[DDSIP_bb->curnode]->ref_scenobj);
            if (DDSIP_param->outlev > 20)
                fprintf (DDSIP_bb->moreoutfile, " lb:  tmpbestbound = %g, DDSIP_RiskLb yields %g\n", tmpbestbound, tmp1);
            if (tmpbestbound < tmp1)
            {
                tmpbestbound = tmp1;
                if (DDSIP_param->outlev > 20)
                    fprintf (DDSIP_bb->moreoutfile, " lb:  new tmpbestbound = %g\n", tmpbestbound);
            }
        }
        else
        {
            double tmp1 = DDSIP_RiskLb (DDSIP_node[DDSIP_bb->curnode]->subbound);
            if (DDSIP_param->outlev > 20)
                fprintf (DDSIP_bb->moreoutfile, " lb:  tmpbestbound = %g, DDSIP_RiskLb using subbound yields %g, diff %g\n", tmpbestbound, tmp1, tmpbestbound - tmp1);
            tmpbestbound = tmp1;
        }
    }
    if (tmpbestbound >= DDSIP_node[DDSIP_bb->curnode]->bound)
    {
        // Set lower bound for the node
        DDSIP_node[DDSIP_bb->curnode]->bound = tmpbestbound;
    }
    else
    {
        if ((tmpbestbound - DDSIP_node[DDSIP_bb->curnode]->bound)/(fabs(tmpbestbound) + 1.e-12) > 1.e-6)
        {
            fprintf(DDSIP_outfile, "*** WARNING: current evaluation gave a worse bound (%16.10g) than the father (%16.10g), difference: %g\n",
                    tmpbestbound, DDSIP_node[DDSIP_bb->curnode]->bound, tmpbestbound - DDSIP_node[DDSIP_bb->curnode]->bound);
            if (DDSIP_param->outlev)
            {
                printf("*** WARNING: current evaluation gave a worse bound (%16.10g) than that of the father (%16.10g), difference: %g\n",
                       tmpbestbound, DDSIP_node[DDSIP_bb->curnode]->bound, tmpbestbound - DDSIP_node[DDSIP_bb->curnode]->bound);
                fprintf(DDSIP_bb->moreoutfile, "*** WARNING: current evaluation gave a worse bound (%16.10g) than the father (%16.10g), difference: %g\n",
                        tmpbestbound, DDSIP_node[DDSIP_bb->curnode]->bound, tmpbestbound - DDSIP_node[DDSIP_bb->curnode]->bound);
            }
        }
    }

    // should upper bounding be skipped in this node?
    if (!DDSIP_bb->skip)
        DDSIP_bb->skip = DDSIP_SkipUB ();
    // If none of the scenario problems has a solution then the dispersion norm (max-min) of each variable
    // is set to DDSIP_infty
    if (DDSIP_bb->skip)
    {
        if (DDSIP_param->outlev)
        {
            fprintf (DDSIP_bb->moreoutfile, "\tLower bound of node %d = %18.16g, skipping upper bound, skip=%d\n", DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound,DDSIP_bb->skip);
        }
        DDSIP_bb->violations = 0;
        maxdispersion = 0.;
        for (j = 0; j < DDSIP_bb->firstvar; j++)
        {
            if (!DDSIP_Equal (maxfirst[j], minfirst[j]))
            {
                DDSIP_bb->violations++;
                maxdispersion = DDSIP_Dmax (maxdispersion, maxfirst[j]-minfirst[j]);
            }
        }
        DDSIP_node[DDSIP_bb->curnode]->violations = DDSIP_bb->violations;
        DDSIP_node[DDSIP_bb->curnode]->dispnorm = maxdispersion;
        if (!DDSIP_bb->violations && DDSIP_node[DDSIP_bb->curnode]->bound < DDSIP_bb->bestvalue + DDSIP_param->accuracy)
        {
            DDSIP_node[DDSIP_bb->curnode]->leaf = 1;
            if (DDSIP_node[DDSIP_bb->curnode]->bound < DDSIP_bb->bestvalue)
            {
                DDSIP_bb->skip = -1;
            }
        }
        else
        {
            for (j = 0; j < DDSIP_bb->firstvar; j++)
            {
                maxfirst[j] = DDSIP_infty;
                minfirst[j] = -DDSIP_infty;
            }
        }
        if ((!(DDSIP_bb->found_optimal_node) && DDSIP_node[DDSIP_bb->curnode]->bound <= DDSIP_bb->bestvalue && DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bestvalue - fabs(DDSIP_bb->bestvalue)*DDSIP_Dmin(0.9*DDSIP_param->relgap, 1.e-8)) ||
                ( (DDSIP_bb->found_optimal_node) && DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bound_optimal_node))
        {
//////////////////////////////////////////////////////////////
            if (DDSIP_param->outlev > 20)
                fprintf (DDSIP_bb->moreoutfile, " ########## setting found_optimal_node: found_optimal_node= %d, DDSIP_node[%d]->bound (%20.15g) - bestvalue (%20.15g) = %.8g, - bestvalue*factor (%20.15g) = %.8g\n",
                         DDSIP_bb->found_optimal_node, DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->bestvalue, DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bestvalue,
                         DDSIP_bb->bestvalue*factor, DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bestvalue*factor);
//////////////////////////////////////////////////////////////
            DDSIP_bb->found_optimal_node = DDSIP_bb->curnode;
            DDSIP_bb->bound_optimal_node = DDSIP_node[DDSIP_bb->curnode]->bound;
        }
//////////////////////////////////////////////////////////////
        else
        {
            if (DDSIP_param->outlev > 20)
                fprintf (DDSIP_bb->moreoutfile, " ########## else: found_optimal_node= %d, DDSIP_node[%d]->bound (%20.15g) - bestvalue (%20.15g) = %.8g, - bestvalue*factor (%20.15g) = %.8g\n",
                         DDSIP_bb->found_optimal_node, DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->bestvalue, DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bestvalue,
                         DDSIP_bb->bestvalue*factor, DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bestvalue*factor);
        }
//////////////////////////////////////////////////////////////
    }
    else
    {
        if ((DDSIP_node[DDSIP_bb->curnode]->bound <  DDSIP_bb->bestvalue + DDSIP_param->accuracy) && (tmpbestvalue > DDSIP_bb->bestvalue + DDSIP_param->accuracy))
        {
            fprintf(DDSIP_outfile, "          WARNING: node %d is possibly not cut off just due to the MIP gaps. max MIP gap = %g%%, upper bound= %16.12g\n",
                    DDSIP_bb->curnode, maxGap, tmpbestvalue);
            if (DDSIP_param->outlev)
            {
                printf("          WARNING: node %d is possibly not cut off just due to the MIP gaps. max MIP gap = %g%%, upper bound= %16.12g\n", DDSIP_bb->curnode, maxGap, tmpbestvalue);
                fprintf(DDSIP_bb->moreoutfile, "          WARNING: node %d is possibly not cut off just due to the MIP gaps. max MIP gap = %g%%, upper bound= %16.12g\n",
                        DDSIP_bb->curnode, maxGap, tmpbestvalue);
            }
            DDSIP_bb->bestBound = 1;
        }
        // Count number of differences within first stage solution in current node
        // DDSIP_bb->violations=k means differences in k components
        if (DDSIP_param->outlev > 34)
        {
            colname = (char **) DDSIP_Alloc (sizeof (char *), (DDSIP_bb->firstvar + DDSIP_bb->secvar), "colname(LowerBound)");
            colstore = (char *) DDSIP_Alloc (sizeof (char), (DDSIP_bb->firstvar + DDSIP_bb->secvar) * DDSIP_ln_varname, "colstore(LowerBound)");
            status =
                CPXgetcolname (DDSIP_env, DDSIP_lp, colname, colstore,
                               (DDSIP_bb->firstvar + DDSIP_bb->secvar) * DDSIP_ln_varname, &j, 0, DDSIP_bb->firstvar + DDSIP_bb->secvar - 1);
            if (status)
            {
                fprintf (stderr, "ERROR: Failed to get column names (UB)\n");
                fprintf (DDSIP_outfile, "ERROR: Failed to get column names (UB)\n");
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to get column names (UB)\n");
                goto TERMINATE;
            }
        }
        DDSIP_bb->violations = 0;
        maxdispersion = 0.;
        for (j = 0; j < DDSIP_bb->firstvar; j++)
        {
            if (!DDSIP_Equal (maxfirst[j], minfirst[j]))
            {
                DDSIP_bb->violations++;
                maxdispersion = DDSIP_Dmax (maxdispersion, maxfirst[j]-minfirst[j]);
                // More Output
                if (DDSIP_param->outlev > 34)
                {
                    fprintf (DDSIP_bb->moreoutfile, "Deviation of variable %d: min=%16.14g\t max=%16.14g\t diff=%16.14g  \t%s\n", j, minfirst[j], maxfirst[j], maxfirst[j] - minfirst[j], colname[DDSIP_bb->firstindex[j]]);
                }
            }
        }
        if (DDSIP_param->outlev > 34)
        {
            DDSIP_Free ((void **) &(colstore));
            DDSIP_Free ((void **) &(colname));
        }
        DDSIP_node[DDSIP_bb->curnode]->dispnorm = maxdispersion;
        DDSIP_node[DDSIP_bb->curnode]->violations = DDSIP_bb->violations;
        DDSIP_bb->meanGapLB = DDSIP_Dmax(DDSIP_bb->meanGapLB, meanGap);
        if (!DDSIP_bb->violations)
        {
            DDSIP_node[DDSIP_bb->curnode]->leaf = 1;
            if (tmpbestvalue < DDSIP_bb->bestvalue)
            {
                DDSIP_bb->heurval = tmpbestvalue;
                DDSIP_bb->skip = -1;
            }
        }
        if ((!(DDSIP_bb->found_optimal_node) && DDSIP_node[DDSIP_bb->curnode]->bound <= DDSIP_bb->bestvalue && DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bestvalue - fabs(DDSIP_bb->bestvalue)*DDSIP_Dmin(0.9*DDSIP_param->relgap, 1.e-8)) ||
                ( (DDSIP_bb->found_optimal_node) && (DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bound_optimal_node)))
        {
//////////////////////////////////////////////////////////////
            if (DDSIP_param->outlev > 20)
                fprintf (DDSIP_bb->moreoutfile, " ########## setting found_optimal_node: found_optimal_node= %d, DDSIP_node[%d]->bound (%20.15g) - bestvalue (%20.15g) = %.8g, - bestvalue*factor (%20.15g) = %.8g\n",
                         DDSIP_bb->found_optimal_node, DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->bestvalue, DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bestvalue,
                         DDSIP_bb->bestvalue*factor, DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bestvalue*factor);
//////////////////////////////////////////////////////////////
            DDSIP_bb->found_optimal_node = DDSIP_bb->curnode;
            DDSIP_bb->bound_optimal_node = DDSIP_node[DDSIP_bb->curnode]->bound;
        }

        // More debugging information
        if (DDSIP_param->outlev)
        {
            if (DDSIP_bb->violations)
            {
                if (DDSIP_param->outlev > 3)
                {
                    fprintf (DDSIP_bb->moreoutfile,
                             "\tNumber of violations of nonanticipativity of first-stage variables:  %d, max: %g\n", DDSIP_bb->violations, maxdispersion);
                    if (DDSIP_param->outlev > 7)
                    {
                        printf ("\tNumber of violations of nonanticipativity of first-stage variables:  %d, max: %g\n", DDSIP_bb->violations, maxdispersion);
                        printf ("\tLower bound of node %d = %18.16g          \t(ub = %.16g,\tmean MIP gap: %g%%)\n", DDSIP_bb->curnode,
                                DDSIP_node[DDSIP_bb->curnode]->bound, tmpbestvalue, meanGap);
                    }
                }
            }
            else
            {
                if (DDSIP_param->outlev > 3)
                {
                    fprintf (DDSIP_bb->moreoutfile,
                             "\tNumber of violations of nonanticipativity of first-stage variables:  %d\n", DDSIP_bb->violations);
                    if (DDSIP_param->outlev > 7)
                    {
                        printf ("\tNumber of violations of nonanticipativity of first-stage variables:  %d\n", DDSIP_bb->violations);
                        printf ("\tLower bound of node %d = %18.16g          \t(ub = %.16g,\tmean MIP gap: %g%%)\n", DDSIP_bb->curnode,
                                DDSIP_node[DDSIP_bb->curnode]->bound, tmpbestvalue, meanGap);
                    }
                }
            }
            fprintf (DDSIP_bb->moreoutfile, "\tLower bound of node %d = %18.16g          \t(ub = %.16g,\tmean MIP gap: %g%%)\n", DDSIP_bb->curnode,
                     DDSIP_node[DDSIP_bb->curnode]->bound, tmpbestvalue, meanGap);
        }

        //Calculate Lagrangean relaxation part of objective
#ifdef CONIC_BUNDLE
        if (DDSIP_param->cb)
        {
            double lagr = 0.0;
            // update BoundNoLag if there were no CB iters before the current node
            if (DDSIP_bb->curnode < 0.5*DDSIP_param->cbBreakIters && (DDSIP_param->cb_depth < 0 || (!DDSIP_param->cb_depth && !DDSIP_param->cbrootitlim)))
            {
                memcpy (DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag, DDSIP_node[DDSIP_bb->curnode]->subbound, DDSIP_param->scenarios*sizeof(double));
                DDSIP_node[DDSIP_bb->curnode]->BoundNoLag = DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag[0] * DDSIP_data->prob[0];
                for (j = 1; j < DDSIP_param->scenarios; j++)
                    DDSIP_node[DDSIP_bb->curnode]->BoundNoLag += DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag[j] * DDSIP_data->prob[j];
            }
            else
            {
                // check whether all Lagrange parameters are zero
                for (i=0; i<DDSIP_bb->dimdual; i++)
                    if (DDSIP_node[DDSIP_bb->curnode]->dual[i])
                        break;
                if (i < DDSIP_bb->dimdual)
                {
                    // there are non-zero multipliers
                    if (DDSIP_param->outlev > 3)
                    {
                        for (scen = 0; scen < DDSIP_param->scenarios; scen++)
                            for (i = 0; i < DDSIP_bb->firstvar; i++)
                                for (j = DDSIP_data->nabeg[scen * DDSIP_bb->firstvar + i];
                                        j < DDSIP_data->nabeg[scen * DDSIP_bb->firstvar + i] + DDSIP_data->nacnt[scen * DDSIP_bb->firstvar + i]; j++)
                                {
                                    lagr += DDSIP_data->naval[j] * DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_data->naind[j]] * (DDSIP_node[DDSIP_bb->curnode]->first_sol[scen])[i];
                                }
                        fprintf (DDSIP_bb->moreoutfile,"        Objective value of solution = %.10g, Lagrangean part = %.10g\n",tmpbestvalue,lagr);
                        if (DDSIP_param->outlev > 29)
                        {
                            printf ("        Objective value of solution = %.10g, Lagrangean part = %.10g\n",tmpbestvalue,lagr);
                        }
                    }
                }
                else
                {
                    // the current bound is without influence of the Lagrangean
                    memcpy (DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag, DDSIP_node[DDSIP_bb->curnode]->subbound, DDSIP_param->scenarios*sizeof(double));
                    DDSIP_node[DDSIP_bb->curnode]->BoundNoLag = DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag[0] * DDSIP_data->prob[0];
                    for (j = 1; j < DDSIP_param->scenarios; j++)
                        DDSIP_node[DDSIP_bb->curnode]->BoundNoLag += DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag[j] * DDSIP_data->prob[j];
                }
            }
        }
#endif

        //DEBUGOUT
#ifdef MDEBUG
        if (DDSIP_param->cb && DDSIP_param->outlev> 49)
        {
            double expsubsol, expsubbound, expscenobj, risksubsol, risksubbound, riskscenobj;
            risksubsol = risksubbound = riskscenobj = expsubsol = expsubbound = expscenobj = 0.;
            fprintf(DDSIP_bb->moreoutfile, "\n############### LowerBound\n scen  cursubsol,   subbound,   ref_scenobj:\n");
            for(j=0; j<DDSIP_param->scenarios; j++)
            {
                fprintf(DDSIP_bb->moreoutfile, "%3d  %19.16g  %19.16g  %19.16g\n",j+1,(DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j],(DDSIP_node[DDSIP_bb->curnode]->subbound)[j],(DDSIP_node[DDSIP_bb->curnode]->ref_scenobj)[j]);
                expsubsol += DDSIP_data->prob[j] *(DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j];
                expsubbound += DDSIP_data->prob[j] *(DDSIP_node[DDSIP_bb->curnode]->subbound)[j];
                expscenobj += DDSIP_data->prob[j] *(DDSIP_node[DDSIP_bb->curnode]->ref_scenobj)[j];
                if (abs(DDSIP_param->riskmod) == 1)
                {
                    if ((DDSIP_node[DDSIP_bb->curnode]->subbound)[j] > DDSIP_param->risktarget)
                        risksubbound += DDSIP_data->prob[j] * (((DDSIP_node[DDSIP_bb->curnode]->subbound)[j] - DDSIP_param->risktarget));
                    if ((DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j] > DDSIP_param->risktarget)
                        risksubsol += DDSIP_data->prob[j] * (((DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j] - DDSIP_param->risktarget));
                    if ((DDSIP_node[DDSIP_bb->curnode]->ref_scenobj)[j] > DDSIP_param->risktarget)
                        riskscenobj += DDSIP_data->prob[j] * (((DDSIP_node[DDSIP_bb->curnode]->ref_scenobj)[j] - DDSIP_param->risktarget));
                }
                else if (abs(DDSIP_param->riskmod) == 2)
                {
                    if ((DDSIP_node[DDSIP_bb->curnode]->subbound)[j] > DDSIP_param->risktarget)
                        risksubbound += DDSIP_data->prob[j];
                    if ((DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j] > DDSIP_param->risktarget)
                        risksubsol += DDSIP_data->prob[j];
                    if ((DDSIP_node[DDSIP_bb->curnode]->ref_scenobj)[j] > DDSIP_param->risktarget)
                        riskscenobj += DDSIP_data->prob[j];
                }
            }
            fprintf(DDSIP_bb->moreoutfile, "exp  %19.16g  %19.16g  %19.16g\n",expsubsol,expsubbound,expscenobj);
            fprintf(DDSIP_bb->moreoutfile, "risk %19.16g  %19.16g  %19.16g\n",risksubsol,risksubbound,riskscenobj);
            fprintf(DDSIP_bb->moreoutfile, "sum  %19.16g  %19.16g  %19.16g\n",expsubsol + DDSIP_param->riskweight * risksubsol,expsubbound + DDSIP_param->riskweight * risksubbound, expscenobj + DDSIP_param->riskweight * riskscenobj);
            fprintf(DDSIP_bb->moreoutfile, "###############\n");
        }
#endif
        //DEBUGOUT

        if (!DDSIP_bb->violations && (!(DDSIP_param->riskmod) || (DDSIP_param->riskmod > 0 && DDSIP_param->riskalg == 1)))
        {
            if (relax < 2)
            {
                DDSIP_node[DDSIP_bb->curnode]->dispnorm = 0.0;
                if (DDSIP_param->scalarization)
                {
                    DDSIP_bb->curexp = 0.0;
                    for (j = 0; j < DDSIP_param->scenarios; j++)
                        DDSIP_bb->curexp += DDSIP_data->prob[j] * (DDSIP_node[DDSIP_bb->curnode]->ref_scenobj)[j];
                    // Calculate risk objective
                    DDSIP_RiskObjective (DDSIP_node[DDSIP_bb->curnode]->ref_scenobj);
                    we = (DDSIP_bb->curexp - DDSIP_param->ref_point[0]) / DDSIP_param->ref_scale[0];
                    wr = (DDSIP_bb->currisk - DDSIP_param->ref_point[1]) / DDSIP_param->ref_scale[1];
                    d = DDSIP_Dmax (we, wr);
                    tmpbestvalue = d + DDSIP_param->ref_eps * (we + wr);
                    if (DDSIP_param->outlev > 29)
                        fprintf (DDSIP_bb->moreoutfile, " Objective %.10g composition: we=%g, wr=%g, max(we,wr)=%.10g, eps*(we+wr)=%.10g\n",
                                 tmpbestvalue, we, wr, d, DDSIP_param->ref_eps * (we + wr));
                    DDSIP_bb->heurval = tmpbestvalue;
                }
                else
                {
                    DDSIP_bb->heurval = tmpbestvalue;

                    if (abs(DDSIP_param->riskmod) < 6 && DDSIP_param->riskalg != 1)
                    {
                        DDSIP_bb->heurval = DDSIP_Dmin (tmpbestvalue, DDSIP_bb->heurval);

                    }
                    if ((abs(DDSIP_param->riskmod) == 6 || abs(DDSIP_param->riskmod) == 7) && DDSIP_param->riskalg != 1)
                        DDSIP_node[DDSIP_bb->curnode]->bound = tmpbestvalue;
                }

                // if branch optimally solved - objective is exact lower (and upper) bound
                if (allopt)
                {
                    if (DDSIP_param->outlev)
                    {
                        fprintf (DDSIP_bb->moreoutfile, "\tLower and upper bound of node %d = %18.16g\n", DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound);
                        if (DDSIP_param->outlev > 7)
                            printf ("\tLower and upper bound of node %d = %18.16g\n", DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound);
                    }
                    DDSIP_bb->heurval = DDSIP_node[DDSIP_bb->curnode]->bound;
                    DDSIP_node[DDSIP_bb->curnode]->leaf = 1;
                    if (!DDSIP_param->cb)
                    {
                        if (DDSIP_bb->heurval >= DDSIP_bb->bestvalue)
                            DDSIP_bb->skip = 4;
                        else
                            DDSIP_bb->skip = -1;
                    }
                }
                // Debugging information
                if (DDSIP_param->outlev)
                {
                    if (DDSIP_param->riskmod && DDSIP_param->riskalg)
                        fprintf (DDSIP_bb->moreoutfile, "\tQ_E = %.10g \t Q_R = %.10g\n", DDSIP_bb->curexp, DDSIP_bb->currisk);
                    fprintf (DDSIP_bb->moreoutfile, "\tNew suggested upper bound = %18.16g\n", DDSIP_bb->heurval);
#ifdef DEBUG
                    if (DDSIP_param->outlev > 40)
                    {
                        char **colname;
                        char *colstore;
                        colname = (char **) DDSIP_Alloc (sizeof (char *), (DDSIP_bb->firstvar + DDSIP_bb->secvar), "colname(LowerBound)");
                        colstore = (char *) DDSIP_Alloc (sizeof (char), (DDSIP_bb->firstvar + DDSIP_bb->secvar) * DDSIP_ln_varname, "colstore(LowerBound)");
                        status =
                            CPXgetcolname (DDSIP_env, DDSIP_lp, colname, colstore,
                                           (DDSIP_bb->firstvar + DDSIP_bb->secvar) * DDSIP_ln_varname, &j, 0, DDSIP_bb->firstvar + DDSIP_bb->secvar - 1);
                        if (status)
                        {
                            fprintf (stderr, "ERROR: Failed to get column names (UB)\n");
                            fprintf (DDSIP_outfile, "ERROR: Failed to get column names (UB)\n");
                            if (DDSIP_param->outlev)
                                fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to get column names (UB)\n");
                            goto TERMINATE;
                        }

                        fprintf (DDSIP_bb->moreoutfile, "\nObjective function composition:\n");
                        for (j = 0; j < DDSIP_bb->firstvar + DDSIP_bb->secvar; j++)
                            if (fabs (DDSIP_bb->objcontrib[j]) > DDSIP_param->accuracy)
                                fprintf (DDSIP_bb->moreoutfile, "\t%s:  %.12g\n", colname[j], DDSIP_bb->objcontrib[j]);
                        DDSIP_Free ((void **) &(colstore));
                        DDSIP_Free ((void **) &(colname));
                    }
#endif
                }
                if (allopt && DDSIP_param->riskmod == 0)
                {
                    DDSIP_node[DDSIP_bb->curnode]->bound = tmpbestvalue;
                }
                if ((!DDSIP_param->riskmod || (DDSIP_param->riskmod && !DDSIP_param->riskalg)) && !DDSIP_param->cb)
                    DDSIP_bb->skip = -1;
                // Update if better solution was found
                if (DDSIP_bb->heurval < DDSIP_bb->bestvalue - 1.e-14)
                {
//////////////////////////////////////////////////////////////////
                    if (DDSIP_param->outlev > 21)
                    {
                        fprintf (DDSIP_bb->moreoutfile, "############ new suggested upper bound %18.12g is less than bestvalue %18.12g by %g, skip= %d #############\n",
                                 DDSIP_bb->heurval, DDSIP_bb->bestvalue, DDSIP_bb->bestvalue -  DDSIP_bb->heurval, DDSIP_bb->skip);
                    }
//////////////////////////////////////////////////////////////////
                    if (DDSIP_bb->skip == -1)
                    {
                        // in case the Lagrangean multipliers are not zero we have to recompute the scenario objective values without the Lagrangean part
                        for (i=0; i<DDSIP_bb->dimdual; i++)
                            if (DDSIP_node[DDSIP_bb->curnode]->dual[i])
                                break;
                        if (i < DDSIP_bb->dimdual)
                        {
                            // there are non-zero multipliers
                            // heuristics will all produce the (identical) first stage solution in this point - use these to produce the correct scenario
                            // objective values for the upper bound
                            //DDSIP_bb->skip=0;
//////////////////////////////////////////////////////////////////
                            if (DDSIP_param->outlev > 21)
                            {
                                fprintf (DDSIP_bb->moreoutfile, "############ new suggested upper bound %18.12g is less than bestvalue %18.12g, non-zero Lagr. mult. -> skip= %d #############\n",
                                         DDSIP_bb->heurval, DDSIP_bb->bestvalue, DDSIP_bb->skip);
                            }
//////////////////////////////////////////////////////////////////
                        }
                        else
                        {
                            // if the node has inherited scenario solutions, we have to solve these for the second stage variables
                            for(scen = 0; scen < DDSIP_param->scenarios; scen++)
                            {
                                if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][DDSIP_bb->firstvar + 2] < DDSIP_bb->curnode)
                                {
                                    nodes_1st = nodes_2nd = -1;
                                    // Output
                                    if (DDSIP_param->outlev)
                                    {
                                        fprintf (DDSIP_bb->moreoutfile,"*****************************************************************************\n");
                                        fprintf (DDSIP_bb->moreoutfile,"Solving scenario problem %d (solution inherited from node %.0f) for second stage variables in node %d....\n", scen + 1, (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][DDSIP_bb->firstvar + 2],  DDSIP_bb->curnode);
                                    }
                                    // Change problem according to scenario
                                    status = DDSIP_ChgProb (scen, DDSIP_bb->multipliers);
                                    if (status)
                                    {
                                        fprintf (stderr, "ERROR: Failed to change problem \n");
                                        fprintf (DDSIP_outfile, "ERROR: Failed to change problem \n");
                                        if (DDSIP_param->outlev)
                                            fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to change problem \n");
                                        goto TERMINATE;
                                    }

                                    // Warm starts
                                    for (k=0; k<DDSIP_param->scenarios; k++)
                                        if (DDSIP_bb->lb_scen_order[k] == scen)
                                            break;
                                    status = DDSIP_Warm (k);
                                    // Error?
                                    if (status)
                                        goto TERMINATE;

                                    status = DDSIP_LowerBoundWriteLp (scen);
                                    // Error?
                                    if (status)
                                        goto TERMINATE;
                                    time_start = DDSIP_GetCpuTime ();
                                    // Optimize
                                    if (relax < 2)
                                    {
#ifdef DEBUG
                                        // query time limit amd mip rel. gap parameters
                                        if (DDSIP_param->cpxscr || DDSIP_param->outlev > 21)
                                        {
                                            status = CPXgetdblparam (DDSIP_env,CPX_PARAM_TILIM,&we);
                                            status = CPXgetdblparam (DDSIP_env,CPX_PARAM_EPGAP,&wr);
                                            printf ("   -- 1st optimization time limit: %gs, rel. gap: %g%% --\n",we,wr*100.0);
                                        }
#endif
                                        // Optimize MIP
                                        DDSIP_bb->scenLBIters++;
                                        optstatus = CPXmipopt (DDSIP_env, DDSIP_lp);
                                        mipstatus = CPXgetstat (DDSIP_env, DDSIP_lp);
                                        if (!optstatus && !DDSIP_Error(optstatus) && !DDSIP_Infeasible (mipstatus))
                                        {
                                            if (CPXgetmiprelgap(DDSIP_env, DDSIP_lp, &mipgap))
                                            {
                                                fprintf (stderr, "ERROR: Query of CPLEX mip gap from 1st optimization failed (LowerBound) \n");
                                                fprintf (stderr, "       CPXgetstat returned: %d\n",status);
                                            }
#ifdef DEBUG
                                            time_lap = DDSIP_GetCpuTime ();
                                            nodes_1st = CPXgetnodecnt (DDSIP_env,DDSIP_lp);
                                            if (DDSIP_param->cpxscr || DDSIP_param->outlev > 11)
                                            {
                                                printf ("      LB: after 1st optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)\n",mipgap*100.0,nodes_1st,time_lap-time_start);
                                                if (DDSIP_param->outlev)
                                                    fprintf (DDSIP_bb->moreoutfile,"      LB: after 1st optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)\n",mipgap*100.0,nodes_1st,time_lap-time_start);
                                            }
#endif
                                            if (DDSIP_param->watchkappa)
                                            {
                                                double maxkappaval, stablekappaval, suspiciouskappaval, unstablekappaval, illposedkappaval;
                                                status = (CPXgetdblquality(DDSIP_env,DDSIP_lp,&maxkappaval,CPX_KAPPA_MAX) ||
                                                          CPXgetdblquality(DDSIP_env,DDSIP_lp,&illposedkappaval,CPX_KAPPA_ILLPOSED) ||
                                                          CPXgetdblquality(DDSIP_env,DDSIP_lp,&stablekappaval,CPX_KAPPA_STABLE) ||
                                                          CPXgetdblquality(DDSIP_env,DDSIP_lp,&suspiciouskappaval,CPX_KAPPA_SUSPICIOUS) ||
                                                          CPXgetdblquality(DDSIP_env,DDSIP_lp,&unstablekappaval,CPX_KAPPA_UNSTABLE));
                                                if (status)
                                                {
                                                    fprintf (stderr, "Failed to obtain Kappa information. \n\n");
                                                }
                                                else
                                                {
                                                    printf("maximal Kappa value        = %g\n",maxkappaval);
                                                    if (stablekappaval < 1.)
                                                    {
                                                        if (stablekappaval)
                                                            printf("numerically stable bases     %7.3f%%\n",stablekappaval*1e2);
                                                        if (suspiciouskappaval)
                                                            printf("numerically suspicious bases %7.3f%%\n",suspiciouskappaval*1e2);
                                                        if (unstablekappaval)
                                                            printf("numerically unstable bases   %7.3f%%\n",unstablekappaval*1e2);
                                                        if (illposedkappaval)
                                                            printf("numerically ill-posed bases  %7.3f%%\n",illposedkappaval*1e2);
                                                    }
                                                    if (DDSIP_param->outlev)
                                                    {
                                                        fprintf(DDSIP_bb->moreoutfile, "             maximal Kappa value        = %g\n",maxkappaval);
                                                        if (stablekappaval < 1.)
                                                        {
                                                            if (stablekappaval)
                                                                fprintf(DDSIP_bb->moreoutfile, "             numerically stable bases     %7.3f%%\n",stablekappaval*1e2);
                                                            if (suspiciouskappaval)
                                                                fprintf(DDSIP_bb->moreoutfile, "             numerically suspicious bases %7.3f%%\n",suspiciouskappaval*1e2);
                                                            if (unstablekappaval)
                                                                fprintf(DDSIP_bb->moreoutfile, "             numerically unstable bases   %7.3f%%\n",unstablekappaval*1e2);
                                                            if (illposedkappaval)
                                                                fprintf(DDSIP_bb->moreoutfile, "             numerically ill-posed bases  %7.3f%%\n",illposedkappaval*1e2);
                                                        }
                                                    }
                                                }
                                            }

                                            if (DDSIP_param->cpxnolb2 && mipstatus != CPXMIP_OPTIMAL)
                                            {
                                                // more iterations with different settings
                                                status = DDSIP_SetCpxPara (DDSIP_param->cpxnolb2, DDSIP_param->cpxlbisdbl2, DDSIP_param->cpxlbwhich2, DDSIP_param->cpxlbwhat2);
                                                if (status)
                                                {
                                                    fprintf (stderr, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
                                                    fprintf (DDSIP_outfile, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
                                                    if (DDSIP_param->outlev)
                                                        fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
                                                    goto TERMINATE;
                                                }
                                                // query time limit amd mip rel. gap parameters
                                                status = CPXgetdblparam (DDSIP_env,CPX_PARAM_EPGAP,&wr);
#ifdef DEBUG
                                                if (DDSIP_param->cpxscr || DDSIP_param->outlev > 21)
                                                {
                                                    status = CPXgetdblparam (DDSIP_env,CPX_PARAM_TILIM,&we);
                                                    printf ("   -- 2nd optimization time limit: %gs, rel. gap: %g%% --\n",we,wr*100.0);
                                                }
#endif
                                                // continue if desired gap is not reached yet
                                                if (mipgap > wr)
                                                {
                                                    if(DDSIP_param->files > 4)
                                                    {
                                                        char fname[DDSIP_ln_fname];
                                                        sprintf (fname, "%s/lb_sc%d_n%d%s", DDSIP_outdir, scen + 1, DDSIP_bb->curnode, "_2.prm");
                                                        status = CPXwriteparam(DDSIP_env,fname);
                                                        if (status)
                                                        {
                                                            fprintf (stderr, "ERROR: Failed to write parameters\n");
                                                            return status;
                                                        }
                                                        else if (DDSIP_param->outlev)
                                                        {
                                                            fprintf (DDSIP_bb->moreoutfile, "      parameter file %s written.\n",fname);
                                                        }
                                                    }
                                                    optstatus = CPXmipopt (DDSIP_env, DDSIP_lp);
                                                    nodes_2nd = CPXgetnodecnt (DDSIP_env,DDSIP_lp);
#ifdef DEBUG
                                                    if (DDSIP_param->cpxscr || DDSIP_param->outlev > 11)
                                                    {
                                                        if (CPXgetmiprelgap(DDSIP_env, DDSIP_lp, &mipgap))
                                                        {
                                                            fprintf (stderr, "ERROR: Query of CPLEX mip gap from 2nd optimization failed (LowerBound) \n");
                                                            fprintf (stderr, "       CPXgetstat returned: %d\n",status);
                                                        }
                                                        else
                                                        {
                                                            time_end = DDSIP_GetCpuTime ();
                                                            printf ("      LB: after 2nd optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)\n",mipgap*100.0,nodes_2nd,time_end-time_lap);
                                                            if (DDSIP_param->outlev)
                                                                fprintf (DDSIP_bb->moreoutfile,"      LB: after 2nd optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)\n",mipgap*100.0,nodes_2nd,time_end-time_lap);
                                                        }
                                                    }
#endif
                                                }
                                                if (DDSIP_param->watchkappa)
                                                {
                                                    double maxkappaval, stablekappaval, suspiciouskappaval, unstablekappaval, illposedkappaval;
                                                    status = (CPXgetdblquality(DDSIP_env,DDSIP_lp,&maxkappaval,CPX_KAPPA_MAX) ||
                                                              CPXgetdblquality(DDSIP_env,DDSIP_lp,&illposedkappaval,CPX_KAPPA_ILLPOSED) ||
                                                              CPXgetdblquality(DDSIP_env,DDSIP_lp,&stablekappaval,CPX_KAPPA_STABLE) ||
                                                              CPXgetdblquality(DDSIP_env,DDSIP_lp,&suspiciouskappaval,CPX_KAPPA_SUSPICIOUS) ||
                                                              CPXgetdblquality(DDSIP_env,DDSIP_lp,&unstablekappaval,CPX_KAPPA_UNSTABLE));
                                                    if (status)
                                                    {
                                                        fprintf (stderr, "Failed to obtain Kappa information. \n\n");
                                                    }
                                                    else
                                                    {
                                                        printf("maximal Kappa value        = %g\n",maxkappaval);
                                                        if (stablekappaval < 1.)
                                                        {
                                                            if (stablekappaval)
                                                                printf("numerically stable bases     %7.3f%%\n",stablekappaval*1e2);
                                                            if (suspiciouskappaval)
                                                                printf("numerically suspicious bases %7.3f%%\n",suspiciouskappaval*1e2);
                                                            if (unstablekappaval)
                                                                printf("numerically unstable bases   %7.3f%%\n",unstablekappaval*1e2);
                                                            if (illposedkappaval)
                                                                printf("numerically ill-posed bases  %7.3f%%\n",illposedkappaval*1e2);
                                                        }
                                                        if (DDSIP_param->outlev)
                                                        {
                                                            fprintf(DDSIP_bb->moreoutfile, "             maximal Kappa value        = %g\n",maxkappaval);
                                                            if (stablekappaval < 1.)
                                                            {
                                                                if (stablekappaval)
                                                                    fprintf(DDSIP_bb->moreoutfile, "             numerically stable bases     %7.3f%%\n",stablekappaval*1e2);
                                                                if (suspiciouskappaval)
                                                                    fprintf(DDSIP_bb->moreoutfile, "             numerically suspicious bases %7.3f%%\n",suspiciouskappaval*1e2);
                                                                if (unstablekappaval)
                                                                    fprintf(DDSIP_bb->moreoutfile, "             numerically unstable bases   %7.3f%%\n",unstablekappaval*1e2);
                                                                if (illposedkappaval)
                                                                    fprintf(DDSIP_bb->moreoutfile, "             numerically ill-posed bases  %7.3f%%\n",illposedkappaval*1e2);
                                                            }
                                                        }
                                                    }
                                                }
                                                // reset CPLEX parameters
                                                //  resetting CPLEX parameters to lb (LowerBound)
                                                status = DDSIP_SetCpxPara (DDSIP_param->cpxnolb, DDSIP_param->cpxlbisdbl, DDSIP_param->cpxlbwhich, DDSIP_param->cpxlbwhat);
                                                if (status)
                                                {
                                                    fprintf (stderr, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
                                                    fprintf (DDSIP_outfile, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
                                                    if (DDSIP_param->outlev)
                                                        fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
                                                    goto TERMINATE;
                                                }
                                            }
                                        }

                                        if ((k = CPXgetnummipstarts(DDSIP_env, DDSIP_lp)) > 2)
                                        {
                                            status    = CPXdelmipstarts (DDSIP_env, DDSIP_lp, 1, k-1);
                                        }
                                    }
                                    else
                                        // Optimize relaxed MIP
                                        optstatus = CPXdualopt (DDSIP_env, DDSIP_lp);
                                    // Infeasible, unbounded ... ?
                                    // We handle some errors separately (blatant infeasible, error in scenario problem)
                                    if (DDSIP_Error (optstatus))
                                    {
                                        fprintf (stderr, "ERROR: Failed to optimize problem for scenario %d.(LowerBound), status=%d\n", scen + 1, optstatus);
                                        fprintf (DDSIP_outfile, "ERROR: Failed to optimize problem for scenario %d.(LowerBound), status=%d\n", scen + 1, optstatus);
                                        if (DDSIP_param->outlev)
                                            fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to optimize problem for scenario %d.(LowerBound) status=%d\n",
                                                     scen + 1, optstatus);
                                        if (optstatus == CPXERR_SUBPROB_SOLVE)
                                        {
                                            optstatus = CPXgetsubstat(DDSIP_env, DDSIP_lp);
                                            fprintf (stderr, "ERROR:                                scenario %d.(LowerBound) subproblem status=%d\n", scen + 1, optstatus);
                                            fprintf (DDSIP_outfile, "ERROR:                                scenario %d.(LowerBound) subproblem status=%d\n", scen + 1, optstatus);
                                            if (DDSIP_param->outlev)
                                                fprintf (DDSIP_bb->moreoutfile, "ERROR:                                scenario %d.(LowerBound) subproblem status=%d\n", scen + 1, optstatus);
                                        }
                                        status = optstatus;
                                        goto TERMINATE;
                                    }

                                    mipstatus = CPXgetstat (DDSIP_env, DDSIP_lp);
                                    time (&DDSIP_bb->cur_time);
                                    // Infeasible, unbounded ... ?
                                    if (mipstatus==CPXMIP_TIME_LIM_INFEAS)
                                    {
                                        if (DDSIP_param->outlev)
                                        {
                                            DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                                            DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                                            fprintf (DDSIP_bb->moreoutfile,
                                                     "%4d Scenario %4.0d:                          \tinfeasible               \t         \tStat.%3.0d\t%3dh %02d:%02.0f / %3dh %02d:%05.2f\n",
                                                     scen + 1, scen + 1, mipstatus, wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs);
                                        }
                                        DDSIP_node[DDSIP_bb->curnode]->bound = DDSIP_infty;
                                        CPXgetdblparam (DDSIP_env,CPX_PARAM_TILIM,&cpu_secs);
                                        if (DDSIP_param->outlev)
                                        {
                                            printf ("ERROR: Problem infeasible within time limit of %g secs for scenario %d (LowerBound)\n", cpu_secs, scen + 1);
                                            fprintf (DDSIP_bb->moreoutfile, "ERROR: Problem infeasible within time limit of %g secs for scenario %d (LowerBound)\n", cpu_secs, scen + 1);
                                        }
                                        DDSIP_bb->skip = -2;
                                        for (k = scen; k < DDSIP_param->scenarios; k++)
                                            DDSIP_bb->solstat[k] = 0;
                                        bobjval = DDSIP_infty;
                                        status = mipstatus;
                                        goto TERMINATE;
                                    }
                                    else if (DDSIP_Infeasible (mipstatus))
                                    {
                                        if (DDSIP_param->outlev)
                                        {
                                            DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                                            DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                                            fprintf (DDSIP_bb->moreoutfile,
                                                     "%4d Scenario %4.0d:                          \tinfeasible               \t         \tStat.%3.0d\t%3dh %02d:%02.0f / %3dh %02d:%05.2f\n",
                                                     scen + 1, scen + 1, mipstatus, wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs);
                                        }
                                        DDSIP_node[DDSIP_bb->curnode]->bound = DDSIP_infty;
                                        DDSIP_bb->skip = 1;
                                        for (k = scen; k < DDSIP_param->scenarios; k++)
                                            DDSIP_bb->solstat[k] = 0;
                                        goto TERMINATE;
                                    }
                                    // We did not detect infeasibility
                                    // If we couldn't find a feasible solution we can at least obtain a lower bound
                                    if (DDSIP_NoSolution (mipstatus))
                                    {
                                        objval = DDSIP_infty;
                                        DDSIP_bb->solstat[scen] = 0;
                                        if (relax != 2)
                                            status = CPXgetbestobjval (DDSIP_env, DDSIP_lp, &bobjval);
                                        else
                                            bobjval = DDSIP_infty;
                                        if (status)
                                        {
                                            fprintf (stderr, "ERROR: Failed to get value of best remaining node (LowerBound)\n");
                                            fprintf (DDSIP_outfile, "ERROR: Failed to get value of best remaining node (LowerBound)\n");
                                            if (DDSIP_param->outlev)
                                                fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to get value of best remaining node (LowerBound)\n");
                                            goto TERMINATE;
                                        }
                                    }
                                    else              // Solution exists
                                    {
                                        DDSIP_bb->solstat[scen] = 1;
                                        status = CPXgetobjval (DDSIP_env, DDSIP_lp, &objval);
                                        if (status)
                                        {
                                            fprintf (stderr, "ERROR*: Failed to get best objective value \n");
                                            fprintf (DDSIP_outfile, "ERROR*: Failed to get best objective value \n");
                                            if (DDSIP_param->outlev)
                                                fprintf (DDSIP_bb->moreoutfile, "ERROR*: Failed to get best objective value \n");
                                            goto TERMINATE;
                                        }

                                        status = CPXgetx (DDSIP_env, DDSIP_lp, mipx, 0, DDSIP_bb->firstvar + DDSIP_bb->secvar - 1);
                                        if (status)
                                        {
                                            fprintf (stderr, "ERROR: Failed to get solution \n");
                                            fprintf (DDSIP_outfile, "ERROR: Failed to get solution \n");
                                            if (DDSIP_param->outlev)
                                                fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to get solution \n");
                                            goto TERMINATE;
                                        }
                                        // Numerical errors?
                                        for (j = 0; j < DDSIP_bb->firstvar; j++)
                                            if (fabs (mipx[DDSIP_bb->firstindex[j]]) < DDSIP_param->accuracy)
                                                mipx[DDSIP_bb->firstindex[j]] = 0.0;
                                        if (relax != 2)
                                        {
                                            //if (mipstatus == CPXMIP_OPTIMAL)
                                            //{
                                            //    bobjval = objval;
                                            //}
                                            //else
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
                                        }
                                        else
                                            bobjval = DDSIP_infty;

                                    }                // end for else

                                    (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen] = objval;
                                    // If all scenario problems were solved to optimality ....
                                    if (!(mipstatus == CPXMIP_OPTIMAL))
                                    {
                                        allopt = 0;
                                        // If the tree was exhausted bobjval is huge
                                        if (fabs (bobjval) > DDSIP_infty)
                                            bobjval = objval;
                                    }
                                    (DDSIP_node[DDSIP_bb->curnode]->mipstatus)[scen] = mipstatus;
                                    gap = 100.0*(objval-bobjval)/(fabs(objval)+1e-4);
                                    if (fabs(gap) < 2.e-13)
                                        gap = 0.;
                                    meanGap += DDSIP_data->prob[scen] * gap;
                                    time (&DDSIP_bb->cur_time);
                                    // Debugging information
                                    if (DDSIP_param->outlev)
                                    {
                                        time_end = DDSIP_GetCpuTime ();
                                        time_help = time_end-time_start;
                                        DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                                        DDSIP_translate_time (time_end,&cpu_hrs,&cpu_mins,&cpu_secs);
                                        if (mipstatus == CPXMIP_OPTIMAL)
                                            fprintf (DDSIP_bb->moreoutfile,
                                                     "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%)     ",
                                                     iscen + 1, scen + 1, objval, bobjval, gap);
                                        else if (mipstatus == CPXMIP_OPTIMAL_TOL)
                                            fprintf (DDSIP_bb->moreoutfile,
                                                     "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) tol.",
                                                     iscen + 1, scen + 1, objval, bobjval, gap);
                                        else if (mipstatus == CPXMIP_NODE_LIM_FEAS)
                                            fprintf (DDSIP_bb->moreoutfile,
                                                     "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) NODE",
                                                     iscen + 1, scen + 1, objval, bobjval, gap);
                                        else if (mipstatus == CPXMIP_TIME_LIM_FEAS)
                                            fprintf (DDSIP_bb->moreoutfile,
                                                     "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) TIME",
                                                     iscen + 1, scen + 1, objval, bobjval, gap);
                                        else
                                            fprintf (DDSIP_bb->moreoutfile,
                                                     "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) %-4.0d ",
                                                     iscen + 1, scen + 1, objval, bobjval, gap, mipstatus);
                                        fprintf (DDSIP_bb->moreoutfile,
                                                 "\t %3dh %02d:%02.0f / %3dh %02d:%05.2f (%7.2fs n:%5d",
                                                 wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs, time_help, nodes_1st);
                                        if (nodes_2nd >= 0)
                                            fprintf (DDSIP_bb->moreoutfile, " +%4d)", nodes_2nd - nodes_1st);
                                        else
                                            fprintf (DDSIP_bb->moreoutfile, ")");
                                    }

                                    (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen] = DDSIP_Dmin (bobjval, objval);

                                    // If a feasible solution exists - compare with the one we had already
                                    if (DDSIP_bb->solstat[scen])
                                    {
                                        for (j = 0; j < DDSIP_bb->firstvar; j++)
                                        {
                                            if (!DDSIP_Equal (mipx[DDSIP_bb->firstindex[j]], (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j]))
                                                break;
                                        }
                                        if (j < DDSIP_bb->firstvar)
                                        {
                                            fprintf (stderr," ERROR: *** resolve of inherited scenario solution %d for second stage variables gave a different first stage.\n",scen);
                                            fprintf (DDSIP_outfile," ERROR: *** resolve of inherited scenario solution %d for second stage variables gave a different first stage.\n",scen);
                                            if (DDSIP_param->outlev)
                                                fprintf (DDSIP_bb->moreoutfile," ERROR: *** resolve of inherited scenario solution %d for second stage variables gave a different first stage.\n",scen);
                                            goto TERMINATE;
                                        }
                                    }
                                    if (DDSIP_param->outlev >= DDSIP_second_stage_outlev)
                                        fprintf (DDSIP_bb->moreoutfile, "    Second-stage solution:\n");
                                    for (j = 0; j < DDSIP_bb->secvar; j++)
                                    {
                                        DDSIP_bb->cur_secstage[scen][j] = mipx[DDSIP_bb->secondindex[j]];
                                        if (DDSIP_param->outlev >= DDSIP_second_stage_outlev)
                                        {
                                            fprintf (DDSIP_bb->moreoutfile, " %21.14g",  mipx[DDSIP_bb->secondindex[j]]);
                                            if (!((j + 1) % 5))
                                                fprintf (DDSIP_bb->moreoutfile, "\n");
                                        }
                                    }
                                    if (DDSIP_param->outlev >= DDSIP_second_stage_outlev)
                                        fprintf (DDSIP_bb->moreoutfile, "\n");
                                }
                            }
                            DDSIP_bb->bestvalue = tmpbestvalue;
                            DDSIP_bb->feasbound = tmpfeasbound;
                            if (!DDSIP_param->riskmod)
                            {
                                DDSIP_bb->bestexp = DDSIP_bb->bestvalue;
                            }
                            else if (DDSIP_param->riskmod < 0)
                            {
                                DDSIP_bb->bestrisk = DDSIP_bb->bestvalue;
                            }
                            // We store all risk measures
                            if (DDSIP_param->riskmod)
                                for (j = 0; j < DDSIP_maxrisk; j++)
                                    DDSIP_bb->bestriskval[j] = DDSIP_bb->curriskval[j];
                            for (j = 0; j < DDSIP_bb->firstvar; j++)
                                DDSIP_bb->bestsol[j] = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[0][j];
                            if (DDSIP_param->outlev)
                                fprintf (DDSIP_bb->moreoutfile, "\t(Current best bound, improvement %16.12g, %g%%)\n",DDSIP_bb->bestvalue-tmpbestvalue,1e+2*(DDSIP_bb->bestvalue-tmpbestvalue)/(fabs(DDSIP_bb->bestvalue) + 1.e-16));
                            // Print a line of output
                            DDSIP_PrintStateUB (0);
                            // copy second stage solutions of current best
                            for (k = 0; k < DDSIP_param->scenarios; k++)
                            {
                                memcpy(DDSIP_bb->secstage[k], DDSIP_bb->cur_secstage[k], DDSIP_bb->secvar*sizeof(double));
                                DDSIP_bb->subsol[k] = (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[k];
                            }
                            // heuristics will all produce the (identical) first stage solution in this point - skip upper bound
                            DDSIP_bb->skip = 4;
                            DDSIP_bb->bestsol_in_curnode = 1;
                        }
                    }
                }
            }
            //if all scenario solutions were inherited (should not happen) we would get the same heuristics again, thus we can skip upper bound
            if (!probs_solved)
            {
                DDSIP_bb->skip = 1;
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile, "\tLower bound of node %d = %18.16g, skipping upper bound, skip=%d\n", DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound,DDSIP_bb->skip);
            }
        }
    }

    // Retrieve branching index and branching value
    if (DDSIP_node[DDSIP_bb->curnode]->dispnorm < DDSIP_param->nulldisp)
    {
        DDSIP_node[DDSIP_bb->curnode]->leaf = 1;
        //
#ifdef DEBUG
        if (DDSIP_param->outlev > 29)
        {
            fprintf(DDSIP_bb->moreoutfile," dispersion in node %d less than nulldisp, make it a leaf\n", DDSIP_bb->curnode);
        }
#endif
        //
    }
    else if (!DDSIP_node[DDSIP_bb->curnode]->leaf && DDSIP_bb->skip != 2)
    {
        for (j = 0; j < DDSIP_bb->firstvar; j++)
            maxfirst[j] -= minfirst[j];
        status = DDSIP_GetBranchIndex (maxfirst);
    }

TERMINATE:

    // Only if no errors, yet.
    if (!status)
    {
        status = DDSIP_RestoreBoundAndType ();
        if (status)
            fprintf (stderr, "ERROR: Failed to restore bounds \n");
    }
    DDSIP_Free ((void **) &(indices));
    DDSIP_Free ((void **) &(mipx));
    DDSIP_Free ((void **) &(minfirst));
    DDSIP_Free ((void **) &(maxfirst));
    DDSIP_Free ((void **) &(type));
    return status;
} // DDSIP_LowerBound




//==========================================================================
// Function solves scenario problems in bundle method for scalarization 0 (weighted sum)
int
DDSIP_CBLowerBound (double *objective_val, double relprec)
{
#ifdef CONIC_BUNDLE
    double objval, bobjval, tmpbestbound = 0.0, tmpupper = 0.0, maxdispersion = 0.;
    double wr, mipgap, time_start, time_end, time_help, wall_secs, cpu_secs, gap, meanGap, maxGap;
#ifdef DEBUG
    double we, time_lap;
#endif
    //double cplexRelGap = 1.e-6;

    int cnt, iscen, i_scen, j, k, k1, status = 0, optstatus, mipstatus, scen, relax = 0, increase = 9;
    int wall_hrs, wall_mins,cpu_hrs, cpu_mins;
    int nodes_1st, nodes_2nd;
    static int nearly_constant = 0, totally_constant = 0;
    static double original_weight = 1., old_value = -1.e30;
    static double DualObj;
    double weight_reset_factor = 1.8;
    static int use_LB_params = 0, no_viol_worse = 0;

    int *indices = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_bb->firstvar + DDSIP_bb->secvar),
                                        "indices (CBLowerBound)");
    double *mipx = (double *) DDSIP_Alloc (sizeof (double), (DDSIP_bb->firstvar + DDSIP_bb->secvar),
                                           "mipx (CBLowerBound)");
    double *minfirst = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar,
                       "minfirst (CBLowerBound)");
    double *maxfirst = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar,
                       "maxfirst (CBLowerBound)");
    //    double *cb_scenrisk = (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios,
    //                                          "cb_scenrisk (CBLowerBound)");
    char *type = (char *) DDSIP_Alloc (sizeof (char), DDSIP_bb->firstvar, "type(CBLowerBound)");
    char **colname;
    char *colstore;
    int shift_in_cb = 0;
    int nrFeasCheck = 0;
    sug_t *tmp;
#ifdef CBHOTSTART
#ifdef SHIFT
    time_t   startTime;
    shift_in_cb = (DDSIP_param->hot == 4 && DDSIP_bb->dualdescitcnt) ||
                  (DDSIP_param->cbhot == 4 && DDSIP_bb->dualdescitcnt) ||
                  (DDSIP_param->cbhot == 3 && !DDSIP_bb->curnode && DDSIP_bb->dualdescitcnt) ||
                  (DDSIP_param->cbhot == 2 && !DDSIP_bb->curnode && (DDSIP_bb->dualdescitcnt == 1 || DDSIP_bb->dualdescitcnt > 5)) ||
                  (DDSIP_param->cbhot == 1 && !DDSIP_bb->curnode && DDSIP_bb->dualdescitcnt <= 3 && DDSIP_bb->dualdescitcnt >= 1);
#endif
#endif
#ifdef SHIFT
    if (!(DDSIP_bb->dualitcnt))
    {
        memset (DDSIP_bb->aggregate_time, '\0', sizeof (double) * (DDSIP_param->scenarios));
    }
#endif

    DDSIP_bb->DDSIP_step = dual;
    DDSIP_bb->skip = 0;
    if (DDSIP_bb->dualitcnt == 1)
    {
        DDSIP_bb->weight_reset = 0;
    }

    // ConicBundle somtimes gave negative relprec
    relprec = fabs(relprec);

    // Insurance in case CB doesn't stop as intended
    if (DDSIP_bb->dualitcnt > DDSIP_param->cbtotalitlim + 100)
    {
        if (DDSIP_param->outlev > 11)
            fprintf (DDSIP_bb->moreoutfile, "######### dualitcnt= %d > %d (cbtotalitlim + 100)\n", DDSIP_bb->dualitcnt, DDSIP_param->cbtotalitlim + 100);
        return 1;
    }
    // Output
    if (DDSIP_param->outlev)
    {
        printf ("\n----------------------\nSolving (dual opt. function eval) in node %d, desc. it. %d (limit: %d), tot. it. %3d (curr. desc.:%2d) weight = %g :\n", DDSIP_bb->curnode, DDSIP_bb->dualdescitcnt,DDSIP_bb->current_itlim,DDSIP_bb->dualitcnt+1, DDSIP_bb->dualitcnt - DDSIP_bb->last_dualitcnt + 1, cb_get_last_weight(DDSIP_bb->dualProblem));
        fprintf (DDSIP_bb->moreoutfile, "\n----------------------\nSolving (dual opt. function eval) in node %d, desc. it. %d (limit: %d), tot. it. %3d (curr. desc.:%2d) weight = %g :\n", DDSIP_bb->curnode, DDSIP_bb->dualdescitcnt,DDSIP_bb->current_itlim,DDSIP_bb->dualitcnt+1, DDSIP_bb->dualitcnt - DDSIP_bb->last_dualitcnt + 1, cb_get_last_weight(DDSIP_bb->dualProblem));
    }
    // Initialization of indices, minfirst, and maxfirst
    for (j = 0; j < DDSIP_bb->firstvar + DDSIP_bb->secvar; j++)
        indices[j] = j;
    for (j = 0; j < DDSIP_bb->firstvar; j++)
    {
        minfirst[j] = DDSIP_infty;
        maxfirst[j] = -DDSIP_infty;
    }

    // CPLEX parameters
    if (use_LB_params)
    {
        status = DDSIP_SetCpxPara (DDSIP_param->cpxnolb, DDSIP_param->cpxlbisdbl, DDSIP_param->cpxlbwhich, DDSIP_param->cpxlbwhat);
    }
    else
    {
        status = DDSIP_SetCpxPara (DDSIP_param->cpxnodual, DDSIP_param->cpxdualisdbl, DDSIP_param->cpxdualwhich, DDSIP_param->cpxdualwhat);
    }
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to set CPLEX parameters (CBLowerBound) \n");
        fprintf (DDSIP_outfile, "ERROR: Failed to set CPLEX parameters (CBLowerBound) \n");
        if (DDSIP_param->outlev)
            fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to set CPLEX parameters (CBLowerBound) \n");
        goto TERMINATE;
    }

    // assure that CPLEX relative gap is less or equal the relprec requested by ConicBundle
    //   if ((status = CPXgetdblparam (DDSIP_env, CPX_PARAM_EPGAP, &cplexRelGap))) {
    //     fprintf (stderr, "ERROR: Failed to get CPLEX parameter CPX_PARAM_EPGAP (CBLowerBound) \n");
    //     fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to get CPLEX parameter CPX_PARAM_EPGAP (CBLowerBound) \n");
    //   }
    //   else {
    //     if (cplexRelGap > relprec) {
    //       if ((status = CPXsetdblparam (DDSIP_env, CPX_PARAM_EPGAP, relprec))) {
    //         fprintf (stderr, "ERROR: Failed to set CPLEX parameter CPX_PARAM_EPGAP to %g (CBLowerBound) \n", relprec);
    //         if (DDSIP_param->outlev)
    //           fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to set CPLEX parameter CPX_PARAM_EPGAP to %g (CBLowerBound) \n", relprec);
    //       }
    //       else {
    //         fprintf (stderr, " reset CPLEX parameter CPX_PARAM_EPGAP from %g to relprec=%g (CBLowerBound) \n", cplexRelGap, relprec);
    //         if (DDSIP_param->outlev)
    //           fprintf (DDSIP_bb->moreoutfile, " reset CPLEX parameter CPX_PARAM_EPGAP from %g to relprec=%g (CBLowerBound) \n", cplexRelGap, relprec);
    //       }
    //     }
    //   }

    // Relax
    relax = 0;
    if (DDSIP_param->relax > 2 && DDSIP_bb->curnode % DDSIP_param->relax)
        relax = 2;
    else if (DDSIP_param->relax == 2)
        relax = 2;
    else if (DDSIP_param->relax == 1)
        relax = 1;

    // Relax first stage
    if (relax == 1)
    {
        for (j = 0; j < DDSIP_bb->firstvar; j++)
        {
            if (DDSIP_bb->firsttype[j] == 'B' || DDSIP_bb->firsttype[j] == 'I')
                type[j] = 'C';
        }
        status = CPXchgctype (DDSIP_env, DDSIP_lp, DDSIP_bb->firstvar, DDSIP_bb->firstindex, type);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change types of first stage variables (CBLowerBound) \n");
            fprintf (DDSIP_outfile, "ERROR: Failed to change types of first stage variables (CBLowerBound) \n");
            if (DDSIP_param->outlev)
                fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to change types of first stage variables (CBLowerBound) \n");
            goto TERMINATE;
        }
    }
    // Relax integrality of all variables
    if (relax == 2)
        status = CPXchgprobtype (DDSIP_env, DDSIP_lp, CPXPROB_LP);

    // Change bounds according to node in bb tree
    if (DDSIP_bb->curbdcnt)
    {
        if (!(DDSIP_bb->dualitcnt))
            status = DDSIP_ChgBounds (1);
        else
            status = DDSIP_ChgBounds (0);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change problem \n");

            if (DDSIP_param->outlev > 11)
                fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to change problem \n");
            return status;
        }
    }
NEXT_TRY:
    // first take care of remains from former calls
    if (!DDSIP_param->cb_inherit || DDSIP_bb->dualitcnt || DDSIP_param->scalarization || DDSIP_bb->skip == 2)
    {
        int keepSolution, solKept = 0;
        double lhs;
        cutpool_t *currentCut;
        for (iscen = 0; iscen < DDSIP_param->scenarios; iscen++)
        {
            if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[iscen] != NULL)
            {
                if (DDSIP_bb->keepSols)
                {
                    // we are repeating the initial dual evaluation
                    // check whether the scenario solution violates a cut
                    keepSolution = 1;
                    if (DDSIP_bb->cutpool)
                    {
                        currentCut = DDSIP_bb->cutpool;
                        while (currentCut)
                        {
                            lhs = 0.;
                            for (int i = 0; i < DDSIP_bb->firstvar; i++)
                            {
                                lhs += (DDSIP_node[DDSIP_bb->curnode]->first_sol)[iscen][i] * currentCut->matval[i];
                            }
                            if (lhs < currentCut->rhs - 1.e-7)
                            {

                                if (DDSIP_param->outlev > 30)
                                    fprintf (DDSIP_bb->moreoutfile, "sol. of scen. %d violates cut %d, violation %g.\n", iscen + 1, currentCut->number, currentCut->rhs - lhs);
                                keepSolution = 0;
                                break;
                            }
                            currentCut = currentCut->prev;
                        }
                    }
                    if (keepSolution)
                    {
                        if (DDSIP_param->outlev > 20)
                        {
                            if (solKept)
                                fprintf (DDSIP_bb->moreoutfile, ", %d", iscen + 1);
                            else
                                fprintf (DDSIP_bb->moreoutfile, "   ## keeping sol. of scen. %d", iscen + 1);
                        }
                        solKept++;
                        continue;
                    }
                }
                if ((cnt = (int) ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[iscen])[DDSIP_bb->firstvar] - 0.9)))
                    for (k1 = iscen + 1; cnt && k1 < DDSIP_param->scenarios; k1++)
                    {
                        if ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[iscen]) == (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[k1]))
                        {
                            ((DDSIP_node[DDSIP_bb->curnode])->first_sol)[k1] = NULL;
                            cnt--;
                        }
                    }
                DDSIP_Free ((void **) &(((DDSIP_node[DDSIP_bb->curnode])->first_sol)[iscen]));
            }
        }
        if (solKept && DDSIP_param->outlev > 20)
            fprintf (DDSIP_bb->moreoutfile, ".\n");
    }
    // This node has been solved (as soon as we enter the loop below)
    DDSIP_node[DDSIP_bb->curnode]->solved = 1;
    DDSIP_bb->violations = DDSIP_bb->firstvar;
    tmpbestbound = 0.0;
    tmpupper     = 0.0;
    maxGap = meanGap = 0.0;

    // LowerBound problem for each scenario
    //****************************************************************************
    for (iscen = 0; iscen < DDSIP_param->scenarios; iscen++)
    {
        scen=DDSIP_bb->lb_scen_order[iscen];

        // User termination
        if (DDSIP_killsignal)
        {
            DDSIP_bb->skip = 2;
            goto TERMINATE;
        }
#ifdef SHIFT
        time (&startTime);
#endif
        nodes_1st = nodes_2nd = -1;
        if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen] == NULL)
        {
            // Output
            if (DDSIP_param->outlev || DDSIP_param->cpxscr)
            {
                if (DDSIP_param->cpxscr)
                    printf ("*****************************************************************************\n");
                printf ("Solving scenario problem %d (for dual optimization) in node %d...\n", scen + 1, DDSIP_bb->curnode);
            }
            nrFeasCheck = 0;
NEXT_SCEN:
            // Change problem according to scenario
            status = DDSIP_ChgProb (scen, DDSIP_bb->multipliers);
            if (status)
            {
                fprintf (stderr, "ERROR: Failed to change problem \n");
                fprintf (DDSIP_outfile, "ERROR: Failed to change problem \n");
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to change problem \n");
                goto TERMINATE;
            }
            // Warm starts
            // copy previous solution of the same scenario, if it exists
            if (DDSIP_param->hot)
            {
                status = DDSIP_Warm (iscen);
                if (status)
                    goto TERMINATE;
            }

            status = DDSIP_LowerBoundWriteLp (scen);
            // Error?
            if (status)
                goto TERMINATE;
            time_start = DDSIP_GetCpuTime ();
            // Optimize
            if (relax < 2)
            {
                // query time limit amd mip rel. gap parameters
#ifdef DEBUG
                if (DDSIP_param->cpxscr || DDSIP_param->outlev > 21)
                {
                    status = CPXgetdblparam (DDSIP_env,CPX_PARAM_TILIM,&we);
                    status = CPXgetdblparam (DDSIP_env,CPX_PARAM_EPGAP,&wr);
                    printf ("   -- 1st optimization time limit: %gs, rel. gap: %g%% --\n",we,wr*100.0);
                }
#endif
                // Optimize MIP
                DDSIP_bb->scenCBIters++;
                optstatus = CPXmipopt (DDSIP_env, DDSIP_lp);
                mipstatus = CPXgetstat (DDSIP_env, DDSIP_lp);
                nodes_1st = CPXgetnodecnt (DDSIP_env,DDSIP_lp);
                if (!optstatus && !DDSIP_Error(optstatus) && !DDSIP_Infeasible (mipstatus))
                {
                    if (CPXgetmiprelgap(DDSIP_env, DDSIP_lp, &mipgap))
                    {
                        fprintf (stderr, "ERROR: Query of CPLEX mip gap from 1st optimization failed (LowerBound) \n");
                        fprintf (stderr, "       CPXgetstat returned: %d\n",mipstatus);
                        mipgap = 1.e+30;
                    }
#ifdef DEBUG
                    time_lap = DDSIP_GetCpuTime ();
                    if (DDSIP_param->cpxscr || DDSIP_param->outlev > 11)
                    {
                        printf ("      CBLB: after 1st optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)\n",mipgap*100.0,nodes_1st,time_lap-time_start);
                        if (DDSIP_param->outlev)
                        {
                            if (use_LB_params)
                                fprintf (DDSIP_bb->moreoutfile,"      CBLB: after 1st optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)   - using CPLEX parameters for LB\n",mipgap*100.0,nodes_1st,time_lap-time_start);
                            else
                                fprintf (DDSIP_bb->moreoutfile,"      CBLB: after 1st optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)\n",mipgap*100.0,nodes_1st,time_lap-time_start);
                        }
                    }
#endif
                    if (DDSIP_param->watchkappa)
                    {
                        double maxkappaval, stablekappaval, suspiciouskappaval, unstablekappaval, illposedkappaval;
                        status = (CPXgetdblquality(DDSIP_env,DDSIP_lp,&maxkappaval,CPX_KAPPA_MAX) ||
                                  CPXgetdblquality(DDSIP_env,DDSIP_lp,&illposedkappaval,CPX_KAPPA_ILLPOSED) ||
                                  CPXgetdblquality(DDSIP_env,DDSIP_lp,&stablekappaval,CPX_KAPPA_STABLE) ||
                                  CPXgetdblquality(DDSIP_env,DDSIP_lp,&suspiciouskappaval,CPX_KAPPA_SUSPICIOUS) ||
                                  CPXgetdblquality(DDSIP_env,DDSIP_lp,&unstablekappaval,CPX_KAPPA_UNSTABLE));
                        if (status)
                        {
                            fprintf (stderr, "Failed to obtain Kappa information. \n\n");
                        }
                        else
                        {
                            printf("maximal Kappa value        = %g\n",maxkappaval);
                            if (stablekappaval < 1.)
                            {
                                if (stablekappaval)
                                    printf("numerically stable bases     %7.3f%%\n",stablekappaval*1e2);
                                if (suspiciouskappaval)
                                    printf("numerically suspicious bases %7.3f%%\n",suspiciouskappaval*1e2);
                                if (unstablekappaval)
                                    printf("numerically unstable bases   %7.3f%%\n",unstablekappaval*1e2);
                                if (illposedkappaval)
                                    printf("numerically ill-posed bases  %7.3f%%\n",illposedkappaval*1e2);
                            }
                            if (DDSIP_param->outlev)
                            {
                                fprintf(DDSIP_bb->moreoutfile, "             maximal Kappa value        = %g\n",maxkappaval);
                                if (stablekappaval < 1.)
                                {
                                    if (stablekappaval)
                                        fprintf(DDSIP_bb->moreoutfile, "             numerically stable bases     %7.3f%%\n",stablekappaval*1e2);
                                    if (suspiciouskappaval)
                                        fprintf(DDSIP_bb->moreoutfile, "             numerically suspicious bases %7.3f%%\n",suspiciouskappaval*1e2);
                                    if (unstablekappaval)
                                        fprintf(DDSIP_bb->moreoutfile, "             numerically unstable bases   %7.3f%%\n",unstablekappaval*1e2);
                                    if (illposedkappaval)
                                        fprintf(DDSIP_bb->moreoutfile, "             numerically ill-posed bases  %7.3f%%\n",illposedkappaval*1e2);
                                }
                            }
                        }
                    }

                    if (DDSIP_param->cpxnodual2 && mipstatus != CPXMIP_OPTIMAL)
                    {
                        // more iterations with different settings
                        if (use_LB_params && DDSIP_param->cpxnolb2)
                        {
                            status = DDSIP_SetCpxPara (DDSIP_param->cpxnolb2, DDSIP_param->cpxlbisdbl2, DDSIP_param->cpxlbwhich2, DDSIP_param->cpxlbwhat2);
                        }
                        else
                        {
                            status = DDSIP_SetCpxPara (DDSIP_param->cpxnodual2, DDSIP_param->cpxdualisdbl2, DDSIP_param->cpxdualwhich2, DDSIP_param->cpxdualwhat2);
                        }
                        if (status)
                        {
                            fprintf (stderr, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
                            fprintf (DDSIP_outfile, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
                            if (DDSIP_param->outlev)
                                fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
                            goto TERMINATE;
                        }
                        // query time limit amd mip rel. gap parameters
                        status = CPXgetdblparam (DDSIP_env,CPX_PARAM_EPGAP,&wr);
#ifdef DEBUG
                        if (DDSIP_param->cpxscr || DDSIP_param->outlev > 21)
                        {
                            status = CPXgetdblparam (DDSIP_env,CPX_PARAM_TILIM,&we);
                            printf ("   -- 2nd optimization time limit: %gs, rel. gap: %g%% --\n",we,wr*100.0);
                        }
#endif
                        // continue if desired gap is not reached yet
                        if (mipgap > wr)
                        {
                            if(DDSIP_param->files > 4)
                            {
                                char fname[DDSIP_ln_fname];
                                sprintf (fname, "%s/lb_sc%d_n%d%s", DDSIP_outdir, scen + 1, DDSIP_bb->curnode, "_2.prm");
                                status = CPXwriteparam(DDSIP_env,fname);
                                if (status)
                                {
                                    fprintf (stderr, "ERROR: Failed to write parameters\n");
                                    fprintf (DDSIP_outfile, "ERROR: Failed to write parameters\n");
                                    if (DDSIP_param->outlev)
                                        fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to write parameters\n");
                                    return status;
                                }
                                else if (DDSIP_param->outlev)
                                {
                                    fprintf (DDSIP_bb->moreoutfile, "      parameter file %s written.\n",fname);
                                }
                            }
                            optstatus = CPXmipopt (DDSIP_env, DDSIP_lp);
                            mipstatus = CPXgetstat (DDSIP_env, DDSIP_lp);
                            nodes_2nd = CPXgetnodecnt (DDSIP_env,DDSIP_lp);
#ifdef DEBUG
                            if (DDSIP_param->cpxscr || DDSIP_param->outlev > 11)
                            {
                                if (CPXgetmiprelgap(DDSIP_env, DDSIP_lp, &mipgap))
                                {
                                    fprintf (stderr, "ERROR: Query of CPLEX mip gap from 2nd optimization failed (LowerBound) \n");
                                    fprintf (stderr, "       CPXgetstat returned: %d\n",mipstatus);
                                    mipgap = 1.e+30;
                                }
                                else
                                {
                                    time_end = DDSIP_GetCpuTime ();
                                    printf ("      CBLB: after 2nd optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)\n",mipgap*100.0,nodes_2nd,time_end-time_lap);
                                    if (DDSIP_param->outlev)
                                    {
                                        if (use_LB_params)
                                            fprintf (DDSIP_bb->moreoutfile,"      CBLB: after 2nd optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)   - using CPLEX parameters for LB\n",mipgap*100.0,nodes_2nd,time_end-time_lap);
                                        else
                                            fprintf (DDSIP_bb->moreoutfile,"      CBLB: after 2nd optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)\n",mipgap*100.0,nodes_2nd,time_end-time_lap);
                                    }
                                }
                            }
#endif
                            if (DDSIP_param->watchkappa)
                            {
                                double maxkappaval, stablekappaval, suspiciouskappaval, unstablekappaval, illposedkappaval;
                                status = (CPXgetdblquality(DDSIP_env,DDSIP_lp,&maxkappaval,CPX_KAPPA_MAX) ||
                                          CPXgetdblquality(DDSIP_env,DDSIP_lp,&illposedkappaval,CPX_KAPPA_ILLPOSED) ||
                                          CPXgetdblquality(DDSIP_env,DDSIP_lp,&stablekappaval,CPX_KAPPA_STABLE) ||
                                          CPXgetdblquality(DDSIP_env,DDSIP_lp,&suspiciouskappaval,CPX_KAPPA_SUSPICIOUS) ||
                                          CPXgetdblquality(DDSIP_env,DDSIP_lp,&unstablekappaval,CPX_KAPPA_UNSTABLE));
                                if (status)
                                {
                                    fprintf (stderr, "Failed to obtain Kappa information. \n\n");
                                }
                                else
                                {
                                    printf("maximal Kappa value        = %g\n",maxkappaval);
                                    if (stablekappaval < 1.)
                                    {
                                        if (stablekappaval)
                                            printf("numerically stable bases     %7.3f%%\n",stablekappaval*1e2);
                                        if (suspiciouskappaval)
                                            printf("numerically suspicious bases %7.3f%%\n",suspiciouskappaval*1e2);
                                        if (unstablekappaval)
                                            printf("numerically unstable bases   %7.3f%%\n",unstablekappaval*1e2);
                                        if (illposedkappaval)
                                            printf("numerically ill-posed bases  %7.3f%%\n",illposedkappaval*1e2);
                                    }
                                    if (DDSIP_param->outlev)
                                    {
                                        fprintf(DDSIP_bb->moreoutfile, "             maximal Kappa value        = %g\n",maxkappaval);
                                        if (stablekappaval < 1.)
                                        {
                                            if (stablekappaval)
                                                fprintf(DDSIP_bb->moreoutfile, "             numerically stable bases     %7.3f%%\n",stablekappaval*1e2);
                                            if (suspiciouskappaval)
                                                fprintf(DDSIP_bb->moreoutfile, "             numerically suspicious bases %7.3f%%\n",suspiciouskappaval*1e2);
                                            if (unstablekappaval)
                                                fprintf(DDSIP_bb->moreoutfile, "             numerically unstable bases   %7.3f%%\n",unstablekappaval*1e2);
                                            if (illposedkappaval)
                                                fprintf(DDSIP_bb->moreoutfile, "             numerically ill-posed bases  %7.3f%%\n",illposedkappaval*1e2);
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            if (mipstatus == CPXMIP_TIME_LIM_FEAS || mipstatus == CPXMIP_NODE_LIM_FEAS)
                                mipstatus = CPXMIP_OPTIMAL_TOL;
                        }
                        // reset CPLEX parameters to lb
                        status = DDSIP_SetCpxPara (DDSIP_param->cpxnodual, DDSIP_param->cpxdualisdbl, DDSIP_param->cpxdualwhich, DDSIP_param->cpxdualwhat);
                        if (status)
                        {
                            fprintf (stderr, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
                            fprintf (DDSIP_outfile, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
                            if (DDSIP_param->outlev)
                                fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to set CPLEX parameters (LowerBound) \n");
                            goto TERMINATE;
                        }
                    }
                }

                if ((k = CPXgetnummipstarts(DDSIP_env, DDSIP_lp)) > 2)
                {
                    status    = CPXdelmipstarts (DDSIP_env, DDSIP_lp, 1, k-1);
                }
            }
            else
            {
                // Optimize relaxed MIP
                optstatus = CPXdualopt (DDSIP_env, DDSIP_lp);
                mipstatus = CPXgetstat (DDSIP_env, DDSIP_lp);
            }

            // We handle some errors separately (blatant infeasible, error in scenario problem)
            if (DDSIP_Error (optstatus))
            {
                fprintf (stderr, "WARNING: Failed to optimize problem for scenario %d.(CBLowerBound) in node %d, status = %d\n",
                         scen + 1, DDSIP_bb->curnode, optstatus);
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile, "WARNING: Failed to optimize problem for scenario %d in node %d.(CBLowerBound) status=%d\n",
                             scen + 1, DDSIP_bb->curnode, optstatus);
                if (optstatus == CPXERR_SUBPROB_SOLVE)
                {
                    optstatus = CPXgetsubstat(DDSIP_env, DDSIP_lp);
                    fprintf (stderr, "ERROR:                                scenario %d.(CBLowerBound) subproblem status=%d\n", scen + 1, optstatus);
                    fprintf (DDSIP_outfile, "ERROR:                                scenario %d.(CBLowerBound) subproblem status=%d\n", scen + 1, optstatus);
                    if (DDSIP_param->outlev)
                        fprintf (DDSIP_bb->moreoutfile, "ERROR:                                scenario %d.(CBLowerBound) subproblem status=%d\n", scen + 1, optstatus);
                }
                status = optstatus;
                goto TERMINATE;
            }
            // Infeasible, unbounded ... ?
            if (mipstatus==CPXMIP_TIME_LIM_INFEAS)
            {
                DDSIP_node[DDSIP_bb->curnode]->bound = DDSIP_infty;
                DDSIP_bb->skip = 1;
                DDSIP_bb->solstat[scen] = 0;
                status = mipstatus;
                CPXgetdblparam (DDSIP_env,CPX_PARAM_TILIM,&cpu_secs);
                fprintf (DDSIP_outfile, "ERROR: Problem infeasible within time limit of %g secs for scenario %d (CBLowerBound)\n", cpu_secs, scen + 1);
                if (DDSIP_param->outlev)
                {
                    printf ("ERROR: Problem infeasible within time limit of %g secs for scenario %d (CBLowerBound)\n", cpu_secs, scen + 1);
                    fprintf (DDSIP_bb->moreoutfile, "ERROR: Problem infeasible within time limit of %g secs for scenario %d (CBLowerBound)\n", cpu_secs, scen + 1);
                }
                goto TERMINATE;
            }
            else if (DDSIP_Infeasible (mipstatus))
            {
                DDSIP_bb->skip = 1;
                DDSIP_bb->solstat[scen] = 0;
                status = mipstatus;
                if (DDSIP_param->outlev)
                {
                    printf ("WARNING: Problem infeasible or unbounded for scenario %d in node %d (CBLowerBound)\n", scen + 1, DDSIP_bb->curnode);
                    fprintf (DDSIP_bb->moreoutfile, "WARNING: Problem infeasible or unbounded for scenario %d in node %d (CBLowerBound)\n", scen + 1, DDSIP_bb->curnode);
                    if (!DDSIP_bb->curnode)
                        fprintf (DDSIP_outfile, "ERROR: Problem infeasible or unbounded for scenario %d in node %d (CBLowerBound)\n", scen + 1, DDSIP_bb->curnode);
                }
                if ((DDSIP_bb->dualitcnt && DDSIP_bb->newTry < 6) || (!DDSIP_bb->dualitcnt && DDSIP_bb->newTry < 3))
                {
                    DDSIP_bb->newTry++;
                    status = -111;
                    if (DDSIP_param->outlev)
                        fprintf (DDSIP_bb->moreoutfile, "         newTry = %d, return value= %d\n", DDSIP_bb->newTry, status);
                }
                else
                {
                    DDSIP_bb->newTry = 0;
                    DDSIP_node[DDSIP_bb->curnode]->bound = DDSIP_infty;
                    fprintf (DDSIP_outfile, "         Problem infeasible or unbounded for scenario %d in node %d (CBLowerBound)\n", scen + 1, DDSIP_bb->curnode);
                }
                goto TERMINATE;
            }
            else
            {
                // We did not detect infeasibility
                // If we couldn't find a feasible solution we can at least obtain a lower bound
                if (DDSIP_NoSolution (mipstatus))
                {
                    objval = DDSIP_infty;
                    DDSIP_bb->solstat[scen] = 0;
                    if (relax != 2)
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
                    else
                        bobjval = DDSIP_infty;
                    if (status)
                    {
                        fprintf (stderr, "ERROR: Failed to get value of best remaining node (CBLowerBound)\n");
                        fprintf (DDSIP_outfile, "ERROR: Failed to get value of best remaining node (CBLowerBound)\n");
                        if (DDSIP_param->outlev)
                            fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to get value of best remaining node (CBLowerBound)\n");
                        goto TERMINATE;
                    }
                }
                else            // Solution exists
                {
                    DDSIP_bb->solstat[scen] = 1;
                    status = CPXgetobjval (DDSIP_env, DDSIP_lp, &objval);
                    if (status)
                    {
                        fprintf (stderr, "ERROR*: Failed to get best objective value \n");
                        fprintf (DDSIP_outfile, "ERROR*: Failed to get best objective value \n");
                        if (DDSIP_param->outlev)
                            fprintf (DDSIP_bb->moreoutfile, "ERROR*: Failed to get best objective value \n");
                        goto TERMINATE;
                    }
                    if (relax != 2)
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
                    else
                        bobjval = DDSIP_infty;
                    if (status)
                    {
                        fprintf (stderr, "ERROR: Failed to get value of best remaining node (CBLowerBound)\n");
                        fprintf (DDSIP_outfile, "ERROR: Failed to get value of best remaining node (CBLowerBound)\n");
                        if (DDSIP_param->outlev)
                            fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to get value of best remaining node (CBLowerBound)\n");
                        goto TERMINATE;
                    }

                    status = CPXgetx (DDSIP_env, DDSIP_lp, mipx, 0, DDSIP_bb->firstvar + DDSIP_bb->secvar - 1);
                    if (status)
                    {
                        fprintf (stderr, "ERROR: Failed to get solution \n");
                        fprintf (DDSIP_outfile, "ERROR: Failed to get solution \n");
                        if (DDSIP_param->outlev)
                            fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to get solution \n");
                        goto TERMINATE;
                    }
                    // Numerical errors?
                    for (j = 0; j < DDSIP_bb->firstvar; j++)
                    {
#ifdef CHECK_BOUNDS
                        // outside of bounds?
                        // loose - restrict the first-stage variables inside of original bounds (disregarding branching)
                        if (DDSIP_bb->lborg[j] - mipx[DDSIP_bb->firstindex[j]] > (DDSIP_bb->lborg[j] + 1.) * 1.e-9)
                        {
//#ifdef DEBUG
                            if (DDSIP_param->outlev > 21)
                                fprintf (DDSIP_bb->moreoutfile, "   ## First-stage var. %6d (%c):  lb= %.16g > %.16g =mipx, diff= %g\n", DDSIP_bb->firstindex[j], DDSIP_bb->firsttype[j], DDSIP_bb->lborg[j], mipx[DDSIP_bb->firstindex[j]], DDSIP_bb->lborg[j] - mipx[DDSIP_bb->firstindex[j]]);
//#endif
                            mipx[DDSIP_bb->firstindex[j]] = DDSIP_bb->lborg[j];
                        }
                        else
                        {
                            if (mipx[DDSIP_bb->firstindex[j]] - DDSIP_bb->uborg[j] > (DDSIP_bb->uborg[j] +1.) * 1.e-9)
                            {
//#ifdef DEBUG
                                if (DDSIP_param->outlev > 21)
                                    fprintf (DDSIP_bb->moreoutfile, "   ## First-stage var. %6d (%c):  ub= %.16g < %.16g =mipx, diff= %g\n", DDSIP_bb->firstindex[j], DDSIP_bb->firsttype[j], DDSIP_bb->uborg[j], mipx[DDSIP_bb->firstindex[j]], mipx[DDSIP_bb->firstindex[j]] - DDSIP_bb->uborg[j]);
//#endif
                                mipx[DDSIP_bb->firstindex[j]] = DDSIP_bb->uborg[j];
                            }
                        }
#endif
                        // tiny value, probably zero?
                        if (fabs (mipx[DDSIP_bb->firstindex[j]]) < 1e-18)
                        {
//#ifdef DEBUG
                            if (DDSIP_param->outlev > 21 && mipx[DDSIP_bb->firstindex[j]])
                                fprintf (DDSIP_bb->moreoutfile, "   ## First-stage var. %6d (%c):  tiny mipx %.6g -> 0.\n", DDSIP_bb->firstindex[j], DDSIP_bb->firsttype[j], mipx[DDSIP_bb->firstindex[j]]);
//#endif
                            mipx[DDSIP_bb->firstindex[j]] = 0.0;
                        }
                    }
                    if (relax != 2)
                    {
                        //if (mipstatus == CPXMIP_OPTIMAL)
                        //{
                        //    bobjval = objval;
                        //}
                        //else
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
                    }
                    else
                        bobjval = DDSIP_infty;
                    if (status)
                    {
                        fprintf (stderr, "ERROR: Failed to get value of best remaining node\n");
                        fprintf (DDSIP_outfile, "ERROR: Failed to get value of best remaining node\n");
                        if (DDSIP_param->outlev)
                            fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to get value of best remaining node\n");
                        goto TERMINATE;
                    }
                }
            }                // end for else

            (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen] = objval;
            // If all scenario problems were solved to optimality ....
            if (!(mipstatus == CPXMIP_OPTIMAL))
            {
                // If the tree was exhausted bobjval is huge
                if (fabs (bobjval) > DDSIP_infty)
                    bobjval = objval;
            }
            (DDSIP_node[DDSIP_bb->curnode]->mipstatus)[scen] = mipstatus;
            gap = 100.0*(objval-bobjval)/(fabs(objval)+1e-4);
            if (fabs(gap) < 2.e-13)
                gap = 0.;
            meanGap += DDSIP_data->prob[scen] * gap;
            maxGap   = DDSIP_Dmax (maxGap, gap);
            time_end = DDSIP_GetCpuTime ();
            time_help = time_end-time_start;
            time (&DDSIP_bb->cur_time);
#ifdef SHIFT
            if (shift_in_cb)
            {
                // prepare for shifting difficult scenarios to the end such that they can make use of more mipstarts
                DDSIP_bb->aggregate_time[scen] += difftime (DDSIP_bb->cur_time, startTime) + 0.1*time_help + 1.e3*gap;
            }
#endif
            // Debugging information
            if (DDSIP_param->outlev)
            {
                DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                DDSIP_translate_time (time_end,&cpu_hrs,&cpu_mins,&cpu_secs);
                if (mipstatus == CPXMIP_OPTIMAL)
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%)     ",
                             iscen + 1, scen + 1, objval, bobjval, gap);
                else if (mipstatus == CPXMIP_OPTIMAL_TOL)
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) tol.",
                             iscen + 1, scen + 1, objval, bobjval, gap);
                else if (mipstatus == CPXMIP_NODE_LIM_FEAS)
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) NODE",
                             iscen + 1, scen + 1, objval, bobjval, gap);
                else if (mipstatus == CPXMIP_TIME_LIM_FEAS)
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) TIME",
                             iscen + 1, scen + 1, objval, bobjval, gap);
                else
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) %-4.0d ",
                             iscen + 1, scen + 1, objval, bobjval, gap, mipstatus);
                fprintf (DDSIP_bb->moreoutfile,
                         "\t %3dh %02d:%02.0f / %3dh %02d:%05.2f (%7.2fs n:%5d",
                         wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs, time_help, nodes_1st);
                if (nodes_2nd >= 0)
                    fprintf (DDSIP_bb->moreoutfile, " +%4d)", nodes_2nd - nodes_1st);
                else
                    fprintf (DDSIP_bb->moreoutfile, ")");
            }
            // Remember values of integer varables for warm starts
            if (DDSIP_bb->solstat[scen] && DDSIP_param->hot)
            {
                for (j = 0; j < DDSIP_bb->total_int; j++)
                    DDSIP_node[DDSIP_bb->curnode]->solut[scen * DDSIP_bb->total_int + j] = mipx[DDSIP_bb->intind[j]];
            }
            wr = DDSIP_Dmin (bobjval, objval);
            (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen] =  wr;

            // Wait-and-see lower bounds, should be set only once
            // Test if btlb contain still the initial value -DDSIP_infty
            if (!(DDSIP_bb->btlb[scen] > -DDSIP_infty))
                DDSIP_bb->btlb[scen] = wr;
            // If a feasible solution exists
            if (DDSIP_bb->solstat[scen])
            {
                if (DDSIP_param->hot)
                {
                    for (j = 0; j < DDSIP_bb->total_int; j++)
                    {
                        DDSIP_bb->intsolvals[scen][j] = floor (mipx[DDSIP_bb->intind[j]] + 0.1);
                    }
                }
                for (i_scen = 0; i_scen < DDSIP_param->scenarios; i_scen++)
                {
                    if (i_scen == scen)
                        continue;
                    if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[i_scen])
                    {
                        for (j = 0; j < DDSIP_bb->firstvar; j++)
                        {
                            if (!DDSIP_Equal (mipx[DDSIP_bb->firstindex[j]], (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i_scen][j]))
                            {
                                break;
                            }
                        }
                        if (j >= DDSIP_bb->firstvar)
                        {
                            break;
                        }
                    }
                }
                if (i_scen < DDSIP_param->scenarios)
                {
                    (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen] = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i_scen];
                    //first_sol[DDSIP_bb->firstvar]  equals the number of identical first stage scenario solutions
                    (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][DDSIP_bb->firstvar] += 1.0;
                    if (DDSIP_param->outlev)
                    {
                        if (nodes_2nd < 0)
                            fprintf (DDSIP_bb->moreoutfile," \t\t->scen %4d (%3g ident.)\n", i_scen + 1, DDSIP_node[DDSIP_bb->curnode]->first_sol[scen][DDSIP_bb->firstvar]);
                        else
                            fprintf (DDSIP_bb->moreoutfile," \t->scen %4d (%3g ident.)\n", i_scen + 1, DDSIP_node[DDSIP_bb->curnode]->first_sol[scen][DDSIP_bb->firstvar]);
                    }
                }
                else
                {
                    if (DDSIP_param->outlev)
                    {
                        fprintf (DDSIP_bb->moreoutfile,"\n");
                    }
                    //DDSIP_Allocate space for the first stage variables
                    if (!(DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen])
                        (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen] =
                            (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar + 3, "DDSIP_node[curode]->first_sol[scen](CBLowerBound)");
                    //first_sol[DDSIP_bb->firstvar]  equals the number of identical first stage scenario solutions
                    (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][DDSIP_bb->firstvar] = 1.0;
                    (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][DDSIP_bb->firstvar + 1] = 0.0;
                    (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][DDSIP_bb->firstvar + 2] = DDSIP_bb->curnode;
                    if (DDSIP_param->outlev >=  DDSIP_first_stage_outlev)
                        fprintf (DDSIP_bb->moreoutfile, "    First-stage solution:\n");
                    for (j = 0; j < DDSIP_bb->firstvar; j++)
                    {
                        // Collect solution vector
                        if (DDSIP_bb->firsttype[j] == 'B' || DDSIP_bb->firsttype[j] == 'I' || DDSIP_bb->firsttype[j] == 'N')
                            (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j] = floor (mipx[DDSIP_bb->firstindex[j]] + 0.1);
                        else
                            (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j] = mipx[DDSIP_bb->firstindex[j]];
                        // Print something
                        if (DDSIP_param->outlev >= DDSIP_first_stage_outlev)
                        {
                            fprintf (DDSIP_bb->moreoutfile," %21.14g", (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j]);
                            if (!((j + 1) % 5))
                                fprintf (DDSIP_bb->moreoutfile, "\n");
                        }
                        // Calculate minimum and maximum of each component
                        minfirst[j] = DDSIP_Dmin (mipx[DDSIP_bb->firstindex[j]], minfirst[j]);
                        maxfirst[j] = DDSIP_Dmax (mipx[DDSIP_bb->firstindex[j]], maxfirst[j]);
                    }
                    if (DDSIP_param->outlev >= DDSIP_first_stage_outlev)
                        fprintf (DDSIP_bb->moreoutfile, "\n");
                    //
                    if ((DDSIP_param->addBendersCuts || DDSIP_param->addIntegerCuts) && (DDSIP_param->numberScenReeval > -1) && ((!(DDSIP_bb->curnode) && (DDSIP_bb->dualdescitcnt < 3)) || ((DDSIP_bb->curnode < 11) && !DDSIP_bb->dualdescitcnt)) && (j = DDSIP_bb->cutCntr))
                    {
                        // check feasibility for other scenarios and possibly add cuts
                        tmp = DDSIP_bb->sug[DDSIP_param->nodelim + 2];
                        if (!tmp)
                        {
                            // DDSIP_Allocate memory for heuristics suggestion if there is none yet
                            tmp = (sug_t *)DDSIP_Alloc (sizeof (sug_t), 1, "sug(LowerBound)");
                            tmp->firstval = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar, "sug->firstval(LowerBound)");
                            tmp->next = NULL;
                            DDSIP_bb->sug[DDSIP_param->nodelim + 2] = tmp;
                        }
                        else
                        {
                            if (tmp->next != NULL)
                            {
                                printf("\n\n FEHLER: !!!!!!!!!!!!!!!!!!!!! DDSIP_bb->sug[DDSIP_param->nodelim + 2]->next = %p\n\n",DDSIP_bb->sug[DDSIP_param->nodelim + 2]->next);
                                exit (1);
                            }
                        }

                        memcpy (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval, (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen], sizeof (double) * DDSIP_bb->firstvar);

                        if ((status = DDSIP_UpperBound (0, 2)))
                        {
                            fprintf (stderr, "ERROR: Failed to perform feasibility check\n");
                            return status;
                        }
                        nrFeasCheck++;
                        // Change bounds according to node in bb tree
                        if (DDSIP_bb->curbdcnt)
                        {
                            status = DDSIP_ChgBounds (0);
                            if (status)
                            {
                                fprintf (stderr, "ERROR: Failed to change problem \n");

                                if (DDSIP_param->outlev > 11)
                                    fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to change problem \n");
                                return status;
                            }
                        }
                        if (j < DDSIP_bb->cutCntr && nrFeasCheck <= DDSIP_param->numberScenReeval)
                        {
                            if(DDSIP_param->outlev)
                                fprintf (DDSIP_bb->moreoutfile, " --> scenario reevaluation %d\n", nrFeasCheck);
                            goto NEXT_SCEN;
                        }
                    }
                    //
                    if (DDSIP_param->outlev >= DDSIP_second_stage_outlev)
                    {
                        fprintf (DDSIP_bb->moreoutfile, "    Second-stage solution:\n");
                        for (j = 0; j < DDSIP_bb->secvar; j++)
                        {
                            fprintf (DDSIP_bb->moreoutfile, " %21.14g",  mipx[DDSIP_bb->secondindex[j]]);
                            if (!((j + 1) % 5))
                                fprintf (DDSIP_bb->moreoutfile, "\n");
                        }
                        fprintf (DDSIP_bb->moreoutfile, "\n");
                    }
                }
                // Remember objective function contribution
                DDSIP_Contrib_LB (mipx, scen);
                if (DDSIP_param->scalarization && abs(DDSIP_param->riskmod) < 2)
                {
                    if ((DDSIP_node[DDSIP_bb->curnode]->ref_scenobj)[scen] > DDSIP_param->risktarget)
                    {
                        if (DDSIP_param->riskmod == 1)
                            DDSIP_bb->ref_risk[scen] = (DDSIP_node[DDSIP_bb->curnode]->ref_scenobj)[scen] - DDSIP_param->risktarget;
                        else if (DDSIP_param->riskmod == 2)
                            DDSIP_bb->ref_risk[scen] = 1.;
                    }
                    else
                        DDSIP_bb->ref_risk[scen] = 0.;
                }
            }                // end if solstat
            wr = DDSIP_Dmin (bobjval, objval);
            // Temporary bound is infinity if a scenario problem is infeasible
            if (DDSIP_Equal (wr, DDSIP_infty))
            {
                tmpbestbound = DDSIP_infty;
                tmpupper     = DDSIP_infty;
            }
            else
            {
                tmpbestbound += DDSIP_data->prob[scen] * wr;
                tmpupper     += DDSIP_data->prob[scen] * objval;
            }
        }
        else
        {
            tmpbestbound += DDSIP_data->prob[scen] *  (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen];
            tmpupper     += DDSIP_data->prob[scen] *  (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen];
            gap = 100.0*((DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen]-(DDSIP_node[DDSIP_bb->curnode]->subbound)[scen])/
                  (fabs((DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen])+1e-4);
            if (fabs(gap) < 2.e-13)
                gap = 0.;
            meanGap += DDSIP_data->prob[scen] * gap;
            maxGap   = DDSIP_Dmax (maxGap, gap);
            time (&DDSIP_bb->cur_time);
            // Debugging information
            if (DDSIP_param->outlev)
            {
                time_end = DDSIP_GetCpuTime ();
                DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                DDSIP_translate_time (time_end,&cpu_hrs,&cpu_mins,&cpu_secs);
                if ((DDSIP_node[DDSIP_bb->curnode]->mipstatus)[scen] == CPXMIP_OPTIMAL)
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%)     ",
                             iscen + 1,scen + 1, (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen], (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen], gap);
                else if ((DDSIP_node[DDSIP_bb->curnode]->mipstatus)[scen] == CPXMIP_OPTIMAL_TOL)
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) tol.",
                             iscen + 1,scen + 1, (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen], (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen], gap);
                else if ((DDSIP_node[DDSIP_bb->curnode]->mipstatus)[scen] == CPXMIP_NODE_LIM_FEAS)
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) NODE",
                             iscen + 1,scen + 1, (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen], (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen], gap);
                else if ((DDSIP_node[DDSIP_bb->curnode]->mipstatus)[scen] == CPXMIP_TIME_LIM_FEAS)
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) TIME",
                             iscen + 1,scen + 1, (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen], (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen], gap);
                else
                    fprintf (DDSIP_bb->moreoutfile,
                             "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g (%10.4g%%) %-4.0d ",
                             iscen + 1,scen + 1, (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[scen], (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen], gap,
                             (DDSIP_node[DDSIP_bb->curnode]->mipstatus)[scen]);
                fprintf (DDSIP_bb->moreoutfile,
                         "\t %3dh %02d:%02.0f / %3dh %02d:%05.2f  > node %3g",
                         wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs, (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][DDSIP_bb->firstvar + 2]);
                for (i_scen = 0; i_scen < scen; i_scen++)
                {
                    // if the solutions have the same address, they were equal
                    if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[i_scen] == (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen])
                        break;
                }
                if (i_scen < scen)
                {
                    //first_sol[DDSIP_bb->firstvar]   equals the number of identical first stage scenario solutions
                    if (nodes_2nd < 0)
                        fprintf (DDSIP_bb->moreoutfile," \t\t->scen %4d (%3g ident.)\n", i_scen + 1, DDSIP_node[DDSIP_bb->curnode]->first_sol[scen][DDSIP_bb->firstvar]);
                    else
                        fprintf (DDSIP_bb->moreoutfile," \t->scen %4d (%3g ident.)\n", i_scen + 1, DDSIP_node[DDSIP_bb->curnode]->first_sol[scen][DDSIP_bb->firstvar]);
                }
                else
                    fprintf (DDSIP_bb->moreoutfile,"\n");
            }
            for (j = 0; j < DDSIP_bb->firstvar; j++)
            {
                // Calculate minimum and maximum of each component
                minfirst[j] = DDSIP_Dmin ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j], minfirst[j]);
                maxfirst[j] = DDSIP_Dmax ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][j], maxfirst[j]);
            }
            wr = (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen];
        }

    }                // end for iscen
#ifdef SHIFT
    if (shift_in_cb)
    {
        DDSIP_qsort_ins_A (DDSIP_bb->aggregate_time, DDSIP_bb->lb_scen_order, 0, DDSIP_param->scenarios-1);
        DDSIP_bb->lb_sorted = 0;
        if (DDSIP_param->outlev > 21)
        {
            // debug output
            fprintf (DDSIP_bb->moreoutfile,"order of scenarios after sorting for lb according to time and gap\n");
            for (iscen = 0; iscen < DDSIP_param->scenarios; iscen++)
            {
                fprintf(DDSIP_bb->moreoutfile," %3d: Scen %3d  %20.14g\n",
                        iscen+1, DDSIP_bb->lb_scen_order[iscen]+1, DDSIP_bb->aggregate_time[DDSIP_bb->lb_scen_order[iscen]]);
            }
        }
    }
#endif

    // Count number of differences within first stage solution in current node
    // DDSIP_bb->violations=k means differences in k components
    DDSIP_bb->violations = 0;
    maxdispersion = 0.;
    for (j = 0; j < DDSIP_bb->firstvar; j++)
        if (!DDSIP_Equal (maxfirst[j], minfirst[j]))
        {
            DDSIP_bb->violations++;
            maxdispersion = DDSIP_Dmax (maxdispersion, maxfirst[j]-minfirst[j]);
        }
    DDSIP_node[DDSIP_bb->curnode]->violations = DDSIP_bb->violations;
    DDSIP_node[DDSIP_bb->curnode]->dispnorm = maxdispersion;
    *objective_val = -tmpbestbound;
    DDSIP_bb->currentDualObjVal = tmpbestbound;

    if ((DDSIP_bb->currentDualObjVal < DDSIP_bb->bestvalue) && !DDSIP_bb->violations)
    {
        int comb = 22;
        j = DDSIP_param->heuristic;
        DDSIP_param->heuristic = 100;
        if (DDSIP_param->outlev)
        {
            fprintf (DDSIP_bb->moreoutfile, "\tCurrent dual objective value             = \t%-18.16g          \t(mean MIP gap: %g%%, max gap: %g%%)\n", -(*objective_val), meanGap, maxGap);
            fprintf (DDSIP_bb->moreoutfile, "\tCurrent best objective value             = \t%-18.16g\n", tmpupper);
            if ((tmpupper < DDSIP_bb->bestvalue) && (!DDSIP_param->riskmod || !DDSIP_param->riskalg))
            {
                fprintf (DDSIP_bb->moreoutfile, "\t(Current best bound, improvement %16.12g, %g%%)\n",DDSIP_bb->bestvalue-tmpupper,1e+2*(DDSIP_bb->bestvalue-tmpupper)/(fabs(DDSIP_bb->bestvalue) + 1.e-16));
            }
            if (!use_LB_params && maxGap > 0.5*DDSIP_param->relgap)
            {
                if (DDSIP_param->outlev)
                {
                    printf ("\nNo violations, but max MIP gap too big - try again with CPLEX parameters for LB\n");
                    fprintf (DDSIP_bb->moreoutfile, "\n### No violations, but max MIP gap too big - try again with CPLEX parameters for LB\n");
                }
                use_LB_params = 1;
                DDSIP_bb->violations = 999;
                goto NEXT_TRY;
            }
            printf ("\nUpper bounds for solution with no violations\n");
            fprintf (DDSIP_bb->moreoutfile, "\nUpper bounds for solution with no violations\n");
        }
        if (!DDSIP_Heuristics (&comb, DDSIP_param->scenarios, 0))
        {
            int jj = -1;
            // Evaluate the proposed first-stage solution
            DDSIP_EvaluateScenarioSolutions (&jj);
        }
        DDSIP_param->heuristic = j;
        // in case the cplex parameters for cb produce too big gaps, try again with LB parameters
        if (!use_LB_params && maxGap > 0.5*DDSIP_param->relgap)
        {
            if (DDSIP_param->outlev)
            {
                printf ("\nNo violations, but max MIP gap too big - try again with CPLEX parameters for LB\n");
                fprintf (DDSIP_bb->moreoutfile, "\n### No violations, but max MIP gap too big - try again with CPLEX parameters for LB\n");
            }
            use_LB_params = 1;
            DDSIP_bb->violations = 999;
            goto NEXT_TRY;
        }
        else
        {
            DDSIP_node[DDSIP_bb->curnode]->leaf = 1;
        }
    }

    DDSIP_bb->ref_max = -1;
    if (DDSIP_param->riskalg || DDSIP_param->scalarization)
    {
        double tmp1 = DDSIP_RiskLb (DDSIP_node[DDSIP_bb->curnode]->ref_scenobj);
        if (DDSIP_param->outlev > 20)
            fprintf (DDSIP_bb->moreoutfile, " cb:  tmpbestbound = %g, DDSIP_RiskLb yields %g\n", tmpbestbound, tmp1);
        if (tmpbestbound < tmp1)
        {
            tmpbestbound = tmp1;
            if (DDSIP_param->outlev > 20)
                fprintf (DDSIP_bb->moreoutfile, " cb:  new tmpbestbound = %g\n", tmpbestbound);
        }
    }
    *objective_val = -tmpbestbound;
    DDSIP_bb->currentDualObjVal = tmpbestbound;
    DDSIP_bb->CBIters++;
    if (tmpbestbound > DDSIP_node[DDSIP_bb->curnode]->bound /*- fabs(DDSIP_node[DDSIP_bb->curnode]->bound)*2.e-15*/)
    {
        DDSIP_bb->dualObjVal = tmpbestbound;
        increase = 1;
        no_viol_worse = use_LB_params = 0;
        if (DDSIP_param->outlev > 6)
            fprintf (DDSIP_bb->moreoutfile,
                     " +++++ dual step     increasing  bound for node %3d, new val: %-18.16g, old value: %-18.16g  incr.  %g, rel %g%%,  weight = %g +++++\n",
                     DDSIP_bb->curnode, tmpbestbound,  DDSIP_node[DDSIP_bb->curnode]->bound, tmpbestbound - DDSIP_node[DDSIP_bb->curnode]->bound,
                     1.e2*(tmpbestbound - DDSIP_node[DDSIP_bb->curnode]->bound)/(fabs(DDSIP_node[DDSIP_bb->curnode]->bound)+1e-16), cb_get_last_weight(DDSIP_bb->dualProblem));
        DDSIP_node[DDSIP_bb->curnode]->violations = DDSIP_bb->violations;
        if (!DDSIP_bb->curnode && !DDSIP_bb->dualdescitcnt && !DDSIP_bb->initial_multiplier)
        {
            memcpy (DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag, DDSIP_node[DDSIP_bb->curnode]->subbound, DDSIP_param->scenarios*sizeof(double));
            DDSIP_node[DDSIP_bb->curnode]->BoundNoLag = DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag[0] * DDSIP_data->prob[0];
            for (j = 1; j < DDSIP_param->scenarios; j++)
                DDSIP_node[DDSIP_bb->curnode]->BoundNoLag += DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag[j] * DDSIP_data->prob[j];
        }
        if (DDSIP_bb->dualdescitcnt)
        {
            // Store the current best multipliers in bestdual
            memcpy (DDSIP_bb->local_bestdual, DDSIP_node[DDSIP_bb->curnode]->dual, sizeof (double) * (DDSIP_bb->dimdual));
            // the current weight and the node number, too
            DDSIP_bb->local_bestdual[DDSIP_bb->dimdual] = cb_get_last_weight(DDSIP_bb->dualProblem);
            DDSIP_bb->local_bestdual[DDSIP_bb->dimdual + 1] = DDSIP_bb->curnode;
            DDSIP_bb->local_bestdual[DDSIP_bb->dimdual + 2] = DDSIP_bb->dualitcnt + 1;
            // for the case of no bound increase, when the initial evaluation (with weight -1) was best, make the weight 0.1
            if (DDSIP_bb->local_bestdual[DDSIP_bb->dimdual] < 0.)
                DDSIP_bb->local_bestdual[DDSIP_bb->dimdual] = 0.1;
            if ((tmpbestbound - DDSIP_node[DDSIP_bb->curnode]->bound)/(fabs(DDSIP_node[DDSIP_bb->curnode]->bound)+1e-16) > 5.e-8)
            {
                if (DDSIP_param->outlev && nearly_constant)
                {
                    fprintf (DDSIP_bb->moreoutfile, " -**** bound: %15.8g, tmpbestbound: %15.8g, rel. change: %g ->nearly_constant=0\n", DDSIP_node[DDSIP_bb->curnode]->bound, tmpbestbound, (DDSIP_node[DDSIP_bb->curnode]->bound - tmpbestbound)/(fabs(DDSIP_node[DDSIP_bb->curnode]->bound)+1e-16));
                }
                nearly_constant = totally_constant = 0;
            }
            if (DDSIP_param->outlev > 10)
                fprintf (DDSIP_bb->moreoutfile," +++-> bestdual for node updated in desc. it. %d, tot. it. %g  (weight %g)\n", DDSIP_bb->dualdescitcnt, DDSIP_bb->local_bestdual[DDSIP_bb->dimdual + 2],  DDSIP_bb->local_bestdual[DDSIP_bb->dimdual]);
        }
        if (DDSIP_param->hot)
        {
            for (j = 0; j < DDSIP_param->scenarios; j++)
                memcpy (DDSIP_bb->boundIncrease_intsolvals[j], DDSIP_bb->intsolvals[j], DDSIP_bb->total_int*sizeof(double));
        }
        // Set lower bound for the node
        DDSIP_node[DDSIP_bb->curnode]->bound = tmpbestbound;
        // Free previous first stage in best point
        for (j = 0; j < DDSIP_param->scenarios; j++)
        {
            if (DDSIP_bb->bestfirst[j].first_sol)
            {
                if ((cnt = (int) (DDSIP_bb->bestfirst[j].first_sol[DDSIP_bb->firstvar] - 0.9)))
                {
                    for (k1 = j + 1; cnt && k1 < DDSIP_param->scenarios; k1++)
                    {
                        if (DDSIP_bb->bestfirst[j].first_sol == DDSIP_bb->bestfirst[k1].first_sol)
                        {
                            DDSIP_bb->bestfirst[k1].first_sol = NULL;
                            cnt--;
                        }
                    }
                }
                DDSIP_Free ((void **) &(DDSIP_bb->bestfirst[j].first_sol));
            }
        }
        // Save first stage in best point
        for (j = 0; j < DDSIP_param->scenarios; j++)
        {
            if (!DDSIP_node[DDSIP_bb->curnode]->first_sol[j])
            {
                fprintf (DDSIP_outfile, "ERROR: DDSIP_node[%d]->first_sol[%d] = %p\n", DDSIP_bb->curnode, j,  DDSIP_node[DDSIP_bb->curnode]->first_sol[j]);
                exit(99);
            }
            if (!DDSIP_bb->bestfirst[j].first_sol)
            {
                DDSIP_bb->bestfirst[j].first_sol = (double *) DDSIP_Alloc(sizeof(double),DDSIP_bb->firstvar+3,"bestfirst.first_sol(CBLowerBound)");
                memcpy (DDSIP_bb->bestfirst[j].first_sol, DDSIP_node[DDSIP_bb->curnode]->first_sol[j], sizeof (double) * (DDSIP_bb->firstvar + 3));
                if (DDSIP_param->outlev > 69)
                {
                    fprintf (DDSIP_bb->moreoutfile,
                             " **** saved solution for scenario %d in bestfirst, %g identical sols.\n", j, (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j])[DDSIP_bb->firstvar]);
                    for (k1 = 0; k1 < DDSIP_bb->firstvar; k1++)
                    {
                        fprintf (DDSIP_bb->moreoutfile," %20.14g,", DDSIP_bb->bestfirst[j].first_sol[k1]);
                        if (!((k1 + 1) % 5))
                            fprintf (DDSIP_bb->moreoutfile, "\n");
                    }
                    fprintf (DDSIP_bb->moreoutfile, "\n");
                }
                if ((cnt = (int) ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j])[DDSIP_bb->firstvar] - 0.9)))
                    for (k1 = j + 1; cnt && k1 < DDSIP_param->scenarios; k1++)
                    {
                        if ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j]) == (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[k1]))
                        {
                            DDSIP_bb->bestfirst[k1].first_sol = DDSIP_bb->bestfirst[j].first_sol;
                            cnt--;
                            if (DDSIP_param->outlev > 69)
                                fprintf (DDSIP_bb->moreoutfile,
                                         "       same solution for scenarios %d and %d, %d identical remaining.\n", j, k1,cnt);
                        }
                    }
            }
            else
            {
                if (DDSIP_param->outlev > 79)
                    fprintf (DDSIP_bb->moreoutfile,
                             " ++++ bestfirst->first_sol for scenario %d already present.\n", j);
            }
            DDSIP_bb->bestfirst[j].cursubsol = (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j];
            DDSIP_bb->bestfirst[j].subbound  = (DDSIP_node[DDSIP_bb->curnode]->subbound)[j];
        }
    }
    else
    {
        increase = 0;
        if (DDSIP_param->outlev > 6)
            fprintf (DDSIP_bb->moreoutfile,
                     " -**** dual step not increasing  bound for node %3d, new val: %-18.16g, old bound: %-18.16g (diff: %.6g, rel %g%%) weight = %g ****-\n",
                     DDSIP_bb->curnode, tmpbestbound, DDSIP_node[DDSIP_bb->curnode]->bound, tmpbestbound -  DDSIP_node[DDSIP_bb->curnode]->bound, 1e2*(tmpbestbound -  DDSIP_node[DDSIP_bb->curnode]->bound)/(fabs(DDSIP_node[DDSIP_bb->curnode]->bound)+1e-16), cb_get_last_weight(DDSIP_bb->dualProblem));
        if ((DDSIP_node[DDSIP_bb->curnode]->bound - tmpbestbound)/(fabs(DDSIP_node[DDSIP_bb->curnode]->bound)+1e-16) > 1.e-7)
        {
            if (DDSIP_param->outlev && nearly_constant)
            {
                fprintf (DDSIP_bb->moreoutfile, " -**** bound: %15.8g, tmpbestbound: %15.8g, rel. change: %g ->nearly_constant=0\n", DDSIP_node[DDSIP_bb->curnode]->bound, tmpbestbound, (DDSIP_node[DDSIP_bb->curnode]->bound - tmpbestbound)/(fabs(DDSIP_node[DDSIP_bb->curnode]->bound)+1e-16));
            }
            nearly_constant = totally_constant = 0;
        }
        // even when bound was not increased, store the first-stage solution in initial evaluation
        if (cb_get_last_weight(DDSIP_bb->dualProblem) < 0.)
        {
            // Free previous first stage in best point
            for (j = 0; j < DDSIP_param->scenarios; j++)
            {
                if (DDSIP_bb->bestfirst[j].first_sol)
                {
                    if ((cnt = (int) (DDSIP_bb->bestfirst[j].first_sol[DDSIP_bb->firstvar] - 0.9)))
                    {
                        for (k1 = j + 1; cnt && k1 < DDSIP_param->scenarios; k1++)
                        {
                            if (DDSIP_bb->bestfirst[j].first_sol == DDSIP_bb->bestfirst[k1].first_sol)
                            {
                                DDSIP_bb->bestfirst[k1].first_sol = NULL;
                                cnt--;
                            }
                        }
                    }
                    DDSIP_Free ((void **) &(DDSIP_bb->bestfirst[j].first_sol));
                }
            }
            // Save first stage in best point
            for (j = 0; j < DDSIP_param->scenarios; j++)
            {
                if (!DDSIP_bb->bestfirst[j].first_sol)
                {
                    DDSIP_bb->bestfirst[j].first_sol = (double *) DDSIP_Alloc(sizeof(double),DDSIP_bb->firstvar+3,"bestfirst.first_sol(CBLowerBound)");
                    memcpy (DDSIP_bb->bestfirst[j].first_sol, DDSIP_node[DDSIP_bb->curnode]->first_sol[j], sizeof (double) * (DDSIP_bb->firstvar + 3));
                    if (DDSIP_param->outlev > 69)
                    {
                        fprintf (DDSIP_bb->moreoutfile,
                                 " **** saved solution for scenario %d in bestfirst, %g identical sols.\n", j, (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j])[DDSIP_bb->firstvar]);
                        for (k1 = 0; k1 < DDSIP_bb->firstvar; k1++)
                        {
                            fprintf (DDSIP_bb->moreoutfile," %20.14g,", DDSIP_bb->bestfirst[j].first_sol[k1]);
                            if (!((k1 + 1) % 5))
                                fprintf (DDSIP_bb->moreoutfile, "\n");
                        }
                        fprintf (DDSIP_bb->moreoutfile, "\n");
                    }
                    if ((cnt = (int) ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j])[DDSIP_bb->firstvar] - 0.9)))
                        for (k1 = j + 1; cnt && k1 < DDSIP_param->scenarios; k1++)
                        {
                            if ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j]) == (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[k1]))
                            {
                                DDSIP_bb->bestfirst[k1].first_sol = DDSIP_bb->bestfirst[j].first_sol;
                                cnt--;
                                if (DDSIP_param->outlev > 69)
                                    fprintf (DDSIP_bb->moreoutfile,
                                             "       same solution for scenarios %d and %d, %d identical remaining.\n", j, k1,cnt);
                            }
                        }
                }
                else
                {
                    if (DDSIP_param->outlev > 79)
                        fprintf (DDSIP_bb->moreoutfile,
                                 " ++++ bestfirst->first_sol for scenario %d already present.\n", j);
                }
                DDSIP_bb->bestfirst[j].cursubsol = (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j];
                DDSIP_bb->bestfirst[j].subbound  = (DDSIP_node[DDSIP_bb->curnode]->subbound)[j];
            }
        }
        else if (DDSIP_bb->dualdescitcnt == 1)
        {
            double weight;
            weight = cb_get_last_weight(DDSIP_bb->dualProblem);
            DDSIP_bb->last_weight = weight;
            if (DDSIP_param->cb_increaseWeight && !DDSIP_bb->weight_reset &&
                    ((weight < 4. &&
                      (DDSIP_node[DDSIP_bb->curnode]->bound - tmpbestbound) > 5.e-3 * (fabs(DDSIP_node[DDSIP_bb->curnode]->bound) + 1.)) ||
                     ((DDSIP_node[DDSIP_bb->curnode]->bound - tmpbestbound) > 0.1 * (fabs(DDSIP_node[DDSIP_bb->curnode]->bound) + 1.)))
               )
            {
                original_weight = weight;
                if ((DDSIP_node[DDSIP_bb->curnode]->bound - tmpbestbound) > (3.*fabs(DDSIP_node[DDSIP_bb->curnode]->bound) + 1.))
                {
                    DDSIP_bb->weight_reset = 6;
                    cb_set_next_weight (DDSIP_bb->dualProblem, 6.*weight_reset_factor*original_weight);
                    if (DDSIP_param->outlev > 20)
                        fprintf (DDSIP_bb->moreoutfile, " ### increased weight to %g\n", 6.*weight_reset_factor*original_weight);
                }
                else if ((DDSIP_node[DDSIP_bb->curnode]->bound - tmpbestbound) > (fabs(DDSIP_node[DDSIP_bb->curnode]->bound) + 1.))
                {
                    DDSIP_bb->weight_reset = 2;
                    cb_set_next_weight (DDSIP_bb->dualProblem, 2.*weight_reset_factor*original_weight);
                    if (DDSIP_param->outlev > 20)
                        fprintf (DDSIP_bb->moreoutfile, " ### increased weight to %g\n", 2.*weight_reset_factor*original_weight);
                }
                else
                {
                    DDSIP_bb->weight_reset = 1;
                    cb_set_next_weight (DDSIP_bb->dualProblem, weight_reset_factor*original_weight);
                    if (DDSIP_param->outlev > 20)
                        fprintf (DDSIP_bb->moreoutfile, " ### increased weight to %g\n", weight_reset_factor*original_weight);
                }
            }
        }
    }
    if ((DDSIP_bb->dualdescitcnt == 1) && (DDSIP_bb->weight_reset) && ((DDSIP_node[DDSIP_bb->curnode]->bound - tmpbestbound) < 1.e-8 * (fabs(DDSIP_node[DDSIP_bb->curnode]->bound) + 100.)))
    {
        if (((DDSIP_bb->weight_reset == 1) &&
                (fabs(   weight_reset_factor*original_weight - cb_get_last_weight(DDSIP_bb->dualProblem)) < DDSIP_param->accuracy)) ||
                ((DDSIP_bb->weight_reset == 2) &&
                 (fabs(2.*weight_reset_factor*original_weight - cb_get_last_weight(DDSIP_bb->dualProblem)) < DDSIP_param->accuracy)) ||
                ((DDSIP_bb->weight_reset == 6) &&
                 (fabs(6.*weight_reset_factor*original_weight - cb_get_last_weight(DDSIP_bb->dualProblem)) < DDSIP_param->accuracy)))
        {
            DDSIP_bb->last_weight = DDSIP_bb->weight_reset * original_weight;
            cb_set_next_weight (DDSIP_bb->dualProblem, DDSIP_bb->last_weight);
            if (DDSIP_param->outlev > 20)
                fprintf (DDSIP_bb->moreoutfile, " ### reset weight to %g\n", DDSIP_bb->last_weight);
            DDSIP_bb->weight_reset = -1;
        }
    }
    if (tmpbestbound > DDSIP_bb->bestvalue + 1.e-10 * fabs (DDSIP_bb->bestvalue))
    {
        if (DDSIP_param->outlev)
            fprintf (DDSIP_bb->moreoutfile, "#### Bound reached %20.14g > bestvalue (%20.14g), diff= %g.\n", tmpbestbound, DDSIP_bb->bestvalue, tmpbestbound - DDSIP_bb->bestvalue);
        DDSIP_bb->skip = 2;
        DDSIP_bb->keepSols = 1;
        goto TERMINATE;
    }
    // Evaluate an upper bound if we found a solution for at least one scenario
    // otherwise it will be skipped
    DDSIP_bb->skip = DDSIP_SkipUB ();
    // If none of the scenario problems has a solution then the dispersion norm (max-min) of each variable
    // is set to DDSIP_infty
    if (DDSIP_bb->skip)
    {
        if (DDSIP_param->outlev)
        {
            if (DDSIP_bb->violations)
                fprintf (DDSIP_bb->moreoutfile, "\tLower bound of node %d = %18.16g, skipping upper bound, skip=%d\n", DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound,DDSIP_bb->skip);
            else
            {
                // for non-zero multipliers (as expected in a Conic Bundle iteration the single scenario objective values here are for the
                // changed objective - thus let UpperBound determine them without the Lagrangean terms
                fprintf (DDSIP_bb->moreoutfile, "\tLower and upper bound of node %d = %18.16g, evaluate scenario solutions in upper bound.\n", DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound);
            }
        }
        if (DDSIP_bb->violations && DDSIP_bb->skip > 0 && DDSIP_bb->skip < 3)
        {
            for (j = 0; j < DDSIP_bb->firstvar; j++)
            {
                maxfirst[j] = DDSIP_infty;
                minfirst[j] = -DDSIP_infty;
            }
        }
        else if (!DDSIP_bb->violations && DDSIP_node[DDSIP_bb->curnode]->bound <  DDSIP_bb->bestvalue)
        {
            DDSIP_bb->skip = -1;
        }
    }

#ifdef DEBUG
    //DEBUGOUT
    if (DDSIP_param->outlev > 49)
    {
        double expsubsol, expsubbound, expscenobj, risksubsol, risksubbound, riskscenobj;
        risksubsol = risksubbound = riskscenobj = expsubsol = expsubbound = expscenobj = 0.;
        fprintf(DDSIP_bb->moreoutfile, "\n############### CBLowerBound\n scen  cursubsol,   subbound,   ref_scenobj:\n");
        for(j=0; j<DDSIP_param->scenarios; j++)
        {
            fprintf(DDSIP_bb->moreoutfile, "%3d  %19.16g  %19.16g  %19.16g\n",j+1,(DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j],(DDSIP_node[DDSIP_bb->curnode]->subbound)[j],(DDSIP_node[DDSIP_bb->curnode]->ref_scenobj)[j]);
            expsubsol += DDSIP_data->prob[j] *(DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j];
            expsubbound += DDSIP_data->prob[j] *(DDSIP_node[DDSIP_bb->curnode]->subbound)[j];
            expscenobj += DDSIP_data->prob[j] *(DDSIP_node[DDSIP_bb->curnode]->ref_scenobj)[j];
            if (abs(DDSIP_param->riskmod) == 1)
            {
                if ((DDSIP_node[DDSIP_bb->curnode]->subbound)[j] > DDSIP_param->risktarget)
                    risksubbound += DDSIP_data->prob[j] * (((DDSIP_node[DDSIP_bb->curnode]->subbound)[j] - DDSIP_param->risktarget));
                if ((DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j] > DDSIP_param->risktarget)
                    risksubsol += DDSIP_data->prob[j] * (((DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j] - DDSIP_param->risktarget));
                if ((DDSIP_node[DDSIP_bb->curnode]->ref_scenobj)[j] > DDSIP_param->risktarget)
                    riskscenobj += DDSIP_data->prob[j] * (((DDSIP_node[DDSIP_bb->curnode]->ref_scenobj)[j] - DDSIP_param->risktarget));
            }
            else if (abs(DDSIP_param->riskmod) == 2)
            {
                if ((DDSIP_node[DDSIP_bb->curnode]->subbound)[j] > DDSIP_param->risktarget)
                    risksubbound += DDSIP_data->prob[j];
                if ((DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j] > DDSIP_param->risktarget)
                    risksubsol += DDSIP_data->prob[j];
                if ((DDSIP_node[DDSIP_bb->curnode]->ref_scenobj)[j] > DDSIP_param->risktarget)
                    riskscenobj += DDSIP_data->prob[j];
            }
        }
        fprintf(DDSIP_bb->moreoutfile, "exp  %19.16g  %19.16g  %19.16g\n",expsubsol,expsubbound,expscenobj);
        fprintf(DDSIP_bb->moreoutfile, "risk %19.16g  %19.16g  %19.16g\n",risksubsol,risksubbound,riskscenobj);
        fprintf(DDSIP_bb->moreoutfile, "sum  %19.16g  %19.16g  %19.16g\n",expsubsol + DDSIP_param->riskweight * risksubsol,expsubbound + DDSIP_param->riskweight * risksubbound, expscenobj + DDSIP_param->riskweight * riskscenobj);
        fprintf(DDSIP_bb->moreoutfile, "###############\n");
    }
    //DEBUGOUT
#endif

    if (DDSIP_bb->skip)
    {
        if (DDSIP_param->outlev)
        {
            if (DDSIP_bb->violations)
                fprintf (DDSIP_bb->moreoutfile, "\tLower bound of node %d = %18.16g, skipping upper bound, skip=%d\n", DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound,DDSIP_bb->skip);
            else
            {
                if (DDSIP_bb->skip == 3)
                {
                    if (DDSIP_node[DDSIP_bb->curnode]->bound < DDSIP_bb->bestvalue)
                        DDSIP_bb->skip = -1;
                    else
                        DDSIP_bb->skip = 0;
                }
                else
                    fprintf (DDSIP_bb->moreoutfile, "\tLower and upper bound of node %d = %18.16g, skipping upper bound, skip=%d\n", DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound,DDSIP_bb->skip);
            }
        }
        if (DDSIP_bb->violations && DDSIP_bb->skip < 3)
        {
            for (j = 0; j < DDSIP_bb->firstvar; j++)
            {
                maxfirst[j] = DDSIP_infty;
                minfirst[j] = 0;
            }
        }
    }
    if (!DDSIP_bb->violations && !increase && DDSIP_bb->bestfirst[0].first_sol && !no_viol_worse)
    {
        // We have no violations, but a worse objective value - a zero subgradient would erroneously(?) be the result
        // Reset first stage solutions to the ones that gave the best bound
        // use the CPLEX settings for LowerBound (probably tighter tols)
        no_viol_worse = 1;
        if (meanGap)
            use_LB_params = 1;
        if (DDSIP_param->outlev)
        {
            printf ("***  no violations, but a worse objective value");
            fprintf (DDSIP_bb->moreoutfile, "***  no violations, but a worse objective value");
        }
        if (tmpupper > DDSIP_node[DDSIP_bb->curnode]->bound)
        {
            if (DDSIP_param->outlev)
            {
                printf (" - probably due to mipgaps\n");
                fprintf (DDSIP_bb->moreoutfile, " - probably due to mipgaps\n");
            }
        }
        else
        {
            if (DDSIP_param->outlev)
            {
                printf (" - reset first stage solutions to the ones that gave the best bound\n");
                fprintf (DDSIP_bb->moreoutfile, " - reset first stage solutions to the ones that gave the best bound\n");
            }
            // Free previous first stage
            for (j = 0; j < DDSIP_param->scenarios; j++)
            {
                if (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j])
                {
                    if ((cnt = (int) (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j][DDSIP_bb->firstvar] - 0.9)))
                    {
                        for (k1 = j + 1; cnt && k1 < DDSIP_param->scenarios; k1++)
                        {
                            if (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j] == ((DDSIP_node[DDSIP_bb->curnode])->first_sol)[k1])
                            {
                                ((DDSIP_node[DDSIP_bb->curnode])->first_sol)[k1] = NULL;
                                cnt--;
                            }
                        }
                    }
                    DDSIP_Free ((void **) &(((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j]));
                }
            }
            // Copy first stage in best point to first_sol
            for (k1 = 0; k1 < DDSIP_bb->firstvar; k1++)
            {
                // Initialize minimum and maximum of each component
                minfirst[k1] = maxfirst[k1] = DDSIP_bb->bestfirst[0].first_sol[k1];
            }
            for (j = 0; j < DDSIP_param->scenarios; j++)
            {
                if (!((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j])
                {
                    ((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j] = (double *) DDSIP_Alloc(sizeof(double),DDSIP_bb->firstvar+3,"first_sol(CBLowerBound)");
                    memcpy (DDSIP_node[DDSIP_bb->curnode]->first_sol[j], DDSIP_bb->bestfirst[j].first_sol, sizeof (double) * (DDSIP_bb->firstvar + 3));
#ifdef DEBUG
                    if (DDSIP_param->outlev > 90)
                        fprintf (DDSIP_bb->moreoutfile,
                                 "       copied solution for scenario %d from bestfirst.\n", j);
#endif
                    if ((cnt = (int) ((DDSIP_bb->bestfirst[j].first_sol)[DDSIP_bb->firstvar] - 0.9)))
                        for (k1 = j + 1; cnt && k1 < DDSIP_param->scenarios; k1++)
                        {
                            if (DDSIP_bb->bestfirst[j].first_sol == DDSIP_bb->bestfirst[k1].first_sol)
                            {
                                ((DDSIP_node[DDSIP_bb->curnode])->first_sol)[k1] = ((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j];
                                cnt--;
#ifdef DEBUG
                                if (DDSIP_param->outlev > 90)
                                    fprintf (DDSIP_bb->moreoutfile,
                                             "       same address for scenarios %d and %d, %d identical remaining.\n", j+1, k1+1,cnt);
#endif
                            }
                        }
                }
                DDSIP_bb->bestfirst[j].cursubsol = (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j];
                DDSIP_bb->bestfirst[j].subbound  = (DDSIP_node[DDSIP_bb->curnode]->subbound)[j];
                for (k1 = 0; k1 < DDSIP_bb->firstvar; k1++)
                {
                    // Calculate minimum and maximum of each component
                    minfirst[k1] = DDSIP_Dmin ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[j][k1], minfirst[k1]);
                    maxfirst[k1] = DDSIP_Dmax ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[j][k1], maxfirst[k1]);
                }
            }
        }
    }
    // Count again the number of differences within first stage solution in current node
    // do this for the bestbound solution in case it was a step with zero multipliers
    // DDSIP_bb->violations=k means differences in k components
    if (DDSIP_param->outlev > 34)
    {
        colname = (char **) DDSIP_Alloc (sizeof (char *), (DDSIP_bb->firstvar + DDSIP_bb->secvar), "colname(LowerBound)");
        colstore = (char *) DDSIP_Alloc (sizeof (char), (DDSIP_bb->firstvar + DDSIP_bb->secvar) * DDSIP_ln_varname, "colstore(LowerBound)");
        status =
            CPXgetcolname (DDSIP_env, DDSIP_lp, colname, colstore,
                           (DDSIP_bb->firstvar + DDSIP_bb->secvar) * DDSIP_ln_varname, &j, 0, DDSIP_bb->firstvar + DDSIP_bb->secvar - 1);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to get column names (UB)\n");
            fprintf (DDSIP_outfile, "ERROR: Failed to get column names (UB)\n");
            if (DDSIP_param->outlev)
                fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to get column names (UB)\n");
            goto TERMINATE;
        }
    }
    DDSIP_bb->violations = 0;
    maxdispersion = 0.;
    for (j = 0; j < DDSIP_bb->firstvar; j++)
    {
        if (!DDSIP_Equal (maxfirst[j], minfirst[j]))
        {
            DDSIP_bb->violations++;
            maxdispersion = DDSIP_Dmax (maxdispersion, maxfirst[j]-minfirst[j]);
            // More Output
            if (DDSIP_param->outlev > 34)
                fprintf (DDSIP_bb->moreoutfile, "Deviation of variable %5d: min=%16.14g\t max=%16.14g\t diff=%16.14g  \t%10.10s\tmaxdisp: %g\n", DDSIP_bb->firstindex[j], minfirst[j], maxfirst[j], maxfirst[j] - minfirst[j], colname[DDSIP_bb->firstindex[j]], maxdispersion);
        }
    }
    if (DDSIP_param->outlev > 34)
    {
        DDSIP_Free ((void **) &(colstore));
        DDSIP_Free ((void **) &(colname));
    }
    DDSIP_node[DDSIP_bb->curnode]->dispnorm = maxdispersion;

    if (DDSIP_bb->dualObjVal < -(*objective_val))
    {
        DDSIP_bb->objIncrease = 1;
        if (DDSIP_param->outlev > 6 && !increase && DDSIP_bb->dualObjVal > -DDSIP_infty)
        {
            fprintf (DDSIP_bb->moreoutfile,
                     " ***++ dual step     increasing objval for node %3d, new val: %-18.16g, old value: %-18.16g (incr:  %g) ++***\n",
                     DDSIP_bb->curnode, -(*objective_val), DDSIP_bb->dualObjVal,  -(*objective_val)-DDSIP_bb->dualObjVal);
        }
        if (!DDSIP_bb->dualdescitcnt && (DualObj <= tmpbestbound) && DDSIP_param->hot)
        {
            for (j = 0; j < DDSIP_param->scenarios; j++)
                memcpy (DDSIP_bb->boundIncrease_intsolvals[j], DDSIP_bb->intsolvals[j], DDSIP_bb->total_int*sizeof(double));
            if (DDSIP_param->outlev > 21)
                fprintf (DDSIP_bb->moreoutfile, "##  Saved mip start info of objective increase\n");
        }
    }
    else
    {
        DDSIP_bb->objIncrease = 0;
        if (DDSIP_param->outlev > 6 && !increase && DDSIP_bb->dualObjVal > -DDSIP_infty)
        {
            fprintf (DDSIP_bb->moreoutfile,
                     " ---** dual step not increasing objval for node %3d, new val: %-18.16g, old value: %-18.16g (diff: %g) **---\n",
                     DDSIP_bb->curnode, -(*objective_val), DDSIP_bb->dualObjVal,  -(*objective_val)-DDSIP_bb->dualObjVal);
        }
    }
    if (!(DDSIP_bb->dualitcnt))
        DualObj = tmpbestbound;
    else
        DualObj = DDSIP_Dmax (tmpbestbound, DualObj);
    DDSIP_bb->dualObjVal =  DDSIP_Dmax(-(*objective_val), DDSIP_bb->dualObjVal);
    DDSIP_bb->meanGapCBLB = DDSIP_Dmax(DDSIP_bb->meanGapCBLB, meanGap);
    // More debugging information
    if (DDSIP_param->outlev)
    {
        if (DDSIP_bb->violations)
        {
            if (DDSIP_param->outlev > 3)
                fprintf (DDSIP_bb->moreoutfile,
                         "\tNumber of violations of nonanticipativity of first-stage variables:  %d, max: %g\n", DDSIP_bb->violations, maxdispersion);
            if (DDSIP_param->outlev > 7)
            {
                printf ("\tNumber of violations of nonanticipativity of first-stage variables:  %d, max: %g\n", DDSIP_bb->violations, maxdispersion);
                printf ("\tDual bound of node  %-10d           = \t%-18.16g\n",
                        DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound);
                printf ("\tCurrent dual objective value             = \t%-18.16g          \t(mean MIP gap: %g%%, max MIP gap: %g%%)\n", -(*objective_val), meanGap, maxGap);
            }
        }
        else
        {
            if (DDSIP_param->outlev > 3)
                fprintf (DDSIP_bb->moreoutfile,
                         "\tNumber of violations of nonanticipativity of first-stage variables:  %d\n", DDSIP_bb->violations);
            if (DDSIP_param->outlev > 7)
            {
                printf ("\tNumber of violations of nonanticipativity of first-stage variables:  %d\n", DDSIP_bb->violations);
                printf ("\tDual bound of node %-10d            = \t%-18.16g\n",
                        DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound);
                printf ("\tCurrent dual objective value             = \t%-18.16g          \t(mean MIP gap: %g%%)\n", -(*objective_val), meanGap);
                printf ("\tBest dual objective value in descent step= \t%-18.16g\n", DDSIP_bb->dualObjVal);
            }
        }
        fprintf (DDSIP_bb->moreoutfile, "\tDual bound of node %-10d            = \t%-18.16g\n",
                 DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound);
        fprintf (DDSIP_bb->moreoutfile, "\tCurrent dual objective value             = \t%-18.16g          \t(mean MIP gap: %g%%, max gap: %g%%)\n", -(*objective_val), meanGap, maxGap);
        fprintf (DDSIP_bb->moreoutfile, "\tBest dual objective value in descent step= \t%-18.16g\n", DDSIP_bb->dualObjVal);
    }
    if (DDSIP_param->cb_increaseWeight)
    {
        if (DDSIP_bb->dualdescitcnt && fabs (-(*objective_val) - old_value)/(fabs (old_value) + 1.e-16) < 1.e-17)
        {
            totally_constant++;
            if (DDSIP_param->outlev > 20)
                fprintf (DDSIP_bb->moreoutfile, "\tobj change: %g, totally_constant: %d\n", fabs (-(*objective_val) - old_value)/(fabs (old_value) + 1.e-16), totally_constant);
        }
        if (DDSIP_bb->dualdescitcnt && (fabs (-(*objective_val) - old_value)/(fabs(old_value) + 1.e-16) < 1.e-12 || fabs (-(*objective_val) - old_value) < 1.e-14))
        {
            nearly_constant++;
            if (DDSIP_param->outlev > 20)
                fprintf (DDSIP_bb->moreoutfile, "##\tobj change: %g, nearly_constant: %d\n", fabs (-(*objective_val) - old_value)/(fabs(old_value) + 1.e-16), nearly_constant);
        }
        if (totally_constant > 2)
        {
            wr =  DDSIP_Dmax (1.0e4, 10.*cb_get_last_weight(DDSIP_bb->dualProblem));
            if (DDSIP_param->outlev > 20)
                fprintf (DDSIP_bb->moreoutfile, "####\tset new weight to %g, totally_constant:%d (nearly_constant: %d)\n", wr, totally_constant, nearly_constant);
            cb_set_next_weight (DDSIP_bb->dualProblem, wr);
            nearly_constant = totally_constant = 0;
        }
        else if (nearly_constant > 4)
        {
            if (cb_get_last_weight(DDSIP_bb->dualProblem) < 1e-3)
                wr =  10.*cb_get_last_weight(DDSIP_bb->dualProblem);
            else if (cb_get_last_weight(DDSIP_bb->dualProblem) < 1e-2)
                wr =  5.*cb_get_last_weight(DDSIP_bb->dualProblem);
            else
                wr =  0.4*cb_get_last_weight(DDSIP_bb->dualProblem);
            if (DDSIP_param->outlev > 20)
                fprintf (DDSIP_bb->moreoutfile, "####\tset new weight to %g, nearly_constant: %d (totally_constant:%d)\n", wr, nearly_constant, totally_constant);
            cb_set_next_weight (DDSIP_bb->dualProblem, wr);
            nearly_constant = totally_constant = 0;
        }
    }
    old_value = -(*objective_val);
TERMINATE:

    // Only if no errors, yet.
    if (!status || status == -111)
    {
        relax = DDSIP_RestoreBoundAndType ();
        if (relax)
            fprintf (stderr, "ERROR: Failed to restore bounds \n");
    }
    if (!status)
        DDSIP_bb->newTry = 0;

    DDSIP_Free ((void **) &(indices));
    DDSIP_Free ((void **) &(mipx));
    DDSIP_Free ((void **) &(minfirst));
    DDSIP_Free ((void **) &(maxfirst));
    DDSIP_Free ((void **) &(type));
    if (DDSIP_param->outlev > 11)
    {
        if (status)
            fprintf (DDSIP_bb->moreoutfile, "#########§§§§§§§§§  CBLB return value = %d  §§§§§§§§§#########\n", status);
    }
    return status;
#else
    return 0;
#endif
} // DDSIP_CBLowerBound
