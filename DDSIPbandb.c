/*  Authors:           Andreas M"arkert, Ralf Gollmer
    Copyright to:      University of Duisburg-Essen
    Language:          C
    Description:
    The main procedures in this file are called `Bound' and `Branch'.
    The first one is invoked to delete nodes in the branch-and-bound tree
    after a new upper bound has been obtained. The second procedure creates
    new nodes in the tree and prepares the evaluation of these nodes.

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

//#define MDEBUG

#include <DDSIP.h>
#include <DDSIPconst.h>


static int DDSIP_Leaf (void);
static int DDSIP_GetCurNode (void);
static int DDSIP_InitNewNodes (void);
static int DDSIP_SetBounds (void);

//==========================================================================
// Select node and value for branching
int
DDSIP_GetCurNode (void)
{
    int i, change = 0;

    // Unsolved nodes have priority (width first)
    if (DDSIP_param->boundstrat <= 1 || DDSIP_param->boundstrat == 3 || DDSIP_param->boundstrat > 4)
    {
        for (i = 0; i < DDSIP_bb->nofront; i++)
            if (!DDSIP_node[DDSIP_bb->front[i]]->solved)
            {
                DDSIP_bb->curnode = DDSIP_bb->front[i];
                change = 1;
                i = DDSIP_bb->nofront;
            }
    }
    if (!change)
    {
        for (i = 0; i < DDSIP_bb->nofront; i++)
        {
            if (!DDSIP_Equal (DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->dispnorm, DDSIP_param->nulldisp) && !DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->leaf)
            {
                DDSIP_bb->curnode = DDSIP_bb->front_nodes_sorted[i];
                change = 1;
                break;
            }
        }
    }

    if (!change)
        return 113;

    return 0;
}

//==========================================================================
// Function initializes new nodes in b&b-tree
int
DDSIP_InitNewNodes (void)
{
    int i, j, status, scen, cnt;
    double branchval, lhs;
    cutpool_t *currentCut;

    if (DDSIP_param->outlev > 2)
    {
        fprintf (DDSIP_bb->moreoutfile, "\n----------------------\n");
        fprintf (DDSIP_bb->moreoutfile, "Branching node %d...\n", DDSIP_bb->curnode);
    }
    DDSIP_node[DDSIP_bb->nonode] = (node_t *) DDSIP_Alloc (sizeof (node_t), 1, "DDSIP_node[nonode](InitNewNodes)");
    DDSIP_node[DDSIP_bb->nonode]->first_sol =
        (double **) DDSIP_Alloc (sizeof (double *), DDSIP_param->scenarios, "DDSIP_node[nonode]->first_sol(InitNewNodes)");
    for (i=0; i<DDSIP_param->scenarios;  DDSIP_node[DDSIP_bb->nonode]->first_sol[i++]=NULL);
    DDSIP_node[DDSIP_bb->nonode]->cursubsol =
        (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "DDSIP_node[nonode]->cursubsol(InitNewNodes)");
    DDSIP_node[DDSIP_bb->nonode]->mipstatus = (int *) DDSIP_Alloc (sizeof (int), DDSIP_param->scenarios, "DDSIP_node[nonode]->mipstatus(InitNewNodes)");
    DDSIP_node[DDSIP_bb->nonode]->ref_scenobj =
        (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "DDSIP_node[nonode]->ref_scenobj(InitNewNodes)");
    DDSIP_node[DDSIP_bb->nonode + 1] = (node_t *) DDSIP_Alloc (sizeof (node_t), 1, "DDSIP_node[nonode](InitNewNodes)");
    DDSIP_node[DDSIP_bb->nonode + 1]->first_sol =
        (double **) DDSIP_Alloc (sizeof (double *), DDSIP_param->scenarios, "DDSIP_node[nonode+1]->first_sol(InitNewNodes)");
    for (i=0; i<DDSIP_param->scenarios;  DDSIP_node[DDSIP_bb->nonode + 1]->first_sol[i++]=NULL);
    DDSIP_node[DDSIP_bb->nonode + 1]->cursubsol =
        (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "DDSIP_node[nonode+1]->cursubsol(InitNewNodes)");
    DDSIP_node[DDSIP_bb->nonode + 1]->mipstatus =
        (int *) DDSIP_Alloc (sizeof (int), DDSIP_param->scenarios, "DDSIP_node[nonode]->mipstatus(InitNewNodes)");
    DDSIP_node[DDSIP_bb->nonode + 1]->ref_scenobj =
        (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "DDSIP_node[nonode+1]->ref_scenobj(InitNewNodes)");
    for (i = 0; i < DDSIP_param->scenarios; i++)
    {
        (DDSIP_node[DDSIP_bb->nonode]->cursubsol)[i] = DDSIP_infty;
        (DDSIP_node[DDSIP_bb->nonode + 1]->cursubsol)[i] = DDSIP_infty;
    }
    // Set father and depth
    DDSIP_node[DDSIP_bb->nonode]->father = DDSIP_bb->curnode;
    DDSIP_node[DDSIP_bb->nonode + 1]->father = DDSIP_bb->curnode;
    DDSIP_node[DDSIP_bb->nonode]->depth = DDSIP_node[DDSIP_bb->curnode]->depth +1;
    DDSIP_node[DDSIP_bb->nonode + 1]->depth = DDSIP_node[DDSIP_bb->curnode]->depth +1;

    // Initialize dispersion norm
    DDSIP_node[DDSIP_bb->nonode]->dispnorm = DDSIP_infty;
    DDSIP_node[DDSIP_bb->nonode + 1]->dispnorm = DDSIP_infty;
    // initailize violations
    DDSIP_node[DDSIP_bb->nonode]->violations = DDSIP_node[DDSIP_bb->curnode]->violations;
    DDSIP_node[DDSIP_bb->nonode + 1]->violations = DDSIP_node[DDSIP_bb->curnode]->violations;

    DDSIP_node[DDSIP_bb->nonode]->solved = 0;
    DDSIP_node[DDSIP_bb->nonode + 1]->solved = 0;

    DDSIP_node[DDSIP_bb->nonode]->cutAdded = 0;
    DDSIP_node[DDSIP_bb->nonode + 1]->cutAdded = 0;
    if (DDSIP_node[DDSIP_bb->curnode]->cbReturn32 && DDSIP_node[DDSIP_bb->curnode]->depth < DDSIP_node[DDSIP_bb->curnode]->cbReturn32 + 3)
    {
        DDSIP_node[DDSIP_bb->nonode]->cbReturn32 = DDSIP_node[DDSIP_bb->curnode]->cbReturn32;
        DDSIP_node[DDSIP_bb->nonode + 1]->cbReturn32 = DDSIP_node[DDSIP_bb->curnode]->cbReturn32;
        if (DDSIP_param->outlev > 20)
            fprintf (DDSIP_bb->moreoutfile, "##########  DDSIP_node[%d]->depth= %d < ->cbReturn32 = %d + 3 -> nodes %d and %d inherit cbReturn\n",
                     DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->depth, DDSIP_node[DDSIP_bb->curnode]->cbReturn32, DDSIP_bb->nonode, DDSIP_bb->nonode+1);
    }
    else
    {
        DDSIP_node[DDSIP_bb->nonode]->cbReturn32 = 0;
        DDSIP_node[DDSIP_bb->nonode + 1]->cbReturn32 = 0;
    }

    if (DDSIP_param->hot)
    {
        // DDSIP_Allocate memory for subsolutions and bound
        DDSIP_node[DDSIP_bb->nonode]->solut =
            (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios * DDSIP_bb->total_int, "solut(InitNewNodes)");
        // Initialization
        memcpy (DDSIP_node[DDSIP_bb->nonode]->solut, DDSIP_node[DDSIP_bb->curnode]->solut, sizeof (double) * DDSIP_param->scenarios * DDSIP_bb->total_int);
        DDSIP_node[DDSIP_bb->nonode + 1]->solut = DDSIP_node[DDSIP_bb->curnode]->solut;
        DDSIP_node[DDSIP_bb->curnode]->solut = NULL;
    }
    // Initialization of lower bounds
    DDSIP_node[DDSIP_bb->nonode]->subbound =
        (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "DDSIP_node[nonode]->subbound(InitNewNodes)");
    memcpy (DDSIP_node[DDSIP_bb->nonode]->subbound, DDSIP_node[DDSIP_bb->curnode]->subbound, sizeof (double) * DDSIP_param->scenarios);
    //DDSIP_node[DDSIP_bb->nonode + 1]->subbound =
    //  (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "DDSIP_node[nonode+1]->subbound(InitNewNodes)");
    //memcpy (DDSIP_node[DDSIP_bb->nonode + 1]->subbound, DDSIP_node[DDSIP_bb->curnode]->subbound, sizeof (double) * DDSIP_param->scenarios);
    DDSIP_node[DDSIP_bb->nonode + 1]->subbound = DDSIP_node[DDSIP_bb->curnode]->subbound;
    DDSIP_node[DDSIP_bb->curnode]->subbound = NULL;
#ifdef CONIC_BUNDLE
    // Initialize multiplier in node
    if (DDSIP_param->cb)
    {
        DDSIP_node[DDSIP_bb->nonode]->dual = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->dimdual + 3, "dual(InitNewNodes)");
        memcpy (DDSIP_node[DDSIP_bb->nonode]->dual, DDSIP_node[DDSIP_bb->curnode]->dual, sizeof (double) * (DDSIP_bb->dimdual + 3));
        DDSIP_node[DDSIP_bb->nonode + 1]->dual = DDSIP_node[DDSIP_bb->curnode]->dual;
        DDSIP_node[DDSIP_bb->curnode]->dual = NULL;

        DDSIP_node[DDSIP_bb->nonode]->scenBoundsNoLag =
            (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "DDSIP_node[nonode]->scenBoundsNoLag(InitNewNodes)");
        memcpy (DDSIP_node[DDSIP_bb->nonode]->scenBoundsNoLag, DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag, sizeof (double) * DDSIP_param->scenarios);
        DDSIP_node[DDSIP_bb->nonode + 1]->scenBoundsNoLag = DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag;
        DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag = NULL;
        DDSIP_node[DDSIP_bb->nonode]->BoundNoLag = DDSIP_node[DDSIP_bb->nonode + 1]->BoundNoLag = DDSIP_node[DDSIP_bb->curnode]->BoundNoLag;
    }
#endif
    // Absolute semideviation
    if (DDSIP_param->riskmod == 3)
    {
        DDSIP_node[DDSIP_bb->nonode]->target = DDSIP_node[DDSIP_bb->curnode]->target;
        DDSIP_node[DDSIP_bb->nonode + 1]->target = DDSIP_node[DDSIP_bb->curnode]->target;
    }
    // Inherit bounds
    DDSIP_node[DDSIP_bb->nonode]->bound = DDSIP_node[DDSIP_bb->curnode]->bound;
    DDSIP_node[DDSIP_bb->nonode + 1]->bound = DDSIP_node[DDSIP_bb->curnode]->bound;

    DDSIP_node[DDSIP_bb->nonode]->branchind = DDSIP_node[DDSIP_bb->curnode]->branchind;
    DDSIP_node[DDSIP_bb->nonode + 1]->branchind = DDSIP_node[DDSIP_bb->curnode]->branchind;

    // Initialize leaf DDSIP_indicator (no leaf)
    DDSIP_node[DDSIP_bb->nonode]->leaf = 0;
    DDSIP_node[DDSIP_bb->nonode + 1]->leaf = 0;

    DDSIP_node[DDSIP_bb->nonode]->numInheritedSols = 0;
    DDSIP_node[DDSIP_bb->nonode + 1]->numInheritedSols = 0;

    // Add new constraints
    // The branch index 'branchind' was calculated
    // in LowerBound based on the maximal dispersion of variables etc.
    // The index of the new bound in the new nodes thus becomes branchind
    DDSIP_node[DDSIP_bb->nonode]->neoind = DDSIP_node[DDSIP_bb->curnode]->branchind;
    DDSIP_node[DDSIP_bb->nonode + 1]->neoind = DDSIP_node[DDSIP_bb->curnode]->branchind;

    // The lower bound of the 'left' node and the upper bound of the 'right' are
    // initialized by the original lower bounds of the according component
    DDSIP_node[DDSIP_bb->nonode]->neolb = DDSIP_bb->lborg[DDSIP_node[DDSIP_bb->curnode]->branchind];
    DDSIP_node[DDSIP_bb->nonode + 1]->neoub = DDSIP_bb->uborg[DDSIP_node[DDSIP_bb->curnode]->branchind];
    //fprintf (stderr," bandb: new nodes variable %d original bounds  lb= %g, ub= %g  branchval= %g\n",DDSIP_node[DDSIP_bb->curnode]->branchind,DDSIP_node[DDSIP_bb->nonode]->neolb,DDSIP_node[DDSIP_bb->nonode + 1]->neoub, DDSIP_node[DDSIP_bb->curnode]->branchval);

    // If the variable has already been branched on, the lower bound of the 'left'
    // node and the upper bound of the 'right' need to be updated.
    i = DDSIP_bb->curnode;
    while (i > 0)
    {
        if (DDSIP_node[i]->neoind == DDSIP_node[DDSIP_bb->curnode]->branchind)
        {
            DDSIP_node[DDSIP_bb->nonode]->neolb = DDSIP_Dmax (DDSIP_node[i]->neolb, DDSIP_node[DDSIP_bb->nonode]->neolb);
            DDSIP_node[DDSIP_bb->nonode + 1]->neoub = DDSIP_Dmin (DDSIP_node[i]->neoub, DDSIP_node[DDSIP_bb->nonode + 1]->neoub);
            if (DDSIP_param->outlev)
                fprintf (DDSIP_bb->moreoutfile," bandb: new nodes variable %d update in node %d lb= %.16g, ub= %.16g\n",DDSIP_node[DDSIP_bb->curnode]->branchind,i,DDSIP_node[DDSIP_bb->nonode]->neolb,DDSIP_node[DDSIP_bb->nonode + 1]->neoub);
        }
        i = DDSIP_node[i]->father;
    }
    if (DDSIP_node[DDSIP_bb->nonode]->neolb > DDSIP_node[DDSIP_bb->nonode + 1]->neoub)
    {
        if (DDSIP_param->outlev > 2)
            fprintf (DDSIP_bb->moreoutfile,
                     "ERROR: Correcting bounds of new nodes (bandb): was: lb=%.16g, ub=%.16g, diff=%g\n",
                     DDSIP_node[DDSIP_bb->nonode]->neolb, DDSIP_node[DDSIP_bb->nonode + 1]->neoub,
                     (DDSIP_node[DDSIP_bb->nonode]->neolb) - (DDSIP_node[DDSIP_bb->nonode + 1]->neoub));
        DDSIP_node[DDSIP_bb->nonode]->neolb = DDSIP_node[DDSIP_bb->nonode + 1]->neoub;
    }

    // branchval = the value chosen in DDSIPlb according to branchstrat
    branchval = DDSIP_node[DDSIP_bb->curnode]->branchval;
 
    // Correction needed ?
    if (branchval < DDSIP_node[DDSIP_bb->nonode]->neolb || branchval > DDSIP_node[DDSIP_bb->nonode + 1]->neoub)
    {
        if (DDSIP_param->outlev > 10)
            fprintf (DDSIP_bb->moreoutfile, "Correcting branchval (bandb): %.16g : neolb= %.16g, neoub= %.16g\n", branchval, DDSIP_node[DDSIP_bb->nonode]->neolb, DDSIP_node[DDSIP_bb->nonode + 1]->neoub);
        branchval = DDSIP_Dmax (DDSIP_Dmin (branchval, DDSIP_node[DDSIP_bb->nonode + 1]->neoub - 1.e-15), DDSIP_node[DDSIP_bb->nonode]->neolb + 1.e-15);
        if (DDSIP_param->outlev > 10)
            fprintf (DDSIP_bb->moreoutfile, "Corrected  branchval (bandb): %.16g\n", branchval);
    }
    // The remaining bounds (upper bd. of left node and lower bd. of right node)
    // depend on branchval
    // Integer or binary variable
    if (DDSIP_bb->firsttype[DDSIP_node[DDSIP_bb->curnode]->branchind] == 'B' || DDSIP_bb->firsttype[DDSIP_node[DDSIP_bb->curnode]->branchind] == 'I' || DDSIP_bb->firsttype[DDSIP_node[DDSIP_bb->curnode]->branchind] == 'N')
    {
        DDSIP_node[DDSIP_bb->nonode]->neoub = DDSIP_Dmax (ceil (branchval) - 1, DDSIP_node[DDSIP_bb->nonode]->neolb);
        DDSIP_node[DDSIP_bb->nonode + 1]->neolb = DDSIP_node[DDSIP_bb->nonode]->neoub + 1;
    }
    // Continuous variable
    // Recall that ub (i+1)-lb(i) > eps and lb(i) <= branchval <= ub(i+1)
    else
    {
        double minus = branchval - DDSIP_Dmax (DDSIP_param->brancheps, fabs(branchval)*1.5e-16);
        double plus  = branchval + DDSIP_Dmax (DDSIP_param->brancheps, fabs(branchval)*1.5e-16);
        // Case: lb (i) <= branchval-eps < branchval <= ub(i+1)
        // -->   ub (i):=branchval-eps, lb(i+1)=branchval
        if (!(minus < DDSIP_node[DDSIP_bb->nonode]->neolb))
        {
            DDSIP_node[DDSIP_bb->nonode]->neoub = minus;
            //take care for numerical errors with very small DDSIP_param->brancheps
            if (branchval < DDSIP_node[DDSIP_bb->nonode + 1]->neoub)
                DDSIP_node[DDSIP_bb->nonode + 1]->neolb = branchval;
            else
                DDSIP_node[DDSIP_bb->nonode + 1]->neolb = DDSIP_node[DDSIP_bb->nonode + 1]->neoub;
        }
        // Case: lb (i) <= branchval < branchval+eps <= ub(i+1)
        // -->   ub (i):=branchval, lb(i+1)=branchval+eps
        else if (!(plus > DDSIP_node[DDSIP_bb->nonode + 1]->neoub))
        {
            //take care for numerical errors with very small DDSIP_param->brancheps
            if (branchval > DDSIP_node[DDSIP_bb->nonode]->neolb)
                DDSIP_node[DDSIP_bb->nonode]->neoub = branchval;
            else
                DDSIP_node[DDSIP_bb->nonode]->neoub = DDSIP_node[DDSIP_bb->nonode]->neolb;
            DDSIP_node[DDSIP_bb->nonode + 1]->neolb = plus;
        }
        // This is the case when ub-lb < 2*eps
        // Last branch on this variable.
        //  --> ub (i):=branchval, lb(i+1)=branchval+.5*(ub(i+1)-branchval)
        // Together with the two cases above we hopefully obtain 0<lb (i+1)-ub(i)<eps
        else
        {
            DDSIP_node[DDSIP_bb->nonode]->neoub = branchval;
            DDSIP_node[DDSIP_bb->nonode + 1]->neolb = branchval + .5 * (DDSIP_node[DDSIP_bb->nonode + 1]->neoub - branchval);
        }

        // We don't want to `lose' the current best value
        if (DDSIP_node[DDSIP_bb->nonode]->neoub < DDSIP_bb->bestsol[DDSIP_node[DDSIP_bb->nonode]->neoind]
                && DDSIP_node[DDSIP_bb->nonode + 1]->neolb > DDSIP_bb->bestsol[DDSIP_node[DDSIP_bb->nonode]->neoind])
            DDSIP_node[DDSIP_bb->nonode]->neoub = DDSIP_bb->bestsol[DDSIP_node[DDSIP_bb->nonode]->neoind];
    }

    if (DDSIP_node[DDSIP_bb->nonode]->neoub > DDSIP_node[DDSIP_bb->nonode + 1]->neolb)
        return 119;

    // Inherit scenario solutions of father node
    for (i = 0; i < DDSIP_param->scenarios; i++)
    {
        if (!(((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i]))
            continue;
        cnt = 0;
        // if cuts have been inserted test the scenario solution for violation
        if (DDSIP_bb->cutpool && (DDSIP_node[DDSIP_bb->curnode])->cutAdded)
        {
            currentCut = DDSIP_bb->cutpool;
            while (currentCut)
            {
                lhs = 0.;
                for (j = 0; j < DDSIP_bb->firstvar; j++)
                {
                    lhs += (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i])[j] * currentCut->matval[j];
                }
                if ((lhs -= currentCut->rhs) < - 1e-7)
                {
                    if (DDSIP_param->outlev > 23)
                        fprintf (DDSIP_bb->moreoutfile,"  nodes %d and %d did not inherit solution of scenario %d from node %d due to added cut %d, violation %g.\n",
                                 DDSIP_bb->nonode, DDSIP_bb->nonode + 1, i+1, DDSIP_bb->curnode, currentCut->number, -lhs);
                    if ((cnt = (int) ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i])[DDSIP_bb->firstvar] - 0.9)))
                        for (j = i + 1; cnt && j < DDSIP_param->scenarios; j++)
                        {
                            if (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j]
                                    && ((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i] == ((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j])
                            {
                                ((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j] = NULL;
                                cnt--;
                            }
                        }
                    DDSIP_Free ((void **) &(((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i]));
                    cnt = 1;
                    break;
                }
                currentCut = currentCut->prev;
            }
        }
        // skip to next scenario if the current solution violated one of the cuts
        if (cnt)
            continue;
        if (DDSIP_param->outlev > 79)
        {
            fprintf (DDSIP_bb->moreoutfile,"  inherit choice: DDSIP_node[%d]->first_sol[%d][%d]= %.15g, DDSIP_node[%d]->neoub: %.15g  DDSIP_node[%d]->neolb: %.15g inherit level %g mipstatus %d\n",DDSIP_bb->curnode,i,DDSIP_node[DDSIP_bb->nonode]->neoind, ((DDSIP_node[DDSIP_bb->curnode])->first_sol[i])[DDSIP_node[DDSIP_bb->nonode]->neoind], DDSIP_bb->nonode, DDSIP_node[DDSIP_bb->nonode]->neoub, DDSIP_bb->nonode + 1, DDSIP_node[DDSIP_bb->nonode + 1]->neolb,(DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][DDSIP_bb->firstvar + 1],(DDSIP_node[DDSIP_bb->curnode]->mipstatus)[i]);
        }
        // due to the change of the lower bound of additional variable in root node for risk model DDSIP_4 (worst case cost) and 5 (TVaR):
        //   do not pass on root node solutions in this case
        // for asd model the same applies in all nodes due to changes of target
        if ((abs(DDSIP_param->riskmod) != 3 && (DDSIP_bb->curnode || (abs(DDSIP_param->riskmod) != 4 && abs(DDSIP_param->riskmod) != 5)) && (DDSIP_node[DDSIP_bb->curnode])->step != dual))
        {
            cnt = (int) ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i])[DDSIP_bb->firstvar]);
            ((DDSIP_node[DDSIP_bb->curnode])->first_sol[i])[DDSIP_bb->firstvar + 1] += 1.;
            if (((DDSIP_node[DDSIP_bb->curnode])->first_sol[i])[DDSIP_node[DDSIP_bb->nonode]->neoind] <= DDSIP_node[DDSIP_bb->nonode]->neoub)
            {
                if (DDSIP_param->outlev > 79)
                {
                    fprintf (DDSIP_bb->moreoutfile,"  DDSIP_node[%d]->first_sol[%d][%d]= %.15g <= %.15g = DDSIP_node[%d]->neoub\n",DDSIP_bb->curnode,i,DDSIP_node[DDSIP_bb->nonode]->neoind, ((DDSIP_node[DDSIP_bb->curnode])->first_sol[i])[DDSIP_node[DDSIP_bb->nonode]->neoind],DDSIP_node[DDSIP_bb->nonode]->neoub, DDSIP_bb->nonode);
                }
                //do check for inheritance conditions for all identical solutions
                for (j = i; cnt && j < DDSIP_param->scenarios; j++)
                {
                    if (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j]
                            && ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i]) == (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j])))
                    {
                        if ((((DDSIP_node[DDSIP_bb->curnode])->first_sol[j])[DDSIP_bb->firstvar + 1]) <= DDSIP_param->maxinherit ||
                                (DDSIP_node[DDSIP_bb->curnode]->mipstatus)[j] == CPXMIP_OPTIMAL ||
                                ((DDSIP_node[DDSIP_bb->curnode]->mipstatus)[j] == CPXMIP_OPTIMAL_TOL && (((DDSIP_node[DDSIP_bb->curnode])->first_sol[j])[DDSIP_bb->firstvar + 1]) < 2.*DDSIP_param->maxinherit))
                        {
                            (DDSIP_node[DDSIP_bb->nonode]->first_sol)[j] = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i];
                            (DDSIP_node[DDSIP_bb->nonode]->cursubsol)[j] = (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j];
                            (DDSIP_node[DDSIP_bb->nonode]->mipstatus)[j] = (DDSIP_node[DDSIP_bb->curnode]->mipstatus)[j];
                            (DDSIP_node[DDSIP_bb->nonode]->ref_scenobj)[j] = (DDSIP_node[DDSIP_bb->curnode]->ref_scenobj)[j];
			    if (j > i)
                                (DDSIP_node[DDSIP_bb->curnode]->first_sol)[j] = NULL;
                            if (DDSIP_param->outlev > 23)
                            {
                                fprintf (DDSIP_bb->moreoutfile,"  node %d inherited solution of scenario %d from node %d (%g identical scen. solutions), mipstatus %d, inh_level %g\n",
                                     DDSIP_bb->nonode, j + 1, DDSIP_bb->curnode, (DDSIP_node[DDSIP_bb->nonode]->first_sol)[j][DDSIP_bb->firstvar],
                                     (DDSIP_node[DDSIP_bb->nonode]->mipstatus)[j], ((DDSIP_node[DDSIP_bb->curnode])->first_sol[i])[DDSIP_bb->firstvar + 1]);
                            }
			    if (j > i)
                                (DDSIP_node[DDSIP_bb->curnode]->first_sol)[j] = NULL;
                            cnt--;
                        }
                        else
                        {
                            if (DDSIP_param->outlev > 21)
                            {
                                fprintf (DDSIP_bb->moreoutfile,"  node %d NOT inherited solution of scenario %d from node %d due to status (%g identical scen. solutions), mipstatus %d, inh_level %g\n",
                                     DDSIP_bb->nonode, j + 1, DDSIP_bb->curnode, (DDSIP_node[DDSIP_bb->curnode]->first_sol)[j][DDSIP_bb->firstvar],
                                     (DDSIP_node[DDSIP_bb->curnode]->mipstatus)[j], ((DDSIP_node[DDSIP_bb->curnode])->first_sol[i])[DDSIP_bb->firstvar + 1]);
                            }
                            ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[i])[DDSIP_bb->firstvar] -= 1.;
			    if (j > i)
                            {
                                (DDSIP_node[DDSIP_bb->curnode]->first_sol)[j] = NULL;
                            }
                            cnt--;
                        }
                    }
                }
                DDSIP_node[DDSIP_bb->nonode]->numInheritedSols += (int) (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i])[DDSIP_bb->firstvar];
                if(fabs((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i])[DDSIP_bb->firstvar]) < 0.2)
                    DDSIP_Free ((void **) &(((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i]));
                else
                    (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i] = NULL;
            }
            else if (((DDSIP_node[DDSIP_bb->curnode])->first_sol[i])[DDSIP_node[DDSIP_bb->nonode]->neoind] >= DDSIP_node[DDSIP_bb->nonode + 1]->neolb)
            {
                if (DDSIP_param->outlev > 79)
                {
                    fprintf (DDSIP_bb->moreoutfile,"  DDSIP_node[%d]->first_sol[%d][%d]= %.15g >= %.15g = DDSIP_node[%d]->neolb\n",DDSIP_bb->curnode,i,DDSIP_node[DDSIP_bb->nonode + 1]->neoind, ((DDSIP_node[DDSIP_bb->curnode])->first_sol[i])[DDSIP_node[DDSIP_bb->nonode + 1]->neoind],DDSIP_node[DDSIP_bb->nonode + 1]->neolb, DDSIP_bb->nonode + 1);
                }
                //do check for inheritance conditions for all identical solutions
                for (j = i; cnt && j < DDSIP_param->scenarios; j++)
                {
                    if (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j]
                            && ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i]) == (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j])))
                    {
                        if ((((DDSIP_node[DDSIP_bb->curnode])->first_sol[j])[DDSIP_bb->firstvar + 1]) <= DDSIP_param->maxinherit ||
                                (DDSIP_node[DDSIP_bb->curnode]->mipstatus)[j] == CPXMIP_OPTIMAL ||
                                ((DDSIP_node[DDSIP_bb->curnode]->mipstatus)[j] == CPXMIP_OPTIMAL_TOL && (((DDSIP_node[DDSIP_bb->curnode])->first_sol[j])[DDSIP_bb->firstvar + 1]) < 2.*DDSIP_param->maxinherit))
                        {
                            (DDSIP_node[DDSIP_bb->nonode + 1]->first_sol)[j] = (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i];
                            (DDSIP_node[DDSIP_bb->nonode + 1]->cursubsol)[j] = (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j];
                            (DDSIP_node[DDSIP_bb->nonode + 1]->mipstatus)[j] = (DDSIP_node[DDSIP_bb->curnode]->mipstatus)[j];
                            (DDSIP_node[DDSIP_bb->nonode + 1]->ref_scenobj)[j] = (DDSIP_node[DDSIP_bb->curnode]->ref_scenobj)[j];
                            if (DDSIP_param->outlev > 23)
                            {
                                fprintf (DDSIP_bb->moreoutfile,"  node %d inherited solution of scenario %d from node %d (%g identical scen. solutions), mipstatus %d, inh_level %g\n",
                                     DDSIP_bb->nonode + 1, j + 1, DDSIP_bb->curnode, (DDSIP_node[DDSIP_bb->curnode]->first_sol)[j][DDSIP_bb->firstvar],
                                     (DDSIP_node[DDSIP_bb->curnode]->mipstatus)[j], ((DDSIP_node[DDSIP_bb->curnode])->first_sol[i])[DDSIP_bb->firstvar + 1]);
                            }
			    if (j > i)
                                (DDSIP_node[DDSIP_bb->curnode]->first_sol)[j] = NULL;
                            cnt--;
                        }
			    else
                        {
                            if (DDSIP_param->outlev > 21)
                            {
                                fprintf (DDSIP_bb->moreoutfile,"  node %d NOT inherited solution of scenario %d from node %d due to status (%g identical scen. solutions), mipstatus %d, inh_level %g\n",
                                     DDSIP_bb->nonode + 1, j + 1, DDSIP_bb->curnode, (DDSIP_node[DDSIP_bb->curnode]->first_sol)[j][DDSIP_bb->firstvar],
                                     (DDSIP_node[DDSIP_bb->curnode]->mipstatus)[j], ((DDSIP_node[DDSIP_bb->curnode])->first_sol[i])[DDSIP_bb->firstvar + 1]);
                            }
                            ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[i])[DDSIP_bb->firstvar] -= 1.;
                            if (j > i)
                            {
                                (DDSIP_node[DDSIP_bb->curnode]->first_sol)[j] = NULL;
                            }
                            cnt--;
                        }
                    }
                }
                DDSIP_node[DDSIP_bb->nonode + 1]->numInheritedSols += (int) (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i])[DDSIP_bb->firstvar];
                if(fabs((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i])[DDSIP_bb->firstvar]) < 0.2)
                    DDSIP_Free ((void **) &(((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i]));
                else
                    (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i] = NULL;
            }
            else
            {
                if (DDSIP_param->outlev > 21)
                {
                    fprintf (DDSIP_bb->moreoutfile,"##scenario %d solution not passed on (%g identical scen. solutions)\n",
                             i+1, (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i])[DDSIP_bb->firstvar]);
                }
                if ((cnt = (int) ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i])[DDSIP_bb->firstvar] - 0.9)))
                    for (j = i + 1; cnt && j < DDSIP_param->scenarios; j++)
                    {
                        if (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j]
                                && ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i]) == (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j])))
                        {
                            if (DDSIP_param->outlev > 21)
                            {
                                fprintf (DDSIP_bb->moreoutfile,"##scenario %d solution not passed on (%g identical scen. solutions)\n",
                                         i+1, (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i])[DDSIP_bb->firstvar]);
                            }
                            ((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j] = NULL;
                            cnt--;
                        }
                    }
                DDSIP_Free ((void **) &(((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i]));
            }
        }
        else
        {
            if (DDSIP_param->outlev > 23)
            {
                if ((abs(DDSIP_param->riskmod) == 3 || abs(DDSIP_param->riskmod) == 4 || abs(DDSIP_param->riskmod) == 5))
                    fprintf (DDSIP_bb->moreoutfile,"  nodes %d and %d did not inherit solution of scenario %d from node %d due to risk model.\n",
                             DDSIP_bb->nonode, DDSIP_bb->nonode + 1, i + 1, DDSIP_bb->curnode);
                else if ((DDSIP_node[DDSIP_bb->curnode])->step == dual)
                    fprintf (DDSIP_bb->moreoutfile,"  nodes %d and %d did not inherit solution of scenario %d from node %d due to ConicBundle step for father.\n",
                             DDSIP_bb->nonode, DDSIP_bb->nonode + 1, i + 1, DDSIP_bb->curnode);
                else
                    fprintf (DDSIP_bb->moreoutfile,"  !!! UNEXPECTED: nodes %d and %d did not inherit solution of scenario %d from node %d.\n",
                             DDSIP_bb->nonode, DDSIP_bb->nonode + 1, i + 1, DDSIP_bb->curnode);
            }
            if ((cnt = (int) ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i])[DDSIP_bb->firstvar] - 0.9)))
                for (j = i + 1; cnt && j < DDSIP_param->scenarios; j++)
                {
                    if (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j]
                            && ((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i] == ((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j])
                    {
                        ((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j] = NULL;
                        cnt--;
                    }
                }
            DDSIP_Free ((void **) &(((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i]));
        }
    }
    if (DDSIP_param->outlev > 4)
    {
        fprintf (DDSIP_bb->moreoutfile,"##**node %d inherited %d solutions from node %d\n",
                 DDSIP_bb->nonode, (DDSIP_node[DDSIP_bb->nonode])->numInheritedSols, DDSIP_bb->curnode);
        fprintf (DDSIP_bb->moreoutfile,"##**node %d inherited %d solutions from node %d\n",
                 DDSIP_bb->nonode + 1, (DDSIP_node[DDSIP_bb->nonode + 1])->numInheritedSols, DDSIP_bb->curnode);
    }
    DDSIP_bb->cutAdded = 0;

    // Debugging information
    if (DDSIP_param->outlev > 2)
    {
        cnt = DDSIP_bb->firstindex[DDSIP_node[DDSIP_bb->nonode]->neoind];
        status = CPXgetcolname (DDSIP_env, DDSIP_lp, DDSIP_bb->name_buffer, DDSIP_bb->n_buffer, DDSIP_bb->n_buffer_len, &i, cnt, cnt);
        if (status)
            fprintf (stderr," Error when querying name of variable %d: %d\n",cnt,status);
        fprintf (DDSIP_bb->moreoutfile, "New nodes (no. of node, no. of variable branched on, lb, ub, name of var.)\n");
        fprintf (DDSIP_bb->moreoutfile, "%3d  %6d %18.16g %18.16g  %s\n", DDSIP_bb->nonode,
                 cnt, DDSIP_node[DDSIP_bb->nonode]->neolb, DDSIP_node[DDSIP_bb->nonode]->neoub, DDSIP_bb->n_buffer);
        fprintf (DDSIP_bb->moreoutfile, "%3d  %6d %18.16g %18.16g  %s\n", DDSIP_bb->nonode + 1,
                 DDSIP_bb->firstindex[DDSIP_node[DDSIP_bb->nonode + 1]->neoind], DDSIP_node[DDSIP_bb->nonode + 1]->neolb, DDSIP_node[DDSIP_bb->nonode + 1]->neoub, DDSIP_bb->n_buffer);
    }
    // Update front tree
    cnt = 0;
    while (!(DDSIP_bb->front[cnt] == DDSIP_bb->curnode) && (cnt < DDSIP_bb->nofront))
        cnt++;

    DDSIP_bb->front[cnt] = DDSIP_bb->nonode;
    DDSIP_bb->front[DDSIP_bb->nofront] = DDSIP_bb->nonode + 1;


    // Update current node numbers in tree and front tree
    DDSIP_bb->nofront++;
    DDSIP_bb->no_reduced_front++;
    DDSIP_bb->nonode += 2;

    // Free current node
    for (scen = 0; scen < DDSIP_param->scenarios; scen++)
    {
        if ((DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen])
        {
            fprintf (stderr, "xxx ERROR in InitNewNodes: (DDSIP_node[%d]->first_sol)[%d]=%p xxx\n", DDSIP_bb->curnode, scen,
                     (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen]);
            exit (1);
        }
    }
    DDSIP_FreeNode (DDSIP_bb->curnode);

    // Select node for branching
    // Select right or left node according to branchdir for solving
    if (DDSIP_param->branchdir == -1)
        DDSIP_bb->curnode = DDSIP_bb->front[cnt];
    else
        DDSIP_bb->curnode = DDSIP_bb->front[DDSIP_bb->nofront - 1];

    return 0;
}

//==========================================================================
// According to the node of the b&b tree, the function adds new bounds to problem.
int
DDSIP_SetBounds (void)
{
    int i, k, change;

    if (DDSIP_param->outlev > 2)
    {
        fprintf (DDSIP_bb->moreoutfile, "\n----------------------\n");
        fprintf (DDSIP_bb->moreoutfile, "Prepare solving of node %d...\n", DDSIP_bb->curnode);
    }
    // Initialize DDSIP_bb->curlb, DDSIP_bb->curub
    for (i = 0; i < DDSIP_bb->firstvar; i++)
    {
        DDSIP_bb->curlb[i] = DDSIP_bb->lborg[i];
        DDSIP_bb->curub[i] = DDSIP_bb->uborg[i];
    }

    // Get all constraints in the upper tree
    // First, the new constraint of the current node
    DDSIP_bb->curind[0] = DDSIP_node[DDSIP_bb->curnode]->neoind;
    DDSIP_bb->curlb[0] = DDSIP_node[DDSIP_bb->curnode]->neolb;
    DDSIP_bb->curub[0] = DDSIP_node[DDSIP_bb->curnode]->neoub;

    // Initialize
    i = DDSIP_node[DDSIP_bb->curnode]->father;
    DDSIP_bb->curbdcnt = 1;
    if (DDSIP_param->outlev)
        fprintf (DDSIP_bb->moreoutfile, "depth %2d:  %d", DDSIP_node[DDSIP_bb->curnode]->depth, DDSIP_bb->curnode);

    // While the father of the current node is greater than the root
    while (i > 0)
    {
        if (DDSIP_param->outlev)
            fprintf (DDSIP_bb->moreoutfile, " -> %d", i);
        k = 0;
        change = 0;

        // If the index has already occured, minima and maxima have to be calculated
        while (k < DDSIP_bb->curbdcnt)
        {
            if (DDSIP_bb->curind[k] == DDSIP_node[i]->neoind)
            {
                DDSIP_bb->curlb[k] = DDSIP_Dmax (DDSIP_node[i]->neolb, DDSIP_bb->curlb[k]);
                DDSIP_bb->curub[k] = DDSIP_Dmin (DDSIP_node[i]->neoub, DDSIP_bb->curub[k]);
                // Actually, the next line should be superfluous
                DDSIP_bb->curlb[k] = DDSIP_Dmin (DDSIP_bb->curlb[k], DDSIP_bb->curub[k]);
                change = 1;
            }
            k++;
        }

        // If the index did not occur yet, a new one is introduced
        if (!(change) && DDSIP_bb->curbdcnt < DDSIP_bb->firstvar)
        {
            DDSIP_bb->curind[DDSIP_bb->curbdcnt] = DDSIP_node[i]->neoind;
            DDSIP_bb->curlb[DDSIP_bb->curbdcnt] = DDSIP_node[i]->neolb;
            DDSIP_bb->curub[DDSIP_bb->curbdcnt] = DDSIP_node[i]->neoub;
            DDSIP_bb->curbdcnt++;
        }

        i = DDSIP_node[i]->father;
    }

    // the depth of the current node may become the depth of the tree
    DDSIP_bb->depth = DDSIP_Imax (DDSIP_node[DDSIP_bb->curnode]->depth, DDSIP_bb->depth);

    // A test on some errors
    for (i = 0; i < DDSIP_bb->curbdcnt; i++)
        if (DDSIP_bb->curlb[i] > DDSIP_bb->curub[i])
            return 115;

    if (DDSIP_bb->curbdcnt > DDSIP_bb->firstvar)
        return 117;

    return 0;
}

//==========================================================================
// Function to set the leaf indicator for the current node
// This indicator is 1 if the current node is a front node of the complete enumeration tree
int
DDSIP_Leaf (void)
{
    int i;
    int fs = DDSIP_bb->firstvar;

    // Branch on additional first-stage variables in risk models WC and TVaR?
    if ((abs (DDSIP_param->riskmod) == 4 || abs (DDSIP_param->riskmod) == 5) && !DDSIP_param->brancheta)
        fs--;

    // Is the selected node a front node of the complete enumeration tree
    if (DDSIP_bb->curbdcnt == fs)
    {
        DDSIP_node[DDSIP_bb->curnode]->leaf = 1;


        for (i = 0; i < DDSIP_bb->curbdcnt; i++)
        {


            // For continuous variables the difference has to be at least DDSIP_param->brancheps
            // For binary or integer variables the difference has to be at least 1
            if (   (DDSIP_bb->firsttype[DDSIP_bb->curind[i]] == 'B' && !(DDSIP_bb->curlb[i] == DDSIP_bb->curub[i]))
                    || (DDSIP_bb->firsttype[DDSIP_bb->curind[i]] == 'I' && DDSIP_bb->curlb[i] < DDSIP_bb->curub[i])
                    || (DDSIP_bb->firsttype[DDSIP_bb->curind[i]] == 'N' && DDSIP_bb->curlb[i] < DDSIP_bb->curub[i])
                    || (DDSIP_bb->firsttype[DDSIP_bb->curind[i]] == 'C' && (DDSIP_bb->curub[i] - DDSIP_bb->curlb[i])/(fabs(DDSIP_bb->curub[i])+1.e-16) > DDSIP_param->brancheps)
                    || (DDSIP_bb->firsttype[DDSIP_bb->curind[i]] == 'S' && (DDSIP_bb->curub[i] - DDSIP_bb->curlb[i])/(fabs(DDSIP_bb->curub[i])+1.e-16) > DDSIP_param->brancheps))
            {
                DDSIP_node[DDSIP_bb->curnode]->leaf = 0;
                break;
            }
        }

        if (DDSIP_node[DDSIP_bb->curnode]->leaf)
            DDSIP_bb->no_reduced_front--;
    }

    return 0;
}

//==========================================================================
// Function gets new bound and fathoms nodes with higher nodevalue than
// best solution
// Bound returns an integer indicating whether maximal dispersion is less than
// DDSIP_nulldisp (1), the whole branching tree was backtracked (2), the bestbound
// in the front tree equals DDSIP_infinity (3) or 0 for continue.

int
DDSIP_Bound (void)
{
    int i, j = 1, k, scen, cnt, status = 1;
    static int callcnt = 0;
    double * front_node_bound, rgap, factor, worstBound;
    double bestAmongTheLast;
    callcnt++;
    factor = (DDSIP_bb->bestvalue < 0.)? 1.-3.e-15 :  1.+3.e-15;

    if (DDSIP_param->outlev > 2)
    {
        fprintf (DDSIP_bb->moreoutfile, "\n----------------------\n");
        fprintf (DDSIP_bb->moreoutfile, "Bounding...\n");
    }
    // Fathom nodes in front tree
    for (i = DDSIP_bb->nofront - 1; i >= 0; i--)
    {
        // if a cut was added in the current node, mark all nodes in the front of the tree (to check cuts before solution passing)
        if (DDSIP_bb->cutAdded)
            DDSIP_node[DDSIP_bb->front[i]]->cutAdded = 1;
        if ((!(DDSIP_bb->found_optimal_node) && (DDSIP_node[DDSIP_bb->front[i]]->bound > DDSIP_bb->bestvalue*factor + DDSIP_bb->correct_bounding)
            ) ||
                ( (DDSIP_bb->found_optimal_node) && (DDSIP_bb->found_optimal_node != DDSIP_bb->front[i] || DDSIP_node[DDSIP_bb->front[i]]->bound > DDSIP_bb->bestvalue + 2.*DDSIP_bb->correct_bounding) &&
                  ((DDSIP_node[DDSIP_bb->front[i]]->violations && (DDSIP_node[DDSIP_bb->front[i]]->bound > DDSIP_bb->bound_optimal_node + DDSIP_bb->correct_bounding)) ||
                   (!(DDSIP_node[DDSIP_bb->front[i]]->violations) && (DDSIP_node[DDSIP_bb->front[i]]->bound > DDSIP_bb->bestvalue))
                  )
                ) ||
                DDSIP_Equal (DDSIP_node[DDSIP_bb->front[i]]->bound, DDSIP_infty)
           )
        {
            // debug info
            if (DDSIP_param->outlev > 19)
            {
#ifdef DEBUG
                printf (" Bounding: delete node %d, bound: %.16g, bestvalue: %.16g, bestvalue*factor= %.16g, previous bestbound: %.16g\n", DDSIP_bb->front[i],
                        DDSIP_node[DDSIP_bb->front[i]]->bound, DDSIP_bb->bestvalue, DDSIP_bb->bestvalue*factor, DDSIP_bb->bestbound);
                printf ("                  node %d, bound - bestvalue = %.16g, bound - (bestvalue*factor(%.17g) + DDSIP_bb->correct_bounding(%g)) = %.16g, accuracy= %g\n", DDSIP_bb->front[i],
                        DDSIP_node[DDSIP_bb->front[i]]->bound- DDSIP_bb->bestvalue, DDSIP_bb->bestvalue*factor,
                        DDSIP_node[DDSIP_bb->front[i]]->bound- (DDSIP_bb->bestvalue*factor), DDSIP_bb->correct_bounding, DDSIP_param->accuracy);
#endif
                fprintf (DDSIP_bb->moreoutfile, " Bounding: delete node %d, bound: %.16g, bestvalue: %.16g, previous bestbound: %.16g\n", DDSIP_bb->front[i],
                         DDSIP_node[DDSIP_bb->front[i]]->bound, DDSIP_bb->bestvalue, DDSIP_bb->bestbound);
                fprintf (DDSIP_bb->moreoutfile, "                  node %d, bound - bestvalue = %.16g, bound - (bestvalue*factor) = %.16g, accuracy= %g\n", DDSIP_bb->front[i],
                         DDSIP_node[DDSIP_bb->front[i]]->bound- DDSIP_bb->bestvalue,
                         DDSIP_node[DDSIP_bb->front[i]]->bound- (DDSIP_bb->bestvalue*factor), DDSIP_param->accuracy);
            }
            // mark ist as a leaf
            DDSIP_node[DDSIP_bb->front[i]]->leaf = 1;
            // if a node previously seeming to be kept is now deleted, remove found_optimal_node info
            if (DDSIP_bb->found_optimal_node == DDSIP_bb->front[i])
                DDSIP_bb->found_optimal_node = 0;
            // Free the node's allocated arrays
            for (scen = 0; scen < DDSIP_param->scenarios; scen++)
            {
                if (((DDSIP_node[DDSIP_bb->front[i]])->first_sol)[scen])
                {
                    if ((cnt = (int) ((((DDSIP_node[DDSIP_bb->front[i]])->first_sol)[scen])[DDSIP_bb->firstvar] - 0.9)))
                    {
                        for (j = scen + 1; cnt && j < DDSIP_param->scenarios; j++)
                        {
                            if (((DDSIP_node[DDSIP_bb->front[i]])->first_sol)[j]
                                    && (((DDSIP_node[DDSIP_bb->front[i]])->first_sol)[scen] == ((DDSIP_node[DDSIP_bb->front[i]])->first_sol)[j]))
                            {
                                ((DDSIP_node[DDSIP_bb->front[i]])->first_sol)[j] = NULL;
                                cnt--;
                            }
                        }
                    }
                }
                DDSIP_Free ((void **) & ((DDSIP_node[DDSIP_bb->front[i]]->first_sol)[scen]));
            }
            DDSIP_FreeNode (DDSIP_bb->front[i]);

            // Delete node from the front
            DDSIP_bb->front[i] = DDSIP_bb->front[DDSIP_bb->nofront - 1];
            DDSIP_bb->nofront--;
        }
    }

    j = 1;
    // If front is not empty
    if (DDSIP_bb->nofront)
    {

        DDSIP_bb->no_reduced_front = DDSIP_bb->nofront;
        for (i = 0; i < DDSIP_bb->nofront; i++)
            if (DDSIP_node[DDSIP_bb->front[i]]->leaf)
                DDSIP_bb->no_reduced_front--;

        // sort front nodes ascending wrt. bound
        DDSIP_Free ((void **) &DDSIP_bb->front_nodes_sorted);
        DDSIP_bb->front_nodes_sorted = (int *) DDSIP_Alloc(sizeof(int), DDSIP_bb->nofront, "DDSIP_bb->front_nodes_sorted(Bound)");
        front_node_bound = (double *) DDSIP_Alloc(sizeof(double), DDSIP_bb->nonode, "front_node_bound(Bound)");
        for  (i = 0; i < DDSIP_bb->nofront; i++)
        {
            DDSIP_bb->front_nodes_sorted[i] = DDSIP_bb->front[i];
            front_node_bound[DDSIP_bb->front_nodes_sorted[i]] =  DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->bound;
        }
        DDSIP_qsort_ins_A (front_node_bound, DDSIP_bb->front_nodes_sorted, 0, DDSIP_bb->nofront-1);
        // Update bestbound = lowest bound within front nodes
        DDSIP_bb->bestbound = DDSIP_node[DDSIP_bb->front_nodes_sorted[0]]->bound;
        worstBound = DDSIP_node[DDSIP_bb->front_nodes_sorted[DDSIP_bb->nofront-1]]->bound;
        // Now sort leaf nodes to the end
        for  (i = 0; i < DDSIP_bb->nofront; i++)
        {
            DDSIP_bb->front_nodes_sorted[i] = DDSIP_bb->front[i];
            front_node_bound[DDSIP_bb->front_nodes_sorted[i]] =  (DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->leaf) ? DDSIP_infty : DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->bound;
        }
        DDSIP_qsort_ins_A (front_node_bound, DDSIP_bb->front_nodes_sorted, 0, DDSIP_bb->nofront-1);
        cnt = 1;

        // sort the least bound nodes according to the violations and depth
        for  (i = 0; i < DDSIP_bb->nofront; i++)
        {
            front_node_bound[DDSIP_bb->front_nodes_sorted[i]] = !(DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->solved) ? -DDSIP_infty : 2*DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->violations + DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->dispnorm - DDSIP_node[DDSIP_bb->front[i]]->depth;
            if ((DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->bound - DDSIP_bb->bestbound) / (fabs(DDSIP_bb->bestbound) + DDSIP_param->accuracy) > 1.e-12)
                break;
        }
        DDSIP_qsort_ins_A (front_node_bound, DDSIP_bb->front_nodes_sorted, 0, i-1);

        if (DDSIP_param->boundstrat > 4 && DDSIP_param->boundstrat < 10 && DDSIP_bb->bestvalue < DDSIP_infty)
        {
            // feasible point found - switch to other bounding strategy
            DDSIP_param->boundstrat -= 5;
        }
        if (DDSIP_param->boundstrat > 4 && !(DDSIP_param->cb < 0 && ((DDSIP_bb->cutoff > 5) && (DDSIP_bb->no_reduced_front < 51) && (DDSIP_bb->noiter % -DDSIP_param->cb) < DDSIP_param->cbContinuous)))
        {
            // depth first until a feasible point is found for boustrat 5..9
            // choose among the most recent nodes according to best bound
            if (DDSIP_bb->no_reduced_front > 6)
            {
                int depth_first_nodes = 6;
                if ((DDSIP_bb->Dive || callcnt%DDSIP_param->bestboundfreq < DDSIP_param->bestboundfreq-2) && (DDSIP_bb->curnode%500 < 439))
                {
                    if (DDSIP_param->outlev > 5)
                    {
                        if (DDSIP_bb->Dive)
                            fprintf (DDSIP_bb->moreoutfile, " - selection of next node: dive\n");
                        else
                            fprintf (DDSIP_bb->moreoutfile, " - selection of next node: depth first\n");
                        if (DDSIP_bb->bestBound)
                        {
                            fprintf (DDSIP_bb->moreoutfile, " DDSIP_bb->bestBound = %d -> select best bound node among all.\n", DDSIP_bb->bestBound);
                        }
                    }
                    // use depth first
                    // sort front nodes wrt. their number
                    for (i = 0; i < DDSIP_bb->nofront; i++)
                    {
                        front_node_bound[DDSIP_bb->front_nodes_sorted[i]] =  (DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->leaf) ? -1 : DDSIP_bb->front_nodes_sorted[i];
                    }
                    DDSIP_qsort_ins_D (front_node_bound, DDSIP_bb->front_nodes_sorted, 0, DDSIP_bb->nofront-1);
                    if (DDSIP_param->boundstrat < 6 || DDSIP_param->boundstrat == 10)
                    {
                        double threshold = -1e6;
                        DDSIP_bb->Dive = 0;
                        bestAmongTheLast = DDSIP_infty;
                        // branch the one with best bound among the last generated nodes
                        if (DDSIP_bb->bestBound)
                        {
                            //depth_first_nodes = DDSIP_Imin(20, DDSIP_bb->nofront);
                            depth_first_nodes = DDSIP_bb->nofront;
                        }
                        else
                        {
                            DDSIP_bb->Dive = (DDSIP_bb->curnode > 40 && DDSIP_bb->curnode%500 < 61) ? 1: 0;
                            if (DDSIP_bb->Dive)
                                depth_first_nodes = 2;
                            else
                                depth_first_nodes = 4;
                        }
                        for (i = 0; i < depth_first_nodes; i++)
                        {
                            front_node_bound[DDSIP_bb->front_nodes_sorted[i]] =  (DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->leaf) ? DDSIP_infty : DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->bound;
                            bestAmongTheLast = DDSIP_Dmin(front_node_bound[DDSIP_bb->front_nodes_sorted[i]], bestAmongTheLast);
                        }
                        DDSIP_qsort_ins_A (front_node_bound, DDSIP_bb->front_nodes_sorted, 0, depth_first_nodes-1);
                        threshold = DDSIP_bb->Dive ? 0.20*DDSIP_bb->bestbound + 0.80*worstBound:(1. - DDSIP_param->btTolerance)*DDSIP_bb->bestbound + DDSIP_param->btTolerance*worstBound;
                        for (i = depth_first_nodes; i < DDSIP_bb->nofront; i++)
                        {
                            front_node_bound[DDSIP_bb->front_nodes_sorted[i]] =  (DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->leaf) ? DDSIP_infty : DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->bound;
                        }
                        DDSIP_qsort_ins_A (front_node_bound, DDSIP_bb->front_nodes_sorted, depth_first_nodes, DDSIP_bb->nofront-1);
                        if (DDSIP_param->outlev > 5 && !DDSIP_bb->bestBound)
                        {
                            fprintf (DDSIP_bb->moreoutfile, "                                least bound among all:  %-16.12g, greatest bound: %-16.12g\n", DDSIP_bb->bestbound, worstBound);
                            fprintf (DDSIP_bb->moreoutfile, "                                least bound among last: %-16.12g, threshold:      %-16.12g", bestAmongTheLast, threshold);
                        }
                        if (DDSIP_param->boundstrat == 10)
                        {
                            // use best bound of all nodes if the best bound within the last nodes is too big and more often if
                            // we have at least one feasible point
                            if ((!DDSIP_bb->Dive &&
                                    (bestAmongTheLast > threshold || DDSIP_bb->bestBound ||
                                     ((fabs (DDSIP_bb->bestvalue) < DDSIP_infty) && !(DDSIP_bb->noiter % 50)))) ||
                                    (DDSIP_bb->Dive && bestAmongTheLast > threshold)
                               )
                            {
                                for  (i = 0; i < DDSIP_bb->nofront; i++)
                                {
                                    front_node_bound[DDSIP_bb->front_nodes_sorted[i]] =   (DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->leaf) ? DDSIP_infty : DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->bound;
                                }
                                DDSIP_qsort_ins_A (front_node_bound, DDSIP_bb->front_nodes_sorted, 0, DDSIP_bb->nofront-1);
                                // sort the least bound nodes according to the violations and depth
                                for  (i = 0; i < DDSIP_bb->nofront; i++)
                                {
                                    if ((DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->bound - DDSIP_bb->bestbound) / (fabs(DDSIP_bb->bestbound) + DDSIP_param->accuracy) > 1.e-12)
                                        break;
                                    front_node_bound[DDSIP_bb->front_nodes_sorted[i]] = !(DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->solved) ? -DDSIP_infty : 2*DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->violations + DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->dispnorm - DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->depth;
                                }
                                DDSIP_qsort_ins_A (front_node_bound, DDSIP_bb->front_nodes_sorted, 0, i-1);
                                if (DDSIP_bb->bestBound > 2)
                                {
                                    DDSIP_bb->bestBound = 0;
                                }
                                else
                                    DDSIP_bb->bestBound++;
                                if (DDSIP_param->outlev > 5 && DDSIP_bb->bestBound == 1)
                                {
                                    if ((fabs (DDSIP_bb->bestvalue) < DDSIP_infty) && !(DDSIP_bb->noiter % 50))
                                        fprintf (DDSIP_bb->moreoutfile, " freq: force bestbound ");
                                    fprintf (DDSIP_bb->moreoutfile, " -> select best bound node among all.");
                                }
                            }
                            if (!DDSIP_bb->Dive)
                            {
                                for  (i = 0; i < DDSIP_bb->nofront; i++)
                                {
                                    if ((DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->bound - DDSIP_bb->bestbound) / (fabs(DDSIP_bb->bestbound) + DDSIP_param->accuracy) > 1.e-12)
                                        break;
                                    front_node_bound[DDSIP_bb->front_nodes_sorted[i]] = !(DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->solved) ?
                                            -DDSIP_infty : 2*DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->violations - DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->dispnorm - DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->depth;
                                }
                            }
                            else
                            {
                                for (i = 0; i < depth_first_nodes; i++)
                                {
                                    if (DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->bound > threshold)
                                        front_node_bound[DDSIP_bb->front_nodes_sorted[i]] =  DDSIP_infty;
                                    else
                                        front_node_bound[DDSIP_bb->front_nodes_sorted[i]] =  (DDSIP_node[DDSIP_bb->front[i]]->leaf) ? DDSIP_infty : DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->bound + 10.*DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->dispnorm + DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->violations;
                                }
                            }
                            DDSIP_qsort_ins_A (front_node_bound, DDSIP_bb->front_nodes_sorted, 0, i-1);
                            // reset the sorting criterion to bound
                            //for  (k = 0; k < DDSIP_Imin(i, DDSIP_bb->nofront); k++)
                            //{
                            //    front_node_bound[DDSIP_bb->front_nodes_sorted[k]] =  DDSIP_node[DDSIP_bb->front_nodes_sorted[k]]->bound;
                            //}
                        }
                        if (DDSIP_param->outlev > 5)
                            fprintf (DDSIP_bb->moreoutfile, "\n");
                    }
                    else if (DDSIP_param->boundstrat < 8)
                    {
                        // branch the one with least dispnorm among the last generated nodes
                        if (DDSIP_param->outlev > 5)
                            fprintf (DDSIP_bb->moreoutfile, "                                least dispnorm\n");
                        for (i = 0; i < depth_first_nodes; i++)
                        {
                            front_node_bound[DDSIP_bb->front_nodes_sorted[i]] =  (DDSIP_node[DDSIP_bb->front[i]]->leaf) ? DDSIP_infty : DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->dispnorm;
                        }
                        DDSIP_qsort_ins_A (front_node_bound, DDSIP_bb->front_nodes_sorted, 0, depth_first_nodes - 1);
                        for (i = depth_first_nodes; i < DDSIP_bb->nofront; i++)
                        {
                            front_node_bound[DDSIP_bb->front_nodes_sorted[i]] =  (DDSIP_node[DDSIP_bb->front[i]]->leaf) ? DDSIP_infty : DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->bound;
                        }
                        DDSIP_qsort_ins_A (front_node_bound, DDSIP_bb->front_nodes_sorted, depth_first_nodes, DDSIP_bb->nofront-1);
                    }
                    else
                    {
                        // branch the one with least violations among the last generated nodes
                        if (DDSIP_param->outlev > 5)
                            fprintf (DDSIP_bb->moreoutfile, "                                least violations\n");
                        for (i = 0; i < depth_first_nodes; i++)
                        {
                            front_node_bound[DDSIP_bb->front_nodes_sorted[i]] =  (DDSIP_node[DDSIP_bb->front[i]]->leaf) ? DDSIP_infty : DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->violations;
                        }
                        DDSIP_qsort_ins_A (front_node_bound, DDSIP_bb->front_nodes_sorted, 0, depth_first_nodes -1);
                        for (i = depth_first_nodes; i < DDSIP_bb->nofront; i++)
                        {
                            front_node_bound[DDSIP_bb->front_nodes_sorted[i]] =  (DDSIP_node[DDSIP_bb->front[i]]->leaf) ? DDSIP_infty : DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->bound;
                        }
                        DDSIP_qsort_ins_A (front_node_bound, DDSIP_bb->front_nodes_sorted, depth_first_nodes, DDSIP_bb->nofront-1);
                    }
                }
                else //if (callcnt%DDSIP_param->bestboundfreq > DDSIP_param->bestboundfreq-3)
                {
                    // choose a best bound node
                    if (DDSIP_param->outlev > 5)
                        fprintf (DDSIP_bb->moreoutfile, " - selection of next node: best bound\n");
                    // reset the sorting criterion to bound
                    for  (k = 0; k < DDSIP_bb->nofront; k++)
                    {
                        front_node_bound[DDSIP_bb->front_nodes_sorted[k]] =  (DDSIP_node[DDSIP_bb->front[k]]->leaf) ? DDSIP_infty : DDSIP_node[DDSIP_bb->front_nodes_sorted[k]]->bound;
                    }
                    DDSIP_qsort_ins_A (front_node_bound, DDSIP_bb->front_nodes_sorted, 0, DDSIP_bb->nofront-1);
                    for  (i = 0; i < DDSIP_bb->nofront; i++)
                    {
                        if (DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->bound > DDSIP_bb->bestbound + fabs(DDSIP_bb->bestbound)*1.e-15 || DDSIP_node[DDSIP_bb->front[i]]->leaf)
                            break;
                        front_node_bound[DDSIP_bb->front_nodes_sorted[i]] = !(DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->solved) ?
                                -DDSIP_infty : 2*DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->violations + DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->dispnorm - DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->depth;
                    }
                    DDSIP_qsort_ins_A (front_node_bound, DDSIP_bb->front_nodes_sorted, 0, i-1);
                }
            }
        }
        else
        {
            double threshold;
            if (DDSIP_param->boundstrat)
            {
                if ((callcnt%DDSIP_param->period)>DDSIP_param->rgapsmall)
                {
                    // DDSIP_param->rgapsmall  strategy - but half of the time branch the least-bound nodes
                    if (!DDSIP_Equal (fabs (DDSIP_bb->bestvalue), 0.0))
                    {
                        threshold = fabs ((DDSIP_bb->bestvalue - DDSIP_bb->bestbound) / DDSIP_bb->bestvalue);
                        if (threshold < 1e-4)
                            rgap = 0.1*threshold;
                        else
                            rgap = 0.02*threshold;
                    }
                    else
                    {
                        threshold = fabs ((DDSIP_bb->bestvalue - DDSIP_bb->bestbound) / (fabs(DDSIP_bb->bestvalue) + DDSIP_param->accuracy));
                        if (threshold < 1e-4)
                            rgap = 0.1*threshold;
                        else
                            rgap = 0.02*threshold;
                    }
                }
                else
                {
                    if ((callcnt%DDSIP_param->period)<DDSIP_param->rgapsmall/2)
                        rgap = 1.e-10;
                    else
                        rgap = 1.e-07;
                    rgap = DDSIP_Dmin (rgap, 0.1*fabs ((DDSIP_bb->bestvalue - DDSIP_bb->bestbound) / (fabs(DDSIP_bb->bestvalue) + DDSIP_param->accuracy)));
                }
            }
            else
            {
                // best bound strategy
                //rgap = 0.01*DDSIP_param->relgap;
                rgap = 1.e-14;
            }
            threshold = DDSIP_bb->bestbound + fabs(DDSIP_bb->bestbound)*rgap;
            for (i = 1; i < DDSIP_bb->nofront; i++)
            {
                if (DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->bound > threshold)
                    break;
                cnt++;
            }

            if (DDSIP_param->outlev > 5)
                fprintf (DDSIP_bb->moreoutfile, " - selection of next node: rgap=%g, threshold value = %18.13g, cnt = %d  (%d%%%d)=%d (comp. with %d)\n", rgap, threshold,cnt,callcnt,DDSIP_param->period,(callcnt%DDSIP_param->period), DDSIP_param->rgapsmall);

            if (cnt > 1)
            {
                if (DDSIP_param->boundstrat && DDSIP_param->boundstrat < 3)
                {
                    // sort front nodes with least bound wrt. dispnorm (in order to branch nodes with bigger dispnorm first)
                    for (i = 0; i < cnt; i++)
                    {
                        front_node_bound[DDSIP_bb->front_nodes_sorted[i]] =  DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->dispnorm;
                    }
                    if (DDSIP_param->boundstrat && (callcnt%DDSIP_param->period) < DDSIP_param->rgapsmall+(DDSIP_param->period - DDSIP_param->rgapsmall)*0.3)
                    {
                        // branch nodes with smaller dispersion/violations norm first (in the hope to get a better feasible heuristics sooner)
                        DDSIP_qsort_ins_A (front_node_bound, DDSIP_bb->front_nodes_sorted, 0, cnt-1);
                        rgap =  DDSIP_node[DDSIP_bb->front_nodes_sorted[0]]->dispnorm;
                        for (i = 1; i < cnt; i++)
                        {
                            if (DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->dispnorm > rgap + 1.e-12)
                                break;
                        }
                        j = DDSIP_Imin(i,cnt);
                        // DEBUGOUT
                        if (DDSIP_param->outlev > 5)
                            fprintf (DDSIP_bb->moreoutfile, " - small dispersion norm, j= %d\n",j);
                        // DEBUGOUT
                    }
                    else
                    {
                        // branch nodes with greater dispersion/violations norm first (in the hope to increase the lower bound sooner)
                        DDSIP_qsort_ins_D (front_node_bound, DDSIP_bb->front_nodes_sorted, 0, cnt-1);
                        rgap =  DDSIP_node[DDSIP_bb->front_nodes_sorted[0]]->dispnorm;
                        for (i = 1; i < cnt; i++)
                        {
                            if (DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->dispnorm < rgap - 1.e-12)
                                break;
                        }
                        j = DDSIP_Imin(i,cnt);
                        // DEBUGOUT
                        if (DDSIP_param->outlev > 5)
                            fprintf (DDSIP_bb->moreoutfile, " - big   dispersion norm, j= %d\n",j);
                        // DEBUGOUT
                    }
                }
                else if (!DDSIP_param->boundstrat || DDSIP_param->boundstrat == 3 || DDSIP_param->boundstrat == 4)
                {
                    // sort front nodes with least bound wrt. violations
                    for (i = 0; i < cnt; i++)
                    {
                        front_node_bound[DDSIP_bb->front_nodes_sorted[i]] =  DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->violations;
                    }
                    // branch nodes with fewer violations first (in the hope to get a better feasible heuristics sooner)
                    DDSIP_qsort_ins_A (front_node_bound, DDSIP_bb->front_nodes_sorted, 0, cnt-1);
                    rgap =  DDSIP_node[DDSIP_bb->front_nodes_sorted[0]]->violations;
                    for (i = 1; i < cnt; i++)
                    {
                        if (DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->violations > rgap + 1)
                            break;
                    }
                    j = DDSIP_Imin(i,cnt);
                    // DEBUGOUT
                    if (DDSIP_param->outlev > 5)
                        fprintf (DDSIP_bb->moreoutfile, " - few   violations, j= %d\n",j);
                    // DEBUGOUT
                }
            }

            if (j > 1)
            {
                for (i = 0; i < j; i++)
                {
                    front_node_bound[DDSIP_bb->front_nodes_sorted[i]] =  DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->bound;
                }
                DDSIP_qsort_ins_A (front_node_bound, DDSIP_bb->front_nodes_sorted, 0, j-1);
            }
        }

        DDSIP_Free ((void**) &front_node_bound);
        // Debugging information
        if (DDSIP_param->outlev > 4)
        {
            fprintf (DDSIP_bb->moreoutfile,
                     "No of front nodes: %d (including %d leaves)     found_optimal_node: %d, bestbound: %18.12g\n", DDSIP_bb->nofront, DDSIP_bb->nofront - DDSIP_bb->no_reduced_front, DDSIP_bb->found_optimal_node, DDSIP_bb->bestbound);
            fprintf (DDSIP_bb->moreoutfile, "     No.   bound             violations dispnorm  branchvar  lower bound   upper        range       depth isleaf solved cutAdded\n");
            j = (DDSIP_param->outlev > 21 || !(DDSIP_bb->curnode % 200)) ? DDSIP_bb->nofront : DDSIP_Imin(DDSIP_bb->nofront,25);
            for (i = 0; i < j; i++)
            {
                if (DDSIP_Equal (DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->dispnorm, DDSIP_infty))
                {
                    if (DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->solved)
                        fprintf (DDSIP_bb->moreoutfile,
                                 "%-4d %-4d %-20.16g  %-9d  inf         %-6d %-12.6g %-12.6g %-11.5g %5d   %-6d %-6d %d\n",
                                 i+1, DDSIP_bb->front_nodes_sorted[i], DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->bound,
                                 DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->violations,
                                 DDSIP_bb->firstindex[DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->neoind],
                                 DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->neolb, DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->neoub,
                                 DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->neoub-DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->neolb,
                                 DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->depth,
                                 DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->leaf, DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->solved,
                                 DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->cutAdded);
                    else
                        fprintf (DDSIP_bb->moreoutfile,
                                 "%-4d %-4d %-20.16g             inf         %-6d %-12.6g %-12.6g %-11.5g %5d   %-6d %-6d\n",
                                 i+1, DDSIP_bb->front_nodes_sorted[i], DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->bound,
                                 DDSIP_bb->firstindex[DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->neoind],
                                 DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->neolb, DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->neoub,
                                 DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->neoub-DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->neolb,
                                 DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->depth,
                                 DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->leaf, DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->solved);
                }
                else
                    fprintf (DDSIP_bb->moreoutfile,
                             "%-4d %-4d %-20.16g  %-9d %-11.5g  %-6d %-12.6g %-12.6g %-11.5g %5d   %-6d %-6d %d\n",
                             i+1, DDSIP_bb->front_nodes_sorted[i], DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->bound,
                             DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->violations,
                             DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->dispnorm,
                             DDSIP_bb->firstindex[DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->neoind],
                             DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->neolb, DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->neoub,
                             DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->neoub-DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->neolb,
                             DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->depth,
                             DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->leaf, DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->solved,
                             DDSIP_node[DDSIP_bb->front_nodes_sorted[i]]->cutAdded);
            }
        }
        // No nodes left ?
        if (!DDSIP_bb->no_reduced_front)
        {
            if (DDSIP_param->outlev)
                fprintf (DDSIP_bb->moreoutfile, "-- No nodes left in the front of the tree. --\n");
            return 2;
        }

        // If there is one node in the front tree with a dispersion norm
        // larger then nulldisp continue branching
        i = 0;
        while (i < DDSIP_bb->nofront && status)
        {
            /*      if (!(DDSIP_node[DDSIP_bb->front[i]]->dispnorm <= DDSIP_param->nulldisp) && !DDSIP_node[DDSIP_bb->front[i]]->leaf)  */
            /* 13.9.2005: Knoten mit kleiner dispnorm werden nicht abgeschnitten, aber auch nicht ausgewaehlt, wenn
             *            dispnorm > nulldisp, aber < accuracy
             *            Loesung: schneide auch solche Knoten weg.                                                     */
            if (!(DDSIP_node[DDSIP_bb->front[i]]->dispnorm <= DDSIP_param->nulldisp + DDSIP_param->accuracy) && !DDSIP_node[DDSIP_bb->front[i]]->leaf)
                status = 0;
            i++;
        }
        return status;
    }
    // If front is empty
    else
    {
        DDSIP_bb->bestbound = DDSIP_Dmax (DDSIP_bb->bestbound, DDSIP_node[DDSIP_bb->curnode]->bound);
        if (DDSIP_param->outlev)
            fprintf (DDSIP_bb->moreoutfile, "-- No nodes left in the front of the tree. --\n");
        return 2;
    }
}

//==========================================================================
// Function generates new nodes in b&b-tree and chooses a node to be solved
// If curnode has been solved --> branch
// Old version: if (!DDSIP_Equal(DDSIP_node[DDSIP_bb->curnode]->dispnorm,DDSIP_infty)) --> invalid if all scenario problems are infeasible
// Branchval is used in the new constraints on some component
int
DDSIP_Branch (void)
{
    int status;

    // Retrieve node to be branched and/or solved
    if ((status = DDSIP_GetCurNode ()))
        return status;

    if (DDSIP_node[DDSIP_bb->curnode]->solved)
    {
        if ((status = DDSIP_InitNewNodes ()))
            return status;
    }

    status = DDSIP_SetBounds ();
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to add new bounds (Branch).\n");
        return status;
    }

    status = DDSIP_Leaf ();
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to perform Leaf() (Branch).\n");
        return status;
    }

    return status;
}
