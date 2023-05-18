/*  Authors:            Andreas M"arkert, Ralf Gollmer
    Copyright to:      University of Duisburg-Essen
    Language:          C

    Description:
    This file contains procedures to manage the memory use of the program.

    License:
    This file is part of DDSIP.

    DDSIP is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the DDSIP_Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    DDSIP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with DDSIP; if not, write to the DDSIP_Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

#include <DDSIP.h>
#include <DDSIPconst.h>

//==========================================================================
void
DDSIP_Free (void **ptr)
{
    if (!(*ptr == NULL))
        free (*ptr);
    *ptr = NULL;
}

//==========================================================================
void *
DDSIP_Alloc (int elsize, int nelem, const char *name)
{
    int status = 0;
    void *ptr = NULL;

    errno = 0;

    if (!nelem)
        return ptr;
    if (elsize)
        ptr = calloc (nelem, elsize);
    else
        errno = 1;

    if (ptr && !errno)
        return ptr;

    printf ("ERROR: Failed to allocate memory for %s.\n", name);
    printf ("\nOutput files in directory `%s'.\n", DDSIP_outdir);

    // DDSIP_Free up problem as allocated above
    if (DDSIP_node != NULL)
    {
        DDSIP_FreeFrontNodes ();
        DDSIP_Free ((void **) &DDSIP_node);
    }

    if (DDSIP_data != NULL)
    {
        DDSIP_FreeData ();
        DDSIP_Free ((void **) &DDSIP_data);
    }

    // DDSIP_Free the problem as allocated by CPXcreateprob, if necessary
    if (DDSIP_lp != NULL)
    {
        status = CPXfreeprob (DDSIP_env, &DDSIP_lp);
        if (status)
            printf ("ERROR: CPXfreeprob failed, error code %d\n", status);
    }
    // DDSIP_Free the CPLEX environment, if necessary
    if (DDSIP_env != NULL)
    {
        status = CPXcloseCPLEX (&DDSIP_env);
        if (status)
        {
            char errmsg[1024];
            fprintf (stderr, "ERROR: Failed to close CPLEX environment.\n");
            CPXgeterrorstring (DDSIP_env, status, errmsg);
            printf ("%s\n", errmsg);
        }
    }

    if (DDSIP_bb != NULL)
    {
        DDSIP_FreeBb ();
        DDSIP_Free ((void **) &DDSIP_bb);
    }

    if (DDSIP_param != NULL)
    {
        DDSIP_FreeParam ();
        DDSIP_Free ((void **) &DDSIP_param);
    }

    if (DDSIP_outfile != NULL)
    {
        fclose (DDSIP_outfile);
        DDSIP_outfile = NULL;
    }

    printf ("Terminating DDSIP.\n");
    exit (1);
    // Never reached but ..
    return NULL;
}

//==========================================================================
void
DDSIP_FreeNode (int nono)
{
    int scen, j, cnt;
    if (DDSIP_param->hot)
        DDSIP_Free ((void **) &(DDSIP_node[nono]->solut));
    if ((DDSIP_node[nono])->first_sol)
    {
        for (scen = 0; scen < DDSIP_param->scenarios; scen++)
        {
            if (((DDSIP_node[nono])->first_sol)[scen]
                    && (cnt = (int) ((((DDSIP_node[nono])->first_sol)[scen])[DDSIP_data->firstvar] - 0.9)))
                for (j = scen + 1; cnt && j < DDSIP_param->scenarios; j++)
                {
                    if (((DDSIP_node[nono])->first_sol)[j]
                            && (((DDSIP_node[nono])->first_sol)[scen] == ((DDSIP_node[nono])->first_sol)[j]))
                    {
                        ((DDSIP_node[nono])->first_sol)[j] = NULL;
                        cnt--;
                    }
                }
            DDSIP_Free ((void **) &((DDSIP_node[nono]->first_sol)[scen]));
        }
    }
    DDSIP_Free ((void **) &(DDSIP_node[nono]->first_sol));
    DDSIP_Free ((void **) &(DDSIP_node[nono]->cursubsol));
    DDSIP_Free ((void **) &(DDSIP_node[nono]->subbound));
    DDSIP_Free ((void **) &(DDSIP_node[nono]->mipstatus));
    DDSIP_Free ((void **) &(DDSIP_node[nono]->ref_scenobj));
    if (DDSIP_param->cb)
    {
        DDSIP_Free ((void **) &(DDSIP_node[nono]->dual));
        DDSIP_Free ((void **) &(DDSIP_node[nono]->scenBoundsNoLag));
    }
    // DDSIP_node itself is needed till the end for neoind, neolb, neoub !
}

//==========================================================================
void
DDSIP_FreeFrontNodes ()
{
    int i, j, scen, cnt;
    if (DDSIP_bb->nofront)
    {
        for (i = 0; i < DDSIP_bb->nofront; i++)
        {
            for (scen = 0; scen < DDSIP_param->scenarios; scen++)
            {
                if (((DDSIP_node[DDSIP_bb->front[i]])->first_sol)[scen]
                        && (cnt = (int) ((((DDSIP_node[DDSIP_bb->front[i]])->first_sol)[scen])[DDSIP_data->firstvar] - 0.9)))
                    for (j = scen + 1; cnt && j < DDSIP_param->scenarios; j++)
                    {
                        if (((DDSIP_node[DDSIP_bb->front[i]])->first_sol)[j]
                                && (((DDSIP_node[DDSIP_bb->front[i]])->first_sol)[scen] == ((DDSIP_node[DDSIP_bb->front[i]])->first_sol)[j]))
                        {
                            ((DDSIP_node[DDSIP_bb->front[i]])->first_sol)[j] = NULL;
                            cnt--;
                        }
                    }
                DDSIP_Free ((void **) &((DDSIP_node[DDSIP_bb->front[i]]->first_sol)[scen]));
            }
            DDSIP_FreeNode (DDSIP_bb->front[i]);
            DDSIP_Free ((void **) &DDSIP_node[DDSIP_bb->front[i]]);
        }
    }
    else
    {
        DDSIP_FreeNode (0);
        DDSIP_Free ((void **) &DDSIP_node[0]);
    }
}

//==========================================================================
// Delete structure components

void
DDSIP_FreeData ()
{
    DDSIP_Free ((void **) &(DDSIP_param->heuristic_vector));
    DDSIP_Free ((void **) &(DDSIP_param->cpxeevwhich));
    DDSIP_Free ((void **) &(DDSIP_param->cpxeevwhat));
    DDSIP_Free ((void **) &(DDSIP_param->cpxwhich));
    DDSIP_Free ((void **) &(DDSIP_param->cpxwhat));
    DDSIP_Free ((void **) &(DDSIP_param->cpxlbwhich));
    DDSIP_Free ((void **) &(DDSIP_param->cpxlbwhat));
    DDSIP_Free ((void **) &(DDSIP_param->cpxlbwhich2));
    DDSIP_Free ((void **) &(DDSIP_param->cpxlbwhat2));
    DDSIP_Free ((void **) &(DDSIP_param->cpxubwhich));
    DDSIP_Free ((void **) &(DDSIP_param->cpxubwhat));
    DDSIP_Free ((void **) &(DDSIP_param->cpxubwhich2));
    DDSIP_Free ((void **) &(DDSIP_param->cpxubwhat2));
    if (DDSIP_param->stocrhs)
    {
        DDSIP_Free ((void **) &(DDSIP_data->rhsind));
        DDSIP_Free ((void **) &(DDSIP_data->rhs));
    }
    if (DDSIP_param->stocmat)
    {
        DDSIP_Free ((void **) &(DDSIP_data->matval));
        DDSIP_Free ((void **) &(DDSIP_data->matcol));
        DDSIP_Free ((void **) &(DDSIP_data->matrow));
    }
    DDSIP_Free ((void **) &(DDSIP_data->cost));
    if (DDSIP_param->stoccost)
    {
        DDSIP_Free ((void **) &(DDSIP_data->costind));
    }
    if (DDSIP_param->cb)
    {
        DDSIP_Free ((void **) &(DDSIP_data->nabeg));
        DDSIP_Free ((void **) &(DDSIP_data->nacnt));
        DDSIP_Free ((void **) &(DDSIP_data->naind));
        DDSIP_Free ((void **) &(DDSIP_data->naval));
    }

    DDSIP_Free ((void **) &(DDSIP_data->prob));
    DDSIP_Free ((void **) &(DDSIP_data->obj_coef));
}

//==========================================================================

void
DDSIP_FreeBb ()
{
    int i, k1, cnt;
    sug_t *tmp, *next;

    if (DDSIP_bb != NULL)
    {
        if (DDSIP_bb->sug != NULL)
        {
            for (i = 0; i < DDSIP_param->nodelim + 3; i++)
            {
                tmp = (DDSIP_bb->sug)[i];
                while (tmp)
                {
                    next = tmp->next;
                    DDSIP_Free ((void **) &((*tmp).firstval));
                    DDSIP_Free ((void **) &(tmp));
                    tmp = next;
                }
            }
        }
        DDSIP_Free ((void **) &(DDSIP_bb->sug));

        DDSIP_Free ((void **) &(DDSIP_bb->subsol));

        if (DDSIP_bb->secstage != NULL)
            for (i = 0; i < DDSIP_param->scenarios; i++)
                DDSIP_Free ((void **) &(DDSIP_bb->secstage[i]));
        if (DDSIP_bb->cur_secstage != NULL)
            for (i = 0; i < DDSIP_param->scenarios; i++)
                DDSIP_Free ((void **) &(DDSIP_bb->cur_secstage[i]));
        DDSIP_Free ((void **) &(DDSIP_bb->secstage));
        DDSIP_Free ((void **) &(DDSIP_bb->cur_secstage));
        DDSIP_Free ((void **) &(DDSIP_bb->cost));
        DDSIP_Free ((void **) &(DDSIP_bb->lborg));
        DDSIP_Free ((void **) &(DDSIP_bb->uborg));
        DDSIP_Free ((void **) &(DDSIP_bb->bestsol));
        DDSIP_Free ((void **) &(DDSIP_bb->front));
        DDSIP_Free ((void **) &(DDSIP_bb->lbident));
        DDSIP_Free ((void **) &(DDSIP_bb->ubident));
        if (DDSIP_param->hot && DDSIP_param->hot != 2)
        {
            DDSIP_Free ((void **) &(DDSIP_bb->values));
            DDSIP_Free ((void **) &(DDSIP_bb->effort));
            DDSIP_Free ((void **) &(DDSIP_bb->beg));
            if (DDSIP_bb->Names)
            {
                for (i = 0; i < 2*DDSIP_param->scenarios + 5; i++)
                    DDSIP_Free ((void **) &(DDSIP_bb->Names[i]));
                DDSIP_Free ((void **) &(DDSIP_bb->Names));
            }
        }
        if (DDSIP_param->cb)
        {
            if (DDSIP_bb->bestfirst)
            {
                for (i = 0; i < DDSIP_param->scenarios; i++)
                {
                    if ( DDSIP_bb->bestfirst[i].first_sol)
                    {
                        if ((cnt = (int) (DDSIP_bb->bestfirst[i].first_sol[DDSIP_bb->firstvar] - 0.9)))
                        {
                            for (k1 = i + 1; cnt && k1 < DDSIP_param->scenarios; k1++)
                            {
                                if (DDSIP_bb->bestfirst[i].first_sol == DDSIP_bb->bestfirst[k1].first_sol)
                                {
                                    DDSIP_bb->bestfirst[k1].first_sol = NULL;
                                    cnt--;
                                }
                            }
                        }
                        DDSIP_Free ((void **) &(DDSIP_bb->bestfirst[i].first_sol));
                    }
                }
                DDSIP_Free ((void **) &(DDSIP_bb->bestfirst));
            }
            DDSIP_Free ((void **) &(DDSIP_bb->startinfo_multipliers));
        }
 
        DDSIP_Free ((void **) &(DDSIP_bb->aggregate_time));

        DDSIP_Free ((void **) &(DDSIP_bb->curind));
        DDSIP_Free ((void **) &(DDSIP_bb->curub));
        DDSIP_Free ((void **) &(DDSIP_bb->curlb));
        DDSIP_Free ((void **) &(DDSIP_bb->firstindex));
        DDSIP_Free ((void **) &(DDSIP_bb->firstindex_reverse));
        DDSIP_Free ((void **) &(DDSIP_bb->secondindex));
        DDSIP_Free ((void **) &(DDSIP_bb->secondindex_reverse));
        DDSIP_Free ((void **) &(DDSIP_bb->firsttype));
        DDSIP_Free ((void **) &(DDSIP_bb->sectype));
        DDSIP_Free ((void **) &(DDSIP_bb->objcontrib));
        DDSIP_Free ((void **) &(DDSIP_bb->solstat));
        DDSIP_Free ((void **) &(DDSIP_bb->btlb));
        if (DDSIP_param->hot)
        {
            DDSIP_Free ((void **) &(DDSIP_bb->intind));
            if (DDSIP_bb->intsolvals)
            {
                DDSIP_Free ((void **) &(DDSIP_bb->intsolvals));
            }
            if (DDSIP_bb->boundIncrease_intsolvals)
            {
                DDSIP_Free ((void **) &(DDSIP_bb->boundIncrease_intsolvals));
            }
        }
        if (DDSIP_param->riskmod)
        {
            DDSIP_Free ((void **) &(DDSIP_bb->bestriskval));
            DDSIP_Free ((void **) &(DDSIP_bb->curriskval));
        }
        if ((DDSIP_param->scalarization || DDSIP_param->cb))
        {
            DDSIP_Free ((void **) &(DDSIP_bb->ref_risk));
            DDSIP_Free ((void **) &(DDSIP_bb->local_bestdual));
        }
        DDSIP_Free ((void **) &(DDSIP_bb->ref_scenobj));
        DDSIP_Free ((void **) &(DDSIP_bb->adv_sol));
        DDSIP_Free ((void **) &(DDSIP_bb->lb_scen_order));
        DDSIP_Free ((void **) &(DDSIP_bb->ub_scen_order));
        DDSIP_Free ((void **) &(DDSIP_bb->n_buffer));
        DDSIP_Free ((void **) &(DDSIP_bb->front_nodes_sorted));
        DDSIP_Free ((void **) &(DDSIP_bb->firstrowind));
        DDSIP_Free ((void **) &(DDSIP_bb->secondrowind));
        DDSIP_Free ((void **) &(DDSIP_bb->firstrowind_reverse));
        DDSIP_Free ((void **) &(DDSIP_bb->secondrowind_reverse));
        DDSIP_FreeCutpool();
        if (DDSIP_bb->bestdual)
        {
            bbest_t * tmp_bestdual, * tmp_next;
            tmp_bestdual = DDSIP_bb->bestdual;
            while (tmp_bestdual)
            {
                tmp_next = tmp_bestdual->next;
                DDSIP_Free ((void **) &(tmp_bestdual->dual));
                DDSIP_Free ((void **) &(tmp_bestdual));
                tmp_bestdual = tmp_next;
            }
        }
    }
#ifndef NEOS
#ifndef _WIN32
    {
        char command[256];
        sprintf(command,"ls -l %s*; gzip -f9 %s; ls -l %s*", DDSIP_solfname, DDSIP_solfname, DDSIP_solfname);
        i = system (command);
        if (DDSIP_param->outlev)
        {
            if (i)
            {
                fprintf (DDSIP_bb->moreoutfile, "XXX  return code of 'gzip -9 %s': %d\n", DDSIP_solfname, i);
                printf ("XXX  return code of 'gzip -9 %s': %d\n", DDSIP_solfname, i);
            }
            if (DDSIP_bb->moreoutfile != NULL)
            {
                fclose (DDSIP_bb->moreoutfile);
                DDSIP_bb->moreoutfile = NULL;
            }
            sprintf(command,"ls -l %s*; gzip -f9 %s; ls -l %s*", DDSIP_moreoutfname, DDSIP_moreoutfname, DDSIP_moreoutfname);
            i = system (command);
            if (i)
                printf ("XXX  return code of 'gzip -9 %s': %d\n", DDSIP_moreoutfname, i);
        }
    }
#endif
#endif
    if (DDSIP_bb->moreoutfile != NULL)
    {
        fclose (DDSIP_bb->moreoutfile);
        DDSIP_bb->moreoutfile = NULL;
    }
}

//==========================================================================
void
DDSIP_FreeParam ()
{
    DDSIP_Free ((void **) &(DDSIP_param->coretype));
    DDSIP_Free ((void **) &(DDSIP_param->prefix));
    DDSIP_Free ((void **) &(DDSIP_param->postfix));

    if (DDSIP_param->cpxno)
    {
        DDSIP_Free ((void **) &(DDSIP_param->cpxwhich));
        DDSIP_Free ((void **) &(DDSIP_param->cpxwhat));
        DDSIP_Free ((void **) &(DDSIP_param->cpxisdbl));
    }
    if (DDSIP_param->cpxnolb)
    {
        DDSIP_Free ((void **) &(DDSIP_param->cpxlbwhich));
        DDSIP_Free ((void **) &(DDSIP_param->cpxlbwhat));
        DDSIP_Free ((void **) &(DDSIP_param->cpxlbisdbl));
    }
    if (DDSIP_param->cpxnolb2)
    {
        DDSIP_Free ((void **) &(DDSIP_param->cpxlbwhich2));
        DDSIP_Free ((void **) &(DDSIP_param->cpxlbwhat2));
        DDSIP_Free ((void **) &(DDSIP_param->cpxlbisdbl2));
    }
    if (DDSIP_param->cpxnoub)
    {
        DDSIP_Free ((void **) &(DDSIP_param->cpxubwhich));
        DDSIP_Free ((void **) &(DDSIP_param->cpxubwhat));
        DDSIP_Free ((void **) &(DDSIP_param->cpxubisdbl));
    }
    if (DDSIP_param->cpxnoub2)
    {
        DDSIP_Free ((void **) &(DDSIP_param->cpxubwhich2));
        DDSIP_Free ((void **) &(DDSIP_param->cpxubwhat2));
        DDSIP_Free ((void **) &(DDSIP_param->cpxubisdbl2));
    }
    if (DDSIP_param->cpxnoeev)
    {
        DDSIP_Free ((void **) &(DDSIP_param->cpxeevwhich));
        DDSIP_Free ((void **) &(DDSIP_param->cpxeevwhat));
        DDSIP_Free ((void **) &(DDSIP_param->cpxeevisdbl));
    }
    if (DDSIP_param->cpxnodual)
    {
        DDSIP_Free ((void **) &(DDSIP_param->cpxdualwhich));
        DDSIP_Free ((void **) &(DDSIP_param->cpxdualwhat));
        DDSIP_Free ((void **) &(DDSIP_param->cpxdualisdbl));
    }
    if (DDSIP_param->cpxnodual2)
    {
        DDSIP_Free ((void **) &(DDSIP_param->cpxdualwhich2));
        DDSIP_Free ((void **) &(DDSIP_param->cpxdualwhat2));
        DDSIP_Free ((void **) &(DDSIP_param->cpxdualisdbl2));
    }
    if (DDSIP_param->scalarization)
    {
        DDSIP_Free ((void **) &(DDSIP_param->ref_point));
        DDSIP_Free ((void **) &(DDSIP_param->ref_scale));
    }
}

void
DDSIP_FreeCutpool (void)
{
    cutpool_t *prev, *current;
    current = DDSIP_bb->cutpool;
    while (current)
    {
        prev = current->prev;
        DDSIP_Free ((void **) &(current->matval));
        DDSIP_Free ((void **) &(current));
        current = prev;
    }
    DDSIP_bb->cutpool = NULL;
}
