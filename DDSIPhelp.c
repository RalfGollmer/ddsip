/*  Authors:           Andreas M"arkert, Ralf Gollmer
    Copyright:         University of Duisburg-Essen
    Language:          C
    Last modification: 06.02.2016
	Description:
	This file contains an auxiliary program for DDSIP.
        It's sole purpose is to read the model file and to write a list of variables
        and constraints with their indices in the internal CPLEX representation.
        These indices are needed for the specification of the locations of
        stochastic matrix coefficients in DDSIP's (simple) input format.
        You don't need this for problems with stochastic rhs and costs only.

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

#include <stdio.h>
#include <stdlib.h>
#include <cplex.h>

int main (int argc, char *argv[])
{
    int          maxvarln = 30;

    CPXENVptr    env = NULL;
    CPXLPptr     lp = NULL;

    FILE         *names;

    char         **colname;
    char         *colstore;
    char         **rowname;
    char         *rowstore;
    char         errmsg[1024];

    int          status, i, numcols, numrows;
    int          *surplus;

    // No. of args ?
    if (argc<2)
    {
        fprintf (stderr,"Usage: siphelp model-file [maxlength [w]]\n   Optional parameter maxlength: maximal length of names in the model, default is 30.\n");
        fprintf (stderr,"   Optional third parameter: if present, write the model as check.lp.gz.\n");
        return 1;
    }

    if (argc>2) maxvarln= atoi(argv[2]);

    // Open CPLEX environment
    env = CPXopenCPLEX (&status);
    if ( env == NULL )
    {
        fprintf (stderr, "Could not open CPLEX environment.\n");
        CPXgeterrorstring (env, status, errmsg);
        fprintf (stderr, "%s", errmsg);
        goto TERMINATE;
    }

    // Create problem
    lp = CPXcreateprob (env, &status, "Unnamed");

    // Read lp file
    status = CPXreadcopyprob (env,lp,argv[1],NULL);
    if (status)
    {
        fprintf (stderr, "Failed to create LP.\n");
        CPXgeterrorstring (env, status, errmsg);
        fprintf (stderr, "%s", errmsg);
        goto TERMINATE;
    }

    // Write lp file
    if (argc>3)
    {
        status = CPXwriteprob (env, lp, "check.lp.gz", NULL);
        if ( status )
        {
            fprintf (stderr, "Failed to write LP to disk.\n");
            goto TERMINATE;
        }
    }

    numrows = CPXgetnumrows (env,lp);
    numcols = CPXgetnumcols (env,lp);
    fprintf (stderr, "%d columns, %d rows in model.\n", numcols, numrows);
    colname = (char **) malloc (sizeof(char *) * numcols);
    colstore = (char *) malloc (sizeof(char) * numcols*(maxvarln+1));
    rowname = (char **) malloc (sizeof(char *) *  numrows);
    rowstore = (char *) malloc (sizeof(char) * numrows*(maxvarln+1));
    surplus = (int *) malloc(sizeof(int));

    status = CPXgetcolname (env, lp, colname, colstore, numcols*(maxvarln+1), surplus, 0, numcols-1);
    if ( status )
    {
        fprintf (stderr, "Failed to obtain column names.\n");
        CPXgeterrorstring (env, status, errmsg);
        fprintf (stderr, "%s", errmsg);
        goto TERMINATE;
    }

    status = CPXgetrowname (env, lp, rowname, rowstore, numrows*(maxvarln+1), surplus, 0, numrows-1);
    if ( status )
    {
        fprintf (stderr, "Failed to obtain row names.\n");
        CPXgeterrorstring (env, status, errmsg);
        fprintf (stderr, "%s", errmsg);
        goto TERMINATE;
    }

    if ((names = fopen("rows+cols","w")) == NULL)
    {
        fprintf (stderr,"Failed to open 'rows+cols' for writing.\n");
        return 1;
    }

    fprintf (names,"Indices and names of columns and rows in file %s\n",argv[1]);
    fprintf (names, "Model contains %d rows and %d columns.\n__________________________\nColumns\n", numrows, numcols);
    for (i=0; i<numrows; i++)
        fprintf (names,"%d\t%s\n", i, rowname[i]);
    fprintf (names, "__________________________\n");
    for (i=0; i<numcols; i++)
        fprintf (names,"%d\t%s\n", i, colname[i]);
    fprintf (names, "__________________________\nRows\n");

    fclose (names);

    printf ("Output written on 'rows+cols'\n");

    free (surplus);
    free (rowstore);
    free (rowname);
    free (colname);
    free (colstore);

    if (!(i=system("gzip -9 rows+cols")))
        printf ("compressed file to 'rows+cols.gz'\n");
    else
        printf ("failed to compress file 'rows+cols.out', return code is %d.\n",i);
//==========================================================================

TERMINATE:

// Free up the problem as allocated by CPXcreateprob, if necessary

    if ( lp != NULL )
    {
        status = CPXfreeprob (env, &lp);
        if ( status )
        {
            fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
        }
    }

// Free up the CPLEX environment, if necessary

    if ( env != NULL )
    {
        status = CPXcloseCPLEX (&env);

        /* Note that CPXcloseCPLEX produces no output,
           so the only way to see the cause of the error is to use
           CPXgeterrorstring.  For other CPLEX routines, the errors will
           be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON. */

        if ( status )
        {
            fprintf (stderr, "Could not close CPLEX environment.\n");
            CPXgeterrorstring (env, status, errmsg);
            fprintf (stderr, "%s", errmsg);
        }
    }
    return (status);

}  // END main */

