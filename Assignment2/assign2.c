////////////////////////////////////////////////////////////////////////
// Solution to assignment #1 for ECE1733.
// This program implements the Quine-McCluskey method for 2-level
// minimization. 
////////////////////////////////////////////////////////////////////////

/**********************************************************************/
/*** HEADER FILES *****************************************************/
/**********************************************************************/

#include <stdlib.h>
#include <conio.h>
#include <stdio.h>
#include <assert.h>
#include <malloc.h>
#include <string.h>
#include "common_types.h"
#include "blif_common.h"
#include "cubical_function_representation.h"

/**********************************************************************/
/*** DATA STRUCTURES DECLARATIONS *************************************/
/**********************************************************************/

/**********************************************************************/
/*** DEFINE STATEMENTS ************************************************/
/**********************************************************************/

/**********************************************************************/
/*** GLOBAL VARIABLES *************************************************/
/**********************************************************************/

/**********************************************************************/
/*** FUNCTION DECLARATIONS ********************************************/
/**********************************************************************/


int cube_cost(t_blif_cube *cube, int num_inputs);
int function_cost(t_blif_cubical_function *f);
int cover_cost(t_blif_cube **cover, int num_cubes, int num_inputs);

void simplify_function(t_blif_cubical_function *f);


/**********************************************************************/
/*** BODY *************************************************************/
/**********************************************************************/


/**********************************************************************/
/*** COST FUNCTIONS ***************************************************/
/**********************************************************************/


int cube_cost(t_blif_cube *cube, int num_inputs)
/* Wires and inverters are free, everything else is #inputs+1*/
{
	int index;
	int cost = 0;

	for (index = 0; index < num_inputs; index++)
	{
		if (read_cube_variable(cube->signal_status, index) != LITERAL_DC)
		{
			cost++;
		}
	}
	if (cost > 1)
	{
		cost++;
	}
	return cost;
}


int function_cost(t_blif_cubical_function *f)
{
	int cost = 0;
	int index;
	
	if (f->cube_count > 0)
	{
		for(index = 0; index < f->cube_count; index++)
		{
			cost += cube_cost(f->set_of_cubes[index], f->input_count);
		}
		if (f->cube_count > 1)
		{
			cost += (f->cube_count+1);
		}
	}

	return cost;
}


int cover_cost(t_blif_cube **cover, int num_cubes, int num_inputs)
{
	int result = 0;
	int index;

	for (index = 0; index < num_cubes; index++)
	{
		result += cube_cost(cover[index], num_inputs);
	}
	if (num_cubes > 1)
	{
		result += num_cubes+1;
	}
	return result;
}


/**********************************************************************/
/*** MINIMIZATION CODE ************************************************/
/**********************************************************************/


void simplify_function(t_blif_cubical_function *f)
/* This function simplifies the function f. The minimized set of cubes is
 * returned though a field in the input structure called set_of_cubes.
 * The number of cubes is stored in the field cube_count.
 */
{
	/* PUT YOUR CODE HERE */
}


/**********************************************************************/
/*** MAIN FUNCTION ****************************************************/
/**********************************************************************/


int main(int argc, char* argv[])
{
	t_blif_logic_circuit *circuit = NULL;

	if (argc != 2)
	{
		printf("Usage: %s <source BLIF file>\r\n", argv[0]);
		return 0;
	}
	printf("\t\tECE1733 - Switching Theory - Assignment 2 \n\t\t Reduced Ordered Binary Decision Diagram.\r\n");

	/* Read BLIF circuit. */
	printf("Reading file %s...\n",argv[1]);
	circuit = ReadBLIFCircuit(argv[1]);

	if (circuit != NULL)
	{
		int index;

		/* Minimize each function, one at a time. */
		printf("Minimizing logic functions\n");
		for (index = 0; index < circuit->function_count; index++)
		{
			t_blif_cubical_function *function = circuit->list_of_functions[index];

			simplify_function(function);
		}

		/* Print out synthesis report. */
		printf("Report:\r\n");
		for (index = 0; index < circuit->function_count; index++)
		{
			t_blif_cubical_function *function = circuit->list_of_functions[index];

			/* Print function information. */
			printf("Function %i: #inputs = %i; #cubes = %i; cost = %i\n", index+1, function->input_count, function->cube_count, function_cost(function)); 
		}

		/* Finish. */
		printf("Done.\r\n");
		DeleteBLIFCircuit(blif_circuit);
	}
	else
	{
		printf("Error reading BLIF file. Terminating.\n");
	}
	return 0;
}

