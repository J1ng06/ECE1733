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
#include "BDD.h"
#include <time.h>

/**********************************************************************/
/*** DATA STRUCTURES DECLARATIONS *************************************/
/**********************************************************************/

/**********************************************************************/
/*** DEFINE STATEMENTS ************************************************/
/**********************************************************************/

/**********************************************************************/
/*** GLOBAL VARIABLES *************************************************/
/**********************************************************************/
t_blif_cubical_function *curFunc;
int input_count = 0;
int numPrimeInput = 0;
bool DEBUG = true;
int *rootNode;

/**********************************************************************/
/*** FUNCTION DECLARATIONS ********************************************/
/**********************************************************************/


int cube_cost(t_blif_cube *cube, int num_inputs);
int function_cost(t_blif_cubical_function *f);
int cover_cost(t_blif_cube **cover, int num_cubes, int num_inputs);



/**********************************************************************/
/*** BODY *************************************************************/
/**********************************************************************/

int op_and(int op1, int op2) {
	return op1 * op2;
}
int op_or (int op1, int op2) {
	return (op1 >= op2) ? op1 : op2;
}
int op_xor(int op1, int op2) {
	return (op1 + op2) % 2;
}

void freeCover(t_blif_cube **cubes, int cube_count)
{
	int i;
	for (i = 0; i < cube_count; i++) {
		free(cubes[i]);
	}
	free(cubes);
}

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
		for (index = 0; index < f->cube_count; index++)
		{
			cost += cube_cost(f->set_of_cubes[index], f->input_count);
		}
		if (f->cube_count > 1)
		{
			cost += (f->cube_count + 1);
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
		result += num_cubes + 1;
	}
	return result;
}


/**********************************************************************/
/*** UTILITY CODE *****************************************************/
/**********************************************************************/

void print_function_info(t_blif_cubical_function *function, int function_index)
{
	int cube_index;

	printf("Function %i: #inputs = %i; #cubes = %i;\n",
		function_index + 1, function->input_count, function->cube_count);

	for (cube_index = 0; cube_index < function->cube_count; cube_index++)
	{
		printf("Cube %d: ", cube_index);
		print_cube(function->input_count, function->set_of_cubes[cube_index]);
	}
}

// create Cofactor
int cofactor(t_blif_cube **oldCover, int oldCubeCount, int varIndex, int val,
	t_blif_cube **newCover, int *newCubeCount)
{
	int i;
	int count = 0;
	for (i = 0; i < oldCubeCount; i++)
	{
		if ((read_cube_variable(oldCover[i]->signal_status, varIndex) == val) ||
			(read_cube_variable(oldCover[i]->signal_status, varIndex) == LITERAL_DC))
		{
			count++;
		}
	}
	*newCubeCount = count;

	int j = 0;
	for (i = 0; (i<oldCubeCount) && (j<count); i++)
	{
		if ((read_cube_variable(oldCover[i]->signal_status, varIndex) == val) ||
			(read_cube_variable(oldCover[i]->signal_status, varIndex) == LITERAL_DC))
		{
			newCover[j] = (t_blif_cube *)malloc(sizeof(t_blif_cube));
			memcpy(newCover[j++], oldCover[i], sizeof(t_blif_cube));
		}
	}

	if (varIndex == input_count - 1)
	{
		if (count == 0)
		{
			return LITERAL_0;
		}
		else
		{
			return LITERAL_1;
		}
	}
	return -1;
}
//---------------------------------------------------------------------
// Build BDD 
//---------------------------------------------------------------------
int build(bdd B, t_blif_cube **inputCover, int cubeCount, int i, int literalVal)
{
	if (i == input_count)
	{
		if (literalVal == LITERAL_0)
			return 0;
		else
			return 1;
	}
	else
	{
		t_blif_cube **v0NewCover = (cubeCount == 0) ? (NULL) :
			((t_blif_cube **)malloc(cubeCount * sizeof(t_blif_cube *)));
		int v0NewCubeCount = 0;
		int v0Literal = cofactor(inputCover, cubeCount, i, LITERAL_0,
			v0NewCover, &v0NewCubeCount);

		int v0 = build(B, v0NewCover, v0NewCubeCount, i + 1, v0Literal);

		if (v0NewCover)
			freeCover(v0NewCover, v0NewCubeCount);

		t_blif_cube **v1NewCover = (cubeCount == 0) ? (NULL) :
			((t_blif_cube **)malloc(cubeCount * sizeof(t_blif_cube *)));
		int v1NewCubeCount = 0;
		int v1Literal = cofactor(inputCover, cubeCount, i, LITERAL_1,
			v1NewCover, &v1NewCubeCount);
		int v1 = build(B, v1NewCover, v1NewCubeCount, i + 1, v1Literal);

		if (v1NewCover)
			freeCover(v1NewCover, v1NewCubeCount);
		int varIdx = curFunc->inputs[i]->data.index;
		return MK(B, varIdx + 1, v0, v1);
	}
}


struct apnode {
	bdd_node u1;
	bdd_node u2;
	bdd_node u12;			
};
typedef struct apnode* apnode;

void apfun_free(ht_elem ap) {
	free((apnode)ap);
}

ht_key ap_key(ht_elem ap) {
	return ap;
}

bool ap_equal(ht_key ap1, ht_key ap2) {
	return ((apnode)ap1)->u1 == ((apnode)ap2)->u1
		&& ((apnode)ap1)->u2 == ((apnode)ap2)->u2;
}

int ap_hash(ht_key ap, int m) {
	unsigned int x = 162885;
	unsigned int y = 1990223;
	unsigned int r = 0xdeadbeef;
	unsigned int h = (unsigned)((apnode)ap)->u1;
	r = r*x + y;
	h = r*h + (unsigned)((apnode)ap)->u2;
	h = h % (unsigned)m;
	return (int)h;
}
bdd_node apply_rec(bdd *B, int(*func)(int op1, int op2),
	bdd_node u1, bdd_node u2, table G)
{
	int u;
	if (u1 <= 1 && u2 <= 1)
		return (*func)(u1, u2);
	else {
		apnode temp = xmalloc(sizeof(struct apnode));
		temp->u1 = u1; temp->u2 = u2;
		apnode apNode = table_search(G, temp);
		if (apNode != NULL) {
			free(temp);
			return apNode->u12;
		}
		int v1 = (B[0]->T[u1])->var;
		int v2 = (B[1]->T[u2])->var;
		if (v1 == v2)
			u = MK(B[2], v1,
				apply_rec(B, func, (B[0]->T[u1])->low, (B[1]->T[u2])->low, G),
				apply_rec(B, func, (B[0]->T[u1])->high, (B[1]->T[u2])->high, G));
		else if (v1 < v2)
			u = MK(B[2], v1,
				apply_rec(B, func, (B[0]->T[u1])->low, u2, G),
				apply_rec(B, func, (B[0]->T[u1])->high, u2, G));
		else
			u = MK(B[2], v2,
				apply_rec(B, func, u1, (B[1]->T[u2])->low, G),
				apply_rec(B, func, u1, (B[1]->T[u2])->high, G));
		temp->u12 = u;
		table_insert(G, temp);
		return u;
	}
}

/* apply(B, f, u1, u2) = f on (u1, u2) in B
* using a table for memoization.  f should be
* a boolean function, where bools are encoded as ints 0 and 1
*/
bdd_node apply(bdd *B, int(*op)(int op1, int op2),
	bdd_node u1, bdd_node u2)
{
	table G = table_new(APPLY_HASHTABLE_SIZE,
		&ap_key, &ap_equal, &ap_hash);
	int u = apply_rec(B, op, u1, u2, G);
	table_free(G, &apfun_free);
	return u;
}


/**********************************************************************/
/*** MAIN FUNCTION ****************************************************/
/**********************************************************************/
void Testing() {
	//int i;
	bdd *B = (bdd*)malloc(3*sizeof(bdd*));
	B[0] = bdd_init(5);
	B[1] = bdd_init(5);
	B[2] = bdd_init(5);
	int* u1 = xcalloc(37, sizeof(int));
	int* u2 = xcalloc(37, sizeof(int));
	u1[1] = MK(B[0], 5, 1, 0);
	u1[2] = MK(B[0], 4, u1[1], 0);
	u1[3] = MK(B[0], 4, 0, u1[1]);
	u1[4] = MK(B[0], 3, u1[2], u1[3]);
	u1[5] = MK(B[0], 2, u1[4], 0);
	u1[6] = MK(B[0], 2, 0, u1[4]);
	u1[7] = MK(B[0], 1, u1[5], u1[6]);

	printTtable(B[0]);
	printf("*****************\n");
	u2[1] = MK(B[1], 5, 1, 0);
	u2[2] = MK(B[1], 3, u2[1], 0);
	u2[3] = MK(B[1], 3, 0, u2[1]);
	u2[4] = MK(B[1], 1, u2[2], u2[3]);
	printTtable(B[1]);
	apply(B, &op_and, u1[7], u2[4]);
	printf("*****************\n");
	printTtable(B[2]);
}

int main(int argc, char* argv[])
{
	t_blif_logic_circuit *circuit = NULL;
	clock_t time1, time2;
	time1 = clock();
	int var_count;
	char* opt;
	if (argc != 3)
	{
		printf("Usage: %s <source BLIF file> <Operator>, where Operator could be AND, OR or XOR\r\n", argv[0]);
		return 0;
	}
	else {
		opt = argv[2];
		if (!(!strcmp(opt, "AND") || !strcmp(opt, "OR") || !strcmp(opt, "XOR"))) {
			printf("Operator should be AND, OR or XOR\r\n");
			return 0;
		}
	}
	//fp = fopen("output.txt", "w");  // use "a" to append text to the file
	//if (fp == NULL) {
	//	perror("Error opening file.");
	//}
	printf("\t\tECE1733 - Switching Theory - Assignment 2 \n\t\t Reduced Ordered Binary Decision Diagram.\r\n");

	/* Read BLIF circuit. */
	printf("Reading file %s...\n",argv[1]);
	circuit = ReadBLIFCircuit(argv[1]);
	if (circuit->function_count != 2) {
		printf("\n[ERROR] Please make sure the BLIF file has exactly TWO functions!\n");
		return 0;
	}
	if (circuit != NULL)
	{
		int index;
		numPrimeInput = circuit->primary_input_count;
		rootNode = xcalloc(circuit->function_count, sizeof(int));
		bdd *B = (bdd*)malloc(3 * sizeof(bdd*));
		B[0] = bdd_init(numPrimeInput);
		B[1] = bdd_init(numPrimeInput);
		B[2] = bdd_init(numPrimeInput);
		for (index = 0; index < circuit->function_count; index++)
		{
			t_blif_cubical_function *function = circuit->list_of_functions[index];
			//curFunc = (t_blif_cubical_function *)malloc(sizeof(t_blif_cubical_function));
			//curFunc->set_of_cubes = (t_blif_cube*)malloc(numInputs * sizeof(t_blif_cube));
			
			input_count = function->input_count;
			
			
			//B[index] = bdd_new(numInputs);
			curFunc = function;
			//memcpy(curFunc->set_of_cubes, function->set_of_cubes, numInputs * sizeof(t_blif_cube));

			print_function_info(function, index);
			
			//printf("\n");
			//build(B[index], function->set_of_cubes,
			rootNode[index] = build(B[index], function->set_of_cubes,
				function->cube_count, 0, function->value);
			printTtable(B[index]);

		}
		if (!strcmp(opt, "AND")){
			apply(B, &op_and, rootNode[0], rootNode[1]);
		}
		else if (!strcmp(opt, "OR")) {
			apply(B, &op_or, rootNode[0], rootNode[1]);
		}
		else if (!strcmp(opt, "XOR")) {
			apply(B, &op_xor , rootNode[0], rootNode[1]);
		}
		//Testing();
		printf("\nAPPLY RESULT:\n");
		printTtable(B[2]);
		time2 = clock();
		printf("Time for Apply = %f ms\n", (((double)time2 - (double)time1) / CLOCKS_PER_SEC) * 1000);
		printf("Done.\r\n");
		DeleteBLIFCircuit(blif_circuit);
	}
	else
	{
		printf("Error reading BLIF file. Terminating.\n");
	}
	return 0;
}

