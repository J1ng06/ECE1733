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
#include "contracts.h"

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
int numInputs = 0;
int numPrimeInput = 0;
bool DEBUG = true;

/**********************************************************************/
/*** FUNCTION DECLARATIONS ********************************************/
/**********************************************************************/


int cube_cost(t_blif_cube *cube, int num_inputs);
int function_cost(t_blif_cubical_function *f);
int cover_cost(t_blif_cube **cover, int num_cubes, int num_inputs);



/**********************************************************************/
/*** BODY *************************************************************/
/**********************************************************************/
char translateLiterals(int literal)
{
	char ret = ' ';
	if (literal == LITERAL_0) ret = '0';
	else if (literal == LITERAL_1) ret = '1';
	else if (literal == LITERAL_DC) ret = 'X';
	return ret;
}

void printCube(t_blif_cube *cube, int numInputs)
{
	int j;
	for (j = 0; j < numInputs; j++) {
		printf("%c ", translateLiterals(read_cube_variable(cube->signal_status, j)));
	}
}
void freeSetOfCubes(t_blif_cube **cubes, int cube_count)
{
	int i;
	for (i = 0; i < cube_count; i++) {
		free(cubes[i]);
	}
	free(cubes);
}
void printSetOfCubes(t_blif_cube **cubes, int numInputs, int numCubes)
{
	int i, j;
	// Print a border
	for (i = 0; i < numInputs; i++) {
		printf("==");
	}
	printf("\n");
	for (i = 0; i < numCubes; i++) {
		for (j = 0; j < numInputs; j++) {
			printf("%c ", translateLiterals(read_cube_variable(cubes[i]->signal_status, j)));
		}
		printf("\n");
	}
	// Print a border
	for (i = 0; i < numInputs; i++) {
		printf("==");
	}
	printf("\n");
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

int and(int b1, int b2) {
	return b1 * b2;
}
int or (int b1, int b2) {
	return (b1 >= b2) ? b1 : b2;
}
int xor(int b1, int b2) {
	return (b1 + b2) % 2;
}
//---------------------------------------------------------------------
// create an updated set_of_cubes based on setting an input (idx) to values val
// returns a newly allocated set_of_cubes
//      or a LITERAL 1 or 0 vslue if it has no cofactors
//---------------------------------------------------------------------
int createCofactor(t_blif_cube **oldCubes, int oldCubeCount, int idx, int val,
	t_blif_cube **newCubes, int *newCubeCount)
{
	int i;
	int cnt = 0;
	// find out how many cubes should be in newCubes
	for (i = 0; i < oldCubeCount; i++)
	{
		if ((read_cube_variable(oldCubes[i]->signal_status, idx) == val) ||
			(read_cube_variable(oldCubes[i]->signal_status, idx) == LITERAL_DC))
		{
			cnt++;
		}
	}
	*newCubeCount = cnt;

	// allocate the newCubes
	int j = 0;
	for (i = 0; (i<oldCubeCount) && (j<cnt); i++)
	{
		if ((read_cube_variable(oldCubes[i]->signal_status, idx) == val) ||
			(read_cube_variable(oldCubes[i]->signal_status, idx) == LITERAL_DC))
		{
			newCubes[j] = (t_blif_cube *)malloc(sizeof(t_blif_cube));
			assert(newCubes[j]);
			memcpy(newCubes[j++], oldCubes[i], sizeof(t_blif_cube));
			//newCubes[j++][0] = oldCubes[i][0]; 
		}
	}
	assert(j == cnt);

	if (idx == numInputs - 1)
	{
		// after setting xi to val, if there are no more matching, then f=false
		if (cnt == 0)
		{
			printf("no matching cubes for x%d\n", idx);
			return LITERAL_0;
		}
		// after setting the last variable x_(n-1), if there are matching ones,
		// then f=true
		else
		{
			printf("assigning last input and still found a match\n");
			return LITERAL_1;
		}
	}


	// if there are more xi to assign, then this should not be examined.
	// just assign a garbage value for now
	return -1;
}

//---------------------------------------------------------------------
// return the node index
//---------------------------------------------------------------------
int build(bdd B, t_blif_cube **setOfCubes, int cubeCount, int i, int literalVal)
{
	if (i == numInputs)
	{
		if (literalVal == LITERAL_0)
			return 0;
		else
			return 1;
	}
	else
	{
		// this allocs more space, but that's ok
		t_blif_cube **v0NewCubes = (cubeCount == 0) ? (NULL) :
			((t_blif_cube **)malloc(cubeCount * sizeof(t_blif_cube *)));
		int v0NewCubeCount = 0;
		// 
		printf("\n***** 0-cofactor on x%d\n", i);
		printf("old cube count=%d\n", cubeCount);
		printf("old cube addr=%p\n", setOfCubes);
		if(cubeCount > 0)
		    printf("old cube[0] addr=%p\n", setOfCubes[0]);
		printSetOfCubes(setOfCubes, numInputs, cubeCount);
		//
		int v0Literal = createCofactor(setOfCubes, cubeCount, i, LITERAL_0,
			v0NewCubes, &v0NewCubeCount);
		//
		printf("new cube count=%d\n", v0NewCubeCount);
		if(v0Literal==-1)
		{
		    printSetOfCubes(v0NewCubes, numInputs, v0NewCubeCount);
		}
		//
		int v0 = build(B, v0NewCubes, v0NewCubeCount, i + 1, v0Literal);

		if (v0NewCubes)
			freeSetOfCubes(v0NewCubes, v0NewCubeCount);



		// this allocs more space, but that's ok
		t_blif_cube **v1NewCubes = (cubeCount == 0) ? (NULL) :
			((t_blif_cube **)malloc(cubeCount * sizeof(t_blif_cube *)));
		int v1NewCubeCount = 0;
		//
		printf("\n***** 1-cofactor on x%d\n", i);
		printf("old cube count=%d\n", cubeCount);
		printf("old cube addr=%p\n", setOfCubes);
		if(cubeCount > 0)
		    printf("old cube[0] addr=%p\n", setOfCubes[0]);
		printSetOfCubes(setOfCubes, numInputs, cubeCount);
		//
		int v1Literal = createCofactor(setOfCubes, cubeCount, i, LITERAL_1,
			v1NewCubes, &v1NewCubeCount);
		//
		printf("new cube count=%d\n", v1NewCubeCount);
		if(v1Literal==-1)
		{
		    printSetOfCubes(v1NewCubes, numInputs, v1NewCubeCount);
		}
		//
		int v1 = build(B, v1NewCubes, v1NewCubeCount, i + 1, v1Literal);

		if (v1NewCubes)
			freeSetOfCubes(v1NewCubes, v1NewCubeCount);

		printf("v0Literal=%d \t v1Literal=%d\n", v0Literal, v1Literal);
		int varIdx = curFunc->inputs[i]->data.index + 1;
		printf("varIdx = %d\n", varIdx);
		return make(B, varIdx, v0, v1);
	}
}

typedef struct varOccur* varOccur;
int* varibableOccurence;
int* variableOrderDesc;
int* initVariableOrder;
int findMaxIndex() {
	int max = 0;
	for (int i = 1; i < numPrimeInput; i++) {
		if (varibableOccurence[max] < varibableOccurence[i]) {
			max = i;
		}
	}
	return max;
}
/**
* Init the variable order in descending order according to their
* associated number of vertices in the BDD graph. That is,
* Beeing x and y varibles,
* x precedes y if  NumVerticesWhoseVariableIs(x) is less than NumVerticesWhoseVariableIs(y).
*/
void initVariableOrderDesc() {
	int* tempVariableOccurence = (int*)malloc(numPrimeInput* sizeof(int*));
	int maxIndex;
	memcpy(tempVariableOccurence, varibableOccurence, numPrimeInput* sizeof(int*));

	for (int i = 0; i < numPrimeInput; i++) {
		maxIndex = findMaxIndex();
		variableOrderDesc[i] = maxIndex;
		varibableOccurence[maxIndex] = INT_MIN;
	}
	memcpy(varibableOccurence, tempVariableOccurence, numPrimeInput* sizeof(int*));
}


/**
* Creates a count of how many vertices have a particular variable.
* Stores this count in variableOccurence class member.
*/
void initVariableOcurrence(t_blif_logic_circuit *circuit, int function_index, bdd B) {
	int varSize = 0;
	for (int i = 0; i < circuit->primary_input_count; i++) {
		for (int j = 2; j < B->size; j++) {
			if (B->T[j]->var == (i + 1)) {
				varSize++;
			}
		}
		varibableOccurence[i] = varSize;
		varSize = 0;
	}	
}
void initVarOrder(t_blif_cubical_function *f) {
	for (int i = 0; i < f->input_count; i++) {
		initVariableOrder[i] = f->inputs[i]->data.index;
	}
}
int getPositionOfVariable(int varIndex) {
	for (int i = 0; i < numInputs; i++) {
		if (varIndex == initVariableOrder[i] + 1) {
			return i;
		}
	}
}

void remove_entry(int idx, int newIdx, node *T, int size) {
	int i;
	for (i = 2; i < size; i++) {
		if (T[i]->low == idx)
			T[i]->low = newIdx;
		if (T[i]->high == idx)
			T[i]->high = newIdx;
	}
	return;
}

void simplify_T(node *T, int size) {
	int i, j;
	for (i = 2; i < size; i++) {
		if (T[i]->low < 0 || T[i]->low > size || T[i]->high < 0 || T[i]->high > size) continue;
		if (T[i]->low == T[i]->high && T[i]->low >= 0) {
			
			remove_entry(i, T[i]->low, T, size);

			T[i]->low = -1;
			T[i]->high = -1;
		}
	}

	for (i = 2; i < size; i++) {
		for (j = i + 1; j < size; j++) {
			if ((T[i]->var == T[j]->var) && (T[i]->low == T[j]->low) && (T[i]->high == T[j]->high)) {
				remove_entry(j, i, T, size);
				T[j]->low = -1;
				T[j]->high = -1;
			}
		}
	}
	//printTTtable(T, size);
	return;

}
int validT_count(node *T, int size)
{
	int i;
	int res = 2;
	for (i = 2; i < size; i++) {
		if (T[i]->var >= 0 && T[i]->var < numPrimeInput && T[i]->low >= 0 && T[i]->low < size
			&& T[i]->high >= 0 && T[i]->high < size) res++;
	}
	return res;
}
int get_pos(int *array, int val, int size) {
	int i;
	for (i = 0; i < size; i++) {
		if (array[i] > val) break;
	}
	return i;
}
void clean_T(bdd B, int size) {
	int i, new_Tnum, num_inval_row, invalid_num, *inval_row;
	node* T = xcalloc(BDD_LIMIT, sizeof(node));
	node* tempT = xcalloc(BDD_LIMIT, sizeof(node));
	T = B->T;
	
	node zero = xmalloc(sizeof(struct node));
	node one = xmalloc(sizeof(struct node));
	zero->var = numPrimeInput + 1;
	one->var = numPrimeInput + 1;
	zero->low = -1;
	zero->high = -1;
	one->low = -1;
	one->high = -1;
	new_Tnum = 2;
	/* low and high for zero and one are irrelevant */
	tempT[0] = zero;
	tempT[1] = one;
	num_inval_row = size - validT_count(T, size);

	inval_row = (int *)malloc(num_inval_row * sizeof(int));

	invalid_num = 0;
	for (i = 2; i < size; i++) {
		if (T[i]->low == -1 || T[i]->high == -1) {
			inval_row[invalid_num++] = i;
		}
	}

	for (i = 2; i < size; i++) {
		if (T[i]->low != -1 && T[i]->high != -1) {
			node a = xmalloc(sizeof(struct node));
			a->var = T[i]->var;
			a->low = T[i]->low - get_pos(inval_row, T[i]->low, invalid_num);
			a->high = T[i]->high - get_pos(inval_row, T[i]->high, invalid_num);
			a->swapped = false;
			tempT[new_Tnum] = a;
			new_Tnum++;
		}

	}

	free(B->T);
	B->T = tempT;
	B->size = new_Tnum;
	return;
}


int swap(t_blif_cubical_function *f, int var1, int var2, bdd B) {

	int children[4] = { -1, -1, -1, -1 };
	node a; elem e; int u;
	a = xmalloc(sizeof(struct node));
	e = xmalloc(sizeof(struct elem));
	for (int i = 0; i < B->size; i++) {
		B->T[i]->swapped = false;
	}
	for (int i = 0; i < B->size; i++) {
		if (B->T[i]->var == var1 && !(B->T[i]->swapped) 
			&& B->T[i]->low >=0 && B->T[i]->low < B->size 
			&& B->T[i]->high >= 0 && B->T[i]->high < B->size) {
		

			if (B->T[B->T[i]->low]->var == var2 && B->T[i]->low > 1) {
				children[0] = B->T[B->T[i]->low]->low;
				children[1] = B->T[B->T[i]->low]->high;
			}
			else {
				children[0] = B->T[i]->low;
				children[1] = B->T[i]->low;
				a->var = var1;
				a->low = children[0];
				a->high = children[2];
				a->swapped = false;
				u = B->size;
				B->size++;
				B->T[u] = a;
				B->T[B->T[i]->low]->swapped = true;
			}

			if (B->T[B->T[i]->high]->var == var2 && B->T[i]->high > 1) {
				children[2] = B->T[B->T[i]->high]->low;
				children[3] = B->T[B->T[i]->high]->high;
			}
			else {
				children[2] = B->T[i]->high;
				children[3] = B->T[i]->high;
				a->var = var1;
				a->low = children[1];
				a->high = children[3];
				a->swapped = false;
				u = B->size;
				B->size++;
				B->T[u] = a;
				B->T[B->T[i]->high]->swapped = true;
			}

			if (B->T[B->T[i]->low]->var == var2 && B->T[i]->low > 1) {
				B->T[B->T[i]->low]->var = var1;
				B->T[B->T[i]->low]->high = children[2];
				B->T[B->T[i]->low]->swapped = true;
			}

			if (B->T[B->T[i]->high]->var == var2 && B->T[i]->high > 1) {
				B->T[B->T[i]->high]->var = var1;
				B->T[B->T[i]->high]->low = children[1];
				B->T[B->T[i]->high]->swapped = true;
			}
			B->T[i]->var = var2;
			B->T[i]->swapped = true;

			simplify_T(B->T, B->size);
			
			clean_T(B, B->size);
			
		}
	}
	
	printTtable(B);
	return B->size;
}

int findBestBackwardPosition(int varIndex, int varIndexBestPosition) {


}
int findBestForwardPosition(t_blif_cubical_function *f, int varIndex, int varIndexBestPosition, bdd B){
	int varIPostion, varJPostion;
	varIPostion = getPositionOfVariable(varIndex);
	for (varJPostion = varIPostion + 1; varJPostion < numInputs; varJPostion++) {
		int sizeAfter = swap(f, initVariableOrder[varIPostion] + 1, initVariableOrder[varJPostion] + 1, B);

		int tempIndex = initVariableOrder[varIPostion];
		initVariableOrder[varIPostion] = initVariableOrder[varJPostion];
		initVariableOrder[varJPostion] = tempIndex;

		varIPostion++;

		

	}
}

void findBestPositionForVariable(int varIndex, t_blif_cubical_function *f, bdd B) {
	int varIndexPosition = getPositionOfVariable(varIndex);
	int varIndexBestPosition = varIndexPosition;

	//if (varIndexPosition <= numInputs / 2) {
	//	varIndexBestPosition = findBestBackwardPosition(varIndex, varIndexBestPosition);
		varIndexBestPosition = findBestForwardPosition(f, varIndex, varIndexBestPosition, B);
	//}
}

void sifting(t_blif_logic_circuit *circuit, int function_index, bdd B) {
	int sizeBefore, sizeAfter;
	bdd tempB = xmalloc(sizeof(struct bdd));
	tempB->T = xcalloc(B->limit, sizeof(node));
	
	varibableOccurence = xcalloc(circuit->primary_input_count, sizeof(int));
	variableOrderDesc = xcalloc(circuit->primary_input_count, sizeof(int));
	initVariableOrder = xcalloc(circuit->list_of_functions[function_index]->input_count, sizeof(int));
	initVariableOcurrence(circuit, function_index, B);
	initVarOrder(circuit->list_of_functions[function_index]);
	initVariableOrderDesc();
	if (DEBUG) {
		for (int i = 0; i < circuit->primary_input_count; i++) {
			printf("x%d has %d occurence\n", i + 1, varibableOccurence[i]);
		}
		printf("Variable order:\n");
		for (int i = 0; i < circuit->list_of_functions[function_index]->input_count; i++) {
			printf("%d: x%d\n", i + 1, initVariableOrder[i] + 1);
		}
		printf("Variable order Desc:\n");
		for (int i = 0; i < circuit->primary_input_count; i++) {
			printf("%d: x%d has %d occurence\n", i + 1, variableOrderDesc[i] + 1, varibableOccurence[variableOrderDesc[i]]);
		}
	}
	sizeBefore = B->size;
	for (int i = 0; i < numPrimeInput; i++) {
		if (varibableOccurence[variableOrderDesc[i]] != 0) {
			if (DEBUG) {
				printf("Start searching for best position for variable x%d\n" , variableOrderDesc[i]+1);
			}
			printf("Before sifting for x%d!\n", variableOrderDesc[i] + 1);
			printTtable(B);
			memcpy(tempB->T, B->T, B->limit*sizeof(node));
			
			
			tempB->size = B->size;
			printf("TempB!\n");
			printTtable(tempB);
			findBestPositionForVariable(variableOrderDesc[i]+1, circuit->list_of_functions[function_index], B);
			memcpy(B->T, tempB->T,B->limit*sizeof(node));
			B->size = tempB->size;
			if (DEBUG) {
				printf("One Iteration Done!\n");
				printTtable(B);
				printf("TempB!\n");
				printTtable(tempB);
			}
		}
	}
	sizeAfter = B->size;
}

/**********************************************************************/
/*** MAIN FUNCTION ****************************************************/
/**********************************************************************/
void Testing() {
	//int i;
	bdd *B = (bdd*)malloc(3*sizeof(bdd*));
	B[0] = bdd_new(5);
	B[1] = bdd_new(5);
	B[2] = bdd_new(5);
	int* u1 = xcalloc(37, sizeof(int));
	int* u2 = xcalloc(37, sizeof(int));
	u1[1] = make(B[0], 5, 1, 0);
	u1[2] = make(B[0], 4, u1[1], 0);
	u1[3] = make(B[0], 4, 0, u1[1]);
	u1[4] = make(B[0], 3, u1[2], u1[3]);
	u1[5] = make(B[0], 2, u1[4], 0);
	u1[6] = make(B[0], 2, 0, u1[4]);
	u1[7] = make(B[0], 1, u1[5], u1[6]);

	printTtable(B[0]);
	printf("*****************\n");
	u2[1] = make(B[1], 5, 1, 0);
	u2[2] = make(B[1], 3, u2[1], 0);
	u2[3] = make(B[1], 3, 0, u2[1]);
	u2[4] = make(B[1], 1, u2[2], u2[3]);
	printTtable(B[1]);
	apply(B, &and, u1[7], u2[4]);
	printf("*****************\n");
	printTtable(B[2]);
}

int main(int argc, char* argv[])
{
	t_blif_logic_circuit *circuit = NULL;
	int var_count;
	if (argc != 2)
	{
		printf("Usage: %s <source BLIF file>\r\n", argv[0]);
		return 0;
	}
	//fp = fopen("output.txt", "w");  // use "a" to append text to the file
	//if (fp == NULL) {
	//	perror("Error opening file.");
	//}
	printf("\t\tECE1733 - Switching Theory - Assignment 2 \n\t\t Reduced Ordered Binary Decision Diagram.\r\n");

	/* Read BLIF circuit. */
	printf("Reading file %s...\n",argv[1]);
	circuit = ReadBLIFCircuit(argv[1]);

	if (circuit != NULL)
	{
		int index;
		bdd B;
		B = bdd_new(circuit->primary_input_count);
		numPrimeInput = circuit->primary_input_count;
		for (index = 0; index < circuit->function_count; index++)
		{
			t_blif_cubical_function *function = circuit->list_of_functions[index];
			//curFunc = (t_blif_cubical_function *)malloc(sizeof(t_blif_cubical_function));
			//curFunc->set_of_cubes = (t_blif_cube*)malloc(numInputs * sizeof(t_blif_cube));
			
			numInputs = function->input_count;
			
			//B[index] = bdd_new(numInputs);
			curFunc = function;
			//memcpy(curFunc->set_of_cubes, function->set_of_cubes, numInputs * sizeof(t_blif_cube));

			print_function_info(function, index);
			
			//initH();
			printf("\n");
			//build(B[index], function->set_of_cubes,
			build(B, function->set_of_cubes,
				function->cube_count, 0, function->value);
			printTtable(B);
			sifting(circuit, index, B);
		}
		//Testing();

		//printTtable(B);
		///* Print out synthesis report. */
		//printf("Report:\r\n");
		//for (index = 0; index < circuit->function_count; index++)
		//{
		//	t_blif_cubical_function *function = circuit->list_of_functions[index];

		//	/* Print function information. */
		//	printf("Function %i: #inputs = %i; #cubes = %i; cost = %i\n", index+1, function->input_count, function->cube_count, function_cost(function)); 
		//}

		//int i;
		//int* u = xcalloc(37, sizeof(int));

		//u[1] = make(B, 2, 0, 1);	/* u[1] = x2 */
		//u[2] = make(B, 1, 0, u[1]);	/* u[2] = x1 /\ x2 */
		//onesat(B, u[2]);
		//assert(satcount(B, u[2]) == 1);
		//u[3] = apply(B, &xor, u[2], 1); /* u[3] = ~(x1 /\ x2) */
		//u[4] = apply(B, &and, u[2], u[3]);	/* u[4] = u[2] /\ u[3] = 0 */
		//assert(satcount(B, u[4]) == 0);
		//assert(u[4] == 0);
		//u[5] = make(B, 2, 0, 1);	/* u[3] = x2 */
		//u[6] = make(B, 1, u[5], 1);	/* u[4] = x1 \/ x2 */
		//onesat(B, u[6]);
		//assert(satcount(B, u[6]) == 3);
		//u[7] = make(B, 1, 0, 1);	/* u[7] = x1 */
		//u[8] = make(B, 2, 0, 1);	/* u[8] = x2 */
		//u[9] = apply(B, &and, u[7], u[8]); 	/* u[9] = x1 /\ x2 */
		//for (i = 1; i < 10; i++) {
		//	printf("u[%d] = %d\n", i, u[i]);
		//}
		//assert(u[9] == u[2]);
		//// print("ROBDD of %d", B->num_vars); //printint();
		// //print(" vars has "); printint(B->size); print(" nodes\n");
		//printf("passed all tests!\n");
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

