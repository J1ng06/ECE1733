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
#include "assign2.h"
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
bool DEBUG = false;
bestOptions BO;
/**********************************************************************/
/*** FUNCTION DECLARATIONS ********************************************/
/**********************************************************************/


int cube_cost(t_blif_cube *cube, int num_inputs);
int function_cost(t_blif_cubical_function *f);
int cover_cost(t_blif_cube **cover, int num_cubes, int num_inputs);



/**********************************************************************/
/*** BODY *************************************************************/
/**********************************************************************/

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
int* varibableOccurence;
int* variableOrderDesc;
int* initVariableOrder;
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
		return MK(B, varIdx, v0, v1);
	}
}


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
			if (B->T[j]->var == i) {
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
	for (int i = 0; i < input_count; i++) {
		if (varIndex == initVariableOrder[i]) {
			return i;
		}
	}
}
int getVariableInPosition(int level) {
	return initVariableOrder[level];
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

int getVerticesWhoseVariableIs(node *vertices, int var, bdd B, int *u) {
	int node_count = 0;
	for (int i = 0; i < B->size; i++) {
		if (B->T[i]->var == initVariableOrder[var]) {
			vertices[node_count] = xmalloc(sizeof(struct node));
			u[node_count] = i;
			memcpy(vertices[node_count], B->T[i], sizeof(struct node));
			vertices[node_count]->low = B->T[i]->low;
			vertices[node_count]->high = B->T[i]->high;
			vertices[node_count]->var = B->T[i]->var;
			vertices[node_count]->level = B->T[i]->level;
			vertices[node_count]->deletion_flag = B->T[i]->deletion_flag;
			node_count++;
		}
	}
	return node_count;
}
bool nodeValidCheck(node node, int size) {
	if (node->low < size && node->high < size && node->low >= 0 && node->high >= 0) {
		return true;
	}
	return false;
}
int addWithoutRedundant(int var, int low, int high, bdd B, int level) {
	if (low == high) {
		return low;
	}
	node a = xmalloc(sizeof(struct node));
	int u;
	a->var = var;
	a->low = low;
	a->high = high;
	a->deletion_flag = false;
	a->level = getPositionOfVariable(var) + 1;
	elem e = table_search(B->H, a);
	if (e != NULL) {
		free(a);			
		return e->u;		
	}
	u = B->size;
	B->size++;
	B->T[u] = a;
	return u;
}
int findVertexParent(bdd B, node v, int parentNodeIndex) {
	elem e = table_search(B->H, v);
	int count = 0;
	if (e != NULL) {
		for (int i = 2; i < B->size; i++) {
			if (B->T[i]->low == e->u || B->T[i]->high == e->u) {
				if (B->T[i]->var != parentNodeIndex) {
					return true;
					break;
				}
			}
		}
	}
	return false;

}
bool swapVertexWithDescendantsWithVariable(node v, int varJ, bdd B, int u) {
	bool swapWasMade = false;
	int varI = v->var;

	node low = xmalloc(sizeof(struct node));
	node high = xmalloc(sizeof(struct node));
	int a, b, c, d;
	int varJLow, varJHigh;
	memcpy(low, B->T[v->low], sizeof(struct node));
	low->var = B->T[v->low]->var;
	low->low = B->T[v->low]->low;
	low->high = B->T[v->low]->high;
	low->level = B->T[v->low]->level;
	low->deletion_flag = B->T[v->low]->deletion_flag;
	memcpy(high, B->T[v->high], sizeof(struct node));
	high->low = B->T[v->high]->low;
	high->high = B->T[v->high]->high;
	high->var = B->T[v->high]->var;
	high->level = B->T[v->high]->level;
	high->deletion_flag = B->T[v->high]->deletion_flag;
	bool keepLow = findVertexParent(B, low, v->var);
	bool keepHigh = findVertexParent(B, high, v->var);
	B->T[u]->level = getPositionOfVariable(v->var) + 1;
	if (nodeValidCheck(low, B->size) && low->var == varJ 
		&& (!nodeValidCheck(high, B->size) || high->var != varJ)) {
		//printf("Case A\n");
		a = B->T[v->low]->low;
		b = B->T[v->low]->high;
		c = v->high;
		varJLow =addWithoutRedundant(varI, a, c, B, v->var);
		varJHigh = addWithoutRedundant(varI, b, c, B, v->var);
		if (!keepLow) {
			B->T[v->low]->deletion_flag = true;
		}
		else {
			B->T[v->low]->level = B->T[v->low]->level - 1;
			B->T[v->low]->deletion_flag = false;
		}
		B->T[u]->var = varJ;
		B->T[u]->low = varJLow;
		B->T[u]->high = varJHigh;
		B->T[u]->level = getPositionOfVariable(v->var);
		//printf("Case A Debug\n");
		//printTtable(B);
		swapWasMade = true;
	}
	else if ((!nodeValidCheck(low, B->size) || low->var != varJ) 
		&& (nodeValidCheck(high, B->size) && high->var == varJ)) {
		//printf("Case B\n");
		a = v->low;
		b = B->T[v->high]->low;
		c = B->T[v->high]->high;
		varJLow = addWithoutRedundant(varI, a, b, B, v->level);
		varJHigh = addWithoutRedundant(varI, a, c, B, v->level);
		if (!keepHigh) {
			B->T[v->high]->deletion_flag = true;
			
		}
		else {
			B->T[v->high]->level = B->T[v->high]->level - 1;
			B->T[v->high]->deletion_flag = false;
		}
		B->T[u]->var = varJ;
		B->T[u]->low = varJLow;
		B->T[u]->high = varJHigh;
		B->T[u]->level = getPositionOfVariable(v->var);
		//printTtable(B);
		swapWasMade = true;
	}
	else if ((nodeValidCheck(low, B->size) && low->var == varJ)
		&& (nodeValidCheck(high, B->size) && high->var == varJ)) {
		//printf("Case C\n");
		a = B->T[v->low]->low;
		b = B->T[v->low]->high;
		c = B->T[v->high]->low;
		d = B->T[v->high]->high;
		varJLow = addWithoutRedundant(varI, a, c, B, v->level);
		varJHigh = addWithoutRedundant(varI, b, d, B, v->level);
		// update varJ nodes
		//addWithoutRedundant(varJ, varJLow, varJHigh, B);
		if (!keepLow) {
			B->T[v->low]->deletion_flag = true;
			
		}
		else {
			B->T[v->low]->level = B->T[v->low]->level - 1;
			B->T[v->low]->deletion_flag = false;
		}
		if (!keepHigh) {
			B->T[v->high]->deletion_flag = true;
			
		}
		else {
			B->T[v->high]->level = B->T[v->high]->level - 1;
			B->T[v->high]->deletion_flag = false;
		}
		B->T[u]->var = varJ;
		B->T[u]->low = varJLow;
		B->T[u]->high = varJHigh;
		B->T[u]->level = getPositionOfVariable(v->var);
		//B->T[u]->deletion_flag = true;
		//printTtable(B);
		swapWasMade = true;
	}
	else if ((!nodeValidCheck(low, B->size) || low->var != varJ)
		&& (!nodeValidCheck(high, B->size) || high->var != varJ)) {
		swapWasMade = false;
		//printf("Case D\n");
	}
	else if ((!nodeValidCheck(low, B->size) || low->var != varJ) 
		&& !nodeValidCheck(high, B->size)) {
		swapWasMade = false;
		//printf("Case E\n");
	}
	//else { //printf("Case F\n"); }
	
	return swapWasMade;
}

bool swap(t_blif_cubical_function *f, int level, bdd B) {
	bool vertexSwapWasMade = false;
	int *u = xmalloc(numPrimeInput, sizeof(int));
	node *verticesOfLevel = xcalloc(B->size, sizeof(node));
	int varI = getVariableInPosition(level);
	int varJ = getVariableInPosition(level + 1);
	int node_count_for_this_level = 
		getVerticesWhoseVariableIs(verticesOfLevel, level, B, u);
	for (int i = 0; i < node_count_for_this_level; i++) {
		bool swapWasMadeVertexV = swapVertexWithDescendantsWithVariable(verticesOfLevel[i], varJ, B, u[i]);
		vertexSwapWasMade = vertexSwapWasMade || swapWasMadeVertexV;

	}
    
	//printf("swap done for x%d at level %d\n", varI, level);
	//printTtable(B);
	return true;
}

void cleanTable(bdd B) {
	int i;
	elem e;
    table_free(B->H, &elem_free);

	i = 0;
	while (i < B->size) {
		if (B->T[i]->deletion_flag == true) {
			for (int j = i; j < B->size - 1; j++) {
				B->T[j]->var = B->T[j + 1]->var;
				B->T[j]->low = (B->T[j + 1]->low > i) ? B->T[j + 1]->low - 1: B->T[j + 1]->low;
				B->T[j]->high = (B->T[j + 1]->high > i) ? B->T[j + 1]->high-1 : B->T[j + 1]->high;
				B->T[j]->deletion_flag = B->T[j + 1]->deletion_flag;
				B->T[j]->level = B->T[j+1]->level;
			}
			B->T[B->size] = NULL;
			B->size--;
		}
		else {
			i++;
		}
	}	
	B->H = table_new(BDD_HASHTABLE_SIZE,
		&elem_key, &node_equal, &hash);
	for (i = 2; i < B->size; i++) {
		/* enter into hash table */
		e = xmalloc(sizeof(struct elem));
		e->node = B->T[i];
		e->u = i;
		table_insert(B->H, e);
	}
	
}
void reordingTable(bdd B) {
	node *newT = xcalloc(B->limit, sizeof(node));
	elem e;
	int i = 0;
	for (i = 2; i < B->size; i++) {
		B->T[i]->newLow = B->T[i]->low;
		B->T[i]->newHigh = B->T[i]->high;
	}
	for (i = 0; i < 2; i++) {
		newT[i] = xmalloc(sizeof(struct node));
		memcpy(newT[i], B->T[i], sizeof(struct node));
		newT[i]->low = B->T[i]->low;
		newT[i]->high = B->T[i]->high;
		newT[i]->var = B->T[i]->var;
		newT[i]->level = B->T[i]->level;
		newT[i]->deletion_flag = B->T[i]->deletion_flag;
	}
	int newTcount = 2;
	int level = numPrimeInput - 1;
	i = 2;
	while (i < B->size) {
		for (int j = 0; j < B->size; j++) {
			if (B->T[j]->level == level && level == numPrimeInput - 1) {
				newT[newTcount] = xmalloc(sizeof(struct node));
				memcpy(newT[newTcount], B->T[j], sizeof(struct node));
				newT[newTcount]->low = B->T[j]->low;
				newT[newTcount]->high = B->T[j]->high;
				newT[newTcount]->newLow = B->T[j]->low;
				newT[newTcount]->newHigh = B->T[j]->high;
				newT[newTcount]->var = B->T[j]->var;
				newT[newTcount]->level = B->T[j]->level;
				newT[newTcount]->deletion_flag = B->T[j]->deletion_flag;
				newTcount++; i++;
			}
			if (B->T[j]->level == level && level != numPrimeInput - 1) {
				//search with old low and high
				for (int k = 2; k < newTcount; k++) {
					 e = table_search(B->H, newT[k]);
					 if (e != NULL) {
						 for (int h = 2; h < B->size; h++) {
							 //update the lows and highs of this level
							 if (B->T[h]->level == level) {
								 if (B->T[h]->low == e->u) {
									 B->T[h]->newLow = k;
								 }
								 else if (B->T[h]->high == e->u) {
									 B->T[h]->newHigh = k;
								 }
							 }
						 }
					 }
				}
				newT[newTcount] = xmalloc(sizeof(struct node));
				memcpy(newT[newTcount], B->T[j], sizeof(struct node));
				newT[newTcount]->low = B->T[j]->low;
				newT[newTcount]->high = B->T[j]->high;
				newT[newTcount]->newLow = B->T[j]->newLow;
				newT[newTcount]->newHigh = B->T[j]->newHigh;
				newT[newTcount]->var = B->T[j]->var;
				newT[newTcount]->level = B->T[j]->level;
				newT[newTcount]->deletion_flag = B->T[j]->deletion_flag;
				newTcount++; i++;
				
			}
		}
		//printf("Debug print:\n");
		//printTTtable(newT, newTcount);pr
		level--;
	}
	for (i = 2; i < B->size; i++) {
		newT[i]->low = newT[i]->newLow;
		newT[i]->high = newT[i]->newHigh;
	}
	free(B->T);
	B->T = newT;
	table_free(B->H, &elem_free);
	B->H = table_new(BDD_HASHTABLE_SIZE,
		&elem_key, &node_equal, &hash);
	for (i = 2; i < B->size; i++) {
		/* enter into hash table */
		e = xmalloc(sizeof(struct elem));
		e->node = B->T[i];
		e->u = i;
		table_insert(B->H, e);
	}
}
int findBestBackwardPosition(t_blif_cubical_function *f, int varIndex, int varIndexBestPosition, bdd B) {
	int varIPosition, varJPosition;
	bool swapWasMade = false;
	int bestPosition;
	int newSize = B->size;
	varIPosition = getPositionOfVariable(varIndex);
	bestPosition = varIndexBestPosition;
	for (varJPosition = varIPosition - 1; varJPosition >= 0; varJPosition--) {
		if (varIPosition == 0) {
			return bestPosition;
		}
		swapWasMade = swap(f, varJPosition, B);

		int tempIndex = initVariableOrder[varIPosition];
		initVariableOrder[varIPosition] = initVariableOrder[varJPosition];
		initVariableOrder[varJPosition] = tempIndex;
		if (swapWasMade) {
			varIPosition--;
		}
		cleanTable(B);
		reordingTable(B);
		if (DEBUG) {
			printf("Debug - from backwards swap\n");
			printTtable(B);
		}
		if (newSize > B->size) {
			newSize = B->size;
			bestPosition = varJPosition;
		}
		if (BO->optSize > B->size) {
			BO->var = varIndex;
			BO->optSize = B->size;
			BO->optPosition = getPositionOfVariable(varIndex);
			for (int i = 0; i < input_count; i++) {
				BO->optOrder[i] = initVariableOrder[i];
			}
		}
	}
	return bestPosition;

}
int findBestForwardPosition(t_blif_cubical_function *f, int varIndex, int varIndexBestPosition, bdd B){
	int varIPosition, varJPosition;
	bool swapWasMade = false;
	int newSize = B->size;
	int bestPosition = varIndexBestPosition;
	varIPosition = getPositionOfVariable(varIndex);
	for (varJPosition = varIPosition + 1; varJPosition < input_count; varJPosition++) {
		
		swapWasMade = swap(f, varIPosition, B);

		int tempIndex = initVariableOrder[varIPosition];
		initVariableOrder[varIPosition] = initVariableOrder[varJPosition];
		initVariableOrder[varJPosition] = tempIndex;

		if (swapWasMade) {
			varIPosition++;
		}
		cleanTable(B);
		//printf("before reordering\n");
		//printTtable(B);
		//printf("Debug - forward swap- varI %d -> varJ %d\n", varIPosition - 1, varJPosition);
		reordingTable(B);
		//printTtable(B);
		if (newSize > B->size) { 
			newSize = B->size; 
			bestPosition = varJPosition;

		}
		if (BO->optSize > B->size) {
			BO->var = varIndex;
			BO->optSize = B->size;
			BO->optPosition = getPositionOfVariable(varIndex);
			for (int i = 0; i < input_count; i++) {
				BO->optOrder[i] = initVariableOrder[i];
			}
		}
	}
	
	return bestPosition;
}

int findBestPositionForVariable(int varIndex, t_blif_cubical_function *f, bdd B) {
	int varIndexPosition = getPositionOfVariable(varIndex);
	int varIndexBestPosition = varIndexPosition;

	
		varIndexBestPosition = findBestForwardPosition(f, varIndex, varIndexBestPosition, B);
		//printf("Variable order after forward swap:\n");
		//for (int i = 0; i < input_count; i++) {
		//	printf("%d: x%d\n", i + 1, initVariableOrder[i]);
		//}
		varIndexBestPosition = findBestBackwardPosition(f, varIndex, varIndexBestPosition, B);
	 //   printf("Variable order after backward swap:\n");
		//for (int i = 0; i < input_count; i++) {
		//	printf("%d: x%d\n", i + 1, initVariableOrder[i]);
		//}
		return varIndexBestPosition;
}
void moveVariable(int varIndex, int varNewPosition, t_blif_cubical_function *f, bdd B) {
	int varIndexPosition = getPositionOfVariable(varIndex);
	bool swapWasMade = true;

	if (varIndexPosition < varNewPosition) {
		while (swapWasMade && varIndexPosition < varNewPosition) {
			swapWasMade = swap(f, varIndexPosition, B);
			//printf("from moveVariable:\n");
			cleanTable(B);
			reordingTable(B);
			//printTtable(B);
			int tempIndex = initVariableOrder[varIndexPosition + 1];
			initVariableOrder[varIndexPosition + 1] = initVariableOrder[varIndexPosition];
			initVariableOrder[varIndexPosition] = tempIndex;
			if (swapWasMade) varIndexPosition++;
		}
	}	else if (varIndexPosition > varNewPosition) {
		while (swapWasMade && varIndexPosition > varNewPosition) {
			swapWasMade = swap(f, varIndexPosition - 1, B);
			cleanTable(B);
			reordingTable(B);
			int tempIndex = initVariableOrder[varIndexPosition - 1];
			initVariableOrder[varIndexPosition - 1] = initVariableOrder[varIndexPosition];
			initVariableOrder[varIndexPosition] = tempIndex;
			if (swapWasMade) varIndexPosition--;
		}
	}
}
void sifting(t_blif_logic_circuit *circuit, int function_index, bdd B) {
	int sizeBefore, sizeAfter;
	elem e;
	int bestPosition;
	node* tempT;
	varibableOccurence = xcalloc(circuit->primary_input_count, sizeof(int));
	variableOrderDesc = xcalloc(circuit->primary_input_count, sizeof(int));
	//initVariableOrder = xcalloc(circuit->list_of_functions[function_index]->input_count, sizeof(int));
	initVariableOcurrence(circuit, function_index, B);
	//initVarOrder(circuit->list_of_functions[function_index]);
	initVariableOrderDesc();
	if (DEBUG) {
		for (int i = 0; i < circuit->primary_input_count; i++) {
			printf("x%d has %d occurence\n", i, varibableOccurence[i]);
		}
		printf("Variable order:\n");
		for (int i = 0; i < circuit->list_of_functions[function_index]->input_count; i++) {
			printf("%d: x%d\n", i + 1, initVariableOrder[i]);
		}
		printf("Variable order Desc:\n");
		for (int i = 0; i < circuit->primary_input_count; i++) {
			printf("%d: x%d has %d occurence\n", i, variableOrderDesc[i], varibableOccurence[variableOrderDesc[i]]);
		}
	}
	sizeBefore = B->size;
	printf("Before sifting, BDD size is %d\n", BO->optSize);
	printTtable(B);
	printf("The variable order is: ");
	for (int i = 0; i < circuit->list_of_functions[function_index]->input_count; i++) {
		printf("X%d ", BO->optOrder[i] + 1);
	}
	for (int i = 0; i < numPrimeInput; i++) {
		if (varibableOccurence[variableOrderDesc[i]] != 0) {
			if (DEBUG) {
				printf("Start searching for best position for variable x%d\n" , variableOrderDesc[i]);
			}
			//printf("Before sifting for x%d!\n", variableOrderDesc[i]);
			//printTtable(B);
			tempT = xcalloc(B->limit, sizeof(node));
			for (int j = 0; j < B->size; j++) {
				tempT[j] = xmalloc(sizeof(struct node));
				memcpy(tempT[j], B->T[j], sizeof(struct node));
				tempT[j]->low = B->T[j]->low;
				tempT[j]->high = B->T[j]->high;
				tempT[j]->newLow = B->T[j]->low;
				tempT[j]->newHigh = B->T[j]->high;
				tempT[j]->var = B->T[j]->var;
				tempT[j]->level = B->T[j]->level;
				tempT[j]->deletion_flag = B->T[j]->deletion_flag;
			}

			//bestPositionForEachVar[variableOrderDesc[i]] = findBestPositionForVariable(0, circuit->list_of_functions[function_index], B);
			bestPosition = findBestPositionForVariable(variableOrderDesc[i], circuit->list_of_functions[function_index], B);
			//printf("Best Position for x%d is at %d!\n", variableOrderDesc[i], bestPosition);
			//printTtable(B);
			free(B->T);
			B->T = tempT;
			//printTTtable(tempT, B->size);
			B->size = sizeBefore;
			table_free(B->H, &elem_free);
			B->H = table_new(BDD_HASHTABLE_SIZE,
				&elem_key, &node_equal, &hash);
			for (int k = 2; k < B->size; k++) {
				/* enter into hash table */
				e = xmalloc(sizeof(struct elem));
				e->node = B->T[k];
				e->u = k;
				table_insert(B->H, e);
			}
			initVarOrder(circuit->list_of_functions[function_index]);
			if (false) {
				printf("Variable order:\n");
				for (int j = 0; j < circuit->list_of_functions[function_index]->input_count; j++) {
					printf("%d: x%d\n", j + 1, initVariableOrder[j]);
				}
				printf("debug print \n");
				printTtable(B);
				printf("TempT!\n");
				printTTtable(tempT, sizeBefore);
			}
		}
	}

	moveVariable(BO->var, BO->optPosition, circuit->list_of_functions[function_index], B);

}

/**********************************************************************/
/*** MAIN FUNCTION ****************************************************/
/**********************************************************************/

int main(int argc, char* argv[])
{
	clock_t time1, time2;
	time1 = clock();
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
		BO = xmalloc(sizeof(struct bestOptions));
		B = bdd_init(circuit->primary_input_count);
		numPrimeInput = circuit->primary_input_count;
		for (index = 0; index < circuit->function_count; index++)
		{
			t_blif_cubical_function *function = circuit->list_of_functions[index];
			//curFunc = (t_blif_cubical_function *)malloc(sizeof(t_blif_cubical_function));
			//curFunc->set_of_cubes = (t_blif_cube*)malloc(input_count * sizeof(t_blif_cube));
			BO->optOrder = xcalloc(function->input_count, sizeof(int));
			for (int i = 0; i < function->input_count; i++) {
				BO->optOrder[i] = function->inputs[i]->data.index;
			}
			initVariableOrder = xcalloc(function->input_count, sizeof(int));
			input_count = function->input_count;
			
			//B[index] = bdd_new(input_count);
			curFunc = function;
			//memcpy(curFunc->set_of_cubes, function->set_of_cubes, input_count * sizeof(t_blif_cube));

			print_function_info(function, index);
			
			//initH();
			printf("\n");
			//build(B[index], function->set_of_cubes,
			build(B, function->set_of_cubes,
				function->cube_count, 0, function->value);
			//reprintTtable(B);
			initVarOrder(circuit->list_of_functions[index]);
			//bestSizeForEachVar = xcalloc(numPrimeInput, sizeof(int));
			BO->var = function->inputs[0]->data.index;
			BO->optSize = B->size;
			BO->optPosition = getPositionOfVariable(BO->var);
			sifting(circuit, index, B);
			printf("\n\nSUMMARY: \nOne of the Best Options: X%d at No.%d position, and BDD size is %d\n", BO->var + 1, BO->optPosition + 1, BO->optSize);
			printf("One of the Best Orders: ");
			for (int i = 0; i < function->input_count; i++) {
				printf("X%d ", BO->optOrder[i] + 1);
			}
			printf("\n");
			printf("After Sifting\n");
			printTtable(B);
		}
		time2 = clock();
		DeleteBLIFCircuit(blif_circuit);
		printf("The time for sifting is %f ms\n", ((double)time2 - (double)time1) * 1000 / CLOCKS_PER_SEC);
		printf("Done.\r\n");
	}
	else
	{
		printf("Error reading BLIF file. Terminating.\n");
	}
	return 0;
}

