#include <malloc.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include "bdd.h"

void elem_free(ht_elem e) {
	free((elem)e);
}
bool node_equal(ht_key k1, ht_key k2) {
	node a1 = (node)k1;		
	node a2 = (node)k2;		
	return a1->var == a2->var
		&& a1->low == a2->low
		&& a1->high == a2->high;
}
ht_key elem_key(ht_elem e)
{
	return ((elem)e)->node;
}
int hash(ht_key k, int m)
{
	node a = (node)k;		
	unsigned int x = 1234;
	unsigned int y = 5678;	
	unsigned int r = 0xdeadbeef;	
	unsigned int h = (unsigned)a->var;
	r = r*x + y;		
	h = r*h + (unsigned)a->low;
	r = r*x + y;		
	h = r*h + (unsigned)a->high;
	h = h % (unsigned)m;		
	return (int)h;
}

void bdd_free(bdd B) {
	for (int i = 0; i < B->size; i++)
		free(B->T[i]);		
	free(B->T);			
	table_free(B->H, &elem_free);	
	free(B);			
}
bdd bdd_init(int num_var)
{
	bdd B = xmalloc(sizeof(struct bdd));
	B->num_vars = num_var;
	B->limit = BDD_LIMIT;
	B->size = 2;
	{
		node* T = xcalloc(B->limit, sizeof(node));
		node zero = xmalloc(sizeof(struct node));
		node one = xmalloc(sizeof(struct node));
		zero->var = num_var + 1;
		one->var = num_var + 1;
		zero->low = -1;
		zero->high = -1;
		zero->deletion_flag = false;
		one->deletion_flag = false;
		one->low = -1;
		one->high = -1;
		zero->level = num_var + 1;
		one->level = num_var + 1;
		T[0] = zero;
		T[1] = one;
		B->T = T;
	}
	B->H = table_new(BDD_HASHTABLE_SIZE,
		&elem_key, &node_equal, &hash);
	return B;
}

bdd_node MK(bdd B, int var, int low, int high) {
	node a; elem e; int u;
	if (low == high) return low;	
	a = xmalloc(sizeof(struct node));
	a->var = var;
	a->low = low;
	a->high = high;
	a->deletion_flag = false;
	a->level = var;
	e = table_search(B->H, a);
	if (e != NULL) {
		free(a);			
		return e->u;		
	}
	u = B->size;
	B->T[u] = a;
	B->size++;
	e = xmalloc(sizeof(struct elem));
	e->node = a;
	e->u = u;
	table_insert(B->H, e);
	return u;
}


/**********************************************************************/
/*** UTILITIES ********************************************************/
/**********************************************************************/

void print_cube(int input_count, t_blif_cube *cube)
{
	int input_index;

	for (input_index = 0; input_index < input_count; input_index++)
		switch (read_cube_variable(cube->signal_status, input_index))
		{
		case LITERAL_DC:
			printf("-");
			break;
		case LITERAL_1:
			printf("1");
			break;
		case LITERAL_0:
			printf("0");
			break;
		case LITERAL_MARKER:
			printf("A");
			break;
		}
	printf("\n");
}



void printTtable(bdd B) {
	printf("===========================================\n");
	printf("|u |\tvar |\tlow |\thigh |\tflag |\tlevel |\n");
	for (int i = 0; i < B->size; i++) {
		if (B->T[i]->low != -1)
			printf("|%d |\t  %d |\t  %d |\t  %d  |\t  %d  |\t  %d  |\n", i, B->T[i]->var , B->T[i]->low, B->T[i]->high,B->T[i]->deletion_flag, B->T[i]->level);
		else 
			printf("|%d |\t  %d |\t    |\t     |\t  %d  |\t  %d  |\n", i, B->T[i]->var , B->T[i]->deletion_flag, B->T[i]->level);
	}
	printf("============================================\n");
}
void printTTtable(node *T, int size) {
	printf("============================================\n");
	printf("|u |\tvar |\tlow |\thigh |\tflag |\tlevel |\n");
	for (int i = 0; i < size; i++) {
		if (T[i]->low != -1)
			printf("|%d |\t  %d |\t  %d |\t  %d  |\t  %d  |\t  %d  |\n", i, T[i]->var, T[i]->low, T[i]->high, T[i]->deletion_flag, T[i]->level);
		else
			printf("|%d |\t  %d |\t    |\t     |\t  %d  |\t  %d  |\n", i, T[i]->var , T[i]->deletion_flag, T[i]->level);
	}
	printf("============================================\n");
}