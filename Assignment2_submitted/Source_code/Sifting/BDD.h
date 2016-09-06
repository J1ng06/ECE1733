#include "cubical_function_representation.h"
#include "hashtable.h"
#include "xalloc.h"

#define BDD_LIMIT 1024
#define BDD_HASHTABLE_SIZE 1024
#define APPLY_HASHTABLE_SIZE 1024
typedef int bdd_node;		
typedef struct node* node;	
struct node {
	int var;			
	bdd_node low;		
	bdd_node high;	
	bdd_node newLow;
	bdd_node newHigh;
	bool deletion_flag;
	int level;
};

typedef struct bdd* bdd;
struct bdd {
	int num_vars;			
	int limit;			
	int size;			
	node* T;			
	table H;			
};
void elem_free(ht_elem e);
bool node_equal(ht_key k1, ht_key k2);
ht_key elem_key(ht_elem e);
int hash(ht_key k, int m);
bdd bdd_init(int k);

typedef struct elem* elem;
struct elem {
	node node;			
	bdd_node u;			
};

bdd_node MK(bdd B, int var, bdd_node low, bdd_node high);
void printTtable(bdd B);
