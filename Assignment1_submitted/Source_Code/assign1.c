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
//#include <curses.h> 
#include <stdio.h>
#include <assert.h>
#include <malloc.h>
#include <string.h>
#include "common_types.h"
#include "blif_common.h"
#include "cubical_function_representation.h"


#define FULLY_COVERED  -1
#define PARTIALLY_COVERED  1
#define NOT_COVERED_AT_ALL  0
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
/*** MINIMIZATION CODE ************************************************/
/**********************************************************************/

void copy_cube_value(t_blif_cube *result_cube, t_blif_cube *source) {
	//int i;
	//result_cube->data_size = &source->data_size;
	//result_cube->is_DC = &source->is_DC;
	//for (i = 0; i < 4; i++){
	//write_cube_variable(result_cube->signal_status, i, read_cube_variable(source->signal_status, i));
	//result_cube->signal_status[i] = source->signal_status[i];}
	//}
	result_cube = source;
}

// combine 2 covers and remove duplicate cubes
int combine_cover(t_blif_cube **result_cover, t_blif_cube **C, t_blif_cube **G, int input_count, int C_cube_count, int G_cube_count) {
	int i, j, k, result_cube_count, same_literal_count;
	int duplicate_flag;
	j = 0;
	same_literal_count = 0;
	result_cube_count = 0;
	duplicate_flag = 0;
	memcpy(result_cover, C, C_cube_count * sizeof(t_blif_cube));
	result_cube_count = C_cube_count;

	i = 0;
	while (i < G_cube_count) {
		while (j < result_cube_count) {
			for (k = 0; k < input_count; k++) {
				if (read_cube_variable(result_cover[j]->signal_status, k) == read_cube_variable(G[i]->signal_status, k)) {
					//					printf("result_cube: G_cube -- %d: %d\n",read_cube_variable(result_cover[j]->signal_status, k), read_cube_variable(G[i]->signal_status,k));
					same_literal_count++;
				}
			}
			if (same_literal_count == input_count) {
				duplicate_flag = 1;
			}
			same_literal_count = 0;
			j++;
		}
		if (duplicate_flag == 0) {
			//copy_cube_value(&result_cover[result_cube_count++], &G[i]);
			result_cover[result_cube_count] = G[i];
			result_cube_count++;
		}
		duplicate_flag = 0;
		i++;
		j = 0;
	}
	return result_cube_count;
}
void sharpRule3(t_blif_cube **sharp_result_cover, t_blif_cube *result_cube, int i) {
	t_blif_cube *cube = (t_blif_cube*)malloc(sizeof(t_blif_cube));
	memcpy(cube, result_cube, sizeof(t_blif_cube));
	sharp_result_cover[i] = cube;
}
int sharp_operation(t_blif_cube **sharp_result_cover, t_blif_cube *A, t_blif_cube *B, t_blif_cube *result_cube, int input_count) {
	int i, result_cover_cube_count, count;
	char result[5];
	//t_blif_cube *result_cube = (t_blif_cube*)malloc(sizeof(t_blif_cube));
	count = 0;
	result[4] = '\0';
	result_cover_cube_count = 0;
	for (i = 0; i < input_count; i++) {
		if (read_cube_variable(A->signal_status, i) == 1 && read_cube_variable(B->signal_status, i) == 1) {
			result[i] = 'e';
			count++;
		}
		else if (read_cube_variable(A->signal_status, i) == 1 && read_cube_variable(B->signal_status, i) == 2) {
			result[i] = 'n';
			//copy_cube_value(sharp_result_cover[0], A);
			sharp_result_cover[0] = A;
			return 0;
		}
		else if (read_cube_variable(A->signal_status, i) == 1 && read_cube_variable(B->signal_status, i) == 3) {
			result[i] = 'e';
			count++;
		}
		else if (read_cube_variable(A->signal_status, i) == 2 && read_cube_variable(B->signal_status, i) == 1) {
			result[i] = 'n';
			sharp_result_cover[0] = A;
			//copy_cube_value(sharp_result_cover[0], A);
			//memcpy(sharp_result_cover[0], A, sizeof(t_blif_cube));
			return 0;
		}
		else if (read_cube_variable(A->signal_status, i) == 2 && read_cube_variable(B->signal_status, i) == 2) {
			result[i] = 'e';
			count++;
		}
		else if (read_cube_variable(A->signal_status, i) == 2 && read_cube_variable(B->signal_status, i) == 3) {
			result[i] = 'e';
			count++;
		}
		else if (read_cube_variable(A->signal_status, i) == 3 && read_cube_variable(B->signal_status, i) == 1) {
			result[i] = '1';
		}
		else if (read_cube_variable(A->signal_status, i) == 3 && read_cube_variable(B->signal_status, i) == 2) {
			result[i] = '0';
		}
		else if (read_cube_variable(A->signal_status, i) == 3 && read_cube_variable(B->signal_status, i) == 3) {
			result[i] = 'e';
			count++;
		}
	}
	if (count == input_count) {
		return -1;
	}
	else {
		for (i = 0; i < input_count; i++) {
			if (read_cube_variable(A->signal_status, i) == LITERAL_DC) {
				if (read_cube_variable(B->signal_status, i) != LITERAL_DC) {
					//copy_cube_value(result_cube, A);
					memcpy(result_cube, A, sizeof(t_blif_cube));
					//result_cube = A;
					//(sharp_result_cover[result_cover_cube_count], A);
					if (read_cube_variable(B->signal_status, i) == LITERAL_0) {
						write_cube_variable(result_cube->signal_status, i, LITERAL_1);
						//memcpy(sharp_result_cover[j], result_cube, sizeof(t_blif_cube)); 
						//copy_cube_value(sharp_result_cover[result_cover_cube_count], result_cube);
						//sharp_result_cover[result_cover_cube_count] = result_cube;
						sharpRule3(sharp_result_cover, result_cube, result_cover_cube_count);
						result_cover_cube_count++;
					}
					else {
						write_cube_variable(result_cube->signal_status, i, LITERAL_0);
						//memcpy(sharp_result_cover[j++], result_cube, sizeof(t_blif_cube)); 
						//sharp_result_cover[j++] = result_cube;
						//copy_cube_value(sharp_result_cover[result_cover_cube_count], result_cube);
						//sharp_result_cover[result_cover_cube_count] = result_cube;
						sharpRule3(sharp_result_cover, result_cube, result_cover_cube_count);
						result_cover_cube_count++;
					}
				}
			}
		}
	}
	return result_cover_cube_count;
}
void star_operation(t_blif_cubical_function *f_G_cover, t_blif_cube *A, t_blif_cube *B, int input_count) {

	int i;
	int count = 0;
	char null = 'n';
	t_blif_cube *result_cube = (t_blif_cube*)malloc(sizeof(t_blif_cube*));
	char* result = (char*)malloc((input_count + 1)*sizeof(char*));
	result[input_count] = '\0';

	// Star-Operations
	for (i = 0; i < input_count; i++) {
		if (read_cube_variable(A->signal_status, i) == 1 && read_cube_variable(B->signal_status, i) == 1) {
			result[i] = '0';
		}
		else if (read_cube_variable(A->signal_status, i) == 1 && read_cube_variable(B->signal_status, i) == 2) {
			result[i] = 'n';
		}
		else if (read_cube_variable(A->signal_status, i) == 1 && read_cube_variable(B->signal_status, i) == 3) {
			result[i] = '0';
		}
		else if (read_cube_variable(A->signal_status, i) == 2 && read_cube_variable(B->signal_status, i) == 1) {
			result[i] = 'n';
		}
		else if (read_cube_variable(A->signal_status, i) == 2 && read_cube_variable(B->signal_status, i) == 2) {
			result[i] = '1';
		}
		else if (read_cube_variable(A->signal_status, i) == 2 && read_cube_variable(B->signal_status, i) == 3) {
			result[i] = '1';
		}
		else if (read_cube_variable(A->signal_status, i) == 3 && read_cube_variable(B->signal_status, i) == 1) {
			result[i] = '0';
		}
		else if (read_cube_variable(A->signal_status, i) == 3 && read_cube_variable(B->signal_status, i) == 2) {
			result[i] = '1';
		}
		else if (read_cube_variable(A->signal_status, i) == 3 && read_cube_variable(B->signal_status, i) == 3) {
			result[i] = 'x';
		}
	}

	// checking for not useful cubes (contains more than 1 'n')
	for (i = 0; result[i] != '\0'; ++i)
	{
		if (result[i] == null)
			count++;
	}
	//printf("%s\n", result);
	// writing a new cube
	if (count <= 1) {
		for (i = 0; i < input_count; i++) {
			switch (result[i])
			{
			case '0': write_cube_variable(result_cube->signal_status, i, LITERAL_0); break;
			case '1': write_cube_variable(result_cube->signal_status, i, LITERAL_1); break;
			case 'n': write_cube_variable(result_cube->signal_status, i, LITERAL_DC); break;
			case 'x': write_cube_variable(result_cube->signal_status, i, LITERAL_DC); break;
			default:
				break;
			}
		}
		f_G_cover->set_of_cubes[f_G_cover->cube_count] = result_cube;
		f_G_cover->cube_count++;
		//copy_cube_value(&result_cover[cubes_count], &result_cube);
		//memcpy(result_cover[cubes_count], result_cube, sizeof(t_blif_cube));
	}
	//printf("%s  %d\n", result, sizeof(result));

}

void enumerateAllMinterms(t_blif_cubical_function *f_all_minterms, t_blif_cubical_function *f) {
	int i, j, k, l, size, count, index, temp_index, temp_value, temp_size;
	char b[64];
	size = 0;
	temp_size = 0;
	// count how many dont care inputs in all the {ON set} cubes
	// and count how many minterms in total 
	for (i = 0; i < f->cube_count; i++) {
		temp_size = 0;
		for (j = 0; j < f->input_count; j++) {
			if (read_cube_variable(f->set_of_cubes[i]->signal_status, j) == LITERAL_DC
				&& f->set_of_cubes[i]->is_DC == T_FALSE) {
				temp_size++;
			}
		}
		if (temp_size == 0 && f->set_of_cubes[i]->is_DC == T_FALSE) {
			size++;
		}
		else if (temp_size != 0 && f->set_of_cubes[i]->is_DC == T_FALSE)
		{
			size += 1 << temp_size;
		}
	}

	index = 0;
	count = 0;
	for (i = 0; i < f->cube_count; i++) {
		for (j = 0; j < f->input_count; j++) {
			if (read_cube_variable(f->set_of_cubes[i]->signal_status, j) == LITERAL_DC
				&& f->set_of_cubes[i]->is_DC == T_FALSE) {
				count++;
			}
		}
		if (count > 0) {
			for (k = 0; k < 1 << count; k++) {
				//for (j = 0; j < f->input_count; j++) {
				//	write_cube_variable(f_all_minterms->set_of_cubes[index]->signal_status, j, read_cube_variable(f->set_of_cubes[i]->signal_status, j));
				//	
				//}
				//copy_cube_value(f_all_minterms->set_of_cubes[index], f->set_of_cubes[i]);
				f_all_minterms->set_of_cubes[index] = f->set_of_cubes[i];
				index++;
			}
			temp_index = index - (1 << count);
			for (k = 0; k < 1 << count; k++) {
				temp_value = k;
				for (j = f->input_count - 1; j >= 0; j--) {
					if (read_cube_variable(f_all_minterms->set_of_cubes[temp_index]->signal_status, j) == LITERAL_DC) {
						write_cube_variable(f_all_minterms->set_of_cubes[temp_index]->signal_status, j, (temp_value % 2 == 0) ? LITERAL_0 : LITERAL_1);
						temp_value = temp_value / 2;
					}
				}
				temp_index++;
			}
		}
		else {
			// for cubes that can be treated as minterms
			for (j = 0; j < f->input_count; j++) {
				write_cube_variable(f_all_minterms->set_of_cubes[index]->signal_status, j, read_cube_variable(f->set_of_cubes[i]->signal_status, j));
			}
			index++;
		}
		count = 0;
	}
}

int is_two_covers_equal(t_blif_cubical_function *prev_cover, t_blif_cubical_function *next_cover, int input_count) {
	int i, j;
	if (prev_cover->cube_count != next_cover->cube_count) {
		return 0;
	}

	else {
		for (i = 0; i < prev_cover->cube_count; ++i) {
			for (j = 0; j < input_count; j++) {
				if (read_cube_variable(prev_cover->set_of_cubes[i]->signal_status, j) !=
					read_cube_variable(next_cover->set_of_cubes[i]->signal_status, j))
				{
					return 0;
				}
			}
		}
	}
	return 1;
}

void get_non_EPIs(t_blif_cube *Non_EPI, t_blif_cubical_function *f_non_EPI) {
	t_blif_cube *non_EPI = (t_blif_cube*)malloc(sizeof(t_blif_cube));
	memcpy(non_EPI, Non_EPI, sizeof(t_blif_cube));
	f_non_EPI->set_of_cubes[f_non_EPI->cube_count++] = non_EPI;
}
void get_EPIs(t_blif_cube *_EPI, t_blif_cubical_function *f_EPI) {
	t_blif_cube *EPI = (t_blif_cube*)malloc(sizeof(t_blif_cube));
	memcpy(EPI, _EPI, sizeof(t_blif_cube));
	f_EPI->set_of_cubes[f_EPI->cube_count++] = EPI;
}

void remove_this_cube(t_blif_cubical_function *f, int i) {
	for (int w = i; w < f->cube_count - 1; w++) {
		f->set_of_cubes[w] = f->set_of_cubes[w + 1];
	}
	// now the cover set is one cube less than original cover set
	f->cube_count = f->cube_count - 1;
}
int combine_cover_EPI(t_blif_cube **result_cover, t_blif_cube **ON, t_blif_cube **DC, int input_count, int ON_cube_count, int DC_cube_count) {
	int i, j, k, result_cube_count, same_literal_count;
	int duplicate_flag;
	j = 0;
	result_cube_count = 0;
	memcpy(result_cover, ON, ON_cube_count* sizeof(t_blif_cube));
	result_cube_count = ON_cube_count;

	i = 0;
	while (i < DC_cube_count) {
		//copy_cube_value(&result_cover[result_cube_count++], &G[i]);
		result_cover[result_cube_count] = DC[i];
		result_cube_count++;
		i++;
	}
	return result_cube_count;
}
int sharp_operation_EPI(t_blif_cube **sharp_result_cover, t_blif_cube *A, t_blif_cube *B, t_blif_cube *result_cube, int input_count) {
	int i, result_cover_cube_count, count;
	char result[5];
	//t_blif_cube *result_cube = (t_blif_cube*)malloc(sizeof(t_blif_cube));
	count = 0;
	result[4] = '\0';
	result_cover_cube_count = 0;
	for (i = 0; i < input_count; i++) {
		if (read_cube_variable(A->signal_status, i) == 1 && read_cube_variable(B->signal_status, i) == 1) {
			result[i] = 'e';
			count++;
		}
		else if (read_cube_variable(A->signal_status, i) == 1 && read_cube_variable(B->signal_status, i) == 2) {
			result[i] = 'n';
			//copy_cube_value(sharp_result_cover[0], A);
			sharp_result_cover[0] = A;
			return 0;
		}
		else if (read_cube_variable(A->signal_status, i) == 1 && read_cube_variable(B->signal_status, i) == 3) {
			result[i] = 'e';
			count++;
		}
		else if (read_cube_variable(A->signal_status, i) == 2 && read_cube_variable(B->signal_status, i) == 1) {
			result[i] = 'n';
			sharp_result_cover[0] = A;
			//copy_cube_value(sharp_result_cover[0], A);
			//memcpy(sharp_result_cover[0], A, sizeof(t_blif_cube));
			return 0;
		}
		else if (read_cube_variable(A->signal_status, i) == 2 && read_cube_variable(B->signal_status, i) == 2) {
			result[i] = 'e';
			count++;
		}
		else if (read_cube_variable(A->signal_status, i) == 2 && read_cube_variable(B->signal_status, i) == 3) {
			result[i] = 'e';
			count++;
		}
		else if (read_cube_variable(A->signal_status, i) == 3 && read_cube_variable(B->signal_status, i) == 1) {
			result[i] = '1';
		}
		else if (read_cube_variable(A->signal_status, i) == 3 && read_cube_variable(B->signal_status, i) == 2) {
			result[i] = '0';
		}
		else if (read_cube_variable(A->signal_status, i) == 3 && read_cube_variable(B->signal_status, i) == 3) {
			result[i] = 'e';
			count++;
		}
	}
	if (count == input_count) {
		return -1;
	}
	else {
		for (i = 0; i < input_count; i++) {
			if (read_cube_variable(A->signal_status, i) == LITERAL_DC) {
				if (read_cube_variable(B->signal_status, i) != LITERAL_DC) {
					//copy_cube_value(result_cube, A);
					memcpy(result_cube, A, sizeof(t_blif_cube));
					//result_cube = A;
					//(sharp_result_cover[result_cover_cube_count], A);
					if (read_cube_variable(B->signal_status, i) == LITERAL_0) {
						write_cube_variable(result_cube->signal_status, i, LITERAL_1);
						//memcpy(sharp_result_cover[j], result_cube, sizeof(t_blif_cube)); 
						//copy_cube_value(sharp_result_cover[result_cover_cube_count], result_cube);
						//sharp_result_cover[result_cover_cube_count] = result_cube;
						sharpRule3(sharp_result_cover, result_cube, result_cover_cube_count);
						result_cover_cube_count++;
					}
					else {
						write_cube_variable(result_cube->signal_status, i, LITERAL_0);
						//memcpy(sharp_result_cover[j++], result_cube, sizeof(t_blif_cube)); 
						//sharp_result_cover[j++] = result_cube;
						//copy_cube_value(sharp_result_cover[result_cover_cube_count], result_cube);
						//sharp_result_cover[result_cover_cube_count] = result_cube;
						sharpRule3(sharp_result_cover, result_cube, result_cover_cube_count);
						result_cover_cube_count++;
					}
				}
			}
		}
	}
	return result_cover_cube_count;
}

void find_EPI_cover_set(t_blif_cubical_function *f_ON, t_blif_cubical_function *f_DC, t_blif_cubical_function *f_EPI, t_blif_cubical_function *f_non_EPI) {
	int i, j, k, h, w, loop, num_of_non_epi, number, result_cover_cube_count, sharp_loop, list_cube_index, counter, prev_sharp_result_cube_count, EPI_flag;
	int index[100];
	t_blif_cube **sharp_result_cover_EPI = (t_blif_cube **)malloc(f_ON->input_count*sizeof(t_blif_cube*));
	t_blif_cube **sharp_result_cover_list = (t_blif_cube **)malloc(1000 * sizeof(t_blif_cube*));
	t_blif_cube *sharp_result_cube_EPI = (t_blif_cube*)malloc(sizeof(t_blif_cube));
	t_blif_cubical_function *f_all_cubes = (t_blif_cubical_function*)malloc(sizeof(t_blif_cubical_function));
	f_all_cubes->set_of_cubes = (t_blif_cube*)malloc(1000 * sizeof(t_blif_cube));
	f_all_cubes->cube_count = combine_cover_EPI(f_all_cubes->set_of_cubes, f_ON->set_of_cubes, f_DC->set_of_cubes, f_ON->input_count, f_ON->cube_count, f_DC->cube_count);
	num_of_non_epi = 0;
	list_cube_index = 0;
	counter = 0;
	for (i = 0; i < f_ON->cube_count; i++) {
		loop = 0;
		for (j = 0; j < f_all_cubes->cube_count; j++) {
			if (i != j) {
				if (loop == 0) {
					loop++;
					result_cover_cube_count = sharp_operation_EPI(sharp_result_cover_EPI, f_ON->set_of_cubes[i], f_all_cubes->set_of_cubes[j], sharp_result_cube_EPI, f_ON->input_count);
					if (result_cover_cube_count == -1) {
						// fully covered, C = nothing yields this PI is not EPI
						//TODO
						index[num_of_non_epi] = i;
						num_of_non_epi++;
						break;
					}
					else {// result_cover_cube_count >= 0
						if (result_cover_cube_count == 0) {
							result_cover_cube_count = 1;
						}
						list_cube_index = 0;
						//printf("****** i = %d, j = %d ****\n", i, j);
						for (w = 0; w < result_cover_cube_count; w++) {
							// apend all sharp result cubes into list
							sharp_result_cover_list[list_cube_index] = sharp_result_cover_EPI[w];
							//for (k = 0; k < f_ON->input_count; k++) {
							//	printf("%d ", read_cube_variable(sharp_result_cover_list[w]->signal_status, k));
							//}
							//printf("\n");
							list_cube_index++;
						}
						list_cube_index--;
					}
				}
				// all other sharp operations other than first one
				else {
					// use resulted cubes to sharp next one
					counter = list_cube_index + 1;
					prev_sharp_result_cube_count = (list_cube_index + 1);
					for (h = 0; h < prev_sharp_result_cube_count; h++) {
						result_cover_cube_count = sharp_operation_EPI(sharp_result_cover_EPI, sharp_result_cover_list[h], f_all_cubes->set_of_cubes[j], sharp_result_cube_EPI, f_ON->input_count);
						if (result_cover_cube_count == -1) {
							// decrease #counter by one
							counter--;
						}
						else {// result_cover_cube_count >= 0

							if (result_cover_cube_count == 0) {
								// C = A return 0, cube count equals 1
								result_cover_cube_count = 1;
							}

							// reset 
							list_cube_index = 0;

							//printf("****** i = %d, j = %d ****\n", i, j);
							for (w = 0; w < result_cover_cube_count; w++) {
								// apend all sharp result cubes into list
								sharp_result_cover_list[list_cube_index] = sharp_result_cover_EPI[w];
								//for (k = 0; k < f_ON->input_count; k++) {
								//	printf("%d ", read_cube_variable(sharp_result_cover_list[w]->signal_status, k));
								//}
								//printf("\n");
								list_cube_index++;
							}
							list_cube_index--;
						}
					}
					// if #counter equals 0, then this is Non-EPI
					//TODO
					if (counter == 0) {
						// fully covered, C = nothing yields this PI is not EPI
						//TODO
						index[num_of_non_epi] = i;
						num_of_non_epi++;
						break;
					}
				}
			}
		}
	}
	f_EPI->cube_count = 0;
	f_non_EPI->cube_count = 0;
	// get all the non-epi cubes, and remove them
	for (k = 0; k < num_of_non_epi; k++) {
		get_non_EPIs(f_ON->set_of_cubes[index[k]], f_non_EPI);
	}
	EPI_flag = 1;
	for (i = 0; i < f_ON->cube_count; i++) {
		for (j = 0; j < num_of_non_epi; j++) {
			if (i == index[j]) {
				EPI_flag = 0;
			}
		}
		if (EPI_flag == 1) {
			get_EPIs(f_ON->set_of_cubes[i], f_EPI);
		}
		EPI_flag = 1;
	}


}
void append_cubes(t_blif_cubical_function *f_result_cover, t_blif_cube *new_cube) {
	t_blif_cube *New_cube = (t_blif_cube*)malloc(sizeof(t_blif_cube));
	memcpy(New_cube, new_cube, sizeof(t_blif_cube));
	f_result_cover->set_of_cubes[f_result_cover->cube_count++] = New_cube;
}
int get_not_covered_ON_set(t_blif_cubical_function *f, t_blif_cubical_function *f_EPI, t_blif_cubical_function *f_not_covered_ON) {
	int i, j, w, h, k, sharp_cube_count, loop, list_cube_index, prev_sharp_result_cube_count, counter;
	t_blif_cube **sharp_result_cover = (t_blif_cube **)malloc(f->input_count * sizeof(t_blif_cube*));

	t_blif_cube **sharp_result_cover_list = (t_blif_cube **)malloc(1000 * sizeof(t_blif_cube*));
	t_blif_cube *sharp_result_cube = (t_blif_cube*)malloc(sizeof(t_blif_cube));
	for (i = 0; i < f->cube_count; i++) {
		if (f->set_of_cubes[i]->is_DC == T_FALSE) {
			loop = 0;
			for (j = 0; j < f_EPI->cube_count; j++) {
					if (loop == 0) {
						loop++;
						sharp_cube_count = sharp_operation_EPI(sharp_result_cover, f->set_of_cubes[i], f_EPI->set_of_cubes[j], sharp_result_cube, f->input_count);
						if (sharp_cube_count == -1) {
							// fully covered, C = nothing yields this PI is not EPI
							//TODO
							break;
						}
						else {// result_cover_cube_count >= 0
							if (sharp_cube_count == 0) {
								sharp_cube_count = 1;
							}
							list_cube_index = 0;
							//printf("****** i = %d, j = %d ****\n", i, j);
							for (w = 0; w < sharp_cube_count; w++) {
								// apend all sharp result cubes into list
								sharp_result_cover_list[list_cube_index] = sharp_result_cover[w];
								//for (k = 0; k < f->input_count; k++) {
									//printf("%d ", read_cube_variable(sharp_result_cover_list[w]->signal_status, k));
								//}
								//printf("\n");
								if (j == f_EPI->cube_count - 1) {
									append_cubes(f_not_covered_ON, sharp_result_cover[w]);
								}
								list_cube_index++;
							}
							list_cube_index--;
						}
					}
					// all other sharp operations other than first one
					else {
						// use resulted cubes to sharp next one
						counter = list_cube_index + 1;
						prev_sharp_result_cube_count = (list_cube_index + 1);
						for (h = 0; h < prev_sharp_result_cube_count; h++) {
							sharp_cube_count = sharp_operation_EPI(sharp_result_cover, sharp_result_cover_list[h], f_EPI->set_of_cubes[j], sharp_result_cube, f->input_count);
							if (sharp_cube_count == -1) {
								// decrease #counter by one
								counter--;
							}
							else {// result_cover_cube_count >= 0

								if (sharp_cube_count == 0) {
									// C = A return 0, cube count equals 1
									sharp_cube_count = 1;
								}

								// reset 
								list_cube_index = 0;

								//printf("****** i = %d, j = %d ****\n", i, j);
								for (w = 0; w < sharp_cube_count; w++) {
									// apend all sharp result cubes into list
									sharp_result_cover_list[list_cube_index] = sharp_result_cover[w];
									if (j == f_EPI->cube_count - 1) {
										append_cubes(f_not_covered_ON, sharp_result_cover[w]);
									}
									//for (k = 0; k < f->input_count; k++) {
									//	printf("%d ", read_cube_variable(sharp_result_cover_list[w]->signal_status, k));
									//}
									//printf("\n");
									
									list_cube_index++;
								}
								list_cube_index--;
							}
						}
						// if #counter equals 0, then this is Non-EPI
						//TODO
						if (counter == 0) {
							// fully covered, C = nothing yields this PI is not EPI
							//TODO
							break;
						}
					}
				
			}
		}
	}
	return f_not_covered_ON->cube_count;
}
int get_min_cost_cube(int cost[1000], int cube_count) {
	int i;
	int min = 0;
	for (i = 0; i < cube_count; i++) {
		if (cost[i] < cost[min]) {
			min = i;
		}
	}
	return min;
}

int get_remaining_cover_set(t_blif_cubical_function *f_not_covered_ON, t_blif_cube *non_EPI) {
	int i, w, loop, num_fully_covered_cubes, num_not_fully_covered_cubes, sharp_result_cube_count;

	t_blif_cube **sharp_result_cover = (t_blif_cube **)malloc(f_not_covered_ON->input_count * sizeof(t_blif_cube*));
	t_blif_cube *sharp_result_cube = (t_blif_cube*)malloc(sizeof(t_blif_cube));
	t_blif_cubical_function *temp_f_not_covered_ON = (t_blif_cubical_function *)malloc(sizeof(t_blif_cubical_function));
	temp_f_not_covered_ON->set_of_cubes = (t_blif_cube **)malloc(1000 * sizeof(t_blif_cube*));
	// make a temp function
	memcpy(temp_f_not_covered_ON->set_of_cubes, f_not_covered_ON->set_of_cubes, f_not_covered_ON->cube_count*sizeof(t_blif_cube));
	temp_f_not_covered_ON->cube_count = f_not_covered_ON->cube_count;
	num_fully_covered_cubes = 0;
	num_not_fully_covered_cubes = 0;
	f_not_covered_ON->cube_count = 0;
	for (i = 0; i < temp_f_not_covered_ON->cube_count; i++) {
			sharp_result_cube_count = sharp_operation_EPI(sharp_result_cover, temp_f_not_covered_ON->set_of_cubes[i], non_EPI, sharp_result_cube, f_not_covered_ON->input_count);
			if (sharp_result_cube_count == FULLY_COVERED) {
				// this PI fully covered current not_covered_ON_cube
				// set counter increment
				num_fully_covered_cubes++;
			}
			else if (sharp_result_cube_count == NOT_COVERED_AT_ALL) {
				// this PI does not cover the current not_covered_ON_cube at all

				// need to add the uncovered cube onto our new cover set

					// append cube(s) into not_covered_ON_set (updating remaining cover set)
				append_cubes(f_not_covered_ON, sharp_result_cover[0]);
				

				// keep a record of it (NOT_COVERED_AT_ALL)
				num_not_fully_covered_cubes++;
			}
			else {
				// sharp_cube_count >= 1
				// patiallly is not covered by this PI
				// need to add the uncovered part(s) onto our new cover set
				for (w = 0; w < sharp_result_cube_count; w++) {
					// append cube(s) into not_covered_ON_set (updating remaining cover set)
					append_cubes(f_not_covered_ON, sharp_result_cover[w]);
				}
			}
		}
	
		if (num_fully_covered_cubes == temp_f_not_covered_ON->cube_count) {
			// this PI fully covered the whole not_covered_cover_set
			// return FULLY_COVERED
			return FULLY_COVERED;
		}

		if (num_not_fully_covered_cubes == temp_f_not_covered_ON->cube_count) {
			// this PI does not cover any of the not_covered_cover_set cubes
			// return NOT_COVERED_AT_ALL
			return NOT_COVERED_AT_ALL;
		}

		// return PARTIALLY_COVERED
		return PARTIALLY_COVERED;
	}
void simplify_function(t_blif_cubical_function *f)
/* This function simplifies the function f. The minimized set of cubes is
* returned though a field in the input structure called set_of_cubes.
* The number of cubes is stored in the field cube_count.
*/

{
	/* PUT YOUR CODE HERE */
	// G_cover for storing the result from star op
	// union_cover for combining G_cover and C_cover with no duplicate cubes
	// sharp_result_cover for partial result cover from sharp op
	int i, j, k, h, input_count, union_cube_count, redundant_flag, iteration, num_of_cube_no_DC, num_of_cube_DC, min_cost_cube_index, remaining_set_flag;
	int num_A, lowest_cost_cube_count, second_lowest_cost_cube_count, cost_for_optimization;
	int sharp_result_cube_count, min_cost, second_min_cost, optimization_counter, tracking_counter;
	
	t_blif_cube **C_temp_cover;
	// this one would be problemetic ***************************
	t_blif_cube *sharp_result_cube = (t_blif_cube*)malloc(sizeof(t_blif_cube));
	//t_blif_cubical_function *f_all_minterms = (t_blif_cubical_function*)malloc(sizeof(t_blif_cubical_function));
	t_blif_cubical_function *f_union_cover = (t_blif_cubical_function*)malloc(sizeof(t_blif_cubical_function));
	t_blif_cubical_function *f_no_redundant_cover = (t_blif_cubical_function*)malloc(sizeof(t_blif_cubical_function));
	t_blif_cubical_function *f_C_cover = (t_blif_cubical_function*)malloc(sizeof(t_blif_cubical_function));
	t_blif_cubical_function *f_G_cover = (t_blif_cubical_function*)malloc(sizeof(t_blif_cubical_function));
	t_blif_cubical_function *f_no_DC_cover = (t_blif_cubical_function*)malloc(sizeof(t_blif_cubical_function));
	t_blif_cubical_function *f_DC_cover = (t_blif_cubical_function*)malloc(sizeof(t_blif_cubical_function));
	t_blif_cubical_function *f_EPI = (t_blif_cubical_function*)malloc(sizeof(t_blif_cubical_function));
	t_blif_cubical_function *temp_f_EPI = (t_blif_cubical_function*)malloc(sizeof(t_blif_cubical_function));
	t_blif_cubical_function *f_non_EPI = (t_blif_cubical_function*)malloc(sizeof(t_blif_cubical_function));
	t_blif_cubical_function *f_not_covered_ON = (t_blif_cubical_function*)malloc(sizeof(t_blif_cubical_function));
	t_blif_cubical_function *temp_f_not_covered_ON = (t_blif_cubical_function*)malloc(sizeof(t_blif_cubical_function));
	t_blif_cube **ONSet_cover = (t_blif_cube **)malloc(1000 * sizeof(t_blif_cube*));
	t_blif_cube **sharp_result_cover = (t_blif_cube **)malloc(f->input_count*sizeof(t_blif_cube*));
	int cost[1000] = { INT_MAX };
	int optimization_index[1000] = { INT_MAX };
	//f_C_cover->set_of_cubes = (t_blif_cube **)malloc(f->cube_count * sizeof (t_blif_cube*));
	//f_G_cover->set_of_cubes = (t_blif_cube **)malloc((((f->cube_count - 1) * f->cube_count) / 2) * sizeof (t_blif_cube*));
	//f_union_cover->set_of_cubes = (t_blif_cube **)malloc((f_G_cover->cube_count + f_C_cover->cube_count)* sizeof(t_blif_cube *));
	//f_no_redundant_cover->set_of_cubes = (t_blif_cube **)malloc((f_G_cover->cube_count + f_C_cover->cube_count)* sizeof(t_blif_cube *));
	//C_temp_cover = (t_blif_cube **)malloc(f_C_cover->cube_count*sizeof(t_blif_cube*));
	f_C_cover->set_of_cubes = (t_blif_cube **)malloc(1000 * sizeof(t_blif_cube*));
	f_G_cover->set_of_cubes = (t_blif_cube **)malloc(1000 * sizeof(t_blif_cube*));
	f_union_cover->set_of_cubes = (t_blif_cube **)malloc(1000 * sizeof(t_blif_cube *));
	f_no_redundant_cover->set_of_cubes = (t_blif_cube **)malloc(1000 * sizeof(t_blif_cube *));
	f_no_DC_cover->set_of_cubes = (t_blif_cube **)malloc(1000 * sizeof(t_blif_cube *));
	f_not_covered_ON->set_of_cubes = (t_blif_cube **)malloc(1000 * sizeof(t_blif_cube *));
	f_DC_cover->set_of_cubes = (t_blif_cube **)malloc(100 * sizeof(t_blif_cube *));
	f_EPI->set_of_cubes = (t_blif_cube **)malloc(1000 * sizeof(t_blif_cube *));
	temp_f_EPI->set_of_cubes = (t_blif_cube **)malloc(1000 * sizeof(t_blif_cube *));
	temp_f_not_covered_ON->set_of_cubes = (t_blif_cube **)malloc(1000 * sizeof(t_blif_cube *));
	f_non_EPI->set_of_cubes = (t_blif_cube **)malloc(100 * sizeof(t_blif_cube *));
	C_temp_cover = (t_blif_cube **)malloc(1000 * sizeof(t_blif_cube*));
	input_count = f->input_count;
	f_C_cover->cube_count = f->cube_count;
	memcpy(f_C_cover->set_of_cubes, f->set_of_cubes, f->cube_count*sizeof(t_blif_cube));
	f_C_cover->input_count = f->input_count;
	
	// find all the minterms before * and # operation
	//f_all_minterms->set_of_cubes = (t_blif_cube **)malloc((1 << f->input_count) * sizeof(t_blif_cube*));
	//enumerateAllMinterms(f_all_minterms, f);
	printf("\n******************************\n");
	for (i = 0; i < f->cube_count; i++) {
		printf("Input Cover %i: ", i);

		for (j = 0; j < input_count; j++) {
			//printf("%d ", read_cube_variable(f_C_cover->set_of_cubes[i]->signal_status, j));
			printf("%s ", read_cube_variable(f_C_cover->set_of_cubes[i]->signal_status, j) == LITERAL_0 ? "0" : (read_cube_variable(f_C_cover->set_of_cubes[i]->signal_status, j) == LITERAL_1 ? "1" : "X"));
		}
		printf("\n");
	}
	printf("\n******************************\n");
	f_no_redundant_cover->cube_count = 0;
	iteration = 0;
	// find prime implicants
	// TODO: compare two covers
	while (is_two_covers_equal(f_C_cover, f_no_redundant_cover, input_count) == 0) {
		// copy new c into old c if two cover are not the same
		if (is_two_covers_equal(f_C_cover, f_no_redundant_cover, input_count) == 0 && iteration != 0) {
			//free(f_C_cover->set_of_cubes);
			//f_C_cover->set_of_cubes = (t_blif_cube**)malloc(f_no_redundant_cover->cube_count*sizeof(t_blif_cube*));
			memcpy(f_C_cover->set_of_cubes, f_no_redundant_cover->set_of_cubes, f_no_redundant_cover->cube_count * sizeof(t_blif_cube*));
			f_C_cover->cube_count = f_no_redundant_cover->cube_count;
		}
//		printf("********************************Iteration %d************************************\n", iteration++);
		iteration++;
		//for (i = 0; i < f_C_cover->cube_count; i++) {
		//	//copy_cube_value(&C_cover[i], &f->set_of_cubes[i]);

		//	printf("cube %i: ", i);
		//	for (j = 0; j < input_count; j++) {
		//		printf("%d ", read_cube_variable(f_C_cover->set_of_cubes[i]->signal_status, j));
		//	}
		//	printf("\n");
		//}
		f_G_cover->cube_count = 0;
		memcpy(C_temp_cover, f_C_cover->set_of_cubes, f_C_cover->cube_count * sizeof(t_blif_cube*));
		for (i = 0; i < f_C_cover->cube_count - 1; i++) {
			for (j = i + 1; j < f_C_cover->cube_count; j++) {
				star_operation(f_G_cover, f_C_cover->set_of_cubes[i], f_C_cover->set_of_cubes[j], f->input_count);
			}
		}
		memcpy(f_C_cover->set_of_cubes, C_temp_cover, f_C_cover->cube_count * sizeof(t_blif_cube*));
		//for (i = 0; i < f_C_cover->cube_count; i++) {
		//	//copy_cube_value(&C_cover[i], &f->set_of_cubes[i]);

		//	printf("after star C cover cube %i: ", i);
		//	for (j = 0; j < input_count; j++) {
		//		printf("%d ", read_cube_variable(f_C_cover->set_of_cubes[i]->signal_status, j));
		//	}
		//	printf("\n");
		//}
		//i = 0;
		//while (i < f_G_cover->cube_count) {
		//	printf("G cover %i: ", i);
		//	for (j = 0; j < input_count; j++) {
		//		printf("%d ", read_cube_variable(f_G_cover->set_of_cubes[i]->signal_status, j));
		//	}
		//	printf("\n");
		//	i++;
		//}
		memcpy(C_temp_cover, f_C_cover->set_of_cubes, f_C_cover->cube_count * sizeof(t_blif_cube*));
		union_cube_count = combine_cover(f_union_cover->set_of_cubes, f_C_cover->set_of_cubes, f_G_cover->set_of_cubes, input_count, f_C_cover->cube_count, f_G_cover->cube_count);

		//i = 0;
		//while (i < union_cube_count) {
		//	printf("union cover %i: ", i);	
		//	for (j = 0; j < input_count; j++) {
		//		printf("%d ", read_cube_variable(f_union_cover->set_of_cubes[i]->signal_status, j));
		//	}
		//	printf("\n");
		//	i++;
		//}
		// for redundant
		redundant_flag = 0;
		f_no_redundant_cover->cube_count = 0;

		for (i = 0; i < union_cube_count; i++) {
			for (j = 0; j < union_cube_count; j++) {
				if (i != j) {
					// sharp op returns the size of sharp_result_cover, if it returns -1, that means A is redundant
					if (sharp_operation(sharp_result_cover, f_union_cover->set_of_cubes[i], f_union_cover->set_of_cubes[j], sharp_result_cube, input_count) == -1) {
						redundant_flag = 1;
					}
				}
			}
			if (redundant_flag == 0) {
				//copy_cube_value(&union_cover_without_redundant[union_without_redundant_cube_count++], &union_cover[i]);
				f_no_redundant_cover->set_of_cubes[f_no_redundant_cover->cube_count] = f_union_cover->set_of_cubes[i];
				f_no_redundant_cover->cube_count++;
			}
			redundant_flag = 0;
		}
		memcpy(f_C_cover->set_of_cubes, C_temp_cover, f_C_cover->cube_count * sizeof(t_blif_cube*));
		//i = 0;
		//while (i < f_no_redundant_cover->cube_count) {
		//	printf("New C Cover %i: ", i);
		//	for (j = 0; j < input_count; j++) {
		//		printf("%d ", read_cube_variable(f_no_redundant_cover->set_of_cubes[i]->signal_status, j));
		//	}
		//	printf("\n");
		//	i++;
		//}
	}
	// update to the final C cover (might includes DCs)
	memcpy(f_C_cover->set_of_cubes, C_temp_cover, f_C_cover->cube_count * sizeof(t_blif_cube*));
	f_C_cover->cube_count = f_no_redundant_cover->cube_count;
	/*
	get all ON-set cubes from f
	int num_ON_cubes = 0;
	for (int i = 0; i < f->cube_count; i++) {
	if (f->set_of_cubes[i]->is_DC == T_FALSE) {
	ONSet_cover[num_ON_cubes] = f->set_of_cubes[i];
	num_ON_cubes++;
	}
	}*/

	/****** below codes are try to eliminate DCs within PIs ******/
	num_A = 0;
	sharp_result_cube_count = 0;
	num_of_cube_no_DC = 0;
	num_of_cube_DC = 0;
	// double loops to check any DCs within PIs
	for (i = 0; i < f_C_cover->cube_count; i++) {
		// reset
		num_A = 0;
		for (j = 0; j < f->cube_count; j++) {
			if (f->set_of_cubes[j]->is_DC == T_FALSE) {
				// 0 means C = A
				sharp_result_cube_count = sharp_operation(sharp_result_cover, f->set_of_cubes[j], f_C_cover->set_of_cubes[i], sharp_result_cube, input_count);
				if (sharp_result_cube_count == 0) {
					num_A++;
				}
				//else {
				//	printf("Sharp Result cube count %i: ", sharp_result_cube_count);
				//	for (k = 0; k < sharp_result_cube_count; k++) {
				//		printf("Sharp Rule 3 %i: ", k);
				//		for (h = 0; h < f->input_count; h++) {
				//			printf("%d ", read_cube_variable(sharp_result_cover[k]->signal_status, h));
				//		}
				//		printf("\n");
				//	}
				//}

			}
			else {
				num_A++;
			}
		}
		if (num_A != f->cube_count) {
			//this cube is DC cube 
			f_no_DC_cover->set_of_cubes[num_of_cube_no_DC++] = f_C_cover->set_of_cubes[i];
		}
		else {
			f_DC_cover->set_of_cubes[num_of_cube_DC++] = f_C_cover->set_of_cubes[i];
		}
	}
	f_no_DC_cover->cube_count = num_of_cube_no_DC;
	f_DC_cover->cube_count = num_of_cube_DC;
	// print out the final PIs (eliminated DC set within PIs)
	i = 0;
	while (i < f_no_DC_cover->cube_count) {
		printf("Actual PIs %i: ", i);
		for (j = 0; j < input_count; j++) {
			//printf("%d ", read_cube_variable(f_no_DC_cover->set_of_cubes[i]->signal_status, j));
			printf("%s ", read_cube_variable(f_no_DC_cover->set_of_cubes[i]->signal_status, j) == LITERAL_0 ? "0" : (read_cube_variable(f_no_DC_cover->set_of_cubes[i]->signal_status, j) == LITERAL_1 ? "1" : "X"));
		}
		printf("\n");
		i++;
	}
	printf("\n******************************\n");
	for (i = 0; i < f->cube_count; i++) {
		if (f->set_of_cubes[i]->is_DC == T_TRUE) {
			f_DC_cover->set_of_cubes[f_DC_cover->cube_count] = f->set_of_cubes[i];
			f_DC_cover->cube_count++;
		}
	}
	///****** below codes are trying to find EPIs from PIs ******/
	f_no_DC_cover->input_count = f_C_cover->input_count;
	f_non_EPI->input_count = input_count;
	find_EPI_cover_set(f_no_DC_cover, f_DC_cover, f_EPI, f_non_EPI);

	// print out the EPIs
	for (i = 0; i < f_EPI->cube_count; i++) {

		printf("Essential PIs %i: ", i);
		for (j = 0; j < input_count; j++) {
			//printf("%d ", read_cube_variable(f_EPI->set_of_cubes[i]->signal_status, j));
			printf("%s ", read_cube_variable(f_EPI->set_of_cubes[i]->signal_status, j) == LITERAL_0 ? "0" : (read_cube_variable(f_EPI->set_of_cubes[i]->signal_status, j) == LITERAL_1 ? "1" : "X"));
		}
		printf("\n");
	}
	//printf("*******************************\n");
	//for (i = 0; i < f_non_EPI->cube_count; i++) {

	//	printf("non EPIs %i: ", i);
	//	for (j = 0; j < input_count; j++) {
	//		printf("%d ", read_cube_variable(f_non_EPI->set_of_cubes[i]->signal_status, j));
	//	}
	//	printf("\n");
	//}
	f_not_covered_ON->cube_count = 0;
	f_not_covered_ON->input_count = input_count;
	get_not_covered_ON_set(f, f_EPI, f_not_covered_ON);
	if (f_not_covered_ON->cube_count == 0) {
		printf("\n******************************\n");
		for (i = 0; i < f_EPI->cube_count; i++) {
			printf("Minimal Cover %i: ", i);
			for (j = 0; j < input_count; j++) {
				//printf("%d ", read_cube_variable(f_EPI->set_of_cubes[i]->signal_status, j));
				printf("%s ", read_cube_variable(f_EPI->set_of_cubes[i]->signal_status, j) == LITERAL_0 ? "0" : (read_cube_variable(f_EPI->set_of_cubes[i]->signal_status, j) == LITERAL_1 ? "1" : "X"));
			}
			printf("\n");
		}
		printf("\n");
		memcpy(f->set_of_cubes, f_EPI->set_of_cubes, f_EPI->cube_count*sizeof(t_blif_cube));
		f->cube_count = f_EPI->cube_count;
	}
	else {
		// use the new cover set (selected PIs + EPIs)
		// calculate cover cost		//TODO
		//printf("*******************************\n");

		//for (i = 0; i < f_not_covered_ON->cube_count; i++) {
		//	printf("Uncovered ONs %i: ", i);
		//	for (j = 0; j < input_count; j++) {
		//		printf("%d ", read_cube_variable(f_not_covered_ON->set_of_cubes[i]->signal_status, j));
		//	}
		//	printf("\n");
		//}
		//TODO find good solutions

		for (i = 0; i < f_non_EPI->cube_count; i++) {
			cost[i] = cube_cost(f_non_EPI->set_of_cubes[i], input_count);
		}
		// count how many lowest cost and how many second lowest cost cube
		lowest_cost_cube_count = 0;
		second_lowest_cost_cube_count = 0;
		min_cost = 0;
		second_min_cost = 0;
		min_cost_cube_index = get_min_cost_cube(cost, f_non_EPI->cube_count);
		for (i = 0; i<f_non_EPI->cube_count; i++) {
			if (cost[i] == cost[min_cost_cube_index]) {
				lowest_cost_cube_count++;
				min_cost = cost[min_cost_cube_index];
			}
		}

		if (lowest_cost_cube_count != f_non_EPI->cube_count) {
			for (i = 0; i<lowest_cost_cube_count; i++) {
				min_cost_cube_index = get_min_cost_cube(cost, f_non_EPI->cube_count);
				cost[min_cost_cube_index] = INT_MAX;	
			}
			min_cost_cube_index = get_min_cost_cube(cost, f_non_EPI->cube_count);
			for (i = 0; i < f_non_EPI->cube_count; i++) {
				if (cost[i] == cost[min_cost_cube_index]) {
					second_lowest_cost_cube_count++;
					second_min_cost = cost[min_cost_cube_index];
				}
			}
		}

		for (i = 0; i < f_non_EPI->cube_count; i++) {
			cost[i] = cube_cost(f_non_EPI->set_of_cubes[i], input_count);
		}
		optimization_counter = (lowest_cost_cube_count > second_lowest_cost_cube_count) ? lowest_cost_cube_count : second_lowest_cost_cube_count;
		cost_for_optimization = (lowest_cost_cube_count > second_lowest_cost_cube_count) ? min_cost : second_min_cost;
		for (k = 0; k < optimization_counter; k++) {
			//printf("\n\n************* k = %d *****************", k);
			// f_not_covered_ON stays unchanged
			memcpy(temp_f_not_covered_ON->set_of_cubes, f_not_covered_ON->set_of_cubes, f_not_covered_ON->cube_count* sizeof(t_blif_cube));
			memcpy(temp_f_EPI->set_of_cubes, f_EPI->set_of_cubes, f_EPI->cube_count* sizeof(t_blif_cube));
			temp_f_EPI->cube_count = f_EPI->cube_count;
			temp_f_not_covered_ON->cube_count = f_not_covered_ON->cube_count;
			temp_f_not_covered_ON->input_count = f_not_covered_ON->input_count;
			temp_f_EPI->input_count = f_EPI->input_count;
			for (i = 0; i < f_non_EPI->cube_count; i++) {
				cost[i] = cube_cost(f_non_EPI->set_of_cubes[i], input_count);
			}
			min_cost_cube_index = get_min_cost_cube(cost, f_non_EPI->cube_count);
			while (cost[min_cost_cube_index] != cost_for_optimization) {
				cost[min_cost_cube_index] = INT_MAX;
				min_cost_cube_index = get_min_cost_cube(cost, f_non_EPI->cube_count);
			}
			h = 0;
			if (cost[min_cost_cube_index] == cost_for_optimization) {
				for (j = 0; j < k; j++) {
					cost[min_cost_cube_index] = INT_MAX;
					optimization_index[h++] = min_cost_cube_index;
					min_cost_cube_index = get_min_cost_cube(cost, f_non_EPI->cube_count);
				}
			}
			tracking_counter = 0;
			for (i = 0; i < f_non_EPI->cube_count; i++) {
				//// check temp_f_not_covered_ON->cube_count == 1
				//if (temp_f_not_covered_ON->cube_count == 1) {
				//	for (j = 0; j < f_non_EPI->cube_count; j++) {
				//		if (sharp_operation_EPI(sharp_result_cover, temp_f_not_covered_ON->set_of_cubes[0], f_non_EPI->set_of_cubes[j], sharp_result_cube, input_count) == FULLY_COVERED) {
				//			append_cubes(f_EPI, f_non_EPI->set_of_cubes[j]);
				//			temp_f_not_covered_ON->cube_count--;
				//			break;
				//		}
				//	}
				//}
				//if (temp_f_not_covered_ON->cube_count == 0) { break; }
				min_cost_cube_index = get_min_cost_cube(cost, f_non_EPI->cube_count);
				if (cost[min_cost_cube_index] == cost_for_optimization) {
					tracking_counter++;
				}
				// resume the ignored PIs
				if ((tracking_counter + k) == optimization_counter) {
					for (j = 0; j < k; j++) {
						cost[optimization_index[j]] = cost_for_optimization;
					}
				}
				remaining_set_flag = get_remaining_cover_set(temp_f_not_covered_ON, f_non_EPI->set_of_cubes[min_cost_cube_index]);
				cost[min_cost_cube_index] = INT_MAX;
				if (remaining_set_flag == NOT_COVERED_AT_ALL) {
					//		// skip this PI
					//		// continue to next loop (PI)
				}
				else if (remaining_set_flag == FULLY_COVERED) {
					//		// use this PI 
					//		// add this PI to EPI cover set
					//		// exit for loop
					append_cubes(temp_f_EPI, f_non_EPI->set_of_cubes[min_cost_cube_index]);
					break;
				}
				else {
					//		// PARTIALLY_COVERED
					//		// add this PI to EPI cover set
					//		// continue to next loop (PI)
					append_cubes(temp_f_EPI, f_non_EPI->set_of_cubes[min_cost_cube_index]);
				}
			}

			if (cover_cost(temp_f_EPI->set_of_cubes, temp_f_EPI->cube_count, input_count) < cover_cost(f->set_of_cubes, f->cube_count, input_count)) {
				memcpy(f->set_of_cubes, temp_f_EPI->set_of_cubes, temp_f_EPI->cube_count*sizeof(t_blif_cube));
				f->cube_count = temp_f_EPI->cube_count;
			}

		}

		printf("\n******************************\n");
		for (i = 0; i < f->cube_count; i++) {
			printf("Minimal Cover %i: ", i);
			for (j = 0; j < input_count; j++) {
				printf("%s ", read_cube_variable(f->set_of_cubes[i]->signal_status, j) == LITERAL_0? "0":(read_cube_variable(f->set_of_cubes[i]->signal_status, j) == LITERAL_1? "1": "X"));
			}
			printf("\n");
		}
		printf("\n");
	}

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
	printf("Star-Sharp Operation 2-level logic minimization program.\r\n");

	/* Read BLIF circuit. */
	printf("Reading file %s...\n", argv[1]);
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
			printf("Function %i: #inputs = %i; #cubes = %i; cost = %i\n", index + 1, function->input_count, function->cube_count, function_cost(function));
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