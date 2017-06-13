#include <stdio.h>
#include <math.h>
#include <hpx/hpx.h>
#include <iostream>

// problem description
const double L = 1.0; // length scale 
const int I = 10; // number of iterations 
const int N = 81; // grid size N 
const double dT = 1e-3; // timestep 
 
// globals 
double h; 
 
// grid point type
typedef struct {
  double T[2]; 
} point_t; 
 
 
// handler and action
int hpx_main_handler(int print); 
HPX_ACTION(HPX_DEFAULT, 0, hpx_main, hpx_main_handler, HPX_INT); 

int hpx_clean_handler(hpx_addr_t grid);
HPX_ACTION(HPX_DEFAULT, 0, hpx_clean, hpx_clean_handler, HPX_ADDR);

int spmd_test_handler();
HPX_ACTION(HPX_DEFAULT, 0, spmd_test_, spmd_test_handler);

int clean_broadcast_handler(hpx_addr_t mem);
HPX_ACTION(HPX_DEFAULT, 0, clean_broadcast_, clean_broadcast_handler, HPX_ADDR);



int spmd_test_handler(){
	std::cout<<"In broadcast handler"<<std::endl;
	int num_ranks = hpx_get_num_ranks(); 
	int my_rank = hpx_get_my_rank();
	int* my_array = new int[num_ranks];
	my_array[my_rank] = my_rank;

	std::cout<<"Leaving broadcast handler"<<std::endl;
	hpx_exit(sizeof(int*),&my_array);
}

//creates some local memory on other localities 
int* spmd_test(){
	
	int* my_array;// = new int[num_ranks];

	hpx_run_spmd(&spmd_test_, &my_array);
	return my_array;
}


int do_clean_broadcast_handler(hpx_addr_t mem){
	return 0;
}

int clean_broadcast_handler(hpx_addr_t mem){
	return 0;
}

//frees local memory on other localities
int clean_broadcast(hpx_addr_t mem){
	int e = hpx_run(&clean_broadcast_, NULL, &mem);
	return e;
}


// implementation
int hpx_main_handler(int print) {
	std::cout<<"In hpx_main_handler:"<<std::endl;
  // allocate memory 
  hpx_addr_t grid = hpx_gas_calloc_cyclic(N, sizeof(point_t), 0); 
  assert(grid != HPX_NULL); 
  printf("allocation is done\n"); 
	std::cout<<"Leaving hpx_main_handler:"<<std::endl;

	hpx_exit(sizeof(hpx_addr_t),&grid);
} 

hpx_addr_t initialize(int print){
	std::cout<<"In initialize:"<<std::endl;
	hpx_addr_t grid;

	hpx_run(&hpx_main, &grid, &print); 
	//e = hpx_run(&hpx_clean, NULL, &grid);
	
	std::cout<<"Leaving initialize:"<<std::endl;
	return grid;
}



int hpx_clean_handler(hpx_addr_t grid){
	std::cout<<"In clean:"<<std::endl;
	// free global memory 
  hpx_gas_free(grid, HPX_NULL); 
 
  // finalize 
	std::cout<<"Leaving clean:"<<std::endl;
  hpx_exit(0, NULL); 
}

int clean(hpx_addr_t grid){
	int e = hpx_run(&hpx_clean, NULL, &grid);
	return e;
}
 
int main(int argc, char *argv[]) {
  // initialize hpx 
  if (hpx_init(&argc, &argv)) {
    hpx_print_help(); 
    return -1;
  }
 
  // parse application arguments
  int print = (argc > 1); 
 
  h = L / (N - 1); 
	//create some global memory
	hpx_addr_t grid = initialize(print);

	//create an array of ranks from multiple localities
	int* my_array = spmd_test();
	for(int i=0;i<hpx_get_num_ranks();++i){
		std::cout<<my_array[i]<<" "<<std::endl;
	}

  //free some global memory
	int e = clean(grid);

  hpx_finalize();


	std::cout<<"End"<<std::endl;
  return e;
} 
