#include "includes.h"
#include "lab8_2.cpp"



int main(){

    clock_t start = clock(); //chronometer for execution time

    Fluid fluid;
    int nb_points = 100; 
    int nb_steps = 157;
    fluid.simulate(nb_points, nb_steps); 

    clock_t end = clock();
	double execution_time = ((double)end - start) / CLOCKS_PER_SEC;
	std::cout << "Code ran in: " << execution_time << " seconds." << std::endl;
    return 0; 
}


//g++ main7_2.cpp lbfgs.c -fopenmp -IC:\Users\Garance\year3\CSE306\project2\ -o pg
