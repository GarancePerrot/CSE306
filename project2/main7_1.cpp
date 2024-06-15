#include "includes.h"
#include "lab7_1.cpp"


void draw_powerdiag(int N){
    PowerDiagram diag; 
    diag.points.resize(N);
    diag.weights.resize(N);
    for (int i= 0 ; i<N; i++ ){
        diag.points[i] = Vector(uniform(engine)/ double(RAND_MAX),uniform(engine)/ double(RAND_MAX), 0 );
        diag.weights[i] = uniform(engine)/ double(RAND_MAX);
    }
    diag.compute();
    save_svg(diag.cells, "Power_diagram.svg");

}

int main(){
    clock_t start = clock(); //chronometer for execution time

    // display the power diagram of a random set of N points: 
    draw_powerdiag(20000); 

    clock_t end = clock();
	double execution_time = ((double)end - start) / CLOCKS_PER_SEC;
	std::cout << "Code ran in: " << execution_time << " seconds." << std::endl;
    return 0; 
}
