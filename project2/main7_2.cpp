
#include "includes.h"
#include "lab7_2.cpp"



void draw_result(int N){
    SemiDiscreteOT ot; 
    ot.diagram.points.resize(N);
    ot.diagram.weights.resize(N);
    for (int i= 0 ; i<N; i++ ){
        ot.diagram.points[i] = Vector(uniform(engine)/ double(RAND_MAX),uniform(engine)/ double(RAND_MAX), 0 );
        ot.diagram.weights[i] = uniform(engine)/ double(RAND_MAX);
    }
    ot.optimize();
    save_svg(ot.diagram.cells, "OT.svg");

}

int main(){

    // checking if an old file "OT.svg" exists and removes it if so because 
    // otherwise it messes up with the optimization
    std::string filename = "OT.svg";
    std::filesystem::path current_path = std::filesystem::current_path();
    std::filesystem::path file_path = current_path / filename;
    if (std::filesystem::exists(file_path)) {
        std::filesystem::remove(file_path);
    } 
    clock_t start = clock(); //chronometer for execution time

    // displays cells for a random set of N points: 
    draw_result(200);

    clock_t end = clock();
	double execution_time = ((double)end - start) / CLOCKS_PER_SEC;
	std::cout << "Code ran in: " << execution_time << " seconds." << std::endl;
    return 0; 
}


//g++ main7_2.cpp lbfgs.c -fopenmp -IC:\Users\Garance\year3\CSE306\project2\ -o pg
