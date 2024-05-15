
#include "includes.h"
#include "lab6.cpp"


int main(){
    clock_t start = clock(); //chronometer for execution time

    // std::vector<Polygon> test; 
    // Polygon pol1, pol2; 
    // pol1.vertices.push_back(Vector(0,0,0));
    // pol1.vertices.push_back(Vector(0,1,0));
    // pol1.vertices.push_back(Vector(1,1,0));
    // pol2.vertices.push_back(Vector(0.5,0.5,0));
    // pol2.vertices.push_back(Vector(0.5,1,0));
    // pol2.vertices.push_back(Vector(1,0.3,0));
    // pol2.vertices.push_back(Vector(0.3,0.3,0));
    // test.push_back(pol1);
    // test.push_back(pol2);

    Voronoi vor; 
    int N = 2000; 
    vor.points.resize(N);

    for (int i= 0 ; i<N; i++ ){
        vor.points[i] = Vector(rand()/ double(RAND_MAX),rand()/ double(RAND_MAX), 0 );
    }
    vor.compute();
    save_svg(vor.cells, "test.svg");

    clock_t end = clock();
	double execution_time = ((double)end - start) / CLOCKS_PER_SEC;
	std::cout << "Code ran in: " << execution_time << " seconds." << std::endl;
    return 0; 
}