
#include "includes.h"
#include "lab6.cpp"

void draw_heart(){
    std::vector<Polygon> test;
    Polygon pol1;
    pol1.vertices.push_back(Vector(0.5, 0.2, 0)); 
    pol1.vertices.push_back(Vector(0.3, 0.4, 0));  
    pol1.vertices.push_back(Vector(0.15, 0.55, 0));  
    pol1.vertices.push_back(Vector(0.1, 0.7, 0));  
    pol1.vertices.push_back(Vector(0.2, 0.85, 0));  
    pol1.vertices.push_back(Vector(0.35, 0.85, 0));  
    pol1.vertices.push_back(Vector(0.5, 0.8, 0));  
    pol1.vertices.push_back(Vector(0.65, 0.85, 0));  
    pol1.vertices.push_back(Vector(0.8, 0.85, 0));  
    pol1.vertices.push_back(Vector(0.9, 0.7, 0));  
    pol1.vertices.push_back(Vector(0.85, 0.55, 0));  
    pol1.vertices.push_back(Vector(0.7, 0.4, 0));  
    pol1.vertices.push_back(Vector(0.5, 0.2, 0)); 
    test.push_back(pol1);
    save_svg(test, "heart.svg"); 
}


void draw_voronoi(int N){
    Voronoi vor; 
    vor.points.resize(N);
    for (int i= 0 ; i<N; i++ ){
        vor.points[i] = Vector(rand()/ double(RAND_MAX),rand()/ double(RAND_MAX), 0 );
    }
    vor.compute();
    save_svg(vor.cells, "Voronoi_cells.svg");

}

int main(){
    clock_t start = clock(); //chronometer for execution time

    // //display a simple polygon in heart shape : 
    //draw_heart();

    // display the Voronoi cells of a random set of N points: 
    draw_voronoi(2000); 


    clock_t end = clock();
	double execution_time = ((double)end - start) / CLOCKS_PER_SEC;
	std::cout << "Code ran in: " << execution_time << " seconds." << std::endl;
    return 0; 
}
