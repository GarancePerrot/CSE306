
#define _CRT_SECURE_NO_WARNINGS 1

#include "includes.h"

#include "scene_elements.cpp"
#include "triangle.cpp"


#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


int main() {
	clock_t start = clock(); //chronometer for execution time

	int W = 512; //width of the screen
	int H = 512; //height of the screen

	double alpha = 60 * PI /100; //horizontal field of view aka fov
	
	Scene s;
	s.L = Vector(-20, 20, 60); //light source
	s.I = 2E11; //light intensity
	Vector C(0,0,55); //camera

	int NB_PATHS = 32; 

	s.addSphere(new Sphere(Vector(0.,1000.,0.),   940.,  Vector(1,0.,0.))); //red ceiling 
    s.addSphere(new Sphere(Vector(0.,-1000.,0.),  940., Vector(0,1,0.))); //blue back wall
    s.addSphere(new Sphere(Vector(0.,0., -1000.), 940.,  Vector(0,0.,1))); //green floor
    s.addSphere(new Sphere(Vector(0.,0.,1000.),   940.,  Vector(0.8,0.,0.2))); //pink front wall
    s.addSphere(new Sphere(Vector(0.,-25.,-15.),  25,  Vector(0.1,0.1,0.1), false, false)); //center sphere : isMirror, isTransparent
	//s.addSphere(new Sphere(Vector(-50.,0.,0.), 5, Vector(0., 1, 0.)));
	//s.addSphere(new Sphere(Vector(50.,0.,0.), 20, Vector(0., 0, 1.)));

	double std_dev = 0.5;
	double z = -W/(2*tan(alpha/2)); 
	
	std::vector<unsigned char> image(W * H * 3, 0);
	#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector pixel_color(0.,0.,0.); 

			for (int k = 0 ; k< NB_PATHS; k++){ //antialiasing 
				double x_mul, y_mul; 
				boxMuller(std_dev, x_mul, y_mul);
				Vector u_mul(x_mul,y_mul,0.);  // contribution from boxMuller

				double x = j-W/2+0.5;  
				double y = H/2-i-0.5;
				Vector u(x, y, z); // coordinate of pixel (i,j) in the scene

				u += u_mul; //combining boxMuller and usual 
			
				u.normalize();
				Ray r(C, u); //from the coordinate of each pixel and the camera center, we compute a normalized ray direction
				int ray_depth = 5;  // controls the number of recursive calls allowed in getColor for refractions / reflections
				pixel_color += s.getColor(r, ray_depth); 
			}
			
			//without gamma correction : 
			// image[(i * W + j) * 3 + 0] = color[0];
			// image[(i * W + j) * 3 + 1] = color[1];
			// image[(i * W + j) * 3 + 2] = color[2];

			//not involved within the light simulation process, last step : 
			//applying gamma correction : elevating RGB values at 1/gamma where gamma = 2.2
			//and pixels values are limited to the range {0..255} to avoid overflowing unsigned char

			image[(i * W + j) * 3 + 0] = std::min(255.0, std::pow(pixel_color[0]/NB_PATHS, 1/2.2));
			image[(i * W + j) * 3 + 1] = std::min(255.0, std::pow(pixel_color[1]/NB_PATHS, 1/2.2));
			image[(i * W + j) * 3 + 2] = std::min(255., std::pow(pixel_color[2]/NB_PATHS, 1/2.2));
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	clock_t end = clock();
	double execution_time = ((double)end - start) / CLOCKS_PER_SEC;
	std::cout << "Code ran in: " << execution_time << " seconds." << std::endl;
	return 0;
}

//g++ -fopenmp -O2 main1.cpp -o test -lgomp
//./test


// Questions : 
// how to change the number of paths per pixel ? (fig 2.17 says 1000 vs 32)