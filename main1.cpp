
#define _CRT_SECURE_NO_WARNINGS 1

#include "includes.h"

#include "scene_elements.cpp"


#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


int main() {

	int W = 512; //width of the screen
	int H = 512; //height of the screen

	double alpha = 60 * PI /100; //horizontal field of view aka fov
	
	Scene s;
	s.L = Vector(-10, 20, 40); //light source
	s.I = 2E10; //light intensity
	Vector C(0,0,55); //camera

	s.addSphere(Sphere(Vector(0.,1000.,0.),   940.,  Vector(1,0.,0.))); //red ceiling 
    s.addSphere(Sphere(Vector(0.,-1000.,0.),  940., Vector(0,1,0.))); //blue back wall
    s.addSphere(Sphere(Vector(0.,0., -1000.), 940.,  Vector(0,0.,1))); //green floor
    s.addSphere(Sphere(Vector(0.,0.,1000.),   940.,  Vector(0.8,0.,0.2))); //pink front wall
    s.addSphere(Sphere(Vector(0.,-20.,0.),  20,  Vector(0.1,0.1,0.1), false, false)); //center sphere : isMirror, isTransparent
	//s.addSphere(Sphere(Vector(-50.,0.,0.), 5, Vector(0., 1, 0.)));
	//s.addSphere(Sphere(Vector(50.,0.,0.), 20, Vector(0., 0, 1.)));


	double z = -W/(2*tan(alpha/2)); 
	
	std::vector<unsigned char> image(W * H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			double x = j-W/2+0.5;
			double y = H/2-i-0.5;
			Vector u(x, y, z); //coordinate of pixel (i,j) in the scene
			u.normalize();
			Ray r(C, u); //from the coordinate of each pixel and the camera center, we compute a normalized ray direction
			int ray_depth = 5;  // controls the number of recursive calls allowed in getColor for refractions / reflections
			Vector color = s.getColor(r, ray_depth); 
			
			//without gamma correction : 
			// image[(i * W + j) * 3 + 0] = color[0];
			// image[(i * W + j) * 3 + 1] = color[1];
			// image[(i * W + j) * 3 + 2] = color[2];

			//not involved within the light simulation process, last step : 
			//applying gamma correction : elevating RGB values at 1/gamma where gamma = 2.2
			//and pixels values are limited to the range {0..255} to avoid overflowing unsigned char

			image[(i * W + j) * 3 + 0] = std::min(255.0, std::pow(color[0], 1/2.2));
			image[(i * W + j) * 3 + 1] = std::min(255.0, std::pow(color[1], 1/2.2));
			image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2], 1/2.2));
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}
