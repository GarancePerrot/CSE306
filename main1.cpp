
#define _CRT_SECURE_NO_WARNINGS 1

#include "includes.h"

#include "scene_elements.cpp"


#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


int main() {

	int W = 512;
	int H = 512;

	double alpha = 60 * PI /100;
	
	Scene s;
	s.L = Vector(-10, 20, 40);
	s.I = 2E10;
	Vector C(0,0,55);

	s.addSphere(Sphere(Vector(0.,1000.,0.),   940.,  Vector(1,0.,0.), false, false));
    s.addSphere(Sphere(Vector(0.,-1000.,0.),  940., Vector(0,1,0.)));
    s.addSphere(Sphere(Vector(0.,0., -1000.), 940.,  Vector(0,0.,1)));
    s.addSphere(Sphere(Vector(0.,0.,1000.),   940.,  Vector(0.8,0.,0.2)));
    s.addSphere(Sphere(Vector(0.,-20.,0.),  20,  Vector(0.1,0.1,0.1), false, false));
	//s.addSphere(Sphere(Vector(-50.,0.,0.), 5, Vector(0., 1, 0.), false, false));
	//s.addSphere(Sphere(Vector(50.,0.,0.), 20, Vector(0., 0, 1.), false, false));


	double z = -W/(2*tan(alpha/2));
	


	std::vector<unsigned char> image(W * H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			double x = j-W/2+0.5;
			double y = H/2-i-0.5;
			Vector u(x, y, z);
			u.normalize();
			Ray r(C, u);

			Vector color = s.getColor(r, 5);

			image[(i * W + j) * 3 + 0] = std::min(255.0, std::pow(color[0], 1/2.2));
			image[(i * W + j) * 3 + 1] = std::min(255.0, std::pow(color[1], 1/2.2));
			image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2], 1/2.2));
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}