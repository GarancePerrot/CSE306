#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <algorithm>
#include <random>
#include <stdio.h>
#include <math.h>

double PI = 3.1415926535897932384626433;
static std::default_random_engine engine(10); // random seed = 10
static std::uniform_real_distribution<double> uniform (0, 1);


//color matching  : 
// transfering the color palette from one image to another


void sliced(double* src, double* tgt, int Npix){ //src = input image, tgt = model image, same nb of pixels

    // We store the dot product and pixel index as a pair of values
	std::vector<std::pair<double, int>> projectionSrc(Npix);
	std::vector<double> projectionTgt(Npix);
	
	for (int iter = 0; iter < 60; iter++){ // nb of iterations, use 100

		// we find a uniformly random direction on the sphere
		double r1 = uniform(engine);
		double r2 = uniform(engine);
		double xLine = cos(2*PI*r1) * sqrt(r2*(1-r2));
		double yLine = sin(2*PI*r1) * sqrt(r2*(1-r2));
		double zLine = 1 - 2*r2;
		double norm = sqrt(xLine*xLine + yLine*yLine + zLine*zLine);
		xLine /=norm;
		yLine /= norm;
		zLine /= norm;
		
		// we project the input and model point clouds (i.e., pixel RGB values) onto this direction 
		//simply computing the dot product between this random direction and the pixel coordinate
		for (int i = 0; i<Npix; i++){  //for each pixel
			projectionSrc[i] = std::pair<double, int>(src[i*3+0]*xLine + src[i*3+1]*yLine + src[i*3+2]*zLine , i);
			projectionTgt[i] = tgt[i*3+0]*xLine + tgt[i*3+1]*yLine + tgt[i*3+2]*zLine ;
		}
		
		// Goal : matching the first projected point of the input point cloud to the first projected point 
		// of the model point cloud, so we sort according to the dot product
		std::sort(projectionSrc.begin(),projectionSrc.end());
		std::sort(projectionTgt.begin(),projectionTgt.end());

		// Advect each point of the initial cloud to its matched point in the model image
		for (int i = 0; i<Npix; i++){
			src[projectionSrc[i].second * 3 + 0] += (projectionTgt[i] - projectionSrc[i].first) * xLine;
			src[projectionSrc[i].second * 3 + 1] += (projectionTgt[i] - projectionSrc[i].first) * yLine;
			src[projectionSrc[i].second * 3 + 2] += (projectionTgt[i] - projectionSrc[i].first) * zLine;
		}
	}
}



int main() {

	int W, H, C;
	int W2, H2, C2;

	//stbi_set_flip_vertically_on_load(true);
	unsigned char *image = stbi_load("8733654151_b9422bb2ec_k.jpg",
                                 &W,
                                 &H,
                                 &C,
                                 STBI_rgb);

	unsigned char* target = stbi_load("redim.jpg",
                                 &W2,
                                 &H2,
                                 &C2,
                                 STBI_rgb);


	std::vector<double> image_double(W*H*3);
	std::vector<double> target_double(W2*H2*3);

	for (int i=0; i<W*H*3; i++)
		image_double[i] = image[i];
	for (int i=0; i<W2*H2*3; i++)
		target_double[i] = target[i];

	sliced(&image_double[0], &target_double[0], W * H);

	std::vector<unsigned char> image_result(W*H * 3, 0);
	//pixels values are clamped in the range {0..255} to avoid overflowing unsigned char

	for (int i=0; i<W*H*3; i++){
		image_result[i] = (unsigned char)(std::max(0., std::min(255., image_double[i])));
	}

	stbi_write_png("image.png", W, H, 3, &image_result[0], 0);

	return 0;
}