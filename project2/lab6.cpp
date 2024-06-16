#include "includes.h"

double sqr(double x) {return x*x;}; 


class Vector {  // we use the same class vector but leave the third coordinate to 0 as the points are of two dimension
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {  
		coords[0] = x;
		coords[1] = y;
		coords[2] = z;
	}
	double norm2() const {
		return coords[0] *coords[0] +coords[1] * coords[1] + coords[2] * coords[2];
	}
	double norm() const {
		return sqrt(norm2());
	}
	void normalize() {
		double n = norm();
		coords[0] /= n;
		coords[1] /= n;
		coords[2] /= n;
	}
    Vector& operator +=(const Vector& b) {
        coords[0] += b[0] ;
        coords[1] += b[1] ;
        coords[2] += b[2] ;
        return (*this);
    }
	double operator[](int i) const { return coords[i]; };
	double& operator[](int i) { return coords[i]; };
	double coords[3];

};

Vector operator-(const Vector& a) {
	return Vector(-a[0], -a[1], -a[2]);
}
Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]); 
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

 

class Polygon {  
public:
    std::vector<Vector> vertices;
};  

// computes the intersection point between bisector of [P0,Pi] and given edge [AB] of a polygon
// from slides p.61
Vector intersection_point(const Vector& P0, const Vector& Pi, const Vector& A, const Vector& B ){ //add comment
    Vector M = (P0 + Pi ) * 0.5;
    double t = dot(M-A, Pi-P0) / dot(B-A, Pi-P0);
    Vector P = A + t*(B-A);
    return P; 
}

// checks if X follows the definition of the Voronoi cell (slides p.61)
bool is_inside(const Vector& point,  const Vector& P0, const Vector& Pi) { // add comment 
        if ((point-P0).norm2()<=(point-Pi).norm2()){
            return true; 
        }
        return false;
    }


// Sutherland Hodgman to clip cell by bisector of [P0,Pi]
// overall algorithm to remove 1 half-space (slides p.16)
void clip_by_bisector(Polygon& result, const Polygon& cell, const Vector& P0, const Vector& Pi ){

    int N = cell.vertices.size();

    for (int i=0; i< N; i++){                        // iterate over the edges
        const Vector& A = cell.vertices[i==0?(N-1): i-1];
	const Vector& B = cell.vertices[i];
        if (is_inside(B, P0,Pi)){                //if B is inside
            if ( ! is_inside(A, P0,Pi)){         // if A is outside      
                Vector P = intersection_point(P0,Pi,A,B);   // we compute the point of intersection P
                result.vertices.push_back(P);              // and add it
            }
            result.vertices.push_back(B);   
        }
        else {
            if (is_inside(A, P0,Pi)){       // A is inside  
                Vector P = intersection_point(P0,Pi,A,B);  // compute P
                result.vertices.push_back(P);         // add P
            }
        }
    }
    //return result; // changed Polygon into void
}

class Voronoi{
public: 
    Voronoi(){
    }
    void compute(){  // computes the voronoi cell for each point
        Polygon square;  
        square.vertices.push_back(Vector(0,0,0));  // anti-clockwise 
        square.vertices.push_back(Vector(1,0,0));
        square.vertices.push_back(Vector(1,1,0));
        square.vertices.push_back(Vector(0,1,0));

        // square.vertices.push_back(Vector(0,0,0));  // clockwise -> does not change anything, even runtime
        // square.vertices.push_back(Vector(0,1,0));  
        // square.vertices.push_back(Vector(1,1,0));
        // square.vertices.push_back(Vector(1,0,0));

        cells.resize(points.size());
#pragma omp parallel for  // parallelize the computation of Voronoi cells for each point
        for (int i= 0 ; i<points.size(); i++){
            Polygon cell = square;              // initial cell shape
            Polygon res ; 
            res.vertices.reserve(50);
            for (int j= 0 ; j< points.size(); j++){   // iterating over all other points
                if (i==j) continue;                   // excluding itself
                res.vertices.clear();
                clip_by_bisector(res, cell, points[i], weights[i], points[j], weights[j]); // clip the current cell with the bissector of the two points 
                std::swap(res,cell); 
            }
            cells[i] = cell;   // storing clipped polygon in cells
        }
    }
    std::vector<Vector> points; //within [0,1]
    std::vector<Polygon> cells; 
};




// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
    void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
        FILE* f = fopen(filename.c_str(), "w+"); 
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        for (int i=0; i<polygons.size(); i++) {
            fprintf(f, "<g>\n");
            fprintf(f, "<polygon points = \""); 
            for (int j = 0; j < polygons[i].vertices.size(); j++) {
                fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
            }
            fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
            fprintf(f, "</g>\n");
        }
        fprintf(f, "</svg>\n");
        fclose(f);
    }
 
 
