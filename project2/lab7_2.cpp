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
    double area() { // computes the area of a polygon in 2D (formula p.34 slides)
        double res = 0;  // sum over the vertices
        int n = vertices.size();
        for (int i = 0; i< n; i++){
            int j = (i+1)%n; // avoid being out of bounds 
            res+= (vertices[i][0]* vertices[j][1] - vertices[i][1]* vertices[j][0]);
        }
        return 0.5*std::abs(res);
    }

    double sum_sq(const Vector& P){  // the sum of squared distances from a point to the polygon vertices 
        int n = vertices.size();
        double result = 0;  // sum over the vertices
        for (int t = 1 ; t < n-1 ; t++){
            Vector c[3] = {vertices[0], vertices[t], vertices[t+1]}; // define a triangle T by three vertices : the first of the polygon and two consecutive
            Vector V1 = c[1] - c[0]; 
            Vector V2 = c[2] - c[0];
            double areaT = 0.5 * std::abs(V2[1]*V1[2] - V2[2]*V1[1]);  // area of the triangle
            
            double sum_dots = 0.; 
            for (int k = 0; k < 3; k++){
                for (int l = k; l < 3; l++){
                    sum_dots += dot(c[k] - P, c[l] - P);  // formula p.34 slides
                }
            }
            result +=  ((1/6)*areaT) * sum_dots; 
        }
        return result;
    }

    std::vector<Vector> vertices;
};  


// new definition for power diagram
Vector intersection_point(const Vector& Mprime, const Vector& P0, const Vector& Pi, const Vector& A, const Vector& B ){ 
    double t = dot(Mprime-A, Pi-P0) / dot(B-A, Pi-P0);
    Vector P = A + t*(B-A);
    return P; 
}

// new definition for power diagram
bool is_inside(const Vector& X,  const Vector& P0, double w0, const Vector& Pi, double wi) { 
        if ((X-P0).norm2() - w0 <=(X-Pi).norm2() - wi){
            return true; 
        }
        return false;
    }


// Sutherland Hodgman algo
Polygon clip_by_bisector(const Polygon& cell, const Vector& P0, double w0,  const Vector& Pi , double wi){
    Polygon result;                         // create a new empty polygon
    Vector A(0,0,0);
    Vector B(0,0,0);
    Vector offset = (w0 - wi) / (2 * (P0- Pi).norm2()) * (Pi-P0);
    Vector M = (P0 + Pi ) *0.5;          // midpoint of bisector
    Vector Mprime = M + offset;
    int N = cell.vertices.size();

    for (int i=0; i< N; i++){                        // iterate over the edges
        if (i==0){                                   // we take the index (i-1)%N
            A = cell.vertices[N-1];
        } else {
            A = cell.vertices[i-1];
        } 
        B = cell.vertices[i]; 
        if (is_inside(B, P0,w0,Pi,wi)){                //if B is inside
            if ( ! is_inside(A, P0,w0, Pi, wi)){         // if A is outside      
                Vector P = intersection_point(Mprime, P0,Pi, A,B);   // we compute the point of intersection P
                result.vertices.push_back(P);              // and add it
            }
            result.vertices.push_back(B);   
        }
        else {
            if (is_inside(A, P0,w0, Pi, wi)){       // A is inside  
                Vector P = intersection_point(Mprime,P0,Pi,A,B);  // compute P
                result.vertices.push_back(P);         // add P
            }
        }
    }
    return result;
}

class PowerDiagram{
public: 
    PowerDiagram(){
    }
    void compute(){  // computes the cell for each point
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
#pragma omp parallel for  // parallelize the computation of cells for each point
        for (int i= 0 ; i<points.size(); i++){
            Polygon cell = square;              // initial cell shape
            for (int j= 0 ; j< points.size(); j++){   // iterating over all other points
                if (i==j) continue;                   // excluding itself
                cell = clip_by_bisector(cell, points[i], weights[i], points[j], weights[j]); // clip the current cell with the bissector of the two points 
            }
            cells[i] = cell;   // storing clipped polygon in cells
        }
    }
    std::vector<Vector> points; //within [0,1]
    std::vector<double> weights; 
    std::vector<Polygon> cells; 
};


class SemiDiscreteOT { // sets up and solves a Semi-Discrete Optimal Transport problem using L-BFGS optimization
public : 
    SemiDiscreteOT(){

    }

    static lbfgsfloatval_t _evaluate(  //static member function acting as a wrapper for calling the non-static evaluate method 
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        return reinterpret_cast<SemiDiscreteOT*>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(  // construct a power diagram whose weights are the variables passed in parameter x
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        for (int i=0; i<n; i++){
            diagram.weights[i] = x[i];  // update power diagram weights
        }
        diagram.compute();  // compute power diagram for current iteration 


        lbfgsfloatval_t fx = 0.0;
                                      // formulas from slides p.32
        double S1 = 0. ;              // the sum of squared distances from a point to the polygon vertices 
        double S2 = 0.;               // the sum of area * weights
        double S3 = 0. ;              // the sum of lambda * weights

        for (int i = 0; i < n; i ++) {  // sum over the vertices 
            double lambda = 1. / n;
            double area = diagram.cells[i].area();
            S2 += area *  diagram.weights[i];
            S3 += lambda * diagram.weights[i];
            S1 += diagram.cells[i].sum_sq(diagram.points[i]); 

            g[i] = -(lambda - area);  // reverse sign because function returns -g 
        }
        fx = -(S1 - S2 + S3);  
        // we reverse sign because function return minus the gradient g  of the objective function 
        //(since the objective function is maximized while this codes minimizes functions)
        return fx;
    }

    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        return reinterpret_cast<SemiDiscreteOT*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }

    int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        printf("Iteration %d:\n", k);   // prints info about current iteration
        printf("  fx = %f\n", fx);
        printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
        printf("\n");
        return 0;
    }

    void optimize(){
        N = diagram.points.size();
        diagram.weights.resize(N);
        for (int i=0; i<N; i++){
            diagram.weights[i] = 10;
        }
        double objectivefct = -1;
        lbfgs_parameter_t param;
        lbfgs_parameter_init(&param);
        param.linesearch = LBFGS_LINESEARCH_BACKTRACKING; // line search algorithm to use (default : BFGS_LINESEARCH_MORETHUENTE, much quicker )
        param.max_iterations = 1000; // sets a maximum number of iterations (default : 0 ,runs until convergence)
        param.epsilon = 1e-8; // sets a convergence tolerance for the gradient norm (default : 1e-5)

        std::vector<double> optimized_weights(N);
        lbfgs(N, &optimized_weights[0], &objectivefct, _evaluate, _progress, this, &param);

    }


    PowerDiagram diagram; // diagram that is updated at each iteration
    int N;     // nb of vertices
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
 
 
