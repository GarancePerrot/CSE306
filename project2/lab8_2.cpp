#include "includes.h"

double sqr(double x) {return x*x;}; 
double PI = 3.1415926535;

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

 

class Polygon {  // in 2D
public:
    double area() { // computes the area of a polygon in 2D (formula p.34 slides)
        double res = 0;  // sum over the vertices
        int n = vertices.size();
        for (int i = 0; i< n; i++){
            int j = (i+1)%n; // avoid being out of bounds 
            res+= (vertices[i][0]* vertices[j][1] - vertices[i][1]* vertices[j][0]); // x_i * y_i+1 - x_i+1 y_i
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

    Vector centroid(){   
    // formulas from slides lab 7 p.9, "centroidal voronoi tesselations"
        double A = area();
        Vector sum_c(0,0,0);
        int n = vertices.size();
        for (int i = 0; i < n; i++){ //sum over vertices
            int j = (i+1)%n; // avoid being out of bounds
            double prod = vertices[i][0]*vertices[j][1] - vertices[j][0]*vertices[i][1]; // x_i * y_i+1 - x_i+1 * y_i
            double X = (vertices[i][0] + vertices[j][0]) * prod; // (x_i + x_i+1) * prod
            double Y = (vertices[i][1]+vertices[j][1]) * prod ; // (y_i + y_i+1) * prod
            sum_c += Vector(X, Y, 0); 
        }
        //return (1/(6*A)) * sum_c; // returns a Vector (C_x, C_y, 0)
        return (-1 / (6*A)) * sum_c; // multiplying by -1 for counter-clockwise polygons
    }

    std::vector<Vector> vertices;
};  


// definition for power diagram
Vector intersection_point(const Vector& Mprime, const Vector& P0, const Vector& Pi, const Vector& A, const Vector& B ){ 
    double t = dot(Mprime-A, Pi-P0) / dot(B-A, Pi-P0);
    Vector P = A + t*(B-A);
    return P; 
}

// definition for power diagram
bool is_inside(const Vector& X, const Vector& P0, double w0, const Vector& Pi, double wi) { 
        if ((X-P0).norm2() - w0 <=(X-Pi).norm2() - wi){
            return true; 
        }
        return false;
    }


// Sutherland Hodgman algo
void clip_by_bisector(Polygon& result, const Polygon& cell, const Vector& P0, double w0,  const Vector& Pi , double wi){
    
    Vector offset = (w0 - wi) / (2 * (P0- Pi).norm2()) * (Pi-P0);
    Vector M = (P0 + Pi ) *0.5;          // midpoint of bisector
    Vector Mprime = M + offset;
    int N = cell.vertices.size();

    for (int i=0; i< N; i++){                        // iterate over the edges
        const Vector& A = cell.vertices[i==0?(N-1): i-1];
	const Vector& B = cell.vertices[i];
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
}



 // slides lab 6 p.16 point of intersection for clipping by line
Vector intersection_point_new(const Vector& u, const Vector& A, const Vector& B, const Vector& normal ){ 
    double t =  dot(u-A, normal) / dot(B-A, normal);
    Vector P = A + t * (B-A);
    return P;
}

//same: slides lab 6 p.16
bool is_inside_new(const Vector& u,  const Vector& vec, const Vector& normal) { 
    if (dot(u-vec, normal ) >= 0){
        return true; 
    }
    return false;
}

//clips the current cell by the line defined by points u and v
//similar as clip_by_bisector but using new definitions of intersection_point and is_inside
void clip_by_line( Polygon& result, const Polygon& cell, const Vector& u,  const Vector& v){

    Vector normal(v[1] - u[1],u[0] - v[0], 0);
    int N = cell.vertices.size();

    for (int i=0; i< N; i++){                        // iterate over the edges
        const Vector& A = cell.vertices[i==0?(N-1): i-1];
	const Vector& B = cell.vertices[i];
        if (is_inside_new(u, B, normal)){                //if B is inside
            if ( ! is_inside_new(u , A , normal)){         // if A is outside      
                Vector P = intersection_point_new(u, A , B , normal);   // we compute the point of intersection P
                result.vertices.push_back(P);              // and add it
            }
            result.vertices.push_back(B);   
        }
        else {
            if (is_inside_new(u, A , normal)){       // A is inside  
                Vector P = intersection_point_new(u, A , B, normal);  // compute P
                result.vertices.push_back(P);         // add P
            }
        }
    }
}



class PowerDiagram{
public: 
    PowerDiagram(){
        for (int i = 0 ; i < Ncircle; i++){  // initialize a unit circle used for clipping by line
            double theta = i*2* PI / Ncircle;
            circle[i][0] = cos(theta);
            circle[i][1] = sin(theta);
            circle[i][2] = 0;
        }
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
        for (int i= 0 ; i < points.size(); i++){
            Polygon cell = square;              // initial cell shape
	    Polygon res ; 
            res.vertices.reserve(points.size());
            for (int j= 0 ; j< points.size(); j++){   // iterating over all other points
                if (i==j) continue;                   // excluding itself
                res.vertices.clear();
		clip_by_bisector(res, cell, points[i], weights[i], points[j], weights[j]); // clip the current cell with the bissector of the two points 
            	std::swap(res,cell);
	    }
            
            // performing an additional clipping step using a series of lines defined by points on a circle
            // (refining the bounding of cells can help for incompressibility)
            for (int j= 0 ; j< Ncircle; j++){
                if (i==j) continue;   
		res.vertices.clear();
                double radius = sqrt(weights[i] - w_air) ; //radius of circle intersecting the polygon, proportional to point's weight              
                Vector u = circle[j] * radius + points[i]; // defined on the perimeter of the circle centered at points[i]
                Vector v = circle[(j+1)%Ncircle] * radius + points[i]; // defined on the perimeter of the circle and adjacent to u
                clip_by_line(res, cell, u, v);
		std::swap(res, cell);

            }

            cells[i] = cell;   // storing clipped polygon in cells
        }
    }
    std::vector<Vector> points; //within [0,1]
    std::vector<double> weights; 
    std::vector<Polygon> cells; 


    double w_air;
    static const int Ncircle = 20;
    Vector circle[Ncircle];
    
};


class SemiDiscreteOT { // sets up and solves a Semi-Discrete Optimal Transport problem using lbfgs optimization
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
        for (int i=0; i<n-1; i++){
            diagram.weights[i] = x[i];  // update power diagram weights
        }
        diagram.w_air = x[n-1];
        diagram.compute();  // compute power diagram for current iteration 


        lbfgsfloatval_t fx = 0.0;
                                      // formulas from slides p.32
        double S1 = 0. ;              // the sum of squared distances from a point to the polygon vertices 
        double S2 = 0.;               // the sum of area * weights
        double S3 = 0. ;              // the sum of lambda * weights

        // we add air particles and make the result incompressible : 

        double estimated_fluid = 0.; // the sum of the area covered by fluid particles
        double lambda = 0.4 / (n-1); // desired volume of fluid : 40% of domain
        double desired_air = 1 - lambda; // desired volume of air : the rest of the domain 



        for (int i = 0; i < n-1 ; i ++) {  // sum over the vertices  UNTIL N-1 = N_fluid

            double area = diagram.cells[i].area();
            S1 += diagram.cells[i].sum_sq(diagram.points[i]); 
            S2 += area *  diagram.weights[i];
            S3 += lambda * diagram.weights[i];
            estimated_fluid += area;
            g[i] = -(lambda - area);  // reverse sign because function returns -g 
        }

        // corresponds to slides lab 8 p.23
        double estimated_air = 1 - estimated_fluid; //domain covered by air + fluid
        g[n-1] = estimated_air - desired_air;  // difference that we want to minimize

        fx = -(S1 - S2 + S3 + diagram.w_air *(desired_air - estimated_air) ); // add air particles
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

	for (int i = 0; i< n-1; i++){
		diagram.weights[i] = x[i];
	}
	diagram.w_air = x[n-1];
	diagram.compute();
	    
        printf("Iteration %d:\n", k);   // prints info about current iteration
        printf("  fx = %f\n", fx);
        printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
        printf("\n");
        return 0;
    }

    void optimize(){
        N = diagram.points.size();
        diagram.weights.resize(N);


        double objectivefct = -1;
        lbfgs_parameter_t param;
        lbfgs_parameter_init(&param);
        param.linesearch = LBFGS_LINESEARCH_BACKTRACKING; // line search algorithm to use (default : LBFGS_LINESEARCH_MORETHUENTE, quicker )
        param.max_iterations = 1000; // sets a maximum number of iterations (default : 0 ,runs until convergence)
        //param.epsilon = 1e-8; // sets a convergence tolerance for the gradient norm (default : 1e-5)

        //std::vector<double> optimized_weights(N);

        diagram.w_air = 0.;  // all air power cells have the same weight w_air
        std::vector<double> optimized_weights(N+1, 0.9);
        
        for (int i=0; i < N; i++){
            optimized_weights[i] = 1; 
            
        }
       
        lbfgs(N+1, &optimized_weights[0], &objectivefct, _evaluate, _progress, this, &param);
    }


    PowerDiagram diagram; // diagram that is updated at each iteration
    int N;     // nb of vertices
};




int sgn(double x){  // used in save_frame
    if (x > 0) {
        return 1;
    }
    if (x < 0 ) {
        return -1;
    }
    return 0;
}

void save_frame(const std::vector<Polygon> &cells, std::string filename, int frameid = 0) {
    int W = 500, H = 500;
    std::vector<unsigned char> image(W*H * 3, 255);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < cells.size(); i++) {

        double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
        for (int j = 0; j < cells[i].vertices.size(); j++) {
            bminx = std::min(bminx, cells[i].vertices[j][0]);
            bminy = std::min(bminy, cells[i].vertices[j][1]);
            bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
            bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
        }
        bminx = std::min(W-1., std::max(0., W * bminx));
        bminy = std::min(H-1., std::max(0., H * bminy));
        bmaxx = std::max(W-1., std::max(0., W * bmaxx));
        bmaxy = std::max(H-1., std::max(0., H * bmaxy));

        for (int y = bminy; y < bmaxy; y++) {
            for (int x = bminx; x < bmaxx; x++) {
                int prevSign = 0;
                bool isInside = true;
                double mindistEdge = 1E9;
                for (int j = 0; j < cells[i].vertices.size(); j++) {
                    double x0 = cells[i].vertices[j][0] * W;
                    double y0 = cells[i].vertices[j][1] * H;
                    double x1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][0] * W;
                    double y1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][1] * H;
                    double det = (x - x0)*(y1-y0) - (y - y0)*(x1-x0);
                    int sign = sgn(det);   // to implement : 1 or -1 
                    if (prevSign == 0) prevSign = sign; else
                        if (sign == 0) sign = prevSign; else
                        if (sign != prevSign) {
                            isInside = false;
                            break;
                        }
                    prevSign = sign;
                    double edgeLen = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
                    double distEdge = std::abs(det)/ edgeLen;
                    double dotp = (x - x0)*(x1 - x0) + (y - y0)*(y1 - y0);
                    if (dotp<0 || dotp>edgeLen*edgeLen) distEdge = 1E9;
                    mindistEdge = std::min(mindistEdge, distEdge);
                }
                if (isInside) {  // UNCOMMENT FOR BLUE PARTICLES
                    if (i < cells.size()) {   // the N first particles may represent fluid, displayed in blue
                        image[((H - y - 1)*W + x) * 3] = 0;
                        image[((H - y - 1)*W + x) * 3 + 1] = 0;
                        image[((H - y - 1)*W + x) * 3 + 2] = 255;
                    }
                    if (mindistEdge <= 2) {
                        image[((H - y - 1)*W + x) * 3] = 0;
                        image[((H - y - 1)*W + x) * 3 + 1] = 0;
                        image[((H - y - 1)*W + x) * 3 + 2] = 0;
                    }

                }  
            }
        }
    }
    std::ostringstream os;
    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}



class Fluid {
public : 
    Fluid(){  // we use parameters from slides p. 17 (epsilon, dt, m_i)
        epsilon = 0.004 ;  
        m_i = 200;
        g = Vector(0,-9.81, 0);
    }

    void time_step(double dt){ 
        //  Algorithm at time step t : [Gallouët & Mérigot 2016] Recover incompressibility throught OT

        OT.diagram.points = positions;
        OT.optimize();    // Compute current OT from particles to uniform density

        for (int i = 0; i < positions.size(); i++){   
            // Add force from particle towards power cell centroid
            Vector F_spring = (1. / sqr(epsilon)) * (OT.diagram.cells[i].centroid()- positions[i]);  // spring force to attract to the point of power diagram 
            Vector F_i = F_spring + (m_i*g);    // Add forces as usual : spring and gravity (Lagrangian scheme)
            velocity[i] += (dt/m_i) *F_i;  // notation in slides : (v_i)^(t+1)
            positions[i] += dt*velocity[i];  // notation in slides : (X_i)^(t+1)

            // bounce particles back into the domain : 
            positions[i][0] = std::min(1 - 1E-9, std::max(positions[i][0], 1E-9 ));
            positions[i][1] = std::min(1 - 1E-9, std::max(positions[i][1], 1E-9 ));
        }
    }

    void simulate(int nb_points, int nb_steps){  // combine everything

        int N = nb_points;
        positions.resize(N);
        velocity.resize(N, Vector(0,0,0));
        for (int i = 0; i < N; i++){
            positions[i] = Vector(uniform(engine)/ double(RAND_MAX),uniform(engine)/ double(RAND_MAX), 0);
        }
        OT.diagram.points = positions; // random set of points
        OT.optimize();
        for (int i = 0 ; i < nb_steps ; i++){
            save_frame(OT.diagram.cells, "anim", i);
            time_step(0.002);   // dt = 0.002
        }
    }

    SemiDiscreteOT OT;
    std::vector<Vector> positions;
    std::vector<Vector> velocity;
    double m_i;
    double epsilon;
    Vector g;

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
 
 
