#include "includes.h"

double PI = 3.1415926535897932384626433832795;
double sqr(double x) { return x * x;}

class Vector {
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
Vector random_cos(const Vector &N) { // for indirect lighting

	double r1 = uniform(engine);
	double r2 = uniform(engine);
	double x = cos(2*PI*r1)*sqrt(1-r2);
	double y = sin(2*PI*r1)*sqrt(1-r2);
	double z = sqrt(r2);
	
    
	//generating the two orthogonal tangent vectors T1 and T2 : 
	Vector T1;
	//we detect the smallest component of N (in absolute value) and set T1 accordingly :
	if (std::abs(N[0]) <= std::abs(N[1]) && std::abs(N[0]) <= std::abs(N[2])) { //N[0] min
			T1 = Vector(0, -N[2], N[1]);
		} else {
			if (std::abs(N[1]) <= std::abs(N[0]) && std::abs(N[1]) <= std::abs(N[2])){ //N[1] min
				T1 = Vector(-N[2], 0, N[0]);
			} else {
				T1 = Vector(-N[1], N[0], 0); // N[2] min
			}
		}
	T1.normalize();
	Vector T2 = cross(N, T1);
	return x*T1 + y*T2 + z*N;
}

void boxMuller (double stdev , double &x , double &y ) {  // for antialiasing 
	double r1 = uniform(engine);
	double r2 = uniform(engine);
	x = sqrt(-2 * log(r1)) *cos( 2 * PI*r2)*stdev ;
	y = sqrt(-2 * log(r1))*sin( 2 * PI*r2) *stdev ;
}


class Ray {
public:
	Ray(Vector O, Vector u) : O(O), u(u) {};
    Vector O; //origin
	Vector u; //unit direction
};


class Geometry {  //general Geometry abstract class to add different objects in the scene
public:
    Geometry(const Vector& albedo, bool isMirror=false, bool isTransparent=false) : albedo(albedo), isMirror(isMirror), isTransparent(isTransparent) {};
    bool isMirror, isTransparent;
    Vector albedo;
    virtual std::string getType() const { return "Geometry"; } // base class returns generic type
};

class Sphere : public Geometry { // changed sphere to inherit from Geometry class
public:
	Sphere(const Vector& C, double R) : C(C), R(R) {};

	Vector C; //center
	double R; //radius
	std::string getType() const override { return "Sphere"; }
};

class IntersectionSphere { //the intersection between a ray and a sphere
public:
    IntersectionSphere(Sphere* S) : S(S) {};

	Sphere* S;
	Vector P, N;
	double t; 

	bool intersect(const Ray& r) {
		Vector d = r.O - S->C ; //for simplification
		double delta = sqr(dot(r.u, d)) - ((d).norm2() - sqr(S->R)); // the discriminant of the quadratic equation satisfied by the point of intersection
		if (delta < 0) return false; //no intersection between the line and the sphere is found
		//otherwise one (double) or two intersections are found :
		double tcom = dot(r.u, S->C - r.O);
		double t1 = tcom - sqrt(delta);
		double t2 = tcom + sqrt(delta);

		if (t2 < 0) {
			return false; //no intersection
		}
		//Figure 2.4 for details 
		if (t1 > 0) {
			t = t1;    
		} else {
			t = t2;
		}
		P = r.O + t * r.u; //P is the intersection point 
		N = P - S->C;  //N is the unit normal at P
		N.normalize();
		return true;
	}
};

class BoundingBox{
public:
	Vector min, max;

	// Checks if there is an intersection between the ray and pairs of planes of the bbox
	bool intersect(const Ray& r) const{

		Vector t0 = do_t0(r); // first plane for constant x, y, z
		Vector t1 = do_t1(r); // second plane for constant x, y, z

		double a = std::max(t0[0], std::max(t0[1], t0[2]));
		double b = std::min(t1[0], std::min(t1[1], t1[2]));

		if (b > a){
			return true; // if true, the actual intersection is a
		}
		return false;
	}

	// Computes the intersection along the ray with the first plane of
	// constant x, y, z for repsectively i = 0, 1, 2
	Vector do_t0(const Ray& r){
		Vector t0;
		for (int i=0; i<3; i++){
			t0[i] = (min[i] - r.O[i]) / r.u[i]; 
		} 
		return t0;
	}

	// Computes the intersection along the ray with the first plane of
	// constant x, y, z for repsectively i = 0, 1, 2
	Vector do_t1(const Ray& r){
		Vector t1;
		for (int i=0; i<3; i++){
			t1[i] = (max[i] - r.O[i]) / r.u[i]; 
		} 
		return t1;
	}

};



// triangle class from the lecture notes
// "The .obj file format encodes this structure as an ASCII file. Each line starting with a v defines a
// vertex coordinate (e.g., v 1.0 3.14 0.00, and each line starting with an f defines a face (most often a
// triangle, but it also supports more general polygonal faces – e.g., f 1 2 3 defines a triangle consisting
// of the first 3 vertices, as indexing starts at 1). Negative indices correspond to offsets relative to the end
// of the vertex list. Normal vectors start with a vn, and UV coordinates with vt. The general syntax
// to define a triangle that has normal and UV coordinates is f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3"
 
class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};
 
 
class TriangleMesh : public Geometry {
public:
  ~TriangleMesh() {}
    TriangleMesh(const Vector& albedo, bool isMirror=false, bool isTransparent=false): ::Geometry(albedo, isMirror, isTransparent){};

    std::string getType() const override { return "TriangleMesh"; }

    void readOBJ(const char* obj) {
 
        char matfile[255];
        char grp[255];
 
        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;
 
            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());
 
            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }
 
            if (line[0] == 'v' && line[1] == ' ') {
                Vector vec;
 
                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));
 
                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
 
                } else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;
 
                char* consumedline = line + 1;
                int offset;
 
                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }
 
                consumedline = consumedline + offset;
 
                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;                             
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
 
            }
 
        }
        fclose(f);
 
    }


	void compute_bbox() { //computes the global min and max of the triangle mesh by going through all the vertices to determine the bounding box
		bbox.min = std::accumulate(vertices.begin(), vertices.end(), bbox.min,
									[](const Vector& a, const Vector& b) {
									return Vector(std::min(a[0], b[0]), std::min(a[1], b[1]), std::min(a[2], b[2]));
									});
		bbox.max = std::accumulate(vertices.begin(), vertices.end(), bbox.max,
									[](const Vector& a, const Vector& b) {
									return Vector(std::max(a[0], b[0]), std::max(a[1], b[1]), std::max(a[2], b[2]));
									});
	}



    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
	BoundingBox bbox;
    
};


class IntersectionMesh { //the intersection between a ray and a triangle mesh
public:
    IntersectionMesh(TriangleMesh* mesh) : mesh(mesh) {};

	TriangleMesh* mesh;
	Vector P, N;
	double t; 

	bool intersect(const Ray& r) { //implementation of Möller Trumbore algorithm

		t = 10E9;
		bool has_inter = false;

		for (int i = 0; i < mesh->indices.size(); i++){ // for all the triangles of the mesh

			// we get the three points of the triangle 
			const Vector& A = mesh->vertices[mesh->indices[i].vtxi];
			const Vector& B = mesh->vertices[mesh->indices[i].vtxj];
			const Vector& C = mesh->vertices[mesh->indices[i].vtxk];

			// we get the vectors e1 and e2 
			// and solve the system using Cramer
			Vector e1 = B - A;
			Vector e2 = C - A;
			Vector d1 = A - r.O; 
			Vector d2 = cross(d1, r.u);
			Vector triangle_N = cross(e1, e2); 
			double beta = dot(e2, d2) / dot(r.u, triangle_N);
			double gamma = - dot(e1, d2) / dot(r.u, triangle_N);
			double alpha = 1 - beta - gamma; 
			double triangle_t = dot(d1, triangle_N) / dot(r.u, triangle_N);

			if ((triangle_t < 0) | (beta < 0) | (beta > 1) | (gamma < 0) | (gamma > 1) | (alpha < 0)){ //no intersection for this triangle
				continue;  //start at next iteration = continue checking for intersections for the other triangles of the mesh
			}

			has_inter = true;
			// we retrieve t, the point of intersection P and the unit normal at P
			if (triangle_t < t){ // the relevant intersection is the one w/ the lowest t (closer to the ray origin)
				t = triangle_t;
				N = triangle_N;
				P = r.O + t * r.u;
			}
		}
		N.normalize();
		return has_inter;
	}

};


class Scene{
public:

	std::vector<Geometry*> objects; 
	Vector L; //light source 
	double I;  //light intensity

	bool intersect(const Ray& r, Vector& P, Vector& N, int& objectID, double& bestt ) const {
		bestt = INFINITY;
		bool has_inter = false;

		for (int i=0; i< objects.size(); i++){ //iterating over all the objects of the scene

			std::string t = objects[i]->getType();
			if (t == "Sphere"){
				Sphere* sphere = dynamic_cast<Geometry*>(objects[i]);
				IntersectionSphere Inter(sphere); 
				if (Inter.intersect(r)) {//check if there is an intersection between the ray r and the objects[i]
					if (Inter.t < bestt){ //The relevant intersection is the one w/ the lowest t (closer to the ray origin)
						bestt = Inter.t; 
						P = Inter.P;
						N = Inter.N; //retrieve t, the point of intersection and the unit normal at P
						objectID = i; //retrieve the ID of the object that has been hit by r
						has_inter = true;
					}
				}
			}
			else if (t == "TriangleMesh") {
				TriangleMesh* mesh = dynamic_cast<Geometry*>(objects[i]); //the object is a triangle mesh, we use the class IntersectionMesh
				IntersectionMesh Inter(mesh);	
				if (Inter.intersect(r)) {//check if there is an intersection between the ray r and the objects[i]
					if (Inter.t < bestt){ //The relevant intersection is the one w/ the lowest t (closer to the ray origin)
						bestt = Inter.t; 
						P = Inter.P;
						N = Inter.N; //retrieve t, the point of intersection and the unit normal at P
						objectID = i; //retrieve the ID of the object that has been hit by r
						has_inter = true;
					}
				}
			}
			else {
				std::cout << "scene.intersect : Unknown type in scene " << std::endl;
			}
			
		}
		return has_inter;
	}

	void addSphere(Sphere* s) { //add a sphere to the scene
		objects.push_back(s);
	}

	void addMesh(TriangleMesh* m) {
		objects.push_back(m);
	}

	Vector getColor(const Ray& r, int ray_depth){ //computing the intensity reflected off a surface at point P
		Vector color(0,0,0); //initialize color at black
		if (ray_depth < 0) return color; // decremented at each recursive call to terminate recursion at some point

		double t;
		int objectID;
		Vector P,N; 
		double epsilon = 0.01; //offsetting the starting point of the reflected/refracted/shadow ray off the surface to avoid numerical precision issues

		bool inter = intersect(r, P, N, objectID, t); //check if the ray has intersected an object, if so, retrieve the info
			
		if (inter){

			//Reflections : a mirror transfers light energy from the incident direction to the reflected direction
			if (objects[objectID]->isMirror) { 
				Vector w_reflected = r.u - 2*dot(r.u, N) * N; //reflected direction 
				Ray reflected(P + epsilon * N, w_reflected );
				return getColor(reflected, ray_depth-1);
			}

			//Refractions : rays bounce off surfaces but passing through them
			if (objects[objectID]->isTransparent) {
				
				 // ensuring that n1 < n2 to avoid imaginary results
				double n1 = 1.; 
				double n2 = 1.5;


				Vector correctN = N;
				if (dot(N, r.u) > 0){ //ray exiting the transparent sphere
					// correcting refraction indices and normal sign in this case
					correctN = -N; 
					n1 = 1.5; 
					n2 = 1.; 
				}

				// Fresnel law to randomly launch either a reflection ray, or a refraction ray : 

				double k0 = sqr(n1-n2) / sqr(n1 + n2); //reflection coefficient at normal incidence
				double R = k0 + (1 - k0)*pow((1- std::abs(dot(N, r.u))), 5); // reflection coefficient for incidence r.u
				double t = 1 - R; // transmission coefficient

				double u = uniform(engine); //a uniform random number u between 0 and 1
				if (u < R) {
					//we launch a reflected ray
					Vector w_reflected = r.u - 2*dot(r.u, N) * N; //reflected direction 
					Ray reflected(P + epsilon * N, w_reflected );
					return getColor(reflected, ray_depth-1);
				} //otherwise we launch a refracted ray

				//decomposing the transmitted direction T in tangential Tt and normal component Tn
				Vector Tt = (n1/ n2) * (r.u - dot(r.u, correctN) * correctN);

				double d = 1 - sqr(n1 / n2) * (1 - sqr(dot(r.u, correctN)));
				if (d < 0) {  // recursive call with the internal reflection
					Ray reflected(P + epsilon * correctN, r.u - 2 * dot(r.u, correctN) * correctN);
					return getColor(reflected, ray_depth-1);
				}
				Vector Tn = -correctN * sqrt(d);

				Vector T = Tn + Tt; 

				Ray refracted(P - epsilon * correctN, T); //offsetting the starting point of the refracted ray towards the correct side
				return getColor(refracted, ray_depth-1);
			}


			// Diffuse surfaces

			// direct light : 
			// shading and shadows computation under the Lambertian model

			Vector wlight = L - P;
			double dlight2 = wlight.norm2(); 
			double dlight = sqrt(dlight2);
			wlight.normalize();
			double tshadow;
			Vector Pshadow, Nshadow;
			int objectshadow;
			Ray rShadow(P + epsilon * N, wlight);

			double reach = I/(4 * PI * dlight2) ; //amount of light reaching P 
			double l = reach * std::max(0.0,dot(N, wlight)); // using max for safety measures to avoid noise 
			color = l * objects[objectID]->albedo / PI;
			if (intersect(rShadow, Pshadow, Nshadow, objectshadow, tshadow)){ //setting the visibility term (so the whole color) to 0 if t < dlight
				if (tshadow < dlight){
					color = Vector(0,0,0);
				}
			}

			// indirect light : add a random contribution 
			Vector w_indirect = random_cos(N);
			Ray randomRay(P + epsilon* N , w_indirect);
			color += objects[objectID]->albedo * getColor(randomRay, ray_depth-1);	
		}
		return color;
	}

};

