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

class Sphere {
public:
	Sphere(const Vector& C, double R, const Vector& albedo, bool isMirror = false, bool isTransparent = false) : C(C), R(R), albedo(albedo), isMirror(isMirror), isTransparent(isTransparent) {};

	Vector C; //center
	double R; //radius
	Vector albedo;  //color
	bool isMirror, isTransparent; //options
};

class IntersectionSphere { //the intersection between a ray and a sphere
public:
    IntersectionSphere(const Sphere& S) : S(S) {};

	Sphere S;
	Vector P, N;
	double t; 

	bool intersect(const Ray& r) {
		Vector d = r.O- S.C ; //for simplification
		double delta = sqr(dot(r.u, d)) - ((d).norm2() - sqr(S.R)); // the discriminant of the quadratic equation satisfied by the point of intersection
		if (delta < 0) return false; //no intersection between the line and the sphere is found
		//otherwise one (double) or two intersections are found :
		double tcom = dot(r.u,S.C-r.O);
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
		N = P - S.C;  //N is the unit normal at P
		N.normalize();
		return true;
	}
};



class Scene{
public:

	std::vector<Sphere> objects; 
	Vector L; //light source 
	double I;  //light intensity

	bool intersect(const Ray& r, Vector& P, Vector& N, int& objectID, double& bestt ) const {
		bestt = INFINITY;
		bool has_inter = false;

		for (int i=0; i< objects.size(); i++){ //iterating over all the objects of the scene
			IntersectionSphere Inter(objects[i]); 
			if (Inter.intersect(r)) {//check if there is an intersection beteween r and the sphere at objects[i]
				if (Inter.t < bestt){ //The relevant intersection is the one w/ the lowest t (closer to the ray origin)
					bestt = Inter.t; 
					P = Inter.P;
					N = Inter.N; //retrieve t, the point of intersection and the unit normal at P
					objectID = i; //retrieve the ID of the object that has been hit by r
					has_inter = true;
				}
			}
		}
		return has_inter;
	}

	void addSphere(const Sphere &s) { //add a sphere to the scene
		objects.push_back(s);
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
			if (objects[objectID].isMirror) { 
				Vector w_reflected = r.u - 2*dot(r.u, N) * N; //reflected direction 
				Ray reflected(P + epsilon * N, w_reflected );
				return getColor(reflected, ray_depth-1);
			}

			//Refractions : rays bounce off surfaces but passing through them
			if (objects[objectID].isTransparent) {
						
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
			color = l * objects[objectID].albedo / PI;
			if (intersect(rShadow, Pshadow, Nshadow, objectshadow, tshadow)){ //setting the visibility term (so the whole color) to 0 if t < dlight
				if (tshadow < dlight){
					color = Vector(0,0,0);
				}
			}

			// indirect light : add a random contribution 
			Vector w_indirect = random_cos(N);
			Ray randomRay(P + epsilon* N , w_indirect);
			color += objects[objectID].albedo * getColor(randomRay, ray_depth-1);	
		}
		return color;
	
	}

};

