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
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}



class Ray {
public:
	Ray(Vector O, Vector u) : O(O), u(u) {};
    Vector O, u ; 
};

class Sphere {
public:
	Sphere(const Vector& C, double R, const Vector& albedo, bool isMirror = false, bool isTransparent = false) : C(C), R(R), albedo(albedo), isMirror(isMirror), isTransparent(isTransparent) {};

	Vector C;
	double R;
	Vector albedo;
	bool isMirror, isTransparent;
};

class IntersectionSphere {
public:
    IntersectionSphere(const Sphere& S) : S(S) {};

	Sphere S;
	Vector P, N;
	double t;

	bool intersect(const Ray& r) {
		Vector d = r.O- S.C ; 
		double delta = sqr(dot(r.u, d)) - ((d).norm2() - sqr(S.R));
		if (delta < 0) return false;

		double tcom = dot(r.u,S.C-r.O);
		double t1 = tcom - sqrt(delta);
		double t2 = tcom + sqrt(delta);

		if (t2 < 0) {
			return false;
		}
		
		if (t1 > 0) {
			t = t1;
		} else {
			t = t2;
		}
		P = r.O + t * r.u;
		N = P - S.C;
		N.normalize();
		return true;
	}
};



// class Geometry{
// 	public virtual bool intersection()
// }

class Scene{
public:

	std::vector<Sphere> objects;
	Vector L;
	double I;

	bool intersect(const Ray& r, Vector& P, Vector& N, int& objectID, double& bestt ) const {
		bestt = INFINITY;
		bool has_inter = false;

		for (int i=0; i< objects.size(); i++){
			IntersectionSphere Inter(objects[i]);
			if (Inter.intersect(r)) {
				if (Inter.t < bestt){
					bestt = Inter.t;
					P = Inter.P;
					N = Inter.N;
					objectID = i;
					has_inter = true;
				}
			}
		}
		return has_inter;
	}

	void addSphere(const Sphere &s) {
		objects.push_back(s);
	}

	Vector getColor(const Ray& r, int bounce){
		Vector color(0,0,0);
		if (bounce < 0) return color;

		double t;
		int objectID;
		Vector P,N;
		bool inter = intersect(r, P, N, objectID, t);
			
		if (inter){

			if (objects[objectID].isMirror) {
				Ray reflected(P+0.001 * N, r.u - 2*dot(r.u, N) * N);
				return getColor(reflected, bounce-1);
			}
			if (objects[objectID].isTransparent) {
				double n1 = 1;
				double n2 = 1.5;

				Vector correctN = N;
				if (dot(N, r.u) > 0){
					correctN = -N; 
					std::swap(n1, n2);
				}

				Vector Tt = n1/ n2 * (r.u - dot(r.u, correctN) * correctN);
				double d = 1 - sqr(n1 / n2) * (1 - sqr(dot(r.u, correctN)));
				if (d < 0) {
					Ray reflected(P+0.001 * correctN, r.u - 2 * dot(r.u, correctN) * correctN);
					return getColor(reflected, bounce-1);
				}

				Vector Tn = -sqrt(d) * correctN;
				Vector T = Tn + Tt;

				Ray refracted(P - 0.001 * correctN, T);
				return getColor(refracted, bounce-1);
			}

			Vector wlight = L - P;
			double dlight2 = wlight.norm2();
			wlight.normalize();
			double tshadow;
			Vector Pshadow, Nshadow;
			int objectshadow;
			Ray rShadow(P + 0.001 * N, wlight);

			double l = I/(4 * PI * dlight2) * std::max(0.0,dot(N, wlight));
			color = l * objects[objectID].albedo / PI;
			if (intersect(rShadow, Pshadow, Nshadow, objectshadow, tshadow)){
				if (sqr(tshadow) < dlight2){
					color = Vector(0,0,0);
				}
			}
		}
		return color;
	}

};
