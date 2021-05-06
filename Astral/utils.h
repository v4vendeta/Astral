#pragma  once
#include <cmath>
#include <iostream>
#include <random>
#include <omp.h>

const double PI = 3.141592653589793238463;
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dis(0, 1);	//uniform distribution between 0 and 1
std::uniform_real_distribution<double> dis2(-1, 1); //uniform distribution between 0 and 1

struct vec2
{
	double x, y;
	vec2(double x_ = 0, double y_ = 0)
	{
		x = x_;
		y = y_;
	}
};
struct vec3
{					// Usage: time ./smallpt 5000 && xv image.ppm
	double x, y, z; // position, also color (r,g,b)
	vec3(double x_ = 0, double y_ = 0, double z_ = 0)
	{
		x = x_;
		y = y_;
		z = z_;
	}

};
vec3 operator+(const vec3& lhs, const vec3& rhs) {
	return vec3(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
}
vec3 operator-(const vec3& lhs, const vec3& rhs) {
	return vec3(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
}
vec3 operator*(const vec3& lhs, double rhs) {
	return vec3(rhs * lhs.x, rhs * lhs.y, rhs * lhs.z);
}
vec3 operator*(double rhs, const vec3& lhs) {
	return vec3(rhs * lhs.x, rhs * lhs.y, rhs * lhs.z);
}
double operator*(const vec3& lhs, const vec3& rhs)
{
	return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}
vec3 operator/(const vec3& lhs, double rhs) {
	return vec3(lhs.x / rhs, lhs.y / rhs, lhs.z / rhs);
}
vec3 operator/(double rhs, const vec3& lhs) {
	return vec3(lhs.x / rhs, lhs.y / rhs, lhs.z / rhs);
}
vec3 normalize(vec3 lhs)
{
	return lhs * (1 / sqrt(lhs.x * lhs.x + lhs.y * lhs.y + lhs.z * lhs.z));
}
vec3 cross(const vec3& lhs, const vec3& rhs)
{
	return vec3(lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x);
}
vec3 multiply(const vec3& lhs, const vec3& rhs) {
	return vec3(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z);
}
inline double random()
{
	return dis(gen);
}
// generate a sample point in unit sphere
inline vec2 random2()
{
	vec2 ret;
	do
	{
		ret.x = dis2(gen);
		ret.y = dis2(gen);
	} while (ret.x * ret.x + ret.y * ret.y > 1.0);
	return ret;
}
inline double clamp(double x){return x < 0 ? 0 : x > 1 ? 1 : x;}
// gamma correction
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }


struct Ray
{
	vec3 o, d;
	Ray(vec3 o_, vec3 d_) : o(o_), d(d_) {}
	vec3 At(double t)const {
		vec3 ret = o + t * d;
		return ret;
	}
};
enum class Mat
{
	DIFF,
	SPEC,
	REFR
}; // material types, used in radiance()
struct IntersectInfo {
	double t = 0.0;
	bool intersected = false;
	vec3 pos;
	vec3 normal;
	vec3 color;
	vec3 e;
	Mat material = Mat::DIFF;

};
struct Geo {
	virtual IntersectInfo intersect(const Ray& r) const=0;
};
struct Sphere:public Geo
{

	Sphere(double _r, vec3 _p, vec3 _e, vec3 _c, Mat _mat) : rad(_r), p(_p), e(_e), c(_c), mat(_mat) {}

	IntersectInfo intersect(const Ray& r) const override
	{
		IntersectInfo inter;
		vec3 op = p - r.o;
		double t, eps = 1e-4, b = op * r.d, det = b * b - op * op + rad * rad;
		if (det < 0)
			inter.intersected = false;
		else
		{
			inter.intersected = true;
			det = sqrt(det);
			inter.t = (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
			inter.pos = r.At(inter.t);
			inter.normal = normalize(r.At(inter.t) - p);
			inter.material = mat;
			inter.color = c;
			inter.e = e;
		}

		return inter;
	}
	double rad;	  // radius
	vec3 p, e, c; // position, emission, color
	Mat mat;  // reflection type (DIFFuse, SPECular, REFRactive)
};

struct Triangle:public Geo {

	Triangle(vec3 _v1, vec3 _v2, vec3 _v3, vec3 _e, vec3 _c, Mat _mat) :v1(_v1), v2(_v2), v3(_v3), c(_c), e(_e), mat(_mat) {}

	IntersectInfo intersect(const Ray& r) const override {

		IntersectInfo inter;
		vec3 E1 = v2 - v1;
		vec3 E2 = v3 - v1;
		vec3 S = r.o - v1;

		vec3 S1 = cross(r.d, E2);
		vec3 S2 = cross(S, E1);
		float coeff = 1.0 /(S1* E1) ;
		float t = coeff * (S2*E2);
		float b1 = coeff * (S1* S);
		float b2 = coeff * (S2* r.d);
		float u = 0.0f, v = 0.0f;
		if (t >= 0 && b1 >= 0 && b2 >= 0 && (1 - b1 - b2) >= 0)
		{
			inter.intersected = true;
			inter.pos = r.At(t);
			inter.t = t;
			inter.normal = cross(E1, E2);
			inter.material = mat;
			inter.color = c;
			inter.e = e;
		}
		else {
			inter.intersected = false; 
		}

		return inter;
	}
	vec3 v1, v2, v3;
	vec3 c, e;
	Mat mat;
};
// 
std::vector<Geo*> scene = {
	new Sphere(1e5, vec3(1e5 - 50 + 1, 40.8 - 52, 81.6), vec3(), vec3(.65, .05, .05), Mat::DIFF),	 //Left
	new Sphere(1e5, vec3(-1e5 + 99 - 50, 40.8 - 52, 81.6), vec3(), vec3(.12, .45, .15), Mat::DIFF), //Rght
	new Sphere(1e5, vec3(50 - 50, 40.8 - 52, 1e5), vec3(), vec3(.73, .73, .73), Mat::DIFF),		 //Back
	new Sphere(1e5, vec3(50 - 50, 1e5 - 52, 81.6), vec3(), vec3(.73, .73, .73), Mat::DIFF),		 //Botm
	new Sphere(1e5, vec3(50 - 50, -1e5 - 52 + 81.6, 81.6), vec3(), vec3(.73, .73, .73), Mat::DIFF), //Top
	new Sphere(16.5, vec3(27 - 50, 16.5 - 52, 47), vec3(), vec3(1, 1, 1) * .999, Mat::SPEC),		 //Mirr
	new Sphere(16.5, vec3(73 - 50, 16.5 - 52, 78), vec3(), vec3(1, 1, 1) * .999, Mat::REFR),		 //Glas
	new Sphere(600, vec3(50 - 50, 681.6 - 52 - .27, 81.6), vec3(12, 12, 12), vec3(), Mat::DIFF),	 //Lite
	//new Sphere(20, vec3(50 - 50, 40.8 - 52, 90), vec3(), vec3(0.75), Mat::DIFF),
	new Triangle(vec3(0,-40,100),vec3(0,-30,90),vec3(-30,-40,90),vec3(), vec3(.65, .75, .05) * .999, Mat::DIFF),
	//new Triangle(vec3(0,40,100),vec3(0,50,100),vec3(-50,50,50),vec3(), vec3(.65, .75, .05) * .999, Mat::DIFF)
};


// get the cloeset hit point and hit info
inline IntersectInfo intersect(const Ray& r)
{
	IntersectInfo ret;
	double inf = 1e20;// inf
	for (int i = 0; i < scene.size(); i++) {
		IntersectInfo inter = scene[i]->intersect(r);
		if (inter.intersected == true && inter.t > 0 && inter.t < inf)
		{
			ret.t = inter.t;
			inf = inter.t; // reset the closest object
			ret.material = inter.material;
			ret.normal = inter.normal;
			ret.pos = inter.pos;
			ret.color = inter.color;
			ret.e = inter.e;
			ret.intersected = inter.intersected;
		}
	}
	return ret;

}
vec3 radiance(const Ray& r, int depth)
{
	// max depth
	if (depth > 20) {
		return vec3();
	}

	IntersectInfo info = intersect(r);
	if (info.intersected == false)
		return vec3();				 // if miss, return black7689

	vec3 x = info.pos;			 //hit point
	vec3 n = info.normal; // hit point normal dir
	vec3 nl = n * r.d < 0 ? n : n * -1;
	vec3 f = info.color;

	double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y: f.z; // max refl
	if (++depth > 5)
		if (random() < p)
			f = f * (1 / p);
		else
			return info.e; //R.R.

	if (info.material == Mat::DIFF)
	{ // Ideal DIFFUSE reflection
		double r1 = 2 * PI * random(), r2 = random(), r2s = sqrt(r2);
		vec3 w = nl;
		vec3 u = normalize(cross(fabs(w.x) > .1 ? vec3(0, 1) : vec3(1), w));
		vec3 v = cross(w, u);
		vec3 d = normalize(u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2));
		return info.e + multiply(f, radiance(Ray(x, d), depth));
	}
	else if (info.material == Mat::SPEC) // Ideal SPECULAR reflection
		return info.e + multiply(f, radiance(Ray(x, r.d - n * 2 * (n * r.d)), depth));

	Ray reflRay(x, r.d - n * 2 * n * r.d); // Ideal dielectric REFRACTION
	bool into = n * nl > 0;			   // Ray from outside going in?
	double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc;
	double ddn = r.d * nl, cos2t;
	if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // Total internal reflection
		return info.e + multiply(f, radiance(reflRay, depth));
	vec3 tdir = normalize(r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t))));
	double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir * n);
	double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
	return info.e + multiply(f, depth > 2 ? (random() < P ? // Russian roulette
		radiance(reflRay, depth) * RP
		: radiance(Ray(x, tdir), depth) * TP)
		: radiance(reflRay, depth) * Re + radiance(Ray(x, tdir), depth) * Tr);
}