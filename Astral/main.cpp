#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
#include <glm/glm.hpp>
#include <iostream>
#include <random>
#include <omp.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#include <vector>

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1

inline double random() {
	return dis(gen);
}
struct Vec
{                 // Usage: time ./smallpt 5000 && xv image.ppm
	double x, y, z; // position, also color (r,g,b)
	Vec(double x_ = 0, double y_ = 0, double z_ = 0)
	{
		x = x_;
		y = y_;
		z = z_;
	}
	Vec operator+(const Vec& b) const { return Vec(x + b.x, y + b.y, z + b.z); }
	Vec operator-(const Vec& b) const { return Vec(x - b.x, y - b.y, z - b.z); }
	Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
	Vec mult(const Vec& b) const { return Vec(x * b.x, y * b.y, z * b.z); }
	Vec& norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
	double dot(const Vec& b) const { return x * b.x + y * b.y + z * b.z; } // cross:
	Vec operator%(Vec& b) { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
};
struct Ray
{
	Vec o, d;
	Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};
enum Refl_t
{
	DIFF,
	SPEC,
	REFR
}; // material types, used in radiance()
struct Sphere
{
	double rad;  // radius
	Vec p, e, c; // position, emission, color
	Refl_t refl; // reflection type (DIFFuse, SPECular, REFRactive)
	Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) : rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
	double intersect(const Ray& r) const
	{                   // returns distance, 0 if nohit
		Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
		double t, eps = 1e-4, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
		if (det < 0)
			return 0;
		else
			det = sqrt(det);
		return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
	}
};


Sphere spheres[] = {
	//Scene: radius, position, emission, color, material
	Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.65, .05, .05), DIFF),   //Left
	Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.12, .45, .15), DIFF), //Rght
	Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.73, .73, .73), DIFF),         //Back
	Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(.65, .05, .05), SPEC),               //Frnt
	Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.73, .73, .73), DIFF),         //Botm
	Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.73, .73, .73), DIFF), //Top
	Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1) * .999, SPEC),        //Mirr
	Sphere(16.5, Vec(73, 16.5, 78), Vec(), Vec(1, 1, 1) * .999, REFR),        //Glas
	Sphere(600, Vec(50, 681.6 - .27, 81.6), Vec(12, 12, 12), Vec(), DIFF)     //Lite
};
inline double clamp(double x) {
	return x < 0 ? 0 : x > 1 ? 1: x;
}
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }
inline bool intersect(const Ray& r, double& t, int& id)
{
	double n = sizeof(spheres) / sizeof(Sphere), d, inf = t = 1e20;
	for (int i = int(n); i--;)
		if ((d = spheres[i].intersect(r)) && d < t)
		{
			t = d;
			id = i;
		}
	return t < inf;
}
Vec radiance(const Ray& r, int depth)
{
	double t;   // distance to intersection
	int id = 0; // id of intersected object
	if (!intersect(r, t, id))
		return Vec();                  // if miss, return black
	const Sphere& obj = spheres[id]; // the hit object
	Vec x = r.o + r.d * t, n = (x - obj.p).norm(), nl = n.dot(r.d) < 0 ? n : n * -1, f = obj.c;
	double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y
		: f.z; // max refl
	if (++depth > 5)
		if (random() < p)
			f = f * (1 / p);
		else
			return obj.e; //R.R.
	if (obj.refl == DIFF)
	{ // Ideal DIFFUSE reflection
		double r1 = 2 * 3.1415926 * random(), r2 = random(), r2s = sqrt(r2);
		Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(), v = w % u;
		Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
		return obj.e + f.mult(radiance(Ray(x, d), depth));
	}
	else if (obj.refl == SPEC) // Ideal SPECULAR reflection
		return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth));
	Ray reflRay(x, r.d - n * 2 * n.dot(r.d)); // Ideal dielectric REFRACTION
	bool into = n.dot(nl) > 0;                // Ray from outside going in?
	double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;
	if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // Total internal reflection
		return obj.e + f.mult(radiance(reflRay, depth));
	Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
	double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir.dot(n));
	double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
	return obj.e + f.mult(depth > 2 ? (random() < P ? // Russian roulette
		radiance(reflRay, depth) * RP
		: radiance(Ray(x, tdir), depth) * TP)
		: radiance(reflRay, depth) * Re + radiance(Ray(x, tdir), depth) * Tr);
}
int main()
{


	int w = 300, h = 300;
	int samps = 3;
	Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm());        // cam pos, dir
	Vec cx = Vec(w * .5135 / h), cy = (cx % cam.d).norm() * .5135, r;
	std::vector<Vec> c;
	c.resize(w * h);
//#pragma omp parallel for schedule(dynamic, 1) private(r) // OpenMP
	for (int x = 0; x < w; x++) {

		for (int y = 0; y < h; y++) {
			fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps, 100. * x / (w - 1));
			r = Vec();
			for (int s = 0; s < samps; s++) {
				
					double r1 = 2 * random(), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
					double r2 = 2 * random(), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
					Vec d = cx * (((x + .5 + dx) / 2 + x) / w - .5) + cy * (((y + .5 + dy) / 2 + y) / h - .5) + cam.d;
					float dirx{ 2.0f * static_cast<float>(x) / w - 1.0f };
					float diry{ 2.0f * static_cast<float>(y) / h - 1.0f };
					Ray ray{ Vec {0,0,0},Vec{dirx * w / h,diry,-1.0f}};
					d = Vec{ dirx * w / h,diry,-1.0f };
					r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0) * (1. / samps);
					
			} // Camera rays are pushed ^^^^^ forward to start in interior
			c[y * w + x] = c[y * w + x] + Vec(clamp(r.x), clamp(r.y), clamp(r.z));
				
			
		}
	}
	//for (int y = 0; y < h; y++)
	//{ // Loop over image rows
	//	fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (h - 1));
	//	for (unsigned short x = 0; x < w; x++) // Loop cols
	//		for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++)         // 2x2 subpixel rows
	//			for (int sx = 0; sx < 2; sx++, r = Vec())
	//			{ // 2x2 subpixel cols
	//				for (int s = 0; s < samps; s++)
	//				{
	//					double r1 = 2 * random(), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
	//					double r2 = 2 * random(), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
	//					Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
	//					r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0) * (1. / samps);
	//				} // Camera rays are pushed ^^^^^ forward to start in interior
	//				c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
	//			}
	//}
	std::vector<uint8_t> pixel_color;
	pixel_color.resize(w * h * 3);
	for (int32_t x{ 0 }; x < w; x++) {
		for (int32_t y{ 0 }; y < h; y++) {
			pixel_color[y * w * 3 + x * 3 + 0] = toInt(c[y*w+x].x);
			pixel_color[y * w * 3 + x * 3 + 1] = toInt(c[y*w+x].y);
			pixel_color[y * w * 3 + x * 3 + 2] = toInt(c[y*w+x].z);
		}
	}
	stbi_write_bmp("../output/img.bmp", w, h, 3, pixel_color.data());

}

//for (int y = 0; y < h; y++)
//{ // Loop over image rows
//	fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (h - 1));
//	for (unsigned short x = 0; x < w; x++) // Loop cols
//		for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++)         // 2x2 subpixel rows
//			for (int sx = 0; sx < 2; sx++, r = Vec())
//			{ // 2x2 subpixel cols
//				for (int s = 0; s < samps; s++)
//				{
//					double r1 = 2 * random(), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
//					double r2 = 2 * random(), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
//					Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) + cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
//					r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0) * (1. / samps);
//				} // Camera rays are pushed ^^^^^ forward to start in interior
//				c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
//			}
//}