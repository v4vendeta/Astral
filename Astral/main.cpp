#include <cstdlib> 
#include <cstdio>
#include <omp.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#include <vector>
#include "scene.h"
#include "Intersection.h"
#include "ray.h"
#include<cmath>
#include <iostream>
#include <memory>
// get the cloeset hit point and hit info
inline IntersectInfo intersect(const Ray& r, std::vector<std::shared_ptr<Geo>>)
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


	IntersectInfo info = intersect(r,scene);
	if (info.intersected == false)
		return vec3();				 // if miss, return black7689

	vec3 x = info.pos;			 //hit point
	vec3 n = info.normal;        // hit point normal dir
	vec3 f = info.color;         // hit color
	vec3 nl = n * r.d < 0 ? n : n * -1; //orienting normal

	//Russian roulette
	double p = 0.8; // 
	
	if (depth++ > 5)
		if (random() > p)
			return info.e;
		else
			f = f / p;



	if (info.material == Material::DIFF)
	{
		// shoot a ray to random dir
		double r1 = 2 * PI * random(), r2 = random(), r2s = sqrt(r2);
		vec3 w = nl;
		// tangent space
		vec3 u = normalize(cross(fabs(w.x()) > 0.1 ? vec3(0, 1) : vec3(1), w));
		vec3 v = cross(w, u);
		vec3 d = normalize(u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)); //random reflect ray
		return info.e + multiply(f, radiance(Ray(x, d), depth));
	}
	// reflect
	else if (info.material == Material::SPEC)
	{
		// shoot a ray to reflect dir
		vec3 reflect_dir = r.d - n * 2 * (n * r.d);
		return info.e + multiply(f, radiance(Ray(x, reflect_dir), depth));
	}
}

int main()
{

	int w = 300, h = 300, channels = 3;
	int spp = 100;
	Ray cam(vec3(0, 0, 300), normalize(vec3(0, -0.042612, -1)));
	double fov = deg2rad(30.0);

	vec3 cx = vec3(w / h) * fov;

	vec3 cy = normalize(cross(cx,cam.d)) * fov; // axis of camera coordinate
	vec3 r;
	std::vector<uint8_t> pixel_color(w * h * channels);
#pragma omp parallel for schedule(dynamic, 1) private(r) // OpenMP
	for (int y = 0; y < h; y++)
	{
		fprintf(stderr, "\rRendering (%d spp) %5.2f%%", spp, 100. * y / (h - 1));
		
		for (int x = 0; x < w; x++)
		{
			r = vec3();
			for (int s = 0; s < spp; s++)
			{
				vec2 dxy = random2(); // sample sphere
				vec3 d = normalize(cx * ((dxy.x() + x) / w - .5) + cy * ((dxy.y() + y) / h - .5) + cam.d);
				r = r + radiance(Ray(cam.o + d * 140, d), 0) * (1.0 / spp);
			}
			pixel_color[(h-1-y) * w * channels + x * channels + 0] = toInt(r.x());
			pixel_color[(h-1-y) * w * channels + x * channels + 1] = toInt(r.y());
			pixel_color[(h-1-y) * w * channels + x * channels + 2] = toInt(r.z());
		}
	}
	std::cout << cx.x() << cx.y();
	stbi_write_bmp("../output/img.bmp", w, h, 3, pixel_color.data());
}
