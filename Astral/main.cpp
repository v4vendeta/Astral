#include <cstdlib> 
#include <cstdio>
#include <omp.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#include <vector>
#include "utils.h"
int main()
{

	int w = 300, h = 300, channels = 3;
	int spp = 500;
	//Ray cam(vec3(50, 52, 295.6), normalize(vec3(0, -0.042612, -1)));
	Ray cam(vec3(0, 0, 300), normalize(vec3(0, -0.042612, -1)));
	double fov = 0.5135;
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
				vec3 d = normalize(cx * ((dxy.x + x) / w - .5) + cy * ((dxy.y + y) / h - .5) + cam.d);
				r = r + radiance(Ray(cam.o + d * 140, d), 0) * (1. / spp);
			}
			pixel_color[(h-1-y) * w * channels + x * channels + 0] = toInt(r.x);
			pixel_color[(h-1-y) * w * channels + x * channels + 1] = toInt(r.y);
			pixel_color[(h-1-y) * w * channels + x * channels + 2] = toInt(r.z);
		}
	}
	stbi_write_bmp("../output/img.bmp", w, h, 3, pixel_color.data());
}
