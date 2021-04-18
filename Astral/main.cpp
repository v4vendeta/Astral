#include <iostream>
#include <string>
#include <glm/glm.hpp>
#include <stdexcept>
#include <cstdint>
#include <array>
#include <filesystem>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#include "Scene.h"
#include "Ray.h"
const uint16_t width{ 300 }, height{ 200 }, channels{ 3 }, spp{ 1 };


glm::vec3 CastRay(const Ray& ray,const  Scene& scene) {
	
	glm::vec3 ret{ 0.0 };
	for (const auto& i : scene.geolist) {
		// intersected
		Intersect_Info inter = i->Intersection(ray);
		if (inter.happened == true) {
			ret = glm::vec3{ 1.0, 1.0, 1.0 };
		}
	}
	return ret;
}
int main(int argc, char** argv) {
	try {

		std::filesystem::path outpath{ "../output/ouput.bmp" };
		std::vector<uint8_t> pixel_color{};
		pixel_color.resize(width * height * channels);
		glm::vec3 accumulate_color{0.0f};
		Scene m_scene;
		for (uint16_t x{ 0 }; x < width; x++) {
			fprintf(stderr, "\rRendering (%d spp) %5.2f%%", spp, 100. * x / (width - 1));
			for (uint16_t y{ 0 }; y < height; y++) {
				// clear accumulate_color
				accumulate_color = glm::vec3{ 0.0f };
				for (uint8_t i{ 0 }; i < spp; i++) {
					//calculate the ray from the eye to the pixel
					float dirx{ 2.0f*static_cast<float>(x) / width-1.0f };
					float diry{ 2.0f*static_cast<float>(y) / height-1.0f };
					//std::cout << "(" << dirx << "," << diry << ")\n";
					// ray origin from (0,0,0) 
					Ray ray{ glm::vec3{0,0,0},glm::normalize(glm::vec3{dirx*width/height,diry,-1.0f}),0.0f };
					accumulate_color += CastRay(ray,m_scene);;
				}
				accumulate_color /= spp;
				
				pixel_color[y * width * channels + x * channels + 0] = 255*accumulate_color.r;
				pixel_color[y * width * channels + x * channels + 1] = 255*accumulate_color.g;
				pixel_color[y * width * channels + x * channels + 2] = 255*accumulate_color.b;
			}
		}

		stbi_write_bmp(outpath.string().c_str(), width, height, channels, pixel_color.data());
	}
	catch (std::exception& e) {
		std::cout << e.what() << "\n";
	}
	return 0;

}
