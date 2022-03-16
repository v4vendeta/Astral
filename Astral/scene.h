#pragma  once
#include "geometry.h"
#include <vector>
#include <memory>
std::vector<std::shared_ptr<Geo>> scene = {
	std::make_shared<Sphere>(1e5, vec3(1e5 - 50 + 1, 40.8 - 52, 81.6), vec3(), vec3(.65, .05, .05), Material::DIFF),	 //Left
	std::make_shared<Sphere>(1e5, vec3(-1e5 + 99 - 50, 40.8 - 52, 81.6), vec3(), vec3(.12, .45, .15), Material::DIFF), //Rght
	std::make_shared<Sphere>(1e5, vec3(50 - 50, 40.8 - 52, 1e5), vec3(), vec3(.73, .73, .73), Material::DIFF),		 //Back
	std::make_shared<Sphere>(1e5, vec3(50 - 50, 1e5 - 52, 81.6), vec3(), vec3(.73, .73, .73), Material::DIFF),		 //Botm
	std::make_shared<Sphere>(1e5, vec3(50 - 50, -1e5 - 52 + 81.6, 81.6), vec3(), vec3(.73, .73, .73), Material::DIFF), //Top
	std::make_shared<Sphere>(16.5, vec3(27 - 50, 16.5 - 52, 47), vec3(), vec3(1, 1, 1) * .999, Material::SPEC),		 //Mirr
	std::make_shared<Sphere>(16.5, vec3(73 - 50, 16.5 - 52, 78), vec3(), vec3(1, 1, 1) * .999, Material::DIFF),		 //Glas
	//new Sphere(600, vec3(50 - 50, 681.6 - 52 - .27, 81.6), vec3(12, 12, 12), vec3(), Mat::DIFF),	 //light
	std::make_shared<Triangle>(vec3(-20,28,90),vec3(20,28,90),vec3(20,28,80),vec3(12, 12, 12), vec3(), Material::DIFF),	 //light
	std::make_shared<Triangle>(vec3(-20,28,90),vec3(20,28,80),vec3(-20,28,80),vec3(12, 12, 12), vec3(), Material::DIFF),	 //light
	//new Sphere(20, vec3(50 - 50, 40.8 - 52, 90), vec3(), vec3(0.75), Mat::DIFF),
	std::make_shared<Triangle>(vec3(0,-40,100),vec3(0,-30,90),vec3(-30,-40,90),vec3(), vec3(.65, .75, .05) * .999, Material::DIFF),
	//std::make_shared<Quad>(vec3(-40,-60,40),vec3(0,-60,40),vec3(0,-50,-40),vec3(-40,-50,-40),vec3(), vec3(.65, .75, .05) * .999, Mat::SPEC)
};


