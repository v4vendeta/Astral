#pragma  once
#include <vector>
#include "vector.h"
#include "material.h"
#include "ray.h"
struct IntersectInfo {
	double t = 0.0;
	bool intersected = false;
	vec3 pos;
	vec3 normal;
	vec3 color;
	vec3 e;
	Material material = Material::DIFF;

};

