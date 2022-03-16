#pragma  once
#include "vector.h"
struct Ray
{
	vec3 o, d;
	Ray(vec3 o_, vec3 d_) : o(o_), d(d_) {}
	vec3 At(double t)const {
		vec3 ret = o + t * d;
		return ret;
	}
};