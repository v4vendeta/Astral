#pragma  once
#include <glm/glm.hpp>
#include "Ray.h"

class Intersect_Info {
	
public:
	Intersect_Info() = default;
	bool happened = false;
	glm::vec3 hit_pos{};
	glm::vec3 normal{};
};
class Geometry
{
public:
	Geometry(){}
	~Geometry(){}
	virtual Intersect_Info Intersection(const Ray& r)=0;
private:

};

class Sphere :public Geometry {
public:
	Sphere(glm::vec3 p, float r) :position{ p }, rad{ r } {}
	Intersect_Info Intersection(const Ray&r) override{
		Intersect_Info inter;
		glm::vec3 op = position - r.origin; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
		float t, eps = 1e-4;
		
		float b = glm::dot(op,r.direction);
		float det = b * b - glm::dot(op,op) + rad * rad;
		if (det < 0)
			return inter;
		else
			det = sqrt(det);
		inter.happened = (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
		return inter;
	}

	glm::vec3 position; 
	float rad;
};

class Triangle :public Geometry {
	Triangle(){}
};