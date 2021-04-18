#pragma once
#include <glm/glm.hpp>

class Ray
{
public:
	Ray(glm::vec3 o, glm::vec3 d, float t = 0.0f) :origin{ o }, direction{ d }, time{ t } {}
//	glm::vec3 Origin() const { return origin; }
//	glm::vec3 Direction() const { return direction; }
//	glm::vec3 Pos(float t) const { return  origin + t * direction; }
//	float Time() const { return time; }
	glm::vec3 At(float t) const { return origin + direction * t; }
//private:
	glm::vec3 origin{0.0f}, direction{0.0f};
	float time{0.0f};
};