#include "Ray.h"
#include "Geometry.h"
#include <glm/glm.hpp>
#include <vector>
#include <memory>
class Scene
{
public:
	Scene() {
		std::shared_ptr<Geometry> sp1=std::make_shared<Sphere>(Sphere{ glm::vec3{0.0f,0.0f,-2.0f},1.0f });
		geolist.push_back(sp1);
	}

	std::vector<std::shared_ptr<Geometry>> geolist{};
};

