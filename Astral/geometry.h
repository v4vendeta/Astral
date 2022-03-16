#pragma  once
#include "Intersection.h"
#include "ray.h"
#include "bvh.h"
struct Geo {
	virtual IntersectInfo intersect(const Ray& r) const = 0;
	aabb bbox;
};
struct Sphere :public Geo
{

	Sphere(double _r, vec3 _p, vec3 _e, vec3 _c, Material _mat) : rad(_r), p(_p), e(_e), c(_c), mat(_mat) {}

	IntersectInfo intersect(const Ray& r) const override
	{
		IntersectInfo inter;
		vec3 op = p - r.o;
		double t, eps = 1e-4, b = op * r.d, det = b * b - op * op + rad * rad;
		if (det < 0)
			inter.intersected = false;
		else
		{
			inter.intersected = true;
			det = sqrt(det);
			inter.t = (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
			inter.pos = r.At(inter.t);
			inter.normal = normalize(r.At(inter.t) - p);
			inter.material = mat;
			inter.color = c;
			inter.e = e;
		}

		return inter;
	}
	double rad;	  // radius
	vec3 p, e, c; // position, emission, color
	Material mat;  // reflection type (DIFFuse, SPECular, REFRactive)
};

struct Triangle :public Geo {

	Triangle(vec3 _v1, vec3 _v2, vec3 _v3, vec3 _e, vec3 _c, Material _mat) :v1(_v1), v2(_v2), v3(_v3), c(_c), e(_e), mat(_mat) {}

	IntersectInfo intersect(const Ray& r) const override {

		IntersectInfo inter;
		vec3 E1 = v2 - v1;
		vec3 E2 = v3 - v1;
		vec3 S = r.o - v1;

		vec3 S1 = cross(r.d, E2);
		vec3 S2 = cross(S, E1);
		float coeff = 1.0 / (S1 * E1);
		float t = coeff * (S2 * E2);
		float b1 = coeff * (S1 * S);
		float b2 = coeff * (S2 * r.d);
		float u = 0.0f, v = 0.0f;
		if (t >= 0 && b1 >= 0 && b2 >= 0 && (1 - b1 - b2) >= 0)
		{
			inter.intersected = true;
			inter.pos = r.At(t);
			inter.t = t;
			inter.normal = cross(E1, E2);
			inter.material = mat;
			inter.color = c;
			inter.e = e;
		}
		else {
			inter.intersected = false;
		}

		return inter;
	}
	vec3 v1, v2, v3;
	vec3 c, e;
	Material mat;
};

struct Quad :public Geo {
	Quad(vec3 ld,vec3 rd,vec3 rt,vec3 lt, vec3 _e, vec3 _c, Material _mat):t1(ld,rd,rt,_e,_c,_mat),t2(ld,rt,lt,_e, _c, _mat){}
	//Quad(Triangle _t1,Triangle _t2):t1(_t1),t2(_t2){}
	IntersectInfo intersect(const Ray& r) const override {
		IntersectInfo inter1 = t1.intersect(r);
		IntersectInfo inter2 = t2.intersect(r);
		if (inter1.intersected)return inter1;
		else if (inter2.intersected)return inter2;

		else return inter1;// not intersected
	}

	Triangle t1, t2;
};

struct Mesh :public Geo {
	std::vector<Triangle> meshes;
};