#pragma  once
#include <random>
#include <cmath>
const double PI = 3.141592653589793238463;
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dis(0, 1);	//uniform distribution between 0 and 1
std::uniform_real_distribution<double> dis2(-1, 1); //uniform distribution between 0 and 1

struct vec2
{
	//double e[0], e[1];
	double e[2];
	vec2(double x_ = 0, double y_ = 0)
	{
		e[0] = x_;
		e[1] = y_;
	}
	double x() const {
		return e[0];
	}
	double y() const {
		return e[1];
	}
};
struct vec3
{					// Usage: time ./smallpt 5000 && xv image.ppm
	//double e[0], e[1], e[2]; // position, also color (r,g,b)
	double e[3];
	vec3(double x_ = 0, double y_ = 0, double z_ = 0)
	{
		e[0] = x_;
		e[1] = y_;
		e[2] = z_;
	}
	double x() const { 
		return e[0]; 
	}
	double y() const {
		return e[1];
	}
	double z() const {
		return e[2];
	}
};

vec3 operator+(const vec3& lhs, const vec3& rhs) {
	return vec3(lhs.e[0] + rhs.e[0], lhs.e[1] + rhs.e[1], lhs.e[2] + rhs.e[2]);
}
vec3 operator-(const vec3& lhs, const vec3& rhs) {
	return vec3(lhs.e[0] - rhs.e[0], lhs.e[1] - rhs.e[1], lhs.e[2] - rhs.e[2]);
}
vec3 operator*(const vec3& lhs, double rhs) {
	return vec3(rhs * lhs.e[0], rhs * lhs.e[1], rhs * lhs.e[2]);
}
vec3 operator*(double rhs, const vec3& lhs) {
	return vec3(rhs * lhs.e[0], rhs * lhs.e[1], rhs * lhs.e[2]);
}
double operator*(const vec3& lhs, const vec3& rhs)
{
	return lhs.e[0] * rhs.e[0] + lhs.e[1] * rhs.e[1] + lhs.e[2] * rhs.e[2];
}
vec3 operator/(const vec3& lhs, double rhs) {
	return vec3(lhs.e[0] / rhs, lhs.e[1] / rhs, lhs.e[2] / rhs);
}
vec3 operator/(double rhs, const vec3& lhs) {
	return vec3(lhs.e[0] / rhs, lhs.e[1] / rhs, lhs.e[2] / rhs);
}
vec3 normalize(vec3 lhs)
{
	return lhs * (1.0 / sqrt(lhs.e[0] * lhs.e[0] + lhs.e[1] * lhs.e[1] + lhs.e[2] * lhs.e[2]));
}
vec3 cross(const vec3& lhs, const vec3& rhs)
{
	return vec3(lhs.e[1] * rhs.e[2] - lhs.e[2] * rhs.e[1], lhs.e[2] * rhs.e[0] - lhs.e[0] * rhs.e[2], lhs.e[0] * rhs.e[1] - lhs.e[1] * rhs.e[0]);
}
vec3 multiply(const vec3& lhs, const vec3& rhs) {
	return vec3(lhs.e[0] * rhs.e[0], lhs.e[1] * rhs.e[1], lhs.e[2] * rhs.e[2]);
}
inline double random()
{
	return dis(gen);
}
// generate a sample point in unit sphere
inline vec2 random2()
{
	vec2 ret;
	do
	{
		ret.e[0] = dis2(gen);
		ret.e[1] = dis2(gen);
	} //while (ret.e[0] * ret.e[0] + ret.e[1] * ret.e[1] > 1.0);
	while (0);
	return ret;
}
inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
// gamma correction
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }

inline double deg2rad(double deg) { return deg*PI / 180.0; }