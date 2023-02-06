/*
3D vector class with various implemented definitions and calculations
*/
#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>

#include "utility.h"

// 3D vector class
class vec3 {
public:
	vec3() : e{ 0, 0, 0 } {}
	vec3(double e0, double e1, double e2) : e{ e0, e1, e2 } {}

	// Get functions
	double x() const { return e[0]; }
	double y() const { return e[1]; }
	double z() const { return e[2]; }
	double operator[](int i) const { return e[i]; }
	double& operator[](int i) { return e[i]; }

	// Negation operator
	vec3 operator-() const {
		return vec3(-e[0], -e[1], -e[2]);
	}

	// Addition assignment operator
	vec3& operator+=(const vec3& other) {
		e[0] += other.e[0];
		e[1] += other.e[1];
		e[2] += other.e[2];
		return *this;
	}

	// Multiplication assignment operator (with scalar)
	vec3& operator*=(const double t) {
		e[0] *= t;
		e[1] *= t;
		e[2] *= t;
		return *this;
	}

	// Multiplication assignment operator (with vec3)
	vec3& operator*=(const vec3& other) {
		e[0] *= other[0];
		e[1] *= other[1];
		e[2] *= other[2];
		return *this;
	}

	// Division assignment operator (with scalar)
	vec3& operator/=(const double t) {
		return *this *= 1 / t;
	}

	// Square of the length of the vector
	double length_squared() const {
		return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
	}

	// Length of the vector
	double length() const {
		return std::sqrt(length_squared());
	}

	// Whether all elements of the vector are near-zero (1e-8)
	bool near_zero() const {
		const double s = 1e-8;
		return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
	}

	// Generate a random vector with components between 0.0 and 1.0
	inline static vec3 random() {
		return vec3(random_double(), random_double(), random_double());
	}

	// Generate a random vector with components between a minimum and a maximum value
	inline static vec3 random(double min, double max) {
		return vec3(random_double(min, max), random_double(min, max), random_double(min, max));
	}

public:
	double e[3];
};

// Aliases
using point3 = vec3;
using color = vec3;

// Insertion operator (for printing vector)
inline std::ostream& operator<<(std::ostream& out, const vec3& v) {
	return out << v.e[0] << " " << v.e[1] << " " << v.e[2];
}

// Addition operator (with vec3)
inline vec3 operator+(const vec3& u, const vec3& v) {
	return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

// Addition operator (with scalar)
inline vec3 operator+(const vec3& u, double v) {
	return vec3(u.e[0] + v, u.e[1] + v, u.e[2] + v);
}

// Subtraction operator (with vec3)
inline vec3 operator-(const vec3& u, const vec3& v) {
	return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

// Subtraction operator (with scalar)
inline vec3 operator-(const vec3& u, double v) {
	return vec3(u.e[0] - v, u.e[1] - v, u.e[2] - v);
}

// Multiplication operator (with vec3, element-wise)
inline vec3 operator*(const vec3& u, const vec3& v) {
	return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

// Multiplication operator (with scalar, left-hand)
inline vec3 operator*(double t, const vec3& v) {
	return vec3(t * v.e[0], t * v.e[1], t * v.e[2]);
}

// Multiplication operator (with scalar, right-hand)
inline vec3 operator*(const vec3& v, double t) {
	return t * v;
}

// Division operator (with scalar)
inline vec3 operator/(const vec3& v, double t) {
	return (1 / t) * v;
}

// Dot product of two vectors
inline double dot(const vec3& u, const vec3& v) {
	return u.e[0] * v.e[0] + u.e[1] * v.e[1] + u.e[2] * v.e[2];
}

// Cross product of two vectors
inline vec3 cross(const vec3& u, const vec3& v) {
	return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
		u.e[2] * v.e[0] - u.e[0] * v.e[2],
		u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

// Element-wise square root
inline vec3 sqrt(const vec3& v) {
	return vec3(sqrt(v.e[0]), sqrt(v.e[1]), sqrt(v.e[0]));
}

// Vector scaled to unit length
inline vec3 unit_vector(const vec3& v) {
	return v / v.length();
}

// Random vector in a unit-sized sphere
inline vec3 random_in_unit_sphere() {
	while (true) {
		vec3 p = vec3::random(-1, 1);
		if (p.length_squared() >= 1) continue;
		return p;
	}
}

// Random unit vector in a unit-sized sphere
inline vec3 random_unit_vector() {
	return unit_vector(random_in_unit_sphere());
}

// Random vector in a unit-sized hemisphere (based on normal direction)
inline vec3 random_in_hemisphere(const vec3& normal) {
	vec3 in_unit_sphere = random_in_unit_sphere();
	if (dot(in_unit_sphere, normal) > 0.0)
		return in_unit_sphere;
	else
		return -in_unit_sphere;
}

// Random vector in a unit-sized disk
inline vec3 random_in_unit_disk() {
	while (true) {
		vec3 p = vec3(random_double(-1, 1), random_double(-1, 1), 0.0);
		if (p.length_squared() >= 1) continue;
		return p;
	}
}

// Vector reflection with based on a normal direction
inline vec3 reflect(const vec3& v, const vec3& n) {
	return v - 2 * dot(v, n) * n;
}

// Vector refraction with a normal and material properties
inline vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat) {
	double cos_theta = fmin(dot(-uv, n), 1.0);
	vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
	vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.length_squared())) * n;
	return r_out_perp + r_out_parallel;
}

// Function that remaps a color value to range int[0, 255]
int convert_color_value(double value) {
	return static_cast<int>(256 * clamp(value, 0.0, 0.999));
}

// Function that writets the red, green, and blue value of a color to a stream (after applying gamma correction)
void write_color(std::ostream& out, const vec3& pixel_color) {
	out << convert_color_value(pixel_color.x()) << " "
		<< convert_color_value(pixel_color.y()) << " "
		<< convert_color_value(pixel_color.z()) << "\n";
}

#endif // !VEC3_H
