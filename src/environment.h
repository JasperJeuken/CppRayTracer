/*
Class that handles an environment texture (used when ray hits infinity / misses all objects)
*/
#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include "texture.h"
#include "utility.h"

// Class that stores a texture used as the environment of the scene
class environment {
public:
	environment() {}
	environment(color c) : t(make_shared<solid_color>(c)) {}
	environment(shared_ptr<texture> _t) : t(_t) {}

	// Color value based on where the ray would intersect an infinitely large sphere
	color value(const ray& r) const {
		vec3 p = unit_vector(r.direction());
		auto theta = acos(-p.y());
		auto phi = atan2(-p.z(), p.x()) + pi;
		double u = phi / (2 * pi);
		double v = theta / pi;
		return t->value(u, v, p);
	}

public:
	shared_ptr<texture> t;
};

#endif // !ENVIRONMENT_H
