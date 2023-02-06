/*
Axis-aligned bounding box. Used in BVH tree for optimized hit calculations.
*/
#ifndef AABB_H
#define AABB_H

#include "utility.h"

// Axis-aligned bounding box, defined by two points (opposite corners of the box)
// Used for optimized hit calculations
class aabb {
public:
	aabb() {}
	aabb(const point3& a, const point3& b) { minimum = a; maximum = b; }

	// Get functions
	point3 min() const { return minimum; }
	point3 max() const { return maximum; }

	// Check if box is hit by a ray (only returns true/false)
	bool hit(const ray& r, double t_min, double t_max) const {
		for (int a = 0; a < 3; a++) {
			double inv_d = 1.0 / r.direction()[a];
			double t0 = (minimum[a] - r.origin()[a]) * inv_d;
			double t1 = (maximum[a] - r.origin()[a]) * inv_d;
			if (inv_d < 0.0)
				std::swap(t0, t1);
			t_min = t0 > t_min ? t0 : t_min;
			t_max = t1 < t_max ? t1 : t_max;
			if (t_max <= t_min)
				return false;
		}
		return true;
	}

public:
	point3 minimum;
	point3 maximum;
};

// Function that creates an aabb that surrounds two other aabbs
aabb surrounding_box(aabb box0, aabb box1) {
	point3 small(fmin(box0.min().x(), box1.min().x()), fmin(box0.min().y(), box1.min().y()), fmin(box0.min().z(), box1.min().z()));
	point3 big(fmax(box0.max().x(), box1.max().x()), fmax(box0.max().y(), box1.max().y()), fmax(box0.max().z(), box1.max().z()));
	return aabb(small, big);
}

#endif // !AABB_H
