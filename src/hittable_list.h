/*
Class for a list of hittable objects. For instance used to store all objects in a scene
*/
#ifndef HITTABLE_LIST_H
#define HITTABLE_LIST_H

#include <memory>
#include <vector>

#include "hittable.h"
#include "aabb.h"

using std::shared_ptr;
using std::make_shared;

// List of hittable objects, with a function that checks if any of these objects is hit by a ray
class hittable_list : public hittable {
public:
	hittable_list() {}
	hittable_list(shared_ptr<hittable> object) { add(object); }

	// Interaction functions
	void clear() { objects.clear(); }
	void add(shared_ptr<hittable> object) { objects.push_back(object); }
	size_t size() { return objects.size(); }

	// Determine if any of the objects in the list are hit
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
		hit_record temp_rec;
		bool hit_anything = false;
		double closest = t_max;

		for (const auto& object : objects) {
			if (object->hit(r, t_min, closest, temp_rec)) {
				hit_anything = true;
				closest = temp_rec.t;
				rec = temp_rec;
			}
		}
		return hit_anything;
	}

	// Create a bounding box that surrounds all of the objects in the list
	virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
		if (objects.empty()) return false;

		aabb temp_box;
		bool first_box = true;
		for (const auto& object : objects) {
			if (!object->bounding_box(time0, time1, temp_box)) return false;
			output_box = first_box ? temp_box : surrounding_box(output_box, temp_box);
			first_box = false;
		}
		return true;
	}

public:
	std::vector<shared_ptr<hittable>> objects;
};

#endif // !HITTABLE_LIST_H