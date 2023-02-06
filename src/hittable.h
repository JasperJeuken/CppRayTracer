/*
Abstract classes for hittable objects
*/
#ifndef HITTABLE_H
#define HITTABLE_H

#include "ray.h"
#include "utility.h"
#include "aabb.h"
#include "texture.h"


class material;

// Struct that stores hit information
struct hit_record {
	point3 p;
	vec3 normal;
	double t;
	double u;
	double v;
	bool front_face;
	shared_ptr<material> mat_ptr;

    // Function that determines the hit's normal direction
	inline void set_face_normal(const ray& r, const vec3& outward_normal) {
		front_face = dot(r.direction(), outward_normal) < 0;
		normal = front_face ? outward_normal : -outward_normal;
	}
};

// Base class for hittable objects. Defines virtual functions to be overridden by derived classes
class hittable {
public:
    // Hit function: returns whether object is hit, and stores hit information in the hit_record struct
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;

    // Bounding box function: returns whether the objects has a bounding box, and sets output_box to this box if so
	virtual bool bounding_box(double time0, double time1, aabb& output_box) const = 0;

    // Function that loads an image_texture from a file, to be used as the object's normal map
    virtual void load_normal_map(const char* filename, double strength = 1.0) {
        auto map = make_shared<image_texture>(filename);
        set_normal_map(map, strength);
    }

    // Function that sets an image_texture as the object's normal map
    virtual void set_normal_map(shared_ptr<image_texture> map, double strength = 1.0) {
        normal_map = map;
        has_normal_map = true;
        normal_map_strength = strength;
    }

   // Function that calculate a new normal based on the loaded normal map using TBN matrix
   vec3 apply_normal_map(vec3& normal, vec3& tangent, double u, double v, point3& p) const {
       // Get normal map value
        color map_value = 2.0 * normal_map->value(u, v, p) - 1.0;
        map_value *= vec3(normal_map_strength, normal_map_strength, 1.0);

        // Calculate bitangent
        vec3 bitangent = cross(normal, tangent);

        // Calculate new normal using TBN matrix calculations
        double x0 = map_value[0];
        double x1 = map_value[1];
        double x2 = map_value[2];
        double row1 = tangent[0] * x0 + bitangent[0] * x1 + normal[0] * x2;
        double row2 = tangent[1] * x0 + bitangent[1] * x1 + normal[1] * x2;
        double row3 = tangent[2] * x0 + bitangent[2] * x1 + normal[2] * x2;
        return unit_vector(vec3(row1, row2, row3));
    }

public:
    shared_ptr<image_texture> normal_map = nullptr;
    bool has_normal_map = false;
    double normal_map_strength;
};

// Class that translates a hittable object using a displacement vector
class translate : public hittable {
public:
	translate(shared_ptr<hittable> p, const vec3& displacement) : ptr(p), offset(displacement) {}

    // Translate the incoming ray and perform the hit calculation
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
		ray moved_r(r.origin() - offset, r.direction(), r.time());
		if (!ptr->hit(moved_r, t_min, t_max, rec))
			return false;

		rec.p += offset;
		rec.set_face_normal(moved_r, rec.normal);

		return true;
	}

    // Translate the bounding box of the object by the offset
	virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
		if (!ptr->bounding_box(time0, time1, output_box))
			return false;

		output_box = aabb(
			output_box.min() + offset,
			output_box.max() + offset);

		return true;
	}

public:
	shared_ptr<hittable> ptr;
	vec3 offset;
};

// Class that rotates a hittable object around the y-axis by some angle
class rotate_y : public hittable {
public:
    rotate_y(shared_ptr<hittable> p, double angle);

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
        output_box = bbox;
        return hasbox;
    }

public:
    shared_ptr<hittable> ptr;
    double sin_theta;
    double cos_theta;
    bool hasbox;
    aabb bbox;
};

// Calculate the angles used for the rotation, and create a new rotated bounding box
rotate_y::rotate_y(shared_ptr<hittable> p, double angle) : ptr(p) {
    // Calculate rotation angles
    double radians = degrees_to_radians(angle);
    sin_theta = sin(radians);
    cos_theta = cos(radians);
    hasbox = ptr->bounding_box(0, 1, bbox);

    // Create rotated bounding box
    point3 min(infinity, infinity, infinity);
    point3 max(-infinity, -infinity, -infinity);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                double x = i * bbox.max().x() + (1.0 - i) * bbox.min().x();
                double y = j * bbox.max().y() + (1.0 - j) * bbox.min().y();
                double z = k * bbox.max().z() + (1.0 - k) * bbox.min().z();

                double newx = cos_theta * x + sin_theta * z;
                double newz = -sin_theta * x + cos_theta * z;

                vec3 tester(newx, y, newz);

                for (int c = 0; c < 3; c++) {
                    min[c] = fmin(min[c], tester[c]);
                    max[c] = fmax(max[c], tester[c]);
                }
            }
        }
    }
    bbox = aabb(min, max);
}

// Rotate the incoming ray and perform the hit calculation
bool rotate_y::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    point3 origin = r.origin();
    vec3 direction = r.direction();

    origin[0] = cos_theta * r.origin()[0] - sin_theta * r.origin()[2];
    origin[2] = sin_theta * r.origin()[0] + cos_theta * r.origin()[2];

    direction[0] = cos_theta * r.direction()[0] - sin_theta * r.direction()[2];
    direction[2] = sin_theta * r.direction()[0] + cos_theta * r.direction()[2];

    ray rotated_r(origin, direction, r.time());

    if (!ptr->hit(rotated_r, t_min, t_max, rec))
        return false;

    point3 p = rec.p;
    vec3 normal = rec.normal;

    p[0] = cos_theta * rec.p[0] + sin_theta * rec.p[2];
    p[2] = -sin_theta * rec.p[0] + cos_theta * rec.p[2];

    normal[0] = cos_theta * rec.normal[0] + sin_theta * rec.normal[2];
    normal[2] = -sin_theta * rec.normal[0] + cos_theta * rec.normal[2];

    rec.p = p;
    rec.set_face_normal(rotated_r, normal);

    return true;
}

// Class that rotates a hittable object around the x-axis by some angle
class rotate_x : public hittable {
public:
    rotate_x(shared_ptr<hittable> p, double angle);

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
        output_box = bbox;
        return hasbox;
    }

public:
    shared_ptr<hittable> ptr;
    double sin_theta;
    double cos_theta;
    bool hasbox;
    aabb bbox;
};

// Calculate the angles used for the rotation, and create a new rotated bounding box
rotate_x::rotate_x(shared_ptr<hittable> p, double angle) : ptr(p) {
    // Calculate rotation angles
    double radians = degrees_to_radians(angle);
    sin_theta = sin(radians);
    cos_theta = cos(radians);
    hasbox = ptr->bounding_box(0, 1, bbox);

    // Create rotated bounding box
    point3 min(infinity, infinity, infinity);
    point3 max(-infinity, -infinity, -infinity);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                double x = i * bbox.max().x() + (1.0 - i) * bbox.min().x();
                double y = j * bbox.max().y() + (1.0 - j) * bbox.min().y();
                double z = k * bbox.max().z() + (1.0 - k) * bbox.min().z();

                double newy = cos_theta * y + sin_theta * z;
                double newz = -sin_theta * y + cos_theta * z;

                vec3 tester(x, newy, newz);

                for (int c = 0; c < 3; c++) {
                    min[c] = fmin(min[c], tester[c]);
                    max[c] = fmax(max[c], tester[c]);
                }
            }
        }
    }
    bbox = aabb(min, max);
}

// Rotate the incoming ray and perform the hit calculation
bool rotate_x::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    point3 origin = r.origin();
    vec3 direction = r.direction();

    origin[1] = cos_theta * r.origin()[1] - sin_theta * r.origin()[2];
    origin[2] = sin_theta * r.origin()[1] + cos_theta * r.origin()[2];

    direction[1] = cos_theta * r.direction()[1] - sin_theta * r.direction()[2];
    direction[2] = sin_theta * r.direction()[1] + cos_theta * r.direction()[2];

    ray rotated_r(origin, direction, r.time());

    if (!ptr->hit(rotated_r, t_min, t_max, rec))
        return false;

    point3 p = rec.p;
    vec3 normal = rec.normal;

    p[1] = cos_theta * rec.p[1] + sin_theta * rec.p[2];
    p[2] = -sin_theta * rec.p[1] + cos_theta * rec.p[2];

    normal[1] = cos_theta * rec.normal[1] + sin_theta * rec.normal[2];
    normal[2] = -sin_theta * rec.normal[1] + cos_theta * rec.normal[2];

    rec.p = p;
    rec.set_face_normal(rotated_r, normal);

    return true;
}

// Class that rotates a hittable object around the z-axis by some angle
class rotate_z : public hittable {
public:
    rotate_z(shared_ptr<hittable> p, double angle);

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
        output_box = bbox;
        return hasbox;
    }

public:
    shared_ptr<hittable> ptr;
    double sin_theta;
    double cos_theta;
    bool hasbox;
    aabb bbox;
};

// Calculate the angles used for the rotation, and create a new rotated bounding box
rotate_z::rotate_z(shared_ptr<hittable> p, double angle) : ptr(p) {
    // Calculate rotation angles
    double radians = degrees_to_radians(angle);
    sin_theta = sin(radians);
    cos_theta = cos(radians);
    hasbox = ptr->bounding_box(0, 1, bbox);

    // Create rotated bounding box
    point3 min(infinity, infinity, infinity);
    point3 max(-infinity, -infinity, -infinity);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                double x = i * bbox.max().x() + (1.0 - i) * bbox.min().x();
                double y = j * bbox.max().y() + (1.0 - j) * bbox.min().y();
                double z = k * bbox.max().z() + (1.0 - k) * bbox.min().z();

                double newx = cos_theta * x + sin_theta * y;
                double newy = -sin_theta * x + cos_theta * y;

                vec3 tester(newx, newy, z);

                for (int c = 0; c < 3; c++) {
                    min[c] = fmin(min[c], tester[c]);
                    max[c] = fmax(max[c], tester[c]);
                }
            }
        }
    }
    bbox = aabb(min, max);
}

// Rotate the incoming ray and perform the hit calculation
bool rotate_z::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    point3 origin = r.origin();
    vec3 direction = r.direction();

    origin[0] = cos_theta * r.origin()[0] - sin_theta * r.origin()[1];
    origin[1] = sin_theta * r.origin()[0] + cos_theta * r.origin()[1];

    direction[0] = cos_theta * r.direction()[0] - sin_theta * r.direction()[1];
    direction[1] = sin_theta * r.direction()[0] + cos_theta * r.direction()[1];

    ray rotated_r(origin, direction, r.time());

    if (!ptr->hit(rotated_r, t_min, t_max, rec))
        return false;

    point3 p = rec.p;
    vec3 normal = rec.normal;

    p[0] = cos_theta * rec.p[0] + sin_theta * rec.p[1];
    p[1] = -sin_theta * rec.p[0] + cos_theta * rec.p[1];

    normal[0] = cos_theta * rec.normal[0] + sin_theta * rec.normal[1];
    normal[1] = -sin_theta * rec.normal[0] + cos_theta * rec.normal[1];

    rec.p = p;
    rec.set_face_normal(rotated_r, normal);

    return true;
}

// Class that rotates a hittable about the x-, y-, and z-axes, in that order
class rotate : public hittable {
public:
    rotate(shared_ptr<hittable> p, vec3 rot) : rotations(rot) {
        auto rotated_x = make_shared<rotate_x>(p, rotations[0]);
        auto rotated_y = make_shared<rotate_y>(rotated_x, rotations[1]);
        ptr = make_shared<rotate_z>(rotated_y, rotations[2]);
    }

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
        return ptr->hit(r, t_min, t_max, rec);
    }

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
        return ptr->bounding_box(time0, time1, output_box);
    }

public:
    shared_ptr<hittable> ptr;
    vec3 rotations;
};

#endif // !HITTABLE_H
