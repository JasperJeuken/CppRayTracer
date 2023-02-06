/*
Collection of hittable objects that can be placed into a scene (such as spheres, rectangles, boxes, triangles, ...)
*/
#ifndef OBJECTS_H
#define OBJECTS_H

#include <cmath>

#include "hittable.h"
#include "hittable_list.h"
#include "vec3.h"
#include "material.h"
#include "bvh.h"

#include "external/OBJ_Loader.h"
#include "external/stl_reader.h"

// Sphere defined by a center point, a radius, and a certain material
class sphere : public hittable {
public:
	sphere() {}
	sphere(point3 cen, double r, shared_ptr<material> m) : center(cen), radius(r), mat_ptr(m) {};

	// Determine if and where the sphere is hit by a ray
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
		// Calculate hit properties
		vec3 oc = r.origin() - center;
		double a = r.direction().length_squared();
		double half_b = dot(oc, r.direction());
		double c = oc.length_squared() - radius * radius;

		double discriminant = half_b * half_b - a * c;
		if (discriminant < 0) return false;

		double sqrt_d = std::sqrt(discriminant);
		double root = (-half_b - sqrt_d) / a;
		if (root < t_min || t_max < root) {
			root = (-half_b + sqrt_d) / a;
			if (root < t_min || t_max < root)
				return false;
		}

		// Store hit information in hit_record
		rec.t = root;
		rec.p = r.at(rec.t);
		vec3 outward_normal = (rec.p - center) / radius;
		get_sphere_uv(outward_normal, rec.u, rec.v);

		// Recalculate normal if the sphere has a normal map
		if (!has_normal_map) {
			rec.set_face_normal(r, outward_normal);
		}
		else {
			vec3 tangent = unit_vector(cross(vec3(0, -1, 0), rec.p - center));
			vec3 new_normal = apply_normal_map(outward_normal, tangent, rec.u, rec.v, rec.p);
			rec.set_face_normal(r, new_normal);
		}

		rec.mat_ptr = mat_ptr;
		return true;
	};

	// Create a bounding box using the sphere's center and radius
	virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
		output_box = aabb(center - vec3(radius, radius, radius), center + vec3(radius, radius, radius));
		return true;
	}

private:
	// Get the UV-coordinates of a point on the sphere
	static void get_sphere_uv(const point3& p, double& u, double& v) {
		double theta = acos(-p.y());
		double phi = atan2(-p.z(), p.x()) + pi;
		u = phi / (2 * pi);
		v = theta / pi;
	}

public:
	point3 center;
	double radius;
	shared_ptr<material> mat_ptr;
};

// Linearly moving sphere defined by two center points (one at time 0, one at time 1), a radius, and a certain material
class moving_sphere : public hittable {
public:
	moving_sphere() {}
	moving_sphere(point3 cen0, point3 cen1, double _time0, double _time1, double r, shared_ptr<material> m) : center0(cen0), center1(cen1), time0(_time0), time1(_time1), radius(r), mat_ptr(m) {}

	// Determine if and where the moving sphere is hit by a ray
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
		// Calculate hit properties
		vec3 oc = r.origin() - center(r.time());
		double a = r.direction().length_squared();
		double half_b = dot(oc, r.direction());
		double c = oc.length_squared() - radius * radius;

		double discriminant = half_b * half_b - a * c;
		if (discriminant < 0) return false;

		double sqrt_d = std::sqrt(discriminant);
		double root = (-half_b - sqrt_d) / a;
		if (root < t_min || t_max < root) {
			root = (-half_b + sqrt_d) / a;
			if (root < t_min || t_max < root)
				return false;
		}
		
		// Store hit information in hit_record
		rec.t = root;
		rec.p = r.at(rec.t);
		vec3 outward_normal = (rec.p - center(r.time())) / radius;
		get_sphere_uv(outward_normal, rec.u, rec.v);

		// Recalculate normal if the moving sphere has a normal map
		if (!has_normal_map) {
			rec.set_face_normal(r, outward_normal);
		}
		else {
			vec3 tangent = unit_vector(cross(vec3(0, 1, 0), rec.p - center(r.time())));
			vec3 new_normal = apply_normal_map(outward_normal, tangent, rec.u, rec.v, rec.p);
			rec.set_face_normal(r, new_normal);
		}

		rec.mat_ptr = mat_ptr;
		return true;
	}

	// Create a bounding box using the sphere's centers and radius
	virtual bool bounding_box(double _time0, double _time1, aabb& output_box) const override {
		aabb box0(center(_time0) - vec3(radius, radius, radius), center(_time0) + vec3(radius, radius, radius));
		aabb box1(center(_time1) - vec3(radius, radius, radius), center(_time1) + vec3(radius, radius, radius));
		output_box = surrounding_box(box0, box1);
		return true;
	}

	// Calculate the center of the moving sphere at a given time
	point3 center(double time) const {
		return center0 + ((time - time0) / (time1 - time0)) * (center1 - center0);
	}

public:
	point3 center0, center1;
	double time0, time1;
	double radius;
	shared_ptr<material> mat_ptr;

private:
	// Get the UV-coordinates of a point on the moving sphere
	static void get_sphere_uv(const point3& p, double& u, double& v) {
		double theta = acos(-p.y());
		double phi = atan2(-p.z(), p.x()) + pi;
		u = phi / (2 * pi);
		v = theta / pi;
	}
};

// Class that creates a rectangle parallel to the xy-plane
class xy_rect : public hittable {
public:
	xy_rect() {}

	xy_rect(double _x0, double _x1, double _y0, double _y1, double _k, shared_ptr<material> mat) : x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k), mp(mat) {};

	// Determine if and where the rectangle is hit by a ray
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
		// Calculate hit properties
		double t = (k - r.origin().z()) / r.direction().z();
		if (t < t_min || t > t_max)
			return false;

		double x = r.origin().x() + t * r.direction().x();
		double y = r.origin().y() + t * r.direction().y();
		if (x < x0 || x > x1 || y < y0 || y > y1)
			return false;

		// Store hit information in hit_record
		rec.u = (x - x0) / (x1 - x0);
		rec.v = (y - y0) / (y1 - y0);
		rec.t = t;
		vec3 outward_normal(0, 0, 1);

		// Recalculate normal if the rectangle has a normal map
		if (!has_normal_map) {
			rec.set_face_normal(r, outward_normal);
		}
		else {
			vec3 tangent(1, 0, 0);
			vec3 new_normal = apply_normal_map(outward_normal, tangent, rec.u, rec.v, rec.p);
			rec.set_face_normal(r, new_normal);
		}

		rec.mat_ptr = mp;
		rec.p = r.at(t);
		return true;
	}

	// Create a bounding box using the rectangle's dimensions (give slight thickness)
	virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
		output_box = aabb(point3(x0, y0, k - 0.0001), point3(x1, y1, k + 0.0001));
		return true;
	}

public:
	shared_ptr<material> mp;
	double x0, x1, y0, y1, k;
};

// Class that creates a rectangle parallel to the xz-plane
class xz_rect : public hittable {
public:
	xz_rect() {}

	xz_rect(double _x0, double _x1, double _z0, double _z1, double _k, shared_ptr<material> mat) : x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), mp(mat) {};

	// Determine if and where the rectangle is hit by a ray
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
		// Calculate hit properties
		double t = (k - r.origin().y()) / r.direction().y();
		if (t < t_min || t > t_max)
			return false;

		double x = r.origin().x() + t * r.direction().x();
		double z = r.origin().z() + t * r.direction().z();
		if (x < x0 || x > x1 || z < z0 || z > z1)
			return false;
		
		// Store hit information in hit_record
		rec.u = (x - x0) / (x1 - x0);
		rec.v = (z - z0) / (z1 - z0);
		rec.t = t;
		vec3 outward_normal(0, 1, 0);
		
		// Recalculate normal if the rectangle has a normal map
		if (!has_normal_map) {
			rec.set_face_normal(r, outward_normal);
		}
		else {
			vec3 tangent(1, 0, 0);
			vec3 new_normal = apply_normal_map(outward_normal, tangent, rec.u, rec.v, rec.p);
			rec.set_face_normal(r, new_normal);
		}

		rec.mat_ptr = mp;
		rec.p = r.at(t);
		return true;
	}

	// Create a bounding box using the rectangle's dimensions (give slight thickness)
	virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
		output_box = aabb(point3(x0, k - 0.0001, z0), point3(x1, k + 0.0001, z1));
		return true;
	}

public:
	shared_ptr<material> mp;
	double x0, x1, z0, z1, k;
};

// Class that creates a rectangle parallel to the yz-plane
class yz_rect : public hittable {
public:
	yz_rect() {}

	yz_rect(double _y0, double _y1, double _z0, double _z1, double _k, shared_ptr<material> mat) : y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), mp(mat) {};

	// Determine if and where the rectangle is hit by a ray
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
		// Calculate hit properties
		double t = (k - r.origin().x()) / r.direction().x();
		if (t < t_min || t > t_max)
			return false;

		double y = r.origin().y() + t * r.direction().y();
		double z = r.origin().z() + t * r.direction().z();
		if (y < y0 || y > y1 || z < z0 || z > z1)
			return false;

		// Store hit information in hit_record
		rec.u = (y - y0) / (y1 - y0);
		rec.v = (z - z0) / (z1 - z0);
		rec.t = t;
		vec3 outward_normal(1, 0, 0);

		// Recalculate normal if the rectangle has a normal map
		if (!has_normal_map) {
			rec.set_face_normal(r, outward_normal);
		}
		else {
			vec3 tangent(0, 1, 0);
			vec3 new_normal = apply_normal_map(outward_normal, tangent, rec.u, rec.v, rec.p);
			rec.set_face_normal(r, new_normal);
		}
		rec.mat_ptr = mp;
		rec.p = r.at(t);
		return true;
	}

	// Create a bounding box using the rectangle's dimensions (give slight thickness)
	virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
		output_box = aabb(point3(k - 0.0001, y0, z0), point3(k + 0.0001, y0, z1));
		return true;
	}

public:
	shared_ptr<material> mp;
	double y0, y1, z0, z1, k;
};

// Class that creates a box with 6 sides (aligned along the xyz-axes), defined by two corner points
class box : public hittable {
public:
	box() {}

	// Create six sides of the box between corner points p0 and p1
	box(const point3& p0, const point3& p1, shared_ptr<material> ptr) {
		box_min = p0;
		box_max = p1;

		sides.add(make_shared<xy_rect>(p0.x(), p1.x(), p0.y(), p1.y(), p1.z(), ptr));
		sides.add(make_shared<xy_rect>(p0.x(), p1.x(), p0.y(), p1.y(), p0.z(), ptr));

		sides.add(make_shared<xz_rect>(p0.x(), p1.x(), p0.z(), p1.z(), p1.y(), ptr));
		sides.add(make_shared<xz_rect>(p0.x(), p1.x(), p0.z(), p1.z(), p0.y(), ptr));

		sides.add(make_shared<yz_rect>(p0.y(), p1.y(), p0.z(), p1.z(), p1.x(), ptr));
		sides.add(make_shared<yz_rect>(p0.y(), p1.y(), p0.z(), p1.z(), p0.x(), ptr));
	}

	// Determine if any of the sides is hit
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
		return sides.hit(r, t_min, t_max, rec);
	}
	
	// Use the corner points to create a bounding box
	virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
		output_box = aabb(box_min, box_max);
		return true;
	}

public:
	point3 box_min;
	point3 box_max;
	hittable_list sides;
};

// Class for a constant-density medium (fog, smoke, ...), only works for convex shapes
class constant_medium : public hittable {
public:
	constant_medium(shared_ptr<hittable> b, double d, shared_ptr<texture> a) : boundary(b), neg_inv_density(-1 / d), phase_function(make_shared<isotropic>(a)) {}

	constant_medium(shared_ptr<hittable> b, double d, color c) : boundary(b), neg_inv_density(-1 / d), phase_function(make_shared<isotropic>(c)) {}

	// Hit function works for rays originating outside and inside the boundary of the object
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
		hit_record rec1, rec2;

		// Calculate hit information
		if (!boundary->hit(r, -infinity, infinity, rec1))
			return false;

		if (!boundary->hit(r, rec1.t + 0.0001, infinity, rec2))
			return false;

		if (rec1.t < t_min) rec1.t = t_min;
		if (rec2.t > t_max) rec2.t = t_max;

		if (rec1.t >= rec2.t)
			return false;

		if (rec1.t < 0)
			rec1.t = 0;

		const double ray_length = r.direction().length();
		const double distance_inside_boundary = (rec2.t - rec1.t) * ray_length;
		const double hit_distance = neg_inv_density * log(random_double());

		if (hit_distance > distance_inside_boundary)
			return false;

		rec.t = rec1.t + hit_distance / ray_length;
		rec.p = r.at(rec.t);
		rec.normal = vec3(1, 0, 0); // arbitrary
		rec.front_face = true; // arbitrary
		rec.mat_ptr = phase_function;
		return true;
	}

	// Return bounding box of boundary
	virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
		return boundary->bounding_box(time0, time1, output_box);
	}

public:
	shared_ptr<hittable> boundary;
	shared_ptr<material> phase_function;
	double neg_inv_density;
};

// Class for a disk (center point, normal, and radius)
class disk : public hittable {
public:
	disk() {}
	disk(point3 _center, vec3 _normal, double _radius, shared_ptr<material> m) : center(_center), normal(_normal), radius(_radius), mat_ptr(m) {
		right = vec3(-normal.y(), normal.x(), normal.z()); // random perpendicular vector
	}
	disk(point3 _center, vec3 _normal, vec3 _right, double _radius, shared_ptr<material> m) : center(_center), normal(_normal), right(_right), radius(_radius), mat_ptr(m) {}

	// Determine if the disk is hit
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
		double denom = dot(normal, r.direction());

		if (fabs(denom) > 0.000001) {
			vec3 q = center - r.origin();
			double t = dot(q, normal) / denom;

			if (t > t_min && t < t_max) {
				point3 p = r.at(t);
				vec3 v = p - center;
				double d2 = dot(v, v);

				if (d2 < radius * radius) {
					rec.t = t;
					rec.p = p;
					rec.mat_ptr = mat_ptr;

					vec3 up = cross(normal, right);
					vec3 i = rec.p - center;
					rec.u = dot(i, right) * 0.5 + 0.5;
					rec.v = dot(i, up) * 0.5 + 0.5;

					// Recalculate normal if the disk has a normal map
					if (!has_normal_map) {
						rec.set_face_normal(r, normal);
					}
					else {
						vec3 tangent = right;
						vec3 outward_normal = normal;
						vec3 new_normal = apply_normal_map(outward_normal, tangent, rec.u, rec.v, rec.p);
						rec.set_face_normal(r, new_normal);
					}
					return true;
				}
			}
		}
		return false;
	}

	// Create bounding box based on the disk center, normal, and radius
	virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
		vec3 e = radius * sqrt(vec3(1.0, 1.0, 1.0) - normal * normal);
		output_box = aabb(center - e, center + e);
		return true;
	}

public:
	point3 center;
	vec3 normal;
	vec3 right;
	double radius;
	shared_ptr<material> mat_ptr;
};

// Class for a 3-sided polygon (triangle) defined by three vertices
class triangle : public hittable {
public:
	triangle() {}
	triangle(const point3& _v0, const point3& _v1, const point3& _v2, shared_ptr<material> m, bool os = false) : v0(_v0), v1(_v1), v2(_v2), mat_ptr(m), one_sided(os) {
		uv0 = vec3(0, 0, 0);
		uv1 = vec3(0, 1, 0);
		uv2 = vec3(1, 0, 0);
	}

	triangle(const point3& _v0, const point3& _v1, const point3& _v2, point3& _uv0, point3& _uv1, point3& _uv2, shared_ptr<material> m, bool os = false) : v0(_v0), v1(_v1), v2(_v2), uv0(_uv0), uv1(_uv1), uv2(_uv2), mat_ptr(m), one_sided(os) {}

	// Determine if a ray hits the triangle
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
		vec3 edge1 = v1 - v0;
		vec3 edge2 = v2 - v0;
		vec3 pvec = cross(r.direction(), edge2);
		double det = dot(edge1, pvec);

		if (one_sided && det < 0.000001)
			return false;

		if (fabs(det) < 0.000001)
			return false;

		double inv_det = 1.0 / det;
		vec3 tvec = r.origin() - v0;
		double u = dot(tvec, pvec) * inv_det;
		if (u < 0.0 || u > 1.0)
			return false;

		vec3 qvec = cross(tvec, edge1);
		double v = dot(r.direction(), qvec) * inv_det;
		if (v < 0.0 || u + v > 1.0)
			return false;

		double t = dot(edge2, qvec) * inv_det;

		if (t < t_min || t > t_max)
			return false;

		rec.t = t;
		rec.p = r.at(rec.t);
		rec.u = u;
		rec.v = v;
		rec.mat_ptr = mat_ptr;
		vec3 normal = unit_vector(cross(edge1, edge2));

		// Recalculate normal if the triangle has a normal map
		if (!has_normal_map) {
			rec.set_face_normal(r, normal);
		}
		else {
			vec3 edge1 = v1 - v0;
			vec3 edge2 = v2 - v0;
			vec3 delta_uv1 = uv1 - uv0;
			vec3 delta_uv2 = uv2 - uv0;

			double f = 1.0 / (delta_uv1.x() * delta_uv2.y() - delta_uv2.x() * delta_uv1.y());
			if (isinf(f)) {
				rec.set_face_normal(r, normal);
			}
			else {
				double tx = f * (delta_uv2.y() * edge1.x() - delta_uv1.y() * edge2.x());
				double ty = f * (delta_uv2.y() * edge1.y() - delta_uv1.y() * edge2.y());
				double tz = f * (delta_uv2.y() * edge1.z() - delta_uv1.y() * edge2.z());
				vec3 tangent(tx, ty, tz);

				vec3 test = cross(normal, tangent);
				double testu = uv0.x() + tangent.x() * rec.u + test.x() * rec.v;
				double testv = uv0.y() + tangent.y() * rec.u + test.y() * rec.v;

				vec3 new_normal = apply_normal_map(normal, tangent, testu, testv, rec.p);
				rec.set_face_normal(r, new_normal);
			}
		}
		return true;
	}

	// Create a bounding box based on the three vertices
	virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
		vec3 minv(fmin(fmin(v0.x(), v1.x()), v2.x()), fmin(fmin(v0.y(), v1.y()), v2.y()), fmin(fmin(v0.z(), v1.z()), v2.z()));
		vec3 maxv(fmax(fmax(v0.x(), v1.x()), v2.x()), fmax(fmax(v0.y(), v1.y()), v2.y()), fmax(fmax(v0.z(), v1.z()), v2.z()));

		vec3 diff = maxv - minv;
		if (diff.x() < 0.000001) maxv[0] += 0.000001;
		if (diff.y() < 0.000001) maxv[1] += 0.000001;
		if (diff.z() < 0.000001) maxv[2] += 0.000001;

		output_box = aabb(minv, maxv);
		return true;
	}

public:
	point3 v0, v1, v2;
	point3 uv0, uv1, uv2;
	shared_ptr<material> mat_ptr;
	bool one_sided;
};

// Class for a cylinder, defined by a radius, length, and start point (one end of cylinder)
class cylinder : public hittable {
public:
	cylinder() {}
	cylinder(shared_ptr<material> m) : radius(0.5), length(1.0), center(point3(0, 0, 0)), mat_ptr(m) { add_caps(); }
	cylinder(double r, double l, shared_ptr<material> m) : radius(r), length(l), mat_ptr(m), center(point3(0, 0, 0)) { add_caps(); }
	cylinder(double r, double l, point3 c, shared_ptr<material> m) : radius(r), length(l), center(c), mat_ptr(m) { add_caps(); }

	// Determine if a ray hits the cylinder
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
		point3 orig = r.origin();
		vec3 dir = r.direction();
		double a = dir.x() * dir.x() + dir.z() * dir.z();
		double b = 2 * (dir.x() * (orig.x() - center.x())) + 2.0 * dir.z() * (orig.z() - center.z());
		double c = (orig.x() - center.x()) * (orig.x() - center.x()) + (orig.z() - center.z()) * (orig.z() - center.z()) - radius * radius;

		double D = b * b - 4 * a * c;
		if (D <= 1e-9) return caps.hit(r, t_min, t_max, rec);

		double t0 = (-b - sqrt(D)) / (2 * a);
		double t1 = (-b + sqrt(D)) / (2 * a);
		double t;

		if (t0 > t1) t = t1;
		else t = t0;

		if (t < t_min || t > t_max) return caps.hit(r, t_min, t_max, rec);

		double y = r.at(t).y();
		if ((y >= center.y()) && (y <= center.y() + length)) {
			rec.t = t;
			rec.p = r.at(rec.t);
			std::pair<double, double> uv = uv_at(rec.p);
			rec.u = uv.first;
			rec.v = uv.second;
			vec3 normal(rec.p.x() - center.x(), 0, rec.p.z() - center.z());
			rec.set_face_normal(r, unit_vector(normal));
			rec.mat_ptr = mat_ptr;
			return true;
		}
		return caps.hit(r, t_min, t_max, rec);
	}

	// Create a bounding box using the center, radius, and length of the cylinder
	virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
		point3 minp(center.x() - radius, center.y(), center.z() - radius);
		point3 maxp(center.x() + radius, center.y() + length, center.z() + radius);
		output_box = aabb(minp, maxp);
		return true;
	}

	// Calculate uv-coordinate of point on the cylinder
	std::pair<double, double> uv_at(const point3& p) const {
		double theta = atan2(p.x(), p.z());
		double u = theta / (2 * pi);
		u = 1 - (u + 0.5);
		double v = std::fmod(p.y(), 1);
		return std::make_pair(u, v);
	}

public:
	double radius;
	double length;
	point3 center;
	shared_ptr<material> mat_ptr;
	hittable_list caps;

private:
	// Function that adds two disks at the ends of the cylinder
	void add_caps() {
		caps.add(make_shared<disk>(center, vec3(0, -1, 0), radius, mat_ptr));
		caps.add(make_shared<disk>(center + vec3(0, length, 0), vec3(0, 1, 0), radius, mat_ptr));
	}
};

// Class for a torus, defined by a center point, major radius, and minor radius
class torus : public hittable {
public:
	torus() {}

	torus(double maj_r, double min_r, shared_ptr<material> m) : major_radius(maj_r), minor_radius(min_r), mat_ptr(m) {
		vec3 minp(-major_radius - minor_radius, -minor_radius, -major_radius - minor_radius);
		vec3 maxp(major_radius + minor_radius, minor_radius, major_radius + minor_radius);
		bbox = aabb(minp, maxp);
	}

	// Determine if a ray hits the torus
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
		// Check if bounding box was hit
		if (!bbox.hit(r, t_min, t_max)) {
			return false;
		}

		// Define arrays for the polynomial equation and its solutions
		double coeffs[5];
		double roots[4];

		// Calculate polynomial coefficients
		double r_dir_sqr = r.direction().length_squared();
		double e = r.origin().length_squared() - major_radius * major_radius - minor_radius * minor_radius;
		double f = dot(r.origin(), r.direction());
		double four_mr_sqr = 4.0 * major_radius * major_radius;
		coeffs[0] = e * e - four_mr_sqr * (minor_radius * minor_radius - r.origin().y() * r.origin().y());
		coeffs[1] = 4.0 * f * e + 2.0 * four_mr_sqr * r.origin().y() * r.direction().y();
		coeffs[2] = 2.0 * r_dir_sqr * e + 4.0 * f * f + four_mr_sqr * r.direction().y() * r.direction().y();
		coeffs[3] = 4.0 * r_dir_sqr * f;
		coeffs[4] = r_dir_sqr * r_dir_sqr;

		// Calculate roots
		int num_real_roots = solve_quartic(coeffs, roots);
		if (num_real_roots == 0) return false;

		// Determine if torus was intersected
		bool intersected = false;
		double t = infinity;
		for (int j = 0; j < num_real_roots; j++) {
			if (roots[j] > 1e-20) {
				intersected = true;
				if (roots[j] < t) {
					t = roots[j];
				}
			}
		}
		if (!intersected) return false;

		// Set hit record data
		rec.t = t;
		rec.p = r.origin() + t * r.direction();
		std::pair<double, double> uv = uv_at(rec.p);
		rec.u = uv.first;
		rec.v = uv.second;
		rec.mat_ptr = mat_ptr;
		vec3 hit_normal = normal_at(rec.p);
		rec.set_face_normal(r, hit_normal);
		
		return true;
	}

	// Return torus bounding box
	virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
		output_box = bbox;
		return true;
	}

	// Calculate the outward normal vector of a point on the torus
	vec3 normal_at(const point3& p) const {
		double param_squared = major_radius * major_radius + minor_radius * minor_radius;
		double sum_sqr = p.x() * p.x() + p.y() * p.y() + p.z() * p.z();

		vec3 result;
		result[0] = 4.0 * p.x() * (sum_sqr - param_squared);
		result[1] = 4.0 * p.y() * (sum_sqr - param_squared + 2.0 * major_radius * major_radius);
		result[2] = 4.0 * p.z() * (sum_sqr - param_squared);

		return unit_vector(result);
	}

	// Calculate uv-coordinate of a point on the torus (handle edge cases)
	std::pair<double, double> uv_at(const point3& p) const {
		double u, v;
		if (p.x() < -1.0 || p.x() > 1.0) {
			u = 0;
		}
		else {
			u = map(asin(p.x() / major_radius), -pi / 2, pi / 2, 0, 1.0);
		}
		if (p.y() < -1.0 || p.y() > 1.0) {
			v = 0;
		} else {
			v = map(acos(p.y() / minor_radius), 0, pi, 0, 1.0);
		}
		return std::make_pair(u, v);
	}

public:
	double major_radius, minor_radius;
	aabb bbox;
	shared_ptr<material> mat_ptr;

private:

	// Utility function that checks if a number is near zero
	bool is_zero(double val, double acc = 1e-20) const {
		return val > -acc && val < acc;
	}

	// Function that solves a quartic equation, returns number of real roots
	// https://github.com/erich666/GraphicsGems/blob/master/gems/Roots3And4.c
	int solve_quartic(double c[5], double s[4]) const {
		double coeffs[4];
		double z, u, v, sub;
		double A, B, C, D;
		double sq_A, p, q, r;
		int i, num;

		// Write to form: X^4 + Ax^3 + Bx^2 + Cx + D = 0
		A = c[3] / c[4];
		B = c[2] / c[4];
		C = c[1] / c[4];
		D = c[0] / c[4];

		// Substitute x = y - A/4: x^4 + px^2 + qx + r = 0
		sq_A = A * A;
		p = -3.0/8 * sq_A + B;
		q =  1.0/8 * sq_A * A - 1.0/2 * A * B + C;
		r = -3.0/256 * sq_A * sq_A + 1.0/16 * sq_A * B - 1.0/4 * A * C + D;

		if (is_zero(r)) {
			// No absolute term: y(y^3 + py + q) = 0
			coeffs[0] = q;
			coeffs[1] = p;
			coeffs[2] = 0;
			coeffs[3] = 1;
			num = solve_cubic(coeffs, s);
			s[num++] = 0;
		}
		else {
			// Solve resolvent cubic
			coeffs[0] = 1.0/2 * r * p - 1.0/8 * q * q;
			coeffs[1] = -r;
			coeffs[2] = -1.0/2 * p;
			coeffs[3] = 1;
			solve_cubic(coeffs, s);

			// Take real solution and build two quadric equations
			z = s[0];
			u = z * z - r;
			v = 2 * z - p;

			if (is_zero(u)) {
				u = 0;
			}
			else if (u > 0) {
				u = sqrt(u);
			}
			else {
				return 0;
			}

			if (is_zero(v)) {
				v = 0;
			}
			else if (v > 0) {
				v = sqrt(v);
			}
			else {
				return 0;
			}

			coeffs[0] = z - u;
			coeffs[1] = q < 0 ? -v : v;
			coeffs[2] = 1;
			num = solve_quadric(coeffs, s);

			coeffs[0] = z + u;
			coeffs[1] = q < 0 ? v : -v;
			coeffs[2] = 1;
			num += solve_quadric(coeffs, s + num);
		}

		// Resubstitute
		sub = 1.0/4 * A;
		for (i = 0; i < num; ++i) {
			s[i] -= sub;
		}

		return num;
	}

	// Function that solves a quadric equation, returns the number of real roots
	int solve_quadric(double c[3], double s[2]) const {
		double p, q, D;

		// Normal form: x^2 + px + q = 0
		p = c[1] / (2 * c[2]);
		q = c[0] / c[2];

		D = p * p - q;

		if (is_zero(D)) {
			s[0] = -p;
			return 1;
		}
		else if (D < 0) {
			return 0;
		}
		else {
			double sqrt_D = sqrt(D);
			s[0] = sqrt_D - p;
			s[1] = -sqrt_D - p;
			return 2;
		}
	}

	// Function that solves a cubic equation, returns the number of real roots
	int solve_cubic(double c[4], double s[3]) const {
		int i, num;
		double sub;
		double A, B, C;
		double sq_A, p, q;
		double cb_p, D;

		// Normal form: x^3 + Ax^2 + Bx + C = 0
		A = c[2] / c[3];
		B = c[1] / c[3];
		C = c[0] / c[3];

		// Substitute x = y - A/3: x^3 + px + q = 0
		sq_A = A * A;
		p = 1.0/3 * (-1.0/3 * sq_A + B);
		q = 1.0/2 * (2.0/27 * A * sq_A - 1.0/3 * A * B + C);

		// Use Cardano's formula
		cb_p = p * p * p;
		D = q * q + cb_p;

		if (is_zero(D)) {
			if (is_zero(q)) {
				s[0] = 0;
				num = 1;
			}
			else {
				double u = cbrt(-q);
				s[0] = 2 * u;
				s[1] = -u;
				num = 2;
			}
		}
		else if (D < 0) {
			double phi = 1.0/3 * acos(-q / sqrt(-cb_p));
			double t = 2 * sqrt(-p);

			s[0] = t * cos(phi);
			s[1] = -t * cos(phi + pi / 3);
			s[2] = -t * cos(phi - pi / 3);
			num = 3;
		}
		else {
			double sqrt_D = sqrt(D);
			double u = cbrt(sqrt_D - q);
			double v = -cbrt(sqrt_D + q);
			s[0] = u + v;
			num = 1;
		}

		// Resubstitute
		sub = 1.0/3 * A;
		for (i = 0; i < num; ++i) {
			s[i] -= sub;
		}

		return num;
	}
};

// Function that loads an OBJ file and adds mesh(es) to the scene with a certain material
shared_ptr<bvh_node> load_obj_file(const char* filename, shared_ptr<material> m, double scale = 1.0) {
	objl::Loader obj_loader;

	// Try opening file
	bool loaded = obj_loader.LoadFile(filename);

	if (!loaded) {
		std::cout << "\n\n" << current_date_time() << "Failed to load model from file '" << filename << "'\n";
		return nullptr;
	}
	std::cout << "\n\n" << current_date_time() << "Loading model from OBJ file '" << filename << "'...\n";

	// Load all meshes
	for (int i = 0; i < obj_loader.LoadedMeshes.size(); i++) {
		objl::Mesh current_mesh = obj_loader.LoadedMeshes[i];
		
		// Collect mesh vertices (in order, index important)
		std::vector<point3> verts(current_mesh.Vertices.size());
		std::vector<vec3> text_coords(current_mesh.Vertices.size());
		for (int j = 0; j < current_mesh.Vertices.size(); j++) {
			verts[j] = point3(current_mesh.Vertices[j].Position.X, current_mesh.Vertices[j].Position.Y, current_mesh.Vertices[j].Position.Z);
			text_coords[j] = vec3(current_mesh.Vertices[j].TextureCoordinate.X, current_mesh.Vertices[j].TextureCoordinate.Y, 0);
		}

		// Create a hittable list with the triangles as specified in the OBJ file
		hittable_list tris;
		for (int j = 0; j < current_mesh.Indices.size(); j += 3) {
			point3 p0 = verts[current_mesh.Indices[j]] * scale;
			point3 p1 = verts[current_mesh.Indices[j+1]] * scale;
			point3 p2 = verts[current_mesh.Indices[j+2]] * scale;
			vec3 uv0 = text_coords[current_mesh.Indices[j]];
			vec3 uv1 = text_coords[current_mesh.Indices[j+1]];
			vec3 uv2 = text_coords[current_mesh.Indices[j+2]];
			tris.add(make_shared<triangle>(p0, p1, p2, uv0, uv1, uv2, m));
		}
		std::cout << "                       - Adding mesh with " << verts.size() << " vertices (" << tris.size() << " triangles)\n";
		
		// Add the triangles, using a BVH tree (speeds up rendering, but takes longer when loading)
		return make_shared<bvh_node>(tris);
	}
}

// Function that loads an STL file and adds a mesh to the scene with a certain material
shared_ptr<bvh_node> load_stl_file(const char* filename, shared_ptr<material> m, double scale = 1.0) {
	try {
		// Try opening file
		stl_reader::StlMesh<double, unsigned int> mesh(filename);

		// Create a hittable list with the triangles as specified in the STL file
		std::cout << "\n\n" << current_date_time() << "Loading model from STL file '" << filename << "'...\n";
		hittable_list tris;
		for (size_t i = 0; i < mesh.num_tris(); ++i) {
			std::vector<point3> tri(3);
			for (size_t j = 0; j < 3; ++j) {
				const double* c = mesh.tri_corner_coords(i, j);
				tri[j] = point3(c[0] * scale, c[1] * scale, c[2] * scale);
			}
			tris.add(make_shared<triangle>(tri[0], tri[1], tri[2], m));
		}
		std::cout << "                       - Adding mesh with " << tris.size() << " triangles\n";

		//Add the triangles, using a BVH tree (speeds up rendering, but takes longer when loading)
		return make_shared<bvh_node>(tris);
	}
	catch (std::exception& e) {
		// Show error when model could not be loaded from file
		std::cout << "\n\n" << current_date_time() << "Failed to load model from file '" << filename << "': " << e.what() << "\n";
	}
}

#endif // !OBJECTS_H
