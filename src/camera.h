/*
Camera class that generates rays based on the camera's properties
*/
#ifndef CAMERA_H
#define CAMERA_H

#include "utility.h"

// Camera base class
class camera {
public:
	camera() {}
	virtual ray get_ray(double s, double t) const = 0;
};

// Perspective camera class
class perspective_camera : public camera {
public:
	perspective_camera() {}

	// Constructor that calculates camera geometry and translates into vectors used for ray generation
	perspective_camera(point3 lookfrom, point3 lookat, double vfov, double aspect_ratio, vec3 vup=vec3(0, 1, 0), double aperture=0.0, double focus_dist=1.0, double _time0 = 0, double _time1 = 0) {
		// Store values
		double theta = degrees_to_radians(vfov);
		double h = tan(theta / 2);
		double viewport_height = 2.0 * h;
		double viewport_width = aspect_ratio * viewport_height;
		origin = lookfrom;
		lens_radius = aperture / 2;
		time0 = _time0;
		time1 = _time1;
		
		// Create camera vectors
		w = unit_vector(lookfrom - lookat);
		u = unit_vector(cross(vup, w));
		v = cross(w, u);
		horizontal = focus_dist * viewport_width * u;
		vertical = focus_dist * viewport_height * v;
		lower_left_corner = origin - horizontal / 2 - vertical / 2 - focus_dist * w;
	}

	// Function that generates a ray from the camera origin to a point (s, t) on the viewport
	virtual ray get_ray(double s, double t) const override {
		vec3 rd = lens_radius * random_in_unit_disk();
		vec3 offset = u * rd.x() + v * rd.y();

		return ray(origin + offset, lower_left_corner + s * horizontal + t * vertical - origin - offset, random_double(time0, time1));
	}

private:
	point3 origin;
	point3 lower_left_corner;
	vec3 horizontal;
	vec3 vertical;
	vec3 u, v, w;
	double lens_radius;
	double time0, time1;
};

// Orthographic camera class
class orthographic_camera : public camera {
public:
	orthographic_camera(){}

	orthographic_camera(point3 lookfrom, point3 lookat, double _width, double _height, vec3 vup = vec3(0, 1, 0), double t0 = 0, double t1 = 0) {
		// Store values
		width = _width;
		height = _height;
		time0 = t0;
		time1 = t1;

		// Create camera vectors
		w = unit_vector(lookat - lookfrom); // look direction
		u = unit_vector(cross(vup, w));
		v = unit_vector(cross(w, u));
		horizontal = width * u;
		vertical = height * v;
		lower_left_corner = lookfrom - horizontal / 2 - vertical / 2;
	}

	// Function that generates a ray from a point on the virtual frame in the orthographic direction
	virtual ray get_ray(double s, double t) const override {
		return ray(lower_left_corner + s * horizontal + t * vertical, w, random_double(time0, time1));
	}

private:
	point3 lower_left_corner;
	vec3 u, v, w;
	vec3 horizontal, vertical;
	double width, height;
	double time0, time1;
};

#endif // !CAMERA_H
