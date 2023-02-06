/*
Various materials to be applied to hittable objects (such as Lambertian, metal, dielectric, ...)
*/
#ifndef MATERIAL_H
#define MATERIAL_H

#include "utility.h"
#include "hittable.h"
#include "texture.h"

struct hit_record;

// Base class for a material. Defines virtual functions to be overridden by derived classes
class material {
public:
	virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const = 0;
	virtual color emitted(double u, double v, const point3& p) const { return color(0, 0, 0); }
};

// Lambertian diffuse material, reflecting rays diffusely (in all directions)
class lambertian : public material {
public:
	lambertian(const color& a) : albedo(make_shared<solid_color>(a)) {}
	lambertian(shared_ptr<texture> a) : albedo(a) {}

	// Random-based scatter
	virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const override {
		vec3 scatter_direction = rec.normal + random_unit_vector();

		// Catch degenerate scatter direction
		if (scatter_direction.near_zero())
			scatter_direction = rec.normal;

		scattered = ray(rec.p, scatter_direction, r_in.time());
		attenuation = albedo->value(rec.u, rec.v, rec.p);
		return true;
	}

public:
	shared_ptr<texture> albedo;
};

// Metal material, reflecting rays (with optional fuzz/randomness)
class metal : public material {
public:
	metal(const color& a) : albedo(a), fuzz(0) {}
	metal(const color& a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}
	metal(const color& a, const char* filename) : albedo(a) {
		roughness_map = make_shared<image_texture>(filename);
		has_roughness_map = true;
	}

	// Scatter with perfect reflection (unless a fuzz has been set, or a roughness map was loaded)
	virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const override {
		vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);

		// Create scattered ray based on roughness map or reflection with fuzz
		if (!has_roughness_map) {
			scattered = ray(rec.p, reflected + fuzz * random_in_unit_sphere(), r_in.time());
		}
		else {
			double roughness_value = roughness_map->value(rec.u, rec.v, rec.p).length() / 3;
			scattered = ray(rec.p, reflected + roughness_value * random_in_unit_sphere(), r_in.time());
		}

		vec3 unit_direction = unit_vector(r_in.direction());
		double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
		attenuation = metal_reflectance(cos_theta, albedo);
		return (dot(scattered.direction(), rec.normal) > 0);
	}
public:
	color albedo;
	double fuzz;
	bool has_roughness_map = false;
	shared_ptr<image_texture> roughness_map;

private:
	// Apply Schlick's approximation for high-incidence rays
	static color metal_reflectance(double cosine, const color& albedo) {
		return albedo + (color(1, 1, 1) - albedo) * pow(1.0 - cosine, 5);
	}
};

// Dielectric material, reflects or refracts rays depending on ray angles (glass, water, ...)
class dielectric : public material {
public:
	dielectric(double index_of_refraction) : ir(index_of_refraction) {}

	// Reflect or refract ray based on ray angle
	virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const override {
		attenuation = color(1.0, 1.0, 1.0);
		double refraction_ratio = rec.front_face ? (1.0 / ir) : ir;

		vec3 unit_direction = unit_vector(r_in.direction());
		double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
		double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

		bool cannot_refract = refraction_ratio * sin_theta > 1.0;
		vec3 direction;

		if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
			direction = reflect(unit_direction, rec.normal);
		else
			direction = refract(unit_direction, rec.normal, refraction_ratio);
				
		scattered = ray(rec.p, direction, r_in.time());
		return true;
	}

public:
	double ir;

private:
	// Apply Schlick's approximation for high-incidence rays
	static double reflectance(double cosine, double ref_idx) {
		double r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
		r0 = r0 * r0;
		return r0 + (1 - r0) * pow(1.0 - cosine, 5);
	}
};

// Diffuse emissive material with a certain color
class diffuse_light : public material {
public:
	diffuse_light(shared_ptr<texture> a) : emit(a) {}
	diffuse_light(color c) : emit(make_shared<solid_color>(c)) {}

	// No scatter
	virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const override {
		return false;
	}

	// Emit color value from light texture
	virtual color emitted(double u, double v, const point3& p) const override {
		return emit->value(u, v, p);
	}

public:
	shared_ptr<texture> emit;
};

// Isotropic material with a certain color
class isotropic : public material {
public:
	isotropic(color c) : albedo(make_shared<solid_color>(c)) {}
	isotropic(shared_ptr<texture> a) : albedo(a) {}

	// Scatter ray in a random direction
	virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const override {
		scattered = ray(rec.p, random_in_unit_sphere(), r_in.time());
		attenuation = albedo->value(rec.u, rec.v, rec.p);
		return true;
	}

public:
	shared_ptr<texture> albedo;
};

#endif // !MATERIAL_H
