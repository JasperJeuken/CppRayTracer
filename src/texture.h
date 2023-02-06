/*
Various textures to be used in materials (solid color, checkered, perlin noise, ...)
*/
#ifndef TEXTURE_H
#define TEXTURE_H

#include "utility.h"
#include "rt_stb_image.h"

// Base class for textures (defines virtual "value" function that gives a color value for certain texture coordinates)
class texture {
public:
	virtual color value(double u, double v, const point3& p) const = 0;
};

// Solid color texture (single color value, returned for all texture coordinates)
class solid_color : public texture {
public:
	solid_color() {}
	solid_color(color c) : color_value(c) {}

	solid_color(double red, double green, double blue) : solid_color(color(red, green, blue)) {}

	// Return color value
	virtual color value(double u, double v, const vec3& p) const override {
		return color_value;
	}

private:
	color color_value;
};

// Checkered texture (two colors, returned in checker pattern based on texture coordinates)
class checker_texture : public texture {
public:
	checker_texture() {}
	checker_texture(shared_ptr<texture> _even, shared_ptr<texture> _odd) : even(_even), odd(_odd) {}

	checker_texture(color c1, color c2) : even(make_shared<solid_color>(c1)), odd(make_shared<solid_color>(c2)) {}

	// Return odd or even color based on point location
	virtual color value(double u, double v, const point3& p) const override {
		double sines = sin(10 * p.x())* sin(10 * p.y())* sin(10 * p.z());
		if (sines < 0)
			return odd->value(u, v, p);
		else
			return even->value(u, v, p);
	}

public:
	shared_ptr<texture> odd;
	shared_ptr<texture> even;
};

// Class that generates perlin noise for 3D points
class perlin {
public:
	// Create perlin noise
	perlin() {
		ranvec = new vec3[point_count];
		for (int i = 0; i < point_count; ++i) {
			ranvec[i] = unit_vector(vec3::random(-1, 1));
		}

		perm_x = perlin_generate_perm();
		perm_y = perlin_generate_perm();
		perm_z = perlin_generate_perm();
	}

	// Destructor
	~perlin() {
		delete[] ranvec;
		delete[] perm_x;
		delete[] perm_y;
		delete[] perm_z;
	}

	// Generate noise at a 3D point p
	double noise(const point3& p) const {
		double u = p.x() - floor(p.x());
		double v = p.y() - floor(p.y());
		double w = p.z() - floor(p.z());

		int i = static_cast<int>(floor(p.x()));
		int j = static_cast<int>(floor(p.y()));
		int k = static_cast<int>(floor(p.z()));
		vec3 c[2][2][2];

		for (int di = 0; di < 2; di++)
			for (int dj = 0; dj < 2; dj++)
				for (int dk = 0; dk < 2; dk++)
					c[di][dj][dk] = ranvec[
						perm_x[(i + di) & 255] ^ perm_y[(j + dj) & 255] ^ perm_z[(k + dk) & 255]
					];

		return perlin_interp(c, u, v, w);
	}

	// Turbulence function
	double turb(const point3& p, int depth = 7) const {
		double accum = 0.0;
		point3 temp_p = p;
		double weight = 1.0;

		for (int i = 0; i < depth; i++) {
			accum += weight * noise(temp_p);
			weight *= 0.5;
			temp_p *= 2;
		}

		return fabs(accum);
	}

private:
	static const int point_count = 256;
	vec3* ranvec;
	int* perm_x;
	int* perm_y;
	int* perm_z;

	// Generate perlin permutations
	static int* perlin_generate_perm() {
		int* p = new int[point_count];

		for (int i = 0; i < perlin::point_count; i++) {
			p[i] = i;
		}

		permute(p, point_count);

		return p;
	}

	// Permute function
	static void permute(int* p, int n) {
		for (int i = n - 1; i > 0; i--) {
			int target = random_int(0, i);
			int temp = p[i];
			p[i] = p[target];
			p[target] = temp;
		}
	}

	// Tri-linear interpolation
	static double trilinear_interp(double c[2][2][2], double u, double v, double w) {
		double accum = 0.0;
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				for (int k = 0; k < 2; k++)
					accum += (i * u + (1.0 - i) * (1 - u)) * (j * v + (1.0 - j) * (1 - v)) * (k * w + (1.0 - k) * (1 - w)) * c[i][j][k];
		return accum;
	}

	// Perlin interpolation
	static double perlin_interp(vec3 c[2][2][2], double u, double v, double w) {
		// Hermitian smoothing
		double uu = u * u * (3 - 2 * u);
		double vv = v * v * (3 - 2 * v);
		double ww = w * w * (3 - 2 * w);
		double accum = 0.0;

		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				for (int k = 0; k < 2; k++) {
					vec3 weight_v(u - i, v - j, w - k);
					accum += (i * uu + (1.0 - i) * (1 - uu))
						   * (j * vv + (1.0 - j) * (1 - vv))
						   * (k * ww + (1.0 - k) * (1 - ww))
						   * dot(c[i][j][k], weight_v);
				}
			}
		}
		return accum;
	}
};

// Perlin noise-based texture (based on texture coordinates)
class noise_texture : public texture {
public:
	noise_texture() {}
	noise_texture(double sc) : scale(sc) {}

	// Return color based on noise value at a point
	virtual color value(double u, double v, const point3& p) const override {
		return color(1, 1, 1) * 0.5 * (1 + sin(scale * p.z() + 10 * noise.turb(p)));
	}

public:
	perlin noise;
	double scale;
};

// Image texture class (selects image pixel based on texture coordinates)
class image_texture : public texture {
public:
	const static int bytes_per_pixel = 3;

	image_texture() : data(nullptr), width(0), height(0), bytes_per_scanline(0) {}

	// Load image from a file
	image_texture(const char* filename) {
		int components_per_pixel = bytes_per_pixel;

		data = stbi_load(filename, &width, &height, &components_per_pixel, components_per_pixel);

		if (!data) {
			std::cout << "Could not load image texture from file '" << filename << "'\n";
			width = height = 0;
		}

		bytes_per_scanline = bytes_per_pixel * width;
	}

	// Destructor
	~image_texture() {
		delete data;
	}

	// Return pixel color value at UV-coordinates
	virtual color value(double u, double v, const vec3& p) const override {
		// Debug color (no image data)
		if (data == nullptr)
			return color(0, 1, 1);

		// Clamp coordinates to [0, 1]
		u = clamp(u, 0.0, 1.0);
		v = 1.0 - clamp(v, 0.0, 1.0);

		// Clamp integer mapping
		int i = static_cast<int>(u * width);
		int j = static_cast<int>(v * height);
		if (i >= width) i = width - 1;
		if (j >= height) j = height - 1;

		const double color_scale = 1.0 / 255.0;
		auto pixel = data + j * bytes_per_scanline + i * bytes_per_pixel;

		return color(color_scale * pixel[0], color_scale * pixel[1], color_scale * pixel[2]);
	}

private:
	unsigned char* data;
	int width, height;
	int bytes_per_scanline;
};

// Gradient texture class (gradient between two colors, using an origin point and a direction vector, and a limit to determine the spread)
class gradient_texture : public texture {
public:
	gradient_texture() {};
	gradient_texture(color _c1, color _c2) : c1(_c1), c2(_c2), orig(point3(0, 0, 0)), dir(vec3(0, 1, 0)), limit(1) {}
	gradient_texture(color _c1, color _c2, point3 _orig, vec3 _dir, double _limit = 1) : c1(_c1), c2(_c2), orig(_orig), dir(unit_vector(_dir)), limit(_limit) {}

	// Return either color 1 or color 2, or a mix of them if in the gradient range
	virtual color value(double u, double v, const vec3& p) const override {
		double dist = dot(dir, p - orig);
		if (dist > limit)
			return c1;

		if (fabs(dist) > limit)
			return c2;

		double t = (dist + limit) / (2 * limit);
		return t * c1 + (1 - t) * c2;
	}

public:
	color c1, c2;
	point3 orig;
	vec3 dir;
	double limit;
};

#endif // !TEXTURE_H
