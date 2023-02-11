/*
Classes used for storing all information needed for rendering a scene
*/
#ifndef RENDER_SCENE
#define RENDER_SCENE

#include "utility.h"
#include "environment.h"
#include "camera.h"
#include "material.h"
#include "texture.h"
#include "objects.h"
#include "hittable_list.h"


// Class that stores render settings (image size, samples, ...)
struct render_settings {
    render_settings() {}

    render_settings(double ar) : aspect_ratio(ar) {}

    render_settings(double ar, int iw, int d, int s) : aspect_ratio(ar), image_width(iw), max_depth(d), samples_per_pixel(s) {}

    // Set image width and recalculate image height based on aspect ratio
    void set_image_width(int iw) {
        image_width = iw;
        image_height = (int)(image_width / aspect_ratio);
    }

    // Set image height and recalculate image width based on aspect ratio
    void set_image_height(int ih) {
        image_height = ih;
        image_width = (int)(image_height * aspect_ratio);
    }

    // Set image width and height and recalculate aspect ratio
    void set_image_size(int iw, int ih) {
        image_width = iw;
        image_height = ih;
        aspect_ratio = (double)iw / (double)ih;
    }

    // Set aspect ratio and recalculate image height based on image width
    void set_aspect_ratio(int ar) {
        aspect_ratio = ar;
        image_height = (int)(image_width * aspect_ratio);
    }

    // Set number of samples per pixel
    void set_samples_per_pixel(int n_samples) {
        samples_per_pixel = n_samples;
    }

    // Set maximum ray depth
    void set_max_ray_bounces(int depth) {
        max_depth = depth;
    }

    // Set filename
    void set_filename(std::string new_filename) {
        filename = new_filename;
    }

public:
    double aspect_ratio;
    int image_width;
    int image_height;
    int max_depth;
    int samples_per_pixel;
    std::string filename;
};

// Class that stores a camera, list of hittable objects (world), and render settings
class scene {
public:
    scene() {
        settings = render_settings(1.0);
        set_camera(point3(0, 0, -1), point3(0, 0, 0), 40.0);
        set_sky_environment();
    }

    scene(double aspect_ratio) {
        settings = render_settings(aspect_ratio);
        set_camera(point3(0, 0, -1), point3(0, 0, 0), 40.0);
        set_sky_environment();
    }

    // Create a new camera using camera properties
    void set_camera(point3 lookfrom, point3 lookat, double vfov, double aperture = 0.0, double focus_dist = 1.0, double _time0 = 0, double _time1 = 0, vec3 vup = vec3(0, 1, 0)) {
        cam = make_shared<perspective_camera>(lookfrom, lookat, vfov, settings.aspect_ratio, vup, aperture, focus_dist, _time0, _time1);
    }

    // Create a new orthographic camera using camera properties
    void set_orthographic_camera(point3 lookfrom, point3 lookat, double width, vec3 vup = vec3(0, 1, 0), double time0 = 0, double time1 = 0) {
        cam = make_shared<orthographic_camera>(lookfrom, lookat, width, width / settings.aspect_ratio, vup, time0, time1);
    }

    // Create a new environment from a texture
    void set_environment(shared_ptr<texture> t) {
        env = environment(t);
    }

    // Create a new environment from a color
    void set_environment(color c) {
        env = environment(c);
    }

    // Set the environment to an existing one
    void set_environment(environment e) {
        env = e;
    }

    // Create an environment that represents a blue sky
    void set_sky_environment() {
        env = environment(make_shared<gradient_texture>(color(0.5, 0.7, 1.0), color(1.0, 1.0, 1.0)));
    }

    // Helper functions for creating textures
    shared_ptr<solid_color> create_solid_color_texture(color c) { return make_shared<solid_color>(c); }
    shared_ptr<checker_texture> create_checker_texture(color even, color odd, double s = 10.0) { return make_shared<checker_texture>(even, odd, s); }
    shared_ptr<checker_texture> create_checker_texture(shared_ptr<texture> even, shared_ptr<texture> odd, double s = 10.0) { return make_shared<checker_texture>(even, odd, s); }
    shared_ptr<noise_texture> create_noise_texture(double scale = 1.0) { return make_shared<noise_texture>(scale); }
    shared_ptr<image_texture> create_image_texture(const char* filename) { return make_shared<image_texture>(filename); }
    shared_ptr<gradient_texture> create_gradient_texture(color c1, color c2) { return make_shared<gradient_texture>(c1, c2); }
    shared_ptr<gradient_texture> create_gradient_texture(color c1, color c2, point3 origin, vec3 direction, double limit = 1.0) { return make_shared<gradient_texture>(c1, c2, origin, direction, limit); }

    // Helper functions for creating materials
    shared_ptr<lambertian> create_lambertian_material(color c) { return make_shared<lambertian>(c); }
    shared_ptr<lambertian> create_lambertian_material(shared_ptr<texture> text) { return make_shared<lambertian>(text); }
    shared_ptr<metal> create_metal_material(color c) { return make_shared<metal>(c); }
    shared_ptr<metal> create_metal_material(color c, double fuzz) { return make_shared<metal>(c, fuzz); }
    shared_ptr<metal> create_metal_material(color c, const char* filename) { return make_shared<metal>(c, filename); }
    shared_ptr<dielectric> create_dielectric_material(double index_of_refraction) { return make_shared<dielectric>(index_of_refraction); }
    shared_ptr<diffuse_light> create_diffuse_light_material(color c, double s = 1.0) { return make_shared<diffuse_light>(c, s); }
    shared_ptr<diffuse_light> create_diffuse_light_material(shared_ptr<texture> text, double s = 1.0) { return make_shared<diffuse_light>(text, s); }
    shared_ptr<isotropic> create_isotropic_material(color c) { return make_shared<isotropic>(c); }
    shared_ptr<isotropic> create_isotropic_material(shared_ptr<texture> text) { return make_shared<isotropic>(text); }

    // Helper functions for creating objects
    shared_ptr<sphere> create_sphere(point3 center, double radius, shared_ptr<material> mat) { return make_shared<sphere>(center, radius, mat); }
    shared_ptr<moving_sphere> create_moving_sphere(point3 center1, point3 center2, double time1, double time2, double radius, shared_ptr<material> mat) { return make_shared<moving_sphere>(center1, center2, time1, time2, radius, mat); }
    shared_ptr<xy_rect> create_xy_rect(double x_min, double x_max, double y_min, double y_max, double z, shared_ptr<material> mat) { return make_shared<xy_rect>(x_min, x_max, y_min, y_max, z, mat); }
    shared_ptr<yz_rect> create_yz_rect(double y_min, double y_max, double z_min, double z_max, double x, shared_ptr<material> mat) { return make_shared<yz_rect>(y_min, y_max, z_min, z_max, x, mat); }
    shared_ptr<xz_rect> create_xz_rect(double x_min, double x_max, double z_min, double z_max, double y, shared_ptr<material> mat) { return make_shared<xz_rect>(x_min, x_max, z_min, z_max, y, mat); }
    shared_ptr<box> create_box(point3 corner1, point3 corner2, shared_ptr<material> mat) { return make_shared<box>(corner1, corner2, mat); }
    shared_ptr<constant_medium> create_constant_medium(shared_ptr<hittable> boundary, double density, shared_ptr<texture> text) { return make_shared<constant_medium>(boundary, density, text); }
    shared_ptr<constant_medium> create_constant_medium(shared_ptr<hittable> boundary, double density, color c) { return make_shared<constant_medium>(boundary, density, c); }
    shared_ptr<disk> create_disk(point3 center, vec3 normal, double radius, shared_ptr<material> mat) { return make_shared<disk>(center, normal, radius, mat); }
    shared_ptr<disk> create_disk(point3 center, vec3 normal, vec3 right, double radius, shared_ptr<material> mat) { return make_shared<disk>(center, normal, right, radius, mat); }
    shared_ptr<triangle> create_triangle(point3 p1, point3 p2, point3 p3, shared_ptr<material> mat, bool one_sided = false) { return make_shared<triangle>(p1, p2, p3, mat, one_sided); }
    shared_ptr<cylinder> create_cylinder(double radius, double length, shared_ptr<material> mat) { return make_shared<cylinder>(radius, length, mat); }
    shared_ptr<cylinder> create_cylinder(double radius, double length, point3 center, shared_ptr<material> mat) { return make_shared<cylinder>(radius, length, center, mat); }
    shared_ptr<torus> create_torus(double major_radius, double minor_radius, shared_ptr<material> mat) { return make_shared<torus>(major_radius, minor_radius, mat); }
    shared_ptr<bvh_node> create_model(const char* filename, shared_ptr<material> mat, double scale=1.0) {
        std::string name = filename;
        if (ends_with(name, ".obj")) return load_obj_file(filename, mat, scale);
        if (ends_with(name, ".stl")) return load_stl_file(filename, mat, scale);
        throw std::runtime_error("This file type is not supported");
    }

    // Helper functions for adding objects
    void add(shared_ptr<hittable> object) { world.add(object); }
    void add_sphere(point3 center, double radius, shared_ptr<material> mat) { world.add(create_sphere(center, radius, mat)); }
    void add_moving_sphere(point3 center1, point3 center2, double time1, double time2, double radius, shared_ptr<material> mat) {world.add(create_moving_sphere(center1, center2, time1, time2, radius, mat)); }
    void add_xy_rect(double x_min, double x_max, double y_min, double y_max, double z, shared_ptr<material> mat) { world.add(create_xy_rect(x_min, x_max, y_min, y_max, z, mat)); }
    void add_yz_rect(double y_min, double y_max, double z_min, double z_max, double x, shared_ptr<material> mat) { world.add(create_yz_rect(y_min, y_max, z_min, z_max, x, mat)); }
    void add_xz_rect(double x_min, double x_max, double z_min, double z_max, double y, shared_ptr<material> mat) { world.add(create_xz_rect(x_min, x_max, z_min, z_max, y, mat)); }
    void add_box(point3 corner1, point3 corner2, shared_ptr<material> mat) { world.add(create_box(corner1, corner2, mat)); }
    void add_constant_medium(shared_ptr<hittable> boundary, double density, shared_ptr<texture> text) { world.add(create_constant_medium(boundary, density, text)); }
    void add_constant_medium(shared_ptr<hittable> boundary, double density, color c) { world.add(create_constant_medium(boundary, density, c)); }
    void add_disk(point3 center, vec3 normal, double radius, shared_ptr<material> mat) { world.add(create_disk(center, normal, radius, mat)); }
    void add_disk(point3 center, vec3 normal, vec3 right, double radius, shared_ptr<material> mat) { world.add(create_disk(center, normal, right, radius, mat)); }
    void add_triangle(point3 p1, point3 p2, point3 p3, shared_ptr<material> mat, bool one_sided = false) { world.add(create_triangle(p1, p2, p3, mat, one_sided)); }
    void add_cylinder(double radius, double length, shared_ptr<material> mat) { world.add(create_cylinder(radius, length, mat)); }
    void add_cylinder(double radius, double length, point3 center, shared_ptr<material> mat) { world.add(create_cylinder(radius, length, center, mat)); }
    void add_torus(double major_radius, double minor_radius, shared_ptr<material> mat) { world.add(create_torus(major_radius, minor_radius, mat)); }
    void add_model(const char* filename, shared_ptr<material> mat, double scale=1.0) { world.add(create_model(filename, mat, scale)); }

    // Helper functions for rotating and translating objects
    shared_ptr<hittable> rotate_object_x(shared_ptr<hittable> obj, double angle) { return make_shared<rotate_x>(obj, angle); }
    shared_ptr<hittable> rotate_object_y(shared_ptr<hittable> obj, double angle) { return make_shared<rotate_y>(obj, angle); }
    shared_ptr<hittable> rotate_object_z(shared_ptr<hittable> obj, double angle) { return make_shared<rotate_z>(obj, angle); }
    shared_ptr<hittable> rotate_object(shared_ptr<hittable> obj, vec3 angles) { return make_shared<rotate>(obj, angles); }
    shared_ptr<hittable> translate_object(shared_ptr<hittable> obj, vec3 displacement) { return make_shared<translate>(obj, displacement); }

public:
    shared_ptr<camera> cam;
    environment env;
    hittable_list world;
    render_settings settings;
};

#endif // !RENDER_SCENE