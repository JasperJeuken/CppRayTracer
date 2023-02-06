/*
File holding pre-defined scenes
*/
#ifndef SCENES_H
#define SCENES_H

#include "render_scene.h"
#include "utility.h"

// Function where a scene is defined and populated with objects
scene your_scene() {
    scene sc(16.0 / 9.0); // Create a scene with an aspect ratio of 16:9

    // Create a camera
    point3 lookfrom(0, 0, -3); // Camera position (x, y, z)
    point3 lookat(0, 0, 0); // Camera target (x, y, z)
    double vfov = 30.0; // Vertical field-of-view (degrees)
    sc.set_camera(lookfrom, lookat, vfov);

    // Create the environment
    sc.set_sky_environment(); // Sky colored gradient

    // Create Lambertian diffuse material (r, g, b)
    auto sphere_material = sc.create_lambertian_material(color(0.9, 0.2, 0.1));

    // Add two spheres next to each other, with a radius of 0.5
    sc.add_sphere(point3(-0.5, 0, 0), 0.5, sphere_material);
    sc.add_sphere(point3(0.5, 0, 0), 0.5, sphere_material);

    return sc;
}

// Example: Cornell box scene
scene cornell_box() {
    // Create scene with aspect ratio
    scene sc(1.0);

    // Set camera and environment
    point3 lookfrom(278, 278, -800);
    point3 lookat(278, 278, 0);
    double vfov = 40.0;
    sc.set_camera(lookfrom, lookat, vfov);
    sc.set_environment(color(0, 0, 0));

    // Create objects in scene
    auto red = sc.create_lambertian_material(color(.65, 0.05, 0.05));
    auto white = sc.create_lambertian_material(color(0.73, 0.73, 0.73));
    auto green = sc.create_lambertian_material(color(0.12, 0.45, 0.15));
    auto light = sc.create_diffuse_light_material(color(15, 15, 15));

    sc.add_yz_rect(0, 555, 0, 555, 555, green);
    sc.add_yz_rect(0, 555, 0, 555, 0, red);
    sc.add_xz_rect(213, 343, 227, 332, 554, light);
    sc.add_xz_rect(0, 555, 0, 555, 0, white);
    sc.add_xz_rect(0, 555, 0, 555, 555, white);
    sc.add_xy_rect(0, 555, 0, 555, 555, white);

    auto box1 = sc.rotate_object_y(sc.create_box(point3(0, 0, 0), point3(165, 330, 165), white), 15);
    box1 = sc.translate_object(box1, vec3(265, 0, 295));
    sc.add(box1);

    auto box2 = sc.rotate_object_y(sc.create_box(point3(0, 0, 0), point3(165, 165, 165), white), -18);
    box2 = sc.translate_object(box2, vec3(130, 0, 65));
    sc.add(box2);

    return sc;
}

#endif // !SCENES_H
