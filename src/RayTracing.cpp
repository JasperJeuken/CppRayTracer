/*
Main function that shows how to render a scene with three passes and denoise the resulting image
*/
#include "scenes.h"
#include "render.h"
#include "denoise.h"

// Main function where a scene is loaded, rendered, denoised, and saved to a file
int main()
{
    // Load scene
    scene sc = cornell_box();

    // Set render settings
    sc.settings.set_image_width(1000);
    sc.settings.set_samples_per_pixel(50);
    sc.settings.set_max_ray_bounces(50);

    // Render three passes (color, albedo, normal)
    auto color_buffer = create_buffer(sc);
    auto albedo_buffer = create_buffer(sc);
    auto normal_buffer = create_buffer(sc);
    render(sc, color_buffer, "color");
    render(sc, albedo_buffer, "albedo");
    render(sc, normal_buffer, "normal");

    // Denoise the image
    auto denoised_buffer = create_buffer(sc);
    denoise(sc, denoised_buffer, color_buffer, albedo_buffer, normal_buffer);

    // Save to JPG
    save_image_to_file(sc, color_buffer, "color");       // "Raw" image
    save_image_to_file(sc, denoised_buffer, "denoised"); // Denoised image

    return 0;
}
