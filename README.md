# C++ Ray Tracer
A ray tracing engine that uses only standard C++ libraries for rendering.
[Click here to see the gallery!](#gallery)

For a detailed guide on how to create scenes and render them, check out the [wiki](https://github.com/JasperJeuken/CppRayTracer/wiki)!

_Includes a live preview window of the multithreaded render_
<img src="images/render_preview.gif" alt="Render preview window" width=80%>

# Features
- **Monte Carlo ray tracing**
- **Integrated denoising using Intel® Open Image Denoise**
- **Improved rendering:**
  - Multi-threaded rendering
  - Live preview window using SDL
  - **Different render passes:**
    - Color, albedo, normal, depth, opacity, emission, ...
  - Object normal maps and roughness maps
  - Simulated motion blur
- **Customizable camera:**
  - Perspective camera (depth of field, aperture, focus distance, field of view, ...)
  - Orthographic camera (position, projection plane size, ...)
- **Easy-to-use scene creation**
- **Customizable environment texture**
- **Object bounding volume hierarchy for improved performance**

### Supports:
- **Objects:**
  - OBJ or STL files
  - **Primitives**:
    - Sphere, planes, torus, triangle, cylinder, disk, constant medium, box, ...
- **Textures:**
  - Solid color, checkered, perlin noise, image, gradient, ...
- **Materials**
  - Diffuse Lambertian, metallic, dielectric, isotropic, ...
- **Image formats:**
  - JPG, PNG, PPM, ...

# Gallery
A collection of some sample scenes rendered using the engine. The following images were all rendered using 10 threads.

### Orthographic scene
*1000x1000px, 6.1 seconds, 50 samples, denoised*
<img src="https://i.imgur.com/Vshhf8i.png" width=100%>

| Raw pass (not denoised)              | Albedo pass                          | Emission pass                        |
|:------------------------------------:|:------------------------------------:|:------------------------------------:|
| ![](https://i.imgur.com/u5GwZtB.png) | ![](https://i.imgur.com/4kIUE8Y.png) | ![](https://i.imgur.com/Do7QjXm.png) |

| Normal pass                          | Depth pass                           | Opacity pass                         |
|:------------------------------------:|:------------------------------------:|:------------------------------------:|
| ![](https://i.imgur.com/Ox5BrlU.png) | ![](https://i.imgur.com/y0gJipm.png) | ![](https://i.imgur.com/7E3iCCR.png) |

### Cornell box
*1000x1000px, 49.6 seconds, 50 samples, denoised*
<img src="https://i.imgur.com/V3rU7Lz.png" width=100%>

### Glass monkey
*1500x1500px, 217.3 seconds, 100 samples, denoised (model: <a href="https://www.dummies.com/article/technology/software/animation-software/blender/meet-suzanne-the-blender-monkey-142918/">Suzanne the monkey</a>)*
<img src="https://i.imgur.com/ZH7Z7T3.png" width=100%>

### Normal map
*1000x1000px, 100 samples, denoised*

| Result with normal map               | Result without normal map            |
|:------------------------------------:|:------------------------------------:|
| 19.0 seconds                         | 15.8 seconds                         |
| ![](https://i.imgur.com/2sCAkkj.png) | ![](https://i.imgur.com/HsM4zb5.png) |
| ![](https://i.imgur.com/cNC9Ehn.png) | ![](https://i.imgur.com/hhRxNob.png) |

### Metal roughness
*2500x1000px, 43.2 seconds, 50 samples, denoised*
<img src="https://i.imgur.com/xwa2WjS.png" width=100%>

### All primitives
*3000x333px, 14.1 seconds, 50 samples, denoised*
<img src="https://i.imgur.com/B0hAhp8.png" width=100%>
Left to right: 3D model, torus, cylinder, triangle, disk, volume, box, x-y-z-rectangles, motion-blurred sphere, sphere

### Roughness map
*1500x1500px, 792 seconds, 1000 samples, not denoised*
<img src="https://i.imgur.com/qdkLm8g.png" width=100%>

# Using the engine
## Install
To install the engine, follow these steps:
1. Clone this repository
2. Install <a href="https://www.libsdl.org/">Simple DirectMedia Layer</a> (for live render preview)
3. Install <a href="https://www.openimagedenoise.org/">Intel® Open Image Denoise</a> (for denoising)

## Creating a scene
Create a scene in `scenes.h` by defining a new function, which specifies the camera properties, the environment, and the objects in the scene. For example:
```c++
scene your_scene() {
  scene sc(16.0 / 9.0); // Create a scene with an aspect ratio of 16:9

  // Create a camera
  point3 lookfrom(0, 0, -5); // Camera position (x, y, z)
  point3 lookat(0, 0, 0); // Camera target (x, y, z)
  double vfov = 30.0; // Vertical field-of-view (degrees)
  sc.set_camera(lookfrom, lookat, vfov);

  // Create the environment
  sc.set_sky_environment(); // Sky colored gradient

  // Create Lambertian diffuse material (r, g, b)
  auto sphere_material = sc.create_lambertian_material(color(0.9, 0.2, 0.1));
  // Add two spheres next to each other, with a radius of 0.5
  sc.add_sphere(point3(-0.5, 0, 0), 0.5, sphere_material);
  sc.add_sphere(point3( 0.5, 0, 0), 0.5, sphere_material);

  return sc;
}
```
This results in a scene with two spheres side-by-side, surrounded by an environment that resembles a blue sky.
<!-- TODO: ADD DOCUMENTATION FOR AVAILABLE FUNCTIONS -->

## Rendering a scene
In your main file, include:
```c++
#include "scenes.h" // Your scenes
#include "render.h" // Render functions
```

We load the scene, and specify some image settings. Note that you can either specify the image width, height, or specify both. The aspect ratio is used to calculate the unknown length if only width or height is specified.
```c++
...
int main() {
  // Load scene
  scene sc = your_scene();

  // Set render settings
  sc.settings.set_image_width(1000);
  sc.settings.set_samples_per_pixel(100);
  sc.settings.set_max_ray_bounces(50);
  ...
```
To perform a render pass, we create an image buffer, and render the scene to it:
```c++
  ...
  // Perform a color render pass
  auto buffer = create_buffer(sc);
  render(sc, buffer, "color");
  ...
```
This starts a render, and will open a window with the live preview. The program automatically detects the number of available threads, and starts workers on each of them.

Instead of a `color` pass, we can also specify one of the following: `normal`, `albedo`, `depth`, `opacity`, or `emission`.

After rendering the scene to the buffer, we can save the image to a file:
```c++
  ...
  // Save to JPG
  save_image_to_file(sc, buffer);
```
This saves the image to a JPG file. The format can alternatively be specified as `ppm` or `png`. This results in the following image:
<img src="https://i.imgur.com/0fQnZi3.jpg" width=100%>


## Denoising an image
The image is noisy, and we can improve the quality by denoising it. To denoise an image using Intel® Open Image Denoise, we need only the `color` pass. However, the results improve if we also supply an `albedo` and `normal` pass. We render these as follows:
```c++
  ...
  // Render three passes (color, normal, albedo)
  auto color_buffer = create_buffer(sc);
  auto albedo_buffer = create_buffer(sc);
  auto normal_buffer = create_buffer(sc);
  render(sc, color_buffer, "color");
  render(sc, albedo_buffer, "albedo");
  render(sc, normal_buffer, "normal");
  ...
```
Next, we simply create another buffer for the denoised result, and run the `denoise` function:
```c++
  // Denoise the image
  auto denoised_buffer = create_buffer(sc);
  denoise(sc, denoised_buffer, color_buffer, albedo_buffer, normal_buffer);
```
Like before, we can save the buffer to a file:
```c++
  // Save to JPG
  save_image_to_file(sc, denoised_buffer, "denoised");
```
Here, the `denoised` string is a suffix for the file name, ensuring we don't overwrite the color pass we saved earlier. This results in the following image:
<img src="https://i.imgur.com/IOnQnoI.jpg" width=100%>

# Attribution
Engine basis from "Ray Tracing in One Weekend Series" (v3.2.0) for C++ by Peter Shirley (Steve Hollasch, Trevor David Black).

Uses <a href="https://www.libsdl.org/">Simple DirectMedia Layer</a> for live render preview.

Uses <a href="https://www.openimagedenoise.org/">Intel® Open Image Denoise</a> for denoising.

Uses <a href="https://github.com/nothings/stb">stb</a> (`stb_image` and `stb_image_write`) for image loading/saving.

Uses <a href="https://github.com/Bly7/OBJ-Loader">OBJ-Loader</a> for loading OBJ files.

Uses <a href="https://github.com/sreiter/stl_reader">stl_reader</a> for loading STL files.
