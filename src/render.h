/*
Various functions used for rendering a scene (multithreaded ray tracing, separating image into render regions, saving image to file)
*/
#ifndef RENDER_H
#define RENDER_H

#include <vector>
#include <mutex>
#include <atomic>
#include <thread>
#include <fstream>
#include <ctime>
#include <vector>
#include <algorithm>
#include <utility>
#include <filesystem>
#include <chrono>

#include "rt_stb_image.h"
#include "utility.h"
#include "scenes.h"


// Variables used for parallelization
std::atomic<int> done_count{ 0 };
std::atomic<int> ray_count{ 0 };
std::vector<bool> threads_done;
static std::mutex m;

static std::vector<int> DEFAULT_INT_VECTOR;

/*
    RENDER REGIONS
*/
// Struct representing a region of the image to be rendered (x- and y-bounds, and number of samples)
struct region {
	region() {}
	region(int x0, int x1, int y0, int y1, int s) : x_start(x0), x_end(x1), y_start(y0), y_end(y1), samples(s) {}

	int x_start, x_end;
	int y_start, y_end;
	int samples;
};

// Create an image buffer based on a scene (image width and height)
std::vector<color> create_buffer(const scene& sc) {
	return std::vector<color>(sc.settings.image_width * sc.settings.image_height);
}

// Create a list of regions by splitting up the image into horizontal lines of pixels
std::vector<region> create_horizontal_regions(const scene& sc) {
	std::vector<region> regions(sc.settings.image_height);
	for (int i = sc.settings.image_height - 1; i >= 0; i--) {
		regions[i] = region(0, sc.settings.image_width, i, i + 1, sc.settings.samples_per_pixel);
	}
	return regions;
}

// Create a list of regions by splitting up the image into vertical lines of pixels
std::vector<region> create_vertical_regions(const scene& sc) {
	std::vector<region> regions(sc.settings.image_width);
	for (int i = sc.settings.image_width - 1; i >= 0; i--) {
		regions[i] = region(i, i + 1, 0, sc.settings.image_height, sc.settings.samples_per_pixel);
	}
	return regions;
}

// Create a list of regions by splitting up the image into sample passes (full image, one sample at a time)
std::vector<region> create_sample_regions(const scene& sc) {
	std::vector<region> regions(sc.settings.samples_per_pixel);
	for (int i = 0; i < sc.settings.samples_per_pixel; i++) {
		regions[i] = region(0, sc.settings.image_width, 0, sc.settings.image_height, 1);
	}
	return regions;
}

// Check if an xy-coordinate is a valid position to place a spiral block
bool is_valid(int x, int y, int w, int h, int size) {
    return (x > -size && x < w) && (y > -size && y < h);
}

// Add a block to a vector of regions (used in spiral generation)
void add_block(std::vector<region>& regions, int x, int y, int size, int samples) {
    regions.push_back(region(x, x + size, y, y + size, samples));
}

// Create a vector of regions with blocks in a spiral pattern
std::vector<region> create_spiral_region(const scene& sc, bool reverse = false) {
    std::vector<region> regions;
    int iw = sc.settings.image_width;
    int ih = sc.settings.image_height;
    int samples = sc.settings.samples_per_pixel;

    // Determine size of blocks and start coordinate
    int block_size = (int)std::ceil(fmin(iw, ih) / 8);
    int x = (int)(iw / 2);
    int y = (int)(ih / 2);

    // Add the center four regions
    add_block(regions, x, y, block_size, samples);
    x -= block_size;
    add_block(regions, x, y, block_size, samples);
    y -= block_size;
    add_block(regions, x, y, block_size, samples);
    x += block_size;
    add_block(regions, x, y, block_size, samples);

    // Try adding blocks in a spiral pattern until image has been covered
    int n = 3;
    bool placed_any = true;
    while (placed_any) {
        placed_any = false;

        // Step 1 right (next spiral layer)
        x += block_size;
        if (is_valid(x, y, iw, ih, block_size)) {
            add_block(regions, x, y, block_size, samples);
            placed_any = true;
        }

        // Step n - 1 down
        for (int i = 0; i < n - 1; i++) {
            y += block_size;
            if (is_valid(x, y, iw, ih, block_size)) {
                add_block(regions, x, y, block_size, samples);
                placed_any = true;
            }
        }

        // Step n left
        for (int i = 0; i < n; i++) {
            x -= block_size;
            if (is_valid(x, y, iw, ih, block_size)) {
                add_block(regions, x, y, block_size, samples);
                placed_any = true;
            }
        }

        // Step n up
        for (int i = 0; i < n; i++) {
            y -= block_size;
            if (is_valid(x, y, iw, ih, block_size)) {
                add_block(regions, x, y, block_size, samples);
                placed_any = true;
            }
        }

        // Step n right
        for (int i = 0; i < n; i++) {
            x += block_size;
            if (is_valid(x, y, iw, ih, block_size)) {
                add_block(regions, x, y, block_size, samples);
                placed_any = true;
            }
        }

        // Increment current spiral arm length by 2
        n += 2;
    }

    if (reverse) {
        std::reverse(regions.begin(), regions.end());
    }
    return regions;
}

/*
    RAY TRACING FUNCTIONS (COLOR, ALBEDO, NORMAL, ...)
*/
// Function that recursively determines the color of a ray shot into a scene with a certain environment and maximum ray depth
color ray_color(const ray& r, const scene* _sc, int depth) {
    if (depth <= 0)
        return color(0, 0, 0);

    ray_count++;

    hit_record rec;
    if (!_sc->world.hit(r, 0.001, infinity, rec))
        return _sc->env.value(r);

    ray scattered;
    color attenuation;
    color emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);

    if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
        return emitted;

    return emitted + attenuation * ray_color(scattered, _sc, depth - 1);
}

// Function that determines the albedo color of whatever object the ray hits first
color ray_albedo(const ray& r, const scene* _sc) {
    ray_count++;
    hit_record rec;
    if (!_sc->world.hit(r, 0.001, infinity, rec)) {
        return _sc->env.value(r);
    }

    ray scattered;
    color attenuation;
    color emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);

    if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
        return emitted;
    }
    return attenuation;
}

// Function that determines the normal of whatever object the ray hits first
vec3 ray_normal(const ray& r, const scene* _sc) {
    hit_record rec;
    ray_count++;
    if (!_sc->world.hit(r, 0.001, infinity, rec)) {
        return vec3(0, 0, 0);
    }

    return unit_vector(rec.normal);
}

// Function that determines the depth (distance from camera) of whatever object the ray hits first
vec3 ray_depth(const ray& r, const scene* _sc) {
    hit_record rec;
    ray_count++;
    if (!_sc->world.hit(r, 0.001, infinity, rec)) {
        return vec3(nan(""), nan(""), nan(""));
    }

    double dist = fabs(rec.t);
    return vec3(dist, dist, dist);
}

// Function that returns a 1-vector if anything is hit, and a 0-vector otherwise
vec3 ray_opacity(const ray& r, const scene* _sc) {
    hit_record rec;
    ray_count++;
    if (!_sc->world.hit(r, 0.001, infinity, rec)) {
        return vec3(0, 0, 0);
    }
    return vec3(1, 1, 1);
}

// Function that returns the color of emissive objects if they are hit directly
vec3 ray_emission(const ray& r, const scene* _sc) {
    hit_record rec;
    ray_count++;
    if (!_sc->world.hit(r, 0.001, infinity, rec)) {
        return vec3(0, 0, 0);
    }

    return rec.mat_ptr->emitted(rec.u, rec.v, rec.p);
}

/*
    WORKER FUNCTION FOR MULTITHREADING
*/
// Worker function that runs on a thread and processes regions of the image
void worker(int id, const scene* sc, std::vector<color>* image_buffer, std::vector<int>* sample_counts, std::vector<region>* q, std::string pass = "color") {
    // Store render settings
    int iw = sc->settings.image_width;
    int ih = sc->settings.image_height;
    int n_samples = sc->settings.samples_per_pixel;
    int depth = sc->settings.max_depth;

    // Grab an item from the queue
    while (done_count < q->size()) {
        int i = done_count++;

        // Get the region that has to be rendered, create a vector for storing the pixel colors in this region
        region r = q->at(i);
        int region_width = r.x_end - r.x_start;
        int region_height = r.y_end - r.y_start;
        std::vector<color> render_area(region_width * region_height, color(0, 0, 0));

        // Render the region
        for (int s = 0; s < r.samples; s++) {
            for (int x = r.x_start; x < r.x_end; x++) {
                for (int y = r.y_start; y < r.y_end; y++) {
                    if (x >= 0 && x < iw && y >= 0 && y < ih) {
                        auto u = ((double)x + random_double()) / (iw - 1.0);
                        auto v = ((double)y + random_double()) / (ih - 1.0);
                        ray from_cam = sc->cam->get_ray(u, v);
                        int idx = region_width * (y - r.y_start) + (x - r.x_start);
                        if (pass == "color") {
                            render_area[idx] = ray_color(from_cam, sc, depth);
                        }
                        else if (pass == "albedo") {
                            render_area[idx] = ray_albedo(from_cam, sc);
                        }
                        else if (pass == "normal") {
                            render_area[idx] = ray_normal(from_cam, sc) * 0.5 + 0.5;
                        }
                        else if (pass == "depth") {
                            render_area[idx] = ray_depth(from_cam, sc);
                        }
                        else if (pass == "opacity") {
                            render_area[idx] = ray_opacity(from_cam, sc);
                        }
                        else if (pass == "emission") {
                            render_area[idx] = ray_emission(from_cam, sc);
                        }
                    }
                }
            }

            {
                // Add the rendered color values to the image buffer
                std::lock_guard<std::mutex> lock(m);
                for (int x = r.x_start; x < r.x_end; x++) {
                    for (int y = r.y_start; y < r.y_end; y++) {
                        if (x >= 0 && x < iw && y >= 0 && y < ih) {
                            int idx = region_width * (y - r.y_start) + (x - r.x_start);
                            image_buffer->at(iw * y + x) += render_area[idx];
                            sample_counts->at(iw * y + x) += 1;
                        }
                    }
                }
            }
        }

        {
            std::lock_guard<std::mutex> lock(m);
            std::cout << "\r   Render progress: " << done_count << "/" << q->size() << " regions (" << (double)done_count / q->size() * 100.0 << "%)            ";
        }
    }

    {
        // Set thread to done before exiting
        std::lock_guard<std::mutex> lock(m);
        threads_done[id] = true;
    }
}

/*
    RENDER UTILITY FUNCTIONS
*/
// Function that starts a number of threads to render the scene
std::vector<std::thread> start_render(const scene& sc, std::vector<color>& image_buffer, std::vector<int>& counts, std::vector<region>& queue, unsigned int n_threads, std::string pass = "color") {
    // Reset tracking variables
    done_count = 0;
    ray_count = 0;
    threads_done = std::vector<bool>(n_threads);
    threads_done.resize(n_threads, false);

    // Create threads using the worker function
    std::vector<std::thread> threads;
    unsigned int new_thread_index = 0;
    while (new_thread_index < n_threads && (int)new_thread_index < sc.settings.image_height) {
        threads.emplace_back(std::thread(worker, new_thread_index, &sc, &image_buffer, &counts, &queue, pass));
        new_thread_index++;
    }

    // Check the actual number of threads that were spawned
    int actual_threads = (int)threads.size();
    if (actual_threads < (int)n_threads) {
        {
            std::lock_guard<std::mutex> lock(m);
            for (int i = 0; i < (int)n_threads - actual_threads; i++) {
                threads_done.pop_back();
            }
        }
    }

    std::cout << current_date_time() << "Successfully spawned " << actual_threads << " threads.\n\n";
    std::cout << "\r   Render progress (lines): 0/" << queue.size() << " regions (0.0%)            ";
    return threads;
}

// Check if all threads have completed rendering
bool check_render_done(std::vector<bool>& done) {
    return std::all_of(done.begin(), done.end(), [](bool v) { return v; });
}

// Join all threads together (waits for all threads to finish executing)
void end_render(std::vector<std::thread>& threads) {
    for (std::thread& thread : threads) {
        if (thread.joinable()) {
            thread.join();
        }
    }
}

// Apply gamma correction and sample scaling to an image buffer
void apply_color_correction(const scene& sc, std::vector<color>& image_buffer, double gamma = 2.0) {
    double scale = 1.0 / double(sc.settings.samples_per_pixel);
    for (int i = 0; i < image_buffer.size(); i++) {
        image_buffer[i][0] = pow(image_buffer[i].x() * scale, 1/gamma);
        image_buffer[i][1] = pow(image_buffer[i].y() * scale, 1/gamma);
        image_buffer[i][2] = pow(image_buffer[i].z() * scale, 1/gamma);
    }
}

// Function that checks if a scene was assigned a filename, if not, it assigns one
void add_filename(scene& sc) {
    if (sc.settings.filename.length() < 1) {
        sc.settings.filename = std::to_string(long(std::time(nullptr)));
    }
}

// Function that saves an image buffer to a .ppm file (specify scene and duration for comments in saved file)
void save_image_to_file(const scene& sc, std::vector<color>& image_buffer, std::string filename_suffix = "", std::string file_format = "jpg", double duration = 0.0) {
    // Parse filename input
    std::string actual_filename = sc.settings.filename;
    if (filename_suffix.length() > 0) {
        actual_filename += "_" + filename_suffix;
    }
    actual_filename += "." + file_format;

    // Create folder for files
    std::string folder_path = "output/" + sc.settings.filename + "/";
    std::string actual_path = folder_path + actual_filename;
    std::filesystem::create_directories(folder_path);

    if (file_format == "ppm") {
        // Open file stream
        std::ofstream file;
        file.open(actual_path);

        // Write image to file (with comments about rendered scene)
        file << "P3\n";
        file << "# Samples: " << sc.settings.samples_per_pixel << ", max. bounces: " << sc.settings.max_depth << ", render time: " << duration << " seconds\n";
        file << sc.settings.image_width << " " << sc.settings.image_height << "\n255\n";
        for (int j = sc.settings.image_height - 1; j >= 0; j--) {
            for (int i = 0; i < sc.settings.image_width; i++) {
                color pixel_color = image_buffer[j * sc.settings.image_width + i];
                write_color(file, pixel_color);
            }
        }

        // Close file
        file.close();
    }
    else if (file_format == "jpg" || file_format == "png") {
        // Write to uint8_t buffer (and flip vertically)
        uint8_t* pixels = new uint8_t[sc.settings.image_width * sc.settings.image_height * 3];
        for (int x = 0; x < sc.settings.image_width; x++) {
            for (int y = 0; y < sc.settings.image_height; y++) {
                for (int k = 0; k < 3; k++) {
                    pixels[(sc.settings.image_width * y + x) * 3 + k] = convert_color_value(image_buffer[sc.settings.image_width * (sc.settings.image_height - 1 - y) + x][k]);
                }
            }
        }

        if (file_format == "jpg") {
            stbi_write_jpg(actual_path.c_str(), sc.settings.image_width, sc.settings.image_height, 3, pixels, 100);
        }
        else {
            stbi_write_png(actual_path.c_str(), sc.settings.image_width, sc.settings.image_height, 3, pixels, sc.settings.image_width * 3);
        }
    }
    
    std::cout << current_date_time() << "Saved to '" << actual_filename << "\n";
}

// Function to remap a depth buffer to values between 0 and 1 (based on 
void rescale_depth_buffer(std::vector<color>& depth_buffer) {
    double min_val = infinity;
    double max_val = -infinity;
    for (int i = 0; i < depth_buffer.size(); i++) {
        if (!isnan(depth_buffer[i][0])) {
            min_val = fmin(min_val, depth_buffer[i].x());
            max_val = fmax(max_val, depth_buffer[i].x());
        }
    }

    for (int i = 0; i < depth_buffer.size(); i++) {
        if (isnan(depth_buffer[i][0])) {
            depth_buffer[i] = color(1.0, 1.0, 1.0);
        }
        else {
            double new_val = map(depth_buffer[i][0], min_val, max_val, 0.0, 1.0);
            depth_buffer[i] = color(new_val, new_val, new_val);
        }
    }
}

#include "gui.h"

/*
    RENDERING FUNCTION
*/
// Function that renders a certain pass (color, albedo, normal) of a scene to an image buffer with a preview window
void render(scene& sc, std::vector<color>& image_buffer, std::string pass = "color", bool preview = true, std::vector<int>& counts = DEFAULT_INT_VECTOR) {
    // Create queue of regions to be rendered
    std::vector<region> queue = create_spiral_region(sc);

    // Reset the counts array
    if (counts.size() != image_buffer.size()) {
        counts = std::vector<int>(image_buffer.size());
    }
    else {
        for (int i = 0; i < counts.size(); i++) {
            counts[i] = 0;
        }
    }

    // Set the gamme for color correction based on the render pass
    double gamma = (pass == "normal" || pass == "opacity") ? 1.0 : 2.0;

    // Add a time-based filename to the scene
    add_filename(sc);

    // Print render settings
    std::cout << "\n" << current_date_time() << "Rendering " << pass << " pass with settings:\n";
    std::cout << "                       - Image size: " << sc.settings.image_width << "x" << sc.settings.image_height << "\n";
    std::cout << "                       - Samples: " << sc.settings.samples_per_pixel << " per pixel\n";
    std::cout << "                       - Max. bounces: " << sc.settings.max_depth << "\n";
    std::cout << "                       - Scene: " << sc.world.objects.size() << " objects\n\n";

    // Determine number of supported concurrent threads
    unsigned int n_threads = std::thread::hardware_concurrency() - 1;
    std::cout << current_date_time() << "Trying to spawn " << n_threads << " threads...\n";

    // Start render by spawning threads
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::thread> threads = start_render(sc, image_buffer, counts, queue, n_threads, pass);

    if (preview) {
        render_preview(threads, threads_done, image_buffer, counts, sc.settings.image_width, sc.settings.image_height, gamma);
    }

    // Wait for all workers to end parsing
    end_render(threads);
    auto end = std::chrono::high_resolution_clock::now();

    // Rescale colors using sample count and optional gamma correction
    if (pass == "depth") {
        rescale_depth_buffer(image_buffer);
    }
    else {
        apply_color_correction(sc, image_buffer, gamma);
    }

    auto duration = (double)std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
    std::cout << "\n\n" << current_date_time() << "Render completed in "<< duration << " seconds, simulated " << ray_count << " rays.\n";
}

#endif // !RENDER_H
