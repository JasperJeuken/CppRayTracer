/*
Functions for creating an SDL window that shows the current image while rendering on multiple threads
*/
#ifndef GUI_H
#define GUI_H


#define SDL_MAIN_HANDLED
#include "SDL.h"
#undef main

#include "render.h"
#include "utility.h"
#include <iostream>
#include <thread>
#include <vector>

// Function that converts a vector of colors to an array of unsigned char (for SDL renderer)
void image_to_texture(const std::vector<color>& image, const std::vector<int>& counts, std::vector<unsigned char>& texture, double gamma) {
	for (int i = 0; i < image.size(); i++) {
		color col = image[i];
		double count = (double)counts[i];

		if (count < 0.000001) {
			texture[i * 4 + 0] = SDL_ALPHA_TRANSPARENT;
			texture[i * 4 + 1] = 0;
			texture[i * 4 + 2] = 0;
			texture[i * 4 + 3] = 0;
		}
		else {
			double r = pow(col[0] / count, 1/gamma);
			double g = pow(col[1] / count, 1/gamma);
			double b = pow(col[2] / count, 1/gamma);

			texture[i * 4 + 0] = SDL_ALPHA_OPAQUE;
			texture[i * 4 + 1] = static_cast<unsigned int>(256 * clamp(r, 0.0, 0.999));  // r
			texture[i * 4 + 2] = static_cast<unsigned int>(256 * clamp(g, 0.0, 0.999));  // g
			texture[i * 4 + 3] = static_cast<unsigned int>(256 * clamp(b, 0.0, 0.999));  // b
		}
	}
}

// Draw a checkered background using SDL renderer
void draw_checker(SDL_Renderer* renderer, int width, int height, int grid_size = 16) {
	for (int y = 0; y < height / grid_size + grid_size; y++) {
		for (int x = 0; x < width / grid_size + grid_size; x++) {
			if (y % 2 == x % 2) {
				SDL_SetRenderDrawColor(renderer, 70, 70, 70, SDL_ALPHA_OPAQUE);
			}
			else {
				SDL_SetRenderDrawColor(renderer, 52, 52, 52, SDL_ALPHA_OPAQUE);
			}
			SDL_Rect checker = { x * grid_size, y * grid_size, grid_size, grid_size };
			SDL_RenderFillRect(renderer, &checker);
		}
	}
}

// Function that creates a window and continuously draws the image buffer
void render_preview(std::vector<std::thread>& threads, std::vector<bool>& threads_done, std::vector<color>& image_buffer, std::vector<int>& counts, int image_width, int image_height, double gamma = 2.0) {
	// Set window variables
	const double ASPECT_RATIO = (double)image_width / (double)image_height;
	const int SCREEN_WIDTH = 800;
	const int SCREEN_HEIGHT = (int)(SCREEN_WIDTH / ASPECT_RATIO);

	// Define window and screen surface
	SDL_Window* window = NULL;	
	SDL_Surface* screen_surace = NULL;

	// Initialise SDL for video, return if error
	if (SDL_Init(SDL_INIT_VIDEO) < 0) {
		std::cout << "SDL could not initialize. Error: " << SDL_GetError() << "\n";
		return;
	}

	// Create window, return if error
	window = SDL_CreateWindow("Raytracing Preview", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
	if (window == NULL) {
		std::cout << "Window could not be created. Error: " << SDL_GetError() << "\n";
		return;
	}

	// Create renderer and texture
	SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
	SDL_Texture* texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_BGRA8888, SDL_TEXTUREACCESS_STREAMING, image_width, image_height);
	
	// Create texture buffer, initialized as fully transparent
	std::vector<unsigned char> texture_buffer((int)image_width * (int)image_height * 4, 0);
	for (int y = 0; y < image_height; y++) {
		for (int x = 0; x < image_width; x++) {
			const unsigned int idx = (image_width * 4 * y) + x * 4;
			texture_buffer[idx + 0] = SDL_ALPHA_TRANSPARENT;
			texture_buffer[idx + 1] = 0;
			texture_buffer[idx + 2] = 0;
			texture_buffer[idx + 3] = 0;
		}
	}

	// GUI loop
	SDL_Event e;
	bool window_open = true;
	while (window_open) {

		// Reset the renderer and draw a checker background
		SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_TRANSPARENT);
		SDL_RenderClear(renderer);
		draw_checker(renderer, SCREEN_WIDTH, SCREEN_HEIGHT);
	
		// Close the window if the user clicks the close button
		while (SDL_PollEvent(&e)) {
			if (e.type == SDL_QUIT) {
				window_open = false;
			}
		}

		// Close the window if the render was completed
		if (check_render_done(threads_done)) {
			window_open = false;
		}

		// Update texture buffer based on current image and sample counts
		image_to_texture(image_buffer, counts, texture_buffer, gamma);

		// Copy texture buffer memory to texture using locking
		unsigned char* locked_pixels = nullptr;
		int pitch = 0;
		SDL_LockTexture(texture, NULL, reinterpret_cast<void**>(&locked_pixels), &pitch);
		SDL_memcpy(locked_pixels, texture_buffer.data(), texture_buffer.size());
		SDL_UnlockTexture(texture);

		// Render texture on screen
		SDL_SetTextureBlendMode(texture, SDL_BLENDMODE_BLEND);
		SDL_RenderCopyEx(renderer, texture, NULL, NULL, 0, NULL, SDL_FLIP_VERTICAL);
		SDL_RenderPresent(renderer);
		SDL_Delay(50);

	}

	// Destroy window and quit SDL
	SDL_DestroyTexture(texture);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();
}

#endif // !GUI_H
