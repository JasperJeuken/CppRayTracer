/*
Functions that use Intel Open Image Denoise to denoise an image (optionally using an albedo and normal pass for better results)
*/
#ifndef DENOISE_H
#define DENOISE_H

#include <vector>
#include <OpenImageDenoise/oidn.h>

#include "scenes.h"
#include "utility.h"

// Function to convert a vector of colors to a float array
float* vector_to_array(const std::vector<color>& vec, bool scale_normal = false) {
	float* arr = new float[vec.size() * 3];
	for (int i = 0; i < vec.size(); i++) {
		arr[i * 3 + 0] = (float)vec[i].x();
		arr[i * 3 + 1] = (float)vec[i].y();
		arr[i * 3 + 2] = (float)vec[i].z();

		if (scale_normal) {
			arr[i * 3 + 0] = arr[i * 3 + 0] * 2.0f - 1.0f;
			arr[i * 3 + 1] = arr[i * 3 + 1] * 2.0f - 1.0f;
			arr[i * 3 + 2] = arr[i * 3 + 2] * 2.0f - 1.0f;
		}
	}
	return arr;
}

// Function to convert a float array to a vector of colors
std::vector<color> array_to_vector(const float* arr, size_t s) {
	std::vector<color> vec(s);
	for (int i = 0; i < vec.size(); i++) {
		vec[i][0] = (double)arr[i * 3 + 0];
		vec[i][1] = (double)arr[i * 3 + 1];
		vec[i][2] = (double)arr[i * 3 + 2];
	}
	return vec;
}

// Function that runs an image (and optional albedo and normal passes) through Open Image Denoise
void denoise(const scene& sc, std::vector<color>& image_output, const std::vector<color>& image_color, const std::vector<color>& image_albedo = std::vector<color>(), const std::vector<color>& image_normal = std::vector<color>()) {
	std::cout << "\n" << current_date_time() << "Denoising image...\n";

	// Convert color vectors to float arrays
	float* color_array = vector_to_array(image_color);
	float* albedo_array = vector_to_array(image_albedo);
	float* normal_array = vector_to_array(image_normal, true);
	float* output_array = vector_to_array(image_output);

	// Create OIDN device
	OIDNDevice device = oidnNewDevice(OIDN_DEVICE_TYPE_DEFAULT);
	oidnCommitDevice(device);

	// Create filter
	OIDNFilter filter = oidnNewFilter(device, "RT");
	oidnSetSharedFilterImage(filter, "color", color_array, OIDN_FORMAT_FLOAT3, sc.settings.image_width, sc.settings.image_height, 0, 0, 0);
	if (image_albedo.size() > 0) oidnSetSharedFilterImage(filter, "albedo", albedo_array, OIDN_FORMAT_FLOAT3, sc.settings.image_width, sc.settings.image_height, 0, 0, 0);
	if (image_normal.size() > 0) oidnSetSharedFilterImage(filter, "normal", normal_array, OIDN_FORMAT_FLOAT3, sc.settings.image_width, sc.settings.image_height, 0, 0, 0);
	oidnSetSharedFilterImage(filter, "output", output_array, OIDN_FORMAT_FLOAT3, sc.settings.image_width, sc.settings.image_height, 0, 0, 0);
	oidnCommitFilter(filter);

	// Run filter and check for errors
	oidnExecuteFilter(filter);
	const char* error;
	if (oidnGetDeviceError(device, &error) != OIDN_ERROR_NONE) {
		std::cout << "Denoise error: " << error << "\n";
	}

	// Write result to output vector
	for (int i = 0; i < image_output.size(); i++) {
		image_output[i][0] = output_array[i * 3 + 0];
		image_output[i][1] = output_array[i * 3 + 1];
		image_output[i][2] = output_array[i * 3 + 2];
	}

	// Cleanup
	delete[] color_array;
	delete[] albedo_array;
	delete[] normal_array;
	delete[] output_array;
	oidnReleaseFilter(filter);
	oidnReleaseDevice(device);

	std::cout << current_date_time() << "Image denoised.\n";
}

#endif // !DENOISE_H
