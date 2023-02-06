/*
Various utility functions and definitions used throughout the program
*/
#ifndef UTILITY_H
#define UTILITY_H

// Common standard headers
#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <chrono>

using std::shared_ptr;
using std::make_shared;

// Constants
const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

// Functions

// Convert degrees value to radians
inline double degrees_to_radians(double degrees) {
	return degrees * pi / 180.0;
}

// Convert radians value to degrees
inline double radians_to_degrees(double radians) {
	return radians * 180.0 / pi;
}

// Generate a random double value between minimum and maximum value (thread-safe)
inline double random_double(double min, double max) {
	thread_local std::mt19937 generator(std::random_device{}());
	std::uniform_real_distribution<double> distribution(min, max);
	return distribution(generator);
}

// Generate a random double value between 0.0 and 1.0 (thread-safe)
inline double random_double() {
	return random_double(0.0, 1.0);
}

// Generate a random integer value between minimum and maximum value (thread-safe)
inline int random_int(int min, int max) {
	return static_cast<int>(random_double(double(min), double(max) + 1));
}

// Clamp a value between a minimum and a maximum value
inline double clamp(double x, double min, double max) {
	if (x < min) return min;
	if (x > max) return max;
	return x;
}

// Map a value from one range to another
inline double map(double value, double min0, double max0, double min1, double max1) {
	return min1 + ((value - min0) * (max1 - min1)) / (max0 - min0);
}

// Return a string with the current date and time (i.e. "[2000-01-01 00:00:00] "
inline std::string current_date_time() {
	time_t now = time(0);
	struct tm tstruct;
	char buf[80];
	localtime_s(&tstruct, &now);
	strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", &tstruct);
	return "[" + std::string(buf) + "] ";
}

template<typename Base, typename T>
inline bool instanceof(const T* ptr) {
	return dynamic_cast<const Base*>(ptr) != nullptr;
}

inline bool ends_with(const std::string& text, const std::string& ending) {
	if (ending.size() > text.size()) return false;
	return std::equal(ending.rbegin(), ending.rend(), text.rbegin());
}

// Common custom headers
#include "vec3.h"
#include "ray.h"

#endif // !UTILITY_H
