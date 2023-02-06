/*
Custom loader of the STB_IMAGE library for loading image files (disables warnings)
*/
#ifndef RT_STB_IMAGE_H
#define RT_STB_IMAGE_H

// Disable warnings
#ifdef _MSC_VER
#pragma warning (push, 0)
#endif

// Include image loading library
#define STB_IMAGE_IMPLEMENTATION
#include "external/stb_image.h"

// Include image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "external/stb_image_write.h"

// Restore warnings
#ifdef _MSC_VER
#pragma warning (pop)
#endif

#endif // !RT_STB_IMAGE_H
