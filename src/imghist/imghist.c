// TODO: Notice stuff, etc.

#include <stdbool.h>

#define STB_IMAGE_IMPLEMENTATION
#include "../../include/imghist/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../../include/imghist/stb_image_write.h"

const int maxPixelValue = 255;

unsigned char *histogramEqualization(unsigned char *input, int width, int height) {
  int inputSize = width * height;
  unsigned char *output = malloc(sizeof(unsigned char) * inputSize);

  // Create a histogram for the pixel values
  int *histogram = (int *) calloc(sizeof(int), maxPixelValue);
  for (int i = 0 ; i < inputSize ; i++) {
    histogram[input[i]]++;
  }

  // Calculate normalized pixel values using CMF
  float *normalizationFunction = (float *) calloc(sizeof(float), maxPixelValue);
  for (int i = 0 ; i < maxPixelValue ; i++) {
    for (int j = 0 ; j < i + 1 ; j++) {
      normalizationFunction[i] += maxPixelValue * ((float) histogram[j]) / (inputSize);
    }
  }

  for (int i = 0 ; i < inputSize ; i++) {
    output[i] = normalizationFunction[input[i]];
  }

  return output;
}

int main(int argc, char *argv[]) {

  int width, height, bpp;
  bool inputFile = false;
  unsigned char *img;

  for (int i = 1 ; i < argc ; ++i) {
    if (strcmp(argv[i], "-r") == 0) {
      img = stbi_load(argv[(i + 1)], &width, &height, &bpp, 1);
      inputFile = true;
    }
  }

  if (!inputFile) {
    printf("Incorrect input file! Using test input: \"lena.bmp\" \n");
    img = stbi_load("lena.bmp", &width, &height, &bpp, 1);
  }

  unsigned char *outputFile = histogramEqualization(img, width, height);

  stbi_write_bmp("out.bmp", width, height, 1, outputFile);
  return 0;
}