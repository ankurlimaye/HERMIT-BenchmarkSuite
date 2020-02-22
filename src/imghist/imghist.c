// TODO: Notice stuff, etc.

//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <stdint.h>
#include <stdbool.h>

#define STB_IMAGE_IMPLEMENTATION
#include "../../include/imghist/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../../include/imghist/stb_image_write.h"

const int color_depth = 255;

int main(int argc, char *argv[]) {

  int width, height, bpp, image_size;
  bool inputFile = false;
  unsigned char * img;

  for (int i = 1 ; i < argc ; ++i) {
    if(strcmp(argv[i], "-r") == 0) {
      img = stbi_load("goldhill.bmp", &width, &height, &bpp, 1);
      inputFile = true;
    }
  }

  if(! inputFile) {
    printf("Incorrect input file! Using test input: \"lena.bmp\" \n");
    img = stbi_load("lena.bmp", &width, &height, &bpp, 1);
  }

  // read from the input bmp image file
  image_size = width * height;
  unsigned char* output_image = malloc(sizeof(unsigned char) * image_size);

  // create a histogram for the pixel values
  int* histogram = (int*)calloc(sizeof(int), color_depth);
  for(int i = 0; i < image_size; i++){
    histogram[img[i]]++;
  }

  // finding the normalised values using cumulative mass function
  float* transfer_function = (float*)calloc(sizeof(float), color_depth);
  for(int i = 0; i < color_depth; i++){
    for(int j = 0; j < i+1; j++){
      transfer_function[i] += color_depth*((float)histogram[j])/(image_size);
    }
  }

  for(int i = 0; i < image_size; i++){
    output_image[i] = transfer_function[img[i]];
  }

  stbi_write_bmp("out.bmp", width, height, 1, output_image);
  return 0;
}