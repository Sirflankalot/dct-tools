#pragma once

#include "matrix.hpp"
#include <utility>

extern Mat8x8f dct_regular;
extern Mat8x8f dct_transpose;
extern Mat8x8f quanttable;

Mat8x8f jpeg_dctii(Mat8x8f& input);
Mat8x8f jpeg_idctii(Mat8x8f& input);
Mat8x8i quantizize(Mat8x8f& input);
Mat8x8f dequantizize(Mat8x8i& input);
std::pair<Mat8x8f, Mat8x8f> create_dct_matrices();
