#pragma once

#include "matrix.hpp"

extern Mat8x8f dct_regular;
extern Mat8x8f dct_transpose;

Mat8x8f jpeg_dctii(Mat8x8f& input);
Mat8x8f jpeg_idctii(Mat8x8f& input);
