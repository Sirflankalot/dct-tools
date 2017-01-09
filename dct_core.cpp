#include "dct_core.hpp"

Mat8x8f dct_regular;
Mat8x8f dct_transpose;

Mat8x8f jpeg_dctii(Mat8x8f& input) {
	auto x = dct_regular * input;
	auto y = x * dct_transpose;

	return y;
}

Mat8x8f jpeg_idctii(Mat8x8f& input) {
	auto x = dct_transpose * input;
	auto y = x * dct_regular;

	return y;
}
