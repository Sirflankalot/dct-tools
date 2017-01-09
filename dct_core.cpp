#include "dct_core.hpp"

#include <cmath>

Mat8x8f dct_regular;
Mat8x8f dct_transpose;

Matrix<float, 8> quanttable{
    {
        //
        {16, 11, 10, 16, 24, 40, 51, 61},     //
        {12, 12, 14, 19, 26, 58, 60, 55},     //
        {14, 13, 16, 24, 40, 57, 69, 56},     //
        {14, 17, 22, 29, 51, 87, 80, 62},     //
        {18, 22, 37, 56, 68, 109, 103, 77},   //
        {24, 35, 55, 64, 81, 104, 113, 92},   //
        {49, 64, 78, 87, 103, 121, 120, 101}, //
        {72, 92, 95, 98, 112, 100, 103, 99},  //
    }                                         //
};

Mat8x8i quantizize(Mat8x8f& input) {
	Mat8x8i ret;

	for (size_t i = 0; i < 8; ++i) {
		for (size_t j = 0; j < 8; ++j) {
			ret[i][j] = std::round(input[i][j] / quanttable[i][j]);
		}
	}

	return ret;
}

Mat8x8f dequantizize(Mat8x8i& input) {
	Mat8x8f ret;

	for (size_t i = 0; i < 8; ++i) {
		for (size_t j = 0; j < 8; ++j) {
			ret[i][j] = input[i][j] * quanttable[i][j];
		}
	}

	return ret;
}

std::pair<Mat8x8f, Mat8x8f> create_dct_matrices() {
	Mat8x8f result;
	Mat8x8f transpose;

	auto c = [](float k) -> float {
		if (k == 0) {
			return std::sqrt(0.5);
		}
		else {
			return 1;
		}
	};

	constexpr float N = 8;
	for (size_t k = 0; k < N; ++k) {
		for (size_t n = 0; n < N; ++n) {
			float leftside = std::sqrt(2 / N);
			float middle = c(k);

			float numerator = (2 * M_PI) * (2 * n + 1) * k;
			float denominator = 4 * N;

			float cosine = std::cos(numerator / denominator);
			float rightside = cosine;

			float pixel = leftside * middle * rightside;

			result[k][n] = pixel;
		}
	}

	transpose = ::transpose(result);

	return std::make_pair(result, transpose);
}

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
