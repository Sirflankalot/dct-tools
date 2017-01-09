#include <array>
#include <chrono>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

#include "dct_core.hpp"

template <class T, size_t size>
auto transpose(Matrix<T, size>& val) {
	Matrix<T, size> transpose;
	for (size_t i = 0; i < size; ++i) {
		for (size_t j = 0; j < size; ++j) {
			transpose[j][i] = val[i][j];
		}
	}

	return transpose;
}

Matrix<float, 8> pixels{{
    //
    {140, 144, 147, 140, 140, 155, 179, 175}, //
    {144, 152, 140, 147, 140, 148, 167, 179}, //
    {152, 155, 136, 167, 163, 162, 152, 172}, //
    {168, 145, 156, 160, 152, 155, 136, 160}, //
    {162, 148, 156, 148, 140, 136, 147, 162}, //
    {147, 167, 140, 155, 155, 140, 136, 162}, //
    {136, 156, 123, 167, 162, 144, 140, 147}, //
    {148, 155, 136, 155, 152, 147, 147, 136}
    //
}};

Matrix<float, 8> pixels2{
    {
        //
        {100, 100, 100, 100, 100, 100, 100, 100}, //
        {100, 100, 100, 100, 100, 100, 100, 100}, //
        {100, 100, 100, 100, 100, 100, 100, 100}, //
        {100, 100, 100, 100, 100, 100, 100, 100}, //
        {100, 100, 100, 100, 100, 100, 100, 100}, //
        {100, 100, 100, 100, 100, 100, 100, 100}, //
        {100, 100, 100, 100, 100, 100, 100, 100}, //
        {100, 100, 100, 100, 100, 100, 100, 100}, //
    }                                             //
};

Matrix<float, 8> pixels3{
    {
        //
        {52, 55, 61, 66, 71, 61, 64, 73},    //
        {63, 59, 66, 90, 109, 585, 69, 72},  //
        {62, 59, 68, 113, 144, 104, 66, 73}, //
        {63, 58, 71, 122, 154, 106, 70, 69}, //
        {67, 61, 68, 104, 126, 88, 68, 70},  //
        {79, 65, 60, 70, 77, 68, 58, 75},    //
        {85, 71, 64, 59, 55, 61, 65, 83},    //
        {87, 79, 69, 58, 65, 76, 78, 94},    //
    }                                        //
};

Matrix<float, 8> pixels4{
    {
        //
        {154, 123, 123, 123, 123, 123, 123, 136}, //
        {192, 180, 136, 154, 154, 154, 136, 110}, //
        {254, 198, 154, 154, 180, 154, 123, 123}, //
        {239, 180, 136, 180, 180, 166, 123, 123}, //
        {180, 154, 136, 167, 166, 149, 136, 136}, //
        {128, 136, 123, 136, 154, 180, 198, 154}, //
        {123, 105, 110, 149, 136, 136, 180, 166}, //
        {110, 136, 123, 123, 123, 136, 154, 136}, //
    }                                             //
};

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

std::array<float, 8> dctii_1d(std::array<float, 8>& input) {
	std::array<float, 8> result; // X

	auto u = [](float in) -> float {
		if (in == 0) {
			return 1.0 / sqrt(2);
		}
		else {
			return 1;
		}
	};

	size_t N = 8;
	for (size_t m = 0; m < N; ++m) {
		float term1 = u(m);
		// float term2 = std::sqrt(2.0 / N);
		float term2 = 1;

		float term3 = 0;
		for (size_t i = 0; i < N; ++i) {
			float leftterm = input[i];

			float numerator = (2 * i + 1) * m * M_PI;
			float denominator = 2 * N;

			float rightterm = std::cos(numerator / denominator);

			term3 += leftterm * rightterm;
		}

		result[m] = term1 * term2 * term3;
	}

	return result;
}

std::array<float, 8> idctii_1d(std::array<float, 8>& input) {
	std::array<float, 8> result; // X

	auto u = [](float in) -> float {
		if (in == 0) {
			return 1.0 / sqrt(2);
		}
		else {
			return 1;
		}
	};

	size_t N = 8;
	for (size_t m = 0; m < N; ++m) {
		float term1 = 2.0 / N;

		float term2 = 0;
		for (size_t i = 0; i < N; ++i) {
			float leftterm = u(i);
			float middleterm = input[i];

			float numerator = (2 * m + 1) * i * M_PI;
			float denominator = 2 * N;

			float rightterm = std::cos(numerator / denominator);

			term2 += leftterm * middleterm * rightterm;
		}

		result[m] = term1 * term2;
	}

	return result;
}

Mat8x8f dctii_1d(Mat8x8f& input) {
	Mat8x8f result;

	for (size_t row = 0; row < 8; ++row) {
		result[row] = dctii_1d(input[row]);
	}

	return result;
}

Mat8x8f idctii_1d(Mat8x8f& input) {
	Mat8x8f result;

	for (size_t row = 0; row < 8; ++row) {
		result[row] = idctii_1d(input[row]);
	}

	return result;
}

Mat8x8f dctii(Mat8x8f& input) {
	Mat8x8f temp, result;

	temp = dctii_1d(input);
	temp = transpose(temp);
	result = dctii_1d(temp);

	// temp = idctii_1d(result);
	// temp = transpose(temp);
	// result = idctii_1d(temp);

	// temp = dctii_1d(input);
	// return idctii_1d(temp);

	return result;
}

bool measure_time(size_t iterations) {
	std::vector<Mat8x8f> input(iterations);
	std::vector<Mat8x8f> output1(iterations);
	std::vector<Mat8x8f> output2(iterations);

	for (size_t i = 0; i < iterations; ++i) {
		input[i] = pixels4;
		input[i][(i % 64) / 8][i % 8] += i / 64;
	}

	auto start = std::chrono::high_resolution_clock::now();

	for (size_t i = 0; i < iterations; ++i) {
		output1[i] = jpeg_dctii(input[i]);
	}

	auto mid = std::chrono::high_resolution_clock::now();

	for (size_t i = 0; i < iterations; ++i) {
		output2[i] = jpeg_idctii(output1[i]);
	}

	auto end = std::chrono::high_resolution_clock::now();

	auto forward = std::chrono::duration_cast<std::chrono::nanoseconds>(mid - start);
	auto backward = std::chrono::duration_cast<std::chrono::nanoseconds>(end - mid);

	float for_time = std::chrono::duration_cast<std::chrono::milliseconds>(forward).count();
	float for_per_element = float(forward.count()) / iterations;

	float back_time = std::chrono::duration_cast<std::chrono::milliseconds>(backward).count();
	float back_per_element = float(backward.count()) / iterations;

	std::cout << "DCT-II: " << for_time << " ms for " << iterations << " iterations\n";
	std::cout << for_per_element << " per iteration\n\n";

	std::cout << "IDCT-II: " << back_time << " ms for " << iterations << " iterations\n";
	std::cout << back_per_element << " per iteration\n\n";

	bool valid = true;
	float err = 0.0f;
	for (size_t i = 0; i < iterations; ++i) {
		for (size_t j = 0; j < 8; ++j) {
			for (size_t k = 0; k < 8; ++k) {
				float diff = std::abs(input[i][j][k] - output2[i][j][k]);
				err = err < diff ? diff : err;
				valid &= diff <= 0.01f;
			}
		}
	}

	std::cout << "Valid " << std::boolalpha << valid << " Error: " << err << '\n';

	return valid;
}

int main(int argc, char**) {
	std::tie(dct_regular, dct_transpose) = create_dct_matrices();

	if (argc == 2) {
		return measure_time(1'000'000);
	}

	std::cout << "Regular dct matrix:\n" << dct_regular << '\n';

	std::cout << "Transposed dct matrix:\n" << dct_transpose << '\n';

	std::cout << "Input Matrix:\n" << pixels4 << '\n';

	Mat8x8f pix;
	for (size_t i = 0; i < 8; ++i) {
		for (size_t j = 0; j < 8; ++j) {
			pix[i][j] = pixels4[i][j] - 128;
		}
	}

	auto dctii = jpeg_dctii(pix);
	std::cout << "Output Matrix:\n" << dctii << '\n';

	auto quant = quantizize(dctii);
	std::cout << "Quantizize:\n" << quant << '\n';

	auto dequant = dequantizize(quant);
	std::cout << "Restored:\n" << dequant << '\n';

	auto idctii = jpeg_idctii(dequant);
	for (size_t i = 0; i < 8; ++i) {
		for (size_t j = 0; j < 8; ++j) {
			idctii[i][j] = std::round(idctii[i][j] + 128);
		}
	}
	std::cout << "Decoded:\n" << idctii << '\n';
}