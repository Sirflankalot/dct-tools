#include <array>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <utility>

#ifdef __FMA__
#include <immintrin.h>
#endif

template <class T, std::size_t size>
class Matrix {
  private:
	using Data_t = std::array<std::array<T, size>, size>;
	Data_t data;

  public:
	using value_type = Data_t;
	Matrix() {
		for (auto&& i : data) {
			for (auto&& j : i) {
				j = 0;
			}
		}
	}
	Matrix(const Data_t& input) : data(input){};
	Matrix(const T (&input)[size][size]) {
		for (size_t i = 0; i < size; ++i) {
			for (size_t j = 0; j < size; ++j) {
				data[i][j] = input[i][j];
			}
		}
	};

	Data_t& get_data() {
		return data;
	};

	typename Data_t::value_type& operator[](int i) {
		return data[i];
	};

	template <class Ti, size_t sizei>
	friend std::ostream& operator<<(std::ostream& os, const Matrix<Ti, sizei>& d);
};

template <class T, std::size_t size>
std::ostream& operator<<(std::ostream& os, const Matrix<T, size>& d) {
	for (auto& r : d.data) {
		for (auto& c : r) {
			std::cout << c << ' ';
		}
		std::cout << '\n';
	}

	return os;
}

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

using Mat8x8i = Matrix<int, 8>;
using Mat8x8f = Matrix<float, 8>;
using Mat8x8d = Matrix<double, 8>;

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

template <class T, class T2, size_t size>
Matrix<decltype(T{} * T2{}), size> operator*(Matrix<T, size>& left, Matrix<T2, size>& right) {
	Matrix<decltype(T{} * T2{}), size> ret;

	for (size_t i = 0; i < size; ++i) {
		for (size_t j = 0; j < size; ++j) {
			float sum = 0;
			for (size_t k = 0; k < size; ++k) {
				sum += left[i][k] * right[k][j];
			}
			ret[i][j] = sum;
		}
	}

	return ret;
}

#ifdef __FMA__
Matrix<float, 8> operator*(Matrix<float, 8>& left, Matrix<float, 8>& right) {
	Matrix<float, 8> ret;

	std::cout << "Blah\n";

	float* a = &left.get_data()[0][0];
	float* b = &right.get_data()[0][0];
	float* r = &ret.get_data()[0][0];

	__m256 a_line, b_line, r_line;
	for (size_t i = 0; i < 64; i += 8) {
		a_line = _mm256_broadcast_ss(a); // a_line = vec8(a[i])
		b_line = _mm256_loadu_ps(b);
		r_line = _mm256_mul_ps(a_line, b_line);
		for (size_t j = 1; j < 8; ++j) {
			a_line = _mm256_broadcast_ss(a + i + j);
			b_line = _mm256_loadu_ps(b + j * 8);
			r_line = _mm256_fmadd_ps(a_line, b_line, r_line);
		}
		_mm256_storeu_ps(r + i, r_line);
	}

	return ret;
}
#endif

Mat8x8f jpeg_dctii(Mat8x8f& input, Mat8x8f& reg, Mat8x8f& trans) {
	auto x = reg * input;
	auto y = x * trans;

	return y;
}

Mat8x8f jpeg_idctii(Mat8x8f& input, Mat8x8f& reg, Mat8x8f& trans) {
	auto x = trans * input;
	auto y = x * reg;

	return y;
}

int main() {
	auto[regular, transpose] = create_dct_matrices();

	std::cout << "Regular dct matrix:\n" << regular << '\n';

	std::cout << "Transposed dct matrix:\n" << transpose << '\n';

	std::cout << "Input Matrix:\n" << pixels4 << '\n';

	Mat8x8f pix;
	for (size_t i = 0; i < 8; ++i) {
		for (size_t j = 0; j < 8; ++j) {
			pix[i][j] = pixels4[i][j] - 128;
		}
	}

	auto dctii = jpeg_dctii(pix, regular, transpose);
	std::cout << "Output Matrix:\n" << dctii << '\n';

	auto quant = quantizize(dctii);
	std::cout << "Quantizize:\n" << quant << '\n';

	auto dequant = dequantizize(quant);
	std::cout << "Restored:\n" << dequant << '\n';

	auto idctii = jpeg_idctii(dequant, regular, transpose);
	for (size_t i = 0; i < 8; ++i) {
		for (size_t j = 0; j < 8; ++j) {
			idctii[i][j] = std::round(idctii[i][j] + 128);
		}
	}
	std::cout << "Decoded:\n" << idctii << '\n';
}