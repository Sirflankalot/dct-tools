#pragma once

#include <array>
#include <iostream>

#ifdef __AVX__
#include <immintrin.h>
#elif defined __SSE__
#include <xmmintrin.h>
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
	}

	Data_t& get_data() {
		return data;
	}

	typename Data_t::value_type& operator[](int i) {
		return data[i];
	}

	template <class Ti, size_t sizei>
	friend std::ostream& operator<<(std::ostream& os, const Matrix<Ti, sizei>& d);
};

using Mat8x8i = Matrix<int, 8>;
using Mat8x8f = Matrix<float, 8>;
using Mat8x8d = Matrix<double, 8>;

template <class T, std::size_t size>
inline std::ostream& operator<<(std::ostream& os, const Matrix<T, size>& d) {
	for (auto& r : d.data) {
		for (auto& c : r) {
			std::cout << c << ' ';
		}
		std::cout << '\n';
	}

	return os;
}
template <class T, class T2, size_t size>
inline Matrix<decltype(T{} * T2{}), size> operator*(Matrix<T, size>& left, Matrix<T2, size>& right) {
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

#if defined __AVX__ && defined SIMD
__attribute__((always_inline)) inline Matrix<float, 8> operator*(Matrix<float, 8>& left, Matrix<float, 8>& right) {
	Matrix<float, 8> ret;

	float* a = &left[0][0];
	float* b = &right[0][0];
	float* r = &ret[0][0];

	__m256 row1 = _mm256_loadu_ps(b + 0);
	__m256 row2 = _mm256_loadu_ps(b + 8);
	__m256 row3 = _mm256_loadu_ps(b + 16);
	__m256 row4 = _mm256_loadu_ps(b + 24);
	__m256 row5 = _mm256_loadu_ps(b + 32);
	__m256 row6 = _mm256_loadu_ps(b + 40);
	__m256 row7 = _mm256_loadu_ps(b + 48);
	__m256 row8 = _mm256_loadu_ps(b + 56);
	for (size_t i = 0; i < 8; ++i) {
		__m256 brod1 = _mm256_broadcast_ss(a + (8 * i + 0));
		__m256 brod2 = _mm256_broadcast_ss(a + (8 * i + 1));
		__m256 brod3 = _mm256_broadcast_ss(a + (8 * i + 2));
		__m256 brod4 = _mm256_broadcast_ss(a + (8 * i + 3));
		__m256 brod5 = _mm256_broadcast_ss(a + (8 * i + 4));
		__m256 brod6 = _mm256_broadcast_ss(a + (8 * i + 5));
		__m256 brod7 = _mm256_broadcast_ss(a + (8 * i + 6));
		__m256 brod8 = _mm256_broadcast_ss(a + (8 * i + 7));

		__m256 row = _mm256_mul_ps(brod1, row1);
#ifdef __FMA__
		row = _mm256_fmadd_ps(brod2, row2, row);
		row = _mm256_fmadd_ps(brod3, row3, row);
		row = _mm256_fmadd_ps(brod4, row4, row);
		row = _mm256_fmadd_ps(brod5, row5, row);
		row = _mm256_fmadd_ps(brod6, row6, row);
		row = _mm256_fmadd_ps(brod7, row7, row);
		row = _mm256_fmadd_ps(brod8, row8, row);
#else
		row = _mm256_add_ps(_mm256_mul_ps(brod2, row2), row);
		row = _mm256_add_ps(_mm256_mul_ps(brod3, row3), row);
		row = _mm256_add_ps(_mm256_mul_ps(brod4, row4), row);
		row = _mm256_add_ps(_mm256_mul_ps(brod5, row5), row);
		row = _mm256_add_ps(_mm256_mul_ps(brod6, row6), row);
		row = _mm256_add_ps(_mm256_mul_ps(brod7, row7), row);
		row = _mm256_add_ps(_mm256_mul_ps(brod8, row8), row);
#endif

		_mm256_storeu_ps(r + (8 * i), row);
	}

	return ret;
}
#elif defined __SSE__ && defined SIMD
__attribute__((always_inline)) inline Matrix<float, 8> operator*(Matrix<float, 8>& left, Matrix<float, 8>& right) {
	Matrix<float, 8> ret;

	float* a = &left[0][0];
	float* b = &right[0][0];
	float* r = &ret[0][0];

	__m128 row1a = _mm_loadu_ps(b + 0);
	__m128 row1b = _mm_loadu_ps(b + 4);
	__m128 row2a = _mm_loadu_ps(b + 8);
	__m128 row2b = _mm_loadu_ps(b + 12);
	__m128 row3a = _mm_loadu_ps(b + 16);
	__m128 row3b = _mm_loadu_ps(b + 20);
	__m128 row4a = _mm_loadu_ps(b + 24);
	__m128 row4b = _mm_loadu_ps(b + 28);

	for (size_t i = 0; i < 8; ++i) {
		__m128 brod1 = _mm_load_ps1(a + (8 * i + 0));
		__m128 brod2 = _mm_load_ps1(a + (8 * i + 1));
		__m128 brod3 = _mm_load_ps1(a + (8 * i + 2));
		__m128 brod4 = _mm_load_ps1(a + (8 * i + 3));

		__m128 rowa = _mm_mul_ps(row1a, brod1);
		__m128 rowb = _mm_mul_ps(row1b, brod1);
		rowa = _mm_add_ps(_mm_mul_ps(row2a, brod2), rowa);
		rowb = _mm_add_ps(_mm_mul_ps(row2b, brod2), rowb);
		rowa = _mm_add_ps(_mm_mul_ps(row3a, brod3), rowa);
		rowb = _mm_add_ps(_mm_mul_ps(row3b, brod3), rowb);
		rowa = _mm_add_ps(_mm_mul_ps(row4a, brod4), rowa);
		rowb = _mm_add_ps(_mm_mul_ps(row4b, brod4), rowb);

		_mm_storeu_ps(r + (8 * i + 0), rowa);
		_mm_storeu_ps(r + (8 * i + 4), rowb);
	}

	__m128 row5a = _mm_loadu_ps(b + 32);
	__m128 row5b = _mm_loadu_ps(b + 36);
	__m128 row6a = _mm_loadu_ps(b + 40);
	__m128 row6b = _mm_loadu_ps(b + 44);
	__m128 row7a = _mm_loadu_ps(b + 48);
	__m128 row7b = _mm_loadu_ps(b + 52);
	__m128 row8a = _mm_loadu_ps(b + 56);
	__m128 row8b = _mm_loadu_ps(b + 60);

	for (size_t i = 0; i < 8; ++i) {
		__m128 brod5 = _mm_load_ps1(a + (8 * i + 4));
		__m128 brod6 = _mm_load_ps1(a + (8 * i + 5));
		__m128 brod7 = _mm_load_ps1(a + (8 * i + 6));
		__m128 brod8 = _mm_load_ps1(a + (8 * i + 7));

		__m128 rowa = _mm_loadu_ps(r + (8 * i + 0));
		__m128 rowb = _mm_loadu_ps(r + (8 * i + 4));

		rowa = _mm_add_ps(_mm_mul_ps(row5a, brod5), rowa);
		rowb = _mm_add_ps(_mm_mul_ps(row5b, brod5), rowb);
		rowa = _mm_add_ps(_mm_mul_ps(row6a, brod6), rowa);
		rowb = _mm_add_ps(_mm_mul_ps(row6b, brod6), rowb);
		rowa = _mm_add_ps(_mm_mul_ps(row7a, brod7), rowa);
		rowb = _mm_add_ps(_mm_mul_ps(row7b, brod7), rowb);
		rowa = _mm_add_ps(_mm_mul_ps(row8a, brod8), rowa);
		rowb = _mm_add_ps(_mm_mul_ps(row8b, brod8), rowb);

		_mm_storeu_ps(r + (8 * i + 0), rowa);
		_mm_storeu_ps(r + (8 * i + 4), rowb);
	}

	return ret;
}
#endif

template <class T, size_t size>
inline auto transpose(Matrix<T, size>& val) {
	Matrix<T, size> transpose;
	for (size_t i = 0; i < size; ++i) {
		for (size_t j = 0; j < size; ++j) {
			transpose[j][i] = val[i][j];
		}
	}

	return transpose;
}
