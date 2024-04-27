#pragma once

#include <initializer_list>
#include <cmath>
#include <cstdint>

template <typename T, uint64_t D>
struct Vector
{
	static constexpr uint64_t size = D;
	using type = T;
	T v[D];
	constexpr Vector();
	constexpr Vector(const T& x);
	template <typename... A>
	constexpr Vector(const T& x, const A&... r);
	constexpr T& operator[](uint64_t i);
	constexpr const T& operator[](uint64_t i) const;

	constexpr Vector<T, D> operator+(const Vector<T, D>& b) const;
	constexpr Vector<T, D> operator-(const Vector<T, D>& b) const;
	constexpr Vector<T, D> operator*(const Vector<T, D>& b) const;
	constexpr Vector<T, D> operator/(const Vector<T, D>& b) const;
	constexpr Vector<T, D>& operator+=(const Vector<T, D>& b);
	constexpr Vector<T, D>& operator-=(const Vector<T, D>& b);
	constexpr Vector<T, D>& operator*=(const Vector<T, D>& b);
	constexpr Vector<T, D>& operator/=(const Vector<T, D>& b);
	constexpr Vector<T, D> operator~() const;
};

namespace vecm
{
	template <typename T, uint64_t D>
	constexpr T length(const Vector<T, D>& v);
	template <typename T, uint64_t D>
	constexpr T length2(const Vector<T, D>& v);
	template <typename T, uint64_t D>
	constexpr T rdist(const Vector<T, D>& v);
	template <typename T, uint64_t D>
	constexpr T dot(const Vector<T, D>& a, const Vector<T, D>& b);
	template <typename T, uint64_t D>
	constexpr Vector<T, D> min(const Vector<T, D>& a, const Vector<T, D>& b);
	template <typename T, uint64_t D>
	constexpr T min(const Vector<T, D>& a);
	template <typename T, uint64_t D>
	constexpr Vector<T, D> max(const Vector<T, D>& a, const Vector<T, D>& b);
	template <typename T, uint64_t D>
	constexpr T max(const Vector<T, D>& a);
	template <typename T, uint64_t D>
	constexpr Vector<T, 2> cat(const T& a, const T& b);
	template <typename T, uint64_t D>
	constexpr Vector<T, D+1> cat(const Vector<T, D>& a, const T& b);
	template <typename T, uint64_t D>
	constexpr Vector<T, D+1> cat(const T& a, const Vector<T, D>& b);
	template <typename T, uint64_t D1, uint64_t D2>
	constexpr Vector<T, D1+D2> cat(const Vector<T, D1>& a, const Vector<T, D2>& b);
	template <typename V1, typename... V, std::enable_if_t<std::greater()(sizeof...(V), 1), bool> = true>
	constexpr auto cat(const V1& a, const V&... b);
	template <uint64_t... I, typename T, uint64_t D>
	constexpr Vector<T, sizeof...(I)> sw(const Vector<T, D>& a);
	template <typename T, uint64_t M, uint64_t D>
	constexpr Vector<T, M/D> mul_mat_vec(const Vector<T, M>& m, const Vector<T, D>& v);
	template <uint64_t K, typename T, uint64_t I, uint64_t J>
	constexpr Vector<T, K> mul_mat_mat(const Vector<T, I>& m1, const Vector<T, J>& m2);
	template <typename T, uint64_t D>
	constexpr Vector<T, D> inv_mat(const Vector<T, D>& m);
};

template <typename T, uint64_t D>
constexpr Vector<T, D>::Vector()
{
	 for (uint64_t i = 0; i < D; i++) v[i] = T();
}

template <typename T, uint64_t D>
constexpr Vector<T, D>::Vector(const T& x)
{
	 for (uint64_t i = 0; i < D; i++) v[i] = x;
}

template <typename T, uint64_t D>
template <typename... A>
constexpr Vector<T, D>::Vector(const T& x, const A&... r) : v{x, static_cast<T>(r)...}
{
	static_assert(D == sizeof...(r) + 1);
}

template <typename T, uint64_t D>
constexpr T& Vector<T, D>::operator[](const uint64_t i)
{
	return v[i];
}

template <typename T, uint64_t D>
constexpr const T& Vector<T, D>::operator[](const uint64_t i) const
{
	return v[i];
}

template <typename T, uint64_t D>
constexpr Vector<T, D> Vector<T, D>::operator+(const Vector<T, D>& b) const
{
	Vector<T, D> out;
	for (uint64_t i = 0; i < D; i++) out.v[i] = v[i] + b.v[i];
	return out;
}

template <typename T, uint64_t D>
constexpr Vector<T, D> Vector<T, D>::operator-(const Vector<T, D>& b) const
{
	Vector<T, D> out;
	for (uint64_t i = 0; i < D; i++) out.v[i] = v[i] - b.v[i];
	return out;
}

template <typename T, uint64_t D>
constexpr Vector<T, D> Vector<T, D>::operator*(const Vector<T, D>& b) const
{
	Vector<T, D> out;
	for (uint64_t i = 0; i < D; i++) out.v[i] = v[i] * b.v[i];
	return out;
}

template <typename T, uint64_t D>
constexpr Vector<T, D> Vector<T, D>::operator/(const Vector<T, D>& b) const
{
	Vector<T, D> out;
	for (uint64_t i = 0; i < D; i++) out.v[i] = v[i] / b.v[i];
	return out;
}

template <typename T, uint64_t D>
constexpr Vector<T, D>& Vector<T, D>::operator+=(const Vector<T, D>& b)
{
	for (uint64_t i = 0; i < D; i++) v[i] += b.v[i];
	return *this;
}

template <typename T, uint64_t D>
constexpr Vector<T, D>& Vector<T, D>::operator-=(const Vector<T, D>& b)
{
	for (uint64_t i = 0; i < D; i++) v[i] -= b.v[i];
	return *this;
}

template <typename T, uint64_t D>
constexpr Vector<T, D>& Vector<T, D>::operator*=(const Vector<T, D>& b)
{
	for (uint64_t i = 0; i < D; i++) v[i] *= b.v[i];
	return *this;
}

template <typename T, uint64_t D>
constexpr Vector<T, D>& Vector<T, D>::operator/=(const Vector<T, D>& b)
{
	for (uint64_t i = 0; i < D; i++) v[i] /= b.v[i];
	return *this;
}

template <typename T, uint64_t D>
constexpr Vector<T, D> Vector<T, D>::operator~() const
{
	return *this * vecm::rdist(*this);
}


namespace vecm
{
	template <typename T, uint64_t D>
	constexpr T length(const Vector<T, D>& v)
	{
		return std::sqrt(dot(v, v));
	}
	template <typename T, uint64_t D>
	constexpr T length2(const Vector<T, D>& v)
	{
		return dot(v, v);
	}
	template <typename T, uint64_t D>
	constexpr T rdist(const Vector<T, D>& v)
	{
		return 1.0 / std::sqrt(dot(v, v));
	}
	template <typename T, uint64_t D>
	constexpr T dot(const Vector<T, D>& a, const Vector<T, D>& b)
	{
		T out{};
		for (uint64_t i = 0; i < D; i++) out += a.v[i] * b.v[i];
		return out;
	}
	template <typename T>
	constexpr Vector<T, 3> cross(const Vector<T, 3>& a, const Vector<T, 3>& b)
	{
		return sw<1,2,0>(a) * sw<2,0,1>(b) - sw<1,2,0>(b) * sw<2,0,1>(a);
	}
	template <typename T, uint64_t D>
	constexpr Vector<T, D> min(const Vector<T, D>& a, const Vector<T, D>& b)
	{
		Vector<T, D> out;
		for (uint64_t i = 0; i < D; i++) out.v[i] = std::min(a.v[i], b.v[i]);
		return out;
	}
	template <typename T, uint64_t D>
	constexpr T min(const Vector<T, D>& a)
	{
		T out = a.v[0];
		for (uint64_t i = 1; i < D; i++) out = std::min(out, a.v[i]);
		return out;
	}
	template <typename T, uint64_t D>
	constexpr Vector<T, D> max(const Vector<T, D>& a, const Vector<T, D>& b)
	{
		Vector<T, D> out;
		for (uint64_t i = 0; i < D; i++) out.v[i] = std::max(a.v[i], b.v[i]);
		return out;
	}
	template <typename T, uint64_t D>
	constexpr T max(const Vector<T, D>& a)
	{
		T out = a.v[0];
		for (uint64_t i = 1; i < D; i++) out = std::max(out, a.v[i]);
		return out;
	}
	template <typename T> requires (!requires (T a) { a.v[0]; })
	constexpr Vector<T, 2> cat(const T& a, const T& b)
	{
		return Vector<T, 2>(a, b);
	}
	template <typename T, uint64_t D>
	constexpr Vector<T, D+1> cat(const Vector<T, D>& a, const T& b)
	{
		Vector<T, D+1> out;
		for (uint64_t i = 0; i < D; i++) out.v[i] = a.v[i];
		out.v[D] = b;
		return out;
	}
	template <typename T, uint64_t D>
	constexpr Vector<T, D+1> cat(const T& a, const Vector<T, D>& b)
	{
		Vector<T, D+1> out;
		out.v[0] = a;
		for (uint64_t i = 0; i < D; i++) out.v[i+1] = a.v[i];
		return out;
	}
	template <typename T, uint64_t D1, uint64_t D2>
	constexpr Vector<T, D1+D2> cat(const Vector<T, D1>& a, const Vector<T, D2>& b)
	{
		Vector<T, D1+D2> out;
		for (uint64_t i = 0; i < D1; i++) out.v[i] = a.v[i];
		for (uint64_t i = 0; i < D2; i++) out.v[i+D1] = b.v[i];
		return out;
	}
	template <typename V1, typename... V, std::enable_if_t<std::greater()(sizeof...(V), 1), bool> = true>
	constexpr auto cat(const V1& a, const V&... b)
	{
		return cat(a, cat(b...));
	}
	template <uint64_t... I, typename T, uint64_t D>
	constexpr Vector<T, sizeof...(I)> sw(const Vector<T, D>& a)
	{
		return Vector<T, sizeof...(I)>(a.v[I]...);
	}
	template <typename T, uint64_t M, uint64_t D>
	constexpr Vector<T, M/D> mul_mat_vec(const Vector<T, M>& m, const Vector<T, D>& v)
	{
		static_assert(M % D == 0);
		Vector<T, M/D> out;
		uint64_t k = 0;
		for (uint64_t i = 0; i < M/D; i++)
			for (uint64_t j = 0; j < D; j++)
				out.v[i] += m.v[k++] * v.v[j];
		return out;
	}
	template <uint64_t K, typename T, uint64_t I, uint64_t J>
	constexpr Vector<T, K> mul_mat_mat(const Vector<T, I>& m1, const Vector<T, J>& m2)
	{
		static_assert(J * K % I == 0);
		constexpr uint64_t Z2 = J * K / I;
		constexpr uint64_t Z = std::sqrt(Z2);
		static_assert(Z * Z == Z2);
		static_assert(I % Z == 0);
		static_assert(J % Z == 0);
		constexpr uint64_t X = J / Z;
		constexpr uint64_t Y = I / Z;

		Vector<T, K> out;
		for (uint64_t i = 0; i < Y; i++)
			for (uint64_t j = 0; j < X; j++)
				for (uint64_t k = 0; k < Z; k++)
					out[i*X+j] += m1[i*Z+k] * m2[k*X+j];
		return out;
	}
	template <typename T, uint64_t D>
	constexpr Vector<T, D> inv_mat(const Vector<T, D>& m)
	{
		constexpr uint64_t S = std::sqrt(D);
		static_assert(S*S == D);
		Vector<T, D> mat = m, out;
		for (uint64_t i = 0; i < D; i += S+1) out.v[i] = 1;
		for (uint64_t i = 0; i < S; i++)
		{
			uint64_t imax = i;
			T max = mat.v[i*S+i];
			if (max < 0) max = -max;
			for (uint64_t j = i+1; j < S; j++)
			{
				T t = mat.v[j*S+i];
				if (t < 0) t = -t;
				if (max < t) imax = j, max = t;
			}
			for (uint64_t j = 0; j < S; j++)
			{
				std::swap(mat.v[i*S+j], mat.v[imax*S+j]);
				std::swap(out.v[i*S+j], out.v[imax*S+j]);
			}
			T pivot = mat.v[i*S+i];
			for (uint64_t j = 0; j < S; j++)
			{
				mat.v[i*S+j] /= pivot;
				out.v[i*S+j] /= pivot;
			}
			for (uint64_t j = i+1; j < S; j++)
			{
				T pivot2 = mat[j*S+i];
				for (uint64_t k = 0; k < S; k++)
				{
					mat.v[j*S+k] -= mat.v[i*S+k] * pivot2;
					out.v[j*S+k] -= out.v[i*S+k] * pivot2;
				}
			}
		}
		for (uint64_t i = S-1; i > 0; i--)
		{
			for (uint64_t j = 0; j < i; j++)
			{
				T pivot = mat[j*S+i];
				mat[j*S+i] = 0;
				for (uint64_t k = 0; k < S; k++)
					out[j*S+k] -= out[i*S+k] * pivot;
			}
		}
		return out;
	}
};
