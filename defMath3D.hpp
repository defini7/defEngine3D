#pragma once

namespace def::math
{
	float PI = 3.1415926f;

	float deg2rad(float a) { return a / 180.0f * PI; }
	float rad2deg(float a) { return a / PI * 180.0f; }

	template <class T>
	struct Vector3
	{
		Vector3(T x = (T)0, T y = (T)0, T z = (T)0, T w = (T)1);

		union
		{
			struct { T x, y, z, w; };
			T xyzw[4];
		};

		Vector3<T>& operator=(const Vector3<T>& v);

		Vector3<T> operator+(const Vector3<T>& v) const;
		Vector3<T> operator-(const Vector3<T>& v) const;
		Vector3<T> operator*(const Vector3<T>& v) const;
		Vector3<T> operator/(const Vector3<T>& v) const;
		Vector3<T> operator+(const T& v) const;
		Vector3<T> operator-(const T& v) const;
		Vector3<T> operator*(const T& v) const;
		Vector3<T> operator/(const T& v) const;

		Vector3<T>& operator+=(const Vector3<T>& v);
		Vector3<T>& operator-=(const Vector3<T>& v);
		Vector3<T>& operator*=(const Vector3<T>& v);
		Vector3<T>& operator/=(const Vector3<T>& v);
		Vector3<T>& operator+=(const T& v);
		Vector3<T>& operator-=(const T& v);
		Vector3<T>& operator*=(const T& v);
		Vector3<T>& operator/=(const T& v);

		bool operator==(const Vector3<T>& v) const;
		bool operator!=(const Vector3<T>& v) const;

		friend Vector3<T> operator*(const float& lhs, const Vector3<T>& rhs)
		{
			return Vector3<T>((T)(lhs * (float)rhs.x), (T)(lhs * (float)rhs.y), (T)(lhs * (float)rhs.z));
		}

		friend Vector3<T> operator*(const double& lhs, const Vector3<T>& rhs)
		{
			return Vector3<T>((T)(lhs * (double)rhs.x), (T)(lhs * (double)rhs.y), (T)(lhs * (double)rhs.z));
		}

		friend Vector3<T> operator*(const int& lhs, const Vector3<T>& rhs)
		{
			return Vector3<T>((T)(lhs * (int)rhs.x), (T)(lhs * (int)rhs.y), (T)(lhs * (int)rhs.z));
		}

		friend Vector3<T> operator/(const float& lhs, const Vector3<T>& rhs)
		{
			return Vector3<T>((T)(lhs / (float)rhs.x), (T)(lhs / (float)rhs.y), (T)(lhs / (float)rhs.z));
		}

		friend Vector3<T> operator/(const double& lhs, const Vector3<T>& rhs)
		{
			return Vector3<T>((T)(lhs / (double)rhs.x), (T)(lhs / (double)rhs.y), (T)(lhs / (double)rhs.z));
		}

		friend Vector3<T> operator/(const int& lhs, const Vector3<T>& rhs)
		{
			return Vector3<T>((T)(lhs / (int)rhs.x), (T)(lhs / (int)rhs.y), (T)(lhs / (int)rhs.z));
		}

		operator Vector3<int>()	const;
		operator Vector3<float>() const;
		operator Vector3<double>() const;

		float dot(const Vector3<T>& v) const;
		Vector3<T> cross(const Vector3<T>& v) const;

		T mag() const;
		T mag2() const;

		void swap(Vector3<T>& v);

		Vector3<T> norm() const;
		Vector3<T> abs() const;
		Vector3<T> perp() const;
		Vector3<T> floor() const;
		Vector3<T> ceil() const;
		Vector3<T> cart() const;
		Vector3<T> polar() const;
		Vector3<T>& ref();

		std::string str() const;

		Vector3<T> intersect_plane(Vector3<T>& end, Vector3<T>& planeP, Vector3<T>& planeN, float& t);
	};

	typedef Vector3<int> Vector3i;
	typedef Vector3<float> Vector3f;
	typedef Vector3<double> Vector3d;

	template <class T>
	struct Tex3
	{
		Tex3(T u = (T)0, T v = (T)0, T w = (T)1);

		union
		{
			struct { T u, v, w; };
			T uvw[3];
		};

		Tex3<T>& operator=(const Tex3<T>& v);

		Tex3<T> operator+(const Tex3<T>& v) const;
		Tex3<T> operator-(const Tex3<T>& v) const;
		Tex3<T> operator*(const Tex3<T>& v) const;
		Tex3<T> operator/(const Tex3<T>& v) const;
		Tex3<T> operator+(const T& v) const;
		Tex3<T> operator-(const T& v) const;
		Tex3<T> operator*(const T& v) const;
		Tex3<T> operator/(const T& v) const;

		Tex3<T>& operator+=(const Tex3<T>& v);
		Tex3<T>& operator-=(const Tex3<T>& v);
		Tex3<T>& operator*=(const Tex3<T>& v);
		Tex3<T>& operator/=(const Tex3<T>& v);
		Tex3<T>& operator+=(const T& v);
		Tex3<T>& operator-=(const T& v);
		Tex3<T>& operator*=(const T& v);
		Tex3<T>& operator/=(const T& v);

		bool operator==(const Tex3<T>& v) const;
		bool operator!=(const Tex3<T>& v) const;

		friend Tex3<T> operator*(const float& lhs, const Tex3<T>& rhs)
		{
			return Tex3<T>((T)(lhs * (float)rhs.u), (T)(lhs * (float)rhs.v), (T)(lhs * (float)rhs.w));
		}

		friend Tex3<T> operator*(const double& lhs, const Tex3<T>& rhs)
		{
			return Tex3<T>((T)(lhs * (double)rhs.u), (T)(lhs * (double)rhs.v), (T)(lhs * (double)rhs.w));
		}

		friend Tex3<T> operator*(const int& lhs, const Tex3<T>& rhs)
		{
			return Tex3<T>((T)(lhs * (int)rhs.u), (T)(lhs * (int)rhs.v), (T)(lhs * (int)rhs.w));
		}

		friend Tex3<T> operator/(const float& lhs, const Tex3<T>& rhs)
		{
			return Tex3<T>((T)(lhs / (float)rhs.u), (T)(lhs / (float)rhs.v), (T)(lhs / (float)rhs.w));
		}

		friend Tex3<T> operator/(const double& lhs, const Tex3<T>& rhs)
		{
			return Tex3<T>((T)(lhs / (double)rhs.u), (T)(lhs / (double)rhs.v), (T)(lhs / (double)rhs.w));
		}

		friend Tex3<T> operator/(const int& lhs, const Tex3<T>& rhs)
		{
			return Tex3<T>((T)(lhs / (int)rhs.u), (T)(lhs / (int)rhs.v), (T)(lhs / (int)rhs.w));
		}

		operator Tex3<int>()	const;
		operator Tex3<float>() const;
		operator Tex3<double>() const;

		void swap(Tex3<T>& v);
		Tex3<T>& ref();
	};

	typedef Tex3<int> Tex3i;
	typedef Tex3<float> Tex3f;
	typedef Tex3<double> Tex3d;

	struct Matrix
	{
		Matrix() = default;
		Matrix(float mat[16]);

		void clear();

		void projection(float fovRad, float aspectRatio, float near, float far);
		void identity();
		void point_at(const Vector3f& pos, const Vector3f& target, const Vector3f& up);
		void quick_inverse(const Matrix& mat);
		void inverse(const Matrix& mat);
		void view(const Vector3f& up, const Vector3f& target, const Vector3f& pos, Vector3f& dir, const Vector3f& rot);

		Matrix operator*(const Matrix& mat);
		Vector3f operator*(const Vector3f& v);

		void translate(float x, float y, float z);
		void rotate_x(float a);
		void rotate_y(float a);
		void rotate_z(float a);
		void scale(float x, float y, float z);

		float* operator[](size_t i) { return m[i]; }
		const float* operator[](size_t i) const { return m[i]; }

	private:
		float m[4][4]{ 0.0f };
	};

#ifdef DEF_MATH_3D
#undef DEF_MATH_3D

	template <class T> Vector3<T>::Vector3(T x, T y, T z, T w)
	{
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = w;
	}

	template <class T> Vector3<T>& Vector3<T>::operator=(const Vector3<T>& v)
	{
		this->x = v.x;
		this->y = v.y;
		this->z = v.z;
		this->w = v.w;
		return ref();
	}

	template <class T> Vector3<T> Vector3<T>::operator+(const Vector3<T>& v) const { return Vector3<T>(this->x + v.x, this->y + v.y, this->z + v.z); }
	template <class T> Vector3<T> Vector3<T>::operator-(const Vector3<T>& v) const { return Vector3<T>(this->x - v.x, this->y - v.y, this->z - v.z); }
	template <class T> Vector3<T> Vector3<T>::operator*(const Vector3<T>& v) const { return Vector3<T>(this->x * v.x, this->y * v.y, this->z * v.z); }
	template <class T> Vector3<T> Vector3<T>::operator/(const Vector3<T>& v) const { return Vector3<T>(this->x / v.x, this->y / v.y, this->z / v.z); }
	template <class T> Vector3<T> Vector3<T>::operator+(const T& v) const { return Vector3<T>(this->x + v, this->y + v, this->z + v); }
	template <class T> Vector3<T> Vector3<T>::operator-(const T& v) const { return Vector3<T>(this->x - v, this->y - v, this->z - v); }
	template <class T> Vector3<T> Vector3<T>::operator*(const T& v) const { return Vector3<T>(this->x * v, this->y * v, this->z * v); }
	template <class T> Vector3<T> Vector3<T>::operator/(const T& v) const { return Vector3<T>(this->x / v, this->y / v, this->z / v); }

	template <class T> Vector3<T>& Vector3<T>::operator+=(const Vector3<T>& v)
	{
		this->x += v.x;
		this->y += v.y;
		this->z += v.z;
		return ref();
	}

	template <class T> Vector3<T>& Vector3<T>::operator-=(const Vector3<T>& v)
	{
		this->x -= v.x;
		this->y -= v.y;
		this->z -= v.z;
		return ref();
	}

	template <class T> Vector3<T>& Vector3<T>::operator*=(const Vector3<T>& v)
	{
		this->x *= v.x;
		this->y *= v.y;
		this->z *= v.z;
		return ref();
	}

	template <class T> Vector3<T>& Vector3<T>::operator/=(const Vector3<T>& v)
	{
		this->x /= v.x;
		this->y /= v.y;
		this->z /= v.z;
		return ref();
	}

	template <class T> Vector3<T>& Vector3<T>::operator+=(const T& v)
	{
		this->x += v;
		this->y += v;
		this->z += v;
		return ref();
	}

	template <class T> Vector3<T>& Vector3<T>::operator-=(const T& v)
	{
		this->x -= v;
		this->y -= v;
		this->z -= v;
		return ref();
	}

	template <class T> Vector3<T>& Vector3<T>::operator*=(const T& v)
	{
		this->x *= v;
		this->y *= v;
		this->z *= v;
		return ref();
	}

	template <class T> Vector3<T>& Vector3<T>::operator/=(const T& v)
	{
		this->x /= v;
		this->y /= v;
		this->z /= v;
		return ref();
	}

	template <class T> bool Vector3<T>::operator==(const Vector3<T>& v) const { return this->x == v.x && this->y == v.y && this->z == v.z; }
	template <class T> bool Vector3<T>::operator!=(const Vector3<T>& v) const { return this->x != v.x || this->y != v.y || this->z != v.z; }

	template <class T> Vector3<T>::operator Vector3<int>()	const { return { static_cast<int32_t>(this->x), static_cast<int32_t>(this->y), static_cast<int32_t>(this->z), static_cast<int32_t>(this->w) }; }
	template <class T> Vector3<T>::operator Vector3<float>() const { return { static_cast<float>(this->x), static_cast<float>(this->y), static_cast<float>(this->z), static_cast<float>(this->w) }; }
	template <class T> Vector3<T>::operator Vector3<double>() const { return { static_cast<double>(this->x), static_cast<double>(this->y), static_cast<double>(this->z), static_cast<double>(this->w) }; }

	template <class T> float Vector3<T>::dot(const Vector3<T>& v) const { return this->x * v.x + this->y * v.y + this->z * v.z; }

	template <class T> Vector3<T> Vector3<T>::cross(const Vector3<T>& v) const { return { this->y * v.z - this->z * v.y, this->z * v.x - this->x * v.z, this->x * v.y - this->y * v.x }; }

	template <class T> T Vector3<T>::mag() const { return static_cast<T>(sqrtf(this->dot(*this))); }
	template <class T> T Vector3<T>::mag2() const { return static_cast<T>(this->dot(*this)); }

	template <class T> void Vector3<T>::swap(Vector3<T>& v)
	{
		std::swap(x, v.x);
		std::swap(y, v.y);
		std::swap(z, v.z);
		std::swap(w, v.w);
	}

	template <class T> Vector3<T> Vector3<T>::norm() const { float l = mag(); return Vector3<T>(this->x / l, this->y / l, this->z / l); }
	template <class T> Vector3<T> Vector3<T>::abs() const { return Vector3<T>(std::abs(this->x), std::abs(this->y), std::abs(this->z)); }
	template <class T> Vector3<T> Vector3<T>::perp() const { return Vector3<T>(-this->y, this->x, this->z); }
	template <class T> Vector3<T> Vector3<T>::floor() const { return Vector3<T>(std::floor(this->x), std::floor(this->y), std::floor(this->z)); }
	template <class T> Vector3<T> Vector3<T>::ceil() const { return Vector3<T>(std::ceil(this->x), std::ceil(this->y), std::ceil(this->z)); }
	template <class T> Vector3<T> Vector3<T>::cart() const { return Vector3<T>(cos(this->y) * this->x, sin(this->y) * this->x, this->z); }
	template <class T> Vector3<T> Vector3<T>::polar() const { return Vector3<T>(mag(), atan2(this->y, this->x), this->z); }
	template <class T> Vector3<T>& Vector3<T>::ref() { return *this; }

	template <class T> std::string Vector3<T>::str() const { return "(" + std::to_string(this->x) + ", " + std::to_string(this->y) + ", " + std::to_string(this->z) + ")"; }

	template <class T> Vector3<T> Vector3<T>::intersect_plane(Vector3<T>& end, Vector3<T>& planeP, Vector3<T>& planeN, float& t)
	{
		planeN = planeN.norm();
		float planeD = -planeN.dot(planeP);
		float ad = dot(planeN);
		float bd = end.dot(planeN);
		t = (-planeD - ad) / (bd - ad);
		return ref() + (end - ref()) * t;
	}

	template <class T> Tex3<T>::Tex3(T u, T v, T w)
	{
		this->u = u;
		this->v = v;
		this->w = w;
	}

	template <class T> Tex3<T>& Tex3<T>::operator=(const Tex3<T>& v)
	{
		this->u = v.u;
		this->v = v.v;
		this->w = v.w;
		return ref();
	}

	template <class T> Tex3<T> Tex3<T>::operator+(const Tex3<T>& v) const { return Tex3<T>(this->u + v.u, this->v + v.v, this->w + v.w); }
	template <class T> Tex3<T> Tex3<T>::operator-(const Tex3<T>& v) const { return Tex3<T>(this->u - v.u, this->v - v.v, this->w - v.w); }
	template <class T> Tex3<T> Tex3<T>::operator*(const Tex3<T>& v) const { return Tex3<T>(this->u * v.u, this->v * v.v, this->w * v.w); }
	template <class T> Tex3<T> Tex3<T>::operator/(const Tex3<T>& v) const { return Tex3<T>(this->u / v.u, this->v / v.v, this->w / v.w); }
	template <class T> Tex3<T> Tex3<T>::operator+(const T& v) const { return Tex3<T>(this->u + v, this->v + v, this->w + v); }
	template <class T> Tex3<T> Tex3<T>::operator-(const T& v) const { return Tex3<T>(this->u - v, this->v - v, this->w - v); }
	template <class T> Tex3<T> Tex3<T>::operator*(const T& v) const { return Tex3<T>(this->u * v, this->v * v, this->w * v); }
	template <class T> Tex3<T> Tex3<T>::operator/(const T& v) const { return Tex3<T>(this->u / v, this->v / v, this->w / v); }

	template <class T> Tex3<T>& Tex3<T>::operator+=(const Tex3<T>& v)
	{
		this->u += v.u;
		this->v += v.v;
		this->w += v.w;
		return ref();
	}

	template <class T> Tex3<T>& Tex3<T>::operator-=(const Tex3<T>& v)
	{
		this->u -= v.u;
		this->v -= v.v;
		this->w -= v.w;
		return ref();
	}

	template <class T> Tex3<T>& Tex3<T>::operator*=(const Tex3<T>& v)
	{
		this->u *= v.u;
		this->v *= v.v;
		this->w *= v.w;
		return ref();
	}

	template <class T> Tex3<T>& Tex3<T>::operator/=(const Tex3<T>& v)
	{
		this->u /= v.u;
		this->v /= v.v;
		this->w /= v.w;
		return ref();
	}

	template <class T> Tex3<T>& Tex3<T>::operator+=(const T& v)
	{
		this->u += v;
		this->v += v;
		this->w += v;
		return ref();
	}

	template <class T> Tex3<T>& Tex3<T>::operator-=(const T& v)
	{
		this->u -= v;
		this->v -= v;
		this->w -= v;
		return ref();
	}

	template <class T> Tex3<T>& Tex3<T>::operator*=(const T& v)
	{
		this->u *= v;
		this->v *= v;
		this->w *= v;
		return ref();
	}

	template <class T> Tex3<T>& Tex3<T>::operator/=(const T& v)
	{
		this->u /= v;
		this->v /= v;
		this->w /= v;
		return ref();
	}

	template <class T> bool Tex3<T>::operator==(const Tex3<T>& v) const { return this->u == v.u && this->v == v.v && this->w == v.w; }
	template <class T> bool Tex3<T>::operator!=(const Tex3<T>& v) const { return this->u != v.u || this->v != v.v || this->w != v.w; }

	template <class T> Tex3<T>::operator Tex3<int>() const { return { static_cast<int32_t>(this->u), static_cast<int32_t>(this->v), static_cast<int32_t>(this->w) }; }
	template <class T> Tex3<T>::operator Tex3<float>() const { return { static_cast<float>(this->u), static_cast<float>(this->v), static_cast<float>(this->w) }; }
	template <class T> Tex3<T>::operator Tex3<double>() const { return { static_cast<double>(this->u), static_cast<double>(this->v), static_cast<double>(this->w) }; }

	template <class T> void Tex3<T>::swap(Tex3<T>& v)
	{
		std::swap(this->u, v.u);
		std::swap(this->v, v.v);
		std::swap(this->w, v.w);
	}

	template <class T> Tex3<T>& Tex3<T>::ref() { return *this; }

	Matrix::Matrix(float mat[16])
	{
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 4; y++)
				m[x][y] = mat[y * 4 + x];
	}

	void Matrix::clear()
	{
		memset(m, 0, sizeof(float) * 16);
	}

	void Matrix::projection(float fovDeg, float aspectRatio, float near, float far)
	{
		float fov = 1.0f / tanf(fovDeg * 0.5f / 180.0f * 3.14159f);
		m[0][0] = aspectRatio * fov;
		m[1][1] = fov;
		m[2][2] = far / (far - near);
		m[3][2] = (-far * near) / (far - near);
		m[2][3] = 1.0f;
		m[3][3] = 0.0f;
	}

	void Matrix::identity()
	{
		m[0][0] = 1.0f;
		m[1][1] = 1.0f;
		m[2][2] = 1.0f;
		m[3][3] = 1.0f;
	}

	void Matrix::translate(float x, float y, float z)
	{
		m[0][0] = 1.0f;
		m[1][1] = 1.0f;
		m[2][2] = 1.0f;
		m[3][3] = 1.0f;
		m[3][0] = x;
		m[3][1] = y;
		m[3][2] = z;
	}

	void Matrix::point_at(const Vector3f& pos, const Vector3f& target, const Vector3f& up)
	{
		Vector3f newForward = (target - pos).norm();
		Vector3f newUp = (up - newForward * up.dot(newForward)).norm();

		Vector3f newRight = newUp.cross(newForward);

		m[0][0] = newRight.x;	m[0][1] = newRight.y;	m[0][2] = newRight.z;	m[0][3] = 0.0f;
		m[1][0] = newUp.x;		m[1][1] = newUp.y;		m[1][2] = newUp.z;		m[1][3] = 0.0f;
		m[2][0] = newForward.x;	m[2][1] = newForward.y;	m[2][2] = newForward.z;	m[2][3] = 0.0f;
		m[3][0] = pos.x;		m[3][1] = pos.y;		m[3][2] = pos.z;		m[3][3] = 1.0f;
	}

	void Matrix::quick_inverse(const Matrix& mat)
	{
		m[0][0] = mat.m[0][0]; m[0][1] = mat.m[1][0]; m[0][2] = mat.m[2][0]; m[0][3] = 0.0f;
		m[1][0] = mat.m[0][1]; m[1][1] = mat.m[1][1]; m[1][2] = mat.m[2][1]; m[1][3] = 0.0f;
		m[2][0] = mat.m[0][2]; m[2][1] = mat.m[1][2]; m[2][2] = mat.m[2][2]; m[2][3] = 0.0f;
		m[3][0] = -(mat.m[3][0] * m[0][0] + mat.m[3][1] * m[1][0] + mat.m[3][2] * m[2][0]);
		m[3][1] = -(mat.m[3][0] * m[0][1] + mat.m[3][1] * m[1][1] + mat.m[3][2] * m[2][1]);
		m[3][2] = -(mat.m[3][0] * m[0][2] + mat.m[3][1] * m[1][2] + mat.m[3][2] * m[2][2]);
		m[3][3] = 1.0f;
	}

	void Matrix::inverse(const Matrix& mat)
	{
		m[0][0] = mat.m[1][1] * mat.m[2][2] * mat.m[3][3] - mat.m[1][1] * mat.m[2][3] * mat.m[3][2] - mat.m[2][1] * mat.m[1][2] * mat.m[3][3] + mat.m[2][1] * mat.m[1][3] * mat.m[3][2] + mat.m[3][1] * mat.m[1][2] * mat.m[2][3] - mat.m[3][1] * mat.m[1][3] * mat.m[2][2];
		m[1][0] = -mat.m[1][0] * mat.m[2][2] * mat.m[3][3] + mat.m[1][0] * mat.m[2][3] * mat.m[3][2] + mat.m[2][0] * mat.m[1][2] * mat.m[3][3] - mat.m[2][0] * mat.m[1][3] * mat.m[3][2] - mat.m[3][0] * mat.m[1][2] * mat.m[2][3] + mat.m[3][0] * mat.m[1][3] * mat.m[2][2];
		m[2][0] = mat.m[1][0] * mat.m[2][1] * mat.m[3][3] - mat.m[1][0] * mat.m[2][3] * mat.m[3][1] - mat.m[2][0] * mat.m[1][1] * mat.m[3][3] + mat.m[2][0] * mat.m[1][3] * mat.m[3][1] + mat.m[3][0] * mat.m[1][1] * mat.m[2][3] - mat.m[3][0] * mat.m[1][3] * mat.m[2][1];
		m[3][0] = -mat.m[1][0] * mat.m[2][1] * mat.m[3][2] + mat.m[1][0] * mat.m[2][2] * mat.m[3][1] + mat.m[2][0] * mat.m[1][1] * mat.m[3][2] - mat.m[2][0] * mat.m[1][2] * mat.m[3][1] - mat.m[3][0] * mat.m[1][1] * mat.m[2][2] + mat.m[3][0] * mat.m[1][2] * mat.m[2][1];
		m[0][1] = -mat.m[0][1] * mat.m[2][2] * mat.m[3][3] + mat.m[0][1] * mat.m[2][3] * mat.m[3][2] + mat.m[2][1] * mat.m[0][2] * mat.m[3][3] - mat.m[2][1] * mat.m[0][3] * mat.m[3][2] - mat.m[3][1] * mat.m[0][2] * mat.m[2][3] + mat.m[3][1] * mat.m[0][3] * mat.m[2][2];
		m[1][1] = mat.m[0][0] * mat.m[2][2] * mat.m[3][3] - mat.m[0][0] * mat.m[2][3] * mat.m[3][2] - mat.m[2][0] * mat.m[0][2] * mat.m[3][3] + mat.m[2][0] * mat.m[0][3] * mat.m[3][2] + mat.m[3][0] * mat.m[0][2] * mat.m[2][3] - mat.m[3][0] * mat.m[0][3] * mat.m[2][2];
		m[2][1] = -mat.m[0][0] * mat.m[2][1] * mat.m[3][3] + mat.m[0][0] * mat.m[2][3] * mat.m[3][1] + mat.m[2][0] * mat.m[0][1] * mat.m[3][3] - mat.m[2][0] * mat.m[0][3] * mat.m[3][1] - mat.m[3][0] * mat.m[0][1] * mat.m[2][3] + mat.m[3][0] * mat.m[0][3] * mat.m[2][1];
		m[3][1] = mat.m[0][0] * mat.m[2][1] * mat.m[3][2] - mat.m[0][0] * mat.m[2][2] * mat.m[3][1] - mat.m[2][0] * mat.m[0][1] * mat.m[3][2] + mat.m[2][0] * mat.m[0][2] * mat.m[3][1] + mat.m[3][0] * mat.m[0][1] * mat.m[2][2] - mat.m[3][0] * mat.m[0][2] * mat.m[2][1];
		m[0][2] = mat.m[0][1] * mat.m[1][2] * mat.m[3][3] - mat.m[0][1] * mat.m[1][3] * mat.m[3][2] - mat.m[1][1] * mat.m[0][2] * mat.m[3][3] + mat.m[1][1] * mat.m[0][3] * mat.m[3][2] + mat.m[3][1] * mat.m[0][2] * mat.m[1][3] - mat.m[3][1] * mat.m[0][3] * mat.m[1][2];
		m[1][2] = -mat.m[0][0] * mat.m[1][2] * mat.m[3][3] + mat.m[0][0] * mat.m[1][3] * mat.m[3][2] + mat.m[1][0] * mat.m[0][2] * mat.m[3][3] - mat.m[1][0] * mat.m[0][3] * mat.m[3][2] - mat.m[3][0] * mat.m[0][2] * mat.m[1][3] + mat.m[3][0] * mat.m[0][3] * mat.m[1][2];
		m[2][2] = mat.m[0][0] * mat.m[1][1] * mat.m[3][3] - mat.m[0][0] * mat.m[1][3] * mat.m[3][1] - mat.m[1][0] * mat.m[0][1] * mat.m[3][3] + mat.m[1][0] * mat.m[0][3] * mat.m[3][1] + mat.m[3][0] * mat.m[0][1] * mat.m[1][3] - mat.m[3][0] * mat.m[0][3] * mat.m[1][1];
		m[3][2] = -mat.m[0][0] * mat.m[1][1] * mat.m[3][2] + mat.m[0][0] * mat.m[1][2] * mat.m[3][1] + mat.m[1][0] * mat.m[0][1] * mat.m[3][2] - mat.m[1][0] * mat.m[0][2] * mat.m[3][1] - mat.m[3][0] * mat.m[0][1] * mat.m[1][2] + mat.m[3][0] * mat.m[0][2] * mat.m[1][1];
		m[0][3] = -mat.m[0][1] * mat.m[1][2] * mat.m[2][3] + mat.m[0][1] * mat.m[1][3] * mat.m[2][2] + mat.m[1][1] * mat.m[0][2] * mat.m[2][3] - mat.m[1][1] * mat.m[0][3] * mat.m[2][2] - mat.m[2][1] * mat.m[0][2] * mat.m[1][3] + mat.m[2][1] * mat.m[0][3] * mat.m[1][2];
		m[1][3] = mat.m[0][0] * mat.m[1][2] * mat.m[2][3] - mat.m[0][0] * mat.m[1][3] * mat.m[2][2] - mat.m[1][0] * mat.m[0][2] * mat.m[2][3] + mat.m[1][0] * mat.m[0][3] * mat.m[2][2] + mat.m[2][0] * mat.m[0][2] * mat.m[1][3] - mat.m[2][0] * mat.m[0][3] * mat.m[1][2];
		m[2][3] = -mat.m[0][0] * mat.m[1][1] * mat.m[2][3] + mat.m[0][0] * mat.m[1][3] * mat.m[2][1] + mat.m[1][0] * mat.m[0][1] * mat.m[2][3] - mat.m[1][0] * mat.m[0][3] * mat.m[2][1] - mat.m[2][0] * mat.m[0][1] * mat.m[1][3] + mat.m[2][0] * mat.m[0][3] * mat.m[1][1];
		m[3][3] = mat.m[0][0] * mat.m[1][1] * mat.m[2][2] - mat.m[0][0] * mat.m[1][2] * mat.m[2][1] - mat.m[1][0] * mat.m[0][1] * mat.m[2][2] + mat.m[1][0] * mat.m[0][2] * mat.m[2][1] + mat.m[2][0] * mat.m[0][1] * mat.m[1][2] - mat.m[2][0] * mat.m[0][2] * mat.m[1][1];

		float det = 1.0f / float(mat.m[0][0] * m[0][0] + mat.m[0][1] * m[1][0] + mat.m[0][2] * m[2][0] + mat.m[0][3] * m[3][0]);

		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				m[i][j] *= (float)det;
	}

	Matrix Matrix::operator*(const Matrix& mat)
	{
		Matrix res;

		for (int c = 0; c < 4; c++)
			for (int r = 0; r < 4; r++)
				res.m[r][c] = m[r][0] * mat.m[0][c] + m[r][1] * mat.m[1][c] + m[r][2] * mat.m[2][c] + m[r][3] * mat.m[3][c];

		return res;
	}


	Vector3f Matrix::operator*(const Vector3f& v)
	{
		Vector3f res;
		res.x = v.x * m[0][0] + v.y * m[1][0] + v.z * m[2][0] + v.w * m[3][0];
		res.y = v.x * m[0][1] + v.y * m[1][1] + v.z * m[2][1] + v.w * m[3][1];
		res.z = v.x * m[0][2] + v.y * m[1][2] + v.z * m[2][2] + v.w * m[3][2];
		res.w = v.x * m[0][3] + v.y * m[1][3] + v.z * m[2][3] + v.w * m[3][3];
		return res;
	}

	void Matrix::rotate_x(float a)
	{
		m[0][0] = 1.0f;
		m[1][1] = cosf(a);
		m[1][2] = sinf(a);
		m[2][1] = -sinf(a);
		m[2][2] = cosf(a);
		m[3][3] = 1.0f;
	}

	void Matrix::rotate_y(float a)
	{
		m[0][0] = cosf(a);
		m[0][2] = -sinf(a);
		m[1][1] = 1.0f;
		m[2][0] = sinf(a);
		m[2][2] = cosf(a);
		m[3][3] = 1.0f;
	}

	void Matrix::rotate_z(float a)
	{
		m[0][0] = cosf(a);
		m[0][1] = sinf(a);
		m[1][0] = -sinf(a);
		m[1][1] = cosf(a);
		m[2][2] = 1.0f;
		m[3][3] = 1.0f;
	}

	void Matrix::scale(float x, float y, float z)
	{
		m[0][0] = x;
		m[1][1] = y;
		m[2][2] = z;
		m[3][3] = 1.0f;
	}

	void Matrix::view(const Vector3f& up, const Vector3f& target, const Vector3f& pos, Vector3f& dir, const Vector3f& rot)
	{
		Matrix rotMatX, rotMatY, rotMatZ;
		rotMatX.rotate_x(rot.x);
		rotMatY.rotate_y(rot.y);
		rotMatZ.rotate_z(rot.z);

		dir = rotMatX * rotMatY * rotMatZ * target;

		Matrix cameraMat;
		cameraMat.point_at(pos, pos + dir, up);

		quick_inverse(cameraMat);
	}

#endif
}
