#ifndef _VECTOR4F_H
#define _VECTOR4F_H

#include "Vector3f.h"

#include <float.h>

class Vector4f
{
public:
	union
	{
		struct
		{
			float x;
			float y;
			float z;
			float w;
		};
	};

	// Constructors -----------------------------------------------------------------------------------------------
	Vector4f()
	{
		x = 0.0f;
		y = 0.0f;
		z = 0.0f;
		w = 0.0f;
	}

	Vector4f(const float x, const float y, const float z, const float w) : x(x), y(y), z(z), w(w) {}

	Vector4f(const Vector4f& vector4) :x(vector4.x), y(vector4.y), z(vector4.z), w(vector4.w) {}

	Vector4f(const Vector3f& vector3, float w) :x(vector3.x), y(vector3.y), z(vector3.z), w(w) {}

	~Vector4f() = default;

	float const Magnitude()
	{
		return sqrt(SqrMagnitude());
	}

	float  const SqrMagnitude()
	{
		return ((x * x) + (y * y) + (z * z) + (w * w));
	}

	void Normalize()
	{
		(*this) *= (1.0f / Magnitude() > FLT_EPSILON ? Magnitude() : FLT_EPSILON);
	}

	//converts vector to a formatted string
	std::string to_string()const
	{
		return "x: " + std::to_string(x) + " y: " + std::to_string(y) + " z: " + std::to_string(z) + " w: " + std::to_string(w);
	}

	bool IsValid() const
	{
		return (!std::isnan(x) && !std::isnan(y) && !std::isnan(z) && !std::isnan(w));
	}

	static Vector4f Lerp(Vector4f v1, Vector4f v2, float alpha)
	{
		return (v1 * alpha) + (v2 * (1 - alpha));
	}

	static float Dot(const Vector4f& v1, const Vector4f& v2)
	{
		return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z) + (v1.w * v2.w);
	}

	static float Dot3(const Vector4f& v1, const Vector4f v2)
	{
		return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
	}

	static Vector4f Cross(const Vector4f& v1, const Vector4f& v2)
	{
		Vector4f result;
		result.x = v1.y * v2.z - v1.z * v2.y;
		result.y = v1.z * v2.x - v1.x * v2.z;
		result.z = v1.x * v2.y - v1.y * v2.x;
		result.w = 0.f;
		return result;
	}

	//Operators-----------------------------------------------------------------------------------------

	Vector4f operator*(float scaler)
	{
		return Vector4f(x * scaler, y * scaler, z * scaler, w * scaler);
	}

	float operator*(const Vector4f other)
	{
		return Dot(*this, other);
	}

	Vector4f operator/(float scaler)
	{
		return Vector4f(x / scaler, y / scaler, z / scaler, w / scaler);
	}

	Vector4f operator+(const Vector4f& other)
	{
		return Vector4f(x + other.x, y + other.y, z + other.z, w + other.w);
	}

	Vector4f operator-(const Vector4f& other)
	{
		return Vector4f(x - other.x, y - other.y, z - other.z, w - other.w);
	}

	Vector4f operator-(void)const
	{
		return Vector4f(-x, -y, -z, -w);
	}

	Vector4f operator+=(const Vector4f& other)
	{
		x += other.x;
		y += other.y;
		z += other.z;
		w += other.w;
		return *this;
	}

	Vector4f operator-=(const Vector4f& other)
	{
		x -= other.x;
		y -= other.y;
		z -= other.z;
		w -= other.w;
		return *this;
	}

	Vector4f operator*=(float other)
	{
		x *= other;
		y *= other;
		z *= other;
		w *= other;
		return *this;
	}

	Vector4f operator/=(float other)
	{
		x /= other;
		y /= other;
		z /= other;
		w /= other;
		return *this;
	}

	bool operator==(const Vector4f& other)
	{
		return (x == other.x && y == other.y && z == other.z && w == other.w);
	}

	bool operator!=(const Vector4f& other)
	{
		return !(x == other.x && y == other.y && z == other.z && w == other.w);
	}

	Vector4f operator=(const Vector4f& other)
	{
		x = other.x;
		y = other.y;
		z = other.z;
		w = other.w;
		return *this;
	}

	const float& operator[](const int i) const
	{
		return i == 0 ? this->x : (i == 1 ? this->y : (i == 2 ? this->z : this->w));
	}

	float& operator[](const int i)
	{
		return i == 0 ? this->x : (i == 1 ? this->y : (i == 2 ? this->z : this->w));
	}

	template<typename Archive>
	void serialize(Archive& archive)
	{
		archive(x, y, z, w);
	}
};

inline std::ostream& operator<<(std::ostream& os, Vector4f& v)
{
	return os << v.to_string();
}
#endif // _VECTOR4F_H
