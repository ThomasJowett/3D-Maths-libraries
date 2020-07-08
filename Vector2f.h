#ifndef _VECTOR2F_H
#define _VECTOR2F_H

#include <string>
#include <cmath>
#include "Vector3f.h"

class Vector2f
{
public:
	union
	{
		struct {
			float x;
			float y;
		};
	};

	//Constructors ----------------------------------------------------------------------------------
	Vector2f()
	{
		x = 0.0f;
		y = 0.0f;
	}

	Vector2f(const float x, const float y) :x(x), y(y) {}

	Vector2f(const Vector2f& vector2) :x(vector2.x), y(vector2.y) {}
	
	Vector2f(const Vector3f& vector3) :x(vector3.x), y(vector3.y) {}

	//Destructor---------------------------------------------------------------------------------------
	~Vector2f() = default;

	//Member Functions --------------------------------------------------------------------------------

	//Length of the vector
	float const Magnitude()
	{
		return sqrt(SqrMagnitude());
	}

	//Square length of the vector
	float const SqrMagnitude()
	{
		return ((x * x) + (y * y));
	}

	//Gets a vector of same direction with magnitude of 1
	Vector2f const GetNormalized()
	{
		float magnitude = Magnitude();
		if (magnitude > 0.0f)
			return Vector2f(x / magnitude, y / magnitude);
		return Vector2f();
	}

	//set this vectors length to 1
	void Normalize()
	{
		Vector2f normalized = GetNormalized();
		x = normalized.x;
		y = normalized.y;
	}

	//clamp this vector to given length
	void Clamp(const float length)
	{
		if (SqrMagnitude() > length * length)
		{
			Normalize();
			*this = *this * length;
		}
	}

	//set x and y to zero
	void Zero()
	{
		x = 0.0f;
		y = 0.0f;
	}

	//Returns Perpendicular either clockwise or anti-clockwise
	Vector2f const Perpendicular(bool clockwise)
	{
		if (clockwise)
			return Vector2f(y, -x);
		else
			return Vector2f(-y, x);
	}

	//converts vector to a formatted string
	std::string to_string() const
	{
		return "x: " + std::to_string(x) + " y: " + std::to_string(y);
	}

	Vector3f to_Vector3D() const
	{
		return Vector3f(x, y, 0);
	}

	bool IsValid() const
	{
		return (!std::isnan(x) && !std::isnan(y));
	}

	//Static-----------------------------------------------------------------------------------------

	//the unsigned angle between v1 and v2
	static float Angle(Vector2f v1, Vector2f v2)
	{
		return acos(Dot(v1, v2));
	}

	//distance between v1 and v2
	static float Distance(Vector2f v1, Vector2f v2)
	{
		return (v1 - v2).Magnitude();
	}

	//square distance between v1 and v2
	static float SqrDistance(Vector2f v1, Vector2f v2)
	{
		return (v1 - v2).SqrMagnitude();
	}

	//returns the sum of the products of v1 and v2
	static float Dot(Vector2f v1, Vector2f v2)
	{
		return v1.x * v2.x + v1.y * v2.y;
	}

	// 2D vector cross product analog.
	// The cross product of 2D vectors results in a 3D vector with only a z component.
	// This function returns the magnitude of the z value.
	static float Cross(Vector2f v1, Vector2f v2)
	{
		return v1.x * v2.y - v1.y * v2.x;
	}

	//linearly interpolate between v1 and v2
	static Vector2f Lerp(Vector2f v1, Vector2f v2, float alpha)
	{
		return (v1 * alpha) + (v2 * (1 - alpha));
	}

	//spherically interpolate between v1 and v2
	//static Vector2f Slerp(Vector2f v1, Vector2f v2, float alpha)
	//{
	//	float theta = Angle(v1, v2);
	//	return Rotated()
	//	float dot = Dot(v1, v2);
	//	float omega = acos(Clamp(dot, -1)
	//}

	//returns a vector that is the reflection v against the normal
	static Vector2f Reflect(Vector2f v, Vector2f normal)
	{
		Vector2f n = normal.GetNormalized();
		return v - (n * (2.0f * Dot(v, n)));
	}

	// return projection v1 on to v2
	static Vector2f Projection(Vector2f v1, Vector2f v2)
	{
		float v2SqrMagnitude = v2.SqrMagnitude();
		if (v2SqrMagnitude > 1.0e-8f)
			return v2 * (Dot(v2, v1) / v2SqrMagnitude);
		else
			return Vector2f();
	}

	//returns unit length vector for given angle (in radians)
	static Vector2f ForAngle(const float angle)
	{
		return Vector2f(cos(angle), sin(angle));
	}

	//returns if v1 is close to v2
	static bool Near(Vector2f v1, Vector2f v2, float distance)
	{
		return SqrDistance(v1, v2) < distance * distance;
	}

	//Operators--------------------------------------------------------------------------------------

	Vector2f operator*(float scaler)
	{
		return Vector2f(x * scaler, y * scaler);
	}

	float operator*(Vector2f other)
	{
		return Dot(*this, other);
	}

	Vector2f operator/(float scaler)
	{
		return Vector2f(x / scaler, y / scaler);
	}

	Vector2f operator+(const Vector2f& other)
	{
		return Vector2f(x + other.x, y + other.y);
	}

	Vector2f operator-(const Vector2f& other)
	{
		return Vector2f(x - other.x, y - other.y);
	}

	Vector2f operator-(void)const
	{
		return Vector2f(-x, -y);
	}

	Vector2f operator+=(const Vector2f& other)
	{
		x += other.x;
		y += other.y;
		return *this;
	}

	Vector2f operator-=(const Vector2f& other)
	{
		x -= other.x;
		y -= other.y;
		return *this;
	}

	bool operator==(const Vector2f& other)
	{
		return (x == other.x && y == other.y);
	}

	bool operator!=(const Vector2f& other)
	{
		return !(x == other.x && y == other.y);
	}

	Vector2f operator=(const Vector2f& other)
	{
		x = other.x;
		y = other.y;
		return *this;
	}

	Vector2f operator=(const Vector3f& other)
	{
		x = other.x;
		y = other.y;
		return *this;
	}

	operator bool() const
	{
		return(IsValid());
	}

	const float& operator[](const int i) const
	{
		return i == 0 ? this->x : this->y;
	}

	float& operator[](const int i)
	{
		return i == 0 ? this->x : this->y;
	}
};

inline Vector2f operator*(float scaler, const Vector2f& v)
{
	return Vector2f(scaler * v.x, scaler * v.y);
}

inline Vector2f operator/(float scaler, const Vector2f& v)
{
	return Vector2f(scaler / v.x, scaler / v.y);
}

inline Vector2f operator+(const Vector2f& v1, const Vector2f& v2)
{
	return Vector2f(v1.x + v2.x, v1.y + v2.y);
}

inline Vector2f operator-(const Vector2f& v1, const Vector2f& v2)
{
	return Vector2f(v1.x - v2.x, v1.y - v2.y);
}

inline std::ostream& operator<<(std::ostream& os, Vector2f& v)
{
	return os << v.to_string();
}

#endif // !_VECTOR2F_H