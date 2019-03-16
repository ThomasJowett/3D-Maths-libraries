#ifndef _VECTOR2F_H
#define _VECTOR2F_H

#include <string>
#include <math.h>
#include "Vector3f.h"

class Vector2f
{
public:
	float x;
	float y;

	//Constructors ----------------------------------------------------------------------------------
	Vector2f()
	{
		x = 0.0f;
		y = 0.0f;
	}

	Vector2f(float x, float y) :x(x), y(y) {}

	Vector2f(Vector3D vector3d) :x(vector3d.x), y(vector3d.y) {}

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
		return ((x*x) + (y*y));
	}

	//Gets a vector of same direction with magnitude of 1
	Vector2f const GetNormalized()
	{
		float magnitude = Magnitude();
		return Vector2f(x / magnitude, y / magnitude);
	}

	//set this vectors length to 1
	void Normalize()
	{
		Vector2f normalized = GetNormalized();
		x = normalized.x;
		y = normalized.y;
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
	std::string to_string()
	{
		return "x: " + std::to_string(x) + " y: " + std::to_string(y);
	}

	Vector3D to_Vector3D()
	{
		return Vector3D(x, y, 0);
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

	//returns the sum of the products of v1 and v2
	static float Dot(Vector2f v1, Vector2f v2)
	{
		return v1.x*v2.x + v1.y*v2.y;
	}

	static float Cross(Vector2f v1, Vector2f v2)
	{
		return v1.x * v2.y - v1.y * v2.x;
	}

	//linearly interpolate between v1 and v2
	static Vector2f Lerp(Vector2f v1, Vector2f v2, float alpha)
	{
		return (v1 * alpha) + (v2 * (1 - alpha));
	}

	//returns a vector that is the reflection v against the normal
	static Vector2f Reflect(Vector2f v, Vector2f normal)
	{
		Vector2f n = normal.GetNormalized();
		return v - (n*(2.0f * Dot(v, n)));
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

	//Operators--------------------------------------------------------------------------------------

	Vector2f operator*(float scaler)
	{
		return Vector2f(x * scaler, y *scaler);
	}

	Vector2f operator*(Vector2f other)
	{
		return Vector2f(x * other.x, y * other.y);
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

	Vector2f operator=(const Vector3D& other)
	{
		x = other.x;
		y = other.y;
		return *this;
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

#endif // !_VECTOR2F_H