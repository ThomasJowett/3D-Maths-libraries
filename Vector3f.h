#ifndef _VECTOR3F_H
#define _VECTOR3F_H

#include <string>
#include <cmath>

class Vector3f
{
public:
	union
	{
		struct {
			float x;
			float y;
			float z;
		};
	};

	//Constructors ----------------------------------------------------------------------------------
	Vector3f()
	{
		x = 0.0f;
		y = 0.0f;
		z = 0.0f;
	}

	Vector3f(const float x, const float y, const float z) : x(x), y(y), z(z) {}

	Vector3f(const Vector3f& vector3) :x(vector3.x), y(vector3.y), z(vector3.z) {}

	//Destructor---------------------------------------------------------------------------------------
	~Vector3f() = default;

	//Member Functions --------------------------------------------------------------------------------

	//Length of the vector
	float const Magnitude()
	{
		return sqrt(SqrMagnitude());
	}

	//Square length of the vector
	float const SqrMagnitude()
	{
		return ((x * x) + (y * y) + (z * z));
	}

	//Gets a vector of same direction with magnitude of 1
	Vector3f const GetNormalized()
	{
		float magnitude = Magnitude();
		return Vector3f(x / magnitude, y / magnitude, z / magnitude);
	}

	//makes the magnitude of this vector 1
	void Normalize()
	{
		Vector3f normalized = GetNormalized();

		x = normalized.x;
		y = normalized.y;
		z = normalized.z;
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

	//set x, y and z to zero
	void Zero()
	{
		x = 0.0f;
		y = 0.0f;
		z = 0.0f;
	}

	//converts vector to a formatted string
	std::string to_string()const
	{
		return "x: " + std::to_string(x) + " y: " + std::to_string(y) + " z: " + std::to_string(z);
	}

	bool IsValid() const
	{
		return (!std::isnan(x) && !std::isnan(y) && !std::isnan(z));
	}

	//Static-------------------------------------------------------------------------------------------

	//returns a vector orthagonal to both v1 and v2
	static Vector3f Cross(Vector3f v1, Vector3f v2)
	{
		Vector3f cross;
		cross.x = (v1.y * v2.z) - (v1.z * v2.y);
		cross.y = -((v1.x * v2.z) - (v1.z * v2.x));
		cross.z = (v1.x * v2.y) - (v1.y * v2.x);
		return cross;
	}

	//returns the sum of the products of v1 and v2
	static float Dot(Vector3f v1, Vector3f v2)
	{
		return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
	}

	//distance between v1 and v2
	static float Distance(Vector3f v1, Vector3f v2)
	{
		return (v1 - v2).Magnitude();
	}

	//returns a vector that is the reflection v against the normal
	static Vector3f Reflect(Vector3f v, Vector3f normal)
	{
		normal.Normalize();
		return v - (normal * (2.0f * Dot(v, normal)));
	}

	//linearly interpolate between v1 and v2
	static Vector3f Lerp(Vector3f v1, Vector3f v2, float alpha)
	{
		return (v1 * alpha) + (v2 * (1 - alpha));
	}

	//spherically interpolate between v1 and v2
	static Vector3f Slerp(Vector3f v1, Vector3f v2, float alpha)
	{
		float dot = Dot(v1, v2);
		float theta = acos(dot) * alpha;

		Vector3f relativeVec = v2 - v1 * dot;
		relativeVec.Normalize();

		return ((v1 * cos(theta)) + (relativeVec * sin(theta)));
	}

	//Operators-----------------------------------------------------------------------------------------

	//multiples each component of the vector by the scaler
	Vector3f operator*(float scaler)
	{
		return Vector3f(x * scaler, y * scaler, z * scaler);
	}

	// dot product
	float operator*(const Vector3f other)
	{
		return Dot(*this, other);
	}

	//divides each component of the vector by the scaler
	Vector3f operator/(float scaler)
	{
		return Vector3f(x / scaler, y / scaler, z / scaler);
	}

	//adds the two vectors together
	Vector3f operator+(const Vector3f& other)
	{
		return Vector3f(x + other.x, y + other.y, z + other.z);
	}

	//subtracts v2 from v1
	Vector3f operator-(const Vector3f& other)
	{
		return Vector3f(x - other.x, y - other.y, z - other.z);
	}

	//the inverse of the vector
	Vector3f operator-(void)const
	{
		return Vector3f(-x, -y, -z);
	}

	Vector3f operator+=(const Vector3f& other)
	{
		x += other.x;
		y += other.y;
		z += other.z;
		return *this;
	}

	Vector3f operator-=(const Vector3f& other)
	{
		x -= other.x;
		y -= other.y;
		z -= other.z;
		return *this;
	}

	bool operator==(const Vector3f& other)
	{
		return (x == other.x && y == other.y && z == other.z);
	}

	bool operator!=(const Vector3f& other)
	{
		return !(x == other.x && y == other.y && z == other.z);
	}

	Vector3f operator=(const Vector3f& other)
	{
		x = other.x;
		y = other.y;
		z = other.z;
		return *this;
	}

	operator bool() const
	{
		return(IsValid());
	}

	const float& operator[](const int i) const
	{
		return i == 0 ? this->x : (i == 1 ? this->y : this->z);
	}

	float& operator[](const int i)
	{
		return i == 0 ? this->x : (i == 1 ? this->y : this->z);
	}
};

//float * Vector3f
inline Vector3f operator*(float scaler, const Vector3f& v)
{
	return Vector3f(scaler * v.x, scaler * v.y, scaler * v.z);
}

//float / Vector3f
inline Vector3f operator/(float scaler, const Vector3f& v)
{
	return Vector3f(scaler / v.x, scaler / v.y, scaler / v.z);
}

//adds the two vectors together
inline Vector3f operator+(const Vector3f& v1, const Vector3f& v2)
{
	return Vector3f(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

inline Vector3f operator-(const Vector3f& v1, const Vector3f& v2)
{
	return Vector3f(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

inline std::ostream& operator<<(std::ostream& os, const Vector3f& v)
{
	return os << v.to_string();
}
#endif // !_VECTOR3F_H