#ifndef _QUATERNION_H
#define _QUATERNION_H

#include <float.h>
#include "Vector3f.h"
#include <string>

#define PI			3.14159265358979323846

class Quaternion
{
public:
	float r;	//Real
	float i;	//First Complex
	float j;	//Second Complex
	float k;	//Thrid Complex

	//Constructors ----------------------------------------------------------------------------------
	Quaternion() : r(1), i(0), j(0), k(0) {}
	Quaternion(const float r, const float i, const float j, const float k) : r(r), i(i), j(j), k(k) {}
	Quaternion(const Quaternion& q) : r(q.r), i(q.i), j(q.j), k(q.k) {}

	//from a vector of 3 euler angles in radians
	Quaternion(const Vector3f angles)
	{
		float cos_x_2 = cosf(0.5f*angles.x);
		float cos_y_2 = cosf(0.5f*angles.y);
		float cos_z_2 = cosf(0.5f*angles.z);

		float sin_x_2 = sinf(0.5f*angles.x);
		float sin_y_2 = sinf(0.5f*angles.y);
		float sin_z_2 = sinf(0.5f*angles.z);

		// and now compute quaternion
		r = cos_z_2*cos_y_2*cos_x_2 + sin_z_2*sin_y_2*sin_x_2;
		i = cos_z_2*cos_y_2*sin_x_2 - sin_z_2*sin_y_2*cos_x_2;
		j = cos_z_2*sin_y_2*cos_x_2 + sin_z_2*cos_y_2*sin_x_2;
		k = sin_z_2*cos_y_2*cos_x_2 - cos_z_2*sin_y_2*sin_x_2;
	}

	//from 3 euler angles in radians
	Quaternion(const float theta_Roll, float theta_Pitch, float theta_Yaw)
	{
		(*this) = Quaternion(Vector3f(theta_Roll, theta_Pitch, theta_Yaw));
	}

	//Destructor---------------------------------------------------------------------------------------
	~Quaternion() = default;

	//Member Functions --------------------------------------------------------------------------------
	float GetSqrMagnitude() const
	{
		return r*r + i*i + j*j + k*k;
	}

	float GetMagnitude() const
	{
		return sqrtf(GetSqrMagnitude());
	}

	//normalize this quaternion
	void Normalize()
	{
		float d = GetSqrMagnitude();

		if (d < FLT_EPSILON)
		{
			r = 1;
		}
		else
		{
			d = static_cast<float>(1.0) / sqrtf(d);
			r *= d;
			i *= d;
			j *= d;
			k *= d;
		}
	}

	//return this quaternion normalized
	Quaternion Normalized() const
	{
		Quaternion quaternion = *this;
		quaternion.Normalize();
		return quaternion;
	}
	Quaternion Conjugate() const
	{
		return Quaternion(r, -i, -j, -k);
	}

	Quaternion Scale(float scaler) const
	{
		return Quaternion(r*scaler, i*scaler, j*scaler, k*scaler);
	}

	Quaternion Inverse()
	{
		return Conjugate().Scale(1 / GetSqrMagnitude());
	}

	Quaternion UnitQuaternion()
	{
		return (*this).Scale(1 / (*this).GetMagnitude());
	}

	//returns the euler angles of the quaternion
	Vector3f EulerAngles(bool homogenous = true)const
	{
		Vector3f euler;

		if (homogenous)
		{
			euler.x = atan2f(2.f * (i*j + k * i), (i*i) - (j*j) - (k*k) + (r*r));
			euler.y = asinf(-2.f * (i*k - j * r));
			euler.z = atan2f(2.f * (j*k + i * r), -(i*i) - (j*j) + (k*k) + (r*r));
		}
		else
		{
			euler.x = atan2f(2.f * (k*j + i * r), 1 - 2 * ((i*i) + (j*j)));
			euler.y = asinf(-2.f * (i*k - j * r));
			euler.z = atan2f(2.f * (i*j + k * r), 1 - 2 * ((j*j) + (k*k)));
		}
		return euler;
	}

	void AxisAngle(float& angle, Vector3f& axis) const
	{
		Quaternion q = *this;
		if (r > 1)
			q.Normalize();
		angle = 2 * (float)acos(q.r);
		float s = sqrtf(1 - q.r*q.r);

		if (s < 0.001)
		{
			axis.x = q.i;
			axis.y = q.j;
			axis.z = q.k;
		}
		else
		{
			axis.x = q.i / s;
			axis.y = q.j / s;
			axis.z = q.k / s;
		}

		//convert to degrees
		angle *= (float)(180 / PI);
	}

	//rotates this quaternion by a vector
	void RotateByVector(const Vector3f& vector)
	{
		Quaternion q(0, vector.x, vector.y, vector.z);
		*this = *this * q;
	}

	//rotates a vector by this quaternion
	Vector3f RotateVectorByQuaternion(Vector3f& vector)
	{
		Quaternion V(0, vector.x, vector.y, vector.z);
		V = (*this * V * this->Conjugate());
		vector.x = V.i;
		vector.z = V.j;
		vector.z = V.k;
		return vector;
	}

	void AddScaledVector(const Vector3f& vector, float scale)
	{
		Quaternion q(0, vector.x * scale, vector.y * scale, vector.z * scale);

		q = *this * q;
		r += q.r * 0.5f;
		i += q.i * 0.5f;
		j += q.j * 0.5f;
		k += q.k * 0.5f;
	}

	//converts quaternion to a formatted string
	std::string to_string()
	{
		return "r: " + std::to_string(r) + "i: " + std::to_string(i) + "j: " + std::to_string(j) + "k: " + std::to_string(k);
	}

	//Static -----------------------------------------------------------------------------------
	float static Dot(Quaternion q1, Quaternion q2)
	{
		return q1.r * q2.r + q1.i * q2.i + q1.j * q2.j + q1.k * q2.k;
	}

	//Spherically interpolates between a and b
	Quaternion static Slerp(Quaternion a, Quaternion b, double t)
	{
		Quaternion result = Quaternion();

		//calculate angle between them
		float cosHalfTheta = Dot(a, b);

		//if a==b or a==-b then theta = 0 and we can return a
		if (abs(cosHalfTheta) >= 1.0)
		{
			result = a;
			return result;
		}

		//calculate temporary values
		float halfTheta = acos(cosHalfTheta);
		float sinHalfTheta = (float)sqrt(1.0 - cosHalfTheta * cosHalfTheta);

		//if theta = 180 degrees then result is not fully defined
		//we could rotate around any axis normal to a or b

		if (fabs(sinHalfTheta) < 0.001)
		{
			result.r = (a.r * 0.5f + b.r * 0.5f);
			result.i = (a.i * 0.5f + b.i * 0.5f);
			result.j = (a.j * 0.5f + b.j * 0.5f);
			result.k = (a.k * 0.5f + b.k * 0.5f);
			return result;
		}

		float ratioA = (float)sin((1 - t)*halfTheta) / sinHalfTheta;
		float ratioB = (float)sin(t *halfTheta) / sinHalfTheta;

		//calculate Quaternion
		result.r = (a.r * ratioA + b.r * ratioB);
		result.i = (a.i * ratioA + b.i * ratioB);
		result.j = (a.j * ratioA + b.j * ratioB);
		result.k = (a.k * ratioA + b.k * ratioB);
		return result;
	}

	//Normally interpolates between a and b
	Quaternion static Nlerp(Quaternion a, Quaternion b, float alpha)
	{
		Quaternion result = Quaternion();

		float dot = Dot(a, b);

		float oneMinusAlpha = 1.0f - alpha;

		if (dot < 0) {
			result.r = oneMinusAlpha * a.r + alpha * -b.r;
			result.i = oneMinusAlpha * a.i + alpha * -b.i;
			result.j = oneMinusAlpha * a.j + alpha * -b.j;
			result.k = oneMinusAlpha * a.k + alpha * -b.k;
		}
		else {
			result.r = oneMinusAlpha * a.r + alpha * b.r;
			result.i = oneMinusAlpha * a.i + alpha * b.i;
			result.j = oneMinusAlpha * a.j + alpha * b.j;
			result.k = oneMinusAlpha * a.k + alpha * b.k;
		}
		result.Normalize();
		return result;
	}

	//Operators---------------------------------------------------------------------------------
	Quaternion operator = (const Quaternion& q)
	{
		r = q.r;
		i = q.i;
		j = q.j;
		k = q.k;

		return (*this);
	}

	Quaternion operator + (const Quaternion& q)
	{
		return Quaternion(r + q.r, i + q.i, j + q.j, k + q.k);
	}

	Quaternion operator - (const Quaternion& q)
	{
		return Quaternion(r - q.r, i - q.i, j - q.j, k - q.k);
	}

	const Quaternion &operator +=(const Quaternion &q)
	{
		r += q.r;
		i += q.i;
		j += q.j;
		k += q.k;
		return *this;
	}
	
	Quaternion operator* (const Quaternion &multiplier)
	{
		Quaternion q;
		q.r = r*multiplier.r - i*multiplier.i - j*multiplier.j - k*multiplier.k;
		q.i = r*multiplier.i + i*multiplier.r + j*multiplier.k - k*multiplier.j;
		q.j = r*multiplier.j + j*multiplier.r + k*multiplier.i - i*multiplier.k;
		q.k = r*multiplier.k + k*multiplier.r + i*multiplier.j - j*multiplier.i;
		return q;
	}

	void operator *=(const Quaternion &multiplier)
	{
		Quaternion q = *this;
		r = q.r*multiplier.r - q.i*multiplier.i -
			q.j*multiplier.j - q.k*multiplier.k;
		i = q.r*multiplier.i + q.i*multiplier.r +
			q.j*multiplier.k - q.k*multiplier.j;
		j = q.r*multiplier.j + q.j*multiplier.r +
			q.k*multiplier.i - q.i*multiplier.k;
		k = q.r*multiplier.k + q.k*multiplier.r +
			q.i*multiplier.j - q.j*multiplier.i;
	}

	Vector3f operator*(const Vector3f& multiplier)
	{
		Vector3f V = multiplier;
		return RotateVectorByQuaternion(V);
	}
};

#endif // !_QUATERNION_H