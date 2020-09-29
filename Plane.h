#ifndef _PLANE_H
#define _PLANE_H

#include "Vector3f.h"
#include "Line3D.h"
#include <vector>

class Plane
{
public:
	float a, b, c, d;

	Plane()
	{
		//Plane with normal pointing up the Y
		a = 0.0f;
		b = 1.0f;
		c = 0.0f;
		d = 0.0f;
	}

	Plane(const Plane &p) :a(p.a), b(p.b), c(p.c), d(p.d) {}
	Plane(float a, float b, float c, float d) : a(a), b(b), c(c), d(d) {}
	Plane(const Vector3f &normal, float distance)
	{
		a = normal.x;
		b = normal.y;
		c = normal.z;
		d = distance;
	}

	// construct from a point and a normal
	Plane(const Vector3f &point, Vector3f &normal)
	{
		Vector3f normalizedNormal = normal.GetNormalized();
		a = normalizedNormal.x;
		b = normalizedNormal.y;
		c = normalizedNormal.z;
		d = -Vector3f::Dot(point, normalizedNormal);
	}

	// construct from three points
	Plane(const Vector3f & p0, const Vector3f &p1, const Vector3f &p2)
	{
		Vector3f normal = Vector3f::Cross(p1 - p0, p2 - p0);
		normal.Normalize();

		a = normal.x;
		b = normal.y;
		c = normal.z;
		d = -Vector3f::Dot(p0, normal);
	}

	float UnsignedDistance(const Vector3f &point)const
	{
		return abs(a * point.x + b * point.y + c * point.z + d);
	}

	float SignedDistance(const Vector3f &point)const
	{
		return (a * point.x + b * point.y + c * point.z + d);
	}

	Vector3f ClosestPoint(const Vector3f &point)
	{
		return (point - Normal() * SignedDistance(point));
	}

	void Normalize()
	{
		float distance = sqrt(a*a + b * b + c * c);
		a = a / distance;
		b = b / distance;
		c = c / distance;
		d = d / distance;
	}

	void Flip()
	{
		a = -a;
		b = -b;
		c = -c;
		d = -d;
	}

	Vector3f Normal()const
	{
		return Vector3f(a, b, c);
	}

	std::string to_string()
	{
		return "a: " + std::to_string(a) + " b: " + std::to_string(b) + " c: " + std::to_string(c) + " d: " + std::to_string(d);
	}

	// Plane intersection
	static bool PlanePlaneIntersection(const Plane &p1, const Plane &p2, Line3D& intersectionLine)
	{
		return false;
	}

	template<typename Archive>
	void serialize(Archive& archive)
	{
		archive(cereal::make_nvp("a", a));
		archive(cereal::make_nvp("b", b));
		archive(cereal::make_nvp("c", c));
		archive(cereal::make_nvp("d", d));
	}
};

inline std::ostream& operator<<(std::ostream& os, Plane& p)
{
	return os << p.to_string();
}
#endif // !_PLANE_H


