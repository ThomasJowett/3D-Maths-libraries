#ifndef _PLANE_H
#define _PLANE_H

#include "Line3D.h"

class Plane
{
public:
	union
	{
		struct
		{
			Vector3f normal;
			float d;
		}nd;

		struct
		{
			float a;	//Normal X
			float b;	//Normal Y
			float c;	//Normal Z
			float d;	//Distance from origin
		};
	};

	Plane()
	{
		//Plane with normal pointing up the Y
		a = 0.0f;
		b = 1.0f;
		c = 0.0f;
		d = 0.0f;
	}

	Plane(const Plane& p) :a(p.a), b(p.b), c(p.c), d(p.d) {}
	Plane(float a, float b, float c, float d) : a(a), b(b), c(c), d(d) {}
	Plane(const Vector3f& normal, float distance)
	{
		a = normal.x;
		b = normal.y;
		c = normal.z;
		d = distance;
	}

	// construct from a point and a normal
	Plane(const Vector3f& point, Vector3f& normal)
	{
		Vector3f normalizedNormal = normal.GetNormalized();
		a = normalizedNormal.x;
		b = normalizedNormal.y;
		c = normalizedNormal.z;
		d = -Vector3f::Dot(point, normalizedNormal);
	}

	// construct from three points
	Plane(const Vector3f& p0, const Vector3f& p1, const Vector3f& p2)
	{
		Vector3f normal = Vector3f::Cross(p1 - p0, p2 - p0);
		normal.Normalize();

		a = normal.x;
		b = normal.y;
		c = normal.z;
		d = -Vector3f::Dot(p0, normal);
	}

	float UnsignedDistance(const Vector3f& point)const
	{
		return std::abs(a * point.x + b * point.y + c * point.z + d);
	}

	float SignedDistance(const Vector3f& point)const
	{
		return (a * point.x + b * point.y + c * point.z + d);
	}

	Vector3f ClosestPoint(const Vector3f& point)
	{
		return (point - nd.normal * SignedDistance(point));
	}

	void Normalize()
	{
		float distance = sqrt(a * a + b * b + c * c);
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

	std::string to_string()
	{
		return "a: " + std::to_string(a) + " b: " + std::to_string(b) + " c: " + std::to_string(c) + " d: " + std::to_string(d);
	}

	// Plane intersection
	static bool PlanePlaneIntersection(const Plane& p1, const Plane& p2, Line3D& intersectionLine)
	{
		float denominator = p1.a * p2.b - p1.b * p2.a;

		if (denominator == 0.0f)
		{
			return false;
		}

		intersectionLine.p = Vector3f((p2.d * p1.b - p1.d * p2.b) / denominator, (p1.d * p2.a - p2.d * p1.a) / denominator, 0.0f);

		intersectionLine.d = Vector3f::Cross(p1.nd.normal, p2.nd.normal);

		if (intersectionLine.d.Magnitude() == 0.0f)
		{
			return false;
		}

		intersectionLine.d.Normalize();

		return false;
	}

	static bool PlaneLineIntersection(const Plane& p, const Line3D& line, const Vector3f& origin, Vector3f& intersectionPoint)
	{
		float dot = Vector3f::Dot(p.nd.normal, origin);

		if (Vector3f::Dot(p.nd.normal, line.d) == 0.0f)
		{
			return false; //No intersection, the line is parallel to the plane
		}
		float x = (dot - Vector3f::Dot(p.nd.normal, line.p)) / Vector3f::Dot(p.nd.normal, line.d);
		intersectionPoint = line.p + line.d.GetNormalized() * x;
		return true;
	}

	template<typename Archive>
	void serialize(Archive& archive)
	{
		archive(a, b, c, d);
	}
};

inline std::ostream& operator<<(std::ostream& os, Plane& p)
{
	return os << p.to_string();
}
#endif // !_PLANE_H
