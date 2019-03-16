#ifndef _PLANE_H
#define _PLANE_H

#include "Vector.h"
#include <vector>

class Plane
{
public:
	float a, b, c, d;

	Plane();
	Plane(const Plane &p);
	Plane(float a, float b, float c, float d);
	Plane(const Vector3D &normal, float d);

	float UnsignedDistance(const Vector3D &point)const;
	float SignedDistance(const Vector3D &point)const;

	void Normalize();

	__forceinline Vector3D Normal()const
	{
		return Vector3D(a, b, c);
	}

};
#endif // !_PLANE_H


