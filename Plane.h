#ifndef _PLANE_H
#define _PLANE_H

#include "Vector3f.h"
#include <vector>

class Plane
{
public:
	float a, b, c, d;

	Plane();
	Plane(const Plane &p);
	Plane(float a, float b, float c, float d);
	Plane(const Vector3f &normal, float d);

	float UnsignedDistance(const Vector3f &point)const;
	float SignedDistance(const Vector3f &point)const;

	void Normalize();

	__forceinline Vector3f Normal()const
	{
		return Vector3f(a, b, c);
	}

};
#endif // !_PLANE_H


