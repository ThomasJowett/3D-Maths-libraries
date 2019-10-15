#include "Plane.h"

//Plane with normal pointing up the Y
Plane::Plane() 
{
	a = 0.0f;
	b = 1.0f;
	c = 0.0f;
	d = 0.0f;
}

Plane::Plane(const Plane & p) :a(p.a), b(p.b), c(p.c), d(p.d) {}

Plane::Plane(float a, float b, float c, float d) : a(a), b(b), c(c), d(d) {}

Plane::Plane(const Vector3f & normal, float distance)
{
	a = normal.x;
	b = normal.y;
	c = normal.z;
	d = distance;
}

void Plane::Normalize()
{
	float distance = sqrt(a*a + b * b + c * c);
	a = a / distance;
	b = b / distance;
	c = c / distance;
	d = d / distance;
}
