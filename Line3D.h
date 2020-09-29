#include "Vector3f.h"
#include "Plane.h"

class Line3D
{
public:
	Vector3f p;
	Vector3f d;

	Line3D() {}
	Line3D(const Vector3f &p, const Vector3f &d) :p(p), d(d) {}

	//Member Functions --------------------------------------------------------------------------------

	float DistanceToPoint(const Vector3f &point) const
	{
		float t0 = Vector3f::Dot(d, point - p) / Vector3f::Dot(d, d);
		float distanceLine = Vector3f::Distance(point, p + t0 * d);
		return distanceLine;
	}

	//Static-------------------------------------------------------------------------------------------

	static float Distance(const Line3D &L1, const Line3D &L2)
	{
		Vector3f Cross = Vector3f::Cross(L1.d, L2.d);
		return abs(Vector3f::Dot(L2.p - L1.p, Cross)) / Cross.Magnitude();
	}

	static float SqrDistance(const Line3D &L1, const Line3D &L2)
	{
		Vector3f cross = Vector3f::Cross(L1.d, L2.d);
		float dot = Vector3f::Dot(L2.p - L1.p, cross);
		return (dot * dot) / cross.SqrMagnitude();
	}

	template<typename Archive>
	void serialize(Archive& archive)
	{
		archive(cereal::make_nvp("Point", p));
		archive(cereal::make_nvp("Distance", d));
	}
};