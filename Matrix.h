#ifndef _MATRIX_H
#define _MATRIX_H

#include "Vector3f.h"
#include "Vector2f.h"
#include "Plane.h"
#include "Quaternion.h"
#include <string>

class Matrix4x4
{
public:
	float m[4][4];

	Matrix4x4();
	Matrix4x4(
		float m00, float m01, float m02, float m03,
		float m10, float m11, float m12, float m13,
		float m20, float m21, float m22, float m23,
		float m30, float m31, float m32, float m33);
	~Matrix4x4();

	static Matrix4x4 Translate(Vector3f translation);
	static Matrix4x4 Scale(Vector3f scale);
	static Matrix4x4 RotateX(float angle);
	static Matrix4x4 RotateY(float angle);
	static Matrix4x4 RotateZ(float angle);
	static Matrix4x4 Rotate(float r, float i, float j, float k);
	static Matrix4x4 Rotate(Quaternion rotation);
	static Matrix4x4 Perspective(float fovY, float aspectRatio, float nearDepth, float farDepth);
	static Matrix4x4 Orthographic(float left, float right, float bottom, float top, float nearDepth, float farDepth);
	static Matrix4x4 LookAt(Vector3f eyePosition, Vector3f lookAtPosition, Vector3f up);

	static void FrustumPlanes(Plane planes[6], Matrix4x4 M);

	static Vector2f MulVec2(Matrix4x4 matrix, Vector2f vector);
	static Vector3f MulVec3(Matrix4x4 matrix, Vector3f vector);

	Vector2f ToVector2f() const;
	Vector3f ToVector3f() const;

	Vector3f ExtractTranslation() const;
	Vector3f ExtractScale()	const;
	Quaternion ExtractRotation() const;
	float ExtractRotationX() const;
	float ExtractRotationY() const;
	float ExtractRotationZ() const;

	std::string to_string()const;

	Matrix4x4 operator*(Matrix4x4 other);
	float operator()(size_t row, size_t column);
};

inline std::ostream& operator<<(std::ostream& os, Matrix4x4& m)
{
	return os << m.to_string();
}

#endif // !_MATRIX_H