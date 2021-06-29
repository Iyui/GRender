#include "GeoMath.h"
#include <iostream>
#include <algorithm>
#include <iomanip>
class GeoMath
{
	template<typename T>//定义模板
	class Vec3
	{
	public:
		// 3种最基本的向量
		Vec3() : x(T(0)), y(T(0)), z(T(0)) {}//初始化时按照类定义中对象成员的顺序分别调用各自对象的构造函数，再执行自己的构造函数
		Vec3(const T& xx) : x(xx), y(xx), z(xx) {}
		Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
		T x, y, z;


		T norm() const
		{
			return x * x + y * y + z * z;
		}

		//向量长度
		T length()
		{
			return sqrt(norm());
		}

		const T& operator [] (uint8_t i) const { return (&x)[i]; }
		T& operator [] (uint8_t i) { return (&x)[i]; }

		//向量归一化
		Vec3<T>& normalize()
		{
			T len = length();
			if (len > 0) {
				T invLen = 1 / len;
				x *= invLen, y *= invLen, z *= invLen;
			}
			return *this;
		}

		//点积（可以看作是一个向量到另一个向量的投影）
		T dot(const Vec3<T>& v) const
		{
			return x * v.x + y * v.y + z * v.z;
		}

		//叉积（两个向量的叉积与这两个向量和垂直）
		Vec3<T> cross(const Vec3<T>& v) const
		{
			return Vec3<T>(
				y * v.z - z * v.y,
				z * v.x - x * v.z,
				x * v.y - y * v.x);
		}

		Vec3<T> operator + (const Vec3<T>& v) const
		{
			return Vec3<T>(x + v.x, y + v.y, z + v.z);
		}
		Vec3<T> operator - (const Vec3<T>& v) const
		{
			return Vec3<T>(x - v.x, y - v.y, z - v.z);
		}
		Vec3<T> operator * (const T& r) const
		{
			return Vec3<T>(x * r, y * r, z * r);
		}
		friend std::ostream& operator << (std::ostream& s, const Vec3<T>& v)
		{
			return s << '(' << v[0] << ' ' << v[1] << ' ' << v[2] << ')';
		}

	};
public:typedef Vec3<float> Vec3f;
public:typedef Vec3<int> Vec3i;
	  //Vec3f a;

	  template<typename T>
	  class Matrix44
	  {
	  public:
		  Matrix44() {}
		
		  // 用单位矩阵的系数初始化矩阵的系数
		  T m[4][4] = { {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1} };

		  Matrix44(T a, T b, T c, T d, T e, T f, T g, T h,
			  T i, T j, T k, T l, T m, T n, T o, T p)
		  {
			  m[0][0] = a;
			  m[0][1] = b;
			  m[0][2] = c;
			  m[0][3] = d;
			  m[1][0] = e;
			  m[1][1] = f;
			  m[1][2] = g;
			  m[1][3] = h;
			  m[2][0] = i;
			  m[2][1] = j;
			  m[2][2] = k;
			  m[2][3] = l;
			  m[3][0] = m;
			  m[3][1] = n;
			  m[3][2] = o;
			  m[3][3] = p;
		  }
		
		  //矩阵相乘
		  Matrix44 operator * (const Matrix44& rhs) const
		  {
			  Matrix44 mult;
			  for (uint8_t i = 0; i < 4; ++i) {
				  for (uint8_t j = 0; j < 4; ++j) {
					  mult[i][j] = m[i][0] * rhs[0][j] +
						  m[i][1] * rhs[1][j] +
						  m[i][2] * rhs[2][j] +
						  m[i][3] * rhs[3][j];
				  }
			  }
			  return mult;
		  }

		  //矩阵转置（交换行和列）
		  Matrix44 transpose() const
		  {
			  Matrix44 transpMat;
			  for (uint8_t i = 0; i < 4; ++i) {
				  for (uint8_t j = 0; j < 4; ++j) {
					  transpMat[i][j] = m[j][i];
				  }
			  }
			  return transpMat;
		  }

		  Matrix44 inverse()
		  {
			  int i, j, k;
			  Matrix44 s;
			  Matrix44 t(*this);

			  // Forward elimination
			  for (i = 0; i < 3; i++) {
				  int pivot = i;

				  T pivotsize = t[i][i];

				  if (pivotsize < 0)
					  pivotsize = -pivotsize;

				  for (j = i + 1; j < 4; j++) {
					  T tmp = t[j][i];

					  if (tmp < 0)
						  tmp = -tmp;

					  if (tmp > pivotsize) {
						  pivot = j;
						  pivotsize = tmp;
					  }
				  }

				  if (pivotsize == 0) {
					  // Cannot invert singular matrix
					  return Matrix44();
				  }

				  if (pivot != i) {
					  for (j = 0; j < 4; j++) {
						  T tmp;

						  tmp = t[i][j];
						  t[i][j] = t[pivot][j];
						  t[pivot][j] = tmp;

						  tmp = s[i][j];
						  s[i][j] = s[pivot][j];
						  s[pivot][j] = tmp;
					  }
				  }

				  for (j = i + 1; j < 4; j++) {
					  T f = t[j][i] / t[i][i];

					  for (k = 0; k < 4; k++) {
						  t[j][k] -= f * t[i][k];
						  s[j][k] -= f * s[i][k];
					  }
				  }
			  }

			  // Backward substitution
			  for (i = 3; i >= 0; --i) {
				  T f;

				  if ((f = t[i][i]) == 0) {
					  // Cannot invert singular matrix
					  return Matrix44();
				  }

				  for (j = 0; j < 4; j++) {
					  t[i][j] /= f;
					  s[i][j] /= f;
				  }

				  for (j = 0; j < i; j++) {
					  f = t[j][i];

					  for (k = 0; k < 4; k++) {
						  t[j][k] -= f * t[i][k];
						  s[j][k] -= f * s[i][k];
					  }
				  }
			  }

			  return s;
		  }

		  // 矩阵求逆
		  const Matrix44<T>& invert()
		  {
			  *this = inverse();
			  return *this;
		  }

		  friend std::ostream& operator << (std::ostream& s, const Matrix44& m)
		  {
			  std::ios_base::fmtflags oldFlags = s.flags();
			  int width = 12; //显示数字的总数
			  s.precision(5); //返回当前的浮点数精度值
			  s.setf(std::ios_base::fixed);

			  s << "(" << std::setw(width) << m[0][0] <<
				  " " << std::setw(width) << m[0][1] <<
				  " " << std::setw(width) << m[0][2] <<
				  " " << std::setw(width) << m[0][3] << "\n" <<

				  " " << std::setw(width) << m[1][0] <<
				  " " << std::setw(width) << m[1][1] <<
				  " " << std::setw(width) << m[1][2] <<
				  " " << std::setw(width) << m[1][3] << "\n" <<

				  " " << std::setw(width) << m[2][0] <<
				  " " << std::setw(width) << m[2][1] <<
				  " " << std::setw(width) << m[2][2] <<
				  " " << std::setw(width) << m[2][3] << "\n" <<

				  " " << std::setw(width) << m[3][0] <<
				  " " << std::setw(width) << m[3][1] <<
				  " " << std::setw(width) << m[3][2] <<
				  " " << std::setw(width) << m[3][3] << ")\n";

			  s.flags(oldFlags);
			  return s;
		  }

		  template<typename S>
		  void multVecMatrix(const Vec3<S>& src, Vec3<S>& dst) const
		  {
			  S a, b, c, w;

			  a = src[0] * m[0][0] + src[1] * m[1][0] + src[2] * m[2][0] + m[3][0];
			  b = src[0] * m[0][1] + src[1] * m[1][1] + src[2] * m[2][1] + m[3][1];
			  c = src[0] * m[0][2] + src[1] * m[1][2] + src[2] * m[2][2] + m[3][2];
			  w = src[0] * m[0][3] + src[1] * m[1][3] + src[2] * m[2][3] + m[3][3];

			  dst.x = a / w;
			  dst.y = b / w;
			  dst.z = c / w;
		  }

		  template<typename S>
		  void multDirMatrix(const Vec3<S>& src, Vec3<S>& dst) const
		  {
			  S a, b, c;

			  a = src[0] * m[0][0] + src[1] * m[1][0] + src[2] * m[2][0];
			  b = src[0] * m[0][1] + src[1] * m[1][1] + src[2] * m[2][1];
			  c = src[0] * m[0][2] + src[1] * m[1][2] + src[2] * m[2][2];

			  dst.x = a;
			  dst.y = b;
			  dst.z = c;
		  }
	  };



	  typedef Matrix44<float> Matrix44f;
	  template<class T>
	  T clamp(T x, T min, T max)
	  {
		  if (x > max)
			  return max;
		  if (x < min)
			  return min;
		  return x;
	  }
	  //从球坐标计算笛卡尔坐标
	  template<typename T>
	  Vec3<T>sphericalToCartesian(const T& theta, const T& phi)
	  {
		  return Vec3<T>(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
	  };

	  template<typename T>
	  inline T sphericalTheta(const Vec3<T>& v)
	  {
		  return acos(clamp<T>(v[2], -1, 1));
	  }

	  //函数 atan 返回范围内的值[ - π: π]. 在范围内重新映射为[ 0 : 2 π]
	  template<typename T>
	  inline T sphericalPhi(const Vec3<T>& v)
	  {
		  T p = atan2(v[1], v[0]);
		  return (p < 0) ? p + 2 * M_PI : p;
	  }

	  template<typename T> inline T cosTheta(const Vec3<T>& w) { return w[2]; }

	  template<typename T>
	  inline T sinTheta2(const Vec3<T>& w)
	  {
		  return std::max(T(0), 1 - cosTheta(w) * cosTheta(w));
	  }

	  template<typename T> inline T sinTheta(const Vec3<T>& w) { return sqrt(sinTheta2(w)); }

	  template<typename T>
	  inline T cosPhi(const Vec3<T>& w)
	  {
		  T sintheta = sinTheta(w);
		  if (sintheta == 0) return 1;
		  return clamp<T>(w[0] / sintheta, -1, 1);
	  }

	  template<typename T>
	  inline T sinPhi(const Vec3<T>& w)
	  {
		  T sintheta = sinTheta(w);
		  if (sintheta == 0) return 0;
		  return clamp<T>(w[1] / sintheta, -1, 1);
	  }
};