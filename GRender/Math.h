#pragma once
#include <iostream>
#include <algorithm>
class Math
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

        //向量长度
        T length()
        {
            return sqrt(x * x + y * y + z * z);
        }

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
        Vec3<T> operator - (const Vec3<T> &v) const
        {
            return Vec3<T>(x - v.x, y - v.y, z - v.z);
        }
        Vec3<T> operator * (const T& r) const
        {
            return Vec3<T>(x * r, y * r, z * r);
        }

        

        
    };

    typedef Vec3<float> Vec3f;

    Vec3f a;
    Vec3f b;

    template<typename T>
    class Matrix44
    {
    public:
        Matrix44() {}
        const T* operator [] (uint8_t i) const { return m[i]; }
        T* operator [] (uint8_t i) { return m[i]; }
        // 用单位矩阵的系数初始化矩阵的系数
        T m[4][4] = { {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1} };

        //向量变换
        void multDirMatrix(const Vec3<T>& src, Vec3<T>& dst) const
        {
            dst.x = src.x * m[0][0] + src.y * m[1][0] + src.z * m[2][0];
            dst.y = src.x * m[0][1] + src.y * m[1][1] + src.z * m[2][1];
            dst.z = src.x * m[0][2] + src.y * m[1][2] + src.z * m[2][2];
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

    template<typename T> inline T sinTheta(const Vec3<T>& w){return sqrt(sinTheta2(w));}

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

