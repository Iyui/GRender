#pragma once
#include <iostream>
#include <algorithm>
class Math
{
    template<typename T>//����ģ��
    class Vec3
    {
    public:
        // 3�������������
        Vec3() : x(T(0)), y(T(0)), z(T(0)) {}//��ʼ��ʱ�����ඨ���ж����Ա��˳��ֱ���ø��Զ���Ĺ��캯������ִ���Լ��Ĺ��캯��
        Vec3(const T& xx) : x(xx), y(xx), z(xx) {}
        Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
        T x, y, z;

        //��������
        T length()
        {
            return sqrt(x * x + y * y + z * z);
        }

        //������һ��
        Vec3<T>& normalize()
        {
            T len = length();
            if (len > 0) {
                T invLen = 1 / len;
                x *= invLen, y *= invLen, z *= invLen;
            }
            return *this;
        }

        //��������Կ�����һ����������һ��������ͶӰ��
        T dot(const Vec3<T>& v) const
        {
            return x * v.x + y * v.y + z * v.z;
        }

        //��������������Ĳ���������������ʹ�ֱ��
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
        // �õ�λ�����ϵ����ʼ�������ϵ��
        T m[4][4] = { {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1} };

        //�����任
        void multDirMatrix(const Vec3<T>& src, Vec3<T>& dst) const
        {
            dst.x = src.x * m[0][0] + src.y * m[1][0] + src.z * m[2][0];
            dst.y = src.x * m[0][1] + src.y * m[1][1] + src.z * m[2][1];
            dst.z = src.x * m[0][2] + src.y * m[1][2] + src.z * m[2][2];
        }

        //�������
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

        //����ת�ã������к��У�
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
    //�����������ѿ�������
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

    //���� atan ���ط�Χ�ڵ�ֵ[ - ��: ��]. �ڷ�Χ������ӳ��Ϊ[ 0 : 2 ��]
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

