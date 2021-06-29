#pragma once
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
    };

    typedef Vec3<float> Vec3f;

    Vec3f a;
    Vec3f b;
};

