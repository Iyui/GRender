#pragma once
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
    };

    typedef Vec3<float> Vec3f;

    Vec3f a;
    Vec3f b;
};

