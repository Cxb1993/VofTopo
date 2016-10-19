#ifndef VECTORS_CUDA_H
#define VECTORS_CUDA_H

#include <cstdlib>

struct int2
{
    int x, y;
};

struct int3
{
    int x, y, z;
};

struct float3
{
    float x, y, z;
};

struct float4
{
    float x, y, z, w;
};

struct double3
{
    double x, y, z;
};

inline int2 make_int2(int x, int y)
{
  int2 t;
  t.x = x; t.y = y;
  return t;
}

inline float3 make_float3(float x, float y, float z)
{
  float3 t;
  t.x = x; t.y = y; t.z = z;
  return t;
}

inline float3 make_float3(float s)
{
  return make_float3(s, s, s);
}

inline float3 make_float3(float4 a)
{
    return make_float3(a.x, a.y, a.z);
}

inline float4 make_float4(float x, float y, float z, float w)
{
  float4 t;
  t.x = x; t.y = y; t.z = z; t.w = w;
  return t;
}

inline float4 make_float4(float3 a)
{
    return make_float4(a.x, a.y, a.z, 0.0f);
}

inline float4 make_float4(float3 a, float w)
{
    return make_float4(a.x, a.y, a.z, w);
}

inline float3 operator+(float3 a, float3 b)
{
    return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline float4 operator+(float4 a, float4 b)
{
    return make_float4(a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w);
}

inline void operator+=(float3 &a, float3 b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
}

inline void operator+=(float4 &a, float4 b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    a.w += b.w;
}

inline float3 operator-(float3 a, float3 b)
{
    return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline float4 operator-(float4 a, float4 b)
{
    return make_float4(a.x - b.x, a.y - b.y, a.z - b.z,  a.w - b.w);
}

inline void operator-=(float4 &a, float4 b)
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    a.w -= b.w;
}

inline float3 operator*(float3 a, float b)
{
    return make_float3(a.x * b, a.y * b, a.z * b);
}

inline float3 operator*(float b, float3 a)
{
    return make_float3(b * a.x, b * a.y, b * a.z);
}

inline float3 operator*(float3 a, float3 b)
{
    return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}

inline float4 operator*(float4 a, float b)
{
    return make_float4(a.x * b, a.y * b, a.z * b,  a.w * b);
}

inline float4 operator*(float b, float4 a)
{
    return make_float4(b * a.x, b * a.y, b * a.z, b * a.w);
}

inline void operator*=(float3 &a, float3 b)
{
    a.x *= b.x;
    a.y *= b.y;
    a.z *= b.z;
}

inline float3 operator/(float3 a, float b)
{
    return make_float3(a.x / b, a.y / b, a.z / b);
}

inline float4 operator/(float4 a, float b)
{
    return make_float4(a.x / b, a.y / b, a.z / b,  a.w / b);
}

inline float3 operator/(float3 a, float3 b)
{
    return make_float3(a.x / b.x, a.y / b.y, a.z / b.z);
}

inline float rsqrtf(float x)
{
    return 1.0f / sqrtf(x);
}

inline float dot(float3 a, float3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline float length(float3 v)
{
    return sqrtf(dot(v, v));
}

inline float3 normalize(float3 v)
{
    float invLen = rsqrtf(dot(v, v));
    return v * invLen;
}

inline float3 cross(float3 a, float3 b)
{
    return make_float3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

inline float3 lerp(float3 a, float3 b, float t)
{
    return a + t*(b-a);
}

inline float4 lerp(float4 a, float4 b, float t)
{
    return a + t*(b-a);
}

inline float3 minf3(float3 a, float3 b)
{
    return make_float3(std::min(a.x,b.x), std::min(a.y,b.y), std::min(a.z,b.z));
}

inline float3 maxf3(float3 a, float3 b)
{
    return make_float3(std::max(a.x,b.x), std::max(a.y,b.y), std::max(a.z,b.z));
}

#endif//VECTORS_CUDA_H
