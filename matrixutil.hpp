// Copyright 2024 Taiga Takano
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef NDTCPP_MATRIX_UTIL_HPP_
#define NDTCPP_MATRIX_UTIL_HPP_

#include "type.hpp"

namespace ndtcpp
{
    auto operator*(const mat3x3& mat1, const mat3x3& mat2)
    {
        mat3x3 result;
        result.a = mat1.a * mat2.a + mat1.b * mat2.d + mat1.c * mat2.g;
        result.b = mat1.a * mat2.b + mat1.b * mat2.e + mat1.c * mat2.h;
        result.c = mat1.a * mat2.c + mat1.b * mat2.f + mat1.c * mat2.i;

        result.d = mat1.d * mat2.a + mat1.e * mat2.d + mat1.f * mat2.g;
        result.e = mat1.d * mat2.b + mat1.e * mat2.e + mat1.f * mat2.h;
        result.f = mat1.d * mat2.c + mat1.e * mat2.f + mat1.f * mat2.i;

        result.g = mat1.g * mat2.a + mat1.h * mat2.d + mat1.i * mat2.g;
        result.h = mat1.g * mat2.b + mat1.h * mat2.e + mat1.i * mat2.h;
        result.i = mat1.g * mat2.c + mat1.h * mat2.f + mat1.i * mat2.i;
        return result;
    }

    auto operator*(const mat3x3& mat, const point3& vec)
    {
        point3 result;
        result.x = mat.a * vec.x + mat.b * vec.y + mat.c * vec.z;
        result.y = mat.d * vec.x + mat.e * vec.y + mat.f * vec.z;
        result.z = mat.g * vec.x + mat.h * vec.y + mat.i * vec.z;
        return result;
    }

    auto operator+=(mat3x3& mat1, const mat3x3& mat2)
    {
        mat1.a += mat2.a;
        mat1.b += mat2.b;
        mat1.c += mat2.c;

        mat1.d += mat2.d;
        mat1.e += mat2.e;
        mat1.f += mat2.f;

        mat1.g += mat2.g;
        mat1.h += mat2.h;
        mat1.i += mat2.i;
    }

    auto operator+=(ndtcpp::point3& point1, const ndtcpp::point3& point2){
        point1.x += point2.x;
        point1.y += point2.y;
        point1.z += point2.z;
    }
} // namespace ndtcpp

#endif // NDTCPP_MATRIX_UTIL_HPP_
