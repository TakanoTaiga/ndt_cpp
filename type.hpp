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

#ifndef NDTCPP_TYPE_H_
#define NDTCPP_TYPE_H_

namespace ndtcpp
{
    struct point2{
        float x, y;
    };
    struct point3{
        float x, y, z;
    };

    struct mat2x2{
        float a, b;
        float c, d;
    };

    struct mat3x3{
        float a, b, c;
        float d, e, f;
        float g, h, i;
    };
} // namespace ndtcpp

#endif // NDTCPP_TYPE_H_
