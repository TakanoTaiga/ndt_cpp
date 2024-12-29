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

#include <cstddef>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>

#include <chrono>

#include "flatkdtree.h"
#include "type.hpp"
#include "ndtcpputil.hpp"
#include "matrixutil.hpp"


template <std::size_t I>
struct kdtree::trait::access<ndtcpp::point2, I> {
    static auto get(const ndtcpp::point2 &p) -> float
    {
        return I == 0 ? p.x : p.y;
    }
};

template <>
struct kdtree::trait::dimension<ndtcpp::point2> {
    static constexpr std::size_t value = 2;
};

struct ndtpoint2 {
    ndtcpp::point2 mean;
    ndtcpp::mat2x2 cov;
};

template <std::size_t I>
struct kdtree::trait::access<ndtpoint2, I> {
    static auto get(const ndtpoint2 &p) -> float
    {
        return I == 0 ? p.mean.x : p.mean.y;
    }
};

template <>
struct kdtree::trait::dimension<ndtpoint2> {
    static constexpr std::size_t value = 2;
};



ndtcpp::mat3x3 makeTransformationMatrix(const float& tx, const float& ty, const float& theta) {
    ndtcpp::mat3x3 mat = {
        cosf(theta), sinf(theta) * -1.0f, tx,
        sinf(theta), cosf(theta)        , ty,
        0.0f, 0.0f, 1.0f
    };
    return mat;
}

void transformPointsZeroCopy(const ndtcpp::mat3x3& mat, std::vector<ndtcpp::point2>& points) {
    ndtcpp::point2 transformedPoint;

    for (auto& point : points) {
        transformedPoint.x = mat.a * point.x + mat.b * point.y + mat.c;
        transformedPoint.y = mat.d * point.x + mat.e * point.y + mat.f;
        point.x = transformedPoint.x;
        point.y = transformedPoint.y;
    }
}

float multtiplyPowPoint3(const ndtcpp::point3& vec){
    return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
}

ndtcpp::point3 solve3x3(const ndtcpp::mat3x3& m, const ndtcpp::point3& p) {
    float A[3][4] = {
        {m.a, m.b, m.c, p.x},
        {m.d, m.e, m.f, p.y},
        {m.g, m.h, m.i, p.z}
    };

    int n = 3;

    for (int i = 0; i < n; i++) {
        // Pivot選択
        float maxEl = std::abs(A[i][i]);
        int maxRow = i;
        for (int k = i+1; k < n; k++) {
            if (std::abs(A[k][i]) > maxEl) {
                maxEl = std::abs(A[k][i]);
                maxRow = k;
            }
        }

        // Pivotのある行を交換
        for (int k = i; k < n+1;k++) {
            float tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }

        // すべての行について消去を行う
        for (int k = i+1; k < n; k++) {
            float c = -A[k][i] / A[i][i];
            for (int j = i; j < n+1; j++) {
                if (i == j) {
                    A[k][j] = 0.0f;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
        }
    }

    // 解の計算 (後退代入)
    ndtcpp::point3 solution;
    solution.z = A[2][3] / A[2][2];
    solution.y = (A[1][3] - A[1][2] * solution.z) / A[1][1];
    solution.x = (A[0][3] - A[0][2] * solution.z - A[0][1] * solution.y) / A[0][0];

    return solution;
}

ndtcpp::mat3x3 expmap(const ndtcpp::point3& point){
    auto t = point.z;
    auto c = cosf(t);
    auto s = sinf(t);

    ndtcpp::mat2x2 R {
        c, s * -1.0f,
        s, c
    };

    ndtcpp::mat3x3 T {
        R.a, R.b, point.x,
        R.c, R.d, point.y,
        0.0f, 0.0f, 1.0f
    };

    return T;
}

ndtcpp::point2 transformPointCopy(const ndtcpp::mat3x3& mat, const ndtcpp::point2& point) {
    ndtcpp::point2 transformedPoint;

    transformedPoint.x = mat.a * point.x + mat.b * point.y + mat.c;
    transformedPoint.y = mat.d * point.x + mat.e * point.y + mat.f;

    return transformedPoint;
}

ndtcpp::mat3x3 inverse3x3Copy(const ndtcpp::mat3x3& mat){
    const auto a = 1.0 / (
        mat.a * mat.e * mat.i +
        mat.b * mat.f * mat.g +
        mat.c * mat.d * mat.h -
        mat.c * mat.e * mat.g -
        mat.b * mat.d * mat.i -
        mat.a * mat.f * mat.h
        );

    ndtcpp::mat3x3 inv_mat;
    inv_mat.a = mat.e * mat.i - mat.f * mat.h;
    inv_mat.b = mat.b * mat.i - mat.c * mat.h;
    inv_mat.c = mat.b * mat.f - mat.c * mat.e;

    inv_mat.d = mat.d * mat.i - mat.f * mat.g;
    inv_mat.e = mat.a * mat.i - mat.c * mat.g;
    inv_mat.f = mat.a * mat.f - mat.c * mat.d;

    inv_mat.g = mat.d * mat.h - mat.e * mat.g;
    inv_mat.h = mat.a * mat.h - mat.b * mat.g;
    inv_mat.i = mat.a * mat.e - mat.b * mat.d;


    inv_mat.a = inv_mat.a * a;
    inv_mat.b = inv_mat.b * a * -1.0f;
    inv_mat.c = inv_mat.c * a;

    inv_mat.d = inv_mat.d * a * -1.0f;
    inv_mat.e = inv_mat.e * a;
    inv_mat.f = inv_mat.f * a * -1.0f;

    inv_mat.g = inv_mat.g * a;
    inv_mat.h = inv_mat.h * a * -1.0f;
    inv_mat.i = inv_mat.i * a;

    return inv_mat;
}

ndtcpp::point2 skewd(const ndtcpp::point2& input_point){
    const ndtcpp::point2 skewd_point {
        input_point.y,
        input_point.x * -1.0f
    };
    return skewd_point;
}

ndtcpp::mat3x3 transpose(const ndtcpp::mat3x3& input_mat){
    const ndtcpp::mat3x3 transpose_mat{
        input_mat.a, input_mat.d, input_mat.g,
        input_mat.b, input_mat.e, input_mat.h,
        input_mat.c, input_mat.f, input_mat.i
    };
    return transpose_mat;
}

ndtcpp::point2 compute_mean(const std::vector<ndtcpp::point2>& points){
    ndtcpp::point2 mean;
    mean.x = 0.0f;
    mean.y = 0.0f;
    for(const auto& point : points){
        mean.x += point.x;
        mean.y += point.y;
    }
    mean.x = mean.x / (float)points.size();
    mean.y = mean.y / (float)points.size();
    return mean;
}

ndtcpp::mat2x2 compute_covariance(const std::vector<ndtcpp::point2>& points, const ndtcpp::point2& mean){
    auto point_size = points.size();
    auto vxx = 0.0f;
    auto vxy = 0.0f;
    auto vyy = 0.0f;

    for(const auto& point : points){
        const auto dx = point.x - mean.x;
        const auto dy = point.y - mean.y;
        vxx += dx * dx;
        vxy += dx * dy;
        vyy += dy * dy;
    }

    ndtcpp::mat2x2 cov;
    cov.a = vxx / point_size;
    cov.b = vxy / point_size;
    cov.c = vxy / point_size;
    cov.d = vyy / point_size;
    return cov;
}

void compute_ndt_points(std::vector<ndtcpp::point2>& points, std::vector<ndtpoint2> &results){
    auto N = 10;

    const auto point_size = points.size();

    kdtree::construct(points.begin(), points.end());
    std::vector<ndtcpp::point2> result_points(N);
    std::vector<float> result_distances(N);

    std::vector<ndtcpp::mat2x2> covs(point_size);
    results.resize(point_size);

    for(std::size_t i = 0; i < point_size; i++) {
        kdtree::search_knn(points.begin(), points.end(), result_points.begin(), result_distances.begin(), N, points[i]);
        const auto mean = compute_mean(result_points);
        const auto cov = compute_covariance(result_points, mean);
        results[i] = {mean, cov};
    }
}

void ndt_scan_matching(ndtcpp::mat3x3& trans_mat, const std::vector<ndtcpp::point2>& source_points, std::vector<ndtpoint2>& target_points){
    const size_t max_iter_num = 20;
    const float max_distance2 = 3.0f * 3.0f;

    const size_t target_points_size = target_points.size();
    const size_t source_points_size = source_points.size();

    kdtree::construct(target_points.begin(), target_points.end());
    for(size_t iter = 0; iter < max_iter_num; iter++){
        ndtcpp::mat3x3 H_Mat {
            0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 0.0f
        };

        ndtcpp::point3 b_Point {
            0.0f, 0.0f, 0.0f
        };

        for(auto point_iter = 0; point_iter < source_points_size; point_iter += 10){
            ndtpoint2 query_point = {transformPointCopy(trans_mat, source_points[point_iter]), {}};
            ndtpoint2 target_point;
            float target_distance;
            kdtree::search_knn(target_points.begin(), target_points.end(), &target_point, &target_distance, 1, query_point);

            if(target_distance > max_distance2){continue;}

            const auto identity_plus_cov = ndtcpp::mat3x3{
                target_point.cov.a, target_point.cov.b, 0.0f,
                target_point.cov.c, target_point.cov.d, 0.0f,
                0.0f, 0.0f, 1.0f
            };

            const ndtcpp::mat3x3 target_cov_inv = inverse3x3Copy(identity_plus_cov); //IM


            const auto error = ndtcpp::point3{
                target_point.mean.x - query_point.mean.x,
                target_point.mean.y - query_point.mean.y,
                0.0f
            };

            const ndtcpp::point2 v_point = transformPointCopy(trans_mat, skewd(source_points[point_iter]));

            const auto mat_J = ndtcpp::mat3x3{
                trans_mat.a * -1.0f, trans_mat.b * -1.0f, v_point.x,
                trans_mat.d * -1.0f, trans_mat.e * -1.0f, v_point.y,
                trans_mat.g * -1.0f, trans_mat.h * -1.0f, trans_mat.i * -1.0f
            };

            const ndtcpp::mat3x3 mat_J_T = transpose(mat_J);

            H_Mat += (mat_J_T * (target_cov_inv * mat_J));

            b_Point += (mat_J_T * (target_cov_inv * error));

        }
        b_Point.x *= -1.0;
        b_Point.y *= -1.0;
        b_Point.z *= -1.0;
        const ndtcpp::point3 delta = solve3x3(H_Mat,b_Point);
        trans_mat = trans_mat * expmap(delta);

        if(multtiplyPowPoint3(delta) < 1e-4){
            std::cout << "END NDT. ITER: " << iter << std::endl;
            break;
        }
    }
}

//debug
void writePointsToSVG(const std::vector<ndtcpp::point2>& point_1, const std::vector<ndtcpp::point2>& point_2, const std::string& file_name) {
    std::ofstream file(file_name);
    if (!file.is_open()) {
        std::cerr << "Cannot open file for writing." << std::endl;
        return;
    }

    file << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"500\" height=\"500\">\n";
    file << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    for (const auto& point : point_1) {
        file << "<circle cx=\"" << point.x * 10.0 + 200.0 << "\" cy=\"" << point.y * 10.0 + 200.0 << "\" r=\"1\" fill=\"red\" />\n";
    }

    for (const auto& point : point_2) {
        file << "<circle cx=\"" << point.x * 10.0 + 200.0 << "\" cy=\"" << point.y * 10.0 + 200.0 << "\" r=\"1\" fill=\"black\" />\n";
    }

    file << "</svg>\n";
    file.close();
}

int main(void){
    auto scan_points1 = ndtcpp::read_scan_points("./data/scan_1.txt");
    auto target_points = ndtcpp::read_scan_points("./data/scan_2.txt");

    auto trans_mat1 = makeTransformationMatrix(1.0f, 0.0f, 0.5f);
    transformPointsZeroCopy(trans_mat1, scan_points1);

    auto ndt_points = std::vector<ndtpoint2>();
    auto start_time = std::chrono::high_resolution_clock::now();

    compute_ndt_points(target_points, ndt_points);
    ndt_scan_matching(trans_mat1, scan_points1, ndt_points);

    auto end_time = std::chrono::high_resolution_clock::now();

    transformPointsZeroCopy(trans_mat1, scan_points1);

    //debug
    auto microsec = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() / 1e6;
    std::cout << (microsec) << " mill sec" << std::endl;

    writePointsToSVG(scan_points1, target_points, "scan_points.svg");
}
