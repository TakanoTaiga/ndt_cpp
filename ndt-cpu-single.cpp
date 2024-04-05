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

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream> 
#include <cmath>
#include <algorithm>
#include <numeric>

#include <chrono>

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

void readScanPoints(const std::string& file_path, std::vector<point2>& points){
    std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "File could not be opened." << std::endl;
        return;
    }

    std::string line_str;
    point2 p;
    while(std::getline(file, line_str)){
        std::istringstream iss(line_str);
        if (!(iss >> p.x >> p.y)) {
            std::cerr << "Failed to parse line: " << line_str << std::endl;
            continue;
        }
        points.push_back(p);
    }
}

mat3x3 makeTransformationMatrix(const float& tx, const float& ty, const float& theta) {
    mat3x3 mat = {
        (float)cos(theta), (float)-sin(theta), tx,
        (float)sin(theta), (float)cos(theta), ty,
        0.0f, 0.0f, 1.0f
    };
    return mat;
}

void transformPointsZeroCopy(const mat3x3& mat, std::vector<point2>& points) {
    point2 transformedPoint;

    for (auto& point : points) {
        transformedPoint.x = mat.a * point.x + mat.b * point.y + mat.c;
        transformedPoint.y = mat.d * point.x + mat.e * point.y + mat.f;
        point.x = transformedPoint.x;
        point.y = transformedPoint.y;
    }
}

mat3x3 multiplyMatrices3x3x2(const mat3x3& mat1, const mat3x3& mat2) {
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

point3 multiplyMatrixPoint3(const mat3x3& mat, const point3& vec) {
    point3 result;
    result.x = mat.a * vec.x + mat.b * vec.y + mat.c * vec.z;
    result.y = mat.d * vec.x + mat.e * vec.y + mat.f * vec.z;
    result.z = mat.g * vec.x + mat.h * vec.y + mat.i * vec.z;
    return result;
}

point3 multtiplyPowPoint3(const point3& vec){
    point3 result;
    result.x = vec.x * vec.x;
    result.y = vec.y * vec.y;
    result.z = vec.z * vec.z;
    return result;
}

mat3x3 addMat3x3(const mat3x3& mat1, const mat3x3& mat2){
    const mat3x3 result{
        mat1.a + mat2.a, mat1.b + mat2.b, mat1.c + mat2.c, 
        mat1.d + mat2.d, mat1.e + mat2.e, mat1.f + mat2.f, 
        mat1.g + mat2.g, mat1.h + mat2.h, mat1.i + mat2.i
    };
    return result;
}

point3 addPoint3(const point3& point1, const point3& point2){
    const point3 result{
        point1.x + point2.x, point1.y + point2.y, point1.z + point2.z 
    };
    return result;
}

point3 solve3x3(const mat3x3& m, const point3& p) {
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
    point3 solution;
    solution.z = A[2][3] / A[2][2];
    solution.y = (A[1][3] - A[1][2] * solution.z) / A[1][1];
    solution.x = (A[0][3] - A[0][2] * solution.z - A[0][1] * solution.y) / A[0][0];

    return solution;
}

mat3x3 expmap(const point3& point){
    auto t = point.z;
    auto c = (float)cos(t);
    auto s = (float)sin(t);

    mat2x2 R {
        c, s * -1.0f,
        s, c
    };

    mat3x3 T {
        R.a, R.b, point.x,
        R.c, R.d, point.y,
        0.0f, 0.0f, 1.0f
    };

    return T;
}

point2 transformPointCopy(const mat3x3& mat, const point2& point) {
    point2 transformedPoint;

    transformedPoint.x = mat.a * point.x + mat.b * point.y + mat.c;
    transformedPoint.y = mat.d * point.x + mat.e * point.y + mat.f;

    return transformedPoint;
}

mat3x3 inverse3x3Copy(const mat3x3& mat){
    const auto a = 1.0 / (
        mat.a * mat.e * mat.i + 
        mat.b * mat.f * mat.g +
        mat.c * mat.d * mat.h -
        mat.c * mat.e * mat.g -
        mat.b * mat.d * mat.i -
        mat.a * mat.f * mat.h
        );

    mat3x3 inv_mat;
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

point2 skewd(const point2& input_point){
    const point2 skewd_point {
        input_point.y,
        input_point.x * -1.0f
    };
    return skewd_point;
}

mat3x3 transpose(const mat3x3& input_mat){
    const mat3x3 transpose_mat{
        input_mat.a, input_mat.d, input_mat.g,
        input_mat.b, input_mat.e, input_mat.h,
        input_mat.c, input_mat.f, input_mat.i
    };
    return transpose_mat;
}

point2 compute_mean(const std::vector<point2>& points){
    point2 mean;
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

mat2x2 compute_covariance(const std::vector<point2>& points, const point2& mean){
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

    mat2x2 cov;
    cov.a = vxx / point_size;
    cov.b = vxy / point_size;
    cov.c = vxy / point_size;
    cov.d = vyy / point_size;
    return cov;
}

std::vector<mat2x2> compute_ndt_points(std::vector<point2>& points){
    auto N = 10;
    
    const auto point_size = points.size();
    std::vector<float> distances(point_size);
    std::vector<float> distances_(point_size);

    std::vector<point2> compute_points(N);

    std::vector<mat2x2> covs;

    for(auto &point : points){
        for(auto i = 0; i < point_size; i++){
            const auto distance = (points[i].x - point.x) * (points[i].x - point.x) + (points[i].y - point.y) * (points[i].y - point.y);
            distances[i] = distance;
            distances_[i] = distance;
        }
        
        std::sort(distances.begin(), distances.end());

        for(auto i = 0; i < N; i++){
            const auto target = distances[i];
            for(auto j = 0; j < point_size; j++){
                if(target == distances_[j]){
                    compute_points[i] = points[j];
                    break;
                }
            }
        }

        const auto mean = compute_mean(compute_points);
        const auto cov = compute_covariance(compute_points, mean);

        point = mean;
        covs.push_back(cov);
    }

    return covs;
}

void ndt_scan_matching(mat3x3& trans_mat, const std::vector<point2>& source_points, const std::vector<point2>& target_points, std::vector<mat2x2> target_covs){
    const size_t max_iter_num = 5;
    const float max_distance2 = 3.0f * 3.0f;
    
    const size_t target_points_size = target_points.size();
    const size_t source_points_size = source_points.size();
    std::vector<float> distances(target_points_size);
    std::vector<float> distances_(target_points_size);

    for(size_t iter = 0; iter < max_iter_num; iter++){
        mat3x3 H_Mat {
            0.0f, 0.0f, 0.0f, 
            0.0f, 0.0f, 0.0f, 
            0.0f, 0.0f, 0.0f
        };

        point3 b_Point {
            0.0f, 0.0f, 0.0f
        };

        for(auto point_iter = 0; point_iter < source_points_size; point_iter += 10){
            point2 query_point = transformPointCopy(trans_mat, source_points[point_iter]);
            for(auto i = 0; i < target_points_size; i++){
                const auto distance = (target_points[i].x - query_point.x) * (target_points[i].x - query_point.x) + (target_points[i].y - query_point.y) * (target_points[i].y - query_point.y);
                distances[i] = distance;
                distances_[i] = distance;
            }
            std::sort(distances.begin(), distances.end());

            size_t target_index = 0;
            for(size_t i = 0; i < target_points_size; i++){
                if(distances[0] == distances_[i]){
                    target_index = i;
                }
            }

            const point2 target_point = target_points[target_index];            
            const float target_distance = distances_[target_index];
            const mat2x2 target_cov = target_covs[target_index];

            if(target_distance > max_distance2){continue;}

            const auto identity_plus_cov = mat3x3{
                target_cov.a, target_cov.b, 0.0f,
                target_cov.c, target_cov.d, 0.0f,
                0.0f, 0.0f, 1.0f
            };

            const mat3x3 target_cov_inv = inverse3x3Copy(identity_plus_cov); //IM


            const auto error = point3{
                target_point.x - query_point.x,
                target_point.y - query_point.y,
                0.0f
            };

            const point2 v_point = transformPointCopy(trans_mat, skewd(source_points[point_iter]));

            const auto mat_J = mat3x3{
                trans_mat.a * -1.0f, trans_mat.b * -1.0f, v_point.x,
                trans_mat.d * -1.0f, trans_mat.e * -1.0f, v_point.y,
                trans_mat.g * -1.0f, trans_mat.h * -1.0f, trans_mat.i * -1.0f
            };

            const mat3x3 mat_J_T = transpose(mat_J);

            const mat3x3 imj = multiplyMatrices3x3x2(target_cov_inv, mat_J);
            const mat3x3 hMat_ = multiplyMatrices3x3x2(mat_J_T, imj);
            H_Mat = addMat3x3(H_Mat, hMat_);

            const point3 imerror = multiplyMatrixPoint3(target_cov_inv, error);
            const point3 bPoint_ = multiplyMatrixPoint3(mat_J_T, imerror);
            b_Point = addPoint3(b_Point, bPoint_);

        } 
        b_Point.x *= -1.0;
        b_Point.y *= -1.0;
        b_Point.z *= -1.0;
        const point3 delta = solve3x3(H_Mat,b_Point);
        trans_mat = multiplyMatrices3x3x2(trans_mat, expmap(delta));
    }
}

//debug
void writePointsToSVG(const std::vector<point2>& point_1, const std::vector<point2>& point_2, const std::string& file_name) {
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
    std::vector<point2> scan_points1;
    std::vector<point2> target_points;
    readScanPoints("./data/scan_1.txt", scan_points1);
    readScanPoints("./data/scan_2.txt", target_points);

    auto trans_mat1 = makeTransformationMatrix(1.0f, 0.0f, 0.5f);
    transformPointsZeroCopy(trans_mat1, scan_points1);

    auto start_time = std::chrono::high_resolution_clock::now();

    const auto covs = compute_ndt_points(target_points);
    ndt_scan_matching(trans_mat1, scan_points1, target_points, covs);

    auto end_time = std::chrono::high_resolution_clock::now(); 

    transformPointsZeroCopy(trans_mat1, scan_points1);

    //debug
    auto microsec = (double)(std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()) / 1000;
    std::cout << (int)(microsec / 1000) << " mill sec" << std::endl;

    writePointsToSVG(scan_points1, target_points, "scan_points.svg");
}