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

#ifndef NDTCPP_NDTCPPUTIL_HPP_
#define NDTCPP_NDTCPPUTIL_HPP_

#include <cstddef>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>

#include "type.hpp"


namespace ndtcpp
{
    auto read_scan_points(const std::string& file_path) -> std::vector<point2>
    {
        std::vector<point2> points;
        std::ifstream file(file_path);
        if (!file.is_open()) {
            std::cerr << "File could not be opened." << std::endl;
            return points;
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

        return points;
    }
}

#endif // NDTCPP_NDTCPPUTIL_HPP_
