#include <numeric>
#include <chrono>
#include "ndt-cpu-single.hpp"


int main(void){

    auto scan_points1 = ndtcpp::read_scan_points("./data/scan_1.txt");
    auto target_points = ndtcpp::read_scan_points("./data/scan_2.txt");

    std::vector<double> durations;
    const size_t N = 10;

    for (size_t i = 0; i < N; ++i) {

        auto source = scan_points1;
        auto target = target_points;
        auto trans_mat1 = ndtcpp::makeTransformationMatrix(1.0f, 0.0f, 0.5f);
        ndtcpp::transformPointsZeroCopy(trans_mat1, source);

        auto ndt_points = std::vector<ndtcpp::ndtpoint2>();
        auto start_time = std::chrono::high_resolution_clock::now();

        const bool verbose = true;
        ndtcpp::compute_ndt_points(target, ndt_points);
        ndtcpp::ndt_scan_matching(trans_mat1, source, ndt_points, verbose);

        auto end_time = std::chrono::high_resolution_clock::now();

        ndtcpp::transformPointsZeroCopy(trans_mat1, source);

        //debug
        auto microsec = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() / 1e6;
        durations.push_back(microsec);
        if (i == N - 1) {
            ndtcpp::writePointsToSVG(source, target, "scan_points.svg");
            ndtcpp::writePointsToSVG(source, ndt_points, "scan_points_ndt.svg");
        }
    }
    const double mean = std::accumulate(durations.begin(), durations.end(), 0.0) / durations.size();
    std::cout << "MEAN: " << mean << " mill sec" << std::endl;

}
