
#include "max_probability_cover.h"
#include "min_telescope_cover.h"
#include "BuildGraph.h"

#include <chrono>
#include <fstream>
#include <filesystem> 
#include <numeric> 
#include <iostream>


int candidate_root_count(double deadline) {
    if (deadline < 19.0) {
        return 1500;
    }
    if (deadline < 49.0) {
        return 400;
    }
    if (deadline < 99.0) {
        return 200;
    }
    if (deadline < 199.0) {
        return 100;
    }
    if (deadline < 399.0) {
        return 50;
    }
    if (deadline < 799.0) {
        return 20;
    }
    return 10;
}

int main_mink(int argc, char** argv) {
    const auto total_start = std::chrono::high_resolution_clock::now();

    // Defaults are used only when the executable is launched without arguments.
    std::string out_file = "../results/min_telescope/out.csv";
    std::string tilefile = "../data/tilings/4.0x2.0_tiling.csv";
    std::string mapfile =
        "../data/large_maps_4.0x2.0_tiling/GW191103_012549_638.txt";

    double w_max = 10.0;
    double w_acc = 10.0;
    double dwell_zenith = 1.0;
    double settle_time = 0.0;
    bool is_deepslow = false;
    double budget = 90.0;

    const double mip_gap = 0.001;
    double time_limit = 7200.0;

    const auto print_usage = [argv]() {
        std::cerr
            << "Usage:\n  " << argv[0]
            << " MAP_FILE TILE_FILE DEADLINE W_MAX W_ACC DWELL_ZENITH"
               " IS_DEEPSLOW OUT_FILE TIME_LIMIT_SECONDS\n\n"
            << "Run with no arguments to use the defaults compiled into the driver.\n";
    };

    if (argc == 2 && std::string(argv[1]) == "--help") {
        print_usage();
        return 0;
    }

    // Nine positional arguments plus argv[0]. Requiring the complete list avoids
    // silently combining user-supplied values with unrelated defaults.
    if (argc != 1 && argc != 10) {
        print_usage();
        return 2;
    }

    try {
        if (argc == 10) {
            mapfile = argv[1];
            tilefile = argv[2];
            budget = std::stod(argv[3]);
            w_max = std::stod(argv[4]);
            w_acc = std::stod(argv[5]);
            dwell_zenith = std::stod(argv[6]);
            is_deepslow = (std::stoi(argv[7]) != 0);
            out_file = argv[8];
            time_limit = std::stod(argv[9]);
        }

        if (budget <= 0.0) {
            throw std::invalid_argument("DEADLINE must be positive");
        }
        if (w_max <= 0.0 || w_acc <= 0.0) {
            throw std::invalid_argument("W_MAX and W_ACC must be positive");
        }
        if (dwell_zenith < 0.0) {
            throw std::invalid_argument("DWELL_ZENITH cannot be negative");
        }
        if (time_limit <= 0.0) {
            throw std::invalid_argument("TIME_LIMIT_SECONDS must be positive");
        }
        if (!std::filesystem::is_regular_file(mapfile)) {
            throw std::runtime_error("Map file not found: " + mapfile);
        }
        if (!std::filesystem::is_regular_file(tilefile)) {
            throw std::runtime_error("Tiling file not found: " + tilefile);
        }

        const std::filesystem::path output_path(out_file);
        if (output_path.has_parent_path()) {
            std::filesystem::create_directories(output_path.parent_path());
        }

        const int candidate_roots = candidate_root_count(budget);

        std::cout << std::boolalpha
                  << "Input parameters:\n"
                  << "  map_file          = " << mapfile << '\n'
                  << "  tile_file         = " << tilefile << '\n'
                  << "  output_file       = " << out_file << '\n'
                  << "  deadline          = " << budget << " s\n"
                  << "  w_max             = " << w_max << '\n'
                  << "  w_acc             = " << w_acc << '\n'
                  << "  dwell_zenith      = " << dwell_zenith << " s\n"
                  << "  settle_time       = " << settle_time << " s\n"
                  << "  is_deepslow       = " << is_deepslow << '\n'
                  << "  candidate_roots   = " << candidate_roots << '\n'
                  << "  mip_gap           = " << mip_gap << '\n'
                  << "  time_limit        = " << time_limit << " s\n\n";

        auto starts = run_mink_rooted(
            tilefile, mapfile, out_file,
            budget, dwell_zenith, w_max, w_acc, settle_time,
            is_deepslow, candidate_roots, mip_gap, time_limit);

        mink_annealing_rooted_random(
            tilefile, mapfile, out_file,
            budget, dwell_zenith, w_max, w_acc, settle_time,
            is_deepslow, candidate_roots, starts);

        const auto total_end = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> elapsed_seconds =
            total_end - total_start;
        std::cout << "Total time (wall clock): " << elapsed_seconds.count()
                  << " seconds\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ts_mink: " << error.what() << '\n';
        return 2;
    }
}


int main_maxp(int argc, char** argv) {
    const auto total_start = std::chrono::high_resolution_clock::now();

    // Defaults are used only when the executable is launched without arguments.
    std::string out_file = "../results/max_probability/out.csv";
    std::string tilefile = "../data/tilings/6.9x6.9_tiling.csv";
    std::string mapfile =
        "../data/small_maps_6.9x6.9_tiling/GW200316_215756_50.txt";

    double w_max = 10.0;
    double w_acc = 10.0;
    double dwell_zenith = 1.0;
    double settle_time = 0.0;
    bool is_deepslow = false;
    double budget = 60.0;
    int n_paths = 4;

    const double mip_gap = 0.001;
    double time_limit = 36000.0;

    const auto print_usage = [argv]() {
        std::cerr
            << "Usage:\n  " << argv[0]
            << " MAP_FILE TILE_FILE DEADLINE W_MAX W_ACC DWELL_ZENITH"
               " IS_DEEPSLOW N_PATHS OUT_FILE TIME_LIMIT_SECONDS\n\n"
            << "Run with no arguments to use the defaults compiled into the driver.\n";
    };

    if (argc == 2 && std::string(argv[1]) == "--help") {
        print_usage();
        return 0;
    }

    // Ten positional arguments plus argv[0]. Requiring the complete list avoids
    // silently combining user-supplied values with unrelated defaults.
    if (argc != 1 && argc != 11) {
        print_usage();
        return 2;
    }

    try {
        if (argc == 11) {
            mapfile = argv[1];
            tilefile = argv[2];
            budget = std::stod(argv[3]);
            w_max = std::stod(argv[4]);
            w_acc = std::stod(argv[5]);
            dwell_zenith = std::stod(argv[6]);
            is_deepslow = (std::stoi(argv[7]) != 0);
            n_paths = std::stoi(argv[8]);
            out_file = argv[9];
            time_limit = std::stod(argv[10]);
        }

        if (budget <= 0.0) {
            throw std::invalid_argument("DEADLINE must be positive");
        }
        if (w_max <= 0.0 || w_acc <= 0.0) {
            throw std::invalid_argument("W_MAX and W_ACC must be positive");
        }
        if (dwell_zenith < 0.0) {
            throw std::invalid_argument("DWELL_ZENITH cannot be negative");
        }
        if (n_paths <= 0) {
            throw std::invalid_argument("N_PATHS must be positive");
        }
        if (time_limit <= 0.0) {
            throw std::invalid_argument("TIME_LIMIT_SECONDS must be positive");
        }
        if (!std::filesystem::is_regular_file(mapfile)) {
            throw std::runtime_error("Map file not found: " + mapfile);
        }
        if (!std::filesystem::is_regular_file(tilefile)) {
            throw std::runtime_error("Tiling file not found: " + tilefile);
        }

        const std::filesystem::path output_path(out_file);
        if (output_path.has_parent_path()) {
            std::filesystem::create_directories(output_path.parent_path());
        }

        std::cout << std::boolalpha
                  << "Input parameters:\n"
                  << "  map_file          = " << mapfile << '\n'
                  << "  tile_file         = " << tilefile << '\n'
                  << "  output_file       = " << out_file << '\n'
                  << "  deadline          = " << budget << " s\n"
                  << "  w_max             = " << w_max << '\n'
                  << "  w_acc             = " << w_acc << '\n'
                  << "  dwell_zenith      = " << dwell_zenith << " s\n"
                  << "  settle_time       = " << settle_time << " s\n"
                  << "  is_deepslow       = " << is_deepslow << '\n'
                  << "  n_paths           = " << n_paths << '\n'
                  << "  mip_gap           = " << mip_gap << '\n'
                  << "  time_limit        = " << time_limit << " s\n\n";

        run_greedy_multistarts_annealing(
            tilefile, mapfile, out_file,
            budget, dwell_zenith, w_max, w_acc, settle_time,
            is_deepslow, n_paths);

        run_greedy_multistarts_ilp(
            tilefile, mapfile, out_file,
            budget, dwell_zenith, w_max, w_acc, settle_time,
            is_deepslow, n_paths, mip_gap, time_limit);

        run_top(
            tilefile, mapfile, out_file,
            budget, dwell_zenith, w_max, w_acc, settle_time,
            is_deepslow, n_paths, mip_gap, time_limit);

        const auto total_end = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> elapsed_seconds =
            total_end - total_start;
        std::cout << "Total time (wall clock): " << elapsed_seconds.count()
                  << " seconds\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ts_maxp: " << error.what() << '\n';
        return 2;
    }
}
