#include "ReadData.h"
#include "BuildGraph.h"
#include "dwell_time.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <random>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#define DEG_TO_RAD (M_PI / 180.0)
#define RAD_TO_DEG (180.0 / M_PI)

double angularDistance(double ra1, double dec1, double ra2, double dec2) {
    double ra1_rad = ra1 * DEG_TO_RAD;
    double dec1_rad = dec1 * DEG_TO_RAD;
    double ra2_rad = ra2 * DEG_TO_RAD;
    double dec2_rad = dec2 * DEG_TO_RAD;

    double delta_ra = fmod(ra1_rad - ra2_rad + M_PI, 2 * M_PI) - M_PI;

    double cos_angle = sin(dec1_rad) * sin(dec2_rad) +
                       cos(dec1_rad) * cos(dec2_rad) * cos(delta_ra);

    // Clamp to [-1, 1]
    cos_angle = std::min(1.0, std::max(-1.0, cos_angle));

    double angle_rad = acos(cos_angle);

    return angle_rad * RAD_TO_DEG;
}


// Normalize angle to [0, 360)
double norm_deg(double x_deg) {
    double y = std::fmod(x_deg, 360.0);
    if (y < 0.0) y += 360.0;
    return y;
}

// Signed delta between two angles (deg) along the track that does not cross 0/360 directly, in (-360, 360)
double non_crossing_delta_deg(double a_deg, double b_deg) {
    return norm_deg(a_deg) - norm_deg(b_deg);
}

// Inputs in degrees; RA can be any real, Dec in [-90, 90]
double two_axes_distance_deg(double ra1_deg, double dec1_deg,
                             double ra2_deg, double dec2_deg) {
    const double dRA = std::abs(non_crossing_delta_deg(ra1_deg, ra2_deg));
    const double dDec = std::abs(dec1_deg - dec2_deg);
    // std::cout << std::max(dRA, dDec) << " ";
    return std::max(dRA, dDec);

}

double compute_slew_time(double theta, double theta_min, double w_max, double w_acc) {
    constexpr double EPS = 1e-15;

    if (theta < 0.0) {
        theta = -theta;
    }

    // if (theta > 180.0) {
    //     theta = 180 - theta;
    // }

    if (theta == 0.0) return 0.0;
    if (w_acc <= 0.0 || w_max <= 0.0) throw std::invalid_argument("alpha must be > 0");

    // Trapezoid: accelerate to w_max, cruise, decelerate.
    if (theta  >= theta_min) {
        const double t_acc    = 2 * w_max / w_acc;
        const double t_cruise = (theta - theta_min) / w_max;
        return t_acc + t_cruise;
    } else { // Triangle
        return 2.0 * std::sqrt(theta / w_acc);
    }
}


// incorporate accerlation, decerlation, settle time
// slew_costs[i][j] = anguler_distance[i][j]/slew_rate + settle_time
void compute_slew_costs(const std::vector<double>& ra, const std::vector<double>& dec, 
                        std::vector<std::vector<double>>& slew_costs, 
                        double w_max, double w_acc, double settle_time, bool is_SlowDeep) {

    double theta_min = w_max * w_max / w_acc;

    int n = ra.size();
    slew_costs.assign(n, std::vector<double>(n, 0.0)); 

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i+1; j < n; ++j) {
            if (i != j) {
                double theta;
                if(is_SlowDeep)
                    theta = two_axes_distance_deg(ra[i], dec[i], ra[j], dec[j]);
                else
                    theta = angularDistance(ra[i], dec[i], ra[j], dec[j]);
                double slew_time = compute_slew_time(theta, theta_min, w_max, w_acc);

                slew_costs[i][j] = slew_time + settle_time;
                slew_costs[j][i] = slew_costs[i][j];

            }
        }

        slew_costs[i][i] = 0.0;
    }

}


//Orienteering Graph incorporate accerlation, decerlation, settle time
// cost[i][j] = anguler_distance[i][j]/slew_rate + 0.5*dwelltime[i] + 0.5*dwelltime[j] + settle_time
void compute_costs_orienteering(const std::vector<double>& ra, const std::vector<double>& dec, 
                        const std::vector<double>& dwell_times,
                        std::vector<std::vector<double>>& costs, 
                        double w_max, double w_acc, double settle_time, 
                        int end_pos, double& _padding, bool is_SlowDeep) {

    size_t n = ra.size();
    if(n != end_pos) {
        std::cerr << "Wrong sizing in compute_costs_orienteering func" << std::endl;
    }

    double theta_min = w_max * w_max / w_acc;

    costs.assign(n + 1, std::vector<double>(n + 1, 0.0)); //to include end_pos t'

    _padding = 0.0;
    // double slew_min = 10000;
    // double slew_max = 0.0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i+1; j < n; ++j) {
            if (i != j) {
                double theta;
                if(is_SlowDeep)
                    theta = two_axes_distance_deg(ra[i], dec[i], ra[j], dec[j]);
                else
                    theta = angularDistance(ra[i], dec[i], ra[j], dec[j]);
                double slew_time = compute_slew_time(theta, theta_min, w_max, w_acc);
                // std::cout << slew_time << " ";
                // slew_min = std::min(slew_min, slew_time);
                // slew_max = std::max(slew_max, slew_time);
                costs[i][j] = slew_time + 0.5 * dwell_times[i] + 0.5 * dwell_times[j] + settle_time;
                costs[j][i] = costs[i][j];
                // std::cout << "theta"<<theta<<":" << costs[j][i] << "\t";
                if(slew_time + settle_time > _padding) _padding = slew_time + settle_time;
            }
        }

        costs[i][i] = dwell_times[i];
    }

    _padding = ceil(_padding); //smallest big enough padding
    for (size_t i = 0; i < n; ++i) {
            costs[i][end_pos] = 0.5 * dwell_times[i] + _padding;
            costs[end_pos][i] = costs[i][end_pos];
    }

}

// dummy start (unrooted)
// Orienteering Graph incorporate accerlation, decerlation, settle time
// cost[i][j] = anguler_distance[i][j]/slew_rate + 0.5*dwelltime[i] + 0.5*dwelltime[j] + settle_time
void compute_costs_orienteering_dummy_roots(const std::vector<double>& ra, const std::vector<double>& dec, 
                        const std::vector<double>& dwell_times,
                        std::vector<std::vector<double>>& costs, 
                        double w_max, double w_acc, double settle_time, int start_pos,
                        int end_pos, double& _padding, bool is_SlowDeep) {

    size_t n = ra.size();
    if(n != start_pos) {
        std::cerr << "Wrong sizing in compute_costs_orienteering func" << std::endl;
    }

    double theta_min = w_max * w_max / w_acc;

    costs.assign(n + 2, std::vector<double>(n + 2, 0.0)); //to include end_pos s' and t'

    _padding = 0.0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i+1; j < n; ++j) {
            if (i != j) {
                double theta;
                if(is_SlowDeep)
                    theta = two_axes_distance_deg(ra[i], dec[i], ra[j], dec[j]);
                else
                    theta = angularDistance(ra[i], dec[i], ra[j], dec[j]);
                double slew_time = compute_slew_time(theta, theta_min, w_max, w_acc);
                costs[i][j] = slew_time + 0.5 * dwell_times[i] + 0.5 * dwell_times[j] + settle_time;
                costs[j][i] = costs[i][j];
                if(slew_time + settle_time > _padding) _padding = slew_time + settle_time;
            }
        }

        costs[i][i] = dwell_times[i];
    }

    _padding = ceil(_padding); //smallest big enough padding
    for (size_t i = 0; i < n; ++i) {
            costs[i][start_pos] = 0.5 * dwell_times[i] + _padding;
            costs[start_pos][i] = costs[i][start_pos];
    }
    for (size_t i = 0; i < n; ++i) {
            costs[i][end_pos] = 0.5 * dwell_times[i] + _padding;
            costs[end_pos][i] = costs[i][end_pos];
    }
    costs[start_pos][end_pos] = _padding * 2;
    costs[end_pos][start_pos] = _padding * 2;
}

// overload, incorporate accerlation decerlation and settle time
// build graph for s-t orienteering from pretiled file
std::tuple<int, int, double> buildGraphOrienteering(const std::string& filename, 
                std::vector<std::vector<double>>& costs, 
                std::vector<double>& probability, std::vector<int>& ranks, 
                std::vector<double>& dwell_times, double dwell_zenith,
                double w_max, double w_acc, double settle_time, 
                std::optional<std::pair<double, double>> zenith_pos, 
                bool is_SlowDeep, int init_pos){

    std::vector<double> ra, dec;
    // get probability, ranks, ra, dec
    read_data_from_file(filename, probability, ranks, ra, dec);
    double init_ra, init_dec;
    int n = (int)ra.size();
    if (init_pos >= 0 && init_pos < n) {
        init_ra = ra[init_pos], init_dec = dec[init_pos];
    } 
    else if (init_pos < 0) {
        unsigned int fixed_seed = 1130;
        std::mt19937 generator(fixed_seed);
        std::uniform_int_distribution<int> distribution(0, n-1);
        init_pos = distribution(generator);
        init_ra = ra[init_pos], init_dec = dec[init_pos];
    }
    else {
        throw std::invalid_argument("Invalid init position.");
    }
    std::cout << "init_pos: " << init_pos << std::endl;
    int end_rank = probability.size();

    dwell_times.resize(probability.size());

    double ra_zenith = ra[0], dec_zenith = dec[0];
    if (zenith_pos.has_value()){
        ra_zenith = zenith_pos->first, dec_zenith = zenith_pos->second;
    }


    std::vector<size_t> valid_indices;
    for (size_t i = 0; i < probability.size(); ++i) {
        double zenithAngle = angularDistance(ra[i], dec[i], ra_zenith, dec_zenith);
        zenithAngle = std::min(zenithAngle, 85.0);
        dwell_times[i] = dwell_time_airmass(zenithAngle, dwell_zenith);
        valid_indices.push_back(i);
    }

    // Filter all vectors based on valid_indices
    std::vector<double> filtered_ra, filtered_dec, filtered_probability, filtered_dwell_times;
    std::vector<int> filtered_ranks;
    for (size_t idx : valid_indices) {
        filtered_ra.push_back(ra[idx]);
        filtered_dec.push_back(dec[idx]);
        filtered_probability.push_back(probability[idx]);
        filtered_ranks.push_back(ranks[idx]);
        filtered_dwell_times.push_back(dwell_times[idx]);
    }
    ra = std::move(filtered_ra);
    dec = std::move(filtered_dec);
    probability = std::move(filtered_probability);
    ranks = std::move(filtered_ranks);
    dwell_times = std::move(filtered_dwell_times);

    int start_pos = probability.size();
    ra.push_back(init_ra);
    dec.push_back(init_dec);
    probability.push_back(0.0);
    ranks.push_back(init_pos);
    dwell_times.push_back(0.0);

    // Push as end node
    int end_pos = probability.size();
    probability.push_back(0.0);
    ranks.push_back(-1);
    dwell_times.push_back(0.0);

    double padding_value;
    compute_costs_orienteering(ra, dec, dwell_times, costs, w_max, w_acc, settle_time, end_pos, padding_value, is_SlowDeep);

    return {start_pos, end_pos, padding_value};
}



// overload, incorporate accerlation decerlation and settle time
// build graph for s-t orienteering from likelihood map and tile file
std::tuple<int, int, double> buildGraphOrienteering(const std::string& tilefile, 
                const std::string& pro_map_file,
                std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, std::vector<int>& ranks, 
               std::vector<double>& dwell_times, double dwell_zenith,
               double w_max, double w_acc, double settle_time, 
               bool is_SlowDeep, int& num_tiles, int init_pos) {

    std::vector<double> ra, dec;
    std::vector<std::vector<int>> member_pixels;
    readTiles(tilefile, ranks, ra, dec, member_pixels);

    // std::vector<int> pixel_ids;
    std::vector<double> pixel_probs;
    // int pixel_count = 12 * 64 * 64;
    // readHEALPixelsTXT(pro_map_file, pixel_probs);
    // readHEALPixelsHDF5(pro_map_file, pixel_probs, 64);
    readHEALPixels(pro_map_file, pixel_probs, 64);


    // std::vector<double> probability;
    computeTileProbs(member_pixels, pixel_probs, probability);

    // double init_ra, init_dec;
    // int n = (int)ra.size();
    // if (init_pos >= 0 && init_pos < n) {
    //     init_ra = ra[init_pos], init_dec = dec[init_pos];
    // } 
    // else if (init_pos < 0) {
    //     unsigned int fixed_seed = 1130;
    //     std::mt19937 generator(fixed_seed);
    //     std::uniform_int_distribution<int> distribution(0, n-1);
    //     init_pos = distribution(generator);
    //     init_ra = ra[init_pos], init_dec = dec[init_pos];
    // }
    // else {
    //     throw std::invalid_argument("Invalid init position.");
    // }
    // std::cout << "init_pos: " << init_pos << std::endl;
    
    double coverage = 0.99;
    int num_tile = getTopTileIndices(probability, ranks, ra, dec, member_pixels, coverage);
    // int num_tile = getNonZeroTileIndices(probability, ranks, ra, dec, member_pixels);
    std::cout << "num tiles after trim: " << num_tile << "\n";

    double init_ra, init_dec;
    int n = (int)ra.size();
    if (init_pos >= 0 && init_pos < n) {
        init_ra = ra[init_pos], init_dec = dec[init_pos];
    } 
    else if (init_pos < 0) {
        unsigned int fixed_seed = 1130;
        std::mt19937 generator(fixed_seed);
        std::uniform_int_distribution<int> distribution(0, n-1);
        init_pos = distribution(generator);
        init_ra = ra[init_pos], init_dec = dec[init_pos];
    }
    else {
        throw std::invalid_argument("Invalid init position.");
    }
    std::cout << "init_pos: " << ranks[init_pos] << std::endl;

    dwell_times.resize(probability.size());
    double dwelltime_ref = dwell_zenith;
    double ra_zenith = ra[0], dec_zenith = dec[0];

    std::vector<size_t> valid_indices;
    for (size_t i = 0; i < probability.size(); ++i) {
        double zenithAngle = angularDistance(ra[i], dec[i], ra_zenith, dec_zenith);
        zenithAngle = std::min(zenithAngle, 85.0);
        dwell_times[i] = dwell_time_airmass(zenithAngle, dwell_zenith);
        valid_indices.push_back(i);
    }
    std::cout << "\n";

    // // Filter all vectors based on valid_indices
    // std::vector<double> filtered_ra, filtered_dec, filtered_probability, filtered_dwell_times;
    // std::vector<int> filtered_ranks;
    // for (size_t idx : valid_indices) {
    //     filtered_ra.push_back(ra[idx]);
    //     filtered_dec.push_back(dec[idx]);
    //     filtered_probability.push_back(probability[idx]);
    //     filtered_ranks.push_back(ranks[idx]);
    //     filtered_dwell_times.push_back(dwell_times[idx]);
    // }
    // ra = std::move(filtered_ra);
    // dec = std::move(filtered_dec);
    // probability = std::move(filtered_probability);
    // ranks = std::move(filtered_ranks);
    // dwell_times = std::move(filtered_dwell_times);
    // std::cout << "num tiles after filter: " << ra.size() << "\n";
    num_tiles = ra.size();

    int start_pos = probability.size();
    ra.push_back(init_ra);
    dec.push_back(init_dec);
    probability.push_back(0.0);
    ranks.push_back(ranks[init_pos]);
    dwell_times.push_back(0.0);

    // Push as end node
    int end_pos = probability.size();
    probability.push_back(0.0);
    ranks.push_back(-1);
    dwell_times.push_back(0.0);

    double padding_value;
    compute_costs_orienteering(ra, dec, dwell_times, costs, w_max, w_acc, settle_time, end_pos, padding_value, is_SlowDeep);

    return {start_pos, end_pos, padding_value};
}

//dummy start (unrooted)
// build graph for s-t orienteering from likelihood map and tile file wirh dummy start
std::tuple<int, int, double> buildGraphOrienteeringDummyRoot(const std::string& tilefile, 
            const std::string& pro_map_file,
            std::vector<std::vector<double>>& costs, 
            std::vector<double>& probability, std::vector<int>& ranks, 
            std::vector<double>& dwell_times, double dwell_zenith,
            double w_max, double w_acc, double settle_time, 
            bool is_SlowDeep, int& num_tiles) {

    std::vector<double> ra, dec;
    std::vector<std::vector<int>> member_pixels;
    readTiles(tilefile, ranks, ra, dec, member_pixels);

    // std::vector<int> pixel_ids;
    std::vector<double> pixel_probs;
    readHEALPixels(pro_map_file, pixel_probs, 64);

    // std::vector<double> probability;
    computeTileProbs(member_pixels, pixel_probs, probability);
    
    double coverage = 0.99;
    num_tiles = getTopTileIndices(probability, ranks, ra, dec, member_pixels, coverage);
    // int num_tile = getNonZeroTileIndices(probability, ranks, ra, dec, member_pixels);
    std::cout << "num tiles after trim: " << num_tiles << "\n";


    dwell_times.resize(probability.size());
    double dwelltime_ref = dwell_zenith;
    double ra_zenith = ra[0], dec_zenith = dec[0];

    std::vector<size_t> valid_indices;
    for (size_t i = 0; i < probability.size(); ++i) {
        double zenithAngle = angularDistance(ra[i], dec[i], ra_zenith, dec_zenith);
        zenithAngle = std::min(zenithAngle, 85.0);
        dwell_times[i] = dwell_time_airmass(zenithAngle, dwell_zenith);
        valid_indices.push_back(i);
    }
    std::cout << "\n";


    // Push as end node
    int start_pos = probability.size();
    probability.push_back(0.0);
    ranks.push_back(-2);
    dwell_times.push_back(0.0);
    int end_pos = probability.size();
    probability.push_back(0.0);
    ranks.push_back(-1);
    dwell_times.push_back(0.0);

    double padding_value;
    // compute_costs_orienteering(ra, dec, dwell_times, costs, w_max, w_acc, settle_time, end_pos, padding_value, is_SlowDeep);
    compute_costs_orienteering_dummy_roots(ra, dec, dwell_times, costs, w_max, w_acc, settle_time, start_pos, end_pos, padding_value, is_SlowDeep);

    return {start_pos, end_pos, padding_value*2};
}

// build team orienteering graph from pre-tiledfile
std::tuple<std::vector<int>, int, double> buildGraphTeamOrienteering(const std::string& filename, 
                std::vector<std::vector<double>>& costs, 
                std::vector<double>& probability, std::vector<int>& ranks, 
                std::vector<double>& dwell_times, double dwell_zenith,
                double w_max, double w_acc, double settle_time, 
                bool is_SlowDeep, int m, 
                std::optional<std::pair<double, double>> zenith_pos, 
                std::vector<std::pair<double, double>> init_poses) {
    std::vector<double> ra, dec;

    // get probability, ranks, ra, dec
    read_data_from_file(filename, probability, ranks, ra, dec);
    double ra_zenith = ra[0], dec_zenith = dec[0];
    if (zenith_pos.has_value()){
        ra_zenith = zenith_pos->first, dec_zenith = zenith_pos->second;
    }

    dwell_times.resize(probability.size());
    double dwelltime_ref = 1.0;

    std::vector<size_t> valid_indices;
    for (size_t i = 0; i < probability.size(); ++i) {
        double zenithAngle = angularDistance(ra[i], dec[i], ra_zenith, dec_zenith);
        zenithAngle = std::min(zenithAngle, 85.0);
        dwell_times[i] = dwell_time_airmass(zenithAngle, dwell_zenith);
        valid_indices.push_back(i);
    }

    // Filter all vectors based on valid_indices
    std::vector<double> filtered_ra, filtered_dec, filtered_probability, filtered_dwell_times;
    std::vector<int> filtered_ranks;
    for (size_t idx : valid_indices) {
        filtered_ra.push_back(ra[idx]);
        filtered_dec.push_back(dec[idx]);
        filtered_probability.push_back(probability[idx]);
        filtered_ranks.push_back(ranks[idx]);
        filtered_dwell_times.push_back(dwell_times[idx]);
    }
    ra = std::move(filtered_ra);
    dec = std::move(filtered_dec);
    probability = std::move(filtered_probability);
    ranks = std::move(filtered_ranks);
    dwell_times = std::move(filtered_dwell_times);

    int n = (int)ra.size();
    if (init_poses.size() < m) {
        unsigned int fixed_seed = 1130;
        std::mt19937 generator(fixed_seed);
        std::uniform_int_distribution<int> distribution(0, n-1);
        for(int i = init_poses.size(); i < m; ++i) {
            int random_index = distribution(generator);
            init_poses.push_back({ra[random_index], dec[random_index]});
            int r_ = ranks[random_index];
            ranks.push_back(r_);
        }
    }

    std::vector<int> start_poses;
    for(auto pos : init_poses) {
        int idx = (int)ra.size();
        start_poses.push_back(idx);
        ra.push_back(pos.first);
        dec.push_back(pos.second);
        probability.push_back(0.0);
        // ranks.push_back(idx);
        dwell_times.push_back(0.0);
    }

    // Push as end node
    int end_rank = probability.size();
    int end_pos = probability.size();
    probability.push_back(0.0);
    ranks.push_back(-1);
    dwell_times.push_back(0.0);

    double padding_value;
    // compute_costs_orienteering(ra, dec, dwell_times, costs, slew_rate, end_pos, padding_value);
    compute_costs_orienteering(ra, dec, dwell_times, costs, w_max, w_acc, 
                                settle_time, end_pos, padding_value, is_SlowDeep);

    return {start_poses, end_pos, padding_value};
}

// build team orienteering graph from likelihood map+tilefile
std::tuple<std::vector<int>, int, double> buildGraphTeamOrienteering(const std::string& tilefile,
                                const std::string& pro_map_file,
                                std::vector<std::vector<double>>& costs,
                                std::vector<double>& probability,
                                std::vector<int>& ranks,
                                std::vector<double>& dwell_times,
                                double dwell_zenith,
                                double w_max, double w_acc, double settle_time,
                                bool is_SlowDeep, int m, int& num_tiles,
                                std::optional<std::pair<double, double>> zenith_pos,
                                std::vector<std::pair<double, double>> init_poses) {

    std::vector<double> ra, dec;
    std::vector<std::vector<int>> member_pixels;
    readTiles(tilefile, ranks, ra, dec, member_pixels);
 
    std::vector<double> pixel_probs;
    readHEALPixels(pro_map_file, pixel_probs, 64);
 
    computeTileProbs(member_pixels, pixel_probs, probability);
 
    double coverage = 0.99;
    int num_tile = getTopTileIndices(probability, ranks, ra, dec,
                                    member_pixels, coverage);
    std::cout << "num tiles after trim: " << num_tile << "\n";
 

    //  3. Fill missing initial positions with random tile positions
    const int n = static_cast<int>(ra.size());  // number of tiles after trim
 
    // Track which tile index each init_pose came from (for rank lookup)
    std::vector<int> init_pose_rank_source;
    init_pose_rank_source.reserve(m);
 
    // For user-supplied positions we don't have a tile index, so store -1
    for (std::size_t i = 0; i < init_poses.size(); ++i)
        init_pose_rank_source.push_back(-2);
 
    if (static_cast<int>(init_poses.size()) < m) {
        unsigned int fixed_seed = 1130;
        std::mt19937 generator(fixed_seed);
        std::uniform_int_distribution<int> distribution(0, n - 1);
 
        for (int i = static_cast<int>(init_poses.size()); i < m; ++i) {
            int random_index = distribution(generator);
            init_poses.push_back({ra[random_index], dec[random_index]});
            init_pose_rank_source.push_back(random_index);
        }
    }
 
    // Compute dwell times for tile nodes
    dwell_times.resize(probability.size());
 
    double ra_zenith  = ra[0];
    double dec_zenith = dec[0];
    if (zenith_pos.has_value()) {
        ra_zenith  = zenith_pos->first;
        dec_zenith = zenith_pos->second;
    }
 
    for (std::size_t i = 0; i < probability.size(); ++i) {
        double zenithAngle = angularDistance(ra[i], dec[i],
                                            ra_zenith, dec_zenith);
        zenithAngle = std::min(zenithAngle, 85.0);
        dwell_times[i] = dwell_time_airmass(zenithAngle, dwell_zenith);
    }
 
    num_tiles = static_cast<int>(ra.size());
 

    //  Append start nodes (one per telescope / init position)
    //     Each start node has prize = 0 and dwell = 0.
    std::vector<int> start_indices;
    start_indices.reserve(m);
 
    for (int i = 0; i < static_cast<int>(init_poses.size()); ++i) {
        int idx = static_cast<int>(ra.size());
        start_indices.push_back(idx);
 
        ra.push_back(init_poses[i].first);
        dec.push_back(init_poses[i].second);
        probability.push_back(0.0);
        dwell_times.push_back(0.0);
 
        // Assign a rank: use the source tile's rank if we have one,
        // otherwise use a sentinel (-1).
        if (init_pose_rank_source[i] >= 0)
            ranks.push_back(ranks[init_pose_rank_source[i]]);
        else
            ranks.push_back(-1);
    }
 

    //  Append end node (dummy node with no real coordinates)
    int end_pos = static_cast<int>(probability.size());
    probability.push_back(0.0);
    ranks.push_back(-1);
    dwell_times.push_back(0.0);
 
  
    //  Build the cost matrix
    double padding_value = 0.0;
    compute_costs_orienteering(ra, dec, dwell_times, costs,
                               w_max, w_acc, settle_time,
                               end_pos, padding_value, is_SlowDeep);
 
    return {start_indices, end_pos, padding_value};
}


// build team orienteering graph from likelihood map+tilefile
// record ra dec
std::tuple<std::vector<int>, int, double> buildGraphTeamOrienteering(const std::string& tilefile,
                                const std::string& pro_map_file,
                                std::vector<std::vector<double>>& costs,
                                std::vector<double>& probability,
                                std::vector<int>& ranks,
                                std::vector<double>& dwell_times,
                                std::vector<double>& ra,
                                std::vector<double>& dec,
                                double dwell_zenith,
                                double w_max, double w_acc, double settle_time,
                                bool is_SlowDeep, int m, int& num_tiles,
                                std::optional<std::pair<double, double>> zenith_pos,
                                std::vector<std::pair<double, double>> init_poses) {

    // std::vector<double> ra, dec;
    std::vector<std::vector<int>> member_pixels;
    readTiles(tilefile, ranks, ra, dec, member_pixels);
 
    std::vector<double> pixel_probs;
    readHEALPixels(pro_map_file, pixel_probs, 64);
 
    computeTileProbs(member_pixels, pixel_probs, probability);

    //  Fill missing initial positions with random tile positions
    const int n = static_cast<int>(ra.size());  // number of tiles after trim
 
    // Track which tile index each init_pose came from (for rank lookup)
    std::vector<int> init_pose_rank_source;
    init_pose_rank_source.reserve(m);
 
    // For user-supplied positions we don't have a tile index, so store -1
    for (std::size_t i = 0; i < init_poses.size(); ++i)
        init_pose_rank_source.push_back(-2);
 
    if (static_cast<int>(init_poses.size()) < m) {
        unsigned int fixed_seed = 1130;
        std::mt19937 generator(fixed_seed);
        std::uniform_int_distribution<int> distribution(0, n - 1);
 
        for (int i = static_cast<int>(init_poses.size()); i < m; ++i) {
            int random_index = distribution(generator);
            init_poses.push_back({ra[random_index], dec[random_index]});
            init_pose_rank_source.push_back(ranks[random_index]);
        }
    }
 
    double coverage = 0.99;
    int num_tile = getTopTileIndices(probability, ranks, ra, dec,
                                    member_pixels, coverage);
    std::cout << "num tiles after trim: " << num_tile << "\n";
 
    // Compute dwell times for tile nodes
    dwell_times.resize(probability.size());
 
    double ra_zenith  = ra[0];
    double dec_zenith = dec[0];
    if (zenith_pos.has_value()) {
        ra_zenith  = zenith_pos->first;
        dec_zenith = zenith_pos->second;
    }
 
    for (std::size_t i = 0; i < probability.size(); ++i) {
        double zenithAngle = angularDistance(ra[i], dec[i],
                                            ra_zenith, dec_zenith);
        zenithAngle = std::min(zenithAngle, 85.0);
        dwell_times[i] = dwell_time_airmass(zenithAngle, dwell_zenith);
    }
 
    num_tiles = static_cast<int>(ra.size());
 

    //  Append start nodes (one per telescope / init position)
    //     Each start node has prize = 0 and dwell = 0.
    std::vector<int> start_indices;
    start_indices.reserve(m);
 
    for (int i = 0; i < static_cast<int>(init_poses.size()); ++i) {
        int idx = static_cast<int>(ra.size());
        start_indices.push_back(idx);
 
        ra.push_back(init_poses[i].first);
        dec.push_back(init_poses[i].second);
        probability.push_back(0.0);
        dwell_times.push_back(0.0);
        ranks.push_back(init_pose_rank_source[i]);
 
        // // Assign a rank: use the source tile's rank if we have one,
        // // otherwise use a sentinel (-2).
        // if (init_pose_rank_source[i] >= 0)
        //     ranks.push_back(ranks[init_pose_rank_source[i]]);
        // else
        //     ranks.push_back(-2);
    }
 

    //  Append end node (dummy node with no real coordinates)
    int end_pos = static_cast<int>(probability.size());
    probability.push_back(0.0);
    ranks.push_back(-1);
    dwell_times.push_back(0.0);
 
  
    //  Build the cost matrix
    double padding_value = 0.0;
    compute_costs_orienteering(ra, dec, dwell_times, costs,
                               w_max, w_acc, settle_time,
                               end_pos, padding_value, is_SlowDeep);
 
    return {start_indices, end_pos, padding_value};
}


// // build team orienteering graph from likelihood map+tilefile with init_pose in index
// // init_pos in index
// std::tuple<std::vector<int>, int, double> buildGraphTeamOrienteering(const std::string& tilefile,
//                                 const std::string& pro_map_file,
//                                 std::vector<std::vector<double>>& costs,
//                                 std::vector<double>& probability,
//                                 std::vector<int>& ranks,
//                                 std::vector<double>& dwell_times,
//                                 double dwell_zenith,
//                                 double w_max, double w_acc, double settle_time,
//                                 bool is_SlowDeep, int m, int& num_tiles,
//                                 std::vector<int> init_poses,
//                                 std::optional<std::pair<double, double>> zenith_pos) {

//     std::vector<double> ra, dec;
//     std::vector<std::vector<int>> member_pixels;
//     readTiles(tilefile, ranks, ra, dec, member_pixels);

//     std::vector<double> pixel_probs;
//     readHEALPixels(pro_map_file, pixel_probs, 64);

//     computeTileProbs(member_pixels, pixel_probs, probability);

//     double coverage = 0.99;
//     int num_tile = getTopTileIndices(probability, ranks, ra, dec,
//                                     member_pixels, coverage);
//     std::cout << "num tiles after trim: " << num_tile << "\n";

//     const int n = static_cast<int>(ra.size());  // number of tiles after trim

//     // Track which tile index each init_pose came from (for rank lookup)
//     std::vector<int> init_pose_rank_source;
//     init_pose_rank_source.reserve(m);

//     for (std::size_t i = 0; i < init_poses.size(); ++i)
//         init_pose_rank_source.push_back(init_poses[i]);

//     // Fill remaining starts with random tile indices if fewer than m provided
//     if (static_cast<int>(init_poses.size()) < m) {
//         unsigned int fixed_seed = 1130;
//         std::mt19937 generator(fixed_seed);
//         std::uniform_int_distribution<int> distribution(0, n - 1);

//         for (int i = static_cast<int>(init_poses.size()); i < m; ++i) {
//             int random_index = distribution(generator);
//             init_poses.push_back(random_index); 
//             init_pose_rank_source.push_back(random_index);
//         }
//     }

//     // Compute dwell times for tile nodes
//     dwell_times.resize(probability.size());

//     double ra_zenith  = ra[0];
//     double dec_zenith = dec[0];
//     if (zenith_pos.has_value()) {
//         ra_zenith  = zenith_pos->first;
//         dec_zenith = zenith_pos->second;
//     }

//     for (std::size_t i = 0; i < probability.size(); ++i) {
//         double zenithAngle = angularDistance(ra[i], dec[i],
//                                             ra_zenith, dec_zenith);
//         zenithAngle = std::min(zenithAngle, 85.0);
//         dwell_times[i] = dwell_time_airmass(zenithAngle, dwell_zenith);
//     }

//     num_tiles = n;  // FIX 2: use n (tile count before appending dummies)

//     std::vector<int> start_indices;
//     start_indices.reserve(m);

//     for (int i = 0; i < static_cast<int>(init_poses.size()); ++i) {
//         int idx = static_cast<int>(ra.size());
//         start_indices.push_back(idx);

//         ra.push_back(ra[init_poses[i]]);
//         dec.push_back(dec[init_poses[i]]);
//         probability.push_back(0.0);
//         dwell_times.push_back(0.0);

//         if (init_pose_rank_source[i] >= 0)
//             ranks.push_back(ranks[init_pose_rank_source[i]]);
//         else
//             ranks.push_back(-1);
//     }

//     // Append end node (dummy node with no real coordinates)
//     int end_pos = static_cast<int>(probability.size());
//     probability.push_back(0.0);
//     ranks.push_back(-1);
//     dwell_times.push_back(0.0);

//     // Build the cost matrix
//     double padding_value = 0.0;
//     compute_costs_orienteering(ra, dec, dwell_times, costs,
//                                w_max, w_acc, settle_time,
//                                end_pos, padding_value, is_SlowDeep);

//     return {start_indices, end_pos, padding_value};
// }


// rebuild graph for s-t orienteering to exlude a existing path 
// and set the last node in the existing path as the start node
std::pair<int, int> removeNodesOrienteering(const std::vector<double>& prizes,
                const std::vector<std::vector<double>>& costs,
                const std::vector<double>& dwell_times,
                const std::vector<int>& nodes_to_remove,
                std::vector<std::vector<double>>& costs_removed,
                std::vector<double>& prizes_removed,
                std::vector<int>& ranks_mapping,
                int start_idx, int end_idx) {

    costs_removed.clear();
    prizes_removed.clear();
    ranks_mapping.clear();

    if (nodes_to_remove.empty()) {
        costs_removed = costs;
        prizes_removed = prizes;
        ranks_mapping.resize(prizes.size());
        std::iota(ranks_mapping.begin(), ranks_mapping.end(), 0);
        return {start_idx, end_idx};
    }
        
    std::unordered_set<int> remove_set(nodes_to_remove.begin(), nodes_to_remove.end());
    remove_set.erase(start_idx);
    remove_set.erase(end_idx);

    std::vector<int> kept;
    kept.reserve(prizes.size() - remove_set.size());

    kept.push_back(start_idx);
    kept.push_back(end_idx);

    for (int old_idx = 0; old_idx < int(prizes.size()); ++old_idx)
        if (remove_set.find(old_idx) == remove_set.end() && old_idx != start_idx && old_idx != end_idx)
            kept.push_back(old_idx);


    const std::size_t S = kept.size();

    costs_removed.assign(S, std::vector<double>(S));
    prizes_removed.resize(S);

    ranks_mapping = kept;

    for (std::size_t i = 0; i < S; ++i) {
        const int old_i = kept[i];
        prizes_removed[i] = prizes[old_i];
        for (std::size_t j = i+1; j < S; ++j) {
            const int old_j = kept[j];
            costs_removed[i][j] = costs[old_i][old_j];
            if(old_i == start_idx || old_j == start_idx) 
                costs_removed[i][j] -= 0.5 * dwell_times[start_idx];
            costs_removed[j][i] = costs_removed[i][j];
        }
    }

    int new_start_idx = 0;
    int new_end_idx = 1;
    return {new_start_idx, new_end_idx}; 
}



std::pair<std::vector<int>, int>
removeVisitedNodes(std::vector<double>& prizes,
                   std::vector<std::vector<double>>& costs,
                   std::vector<int>& ranks,
                   const std::vector<int>& nodes_to_remove,
                   std::vector<int>& start_indices,
                   int end_idx) {
    if (nodes_to_remove.empty()) {
        return {start_indices, end_idx};
    }
 
    // Build the removal set, but never remove end_idx or any start index
    std::unordered_set<int> remove_set(nodes_to_remove.begin(),
                                       nodes_to_remove.end());
    remove_set.erase(end_idx);
    for (int si : start_indices)
        remove_set.erase(si);
 
    const int old_N = static_cast<int>(prizes.size());
 
    // Build the list of kept old indices (preserving original order)
    std::vector<int> kept;
    kept.reserve(old_N - static_cast<int>(remove_set.size()));
 
    for (int old_idx = 0; old_idx < old_N; ++old_idx) {
        if (remove_set.find(old_idx) == remove_set.end())
            kept.push_back(old_idx);
    }
 
    // Build old to new index map
    std::unordered_map<int, int> old_to_new;
    for (int new_idx = 0; new_idx < static_cast<int>(kept.size()); ++new_idx)
        old_to_new[kept[new_idx]] = new_idx;
 
    const std::size_t S = kept.size();
 
    // Build reduced arrays
    std::vector<std::vector<double>> costs_reduced(S, std::vector<double>(S, 0.0));
    std::vector<double> prizes_reduced(S);
    std::vector<int> ranks_reduced(S);
 
    for (std::size_t i = 0; i < S; ++i) {
        const int old_i = kept[i];
        prizes_reduced[i] = prizes[old_i];
        ranks_reduced[i]  = ranks[old_i];
        for (std::size_t j = i + 1; j < S; ++j) {
            const int old_j = kept[j];
            costs_reduced[i][j] = costs[old_i][old_j];
            costs_reduced[j][i] = costs_reduced[i][j];
        }
    }
 
    // Remap start indices and end index
    std::vector<int> new_start_indices;
    new_start_indices.reserve(start_indices.size());
    for (int si : start_indices) {
        auto it = old_to_new.find(si);
        if (it != old_to_new.end())
            new_start_indices.push_back(it->second);
    }
    int new_end_idx = old_to_new.at(end_idx);
 
    // Overwrite in place
    costs  = std::move(costs_reduced);
    prizes = std::move(prizes_reduced);
    ranks  = std::move(ranks_reduced);
 
    return {new_start_indices, new_end_idx};
}


// // build team orienteering graph from likelihood map+tilefile
// std::vector<int> buildGraphMultiStarts(const std::string& tilefile,
//                                 const std::string& pro_map_file,
//                                 std::vector<std::vector<double>>& costs,
//                                 std::vector<double>& probability,
//                                 std::vector<int>& ranks,
//                                 std::vector<double>& dwell_times,
//                                 double dwell_zenith,
//                                 double w_max, double w_acc, double settle_time,
//                                 bool is_SlowDeep, int m, int& num_tiles,
//                                 std::optional<std::pair<double, double>> zenith_pos,
//                                 std::vector<std::pair<double, double>> init_poses) {

//     std::vector<double> ra, dec;
//     std::vector<std::vector<int>> member_pixels;
//     readTiles(tilefile, ranks, ra, dec, member_pixels);
 
//     std::vector<double> pixel_probs;
//     readHEALPixels(pro_map_file, pixel_probs, 64);
 
//     computeTileProbs(member_pixels, pixel_probs, probability);
 
//     double coverage = 0.99;
//     int num_tile = getTopTileIndices(probability, ranks, ra, dec,
//                                     member_pixels, coverage);
//     std::cout << "num tiles after trim: " << num_tile << "\n";
 

//     //  Fill missing initial positions with random tile positions
//     const int n = static_cast<int>(ra.size());  // number of tiles after trim
 
//     // Track which tile index each init_pose came from (for rank lookup)
//     std::vector<int> init_pose_rank_source;
//     init_pose_rank_source.reserve(m);
 
//     // For user-supplied positions we don't have a tile index, so store -2, -1 is t.
//     for (std::size_t i = 0; i < init_poses.size(); ++i)
//         init_pose_rank_source.push_back(-2);
 
//     if (static_cast<int>(init_poses.size()) < m && n > 0) {
//         unsigned int fixed_seed = 1130;
//         std::mt19937 generator(fixed_seed);
//         std::uniform_int_distribution<int> distribution(0, n - 1);
 
//         for (int i = static_cast<int>(init_poses.size()); i < m; ++i) {
//             int random_index = distribution(generator);
//             init_poses.push_back({ra[random_index], dec[random_index]});
//             init_pose_rank_source.push_back(random_index);
//         }
//     }
 
//     // Compute dwell times for tile nodes
//     dwell_times.resize(probability.size());
 
//     double ra_zenith  = ra[0];
//     double dec_zenith = dec[0];
//     if (zenith_pos.has_value()) {
//         ra_zenith  = zenith_pos->first;
//         dec_zenith = zenith_pos->second;
//     }
 
//     for (std::size_t i = 0; i < probability.size(); ++i) {
//         double zenithAngle = angularDistance(ra[i], dec[i],
//                                             ra_zenith, dec_zenith);
//         zenithAngle = std::min(zenithAngle, 85.0);
//         dwell_times[i] = dwell_time_airmass(zenithAngle, dwell_zenith);
//     }
 
//     num_tiles = static_cast<int>(ra.size());
 

//     //  Append start nodes (one per telescope / init position)
//     //     Each start node has prize = 0 and dwell = 0.
//     std::vector<int> start_indices;
//     start_indices.reserve(m);
 
//     for (int i = 0; i < static_cast<int>(init_poses.size()); ++i) {
//         int idx = static_cast<int>(ra.size());
//         start_indices.push_back(idx);
 
//         ra.push_back(init_poses[i].first);
//         dec.push_back(init_poses[i].second);
//         probability.push_back(0.0);
//         dwell_times.push_back(0.0);
 
//         // Assign a rank: use the source tile's rank if we have one,
//         // otherwise use a sentinel (-1).
//         if (init_pose_rank_source[i] >= 0)
//             ranks.push_back(ranks[init_pose_rank_source[i]]);
//         else
//             ranks.push_back(-1);
//     }

//     compute_slew_costs(ra, dec, costs, w_max, w_acc, settle_time, is_SlowDeep);
 
//     return start_indices;
// }



// build team orienteering graph from likelihood map+tilefile
std::vector<int> buildGraphMultiStarts(const std::string& tilefile,
                                const std::string& pro_map_file,
                                std::vector<std::vector<double>>& costs,
                                std::vector<double>& probability,
                                std::vector<int>& ranks,
                                std::vector<double>& dwell_times,
                                std::vector<double>& ra,
                                std::vector<double>& dec,
                                double dwell_zenith,
                                double w_max, double w_acc, double settle_time,
                                bool is_SlowDeep, int m, int& num_tiles,
                                std::optional<std::pair<double, double>> zenith_pos,
                                std::vector<std::pair<double, double>> init_poses) {

    // std::vector<double> ra, dec;
    std::vector<std::vector<int>> member_pixels;
    readTiles(tilefile, ranks, ra, dec, member_pixels);
 
    std::vector<double> pixel_probs;
    readHEALPixels(pro_map_file, pixel_probs, 64);
 
    computeTileProbs(member_pixels, pixel_probs, probability);

    //  Fill missing initial positions with random tile positions
    const int n = static_cast<int>(ra.size());  
 
    // Track which tile index each init_pose came from (for rank lookup)
    std::vector<int> init_pose_rank_source;
    init_pose_rank_source.reserve(m);
    // For user-supplied positions we don't have a tile index, so store -2, -1 is t.
    for (std::size_t i = 0; i < init_poses.size(); ++i)
        init_pose_rank_source.push_back(-2);
 
    if (static_cast<int>(init_poses.size()) < m && n > 0) {
        unsigned int fixed_seed = 1130;
        std::mt19937 generator(fixed_seed);
        std::uniform_int_distribution<int> distribution(0, n - 1);
 
        for (int i = static_cast<int>(init_poses.size()); i < m; ++i) {
            int random_index = distribution(generator);
            init_poses.push_back({ra[random_index], dec[random_index]});
            init_pose_rank_source.push_back(ranks[random_index]);
        }
    }
 
    double coverage = 0.99;
    int num_tile = getTopTileIndices(probability, ranks, ra, dec,
                                    member_pixels, coverage);
    std::cout << "num tiles after trim: " << num_tile << "\n";
 
    // const int n = static_cast<int>(ra.size());  // number of tiles after trim
 
    // // Track which tile index each init_pose came from (for rank lookup)
    // std::vector<int> init_pose_rank_source;
    // init_pose_rank_source.reserve(m);
 
    // Compute dwell times for tile nodes
    dwell_times.resize(probability.size());
 
    double ra_zenith  = ra[0];
    double dec_zenith = dec[0];
    if (zenith_pos.has_value()) {
        ra_zenith  = zenith_pos->first;
        dec_zenith = zenith_pos->second;
    }
 
    for (std::size_t i = 0; i < probability.size(); ++i) {
        double zenithAngle = angularDistance(ra[i], dec[i],
                                            ra_zenith, dec_zenith);
        zenithAngle = std::min(zenithAngle, 85.0);
        dwell_times[i] = dwell_time_airmass(zenithAngle, dwell_zenith);
    }
 
    num_tiles = static_cast<int>(ra.size());
 

    //  Append start nodes (one per telescope / init position)
    //     Each start node has prize = 0 and dwell = 0.
    std::vector<int> start_indices;
    start_indices.reserve(m);
 
    for (int i = 0; i < static_cast<int>(init_poses.size()); ++i) {
        int idx = static_cast<int>(ra.size());
        start_indices.push_back(idx);
 
        ra.push_back(init_poses[i].first);
        dec.push_back(init_poses[i].second);
        probability.push_back(0.0);
        dwell_times.push_back(0.0);
        ranks.push_back(init_pose_rank_source[i]);
 
        // // Assign a rank: use the source tile's rank if we have one,
        // // otherwise use a sentinel (-1).
        // if (init_pose_rank_source[i] >= 0)
        //     ranks.push_back(ranks[init_pose_rank_source[i]]);
        // else
        //     ranks.push_back(-1);
    }

    compute_slew_costs(ra, dec, costs, w_max, w_acc, settle_time, is_SlowDeep);
 
    return start_indices;
}