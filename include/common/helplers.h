#pragma once

#include "ReadData.h"

#include <vector>

//random number generater
double random_double(double min, double max);
int random_int(int min, int max);
std::vector<int> unique_random_ints(int min, int max, int n);

// //2-opt
// void fix_cross(std::vector<int>& path, const std::vector<std::vector<double>>& costs);
// void fix_cross_st_path(std::vector<int>& path, const std::vector<std::vector<double>>& costs);


// //compute cost of a path
// double calculate_path_cost(const std::vector<int>& path, const std::vector<std::vector<double>>& costs);


// double calculate_paths_prize(const std::vector<std::vector<int>>& paths, 
//                              const std::vector<double>& prizes);



/***************************************** */
bool verify_paths(const std::vector<std::vector<double>>& costs,
                  const std::vector<std::vector<int>>& routes,
                  const std::vector<int>& s,   // per-route starts (duplicates allowed)
                  const int t, const int N,
                  double Budget,
                  int m_expected = -1,
                  double tol = 1e-9);

bool verify_paths_any_order(const std::vector<std::vector<double>>& costs,
                            const std::vector<std::vector<int>>& routes,
                            const std::vector<int>& s,   // allowed starts (multiset semantics)
                            int t, int N,
                            double Budget,
                            int m_expected = -1,
                            double tol = 1e-9);

bool verify_routes(const Instance& I,
                   const std::vector<std::vector<int>>& routes,
                   double Budget,
                   int m_expected = -1,                 // −1 don’t enforce
                   double tol = 1e-9);

void remove_st (std::vector<int>& path);

double print_path(const std::vector<std::vector<double>>& costs, const std::vector<double>& probability, 
                    std::vector<int> rank, std::vector<int> best_path, double padding);

double print_paths(const std::vector<std::vector<double>>& costs,
                   const std::vector<double>& probability,
                   std::vector<int> rank,                      // optional; can be empty
                   const std::vector<std::vector<int>>& paths, // multiple paths
                   double padding);



void save_result(const std::string& filename, 
                  const std::string& method, 
                  const std::string& data, 
                  double budget, 
                  double slew_rate, 
                  int npaths,
                  const double prize,
                  double elapsed_time,
                  bool is_valid = true);

void save_result_with_paths(const std::string& filename, 
                            const std::string& method, 
                            const std::string& data, 
                            const std::vector<std::vector<int>>& paths,
                            const std::vector<int>& rank,
                            double budget, 
                            double w_max, double w_acc,
                            int npaths,
                            int nTiles,
                            const double prize,
                            const double best_prize,
                            double tiling_time,
                            double planning_time,
                            bool IsValid=true);


void save_result_with_paths2(const std::string& filename, 
                            const std::string& method, 
                            const std::string& data, 
                            const std::vector<std::vector<int>>& paths,
                            double budget, 
                            double w_max, double w_acc,
                            int npaths,
                            int nTiles,
                            const double prize,
                            const double best_prize,
                            double tiling_time,
                            double planning_time,
                            bool IsValid=true);


void save_result_without_paths(const std::string& filename, 
                            const std::string& method, 
                            const std::string& data, 
                            double budget, 
                            double w_max, double w_acc,
                            double num_tiles,
                            double covered_tiles,
                            double npaths,
                            double npaths_bound,
                            double tiling_time,
                            double planning_time,
                            const std::vector<std::vector<int>>& paths,
                            bool IsValid);
                            
// void save_result_without_paths(const std::string& filename, 
//                             const std::string& method, 
//                             const std::string& data, 
//                             double budget, 
//                             double w_max, double w_acc,
//                             double npaths,
//                             double npaths_bound,
//                             double tiling_time,
//                             double planning_time,
//                             const std::vector<std::vector<int>>& paths,
//                             bool IsValid=true);


// void save_result_without_paths(const std::string& filename, 
//                             const std::string& method, 
//                             const std::string& data, 
//                             double budget, 
//                             double w_max, double w_acc,
//                             double num_tiles,
//                             double covered_tiles,
//                             double npaths,
//                             double tiling_time,
//                             double planning_time,
//                             bool IsValid=true);

double print_path(const std::vector<std::vector<double>>& costs, const std::vector<double>& probability, 
                    std::vector<int> rank, std::vector<int> best_path, double padding);

std::pair<int,int> extractSingleStartSubgraph(
        const std::vector<std::vector<double>>& costs_full,
        const std::vector<int>& start_indices,
        int start_to_keep,
        int end_index,
        std::vector<std::vector<double>>& costs_out);


std::vector<double> extractSingleStartProbability(
        const std::vector<double>& probability_full,
        const std::vector<int>& start_indices);