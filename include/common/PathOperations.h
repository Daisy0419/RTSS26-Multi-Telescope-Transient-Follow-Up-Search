#pragma once
#include <vector>
#include <tuple>
#include <cstddef>  

// **** Single path operation helplers****

double pathCost(const std::vector<std::vector<double>>& costs, const std::vector<int>& path);
double pathPrize(const std::vector<double>& prizes, const std::vector<int>& path);

// Calculate delta cost when removing a node
double removalDeltaCost(const std::vector<std::vector<double>>& costs,
                         const std::vector<int>& path,
                         size_t node_index);

// Calculate delta cost when inserting a node
double insertionDeltaCost(const std::vector<std::vector<double>>& costs,
                         const std::vector<int>& path,
                         size_t node, size_t insert_index);

// Remove duplicate nodes in a path based, preserving the first occurance
void removeDuplicates(const std::vector<std::vector<double>>& costs,
                            const std::vector<double>& prizes,
                            std::vector<int>& path);

// Remove duplicate nodes in a path based on prize/delta_cost ratio
void removeDuplicatesByRatio(const std::vector<std::vector<double>>& costs,
                            const std::vector<double>& prizes,
                            std::vector<int>& path);

// Drop nodes with lowest prize until budget is met
void prunePathByPrize(const std::vector<std::vector<double>>& costs,
                         const std::vector<double>& prizes,
                         std::vector<int>& path,
                         double budget,
                         double& total_cost);

// Drop nodes with lowest prize/delta_cost ratio until budget is met
void prunePathByRatio(const std::vector<std::vector<double>>& costs,
                  const std::vector<double>& prizes,
                  std::vector<int>& path,
                  double budget,
                  double& total_cost);

// Find cheapest insertion 'pos' (between pos and pos+1), pos ∈ [0, n-2] ensure s-t path
std::pair<int, double> cheapestInsertionOnePath(const std::vector<std::vector<double>>& costs,
                                        const std::vector<int>& path, int node);

// Insert nodes with highest possible value in cheapest position until budget exhausts
void insertHighValueCheapest(std::vector<int>& path,
                                const std::vector<std::vector<double>>& cost,
                                const std::vector<double>& prize,
                                double budget,
                                double& current_cost,
                                std::vector<bool>& visited);

// Insert nodes with highest possible value before end node until budget exhausts
void insertHighValueEnd(std::vector<int>& path,
                                  const std::vector<std::vector<double>>& costs,
                                  const std::vector<double>& prizes,
                                  double budget,
                                  double& current_cost,
                                  std::vector<bool>& visited);


// Insert nodes with global best prize/delta_cost ratio until budget exhausts
void insertBestRatio(std::vector<int>& path,
                     const std::vector<std::vector<double>>& cost,
                     const std::vector<double>& prize,
                     double budget,
                     double& current_cost,
                     std::vector<bool>& visited);


void repairPath(const std::vector<std::vector<double>>& costs,
                const std::vector<double>& prizes,
                std::vector<int>& path,
                double budget);

void repairPath(const std::vector<std::vector<double>>& costs,
                const std::vector<double>& prizes,
                std::vector<int>& path,
                double budget, 
                std::vector<bool>& visited);

std::vector<int> repair_path(const std::vector<int>& path,
               const std::vector<std::vector<double>>& costs,
               const std::vector<double>& prizes,
               double budget, std::vector<bool>& visited);




// *****Team paths operation helplers*****
std::vector<double> multiPathsCost(const std::vector<std::vector<double>>& costs,
                                   const std::vector<std::vector<int>>& paths);

std::vector<double> multiPathsPrize(const std::vector<double>& prizes,
                       const std::vector<std::vector<int>>& paths);
                       
double multiPathsMaxCost(const std::vector<std::vector<double>>& costs,
                       const std::vector<std::vector<int>>& paths);

double multiPathsSumPrize(const std::vector<double>& prizes,
                       const std::vector<std::vector<int>>& paths);

std::vector<bool> getVisited(const std::vector<std::vector<int>>& paths, int nNodes);                 

// Remove duplicates across paths; keep the FIRST ownership
void removeInterPathDuplicatesKeepFirst(const std::vector<std::vector<double>>& costs,
                                        std::vector<std::vector<int>>& paths);

// Keep the single best (prize / delta_cost) occurrence across all paths per node
void removeInterPathDuplicatesByRatio(const std::vector<std::vector<double>>& costs,
                                      const std::vector<double>& prizes,
                                      std::vector<std::vector<int>>& paths);

// Greedy global insertion: pick (city, route, pos) with best prize/delta ratio each time.
void insertHighValueCheapestMulti(const std::vector<std::vector<double>>& costs,
                                        const std::vector<double>& prizes,
                                        std::vector<std::vector<int>>& paths,
                                        double budget_per_route,
                                        int s, int t);

void insertHighValueCheapestMulti( const std::vector<std::vector<double>>& costs,
                                    const std::vector<double>& prizes,
                                    std::vector<std::vector<int>>& paths,
                                    double budget_per_route,
                                    const std::vector<int>& start_indices,
                                    int t);


 // multi-route repair
void repairMultiPaths(const std::vector<std::vector<double>>& costs,
                      const std::vector<double>& prizes,
                      std::vector<std::vector<int>>& paths,
                      double budget_per_route,
                      int s, int t);   
void repairMultiPaths(const std::vector<std::vector<double>>& costs,
                      const std::vector<double>& prizes,
                      std::vector<std::vector<int>>& paths,
                      double budget_per_route,
                      const std::vector<int>& start_indices,
                      int t);                                    


std::tuple<int, int, double> cheapestInsertionAll(const std::vector<std::vector<double>>& costs,
                                                       const std::vector<std::vector<int>>& paths,
                                                       const std::vector<double>& pathCosts,
                                                       const std::vector<double>& prizes,
                                                       double budget, 
                                                       int node);



bool verify_routes(const std::vector<std::vector<double>>& costs,
                   const std::vector<double>& prizes,
                   int s, int t,
                   const std::vector<std::vector<int>>& routes,
                   double Budget,
                   int m_expected = -1,  
                   double tol = 1e-10);