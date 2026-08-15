#include "ReadData.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <stdexcept>

// #include <highfive/H5File.hpp>
// #include <highfive/H5DataSet.hpp>
// #include <highfive/H5DataSpace.hpp>


//Read Tiling Data
void readTiles(const std::string& filename, 
             std::vector<int>& ranks, 
             std::vector<double>& ras, 
             std::vector<double>& decs, 
             std::vector<std::vector<int>>& member_pixels) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    // Skip the header line
    if (!std::getline(file, line)) {
        std::cerr << "Error reading header" << std::endl;
        return;
    }

    while (std::getline(file, line)) {
        if (line.empty()) continue; // Skip empty lines

        std::istringstream ss(line);
        std::string token;

        // Read ID (rank)
        if (!std::getline(ss, token, ',')) continue;
        int rank;
        try {
            rank = std::stoi(token);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing rank: " << token << std::endl;
            continue;
        }
        ranks.push_back(rank);

        // Read RA
        if (!std::getline(ss, token, ',')) continue;
        double ra;
        try {
            ra = std::stod(token);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing RA: " << token << std::endl;
            continue;
        }
        ras.push_back(ra);

        // Read DEC
        if (!std::getline(ss, token, ',')) continue;
        double dec;
        try {
            dec = std::stod(token);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing DEC: " << token << std::endl;
            continue;
        }
        decs.push_back(dec);

        // Read the remaining line as HEALPixels (space-separated)
        std::string heal_str;
        std::getline(ss, heal_str);
        std::istringstream heal_ss(heal_str);
        std::vector<int> pixels;
        int pixel;
        while (heal_ss >> pixel) {
            pixels.push_back(pixel);
        }
        member_pixels.push_back(pixels);
    }

    file.close();

}



//Read Probability Map - csv
bool readHEALPixelsCSV(const std::string& filename,
                       std::vector<double>& pixel_probs,
                       int nside) {
    const int pixel_count = 12 * nside * nside;
    pixel_probs.assign(pixel_count, 0.0);

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << "\n";
        return false;
    }

    auto trim = [](const std::string& s) -> std::string {
        size_t b = 0, e = s.size();
        while (b < e && std::isspace((unsigned char)s[b])) ++b;
        while (e > b && std::isspace((unsigned char)s[e - 1])) --e;
        return s.substr(b, e - b);
    };

    std::vector<unsigned char> seen(pixel_count, 0);
    std::string line;
    size_t lineno = 0, loaded = 0, dup = 0, oob = 0, bad = 0, nan_inf = 0;

    bool header_skipped = false;

    while (std::getline(file, line)) {
        ++lineno;
        line = trim(line);
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        // Skip exactly one header line
        if (!header_skipped) {
            header_skipped = true;
            continue;
        }

        // Parse CSV: pix_id,ra_deg,dec_deg,probability
        std::istringstream ss(line);
        std::string token;

        int pix_id;
        double ra, dec, prob;

        try {
            if (!std::getline(ss, token, ',')) throw std::runtime_error("missing pix_id");
            pix_id = std::stoi(trim(token));

            if (!std::getline(ss, token, ',')) throw std::runtime_error("missing ra");
            ra = std::stod(trim(token)); (void)ra;

            if (!std::getline(ss, token, ',')) throw std::runtime_error("missing dec");
            dec = std::stod(trim(token)); (void)dec;

            if (!std::getline(ss, token, ',')) throw std::runtime_error("missing prob");
            prob = std::stod(trim(token));
        } catch (...) {
            std::cerr << "Parse error at line " << lineno << ": " << line << "\n";
            ++bad;
            continue;
        }

        if (pix_id < 0 || pix_id >= pixel_count) {
            std::cerr << "Out-of-bounds pix_id " << pix_id
                      << " at line " << lineno << " (pixel_count=" << pixel_count << ")\n";
            ++oob;
            continue;
        }

        if (!std::isfinite(prob)) {
            std::cerr << "Non-finite prob for pix_id " << pix_id
                      << " at line " << lineno << " -> skipped\n";
            ++nan_inf;
            continue;
        }

        if (seen[pix_id]) ++dup;   // keep last value
        seen[pix_id] = 1;
        pixel_probs[pix_id] = prob;
        ++loaded;
    }

    const size_t seen_count = std::accumulate(seen.begin(), seen.end(), size_t{0});
    const size_t missing = (size_t)pixel_count > seen_count ? ((size_t)pixel_count - seen_count) : 0;
    const double sum_prob = std::accumulate(pixel_probs.begin(), pixel_probs.end(), 0.0);

    std::cerr << "Loaded " << loaded << " entries "
              << "(unique=" << seen_count
              << ", duplicates=" << dup
              << ", missing=" << missing
              << ", out-of-bound=" << oob
              << ", non-finite=" << nan_inf
              << ", parse_errors=" << bad << ")\n";
    std::cerr << "Sum(prob) = " << sum_prob << "\n";

    return true;
}



//Read Probability Map - txt(space seperated)
bool readHEALPixelsTXT(const std::string& filename,
                    std::vector<double>& pixel_probs,
                    int nside) {
    int pixel_count = 12 * nside * nside;
    pixel_probs.assign(pixel_count, 0.0);

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << "\n";
        return false;
    }

    std::vector<unsigned char> seen(pixel_count, 0);
    std::string line;
    size_t lineno = 0, loaded = 0, dup = 0, oob = 0, bad = 0, nan_inf = 0;

    while (std::getline(file, line)) {
        ++lineno;
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        std::istringstream ss(line);
        int pix_id;
        double lat, lon, prob;
        if (!(ss >> pix_id >> lat >> lon >> prob)) {
            std::cerr << "Parse error at line " << lineno << ": " << line << "\n";
            ++bad;
            continue;
        }

        if (pix_id < 0 || pix_id >= pixel_count) {
            std::cerr << "Out-of-bounds pix_id " << pix_id
                      << " at line " << lineno << " (pixel_count=" << pixel_count << ")\n";
            ++oob;
            continue;
        }

        if (!std::isfinite(prob)) {
            std::cerr << "Non-finite prob for pix_id " << pix_id
                      << " at line " << lineno << " -> skipped\n";
            ++nan_inf;
            continue;
        }

        if (seen[pix_id]) ++dup; // keep last value
        seen[pix_id] = 1;
        pixel_probs[pix_id] = prob;
        ++loaded;
    }

    const size_t seen_count = std::accumulate(seen.begin(), seen.end(), size_t{0});
    const size_t missing = pixel_count > seen_count ? (pixel_count - seen_count) : 0;
    const double finite_sum = std::accumulate(pixel_probs.begin(), pixel_probs.end(), 0.0);

    std::cerr << "Loaded " << loaded << " entries "
              << "(unique=" << seen_count
              << ", duplicates=" << dup
              << ", missing=" << missing
              << ", out bound=" << oob
              << ", non-finite=" << nan_inf
              << ", parse_errors=" << bad << ")\n";
    std::cerr << "Sum(prob) over finite entries = " << finite_sum << "\n";

    return true;
}


// bool readHEALPixelsHDF5_Compound(const std::string& filename,
//                                   std::vector<double>& pixel_probs,
//                                   int nside,
//                                   const std::string& dataset_name = "likelihood_map") {
//     int pixel_count = 12 * nside * nside;
//     pixel_probs.assign(pixel_count, 0.0);
//     try {
//         HighFive::File file(filename, HighFive::File::ReadOnly);

//         struct MapEntry {
//             int64_t pixel;
//             double lat;
//             double lon;
//             double probability;
//         };

//         auto dataset = file.getDataSet(dataset_name);
//         size_t n = dataset.getSpace().getDimensions()[0];
//         std::vector<MapEntry> entries(n);

//         hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(MapEntry));
//         H5Tinsert(memtype, "pixel",       offsetof(MapEntry, pixel),       H5T_NATIVE_INT64);
//         H5Tinsert(memtype, "lat",         offsetof(MapEntry, lat),         H5T_NATIVE_DOUBLE);
//         H5Tinsert(memtype, "lon",         offsetof(MapEntry, lon),         H5T_NATIVE_DOUBLE);
//         H5Tinsert(memtype, "probability", offsetof(MapEntry, probability), H5T_NATIVE_DOUBLE);

//         herr_t status = H5Dread(dataset.getId(), memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, entries.data());
//         H5Tclose(memtype);

//         if (status < 0) {
//             std::cerr << "H5Dread failed for dataset " << dataset_name << " in " << filename << "\n";
//             return false;
//         }

//         std::vector<unsigned char> seen(pixel_count, 0);
//         size_t loaded = 0, dup = 0, oob = 0, nan_inf = 0;

//         for (size_t i = 0; i < entries.size(); ++i) {
//             int pix_id = static_cast<int>(entries[i].pixel);
//             double prob = entries[i].probability;

//             if (pix_id < 0 || pix_id >= pixel_count) {
//                 std::cerr << "Out-of-bounds pix_id " << pix_id
//                           << " at index " << i << " (pixel_count=" << pixel_count << ")\n";
//                 ++oob;
//                 continue;
//             }
//             if (!std::isfinite(prob)) {
//                 std::cerr << "Non-finite prob for pix_id " << pix_id
//                           << " at index " << i << " -> skipped\n";
//                 ++nan_inf;
//                 continue;
//             }
//             if (seen[pix_id]) ++dup;
//             seen[pix_id] = 1;
//             pixel_probs[pix_id] = prob;
//             ++loaded;
//         }

//         const size_t seen_count = std::accumulate(seen.begin(), seen.end(), size_t{0});
//         const size_t missing = pixel_count > seen_count ? (pixel_count - seen_count) : 0;
//         const double finite_sum = std::accumulate(pixel_probs.begin(), pixel_probs.end(), 0.0);
//         std::cerr << "Loaded " << loaded << " entries "
//                   << "(unique=" << seen_count
//                   << ", duplicates=" << dup
//                   << ", missing=" << missing
//                   << ", out_of_bounds=" << oob
//                   << ", non-finite=" << nan_inf << ")\n";
//         std::cerr << "Sum(prob) over finite entries = " << finite_sum << "\n";
//         return true;

//     } catch (const HighFive::Exception& e) {
//         std::cerr << "HDF5 error reading file " << filename << ": " << e.what() << "\n";
//         return false;
//     }
// }



// //Read Probability Map - HDF5 (using HighFive)
// bool readHEALPixelsHDF5(const std::string& filename,
//                         std::vector<double>& pixel_probs,
//                         int nside,
//                         const std::string& pix_id_dataset,
//                         const std::string& prob_dataset) {
//     int pixel_count = 12 * nside * nside;
//     pixel_probs.assign(pixel_count, 0.0);

//     try {
//         HighFive::File file(filename, HighFive::File::ReadOnly);

//         // Read datasets
//         std::vector<int> pix_ids;
//         std::vector<double> probs;

//         file.getDataSet(pix_id_dataset).read(pix_ids);
//         file.getDataSet(prob_dataset).read(probs);

//         if (pix_ids.size() != probs.size()) {
//             std::cerr << "Error: pix_id and prob datasets have different sizes ("
//                       << pix_ids.size() << " vs " << probs.size() << ")\n";
//             return false;
//         }

//         std::vector<unsigned char> seen(pixel_count, 0);
//         size_t loaded = 0, dup = 0, oob = 0, nan_inf = 0;

//         for (size_t i = 0; i < pix_ids.size(); ++i) {
//             int pix_id = pix_ids[i];
//             double prob = probs[i];

//             if (pix_id < 0 || pix_id >= pixel_count) {
//                 std::cerr << "Out-of-bounds pix_id " << pix_id
//                           << " at index " << i << " (pixel_count=" << pixel_count << ")\n";
//                 ++oob;
//                 continue;
//             }
//             if (!std::isfinite(prob)) {
//                 std::cerr << "Non-finite prob for pix_id " << pix_id
//                           << " at index " << i << " -> skipped\n";
//                 ++nan_inf;
//                 continue;
//             }
//             if (seen[pix_id]) ++dup; // keep last value
//             seen[pix_id] = 1;
//             pixel_probs[pix_id] = prob;
//             ++loaded;
//         }

//         const size_t seen_count = std::accumulate(seen.begin(), seen.end(), size_t{0});
//         const size_t missing = pixel_count > seen_count ? (pixel_count - seen_count) : 0;
//         const double finite_sum = std::accumulate(pixel_probs.begin(), pixel_probs.end(), 0.0);

//         std::cerr << "Loaded " << loaded << " entries "
//                   << "(unique=" << seen_count
//                   << ", duplicates=" << dup
//                   << ", missing=" << missing
//                   << ", out_of_bounds=" << oob
//                   << ", non-finite=" << nan_inf << ")\n";
//         std::cerr << "Sum(prob) over finite entries = " << finite_sum << "\n";

//         return true;

//     } catch (const HighFive::Exception& e) {
//         std::cerr << "HDF5 error reading file " << filename << ": " << e.what() << "\n";
//         return false;
//     }
// }


// Helper to extract and lowercase the file extension
std::string getExtension(const std::string& filename) {
    size_t dot_pos = filename.rfind('.');
    if (dot_pos == std::string::npos || dot_pos == filename.length() - 1) {
        return "";
    }
    std::string ext = filename.substr(dot_pos + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return ext;
}

// map reading wrapper
bool readHEALPixels(const std::string& filename,
                    std::vector<double>& pixel_probs,
                    int nside) {
    std::string ext = getExtension(filename);

    if (ext == "csv") {
        return readHEALPixelsCSV(filename, pixel_probs, nside);
    } else if (ext == "txt" || ext == "dat") {
        return readHEALPixelsTXT(filename, pixel_probs, nside);
    // } else if (ext == "h5" || ext == "hdf5" || ext == "he5") {
    //     return readHEALPixelsHDF5(filename, pixel_probs, nside);
    //     // return readHEALPixelsHDF5_Compound(filename, pixel_probs,nside);
    } else {
        std::cerr << "Unknown file extension '" << ext
                  << "' for file: " << filename << "\n"
                  << "Supported extensions: .csv, .txt, .dat\n";
                //   << "Supported extensions: .csv, .txt, .dat, .h5, .hdf5, .he5\n";
        return false;
    }
}


// //Read Probability Map
// void readHEALPixels(const std::string& filename, std::vector<int>& pixel_ids, 
//                 std::vector<double>& pixel_probs) {
//     pixel_ids.clear();
//     pixel_probs.clear();

//     std::ifstream file(filename);
//     if (!file.is_open()) {
//         std::cerr << "Error opening file: " << filename << std::endl;
//         return;
//     }

//     std::vector<std::pair<int, double>> temp;
//     std::string line;
//     while (std::getline(file, line)) {
//         if (line.empty() || line[0] == '#') continue; // Skip empty or comment lines

//         std::istringstream ss(line);
//         int pix_id;
//         double lat, lon, prob;
//         if (ss >> pix_id >> lat >> lon >> prob) {
//             temp.push_back({pix_id, prob});
//         } else {
//             std::cerr << "Error parsing line: " << line << std::endl;
//         }
//     }

//     file.close();

//     // Sort by pix_id, ensure idx align with pix_id
//     std::sort(temp.begin(), temp.end());

//     for (const auto& pr : temp) {
//         pixel_ids.push_back(pr.first);
//         pixel_probs.push_back(pr.second);
//     }

//     // Check if pixel_ids align with indices 
//     for (size_t i = 0; i < pixel_ids.size(); ++i) {
//         std::cout << pixel_ids[i] << " " << i << "\n";
//         if (pixel_ids[i] != static_cast<int>(i)) {
//             pixel_ids.clear();
//             pixel_probs.clear();
//             throw std::runtime_error("Error: Pixel IDs do not align with indices at position " + std::to_string(i) +
//                                     " (expected " + std::to_string(i) + ", got " + std::to_string(pixel_ids[i]) + ")");
//         }
//     }
// }


// //Read Probability Map
// bool readHEALPixels(const std::string& filename,
//                     std::vector<int>& pixel_ids,
//                     std::vector<double>& pixel_probs,
//                     int pixel_count) {
//     pixel_ids.resize(pixel_count);
//     std::iota(pixel_ids.begin(), pixel_ids.end(), 0);
//     pixel_probs.assign(pixel_count, 0.0);

//     std::ifstream file(filename);
//     if (!file.is_open()) {
//         std::cerr << "Error opening file: " << filename << "\n";
//         return false;
//     }

//     std::vector<unsigned char> seen(pixel_count, 0);
//     std::string line;
//     size_t lineno = 0, loaded = 0, dup = 0, oob = 0, bad = 0, nan_inf = 0;

//     while (std::getline(file, line)) {
//         ++lineno;
//         if (line.empty()) continue;
//         if (line[0] == '#') continue;

//         std::istringstream ss(line);
//         int pix_id;
//         double lat, lon, prob;
//         if (!(ss >> pix_id >> lat >> lon >> prob)) {
//             std::cerr << "Parse error at line " << lineno << ": " << line << "\n";
//             ++bad;
//             continue;
//         }

//         if (pix_id < 0 || pix_id >= pixel_count) {
//             std::cerr << "Out-of-bounds pix_id " << pix_id
//                       << " at line " << lineno << " (pixel_count=" << pixel_count << ")\n";
//             ++oob;
//             continue;
//         }

//         if (!std::isfinite(prob)) {
//             std::cerr << "Non-finite prob for pix_id " << pix_id
//                       << " at line " << lineno << " -> skipped\n";
//             ++nan_inf;
//             continue;
//         }

//         if (seen[pix_id]) ++dup; // keep last value
//         seen[pix_id] = 1;
//         pixel_probs[pix_id] = prob;
//         ++loaded;
//     }

//     const size_t seen_count = std::accumulate(seen.begin(), seen.end(), size_t{0});
//     const size_t missing = pixel_count > seen_count ? (pixel_count - seen_count) : 0;
//     const double finite_sum = std::accumulate(pixel_probs.begin(), pixel_probs.end(), 0.0);

//     std::cerr << "Loaded " << loaded << " entries "
//               << "(unique=" << seen_count
//               << ", duplicates=" << dup
//               << ", missing=" << missing
//               << ", out bound=" << oob
//               << ", non-finite=" << nan_inf
//               << ", parse_errors=" << bad << ")\n";
//     std::cerr << "Sum(prob) over finite entries = " << finite_sum << "\n";

//     return true;
// }

void computeTileProbs(const std::vector<std::vector<int>>& member_pixels,
                    const std::vector<double>& pixel_probs,
                    std::vector<double>& tile_probs) {
    tile_probs.resize(member_pixels.size(), 0.0);

    double sum = std::accumulate(pixel_probs.begin(), pixel_probs.end(), 0.0);
    std::cout << "init sum: " << sum << "\n";
    std::cout << "num pixels: " << pixel_probs.size() << "\n";


    for (size_t i = 0; i < member_pixels.size(); ++i) {
        tile_probs[i] = 0.0;
        for (int p : member_pixels[i]) {
                if(p > pixel_probs.size()) 
                    throw std::runtime_error("Error: pixel counts do not match tiles!");
                tile_probs[i] += pixel_probs[p];

        }
    }
}


int getTopTileIndices(std::vector<double>& tile_probs,
                      std::vector<int>& ranks,
                      std::vector<double>& ras,
                      std::vector<double>& decs,
                      std::vector<std::vector<int>>& member_pixels,
                      double coverage) {
    std::cout << "tile num: " << member_pixels.size() << "\n";
    std::cout << "total tile prob: " << std::accumulate(tile_probs.begin(), tile_probs.end(), 0.0) << "\n";
    if (tile_probs.empty()) return 0;

    // // Compute total sum to treat coverage as a fraction
    // double total_sum = 0.0;
    // for (double p : tile_probs) total_sum += p;
    // double threshold = coverage * total_sum;

    std::vector<size_t> indices(tile_probs.size());
    std::iota(indices.begin(), indices.end(), 0);

    // Sort indices in descending order of tile_probs
    std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
        return tile_probs[a] > tile_probs[b];
    });

    // Create sorted versions
    std::vector<double> sorted_tile_probs(tile_probs.size());
    std::vector<int> sorted_ranks(ranks.size());
    std::vector<double> sorted_ras(ras.size());
    std::vector<double> sorted_decs(decs.size());
    std::vector<std::vector<int>> sorted_member_pixels(member_pixels.size());

    for (size_t i = 0; i < indices.size(); ++i) {
        sorted_tile_probs[i] = tile_probs[indices[i]];
        sorted_ranks[i] = ranks[indices[i]];
        sorted_ras[i] = ras[indices[i]];
        sorted_decs[i] = decs[indices[i]];
        sorted_member_pixels[i] = std::move(member_pixels[indices[i]]);
    }

    // Assign back to originals
    tile_probs = std::move(sorted_tile_probs);
    ranks = std::move(sorted_ranks);
    ras = std::move(sorted_ras);
    decs = std::move(sorted_decs);
    member_pixels = std::move(sorted_member_pixels);

    // Accumulate until >= threshold using the now sorted tile_probs
    double eps = 10e-5;
    double cum_sum = 0.0;
    size_t n = 0;
    for (; n < tile_probs.size(); ++n) {
        // if (tile_probs[n] < eps) break;
        cum_sum += tile_probs[n];
        if (cum_sum >= coverage) break;
    }

    // Include the last one if we broke out
    if (n < tile_probs.size()) ++n;

    // Resize all to keep only the top N
    tile_probs.resize(n);
    ranks.resize(n);
    ras.resize(n);
    decs.resize(n);
    member_pixels.resize(n);
    std::cout << "total tile prob after trim: " << std::accumulate(tile_probs.begin(), tile_probs.end(), 0.0) << "\n";

    return static_cast<int>(n);
}

int getNonZeroTileIndices(std::vector<double>& tile_probs,
                      std::vector<int>& ranks,
                      std::vector<double>& ras,
                      std::vector<double>& decs,
                      std::vector<std::vector<int>>& member_pixels) {
    std::cout << "tile num: " << member_pixels.size() << "\n";
    std::cout << "total tile prob: " << std::accumulate(tile_probs.begin(), tile_probs.end(), 0.0) << "\n";
    if (tile_probs.empty()) return 0;

    // Find indices of tiles with prob > 0
    std::vector<size_t> keep;
    for (size_t i = 0; i < tile_probs.size(); ++i) {
        if (tile_probs[i] > 0.0) keep.push_back(i);
    }

    // Sort kept indices in descending order of tile_probs
    std::sort(keep.begin(), keep.end(), [&](size_t a, size_t b) {
        return tile_probs[a] > tile_probs[b];
    });

    size_t n = keep.size();

    // Create compacted versions
    std::vector<double> new_tile_probs(n);
    std::vector<int> new_ranks(n);
    std::vector<double> new_ras(n);
    std::vector<double> new_decs(n);
    std::vector<std::vector<int>> new_member_pixels(n);

    for (size_t i = 0; i < n; ++i) {
        new_tile_probs[i] = tile_probs[keep[i]];
        new_ranks[i] = ranks[keep[i]];
        new_ras[i] = ras[keep[i]];
        new_decs[i] = decs[keep[i]];
        new_member_pixels[i] = std::move(member_pixels[keep[i]]);
    }

    tile_probs = std::move(new_tile_probs);
    ranks = std::move(new_ranks);
    ras = std::move(new_ras);
    decs = std::move(new_decs);
    member_pixels = std::move(new_member_pixels);

    std::cout << "tiles with prob > 0: " << n << "\n";
    std::cout << "total tile prob after trim: " << std::accumulate(tile_probs.begin(), tile_probs.end(), 0.0) << "\n";

    return static_cast<int>(n);
}



// read rank, ra, dec and prize from Shaon tiling file
void read_data_from_file(const std::string& filename, std::vector<double>& probability, 
                        std::vector<int>& ranks, std::vector<double>& ra, std::vector<double>& dec) {

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line, token;
    getline(file, line);  // Skip header

    while (getline(file, line)) {
        std::stringstream ss(line);
        std::vector<std::string> row;

        while (getline(ss, token, ',')) {
            row.push_back(token);
        }

        ranks.push_back(std::stoi(row[0])); 
        ra.push_back(std::stod(row[2]));  
        dec.push_back(std::stod(row[3]));  
        probability.push_back(std::stod(row[4]));  
    }
    file.close();

}


