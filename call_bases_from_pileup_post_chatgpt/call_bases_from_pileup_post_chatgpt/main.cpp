#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <regex>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <chrono>

std::vector<std::string> header = {"A", "C", "G", "T", "-"};

std::vector<double> get_q(const std::string& ref, const std::string& match, double e) {
    std::vector<double> qualities(5, e/3);
    std::string updatedMatch = match;
    qualities[4] = 0.0;
    if (match == "," || match == ".") {
        updatedMatch = ref;
    }
    auto it = std::find(header.begin(), header.end(), updatedMatch);
    if (it != header.end()) {
        qualities[std::distance(header.begin(), it)] = 1 - e;
        return qualities;
    }
    if (match == "*") {
        return std::vector<double>(5, 1);
    }
    return std::vector<double>(5, NAN);
}

std::vector<int> get_all_comparisons(const std::vector<std::string>& list1, const std::vector<std::string>& list2) {
    std::vector<int> comparisons;
    comparisons.reserve(list1.size());
    for (size_t i = 0; i < list1.size(); ++i) {
        if (list1[i] == "N" || list2[i] == "N" || list2[i] == "-") {
            comparisons.push_back(-1);
        } else {
            comparisons.push_back(list1[i] == list2[i]);
        }
    }
    return comparisons;
}

std::vector<std::string> get_agreements(const std::vector<std::string>& list1, const std::vector<std::string>& list2, const std::vector<int>& pos_list) {
    std::vector<int> ag_all;
    ag_all.reserve(list1.size());
    for (size_t i = 0; i < list1.size(); ++i) {
        if (list1[i] == "N" || list2[i] == "N" || list2[i] == "-") {
            ag_all.push_back(-1);
        } else {
            ag_all.push_back(list1[i] == list2[i]);
        }
    }
    
    std::vector<int> ag;
    ag.reserve(ag_all.size());
    for (int a : ag_all) {
        if (a != -1) {
            ag.push_back(a);
        }
    }
    
    double mean_agreement = static_cast<double>(std::accumulate(ag.begin(), ag.end(), 0)) / ag.size();
    int sum_agreement = std::accumulate(ag.begin(), ag.end(), 0);
    int agreement_count = ag.size();
    int disagreement_count = static_cast<int>(ag.size()) - sum_agreement;
    
    return {std::to_string(mean_agreement), std::to_string(sum_agreement), std::to_string(agreement_count), std::to_string(disagreement_count)};
}

int main(int argc, char* argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    
    if (argc < 3) {
        std::cerr << "Usage: ./program_name pileup_file pileup_results_path" << std::endl;
        return 1;
    }
    
    std::string pileup_file = argv[1];
    std::string pileup_results_path = argv[2];
    
    std::string home_dir = "/Users/marniella/research/nielsen_lab/phylogenetic-base-calling/";
    std::string msa_dir = "ou_pruned/";
    
    std::ifstream f(home_dir + "justmn.fasta");
    std::string reference;
    int lineCount = 1;
    while (std::getline(f, reference)) {
        if (lineCount == 2) {
            break;
        }
        lineCount++;
    }
    f.close();

    std::ifstream f2(home_dir + "justou.fasta");
    std::string og_reference;
    int lineCount2 = 1;
    while (std::getline(f2, og_reference)) {
        if (lineCount2 == 2) {
            break;
        }
        lineCount2++;
    }
    f2.close();
    
    std::vector<int> msa_og_j;
    for (size_t j = 0; j < reference.size(); ++j) {
        if (reference[j] != '-') {
            msa_og_j.push_back(j);
        }
    }
    
    std::regex pattern1("-[0-9]+[ACGTNacgtn]+");
    std::regex pattern2("\\+[0-9]+[ACGTNacgtn*#]+");
    std::regex pattern3("\\$");
    std::regex pattern4("\\^.{1}");
    
    std::map<int, std::vector<double>> bases_prior_dict;
    std::ifstream f3(home_dir + msa_dir + "OU_PP_indexed.txt");
    std::string line;
    std::getline(f3, line); // Skip the first line
    while (std::getline(f3, line)) {
        std::istringstream iss(line);
        std::string msa_pos, prior_A, prior_C, prior_G, prior_T;
        iss >> msa_pos >> prior_A >> prior_C >> prior_G >> prior_T;
        bases_prior_dict[std::stoi(msa_pos)] = {std::stod(prior_A), std::stod(prior_C), std::stod(prior_G), std::stod(prior_T), 0.0};
    }
    f3.close();
    
    std::map<int, std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>> pileup_dict;
    std::ifstream f4(home_dir + pileup_results_path + pileup_file);
    while (std::getline(f4, line)) {
        std::istringstream iss(line);
        std::string name, pos, ref_base, n_reads, match_str, phred_list;
        iss >> name >> pos >> ref_base >> n_reads >> match_str >> phred_list;
        int position = std::stoi(pos) - 1;
        match_str = std::regex_replace(
                std::regex_replace(
                    std::regex_replace(
                        std::regex_replace(match_str, pattern4, ""), pattern3, ""),pattern2, ""),pattern1, "");
        
        std::vector<std::vector<double>> q_list;
        for (size_t i = 0; i < match_str.size(); ++i) {
            q_list.push_back(get_q(ref_base, std::string(1, match_str[i]), std::pow(10, -(static_cast<double>(phred_list[i]) - 33) / 10.0)));
        }
        
        std::vector<double> ll_list(q_list[0].size());

        for (size_t i = 0; i < q_list[0].size(); ++i) {
            double prod = 1.0;
            for (size_t j = 0; j < q_list.size(); ++j) {
                prod *= q_list[j][i];
            }
            ll_list[i] = std::isnan(prod) ? 0.0 : prod;
        }
        
        int msa_pos = msa_og_j[position];
        std::vector<double> pr_list = bases_prior_dict[msa_pos];
        std::vector<double> pp_list(ll_list.size(), 0.0);
        
        if ((std::accumulate(ll_list.begin(), ll_list.end(), 0.0) == 5.0) && (og_reference[msa_pos] == '-')) {
            pp_list = {0.0, 0.0, 0.0, 0.0, 1.0};
        } else {
            for (size_t i = 0; i < ll_list.size(); ++i) {
                pp_list[i] = pr_list[i] * ll_list[i];
            }
        }
        
        pileup_dict[position] = std::make_tuple(pr_list, ll_list, pp_list);
    }
    f4.close();
    
    std::vector<int> positions_idx;
    for (const auto& entry : pileup_dict) {
        positions_idx.push_back(entry.first);
    }
    std::sort(positions_idx.begin(), positions_idx.end());
    
    std::vector<std::string> pr, pp, ll, og;
    for (int pos : positions_idx) {
        pr.push_back(header[std::distance(std::get<0>(pileup_dict[pos]).begin(), std::max_element(std::get<0>(pileup_dict[pos]).begin(), std::get<0>(pileup_dict[pos]).end()))]);
        ll.push_back(header[std::distance(std::get<1>(pileup_dict[pos]).begin(), std::max_element(std::get<1>(pileup_dict[pos]).begin(), std::get<1>(pileup_dict[pos]).end()))]);
        pp.push_back(header[std::distance(std::get<2>(pileup_dict[pos]).begin(), std::max_element(std::get<2>(pileup_dict[pos]).begin(), std::get<2>(pileup_dict[pos]).end()))]);
        
        og.push_back(std::string(1, og_reference[msa_og_j[pos]]));
    }
    
    std::vector<int> pr_og_all_comp = get_all_comparisons(pr, og);
    std::vector<int> ll_og_all_comp = get_all_comparisons(pp, og);
    std::vector<int> pp_og_all_comp = get_all_comparisons(ll, og);
    
    std::vector<std::string> pr_og_agreements = get_agreements(pr, og, positions_idx);
    std::vector<std::string> pp_og_agreements = get_agreements(pp, og, positions_idx);
    std::vector<std::string> ll_og_agreements = get_agreements(ll, og, positions_idx);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    
    std::cout << pileup_file << "," << pr_og_agreements[0] << "," << pp_og_agreements[0] << "," << ll_og_agreements[0] << "\n" << "Execution time: " << duration/1000.0 << " sec" << std::endl;
    
    return 0;
}


//#include <iostream>
//#include <string>
//#include <fstream>
//#include <sstream>
//#include <vector>
//#include <map>
//#include <regex>
//#include <cmath>
//#include <numeric>
//
//std::vector<std::string> header = {"A", "C", "G", "T", "-"};
//
//std::vector<double> get_q(std::string ref, std::string match, double e) {
//    std::vector<double> qualities(5, e/3);
//    qualities[4] = 0.0;
//    if (match == "," || match == ".") {
//        match = ref;
//    }
//    if (std::find(header.begin(), header.end(), match) != header.end()) {
//        qualities[std::distance(header.begin(), std::find(header.begin(), header.end(), match))] = 1 - e;
//        return qualities;
//    }
//    if (match == "*") {
//        return std::vector<double>(5, 1);
//    }
//    return std::vector<double>(5, NAN);
//}
//
//std::vector<int> get_all_comparisons(std::vector<std::string> list1, std::vector<std::string> list2) {
//    std::vector<int> comparisons;
//    comparisons.reserve(list1.size());
//    for (size_t i = 0; i < list1.size(); ++i) {
//        if (list1[i] == "N" || list2[i] == "N") {
//            comparisons.push_back(-1);
//        } else if (list2[i] == "-") {
//            comparisons.push_back(-1);
//        } else {
//            comparisons.push_back(list1[i] == list2[i]);
//        }
//    }
//    return comparisons;
//}
//
//std::vector<std::string> get_agreements(std::vector<std::string> list1, std::vector<std::string> list2, std::vector<int> pos_list) {
//    std::vector<int> ag_all;
//    ag_all.reserve(list1.size());
//    for (size_t i = 0; i < list1.size(); ++i) {
//        if (list1[i] == "N" || list2[i] == "N") {
//            ag_all.push_back(-1);
//        } else if (list2[i] == "-") {
//            ag_all.push_back(-1);
//        } else {
//            ag_all.push_back(list1[i] == list2[i]);
//        }
//    }
//
//    std::vector<int> ag;
//    ag.reserve(ag_all.size());
//    for (int a : ag_all) {
//        if (a != -1) {
//            ag.push_back(a);
//        }
//    }
//
//    double mean_agreement = static_cast<double>(std::accumulate(ag.begin(), ag.end(), 0)) / ag.size();
//    int sum_agreement = std::accumulate(ag.begin(), ag.end(), 0);
//    int agreement_count = ag.size();
//    int disagreement_count = static_cast<int>(ag.size()) - sum_agreement;
//
//    return {std::to_string(mean_agreement), std::to_string(sum_agreement), std::to_string(agreement_count), std::to_string(disagreement_count)};
//}
//
//int main(int argc, char* argv[]) {
//    if (argc < 3) {
//        std::cerr << "Usage: ./program_name pileup_file pileup_results_path" << std::endl;
//        return 1;
//    }
//
//    std::string pileup_file = argv[1];
//    std::string pileup_results_path = argv[2];
//
//    std::string home_dir = "/Users/marniella/research/nielsen_lab/phylogenetic-base-calling/";
//    std::string msa_dir = "ou_pruned/";
//
//    std::ifstream f(home_dir + "justmn.fasta");
//    std::string reference;
//    int lineCount = 1;
//    while (std::getline(f, reference)) {
//        if (lineCount == 2) {
//            break;
//        }
//        lineCount++;
//    }
//    f.close();
//
//    std::ifstream f2(home_dir + "justou.fasta");
//    std::string og_reference;
//    int lineCount2 = 1;
//    while (std::getline(f2, og_reference)) {
//        if (lineCount2 == 2) {
//            break;
//        }
//        lineCount2++;
//    }
//    f2.close();
//
//
//    std::vector<int> msa_og_j;
//    for (size_t j = 0; j < reference.size(); ++j) {
//        if (reference[j] != '-') {
//            msa_og_j.push_back(j);
//        }
//    }
//
//    std::regex pattern1("-[0-9]+[ACGTNacgtn]+");
//    std::regex pattern2("\\+[0-9]+[ACGTNacgtn*#]+");
//    std::regex pattern3("\\$");
//    std::regex pattern4("\\^.{1}");
//
//    std::map<int, std::vector<double>> bases_prior_dict;
//    std::ifstream f3(home_dir + msa_dir + "OU_PP_indexed.txt");
//    std::string line;
//    std::getline(f3, line); // Skip the first line
//    while (std::getline(f3, line)) {
//        std::istringstream iss(line);
//        std::string msa_pos, prior_A, prior_C, prior_G, prior_T;
//        iss >> msa_pos >> prior_A >> prior_C >> prior_G >> prior_T;
//        bases_prior_dict[std::stoi(msa_pos)] = {std::stod(prior_A), std::stod(prior_C), std::stod(prior_G), std::stod(prior_T), 0.0};
//    }
//    f3.close();
//
//    std::map<int, std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>> pileup_dict;
//    std::ifstream f4(home_dir + pileup_results_path + pileup_file);
//    while (std::getline(f4, line)) {
//        std::istringstream iss(line);
//        std::string name, pos, ref_base, n_reads, match_str, phred_list;
//        iss >> name >> pos >> ref_base >> n_reads >> match_str >> phred_list;
//        int position = std::stoi(pos) - 1;
//        match_str = std::regex_replace(
//                std::regex_replace(
//                    std::regex_replace(
//                        std::regex_replace(match_str, pattern4, ""), pattern3, ""),pattern2, ""),pattern1, "");
//
//        std::vector<std::vector<double>> q_list;
//        for (size_t i = 0; i < match_str.size(); ++i) {
//            q_list.push_back(get_q(ref_base, std::string(1, match_str[i]), std::pow(10, -(static_cast<double>(phred_list[i]) - 33) / 10.0)));
//        }
//
//        std::vector<double> ll_list(q_list[0].size());
//
//        for (size_t i = 0; i < q_list[0].size(); ++i) {
//            double prod = 1.0;
//            for (size_t j = 0; j < q_list.size(); ++j) {
//                prod *= q_list[j][i];
//            }
//            ll_list[i] = std::isnan(prod) ? 0.0 : prod;
//        }
//
//        int msa_pos = msa_og_j[position];
//        std::vector<double> pr_list = bases_prior_dict[msa_pos];
//        std::vector<double> pp_list(ll_list.size(), 0.0);
//
//        if ((std::accumulate(ll_list.begin(), ll_list.end(), 0.0) == 5.0) && (og_reference[msa_pos] == '-')) {
//            pp_list = {0.0, 0.0, 0.0, 0.0, 1.0};
//        } else {
//            for (size_t i = 0; i < ll_list.size(); ++i) {
//                pp_list[i] = pr_list[i] * ll_list[i];
//            }
//        }
//
//        pileup_dict[position] = std::make_tuple(pr_list, ll_list, pp_list);
//    }
//    f4.close();
//
//    std::vector<int> positions_idx;
//    for (const auto& entry : pileup_dict) {
//        positions_idx.push_back(entry.first);
//    }
//    std::sort(positions_idx.begin(), positions_idx.end());
//
//    std::vector<std::string> pr, pp, ll, og;
//    for (int pos : positions_idx) {
//        pr.push_back(header[std::distance(std::get<0>(pileup_dict[pos]).begin(), std::max_element(std::get<0>(pileup_dict[pos]).begin(), std::get<0>(pileup_dict[pos]).end()))]);
//        ll.push_back(header[std::distance(std::get<1>(pileup_dict[pos]).begin(), std::max_element(std::get<1>(pileup_dict[pos]).begin(), std::get<1>(pileup_dict[pos]).end()))]);
//        pp.push_back(header[std::distance(std::get<2>(pileup_dict[pos]).begin(), std::max_element(std::get<2>(pileup_dict[pos]).begin(), std::get<2>(pileup_dict[pos]).end()))]);
//
//        og.push_back(std::string(1, og_reference[msa_og_j[pos]]));
//    }
//
//    std::vector<int> pr_og_all_comp = get_all_comparisons(pr, og);
//    std::vector<int> ll_og_all_comp = get_all_comparisons(pp, og);
//    std::vector<int> pp_og_all_comp = get_all_comparisons(ll, og);
//
//    std::vector<std::string> pr_og_agreements = get_agreements(pr, og, positions_idx);
//    std::vector<std::string> pp_og_agreements = get_agreements(pp, og, positions_idx);
//    std::vector<std::string> ll_og_agreements = get_agreements(ll, og, positions_idx);
//
//    std::cout << pileup_file << "," << pr_og_agreements[0] << "," << pp_og_agreements[0] << "," << ll_og_agreements[0] << std::endl;
//
//    return 0;
//}
