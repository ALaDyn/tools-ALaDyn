
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <boost/algorithm/string.hpp> 

#define SEPARATORS       ",= \t"
#define NPTOT_T_POS      4
#define NPTOT_N_POS      5
#define NPTOT_STR_POS    4
#define NPTOT_STR        "nptot"

std::vector<std::vector<std::string>> Parse_file(std::string, std::string);

typedef struct np_per_t {
  double t;
  size_t np;
} nptot;


int main(int argc, char*argv[]) {
  std::string filename;

  filename = "opic.txt";

  std::vector< std::vector<std::string> > parsed_file = Parse_file(filename, SEPARATORS);
  std::vector<np_per_t> npx;

  size_t nptot_str_pos = NPTOT_STR_POS - 1;
  size_t nptot_n_pos = NPTOT_N_POS - 1;
  size_t nptot_t_pos = NPTOT_T_POS - 1;

  for (size_t i = 1; i < parsed_file.size(); i++) {
    if (parsed_file[i].size() > nptot_str_pos && parsed_file[i].size() > nptot_n_pos &&parsed_file[i - 1].size() > nptot_t_pos) {
      if (parsed_file[i][nptot_str_pos] == NPTOT_STR) {
        //std::cout << "nptot:\"" << parsed_file[i][nptot_str_pos] << "\"\tnptot=" << parsed_file[i][nptot_n_pos] << "\tt=" << parsed_file[i - 1][nptot_t_pos] << std::endl;
        np_per_t npl;
        try {
          npl.np = stoull(parsed_file[i][nptot_n_pos]);
        }
        catch (std::exception e) {
          std::cout << e.what() << std::endl;
        }
        try {
          npl.t = stod(parsed_file[i - 1][nptot_t_pos]);
        }
        catch (std::exception e) {
          std::cout << e.what() << std::endl;
        }
        npx.push_back(npl);
      }
    }
  }

  filename = "nptot.dat";
  std::ofstream nptot(filename);
  double perc = 0.0;
  size_t ref_np = npx.front().np;
  nptot << "#time\tnp\t%" << std::endl;
  for (auto i : npx) {
    perc = double(i.np) / double(ref_np) * 100.0;
    nptot << i.t << "\t" << i.np << "\t" << perc << std::endl;
  }
  nptot.close();

  return 0;
}



std::vector< std::vector<std::string> > Parse_file(std::string file_name, std::string separators) {
  // Safe file opening
  std::ifstream file_to_parse(file_name, std::ios::in);
  if (!file_to_parse) {
    std::cout << "Cannot open " << file_name << ". Quitting..." << std::endl;
    exit(12);
  }
  // Internal variables
  std::string line;
  std::vector<std::string> tokens;
  std::vector< std::vector<std::string> > parsed;
  while (getline(file_to_parse, line)) {
    boost::algorithm::trim(line);  // remove leading/trailing spaces
    boost::algorithm::split(tokens, line, boost::algorithm::is_any_of(separators), boost::token_compress_on);
    //std::transform(tokens[0].begin(), tokens[0].end(), tokens[0].begin(), ::tolower);
    if (tokens.size()) {
      parsed.push_back(tokens);
    }
    line.clear(); tokens.clear();
  }
  file_to_parse.close();
  return parsed;
}


