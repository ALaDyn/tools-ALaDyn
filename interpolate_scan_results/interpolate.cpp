#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS

#include <iostream>
#include <cstdio>
#include <sstream>
#include <vector>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdint>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>



int column_x, column_y, column_E;


bool sortAscendingByTwoColumns(double * riga1, double * riga2) {
  return ((riga1[column_x] < riga2[column_x]) || ((riga1[column_x] == riga2[column_x]) && (riga1[column_y] < riga2[column_y])));
}


int main(int argc, const char* argv[]) {
  size_t ncolumns, nrows;
  size_t interpolation_x, interpolation_y;
  double x1, y1, x2, y2, E11, E12, E21, E22, x, y, E, dx, dy;

  std::string filename;
  std::ifstream infile;
  std::ofstream outfile;
  std::string riga;
  std::vector<std::string> tokens;
  std::vector<double> dtokens;
  std::vector< std::vector<double> > matrix;
  std::vector< std::vector< std::vector<double> > > bigblock;

  bool found;
  std::vector<std::string> unique_x_values, unique_y_values;
  std::vector<double> dunique_x_values, dunique_y_values;
  double ** pmatrix;
  double ** values;

  for (int i = 1; i < argc; i++) {
    /************************************************************************
    We will iterate over argv[] to get the parameters stored inside.
    Note that we're starting on 1 because we don't need to know the
    path of the program, which is stored in argv[0]
    ************************************************************************/
    if (std::string(argv[i]) == "-cx") {
      column_x = boost::lexical_cast<int>(argv[++i]) - 1; // mind the -1 to port the index to C-style
    }
    else if (std::string(argv[i]) == "-cy") {
      column_y = boost::lexical_cast<int>(argv[++i]) - 1; // mind the -1 to port the index to C-style
    }
    else if (std::string(argv[i]) == "-ce") {
      column_E = boost::lexical_cast<int>(argv[++i]) - 1; // mind the -1 to port the index to C-style
    }
    else if (std::string(argv[i]) == "-nx") {
      interpolation_x = boost::lexical_cast<size_t>(argv[++i]);
    }
    else if (std::string(argv[i]) == "-ny") {
      interpolation_y = boost::lexical_cast<size_t>(argv[++i]);
    }
    else if (std::string(argv[i]) == "-file") {
      filename = argv[++i];
    }
  }

  infile.open(filename, std::ifstream::in);
  outfile.open(filename+"_int", std::ofstream::out);

  if (!infile.is_open()) {
    std::cerr << "Unable to open input file " << filename << "!" << std::endl;
    exit(4);
  }

  if (!outfile.is_open()) {
    std::cerr << "Unable to open output file " << filename+"_int" << "!" << std::endl;
    exit(5);
  }


  while (true) {
    riga.clear(), tokens.clear(), dtokens.clear();
    std::getline(infile, riga);
    if (infile.eof()) break;
    boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_off);
    if (tokens[0][0] == '#') continue;
    if (!matrix.size()) ncolumns = tokens.size();
    if (ncolumns != tokens.size()) {
      std::cout << "Warning, row has an unexpected length" << std::endl;
      exit(1);
    }

    found = false;
    for (auto i : unique_x_values) {
      if (i == tokens[column_x]) found = true;
    }
    if (!found) unique_x_values.push_back(tokens[column_x]);

    found = false;
    for (auto i : unique_y_values) {
      if (i == tokens[column_y]) found = true;
    }
    if (!found) unique_y_values.push_back(tokens[column_y]);

    for (auto i : tokens) dtokens.push_back(boost::lexical_cast<double>(i));
    matrix.push_back(dtokens);
  }
  for (auto i : unique_x_values) dunique_x_values.push_back(boost::lexical_cast<double>(unique_x_values));
  for (auto i : unique_y_values) dunique_y_values.push_back(boost::lexical_cast<double>(unique_y_values));


  nrows = matrix.size();
  pmatrix = new double*[nrows];
  for (size_t i = 0; i < nrows; i++) pmatrix[i] = new double[ncolumns];

  for (size_t i = 0; i < nrows; i++) {
    for (size_t j = 0; j < ncolumns; j++) pmatrix[i][j] = matrix[i][j];
  }
  riga.clear(), tokens.clear(), dtokens.clear(), matrix.clear();

  std::sort(pmatrix, pmatrix + nrows, &sortAscendingByTwoColumns);

  values = new double*[nrows];
  for (size_t i = 0; i < nrows; i++) values[i] = new double[3];
  for (size_t i = 0; i < nrows; i++) {
    values[i][0] = pmatrix[i][column_x];
    values[i][1] = pmatrix[i][column_y];
    values[i][2] = pmatrix[i][column_E];
  }

  // manca il riempibuchi


  for (size_t i = 1; i < dunique_x_values.size(); i++) {
    for (size_t j = 1; j < dunique_x_values.size(); j++) {
      x1 = dunique_x_values[i - 1];
      y1 = dunique_y_values[j - 1];
      x2 = dunique_x_values[i];
      y2 = dunique_y_values[j];
      E11 = values[(i - 1)*dunique_x_values.size() + j - 1][2];
      E12 = values[(i - 1)*dunique_x_values.size() + j][2];
      E21 = values[i*dunique_x_values.size() + j - 1][2];
      E22 = values[i*dunique_x_values.size() + j][2];

      for (size_t k = 0; k < interpolation_x; k++) {
        for (size_t l = 0; l < interpolation_y; l++) {

        }
      }
    }
  }

  infile.close();
  outfile.close();
  return 0;
}