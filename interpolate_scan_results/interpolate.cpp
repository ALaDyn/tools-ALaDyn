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
#include <cmath>
#include <iomanip>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#define EPSILON    0.0001
#define MAX(x,y)   (x > y ? x : y)

int column_x = 0, column_y = 0, column_E = 0;

bool checkEqual(double a, double b) {
  return (fabs(a - b) < EPSILON);
}

bool sortAscendingByTwoColumns(std::vector<double>& riga1, std::vector<double>& riga2) {
  return ((riga1[column_x] < riga2[column_x]) || ((checkEqual(riga1[column_x], riga2[column_x])) && (riga1[column_y] < riga2[column_y])));
}


int main(int argc, const char* argv[]) {
  size_t ncolumns, precision = 4;
  size_t interpolation_x, interpolation_y;
  double x1, y1, x2, y2, E11, E12, E21, E22, x, y, E, dx, dy, k0, E_magn = 1.0;
  int column_max;
  std::string filename_in, filename_out, filename_gnuplot_plt, filename_gnuplot_png, column_E_string;
  std::ifstream infile;
  std::ofstream outfile;
  std::string riga, xlabel = "", ylabel = "", cblabel = "", title = "";
  std::vector<std::string> tokens;
  std::vector<double> dtokens;
  std::vector< std::vector<double> > matrix;
  std::vector< std::vector<double> > results;

  bool found;
  bool gnuplot = false;
  std::vector<std::string> unique_x_values, unique_y_values;
  std::vector<double> dunique_x_values, dunique_y_values;

  if (argc < 7) {
    std::cerr << "Usage: " << argv[0] << " -cx N -cy M -ce L -nx A -ny B -file filename.txt\nOptional: [-gnuplot] [-title] [-xlabel] [-ylabel] [-cblabel] [-cb_magn] [-precision]" << std::endl;
    std::cerr << "Press a key to exit..." << std::endl;
    std::cin.get();
    exit(1);
  }

  for (int i = 1; i < argc; i++) {
    /************************************************************************
    We will iterate over argv[] to get the parameters stored inside.
    Note that we're starting on 1 because we don't need to know the
    path of the program, which is stored in argv[0]
    ************************************************************************/
    if (std::string(argv[i]) == "-cx") {
      column_x = boost::lexical_cast<int>(std::string(argv[++i])) - 1; // mind the -1 to port the index to C-style
    }
    else if (std::string(argv[i]) == "-cy") {
      column_y = boost::lexical_cast<int>(std::string(argv[++i])) - 1; // mind the -1 to port the index to C-style
    }
    else if (std::string(argv[i]) == "-ce") {
      column_E_string = std::string(argv[++i]);
      column_E = boost::lexical_cast<int>(column_E_string) - 1; // mind the -1 to port the index to C-style
    }
    else if (std::string(argv[i]) == "-nx") {
      interpolation_x = boost::lexical_cast<size_t>(std::string(argv[++i]));
    }
    else if (std::string(argv[i]) == "-ny") {
      interpolation_y = boost::lexical_cast<size_t>(std::string(argv[++i]));
    }
    else if (std::string(argv[i]) == "-file") {
      filename_in = argv[++i];
    }
    else if (std::string(argv[i]) == "-precision") {
      precision = boost::lexical_cast<size_t>(argv[++i]);
    }
    else if (std::string(argv[i]) == "-gnuplot") {
      gnuplot = true;
    }
    else if (std::string(argv[i]) == "-title") {
      gnuplot = true;
      title = std::string(argv[++i]);
    }
    else if (std::string(argv[i]) == "-xlabel") {
      gnuplot = true;
      xlabel = std::string(argv[++i]);
    }
    else if (std::string(argv[i]) == "-ylabel") {
      gnuplot = true;
      ylabel = std::string(argv[++i]);
    }
    else if (std::string(argv[i]) == "-cblabel") {
      gnuplot = true;
      cblabel = std::string(argv[++i]);
    }
    else if (std::string(argv[i]) == "-cb_magn") {
      E_magn = boost::lexical_cast<double>(argv[++i]);
    }
  }

  filename_out = "int" + column_E_string + '_' + filename_in;
  filename_gnuplot_plt = "plot.plt";
  filename_gnuplot_png = filename_out + ".png";
  column_max = MAX(column_x, column_y);
  column_max = MAX(column_max, column_E);

#ifdef ENABLE_DEBUG
  int linecounter = 0;
  std::cout << "column_x = " << column_x << std::endl;
  std::cout << "column_y = " << column_y << std::endl;
  std::cout << "column_E = " << column_E << ", column_E_string = " << column_E_string << std::endl;
  std::cout << "interpolation_x = " << interpolation_x << std::endl;
  std::cout << "interpolation_y = " << interpolation_y << std::endl;
  std::cout << "filename_in = " << filename_in << std::endl;
  std::cout << "filename_out = " << filename_out << std::endl;
  std::cout << "filename_gnuplot_plt = " << filename_gnuplot_plt << std::endl;
  std::cout << "filename_gnuplot_png = " << filename_gnuplot_png << std::endl;
  std::cout << "gnuplot = " << gnuplot << std::endl;
  std::cout << "title = " << title << std::endl;
  std::cout << "xlabel = " << xlabel << std::endl;
  std::cout << "ylabel = " << ylabel << std::endl;
  std::cout << "cblabel = " << cblabel << std::endl;
  std::cout << "E_magn = " << E_magn << std::endl;
  std::cout << "Press a key to continue... " << std::endl;
  std::cin.get();
#endif


  infile.open(filename_in, std::ifstream::in);
  outfile.open(filename_out, std::ofstream::out);

  if (!infile.is_open()) {
    std::cerr << "Unable to open input file " << filename_in << "!" << std::endl;
    exit(2);
  }

  if (!outfile.is_open()) {
    std::cerr << "Unable to open output file " << filename_out << "!" << std::endl;
    exit(3);
  }


  while (true) {
    riga.clear(), tokens.clear(), dtokens.clear();
    std::getline(infile, riga);
    if (infile.eof()) break;
    if (riga[0] == '#') continue;
    boost::algorithm::trim(riga);
    //boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_off);
    boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
    if (!matrix.size()) ncolumns = tokens.size();
    if (ncolumns != tokens.size()) {
      std::cout << "Warning, row has an unexpected length" << std::endl;
      exit(4);
    }
    if (ncolumns < column_max) {
      std::cout << "Warning, not enough columns in your file" << std::endl;
      exit(5);
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
  for (auto i : unique_x_values) dunique_x_values.push_back(boost::lexical_cast<double>(i));
  for (auto i : unique_y_values) dunique_y_values.push_back(boost::lexical_cast<double>(i));
  infile.close();

  riga.clear(), tokens.clear();
  dtokens.clear();
  dtokens.resize(3, 0);

  std::sort(&dunique_x_values.front(), &dunique_x_values.front() + dunique_x_values.size());
  std::sort(&dunique_y_values.front(), &dunique_y_values.front() + dunique_y_values.size());
  std::sort(&matrix.front(), &matrix.front() + matrix.size(), &sortAscendingByTwoColumns);

#ifdef ENABLE_DEBUG
  for (auto i : matrix) {
    for (auto j : i) std::cout << j << ' ';
    std::cout << std::endl;
  }
  for (auto i : dunique_x_values) std::cout << i << ' ';
  std::cout << std::endl;
  for (auto i : dunique_y_values) std::cout << i << ' ';
  std::cout << std::endl;
#endif


  // manca il riempibuchi

  for (size_t i = 1; i < dunique_x_values.size(); i++) {
    for (size_t j = 1; j < dunique_y_values.size(); j++) {
      x1 = dunique_x_values[i - 1];
      y1 = dunique_y_values[j - 1];
      x2 = dunique_x_values[i];
      y2 = dunique_y_values[j];
      E11 = matrix[(i - 1)*dunique_y_values.size() + j - 1][column_E];
      E12 = matrix[(i - 1)*dunique_y_values.size() + j][column_E];
      E21 = matrix[i*dunique_y_values.size() + j - 1][column_E];
      E22 = matrix[i*dunique_y_values.size() + j][column_E];
      k0 = 1.0 / ((x2 - x1)*(y2 - y1));
      dx = (x2 - x1) / interpolation_x;
      dy = (y2 - y1) / interpolation_y;
      for (size_t k = 0; k < interpolation_x; k++) {
        for (size_t l = 0; l < interpolation_y; l++) {
          x = x1 + k*dx;
          y = y1 + l*dy;
          E = k0*(E11*(x2 - x)*(y2 - y) + E21*(x - x1)*(y2 - y) + E12*(x2 - x)*(y - y1) + E22*(x - x1)*(y - y1));
          dtokens[0] = x;
          dtokens[1] = y;
          dtokens[2] = E*E_magn;
          results.push_back(dtokens);
        }
      }
    }
  }

  std::sort(&results.front(), &results.front() + results.size(), &sortAscendingByTwoColumns);
  for (auto i : results) {
    for (auto j : i) outfile << std::fixed << std::setprecision(precision) << j << ' ';
    outfile << std::endl;
  }
  outfile.close();

  if (gnuplot)
  {
    FILE*  outfile_gnuplot;
    outfile_gnuplot = fopen(filename_gnuplot_plt.c_str(), "w");
    int Xres = 1280;
    int Yres = 720;
    int fontsize = 15;
    char fontname[] = "";
    char image_type[] = "pngcairo";

    fprintf(outfile_gnuplot, "#!/gnuplot\n");
    fprintf(outfile_gnuplot, "set terminal %s size %i,%i font \"%s,%i\"\n", image_type, Xres, Yres, fontname, fontsize);
    fprintf(outfile_gnuplot, "set title '%s'\n", title.c_str());
    fprintf(outfile_gnuplot, "set xlabel '%s'\n", xlabel.c_str());
    fprintf(outfile_gnuplot, "set ylabel '%s'\n", ylabel.c_str());
    fprintf(outfile_gnuplot, "set cblabel '%s'\n", cblabel.c_str());
    fprintf(outfile_gnuplot, "set size ratio - 1\n");
    fprintf(outfile_gnuplot, "set xrange[%f:%f]\n", dunique_x_values.front(), dunique_x_values.back());
    fprintf(outfile_gnuplot, "set yrange[%f:%f]\n", dunique_y_values.front(), dunique_y_values.back());
    fprintf(outfile_gnuplot, "set palette rgbformulae 22, 13, 10\n");
    fprintf(outfile_gnuplot, "FILE_IN = \"%s\"\n", filename_out.c_str());
    fprintf(outfile_gnuplot, "FILE_OUT = \"%s\"\n", filename_gnuplot_png.c_str());
    fprintf(outfile_gnuplot, "set output FILE_OUT\n");
    fprintf(outfile_gnuplot, "plot FILE_IN u 1:2:3 w image notitle\n");

    fclose(outfile_gnuplot);
  }

  return 0;
}



