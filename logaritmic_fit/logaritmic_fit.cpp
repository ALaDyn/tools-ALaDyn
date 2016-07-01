/******************************************************************************
Copyright 2015 Stefano Sinigardi
The program is distributed under the terms of the GNU General Public License
******************************************************************************/

/**************************************************************************
This file is part of "tools-ALaDyn".

tools-ALaDyn is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

tools-ALaDyn is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with tools-ALaDyn.  If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/



#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS



#include <iostream>
#include <cstdio>
#include <sstream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <limits>



bool AreSame(double a, double b) {
  return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
}


// this program tries a logaritmic fit for energy evolution of particle species.
// data is supposed to come from diags converted from aladyn output to columns through another tool in this collection
// you have to specify on the command line the X column and the Y column

int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    std::cerr << "Please write on command line the input file (x,y), in which columns we will find the data and the working mode!" << std::endl;
    std::cerr << "-x N means that the x data will be in the N column of the file" << std::endl;
    std::cerr << "-y M means that the y data will be in the M column of the file" << std::endl;
    std::cerr << "-scan to write on the output, on a single line and without the newline at the end, just the fitting parameters" << std::endl;
    std::cerr << "-func to write on the output the fitting functions" << std::endl;
    std::cerr << "-gnuplot to write on the output the gnuplot script useful to plot the input file including the fitting curve" << std::endl;
    std::cerr << "\nIn all cases, this program works best using output redirection" << std::endl;
    exit(1);
  }

  int Xres = 1280;
  int Yres = 720;
  char image_type[] = "png";
  int colonna_x = 1, colonna_y = 2;
  bool scan = false, func = false, gnuplot = false;
  int inputfile_position = 0;
  std::string riga;
  std::vector<std::string> righe;
  std::ifstream infile;
  std::stringstream ss;
  double sum_y = 0., sum_logx = 0., sum_ylogx = 0., sum_log2x = 0.;
  double x = 0., y = 0., logx = 0., logy = 0.;
  double min_val = std::numeric_limits<double>::epsilon();
  double max_val = std::numeric_limits<double>::infinity();


  for (int i = 1; i < argc; i++) {
    /*******************************************************************
    * We will iterate over argv[] to get the parameters stored inside. *
    * Note that we're starting on 1 because we don't need to know the  *
    * path of the program, which is stored in argv[0]                  *
    *******************************************************************/
    if      (std::string(argv[i]) == "-scan")     scan = true;
    else if (std::string(argv[i]) == "-gnuplot")  gnuplot = true;
    else if (std::string(argv[i]) == "-func")     func = true;
    else if (std::string(argv[i]) == "-x")        colonna_x = atoi(argv[++i]);
    else if (std::string(argv[i]) == "-y")        colonna_y = atoi(argv[++i]);
    else if (std::string(argv[i]) == "-min")      min_val = atof(argv[++i]);
    else if (std::string(argv[i]) == "-max")      max_val = atof(argv[++i]);
    else                                          inputfile_position = i;
  }

  if (inputfile_position < 1) {
    std::cerr << "Unable to find an input file on the command line" << std::endl;
    std::cerr << "   Scan [0:disabled, 1:enabled] --> " << scan << std::endl;
    std::cerr << "Gnuplot [0:disabled, 1:enabled] --> " << gnuplot << std::endl;
    std::cerr << "   Func [0:disabled, 1:enabled] --> " << func << std::endl;
    exit(2);
  }

  infile.open(argv[inputfile_position], std::ifstream::in);
  if (!infile.is_open()) {
    std::cerr << "Unable to open input file " << argv[inputfile_position] << "!" << std::endl;
    exit(3);
  }

  while (true) {
    std::getline(infile, riga);
    if (infile.eof()) break;
    std::size_t found;
    found = riga.find("#");
    if ((found == std::string::npos || found > 1) && !infile.eof())
      righe.push_back(riga);
  }


  int contacolonne = 0;
  char * pch;
  char * str = new char[righe[0].length() + 1];
  std::strcpy(str, righe[0].c_str());
  pch = std::strtok(str, " \t");
  if (pch != NULL) contacolonne++;
  while (pch != NULL) {
    pch = strtok(NULL, " \t");
    if (pch != NULL) contacolonne++;
  }
  if (colonna_x > contacolonne || colonna_x < 1) {
    std::cerr << "You have requested to analyze the x-column at a position (" << colonna_x << ") not available! " << contacolonne << " columns found." << std::endl;
    exit(4);
  }
  if (colonna_y > contacolonne || colonna_y < 1) {
    std::cerr << "You have requested to analyze the y-column at a position (" << colonna_y << ") not available! " << contacolonne << " columns found." << std::endl;
    exit(5);
  }

  colonna_x--;  // to convert from human-form to C-style (first column in a file is represented as column #0 in C)
  colonna_y--;

  infile.close();
  size_t n = righe.size();
  int reduced_n = 0;
  double * energy = new double[n];
  double * time = new double[n];
  double * values = new double[contacolonne];

  for (size_t it = 0; it < n; it++) {
    std::stringstream ss(righe.at(it));
    for (size_t line_it = 0; line_it < contacolonne; line_it++) ss >> values[line_it];
    time[it] = values[colonna_x];
    energy[it] = values[colonna_y];
  }

  for (size_t it = 0; it < righe.size(); it++) {
    if (time[it] < min_val || time[it] > max_val) continue;
    x = time[it];
    y = energy[it];
    x > 0.0 ? logx = log(x) : 0.0;
    y > 0.0 ? logy = log(y) : 0.0;

    sum_y += y;
    sum_logx += logx;
    sum_ylogx += y*logx;
    sum_log2x += (logx*logx);

    reduced_n++;
  }


  double fit_a1 = 0.0, fit_b1 = 0.0;    // see http://mathworld.wolfram.com/LeastSquaresFittingLogarithmic.html
  double denom_1 = reduced_n*sum_log2x - sum_logx*sum_logx;

  if (!AreSame(denom_1, 0.0) && reduced_n != 0) {
    fit_b1 = (reduced_n*sum_ylogx - sum_y*sum_logx) / denom_1;
    fit_a1 = (sum_y - fit_b1*sum_logx) / reduced_n;
  }
  else {
    std::cerr << denom_1 << " is too close to zero to work!" << std::endl;
    exit(4);
  }


  if (gnuplot) {
    FILE*  outfile;
    outfile = fopen("plot.plt", "w");

    fprintf(outfile, "#!/gnuplot\n");
    fprintf(outfile, "FILE_IN='%s'\n", argv[inputfile_position]);
    fprintf(outfile, "FILE_OUT='%s.%s'\n", argv[inputfile_position], image_type);
    fprintf(outfile, "set terminal %s truecolor enhanced size %i,%i\n", image_type, Xres, Yres);
    fprintf(outfile, "set output FILE_OUT\n");
    fprintf(outfile, "a = %g\n", fit_a1);
    fprintf(outfile, "b = %g\n", fit_b1);
    fprintf(outfile, "f(x) = a + b*log(x)\n");
    fprintf(outfile, "set xlabel 't' \n");
    fprintf(outfile, "set ylabel 'E (MeV)'\n");
    fprintf(outfile, "set logscale x\n");
    fprintf(outfile, "plot FILE_IN u %i:%i w points pt 7 lc rgb 'blue' ps 1.5 notitle,\\",colonna_x+1, colonna_y+1);
    fprintf(outfile, "\n");
    fprintf(outfile, "f(x) w lines lt 1 lc rgb 'red' lw 2 notitle");
    fprintf(outfile, "\n");
    fclose(outfile);
  }

  if (func) std::cout << "Fit function: y=" << fit_a1 << " + " << fit_b1 << "*ln(x)" << std::endl;
  if (scan) std::cout << " \t " << fit_a1 << " \t " << fit_b1 << " \t ";

  return 0;
}

