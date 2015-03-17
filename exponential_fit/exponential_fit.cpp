/******************************************************************************
Copyright 2014, 2015 Stefano Sinigardi
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



#define LUNGHEZZA_MAX_RIGA 1024
#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS



#include <iostream>
#include <cstdio>
#include <sstream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>


using namespace std;



int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    cerr << "Please write input file on command line and also working mode!" << endl;
    cerr << "-scan to write on the output, on a single line and without the newline at the end, just the mean energy and the total number of particles (fitting parameters)" << endl;
    cerr << "-func to write on the output the fitting functions" << endl;
    cerr << "-gnuplot to write on the output the gnuplot script useful to plot the input file including the fitting curves" << endl;
    cerr << "\nIn all cases, this program works best using output redirection" << endl;
    exit(1);
  }

  bool scan = false, func = false, gnuplot = false;
  int inputfile_position = 0;

  for (int i = 1; i < argc; i++)
    /************************************************************************
    We will iterate over argv[] to get the parameters stored inside.
    Note that we're starting on 1 because we don't need to know the
    path of the program, which is stored in argv[0]
    ************************************************************************/
  {
    if (std::string(argv[i]) == "-scan")
    {
      scan = true;
    }
    else if (std::string(argv[i]) == "-gnuplot")
    {
      gnuplot = true;
    }
    else if (std::string(argv[i]) == "-func")
    {
      func = true;
    }
    else
    {
      inputfile_position = i;
    }


  }

  if (inputfile_position < 1)
  {
    cerr << "Unable to find an input file on the command line" << endl;
    cerr << "Scan[0:disabled, 1:enabled] --> " << scan << endl;
    cerr << "Gnuplot[0:disabled, 1:enabled] --> " << gnuplot << endl;
    cerr << "Func[0:disabled, 1:enabled] --> " << func << endl;
    exit(2);
  }

  string riga;
  vector<string> righe;

  ifstream infile;
  infile.open(argv[inputfile_position], ifstream::in);
  if (!infile.is_open())
  {
    cerr << "Unable to open input file " << argv[inputfile_position] << "!" << endl;
    exit(3);
  }


  char * riga_letta;
  riga_letta = new char[LUNGHEZZA_MAX_RIGA];

  while (true)
  {
    infile.getline(riga_letta, LUNGHEZZA_MAX_RIGA);
    if (infile.eof()) break;
    riga = riga_letta;
    std::size_t found;
    found = riga.find("#");
    if ((found == std::string::npos || found > 1) && !infile.eof())
      righe.push_back(riga);
  }

  infile.close();

  double * energies = new double[righe.size()];
  double * particles = new double[righe.size()];
  double * particles_selected = new double[righe.size()];
  stringstream ss;

  for (unsigned int it = 0; it < righe.size(); it++)
  {
    stringstream ss(righe.at(it));
    ss >> energies[it] >> particles[it] >> particles_selected[it];
  }


  double sum_x = 0., sum_x2 = 0.;
  double sum_y = 0., sum_y2 = 0., sum_logy = 0., sum_x2y = 0., sum_ylogy = 0., sum_xy = 0., sum_xylogy = 0.;
  double sum_z = 0., sum_z2 = 0., sum_logz = 0., sum_x2z = 0., sum_zlogz = 0., sum_xz = 0., sum_xzlogz = 0.;
  double x = 0., y = 0., z = 0., logx = 0., logy = 0., logz = 0.;

  unsigned int inizio = (int)(righe.size() / 5.);
  unsigned int fine = (int)(4. * righe.size() / 5.);
  //  for (unsigned int it = 0; it < righe.size(); it++)
  for (unsigned int it = inizio; it < fine; it++)
  {
    x = energies[it];
    y = particles[it];
    z = particles_selected[it];
    x > 0.0 ? logx = log(x) : 0.0;
    y > 0.0 ? logy = log(y) : 0.0;
    z > 0.0 ? logz = log(z) : 0.0;

    sum_x += x;
    sum_x2 += x*x;
    sum_y += y;
    sum_y2 += y*y;
    sum_logy += logx;
    sum_x2y += x*x*y;
    sum_ylogy += y*logy;
    sum_xy += x*y;
    sum_xylogy += x*y*logy;
    sum_z += z;
    sum_z2 += z*z;
    sum_logz += logz;
    sum_x2z += x*x*z;
    sum_zlogz += z*logz;
    sum_xz += x*z;
    sum_xzlogz += x*z*logz;
  }


  double fit_a1, fit_b1; // see http://mathworld.wolfram.com/LeastSquaresFittingExponential.html
  double fit_a2, fit_b2;

  double den;

  den = (sum_y*sum_x2y - sum_xy*sum_xy);
  if (den > 0 || den < 0) {
    fit_a1 = (sum_x2y*sum_ylogy - sum_xy*sum_xylogy) / den;
    fit_b1 = (sum_y*sum_xylogy - sum_xy*sum_ylogy) / den;
  }
  else {
    fit_a1 = fit_b1 = 0.0;
  }

  den = (sum_z*sum_x2z - sum_xz*sum_xz);
  if (den > 0 || den < 0) {
    fit_a2 = (sum_x2z*sum_zlogz - sum_xz*sum_xzlogz) / den;
    fit_b2 = (sum_z*sum_xzlogz - sum_xz*sum_zlogz) / den;
  }
  else {
    fit_a2 = fit_b2 = 0.0;
  }


  double aveE1 = -1. / fit_b1;
  double aveE2 = -1. / fit_b2;
  int N0_1 = (int)(exp(fit_a1) * aveE1);
  int N0_2 = (int)(exp(fit_a2) * aveE2);
  int weight = 2; // fix, read it from the first line of the infile
  int subsample_factor = 48230; // fix, read it from the first line of the infile




  if (func)
  {
    cout << "Fit for full spectrum: y=" << exp(fit_a1) << "e^(" << fit_b1 << "x)" << endl;
    cout << "Fit for selected spectrum: y=" << exp(fit_a2) << "e^(" << fit_b2 << "x)" << endl;
  }


  if (gnuplot)
  {
    FILE*  outfile;
    outfile = fopen("plot.plt", "w");

    int Xres = 1280; // fix, make it possible to define on command line
    int Yres = 720; // fix, make it possible to define on command line
    //  int Emin = 0; // fix, read it from the infile
    //  int Emax = 60; // fix, read it from the infile
    char image_type[] = "png";

    fprintf(outfile, "#!/gnuplot\n");
    fprintf(outfile, "FILE_IN='%s'\n", argv[1]);
    fprintf(outfile, "FILE_OUT='%s.%s'\n", argv[1], image_type);
    fprintf(outfile, "set terminal %s truecolor enhanced size %i,%i\n", image_type, Xres, Yres);
    fprintf(outfile, "set output FILE_OUT\n");
    fprintf(outfile, "AVERAGE_E1 = %3.2f\n", aveE1);
    fprintf(outfile, "AVERAGE_E2 = %3.2f\n", aveE2);
    fprintf(outfile, "WEIGHT = %i\n", weight);
    fprintf(outfile, "SUBSAMPLE = %i\n", subsample_factor);
    fprintf(outfile, "N0_1 = %i*WEIGHT*SUBSAMPLE\n", N0_1);
    fprintf(outfile, "N0_2 = %i*WEIGHT*SUBSAMPLE\n", N0_2);
    fprintf(outfile, "f(x) = (N0_1 / AVERAGE_E1)*exp(-x / AVERAGE_E1)\n");
    fprintf(outfile, "g(x) = (N0_2 / AVERAGE_E2)*exp(-x / AVERAGE_E2)\n");
    fprintf(outfile, "set xlabel 'E (MeV)' \n");
    fprintf(outfile, "set ylabel 'dN/dE (MeV^{-1})'\n");
    fprintf(outfile, "set format y '10^{%%L}'\n");
    //  fprintf(outfile, "set xrange[%i:%i]\n", Emin, Emax);
    fprintf(outfile, "set logscale y\n");
    fprintf(outfile, "plot FILE_IN u 1:($2*%i*%i) w histeps lt 1 lc rgb 'blue' lw 3 t 'full spectrum',\\", weight, subsample_factor);
    fprintf(outfile, "\n");
    fprintf(outfile, "FILE_IN u 1:($3*%i*%i) w histeps lt 1 lc rgb 'red' lw 3 t 'selected spectrum',\\", weight, subsample_factor);
    fprintf(outfile, "\n");
    fprintf(outfile, "f(x) w lines lt 1 lc rgb 'purple' lw 3 t 'exponential fit E_01 = %1.1f MeV',\\", aveE1);
    fprintf(outfile, "\n");
    fprintf(outfile, "g(x) w lines lt 1 lc rgb 'dark-green' lw 3 t 'exponential fit E_02 = %1.1f MeV'", aveE2);

    fprintf(outfile, "\n");
    fclose(outfile);
  }


  if (scan)
  {
    printf(" \t %3.2f \t %i \t %3.2f \t %i \t ", aveE1, N0_1, aveE2, N0_2);
  }



  return 0;

}


