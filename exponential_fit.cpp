/******************************************************************************
Copyright 2014 Stefano Sinigardi
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



int main(int args, char* argv[])
{
  if (args < 2)
  {
    cerr << "Specificare su linea comando il file di input!" << endl;
    exit(1);
  }

  string riga;
  vector<string> righe;

  ifstream infile;
  infile.open(argv[1], ifstream::in);

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
  double x, y, z;

  unsigned int inizio = (int)(righe.size() / 5.);
  unsigned int fine = (int)(4. * righe.size() / 5.);
  //  for (unsigned int it = 0; it < righe.size(); it++)
  for (unsigned int it = inizio; it < fine; it++)
  {
    x = energies[it];
    y = particles[it];
    z = particles_selected[it];

    sum_x += x;
    sum_x2 += x*x;
    sum_y += y;
    sum_y2 += y*y;
    sum_logy += log(x);
    sum_x2y += x*x*y;
    sum_ylogy += y*log(y);
    sum_xy += x*y;
    sum_xylogy += x*y*log(y);
    sum_z += z;
    sum_z2 += z*z;
    sum_logz += log(z);
    sum_x2z += x*x*z;
    sum_zlogz += z*log(z);
    sum_xz += x*z;
    sum_xzlogz += x*z*log(z);
  }


  double fit_a1, fit_b1; // see http://mathworld.wolfram.com/LeastSquaresFittingExponential.html
  double fit_a2, fit_b2;

  fit_a1 = (sum_x2y*sum_ylogy - sum_xy*sum_xylogy) / (sum_y*sum_x2y - sum_xy*sum_xy);
  fit_b1 = (sum_y*sum_xylogy - sum_xy*sum_ylogy) / (sum_y*sum_x2y - sum_xy*sum_xy);

  fit_a2 = (sum_x2z*sum_zlogz - sum_xz*sum_xzlogz) / (sum_z*sum_x2z - sum_xz*sum_xz);
  fit_b2 = (sum_z*sum_xzlogz - sum_xz*sum_zlogz) / (sum_z*sum_x2z - sum_xz*sum_xz);

  cout << "Fit for full spectrum: y=" << exp(fit_a1) << "e^(" << fit_b1 << "x)" << endl;
  cout << "Fit for selected spectrum: y=" << exp(fit_a2) << "e^(" << fit_b2 << "x)" << endl;



  FILE*  outfile;
  outfile = fopen("plot.plt", "w");

  double aveE1 = -1. / fit_b1;
  double aveE2 = -1. / fit_b2;
  int N0_1 = (int)(exp(fit_a1) * aveE1);
  int N0_2 = (int)(exp(fit_a2) * aveE2);
  int weight = 2;
  int subsample_factor = 48230;
  int Xres = 1280;
  int Yres = 720;
  //  int Emin = 0;
  //  int Emax = 60;
  char image_type[] = "png";

  fprintf(outfile, "#!/gnuplot\n");
  fprintf(outfile, "FILE_IN='%s'\n", argv[1]);
  fprintf(outfile, "FILE_OUT='%s.%s'\n", argv[1], image_type);
  fprintf(outfile, "set terminal %s truecolor enhanced size %i,%i\n", image_type, Xres, Yres);
  fprintf(outfile, "set output FILE_OUT\n");
  fprintf(outfile, "AVERAGE_E1 = %2.1f\n", aveE1);
  fprintf(outfile, "AVERAGE_E2 = %2.1f\n", aveE2);
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

  return 0;

}


