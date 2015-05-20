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



#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#define MEV_TO_JOULE 1.602176565E-13

#define SKIP_INIZIALE  (1./5.)    // skip the first one fifth of the data
#define SKIP_FINALE    (4./5.)    // skip the last fifth of the data, since they do not fit very well into an exponential


#include <iostream>
#include <cstdio>
#include <sstream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <limits>



bool AreSame(double a, double b) {
  return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
}




int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    std::cerr << "Please write input file on command line and also working mode!" << std::endl;
    std::cerr << "-scan to write on the output, on a single line and without the newline at the end, just the mean energy and the total number of particles (fitting parameters)" << std::endl;
    std::cerr << "-func to write on the output the fitting functions" << std::endl;
    std::cerr << "-gnuplot to write on the output the gnuplot script useful to plot the input file including the fitting curves" << std::endl;
    std::cerr << "-piccante to read spectra made by piccante [in this case a -mass with the particle mass in MeV following is MANDATORY]" << std::endl;
    std::cerr << "\nIn all cases, this program works best using output redirection" << std::endl;
    exit(1);
  }

  bool scan = false, func = false, gnuplot = false, piccante = false;
  int inputfile_position = 0;
  double mass = -1.0;
  double nmacro_to_nphys = 1.0;

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
    else if (std::string(argv[i]) == "-piccante")
    {
      piccante = true;
    }
    else if (std::string(argv[i]) == "-mass")
    {
      mass = atof(argv[i + 1]);
      i++;
    }
    else if (std::string(argv[i]) == "-nm")
    {
      nmacro_to_nphys = atof(argv[i + 1]);
      i++;
    }
    else
    {
      inputfile_position = i;
    }
  }

  if (inputfile_position < 1)
  {
    std::cerr << "Unable to find an input file on the command line" << std::endl;
    std::cerr << "Scan[0:disabled, 1:enabled] --> " << scan << std::endl;
    std::cerr << "Gnuplot[0:disabled, 1:enabled] --> " << gnuplot << std::endl;
    std::cerr << "Func[0:disabled, 1:enabled] --> " << func << std::endl;
    exit(2);
  }

  if (piccante && mass <= 0)
  {
    std::cerr << "With piccante mode, you must give a particle mass in MeV/c^2 in order to calculate the total energy." << std::endl;
    std::cerr << "Please use -mass 1.0 and then disregard the total energy if not required" << std::endl;
    std::cerr << "Useful masses: proton=938.272013, electron=0.51099891" << std::endl;
    exit(3);
  }

  if (!piccante) mass = 1.0;

  std::string riga;
  std::vector<std::string> righe;

  std::ifstream infile;
  infile.open(argv[inputfile_position], std::ifstream::in);
  if (!infile.is_open())
  {
    std::cerr << "Unable to open input file " << argv[inputfile_position] << "!" << std::endl;
    exit(4);
  }


  while (true)
  {
    std::getline(infile, riga);
    if (infile.eof()) break;
    std::size_t found;
    found = riga.find("#");
    if ((found == std::string::npos || found > 1) && !infile.eof())
      righe.push_back(riga);
  }

  infile.close();

  double * energies = new double[righe.size()];
  double min_energy, max_energy, tot_energy = 0.0, tot_energy_front = 0.0, tot_energy_rear = 0.0;
  double * particles = new double[righe.size()];
  double * particles_front = new double[righe.size()];
  double * particles_rear = new double[righe.size()];
  std::stringstream ss;

  if (piccante)
  {
    for (unsigned int it = 0; it < righe.size(); it++)
    {
      std::stringstream ss(righe.at(it));
      ss >> min_energy >> max_energy >> particles[it];
      energies[it] = 0.5*(min_energy + max_energy);
      particles_front[it] = particles_rear[it] = 0.0;
      tot_energy += particles[it] * energies[it];
    }
  }
  else
  {
    for (unsigned int it = 0; it < righe.size(); it++)
    {
      std::stringstream ss(righe.at(it));
      ss >> energies[it] >> particles[it] >> particles_front[it] >> particles_rear[it];
      tot_energy += particles[it] * energies[it];
      tot_energy_front += particles_front[it] * energies[it];
      tot_energy_rear += particles_rear[it] * energies[it];
    }
  }

  max_energy = energies[righe.size() - 1] * mass;
  tot_energy *= mass * nmacro_to_nphys * MEV_TO_JOULE;
  tot_energy_front *= mass * nmacro_to_nphys * MEV_TO_JOULE;
  tot_energy_rear *= mass * nmacro_to_nphys * MEV_TO_JOULE;

  double sum_x = 0., sum_x2 = 0.;
  double sum_y = 0., sum_y2 = 0., sum_logy = 0., sum_x2y = 0., sum_ylogy = 0., sum_xy = 0., sum_xylogy = 0.;

  double x = 0., y = 0., f = 0., r = 0., logx = 0., logy = 0., logf = 0., logr = 0.;

  double sum_f = 0., sum_f2 = 0., sum_logf = 0., sum_x2f = 0., sum_flogf = 0., sum_xf = 0., sum_xflogf = 0.;
  double sum_r = 0., sum_r2 = 0., sum_logr = 0., sum_x2r = 0., sum_rlogr = 0., sum_xr = 0., sum_xrlogr = 0.;

  unsigned int inizio = (int)(SKIP_INIZIALE * righe.size());
  unsigned int fine = (int)(SKIP_FINALE * righe.size());
  for (unsigned int it = inizio; it < fine; it++)
  {
    x = energies[it];
    y = particles[it];
    f = particles_front[it];
    r = particles_rear[it];
    x > 0.0 ? logx = log(x) : 0.0;
    y > 0.0 ? logy = log(y) : 0.0;
    f > 0.0 ? logf = log(f) : 0.0;
    r > 0.0 ? logr = log(r) : 0.0;

    sum_x += x;
    sum_x2 += x*x;
    sum_y += y;
    sum_y2 += y*y;
    sum_logy += logx;
    sum_x2y += x*x*y;
    sum_ylogy += y*logy;
    sum_xy += x*y;
    sum_xylogy += x*y*logy;

    sum_f += f;
    sum_f2 += f*f;
    sum_logf += logf;
    sum_x2f += x*x*f;
    sum_flogf += f*logf;
    sum_xf += x*f;
    sum_xflogf += x*f*logf;

    sum_r += r;
    sum_r2 += r*r;
    sum_logr += logr;
    sum_x2r += x*x*r;
    sum_rlogr += r*logr;
    sum_xr += x*r;
    sum_xrlogr += x*r*logr;
  }


  double fit_a1 = 0.0, fit_b1 = 0.0; // see http://mathworld.wolfram.com/LeastSquaresFittingExponential.html
  double fit_a2 = 0.0, fit_b2 = 0.0;
  double fit_a3 = 0.0, fit_b3 = 0.0;
  double denom_1 = 0.0, denom_2 = 0.0, denom_3 = 0.0;
  double aveE1 = 0.0, aveE2 = 0.0, aveE3 = 0.0;
  int N0_1 = 0, N0_2 = 0, N0_3 = 0;

  denom_1 = (sum_y*sum_x2y - sum_xy*sum_xy);
  denom_2 = (sum_f*sum_x2f - sum_xf*sum_xf);
  denom_3 = (sum_r*sum_x2r - sum_xr*sum_xr);
  if (!AreSame(denom_1, 0.0)) {
    fit_a1 = (sum_x2y*sum_ylogy - sum_xy*sum_xylogy) / denom_1;
    fit_b1 = (sum_y*sum_xylogy - sum_xy*sum_ylogy) / denom_1;
    aveE1 = -1. / fit_b1;
    N0_1 = (int)(exp(fit_a1) * aveE1);
  }
  if (!AreSame(denom_2, 0.0)) {
    fit_a2 = (sum_x2f*sum_flogf - sum_xf*sum_xflogf) / denom_2;
    fit_b2 = (sum_f*sum_xflogf - sum_xf*sum_flogf) / denom_2;
    aveE2 = -1. / fit_b2;
    N0_2 = (int)(exp(fit_a2) * aveE2);
  }
  if (!AreSame(denom_3, 0.0)) {
    fit_a3 = (sum_x2r*sum_rlogr - sum_xr*sum_xrlogr) / denom_3;
    fit_b3 = (sum_r*sum_xrlogr - sum_xr*sum_rlogr) / denom_3;
    aveE3 = -1. / fit_b3;
    N0_3 = (int)(exp(fit_a3) * aveE3);
  }


  int weight = 1; // fix, read it from the infile
  int subsample_factor = 1; // fix, read it from the infile


  if (func)
  {
    std::cout << "Fit for full spectrum: y=" << exp(fit_a1) << "e^(" << fit_b1 << "x)" << std::endl;
    std::cout << "Fit for front spectrum: y=" << exp(fit_a2) << "e^(" << fit_b2 << "x)" << std::endl;
    std::cout << "Fit for rear spectrum: y=" << exp(fit_a3) << "e^(" << fit_b3 << "x)" << std::endl;
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
    fprintf(outfile, "AVERAGE_E3 = %3.2f\n", aveE3);
    fprintf(outfile, "WEIGHT = %i\n", weight);
    fprintf(outfile, "SUBSAMPLE = %i\n", subsample_factor);
    fprintf(outfile, "N0_1 = %i*WEIGHT*SUBSAMPLE\n", N0_1);
    fprintf(outfile, "N0_2 = %i*WEIGHT*SUBSAMPLE\n", N0_2);
    fprintf(outfile, "N0_3 = %i*WEIGHT*SUBSAMPLE\n", N0_3);
    fprintf(outfile, "f(x) = (N0_1 / AVERAGE_E1)*exp(-x / AVERAGE_E1)\n");
    fprintf(outfile, "g(x) = (N0_2 / AVERAGE_E2)*exp(-x / AVERAGE_E2)\n");
    fprintf(outfile, "h(x) = (N0_3 / AVERAGE_E3)*exp(-x / AVERAGE_E3)\n");
    fprintf(outfile, "set xlabel 'E (MeV)' \n");
    fprintf(outfile, "set ylabel 'dN/dE (MeV^{-1})'\n");
    fprintf(outfile, "set format y '10^{%%L}'\n");
    //  fprintf(outfile, "set xrange[%i:%i]\n", Emin, Emax);
    fprintf(outfile, "set logscale y\n");
    fprintf(outfile, "plot FILE_IN u 1:($2*%i*%i) w histeps lt 1 lc rgb 'blue' lw 3 t 'full spectrum',\\", weight, subsample_factor);
    fprintf(outfile, "\n");
    fprintf(outfile, "FILE_IN u 1:($3*%i*%i) w histeps lt 1 lc rgb 'red' lw 3 t 'front spectrum',\\", weight, subsample_factor);
    fprintf(outfile, "\n");
    fprintf(outfile, "FILE_IN u 1:($4*%i*%i) w histeps lt 1 lc rgb 'cyan' lw 3 t 'rear spectrum',\\", weight, subsample_factor);
    fprintf(outfile, "\n");
    fprintf(outfile, "f(x) w lines lt 1 lc rgb 'purple' lw 3 t 'exponential fit E_01 = %1.1f MeV',\\", aveE1);
    fprintf(outfile, "\n");
    fprintf(outfile, "g(x) w lines lt 1 lc rgb 'brown' lw 3 t 'exponential fit E_02 = %1.1f MeV',\\", aveE2);
    fprintf(outfile, "\n");
    fprintf(outfile, "h(x) w lines lt 1 lc rgb 'dark-green' lw 3 t 'exponential fit E_03 = %1.1f MeV'", aveE3);

    fprintf(outfile, "\n");
    fclose(outfile);
  }


  if (scan)
  {
    printf(" \t %3.2g \t %i \t %3.2g \t %i \t %3.2g \t %i \t %3.2g \t %3.2g \t %3.2g \t %3.2g \t ", aveE1, N0_1, aveE2, N0_2, aveE3, N0_3, max_energy, tot_energy, tot_energy_front, tot_energy_rear);
  }



  return 0;

}


