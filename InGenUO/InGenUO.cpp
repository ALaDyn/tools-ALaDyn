//
//	"Input Generator User Oriented" for ALaDyn
//
//	Onofrio Mazzarisi (o.mazzarisi@mazzarisi.it)
//  Stefano Sinigardi (stesinigardi@hotmail.com)
//


#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <string>
#include <cstdlib>
#include <sys/stat.h>
#include <errno.h>
#if defined (_MSC_VER)
#include <boost/filesystem.hpp>
#else
#include <unistd.h>
#endif

using namespace std;

string cartella(double j, double k, double l)
{
  char* path = (char*)malloc(sizeof(char) * 100);
#if defined (_MSC_VER)
  std::string pathstring;
  boost::filesystem::path mydir = boost::filesystem::system_complete(pathstring);
  const char* curr_dir = mydir.string().c_str();
#else
  int ris;
  const char* curr_dir = getcwd(NULL, 0);
#endif
  string new_dir_name = "/par_j" + std::to_string(j) + "par_k" + std::to_string(k) + "par_l" + std::to_string(l) + "/";
  strcpy(path, curr_dir);
  strcat(path, new_dir_name.data());
#if defined (_MSC_VER)
  if (!boost::filesystem::exists(path)){
    boost::filesystem::create_directories(path);
  }
#else
  ris = mkdir(path,S_IRWXU);
  free(curr_dir);
#endif
  free(path);
  return string(path);
}



int main()
{

  //
  //	Choose max values and step values for the increments
  //

  double max_j = 2.5;
  double max_k = 2.5;
  double max_l = 2.5;
  double step_j = 0.5;
  double step_k = 0.5;
  double step_l = 0.5;

  for (double j = 0.0, k = 0.0, l = 0.0; j < max_j&&k < max_k&&l < max_l; j += step_j, k += step_k, l += step_l)
  {
    string percorso = cartella(j, k, l);
    ofstream myfile;
    myfile.open((percorso + "Input.nml").data());
    myfile << "\n";
    myfile << "\n";

    //
    //	Add j,k,l to the parameters you want to vary
    //

    myfile << "&GRID \n";

    myfile << " nx = ";          myfile << 256;      myfile << ",\n";
    myfile << " ny = ";          myfile << 128;      myfile << ",\n";
    myfile << " nz = ";          myfile << 1;        myfile << ",\n";
    myfile << " ny_targ = ";     myfile << 120;      myfile << ",\n";
    myfile << " k0 = ";          myfile << 30.;      myfile << ",\n";
    myfile << " yx_rat = ";      myfile << 2.0;      myfile << ",\n";
    myfile << "/ \n\n\n";

    myfile << "&SIMULATION \n";

    myfile << " LPf_ord = ";      myfile << 2;        myfile << ",\n";
    myfile << " Der_ord = ";      myfile << 2;        myfile << ",\n";
    myfile << " Str_flag = ";     myfile << 2;        myfile << ",\n";
    myfile << " iform = ";        myfile << 1;        myfile << ",\n";
    myfile << " model_id = ";     myfile << 1;        myfile << ",\n";
    myfile << " dmodel_id = ";    myfile << 4;        myfile << ",\n";
    myfile << " ibeam = ";        myfile << 1;        myfile << ",\n";
    myfile << " nsp = ";          myfile << 3;        myfile << ",\n";
    myfile << " nsb = ";          myfile << 0;        myfile << ",\n";
    myfile << " Z_ion = ";        myfile << 9;        myfile << ",\n";
    myfile << " A_ion = ";        myfile << 27;       myfile << ",\n";
    myfile << " np_per_xc(1) = "; myfile << 6;        myfile << ",\n";
    myfile << " np_per_xc(2) = "; myfile << 2;        myfile << ",\n";
    myfile << " np_per_xc(3) = "; myfile << 2;        myfile << ",\n";
    myfile << " np_per_xc(4) = "; myfile << 2;        myfile << ",\n";
    myfile << " np_per_xc(5) = "; myfile << 2;        myfile << ",\n";
    myfile << " np_per_xc(6) = "; myfile << 2;        myfile << ",\n";
    myfile << " np_per_yc(1) = "; myfile << 2;        myfile << ",\n";
    myfile << " np_per_yc(2) = "; myfile << 2;        myfile << ",\n";
    myfile << " np_per_yc(3) = "; myfile << 2;        myfile << ",\n";
    myfile << " np_per_yc(4) = "; myfile << 2;        myfile << ",\n";
    myfile << " np_per_yc(5) = "; myfile << 2;        myfile << ",\n";
    myfile << " np_per_yc(6) = "; myfile << 2;        myfile << ",\n";
    myfile << " t0_lp = ";        myfile << 2.0;      myfile << ",\n";
    myfile << " xc_lp = ";        myfile << 1.0;      myfile << ",\n";
    myfile << " w0_x = ";         myfile << 2.0;      myfile << ",\n";
    myfile << " w0_y = ";         myfile << 1.0;      myfile << ",\n";
    myfile << " a0 = ";           myfile << 5.0;      myfile << ",\n";
    myfile << " lam0 = ";         myfile << 0.8;      myfile << ",\n";
    myfile << " lpx(1) = ";       myfile << 0.0;      myfile << ",\n";
    myfile << " lpx(2) = ";       myfile << 0.0;      myfile << ",\n";
    myfile << " lpx(3) = ";       myfile << 0.2;      myfile << ",\n";
    myfile << " lpx(4) = ";       myfile << 0.0;      myfile << ",\n";
    myfile << " lpx(5) = ";       myfile << 0.0;      myfile << ",\n";
    myfile << " lpx(6) = ";       myfile << 0.0;      myfile << ",\n";
    myfile << " lpx(7) = ";       myfile << 0.01;     myfile << ",\n";
    myfile << " n_over_nc = ";    myfile << 100.;     myfile << ",\n";
    myfile << " n1_over_n = ";    myfile << 10.0;     myfile << ",\n";
    myfile << " n2_over_n = ";    myfile << 10.0;     myfile << ",\n";
    myfile << "w_sh = ";          myfile << 20;       myfile << ",\n";
    myfile << "wi_time = ";       myfile << 120.;     myfile << ",\n";
    myfile << "wf_time = ";       myfile << 120.;     myfile << ",\n";
    myfile << "w_speed = ";       myfile << 1.0;      myfile << ",\n";
    myfile << "/ \n\n";

    myfile << "&OUTPUT \n";

    myfile << " nouts = ";        myfile << 1;        myfile << ",\n";
    myfile << " iene = ";         myfile << 10;       myfile << ",\n";
    myfile << " nvout = ";        myfile << 0;        myfile << ",\n";
    myfile << " nden = ";         myfile << 0;        myfile << ",\n";
    myfile << " npout = ";        myfile << 4;        myfile << ",\n";
    myfile << " nbout = ";        myfile << 0;        myfile << ",\n";
    myfile << " jump = ";         myfile << 1;        myfile << ",\n";
    myfile << " pjump = ";        myfile << 1;        myfile << ",\n";
    myfile << " xp0_out = ";      myfile << 0.0;      myfile << ",\n";
    myfile << " xp1_out = ";      myfile << 100.0, myfile << ",\n";
    myfile << " yp_out = ";       myfile << 20.0;     myfile << ",\n";
    myfile << " tmax = ";         myfile << 1.0;      myfile << ",\n";
    myfile << " cfl = ";          myfile << 0.85;     myfile << ",\n";
    myfile << " new_sim = ";      myfile << 0;        myfile << ",\n";
    myfile << " id_new = ";       myfile << 0;        myfile << ",\n";
    myfile << " dump = ";         myfile << 0;        myfile << ",\n";
    myfile << " npe_yz = ";       myfile << 2;        myfile << ",\n";
    myfile << "/ \n";


    myfile.close();
  }
}
