//
//  "Input Generator User Oriented" for ALaDyn
//  Converts a json input file to a namelist input file
//  It is the basis for a more user-oriented input to ALaDyn
//  For now, the json and the nml are just the representation one of the other, will evolve in the future
//
//  Copyright 2014-2016
//  Stefano Sinigardi (stesinigardi@hotmail.com)
//  Onofrio Mazzarisi (o.mazzarisi@mazzarisi.it)
//


#define _CRT_SECURE_NO_WARNINGS
#define MAX_SIZE_OF_PATH 1024

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
#include "jsoncons/json.hpp"


class ALaDyn_parameters {
public:
  ALaDyn_parameters();
  unsigned int nx, ny, nz, ny_targ;
  double k0, yx_rat, zx_rat;
  unsigned short lpf_ord, der_ord, str_flag, iform, model_id, dmodel_id, ibeam, nsp, nsb;
  unsigned short ibx, iby, ibz;
  unsigned short ionz_lev, ionz_model;
  unsigned short ion_min1, ion_min2, ion_min3;
  unsigned short ion_max1, ion_max2, ion_max3;
  unsigned short atomic_number1, atomic_number2, atomic_number3;
  double mass_number1, mass_number2, mass_number3;
  double t0_pl1, t0_pl2, t0_pl3, t0_pl4;
  unsigned short np_per_xc1, np_per_xc2, np_per_xc3, np_per_xc4, np_per_xc5, np_per_xc6;
  unsigned short np_per_yc1, np_per_yc2, np_per_yc3, np_per_yc4, np_per_yc5, np_per_yc6;
  double t0_lp, xc_lp, w0_x, w0_y, a0, lam0;
  double lpx1, lpx2, lpx3, lpx4, lpx5, lpx6, lpx7;
  double n_over_nc, n1_over_n, n2_over_n;
  unsigned short w_sh;
  double wi_time, wf_time, w_speed;
  unsigned int nouts, iene;
  unsigned short nvout, nden, npout, nbout, jump, pjump;
  double xp0_out, xp1_out, yp_out;
  double tmax;
  double cfl;
  unsigned short new_sim, dump;
  unsigned int id_new;
  unsigned int nprocx, nprocy, nprocz;
};


ALaDyn_parameters::ALaDyn_parameters() {
  nx = 256, ny = 128, nz = 1, ny_targ = 120;
  k0 = 30.0, yx_rat = 2.0, zx_rat = 2.0;
  lpf_ord = 2, der_ord = 2, str_flag = 2, iform = 1, model_id = 1, dmodel_id = 4, ibeam = 1, nsp = 3, nsb = 0;
  ibx = 0, iby = 0, ibz = 0;
  ionz_lev = 0, ionz_model = 1;
  ion_min1 = 9, ion_min2 = 1, ion_min3 = 1;
  ion_max1 = 9, ion_max2 = 1, ion_max3 = 1;
  atomic_number1 = 13, atomic_number2 = 1, atomic_number3 = 1;
  mass_number1 = 27.0, mass_number2 = 1.0, mass_number3 = 1.0;
  t0_pl1 = 0.003, t0_pl2 = 0.0, t0_pl3 = 0.0, t0_pl4 = 0.0;
  np_per_xc1 = 6, np_per_xc2 = 2, np_per_xc3 = 2, np_per_xc4 = 2, np_per_xc5 = 2, np_per_xc6 = 2;
  np_per_yc1 = 2, np_per_yc2 = 2, np_per_yc3 = 2, np_per_yc4 = 2, np_per_yc5 = 2, np_per_yc6 = 2;
  t0_lp = 2.0, xc_lp = 1.0, w0_x = 2.0, w0_y = 1.0, a0 = 5.0, lam0 = 0.8;
  lpx1 = 0.0, lpx2 = 0.0, lpx3 = 0.2, lpx4 = 0.0, lpx5 = 0.0, lpx6 = 0.0, lpx7 = 0.01;
  n_over_nc = 100.0, n1_over_n = 10.0, n2_over_n = 10.0;
  w_sh = 20;
  wi_time = 120.0, wf_time = 120.0, w_speed = 1.0;
  nouts = 1, iene = 10;
  nvout = 0, nden = 0, npout = 4, nbout = 0, jump = 1, pjump = 1;
  xp0_out = 0.0, xp1_out = 100.0, yp_out = 20.0;
  tmax = 1.0;
  cfl = 0.85;
  new_sim = 0, dump = 0;
  id_new = 0;
  nprocx = 1, nprocy = 2, nprocz = 1;
}

void parse_json_file(ALaDyn_parameters& parameters, const char * filename) {
  jsoncons::json json_parameters = jsoncons::json::parse_file(filename);

  if (json_parameters.has_member("nx")) parameters.nx = json_parameters["nx"].as<unsigned int>();
  else std::cout << "Missing nx definition, default to " << parameters.nx << "\n";

  if (json_parameters.has_member("ny")) parameters.ny = json_parameters["ny"].as<unsigned int>();
  else std::cout << "Missing ny definition, default to " << parameters.ny << "\n";

  if (json_parameters.has_member("nz")) parameters.nz = json_parameters["nz"].as<unsigned int>();
  else std::cout << "Missing nz definition, default to " << parameters.nz << "\n";

  if (json_parameters.has_member("ny_targ")) parameters.ny_targ = json_parameters["ny_targ"].as<unsigned int>();
  else std::cout << "Missing ny_targ definition, default to " << parameters.ny_targ << "\n";

  if (json_parameters.has_member("k0")) parameters.k0 = json_parameters["k0"].as<double>();
  else std::cout << "Missing k0 definition, default to " << parameters.k0 << "\n";

  if (json_parameters.has_member("yx_rat")) parameters.yx_rat = json_parameters["yx_rat"].as<double>();
  else std::cout << "Missing yx_rat definition, default to " << parameters.yx_rat << "\n";

  if (json_parameters.has_member("zx_rat")) parameters.zx_rat = json_parameters["zx_rat"].as<double>();
  else std::cout << "Missing zx_rat definition, default to " << parameters.zx_rat << "\n";

  if (json_parameters.has_member("lpf_ord")) parameters.lpf_ord = json_parameters["lpf_ord"].as<unsigned short>();
  else std::cout << "Missing lpf_ord definition, default to " << parameters.lpf_ord << "\n";

  if (json_parameters.has_member("der_ord")) parameters.der_ord = json_parameters["der_ord"].as<unsigned short>();
  else std::cout << "Missing der_ord definition, default to " << parameters.der_ord << "\n";

  if (json_parameters.has_member("str_flag")) parameters.str_flag = json_parameters["str_flag"].as<unsigned short>();
  else std::cout << "Missing str_flag definition, default to " << parameters.str_flag << "\n";

  if (json_parameters.has_member("iform")) parameters.iform = json_parameters["iform"].as<unsigned short>();
  else std::cout << "Missing iform definition, default to " << parameters.iform << "\n";

  if (json_parameters.has_member("model_id")) parameters.model_id = json_parameters["model_id"].as<unsigned short>();
  else std::cout << "Missing model_id definition, default to " << parameters.model_id << "\n";

  if (json_parameters.has_member("dmodel_id")) parameters.dmodel_id = json_parameters["dmodel_id"].as<unsigned short>();
  else std::cout << "Missing dmodel_id definition, default to " << parameters.dmodel_id << "\n";

  if (json_parameters.has_member("ibx")) parameters.ibx = json_parameters["ibx"].as<unsigned short>();
  else std::cout << "Missing ibx definition, default to " << parameters.ibx << "\n";

  if (json_parameters.has_member("iby")) parameters.iby = json_parameters["iby"].as<unsigned short>();
  else std::cout << "Missing iby definition, default to " << parameters.iby << "\n";

  if (json_parameters.has_member("ibz")) parameters.ibz = json_parameters["ibz"].as<unsigned short>();
  else std::cout << "Missing ibz definition, default to " << parameters.ibz << "\n";

  if (json_parameters.has_member("ibeam")) parameters.ibeam = json_parameters["ibeam"].as<unsigned short>();
  else std::cout << "Missing ibeam definition, default to " << parameters.ibeam << "\n";

  if (json_parameters.has_member("nsp")) parameters.nsp = json_parameters["nsp"].as<unsigned short>();
  else std::cout << "Missing nsp definition, default to " << parameters.nsp << "\n";

  if (json_parameters.has_member("nsb")) parameters.nsb = json_parameters["nsb"].as<unsigned short>();
  else std::cout << "Missing nsb definition, default to " << parameters.nsb << "\n";

  if (json_parameters.has_member("ionz_lev")) parameters.ionz_lev = json_parameters["ionz_lev"].as<unsigned short>();
  else std::cout << "Missing ionz_lev definition, default to " << parameters.ionz_lev << "\n";

  if (json_parameters.has_member("ionz_model")) parameters.ionz_model = json_parameters["ionz_model"].as<unsigned short>();
  else std::cout << "Missing ionz_model definition, default to " << parameters.ionz_model << "\n";

  if (json_parameters.has_member("ion_min1")) parameters.ion_min1 = json_parameters["ion_min1"].as<unsigned short>();
  else std::cout << "Missing ion_min1 definition, default to " << parameters.ion_min1 << "\n";

  if (json_parameters.has_member("ion_min2")) parameters.ion_min2 = json_parameters["ion_min2"].as<unsigned short>();
  else std::cout << "Missing ion_min2 definition, default to " << parameters.ion_min2 << "\n";

  if (json_parameters.has_member("ion_min3")) parameters.ion_min3 = json_parameters["ion_min3"].as<unsigned short>();
  else std::cout << "Missing ion_min3 definition, default to " << parameters.ion_min3 << "\n";

  if (json_parameters.has_member("ion_max1")) parameters.ion_max1 = json_parameters["ion_max1"].as<unsigned short>();
  else std::cout << "Missing ion_max1 definition, default to " << parameters.ion_max1 << "\n";

  if (json_parameters.has_member("ion_max2")) parameters.ion_max2 = json_parameters["ion_max2"].as<unsigned short>();
  else std::cout << "Missing ion_max2 definition, default to " << parameters.ion_max2 << "\n";

  if (json_parameters.has_member("ion_max3")) parameters.ion_max3 = json_parameters["ion_max3"].as<unsigned short>();
  else std::cout << "Missing ion_max3 definition, default to " << parameters.ion_max3 << "\n";

  if (json_parameters.has_member("atomic_number1")) parameters.atomic_number1 = json_parameters["atomic_number1"].as<unsigned short>();
  else std::cout << "Missing atomic_number1 definition, default to " << parameters.atomic_number1 << "\n";

  if (json_parameters.has_member("atomic_number2")) parameters.atomic_number2 = json_parameters["atomic_number2"].as<unsigned short>();
  else std::cout << "Missing atomic_number2 definition, default to " << parameters.atomic_number2 << "\n";

  if (json_parameters.has_member("atomic_number3")) parameters.atomic_number3 = json_parameters["atomic_number3"].as<unsigned short>();
  else std::cout << "Missing atomic_number3 definition, default to " << parameters.atomic_number3 << "\n";

  if (json_parameters.has_member("mass_number1")) parameters.mass_number1 = json_parameters["mass_number1"].as<double>();
  else std::cout << "Missing mass_number1 definition, default to " << parameters.mass_number1 << "\n";

  if (json_parameters.has_member("mass_number2")) parameters.mass_number2 = json_parameters["mass_number2"].as<double>();
  else std::cout << "Missing mass_number2 definition, default to " << parameters.mass_number2 << "\n";

  if (json_parameters.has_member("mass_number3")) parameters.mass_number3 = json_parameters["mass_number3"].as<double>();
  else std::cout << "Missing mass_number3 definition, default to " << parameters.mass_number3 << "\n";

  if (json_parameters.has_member("t0_pl1")) parameters.t0_pl1 = json_parameters["t0_pl1"].as<double>();
  else std::cout << "Missing t0_pl1 definition, default to " << parameters.t0_pl1 << "\n";

  if (json_parameters.has_member("t0_pl2")) parameters.t0_pl2 = json_parameters["t0_pl2"].as<double>();
  else std::cout << "Missing t0_pl2 definition, default to " << parameters.t0_pl2 << "\n";

  if (json_parameters.has_member("t0_pl3")) parameters.t0_pl3 = json_parameters["t0_pl3"].as<double>();
  else std::cout << "Missing t0_pl3 definition, default to " << parameters.t0_pl3 << "\n";

  if (json_parameters.has_member("t0_pl4")) parameters.t0_pl4 = json_parameters["t0_pl4"].as<double>();
  else std::cout << "Missing t0_pl4 definition, default to " << parameters.t0_pl4 << "\n";

  if (json_parameters.has_member("np_per_xc1")) parameters.np_per_xc1 = json_parameters["np_per_xc1"].as<unsigned short>();
  else std::cout << "Missing np_per_xc1 definition, default to " << parameters.np_per_xc1 << "\n";

  if (json_parameters.has_member("np_per_xc2")) parameters.np_per_xc2 = json_parameters["np_per_xc2"].as<unsigned short>();
  else std::cout << "Missing np_per_xc2 definition, default to " << parameters.np_per_xc2 << "\n";

  if (json_parameters.has_member("np_per_xc3")) parameters.np_per_xc3 = json_parameters["np_per_xc3"].as<unsigned short>();
  else std::cout << "Missing np_per_xc3 definition, default to " << parameters.np_per_xc3 << "\n";

  if (json_parameters.has_member("np_per_xc4")) parameters.np_per_xc4 = json_parameters["np_per_xc4"].as<unsigned short>();
  else std::cout << "Missing np_per_xc4 definition, default to " << parameters.np_per_xc4 << "\n";

  if (json_parameters.has_member("np_per_xc5")) parameters.np_per_xc5 = json_parameters["np_per_xc5"].as<unsigned short>();
  else std::cout << "Missing np_per_xc5 definition, default to " << parameters.np_per_xc5 << "\n";

  if (json_parameters.has_member("np_per_xc6")) parameters.np_per_xc6 = json_parameters["np_per_xc6"].as<unsigned short>();
  else std::cout << "Missing np_per_xc6 definition, default to " << parameters.np_per_xc6 << "\n";

  if (json_parameters.has_member("np_per_yc1")) parameters.np_per_yc1 = json_parameters["np_per_yc1"].as<unsigned short>();
  else std::cout << "Missing np_per_yc1 definition, default to " << parameters.np_per_yc1 << "\n";

  if (json_parameters.has_member("np_per_yc2")) parameters.np_per_yc2 = json_parameters["np_per_yc2"].as<unsigned short>();
  else std::cout << "Missing np_per_yc2 definition, default to " << parameters.np_per_yc2 << "\n";

  if (json_parameters.has_member("np_per_yc3")) parameters.np_per_yc3 = json_parameters["np_per_yc3"].as<unsigned short>();
  else std::cout << "Missing np_per_yc3 definition, default to " << parameters.np_per_yc3 << "\n";

  if (json_parameters.has_member("np_per_yc4")) parameters.np_per_yc4 = json_parameters["np_per_yc4"].as<unsigned short>();
  else std::cout << "Missing np_per_yc4 definition, default to " << parameters.np_per_yc4 << "\n";

  if (json_parameters.has_member("np_per_yc5")) parameters.np_per_yc5 = json_parameters["np_per_yc5"].as<unsigned short>();
  else std::cout << "Missing np_per_yc5 definition, default to " << parameters.np_per_yc5 << "\n";

  if (json_parameters.has_member("np_per_yc6")) parameters.np_per_yc6 = json_parameters["np_per_yc6"].as<unsigned short>();
  else std::cout << "Missing np_per_yc6 definition, default to " << parameters.np_per_yc6 << "\n";

  if (json_parameters.has_member("t0_lp")) parameters.t0_lp = json_parameters["t0_lp"].as<double>();
  else std::cout << "Missing t0_lp definition, default to " << parameters.t0_lp << "\n";

  if (json_parameters.has_member("xc_lp")) parameters.xc_lp = json_parameters["xc_lp"].as<double>();
  else std::cout << "Missing xc_lp definition, default to " << parameters.xc_lp << "\n";

  if (json_parameters.has_member("w0_x")) parameters.w0_x = json_parameters["w0_x"].as<double>();
  else std::cout << "Missing w0_x definition, default to " << parameters.w0_x << "\n";

  if (json_parameters.has_member("w0_y")) parameters.w0_y = json_parameters["w0_y"].as<double>();
  else std::cout << "Missing w0_y definition, default to " << parameters.w0_y << "\n";

  if (json_parameters.has_member("a0")) parameters.a0 = json_parameters["a0"].as<double>();
  else std::cout << "Missing a0 definition, default to " << parameters.a0 << "\n";

  if (json_parameters.has_member("lam0")) parameters.lam0 = json_parameters["lam0"].as<double>();
  else std::cout << "Missing lam0 definition, default to " << parameters.lam0 << "\n";

  if (json_parameters.has_member("lpx1")) parameters.lpx1 = json_parameters["lpx1"].as<double>();
  else std::cout << "Missing lpx1 definition, default to " << parameters.lpx1 << "\n";

  if (json_parameters.has_member("lpx2")) parameters.lpx2 = json_parameters["lpx2"].as<double>();
  else std::cout << "Missing lpx2 definition, default to " << parameters.lpx2 << "\n";

  if (json_parameters.has_member("lpx3")) parameters.lpx3 = json_parameters["lpx3"].as<double>();
  else std::cout << "Missing lpx3 definition, default to " << parameters.lpx3 << "\n";

  if (json_parameters.has_member("lpx4")) parameters.lpx4 = json_parameters["lpx4"].as<double>();
  else std::cout << "Missing lpx4 definition, default to " << parameters.lpx4 << "\n";

  if (json_parameters.has_member("lpx5")) parameters.lpx5 = json_parameters["lpx5"].as<double>();
  else std::cout << "Missing lpx5 definition, default to " << parameters.lpx5 << "\n";

  if (json_parameters.has_member("lpx6")) parameters.lpx6 = json_parameters["lpx6"].as<double>();
  else std::cout << "Missing lpx6 definition, default to " << parameters.lpx6 << "\n";

  if (json_parameters.has_member("lpx7")) parameters.lpx7 = json_parameters["lpx7"].as<double>();
  else std::cout << "Missing lpx7 definition, default to " << parameters.lpx7 << "\n";

  if (json_parameters.has_member("n_over_nc")) parameters.n_over_nc = json_parameters["n_over_nc"].as<double>();
  else std::cout << "Missing n_over_nc definition, default to " << parameters.n_over_nc << "\n";

  if (json_parameters.has_member("n1_over_n")) parameters.n1_over_n = json_parameters["n1_over_n"].as<double>();
  else std::cout << "Missing n1_over_n definition, default to " << parameters.n1_over_n << "\n";

  if (json_parameters.has_member("n2_over_n")) parameters.n2_over_n = json_parameters["n2_over_n"].as<double>();
  else std::cout << "Missing n2_over_n definition, default to " << parameters.n2_over_n << "\n";

  if (json_parameters.has_member("w_sh")) parameters.w_sh = json_parameters["w_sh"].as<unsigned short>();
  else std::cout << "Missing w_sh definition, default to " << parameters.w_sh << "\n";

  if (json_parameters.has_member("wi_time")) parameters.wi_time = json_parameters["wi_time"].as<double>();
  else std::cout << "Missing wi_time definition, default to " << parameters.wi_time << "\n";

  if (json_parameters.has_member("wf_time")) parameters.wf_time = json_parameters["wf_time"].as<double>();
  else std::cout << "Missing wf_time definition, default to " << parameters.wf_time << "\n";

  if (json_parameters.has_member("w_speed")) parameters.w_speed = json_parameters["w_speed"].as<double>();
  else std::cout << "Missing w_speed definition, default to " << parameters.w_speed << "\n";

  if (json_parameters.has_member("nouts")) parameters.nouts = json_parameters["nouts"].as<unsigned int>();
  else std::cout << "Missing nouts definition, default to " << parameters.nouts << "\n";

  if (json_parameters.has_member("iene")) parameters.iene = json_parameters["iene"].as<unsigned int>();
  else std::cout << "Missing iene definition, default to " << parameters.iene << "\n";

  if (json_parameters.has_member("nvout")) parameters.nvout = json_parameters["nvout"].as<unsigned short>();
  else std::cout << "Missing nvout definition, default to " << parameters.nvout << "\n";

  if (json_parameters.has_member("nden")) parameters.nden = json_parameters["nden"].as<unsigned short>();
  else std::cout << "Missing nden definition, default to " << parameters.nden << "\n";

  if (json_parameters.has_member("npout")) parameters.npout = json_parameters["npout"].as<unsigned short>();
  else std::cout << "Missing npout definition, default to " << parameters.npout << "\n";

  if (json_parameters.has_member("nbout")) parameters.nbout = json_parameters["nbout"].as<unsigned short>();
  else std::cout << "Missing nbout definition, default to " << parameters.nbout << "\n";

  if (json_parameters.has_member("jump")) parameters.jump = json_parameters["jump"].as<unsigned short>();
  else std::cout << "Missing jump definition, default to " << parameters.jump << "\n";

  if (json_parameters.has_member("pjump")) parameters.pjump = json_parameters["pjump"].as<unsigned short>();
  else std::cout << "Missing pjump definition, default to " << parameters.pjump << "\n";

  if (json_parameters.has_member("xp0_out")) parameters.xp0_out = json_parameters["xp0_out"].as<double>();
  else std::cout << "Missing xp0_out definition, default to " << parameters.xp0_out << "\n";

  if (json_parameters.has_member("xp1_out")) parameters.xp1_out = json_parameters["xp1_out"].as<double>();
  else std::cout << "Missing xp1_out definition, default to " << parameters.xp1_out << "\n";

  if (json_parameters.has_member("yp_out")) parameters.yp_out = json_parameters["yp_out"].as<double>();
  else std::cout << "Missing yp_out definition, default to " << parameters.yp_out << "\n";

  if (json_parameters.has_member("tmax")) parameters.tmax = json_parameters["tmax"].as<double>();
  else std::cout << "Missing tmax definition, default to " << parameters.tmax << "\n";

  if (json_parameters.has_member("cfl")) parameters.cfl = json_parameters["cfl"].as<double>();
  else std::cout << "Missing cfl definition, default to " << parameters.cfl << "\n";

  if (json_parameters.has_member("new_sim")) parameters.new_sim = json_parameters["new_sim"].as<unsigned short>();
  else std::cout << "Missing new_sim definition, default to " << parameters.new_sim << "\n";

  if (json_parameters.has_member("id_new")) parameters.id_new = json_parameters["id_new"].as<unsigned int>();
  else std::cout << "Missing id_new definition, default to " << parameters.id_new << "\n";

  if (json_parameters.has_member("dump")) parameters.dump = json_parameters["dump"].as<unsigned short>();
  else std::cout << "Missing dump definition, default to " << parameters.dump << "\n";

  if (json_parameters.has_member("nprocx")) parameters.nprocx = json_parameters["nprocx"].as<unsigned int>();
  else std::cout << "Missing nprocx definition, default to " << parameters.nprocx << "\n";

  if (json_parameters.has_member("nprocy")) parameters.nprocy = json_parameters["nprocy"].as<unsigned int>();
  else std::cout << "Missing nprocy definition, default to " << parameters.nprocy << "\n";

  if (json_parameters.has_member("nprocz")) parameters.nprocz = json_parameters["nprocz"].as<unsigned int>();
  else std::cout << "Missing nprocz definition, default to " << parameters.nprocz << "\n";
}

std::string cartella(std::string prefisso)
{
  char* path = (char*)malloc(sizeof(char) * MAX_SIZE_OF_PATH);
#if defined (_MSC_VER)
  std::string pathstring;
  boost::filesystem::path mydir = boost::filesystem::system_complete(pathstring);
  const char* curr_dir = mydir.string().c_str();
#else
  int ris;
  const char* curr_dir = getcwd(NULL, 0);
#endif
  strcpy(path, curr_dir);
  strcat(path, prefisso.data());
#if defined (_MSC_VER)
  if (!boost::filesystem::exists(path)) {
    boost::filesystem::create_directories(path);
  }
#else
  ris = mkdir(path, S_IRWXU);
#endif
  free(path);
  return std::string(path);
}

void print_input_nml(ALaDyn_parameters &parameters, std::string filename) {
  std::ofstream myfile;
  myfile.open(filename);

  myfile << "&GRID \n";
  myfile << " nx = ";               myfile << parameters.nx;              myfile << ",\n";
  myfile << " ny = ";               myfile << parameters.ny;              myfile << ",\n";
  myfile << " nz = ";               myfile << parameters.nz;              myfile << ",\n";
  myfile << " ny_targ = ";          myfile << parameters.ny_targ;         myfile << ",\n";
  myfile << " k0 = ";               myfile << parameters.k0;              myfile << ",\n";
  myfile << " yx_rat = ";           myfile << parameters.yx_rat;          myfile << ",\n";
  myfile << " zx_rat = ";           myfile << parameters.zx_rat;          myfile << "\n";
  myfile << "/ \n\n\n";

  myfile << "&SIMULATION \n";
  myfile << " lpf_ord = ";          myfile << parameters.lpf_ord;         myfile << ",\n";
  myfile << " der_ord = ";          myfile << parameters.der_ord;         myfile << ",\n";
  myfile << " str_flag = ";         myfile << parameters.str_flag;        myfile << ",\n";
  myfile << " iform = ";            myfile << parameters.iform;           myfile << ",\n";
  myfile << " model_id = ";         myfile << parameters.model_id;        myfile << ",\n";
  myfile << " dmodel_id = ";        myfile << parameters.dmodel_id;       myfile << ",\n";
  myfile << " ibx = ";              myfile << parameters.ibx;             myfile << ",\n";
  myfile << " iby = ";              myfile << parameters.iby;             myfile << ",\n";
  myfile << " ibz = ";              myfile << parameters.ibz;             myfile << ",\n";
  myfile << " ibeam = ";            myfile << parameters.ibeam;           myfile << "\n";
  myfile << "/ \n\n";

  myfile << "&TARGET_DESCRIPTION \n";
  myfile << " nsp = ";              myfile << parameters.nsp;             myfile << ",\n";
  myfile << " nsb = ";              myfile << parameters.nsb;             myfile << ",\n";
  myfile << " ionz_lev = ";         myfile << parameters.ionz_lev;        myfile << ",\n";
  myfile << " ionz_model = ";       myfile << parameters.ionz_model;      myfile << ",\n";
  myfile << " ion_min(1) = ";       myfile << parameters.ion_min1;        myfile << ",\n";
  myfile << " ion_min(2) = ";       myfile << parameters.ion_min2;        myfile << ",\n";
  myfile << " ion_min(3) = ";       myfile << parameters.ion_min3;        myfile << ",\n";
  myfile << " ion_max(1) = ";       myfile << parameters.ion_max1;        myfile << ",\n";
  myfile << " ion_max(2) = ";       myfile << parameters.ion_max2;        myfile << ",\n";
  myfile << " ion_max(3) = ";       myfile << parameters.ion_max3;        myfile << ",\n";
  myfile << " atomic_number(1) = "; myfile << parameters.atomic_number1;  myfile << ",\n";
  myfile << " atomic_number(2) = "; myfile << parameters.atomic_number2;  myfile << ",\n";
  myfile << " atomic_number(3) = "; myfile << parameters.atomic_number3;  myfile << ",\n";
  myfile << " mass_number(1) = ";   myfile << parameters.mass_number1;    myfile << ",\n";
  myfile << " mass_number(2) = ";   myfile << parameters.mass_number2;    myfile << ",\n";
  myfile << " mass_number(3) = ";   myfile << parameters.mass_number3;    myfile << ",\n";
  myfile << " t0_pl(1) = ";         myfile << parameters.t0_pl1;          myfile << ",\n";
  myfile << " t0_pl(2) = ";         myfile << parameters.t0_pl2;          myfile << ",\n";
  myfile << " t0_pl(3) = ";         myfile << parameters.t0_pl3;          myfile << ",\n";
  myfile << " t0_pl(4) = ";         myfile << parameters.t0_pl4;          myfile << ",\n";
  myfile << " np_per_xc(1) = ";     myfile << parameters.np_per_xc1;      myfile << ",\n";
  myfile << " np_per_xc(2) = ";     myfile << parameters.np_per_xc2;      myfile << ",\n";
  myfile << " np_per_xc(3) = ";     myfile << parameters.np_per_xc3;      myfile << ",\n";
  myfile << " np_per_xc(4) = ";     myfile << parameters.np_per_xc4;      myfile << ",\n";
  myfile << " np_per_xc(5) = ";     myfile << parameters.np_per_xc5;      myfile << ",\n";
  myfile << " np_per_xc(6) = ";     myfile << parameters.np_per_xc6;      myfile << ",\n";
  myfile << " np_per_yc(1) = ";     myfile << parameters.np_per_yc1;      myfile << ",\n";
  myfile << " np_per_yc(2) = ";     myfile << parameters.np_per_yc2;      myfile << ",\n";
  myfile << " np_per_yc(3) = ";     myfile << parameters.np_per_yc3;      myfile << ",\n";
  myfile << " np_per_yc(4) = ";     myfile << parameters.np_per_yc4;      myfile << ",\n";
  myfile << " np_per_yc(5) = ";     myfile << parameters.np_per_yc5;      myfile << ",\n";
  myfile << " np_per_yc(6) = ";     myfile << parameters.np_per_yc6;      myfile << ",\n";
  myfile << " lpx(1) = ";           myfile << parameters.lpx1;            myfile << ",\n";
  myfile << " lpx(2) = ";           myfile << parameters.lpx2;            myfile << ",\n";
  myfile << " lpx(3) = ";           myfile << parameters.lpx3;            myfile << ",\n";
  myfile << " lpx(4) = ";           myfile << parameters.lpx4;            myfile << ",\n";
  myfile << " lpx(5) = ";           myfile << parameters.lpx5;            myfile << ",\n";
  myfile << " lpx(6) = ";           myfile << parameters.lpx6;            myfile << ",\n";
  myfile << " lpx(7) = ";           myfile << parameters.lpx7;            myfile << ",\n";
  myfile << " n_over_nc = ";        myfile << parameters.n_over_nc;       myfile << ",\n";
  myfile << " n1_over_n = ";        myfile << parameters.n1_over_n;       myfile << ",\n";
  myfile << " n2_over_n = ";        myfile << parameters.n2_over_n;       myfile << "\n";
  myfile << "/ \n\n";

  myfile << "&LASER \n";
  myfile << " t0_lp = ";            myfile << parameters.t0_lp;           myfile << ",\n";
  myfile << " xc_lp = ";            myfile << parameters.xc_lp;           myfile << ",\n";
  myfile << " w0_x = ";             myfile << parameters.w0_x;            myfile << ",\n";
  myfile << " w0_y = ";             myfile << parameters.w0_y;            myfile << ",\n";
  myfile << " a0 = ";               myfile << parameters.a0;              myfile << ",\n";
  myfile << " lam0 = ";             myfile << parameters.lam0;            myfile << "\n";
  myfile << "/ \n\n";

  myfile << "&MOVING_WINDOW \n";
  myfile << " w_sh = ";             myfile << parameters.w_sh;            myfile << ",\n";
  myfile << " wi_time = ";          myfile << parameters.wi_time;         myfile << ",\n";
  myfile << " wf_time = ";          myfile << parameters.wf_time;         myfile << ",\n";
  myfile << " w_speed = ";          myfile << parameters.w_speed;         myfile << "\n";
  myfile << "/ \n\n";

  myfile << "&OUTPUT \n";
  myfile << " nouts = ";            myfile << parameters.npout;           myfile << ",\n";
  myfile << " iene = ";             myfile << parameters.iene;            myfile << ",\n";
  myfile << " nvout = ";            myfile << parameters.nvout;           myfile << ",\n";
  myfile << " nden = ";             myfile << parameters.nden;            myfile << ",\n";
  myfile << " npout = ";            myfile << parameters.npout;           myfile << ",\n";
  myfile << " nbout = ";            myfile << parameters.nbout;           myfile << ",\n";
  myfile << " jump = ";             myfile << parameters.jump;            myfile << ",\n";
  myfile << " pjump = ";            myfile << parameters.pjump;           myfile << ",\n";
  myfile << " xp0_out = ";          myfile << parameters.xp0_out;         myfile << ",\n";
  myfile << " xp1_out = ";          myfile << parameters.xp1_out;         myfile << ",\n";
  myfile << " yp_out = ";           myfile << parameters.yp_out;          myfile << ",\n";
  myfile << " tmax = ";             myfile << parameters.tmax;            myfile << ",\n";
  myfile << " cfl = ";              myfile << parameters.cfl;             myfile << ",\n";
  myfile << " new_sim = ";          myfile << parameters.new_sim;         myfile << ",\n";
  myfile << " id_new = ";           myfile << parameters.id_new;          myfile << ",\n";
  myfile << " dump = ";             myfile << parameters.dump;            myfile << "\n";
  myfile << "/ \n\n";

  myfile << "&MPIPARAMS \n";
  myfile << " nprocx = ";           myfile << parameters.nprocx;          myfile << ",\n";
  myfile << " nprocy = ";           myfile << parameters.nprocy;          myfile << ",\n";
  myfile << " nprocz = ";           myfile << parameters.nprocz;          myfile << "\n";
  myfile << "/ \n";


  myfile.close();
}



int main()
{
  ALaDyn_parameters parametri;

  parse_json_file(parametri, "input.json");

  std::string percorso = cartella("test");
  std::string inputnml = percorso + "input.nml";

  print_input_nml(parametri, inputnml);
}


