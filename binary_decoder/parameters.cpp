
#include "parameters.h"
#include "swap_tools.h"


std::string histo::get_filename_out() {
  return basename + "_histogram_" + data_to_bin + ".txt";
}

void histo::write_binned_data() {
  std::ofstream file_out;
  file_out.open(get_filename_out());
  aladyn_float minT = min_value - get_bin_size();
  aladyn_float maxT = min_value;
  for (int i = 0; i < nbin; i++)
  {
    file_out << std::setprecision(6) << minT << "\t" << maxT << "\t" << data[i] << std::endl;
    minT += get_bin_size();
    maxT += get_bin_size();
  }
  file_out.close();
}

aladyn_float histo::get_bin_size() {
  return (max_value - min_value) / nbin;
}

size_t histo::get_column_to_bin() {
  size_t bin_on_x = std::numeric_limits<size_t>::signaling_NaN();;
  if (data_to_bin == "x") bin_on_x = COLUMN_X;
  else if (data_to_bin == "y") bin_on_x = COLUMN_Y;
  else if (data_to_bin == "z") bin_on_x = COLUMN_Z;
  else if (data_to_bin == "px") bin_on_x = COLUMN_PX;
  else if (data_to_bin == "py") bin_on_x = COLUMN_PY;
  else if (data_to_bin == "pz") bin_on_x = COLUMN_PZ;
  else if (data_to_bin == "gamma") bin_on_x = COLUMN_GAMMA;
  else if (data_to_bin == "theta") bin_on_x = COLUMN_THETA;
  else if (data_to_bin == "E") bin_on_x = COLUMN_E;
  else if (data_to_bin == "thetat") bin_on_x = COLUMN_THETAT;
  else if (data_to_bin == "ty") bin_on_x = COLUMN_TY;
  else if (data_to_bin == "tz") bin_on_x = COLUMN_TZ;
  else if (data_to_bin == "w") bin_on_x = COLUMN_W;
  else if (data_to_bin == "ch") bin_on_x = COLUMN_CH;

  return bin_on_x;
}

histo::histo() {}

histo::histo(std::string file_basename, std::string axis_name, size_t bin_counter, aladyn_float mint, aladyn_float maxt) {
  basename = file_basename;
  data_to_bin = axis_name;
  nbin = bin_counter;
  min_value = mint;
  max_value = maxt;
  data = new aladyn_float[nbin + 3];
  for (int i = 0; i < nbin + 3; i++) data[i] = 0.0;
}

std::string densityplot::get_filename_out() {
  return basename + "_densityplot_" + data_to_bin_x + data_to_bin_y + ".txt";
}

aladyn_float densityplot::get_bin_size_x() {
  return (max_value_x - min_value_x) / nbin_x;
}

aladyn_float densityplot::get_bin_size_y() {
  return (max_value_y - min_value_y) / nbin_y;
}

void densityplot::write_binned_data() {
  std::ofstream file_out;
  file_out.open(get_filename_out());
  aladyn_float yminT = min_value_y - get_bin_size_y();
  aladyn_float ymaxT = min_value_y;
  aladyn_float xminT = min_value_x - get_bin_size_x();
  aladyn_float xmaxT = min_value_x;
  for (int i = 0; i < nbin_x; i++)
  {
    for (int j = 0; j < nbin_y; j++)
    {
      file_out << std::setprecision(6) << xminT << "\t" << xmaxT << "\t" << yminT << "\t" << ymaxT << "\t" << data[i][j] << std::endl;
      yminT += get_bin_size_y();
      ymaxT += get_bin_size_y();
    }
    xminT += get_bin_size_x();
    xmaxT += get_bin_size_x();
    yminT = min_value_y - get_bin_size_y();
    ymaxT = min_value_y;
  }
  file_out.close();
}

size_t densityplot::get_x_column_to_bin() {
  size_t bin_on_x = std::numeric_limits<size_t>::signaling_NaN();;
  if (data_to_bin_x == "x") bin_on_x = COLUMN_X;
  else if (data_to_bin_x == "y") bin_on_x = COLUMN_Y;
  else if (data_to_bin_x == "z") bin_on_x = COLUMN_Z;
  else if (data_to_bin_x == "px") bin_on_x = COLUMN_PX;
  else if (data_to_bin_x == "py") bin_on_x = COLUMN_PY;
  else if (data_to_bin_x == "pz") bin_on_x = COLUMN_PZ;
  else if (data_to_bin_x == "gamma") bin_on_x = COLUMN_GAMMA;
  else if (data_to_bin_x == "theta") bin_on_x = COLUMN_THETA;
  else if (data_to_bin_x == "E") bin_on_x = COLUMN_E;
  else if (data_to_bin_x == "thetat") bin_on_x = COLUMN_THETAT;
  else if (data_to_bin_x == "ty") bin_on_x = COLUMN_TY;
  else if (data_to_bin_x == "tz") bin_on_x = COLUMN_TZ;
  else if (data_to_bin_x == "w") bin_on_x = COLUMN_W;
  else if (data_to_bin_x == "ch") bin_on_x = COLUMN_CH;

  return bin_on_x;
}

size_t densityplot::get_y_column_to_bin() {
  size_t bin_on_y = std::numeric_limits<size_t>::signaling_NaN();;
  if (data_to_bin_y == "x") bin_on_y = COLUMN_X;
  else if (data_to_bin_y == "y") bin_on_y = COLUMN_Y;
  else if (data_to_bin_y == "z") bin_on_y = COLUMN_Z;
  else if (data_to_bin_y == "px") bin_on_y = COLUMN_PX;
  else if (data_to_bin_y == "py") bin_on_y = COLUMN_PY;
  else if (data_to_bin_y == "pz") bin_on_y = COLUMN_PZ;
  else if (data_to_bin_y == "gamma") bin_on_y = COLUMN_GAMMA;
  else if (data_to_bin_y == "theta") bin_on_y = COLUMN_THETA;
  else if (data_to_bin_y == "E") bin_on_y = COLUMN_E;
  else if (data_to_bin_y == "thetat") bin_on_y = COLUMN_THETAT;
  else if (data_to_bin_y == "ty") bin_on_y = COLUMN_TY;
  else if (data_to_bin_y == "tz") bin_on_y = COLUMN_TZ;
  else if (data_to_bin_y == "w") bin_on_y = COLUMN_W;
  else if (data_to_bin_y == "ch") bin_on_y = COLUMN_CH;

  return bin_on_y;
}

densityplot::densityplot() {}

densityplot::densityplot(std::string file_basename, std::string xaxis_name, std::string yaxis_name, size_t xbin_counter, size_t ybin_counter, aladyn_float minx, aladyn_float maxx, aladyn_float miny, aladyn_float maxy) {
  basename = file_basename;
  data_to_bin_x = xaxis_name;
  data_to_bin_y = yaxis_name;
  nbin_x = xbin_counter;
  nbin_y = ybin_counter;
  min_value_x = minx;
  min_value_y = miny;
  max_value_x = maxx;
  max_value_y = maxy;
  data = new aladyn_float*[nbin_x + 3];
  for (size_t i = 0; i < nbin_x + 3; i++)
  {
    data[i] = new aladyn_float[nbin_y + 3];
    for (size_t j = 0; j < nbin_y + 3; j++) data[i][j] = 0.0;
  }
}



Parameters::Parameters() {
  intpar.resize(20, 0);
  realpar.resize(20, 0.0);
  phasespace_file_labels.push_back("prp");
  phasespace_file_labels.push_back("hip");
  phasespace_file_labels.push_back("h1p");
  phasespace_file_labels.push_back("h2p");
  phasespace_file_labels.push_back("lip");
  phasespace_file_labels.push_back("elp");
  grid_file_labels.push_back("pren");
  grid_file_labels.push_back("pden");
  grid_file_labels.push_back("hidn");
  grid_file_labels.push_back("hien");
  grid_file_labels.push_back("h1dn");
  grid_file_labels.push_back("h1en");
  grid_file_labels.push_back("h2dn");
  grid_file_labels.push_back("h2en");
  grid_file_labels.push_back("lidn");
  grid_file_labels.push_back("lien");
  grid_file_labels.push_back("eden");
  grid_file_labels.push_back("elen");
  grid_file_labels.push_back("bden");
  grid_file_labels.push_back("ex");
  grid_file_labels.push_back("ey");
  grid_file_labels.push_back("ez");
  grid_file_labels.push_back("bx");
  grid_file_labels.push_back("by");
  grid_file_labels.push_back("bz");
  grid_file_labels.push_back("jx");
  grid_file_labels.push_back("jy");
  grid_file_labels.push_back("jz");
  header_size_bytes = 0;
  subsample = 1;
  span = 5;
  ncpu_x = ncpu_y = ncpu_z = ncpu = 0;
  nptot = ndv = 0;
  we_dont_know_file_version = true;
  we_dont_know_if_sim_is_2d = true;
  sim_is_2d = false;
  file_has_weight = false;
  we_dont_know_if_file_has_weight = true;
  file_has_charge = false;
  we_dont_know_if_file_has_charge = true;
  we_have_to_do_swap = false;
  we_dont_know_if_we_have_to_do_swap = true;
  out_ppg = out_json = out_csv = out_xyze = out_cutx = out_cuty = out_cutz = out_grid2d = out_clean_bin = out_lineoutx = out_vtk = out_vtk_nostretch = false;
  npx = npy = npz = npx_per_cpu = npy_per_cpu = npz_per_cpu = 0;
  npx_resampled = npy_resampled = npz_resampled = npx_resampled_per_cpu = npy_resampled_per_cpu = npz_resampled_per_cpu = 0;
  resampling_factor = 0;
  endianness = 0;
  file_version = 0; // initialized to an invalid file version
  fixed_aladyn_version = false;
  multifile = false;
  stretched_grid = true;
  stretched_along_x = 1;
  mass_MeV = 0.;
  tnow = 0.0;
  xmin = ymin = zmin = 0.0;
  xmax = ymax = zmax = 1.0;
  overwrite_weight = false;
  overwrite_charge = false;
  overwrite_weight_value = 1.0;
  overwrite_charge_value = 1.0;
  endian_file = 0;
  endian_machine = is_big_endian();
  phasespace_file = false;
  grid_file = false;
}


void Parameters::man(const char argv[]) {
  std::cerr << argv[0] << " filebasename -arguments ..." << std::endl;
  std::cerr << "Example: " << argv[0] << " Edenout01 -densityplot E 0.0 30.0 60 thetaT 0 0.20 20" << std::endl;
  std::cerr << "----------Full argument list------------------- " << std::endl;
  std::cerr << "-swap/-noswap (force endianess swap/noswap) -force_v1 -force_v2 -force_v3 -force_v4 (force specific file format)" << std::endl;
  std::cerr << "-dump_ppg -dump_csv -dump_clean -dump_xyzE -find_minxmax" << std::endl;
  std::cerr << "-dump_cutx #x -dump_cuty #y -dump_cutz #z  -dump_lineoutx -dump_gnuplot" << std::endl;
  std::cerr << "-dump_vtk -dump_vtk_nostretch (dumps in the vtk just the unstretched part of the grid)" << std::endl;
  std::cerr << "(use -no_stretch_x if the grid is not stretched along x axis)" << std::endl;
  std::cerr << "-[x,y,z,px,py,pz,theta,thetaT,gamma,E,ty,tz,w,ch]min/max #number (to force min/max values)" << std::endl;
  std::cerr << "-densityplot \"A\" minA maxA nbinA \"B\" minB maxB nbinB [A,B={x,y,z,px,py,pz,theta,thetaT,gamma,E,ty,tz,w,ch}] (to obtain a density plot of the two variables)" << std::endl;
  std::cerr << "-histogram \"A\" minA maxA nbinA [A={x,y,z,px,py,pz,theta,thetaT,gamma,E,ty,tz,w,ch}] (to obtain a histogram of the variable)" << std::endl;
  std::cerr << "Filters: \n +[x,y,z,px,py,pz,theta,thetaT,gamma,E,ty,tz,w,ch]min/max #num (to completely exclude particles not inside the filter)" << std::endl;
}



void Parameters::read_params_from_bin_file(const char * filename)
{
  std::FILE * file_in = nullptr;
  int fortran_buff;
  size_t fread_size = 0;
  int nparams;

  file_in = fopen(filename, "rb");
  if (file_in == nullptr) std::cerr << "Unable to open file!" << std::endl;
  else std::cout << "File opened to read parameters!" << std::endl;

  fread_size += std::fread(&fortran_buff, sizeof(int), 1, file_in);
  fread_size += std::fread(&nparams, sizeof(int), 1, file_in);
  fread_size += std::fread(&fortran_buff, sizeof(int), 1, file_in);

  if (we_have_to_do_swap) swap_endian_i(&nparams, 1);
  if (nparams > intpar.size() || nparams > realpar.size()) std::cerr << "Bad number of parameters found in bin file, exiting" << std::endl, exit(-1);

  fread_size += std::fread(&fortran_buff, sizeof(int), 1, file_in);
  fread_size += std::fread(&intpar[0], sizeof(int), nparams, file_in);
  fread_size += std::fread(&fortran_buff, sizeof(int), 1, file_in);
  fread_size += std::fread(&fortran_buff, sizeof(int), 1, file_in);
  fread_size += std::fread(&realpar[0], sizeof(aladyn_float), nparams, file_in);
  fread_size += std::fread(&fortran_buff, sizeof(int), 1, file_in);

  if (we_have_to_do_swap) swap_endian_i(intpar);
  if (we_have_to_do_swap) swap_endian_f(realpar);

  fclose(file_in);

  /* overwrite default with good values */
  ncpu_x = 1;
  ncpu_y = intpar[0];
  ncpu_z = intpar[1];
  npx_resampled_per_cpu = intpar[2];
  npx = intpar[3];
  npy = intpar[4];
  npy_resampled_per_cpu = intpar[5];
  npz = intpar[6];
  npz_resampled_per_cpu = intpar[7];
  nptot = (long long int) intpar[16];
  ndv = intpar[17];
  file_version = intpar[18];
  endianness = intpar[19];
  tnow = realpar[0];
  xmin = realpar[1];
  xmax = realpar[2];
  ymin = realpar[3];
  ymax = realpar[4];
  zmin = realpar[5];
  zmax = realpar[6];

  sim_is_2d = ((ndv == 4 || ndv == 5) && file_version < 3) || (ndv == 6 && file_version >= 3);

  header_size_bytes = (7 + nparams) * sizeof(int) + nparams * sizeof(aladyn_float);
  if (fread_size != header_size_bytes) std::cout << "error: header size is different than expected" << std::endl << std::flush;
}


void Parameters::check_swap() {
  if (endian_file == endian_machine) we_have_to_do_swap = false;
  else we_have_to_do_swap = true;
  we_dont_know_if_we_have_to_do_swap = false;
}


void Parameters::read_params_from_dat_file(std::ifstream& file_dat)
{
  std::string forget_this_line;
  int resampling_factor;
  aladyn_float coord;
  std::getline(file_dat, forget_this_line); // per leggere la riga Integer parameters


  for (size_t i = 0; i < intpar.size(); i++)
  {
    file_dat >> intpar[i];
    if (file_dat.fail())
    {
      file_dat.clear();
      std::cout << "Unable to parse int_par #" << i + 1 << std::endl;
      if (i <= 7 || i >= 16)
      {
        std::cout << "Bad error; please fix the .dat file if possibile and then re-run the program" << std::endl;
        exit(-77);
      }
      file_dat.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
    }
  }
  std::getline(file_dat, forget_this_line); // per pulire i caratteri rimanenti sull'ultima riga degli interi
  std::getline(file_dat, forget_this_line); // per leggere la riga Real parameters
  for (size_t i = 0; i < realpar.size(); i++)
  {
    file_dat >> realpar[i];
    if (file_dat.fail())
    {
      file_dat.clear();
      std::cout << "Unable to parse real_par #" << i + 1 << std::endl;
      if (i <= 6)
      {
        std::cout << "Bad error; please fix the .dat file if possibile and then re-run the program" << std::endl;
        exit(-77);
      }
      file_dat.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
    }
  }

  if (!fixed_aladyn_version) file_version = intpar[18];
  std::cout << "Parsing an ALaDyn file versioned as v" << file_version << std::endl;

  if (file_version == 1) {

    /*
    real_par(1:20) =(/tnow,xmin,xmax,ymin,ymax,zmin,zmax,w0_x,w0_y,&
    n_over_nc,a0,lam0,E0,ompe,targ_in,targ_end,&
    gam0,nb_over_np,b_charge,vbeam/)

    int_par(1:20) = (/npe_loc,npe_zloc,nx,nxh,ny1,loc_nyc_max,nz1,loc_nzc_max,jump,iby,iform,&
    model_id,dmodel_id,nsp,curr_ndim,np_per_cell(1),&
    LPf_ord,der_ord,file_version,i_end/)
    */

    ncpu_y = intpar[0];
    ncpu_z = intpar[1];
    ncpu_x = 1;
    ncpu = ncpu_x * ncpu_y * ncpu_z;

    endianness = intpar[19];
    tnow = realpar[0];
    xmin = realpar[1];
    xmax = realpar[2];
    ymin = realpar[3];
    ymax = realpar[4];
    zmin = realpar[5];
    zmax = realpar[6];

    if (grid_file) {
      npx = intpar[2];
      npx_resampled = intpar[3];
      npx_per_cpu = npx / ncpu_x;
      resampling_factor = (int)(npx / npx_resampled);
      npy_resampled = intpar[4];
      npy_per_cpu = intpar[5];
      npy = npy_resampled * resampling_factor;
      npz_resampled = intpar[6];
      npz_per_cpu = intpar[7];
      if (npz_resampled > 1) npz = npz_resampled * resampling_factor;
      else npz = npz_resampled;
      npx_resampled_per_cpu = npx_resampled / ncpu_x;
      npy_resampled_per_cpu = npy_resampled / ncpu_y;
      npz_resampled_per_cpu = npz_resampled / ncpu_z;
      header_size_bytes = (intpar.size() + 1) * sizeof(int) + realpar.size() * sizeof(aladyn_float) + 6 * sizeof(int); // there are 6 fortran buffers: two around nparams, two around intpars and two around realpars
    }
    else {
      nptot = (long long int) intpar[16];
      ndv = intpar[17];
      header_size_bytes = (intpar.size() + 1) * sizeof(int) + realpar.size() * sizeof(aladyn_float);
    }
  }
  else if (file_version == 2) {
    endianness = intpar[19];

    tnow = realpar[0];
    xmin = realpar[1];
    xmax = realpar[2];
    ymin = realpar[3];
    ymax = realpar[4];
    zmin = realpar[5];
    zmax = realpar[6];

    if (grid_file) {

      /*
      real_par(1:20) = (/ tnow, xmin, xmax, ymin, ymax, zmin, zmax, w0_x, w0_y, &
      n_over_nc, a0, lam0, E0, ompe, targ_in, targ_end, &
      gam0, nb_over_np, b_charge, vbeam / )

      int_par(1:20) = (/ npe_yloc, npe_zloc, npe_xloc, &
      nx1, ny1, loc_nyc_max, nz1, loc_nzc_max, jump, iby, iform, &
      model_id, dmodel_id, nsp, curr_ndim, mp_per_cell(1), &
      LPf_ord, der_ord, file_version, i_end / )
      */

      ncpu_y = intpar[0];
      ncpu_z = intpar[1];
      ncpu_x = intpar[2];
      ncpu = ncpu_x * ncpu_y * ncpu_z;

      npx_resampled = intpar[3];
      npy_resampled = intpar[4];
      npy_per_cpu = intpar[5];
      npz_resampled = intpar[6];
      npz_per_cpu = intpar[7];
      resampling_factor = intpar[8];
      npx_resampled_per_cpu = npx_resampled / ncpu_x;
      npy_resampled_per_cpu = npy_per_cpu / resampling_factor;
      npz_resampled_per_cpu = npz_per_cpu / resampling_factor;
      npx = npx_resampled * resampling_factor;
      npy = npy_resampled * resampling_factor;
      if (npz_resampled > 1) npz = npz_resampled * resampling_factor;
      else npz = npz_resampled;
      npx_per_cpu = npx / ncpu_x;
      header_size_bytes = (intpar.size() + 1 + 6) * sizeof(int) + realpar.size() * sizeof(aladyn_float); // +1 for n_par, +6 for fortran buffers (2 around nparams, 2 around intpars, 2 around realpars)
    }

    else {
      nptot = (long long int) intpar[16];
      ndv = intpar[17];
      header_size_bytes = (intpar.size() + 1) * sizeof(int) + realpar.size() * sizeof(aladyn_float);
    }

  }
  else if (file_version == 3) {
    tnow = realpar[0];
    xmin = realpar[1];
    xmax = realpar[2];
    ymin = realpar[3];
    ymax = realpar[4];
    zmin = realpar[5];
    zmax = realpar[6];

    header_size_bytes = (intpar.size() + 1) * sizeof(int) + realpar.size() * sizeof(aladyn_float);
    endianness = intpar[19];

    if (grid_file) {
      /*
      real_par(1:20) =(/tnow,xmin,xmax,ymin,ymax,zmin,zmax,w0_x,w0_y,&
      n_over_nc,a0,lam0,E0,ompe,targ_in,targ_end,&
      gam0,nb_over_np,b_charge,vbeam/)

      int_par(1:20) = (/npe_loc,npe_zloc,npe_xloc,&
      nx1,ny1,loc_nyc_max,nz1,loc_nzc_max,jump,iby,iform,&
      model_id,dmodel_id,nsp,curr_ndim,mp_per_cell(1),&
      LPf_ord,der_ord,file_version,i_end/)
      */

      ncpu_y = intpar[0];
      ncpu_z = intpar[1];
      ncpu_x = intpar[2];
      ncpu = ncpu_x * ncpu_y * ncpu_z;

      npx_resampled = intpar[3];
      npx_resampled_per_cpu = npx_resampled / ncpu_x;
      npy_resampled = intpar[4];
      npy_resampled_per_cpu = intpar[5];
      npz_resampled = intpar[6];
      npz_resampled_per_cpu = intpar[7];
      resampling_factor = intpar[8];
      npx = npx_resampled * resampling_factor;
      npy = npy_resampled * resampling_factor;
      if (npz_resampled > 1) npz = npz_resampled * resampling_factor;
      else npz = npz_resampled;
      npx_per_cpu = npx / ncpu_x;
      npy_per_cpu = npy / ncpu_y;
      npz_per_cpu = npz / ncpu_z;
    }
    else {
      nptot = (long long int) intpar[16];
      ndv = intpar[17];
    }
  }
  else if (file_version == 4) {
    tnow = realpar[0];
    xmin = realpar[1];
    xmax = realpar[2];
    ymin = realpar[3];
    ymax = realpar[4];
    zmin = realpar[5];
    zmax = realpar[6];

    header_size_bytes = 0;
    /*i_end*/    endianness = intpar[19];

    if (grid_file) {
      /*
      real_par(1:20) =(/tnow,xmin,xmax,ymin,ymax,zmin,zmax,w0_x,w0_y,&
      n_over_nc,a0,lam0,E0,ompe,targ_in,targ_end,&
      gam0,nb_over_np,b_charge,vbeam/)

      int_par(1:20) = (/ npe_loc, npe_zloc, npe_xloc, nx1, ny1, nz1, &
      jump, ibx, iby, iform, pid, &
      model_id, dmodel_id, nsp, curr_ndim, mp_per_cell(1), &
      nptot, ndv, file_version, i_end / )
      */

      /*npe_xloc*/ ncpu_x = intpar[2];
      /*npe_loc*/  ncpu_y = intpar[0];
      /*npe_zloc*/ ncpu_z = intpar[1];
      ncpu = ncpu_x * ncpu_y * ncpu_z;
      /*nx1*/      npx_per_cpu = intpar[3];
      /*jump*/     resampling_factor = intpar[6];
      /**/         npx_resampled_per_cpu = npx_per_cpu / resampling_factor;
      /*ny1*/      npy_per_cpu = intpar[4];
      /**/         npy_resampled_per_cpu = npy_per_cpu / resampling_factor;
      /*nz1*/      npz_per_cpu = intpar[5];
      /**/         npz_resampled_per_cpu = npz_per_cpu / resampling_factor;
      /**/         npx = npx_per_cpu * ncpu_x;
      /**/         npy = npy_per_cpu * ncpu_y;
      /**/         npz = npz_per_cpu * ncpu_z;
      /**/         npx_resampled = npx_resampled_per_cpu * ncpu_x;
      /**/         npy_resampled = npy_resampled_per_cpu * ncpu_y;
      /**/         npz_resampled = npz_resampled_per_cpu * ncpu_z;
    }
    else {
      /*
      int_par(1:20) = (/ npe, nx, ny_loc, nz_loc, jmp, iby, iform, &
        model_id, dmodel_id, nsp, curr_ndim, mp_per_cell(1), &
        LPf_ord, der_ord, iform, pid, nptot, ndv, file_version, i_end / )
      */
      ncpu = intpar[0];
      nptot = (long long int) intpar[16];
      ndv = intpar[17];
    }
  }


#ifdef ENABLE_DEBUG
  debug_read_parameters();
#endif


  if (phasespace_file) {
    if (file_version < 3)
    {
      if (ndv == 4 || ndv == 6) file_has_weight = 0;
      else if (ndv == 5 || ndv == 7) file_has_weight = 1;
      else std::cerr << "Attention: illegal value for ndv" << std::endl, exit(-17);
      if (ndv == 4 || ndv == 5) zmin = 0.0, zmax = 1.0;
      we_dont_know_if_file_has_weight = false;
    }
    else
    {
      file_has_weight = 1;
      we_dont_know_if_file_has_weight = false;
      if (ndv < 7) zmin = 0.0, zmax = 1.0; // a 2D file has 6 floats (columns): x, y, px, py, w, ch
      if (nptot == -1)
      {
        std::getline(file_dat, forget_this_line); // per pulire i caratteri rimanenti sull'ultima riga dei real
        std::getline(file_dat, forget_this_line); // per leggere la riga Number of particles
        file_dat >> nptot;
      }
    }
    sim_is_2d = ((ndv == 4 || ndv == 5) && file_version < 3) || (ndv == 6 && file_version >= 3);
    we_dont_know_if_sim_is_2d = false;
  }
  else {
    if (npz == 1) zmin = 0.0, zmax = 1.0;
    if (file_version > 2) {
      std::getline(file_dat, forget_this_line); // per pulire i caratteri rimanenti sull'ultima riga dei aladyn_float
      std::getline(file_dat, forget_this_line); // per togliere la riga vuota che separa la griglia dai params


      for (unsigned int i = 0; i < npx_resampled; i++)
      {
        file_dat >> coord;
        xcoord.push_back(coord);
      }
      for (unsigned int i = 0; i < npy_resampled; i++)
      {
        file_dat >> coord;
        ycoord.push_back(coord);
      }
      for (unsigned int i = 0; i < npz_resampled; i++)
      {
        file_dat >> coord;
        zcoord.push_back(coord);
      }
    }
    else {     // mettiamo una griglia temporanea fissa, che al limite sara' sovrascritta da quella stretchata se presente nel binario
      aladyn_float dx, dy, dz;
      if (npx_resampled > 1) dx = (xmax - xmin) / (npx_resampled - 1);
      else dx = (xmax - xmin);
      if (npy_resampled > 1) dy = (ymax - ymin) / (npy_resampled - 1);
      else dy = (ymax - ymin);
      if (npz_resampled > 1) dz = (zmax - zmin) / (npz_resampled - 1);
      else dz = (zmax - zmin);

      for (unsigned int i = 0; i < npx_resampled; i++)
      {
        coord = xmin + dx*i;
        xcoord.push_back(coord);
      }
      for (unsigned int i = 0; i < npy_resampled; i++)
      {
        coord = ymin + dy*i;
        ycoord.push_back(coord);
      }
      for (unsigned int i = 0; i < npz_resampled; i++)
      {
        coord = zmin + dz*i;
        zcoord.push_back(coord);
      }
    }
    sim_is_2d = (npz == 1);
    we_dont_know_if_sim_is_2d = false;
  }
  endian_file = (endianness - 1);
}




void Parameters::debug_read_parameters()
{
  std::cout << "Integer parameters" << std::endl;
  for (size_t i = 0; i < intpar.size(); i++)
  {
    std::cout << std::setw(14) << intpar[i];
    if (i > 0 && !((i + 1) % 4)) std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Real parameters" << std::endl;
  for (size_t i = 0; i < realpar.size(); i++)
  {
    std::cout << std::setw(14) << realpar[i];
    if (i > 0 && !((i + 1) % 4)) std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "ncpu_x = " << ncpu_x << std::endl;
  std::cout << "ncpu_y = " << ncpu_y << std::endl;
  std::cout << "ncpu_z = " << ncpu_z << std::endl;
  std::cout << "ncpu = " << ncpu << std::endl;
  std::cout << "npx = " << npx << std::endl;
  std::cout << "npy = " << npy << std::endl;
  std::cout << "npz = " << npz << std::endl;
  std::cout << "npx_resampled = " << npx_resampled << std::endl;
  std::cout << "npy_resampled = " << npy_resampled << std::endl;
  std::cout << "npz_resampled = " << npz_resampled << std::endl;
  std::cout << "resampling_factor = " << resampling_factor << std::endl;
  std::cout << "npx_per_cpu = " << npx_per_cpu << std::endl;
  std::cout << "npy_per_cpu = " << npy_per_cpu << std::endl;
  std::cout << "npz_per_cpu = " << npz_per_cpu << std::endl;
  std::cout << "npx_resampled_per_cpu = " << npx_resampled_per_cpu << std::endl;
  std::cout << "npy_resampled_per_cpu = " << npy_resampled_per_cpu << std::endl;
  std::cout << "npz_resampled_per_cpu = " << npz_resampled_per_cpu << std::endl;
  std::cout << "nptot = " << nptot << std::endl;
  std::cout << "ndv = " << ndv << std::endl;
  std::cout << "endianness = " << endianness << std::endl;
}


void Parameters::check_filename(const char *filename_in) {
  for (auto local_support_label : phasespace_file_labels) {
    std::string filename = filename_in;
    std::transform(filename.begin(), filename.end(), filename.begin(), ::tolower);
    if (boost::starts_with(filename, local_support_label)) {
      mass_MeV = (aladyn_float)MP_MEV; // fix_wrong!!!!
      phasespace_file = true;
      return;
    }

  }
  for (auto local_support_label : grid_file_labels) {
    std::string filename = filename_in;
    std::transform(filename.begin(), filename.end(), filename.begin(), ::tolower);
    if (boost::starts_with(filename, local_support_label)) {
      grid_file = true;
      return;
    }

  }
  std::cerr << "WARNING: unable to understand file type: " << filename_in << std::endl;
  exit(-10);

}


void Parameters::check_forced_version(const int cl_argc, const char ** cl_argv)
{
  argc = cl_argc;
  argv = new std::string[argc];
  for (int i = 0; i < argc; i++) argv[i] = std::string(cl_argv[i]);

  for (int i = 2; i < argc; i++) {
    if (argv[i] == "-force_v1") {
      std::cout << "Forced using routines for file v1" << std::endl;
      file_version = 1;
      fixed_aladyn_version = true;
    }
    else if (argv[i] == "-force_v2") {
      std::cout << "Forced using routines for file v2" << std::endl;
      file_version = 2;
      fixed_aladyn_version = true;
    }
    else if (argv[i] == "-force_v3") {
      std::cout << "Forced using routines for file v3" << std::endl;
      file_version = 3;
      fixed_aladyn_version = true;
    }
    else if (argv[i] == "-force_v4") {
      std::cout << "Forced using routines for file v4" << std::endl;
      file_version = 4;
      fixed_aladyn_version = true;
    }
  }
}




void Parameters::parse_command_line()
{
  for (int i = 2; i < argc; i++)
    /************************************************************************
    We will iterate over argv[] to get the parameters stored inside.
    Note that we're starting on 1 because we don't need to know the
    path of the program, which is stored in argv[0], and the input file,
    which is supposed to be given as the first argument and so is in argv[1]
    ************************************************************************/
  {
    if (argv[i] == "-swap")
    {
      std::cout << "Forcing byte swapping (endianness)" << std::endl;
      we_have_to_do_swap = true;
      we_dont_know_if_we_have_to_do_swap = false;
    }
    else if (argv[i] == "-noswap")
    {
      std::cout << "Forcing byte NON swapping (endianness)" << std::endl;
      we_have_to_do_swap = 0;
      we_dont_know_if_we_have_to_do_swap = false;
    }
    else if (argv[i] == "-span")
    {
      span = atoi(argv[i + 1].c_str());
      if (span < 0) span = 0;
      std::cout << "Span factor for lineout: " << span << std::endl;
      i++;
    }
    else if (argv[i] == "-subsample")
    {
      subsample = atoi(argv[i + 1].c_str());
      if (subsample < 1)
      {
        subsample = 1;
        std::cout << "Value for subsampling not valid, disabled" << std::endl;
      }
      else
      {
        std::cout << "Will subsample with a ratio of 1:" << subsample << " if any ASCII dump will be requested (valid only for phase-space files)" << std::endl;
      }
      i++;
    }
    else if (argv[i] == "-ncol3d")
    {
      sim_is_2d = false;
      int ncolumns = atoi(argv[i + 1].c_str());
      if (ncolumns == 6) file_has_weight = false, file_has_charge = false;
      else if (ncolumns == 7) file_has_weight = true, file_has_charge = false;
      else if (ncolumns == 8) file_has_weight = true, file_has_charge = true;
      std::cout << "Forced number of columns in 3D binary file to " << ncolumns << std::endl;
      we_dont_know_if_file_has_weight = false;
      we_dont_know_if_file_has_charge = false;
      we_dont_know_if_sim_is_2d = false;
      i++;
    }
    else if (argv[i] == "-ncol2d")
    {
      sim_is_2d = true;
      int ncolumns = atoi(argv[i + 1].c_str());
      if (ncolumns == 4) file_has_weight = false, file_has_charge = false;
      else if (ncolumns == 5) file_has_weight = true, file_has_charge = false;
      else if (ncolumns == 6) file_has_weight = true, file_has_charge = true;
      std::cout << "Forced number of columns in 2D binary file to " << ncolumns << std::endl;
      we_dont_know_if_file_has_weight = false;
      we_dont_know_if_file_has_charge = false;
      we_dont_know_if_sim_is_2d = false;
      i++;
    }
    else if (argv[i] == "-dump_vtk")
    {
      std::cout << "You asked to have a VTK dump of the input file" << std::endl;
      out_vtk = true;
    }
    else if (argv[i] == "-dump_vtk_nostretch")
    {
      std::cout << "You asked to have a VTK dump of the non-stretched grid.\n";
      std::cout << "If not explicitly said, the grid will be considered stretched in ALL directions" << std::endl;
      out_vtk_nostretch = true;
    }
    else if (argv[i] == "-no_stretch_x")
    {
      std::cout << "Assuming the grid is NOT stretched along x axis.\n";
      stretched_along_x = 0;
    }
    else if (argv[i] == "-dump_cutx")
    {
      if (argv[i + 1][0] != '-')
      {
        aladyn_float posizione_taglio = (aladyn_float)atof(argv[i + 1].c_str());
        where_to_cut_grid_along_x.push_back(posizione_taglio);
        std::cout << "You asked to cut the grid at x = " << posizione_taglio << std::endl;
        i++;
      }
      else
      {
        std::cout << "You asked to cut the grid at the middle of the x-axis" << std::endl;
      }
      out_cutx = 1;
    }

    else if (argv[i] == "-dump_cuty")
    {
      if (argv[i + 1][0] != '-')
      {
        aladyn_float posizione_taglio = (aladyn_float)atof(argv[i + 1].c_str());
        where_to_cut_grid_along_y.push_back(posizione_taglio);
        std::cout << "You asked to cut the grid at y = " << posizione_taglio << std::endl;
        i++;
      }
      else
      {
        std::cout << "You asked to cut the grid at the middle of the y-axis" << std::endl;
      }
      out_cuty = 1;
    }

    else if (argv[i] == "-dump_cutz")
    {
      if (argv[i + 1][0] != '-')
      {
        aladyn_float posizione_taglio = (aladyn_float)atof(argv[i + 1].c_str());
        where_to_cut_grid_along_z.push_back(posizione_taglio);
        std::cout << "You asked to cut the grid at z = " << posizione_taglio << std::endl;
        i++;
      }
      else
      {
        std::cout << "You asked to cut the grid at the middle of the z-axis" << std::endl;
      }
      out_cutz = true;
    }

    else if (argv[i] == "-dump_lineoutx")
    {
      std::cout << "You asked to have a lineout of the grid along the x-axis" << std::endl;
      out_lineoutx = 1;
    }

    else if (argv[i] == "-dump_gnuplot")
    {
      std::cout << "You asked to rewrite the 2D grid in ASCII format for gnuplot" << std::endl;
      out_grid2d = 1;
    }
    else if (argv[i] == "-dump_propaga")
    {
      std::cout << "You asked to have a .ppg dump of the input phase space" << std::endl;
      out_ppg = 1;
    }
    else if (argv[i] == "-dump_csv")
    {
      std::cout << "You asked to have a .csv dump of the input phase space" << std::endl;
      out_csv = 1;
    }
    else if (argv[i] == "-dump_xyzE")
    {
      std::cout << "You asked to have a xy(z)E dump of the input phase space" << std::endl;
      out_xyze = 1;
    }
    else if (argv[i] == "-dump_clean")
    {
      std::cout << "You asked to have a unique, clean binary file as the output" << std::endl;
      out_clean_bin = 1;
    }
    else if (argv[i] == "-dump_json")
    {
      std::cout << "You asked to write the simulation parameters file" << std::endl;
      out_json = 1;
    }
    else if (argv[i] == "-weight")
    {
      overwrite_weight = true;
      overwrite_weight_value = (aladyn_float)atof(argv[i + 1].c_str());
      i++;
    }
    else if (argv[i] == "-charge")
    {
      overwrite_charge = true;
      overwrite_charge_value = (aladyn_float)atof(argv[i + 1].c_str());
      i++;
    }
    else if (argv[i] == "-densityplot")
    {
      std::string x_axis = argv[++i];
      aladyn_float x_min = (aladyn_float)atof(argv[++i].c_str());
      aladyn_float x_max = (aladyn_float)atof(argv[++i].c_str());
      size_t binx = (size_t)atoll(argv[++i].c_str());
      std::string y_axis = argv[++i];
      aladyn_float y_min = (aladyn_float)atof(argv[++i].c_str());
      aladyn_float y_max = (aladyn_float)atof(argv[++i].c_str());
      size_t biny = (size_t)atoll(argv[++i].c_str());
      densityplots.push_back(densityplot(filebasename, x_axis, y_axis, binx, biny, x_min, x_max, y_min, y_max));
    }
    else if (argv[i] == "-histogram")
    {
      std::string x_axis = argv[++i];
      aladyn_float x_min = (aladyn_float)atof(argv[++i].c_str());
      aladyn_float x_max = (aladyn_float)atof(argv[++i].c_str());
      size_t binx = (size_t)atoll(argv[++i].c_str());
      histograms.push_back(histo(filebasename, x_axis, binx, x_min, x_max));
    }
  }
}


void Parameters::parse_json() {
  std::string json_filename = filebasename + ".json";

  jsoncons::json jsonpar;
  try {
    jsonpar = jsoncons::json::parse_file(json_filename);
  }
  catch (std::exception &e) {
    std::cerr << "Exception found: " << e.what() << std::endl;
  }
  //TODO
}


