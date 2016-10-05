
#include "leggi_binario_ALaDyn_fortran.h"
#include "swap_tools.h"


Parameters::Parameters()
{
  nparams = NUMBER_OF_PARAMS_IN_DAT_FILE;
  intpar.resize(NUMBER_OF_PARAMS_IN_DAT_FILE, 0);
  realpar.resize(NUMBER_OF_PARAMS_IN_DAT_FILE, 0.0);
  header_size_bytes = 0;
  subsample = 1;
  span = 5;
  ncpu_x = ncpu_y = ncpu_z = ncpu = 0;
  nptot = ndv = 0;
  we_dont_know_file_version = true;
  we_dont_know_if_sim_is_2d = true;
  sim_is_2d = false;
  we_have_to_find_minmax = false;
  we_dont_know_if_we_have_to_find_minmax = true;
  we_have_to_do_binning = false;
  we_dont_know_if_we_have_to_do_binning = true;
  file_has_weight = false;
  we_dont_know_if_file_has_weight = true;
  file_has_charge = false;
  we_dont_know_if_file_has_charge = true;
  we_have_to_do_swap = false;
  we_dont_know_if_we_have_to_do_swap = true;
  out_vtk = out_ppg = out_params = out_csv = out_xyze = out_cutx = out_cuty = out_cutz = out_grid2d = out_clean_bin = out_lineoutx = out_vtk_nostretch = false;
  minmax_found = false;
  npx = npy = npz = npx_per_cpu = npy_per_cpu = npz_per_cpu = 0;
  npx_resampled = npy_resampled = npz_resampled = npx_resampled_per_cpu = npy_resampled_per_cpu = npz_resampled_per_cpu = 0;
  resampling_factor = 0;
  endianness = 0;
  file_version = -10;
  fixed_aladyn_version = false;
  multifile = false;
  stretched_grid = true;
  stretched_along_x = 1;
  mass_MeV = 0.;
  nbin = nbin_x = nbin_px = nbin_y = nbin_py = nbin_z = nbin_pz = nbin_w = nbin_ch = nbin_E = nbin_gamma = nbin_theta = nbin_thetaT = nbin_ty = nbin_tz = NUMBER_OF_BIN_BY_DEFAULT;
  tnow = 0.0;
  xmin = pxmin = ymin = pymin = zmin = pzmin = wmin = chmin = thetamin = thetaTmin = Emin = gammamin = 0.0;
  xmax = pxmax = ymax = pymax = zmax = pzmax = wmax = chmax = thetamax = thetaTmax = Emax = gammamax = 1.0;
  tymin = tzmin = -1.0;
  tymax = tzmax = 1.0;
  ymin_b = ymax_b = pymin_b = pymax_b = zmin_b = zmax_b = pzmin_b = pzmax_b = wmin_b = wmax_b = chmin_b = chmax_b = gammamin_b = gammamax_b = true;
  xmin_b = xmax_b = pxmin_b = pxmax_b = Emin_b = Emax_b = thetaTmin_b = thetaTmax_b = thetamin_b = thetamax_b = tymin_b = tymax_b = tzmin_b = tzmax_b = true;
  nbin_b = true;
  nbin_E_b = nbin_theta_b = nbin_thetaT_b = nbin_ty_b = nbin_tz_b = nbin_gamma_b = true;
  nbin_x_b = nbin_px_b = nbin_y_b = nbin_py_b = nbin_z_b = nbin_pz_b = true;
  do_plot_wspec = do_plot_chspec = do_plot_Espec = do_plot_thetaspec = do_plot_thetaTspec = do_plot_Etheta = do_plot_EthetaT = false;
  do_plot_xy = do_plot_xz = do_plot_yz = do_plot_xpx = do_plot_xpy = do_plot_xpz = do_plot_ypx = false;
  do_plot_ypy = do_plot_ypz = do_plot_zpx = do_plot_zpy = do_plot_zpz = do_plot_pxpy = do_plot_pxpz = do_plot_pypz = do_plot_xw = do_plot_rcf = false;
  overwrite_weight = false;
  overwrite_charge = false;
  overwrite_weight_value = 1.0;
  overwrite_charge_value = 1.0;
  do_not_ask_missing = false;
  memset(&support_label[0], 0, sizeof(support_label));
  last_cpu = MAX_NUMBER_OF_CPUS;    // il tool funziona quindi per un ncpu_max, attualmente, pari a 32768
  endian_file = 0;
  endian_machine = is_big_endian();
  phasespace_file = false;
  grid_file = false;
}


void Parameters::man(const char argv[]) {
  std::cout << "Interactive: " << argv << std::endl;
  std::cout << "Batch:       " << argv << " filebasename -arguments" << std::endl;
  std::cout << "----------Argument list------------------- " << std::endl;
  std::cout << "-params (write a .parameters file with params from .bin/.dat files)" << std::endl;
  std::cout << "-swap/-noswap (force endianess swap) -force_v1 -force_v2 -force_v3 -force_v4 (force specific file format)" << std::endl;
  std::cout << "-dump_vtk -dump_cutx #x -dump_cuty #y -dump_cutz #z  -dump_lineoutx -dump_gnuplot" << std::endl;
  std::cout << "-dump_vtk_nostretch (dumps in the vtk just the unstretched part of the grid)" << std::endl;
  std::cout << "(use -no_stretch_x if the grid is not stretched along x axis)" << std::endl;
  std::cout << "-dump_ppg -dump_csv -dump_clean -dump_xyzE -parameters -find_minxmax" << std::endl;
  std::cout << "-[x,y,z,px,py,pz,theta,thetaT,gamma,E,ty,tz,w,ch]min/max #number" << std::endl;
  std::cout << "-plot_AB A,B={x,y,z,px,py,pz}" << std::endl;
  std::cout << "-plot_etheta -plot_ethetaT -plot_rfc -plot_espec -plot_thetaspec -plot_thetaTspec -plot_chspec" << std::endl;
  std::cout << "-nbin #num" << std::endl;
  std::cout << "-nbin[x,y,z,px,py,pz,theta,thetaT,gamma,E,ty,tz,w,ch] #num" << std::endl;
  std::cout << "-dontask [TRY TO RUN IN NON-INTERACTIVE MODE]" << std::endl;
  std::cout << "Filters: \n +[x,y,z,px,py,pz,theta,thetaT,gamma,E,ty,tz,w,ch]min/max #num" << std::endl;
  std::cout << "----------Argument list------------------- " << std::endl;
}


aladyn_float Parameters::get_size_x()
{
  return (xmax - xmin) / static_cast <aladyn_float> (nbin_x);
}
aladyn_float Parameters::get_size_y()
{
  return (ymax - ymin) / static_cast <aladyn_float> (nbin_y);
}
aladyn_float Parameters::get_size_z()
{
  return (zmax - zmin) / static_cast <aladyn_float> (nbin_z);
}
aladyn_float Parameters::get_size_ty()
{
  return (tymax - tymin) / static_cast <aladyn_float> (nbin_ty);
}
aladyn_float Parameters::get_size_tz()
{
  return (tzmax - tzmin) / static_cast <aladyn_float> (nbin_tz);
}
aladyn_float Parameters::get_size_px()
{
  return (pxmax - pxmin) / static_cast <aladyn_float> (nbin_px);
}
aladyn_float Parameters::get_size_py()
{
  return (pymax - pymin) / static_cast <aladyn_float> (nbin_py);
}
aladyn_float Parameters::get_size_pz()
{
  return (pzmax - pzmin) / static_cast <aladyn_float> (nbin_pz);
}
aladyn_float Parameters::get_size_w()
{
  return (wmax - wmin) / static_cast <aladyn_float> (nbin_w);
}
aladyn_float Parameters::get_size_ch()
{
  return (chmax - chmin) / static_cast <aladyn_float> (nbin_ch);
}
aladyn_float Parameters::get_size_gamma()
{
  return (gammamax - gammamin) / static_cast <aladyn_float> (nbin_gamma);
}
aladyn_float Parameters::get_size_theta()
{
  return (thetamax - thetamin) / static_cast <aladyn_float> (nbin_theta);
}
aladyn_float Parameters::get_size_thetaT()
{
  return (thetaTmax - thetaTmin) / static_cast <aladyn_float> (nbin_thetaT);
}
aladyn_float Parameters::get_size_E()
{
  return (Emax - Emin) / static_cast <aladyn_float> (nbin_E);
}


int Parameters::get_nbin(int column)
{
  if (column == COLUMN_X)           return nbin_x;
  else if (column == COLUMN_Y)      return nbin_y;
  else if (column == COLUMN_Z)      return nbin_z;
  else if (column == COLUMN_PX)     return nbin_px;
  else if (column == COLUMN_PY)     return nbin_py;
  else if (column == COLUMN_PZ)     return nbin_pz;
  else if (column == COLUMN_GAMMA)  return nbin_gamma;
  else if (column == COLUMN_THETA)  return nbin_theta;
  else if (column == COLUMN_E)      return nbin_E;
  else if (column == COLUMN_THETAT) return nbin_thetaT;
  else if (column == COLUMN_TY)     return nbin_ty;
  else if (column == COLUMN_TZ)     return nbin_tz;
  else if (column == COLUMN_W)      return nbin_w;
  else if (column == COLUMN_CH)     return nbin_ch;
  else return INVALID_COLUMN;
}


aladyn_float Parameters::get_min(int column)
{
  if (column == COLUMN_X)           return xmin;
  else if (column == COLUMN_Y)      return ymin;
  else if (column == COLUMN_Z)      return zmin;
  else if (column == COLUMN_PX)     return pxmin;
  else if (column == COLUMN_PY)     return pymin;
  else if (column == COLUMN_PZ)     return pzmin;
  else if (column == COLUMN_GAMMA)  return gammamin;
  else if (column == COLUMN_THETA)  return thetamin;
  else if (column == COLUMN_E)      return Emin;
  else if (column == COLUMN_THETAT) return thetaTmin;
  else if (column == COLUMN_TY)     return tymin;
  else if (column == COLUMN_TZ)     return tzmin;
  else if (column == COLUMN_W)      return wmin;
  else if (column == COLUMN_CH)     return chmin;
  else return INVALID_COLUMN;
}


aladyn_float Parameters::get_max(int column)
{
  if (column == COLUMN_X)           return xmax;
  else if (column == COLUMN_Y)      return ymax;
  else if (column == COLUMN_Z)      return zmax;
  else if (column == COLUMN_PX)     return pxmax;
  else if (column == COLUMN_PY)     return pymax;
  else if (column == COLUMN_PZ)     return pzmax;
  else if (column == COLUMN_GAMMA)  return gammamax;
  else if (column == COLUMN_THETA)  return thetamax;
  else if (column == COLUMN_E)      return Emax;
  else if (column == COLUMN_THETAT) return thetaTmax;
  else if (column == COLUMN_TY)     return tymax;
  else if (column == COLUMN_TZ)     return tzmax;
  else if (column == COLUMN_W)      return wmax;
  else if (column == COLUMN_CH)     return chmax;
  else return INVALID_COLUMN;
}



aladyn_float Parameters::get_size_(int column)
{
  if (column == 0)       return get_size_x();
  else if (column == 1)  return get_size_y();
  else if (column == 2)  return get_size_z();
  else if (column == 3)  return get_size_px();
  else if (column == 4)  return get_size_py();
  else if (column == 5)  return get_size_pz();
  else if (column == 6)  return get_size_gamma();
  else if (column == 7)  return get_size_theta();
  else if (column == 8)  return get_size_E();
  else if (column == 9)  return get_size_thetaT();
  else if (column == 10) return get_size_ty();
  else if (column == 11) return get_size_tz();
  else if (column == 12) return get_size_w();
  else if (column == 13) return get_size_ch();
  else return INVALID_COLUMN;
}



void Parameters::read_params_from_bin_file(const char * filename)
{
  std::FILE * file_in = NULL;
  int fortran_buff;
  size_t fread_size = 0;

  file_in = fopen(filename, "rb");
  if (file_in == NULL) std::cout << "Unable to open file!" << std::endl;
  else std::cout << "File opened to read parameters!" << std::endl;

  fread_size += std::fread(&fortran_buff, sizeof(int), 1, file_in);
  fread_size += std::fread(&nparams, sizeof(int), 1, file_in);
  fread_size += std::fread(&fortran_buff, sizeof(int), 1, file_in);

  if (we_have_to_do_swap) swap_endian_i(&nparams, 1);

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
  we_dont_know_if_we_have_to_do_swap = false;
}


void Parameters::read_params_from_dat_file(std::ifstream& file_dat)
{
  std::string forget_this_line;
  int resampling_factor;
  //  int discriminante_versione_file;
  aladyn_float coord;
  std::getline(file_dat, forget_this_line); // per leggere la riga Integer parameters


  for (int i = 0; i < nparams; i++)
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
  for (int i = 0; i < nparams; i++)
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

  if (!fixed_aladyn_version) {
    file_version = intpar[18];
    // compatibility fixes (sometimes aladyn versions were defined as negatives with this convention)
    if (file_version == -1) file_version = 2;
    if (file_version == -2) file_version = 3;
  }

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
    last_cpu = ncpu_x * ncpu_y * ncpu_z;
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
      header_size_bytes = (nparams + 1) * sizeof(int) + nparams * sizeof(aladyn_float) + 6 * sizeof(int); // there are 6 fortran buffers: two around nparams, two around intpars and two around realpars
    }
    else {
      nptot = (long long int) intpar[16];
      ndv = intpar[17];
      header_size_bytes = (nparams + 1) * sizeof(int) + nparams * sizeof(aladyn_float);
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
      last_cpu = ncpu_x * ncpu_y * ncpu_z;

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
      header_size_bytes = (nparams + 1 + 6) * sizeof(int) + nparams * sizeof(aladyn_float); // +1 for n_par, +6 for fortran buffers (2 around nparams, 2 around intpars, 2 around realpars)
    }

    else {
      nptot = (long long int) intpar[16];
      ndv = intpar[17];
      header_size_bytes = (nparams + 1) * sizeof(int) + nparams * sizeof(aladyn_float);
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

    header_size_bytes = (nparams + 1) * sizeof(int) + nparams * sizeof(aladyn_float);
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
      last_cpu = ncpu_x * ncpu_y * ncpu_z;

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
      last_cpu = ncpu_x * ncpu_y * ncpu_z;
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
      last_cpu = intpar[0];
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
      else printf("Attention: illegal value for ndv\n"), exit(-17);
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
  }
  endian_file = (endianness - 1);
  sim_is_2d = ((ndv == 4 || ndv == 5) && file_version < 3) || (ndv == 6 && file_version >= 3);
  we_dont_know_if_sim_is_2d = false;
}




void Parameters::debug_read_parameters()
{
  std::cout << "Integer parameters" << std::endl;
  for (int i = 0; i < NUMBER_OF_PARAMS_IN_DAT_FILE; i++)
  {
    std::cout << std::setw(14) << intpar[i];
    if (i > 0 && !((i + 1) % 4)) std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Real parameters" << std::endl;
  for (int i = 0; i < NUMBER_OF_PARAMS_IN_DAT_FILE; i++)
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

void Parameters::ask_file_dims()
{
  int ncolonne;
  std::cout << "Which file version is this? (1|2|3|4): ";
  std::cin >> file_version;
  if (file_version < 3)
  {
    std::cout << "How many columns are present in the binary file? (4|5|6|7): ";
    std::cin >> ncolonne;
    if (ncolonne == 6 || ncolonne == 4) file_has_weight = 0;
    else if (ncolonne == 7 || ncolonne == 5) file_has_weight = 1;
    else exit(-5);
    we_dont_know_if_file_has_weight = false;
    if (ncolonne == 4 || ncolonne == 5) sim_is_2d = true;
    we_dont_know_if_sim_is_2d = false;
  }
  else if (file_version == 3)
  {
    std::cout << "is the sim 2D (2) or 3D (3)? ";
    std::cin >> ncolonne;
    if (ncolonne == 2) sim_is_2d = true;
    we_dont_know_if_sim_is_2d = false;
    file_has_weight = true;
    we_dont_know_if_file_has_weight = false;
    file_has_charge = false;
    we_dont_know_if_file_has_charge = false;
  }
  else
  {
    std::cout << "is the sim 2D (2) or 3D (3)? ";
    std::cin >> ncolonne;
    if (ncolonne == 2) sim_is_2d = true;
    we_dont_know_if_sim_is_2d = false;
    file_has_weight = true;
    we_dont_know_if_file_has_weight = false;
    file_has_charge = true;
    we_dont_know_if_file_has_charge = false;
  }

}


void Parameters::ask_file_endianness()
{
  std::cout << "Tell me something about endianness: is it little [x86] (0) or big [ppc] (1)? ";
  std::cin >> endian_file;
  if (endian_file != 1 && endian_file != 0) exit(-4);
}


void Parameters::check_filename(const char *nomefile)
{
  if (nomefile[0] == 'P')
  {
    if (nomefile[1] == 'r')
    {
      if (nomefile[2] == 'p')
      {
        mass_MeV = (aladyn_float)MP_MEV;
        phasespace_file = true;
      }
      else if (nomefile[2] == 'e')
      {
        grid_file = true;
        support_label = "pren";
      }
      else
      {
        std::cout << "Unrecognized file" << std::endl;
        exit(-15);
      }
    }
    else if (nomefile[1] == 'd')
    {
      grid_file = true;
      support_label = "pden";
    }
    else
    {
      std::cout << "Unrecognized file" << std::endl;
      exit(-15);
    }
  }
  else if (nomefile[0] == 'H')
  {
    if (nomefile[1] == 'i')
    {
      if (nomefile[2] == 'p')
      {
        mass_MeV = (aladyn_float)MP_MEV; // fix wrong!
        phasespace_file = true;
      }
      else if (nomefile[2] == 'd')
      {
        grid_file = true;
        support_label = "hidn";
      }
      else if (nomefile[2] == 'e')
      {
        grid_file = true;
        support_label = "hien";
      }
      else
      {
        std::cout << "Unrecognized file" << std::endl;
        exit(-15);
      }
    }
    else if (nomefile[1] == '1')
    {
      if (nomefile[2] == 'p')
      {
        mass_MeV = (aladyn_float)MP_MEV; // fix wrong!
        phasespace_file = true;
      }
      else if (nomefile[2] == 'd')
      {
        grid_file = true;
        support_label = "h1dn";
      }
      else if (nomefile[2] == 'e')
      {
        grid_file = true;
        support_label = "h1en";
      }
      else
      {
        std::cout << "Unrecognized file" << std::endl;
        exit(-15);
      }
    }
    else if (nomefile[1] == '2')
    {
      if (nomefile[2] == 'p')
      {
        mass_MeV = (aladyn_float)MP_MEV; // fix wrong!
        phasespace_file = true;
      }
      else if (nomefile[2] == 'd')
      {
        grid_file = true;
        support_label = "h2dn";
      }
      else if (nomefile[2] == 'e')
      {
        grid_file = true;
        support_label = "h2en";
      }
      else
      {
        std::cout << "Unrecognized file" << std::endl;
        exit(-15);
      }
    }
    else
    {
      std::cout << "Unrecognized file" << std::endl;
      exit(-15);
    }
  }
  else if (nomefile[0] == 'L')
  {
    if (nomefile[1] == 'i')
    {
      if (nomefile[2] == 'p')
      {
        mass_MeV = (aladyn_float)MP_MEV; // fix wrong!
        phasespace_file = true;
      }
      else if (nomefile[2] == 'd')
      {
        grid_file = true;
        support_label = "lidn";
      }
      else if (nomefile[2] == 'e')
      {
        grid_file = true;
        support_label = "lien";
      }
      else
      {
        std::cout << "Unrecognized file" << std::endl;
        exit(-15);
      }
    }
    else
    {
      std::cout << "Unrecognized file" << std::endl;
      exit(-15);
    }
  }
  else if (nomefile[0] == 'E')
  {
    if (nomefile[1] == 'l')
    {
      if (nomefile[2] == 'p')
      {
        mass_MeV = (aladyn_float)ME_MEV;
        phasespace_file = true;
      }
      else if (nomefile[2] == 'e')
      {
        grid_file = true;
        support_label = "elen";
      }
      else
      {
        std::cout << "Unrecognized file" << std::endl;
        exit(-15);
      }
    }
    else if (nomefile[1] == 'x')
    {
      if (nomefile[2] == 'f')
      {
        grid_file = true;
        support_label = "Ex";
      }
      else
      {
        std::cout << "Unrecognized file" << std::endl;
        exit(-15);
      }
    }
    else if (nomefile[1] == 'y')
    {
      if (nomefile[2] == 'f')
      {
        grid_file = true;
        support_label = "Ey";
      }
      else
      {
        std::cout << "Unrecognized file" << std::endl;
        exit(-15);
      }
    }
    else if (nomefile[1] == 'z')
    {
      if (nomefile[2] == 'f')
      {
        grid_file = true;
        support_label = "Ez";
      }
      else
      {
        std::cout << "Unrecognized file" << std::endl;
        exit(-15);
      }
    }
    else if (nomefile[1] == 'd')
    {
      grid_file = true;
      support_label = "eden";
    }
    else
    {
      std::cout << "Unrecognized file" << std::endl;
      exit(-15);
    }
  }
  else if (nomefile[0] == 'B')
  {
    if (nomefile[1] == 'x')
    {
      if (nomefile[2] == 'f')
      {
        grid_file = true;
        support_label = "Bx";
      }
      else
      {
        std::cout << "Unrecognized file" << std::endl;
        exit(-15);
      }
    }
    else if (nomefile[1] == 'y')
    {
      if (nomefile[2] == 'f')
      {
        grid_file = true;
        support_label = "By";
      }
      else
      {
        std::cout << "Unrecognized file" << std::endl;
        exit(-15);
      }
    }
    else if (nomefile[1] == 'z')
    {
      if (nomefile[2] == 'f')
      {
        grid_file = true;
        support_label = "Bz";
      }
      else
      {
        std::cout << "Unrecognized file" << std::endl;
        exit(-15);
      }
    }
    else if (nomefile[1] == 'd')
    {
      grid_file = true;
      support_label = "Bd";
    }
    else
    {
      std::cout << "Unrecognized file" << std::endl;
      exit(-15);
    }
  }
  else if (nomefile[0] == 'J')
  {
    if (nomefile[1] == 'x')
    {
      if (nomefile[2] == 'f')
      {
        grid_file = true;
        support_label = "Jx";
      }
      else
      {
        std::cout << "Unrecognized file" << std::endl;
        exit(-15);
      }
    }
    else if (nomefile[1] == 'y')
    {
      if (nomefile[2] == 'f')
      {
        grid_file = true;
        support_label = "Jy";
      }
      else
      {
        std::cout << "Unrecognized file" << std::endl;
        exit(-15);
      }
    }
    else if (nomefile[1] == 'z')
    {
      if (nomefile[2] == 'f')
      {
        grid_file = true;
        support_label = "Jz";
      }
      else
      {
        std::cout << "Unrecognized file" << std::endl;
        exit(-15);
      }
    }
    else
    {
      std::cout << "Unrecognized file" << std::endl;
      exit(-15);
    }
  }
  else
  {
    std::cout << "Unrecognized file" << std::endl;
    exit(-15);
  }
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
  std::ifstream fileParameters;
  std::string nomefile;
  bool usa_file_params = false;
  bool failed_opening_file;

  for (int i = 2; i < argc; i++)
    /************************************************************************
    We will iterate over argv[] to get the parameters stored inside.
    Note that we're starting on 1 because we don't need to know the
    path of the program, which is stored in argv[0], and the input file,
    which is supposed to be given as the first argument and so is in argv[1]
    ************************************************************************/
  {
    if (argv[i] == "-readParamsfromFile" || argv[i] == "-readParamsFromFile" || argv[i] == "-readParams" || argv[i] == "-readparamsfromfile")
    {
      if (i < argc - 1 && argv[i + 1][0] != '-')
      {
        nomefile = std::string(argv[i + 1]);
        usa_file_params = true;
        i++;
        std::cout << "Using " << nomefile << " as the binning parameters file" << std::endl;
      }
      else
      {
        nomefile = filebasename + ".extremes";
        usa_file_params = true;
        std::cout << "Using " << nomefile << " as the binning parameters file" << std::endl;
      }
    }
    else if (argv[i] == "-params")
    {
      std::cout << "Enabled the output of the parameters from .dat/.bin files" << std::endl;
      out_params = true;
    }
    else if (argv[i] == "-swap")
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
    else if (argv[i] == "-stop")
    {
      last_cpu = atoi(argv[i + 1].c_str());
      if (last_cpu < 1) last_cpu = 1;
      std::cout << "Forced stopping reading at CPU #" << last_cpu << std::endl;
      i++;
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
      if (phasespace_file)
      {
        if (subsample < 1)
        {
          subsample = 1;
          std::cout << "Value for subsampling not valid, disabled" << std::endl;
        }
        else
        {
          std::cout << "Will subsample with a ratio of 1:" << subsample << " if any ASCII dump will be requested" << std::endl;
        }
      }
      else
      {
        std::cout << "Subsample factor is valid only for phase space files and will be used only for ASCII output" << std::endl;
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
      if (!sim_is_2d)
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
      else
      {
        std::cout << "Unable to apply a cut on the grid in 2D, please use -dump_gnuplot" << std::endl;
        if (argv[i + 1][0] != '-') i++;
        out_cutx = 0;
      }
    }

    else if (argv[i] == "-dump_cuty")
    {
      if (!sim_is_2d)
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
      else
      {
        std::cout << "Unable to apply a cut on the grid in 2D, please use -dump_gnuplot" << std::endl;
        if (argv[i + 1][0] != '-') i++;
        out_cuty = 0;
      }
    }

    else if (argv[i] == "-dump_cutz")
    {
      if (!sim_is_2d)
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
      else
      {
        std::cout << "Unable to apply a cut on the grid in 2D, please use -dump_gnuplot" << std::endl;
        if (argv[i + 1][0] != '-') i++;
        out_cutz = 0;
      }
    }

    else if (argv[i] == "-dump_lineoutx")
    {
      std::cout << "You asked to have a lineout of the grid along the x-axis" << std::endl;
      out_lineoutx = 1;
    }

    else if (argv[i] == "-dump_gnuplot")
    {
      if (sim_is_2d)
      {
        std::cout << "You asked to rewrite the 2D grid in ASCII format for gnuplot" << std::endl;
        out_grid2d = 1;
      }
      else
      {
        std::cout << "Unable to write a 3D grid for gnuplot without slicing it, please use dump_cutx/y/z" << std::endl;
        out_grid2d = 0;
      }
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
    else if (argv[i] == "-parameters")
    {
      std::cout << "You asked to write the simulation parameters file" << std::endl;
      out_params = 1;
    }
    else if (argv[i] == "-find_minmax")
    {
      std::cout << "You asked to search for minima and maxima" << std::endl;
      we_have_to_find_minmax = true;
      we_dont_know_if_we_have_to_find_minmax = false;
    }
    else if (argv[i] == "-do_binning")
    {
      std::cout << "You asked to enable plotting functions" << std::endl;
      we_have_to_do_binning = true;
      we_dont_know_if_we_have_to_do_binning = false;
    }
    else if (argv[i] == "-xmin")
    {
      xmin = (aladyn_float)atof(argv[i + 1].c_str());
      xmin_b = false;
      i++;
    }
    else if (argv[i] == "-xmax")
    {
      xmax = (aladyn_float)atof(argv[i + 1].c_str());
      xmax_b = false;
      i++;
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
    else if (argv[i] == "-wmin")
    {
      wmin = (aladyn_float)atof(argv[i + 1].c_str());
      wmin_b = false;
      i++;
    }
    else if (argv[i] == "-wmax")
    {
      wmax = (aladyn_float)atof(argv[i + 1].c_str());
      wmax_b = false;
      i++;
    }
    else if (argv[i] == "-chmin")
    {
      chmin = (aladyn_float)atof(argv[i + 1].c_str());
      chmin_b = false;
      i++;
    }
    else if (argv[i] == "-chmax")
    {
      chmax = (aladyn_float)atof(argv[i + 1].c_str());
      chmax_b = false;
      i++;
    }
    else if (argv[i] == "-ymin")
    {
      ymin = (aladyn_float)atof(argv[i + 1].c_str());
      ymin_b = false;
      i++;
    }
    else if (argv[i] == "-ymax")
    {
      ymax = (aladyn_float)atof(argv[i + 1].c_str());
      ymax_b = false;
      i++;
    }
    else if (argv[i] == "-zmin")
    {
      if (sim_is_2d)
      {
        std::cout << "Unable to apply a cut in z-dimension on a 2D file" << std::endl;
      }
      else
      {
        zmin = (aladyn_float)atof(argv[i + 1].c_str());
        zmin_b = false;
        i++;
      }
    }
    else if (argv[i] == "-zmax")
    {
      if (sim_is_2d)
      {
        std::cout << "Unable to apply a cut in z-dimension on a 2D file" << std::endl;
      }
      else
      {
        zmax = (aladyn_float)atof(argv[i + 1].c_str());
        zmax_b = false;
        i++;
      }
    }
    else if (argv[i] == "-tymin")
    {
      tymin = (aladyn_float)atof(argv[i + 1].c_str());
      tymin_b = false;
      i++;
    }
    else if (argv[i] == "-tymax")
    {
      tymax = (aladyn_float)atof(argv[i + 1].c_str());
      tymax_b = false;
      i++;
    }
    else if (argv[i] == "-tzmin")
    {
      if (sim_is_2d)
      {
        std::cout << "Unable to apply a cut in z-dimension on a 2D file" << std::endl;
      }
      else
      {
        tzmin = (aladyn_float)atof(argv[i + 1].c_str());
        tzmin_b = false;
        i++;
      }
    }
    else if (argv[i] == "-tzmax")
    {
      if (sim_is_2d)
      {
        std::cout << "Unable to apply a cut in z-dimension on a 2D file" << std::endl;
      }
      else
      {
        tzmax = (aladyn_float)atof(argv[i + 1].c_str());
        tzmax_b = false;
        i++;
      }
    }
    else if (argv[i] == "-pxmin")
    {
      pxmin = (aladyn_float)atof(argv[i + 1].c_str());
      pxmin_b = false;
      i++;
    }
    else if (argv[i] == "-pxmax")
    {
      pxmax = (aladyn_float)atof(argv[i + 1].c_str());
      pxmax_b = false;
      i++;
    }
    else if (argv[i] == "-pymin")
    {
      pymin = (aladyn_float)atof(argv[i + 1].c_str());
      pymin_b = false;
      i++;
    }
    else if (argv[i] == "-pymax")
    {
      pymax = (aladyn_float)atof(argv[i + 1].c_str());
      pymax_b = false;
      i++;
    }
    else if (argv[i] == "-pzmin")
    {
      if (sim_is_2d)
      {
        std::cout << "Unable to apply a cut in z-dimension on a 2D file" << std::endl;
      }
      else
      {
        pzmin = (aladyn_float)atof(argv[i + 1].c_str());
        pzmin_b = false;
        i++;
      }
    }
    else if (argv[i] == "-pzmax")
    {
      if (sim_is_2d)
      {
        std::cout << "Unable to apply a cut in z-dimension on a 2D file" << std::endl;
      }
      else
      {
        pzmax = (aladyn_float)atof(argv[i + 1].c_str());
        pzmax_b = false;
        i++;
      }
    }
    else if (argv[i] == "-thetamin")
    {
      thetamin = (aladyn_float)atof(argv[i + 1].c_str());
      thetamin_b = false;
      i++;
    }
    else if (argv[i] == "-thetamax")
    {
      thetamax = (aladyn_float)atof(argv[i + 1].c_str());
      thetamax_b = false;
      i++;
    }
    else if (argv[i] == "-thetaTmin")
    {
      thetaTmin = (aladyn_float)atof(argv[i + 1].c_str());
      thetaTmin_b = false;
      i++;
    }
    else if (argv[i] == "-thetaTmax")
    {
      thetaTmax = (aladyn_float)atof(argv[i + 1].c_str());
      thetaTmax_b = false;
      i++;
    }
    else if (argv[i] == "-gammamin")
    {
      gammamin = (aladyn_float)atof(argv[i + 1].c_str());
      gammamin_b = false;
      i++;
    }
    else if (argv[i] == "-gammamax")
    {
      gammamax = (aladyn_float)atof(argv[i + 1].c_str());
      gammamax_b = false;
      i++;
    }
    else if (argv[i] == "-Emin")
    {
      Emin = (aladyn_float)atof(argv[i + 1].c_str());
      Emin_b = false;
      i++;
    }
    else if (argv[i] == "-Emax")
    {
      Emax = (aladyn_float)atof(argv[i + 1].c_str());
      Emax_b = false;
      i++;
    }
    else if (argv[i] == "-plot_xy")
    {
      if (we_dont_know_if_we_have_to_do_binning) {
        std::cout << "You asked to enable plotting functions" << std::endl;
        we_have_to_do_binning = 1;
        we_dont_know_if_we_have_to_do_binning = false;
      }
      do_plot_xy = 1;
    }
    else if (argv[i] == "-plot_xw")
    {
      if (file_has_weight)
      {
        if (we_dont_know_if_we_have_to_do_binning) {
          std::cout << "You asked to enable plotting functions" << std::endl;
          we_have_to_do_binning = 1;
          we_dont_know_if_we_have_to_do_binning = false;
        }
        do_plot_xw = 1;
      }
      else
      {
        std::cout << "Unable to do a plot with weight using a file without weight!" << std::endl;
      }
    }
    else if (argv[i] == "-plot_xz")
    {
      if (sim_is_2d)
      {
        std::cout << "Unable to do a plot with z using a 2D file" << std::endl;
      }
      else
      {
        if (we_dont_know_if_we_have_to_do_binning) {
          std::cout << "You asked to enable plotting functions" << std::endl;
          we_have_to_do_binning = 1;
          we_dont_know_if_we_have_to_do_binning = false;
        }
        do_plot_xz = 1;
      }
    }
    else if (argv[i] == "-plot_yz")
    {
      if (sim_is_2d)
      {
        std::cout << "Unable to do a plot with z using a 2D file" << std::endl;
      }
      else
      {
        if (we_dont_know_if_we_have_to_do_binning) {
          std::cout << "You asked to enable plotting functions" << std::endl;
          we_have_to_do_binning = 1;
          we_dont_know_if_we_have_to_do_binning = false;
        }
        do_plot_yz = 1;
      }
    }
    else if (argv[i] == "-plot_rcf")
    {
      if (sim_is_2d)
      {
        std::cout << "Unable to do a plot with tz using a 2D file" << std::endl;
      }
      else
      {
        if (we_dont_know_if_we_have_to_do_binning) {
          std::cout << "You asked to enable plotting functions" << std::endl;
          we_have_to_do_binning = 1;
          we_dont_know_if_we_have_to_do_binning = false;
        }
        do_plot_rcf = 1;
      }
    }
    else if (argv[i] == "-plot_xpx")
    {
      if (we_dont_know_if_we_have_to_do_binning) {
        std::cout << "You asked to enable plotting functions" << std::endl;
        we_have_to_do_binning = 1;
        we_dont_know_if_we_have_to_do_binning = false;
      }
      do_plot_xpx = 1;
    }
    else if (argv[i] == "-plot_xpy")
    {
      if (we_dont_know_if_we_have_to_do_binning) {
        std::cout << "You asked to enable plotting functions" << std::endl;
        we_have_to_do_binning = 1;
        we_dont_know_if_we_have_to_do_binning = false;
      }
      do_plot_xpy = 1;
    }
    else if (argv[i] == "-plot_xpz")
    {
      if (sim_is_2d)
      {
        std::cout << "Unable to do a plot with z using a 2D file" << std::endl;
      }
      else
      {
        if (we_dont_know_if_we_have_to_do_binning) {
          std::cout << "You asked to enable plotting functions" << std::endl;
          we_have_to_do_binning = 1;
          we_dont_know_if_we_have_to_do_binning = false;
        }
        do_plot_xpz = 1;
      }
    }
    else if (argv[i] == "-plot_ypx")
    {
      if (we_dont_know_if_we_have_to_do_binning) {
        std::cout << "You asked to enable plotting functions" << std::endl;
        we_have_to_do_binning = 1;
        we_dont_know_if_we_have_to_do_binning = false;
      }
      do_plot_ypx = 1;
    }
    else if (argv[i] == "-plot_ypy")
    {
      if (we_dont_know_if_we_have_to_do_binning) {
        std::cout << "You asked to enable plotting functions" << std::endl;
        we_have_to_do_binning = 1;
        we_dont_know_if_we_have_to_do_binning = false;
      }
      do_plot_ypy = 1;
    }
    else if (argv[i] == "-plot_ypz")
    {
      if (sim_is_2d)
      {
        std::cout << "Unable to do a plot with z using a 2D file" << std::endl;
      }
      else
      {
        if (we_dont_know_if_we_have_to_do_binning) {
          std::cout << "You asked to enable plotting functions" << std::endl;
          we_have_to_do_binning = 1;
          we_dont_know_if_we_have_to_do_binning = false;
        }
        do_plot_ypz = 1;
      }
    }
    else if (argv[i] == "-plot_zpx")
    {
      if (sim_is_2d)
      {
        std::cout << "Unable to do a plot with z using a 2D file" << std::endl;
      }
      else
      {
        if (we_dont_know_if_we_have_to_do_binning) {
          std::cout << "You asked to enable plotting functions" << std::endl;
          we_have_to_do_binning = 1;
          we_dont_know_if_we_have_to_do_binning = false;
        }
        do_plot_zpx = 1;
      }
    }
    else if (argv[i] == "-plot_zpy")
    {
      if (sim_is_2d)
      {
        std::cout << "Unable to do a plot with z using a 2D file" << std::endl;
      }
      else
      {
        if (we_dont_know_if_we_have_to_do_binning) {
          std::cout << "You asked to enable plotting functions" << std::endl;
          we_have_to_do_binning = 1;
          we_dont_know_if_we_have_to_do_binning = false;
        }
        do_plot_zpy = 1;
      }
    }
    else if (argv[i] == "-plot_zpz")
    {
      if (sim_is_2d)
      {
        std::cout << "Unable to do a plot with z using a 2D file" << std::endl;
      }
      else
      {
        if (we_dont_know_if_we_have_to_do_binning) {
          std::cout << "You asked to enable plotting functions" << std::endl;
          we_have_to_do_binning = 1;
          we_dont_know_if_we_have_to_do_binning = false;
        }
        do_plot_zpz = 1;
      }
    }
    else if (argv[i] == "-plot_pxpy")
    {
      if (we_dont_know_if_we_have_to_do_binning) {
        std::cout << "You asked to enable plotting functions" << std::endl;
        we_have_to_do_binning = 1;
        we_dont_know_if_we_have_to_do_binning = false;
      }
      do_plot_pxpy = 1;
    }
    else if (argv[i] == "-plot_pxpz")
    {
      if (sim_is_2d)
      {
        std::cout << "Unable to do a plot with z using a 2D file" << std::endl;
      }
      else
      {
        if (we_dont_know_if_we_have_to_do_binning) {
          std::cout << "You asked to enable plotting functions" << std::endl;
          we_have_to_do_binning = 1;
          we_dont_know_if_we_have_to_do_binning = false;
        }
        do_plot_pxpz = 1;
      }
    }
    else if (argv[i] == "-plot_pypz")
    {
      if (sim_is_2d)
      {
        std::cout << "Unable to do a plot with z using a 2D file" << std::endl;
      }
      else
      {
        if (we_dont_know_if_we_have_to_do_binning) {
          std::cout << "You asked to enable plotting functions" << std::endl;
          we_have_to_do_binning = 1;
          we_dont_know_if_we_have_to_do_binning = false;
        }
        do_plot_pypz = 1;
      }
    }
    else if (argv[i] == "-plot_etheta")
    {
      if (we_dont_know_if_we_have_to_do_binning) {
        std::cout << "You asked to enable plotting functions" << std::endl;
        we_have_to_do_binning = 1;
        we_dont_know_if_we_have_to_do_binning = false;
      }
      do_plot_Etheta = 1;
    }
    else if (argv[i] == "-plot_ethetaT")
    {
      if (we_dont_know_if_we_have_to_do_binning) {
        std::cout << "You asked to enable plotting functions" << std::endl;
        we_have_to_do_binning = 1;
        we_dont_know_if_we_have_to_do_binning = false;
      }
      do_plot_EthetaT = 1;
    }
    else if (argv[i] == "-plot_espec")
    {
      if (we_dont_know_if_we_have_to_do_binning) {
        std::cout << "You asked to enable plotting functions" << std::endl;
        we_have_to_do_binning = 1;
        we_dont_know_if_we_have_to_do_binning = false;
      }
      do_plot_Espec = 1;
    }
    else if (argv[i] == "-plot_chspec")
    {
      if (file_has_charge)
      {
        if (we_dont_know_if_we_have_to_do_binning) {
          std::cout << "You asked to enable plotting functions" << std::endl;
          we_have_to_do_binning = 1;
          we_dont_know_if_we_have_to_do_binning = false;
        }
        do_plot_chspec = 1;
      }
      else
      {
        std::cout << "Unable to do a plot with charge using a file without charges!" << std::endl;
      }
    }
    else if (argv[i] == "-plot_wspec")
    {
      if (file_has_weight)
      {
        if (we_dont_know_if_we_have_to_do_binning) {
          std::cout << "You asked to enable plotting functions" << std::endl;
          we_have_to_do_binning = 1;
          we_dont_know_if_we_have_to_do_binning = false;
        }
        do_plot_wspec = 1;
      }
      else
      {
        std::cout << "Unable to do a plot with weight using a file without weights!" << std::endl;
      }
    }
    else if (argv[i] == "-plot_thetaspec")
    {
      if (we_dont_know_if_we_have_to_do_binning) {
        std::cout << "You asked to enable plotting functions" << std::endl;
        we_have_to_do_binning = 1;
        we_dont_know_if_we_have_to_do_binning = false;
      }
      do_plot_thetaspec = 1;
    }
    else if (argv[i] == "-plot_thetaTspec")
    {
      if (we_dont_know_if_we_have_to_do_binning) {
        std::cout << "You asked to enable plotting functions" << std::endl;
        we_have_to_do_binning = 1;
        we_dont_know_if_we_have_to_do_binning = false;
      }
      do_plot_thetaTspec = 1;
    }
    else if (argv[i] == "-nbin")
    {
      nbin = atoi(argv[i + 1].c_str());
      if (nbin_x_b) nbin_x = nbin;
      if (nbin_y_b) nbin_y = nbin;
      if (nbin_z_b) nbin_z = nbin;
      if (nbin_px_b) nbin_px = nbin;
      if (nbin_py_b) nbin_py = nbin;
      if (nbin_pz_b) nbin_pz = nbin;
      if (nbin_ty_b) nbin_ty = nbin;
      if (nbin_tz_b) nbin_tz = nbin;
      if (nbin_gamma_b) nbin_gamma = nbin;
      if (nbin_theta_b) nbin_theta = nbin;
      if (nbin_theta_b) nbin_thetaT = nbin;
      if (nbin_E_b) nbin_E = nbin;
      if (nbin_w_b) nbin_w = nbin;
      if (nbin_ch_b) nbin_ch = nbin;
      nbin_b = false;
      i++;
    }
    else if (argv[i] == "-nbinx")
    {
      nbin_x = atoi(argv[i + 1].c_str());
      nbin_x_b = false;
      nbin_b = false;
      i++;
    }
    else if (argv[i] == "-nbiny")
    {
      nbin_y = atoi(argv[i + 1].c_str());
      nbin_y_b = false;
      nbin_b = false;
      i++;
    }
    else if (argv[i] == "-nbinz")
    {
      nbin_z = atoi(argv[i + 1].c_str());
      nbin_z_b = false;
      nbin_b = false;
      i++;
    }
    else if (argv[i] == "-nbinty")
    {
      nbin_ty = atoi(argv[i + 1].c_str());
      nbin_ty_b = false;
      nbin_b = false;
      i++;
    }
    else if (argv[i] == "-nbintz")
    {
      nbin_tz = atoi(argv[i + 1].c_str());
      nbin_tz_b = false;
      nbin_b = false;
      i++;
    }
    else if (argv[i] == "-nbinpx")
    {
      nbin_px = atoi(argv[i + 1].c_str());
      nbin_px_b = false;
      nbin_b = false;
      i++;
    }
    else if (argv[i] == "-nbinpy")
    {
      nbin_py = atoi(argv[i + 1].c_str());
      nbin_py_b = false;
      nbin_b = false;
      i++;
    }
    else if (argv[i] == "-nbinpz")
    {
      nbin_pz = atoi(argv[i + 1].c_str());
      nbin_pz_b = false;
      nbin_b = false;
      i++;
    }
    else if (argv[i] == "-nbintheta")
    {
      nbin_theta = atoi(argv[i + 1].c_str());
      nbin_theta_b = false;
      nbin_b = false;
      i++;
    }
    else if (argv[i] == "-nbinthetaT")
    {
      nbin_thetaT = atoi(argv[i + 1].c_str());
      nbin_thetaT_b = false;
      nbin_b = false;
      i++;
    }
    else if (argv[i] == "-nbingamma")
    {
      nbin_gamma = atoi(argv[i + 1].c_str());
      nbin_gamma_b = false;
      nbin_b = false;
      i++;
    }
    else if (argv[i] == "-nbinE")
    {
      nbin_E = atoi(argv[i + 1].c_str());
      nbin_E_b = false;
      nbin_b = false;
      i++;
    }
    else if (argv[i] == "-nbinw")
    {
      nbin_w = atoi(argv[i + 1].c_str());
      nbin_w_b = false;
      nbin_b = false;
      i++;
    }
    else if (argv[i] == "-nbinch")
    {
      nbin_ch = atoi(argv[i + 1].c_str());
      nbin_ch_b = false;
      nbin_b = false;
      i++;
    }
    else if (argv[i] == "-dontask")
    {
      do_not_ask_missing = true;
    }
  }


  std::string nomepar, evita, leggi;

  if (phasespace_file)
  {
    if (usa_file_params)
    {
      fileParameters.open(nomefile.c_str());
      failed_opening_file = fileParameters.fail();
      if (failed_opening_file)
      {
        std::cout << "Unable to open file " << nomefile << " which contains the binning parameters" << std::endl;
        exit(50);
      }
      while (!fileParameters.eof())
      {
        fileParameters >> nomepar >> evita >> leggi;
        if ((nomepar == "xmin" || nomepar == "XMIN") && xmin_b)
        {
          xmin = (aladyn_float)std::atof(leggi.c_str());
          xmin_b = false;
        }
        else if ((nomepar == "xmax" || nomepar == "XMAX") && xmax_b)
        {
          xmax = (aladyn_float)std::atof(leggi.c_str());
          xmax_b = false;
        }
        else if ((nomepar == "ymin" || nomepar == "YMIN") && ymin_b)
        {
          ymin = (aladyn_float)std::atof(leggi.c_str());
          ymin_b = false;
        }
        else if ((nomepar == "ymax" || nomepar == "YMAX") && ymax_b)
        {
          ymax = (aladyn_float)std::atof(leggi.c_str());
          ymax_b = false;
        }
        else if ((nomepar == "zmin" || nomepar == "ZMIN") && zmin_b)
        {
          zmin = (aladyn_float)std::atof(leggi.c_str());
          zmin_b = false;
        }
        else if ((nomepar == "zmax" || nomepar == "ZMAX") && zmax_b)
        {
          zmax = (aladyn_float)std::atof(leggi.c_str());
          zmax_b = false;
        }
        else if ((nomepar == "tymin" || nomepar == "TYMIN") && tymin_b)
        {
          tymin = (aladyn_float)std::atof(leggi.c_str());
          tymin_b = false;
        }
        else if ((nomepar == "tymax" || nomepar == "TYMAX") && tymax_b)
        {
          tymax = (aladyn_float)std::atof(leggi.c_str());
          tymax_b = false;
        }
        else if ((nomepar == "tzmin" || nomepar == "TZMIN") && tzmin_b)
        {
          tzmin = (aladyn_float)std::atof(leggi.c_str());
          tzmin_b = false;
        }
        else if ((nomepar == "tzmax" || nomepar == "TZMAX") && tzmax_b)
        {
          tzmax = (aladyn_float)std::atof(leggi.c_str());
          tzmax_b = false;
        }
        else if ((nomepar == "pxmin" || nomepar == "PXMIN") && pxmin_b)
        {
          pxmin = (aladyn_float)std::atof(leggi.c_str());
          pxmin_b = false;
        }
        else if ((nomepar == "pxmax" || nomepar == "PXMAX") && pxmax_b)
        {
          pxmax = (aladyn_float)std::atof(leggi.c_str());
          pxmax_b = false;
        }
        else if ((nomepar == "pymin" || nomepar == "PYMIN") && pymin_b)
        {
          pymin = (aladyn_float)std::atof(leggi.c_str());
          pymin_b = false;
        }
        else if ((nomepar == "pymax" || nomepar == "PYMAX") && pymax_b)
        {
          pymax = (aladyn_float)std::atof(leggi.c_str());
          pymax_b = false;
        }
        else if ((nomepar == "pzmin" || nomepar == "PZMIN") && pzmin_b)
        {
          pzmin = (aladyn_float)std::atof(leggi.c_str());
          pzmin_b = false;
        }
        else if ((nomepar == "pzmax" || nomepar == "PZMAX") && pzmax_b)
        {
          pzmax = (aladyn_float)std::atof(leggi.c_str());
          pzmax_b = false;
        }
        else if ((nomepar == "gammamin" || nomepar == "GAMMAMIN") && gammamin_b)
        {
          gammamin = (aladyn_float)std::atof(leggi.c_str());
          gammamin_b = false;
        }
        else if ((nomepar == "gammamax" || nomepar == "GAMMAMAX") && gammamax_b)
        {
          gammamax = (aladyn_float)std::atof(leggi.c_str());
          gammamax_b = false;
        }
        else if ((nomepar == "thetamin" || nomepar == "THETAMIN") && thetamin_b)
        {
          thetamin = (aladyn_float)std::atof(leggi.c_str());
          thetamin_b = false;
        }
        else if ((nomepar == "thetamax" || nomepar == "THETAMAX") && thetamax_b)
        {
          thetamax = (aladyn_float)std::atof(leggi.c_str());
          thetamax_b = false;
        }
        else if ((nomepar == "thetaradmin" || nomepar == "THETARADMIN") && thetaTmin_b)
        {
          thetaTmin = (aladyn_float)std::atof(leggi.c_str());
          thetaTmin_b = false;
        }
        else if ((nomepar == "thetaradmax" || nomepar == "THETARADMAX") && thetaTmax_b)
        {
          thetaTmax = (aladyn_float)std::atof(leggi.c_str());
          thetaTmax_b = false;
        }
        else if ((nomepar == "emin" || nomepar == "EMIN") && Emin_b)
        {
          Emin = (aladyn_float)std::atof(leggi.c_str());
          Emin_b = false;
        }
        else if ((nomepar == "emax" || nomepar == "EMAX") && Emax_b)
        {
          Emax = (aladyn_float)std::atof(leggi.c_str());
          Emax_b = false;
        }
        else if ((nomepar == "wmin" || nomepar == "WMIN") && wmin_b)
        {
          wmin = (aladyn_float)std::atof(leggi.c_str());
          wmin_b = false;
        }
        else if ((nomepar == "wmax" || nomepar == "WMAX") && wmax_b)
        {
          wmax = (aladyn_float)std::atof(leggi.c_str());
          wmax_b = false;
        }
        else if ((nomepar == "chmin" || nomepar == "CHMIN") && chmin_b)
        {
          chmin = (aladyn_float)std::atof(leggi.c_str());
          chmin_b = false;
        }
        else if ((nomepar == "chmax" || nomepar == "CHMAX") && chmax_b)
        {
          chmax = (aladyn_float)std::atof(leggi.c_str());
          chmax_b = false;
        }
        /*
        else
        {
        std::cout << "Parametro " << nomepar << " non riconosciuto." << std::endl;
        }
        */
      }
      fileParameters.close();
    }
    if (we_dont_know_if_we_have_to_find_minmax && !do_not_ask_missing)
    {
      std::cout << "Do you want to find minimum and maximum values for common parameters? 0 (NO) - 1 (YES): ";
      std::cin >> we_have_to_find_minmax;
      we_dont_know_if_we_have_to_find_minmax = false;
    }
    if (we_dont_know_if_we_have_to_do_binning && !do_not_ask_missing)
    {
      std::cout << "Do you want to do a histogram (data binning)? 0 (NO) - 1 (YES): ";
      std::cin >> we_have_to_do_binning;
      we_dont_know_if_we_have_to_do_binning = false;
    }
    if (we_have_to_do_binning == 1 && nbin_b && !do_not_ask_missing)
    {
      std::cout << "How many bins per axis? (common ALaDyn choice: 120): ";
      std::cin >> nbin;
      nbin_x = nbin_y = nbin_z = nbin_px = nbin_py = nbin_pz = nbin_E = nbin_theta = nbin_thetaT = nbin_ty = nbin_tz = nbin_w = nbin_ch = nbin;
    }

    if (we_have_to_do_binning == 1 && !do_not_ask_missing)
    {
      std::cout << "Do you want to do an x-px plot? 0 (NO) - 1 (YES): ";
      std::cin >> do_plot_xpx;
      std::cout << "E-theta (deg)? 0 (NO) - 1 (YES): ";
      std::cin >> do_plot_Etheta;
      std::cout << "E-theta (rad)? 0 (NO) - 1 (YES): ";
      std::cin >> do_plot_EthetaT;
      std::cout << "Energy spectrum? 0 (NO) - 1 (YES): ";
      std::cin >> do_plot_Espec;
      std::cout << "RFC plot? 0 (NO) - 1 (YES): ";
      std::cin >> do_plot_rcf;
      if (do_plot_xpx)
      {
        if (xmin_b)
        {
          std::cout << "xmin = ";
          std::cin >> xmin;
          xmin_b = false;
        }
        if (xmax_b)
        {
          std::cout << "xmax = ";
          std::cin >> xmax;
          xmax_b = false;
        }
        if (nbin_x_b)
        {
          std::cout << "nbin_x = ";
          std::cin >> nbin_x;
          nbin_x_b = false;
        }
        if (pxmin_b)
        {
          std::cout << "pxmin = ";
          std::cin >> pxmin;
          pxmin_b = false;
        }
        if (pxmax_b)
        {
          std::cout << "pxmax = ";
          std::cin >> pxmax;
          pxmax_b = false;
        }
        if (nbin_px_b)
        {
          std::cout << "nbin_px = ";
          std::cin >> nbin_px;
          nbin_px_b = false;
        }
      }
      if (do_plot_Etheta || do_plot_Espec || do_plot_EthetaT)
      {
        if (Emin_b)
        {
          std::cout << "Emin = ";
          std::cin >> Emin;
          Emin_b = false;
        }
        if (Emax_b)
        {
          std::cout << "Emax = ";
          std::cin >> Emax;
          Emax_b = false;
        }
        if (nbin_E_b)
        {
          std::cout << "nbin_E = ";
          std::cin >> nbin_E;
          nbin_E_b = false;
        }
      }
      if (do_plot_Etheta)
      {
        if (thetamin_b)
        {
          std::cout << "thetamin = ";
          std::cin >> thetamin;
          thetamin_b = false;
        }
        if (thetamax_b)
        {
          std::cout << "thetamax = ";
          std::cin >> thetamax;
          thetamax_b = false;
        }
        if (nbin_theta_b)
        {
          std::cout << "nbin_theta = ";
          std::cin >> nbin_theta;
          nbin_theta_b = false;
        }
      }
      if (do_plot_EthetaT)
      {
        if (thetaTmin_b)
        {
          std::cout << "thetaRADmin = ";
          std::cin >> thetaTmin;
          thetaTmin_b = false;
        }
        if (thetaTmax_b)
        {
          std::cout << "thetaRADmax = ";
          std::cin >> thetaTmax;
          thetaTmax_b = false;
        }
        if (nbin_thetaT_b)
        {
          std::cout << "nbin_thetaRAD = ";
          std::cin >> nbin_thetaT;
          nbin_thetaT_b = false;
        }
      }
      if (do_plot_rcf)
      {
        if (tymin_b)
        {
          std::cout << "tymin = ";
          std::cin >> tymin;
          tymin_b = false;
        }
        if (tymax_b)
        {
          std::cout << "tymax = ";
          std::cin >> tymax;
          tymax_b = false;
        }
        if (nbin_ty_b)
        {
          std::cout << "nbin_ty = ";
          std::cin >> nbin_ty;
          nbin_ty_b = false;
        }
        if (tzmin_b)
        {
          std::cout << "tzmin = ";
          std::cin >> tzmin;
          tzmin_b = false;
        }
        if (tzmax_b)
        {
          std::cout << "tzmax = ";
          std::cin >> tzmax;
          tzmax_b = false;
        }
        if (nbin_tz_b)
        {
          std::cout << "nbin_tz = ";
          std::cin >> nbin_tz;
          nbin_tz_b = false;
        }
      }
    }




#ifdef ENABLE_DEBUG
    std::cout << "In file " << nomefile << " I found these binning parameters:" << std::endl;
    std::cout << "XMIN = " << xmin << std::endl;
    std::cout << "XMAX = " << xmax << std::endl;
    std::cout << "YMIN = " << ymin << std::endl;
    std::cout << "YMAX = " << ymax << std::endl;
    std::cout << "ZMIN = " << zmin << std::endl;
    std::cout << "ZMAX = " << zmax << std::endl;
    std::cout << "PXMIN = " << pxmin << std::endl;
    std::cout << "PXMAX = " << pxmax << std::endl;
    std::cout << "PYMIN = " << pymin << std::endl;
    std::cout << "PYMAX = " << pymax << std::endl;
    std::cout << "PZMIN = " << pzmin << std::endl;
    std::cout << "PZMAX = " << pzmax << std::endl;
    std::cout << "GAMMAMIN = " << gammamin << std::endl;
    std::cout << "GAMMAMAX = " << gammamax << std::endl;
    std::cout << "THETAMIN = " << thetamin << std::endl;
    std::cout << "THETAMAX = " << thetamax << std::endl;
    std::cout << "THETARADMIN = " << thetaTmin << std::endl;
    std::cout << "THETARADMAX = " << thetaTmax << std::endl;
    std::cout << "EMIN = " << Emin << std::endl;
    std::cout << "EMAX = " << Emax << std::endl;
    std::cout << "TYMIN = " << tymin << std::endl;
    std::cout << "TYMAX = " << tymax << std::endl;
    std::cout << "TZMIN = " << tzmin << std::endl;
    std::cout << "TZMAX = " << tzmax << std::endl;
    std::cout << "WMIN = " << wmin << std::endl;
    std::cout << "WMAX = " << wmax << std::endl;
    std::cout << "CHMIN = " << chmin << std::endl;
    std::cout << "CHMAX = " << chmax << std::endl;
#endif
    if (!out_vtk && !do_not_ask_missing)
    {
      std::cout << "Do you want the .vtk binary output?? 0 (NO) - 1 (YES): ";
      std::cin >> out_vtk;
    }
    if (!out_clean_bin && !do_not_ask_missing)
    {
      std::cout << "Do you want the .bin output, cleaned from Fortran metadata? 0 (NO) - 1 (YES): ";
      std::cin >> out_clean_bin;
    }
    if (!out_ppg && !do_not_ask_missing)
    {
      std::cout << "Do you want the .ppg output for Propaga? 0 (NO) - 1 (YES): ";
      std::cin >> out_ppg;
    }
    if (!out_xyze && !do_not_ask_missing)
    {
      std::cout << "Do you want a .txt file with x y (z) and Energy? 0 (NO) - 1 (YES): ";
      std::cin >> out_xyze;
    }
    if (!out_csv && !do_not_ask_missing)
    {
      std::cout << "Do you want a .csv file for Paraview? 0 (NO) - 1 (YES): ";
      std::cin >> out_csv;
    }
    if (!out_params && !do_not_ask_missing)
    {
      std::cout << "Do you want a .dat files with simulation parameters? 0 (NO) - 1 (YES): ";
      std::cin >> out_params;
    }
  }
  else
  {
    if (!out_vtk && !do_not_ask_missing)
    {
      std::cout << "Do you want the .vtk binary output?? 0 (NO) - 1 (YES): ";
      std::cin >> out_vtk;
    }
    if (!out_vtk_nostretch && !do_not_ask_missing)
    {
      std::cout << "Do you want the .vtk binary output only for the unstretched part of the grid? 0 (NO) - 1 (YES): ";
      std::cin >> out_vtk_nostretch;
    }
    if (out_vtk_nostretch && !do_not_ask_missing)
    {
      std::cout << "Is the grid stretched also along x? 0 (NO) - 1 (YES): ";
      std::cin >> stretched_along_x;
    }
    if (!sim_is_2d)
    {
      if (!out_cutx && !do_not_ask_missing)
      {
        std::cout << "Do you want a .txt from a slice along x, for gnuplot? 0 (NO) - 1 (YES): ";
        std::cin >> out_cutx;
        aladyn_float posizione_taglio = 0.0;
        if (out_cutx == 1)
        {
          std::cout << "Please tell me at which position should I cut the slice (in micrometers): ";
          std::cin >> posizione_taglio;
          where_to_cut_grid_along_x.push_back(posizione_taglio);
        }
      }
      if (!out_cuty && !do_not_ask_missing)
      {
        std::cout << "Do you want a .txt from a slice along y, for gnuplot? 0 (NO) - 1 (YES): ";
        std::cin >> out_cuty;
        aladyn_float posizione_taglio = 0.0;
        if (out_cuty == 1)
        {
          std::cout << "Please tell me at which position should I cut the slice (in micrometers): ";
          std::cin >> posizione_taglio;
          where_to_cut_grid_along_y.push_back(posizione_taglio);
        }
      }
      if (!out_cutz && !do_not_ask_missing)
      {
        std::cout << "Do you want a .txt from a slice along z, for gnuplot? 0 (NO) - 1 (YES): ";
        std::cin >> out_cutz;
        aladyn_float posizione_taglio = 0.0;
        if (out_cutz == 1)
        {
          std::cout << "Please tell me at which position should I cut the slice (in micrometers): ";
          std::cin >> posizione_taglio;
          where_to_cut_grid_along_z.push_back(posizione_taglio);
        }
      }
      out_grid2d = 0;
    }
    else
    {
      if (!out_grid2d && !do_not_ask_missing)
      {
        std::cout << "Do you want a .txt file for gnuplot? 0 (NO) - 1 (YES): ";
        std::cin >> out_grid2d;
      }
      out_cutx = 0;
      out_cuty = 0;
      out_cutz = 0;
    }
    if (!out_params && !do_not_ask_missing)
    {
      std::cout << "Do you want a .dat files with simulation parameters? 0 (NO) - 1 (YES): ";
      std::cin >> out_params;
    }
    if (!out_lineoutx && !do_not_ask_missing)
    {
      std::cout << "Do you want a lineout along x? 0 (NO) - 1 (YES): ";
      std::cin >> out_lineoutx;
    }
  }
}



bool Parameters::check_params()
{
  bool test = true;
  if (we_dont_know_if_we_have_to_do_swap && do_not_ask_missing)
  {
    we_have_to_do_swap = 0;
    we_dont_know_if_we_have_to_do_swap = false;
    test = true;
  }
  else test = false;
  if (phasespace_file)
  {
    if (xmin > xmax)
    {
      printf("Warning: xmin > xmax\n");
      test = false;
    }
    if (ymin > ymax)
    {
      printf("Warning: ymin > ymax\n");
      test = false;
    }
    if (zmin > zmax)
    {
      printf("Warning: zmin > zmax\n");
      test = false;
    }
    if (pxmin > pxmax)
    {
      printf("Warning: pxmin > pxmax\n");
      test = false;
    }
    if (pymin > pymax)
    {
      printf("Warning: pymin > pymax\n");
      test = false;
    }
    if (pzmin > pzmax)
    {
      printf("Warning: pzmin > pzmax\n");
      test = false;
    }
    if (Emin > Emax)
    {
      printf("Warning: Emin > Emax\n");
      test = false;
    }
    if (thetamin > thetamax)
    {
      printf("Warning: thetamin > thetamax\n");
      test = false;
    }
    if (thetaTmin > thetaTmax)
    {
      printf("Warning: thetaTmin > thetaTmax\n");
      test = false;
    }
    if (tymin > tymax)
    {
      printf("Warning: tymin > tymax\n");
      test = false;
    }
    if (tzmin > tzmax)
    {
      printf("Warning: tzmin > tzmax\n");
      test = false;
    }
    if (wmin > wmax)
    {
      printf("Warning: wmin > wmax\n");
      test = false;
    }
    if (chmin > chmax)
    {
      printf("Warning: chmin > chmax\n");
      test = false;
    }
    if (nbin_x <= 0)
    {
      printf("Warning: nbin_x < 0\n");
      test = false;
    }
    if (nbin_y <= 0)
    {
      printf("Warning: nbin_y < 0\n");
      test = false;
    }
    if (nbin_z <= 0)
    {
      printf("Warning: nbin_z < 0\n");
      test = false;
    }
    if (nbin_px <= 0)
    {
      printf("Warning: nbin_px < 0\n");
      test = false;
    }
    if (nbin_py <= 0)
    {
      printf("Warning: nbin_py < 0\n");
      test = false;
    }
    if (nbin_pz <= 0)
    {
      printf("Warning: nbin_pz < 0\n");
      test = false;
    }
    if (nbin_E <= 0)
    {
      printf("Warning: nbin_E < 0\n");
      test = false;
    }
    if (nbin_theta <= 0)
    {
      printf("Warning: nbin_theta < 0\n");
      test = false;
    }
    if (nbin_thetaT <= 0)
    {
      printf("Warning: nbin_thetaT < 0\n");
      test = false;
    }
    if (nbin_ty <= 0)
    {
      printf("Warning: nbin_ty < 0\n");
      test = false;
    }
    if (nbin_tz <= 0)
    {
      printf("Warning: nbin_tz < 0\n");
      test = false;
    }
    if (nbin_w <= 0)
    {
      printf("Warning: nbin_w < 0\n");
      test = false;
    }
    if (nbin_ch <= 0)
    {
      printf("Warning: nbin_ch < 0\n");
      test = false;
    }
  }


  return test;
}


