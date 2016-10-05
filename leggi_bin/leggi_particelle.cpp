
#include "leggi_particelle.h"


int read_phase_space_file(Parameters * params)
{
  size_t phasespace_size = 14; // x, y, z, px, py, pz, gamma, theta, thetaT, E, ty, tz, w, ch
  std::FILE *file_in = NULL;
  int multifile_index = 0;
  int contatori[] = { 0, 0, 0 };
  aladyn_float zero = 0.0f;
  long long parts_accumulate = 0;
  double peso_accumulato = 0.0;
  double carica_accumulata = 0.0;
  size_t dim_file_in_bytes = 0, num_of_floats_in_file = 0, num_of_particles_in_file = 0, num_of_passes = 0, num_residual_particles = 0, dimensione_array_parts = 0;
  unsigned int val[2] = { 0, 0 };
  aladyn_float array_supporto8[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  aladyn_float array_supporto6[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  size_t max_number_of_particles_in_memory = MAX_NUM_OF_PARTICLES_IN_MEMORY;

  std::FILE *binary_vtk = NULL;
  std::FILE *binary_clean = NULL;
  std::FILE *ascii_propaga = NULL;
  std::FILE *ascii_xyze = NULL;
  std::FILE *ascii_csv = NULL;
  std::FILE *parameters = NULL;
  std::ofstream Estremi_out;
  int conta_processori = 0;

  int npart_loc = 0;
  int buff = 0;

  aladyn_float x = 0.0, y = 0.0, z = 0.0, px = 0.0, py = 0.0, pz = 0.0, ptot = 0.0;
  aladyn_float ch = 0.0, w = 0.0;
  aladyn_float gamma = 0.0, theta = 0.0, thetaT = 0.0, E = 0.0, ty = 0.0, tz = 0.0;
  double *estremi_min = NULL, *estremi_max = NULL;

  short buffshort[2] = { 0, 0 };
  aladyn_float *parts = NULL;
  std::string bin_filename;
  std::string dat_filename;
  std::string minmax_filename;
  std::string vtk_filename;
  std::string clean_bin_filename;
  std::string ppg_filename;
  std::string xyze_filename;
  std::string csv_filename;
  std::string params_filename;
  std::string binned_filename;

  size_t fread_size = 0;
  double emittance_x = 0.0, emittance_y = 0.0, emittance_z = 0.0;
  double em_x2 = 0.0, em_x = 0.0, em_y2 = 0.0, em_y = 0.0, em_z2 = 0.0, em_z = 0.0;
  double em_px2 = 0.0, em_px = 0.0, em_py2 = 0.0, em_py = 0.0, em_pz2 = 0.0, em_pz = 0.0, em_xpx = 0.0, em_ypy = 0.0, em_zpz = 0.0;

  estremi_min = new double[phasespace_size];
  estremi_max = new double[phasespace_size];
  for (size_t i = 0; i < phasespace_size; i++)
  {
    estremi_min[i] = std::numeric_limits<double>::max();
    estremi_max[i] = std::numeric_limits<double>::lowest();
  }

  aladyn_float **xw = new aladyn_float*[params->nbin_x + 3];
  for (int i = 0; i < params->nbin_x + 3; i++)
  {
    xw[i] = new aladyn_float[params->nbin_w + 3];
    for (int j = 0; j < params->nbin_w + 3; j++) xw[i][j] = 0.0;
  }

  aladyn_float **xy = new aladyn_float*[params->nbin_x + 3];
  for (int i = 0; i < params->nbin_x + 3; i++)
  {
    xy[i] = new aladyn_float[params->nbin_y + 3];
    for (int j = 0; j < params->nbin_y + 3; j++) xy[i][j] = 0.0;
  }

  aladyn_float **xz = new aladyn_float*[params->nbin_x + 3];
  for (int i = 0; i < params->nbin_x + 3; i++)
  {
    xz[i] = new aladyn_float[params->nbin_z + 3];
    for (int j = 0; j < params->nbin_z + 3; j++) xz[i][j] = 0.0;
  }

  aladyn_float **yz = new aladyn_float*[params->nbin_y + 3];
  for (int i = 0; i < params->nbin_y + 3; i++)
  {
    yz[i] = new aladyn_float[params->nbin_z + 3];
    for (int j = 0; j < params->nbin_z + 3; j++) yz[i][j] = 0.0;
  }

  aladyn_float **rcf = new aladyn_float*[params->nbin_ty + 3];
  for (int i = 0; i < params->nbin_ty + 3; i++)
  {
    rcf[i] = new aladyn_float[params->nbin_tz + 3];
    for (int j = 0; j < params->nbin_tz + 3; j++) rcf[i][j] = 0.0;
  }

  aladyn_float **xpx = new aladyn_float*[params->nbin_x + 3];
  for (int i = 0; i < params->nbin_x + 3; i++)
  {
    xpx[i] = new aladyn_float[params->nbin_px + 3];
    for (int j = 0; j < params->nbin_px + 3; j++) xpx[i][j] = 0.0;
  }

  aladyn_float **xpy = new aladyn_float*[params->nbin_x + 3];
  for (int i = 0; i < params->nbin_x + 3; i++)
  {
    xpy[i] = new aladyn_float[params->nbin_py + 3];
    for (int j = 0; j < params->nbin_py + 3; j++) xpy[i][j] = 0.0;
  }

  aladyn_float **xpz = new aladyn_float*[params->nbin_x + 3];
  for (int i = 0; i < params->nbin_x + 3; i++)
  {
    xpz[i] = new aladyn_float[params->nbin_pz + 3];
    for (int j = 0; j < params->nbin_pz + 3; j++) xpz[i][j] = 0.0;
  }

  aladyn_float **ypx = new aladyn_float*[params->nbin_y + 3];
  for (int i = 0; i < params->nbin_y + 3; i++)
  {
    ypx[i] = new aladyn_float[params->nbin_px + 3];
    for (int j = 0; j < params->nbin_px + 3; j++) ypx[i][j] = 0.0;
  }

  aladyn_float **ypy = new aladyn_float*[params->nbin_y + 3];
  for (int i = 0; i < params->nbin_y + 3; i++)
  {
    ypy[i] = new aladyn_float[params->nbin_py + 3];
    for (int j = 0; j < params->nbin_py + 3; j++) ypy[i][j] = 0.0;
  }

  aladyn_float **ypz = new aladyn_float*[params->nbin_y + 3];
  for (int i = 0; i < params->nbin_y + 3; i++)
  {
    ypz[i] = new aladyn_float[params->nbin_pz + 3];
    for (int j = 0; j < params->nbin_pz + 3; j++) ypz[i][j] = 0.0;
  }

  aladyn_float **zpx = new aladyn_float*[params->nbin_z + 3];
  for (int i = 0; i < params->nbin_z + 3; i++)
  {
    zpx[i] = new aladyn_float[params->nbin_px + 3];
    for (int j = 0; j < params->nbin_px + 3; j++) zpx[i][j] = 0.0;
  }

  aladyn_float **zpy = new aladyn_float*[params->nbin_z + 3];
  for (int i = 0; i < params->nbin_z + 3; i++)
  {
    zpy[i] = new aladyn_float[params->nbin_py + 3];
    for (int j = 0; j < params->nbin_py + 3; j++) zpy[i][j] = 0.0;
  }

  aladyn_float **zpz = new aladyn_float*[params->nbin_z + 3];
  for (int i = 0; i < params->nbin_z + 3; i++)
  {
    zpz[i] = new aladyn_float[params->nbin_pz + 3];
    for (int j = 0; j < params->nbin_pz + 3; j++) zpz[i][j] = 0.0;
  }

  aladyn_float **pxpy = new aladyn_float*[params->nbin_px + 3];
  for (int i = 0; i < params->nbin_px + 3; i++)
  {
    pxpy[i] = new aladyn_float[params->nbin_py + 3];
    for (int j = 0; j < params->nbin_py + 3; j++) pxpy[i][j] = 0.0;
  }

  aladyn_float **pxpz = new aladyn_float*[params->nbin_px + 3];
  for (int i = 0; i < params->nbin_px + 3; i++)
  {
    pxpz[i] = new aladyn_float[params->nbin_pz + 3];
    for (int j = 0; j < params->nbin_pz + 3; j++) pxpz[i][j] = 0.0;
  }

  aladyn_float **pypz = new aladyn_float*[params->nbin_py + 3];
  for (int i = 0; i < params->nbin_py + 3; i++)
  {
    pypz[i] = new aladyn_float[params->nbin_pz + 3];
    for (int j = 0; j < params->nbin_pz + 3; j++) pypz[i][j] = 0.0;
  }

  aladyn_float **Etheta = new aladyn_float*[params->nbin_E + 3];
  for (int i = 0; i < params->nbin_E + 3; i++)
  {
    Etheta[i] = new aladyn_float[params->nbin_theta + 3];
    for (int j = 0; j < params->nbin_theta + 3; j++) Etheta[i][j] = 0.0;
  }

  aladyn_float **EthetaT = new aladyn_float*[params->nbin_E + 3];
  for (int i = 0; i < params->nbin_E + 3; i++)
  {
    EthetaT[i] = new aladyn_float[params->nbin_thetaT + 3];
    for (int j = 0; j < params->nbin_thetaT + 3; j++) EthetaT[i][j] = 0.0;
  }

  aladyn_float *wspec = new aladyn_float[params->nbin_w + 3];
  for (int i = 0; i < params->nbin_w + 3; i++) wspec[i] = 0.0;

  aladyn_float *chspec = new aladyn_float[params->nbin_ch + 3];
  for (int i = 0; i < params->nbin_ch + 3; i++) chspec[i] = 0.0;

  aladyn_float *Espec = new aladyn_float[params->nbin_E + 3];
  for (int i = 0; i < params->nbin_E + 3; i++) Espec[i] = 0.0;

  aladyn_float *thetaspec = new aladyn_float[params->nbin_theta + 3];
  for (int i = 0; i < params->nbin_theta + 3; i++) thetaspec[i] = 0.0;

  aladyn_float *thetaTspec = new aladyn_float[params->nbin_thetaT + 3];
  for (int i = 0; i < params->nbin_thetaT + 3; i++) thetaTspec[i] = 0.0;


  ppg_filename = params->filebasename + ".ppg";
  xyze_filename = params->filebasename + "_xyzE.ppg";
  csv_filename = params->filebasename + ".csv";
  vtk_filename = params->filebasename + ".vtk";
  clean_bin_filename = params->filebasename + "_clean.bin";
  bin_filename = params->filebasename + ".bin";
  dat_filename = params->filebasename + ".dat";
  minmax_filename = params->filebasename + ".extremes";
  params_filename = params->filebasename + ".parameters";

  if (!params->multifile) file_in = fopen(bin_filename.c_str(), "rb");


#ifdef ENABLE_DEBUG
  std::cout << "NPTOT = " << params->nptot << std::endl;
  if (params->we_have_to_do_binning)
  {
    std::cout << "XMIN = " << params->xmin << std::endl;
    std::cout << "XMAX = " << params->xmax << std::endl;
    std::cout << "YMIN = " << params->ymin << std::endl;
    std::cout << "YMAX = " << params->ymax << std::endl;
    std::cout << "ZMIN = " << params->zmin << std::endl;
    std::cout << "ZMAX = " << params->zmax << std::endl;
    std::cout << "PXMIN = " << params->pxmin << std::endl;
    std::cout << "PXMAX = " << params->pxmax << std::endl;
    std::cout << "PYMIN = " << params->pymin << std::endl;
    std::cout << "PYMAX = " << params->pymax << std::endl;
    std::cout << "PZMIN = " << params->pzmin << std::endl;
    std::cout << "PZMAX = " << params->pzmax << std::endl;
    std::cout << "TYMIN = " << params->tymin << std::endl;
    std::cout << "TYMAX = " << params->tymax << std::endl;
    std::cout << "TZMIN = " << params->tzmin << std::endl;
    std::cout << "TZMAX = " << params->tzmax << std::endl;
    std::cout << "GAMMAMIN = " << params->gammamin << std::endl;
    std::cout << "GAMMAMAX = " << params->gammamax << std::endl;
    std::cout << "THETAMIN = " << params->thetamin << std::endl;
    std::cout << "THETAMAX = " << params->thetamax << std::endl;
    std::cout << "THETARADMIN = " << params->thetaTmin << std::endl;
    std::cout << "THETARADMAX = " << params->thetaTmax << std::endl;
    std::cout << "EMIN = " << params->Emin << std::endl;
    std::cout << "EMAX = " << params->Emax << std::endl;
    std::cout << "WMIN = " << params->wmin << std::endl;
    std::cout << "WMAX = " << params->wmax << std::endl;
    std::cout << "CHMIN = " << params->chmin << std::endl;
    std::cout << "CHMAX = " << params->chmax << std::endl;
  }
#endif

  if (params->out_vtk)
  {
    printf("\nENABLED .vtk FILE\n");
    binary_vtk = fopen(vtk_filename.c_str(), "wb");

    // Scrittura primo Header VTK e memorizzazione sua dimensione in contatori[0]
    contatori[0] += fprintf(binary_vtk, "# vtk DataFile Version 2.0\n");
    contatori[0] += fprintf(binary_vtk, "titolo nostro\n");
    contatori[0] += fprintf(binary_vtk, "BINARY\n");
    contatori[0] += fprintf(binary_vtk, "DATASET UNSTRUCTURED_GRID\n");
    contatori[0] += fprintf(binary_vtk, "POINTS %llu aladyn_float\n", (unsigned long long int) params->nptot);

    fseeko(binary_vtk, contatori[0] + params->nptot * sizeof(aladyn_float) * 3, SEEK_SET);

    //Scrittura secondo Header VTK e memorizzazione sua dimensione in contatori[1]
    //contatori[1] += fprintf(binary_vtk, "DATASET UNSTRUCTURED_GRID\n");
    contatori[1] += fprintf(binary_vtk, "POINT_DATA %llu\n", (unsigned long long int) params->nptot);
    //contatori[1] += fprintf(binary_vtk, "POINTS %i aladyn_float\n", params->nptot);
    contatori[1] += fprintf(binary_vtk, "VECTORS p aladyn_float\n");
    //contatori[1] += fprintf(binary_vtk, "LOOKUP_TABLE default\n");

    if (params->file_has_weight)
    {
      fseeko(binary_vtk, contatori[0] + params->nptot * sizeof(aladyn_float) * 3 + contatori[1] + params->nptot * sizeof(aladyn_float) * 3, SEEK_SET);

      //Scrittura terzo Header VTK e memorizzazione sua dimensione in contatori[2]
      //contatori[2] += fprintf(binary_vtk,"DATASET STRUCTURED_POINTS\n");
      //contatori[2] += fprintf(binary_vtk,"DIMENSIONS %ui %i %i\n",params->nptot, 1, 1);
      //contatori[2] += fprintf(binary_vtk,"ORIGIN 0 0 0\n");
      //contatori[2] += fprintf(binary_vtk,"SPACING 1 1 1\n");
      //contatori[2] += fprintf(binary_vtk,"POINT_DATA %ui\n",params->nptot);
      contatori[2] += fprintf(binary_vtk, "SCALARS w aladyn_float 1\n");
      contatori[2] += fprintf(binary_vtk, "LOOKUP_TABLE default\n");
    }
  }

  if (params->out_clean_bin)
  {
    printf("\nENABLED CLEAN .bin FILE\n");
    binary_clean = fopen(clean_bin_filename.c_str(), "wb");
  }

  if (params->out_ppg)
  {
    printf("\nENABLED .txt FILE FOR PROPAGA\n");
    ascii_propaga = fopen(ppg_filename.c_str(), "wb");
  }

  if (params->out_xyze)
  {
    printf("\nENABLED .txt FILE WITH x, y, z, E\n");
    ascii_xyze = fopen(xyze_filename.c_str(), "wb");
  }

  if (params->out_csv)
  {
    printf("\nENABLED .csv FILE FOR PARAVIEW\n");
    ascii_csv = fopen(csv_filename.c_str(), "wb");
  }

  fflush(stdout);

  while (1)
  {
    if (!params->multifile)
    {
      if (conta_processori >= params->last_cpu) break;
      if (conta_processori == 0) {
        /*skip header*/
        std::fseek(file_in, (long)params->header_size_bytes, SEEK_SET);
      }

      if (params->file_version == 1)
      {
        fread_size = std::fread(&buff, sizeof(int), 1, file_in);
        fread_size = std::fread(&npart_loc, sizeof(int), 1, file_in);
        fread_size = std::fread(&buff, sizeof(int), 1, file_in);
      }
      else
      {
        fread_size = std::fread(&npart_loc, sizeof(int), 1, file_in);
      }

      if (feof(file_in)) break;
      if (params->we_have_to_do_swap) swap_endian_i(&npart_loc, 1);
      if (npart_loc > (long long int) params->nptot || npart_loc < 0)
      {
        printf("Read a npart=%i, non valid. Exiting!", npart_loc);
        break;
      }
      dimensione_array_parts = npart_loc;
      val[0] = (unsigned int)npart_loc;
      val[1] = (unsigned int)params->ndv;
#ifdef ENABLE_DEBUG
      printf("proc number \t %i \t npart=%i \n", conta_processori, npart_loc);
#else
      printf("proc number \t %i \t npart=%i \r", conta_processori, npart_loc);
#endif
      fflush(stdout);
      num_of_passes = 1;
    }
    else  //we do have multifiles i.e. Prpout00_000.bin
    {
      memset(&bin_filename[0], 0, sizeof(bin_filename));
      bin_filename = params->filebasename + "_" + std::to_string(multifile_index) + ".bin";

      if ((file_in = fopen(bin_filename.c_str(), "rb")) == NULL)
      {
        printf("End of files! \n");
        break;
      }
      fseeko(file_in, 0, SEEK_END);
      dim_file_in_bytes = (int)ftello(file_in);
      rewind(file_in);
      num_of_floats_in_file = (dim_file_in_bytes / sizeof(aladyn_float));
      num_of_particles_in_file = (int)(num_of_floats_in_file / params->ndv);
      printf("File %s_%.3i.bin has %llu particles\n", params->filebasename.c_str(), multifile_index, (unsigned long long int) num_of_particles_in_file);
      fflush(stdout);
      num_of_passes = (int)((aladyn_float)(num_of_particles_in_file) / (aladyn_float)(max_number_of_particles_in_memory)) + 1;
      num_residual_particles = num_of_particles_in_file % max_number_of_particles_in_memory;
      dimensione_array_parts = std::min(max_number_of_particles_in_memory, num_of_particles_in_file);
      if (dimensione_array_parts > params->nptot)
      {
        printf("Read npart=%llu: not valid! Exiting...\n", (unsigned long long int) dimensione_array_parts);
        break;
      }
      val[0] = (unsigned int)dimensione_array_parts;
      val[1] = (unsigned int)params->ndv;
    }

    if (val[0] > 0)
    {
      fflush(stdout);
      for (size_t h = 0; h < num_of_passes; h++)
      {
        if (num_of_passes > 1) printf("File is very big, will be splitted in multiple readings: step %llu of %llu\n", (unsigned long long int) (h + 1), (unsigned long long int) num_of_passes);
        if (!params->multifile)
        {
          parts = new aladyn_float[npart_loc*params->ndv];
          if (params->file_version == 1)
          {
            fread_size = std::fread(buffshort, sizeof(short), 2, file_in);
            fread_size = std::fread(parts, sizeof(aladyn_float), npart_loc*params->ndv, file_in);
            fread_size = std::fread(&buff, sizeof(int), 1, file_in);
          }
          else fread_size = std::fread(parts, sizeof(aladyn_float), npart_loc*params->ndv, file_in);
          if (params->we_have_to_do_swap) swap_endian_f(parts, (size_t)npart_loc*params->ndv);
        }
        else
        {
          if (h == num_of_passes - 1 && num_of_passes > 1) dimensione_array_parts = num_residual_particles;
          parts = new aladyn_float[dimensione_array_parts*params->ndv];
          val[0] = (unsigned int)dimensione_array_parts;
#ifdef ENABLE_DEBUG
          printf("npart_loc = %i\t\t ndv=%i\n", val[0], params->ndv);
#else
          printf("npart_loc = %i\t\t ndv=%i\r", val[0], params->ndv);
#endif
          fflush(stdout);
          fread_size = std::fread(parts, sizeof(aladyn_float), val[0] * params->ndv, file_in);
          if (params->we_have_to_do_swap) swap_endian_f(parts, (size_t)val[0] * params->ndv);
        }

        _Filter(params, parts, val, _Filter::build_filter(params));

        if (params->out_params)
        {
          for (unsigned int i = 0; i < val[0]; i++)
          {
            if (params->sim_is_2d)
            {
              x = *(parts + i*params->ndv);
              y = *(parts + i*params->ndv + 1);
              z = 0.0;
              px = *(parts + i*params->ndv + 2);
              py = *(parts + i*params->ndv + 3);
              pz = 0.0;
              ptot = sqrt(px*px + py*py);
              px /= ptot;
              py /= ptot;
              if (params->file_has_weight && !params->overwrite_weight)
                w = *(parts + i*params->ndv + 4);
              else
                w = params->overwrite_weight_value;
              if (params->file_version >= 3 && !params->overwrite_charge)
                ch = *(parts + i*params->ndv + 5);
              else
                ch = params->overwrite_charge_value;
            }
            else
            {
              x = *(parts + i*params->ndv);
              y = *(parts + i*params->ndv + 1);
              z = *(parts + i*params->ndv + 2);
              px = *(parts + i*params->ndv + 3);
              py = *(parts + i*params->ndv + 4);
              pz = *(parts + i*params->ndv + 5);
              ptot = sqrt(px*px + py*py + pz*pz);
              px /= ptot;
              py /= ptot;
              pz /= ptot;
              if (params->file_has_weight && !params->overwrite_weight)
                w = *(parts + i*params->ndv + 6);
              else
                w = params->overwrite_weight_value;
              if (params->file_version >= 3 && !params->overwrite_charge)
                ch = *(parts + i*params->ndv + 7);
              else
                ch = params->overwrite_charge_value;
            }

            em_x += (double)(x*w);
            em_y += (double)(y*w);
            em_z += (double)(z*w);
            em_x2 += (double)((x*x)*w);
            em_y2 += (double)((y*y)*w);
            em_z2 += (double)((z*z)*w);
            em_px += (double)(px*w);
            em_py += (double)(py*w);
            em_pz += (double)(pz*w);
            em_px2 += (double)((px*px)*w);
            em_py2 += (double)((py*py)*w);
            em_pz2 += (double)((pz*pz)*w);
            em_xpx += (double)((x*px)*w);
            em_ypy += (double)((y*py)*w);
            em_zpz += (double)((z*pz)*w);
            peso_accumulato += w;
            carica_accumulata += ch;
          }
        }

        if (params->we_have_to_find_minmax)
        {
          for (unsigned int i = 0; i < val[0]; i++)
          {
            if (params->sim_is_2d)
            {
              x = *(parts + i*params->ndv);
              y = *(parts + i*params->ndv + 1);
              px = *(parts + i*params->ndv + 2);
              py = *(parts + i*params->ndv + 3);
              if (params->file_has_weight && !params->overwrite_weight)
                w = *(parts + i*params->ndv + 4);
              else
                w = params->overwrite_weight_value;
              if (params->file_version >= 3 && !params->overwrite_charge)
                ch = *(parts + i*params->ndv + 5);
              else
                ch = params->overwrite_charge_value;
              gamma = (aladyn_float)(sqrt(1. + px*px + py*py) - 1.);       //gamma
              theta = (aladyn_float)(atan2(py, px)*180. / M_PI);       //theta
              thetaT = (aladyn_float)atan(sqrt((py*py) / (px*px)));
              E = (aladyn_float)(gamma*params->mass_MeV); //energia
              if (px > 0) ty = py / px;
              else ty = 0.0;
              if (x < estremi_min[0]) estremi_min[0] = x;
              if (x > estremi_max[0]) estremi_max[0] = x;
              if (y < estremi_min[1]) estremi_min[1] = y;
              if (y > estremi_max[1]) estremi_max[1] = y;
              estremi_min[2] = 0., estremi_max[2] = 1.;
              if (px < estremi_min[3]) estremi_min[3] = px;
              if (px > estremi_max[3]) estremi_max[3] = px;
              if (py < estremi_min[4]) estremi_min[4] = py;
              if (py > estremi_max[4]) estremi_max[4] = py;
              estremi_min[5] = 0., estremi_max[5] = 1.;
              if (gamma < estremi_min[6]) estremi_min[6] = gamma;
              if (gamma > estremi_max[6]) estremi_max[6] = gamma;
              if (theta < estremi_min[7]) estremi_min[7] = theta;
              if (theta > estremi_max[7]) estremi_max[7] = theta;
              if (E < estremi_min[8]) estremi_min[8] = E;
              if (E > estremi_max[8]) estremi_max[8] = E;
              if (thetaT < estremi_min[9]) estremi_min[9] = thetaT;
              if (thetaT > estremi_max[9]) estremi_max[9] = thetaT;
              if (ty < estremi_min[10]) estremi_min[10] = ty;
              if (ty > estremi_max[10]) estremi_max[10] = ty;
              estremi_min[11] = -1., estremi_max[11] = 1.;
              if (w < estremi_min[12]) estremi_min[12] = w;
              if (w > estremi_max[12]) estremi_max[12] = w;
              if (ch < estremi_min[13]) estremi_min[13] = ch;
              if (ch > estremi_max[13]) estremi_max[13] = ch;
            }
            else
            {
              x = *(parts + i*params->ndv);
              y = *(parts + i*params->ndv + 1);
              z = *(parts + i*params->ndv + 2);
              px = *(parts + i*params->ndv + 3);
              py = *(parts + i*params->ndv + 4);
              pz = *(parts + i*params->ndv + 5);
              if (params->file_has_weight && !params->overwrite_weight)
                w = *(parts + i*params->ndv + 6);
              else
                w = params->overwrite_weight_value;
              if (params->file_version >= 3 && !params->overwrite_charge)
                ch = *(parts + i*params->ndv + 7);
              else
                ch = params->overwrite_charge_value;
              gamma = (aladyn_float)(sqrt(1. + px*px + py*py + pz*pz) - 1.);     //gamma
              theta = (aladyn_float)(atan2(sqrt(py*py + pz*pz), px)*180. / M_PI);  //theta nb: py e pz sono quelli trasversi in ALaDyn!
              thetaT = (aladyn_float)atan(sqrt((py*py + pz*pz) / (px*px)));
              E = (aladyn_float)(gamma*params->mass_MeV);   //energia
              if (px > 0)
              {
                ty = py / px;
                tz = pz / px;
              }
              else
              {
                ty = 0.0;
                tz = 0.0;
              }
              if (x < estremi_min[0]) estremi_min[0] = x;
              if (x > estremi_max[0]) estremi_max[0] = x;
              if (y < estremi_min[1]) estremi_min[1] = y;
              if (y > estremi_max[1]) estremi_max[1] = y;
              if (z < estremi_min[2]) estremi_min[2] = z;
              if (z > estremi_max[2]) estremi_max[2] = z;
              if (px < estremi_min[3]) estremi_min[3] = px;
              if (px > estremi_max[3]) estremi_max[3] = px;
              if (py < estremi_min[4]) estremi_min[4] = py;
              if (py > estremi_max[4]) estremi_max[4] = py;
              if (pz < estremi_min[5]) estremi_min[5] = pz;
              if (pz > estremi_max[5]) estremi_max[5] = pz;
              if (gamma < estremi_min[6]) estremi_min[6] = gamma;
              if (gamma > estremi_max[6]) estremi_max[6] = gamma;
              if (theta < estremi_min[7]) estremi_min[7] = theta;
              if (theta > estremi_max[7]) estremi_max[7] = theta;
              if (E < estremi_min[8]) estremi_min[8] = E;
              if (E > estremi_max[8]) estremi_max[8] = E;
              if (thetaT < estremi_min[9]) estremi_min[9] = thetaT;
              if (thetaT > estremi_max[9]) estremi_max[9] = thetaT;
              if (ty < estremi_min[10]) estremi_min[10] = ty;
              if (ty > estremi_max[10]) estremi_max[10] = ty;
              if (tz < estremi_min[11]) estremi_min[11] = tz;
              if (tz > estremi_max[11]) estremi_max[11] = tz;
              if (w < estremi_min[12]) estremi_min[12] = w;
              if (w > estremi_max[12]) estremi_max[12] = w;
              if (ch < estremi_min[13]) estremi_min[13] = ch;
              if (ch > estremi_max[13]) estremi_max[13] = ch;
            }
          }
        }
        if (params->we_have_to_do_binning)
        {
          if (params->do_plot_xy)         _Binning(parts, params, xy, "x", "y");
          if (params->do_plot_xz)         _Binning(parts, params, xz, "x", "z");
          if (params->do_plot_yz)         _Binning(parts, params, yz, "y", "z");
          if (params->do_plot_rcf)        _Binning(parts, params, rcf, "ty", "tz");
          if (params->do_plot_xpx)        _Binning(parts, params, xpx, "x", "px");
          if (params->do_plot_xpy)        _Binning(parts, params, xpy, "x", "py");
          if (params->do_plot_xpz)        _Binning(parts, params, xpz, "x", "pz");
          if (params->do_plot_ypx)        _Binning(parts, params, ypx, "y", "px");
          if (params->do_plot_ypy)        _Binning(parts, params, ypy, "y", "py");
          if (params->do_plot_ypz)        _Binning(parts, params, ypz, "y", "pz");
          if (params->do_plot_zpx)        _Binning(parts, params, zpx, "z", "px");
          if (params->do_plot_zpy)        _Binning(parts, params, zpy, "z", "py");
          if (params->do_plot_zpz)        _Binning(parts, params, zpz, "z", "pz");
          if (params->do_plot_pxpy)       _Binning(parts, params, pxpy, "px", "py");
          if (params->do_plot_pxpz)       _Binning(parts, params, pxpz, "px", "pz");
          if (params->do_plot_pypz)       _Binning(parts, params, pypz, "py", "pz");
          if (params->do_plot_xw)         _Binning(parts, params, xw, "x", "w");
          if (params->do_plot_Etheta)     _Binning(parts, params, Etheta, "E", "theta");
          if (params->do_plot_EthetaT)    _Binning(parts, params, EthetaT, "E", "thetaT");
          if (params->do_plot_Espec)      _Binning(parts, params, Espec, "E");
          if (params->do_plot_wspec)      _Binning(parts, params, wspec, "w");
          if (params->do_plot_chspec)     _Binning(parts, params, chspec, "ch");
          if (params->do_plot_thetaspec)  _Binning(parts, params, thetaspec, "theta");
          if (params->do_plot_thetaTspec) _Binning(parts, params, thetaTspec, "thetaT");
        }

        if (params->out_ppg)
        {
          if (params->sim_is_2d)
          {
            for (unsigned int i = 0; i < val[0]; i++)
            {
              x = parts[i*params->ndv + 1] * ((aladyn_float)1.e-4);
              z = parts[i*params->ndv + 0] * ((aladyn_float)1.e-4);
              px = parts[i*params->ndv + 3];
              pz = parts[i*params->ndv + 2];
              if (params->file_has_weight && !params->overwrite_weight)
                w = parts[i*params->ndv + 4];
              else
                w = params->overwrite_weight_value;

              if (i % params->subsample == 0) {
                fprintf(ascii_propaga, "%e 0 %e %e 0 %e %d %e 0 %d\n", x, z, px, pz, 1, w, i + 1); // fix: no distinction in the particle type column, always written as protons!
              }
            }
          }
          else
          {
            for (unsigned int i = 0; i < val[0]; i++)
            {
              x = parts[i*params->ndv + 1] * ((aladyn_float)1.e-4);
              y = parts[i*params->ndv + 2] * ((aladyn_float)1.e-4);
              z = parts[i*params->ndv + 0] * ((aladyn_float)1.e-4);
              px = parts[i*params->ndv + 4];
              py = parts[i*params->ndv + 5];
              pz = parts[i*params->ndv + 3];
              if (params->file_has_weight && !params->overwrite_weight)
                w = parts[i*params->ndv + 6];
              else
                w = params->overwrite_weight_value;

              if (i % params->subsample == 0) {
                fprintf(ascii_propaga, "%e %e %e %e %e %e %d %e 0 %d\n", x, y, z, px, py, pz, 1, w, i + 1);  // fix: no distinction in the particle type column, always written as protons!
              }
            }
          }
        }



        if (params->out_xyze)
        {
          if (params->sim_is_2d)
          {
            for (unsigned int i = 0; i < val[0]; i++)
            {
              x = parts[i*params->ndv + 1] * ((aladyn_float)1.e-4);
              z = parts[i*params->ndv + 0] * ((aladyn_float)1.e-4);
              px = parts[i*params->ndv + 3];
              pz = parts[i*params->ndv + 2];
              gamma = (aladyn_float)(sqrt(1. + px*px + py*py) - 1.);
              E = (aladyn_float)(gamma*params->mass_MeV);
              if (i % params->subsample == 0) fprintf(ascii_xyze, "%e %e %e\n", x, z, E);
            }
          }
          else
          {
            for (unsigned int i = 0; i < val[0]; i++)
            {
              x = parts[i*params->ndv + 1] * ((aladyn_float)1.e-4);
              y = parts[i*params->ndv + 2] * ((aladyn_float)1.e-4);
              z = parts[i*params->ndv + 0] * ((aladyn_float)1.e-4);
              px = parts[i*params->ndv + 4];
              py = parts[i*params->ndv + 5];
              pz = parts[i*params->ndv + 3];
              gamma = (aladyn_float)(sqrt(1. + px*px + py*py + pz*pz) - 1.);
              E = (aladyn_float)(gamma*params->mass_MeV);
              if (i % params->subsample == 0) fprintf(ascii_xyze, "%e %e %e %e\n", x, y, z, E);
            }
          }
        }



        if (params->out_csv)
        {
          for (unsigned int i = 0; i < val[0]; i++)
          {
            x = parts[i*params->ndv + 0];
            y = parts[i*params->ndv + 1];
            if (params->sim_is_2d) z = 0;
            else z = parts[i*params->ndv + 2];
            px = parts[i*params->ndv + 3];
            py = parts[i*params->ndv + 4];
            if (params->sim_is_2d) pz = 0;
            else pz = parts[i*params->ndv + 5];

            if (params->sim_is_2d) {
              if (params->file_has_weight && !params->overwrite_weight) w = parts[i*params->ndv + 4];
              else w = params->overwrite_weight_value;
            }
            else {
              if (params->file_has_weight && !params->overwrite_weight) w = parts[i*params->ndv + 6];
              else w = params->overwrite_weight_value;
            }

            if (params->sim_is_2d) {
              if (params->file_has_charge && !params->overwrite_charge) ch = parts[i*params->ndv + 5];
              else ch = params->overwrite_charge_value;
            }
            else {
              if (params->file_has_charge && !params->overwrite_charge) ch = parts[i*params->ndv + 7];
              else ch = params->overwrite_charge_value;
            }


            if (i % params->subsample == 0) fprintf(ascii_csv, "%e, %e, %e, %e, %e, %e, %e, %e\n", x, y, z, px, py, pz, w, ch);
          }
        }



        if (params->out_clean_bin)
        {
          for (unsigned int i = 0; i < val[0]; i++)
          {
            array_supporto8[0] = parts[i*params->ndv + 0];
            array_supporto8[1] = parts[i*params->ndv + 1];
            if (params->sim_is_2d) array_supporto8[2] = 0;
            else array_supporto8[2] = parts[i*params->ndv + 2];
            array_supporto8[3] = parts[i*params->ndv + 3];
            array_supporto8[4] = parts[i*params->ndv + 4];
            if (params->sim_is_2d) array_supporto8[5] = 0;
            else array_supporto8[5] = parts[i*params->ndv + 5];

            if (params->sim_is_2d) {
              if (params->file_has_weight && !params->overwrite_weight) array_supporto8[6] = parts[i*params->ndv + 4];
              else array_supporto8[6] = params->overwrite_weight_value;
            }
            else {
              if (params->file_has_weight && !params->overwrite_weight) array_supporto8[6] = parts[i*params->ndv + 6];
              else array_supporto8[6] = params->overwrite_weight_value;
            }

            if (params->sim_is_2d) {
              if (params->file_has_charge && !params->overwrite_charge) array_supporto8[7] = parts[i*params->ndv + 5];
              else array_supporto8[7] = params->overwrite_charge_value;
            }
            else {
              if (params->file_has_charge && !params->overwrite_charge) array_supporto8[7] = parts[i*params->ndv + 7];
              else array_supporto8[7] = params->overwrite_charge_value;
            }

            if (i % params->subsample == 0) fwrite((void*)(array_supporto8), sizeof(aladyn_float), 8, binary_clean);
          }
        }


        if (params->out_vtk)
        {
          if (params->endian_machine == 0)
          {
            swap_endian_f(parts, (size_t)val[0] * params->ndv);
          }
          int weight_index, charge_index;

          if (params->sim_is_2d)
          {
            weight_index = 4;
            if (params->file_has_charge) charge_index = 5;
            else charge_index = 0;
            // scrittura coordinate x, y, z
            fseeko(binary_vtk, contatori[0] + parts_accumulate * sizeof(aladyn_float) * 3, SEEK_SET);
            for (unsigned int i = 0; i < val[0]; i += params->ndv)
            {
              fwrite((void*)(parts + i), sizeof(aladyn_float), 2, binary_vtk);
              fwrite((void*)&zero, sizeof(aladyn_float), 1, binary_vtk);
            }

            // scrittura momenti px, py, pz
            fseeko(binary_vtk, contatori[0] + params->nptot * sizeof(aladyn_float) * 3 + contatori[1] + parts_accumulate * sizeof(aladyn_float) * 3, SEEK_SET);
            for (unsigned int i = 2; i < val[0]; i += params->ndv)
            {
              fwrite((void*)(parts + i), sizeof(aladyn_float), 2, binary_vtk);
              fwrite((void*)&zero, sizeof(aladyn_float), 1, binary_vtk);
            }
          }
          else
          {
            weight_index = 6;
            if (params->file_has_charge) charge_index = 7;
            else charge_index = 0;
            // scrittura coordinate x, y, z
            fseeko(binary_vtk, contatori[0] + parts_accumulate * sizeof(aladyn_float) * 3, SEEK_SET);
            for (unsigned int i = 0; i < val[0]; i += params->ndv)
              fwrite((void*)(parts + i), sizeof(aladyn_float), 3, binary_vtk);

            // scrittura momenti px, py, pz
            fseeko(binary_vtk, contatori[0] + params->nptot * sizeof(aladyn_float) * 3 + contatori[1] + parts_accumulate * sizeof(aladyn_float) * 3, SEEK_SET);
            for (unsigned int i = 3; i < val[0]; i += params->ndv)
              fwrite((void*)(parts + i), sizeof(aladyn_float), 3, binary_vtk);
          }



          if (params->file_has_weight && !params->overwrite_weight)
          {
            // scrittura pesi
            fseeko(binary_vtk, contatori[0] + params->nptot * sizeof(aladyn_float) * 3 + contatori[1] + params->nptot * sizeof(aladyn_float) * 3 + contatori[2] + parts_accumulate * sizeof(aladyn_float), SEEK_SET);
            for (unsigned int i = weight_index; i < val[0]; i += params->ndv)
              fwrite((void*)(parts + i), sizeof(aladyn_float), 1, binary_vtk);
          }
          else if (params->file_has_weight && params->overwrite_weight)
          {
            // scrittura pesi sovrascritti da linea comando
            fseeko(binary_vtk, contatori[0] + params->nptot * sizeof(aladyn_float) * 3 + contatori[1] + params->nptot * sizeof(aladyn_float) * 3 + contatori[2] + parts_accumulate * sizeof(aladyn_float), SEEK_SET);
            for (unsigned int i = weight_index; i < val[0]; i += params->ndv)
              fwrite((void*)(&(params->overwrite_weight_value)), sizeof(aladyn_float), 1, binary_vtk);
          }

          if (params->file_has_charge && !params->overwrite_charge)
          {
            // scrittura pesi
            fseeko(binary_vtk, contatori[0] + params->nptot * sizeof(aladyn_float) * 3 + contatori[1] + params->nptot * sizeof(aladyn_float) * 3 + contatori[2] + parts_accumulate * sizeof(aladyn_float), SEEK_SET);
            for (unsigned int i = charge_index; i < val[0]; i += params->ndv)
              fwrite((void*)(parts + i), sizeof(aladyn_float), 1, binary_vtk);
          }
          else if (params->file_has_charge && params->overwrite_charge)
          {
            // scrittura pesi sovrascritti da linea comando
            fseeko(binary_vtk, contatori[0] + params->nptot * sizeof(aladyn_float) * 3 + contatori[1] + params->nptot * sizeof(aladyn_float) * 3 + contatori[2] + parts_accumulate * sizeof(aladyn_float), SEEK_SET);
            for (unsigned int i = charge_index; i < val[0]; i += params->ndv)
              fwrite((void*)(&(params->overwrite_weight_value)), sizeof(aladyn_float), 1, binary_vtk);
          }
        }
        parts_accumulate += val[0];
        delete[] parts;
        parts = NULL;
      }
    }
    multifile_index++;
    conta_processori++;
  }


  if (params->out_params)
  {
    em_x /= (double)peso_accumulato;
    em_y /= (double)peso_accumulato;
    em_z /= (double)peso_accumulato;
    em_px /= (double)peso_accumulato;
    em_py /= (double)peso_accumulato;
    em_pz /= (double)peso_accumulato;
    em_x2 /= (double)peso_accumulato;
    em_y2 /= (double)peso_accumulato;
    em_z2 /= (double)peso_accumulato;
    em_px2 /= (double)peso_accumulato;
    em_py2 /= (double)peso_accumulato;
    em_pz2 /= (double)peso_accumulato;
    em_xpx /= (double)peso_accumulato;
    em_ypy /= (double)peso_accumulato;
    em_zpz /= (double)peso_accumulato;

    emittance_x = sqrt((em_x2 - em_x*em_x)*(em_px2 - em_px*em_px) - (em_xpx - em_x*em_px)*(em_xpx - em_x*em_px));
    emittance_y = sqrt((em_y2 - em_y*em_y)*(em_py2 - em_py*em_py) - (em_ypy - em_y*em_py)*(em_ypy - em_y*em_py));
    emittance_z = sqrt((em_z2 - em_z*em_z)*(em_pz2 - em_pz*em_pz) - (em_zpz - em_z*em_pz)*(em_zpz - em_z*em_pz));

    parameters = fopen(params_filename.c_str(), "wb");
    printf("\nWriting the parameters file\n");
    fprintf(parameters, "ncpu_x=%i\n", params->ncpu_x);
    fprintf(parameters, "ncpu_y=%i\n", params->ncpu_y);
    fprintf(parameters, "ncpu_z=%i\n", params->ncpu_z);
    fprintf(parameters, "npx_resampled=%llu\n", (unsigned long long int) params->npx_resampled);
    fprintf(parameters, "npy_resampled=%llu\n", (unsigned long long int) params->npy_resampled);
    fprintf(parameters, "npz_resampled=%llu\n", (unsigned long long int) params->npz_resampled);
    fprintf(parameters, "npx_resampled_per_cpu=%llu\n", (unsigned long long int) params->npx_resampled_per_cpu);
    fprintf(parameters, "npy_resampled_per_cpu=%llu\n", (unsigned long long int) params->npy_resampled_per_cpu);
    fprintf(parameters, "npz_resampled_per_cpu=%llu\n", (unsigned long long int) params->npz_resampled_per_cpu);
    fprintf(parameters, "tnow=%f\n", params->tnow);
    fprintf(parameters, "xmin=%f\n", params->xmin);
    fprintf(parameters, "xmax=%f\n", params->xmax);
    fprintf(parameters, "ymin=%f\n", params->ymin);
    fprintf(parameters, "ymax=%f\n", params->ymax);
    fprintf(parameters, "zmin=%f\n", params->zmin);
    fprintf(parameters, "zmax=%f\n", params->zmax);
    ////////////////////////////////////////////////////////////////////////////
    fprintf(parameters, "rms_emittance_x=%g\n", emittance_x);    //massa parts su massa elettrone
    fprintf(parameters, "rms_emittance_y=%g\n", emittance_y);    //massa parts su massa elettrone
    fprintf(parameters, "rms_emittance_z=%g\n", emittance_z);    //massa parts su massa elettrone
    ////////////////////////////////////////////////////////////////////////////
    fclose(parameters);
  }

  if (params->we_have_to_find_minmax)
  {
    Estremi_out.open(minmax_filename);
    if (Estremi_out.fail()) printf("\nunable to create .extremes file");
    Estremi_out << "XMIN = " << estremi_min[0] << std::endl;
    Estremi_out << "XMAX = " << estremi_max[0] << std::endl;
    Estremi_out << "YMIN = " << estremi_min[1] << std::endl;
    Estremi_out << "YMAX = " << estremi_max[1] << std::endl;
    Estremi_out << "ZMIN = " << estremi_min[2] << std::endl;
    Estremi_out << "ZMAX = " << estremi_max[2] << std::endl;
    Estremi_out << "PXMIN = " << estremi_min[3] << std::endl;
    Estremi_out << "PXMAX = " << estremi_max[3] << std::endl;
    Estremi_out << "PYMIN = " << estremi_min[4] << std::endl;
    Estremi_out << "PYMAX = " << estremi_max[4] << std::endl;
    Estremi_out << "PZMIN = " << estremi_min[5] << std::endl;
    Estremi_out << "PZMAX = " << estremi_max[5] << std::endl;
    Estremi_out << "GAMMAMIN = " << estremi_min[6] << std::endl;
    Estremi_out << "GAMMAMAX = " << estremi_max[6] << std::endl;
    Estremi_out << "THETAMIN = " << estremi_min[7] << std::endl;
    Estremi_out << "THETAMAX = " << estremi_max[7] << std::endl;
    Estremi_out << "THETARADMIN = " << estremi_min[9] << std::endl;
    Estremi_out << "THETARADMAX = " << estremi_max[9] << std::endl;
    Estremi_out << "EMIN = " << estremi_min[8] << std::endl;
    Estremi_out << "EMAX = " << estremi_max[8] << std::endl;
    Estremi_out << "TYMIN = " << estremi_min[10] << std::endl;
    Estremi_out << "TYMAX = " << estremi_max[10] << std::endl;
    Estremi_out << "TZMIN = " << estremi_min[11] << std::endl;
    Estremi_out << "TZMAX = " << estremi_max[11] << std::endl;
    Estremi_out << "WMIN = " << estremi_min[12] << std::endl;
    Estremi_out << "WMAX = " << estremi_max[12] << std::endl;
    Estremi_out << "CHMIN = " << estremi_min[12] << std::endl;
    Estremi_out << "CHMAX = " << estremi_max[12] << std::endl;
    if (params->out_params) Estremi_out << "collected_weight = " << peso_accumulato << std::endl;
    if (params->out_params) Estremi_out << "collected_charge = " << carica_accumulata << std::endl;
    Estremi_out.close();
  }


  if (params->we_have_to_do_binning)
  {
    if (params->do_plot_xy)
    {
      binned_filename = params->filebasename + "_xy.txt";
      _Write(params, xy, "x", "y", std::string(binned_filename));
    }
    if (params->do_plot_xz)
    {
      binned_filename = params->filebasename + "_xz.txt";
      _Write(params, xz, "x", "z", std::string(binned_filename));
    }
    if (params->do_plot_yz)
    {
      binned_filename = params->filebasename + "_yz.txt";
      _Write(params, yz, "y", "z", std::string(binned_filename));
    }
    if (params->do_plot_rcf)
    {
      binned_filename = params->filebasename + "_rcf.txt";
      _Write(params, rcf, "ty", "tz", std::string(binned_filename));
    }
    if (params->do_plot_xpx)
    {
      binned_filename = params->filebasename + "_xpx.txt";
      _Write(params, xpx, "x", "px", std::string(binned_filename));
    }
    if (params->do_plot_xpy)
    {
      binned_filename = params->filebasename + "_xpy.txt";
      _Write(params, xpy, "x", "py", std::string(binned_filename));
    }
    if (params->do_plot_xpz)
    {
      binned_filename = params->filebasename + "_xpz.txt";
      _Write(params, xpz, "x", "pz", std::string(binned_filename));
    }
    if (params->do_plot_ypx)
    {
      binned_filename = params->filebasename + "_ypx.txt";
      _Write(params, ypx, "y", "px", std::string(binned_filename));
    }
    if (params->do_plot_ypy)
    {
      binned_filename = params->filebasename + "_ypy.txt";
      _Write(params, ypy, "y", "py", std::string(binned_filename));
    }
    if (params->do_plot_ypz)
    {
      binned_filename = params->filebasename + "_ypz.txt";
      _Write(params, ypz, "y", "pz", std::string(binned_filename));
    }
    if (params->do_plot_zpx)
    {
      binned_filename = params->filebasename + "_zpx.txt";
      _Write(params, zpx, "z", "px", std::string(binned_filename));
    }
    if (params->do_plot_zpy)
    {
      binned_filename = params->filebasename + "_zpy.txt";
      _Write(params, zpy, "z", "py", std::string(binned_filename));
    }
    if (params->do_plot_xpz)
    {
      binned_filename = params->filebasename + "_zpz.txt";
      _Write(params, zpz, "z", "pz", std::string(binned_filename));
    }
    if (params->do_plot_pxpy)
    {
      binned_filename = params->filebasename + "_pxpy.txt";
      _Write(params, pxpy, "px", "py", std::string(binned_filename));
    }
    if (params->do_plot_pxpz)
    {
      binned_filename = params->filebasename + "_pxpz.txt";
      _Write(params, pxpz, "px", "pz", std::string(binned_filename));
    }
    if (params->do_plot_pypz)
    {
      binned_filename = params->filebasename + "_pypz.txt";
      _Write(params, pypz, "py", "pz", std::string(binned_filename));
    }
    if (params->do_plot_xw)
    {
      binned_filename = params->filebasename + "_xw.txt";
      _Write(params, xw, "x", "w", std::string(binned_filename));
    }
    if (params->do_plot_Etheta)
    {
      binned_filename = params->filebasename + "_Etheta.txt";
      _Write(params, Etheta, "E", "theta", std::string(binned_filename));
    }
    if (params->do_plot_EthetaT)
    {
      binned_filename = params->filebasename + "_EthetaT.txt";
      _Write(params, EthetaT, "E", "thetaT", std::string(binned_filename));
    }
    if (params->do_plot_Espec)
    {
      binned_filename = params->filebasename + "_Espec.txt";
      _Write(params, Espec, "E", std::string(binned_filename));
    }
    if (params->do_plot_wspec)
    {
      binned_filename = params->filebasename + "_wspec.txt";
      _Write(params, wspec, "w", std::string(binned_filename));
    }
    if (params->do_plot_chspec)
    {
      binned_filename = params->filebasename + "_chspec.txt";
      _Write(params, chspec, "ch", std::string(binned_filename));
    }
    if (params->do_plot_thetaspec)
    {
      binned_filename = params->filebasename + "_thetaspec.txt";
      _Write(params, thetaspec, "theta", std::string(binned_filename));
    }
    if (params->do_plot_thetaTspec)
    {
      binned_filename = params->filebasename + "_thetaTspec.txt";
      _Write(params, thetaTspec, "thetaT", std::string(binned_filename));
    }
  }

  printf("\nfread_size=%lu\nEND\n\n", (unsigned long)fread_size);


  if (params->out_vtk)
  {
    fflush(binary_vtk);
    fclose(binary_vtk);
  }

  if (params->out_clean_bin)
  {
    fflush(binary_clean);
    fclose(binary_clean);
  }

  if (params->out_ppg)
  {
    fflush(ascii_propaga);
    fclose(ascii_propaga);
  }

  if (params->out_xyze)
  {
    fflush(ascii_xyze);
    fclose(ascii_xyze);
  }

  if (params->out_csv)
  {
    fflush(ascii_csv);
    fclose(ascii_csv);
  }

  if (!params->multifile) fclose(file_in);

  return 0;
}

