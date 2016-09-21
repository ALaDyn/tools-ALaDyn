
#include "leggi_particelle.h"


int leggi_particelle(Parametri * parametri)
{
  std::FILE *file_in = NULL;
  int indice_multifile = 0;
  int contatori[] = { 0, 0, 0 };
  float zero = 0.0f;
  long long particelle_accumulate = 0;
  double peso_accumulato = 0.0;
  double carica_accumulata = 0.0;
  size_t dim_file_in_bytes = 0, num_of_floats_in_file = 0, num_of_particles_in_file = 0, num_of_passes = 0, num_residual_particles = 0, dimensione_array_particelle = 0;
  unsigned int val[2] = { 0, 0 };
  float array_supporto8[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  float array_supporto6[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };


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

  float x = 0.0, y = 0.0, z = 0.0, px = 0.0, py = 0.0, pz = 0.0, ptot = 0.0;
  float ch = 0.0, w = 0.0;
  float gamma = 0.0, theta = 0.0, thetaT = 0.0, E = 0.0, ty = 0.0, tz = 0.0;
  float *estremi_min = NULL, *estremi_max = NULL;

  short buffshort[2] = { 0, 0 };
  float *particelle = NULL;
  char nomefile_bin[MAX_LENGTH_FILENAME];
  char nomefile_dat[MAX_LENGTH_FILENAME];
  char nomefile_extremes[MAX_LENGTH_FILENAME];
  char nomefile_vtk[MAX_LENGTH_FILENAME];
  char nomefile_bin_clean[MAX_LENGTH_FILENAME];
  char nomefile_propaga[MAX_LENGTH_FILENAME];
  char nomefile_xyze[MAX_LENGTH_FILENAME];
  char nomefile_csv[MAX_LENGTH_FILENAME];
  char nomefile_parametri[MAX_LENGTH_FILENAME];
  char nomefile_binnato[MAX_LENGTH_FILENAME];

  size_t fread_size = 0;
  double emittance_x = 0.0, emittance_y = 0.0, emittance_z = 0.0;
  double em_x2 = 0.0, em_x = 0.0, em_y2 = 0.0, em_y = 0.0, em_z2 = 0.0, em_z = 0.0;
  double em_px2 = 0.0, em_px = 0.0, em_py2 = 0.0, em_py = 0.0, em_pz2 = 0.0, em_pz = 0.0, em_xpx = 0.0, em_ypy = 0.0, em_zpz = 0.0;

  estremi_min = new float[SEI_DIMENSIONI + ALTRI_PARAMETRI];
  estremi_max = new float[SEI_DIMENSIONI + ALTRI_PARAMETRI];
  for (int i = 0; i < (SEI_DIMENSIONI + ALTRI_PARAMETRI); i++)
  {
    estremi_min[i] = (float)NUMERO_MASSIMO;
    estremi_max[i] = (float)-NUMERO_MASSIMO;
  }

  float **xw = new float*[parametri->nbin_x + 3];
  for (int i = 0; i < parametri->nbin_x + 3; i++)
  {
    xw[i] = new float[parametri->nbin_w + 3];
    for (int j = 0; j < parametri->nbin_w + 3; j++) xw[i][j] = 0.0;
  }

  float **xy = new float*[parametri->nbin_x + 3];
  for (int i = 0; i < parametri->nbin_x + 3; i++)
  {
    xy[i] = new float[parametri->nbin_y + 3];
    for (int j = 0; j < parametri->nbin_y + 3; j++) xy[i][j] = 0.0;
  }

  float **xz = new float*[parametri->nbin_x + 3];
  for (int i = 0; i < parametri->nbin_x + 3; i++)
  {
    xz[i] = new float[parametri->nbin_z + 3];
    for (int j = 0; j < parametri->nbin_z + 3; j++) xz[i][j] = 0.0;
  }

  float **yz = new float*[parametri->nbin_y + 3];
  for (int i = 0; i < parametri->nbin_y + 3; i++)
  {
    yz[i] = new float[parametri->nbin_z + 3];
    for (int j = 0; j < parametri->nbin_z + 3; j++) yz[i][j] = 0.0;
  }

  float **rcf = new float*[parametri->nbin_ty + 3];
  for (int i = 0; i < parametri->nbin_ty + 3; i++)
  {
    rcf[i] = new float[parametri->nbin_tz + 3];
    for (int j = 0; j < parametri->nbin_tz + 3; j++) rcf[i][j] = 0.0;
  }

  float **xpx = new float*[parametri->nbin_x + 3];
  for (int i = 0; i < parametri->nbin_x + 3; i++)
  {
    xpx[i] = new float[parametri->nbin_px + 3];
    for (int j = 0; j < parametri->nbin_px + 3; j++) xpx[i][j] = 0.0;
  }

  float **xpy = new float*[parametri->nbin_x + 3];
  for (int i = 0; i < parametri->nbin_x + 3; i++)
  {
    xpy[i] = new float[parametri->nbin_py + 3];
    for (int j = 0; j < parametri->nbin_py + 3; j++) xpy[i][j] = 0.0;
  }

  float **xpz = new float*[parametri->nbin_x + 3];
  for (int i = 0; i < parametri->nbin_x + 3; i++)
  {
    xpz[i] = new float[parametri->nbin_pz + 3];
    for (int j = 0; j < parametri->nbin_pz + 3; j++) xpz[i][j] = 0.0;
  }

  float **ypx = new float*[parametri->nbin_y + 3];
  for (int i = 0; i < parametri->nbin_y + 3; i++)
  {
    ypx[i] = new float[parametri->nbin_px + 3];
    for (int j = 0; j < parametri->nbin_px + 3; j++) ypx[i][j] = 0.0;
  }

  float **ypy = new float*[parametri->nbin_y + 3];
  for (int i = 0; i < parametri->nbin_y + 3; i++)
  {
    ypy[i] = new float[parametri->nbin_py + 3];
    for (int j = 0; j < parametri->nbin_py + 3; j++) ypy[i][j] = 0.0;
  }

  float **ypz = new float*[parametri->nbin_y + 3];
  for (int i = 0; i < parametri->nbin_y + 3; i++)
  {
    ypz[i] = new float[parametri->nbin_pz + 3];
    for (int j = 0; j < parametri->nbin_pz + 3; j++) ypz[i][j] = 0.0;
  }

  float **zpx = new float*[parametri->nbin_z + 3];
  for (int i = 0; i < parametri->nbin_z + 3; i++)
  {
    zpx[i] = new float[parametri->nbin_px + 3];
    for (int j = 0; j < parametri->nbin_px + 3; j++) zpx[i][j] = 0.0;
  }

  float **zpy = new float*[parametri->nbin_z + 3];
  for (int i = 0; i < parametri->nbin_z + 3; i++)
  {
    zpy[i] = new float[parametri->nbin_py + 3];
    for (int j = 0; j < parametri->nbin_py + 3; j++) zpy[i][j] = 0.0;
  }

  float **zpz = new float*[parametri->nbin_z + 3];
  for (int i = 0; i < parametri->nbin_z + 3; i++)
  {
    zpz[i] = new float[parametri->nbin_pz + 3];
    for (int j = 0; j < parametri->nbin_pz + 3; j++) zpz[i][j] = 0.0;
  }

  float **pxpy = new float*[parametri->nbin_px + 3];
  for (int i = 0; i < parametri->nbin_px + 3; i++)
  {
    pxpy[i] = new float[parametri->nbin_py + 3];
    for (int j = 0; j < parametri->nbin_py + 3; j++) pxpy[i][j] = 0.0;
  }

  float **pxpz = new float*[parametri->nbin_px + 3];
  for (int i = 0; i < parametri->nbin_px + 3; i++)
  {
    pxpz[i] = new float[parametri->nbin_pz + 3];
    for (int j = 0; j < parametri->nbin_pz + 3; j++) pxpz[i][j] = 0.0;
  }

  float **pypz = new float*[parametri->nbin_py + 3];
  for (int i = 0; i < parametri->nbin_py + 3; i++)
  {
    pypz[i] = new float[parametri->nbin_pz + 3];
    for (int j = 0; j < parametri->nbin_pz + 3; j++) pypz[i][j] = 0.0;
  }

  float **Etheta = new float*[parametri->nbin_E + 3];
  for (int i = 0; i < parametri->nbin_E + 3; i++)
  {
    Etheta[i] = new float[parametri->nbin_theta + 3];
    for (int j = 0; j < parametri->nbin_theta + 3; j++) Etheta[i][j] = 0.0;
  }

  float **EthetaT = new float*[parametri->nbin_E + 3];
  for (int i = 0; i < parametri->nbin_E + 3; i++)
  {
    EthetaT[i] = new float[parametri->nbin_thetaT + 3];
    for (int j = 0; j < parametri->nbin_thetaT + 3; j++) EthetaT[i][j] = 0.0;
  }

  float *wspec = new float[parametri->nbin_w + 3];
  for (int i = 0; i < parametri->nbin_w + 3; i++) wspec[i] = 0.0;

  float *chspec = new float[parametri->nbin_ch + 3];
  for (int i = 0; i < parametri->nbin_ch + 3; i++) chspec[i] = 0.0;

  float *Espec = new float[parametri->nbin_E + 3];
  for (int i = 0; i < parametri->nbin_E + 3; i++) Espec[i] = 0.0;

  float *thetaspec = new float[parametri->nbin_theta + 3];
  for (int i = 0; i < parametri->nbin_theta + 3; i++) thetaspec[i] = 0.0;

  float *thetaTspec = new float[parametri->nbin_thetaT + 3];
  for (int i = 0; i < parametri->nbin_thetaT + 3; i++) thetaTspec[i] = 0.0;

  sprintf(nomefile_propaga, "%s.ppg", parametri->filebasename.c_str());
  sprintf(nomefile_xyze, "%s_xyzE.ppg", parametri->filebasename.c_str());
  sprintf(nomefile_csv, "%s.csv", parametri->filebasename.c_str());
  sprintf(nomefile_vtk, "%s.vtk", parametri->filebasename.c_str());
  sprintf(nomefile_bin_clean, "%s_clean.bin", parametri->filebasename.c_str());

  sprintf(nomefile_bin, "%s.bin", parametri->filebasename.c_str());
  sprintf(nomefile_dat, "%s.dat", parametri->filebasename.c_str());
  sprintf(nomefile_extremes, "%s.extremes", parametri->filebasename.c_str());
  sprintf(nomefile_parametri, "%s.parameters", parametri->filebasename.c_str());

  memset(&nomefile_binnato[0], 0, sizeof(nomefile_binnato));

  if (!parametri->multifile) file_in = fopen(nomefile_bin, "rb");


#ifdef ENABLE_DEBUG
  std::cout << "NPTOT = " << parametri->nptot << std::endl;
  if (parametri->p[DO_BINNING])
  {
    std::cout << "XMIN = " << parametri->xmin << std::endl;
    std::cout << "XMAX = " << parametri->xmax << std::endl;
    std::cout << "YMIN = " << parametri->ymin << std::endl;
    std::cout << "YMAX = " << parametri->ymax << std::endl;
    std::cout << "ZMIN = " << parametri->zmin << std::endl;
    std::cout << "ZMAX = " << parametri->zmax << std::endl;
    std::cout << "PXMIN = " << parametri->pxmin << std::endl;
    std::cout << "PXMAX = " << parametri->pxmax << std::endl;
    std::cout << "PYMIN = " << parametri->pymin << std::endl;
    std::cout << "PYMAX = " << parametri->pymax << std::endl;
    std::cout << "PZMIN = " << parametri->pzmin << std::endl;
    std::cout << "PZMAX = " << parametri->pzmax << std::endl;
    std::cout << "TYMIN = " << parametri->tymin << std::endl;
    std::cout << "TYMAX = " << parametri->tymax << std::endl;
    std::cout << "TZMIN = " << parametri->tzmin << std::endl;
    std::cout << "TZMAX = " << parametri->tzmax << std::endl;
    std::cout << "GAMMAMIN = " << parametri->gammamin << std::endl;
    std::cout << "GAMMAMAX = " << parametri->gammamax << std::endl;
    std::cout << "THETAMIN = " << parametri->thetamin << std::endl;
    std::cout << "THETAMAX = " << parametri->thetamax << std::endl;
    std::cout << "THETARADMIN = " << parametri->thetaTmin << std::endl;
    std::cout << "THETARADMAX = " << parametri->thetaTmax << std::endl;
    std::cout << "EMIN = " << parametri->Emin << std::endl;
    std::cout << "EMAX = " << parametri->Emax << std::endl;
    std::cout << "WMIN = " << parametri->wmin << std::endl;
    std::cout << "WMAX = " << parametri->wmax << std::endl;
    std::cout << "CHMIN = " << parametri->chmin << std::endl;
    std::cout << "CHMAX = " << parametri->chmax << std::endl;
  }
#endif

  if (parametri->p[OUT_VTK])
  {
    printf("\nENABLED .vtk FILE\n");
    binary_vtk = fopen(nomefile_vtk, "wb");

    // Scrittura primo Header VTK e memorizzazione sua dimensione in contatori[0]
    contatori[0] += fprintf(binary_vtk, "# vtk DataFile Version 2.0\n");
    contatori[0] += fprintf(binary_vtk, "titolo nostro\n");
    contatori[0] += fprintf(binary_vtk, "BINARY\n");
    contatori[0] += fprintf(binary_vtk, "DATASET UNSTRUCTURED_GRID\n");
    contatori[0] += fprintf(binary_vtk, "POINTS %llu float\n", (unsigned long long int) parametri->nptot);

    fseeko(binary_vtk, contatori[0] + parametri->nptot * sizeof(float) * 3, SEEK_SET);

    //Scrittura secondo Header VTK e memorizzazione sua dimensione in contatori[1]
    //contatori[1] += fprintf(binary_vtk, "DATASET UNSTRUCTURED_GRID\n");
    contatori[1] += fprintf(binary_vtk, "POINT_DATA %llu\n", (unsigned long long int) parametri->nptot);
    //contatori[1] += fprintf(binary_vtk, "POINTS %i float\n", parametri->nptot);
    contatori[1] += fprintf(binary_vtk, "VECTORS p float\n");
    //contatori[1] += fprintf(binary_vtk, "LOOKUP_TABLE default\n");

    if (parametri->p[WEIGHT])
    {
      fseeko(binary_vtk, contatori[0] + parametri->nptot * sizeof(float) * 3 + contatori[1] + parametri->nptot * sizeof(float) * 3, SEEK_SET);

      //Scrittura terzo Header VTK e memorizzazione sua dimensione in contatori[2]
      //contatori[2] += fprintf(binary_vtk,"DATASET STRUCTURED_POINTS\n");
      //contatori[2] += fprintf(binary_vtk,"DIMENSIONS %ui %i %i\n",parametri->nptot, 1, 1);
      //contatori[2] += fprintf(binary_vtk,"ORIGIN 0 0 0\n");
      //contatori[2] += fprintf(binary_vtk,"SPACING 1 1 1\n");
      //contatori[2] += fprintf(binary_vtk,"POINT_DATA %ui\n",parametri->nptot);
      contatori[2] += fprintf(binary_vtk, "SCALARS w float 1\n");
      contatori[2] += fprintf(binary_vtk, "LOOKUP_TABLE default\n");
    }
  }

  if (parametri->p[OUT_CLEAN_BINARY])
  {
    printf("\nENABLED CLEAN .bin FILE\n");
    binary_clean = fopen(nomefile_bin_clean, "wb");
  }

  if (parametri->p[OUT_PROPAGA])
  {
    printf("\nENABLED .txt FILE FOR PROPAGA\n");
    ascii_propaga = fopen(nomefile_propaga, "w");
  }

  if (parametri->p[OUT_XYZE])
  {
    printf("\nENABLED .txt FILE WITH x, y, z, E\n");
    ascii_xyze = fopen(nomefile_xyze, "w");
  }

  if (parametri->p[OUT_CSV])
  {
    printf("\nENABLED .csv FILE FOR PARAVIEW\n");
    ascii_csv = fopen(nomefile_csv, "w");
  }

  fflush(stdout);

  while (1)
  {
    if (!parametri->multifile)
    {
      if (conta_processori >= parametri->last_cpu) break;
      if (conta_processori == 0) {
        /*skip header*/
        std::fseek(file_in, (long)parametri->header_size_bytes, SEEK_SET);
      }

      if (parametri->file_version == 1)
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
      if (parametri->p[SWAP]) swap_endian_i(&npart_loc, 1);
      if (npart_loc > (long long int) parametri->nptot || npart_loc < 0)
      {
        printf("Read a npart=%i, non valid. Exiting!", npart_loc);
        break;
      }
      dimensione_array_particelle = npart_loc;
      val[0] = (unsigned int)npart_loc;
      val[1] = (unsigned int)parametri->ndv;
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
      memset(&nomefile_bin[0], 0, sizeof(nomefile_bin));
      sprintf(nomefile_bin, "%s_%03d.bin", parametri->filebasename.c_str(), indice_multifile);

      if ((file_in = fopen(nomefile_bin, "rb")) == NULL)
      {
        printf("End of files! \n");
        break;
      }
      fseeko(file_in, 0, SEEK_END);
      dim_file_in_bytes = (int)ftello(file_in);
      rewind(file_in);
      num_of_floats_in_file = (dim_file_in_bytes / sizeof(float));
      num_of_particles_in_file = (int)(num_of_floats_in_file / parametri->ndv);
      printf("File %s_%.3i.bin has %llu particles\n", parametri->filebasename.c_str(), indice_multifile, (unsigned long long int) num_of_particles_in_file);
      fflush(stdout);
      num_of_passes = (int)((float)(num_of_particles_in_file) / (float)(MAX_NUM_OF_PARTICLES_PER_SHOT)) + 1;
      num_residual_particles = num_of_particles_in_file % MAX_NUM_OF_PARTICLES_PER_SHOT;
      dimensione_array_particelle = MIN(MAX_NUM_OF_PARTICLES_PER_SHOT, num_of_particles_in_file);
      if (dimensione_array_particelle > parametri->nptot)
      {
        printf("Read npart=%llu: not valid! Exiting...\n", (unsigned long long int) dimensione_array_particelle);
        break;
      }
      val[0] = (unsigned int)dimensione_array_particelle;
      val[1] = (unsigned int)parametri->ndv;
    }

    if (val[0] > 0)
    {
      fflush(stdout);
      for (size_t h = 0; h < num_of_passes; h++)
      {
        if (num_of_passes > 1) printf("File is very big, will be splitted in multiple readings: step %llu of %llu\n", (unsigned long long int) (h + 1), (unsigned long long int) num_of_passes);
        if (!parametri->multifile)
        {
          particelle = new float[npart_loc*parametri->ndv];
          if (parametri->file_version == 1)
          {
            fread_size = std::fread(buffshort, sizeof(short), 2, file_in);
            fread_size = std::fread(particelle, sizeof(float), npart_loc*parametri->ndv, file_in);
            fread_size = std::fread(&buff, sizeof(int), 1, file_in);
          }
          else fread_size = std::fread(particelle, sizeof(float), npart_loc*parametri->ndv, file_in);
          if (parametri->p[SWAP]) swap_endian_f(particelle, (size_t)npart_loc*parametri->ndv);
        }
        else
        {
          if (h == num_of_passes - 1 && num_of_passes > 1) dimensione_array_particelle = num_residual_particles;
          particelle = new float[dimensione_array_particelle*parametri->ndv];
          val[0] = (unsigned int)dimensione_array_particelle;
#ifdef ENABLE_DEBUG
          printf("npart_loc = %i\t\t ndv=%i\n", val[0], parametri->ndv);
#else
          printf("npart_loc = %i\t\t ndv=%i\r", val[0], parametri->ndv);
#endif
          fflush(stdout);
          fread_size = std::fread(particelle, sizeof(float), val[0] * parametri->ndv, file_in);
          if (parametri->p[SWAP]) swap_endian_f(particelle, (size_t)val[0] * parametri->ndv);
        }

        _Filtro(parametri, particelle, val, _Filtro::costruisci_filtro(parametri));

        if (parametri->p[OUT_PARAMS])
        {
          for (unsigned int i = 0; i < val[0]; i++)
          {
            if (((parametri->ndv == 6 || parametri->ndv == 7) && parametri->file_version < 3) || (parametri->ndv == 8 && parametri->file_version >= 3))
            {
              x = *(particelle + i*parametri->ndv);
              y = *(particelle + i*parametri->ndv + 1);
              z = *(particelle + i*parametri->ndv + 2);
              px = *(particelle + i*parametri->ndv + 3);
              py = *(particelle + i*parametri->ndv + 4);
              pz = *(particelle + i*parametri->ndv + 5);
              ptot = sqrt(px*px + py*py + pz*pz);
              px /= ptot;
              py /= ptot;
              pz /= ptot;
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight)
                w = *(particelle + i*parametri->ndv + 6);
              else
                w = parametri->overwrite_weight_value;
              if (parametri->file_version >= 3 && !parametri->overwrite_charge)
                ch = *(particelle + i*parametri->ndv + 7);
              else
                ch = parametri->overwrite_charge_value;
            }
            else if (((parametri->ndv == 4 || parametri->ndv == 5) && parametri->file_version < 3) || (parametri->ndv == 6 && parametri->file_version >= 3))
            {
              x = *(particelle + i*parametri->ndv);
              y = *(particelle + i*parametri->ndv + 1);
              z = 0.0;
              px = *(particelle + i*parametri->ndv + 2);
              py = *(particelle + i*parametri->ndv + 3);
              pz = 0.0;
              ptot = sqrt(px*px + py*py);
              px /= ptot;
              py /= ptot;
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight)
                w = *(particelle + i*parametri->ndv + 4);
              else
                w = parametri->overwrite_weight_value;
              if (parametri->file_version >= 3 && !parametri->overwrite_charge)
                ch = *(particelle + i*parametri->ndv + 5);
              else
                ch = parametri->overwrite_charge_value;
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

        if (parametri->p[FIND_MINMAX])
        {
          for (unsigned int i = 0; i < val[0]; i++)
          {
            if (((parametri->ndv == 6 || parametri->ndv == 7) && parametri->file_version < 3) || (parametri->ndv == 8 && parametri->file_version >= 3))
            {
              x = *(particelle + i*parametri->ndv);
              y = *(particelle + i*parametri->ndv + 1);
              z = *(particelle + i*parametri->ndv + 2);
              px = *(particelle + i*parametri->ndv + 3);
              py = *(particelle + i*parametri->ndv + 4);
              pz = *(particelle + i*parametri->ndv + 5);
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight)
                w = *(particelle + i*parametri->ndv + 6);
              else
                w = parametri->overwrite_weight_value;
              if (parametri->file_version >= 3 && !parametri->overwrite_charge)
                ch = *(particelle + i*parametri->ndv + 7);
              else
                ch = parametri->overwrite_charge_value;
              gamma = (float)(sqrt(1. + px*px + py*py + pz*pz) - 1.);     //gamma
              theta = (float)(atan2(sqrt(py*py + pz*pz), px)*180. / M_PI);  //theta nb: py e pz sono quelli trasversi in ALaDyn!
              thetaT = (float)atan(sqrt((py*py + pz*pz) / (px*px)));
              E = (float)(gamma*parametri->massa_particella_MeV);   //energia
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
            else if (((parametri->ndv == 4 || parametri->ndv == 5) && parametri->file_version < 3) || (parametri->ndv == 6 && parametri->file_version >= 3))
            {
              x = *(particelle + i*parametri->ndv);
              y = *(particelle + i*parametri->ndv + 1);
              px = *(particelle + i*parametri->ndv + 2);
              py = *(particelle + i*parametri->ndv + 3);
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight)
                w = *(particelle + i*parametri->ndv + 4);
              else
                w = parametri->overwrite_weight_value;
              if (parametri->file_version >= 3 && !parametri->overwrite_charge)
                ch = *(particelle + i*parametri->ndv + 5);
              else
                ch = parametri->overwrite_charge_value;
              gamma = (float)(sqrt(1. + px*px + py*py) - 1.);       //gamma
              theta = (float)(atan2(py, px)*180. / M_PI);       //theta
              thetaT = (float)atan(sqrt((py*py) / (px*px)));
              E = (float)(gamma*parametri->massa_particella_MeV); //energia
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
          }
        }
        if (parametri->p[DO_BINNING])
        {
          if (parametri->fai_plot_xy)     _Binnaggio(particelle, val[0], parametri->ndv, parametri, xy, "x", "y");
          if (parametri->fai_plot_xz)     _Binnaggio(particelle, val[0], parametri->ndv, parametri, xz, "x", "z");
          if (parametri->fai_plot_yz)     _Binnaggio(particelle, val[0], parametri->ndv, parametri, yz, "y", "z");
          if (parametri->fai_plot_rcf)    _Binnaggio(particelle, val[0], parametri->ndv, parametri, rcf, "ty", "tz");
          if (parametri->fai_plot_xpx)    _Binnaggio(particelle, val[0], parametri->ndv, parametri, xpx, "x", "px");
          if (parametri->fai_plot_xpy)    _Binnaggio(particelle, val[0], parametri->ndv, parametri, xpy, "x", "py");
          if (parametri->fai_plot_xpz)    _Binnaggio(particelle, val[0], parametri->ndv, parametri, xpz, "x", "pz");
          if (parametri->fai_plot_ypx)    _Binnaggio(particelle, val[0], parametri->ndv, parametri, ypx, "y", "px");
          if (parametri->fai_plot_ypy)    _Binnaggio(particelle, val[0], parametri->ndv, parametri, ypy, "y", "py");
          if (parametri->fai_plot_ypz)    _Binnaggio(particelle, val[0], parametri->ndv, parametri, ypz, "y", "pz");
          if (parametri->fai_plot_zpx)    _Binnaggio(particelle, val[0], parametri->ndv, parametri, zpx, "z", "px");
          if (parametri->fai_plot_zpy)    _Binnaggio(particelle, val[0], parametri->ndv, parametri, zpy, "z", "py");
          if (parametri->fai_plot_zpz)    _Binnaggio(particelle, val[0], parametri->ndv, parametri, zpz, "z", "pz");
          if (parametri->fai_plot_pxpy)   _Binnaggio(particelle, val[0], parametri->ndv, parametri, pxpy, "px", "py");
          if (parametri->fai_plot_pxpz)   _Binnaggio(particelle, val[0], parametri->ndv, parametri, pxpz, "px", "pz");
          if (parametri->fai_plot_pypz)   _Binnaggio(particelle, val[0], parametri->ndv, parametri, pypz, "py", "pz");
          if (parametri->fai_plot_xw)     _Binnaggio(particelle, val[0], parametri->ndv, parametri, xw, "x", "w");
          if (parametri->fai_plot_Etheta)   _Binnaggio(particelle, val[0], parametri->ndv, parametri, Etheta, "E", "theta");
          if (parametri->fai_plot_EthetaT)  _Binnaggio(particelle, val[0], parametri->ndv, parametri, EthetaT, "E", "thetaT");
          if (parametri->fai_plot_Espec)    _Binnaggio(particelle, val[0], parametri->ndv, parametri, Espec, "E");
          if (parametri->fai_plot_wspec)    _Binnaggio(particelle, val[0], parametri->ndv, parametri, wspec, "w");
          if (parametri->fai_plot_chspec)   _Binnaggio(particelle, val[0], parametri->ndv, parametri, chspec, "ch");
          if (parametri->fai_plot_thetaspec)  _Binnaggio(particelle, val[0], parametri->ndv, parametri, thetaspec, "theta");
          if (parametri->fai_plot_thetaTspec) _Binnaggio(particelle, val[0], parametri->ndv, parametri, thetaTspec, "thetaT");
        }

        if (parametri->p[OUT_PROPAGA])
        {
          if (((parametri->ndv == 6 || parametri->ndv == 7) && parametri->file_version < 3) || (parametri->ndv == 8 && parametri->file_version >= 3))
          {
            for (unsigned int i = 0; i < val[0]; i++)
            {
              x = particelle[i*parametri->ndv + 1] * ((float)1.e-4);
              y = particelle[i*parametri->ndv + 2] * ((float)1.e-4);
              z = particelle[i*parametri->ndv + 0] * ((float)1.e-4);
              px = particelle[i*parametri->ndv + 4];
              py = particelle[i*parametri->ndv + 5];
              pz = particelle[i*parametri->ndv + 3];
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight)
                w = particelle[i*parametri->ndv + 6];
              else
                w = parametri->overwrite_weight_value;
              //if (parametri->file_version >= 3 && !parametri->overwrite_charge)
              //  ch = particelle[i*parametri->ndv + 7];
              //else
              //  ch = parametri->overwrite_charge_value;

              if (i % parametri->subsample == 0) {
                if (parametri->file_particelle_E) fprintf(ascii_propaga, "%e %e %e %e %e %e %d %e 0 %d\n", x, y, z, px, py, pz, 3, w, i + 1);
                else if (parametri->file_particelle_P || parametri->file_particelle_LI || parametri->file_particelle_HI || parametri->file_particelle_generic_ion) fprintf(ascii_propaga, "%e %e %e %e %e %e %d %e 0 %d\n", x, y, z, px, py, pz, 1, w, i + 1);
              }
            }
          }
          else if (((parametri->ndv == 4 || parametri->ndv == 5) && parametri->file_version < 3) || (parametri->ndv == 6 && parametri->file_version >= 3))
          {
            for (unsigned int i = 0; i < val[0]; i++)
            {
              x = particelle[i*parametri->ndv + 1] * ((float)1.e-4);
              z = particelle[i*parametri->ndv + 0] * ((float)1.e-4);
              px = particelle[i*parametri->ndv + 3];
              pz = particelle[i*parametri->ndv + 2];
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight)
                w = particelle[i*parametri->ndv + 4];
              else
                w = parametri->overwrite_weight_value;
              //if (parametri->file_version >= 3 && !parametri->overwrite_charge)
              //  ch = particelle[i*parametri->ndv + 5];
              //else
              //  ch = parametri->overwrite_charge_value;

              if (i % parametri->subsample == 0) {
                if (parametri->file_particelle_E) fprintf(ascii_propaga, "%e 0 %e %e 0 %e %d %e 0 %d\n", x, z, px, pz, 3, w, i + 1);
                else if (parametri->file_particelle_P || parametri->file_particelle_LI || parametri->file_particelle_HI || parametri->file_particelle_generic_ion) fprintf(ascii_propaga, "%e 0 %e %e 0 %e %d %e 0 %d\n", x, z, px, pz, 1, w, i + 1);
              }
            }
          }
        }



        if (parametri->p[OUT_XYZE])
        {
          if (((parametri->ndv == 6 || parametri->ndv == 7) && parametri->file_version < 3) || (parametri->ndv == 8 && parametri->file_version >= 3))
          {
            for (unsigned int i = 0; i < val[0]; i++)
            {
              x = particelle[i*parametri->ndv + 1] * ((float)1.e-4);
              y = particelle[i*parametri->ndv + 2] * ((float)1.e-4);
              z = particelle[i*parametri->ndv + 0] * ((float)1.e-4);
              px = particelle[i*parametri->ndv + 4];
              py = particelle[i*parametri->ndv + 5];
              pz = particelle[i*parametri->ndv + 3];
              gamma = (float)(sqrt(1. + px*px + py*py + pz*pz) - 1.);     //gamma
              E = (float)(gamma*parametri->massa_particella_MeV);         //energia
              if (i % parametri->subsample == 0) fprintf(ascii_xyze, "%e %e %e %e\n", x, y, z, E);
            }
          }
          else if (((parametri->ndv == 4 || parametri->ndv == 5) && parametri->file_version < 3) || (parametri->ndv == 6 && parametri->file_version >= 3))
          {
            for (unsigned int i = 0; i < val[0]; i++)
            {
              x = particelle[i*parametri->ndv + 1] * ((float)1.e-4);
              z = particelle[i*parametri->ndv + 0] * ((float)1.e-4);
              px = particelle[i*parametri->ndv + 3];
              pz = particelle[i*parametri->ndv + 2];
              gamma = (float)(sqrt(1. + px*px + py*py) - 1.);       //gamma
              E = (float)(gamma*parametri->massa_particella_MeV);   //energia
              if (i % parametri->subsample == 0) fprintf(ascii_xyze, "%e %e %e\n", x, z, E);
            }
          }
        }



        if (parametri->p[OUT_CSV])
        {
          for (unsigned int i = 0; i < val[0]; i++)
          {
            x = particelle[i*parametri->ndv + 0];
            y = particelle[i*parametri->ndv + 1];
            if (((parametri->ndv == 4 || parametri->ndv == 5) && parametri->file_version < 3) || (parametri->ndv == 6 && parametri->file_version >= 3)) z = 0;
            else z = particelle[i*parametri->ndv + 2];
            px = particelle[i*parametri->ndv + 3];
            py = particelle[i*parametri->ndv + 4];
            if (((parametri->ndv == 4 || parametri->ndv == 5) && parametri->file_version < 3) || (parametri->ndv == 6 && parametri->file_version >= 3)) pz = 0;
            else pz = particelle[i*parametri->ndv + 5];

            if (((parametri->ndv == 4 || parametri->ndv == 5) && parametri->file_version < 3) || (parametri->ndv == 6 && parametri->file_version >= 3)) {
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight) w = particelle[i*parametri->ndv + 4];
              else w = parametri->overwrite_weight_value;
            }
            else {
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight) w = particelle[i*parametri->ndv + 6];
              else w = parametri->overwrite_weight_value;
            }

            if (((parametri->ndv == 4 || parametri->ndv == 5) && parametri->file_version < 3) || (parametri->ndv == 6 && parametri->file_version >= 3)) {
              if (parametri->file_version >= 3 && !parametri->overwrite_charge) ch = particelle[i*parametri->ndv + 5];
              else ch = parametri->overwrite_charge_value;
            }
            else {
              if (parametri->file_version >= 3 && !parametri->overwrite_charge) ch = particelle[i*parametri->ndv + 7];
              else ch = parametri->overwrite_charge_value;
            }


            if (i % parametri->subsample == 0) fprintf(ascii_csv, "%e, %e, %e, %e, %e, %e, %e, %e\n", x, y, z, px, py, pz, w, ch);
          }
        }



        if (parametri->p[OUT_CLEAN_BINARY])
        {
          for (unsigned int i = 0; i < val[0]; i++)
          {
            array_supporto8[0] = particelle[i*parametri->ndv + 0];
            array_supporto8[1] = particelle[i*parametri->ndv + 1];
            if (((parametri->ndv == 4 || parametri->ndv == 5) && parametri->file_version < 3) || (parametri->ndv == 6 && parametri->file_version >= 3)) array_supporto8[2] = 0;
            else array_supporto8[2] = particelle[i*parametri->ndv + 2];
            array_supporto8[3] = particelle[i*parametri->ndv + 3];
            array_supporto8[4] = particelle[i*parametri->ndv + 4];
            if (((parametri->ndv == 4 || parametri->ndv == 5) && parametri->file_version < 3) || (parametri->ndv == 6 && parametri->file_version >= 3)) array_supporto8[5] = 0;
            else array_supporto8[5] = particelle[i*parametri->ndv + 5];

            if (((parametri->ndv == 4 || parametri->ndv == 5) && parametri->file_version < 3) || (parametri->ndv == 6 && parametri->file_version >= 3)) {
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight) array_supporto8[6] = particelle[i*parametri->ndv + 4];
              else array_supporto8[6] = parametri->overwrite_weight_value;
            }
            else {
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight) array_supporto8[6] = particelle[i*parametri->ndv + 6];
              else array_supporto8[6] = parametri->overwrite_weight_value;
            }

            if (((parametri->ndv == 4 || parametri->ndv == 5) && parametri->file_version < 3) || (parametri->ndv == 6 && parametri->file_version >= 3)) {
              if (parametri->file_version >= 3 && !parametri->overwrite_charge) array_supporto8[7] = particelle[i*parametri->ndv + 5];
              else array_supporto8[7] = parametri->overwrite_charge_value;
            }
            else {
              if (parametri->file_version >= 3 && !parametri->overwrite_charge) array_supporto8[7] = particelle[i*parametri->ndv + 7];
              else array_supporto8[7] = parametri->overwrite_charge_value;
            }

            if (i % parametri->subsample == 0) fwrite((void*)(array_supporto8), sizeof(float), 8, binary_clean);
          }
        }


        if (parametri->p[OUT_VTK])
        {
          if (parametri->endian_machine == 0)
          {
            swap_endian_f(particelle, (size_t)val[0] * parametri->ndv);
          }
          int weight_index, charge_index;

          if (((parametri->ndv == 4 || parametri->ndv == 5) && parametri->file_version < 3) || (parametri->ndv == 6 && parametri->file_version >= 3))
          {
            weight_index = 4;
            if (parametri->file_version >= 3) charge_index = 5;
            else charge_index = 0;
            // scrittura coordinate x, y, z
            fseeko(binary_vtk, contatori[0] + particelle_accumulate * sizeof(float) * 3, SEEK_SET);
            for (unsigned int i = 0; i < val[0]; i += parametri->ndv)
            {
              fwrite((void*)(particelle + i), sizeof(float), 2, binary_vtk);
              fwrite((void*)&zero, sizeof(float), 1, binary_vtk);
            }

            // scrittura momenti px, py, pz
            fseeko(binary_vtk, contatori[0] + parametri->nptot * sizeof(float) * 3 + contatori[1] + particelle_accumulate * sizeof(float) * 3, SEEK_SET);
            for (unsigned int i = 2; i < val[0]; i += parametri->ndv)
            {
              fwrite((void*)(particelle + i), sizeof(float), 2, binary_vtk);
              fwrite((void*)&zero, sizeof(float), 1, binary_vtk);
            }
          }

          else if (((parametri->ndv == 6 || parametri->ndv == 7) && parametri->file_version < 3) || (parametri->ndv == 8 && parametri->file_version >= 3))
          {
            weight_index = 6;
            if (parametri->file_version >= 3) charge_index = 7;
            else charge_index = 0;
            // scrittura coordinate x, y, z
            fseeko(binary_vtk, contatori[0] + particelle_accumulate * sizeof(float) * 3, SEEK_SET);
            for (unsigned int i = 0; i < val[0]; i += parametri->ndv)
              fwrite((void*)(particelle + i), sizeof(float), 3, binary_vtk);

            // scrittura momenti px, py, pz
            fseeko(binary_vtk, contatori[0] + parametri->nptot * sizeof(float) * 3 + contatori[1] + particelle_accumulate * sizeof(float) * 3, SEEK_SET);
            for (unsigned int i = 3; i < val[0]; i += parametri->ndv)
              fwrite((void*)(particelle + i), sizeof(float), 3, binary_vtk);
          }

          else printf("\nparametri->ndv imprevisto: %i\n", parametri->ndv);


          if (parametri->p[WEIGHT] && !parametri->overwrite_weight)
          {
            // scrittura pesi
            fseeko(binary_vtk, contatori[0] + parametri->nptot * sizeof(float) * 3 + contatori[1] + parametri->nptot * sizeof(float) * 3 + contatori[2] + particelle_accumulate * sizeof(float), SEEK_SET);
            for (unsigned int i = weight_index; i < val[0]; i += parametri->ndv)
              fwrite((void*)(particelle + i), sizeof(float), 1, binary_vtk);
          }
          else if (parametri->p[WEIGHT] && parametri->overwrite_weight)
          {
            // scrittura pesi sovrascritti da linea comando
            fseeko(binary_vtk, contatori[0] + parametri->nptot * sizeof(float) * 3 + contatori[1] + parametri->nptot * sizeof(float) * 3 + contatori[2] + particelle_accumulate * sizeof(float), SEEK_SET);
            for (unsigned int i = weight_index; i < val[0]; i += parametri->ndv)
              fwrite((void*)(&(parametri->overwrite_weight_value)), sizeof(float), 1, binary_vtk);
          }

          if (parametri->file_version >= 3 && !parametri->overwrite_charge)
          {
            // scrittura pesi
            fseeko(binary_vtk, contatori[0] + parametri->nptot * sizeof(float) * 3 + contatori[1] + parametri->nptot * sizeof(float) * 3 + contatori[2] + particelle_accumulate * sizeof(float), SEEK_SET);
            for (unsigned int i = charge_index; i < val[0]; i += parametri->ndv)
              fwrite((void*)(particelle + i), sizeof(float), 1, binary_vtk);
          }
          else if (parametri->file_version >= 3 && parametri->overwrite_charge)
          {
            // scrittura pesi sovrascritti da linea comando
            fseeko(binary_vtk, contatori[0] + parametri->nptot * sizeof(float) * 3 + contatori[1] + parametri->nptot * sizeof(float) * 3 + contatori[2] + particelle_accumulate * sizeof(float), SEEK_SET);
            for (unsigned int i = charge_index; i < val[0]; i += parametri->ndv)
              fwrite((void*)(&(parametri->overwrite_weight_value)), sizeof(float), 1, binary_vtk);
          }
        }
        particelle_accumulate += val[0];
        delete[] particelle;
        particelle = NULL;
      }
    }
    indice_multifile++;
    conta_processori++;
  }


  if (parametri->p[OUT_PARAMS])
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

    parameters = fopen(nomefile_parametri, "w");
    printf("\nWriting the parameters file\n");
    fprintf(parameters, "ncpu_x=%i\n", parametri->ncpu_x);
    fprintf(parameters, "ncpu_y=%i\n", parametri->ncpu_y);
    fprintf(parameters, "ncpu_z=%i\n", parametri->ncpu_z);
    fprintf(parameters, "npx_resampled=%llu\n", (unsigned long long int) parametri->npx_ricampionati);
    fprintf(parameters, "npy_resampled=%llu\n", (unsigned long long int) parametri->npy_ricampionati);
    fprintf(parameters, "npz_resampled=%llu\n", (unsigned long long int) parametri->npz_ricampionati);
    fprintf(parameters, "npx_resampled_per_cpu=%llu\n", (unsigned long long int) parametri->npx_ricampionati_per_cpu);
    fprintf(parameters, "npy_resampled_per_cpu=%llu\n", (unsigned long long int) parametri->npy_ricampionati_per_cpu);
    fprintf(parameters, "npz_resampled_per_cpu=%llu\n", (unsigned long long int) parametri->npz_ricampionati_per_cpu);
    fprintf(parameters, "tnow=%f\n", parametri->tnow);
    fprintf(parameters, "xmin=%f\n", parametri->xmin);
    fprintf(parameters, "xmax=%f\n", parametri->xmax);
    fprintf(parameters, "ymin=%f\n", parametri->ymin);
    fprintf(parameters, "ymax=%f\n", parametri->ymax);
    fprintf(parameters, "zmin=%f\n", parametri->zmin);
    fprintf(parameters, "zmax=%f\n", parametri->zmax);
    ////////////////////////////////////////////////////////////////////////////
    fprintf(parameters, "rms_emittance_x=%g\n", emittance_x);    //massa particelle su massa elettrone
    fprintf(parameters, "rms_emittance_y=%g\n", emittance_y);    //massa particelle su massa elettrone
    fprintf(parameters, "rms_emittance_z=%g\n", emittance_z);    //massa particelle su massa elettrone
    ////////////////////////////////////////////////////////////////////////////
    fclose(parameters);
  }

  if (parametri->p[FIND_MINMAX])
  {
    Estremi_out.open(nomefile_extremes);
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
    if (parametri->p[OUT_PARAMS]) Estremi_out << "collected_weight = " << peso_accumulato << std::endl;
    if (parametri->p[OUT_PARAMS]) Estremi_out << "collected_charge = " << carica_accumulata << std::endl;
    Estremi_out.close();
  }


  if (parametri->p[DO_BINNING])
  {
    if (parametri->fai_plot_xy)
    {
      sprintf(nomefile_binnato, "%s_xy.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, xy, "x", "y", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_xz)
    {
      sprintf(nomefile_binnato, "%s_xz.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, xz, "x", "z", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_yz)
    {
      sprintf(nomefile_binnato, "%s_yz.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, yz, "y", "z", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_rcf)
    {
      sprintf(nomefile_binnato, "%s_rcf.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, rcf, "ty", "tz", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_xpx)
    {
      sprintf(nomefile_binnato, "%s_xpx.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, xpx, "x", "px", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_xpy)
    {
      sprintf(nomefile_binnato, "%s_xpy.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, xpy, "x", "py", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_xpz)
    {
      sprintf(nomefile_binnato, "%s_xpz.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, xpz, "x", "pz", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_ypx)
    {
      sprintf(nomefile_binnato, "%s_ypx.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, ypx, "y", "px", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_ypy)
    {
      sprintf(nomefile_binnato, "%s_ypy.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, ypy, "y", "py", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_ypz)
    {
      sprintf(nomefile_binnato, "%s_ypz.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, ypz, "y", "pz", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_zpx)
    {
      sprintf(nomefile_binnato, "%s_zpx.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, zpx, "z", "px", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_zpy)
    {
      sprintf(nomefile_binnato, "%s_zpy.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, zpy, "z", "py", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_xpz)
    {
      sprintf(nomefile_binnato, "%s_zpz.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, zpz, "z", "pz", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_pxpy)
    {
      sprintf(nomefile_binnato, "%s_pxpy.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, pxpy, "px", "py", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_pxpz)
    {
      sprintf(nomefile_binnato, "%s_pxpz.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, pxpz, "px", "pz", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_pypz)
    {
      sprintf(nomefile_binnato, "%s_pypz.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, pypz, "py", "pz", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_xw)
    {
      sprintf(nomefile_binnato, "%s_xw.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, xw, "x", "w", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_Etheta)
    {
      sprintf(nomefile_binnato, "%s_Etheta.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, Etheta, "E", "theta", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_EthetaT)
    {
      sprintf(nomefile_binnato, "%s_EthetaT.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, EthetaT, "E", "thetaT", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_Espec)
    {
      sprintf(nomefile_binnato, "%s_Espec.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, Espec, "E", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_wspec)
    {
      sprintf(nomefile_binnato, "%s_wspec.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, wspec, "w", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_chspec)
    {
      sprintf(nomefile_binnato, "%s_chspec.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, chspec, "ch", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_thetaspec)
    {
      sprintf(nomefile_binnato, "%s_thetaspec.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, thetaspec, "theta", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_thetaTspec)
    {
      sprintf(nomefile_binnato, "%s_thetaTspec.txt", parametri->filebasename.c_str());
      _Scrittura(parametri, thetaTspec, "thetaT", std::string(nomefile_binnato));
    }
  }

  printf("\nfread_size=%lu\nEND\n\n", (unsigned long)fread_size);


  if (parametri->p[OUT_VTK])
  {
    fflush(binary_vtk);
    fclose(binary_vtk);
  }

  if (parametri->p[OUT_CLEAN_BINARY])
  {
    fflush(binary_clean);
    fclose(binary_clean);
  }

  if (parametri->p[OUT_PROPAGA])
  {
    fflush(ascii_propaga);
    fclose(ascii_propaga);
  }

  if (parametri->p[OUT_XYZE])
  {
    fflush(ascii_xyze);
    fclose(ascii_xyze);
  }

  if (parametri->p[OUT_CSV])
  {
    fflush(ascii_csv);
    fclose(ascii_csv);
  }

  if (!parametri->multifile) fclose(file_in);

  return 0;
}

