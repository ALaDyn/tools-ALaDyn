#ifndef __LEGGI_PARTICELLE_C
#define __LEGGI_PARTICELLE_C

#include "leggi_binario_ALaDyn_fortran.h"

#define MAX_NUM_OF_PARTICLES_PER_SHOT 10000000			// in realta' il numero massimo che viene caricato in memoria e' il doppio di questo -1


class Particella_v1 {
  double x, y, z;
  double px, py, pz;
};
class Particella_v2 {
  double x, y, z;
  double px, py, pz;
  double weight;
};
class Particella_v3 {
  double x, y, z;
  double px, py, pz;
  float weight;
  int charge;
};
union Particella_vX {
  Particella_v1 particella_v1;
  Particella_v2 particella_v2;
  Particella_v3 particella_v3;
};



int leggi_particelle(int argc, const char ** argv, Parametri * parametri)
{
  std::ostringstream nomefile_bin, nomefile_dat, nomefile_Estremi;
  nomefile_bin << std::string(argv[1]) << ".bin";
  nomefile_dat << std::string(argv[1]) << ".dat";
  std::string riga_persa;
  char* trascura;
  trascura = new char[MAX_LENGTH_FILENAME];
  std::FILE *file_in = NULL;
  int indice_multifile = 0;
  int contatori[] = { 0, 0, 0 };
  float zero = 0.0f;
  long long particelle_accumulate = 0;
  double peso_accumulato = 0.0;
  int dim_file_in_bytes, num_of_floats_in_file, num_of_particles_in_file, num_of_passes, num_residual_particles;
  int dimensione_array_particelle;
  unsigned int val[2];
  float array_supporto7[7];
  float array_supporto5[5];
  double_as_two_float acc;

  if (!parametri->multifile) file_in = fopen(nomefile_bin.str().c_str(), "rb");

  std::FILE *binary_vtk;
  std::FILE *binary_clean;
  std::FILE *ascii_propaga;
  std::FILE *ascii_xyze;
  std::FILE *ascii_csv;
  std::FILE *parameters;
  std::ofstream Estremi_out;
  int conta_processori = 0;

  const int out_swap = parametri->p[SWAP];
  const int out_vtk = parametri->p[OUT_VTK];
  const int out_ascii_propaga = parametri->p[OUT_PROPAGA];
  const int out_ascii_csv = parametri->p[OUT_CSV];
  const int out_parameters = parametri->p[OUT_PARAMS];
  const int fai_binning = parametri->p[DO_BINNING];
  const int cerca_minmax = parametri->p[FIND_MINMAX];
  const int out_xyzE = parametri->p[OUT_XYZE];
  const int out_clean_binary = parametri->p[OUT_CLEAN_BINARY];

  const int stop_at_cpu_number = parametri->last_cpu;
  const int sottocampionamento = parametri->subsample;


  int N_param, *int_param, npart_loc;
  int buff, pID;

  float x, y, z, px, py, pz, w, ptot;
  float gamma, theta, thetaT, E, ty, tz;
  float rx, ry, rz, ux, uy, uz, wgh;
  float *estremi_min, *estremi_max;

  short buffshort[2];
  float *particelle = NULL;
  float *real_param = NULL;
  char nomefile_vtk[MAX_LENGTH_FILENAME];
  char nomefile_bin_clean[MAX_LENGTH_FILENAME];
  char nomefile_propaga[MAX_LENGTH_FILENAME];
  char nomefile_xyze[MAX_LENGTH_FILENAME];
  char nomefile_csv[MAX_LENGTH_FILENAME];
  char nomefile_parametri[MAX_LENGTH_FILENAME];
  char nomefile_binnato[MAX_LENGTH_FILENAME];

  size_t fread_size = 0;
  int npe, nx, nz, ibx, iby, ibz, model, dmodel, nsp, ndimen, lpord, deord, nptot, np_loc, ny_loc, ndv;
  ndv = parametri->p[NCOLONNE];
  float tnow, xmin, xmax, ymin, ymax, zmin, zmax, w0x, w0y, nrat, a0, lam0, E0, ompe, xt_in, xt_end, charge, mass, np_over_nm;
  double emittance_x = 0.0, emittance_y = 0.0, emittance_z = 0.0;
  double em_weight = 0.0, em_x2 = 0.0, em_x = 0.0, em_y2 = 0.0, em_y = 0.0, em_z2 = 0.0, em_z = 0.0;
  double em_px2 = 0.0, em_px = 0.0, em_py2 = 0.0, em_py = 0.0, em_pz2 = 0.0, em_pz = 0.0, em_xpx = 0.0, em_ypy = 0.0, em_zpz = 0.0;
  int tipo = 0;
  if (argv[1][0] == 'E') tipo = 3;
  else if (argv[1][0] == 'P') tipo = 1;
  else printf("Tipo non riconosciuto!\n");

  if (parametri->aladyn_version == 1)
  {
    fread_size = std::fread((void*)&buff, sizeof(int), 1, file_in);
    fread_size = std::fread((void*)&N_param, sizeof(int), 1, file_in);
    fread_size = std::fread((void*)&buff, sizeof(int), 1, file_in);
    if (out_swap) swap_endian_i(&N_param, 1);
    printf("numero parametri %i\n", N_param);
    fflush(stdout);
    int_param = new int[N_param];
    real_param = new float[N_param];
    //	int_param=(int*)malloc(N_param*sizeof(int));
    //	real_param=(float*)malloc(N_param*sizeof(float));
    fread_size = std::fread(&buff, sizeof(int), 1, file_in);
    fread_size = std::fread(int_param, sizeof(int), N_param, file_in);
    fread_size = std::fread(&buff, sizeof(int), 1, file_in);
    fread_size = std::fread(&buff, sizeof(int), 1, file_in);
    fread_size = std::fread(real_param, sizeof(float), N_param, file_in);
    fread_size = std::fread(&buff, sizeof(int), 1, file_in);
    if (out_swap) swap_endian_i(int_param, N_param);
    if (out_swap) swap_endian_f(real_param, N_param);
    npe = int_param[0];				//numero processori
    nx = int_param[1];
    ny_loc = int_param[2];
    nz = int_param[3];
    ibx = int_param[4];
    iby = int_param[5];
    ibz = int_param[6];
    model = int_param[7];				//modello di laser utilizzato
    dmodel = int_param[8];			//modello di condizioni iniziali
    nsp = int_param[9];				//numero di speci
    ndimen = int_param[10];			//numero di componenti dello spazio dei momenti
    np_loc = int_param[11];
    lpord = int_param[12];			//ordine dello schema leapfrog
    deord = int_param[13];			//ordine derivate
    //		int_param[14];					// current type
    pID = int_param[15];				// tipo delle particelle
    nptot = int_param[16];			// numero di particelle da leggere nella specie
    ndv = int_param[17];				// numero di particelle da leggere nella specie
    tnow = real_param[0];				// tempo dell'output
    xmin = real_param[1];				// estremi della griglia
    xmax = real_param[2];				// estremi della griglia
    ymin = real_param[3];				// estremi della griglia
    ymax = real_param[4];				// estremi della griglia
    zmin = real_param[5];				// estremi della griglia
    zmax = real_param[6];				// estremi della griglia
    w0x = real_param[7];				// waist del laser in x
    w0y = real_param[8];				// waist del laser in y
    nrat = real_param[9];				// n over n critical
    a0 = real_param[10];				// a0 laser
    lam0 = real_param[11];			// lambda
    E0 = real_param[12];				// conversione da campi numerici a TV/m
    ompe = real_param[13];			// costante accoppiamento correnti campi
    np_over_nm = real_param[14];		// numerical2physical particles 14 
    xt_in = real_param[15];			// inizio plasma
    xt_end = real_param[16];
    charge = real_param[17];			// carica particella su carica elettrone
    mass = real_param[18];			// massa particelle su massa elettrone
  }
  else
  {
    npe = parametri->intpar[0];				//numero processori
    nx = parametri->intpar[1];
    ny_loc = parametri->intpar[2];
    nz = parametri->intpar[3];
    ibx = parametri->intpar[4];
    iby = parametri->intpar[5];
    ibz = parametri->intpar[6];
    model = parametri->intpar[7];				//modello di laser utilizzato
    dmodel = parametri->intpar[8];			//modello di condizioni iniziali
    nsp = parametri->intpar[9];				//numero di speci
    ndimen = parametri->intpar[10];			//numero di componenti dello spazio dei momenti
    np_loc = parametri->intpar[11];
    lpord = parametri->intpar[12];			//ordine dello schema leapfrog
    deord = parametri->intpar[13];			//ordine derivate
    pID = parametri->intpar[15];				// tipo delle particelle
    nptot = parametri->intpar[16];			// numero di particelle da leggere nella specie
    ndv = parametri->intpar[17];				// numero di particelle da leggere nella specie
    tnow = parametri->realpar[0];				// tempo dell'output
    xmin = parametri->realpar[1];				// estremi della griglia
    xmax = parametri->realpar[2];				// estremi della griglia
    ymin = parametri->realpar[3];				// estremi della griglia
    ymax = parametri->realpar[4];				// estremi della griglia
    zmin = parametri->realpar[5];				// estremi della griglia
    zmax = parametri->realpar[6];				// estremi della griglia
    w0x = parametri->realpar[7];				// waist del laser in x
    w0y = parametri->realpar[8];				// waist del laser in y
    nrat = parametri->realpar[9];				// n over n critical
    a0 = parametri->realpar[10];				// a0 laser
    lam0 = parametri->realpar[11];			// lambda
    E0 = parametri->realpar[12];				// conversione da campi numerici a TV/m
    ompe = parametri->realpar[13];			// costante accoppiamento correnti campi
    np_over_nm = parametri->realpar[14];		// numerical2physical particles 14 
    xt_in = parametri->realpar[15];			// inizio plasma
    xt_end = parametri->realpar[16];
    charge = parametri->realpar[17];			// carica particella su carica elettrone
    mass = parametri->realpar[18];			// massa particelle su massa elettrone
  }



  estremi_min = new float[SEI_DIMENSIONI + ALTRI_PARAMETRI];
  estremi_max = new float[SEI_DIMENSIONI + ALTRI_PARAMETRI];
  for (int i = 0; i < (SEI_DIMENSIONI + ALTRI_PARAMETRI); i++)
  {
    estremi_min[i] = (float)NUMERO_MASSIMO;
    estremi_max[i] = (float)-NUMERO_MASSIMO;
  }

  printf("nptot=%i\n", nptot);
  fflush(stdout);


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

  float *Espec = new float[parametri->nbin_E + 3];
  for (int i = 0; i < parametri->nbin_E + 3; i++) Espec[i] = 0.0;

  float *thetaspec = new float[parametri->nbin_theta + 3];
  for (int i = 0; i < parametri->nbin_theta + 3; i++) thetaspec[i] = 0.0;

  float *thetaTspec = new float[parametri->nbin_thetaT + 3];
  for (int i = 0; i < parametri->nbin_thetaT + 3; i++) thetaTspec[i] = 0.0;



  sprintf(nomefile_propaga, "%s.ppg", argv[1]);
  sprintf(nomefile_xyze, "%s_xyzE.ppg", argv[1]);
  sprintf(nomefile_csv, "%s.csv", argv[1]);
  sprintf(nomefile_vtk, "%s.vtk", argv[1]);
  sprintf(nomefile_bin_clean, "%s_clean.bin", argv[1]);


#ifdef ENABLE_DEBUG
  if (fai_binning)
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
  }
#endif

  if (out_vtk)
  {
    printf("\nRichiesta scrittura file .vtk\n");
    binary_vtk = fopen(nomefile_vtk, "wb");

    // Scrittura primo Header VTK e memorizzazione sua dimensione in contatori[0]
    contatori[0] += fprintf(binary_vtk, "# vtk DataFile Version 2.0\n");
    contatori[0] += fprintf(binary_vtk, "titolo nostro\n");
    contatori[0] += fprintf(binary_vtk, "BINARY\n");
    contatori[0] += fprintf(binary_vtk, "DATASET UNSTRUCTURED_GRID\n");
    contatori[0] += fprintf(binary_vtk, "POINTS %i float\n", nptot);

    fseeko(binary_vtk, contatori[0] + nptot*sizeof(float) * 3, SEEK_SET);

    // Scrittura secondo Header VTK e memorizzazione sua dimensione in contatori[1]
    //		contatori[1] += fprintf(binary_vtk, "DATASET UNSTRUCTURED_GRID\n");
    contatori[1] += fprintf(binary_vtk, "POINT_DATA %i\n", nptot);
    //		contatori[1] += fprintf(binary_vtk, "POINTS %i float\n", nptot);
    contatori[1] += fprintf(binary_vtk, "VECTORS p float\n");
    //		contatori[1] += fprintf(binary_vtk, "LOOKUP_TABLE default\n");

    if (parametri->p[WEIGHT])
    {
      fseeko(binary_vtk, contatori[0] + nptot*sizeof(float) * 3 + contatori[1] + nptot*sizeof(float) * 3, SEEK_SET);

      // Scrittura terzo Header VTK e memorizzazione sua dimensione in contatori[2]
      //			contatori[2] += fprintf(binary_vtk,"DATASET STRUCTURED_POINTS\n");
      //			contatori[2] += fprintf(binary_vtk,"DIMENSIONS %i %i %i\n",nptot, 1, 1);
      //			contatori[2] += fprintf(binary_vtk,"ORIGIN 0 0 0\n");
      //			contatori[2] += fprintf(binary_vtk,"SPACING 1 1 1\n");
      //			contatori[2] += fprintf(binary_vtk,"POINT_DATA %i\n",nptot);
      contatori[2] += fprintf(binary_vtk, "SCALARS w float 1\n");
      contatori[2] += fprintf(binary_vtk, "LOOKUP_TABLE default\n");
    }
  }

  if (out_clean_binary)
  {
    printf("\nRichiesta scrittura file binario unico e pulito\n");
    binary_clean = fopen(nomefile_bin_clean, "wb");
  }


  if (out_ascii_propaga)
  {
    printf("\nRichiesta scrittura file ASCII per Propaga\n");
    ascii_propaga = fopen(nomefile_propaga, "w");
  }

  if (out_xyzE)
  {
    printf("\nRichiesta scrittura file ASCII con x, y, z, E\n");
    ascii_xyze = fopen(nomefile_xyze, "w");
  }

  if (out_ascii_csv)
  {
    printf("\nRichiesta scrittura file ASCII CSV per Paraview\n");
    ascii_csv = fopen(nomefile_csv, "w");
  }

  fflush(stdout);

  while (1)
  {
    if (!parametri->multifile)
    {
      if (conta_processori >= stop_at_cpu_number) break;
      if (parametri->aladyn_version == 1)
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
      if (out_swap) swap_endian_i(&npart_loc, 1);
      if (npart_loc > nptot || npart_loc < 0)
      {
        printf("Read a npart=%i, non valid. Exiting!", npart_loc);
        break;
      }
      dimensione_array_particelle = npart_loc;
      val[0] = (unsigned int)npart_loc;
      val[1] = (unsigned int)ndv;
#ifdef ENABLE_DEBUG
      printf("proc number \t %i \t npart=%i \n", conta_processori, npart_loc);
#else
      printf("proc number \t %i \t npart=%i \r", conta_processori, npart_loc);
#endif
      printf("\n");
      fflush(stdout);
      num_of_passes = 1;
    }
    else  //we do have multifiles i.e. Prpout00_000.bin
    {
      nomefile_bin.str("");
      nomefile_bin << std::string(argv[1]) << "_" << std::setfill('0') << std::setw(3) << indice_multifile << ".bin";
      if ((file_in = fopen(nomefile_bin.str().c_str(), "rb")) == NULL)
      {
        printf("Sono finiti i files! \n");
        break;
      }
      fseeko(file_in, 0, SEEK_END);
      dim_file_in_bytes = (int)ftello(file_in);
      rewind(file_in);
      num_of_floats_in_file = (dim_file_in_bytes / sizeof(float));
      num_of_particles_in_file = (int)(num_of_floats_in_file / ndv);
      printf("Il file %s_%.3i.bin contiene %i particelle\n", argv[1], indice_multifile, num_of_particles_in_file);
      fflush(stdout);
      num_of_passes = (int)((float)(num_of_particles_in_file) / (float)(MAX_NUM_OF_PARTICLES_PER_SHOT)) + 1;
      num_residual_particles = num_of_particles_in_file % MAX_NUM_OF_PARTICLES_PER_SHOT;
      dimensione_array_particelle = MIN(MAX_NUM_OF_PARTICLES_PER_SHOT, num_of_particles_in_file);
      if (dimensione_array_particelle > nptot || dimensione_array_particelle < 0)
      {
        printf("Read a npart=%i, non valid. Exiting!", dimensione_array_particelle);
        break;
      }
      val[0] = (unsigned int)dimensione_array_particelle;
      val[1] = (unsigned int)ndv;
    }

    if (val[0] > 0)
    {
      fflush(stdout);
      for (int h = 0; h < num_of_passes; h++)
      {
        if (num_of_passes > 1) printf("File is very big, will be splitted in multiple readings: step %i of %i\n", h + 1, num_of_passes);
        if (!parametri->multifile)
        {
          particelle = new float[npart_loc*ndv];
          //	particelle=(float*)malloc(npart_loc*ndv*sizeof(float));
          //	printf("Reading file %s.bin \n",argv[1]);
          if (parametri->aladyn_version == 1)
          {
            fread_size = std::fread(buffshort, sizeof(short), 2, file_in);
            fread_size = std::fread(particelle, sizeof(float), npart_loc*ndv, file_in);
            fread_size = std::fread(&buff, sizeof(int), 1, file_in);
          }
          else fread_size = std::fread(particelle, sizeof(float), npart_loc*ndv, file_in);
          if (out_swap) swap_endian_f(particelle, (size_t)npart_loc*ndv);
        }
        else
        {
          if (h == num_of_passes - 1 && num_of_passes > 1) dimensione_array_particelle = num_residual_particles;
          particelle = new float[dimensione_array_particelle*ndv];
          //	particelle=(float*)malloc(dimensione_array_particelle*ndv*sizeof(float));
          //	printf("File %s has been splitted, reading %s_%.3i.bin\n",argv[1],argv[1],indice_multifile);
          val[0] = (unsigned int)dimensione_array_particelle;
#ifdef ENABLE_DEBUG
          printf("npart_loc = %i\t\t ndv=%i\n", val[0], ndv);
#else
          printf("npart_loc = %i\t\t ndv=%i\r", val[0], ndv);
#endif
          printf("\n");
          fflush(stdout);
          fread_size = std::fread(particelle, sizeof(float), val[0] * ndv, file_in);
          if (out_swap) swap_endian_f((size_t)particelle, val[0] * ndv);
        }

        _Filtro(parametri, particelle, val, _Filtro::costruisci_filtro(argc, argv));

        if (out_parameters)
        {
          for (unsigned int i = 0; i < val[0]; i++)
          {
            if (ndv == 6 || ndv == 7)
            {
              x = *(particelle + i*ndv);
              y = *(particelle + i*ndv + 1);
              z = *(particelle + i*ndv + 2);
              px = *(particelle + i*ndv + 3);
              py = *(particelle + i*ndv + 4);
              pz = *(particelle + i*ndv + 5);
              ptot = sqrt(px*px + py*py + pz*pz);
              px /= ptot;
              py /= ptot;
              pz /= ptot;
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight && parametri->aladyn_version == 3)
              {
                acc.d = *(particelle + i*ndv + 6);
                em_weight = (double)acc.f[0];
              }
              else if (parametri->p[WEIGHT] && !parametri->overwrite_weight && parametri->aladyn_version < 3)
              {
                em_weight = *(particelle + i*ndv + 6);
              }
              else
                em_weight = parametri->overwrite_weight_value;
            }
            else if (ndv == 4 || ndv == 5)
            {
              x = *(particelle + i*ndv);
              y = *(particelle + i*ndv + 1);
              z = 0.0;
              px = *(particelle + i*ndv + 2);
              py = *(particelle + i*ndv + 3);
              pz = 0.0;
              ptot = sqrt(px*px + py*py);
              px /= ptot;
              py /= ptot;
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight && parametri->aladyn_version == 3)
              {
                acc.d = *(particelle + i*ndv + 4);
                em_weight = (double)acc.f[0];
              }
              else if (parametri->p[WEIGHT] && !parametri->overwrite_weight && parametri->aladyn_version < 3)
                em_weight = *(particelle + i*ndv + 4);
              else
                em_weight = parametri->overwrite_weight_value;
            }

            em_x += (double)(x*em_weight);
            em_y += (double)(y*em_weight);
            em_z += (double)(z*em_weight);
            em_x2 += (double)((x*x)*em_weight);
            em_y2 += (double)((y*y)*em_weight);
            em_z2 += (double)((z*z)*em_weight);
            em_px += (double)(px*em_weight);
            em_py += (double)(py*em_weight);
            em_pz += (double)(pz*em_weight);
            em_px2 += (double)((px*px)*em_weight);
            em_py2 += (double)((py*py)*em_weight);
            em_pz2 += (double)((pz*pz)*em_weight);
            em_xpx += (double)((x*px)*em_weight);
            em_ypy += (double)((y*py)*em_weight);
            em_zpz += (double)((z*pz)*em_weight);
            peso_accumulato += em_weight;
          }
        }

        if (cerca_minmax)
        {
          for (unsigned int i = 0; i < val[0]; i++)
          {
            if (ndv == 6 || ndv == 7)
            {
              x = *(particelle + i*ndv);
              y = *(particelle + i*ndv + 1);
              z = *(particelle + i*ndv + 2);
              px = *(particelle + i*ndv + 3);
              py = *(particelle + i*ndv + 4);
              pz = *(particelle + i*ndv + 5);
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight && parametri->aladyn_version == 3)
              {
                acc.d = *(particelle + i*ndv + 6);
                w = (double)acc.f[0];
              }
              else if (parametri->p[WEIGHT] && !parametri->overwrite_weight && parametri->aladyn_version < 3)
              {
                w = *(particelle + i*ndv + 6);
              }
              else
                w = parametri->overwrite_weight_value;
              gamma = (float)(sqrt(1. + px*px + py*py + pz*pz) - 1.);			//gamma
              theta = (float)(atan2(sqrt(py*py + pz*pz), px)*180. / M_PI);	//theta nb: py e pz sono quelli trasversi in ALaDyn!
              thetaT = (float)atan(sqrt((py*py + pz*pz) / (px*px)));
              E = (float)(gamma*parametri->massa_particella_MeV);		//energia
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
            }
            else if (ndv == 4 || ndv == 5)
            {
              x = *(particelle + i*ndv);
              y = *(particelle + i*ndv + 1);
              px = *(particelle + i*ndv + 2);
              py = *(particelle + i*ndv + 3);
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight && parametri->aladyn_version == 3)
              {
                acc.d = *(particelle + i*ndv + 4);
                w = (double)acc.f[0];
              }
              else if (parametri->p[WEIGHT] && !parametri->overwrite_weight && parametri->aladyn_version < 3)
                w = *(particelle + i*ndv + 4);
              else
                w = parametri->overwrite_weight_value;
              gamma = (float)(sqrt(1. + px*px + py*py) - 1.);				//gamma
              theta = (float)(atan2(py, px)*180. / M_PI);				//theta
              thetaT = (float)atan(sqrt((py*py) / (px*px)));
              E = (float)(gamma*parametri->massa_particella_MeV);	//energia
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
            }
          }
        }
        if (fai_binning)
        {
          if (parametri->fai_plot_xy)			_Binnaggio(particelle, val[0], ndv, parametri, xy, "x", "y");
          if (parametri->fai_plot_xz)			_Binnaggio(particelle, val[0], ndv, parametri, xz, "x", "z");
          if (parametri->fai_plot_yz)			_Binnaggio(particelle, val[0], ndv, parametri, yz, "y", "z");
          if (parametri->fai_plot_rcf)		_Binnaggio(particelle, val[0], ndv, parametri, rcf, "ty", "tz");
          if (parametri->fai_plot_xpx)		_Binnaggio(particelle, val[0], ndv, parametri, xpx, "x", "px");
          if (parametri->fai_plot_xpy)		_Binnaggio(particelle, val[0], ndv, parametri, xpy, "x", "py");
          if (parametri->fai_plot_xpz)		_Binnaggio(particelle, val[0], ndv, parametri, xpz, "x", "pz");
          if (parametri->fai_plot_ypx)		_Binnaggio(particelle, val[0], ndv, parametri, ypx, "y", "px");
          if (parametri->fai_plot_ypy)		_Binnaggio(particelle, val[0], ndv, parametri, ypy, "y", "py");
          if (parametri->fai_plot_ypz)		_Binnaggio(particelle, val[0], ndv, parametri, ypz, "y", "pz");
          if (parametri->fai_plot_zpx)		_Binnaggio(particelle, val[0], ndv, parametri, zpx, "z", "px");
          if (parametri->fai_plot_zpy)		_Binnaggio(particelle, val[0], ndv, parametri, zpy, "z", "py");
          if (parametri->fai_plot_zpz)		_Binnaggio(particelle, val[0], ndv, parametri, zpz, "z", "pz");
          if (parametri->fai_plot_pxpy)		_Binnaggio(particelle, val[0], ndv, parametri, pxpy, "px", "py");
          if (parametri->fai_plot_pxpz)		_Binnaggio(particelle, val[0], ndv, parametri, pxpz, "px", "pz");
          if (parametri->fai_plot_pypz)		_Binnaggio(particelle, val[0], ndv, parametri, pypz, "py", "pz");
          if (parametri->fai_plot_xw)			_Binnaggio(particelle, val[0], ndv, parametri, xw, "x", "w");
          if (parametri->fai_plot_Etheta)		_Binnaggio(particelle, val[0], ndv, parametri, Etheta, "E", "theta");
          if (parametri->fai_plot_EthetaT)	_Binnaggio(particelle, val[0], ndv, parametri, EthetaT, "E", "thetaT");
          if (parametri->fai_plot_Espec)		_Binnaggio(particelle, val[0], ndv, parametri, Espec, "E");
          if (parametri->fai_plot_wspec)		_Binnaggio(particelle, val[0], ndv, parametri, wspec, "w");
          if (parametri->fai_plot_thetaspec)	_Binnaggio(particelle, val[0], ndv, parametri, thetaspec, "theta");
          if (parametri->fai_plot_thetaTspec)	_Binnaggio(particelle, val[0], ndv, parametri, thetaTspec, "thetaT");
        }

        if (out_ascii_propaga)
        {
          if (ndv == 6 || ndv == 7)
          {
            for (unsigned int i = 0; i < val[0]; i++)
            {
              rx = particelle[i*ndv + 1] * ((float)1.e-4);
              ry = particelle[i*ndv + 2] * ((float)1.e-4);
              rz = particelle[i*ndv + 0] * ((float)1.e-4);
              ux = particelle[i*ndv + 4];
              uy = particelle[i*ndv + 5];
              uz = particelle[i*ndv + 3];
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight && parametri->aladyn_version == 3)
              {
                acc.d = *(particelle + i*ndv + 6);
                wgh = (double)acc.f[0];
                if (i % sottocampionamento == 0) fprintf(ascii_propaga, "%e %e %e %e %e %e %d %e 0 %d\n", rx, ry, rz, ux, uy, uz, tipo, wgh, i + 1);
              }
              else if (parametri->p[WEIGHT] && !parametri->overwrite_weight && parametri->aladyn_version < 3)
              {
                wgh = *(particelle + i*ndv + 6);
                if (i % sottocampionamento == 0) fprintf(ascii_propaga, "%e %e %e %e %e %e %d %e 0 %d\n", rx, ry, rz, ux, uy, uz, tipo, wgh, i + 1);
              }
              else
                if (i % sottocampionamento == 0) fprintf(ascii_propaga, "%e %e %e %e %e %e %d %e 0 %d\n", rx, ry, rz, ux, uy, uz, tipo, parametri->overwrite_weight_value, i + 1);

            }
          }
          else if (ndv == 4 || ndv == 5)
          {
            for (unsigned int i = 0; i < val[0]; i++)
            {
              rx = particelle[i*ndv + 1] * ((float)1.e-4);
              rz = particelle[i*ndv + 0] * ((float)1.e-4);
              ux = particelle[i*ndv + 3];
              uz = particelle[i*ndv + 2];
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight && parametri->aladyn_version == 3)
              {
                acc.d = *(particelle + i*ndv + 4);
                wgh = (double)acc.f[0];
                if (i % sottocampionamento == 0) fprintf(ascii_propaga, "%e 0 %e %e 0 %e %d %e 0 %d\n", rx, rz, ux, uz, tipo, wgh, i + 1);
              }
              else if (parametri->p[WEIGHT] && !parametri->overwrite_weight && parametri->aladyn_version < 3)
              {
                wgh = *(particelle + i*ndv + 4);
                if (i % sottocampionamento == 0) fprintf(ascii_propaga, "%e 0 %e %e 0 %e %d %e 0 %d\n", rx, rz, ux, uz, tipo, wgh, i + 1);
              }
              else
              {
                if (i % sottocampionamento == 0) fprintf(ascii_propaga, "%e 0 %e %e 0 %e %d %e 0 %d\n", rx, rz, ux, uz, tipo, parametri->overwrite_weight_value, i + 1);
              }
            }
          }
        }



        if (out_xyzE)
        {
          if (ndv == 6 || ndv == 7)
          {
            for (unsigned int i = 0; i < val[0]; i++)
            {
              rx = particelle[i*ndv + 1] * ((float)1.e-4);
              ry = particelle[i*ndv + 2] * ((float)1.e-4);
              rz = particelle[i*ndv + 0] * ((float)1.e-4);
              ux = particelle[i*ndv + 4];
              uy = particelle[i*ndv + 5];
              uz = particelle[i*ndv + 3];
              gamma = (float)(sqrt(1. + ux*ux + uy*uy + uz*uz) - 1.);			//gamma
              E = (float)(gamma*parametri->massa_particella_MeV);		//energia
              if (i % sottocampionamento == 0) fprintf(ascii_xyze, "%e %e %e %e\n", rx, ry, rz, E);
            }
          }
          else if (ndv == 4 || ndv == 5)
          {
            for (unsigned int i = 0; i < val[0]; i++)
            {
              rx = particelle[i*ndv + 1] * ((float)1.e-4);
              rz = particelle[i*ndv + 0] * ((float)1.e-4);
              ux = particelle[i*ndv + 3];
              uz = particelle[i*ndv + 2];
              gamma = (float)(sqrt(1. + ux*ux + uy*uy) - 1.);				//gamma
              E = (float)(gamma*parametri->massa_particella_MeV);	//energia
              if (i % sottocampionamento == 0) fprintf(ascii_xyze, "%e %e %e\n", rx, rz, E);
            }
          }
        }



        if (out_ascii_csv)
        {
          if (ndv == 6 || ndv == 7)
          {
            for (unsigned int i = 0; i < val[0]; i++)
            {
              rx = particelle[i*ndv + 0];
              ry = particelle[i*ndv + 1];
              rz = particelle[i*ndv + 2];
              ux = particelle[i*ndv + 3];
              uy = particelle[i*ndv + 4];
              uz = particelle[i*ndv + 5];
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight) wgh = particelle[i*ndv + 6];
              else wgh = parametri->overwrite_weight_value;
              if (i % sottocampionamento == 0) fprintf(ascii_csv, "%e, %e, %e, %e, %e, %e, %e\n", rx, ry, rz, ux, uy, uz, wgh);
            }
          }
          else if (ndv == 4 || ndv == 5)
          {
            for (unsigned int i = 0; i < val[0]; i++)
            {
              rx = particelle[i*ndv + 0];
              rz = particelle[i*ndv + 1];
              ux = particelle[i*ndv + 2];
              uz = particelle[i*ndv + 3];
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight) wgh = particelle[i*ndv + 4];
              else wgh = parametri->overwrite_weight_value;
              if (i % sottocampionamento == 0) fprintf(ascii_csv, "%e, 0, %e, %e, 0, %e, %e\n", rx, rz, ux, uz, wgh);
            }
          }
        }



        if (out_clean_binary)
        {
          if (ndv == 6 || ndv == 7)
          {
            for (unsigned int i = 0; i < val[0]; i++)
            {
              array_supporto7[0] = particelle[i*ndv + 0];
              array_supporto7[1] = particelle[i*ndv + 1];
              array_supporto7[2] = particelle[i*ndv + 2];
              array_supporto7[3] = particelle[i*ndv + 3];
              array_supporto7[4] = particelle[i*ndv + 4];
              array_supporto7[5] = particelle[i*ndv + 5];
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight) array_supporto7[6] = particelle[i*ndv + 6];
              else array_supporto7[6] = parametri->overwrite_weight_value;
              if (i % sottocampionamento == 0) fwrite((void*)(array_supporto7), sizeof(float), 7, binary_clean);
            }
          }
          else if (ndv == 4 || ndv == 5)
          {
            for (unsigned int i = 0; i < val[0]; i++)
            {
              array_supporto5[0] = particelle[i*ndv + 0];
              array_supporto5[1] = particelle[i*ndv + 1];
              array_supporto5[2] = particelle[i*ndv + 2];
              array_supporto5[3] = particelle[i*ndv + 3];
              if (parametri->p[WEIGHT] && !parametri->overwrite_weight) array_supporto5[4] = particelle[i*ndv + 4];
              else array_supporto5[4] = parametri->overwrite_weight_value;
              if (i % sottocampionamento == 0) fwrite((void*)(array_supporto5), sizeof(float), 5, binary_clean);
            }
          }
        }


        if (out_vtk)
        {
          if (parametri->endian_machine == 0)
          {
            swap_endian_f(particelle, (size_t)val[0] * ndv);
          }
          switch (ndv)
          {
          case 4:
          case 5:
            // scrittura coordinate x, y, z
            fseeko(binary_vtk, contatori[0] + particelle_accumulate*sizeof(float) * 3, SEEK_SET);
            for (unsigned int i = 0; i < val[0]; i += ndv)
              fwrite((void*)(particelle + i), sizeof(float), 2, binary_vtk),
              fwrite((void*)&zero, sizeof(float), 1, binary_vtk);

            // scrittura momenti px, py, pz
            fseeko(binary_vtk, contatori[0] + nptot*sizeof(float) * 3 + contatori[1] + particelle_accumulate*sizeof(float) * 3, SEEK_SET);
            for (unsigned int i = 2; i < val[0]; i += ndv)
              fwrite((void*)(particelle + i), sizeof(float), 2, binary_vtk),
              fwrite((void*)&zero, sizeof(float), 1, binary_vtk);

            break;
          case 6:
          case 7:
            // scrittura coordinate x, y, z
            fseeko(binary_vtk, contatori[0] + particelle_accumulate*sizeof(float) * 3, SEEK_SET);
            for (unsigned int i = 0; i < val[0]; i += ndv)
              fwrite((void*)(particelle + i), sizeof(float), 3, binary_vtk);

            // scrittura momenti px, py, pz
            fseeko(binary_vtk, contatori[0] + nptot*sizeof(float) * 3 + contatori[1] + particelle_accumulate*sizeof(float) * 3, SEEK_SET);
            for (unsigned int i = 3; i < val[0]; i += ndv)
              fwrite((void*)(particelle + i), sizeof(float), 3, binary_vtk);

            break;
          default:
            printf("ndv imprevisto: %i\n", ndv);
          }


          if (parametri->p[WEIGHT] && !parametri->overwrite_weight)
          {
            // scrittura pesi
            fseeko(binary_vtk, contatori[0] + nptot*sizeof(float) * 3 + contatori[1] + nptot*sizeof(float) * 3 + contatori[2] + particelle_accumulate*sizeof(float), SEEK_SET);
            for (unsigned int i = ndv - 1; i < val[0]; i += ndv)
              fwrite((void*)(particelle + i), sizeof(float), 1, binary_vtk);
          }
          else if (parametri->p[WEIGHT] && parametri->overwrite_weight)
          {
            // scrittura pesi sovrascritti da linea comando
            fseeko(binary_vtk, contatori[0] + nptot*sizeof(float) * 3 + contatori[1] + nptot*sizeof(float) * 3 + contatori[2] + particelle_accumulate*sizeof(float), SEEK_SET);
            for (unsigned int i = ndv - 1; i < val[0]; i += ndv)
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


  if (out_parameters)
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

    sprintf(nomefile_parametri, "%s.parameters", argv[1]);
    parameters = fopen(nomefile_parametri, "w");
    printf("\nWriting the parameters file\n");
    fprintf(parameters, "interi\n");
    fprintf(parameters, "npe=%i\n", npe);     //numero processori
    fprintf(parameters, "nx=%i\n", nx);
    fprintf(parameters, "ny=%i\n", ny_loc);
    fprintf(parameters, "nz=%i\n", nz);
    fprintf(parameters, "ibx=%i\n", ibx);
    fprintf(parameters, "iby=%i\n", iby);
    fprintf(parameters, "ibz=%i\n", ibz);
    fprintf(parameters, "model=%i\n", model);  //modello di laser utilizzato
    fprintf(parameters, "dmodel=%i\n", dmodel); //modello di condizioni iniziali
    fprintf(parameters, "nsp=%i\n", nsp);    //numero di speci
    fprintf(parameters, "ndim=%i\n", ndimen);
    fprintf(parameters, "np_loc=%i\n", np_loc);
    fprintf(parameters, "lpord=%i\n", lpord); //ordine dello schema leapfrog
    fprintf(parameters, "deord=%i\n", deord); //ordine derivate
    fprintf(parameters, "pID=%i\n", pID);
    fprintf(parameters, "nptot=%i\n", nptot);
    fprintf(parameters, "ndv=%i\n", ndv);
    fprintf(parameters, "========= fine interi\n");
    fprintf(parameters, "\n floating\n");
    fprintf(parameters, "tnow=%f\n", tnow);  //tempo dell'output
    fprintf(parameters, "xmin=%f\n", xmin);  //estremi della griglia
    fprintf(parameters, "xmax=%f\n", xmax);  //estremi della griglia
    fprintf(parameters, "ymin=%f\n", ymin);  //estremi della griglia
    fprintf(parameters, "ymax=%f\n", ymax);  //estremi della griglia
    fprintf(parameters, "zmin=%f\n", zmin);  //estremi della griglia
    fprintf(parameters, "zmax=%f\n", zmax);  //estremi della griglia
    fprintf(parameters, "w0x=%f\n", w0x);      //waist del laser in x
    fprintf(parameters, "w0y=%f\n", w0y);      //waist del laser in y
    fprintf(parameters, "nrat=%f\n", nrat);     //n over n critical
    fprintf(parameters, "a0=%f\n", a0);      // a0 laser
    fprintf(parameters, "lam0=%f\n", lam0);    // lambda
    fprintf(parameters, "E0=%f\n", E0);      //conversione da campi numerici a TV/m
    fprintf(parameters, "ompe=%f\n", ompe);    //costante accoppiamento correnti campi
    fprintf(parameters, "np_over_nm=%f\n", np_over_nm);   //numerical2physical particles 14 
    fprintf(parameters, "xt_in=%f\n", xt_in);
    fprintf(parameters, "xt_end=%f\n", xt_end);
    fprintf(parameters, "charge=%f\n", charge);  //carica particella su carica elettrone
    fprintf(parameters, "mass=%f\n", mass);    //massa particelle su massa elettrone
    //		if(WEIGHT) fprintf(parameters,"weight=%f\n",particelle[6]);    //massa particelle su massa elettrone
    ////////////////////////////////////////////////////////////////////////////
    fprintf(parameters, "rms_emittance_x=%g\n", emittance_x);    //massa particelle su massa elettrone
    fprintf(parameters, "rms_emittance_y=%g\n", emittance_y);    //massa particelle su massa elettrone
    fprintf(parameters, "rms_emittance_z=%g\n", emittance_z);    //massa particelle su massa elettrone
    fclose(parameters);
  }

  if (cerca_minmax)
  {
    nomefile_Estremi << argv[1] << ".extremes";
    Estremi_out.open(nomefile_Estremi.str().c_str());
    if (Estremi_out.fail()) printf("unable to create .extremes file");
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
    Estremi_out.close();
  }


  if (fai_binning)
  {
    if (parametri->fai_plot_xy)
    {
      sprintf(nomefile_binnato, "%s_xy.txt", argv[1]);
      _Scrittura(parametri, xy, "x", "y", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_xz)
    {
      sprintf(nomefile_binnato, "%s_xz.txt", argv[1]);
      _Scrittura(parametri, xz, "x", "z", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_yz)
    {
      sprintf(nomefile_binnato, "%s_yz.txt", argv[1]);
      _Scrittura(parametri, yz, "y", "z", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_rcf)
    {
      sprintf(nomefile_binnato, "%s_rcf.txt", argv[1]);
      _Scrittura(parametri, rcf, "ty", "tz", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_xpx)
    {
      sprintf(nomefile_binnato, "%s_xpx.txt", argv[1]);
      _Scrittura(parametri, xpx, "x", "px", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_xpy)
    {
      sprintf(nomefile_binnato, "%s_xpy.txt", argv[1]);
      _Scrittura(parametri, xpy, "x", "py", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_xpz)
    {
      sprintf(nomefile_binnato, "%s_xpz.txt", argv[1]);
      _Scrittura(parametri, xpz, "x", "pz", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_ypx)
    {
      sprintf(nomefile_binnato, "%s_ypx.txt", argv[1]);
      _Scrittura(parametri, ypx, "y", "px", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_ypy)
    {
      sprintf(nomefile_binnato, "%s_ypy.txt", argv[1]);
      _Scrittura(parametri, ypy, "y", "py", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_ypz)
    {
      sprintf(nomefile_binnato, "%s_ypz.txt", argv[1]);
      _Scrittura(parametri, ypz, "y", "pz", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_zpx)
    {
      sprintf(nomefile_binnato, "%s_zpx.txt", argv[1]);
      _Scrittura(parametri, zpx, "z", "px", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_zpy)
    {
      sprintf(nomefile_binnato, "%s_zpy.txt", argv[1]);
      _Scrittura(parametri, zpy, "z", "py", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_xpz)
    {
      sprintf(nomefile_binnato, "%s_zpz.txt", argv[1]);
      _Scrittura(parametri, zpz, "z", "pz", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_pxpy)
    {
      sprintf(nomefile_binnato, "%s_pxpy.txt", argv[1]);
      _Scrittura(parametri, pxpy, "px", "py", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_pxpz)
    {
      sprintf(nomefile_binnato, "%s_pxpz.txt", argv[1]);
      _Scrittura(parametri, pxpz, "px", "pz", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_pypz)
    {
      sprintf(nomefile_binnato, "%s_pypz.txt", argv[1]);
      _Scrittura(parametri, pypz, "py", "pz", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_xw)
    {
      sprintf(nomefile_binnato, "%s_xw.txt", argv[1]);
      _Scrittura(parametri, xw, "x", "w", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_Etheta)
    {
      sprintf(nomefile_binnato, "%s_Etheta.txt", argv[1]);
      _Scrittura(parametri, Etheta, "E", "theta", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_EthetaT)
    {
      sprintf(nomefile_binnato, "%s_EthetaT.txt", argv[1]);
      _Scrittura(parametri, EthetaT, "E", "thetaT", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_Espec)
    {
      sprintf(nomefile_binnato, "%s_Espec.txt", argv[1]);
      _Scrittura(parametri, Espec, "E", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_wspec)
    {
      sprintf(nomefile_binnato, "%s_wspec.txt", argv[1]);
      _Scrittura(parametri, wspec, "w", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_thetaspec)
    {
      sprintf(nomefile_binnato, "%s_thetaspec.txt", argv[1]);
      _Scrittura(parametri, thetaspec, "theta", std::string(nomefile_binnato));
    }
    if (parametri->fai_plot_thetaTspec)
    {
      sprintf(nomefile_binnato, "%s_thetaTspec.txt", argv[1]);
      _Scrittura(parametri, thetaTspec, "thetaT", std::string(nomefile_binnato));
    }
  }

  printf("fread_size=%lu\nFine\n\n", (unsigned long)fread_size);


  if (out_vtk)
  {
    fflush(binary_vtk);
    fclose(binary_vtk);
  }

  if (out_clean_binary)
  {
    fflush(binary_clean);
    fclose(binary_clean);
  }

  if (out_ascii_propaga)
  {
    fflush(ascii_propaga);
    fclose(ascii_propaga);
  }

  if (out_xyzE)
  {
    fflush(ascii_xyze);
    fclose(ascii_xyze);
  }

  if (out_ascii_csv)
  {
    fflush(ascii_csv);
    fclose(ascii_csv);
  }

  if (!parametri->multifile) fclose(file_in);

  return 0;
}


#endif
