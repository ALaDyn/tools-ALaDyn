#ifndef __LEGGI_CAMPI_C
#define __LEGGI_CAMPI_C
#include "leggi_binario_ALaDyn_fortran.h"


int leggi_campi(Parametri * parametri)
{
  std::ostringstream nomefile_bin;
  float taglio;
  std::vector<float> cutx, cuty, cutz;
  std::vector<size_t> gridIndex_cutx, gridIndex_cuty, gridIndex_cutz;
  int indice_multifile = 0;
  int fortran_buff;
  int header_size = 3;
  int * header = new int[header_size];
  float *buffer = NULL;
  float *x_lineout = NULL;
  char nomefile_parametri[MAX_LENGTH_FILENAME];
  char nomefile_campi[MAX_LENGTH_FILENAME];
  std::FILE *file_in = NULL;
  std::FILE *parameters = NULL;
  std::FILE *clean_fields = NULL;
  size_t fread_size = 0;
  size_t allocated_size = 0;

  printf("Expected at least %f MB of RAM occupancy\n", ((parametri->npx_ricampionati*parametri->npy_ricampionati*parametri->npz_ricampionati + parametri->npx_ricampionati + parametri->npx_ricampionati_per_cpu*parametri->npy_ricampionati_per_cpu*parametri->npz_ricampionati_per_cpu)*sizeof(float) + ((parametri->npy_ricampionati*parametri->npz_ricampionati + 1)*sizeof(buffer))) / 1024. / 1024.);
  printf("ALLOCATING MEMORY\n");
  fflush(stdout);

  float *** field = new float**[parametri->npz_ricampionati];
  for (size_t i = 0; i < parametri->npz_ricampionati; i++)
  {
    field[i] = new float*[parametri->npy_ricampionati];
    for (size_t j = 0; j < parametri->npy_ricampionati; j++)
    {
      field[i][j] = new float[parametri->npx_ricampionati];
      allocated_size += (parametri->npx_ricampionati)*sizeof(float);
    }
#ifdef ENABLE_DEBUG
    printf("Allocated %llu bytes for fields\r", allocated_size);
#endif
  }
  x_lineout = new float[parametri->npx_ricampionati];

  printf("READING\n");
  if (!parametri->multifile)
  {
    nomefile_bin.str("");
    nomefile_bin << parametri->filebasename << ".bin";
    file_in = fopen(nomefile_bin.str().c_str(), "rb");
    if (file_in == NULL) std::cout << "Unable to open file!" << std::endl;
    else std::cout << "File opened to read data!" << std::endl;

    /*skip header*/
    std::fseek(file_in, (long)parametri->header_size_bytes, SEEK_SET);
    std::cout << "Fseek of " << parametri->header_size_bytes << " bytes from beginning of file done" << std::endl << std::flush;

    for (unsigned int ipx = 0; ipx < parametri->ncpu_x; ipx++)
    {
      for (unsigned int ipz = 0; ipz < parametri->ncpu_z; ipz++)
      {
        for (unsigned int ipy = 0; ipy < parametri->ncpu_y; ipy++)
        {
          if (parametri->aladyn_version == 1) fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);
          fread_size += sizeof(int)*std::fread(header, sizeof(int), header_size, file_in);
          if (parametri->aladyn_version == 1) fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);

          if (parametri->p[SWAP]) swap_endian_i(header, header_size);

#ifdef ENABLE_DEBUG
          if (header[0] != parametri->npx_ricampionati_per_cpu ||
            header[1] != parametri->npy_ricampionati_per_cpu ||
            header[2] != parametri->npz_ricampionati_per_cpu)
            std::cout << "WARNING: unexpected number of points in this chunk!" << std::endl << std::flush;
#endif
          printf("header[] = {%i/%llu, %i/%llu, %i/%llu}, cpu[] = {%u/%u, %u/%u, %u/%u}\r", header[0], parametri->npx_ricampionati_per_cpu, header[1], parametri->npy_ricampionati_per_cpu, header[2], parametri->npz_ricampionati_per_cpu, ipx + 1, parametri->ncpu_x, ipy + 1, parametri->ncpu_y, ipz + 1, parametri->ncpu_z);
          fflush(stdout);

          buffer = new float[header[0] * header[1] * header[2]];
          if (parametri->aladyn_version == 1) fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);
          fread_size += sizeof(float)*std::fread(buffer, sizeof(float), header[0] * header[1] * header[2], file_in);
          if (parametri->aladyn_version == 1) fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);

          if (parametri->p[SWAP]) swap_endian_f(buffer, parametri->npx_ricampionati_per_cpu*parametri->npy_ricampionati_per_cpu*parametri->npz_ricampionati_per_cpu);

          for (size_t k = 0; k < header[2]; k++)
            for (size_t j = 0; j < header[1]; j++)
              for (size_t i = 0; i < header[0]; i++)
                field[k + (ipz * parametri->npz_ricampionati_per_cpu)]
                /* */[j + (ipy * parametri->npy_ricampionati_per_cpu)]
          /*       */[i + (ipx * parametri->npx_ricampionati_per_cpu)] =
            /*     */ buffer[i + j*header[0] + k*header[0] * header[1]];

          delete[] buffer;
          buffer = NULL;
        }
      }
    }

    // leggiamo ora le coordinate dei punti di griglia, presenti solo nelle versioni che possono prevedere griglia stretchata e che ancora non la scrivevano nel .dat
    // se presenti, sovrascrivono quelle lette o precostruite (se non trovate nel file .dat) dalle routine dei parametri

    fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);  // facciamo il test sul buffer Fortran della prima coordinata;
    // se esiste, non e' necessario tornare indietro perche' il buffer fortran che precede i dati non e' di alcun interesse

    if (!std::feof(file_in))
    {
      float *x_coordinates, *y_coordinates, *z_coordinates;
      x_coordinates = new float[parametri->npx_ricampionati];
      y_coordinates = new float[parametri->npy_ricampionati];
      z_coordinates = new float[parametri->npz_ricampionati];
      fread_size += sizeof(float)*std::fread(x_coordinates, sizeof(float), parametri->npx_ricampionati, file_in);
      fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);
      fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);
      fread_size += sizeof(float)*std::fread(y_coordinates, sizeof(float), parametri->npy_ricampionati, file_in);
      fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);
      fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);
      fread_size += sizeof(float)*std::fread(z_coordinates, sizeof(float), parametri->npz_ricampionati, file_in);
      fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);

      if (parametri->p[SWAP])
      {
        swap_endian_f(x_coordinates, parametri->npx_ricampionati);
        swap_endian_f(y_coordinates, parametri->npy_ricampionati);
        swap_endian_f(z_coordinates, parametri->npz_ricampionati);
      }

      parametri->xcoord.resize(parametri->npx_ricampionati, 0);
      parametri->ycoord.resize(parametri->npy_ricampionati, 0);
      parametri->zcoord.resize(parametri->npz_ricampionati, 0);

      for (size_t i = 0; i < parametri->npx_ricampionati; i++)
        parametri->xcoord[i] = x_coordinates[i];
      for (size_t i = 0; i < parametri->npy_ricampionati; i++)
        parametri->ycoord[i] = y_coordinates[i];
      for (size_t i = 0; i < parametri->npz_ricampionati; i++)
        parametri->zcoord[i] = z_coordinates[i];
    }
    else parametri->stretched_grid = false;
  }
  else
  {
    int header_size = 3;
    int * header = new int[header_size];
    while (1)
    {
      nomefile_bin.str("");
      nomefile_bin << parametri->filebasename << "_" << std::setfill('0') << std::setw(3) << indice_multifile << ".bin";
      file_in = fopen(nomefile_bin.str().c_str(), "rb");
      if (file_in == NULL)
      {
        std::cout << "End of files!" << std::endl;
        break;
      }
      else std::cout << "Opened file #" << indice_multifile << " to read data!" << std::endl;

      /*skip header*/
      std::fseek(file_in, (long)parametri->header_size_bytes, SEEK_SET);
      std::cout << "Fseek of " << parametri->header_size_bytes << " bytes from beginning of file done" << std::endl << std::flush;

      for (unsigned int ipx = 0; ipx < parametri->ncpu_x; ipx++)
      {
        for (unsigned int ipz = 0; ipz < parametri->ncpu_z; ipz++)
        {
          for (unsigned int ipy = 0; ipy < parametri->ncpu_y; ipy++)
          {
            fread_size += sizeof(int)*std::fread(header, sizeof(int), header_size, file_in);
            if (parametri->p[SWAP]) swap_endian_i(header, header_size);

#ifdef ENABLE_DEBUG
            if (header[0] != parametri->npx_ricampionati_per_cpu ||
              header[1] != parametri->npy_ricampionati_per_cpu ||
              header[2] != parametri->npz_ricampionati_per_cpu)
              std::cout << "WARNING: unexpected number of points in this chunk!" << std::endl << std::flush;
#endif

            printf("file %i, header[] = {%i/%llu, %i/%llu, %i/%llu}, cpu[] = {%u/%u, %u/%u, %u/%u}\r", indice_multifile, header[0], parametri->npx_ricampionati_per_cpu, header[1], parametri->npy_ricampionati_per_cpu, header[2], parametri->npz_ricampionati_per_cpu, ipx + 1, parametri->ncpu_x, ipy + 1, parametri->ncpu_y, ipz + 1, parametri->ncpu_z);
            fflush(stdout);

            buffer = new float[header[0] * header[1] * header[2]];
            fread_size += sizeof(float)*std::fread(buffer, sizeof(float), header[0] * header[1] * header[2], file_in);

            if (parametri->p[SWAP]) swap_endian_f(buffer, parametri->npx_ricampionati_per_cpu*parametri->npy_ricampionati_per_cpu*parametri->npz_ricampionati_per_cpu);

            for (size_t k = 0; k < header[2]; k++)
              for (size_t j = 0; j < header[1]; j++)
                for (size_t i = 0; i < header[0]; i++)
                  field[k + (ipz * parametri->npz_ricampionati_per_cpu)]
                  /* */[j + (ipy * parametri->npy_ricampionati_per_cpu)]
            /*       */[i + (ipx * parametri->npx_ricampionati_per_cpu)] =
              /*     */ buffer[i + j*header[0] + k*header[0] * header[1]];
            delete[] buffer;
            buffer = NULL;
          }
        }
      }
      indice_multifile++;
      fclose(file_in);
    }
  }

  std::cout << std::endl << "END READING: TOTAL " << fread_size / (1024. * 1024.) << " MB READ" << std::endl << std::flush;

  if (parametri->npz_ricampionati == 1 && parametri->p[OUT_GRID2D])
  {
    printf("Writing the ASCII 2D field file\n");
    sprintf(nomefile_campi, "%s.txt", parametri->filebasename.c_str());
    clean_fields = std::fopen(nomefile_campi, "wb"); // binary has just the meaning here not to convert \n to \r\n in Windows. The file is ASCII anyway ;)

#if defined(FORCE_PRINTF_BUFFER_SIZE) && (FORCE_PRINTF_BUFFER_SIZE > 0)
    // create a buffer to optimize writing to output_file, here the buffer size is 1k
    const int LEN = FORCE_PRINTF_BUFFER_SIZE;
    char * buffer_out = new char[LEN];
    if (setvbuf(clean_fields, buffer_out, _IOFBF, LEN) != 0)
      printf("Incorrect type or size of buffer for stream\n");
#ifdef ENABLE_DEBUG
    else printf("FILE * now has a buffer of 1024 bytes\n");
#endif
#endif

    fprintf(clean_fields, "# %llu \n # %llu \n # %llu\n# %g  %g \n # %g  %g\n", parametri->npx_ricampionati,
      parametri->npy_ricampionati, parametri->npz_ricampionati, parametri->xmin, parametri->ymin, parametri->xmax, parametri->ymax);
    for (size_t j = 0; j < parametri->npy_ricampionati; j++)
    {
      for (size_t i = 0; i < parametri->npx_ricampionati; i++)
      {
        fprintf(clean_fields, "%.4g %.4g %.4g\n", parametri->xcoord[i], parametri->ycoord[j], field[0][j][i]);
      }
      printf("grid[y] = {%llu/%llu}\r", j, parametri->npy_ricampionati);
    }
    fclose(clean_fields);

#if defined(FORCE_PRINTF_BUFFER_SIZE) && (FORCE_PRINTF_BUFFER_SIZE > 0)
    delete[] buffer_out;
#endif
  }


  if (parametri->npz_ricampionati == 1 && parametri->p[OUT_LINEOUT_X])
  {
    size_t myj;
    printf("Writing 1D lineout along x\n");
    sprintf(nomefile_campi, "%s_lineout.txt", parametri->filebasename.c_str());
    clean_fields = fopen(nomefile_campi, "w");

    for (size_t i = 0; i < parametri->npx_ricampionati; i++) x_lineout[i] = 0;

    for (size_t j = 0; j < parametri->npy_ricampionati; j++)
    {
      if (parametri->ycoord[j] >= 0)
      {
        myj = j;
        break;
      }
    }

    fprintf(clean_fields, "#");
    for (size_t j = myj - parametri->span; j < (myj + parametri->span + 1); j++)
    {
      fprintf(clean_fields, "%.4g\t", parametri->ycoord[j]);
      for (size_t i = 0; i < parametri->npx_ricampionati; i++)
      {
        x_lineout[i] += field[0][j][i] / (2.0f * parametri->span + 1.0f);
      }
    }
    fprintf(clean_fields, "\n");
    for (size_t i = 0; i < parametri->npx_ricampionati; i++)
    {
      fprintf(clean_fields, "%.4g %.4g\n", parametri->xcoord[i], x_lineout[i]);
    }
    fclose(clean_fields);
  }

  if (parametri->npz_ricampionati > 1 && parametri->p[OUT_LINEOUT_X])
  {
    size_t myj;
    printf("Writing 1D lineout along x\n");
    sprintf(nomefile_campi, "%s_lineout.txt", parametri->filebasename.c_str());
    clean_fields = fopen(nomefile_campi, "w");

    for (size_t i = 0; i < parametri->npx_ricampionati; i++)
      x_lineout[i] = 0;
    for (size_t j = 0; j < parametri->npy_ricampionati; j++)
    {
      if (parametri->ycoord[j] >= 0)
      {
        myj = j;
        break;
      }
    }
    fprintf(clean_fields, "#");
    for (size_t k = myj - parametri->span; k < (myj + parametri->span + 1); k++)
    {
      fprintf(clean_fields, "%.4g\t", parametri->zcoord[k]);

      for (size_t j = myj - parametri->span; j < (myj + parametri->span + 1); j++)
      {
        for (size_t i = 0; i < parametri->npx_ricampionati; i++)
        {
          x_lineout[i] += field[k][j][i] / ((2.0f * parametri->span + 1.0f)*(2.0f * parametri->span + 1.0f));
        }
      }
    }
    fprintf(clean_fields, "\n");
    for (size_t i = 0; i < parametri->npx_ricampionati; i++)
    {
      fprintf(clean_fields, "%.4g %.4g\n", parametri->xcoord[i], x_lineout[i]);
    }
    fclose(clean_fields);
  }


  if (parametri->npz_ricampionati > 1 && parametri->p[OUT_CUTZ])
  {
    if (parametri->posizioni_taglio_griglia_z.size() == 0)
    {
      taglio = parametri->zcoord[parametri->npz_ricampionati / 2];
      cutz.push_back(taglio);
      gridIndex_cutz.push_back(parametri->npz_ricampionati / 2);
    }
    else
    {
      taglio = parametri->zmin;
      int i = 0;
      for (size_t j = 0; j < parametri->posizioni_taglio_griglia_z.size(); j++)
      {
        while (taglio < parametri->posizioni_taglio_griglia_z.at(j)) taglio = parametri->zcoord[i], i++;
        cutz.push_back(taglio);
        gridIndex_cutz.push_back(i - 1);
        i = 0;
      }
    }
    for (size_t n = 0; n < cutz.size(); n++)
    {
      sprintf(nomefile_campi, "%s_cutz_%g.txt", parametri->filebasename.c_str(), cutz[n]);
      printf("Writing the 2D ASCII field file at z=%g\n", cutz[n]);
      clean_fields = fopen(nomefile_campi, "wb");

#if defined(FORCE_PRINTF_BUFFER_SIZE) && (FORCE_PRINTF_BUFFER_SIZE > 0)
      // create a buffer to optimize writing to output_file, here the buffer size is 1k
      const int LEN = FORCE_PRINTF_BUFFER_SIZE;
      char * buffer_out = new char[LEN];
      if (setvbuf(clean_fields, buffer_out, _IOFBF, LEN) != 0)
        printf("Incorrect type or size of buffer for stream\n");
#ifdef ENABLE_DEBUG
      else printf("FILE * now has a buffer of 1024 bytes\n");
#endif
#endif

      //output per gnuplot (x:y:valore) compatibile con programmino passe_par_tout togliendo i #
      fprintf(clean_fields, "# 2D cut at z=%g\n", cutz[n]);
      fprintf(clean_fields, "# %lu\n#%lu\n#%i\n", (long)parametri->npx_ricampionati, (long)parametri->npy_ricampionati, 1);
      fprintf(clean_fields, "#%f %f\n#%f %f\n", parametri->xmin, parametri->ymin, parametri->xmax, parametri->ymax);
      size_t k = gridIndex_cutz[n];
      for (size_t j = 0; j < parametri->npy_ricampionati; j++)
      {
        for (size_t i = 0; i < parametri->npx_ricampionati; i++)
        {
          fprintf(clean_fields, "%.4g %.4g %.4g\n", parametri->xcoord[i], parametri->ycoord[j], field[k][j][i]);
        }
        printf("grid[y] = {%llu/%llu}\r", j, parametri->npy_ricampionati);
      }
      fclose(clean_fields);

#if defined(FORCE_PRINTF_BUFFER_SIZE) && (FORCE_PRINTF_BUFFER_SIZE > 0)
      delete[] buffer_out;
#endif
    }
  }

  if (parametri->npz_ricampionati > 1 && parametri->p[OUT_CUTY])
  {
    if (parametri->posizioni_taglio_griglia_y.size() == 0)
    {
      taglio = parametri->ycoord[parametri->npy_ricampionati / 2];
      cuty.push_back(taglio);
      gridIndex_cuty.push_back(parametri->npy_ricampionati / 2);
    }
    else
    {
      taglio = parametri->ymin;
      int i = 0;
      for (size_t j = 0; j < parametri->posizioni_taglio_griglia_y.size(); j++)
      {
        while (taglio < parametri->posizioni_taglio_griglia_y.at(j)) taglio = parametri->ycoord[i], i++;
        cuty.push_back(taglio);
        gridIndex_cuty.push_back(i - 1);
        i = 0;
      }
    }
    for (size_t n = 0; n < cuty.size(); n++)
    {
      sprintf(nomefile_campi, "%s_cuty_%g.txt", parametri->filebasename.c_str(), cuty[n]);
      printf("Writing the 2D ASCII field file at y=%g\n", cuty[n]);
      clean_fields = fopen(nomefile_campi, "wb");

#if defined(FORCE_PRINTF_BUFFER_SIZE) && (FORCE_PRINTF_BUFFER_SIZE > 0)
      // create a buffer to optimize writing to output_file, here the buffer size is 1k
      const int LEN = FORCE_PRINTF_BUFFER_SIZE;
      char * buffer_out = new char[LEN];
      if (setvbuf(clean_fields, buffer_out, _IOFBF, LEN) != 0)
        printf("Incorrect type or size of buffer for stream\n");
#ifdef ENABLE_DEBUG
      else printf("FILE * now has a buffer of 1024 bytes\n");
#endif
#endif

      //output per gnuplot (x:y:valore) compatibile con programmino passe_par_tout togliendo i #
      fprintf(clean_fields, "# 2D cut at y=%g\n", cuty[n]);
      fprintf(clean_fields, "# %lu\n#%lu\n#%i\n", (long)parametri->npx_ricampionati, (long)parametri->npz_ricampionati, 1);
      fprintf(clean_fields, "#%f %f\n#%f %f\n", parametri->xmin, parametri->zmin, parametri->xmax, parametri->zmax);
      size_t j = gridIndex_cuty[n];
      for (size_t k = 0; k < parametri->npz_ricampionati; k++)
      {
        for (size_t i = 0; i < parametri->npx_ricampionati; i++)
        {
          fprintf(clean_fields, "%.4g %.4g %.4g\n", parametri->xcoord[i], parametri->zcoord[k], field[k][j][i]);
        }
        printf("grid[z] = {%llu/%llu}\r", k, parametri->npz_ricampionati);
      }
      fclose(clean_fields);

#if defined(FORCE_PRINTF_BUFFER_SIZE) && (FORCE_PRINTF_BUFFER_SIZE > 0)
      delete[] buffer_out;
#endif
    }
  }

  if (parametri->npz_ricampionati > 1 && parametri->p[OUT_CUTX])
  {
    if (parametri->posizioni_taglio_griglia_x.size() == 0)
    {
      taglio = parametri->xcoord[parametri->npx_ricampionati / 2];
      cutx.push_back(taglio);
      gridIndex_cutx.push_back(parametri->npx_ricampionati / 2);
    }
    else
    {
      taglio = parametri->xmin;
      int i = 0;
      for (size_t j = 0; j < parametri->posizioni_taglio_griglia_x.size(); j++)
      {
        while (taglio < parametri->posizioni_taglio_griglia_x.at(j)) taglio = parametri->xcoord[i], i++;
        cutx.push_back(taglio);
        gridIndex_cutx.push_back(i - 1);
        i = 0;
      }
    }
    for (size_t n = 0; n < cutx.size(); n++)
    {
      sprintf(nomefile_campi, "%s_cutx_%g.txt", parametri->filebasename.c_str(), cutx[n]);
      printf("Writing the 2D ASCII field file at x=%g\n", cutx[n]);
      clean_fields = fopen(nomefile_campi, "wb");

#if defined(FORCE_PRINTF_BUFFER_SIZE) && (FORCE_PRINTF_BUFFER_SIZE > 0)
      // create a buffer to optimize writing to output_file, here the buffer size is 1k
      const int LEN = FORCE_PRINTF_BUFFER_SIZE;
      char * buffer_out = new char[LEN];
      if (setvbuf(clean_fields, buffer_out, _IOFBF, LEN) != 0)
        printf("Incorrect type or size of buffer for stream\n");
#ifdef ENABLE_DEBUG
      else printf("FILE * now has a buffer of 1024 bytes\n");
#endif
#endif

      //output per gnuplot (x:y:valore) compatibile con programmino passe_par_tout togliendo i #
      fprintf(clean_fields, "# 2D cut at x=%g\n", cutx[n]);
      fprintf(clean_fields, "# %lu\n#%lu\n#%i\n", (long)parametri->npy_ricampionati, (long)parametri->npz_ricampionati, 1);
      fprintf(clean_fields, "#%f %f\n#%f %f\n", parametri->ymin, parametri->zmin, parametri->ymax, parametri->zmax);
      size_t i = gridIndex_cutx[n];
      for (size_t k = 0; k < parametri->npz_ricampionati; k++)
      {
        for (size_t j = 0; j < parametri->npy_ricampionati; j++)
        {
          fprintf(clean_fields, "%.4g %.4g %.4g\n", parametri->ycoord[j], parametri->zcoord[k], field[k][j][i]);
        }
        printf("grid[z] = {%llu/%llu}\r", k, parametri->npz_ricampionati);
      }
      fclose(clean_fields);

#if defined(FORCE_PRINTF_BUFFER_SIZE) && (FORCE_PRINTF_BUFFER_SIZE > 0)
      delete[] buffer_out;
#endif
    }
  }


  if (parametri->p[OUT_VTK])
  {
    printf("Writing the vtk fields file\n");

    float *x_coordinates, *y_coordinates, *z_coordinates;
    x_coordinates = new float[parametri->npx_ricampionati];
    y_coordinates = new float[parametri->npy_ricampionati];
    z_coordinates = new float[parametri->npz_ricampionati];

    for (size_t i = 0; i < parametri->npx_ricampionati; i++)
      x_coordinates[i] = parametri->xcoord[i];
    for (size_t i = 0; i < parametri->npy_ricampionati; i++)
      y_coordinates[i] = parametri->ycoord[i];
    for (size_t i = 0; i < parametri->npz_ricampionati; i++)
      z_coordinates[i] = parametri->zcoord[i];

    if (parametri->endian_machine == 0)
    {
      swap_endian_f(field, parametri->npz_ricampionati, parametri->npy_ricampionati, parametri->npx_ricampionati);
      swap_endian_f(x_coordinates, parametri->npx_ricampionati);
      swap_endian_f(y_coordinates, parametri->npy_ricampionati);
      swap_endian_f(z_coordinates, parametri->npz_ricampionati);
    }

    //////// DATASET UNSTRUCTURED_GRID VERSION    ////////

    sprintf(nomefile_campi, "%s.vtk", parametri->filebasename.c_str());
    clean_fields = fopen(nomefile_campi, "w");
    fprintf(clean_fields, "# vtk DataFile Version 2.0\n");
    fprintf(clean_fields, "titolo mio\n");
    fprintf(clean_fields, "BINARY\n");
    fprintf(clean_fields, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(clean_fields, "POINTS %lu float\n", (long)(parametri->npx_ricampionati*parametri->npy_ricampionati*parametri->npz_ricampionati));
    float rr[3];
    for (size_t k = 0; k < parametri->npz_ricampionati; k++)
    {
      rr[2] = z_coordinates[k];
      for (size_t j = 0; j < parametri->npy_ricampionati; j++)
      {
        rr[1] = y_coordinates[j];
        for (size_t i = 0; i < parametri->npx_ricampionati; i++)
        {
          rr[0] = x_coordinates[i];
          fwrite((void*)rr, sizeof(float), 3, clean_fields);
        }
      }
    }

    fprintf(clean_fields, "POINT_DATA %lu\n", (long)(parametri->npx_ricampionati*parametri->npy_ricampionati*parametri->npz_ricampionati));
    fprintf(clean_fields, "SCALARS %s float 1\n", parametri->support_label);
    fprintf(clean_fields, "LOOKUP_TABLE default\n");
    fwrite((void*)field, sizeof(float), parametri->npz_ricampionati*parametri->npy_ricampionati*parametri->npx_ricampionati, clean_fields);
    fclose(clean_fields);
  }



  if (parametri->p[OUT_VTK_NOSTRETCH])
  {
    printf("Writing the vtk file related to the unstretched part of the grid\n");

    int inizio_punti_non_stretchati_x, inizio_punti_non_stretchati_y, inizio_punti_non_stretchati_z;
    int fine_punti_non_stretchati_x, fine_punti_non_stretchati_y, fine_punti_non_stretchati_z;
    size_t npunti_non_stretchati_x, npunti_non_stretchati_y, npunti_non_stretchati_z;
    if (parametri->stretched_grid)
    {
      if (parametri->stretched_along_x) inizio_punti_non_stretchati_x = (int)(parametri->npx_ricampionati / 6.0);
      else inizio_punti_non_stretchati_x = 0;
      inizio_punti_non_stretchati_y = (int)(parametri->npy_ricampionati / 6.0);
      inizio_punti_non_stretchati_z = (int)(parametri->npz_ricampionati / 6.0);
      if (parametri->stretched_along_x) fine_punti_non_stretchati_x = (int)(parametri->npx_ricampionati*5.0 / 6.0);
      else fine_punti_non_stretchati_x = (int)(parametri->npx_ricampionati);
      fine_punti_non_stretchati_y = (int)(parametri->npy_ricampionati*5.0 / 6.0);
      fine_punti_non_stretchati_z = (int)(parametri->npz_ricampionati*5.0 / 6.0);
      npunti_non_stretchati_x = (size_t)(fine_punti_non_stretchati_x - inizio_punti_non_stretchati_x);
      npunti_non_stretchati_y = (size_t)(fine_punti_non_stretchati_y - inizio_punti_non_stretchati_y);
      npunti_non_stretchati_z = (size_t)(fine_punti_non_stretchati_z - inizio_punti_non_stretchati_z);
    }
    else
    {
      fine_punti_non_stretchati_x = (int)parametri->npx_ricampionati;
      npunti_non_stretchati_x = parametri->npx_ricampionati;
      fine_punti_non_stretchati_y = (int)parametri->npy_ricampionati;
      npunti_non_stretchati_y = parametri->npy_ricampionati;
      fine_punti_non_stretchati_z = (int)parametri->npz_ricampionati;
      npunti_non_stretchati_z = parametri->npz_ricampionati;
      inizio_punti_non_stretchati_x = inizio_punti_non_stretchati_y = inizio_punti_non_stretchati_z = 0;
    }

    float * field_non_stretchato = new float[npunti_non_stretchati_x*npunti_non_stretchati_y*npunti_non_stretchati_z];
    int a = 0, b = 0, c = 0;

    for (int k = inizio_punti_non_stretchati_z; k < fine_punti_non_stretchati_z; k++)
    {
      c = k - inizio_punti_non_stretchati_z;
      for (int j = inizio_punti_non_stretchati_y; j < fine_punti_non_stretchati_y; j++)
      {
        b = j - inizio_punti_non_stretchati_y;
        for (int i = inizio_punti_non_stretchati_x; i < fine_punti_non_stretchati_x; i++)
        {
          a = i - inizio_punti_non_stretchati_x;
          field_non_stretchato[a + b*npunti_non_stretchati_x + c*npunti_non_stretchati_x*npunti_non_stretchati_y] = field[k][j][i];
        }
      }
    }

    if (parametri->endian_machine == 0) swap_endian_f(field_non_stretchato, npunti_non_stretchati_x*npunti_non_stretchati_y*npunti_non_stretchati_z);



    //////// DATASET STRUCTURED_POINTS VERSION    ////////

    float xmin_non_stretchato = parametri->xcoord[inizio_punti_non_stretchati_x];
    //    float xmax_non_stretchato = parametri->xcoord[fine_punti_non_stretchati_x];
    float ymin_non_stretchato = parametri->ycoord[inizio_punti_non_stretchati_y];
    //    float ymax_non_stretchato = parametri->ycoord[fine_punti_non_stretchati_y];
    float zmin_non_stretchato = parametri->zcoord[inizio_punti_non_stretchati_z];
    //    float zmax_non_stretchato = parametri->zcoord[fine_punti_non_stretchati_z];

    sprintf(nomefile_campi, "%s_nostretch.vtk", parametri->filebasename.c_str());
    clean_fields = fopen(nomefile_campi, "wb");
    fprintf(clean_fields, "# vtk DataFile Version 2.0\n");
    fprintf(clean_fields, "titolo mio\n");
    fprintf(clean_fields, "BINARY\n");
    fprintf(clean_fields, "DATASET STRUCTURED_POINTS\n");
    fprintf(clean_fields, "DIMENSIONS %lu %lu %lu\n", (long)npunti_non_stretchati_x, (long)npunti_non_stretchati_y, (long)npunti_non_stretchati_z);
    fprintf(clean_fields, "ORIGIN %f %f %f\n", xmin_non_stretchato, ymin_non_stretchato, zmin_non_stretchato);
    fprintf(clean_fields, "SPACING %f %f %f\n",
      parametri->xcoord[parametri->npx_ricampionati / 2] - parametri->xcoord[parametri->npx_ricampionati / 2 - 1],
      parametri->ycoord[parametri->npy_ricampionati / 2] - parametri->ycoord[parametri->npy_ricampionati / 2 - 1],
      parametri->zcoord[parametri->npz_ricampionati / 2] - parametri->zcoord[parametri->npz_ricampionati / 2 - 1]);
    fprintf(clean_fields, "POINT_DATA %lu\n", (long)(npunti_non_stretchati_x*npunti_non_stretchati_y*npunti_non_stretchati_z));
    fprintf(clean_fields, "SCALARS %s float 1\n", parametri->support_label);
    fprintf(clean_fields, "LOOKUP_TABLE default\n");
    fwrite((void*)field_non_stretchato, sizeof(float), npunti_non_stretchati_x*npunti_non_stretchati_y*npunti_non_stretchati_z, clean_fields);




    /******************************************************************************
    //////// DATASET RECTILINEAR_GRID VERSION    ////////

    float *x_coordinates, *y_coordinates, *z_coordinates;
    x_coordinates=new float[npunti_non_stretchati_x];
    y_coordinates=new float[npunti_non_stretchati_y];
    z_coordinates=new float[npunti_non_stretchati_z];
    for (int i = inizio_punti_non_stretchati_x; i < fine_punti_non_stretchati_x; i++) x_coordinates[i] = parametri->xcoord[i];
    for (int i = inizio_punti_non_stretchati_y; i < fine_punti_non_stretchati_y; i++) y_coordinates[i] = parametri->ycoord[i];
    for (int i = inizio_punti_non_stretchati_z; i < fine_punti_non_stretchati_z; i++) z_coordinates[i] = parametri->zcoord[i];
    if(parametri->endian_machine == 0)
    {
    swap_endian_f(x_coordinates,npunti_non_stretchati_x);
    swap_endian_f(y_coordinates,npunti_non_stretchati_y);
    swap_endian_f(z_coordinates,npunti_non_stretchati_z);
    }

    sprintf(nomefile_campi,"%s_out.vtk",parametri->filebasename);
    clean_fields=fopen(nomefile_campi, "wb");
    printf("\nWriting the fields file\n");
    fprintf(clean_fields,"# vtk DataFile Version 2.0\n");
    fprintf(clean_fields,"titolo mio\n");
    fprintf(clean_fields,"BINARY\n");
    fprintf(clean_fields,"DATASET RECTILINEAR_GRID\n");
    fprintf(clean_fields,"DIMENSIONS %i %i %i\n",npunti_non_stretchati_x, npunti_non_stretchati_y, npunti_non_stretchati_z);
    fprintf(clean_fields,"X_COORDINATES %i float\n",npunti_non_stretchati_x);
    fwrite((void*)x_coordinates,sizeof(float),npunti_non_stretchati_x,clean_fields);
    fprintf(clean_fields,"Y_COORDINATES %i float\n",npunti_non_stretchati_y);
    fwrite((void*)y_coordinates,sizeof(float),npunti_non_stretchati_y,clean_fields);
    fprintf(clean_fields,"Z_COORDINATES %i float\n",npunti_non_stretchati_z);
    fwrite((void*)z_coordinates,sizeof(float),npunti_non_stretchati_z,clean_fields);
    fprintf(clean_fields,"POINT_DATA %i\n",npunti_non_stretchati_x*npunti_non_stretchati_y*npunti_non_stretchati_z);
    fprintf(clean_fields,"SCALARS %s float 1\n",parametri->support_label);
    fprintf(clean_fields,"LOOKUP_TABLE default\n");
    fwrite((void*)field_non_stretchato,sizeof(float),npunti_non_stretchati_x*npunti_non_stretchati_y*npunti_non_stretchati_z,clean_fields);
    ******************************************************************************/

    fclose(clean_fields);
  }

  if (parametri->p[OUT_PARAMS])
  {
    sprintf(nomefile_parametri, "%s.parameters", parametri->filebasename.c_str());
    parameters = fopen(nomefile_parametri, "w");
    printf("Writing parameters to file\n");
    fprintf(parameters, "ncpu_x=%u\n", parametri->ncpu_x);
    fprintf(parameters, "ncpu_y=%u\n", parametri->ncpu_y);
    fprintf(parameters, "ncpu_z=%u\n", parametri->ncpu_z);
    fprintf(parameters, "npx_ricampionati=%llu\n", parametri->npx_ricampionati);
    fprintf(parameters, "npy_ricampionati=%llu\n", parametri->npy_ricampionati);
    fprintf(parameters, "npz_ricampionati=%llu\n", parametri->npz_ricampionati);
    fprintf(parameters, "npx_ricampionati_per_cpu=%llu\n", parametri->npx_ricampionati_per_cpu);
    fprintf(parameters, "npy_ricampionati_per_cpu=%llu\n", parametri->npy_ricampionati_per_cpu);
    fprintf(parameters, "npz_ricampionati_per_cpu=%llu\n", parametri->npz_ricampionati_per_cpu);
    fprintf(parameters, "tnow=%f\n", parametri->tnow);
    fprintf(parameters, "xmin=%f\n", parametri->xmin);
    fprintf(parameters, "xmax=%f\n", parametri->xmax);
    fprintf(parameters, "ymin=%f\n", parametri->ymin);
    fprintf(parameters, "ymax=%f\n", parametri->ymax);
    fprintf(parameters, "zmin=%f\n", parametri->zmin);
    fprintf(parameters, "zmax=%f\n", parametri->zmax);

    fprintf(parameters, "\n\nGrid along x axis\n");
    for (unsigned int i = 0; i < parametri->npx_ricampionati; i++)
    {
      fprintf(parameters, "%.4g  ", parametri->xcoord[i]);
      if (i > 0 && i % 10 == 0) fprintf(parameters, "\n");
    }

    fprintf(parameters, "\n\nGrid along y axis\n");
    for (unsigned int i = 0; i < parametri->npy_ricampionati; i++)
    {
      fprintf(parameters, "%.4g  ", parametri->ycoord[i]);
      if (i > 0 && i % 10 == 0) fprintf(parameters, "\n");
    }

    fprintf(parameters, "\n\nGrid along z axis\n");
    for (unsigned int i = 0; i < parametri->npz_ricampionati; i++)
    {
      fprintf(parameters, "%.4g  ", parametri->zcoord[i]);
      if (i > 0 && i % 10 == 0) fprintf(parameters, "\n");
    }

    fclose(parameters);
  }


  return 0;

}

#endif

