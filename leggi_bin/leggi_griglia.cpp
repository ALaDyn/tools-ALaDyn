
#include "leggi_griglia.h"


int read_grid_file(Parameters * params)
{
  std::ostringstream bin_filename;
  aladyn_float cut;
  std::vector<aladyn_float> cutx, cuty, cutz;
  std::vector<size_t> gridIndex_cutx, gridIndex_cuty, gridIndex_cutz;
  int multifile_index = 0;
  int fortran_buff;
  int header_size = 3;
  int * header = new int[header_size];
  aladyn_float *buffer = NULL;
  aladyn_float *x_lineout = NULL;
  std::string params_filename;
  std::string grid_filename;
  std::FILE *file_in = NULL;
  std::FILE *parameters = NULL;
  std::FILE *clean_fields = NULL;
  size_t fread_size = 0;
  size_t allocated_size = 0;

  printf("Expected at least %f MB of RAM occupancy\n", ((params->npx_resampled*params->npy_resampled*params->npz_resampled + params->npx_resampled + params->npx_resampled_per_cpu*params->npy_resampled_per_cpu*params->npz_resampled_per_cpu)*sizeof(aladyn_float) + ((params->npy_resampled*params->npz_resampled + 1)*sizeof(buffer))) / 1024. / 1024.);
  printf("ALLOCATING MEMORY\n");
  fflush(stdout);

  aladyn_float *** field = new aladyn_float**[params->npz_resampled];
  for (size_t i = 0; i < params->npz_resampled; i++)
  {
    field[i] = new aladyn_float*[params->npy_resampled];
    for (size_t j = 0; j < params->npy_resampled; j++)
    {
      field[i][j] = new aladyn_float[params->npx_resampled];
      allocated_size += (params->npx_resampled)*sizeof(aladyn_float);
    }
#ifdef ENABLE_DEBUG
    printf("Allocated %llu bytes for fields\r", (unsigned long long int) allocated_size);
#endif
  }
  x_lineout = new aladyn_float[params->npx_resampled];

  printf("READING\n");
  if (!params->multifile)
  {
    bin_filename.str("");
    bin_filename << params->filebasename << ".bin";
    file_in = fopen(bin_filename.str().c_str(), "rb");
    if (file_in == NULL) std::cout << "Unable to open file!" << std::endl;
    else std::cout << "File opened to read data!" << std::endl;

    /*skip header*/
    std::fseek(file_in, (long)params->header_size_bytes, SEEK_SET);
    std::cout << "Fseek of " << params->header_size_bytes << " bytes from beginning of file done" << std::endl << std::flush;

    for (unsigned int ipx = 0; ipx < params->ncpu_x; ipx++)
    {
      for (unsigned int ipz = 0; ipz < params->ncpu_z; ipz++)
      {
        for (unsigned int ipy = 0; ipy < params->ncpu_y; ipy++)
        {
          if (params->file_version == 1 || params->file_version == 2) fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);
          fread_size += sizeof(int)*std::fread(header, sizeof(int), header_size, file_in);
          if (params->file_version == 1 || params->file_version == 2) fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);

          if (params->we_have_to_do_swap) swap_endian_i(header, header_size);

          if (header[0] != params->npx_resampled_per_cpu ||
            header[1] != params->npy_resampled_per_cpu ||
            (header[2] != params->npz_resampled_per_cpu && header[2] != 1 && params->npz_resampled_per_cpu != 0)) // fix for 2D files, which have 1 point in z for parameters but none for the binary dump
          {
            printf("\nWARNING: unexpected number of points in this chunk!\n");
            printf("header[] = [%i,%i,%i], parameters[] = [%llu,%llu,%llu]\n", header[0], header[1], header[2], (unsigned long long int) params->npx_resampled_per_cpu, (unsigned long long int) params->npy_resampled_per_cpu, (unsigned long long int) params->npz_resampled_per_cpu);
          }

          printf("header[] = {%i/%llu, %i/%llu, %i/%llu}, cpu[] = {%u/%u, %u/%u, %u/%u}\r", header[0], (unsigned long long int) params->npx_resampled_per_cpu, header[1], (unsigned long long int) params->npy_resampled_per_cpu, header[2], (unsigned long long int) params->npz_resampled_per_cpu, ipx + 1, params->ncpu_x, ipy + 1, params->ncpu_y, ipz + 1, params->ncpu_z);
          fflush(stdout);

          buffer = new aladyn_float[header[0] * header[1] * header[2]];
          if (params->file_version == 1 || params->file_version == 2) fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);
          fread_size += sizeof(aladyn_float)*std::fread(buffer, sizeof(aladyn_float), header[0] * header[1] * header[2], file_in);
          if (params->file_version == 1 || params->file_version == 2) fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);

          if (params->we_have_to_do_swap) swap_endian_f(buffer, header[0] * header[1] * header[2]);

          for (size_t k = 0; k < header[2]; k++)
            for (size_t j = 0; j < header[1]; j++)
              for (size_t i = 0; i < header[0]; i++)
                field[k + (ipz * params->npz_resampled_per_cpu)]
                /* */[j + (ipy * params->npy_resampled_per_cpu)]
          /*       */[i + (ipx * params->npx_resampled_per_cpu)] =
            /*     */ buffer[i + j*header[0] + k*header[0] * header[1]];

          delete[] buffer;
          buffer = NULL;
        }
      }
    }

    // leggiamo ora le coordinate dei punti di griglia, presenti solo nelle versioni che possono prevedere griglia stretchata e che ancora non la scrivevano nel .dat
    // se presenti, sovrascrivono quelle lette o precostruite (se non trovate nel file .dat) dalle routine dei params

    fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);  // facciamo il test sul buffer Fortran della prima coordinata;
    // se esiste, non e' necessario tornare indietro perche' il buffer fortran che precede i data non e' di alcun interesse

    if (!std::feof(file_in))
    {
      aladyn_float *x_coordinates, *y_coordinates, *z_coordinates;
      x_coordinates = new aladyn_float[params->npx_resampled];
      y_coordinates = new aladyn_float[params->npy_resampled];
      z_coordinates = new aladyn_float[params->npz_resampled];
      fread_size += sizeof(aladyn_float)*std::fread(x_coordinates, sizeof(aladyn_float), params->npx_resampled, file_in);
      fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);
      fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);
      fread_size += sizeof(aladyn_float)*std::fread(y_coordinates, sizeof(aladyn_float), params->npy_resampled, file_in);
      fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);
      fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);
      fread_size += sizeof(aladyn_float)*std::fread(z_coordinates, sizeof(aladyn_float), params->npz_resampled, file_in);
      fread_size += sizeof(int)*std::fread(&fortran_buff, sizeof(int), 1, file_in);

      if (params->we_have_to_do_swap)
      {
        swap_endian_f(x_coordinates, params->npx_resampled);
        swap_endian_f(y_coordinates, params->npy_resampled);
        swap_endian_f(z_coordinates, params->npz_resampled);
      }

      params->xcoord.resize(params->npx_resampled, 0);
      params->ycoord.resize(params->npy_resampled, 0);
      params->zcoord.resize(params->npz_resampled, 0);

      for (size_t i = 0; i < params->npx_resampled; i++)
        params->xcoord[i] = x_coordinates[i];
      for (size_t i = 0; i < params->npy_resampled; i++)
        params->ycoord[i] = y_coordinates[i];
      for (size_t i = 0; i < params->npz_resampled; i++)
        params->zcoord[i] = z_coordinates[i];
    }
    else params->stretched_grid = false;
  }
  else
  {
    int header_size = 3;
    int * header = new int[header_size];
    while (1)
    {
      bin_filename.str("");
      bin_filename << params->filebasename << "_" << std::setfill('0') << std::setw(3) << multifile_index << ".bin";
      file_in = fopen(bin_filename.str().c_str(), "rb");
      if (file_in == NULL)
      {
        std::cout << "End of files!" << std::endl;
        break;
      }
      else std::cout << "Opened file #" << multifile_index << " to read data!" << std::endl;

      /*skip header*/
      std::fseek(file_in, (long)params->header_size_bytes, SEEK_SET);
      std::cout << "Fseek of " << params->header_size_bytes << " bytes from beginning of file done" << std::endl << std::flush;

      for (unsigned int ipx = 0; ipx < params->ncpu_x; ipx++)
      {
        for (unsigned int ipz = 0; ipz < params->ncpu_z; ipz++)
        {
          for (unsigned int ipy = 0; ipy < params->ncpu_y; ipy++)
          {
            fread_size += sizeof(int)*std::fread(header, sizeof(int), header_size, file_in);
            if (params->we_have_to_do_swap) swap_endian_i(header, header_size);

#ifdef ENABLE_DEBUG
            if (header[0] != params->npx_resampled_per_cpu ||
              header[1] != params->npy_resampled_per_cpu ||
              header[2] != params->npz_resampled_per_cpu)
              std::cout << "WARNING: unexpected number of points in this chunk!" << std::endl << std::flush;
#endif

            printf("file %i, header[] = {%i/%llu, %i/%llu, %i/%llu}, cpu[] = {%u/%u, %u/%u, %u/%u}\r", multifile_index, header[0], (unsigned long long int) params->npx_resampled_per_cpu, header[1], (unsigned long long int) params->npy_resampled_per_cpu, header[2], (unsigned long long int) params->npz_resampled_per_cpu, ipx + 1, params->ncpu_x, ipy + 1, params->ncpu_y, ipz + 1, params->ncpu_z);
            fflush(stdout);

            buffer = new aladyn_float[header[0] * header[1] * header[2]];
            fread_size += sizeof(aladyn_float)*std::fread(buffer, sizeof(aladyn_float), header[0] * header[1] * header[2], file_in);

            if (params->we_have_to_do_swap) swap_endian_f(buffer, params->npx_resampled_per_cpu*params->npy_resampled_per_cpu*params->npz_resampled_per_cpu);

            for (size_t k = 0; k < header[2]; k++)
              for (size_t j = 0; j < header[1]; j++)
                for (size_t i = 0; i < header[0]; i++)
                  field[k + (ipz * params->npz_resampled_per_cpu)]
                  /* */[j + (ipy * params->npy_resampled_per_cpu)]
            /*       */[i + (ipx * params->npx_resampled_per_cpu)] =
              /*     */ buffer[i + j*header[0] + k*header[0] * header[1]];
            delete[] buffer;
            buffer = NULL;
          }
        }
      }
      multifile_index++;
      fclose(file_in);
    }
  }

  std::cout << std::endl << "END READING: TOTAL " << fread_size / (1024. * 1024.) << " MB READ" << std::endl << std::flush;

  if (params->npz_resampled == 1 && params->out_grid2d)
  {
    printf("Writing the ASCII 2D field file\n");
    grid_filename = params->filebasename + ".txt";
    clean_fields = std::fopen(grid_filename.c_str(), "wb"); // binary has just the meaning here not to convert \n to \r\n in Windows. The file is ASCII anyway ;)

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

    fprintf(clean_fields, "# %llu \n # %llu \n # %llu\n# %g  %g \n # %g  %g\n", (unsigned long long int) params->npx_resampled,
      (unsigned long long int) params->npy_resampled, (unsigned long long int) params->npz_resampled, params->xmin, params->ymin, params->xmax, params->ymax);
    for (size_t j = 0; j < params->npy_resampled; j++)
    {
      for (size_t i = 0; i < params->npx_resampled; i++)
      {
        fprintf(clean_fields, "%.4g %.4g %.4g\n", params->xcoord[i], params->ycoord[j], field[0][j][i]);
      }
      printf("grid[y] = {%llu/%llu}\r", (unsigned long long int) j, (unsigned long long int) params->npy_resampled);
    }
    fclose(clean_fields);
    printf("\nASCII FILE COMPLETED\n");

#if defined(FORCE_PRINTF_BUFFER_SIZE) && (FORCE_PRINTF_BUFFER_SIZE > 0)
    delete[] buffer_out;
#endif
  }


  if (params->npz_resampled == 1 && params->out_lineoutx)
  {
    size_t myj;
    printf("Writing 1D lineout along x\n");
    grid_filename = params->filebasename + "_lineout.txt";
    clean_fields = std::fopen(grid_filename.c_str(), "wb"); // binary has just the meaning here not to convert \n to \r\n in Windows. The file is ASCII anyway ;)

    for (size_t i = 0; i < params->npx_resampled; i++) x_lineout[i] = 0;

    for (size_t j = 0; j < params->npy_resampled; j++)
    {
      if (params->ycoord[j] >= 0)
      {
        myj = j;
        break;
      }
    }

    fprintf(clean_fields, "#");
    for (size_t j = myj - params->span; j < (myj + params->span + 1); j++)
    {
      fprintf(clean_fields, "%.4g\t", params->ycoord[j]);
      for (size_t i = 0; i < params->npx_resampled; i++)
      {
        x_lineout[i] += field[0][j][i] / (2.0f * params->span + 1.0f);
      }
    }
    fprintf(clean_fields, "\n");
    for (size_t i = 0; i < params->npx_resampled; i++)
    {
      fprintf(clean_fields, "%.4g %.4g\n", params->xcoord[i], x_lineout[i]);
    }
    fclose(clean_fields);
  }

  if (params->npz_resampled > 1 && params->out_lineoutx)
  {
    size_t myj;
    printf("Writing 1D lineout along x\n");
    grid_filename = params->filebasename + "_lineout.txt";
    clean_fields = std::fopen(grid_filename.c_str(), "wb"); // binary has just the meaning here not to convert \n to \r\n in Windows. The file is ASCII anyway ;)

    for (size_t i = 0; i < params->npx_resampled; i++)
      x_lineout[i] = 0;
    for (size_t j = 0; j < params->npy_resampled; j++)
    {
      if (params->ycoord[j] >= 0)
      {
        myj = j;
        break;
      }
    }
    fprintf(clean_fields, "#");
    for (size_t k = myj - params->span; k < (myj + params->span + 1); k++)
    {
      fprintf(clean_fields, "%.4g\t", params->zcoord[k]);

      for (size_t j = myj - params->span; j < (myj + params->span + 1); j++)
      {
        for (size_t i = 0; i < params->npx_resampled; i++)
        {
          x_lineout[i] += field[k][j][i] / ((2.0f * params->span + 1.0f)*(2.0f * params->span + 1.0f));
        }
      }
    }
    fprintf(clean_fields, "\n");
    for (size_t i = 0; i < params->npx_resampled; i++)
    {
      fprintf(clean_fields, "%.4g %.4g\n", params->xcoord[i], x_lineout[i]);
    }
    fclose(clean_fields);
  }


  if (params->npz_resampled > 1 && params->out_cutz)
  {
    if (params->where_to_cut_grid_along_z.size() == 0)
    {
      cut = params->zcoord[params->npz_resampled / 2];
      cutz.push_back(cut);
      gridIndex_cutz.push_back(params->npz_resampled / 2);
    }
    else
    {
      cut = params->zmin;
      int i = 0;
      for (size_t j = 0; j < params->where_to_cut_grid_along_z.size(); j++)
      {
        while (cut < params->where_to_cut_grid_along_z.at(j)) cut = params->zcoord[i], i++;
        cutz.push_back(cut);
        gridIndex_cutz.push_back(i - 1);
        i = 0;
      }
    }
    for (size_t n = 0; n < cutz.size(); n++)
    {
      grid_filename = params->filebasename + "_cutz_" + std::to_string(cutz[n]) + ".txt";
      printf("Writing the 2D ASCII field file at z=%g\n", cutz[n]);
      clean_fields = std::fopen(grid_filename.c_str(), "wb"); // binary has just the meaning here not to convert \n to \r\n in Windows. The file is ASCII anyway ;)

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
      fprintf(clean_fields, "# %lu\n#%lu\n#%i\n", (long)params->npx_resampled, (long)params->npy_resampled, 1);
      fprintf(clean_fields, "#%f %f\n#%f %f\n", params->xmin, params->ymin, params->xmax, params->ymax);
      size_t k = gridIndex_cutz[n];
      for (size_t j = 0; j < params->npy_resampled; j++)
      {
        for (size_t i = 0; i < params->npx_resampled; i++)
        {
          fprintf(clean_fields, "%.4g %.4g %.4g\n", params->xcoord[i], params->ycoord[j], field[k][j][i]);
        }
        printf("grid[y] = {%llu/%llu}\r", (unsigned long long int) j, (unsigned long long int) params->npy_resampled);
      }
      fclose(clean_fields);
      printf("\nASCII FILE COMPLETED\n");

#if defined(FORCE_PRINTF_BUFFER_SIZE) && (FORCE_PRINTF_BUFFER_SIZE > 0)
      delete[] buffer_out;
#endif
    }
  }

  if (params->npz_resampled > 1 && params->out_cuty)
  {
    if (params->where_to_cut_grid_along_y.size() == 0)
    {
      cut = params->ycoord[params->npy_resampled / 2];
      cuty.push_back(cut);
      gridIndex_cuty.push_back(params->npy_resampled / 2);
    }
    else
    {
      cut = params->ymin;
      int i = 0;
      for (size_t j = 0; j < params->where_to_cut_grid_along_y.size(); j++)
      {
        while (cut < params->where_to_cut_grid_along_y.at(j)) cut = params->ycoord[i], i++;
        cuty.push_back(cut);
        gridIndex_cuty.push_back(i - 1);
        i = 0;
      }
    }
    for (size_t n = 0; n < cuty.size(); n++)
    {
      grid_filename = params->filebasename + "_cuty_" + std::to_string(cuty[n]) + ".txt";
      printf("Writing the 2D ASCII field file at y=%g\n", cuty[n]);
      clean_fields = std::fopen(grid_filename.c_str(), "wb"); // binary has just the meaning here not to convert \n to \r\n in Windows. The file is ASCII anyway ;)

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
      fprintf(clean_fields, "# %lu\n#%lu\n#%i\n", (long)params->npx_resampled, (long)params->npz_resampled, 1);
      fprintf(clean_fields, "#%f %f\n#%f %f\n", params->xmin, params->zmin, params->xmax, params->zmax);
      size_t j = gridIndex_cuty[n];
      for (size_t k = 0; k < params->npz_resampled; k++)
      {
        for (size_t i = 0; i < params->npx_resampled; i++)
        {
          fprintf(clean_fields, "%.4g %.4g %.4g\n", params->xcoord[i], params->zcoord[k], field[k][j][i]);
        }
        printf("grid[z] = {%llu/%llu}\r", (unsigned long long int) k, (unsigned long long int) params->npz_resampled);
      }
      fclose(clean_fields);
      printf("\nASCII FILE COMPLETED\n");

#if defined(FORCE_PRINTF_BUFFER_SIZE) && (FORCE_PRINTF_BUFFER_SIZE > 0)
      delete[] buffer_out;
#endif
    }
  }

  if (params->npz_resampled > 1 && params->out_cutx)
  {
    if (params->where_to_cut_grid_along_x.size() == 0)
    {
      cut = params->xcoord[params->npx_resampled / 2];
      cutx.push_back(cut);
      gridIndex_cutx.push_back(params->npx_resampled / 2);
    }
    else
    {
      cut = params->xmin;
      int i = 0;
      for (size_t j = 0; j < params->where_to_cut_grid_along_x.size(); j++)
      {
        while (cut < params->where_to_cut_grid_along_x.at(j)) cut = params->xcoord[i], i++;
        cutx.push_back(cut);
        gridIndex_cutx.push_back(i - 1);
        i = 0;
      }
    }
    for (size_t n = 0; n < cutx.size(); n++)
    {
      grid_filename = params->filebasename + "_cutx_" + std::to_string(cutx[n]) + ".txt";
      printf("Writing the 2D ASCII field file at x=%g\n", cutx[n]);
      clean_fields = std::fopen(grid_filename.c_str(), "wb"); // binary has just the meaning here not to convert \n to \r\n in Windows. The file is ASCII anyway ;)

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
      fprintf(clean_fields, "# %lu\n#%lu\n#%i\n", (long)params->npy_resampled, (long)params->npz_resampled, 1);
      fprintf(clean_fields, "#%f %f\n#%f %f\n", params->ymin, params->zmin, params->ymax, params->zmax);
      size_t i = gridIndex_cutx[n];
      for (size_t k = 0; k < params->npz_resampled; k++)
      {
        for (size_t j = 0; j < params->npy_resampled; j++)
        {
          fprintf(clean_fields, "%.4g %.4g %.4g\n", params->ycoord[j], params->zcoord[k], field[k][j][i]);
        }
        printf("grid[z] = {%llu/%llu}\r", (unsigned long long int) k, (unsigned long long int) params->npz_resampled);
      }
      fclose(clean_fields);
      printf("\nASCII FILE COMPLETED\n");

#if defined(FORCE_PRINTF_BUFFER_SIZE) && (FORCE_PRINTF_BUFFER_SIZE > 0)
      delete[] buffer_out;
#endif
    }
  }


  if (params->out_vtk)
  {
    printf("Writing the vtk fields file\n");

    aladyn_float *x_coordinates, *y_coordinates, *z_coordinates;
    x_coordinates = new aladyn_float[params->npx_resampled];
    y_coordinates = new aladyn_float[params->npy_resampled];
    z_coordinates = new aladyn_float[params->npz_resampled];

    for (size_t i = 0; i < params->npx_resampled; i++)
      x_coordinates[i] = params->xcoord[i];
    for (size_t i = 0; i < params->npy_resampled; i++)
      y_coordinates[i] = params->ycoord[i];
    for (size_t i = 0; i < params->npz_resampled; i++)
      z_coordinates[i] = params->zcoord[i];

    if (params->endian_machine == 0)
    {
      swap_endian_f(field, params->npz_resampled, params->npy_resampled, params->npx_resampled);
      swap_endian_f(x_coordinates, params->npx_resampled);
      swap_endian_f(y_coordinates, params->npy_resampled);
      swap_endian_f(z_coordinates, params->npz_resampled);
    }

    //////// DATASET UNSTRUCTURED_GRID VERSION    ////////

    grid_filename = params->filebasename + ".vtk";
    clean_fields = std::fopen(grid_filename.c_str(), "wb");
    fprintf(clean_fields, "# vtk DataFile Version 2.0\n");
    fprintf(clean_fields, "titolo mio\n");
    fprintf(clean_fields, "BINARY\n");
    fprintf(clean_fields, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(clean_fields, "POINTS %lu aladyn_float\n", (long)(params->npx_resampled*params->npy_resampled*params->npz_resampled));
    aladyn_float rr[3];
    for (size_t k = 0; k < params->npz_resampled; k++)
    {
      rr[2] = z_coordinates[k];
      for (size_t j = 0; j < params->npy_resampled; j++)
      {
        rr[1] = y_coordinates[j];
        for (size_t i = 0; i < params->npx_resampled; i++)
        {
          rr[0] = x_coordinates[i];
          fwrite((void*)rr, sizeof(aladyn_float), 3, clean_fields);
        }
      }
    }

    fprintf(clean_fields, "POINT_DATA %lu\n", (long)(params->npx_resampled*params->npy_resampled*params->npz_resampled));
    fprintf(clean_fields, "SCALARS %s aladyn_float 1\n", params->support_label.c_str());
    fprintf(clean_fields, "LOOKUP_TABLE default\n");
    fwrite((void*)field, sizeof(aladyn_float), params->npz_resampled*params->npy_resampled*params->npx_resampled, clean_fields);
    fclose(clean_fields);
  }



  if (params->out_vtk_nostretch)
  {
    printf("Writing the vtk file related to the unstretched part of the grid\n");

    int inizio_punti_non_stretchati_x, inizio_punti_non_stretchati_y, inizio_punti_non_stretchati_z;
    int fine_punti_non_stretchati_x, fine_punti_non_stretchati_y, fine_punti_non_stretchati_z;
    size_t npunti_non_stretchati_x, npunti_non_stretchati_y, npunti_non_stretchati_z;
    if (params->stretched_grid)
    {
      if (params->stretched_along_x) inizio_punti_non_stretchati_x = (int)(params->npx_resampled / 6.0);
      else inizio_punti_non_stretchati_x = 0;
      inizio_punti_non_stretchati_y = (int)(params->npy_resampled / 6.0);
      inizio_punti_non_stretchati_z = (int)(params->npz_resampled / 6.0);
      if (params->stretched_along_x) fine_punti_non_stretchati_x = (int)(params->npx_resampled*5.0 / 6.0);
      else fine_punti_non_stretchati_x = (int)(params->npx_resampled);
      fine_punti_non_stretchati_y = (int)(params->npy_resampled*5.0 / 6.0);
      fine_punti_non_stretchati_z = (int)(params->npz_resampled*5.0 / 6.0);
      npunti_non_stretchati_x = (size_t)(fine_punti_non_stretchati_x - inizio_punti_non_stretchati_x);
      npunti_non_stretchati_y = (size_t)(fine_punti_non_stretchati_y - inizio_punti_non_stretchati_y);
      npunti_non_stretchati_z = (size_t)(fine_punti_non_stretchati_z - inizio_punti_non_stretchati_z);
    }
    else
    {
      fine_punti_non_stretchati_x = (int)params->npx_resampled;
      npunti_non_stretchati_x = params->npx_resampled;
      fine_punti_non_stretchati_y = (int)params->npy_resampled;
      npunti_non_stretchati_y = params->npy_resampled;
      fine_punti_non_stretchati_z = (int)params->npz_resampled;
      npunti_non_stretchati_z = params->npz_resampled;
      inizio_punti_non_stretchati_x = inizio_punti_non_stretchati_y = inizio_punti_non_stretchati_z = 0;
    }

    aladyn_float * field_non_stretchato = new aladyn_float[npunti_non_stretchati_x*npunti_non_stretchati_y*npunti_non_stretchati_z];
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

    if (params->endian_machine == 0) swap_endian_f(field_non_stretchato, npunti_non_stretchati_x*npunti_non_stretchati_y*npunti_non_stretchati_z);



    //////// DATASET STRUCTURED_POINTS VERSION    ////////

    aladyn_float xmin_non_stretchato = params->xcoord[inizio_punti_non_stretchati_x];
    //    aladyn_float xmax_non_stretchato = params->xcoord[fine_punti_non_stretchati_x];
    aladyn_float ymin_non_stretchato = params->ycoord[inizio_punti_non_stretchati_y];
    //    aladyn_float ymax_non_stretchato = params->ycoord[fine_punti_non_stretchati_y];
    aladyn_float zmin_non_stretchato = params->zcoord[inizio_punti_non_stretchati_z];
    //    aladyn_float zmax_non_stretchato = params->zcoord[fine_punti_non_stretchati_z];

    grid_filename = params->filebasename + "_nostretch.vtk";
    clean_fields = std::fopen(grid_filename.c_str(), "wb");
    fprintf(clean_fields, "# vtk DataFile Version 2.0\n");
    fprintf(clean_fields, "titolo mio\n");
    fprintf(clean_fields, "BINARY\n");
    fprintf(clean_fields, "DATASET STRUCTURED_POINTS\n");
    fprintf(clean_fields, "DIMENSIONS %lu %lu %lu\n", (long)npunti_non_stretchati_x, (long)npunti_non_stretchati_y, (long)npunti_non_stretchati_z);
    fprintf(clean_fields, "ORIGIN %f %f %f\n", xmin_non_stretchato, ymin_non_stretchato, zmin_non_stretchato);
    fprintf(clean_fields, "SPACING %f %f %f\n",
      params->xcoord[params->npx_resampled / 2] - params->xcoord[params->npx_resampled / 2 - 1],
      params->ycoord[params->npy_resampled / 2] - params->ycoord[params->npy_resampled / 2 - 1],
      params->zcoord[params->npz_resampled / 2] - params->zcoord[params->npz_resampled / 2 - 1]);
    fprintf(clean_fields, "POINT_DATA %lu\n", (long)(npunti_non_stretchati_x*npunti_non_stretchati_y*npunti_non_stretchati_z));
    fprintf(clean_fields, "SCALARS %s aladyn_float 1\n", params->support_label.c_str());
    fprintf(clean_fields, "LOOKUP_TABLE default\n");
    fwrite((void*)field_non_stretchato, sizeof(aladyn_float), npunti_non_stretchati_x*npunti_non_stretchati_y*npunti_non_stretchati_z, clean_fields);




    /******************************************************************************
    //////// DATASET RECTILINEAR_GRID VERSION    ////////

    aladyn_float *x_coordinates, *y_coordinates, *z_coordinates;
    x_coordinates=new aladyn_float[npunti_non_stretchati_x];
    y_coordinates=new aladyn_float[npunti_non_stretchati_y];
    z_coordinates=new aladyn_float[npunti_non_stretchati_z];
    for (int i = inizio_punti_non_stretchati_x; i < fine_punti_non_stretchati_x; i++) x_coordinates[i] = params->xcoord[i];
    for (int i = inizio_punti_non_stretchati_y; i < fine_punti_non_stretchati_y; i++) y_coordinates[i] = params->ycoord[i];
    for (int i = inizio_punti_non_stretchati_z; i < fine_punti_non_stretchati_z; i++) z_coordinates[i] = params->zcoord[i];
    if(params->endian_machine == 0)
    {
    swap_endian_f(x_coordinates,npunti_non_stretchati_x);
    swap_endian_f(y_coordinates,npunti_non_stretchati_y);
    swap_endian_f(z_coordinates,npunti_non_stretchati_z);
    }

    sprintf(grid_filename,"%s_out.vtk",params->filebasename);
    clean_fields=fopen(grid_filename, "wb");
    printf("\nWriting the fields file\n");
    fprintf(clean_fields,"# vtk DataFile Version 2.0\n");
    fprintf(clean_fields,"titolo mio\n");
    fprintf(clean_fields,"BINARY\n");
    fprintf(clean_fields,"DATASET RECTILINEAR_GRID\n");
    fprintf(clean_fields,"DIMENSIONS %i %i %i\n",npunti_non_stretchati_x, npunti_non_stretchati_y, npunti_non_stretchati_z);
    fprintf(clean_fields,"X_COORDINATES %i aladyn_float\n",npunti_non_stretchati_x);
    fwrite((void*)x_coordinates,sizeof(aladyn_float),npunti_non_stretchati_x,clean_fields);
    fprintf(clean_fields,"Y_COORDINATES %i aladyn_float\n",npunti_non_stretchati_y);
    fwrite((void*)y_coordinates,sizeof(aladyn_float),npunti_non_stretchati_y,clean_fields);
    fprintf(clean_fields,"Z_COORDINATES %i aladyn_float\n",npunti_non_stretchati_z);
    fwrite((void*)z_coordinates,sizeof(aladyn_float),npunti_non_stretchati_z,clean_fields);
    fprintf(clean_fields,"POINT_DATA %i\n",npunti_non_stretchati_x*npunti_non_stretchati_y*npunti_non_stretchati_z);
    fprintf(clean_fields,"SCALARS %s aladyn_float 1\n",params->support_label);
    fprintf(clean_fields,"LOOKUP_TABLE default\n");
    fwrite((void*)field_non_stretchato,sizeof(aladyn_float),npunti_non_stretchati_x*npunti_non_stretchati_y*npunti_non_stretchati_z,clean_fields);
    ******************************************************************************/

    fclose(clean_fields);
  }

  if (params->out_params)
  {
    params_filename = params->filebasename + ".parameters";
    parameters = std::fopen(params_filename.c_str(), "wb");
    printf("Writing parameters to file\n");
    fprintf(parameters, "ncpu_x=%u\n", params->ncpu_x);
    fprintf(parameters, "ncpu_y=%u\n", params->ncpu_y);
    fprintf(parameters, "ncpu_z=%u\n", params->ncpu_z);
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

    fprintf(parameters, "\n\nGrid along x axis\n");
    for (unsigned int i = 0; i < params->npx_resampled; i++)
    {
      fprintf(parameters, "%.4g  ", params->xcoord[i]);
      if (i > 0 && i % 10 == 0) fprintf(parameters, "\n");
    }

    fprintf(parameters, "\n\nGrid along y axis\n");
    for (unsigned int i = 0; i < params->npy_resampled; i++)
    {
      fprintf(parameters, "%.4g  ", params->ycoord[i]);
      if (i > 0 && i % 10 == 0) fprintf(parameters, "\n");
    }

    fprintf(parameters, "\n\nGrid along z axis\n");
    for (unsigned int i = 0; i < params->npz_resampled; i++)
    {
      fprintf(parameters, "%.4g  ", params->zcoord[i]);
      if (i > 0 && i % 10 == 0) fprintf(parameters, "\n");
    }

    fclose(parameters);
  }


  return 0;

}

