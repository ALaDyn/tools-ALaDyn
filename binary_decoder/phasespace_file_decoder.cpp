
#include "phasespace_file_decoder.h"


int read_phasespace_file(Parameters * params)
{
  size_t phasespace_size = 14; // x, y, z, px, py, pz, gamma, theta, thetaT, E, ty, tz, w, ch
  std::FILE *file_in = nullptr;
  int multifile_index = 0;
  int contatori[] = { 0, 0, 0 };
  aladyn_float zero = 0.0f;
  long long parts_accumulate = 0;
  double peso_accumulato = 0.0;
  double carica_accumulata = 0.0;
  size_t dim_file_in_bytes = 0, num_of_floats_in_file = 0, num_of_particles_in_file = 0, num_of_passes = 0, num_residual_particles = 0;
  long long dimensione_array_parts = 0;
  unsigned int val[2] = { 0, 0 };
  aladyn_float array_supporto8[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  aladyn_float array_supporto6[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  const size_t max_number_of_particles_in_memory = 10000000;          // in reality we store double this number -1

  std::FILE *binary_vtk = nullptr;
  std::FILE *binary_clean = nullptr;
  std::FILE *ascii_propaga = nullptr;
  std::FILE *ascii_xyze = nullptr;
  std::FILE *ascii_csv = nullptr;
  std::FILE *parameters = nullptr;
  std::ofstream Estremi_out;
  unsigned int conta_processori = 0;

  int npart_loc = 0;
  int buff = 0;

  aladyn_float x = 0.0, y = 0.0, z = 0.0, px = 0.0, py = 0.0, pz = 0.0, ptot = 0.0;
  aladyn_float ch = 0.0, w = 0.0;
  aladyn_float gamma = 0.0, theta = 0.0, thetaT = 0.0, E = 0.0, ty = 0.0, tz = 0.0;
  double *estremi_min = nullptr, *estremi_max = nullptr;

  short buffshort[2] = { 0, 0 };
  aladyn_float *parts = nullptr;
  std::string bin_filename;
  std::string dat_filename;
  std::string vtk_filename;
  std::string clean_bin_filename;
  std::string ppg_filename;
  std::string xyze_filename;
  std::string csv_filename;
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

  ppg_filename = params->filebasename + ".ppg";
  xyze_filename = params->filebasename + "_xyzE.ppg";
  csv_filename = params->filebasename + ".csv";
  vtk_filename = params->filebasename + ".vtk";
  clean_bin_filename = params->filebasename + "_clean.bin";
  bin_filename = params->filebasename + ".bin";
  dat_filename = params->filebasename + ".dat";

  if (!params->multifile) file_in = fopen(bin_filename.c_str(), "rb");

  if (params->out_vtk)
  {
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

  if (params->out_clean_bin) {
    binary_clean = fopen(clean_bin_filename.c_str(), "wb");
  }

  if (params->out_ppg) {
    ascii_propaga = fopen(ppg_filename.c_str(), "wb");
  }

  if (params->out_xyze) {
    ascii_xyze = fopen(xyze_filename.c_str(), "wb");
  }

  if (params->out_csv) {
    ascii_csv = fopen(csv_filename.c_str(), "wb");
  }

  while (1)
  {
    if (!params->multifile)
    {
      if (conta_processori >= params->ncpu) break;
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
      printf("proc number \t %i/%u \t npart=%i \n", conta_processori + 1, params->ncpu, npart_loc);
#else
      printf("proc number \t %i/%u \t npart=%i \r", conta_processori + 1, params->ncpu, npart_loc);
#endif
      fflush(stdout);
      num_of_passes = 1;
    }
    else  //we do have multifiles i.e. Prpout00_000.bin
    {
      memset(&bin_filename[0], 0, sizeof(bin_filename));
      bin_filename = params->filebasename + "_" + std::to_string(multifile_index) + ".bin";

      if ((file_in = fopen(bin_filename.c_str(), "rb")) == nullptr)
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
          fread_size = std::fread(parts, sizeof(aladyn_float), val[0] * params->ndv, file_in);
          if (params->we_have_to_do_swap) swap_endian_f(parts, (size_t)val[0] * params->ndv);
        }

        _Filter(params, parts, val, _Filter::build_filter(params));

        for (auto histogram : params->histograms) {
          if (histogram.enabled) _Binning(parts, val[0], params, &histogram);
        }
        for (auto densityplot : params->densityplots) {
          if (densityplot.enabled) _Binning(parts, val[0], params, &densityplot);
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
              y = 0.0;
              z = parts[i*params->ndv + 0] * ((aladyn_float)1.e-4);
              px = parts[i*params->ndv + 3];
              py = 0.0;
              pz = parts[i*params->ndv + 2];
              gamma = (aladyn_float)(sqrt(1. + px*px + py*py + pz*pz) - 1.);
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
        parts = nullptr;
      }
    }
    multifile_index++;
    conta_processori++;
  }

  for (auto histogram : params->histograms) {
    if (histogram.enabled) histogram.write_binned_data();
  }

  for (auto densityplot : params->densityplots) {
    if (densityplot.enabled) densityplot.write_binned_data();
  }

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

  std::cout << std::endl << "fread_size = " << fread_size << std::endl << "END!" << std::endl;

  return 0;
}



int create_json_from_phasespace_file(Parameters * params)
{
  std::FILE *file_in = nullptr;
  int multifile_index = 0;
  int contatori[] = { 0, 0, 0 };
  aladyn_float zero = 0.0f;
  long long parts_accumulate = 0;
  double peso_accumulato = 0.0;
  double carica_accumulata = 0.0;
  size_t dim_file_in_bytes = 0, num_of_floats_in_file = 0, num_of_particles_in_file = 0, num_of_passes = 0, num_residual_particles = 0;
  long long dimensione_array_parts = 0;
  unsigned int val[2] = { 0, 0 };
  const size_t max_number_of_particles_in_memory = 10000000;          // in reality we store double this number -1

  unsigned int conta_processori = 0;

  int npart_loc = 0;
  int buff = 0;
  aladyn_float x = 0.0, y = 0.0, z = 0.0, px = 0.0, py = 0.0, pz = 0.0, ptot = 0.0;
  aladyn_float ch = 0.0, w = 0.0;
  aladyn_float gamma = 0.0, theta = 0.0, thetaT = 0.0, E = 0.0, ty = 0.0, tz = 0.0;
  double emittance_x = 0.0, emittance_y = 0.0, emittance_z = 0.0;
  double em_x2 = 0.0, em_x = 0.0, em_y2 = 0.0, em_y = 0.0, em_z2 = 0.0, em_z = 0.0;
  double em_px2 = 0.0, em_px = 0.0, em_py2 = 0.0, em_py = 0.0, em_pz2 = 0.0, em_pz = 0.0, em_xpx = 0.0, em_ypy = 0.0, em_zpz = 0.0;
  short buffshort[2] = { 0, 0 };
  aladyn_float *parts = nullptr;

  size_t fread_size = 0;
  double *estremi_min, *estremi_max;
  estremi_min = new double[14];
  estremi_max = new double[14];
  for (size_t i = 0; i < 14; i++)
  {
    estremi_min[i] = std::numeric_limits<double>::max();
    estremi_max[i] = std::numeric_limits<double>::lowest();
  }

  std::string bin_filename, json_filename;
  json_filename = params->filebasename + ".json";
  bin_filename = params->filebasename + ".bin";

  if (!params->multifile) file_in = fopen(bin_filename.c_str(), "rb");

  while (1)
  {
    if (!params->multifile)
    {
      if (conta_processori >= params->ncpu) break;
      if (conta_processori == 0) {
        /*we are at the beginning of the file, so skip header*/
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
      printf("proc number \t %i/%u \t npart=%i \n", conta_processori + 1, params->ncpu, npart_loc);
#else
      printf("proc number \t %i/%u \t npart=%i \r", conta_processori + 1, params->ncpu, npart_loc);
#endif
      num_of_passes = 1;
    }
    else  //we do have multifiles i.e. Prpout00_000.bin
    {
      memset(&bin_filename[0], 0, sizeof(bin_filename));
      bin_filename = params->filebasename + "_" + std::to_string(multifile_index) + ".bin";

      if ((file_in = fopen(bin_filename.c_str(), "rb")) == nullptr)
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
          fread_size = std::fread(parts, sizeof(aladyn_float), val[0] * params->ndv, file_in);
          if (params->we_have_to_do_swap) swap_endian_f(parts, (size_t)val[0] * params->ndv);
        }

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
        parts_accumulate += val[0];
        delete[] parts;
        parts = nullptr;
      }
    }
    multifile_index++;
    conta_processori++;
  }


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

  std::ofstream json_file;
  json_file.open(json_filename);
  if (json_file.fail()) std::cerr << "Unable to create .json file!" << std::endl;

  jsoncons::json parameters;

  parameters["ncpu_x"] = params->ncpu_x;
  parameters["ncpu_y"] = params->ncpu_y;
  parameters["ncpu_z"] = params->ncpu_z;
  parameters["npx_resampled"] = params->npx_resampled;
  parameters["npy_resampled"] = params->npy_resampled;
  parameters["npz_resampled"] = params->npz_resampled;
  parameters["npx_resampled_per_cpu"] = params->npx_resampled_per_cpu;
  parameters["npy_resampled_per_cpu"] = params->npy_resampled_per_cpu;
  parameters["npz_resampled_per_cpu"] = params->npz_resampled_per_cpu;
  parameters["tnow"] = params->tnow;
  parameters["xmin"] = params->xmin;
  parameters["xmax"] = params->xmax;
  parameters["ymin"] = params->ymin;
  parameters["ymax"] = params->ymax;
  parameters["zmin"] = params->zmin;
  parameters["zmax"] = params->zmax;

  jsoncons::json dump_analysis;

  dump_analysis["rms_emittance_x"] = emittance_x;
  dump_analysis["rms_emittance_y"] = emittance_y;
  dump_analysis["rms_emittance_z"] = emittance_z;
  dump_analysis["XMIN"] = estremi_min[0];
  dump_analysis["XMAX"] = estremi_max[0];
  dump_analysis["YMIN"] = estremi_min[1];
  dump_analysis["YMAX"] = estremi_max[1];
  dump_analysis["ZMIN"] = estremi_min[2];
  dump_analysis["ZMAX"] = estremi_max[2];
  dump_analysis["PXMIN"] = estremi_min[3];
  dump_analysis["PXMAX"] = estremi_max[3];
  dump_analysis["PYMIN"] = estremi_min[4];
  dump_analysis["PYMAX"] = estremi_max[4];
  dump_analysis["PZMIN"] = estremi_min[5];
  dump_analysis["PZMAX"] = estremi_max[5];
  dump_analysis["GAMMAMIN"] = estremi_min[6];
  dump_analysis["GAMMAMAX"] = estremi_max[6];
  dump_analysis["THETAMIN"] = estremi_min[7];
  dump_analysis["THETAMAX"] = estremi_max[7];
  dump_analysis["THETARADMIN"] = estremi_min[9];
  dump_analysis["THETARADMAX"] = estremi_max[9];
  dump_analysis["EMIN"] = estremi_min[8];
  dump_analysis["EMAX"] = estremi_max[8];
  dump_analysis["TYMIN"] = estremi_min[10];
  dump_analysis["TYMAX"] = estremi_max[10];
  dump_analysis["TZMIN"] = estremi_min[11];
  dump_analysis["TZMAX"] = estremi_max[11];
  dump_analysis["WMIN"] = estremi_min[12];
  dump_analysis["WMAX"] = estremi_max[12];
  dump_analysis["CHMIN"] = estremi_min[12];
  dump_analysis["CHMAX"] = estremi_max[12];
  dump_analysis["collected_weight"] = peso_accumulato;
  dump_analysis["collected_charge"] = carica_accumulata;

  jsoncons::json report;
  report["parameters"] = parameters;
  report["dump_analysis"] = dump_analysis;

  json_file << jsoncons::pretty_print(report) << std::endl;

  if (!params->multifile) fclose(file_in);

  std::cout << std::endl << "fread_size = " << fread_size << std::endl << "json parameters file produced!" << std::endl;

  delete[] estremi_min;
  estremi_min = nullptr;
  delete[] estremi_max;
  estremi_max = nullptr;

  return 0;
}

