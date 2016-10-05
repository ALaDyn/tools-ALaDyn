
#include "scrittura.h"


_Write::_Write(Parameters * params, aladyn_float ** data_binned, std::string x, std::string y, std::string filename_out)
{
  aladyn_float xmin, xmax, dimx, ymin, ymax, dimy;
  int nbinx, nbiny;
  std::ofstream file_out;
  file_out.open(filename_out.c_str(), std::ofstream::out);

  if (x == "x")
  {
    dimx = (params->get_size_x());
    xmin = (params->xmin) - dimx;
    xmax = (params->xmin);
    nbinx = (params->nbin_x + 3);
  }
  else if (x == "y")
  {
    dimx = (params->get_size_y());
    xmin = (params->ymin) - dimx;
    xmax = (params->ymin);
    nbinx = (params->nbin_y + 3);
  }
  else if (x == "z")
  {
    dimx = (params->get_size_z());
    xmin = (params->zmin) - dimx;
    xmax = (params->zmin);
    nbinx = (params->nbin_z + 3);
  }
  else if (x == "px")
  {
    dimx = (params->get_size_px());
    xmin = (params->pxmin) - dimx;
    xmax = (params->pxmin);
    nbinx = (params->nbin_px + 3);
  }
  else if (x == "py")
  {
    dimx = (params->get_size_py());
    xmin = (params->pymin) - dimx;
    xmax = (params->pymin);
    nbinx = (params->nbin_py + 3);
  }
  else if (x == "pz")
  {
    dimx = (params->get_size_pz());
    xmin = (params->pzmin) - dimx;
    xmax = (params->pzmin);
    nbinx = (params->nbin_pz + 3);
  }
  else if (x == "gamma")
  {
    dimx = (params->get_size_gamma());
    xmin = (params->gammamin) - dimx;
    xmax = (params->gammamin);
    nbinx = (params->nbin_gamma + 3);
  }
  else if (x == "theta")
  {
    dimx = (params->get_size_theta());
    xmin = (params->thetamin) - dimx;
    xmax = (params->thetamin);
    nbinx = (params->nbin_theta + 3);
  }
  else if (x == "E")
  {
    dimx = (params->get_size_E());
    xmin = (params->Emin) - dimx;
    xmax = (params->Emin);
    nbinx = (params->nbin_E + 3);
  }
  else if (x == "thetaT")
  {
    dimx = (params->get_size_thetaT());
    xmin = (params->thetaTmin) - dimx;
    xmax = (params->thetaTmin);
    nbinx = (params->nbin_thetaT + 3);
  }
  else if (x == "ty")
  {
    dimx = (params->get_size_ty());
    xmin = (params->tymin) - dimx;
    xmax = (params->tymin);
    nbinx = (params->nbin_ty + 3);
  }
  else if (x == "tz")
  {
    dimx = (params->get_size_tz());
    xmin = (params->tzmin) - dimx;
    xmax = (params->tzmin);
    nbinx = (params->nbin_tz + 3);
  }
  else if (x == "w")
  {
    dimx = (params->get_size_w());
    xmin = (params->wmin) - dimx;
    xmax = (params->wmin);
    nbinx = (params->nbin_w + 3);
  }
  else if (x == "ch")
  {
    dimx = (params->get_size_ch());
    xmin = (params->chmin) - dimx;
    xmax = (params->chmin);
    nbinx = (params->nbin_ch + 3);
  }
  else std::cerr << "Unrecognized x variable" << std::endl;


  if (y == "x")
  {
    dimy = (params->get_size_x());
    ymin = (params->xmin) - dimy;
    ymax = (params->xmin);
    nbiny = (params->nbin_x + 3);
  }
  else if (y == "y")
  {
    dimy = (params->get_size_y());
    ymin = (params->ymin) - dimy;
    ymax = (params->ymin);
    nbiny = (params->nbin_y + 3);
  }
  else if (y == "z")
  {
    dimy = (params->get_size_z());
    ymin = (params->zmin) - dimy;
    ymax = (params->zmin);
    nbiny = (params->nbin_z + 3);
  }
  else if (y == "px")
  {
    dimy = (params->get_size_px());
    ymin = (params->pxmin) - dimy;
    ymax = (params->pxmin);
    nbiny = (params->nbin_px + 3);
  }
  else if (y == "py")
  {
    dimy = (params->get_size_py());
    ymin = (params->pymin) - dimy;
    ymax = (params->pymin);
    nbiny = (params->nbin_py + 3);
  }
  else if (y == "pz")
  {
    dimy = (params->get_size_pz());
    ymin = (params->pzmin) - dimy;
    ymax = (params->pzmin);
    nbiny = (params->nbin_pz + 3);
  }
  else if (y == "gamma")
  {
    dimy = (params->get_size_gamma());
    ymin = (params->gammamin) - dimy;
    ymax = (params->gammamin);
    nbiny = (params->nbin_gamma + 3);
  }
  else if (y == "theta")
  {
    dimy = (params->get_size_theta());
    ymin = (params->thetamin) - dimy;
    ymax = (params->thetamin);
    nbiny = (params->nbin_theta + 3);
  }
  else if (y == "E")
  {
    dimy = (params->get_size_E());
    ymin = (params->Emin) - dimy;
    ymax = (params->Emin);
    nbiny = (params->nbin_E + 3);
  }
  else if (y == "thetaT")
  {
    dimy = (params->get_size_thetaT());
    ymin = (params->thetaTmin) - dimy;
    ymax = (params->thetaTmin);
    nbiny = (params->nbin_thetaT + 3);
  }
  else if (y == "ty")
  {
    dimy = (params->get_size_ty());
    ymin = (params->tymin) - dimy;
    ymax = (params->tymin);
    nbiny = (params->nbin_ty + 3);
  }
  else if (y == "tz")
  {
    dimy = (params->get_size_tz());
    ymin = (params->tzmin) - dimy;
    ymax = (params->tzmin);
    nbiny = (params->nbin_tz + 3);
  }
  else if (y == "w")
  {
    dimy = (params->get_size_w());
    ymin = (params->wmin) - dimy;
    ymax = (params->wmin);
    nbiny = (params->nbin_w + 3);
  }
  else if (y == "ch")
  {
    dimy = (params->get_size_ch());
    ymin = (params->chmin) - dimy;
    ymax = (params->chmin);
    nbiny = (params->nbin_ch + 3);
  }
  else std::cerr << "Unrecognized y variable" << std::endl;


  aladyn_float yminT = ymin;
  aladyn_float ymaxT = ymax;

  for (int i = 0; i < nbinx; i++)
  {
    for (int j = 0; j < nbiny; j++)
    {
      file_out << std::setprecision(6) << xmin << "\t" << xmax << "\t" << yminT << "\t" << ymaxT << "\t" << data_binned[i][j] << std::endl;
      yminT += dimy;
      ymaxT += dimy;
    }
    xmin += dimx;
    xmax += dimx;
    yminT = ymin;
    ymaxT = ymax;
  }
  file_out.close();
}



_Write::_Write(Parameters * params, aladyn_float * data_binned, std::string x, std::string filename_out)
{
  aladyn_float xmin, xmax, dimx;
  int nbinx;
  std::ofstream file_out;
  file_out.open(filename_out.c_str(), std::ofstream::out);

  if (x == "x")
  {
    dimx = (params->get_size_x());
    xmin = (params->xmin) - dimx;
    xmax = (params->xmin);
    nbinx = (params->nbin_x + 3);
  }
  else if (x == "y")
  {
    dimx = (params->get_size_y());
    xmin = (params->ymin) - dimx;
    xmax = (params->ymin);
    nbinx = (params->nbin_y + 3);
  }
  else if (x == "z")
  {
    dimx = (params->get_size_z());
    xmin = (params->zmin) - dimx;
    xmax = (params->zmin);
    nbinx = (params->nbin_z + 3);
  }
  else if (x == "px")
  {
    dimx = (params->get_size_px());
    xmin = (params->pxmin) - dimx;
    xmax = (params->pxmin);
    nbinx = (params->nbin_px + 3);
  }
  else if (x == "py")
  {
    dimx = (params->get_size_py());
    xmin = (params->pymin) - dimx;
    xmax = (params->pymin);
    nbinx = (params->nbin_py + 3);
  }
  else if (x == "pz")
  {
    dimx = (params->get_size_pz());
    xmin = (params->pzmin) - dimx;
    xmax = (params->pzmin);
    nbinx = (params->nbin_pz + 3);
  }
  else if (x == "gamma")
  {
    dimx = (params->get_size_gamma());
    xmin = (params->gammamin) - dimx;
    xmax = (params->gammamin);
    nbinx = (params->nbin_gamma + 3);
  }
  else if (x == "theta")
  {
    dimx = (params->get_size_theta());
    xmin = (params->thetamin) - dimx;
    xmax = (params->thetamin);
    nbinx = (params->nbin_theta + 3);
  }
  else if (x == "E")
  {
    dimx = (params->get_size_E());
    xmin = (params->Emin) - dimx;
    xmax = (params->Emin);
    nbinx = (params->nbin_E + 3);
  }
  else if (x == "thetaT")
  {
    dimx = (params->get_size_thetaT());
    xmin = (params->thetaTmin) - dimx;
    xmax = (params->thetaTmin);
    nbinx = (params->nbin_thetaT + 3);
  }
  else if (x == "ty")
  {
    dimx = (params->get_size_ty());
    xmin = (params->tymin) - dimx;
    xmax = (params->tymin);
    nbinx = (params->nbin_ty + 3);
  }
  else if (x == "tz")
  {
    dimx = (params->get_size_tz());
    xmin = (params->tzmin) - dimx;
    xmax = (params->tzmin);
    nbinx = (params->nbin_tz + 3);
  }
  else if (x == "w")
  {
    dimx = (params->get_size_w());
    xmin = (params->wmin) - dimx;
    xmax = (params->wmin);
    nbinx = (params->nbin_w + 3);
  }
  else if (x == "ch")
  {
    dimx = (params->get_size_ch());
    xmin = (params->chmin) - dimx;
    xmax = (params->chmin);
    nbinx = (params->nbin_ch + 3);
  }
  else std::cerr << "Unrecognized x variable" << std::endl;


  for (int i = 0; i < nbinx; i++)
  {
    file_out << std::setprecision(6) << xmin << "\t" << xmax << "\t" << data_binned[i] << std::endl;
    xmin += dimx;
    xmax += dimx;
  }
  file_out.close();
}

