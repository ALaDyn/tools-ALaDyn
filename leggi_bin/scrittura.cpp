
#include "scrittura.h"


_Scrittura::_Scrittura(Parametri * parametri, float ** data_binned, std::string x, std::string y, std::string nomefile_out)
{
  float xmin, xmax, dimx, ymin, ymax, dimy;
  int nbinx, nbiny;
  std::ofstream file_out;
  file_out.open(nomefile_out.c_str(), std::ofstream::out);

  if (x == "x")
  {
    dimx = (parametri->dimmi_dimx());
    xmin = (parametri->xmin) - dimx;
    xmax = (parametri->xmin);
    nbinx = (parametri->nbin_x + 3);
  }
  else if (x == "y")
  {
    dimx = (parametri->dimmi_dimy());
    xmin = (parametri->ymin) - dimx;
    xmax = (parametri->ymin);
    nbinx = (parametri->nbin_y + 3);
  }
  else if (x == "z")
  {
    dimx = (parametri->dimmi_dimz());
    xmin = (parametri->zmin) - dimx;
    xmax = (parametri->zmin);
    nbinx = (parametri->nbin_z + 3);
  }
  else if (x == "px")
  {
    dimx = (parametri->dimmi_dimpx());
    xmin = (parametri->pxmin) - dimx;
    xmax = (parametri->pxmin);
    nbinx = (parametri->nbin_px + 3);
  }
  else if (x == "py")
  {
    dimx = (parametri->dimmi_dimpy());
    xmin = (parametri->pymin) - dimx;
    xmax = (parametri->pymin);
    nbinx = (parametri->nbin_py + 3);
  }
  else if (x == "pz")
  {
    dimx = (parametri->dimmi_dimpz());
    xmin = (parametri->pzmin) - dimx;
    xmax = (parametri->pzmin);
    nbinx = (parametri->nbin_pz + 3);
  }
  else if (x == "gamma")
  {
    dimx = (parametri->dimmi_dimgamma());
    xmin = (parametri->gammamin) - dimx;
    xmax = (parametri->gammamin);
    nbinx = (parametri->nbin_gamma + 3);
  }
  else if (x == "theta")
  {
    dimx = (parametri->dimmi_dimtheta());
    xmin = (parametri->thetamin) - dimx;
    xmax = (parametri->thetamin);
    nbinx = (parametri->nbin_theta + 3);
  }
  else if (x == "E")
  {
    dimx = (parametri->dimmi_dimE());
    xmin = (parametri->Emin) - dimx;
    xmax = (parametri->Emin);
    nbinx = (parametri->nbin_E + 3);
  }
  else if (x == "thetaT")
  {
    dimx = (parametri->dimmi_dimthetaT());
    xmin = (parametri->thetaTmin) - dimx;
    xmax = (parametri->thetaTmin);
    nbinx = (parametri->nbin_thetaT + 3);
  }
  else if (x == "ty")
  {
    dimx = (parametri->dimmi_dimty());
    xmin = (parametri->tymin) - dimx;
    xmax = (parametri->tymin);
    nbinx = (parametri->nbin_ty + 3);
  }
  else if (x == "tz")
  {
    dimx = (parametri->dimmi_dimtz());
    xmin = (parametri->tzmin) - dimx;
    xmax = (parametri->tzmin);
    nbinx = (parametri->nbin_tz + 3);
  }
  else if (x == "w")
  {
    dimx = (parametri->dimmi_dimw());
    xmin = (parametri->wmin) - dimx;
    xmax = (parametri->wmin);
    nbinx = (parametri->nbin_w + 3);
  }
  else if (x == "ch")
  {
    dimx = (parametri->dimmi_dimch());
    xmin = (parametri->chmin) - dimx;
    xmax = (parametri->chmin);
    nbinx = (parametri->nbin_ch + 3);
  }
  else printf("variabile x non riconosciuta\n");


  if (y == "x")
  {
    dimy = (parametri->dimmi_dimx());
    ymin = (parametri->xmin) - dimy;
    ymax = (parametri->xmin);
    nbiny = (parametri->nbin_x + 3);
  }
  else if (y == "y")
  {
    dimy = (parametri->dimmi_dimy());
    ymin = (parametri->ymin) - dimy;
    ymax = (parametri->ymin);
    nbiny = (parametri->nbin_y + 3);
  }
  else if (y == "z")
  {
    dimy = (parametri->dimmi_dimz());
    ymin = (parametri->zmin) - dimy;
    ymax = (parametri->zmin);
    nbiny = (parametri->nbin_z + 3);
  }
  else if (y == "px")
  {
    dimy = (parametri->dimmi_dimpx());
    ymin = (parametri->pxmin) - dimy;
    ymax = (parametri->pxmin);
    nbiny = (parametri->nbin_px + 3);
  }
  else if (y == "py")
  {
    dimy = (parametri->dimmi_dimpy());
    ymin = (parametri->pymin) - dimy;
    ymax = (parametri->pymin);
    nbiny = (parametri->nbin_py + 3);
  }
  else if (y == "pz")
  {
    dimy = (parametri->dimmi_dimpz());
    ymin = (parametri->pzmin) - dimy;
    ymax = (parametri->pzmin);
    nbiny = (parametri->nbin_pz + 3);
  }
  else if (y == "gamma")
  {
    dimy = (parametri->dimmi_dimgamma());
    ymin = (parametri->gammamin) - dimy;
    ymax = (parametri->gammamin);
    nbiny = (parametri->nbin_gamma + 3);
  }
  else if (y == "theta")
  {
    dimy = (parametri->dimmi_dimtheta());
    ymin = (parametri->thetamin) - dimy;
    ymax = (parametri->thetamin);
    nbiny = (parametri->nbin_theta + 3);
  }
  else if (y == "E")
  {
    dimy = (parametri->dimmi_dimE());
    ymin = (parametri->Emin) - dimy;
    ymax = (parametri->Emin);
    nbiny = (parametri->nbin_E + 3);
  }
  else if (y == "thetaT")
  {
    dimy = (parametri->dimmi_dimthetaT());
    ymin = (parametri->thetaTmin) - dimy;
    ymax = (parametri->thetaTmin);
    nbiny = (parametri->nbin_thetaT + 3);
  }
  else if (y == "ty")
  {
    dimy = (parametri->dimmi_dimty());
    ymin = (parametri->tymin) - dimy;
    ymax = (parametri->tymin);
    nbiny = (parametri->nbin_ty + 3);
  }
  else if (y == "tz")
  {
    dimy = (parametri->dimmi_dimtz());
    ymin = (parametri->tzmin) - dimy;
    ymax = (parametri->tzmin);
    nbiny = (parametri->nbin_tz + 3);
  }
  else if (y == "w")
  {
    dimy = (parametri->dimmi_dimw());
    ymin = (parametri->wmin) - dimy;
    ymax = (parametri->wmin);
    nbiny = (parametri->nbin_w + 3);
  }
  else if (y == "ch")
  {
    dimy = (parametri->dimmi_dimch());
    ymin = (parametri->chmin) - dimy;
    ymax = (parametri->chmin);
    nbiny = (parametri->nbin_ch + 3);
  }
  else printf("variabile y non riconosciuta\n");


  float yminT = ymin;
  float ymaxT = ymax;

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



_Scrittura::_Scrittura(Parametri * parametri, float * data_binned, std::string x, std::string nomefile_out)
{
  float xmin, xmax, dimx;
  int nbinx;
  std::ofstream file_out;
  file_out.open(nomefile_out.c_str(), std::ofstream::out);

  if (x == "x")
  {
    dimx = (parametri->dimmi_dimx());
    xmin = (parametri->xmin) - dimx;
    xmax = (parametri->xmin);
    nbinx = (parametri->nbin_x + 3);
  }
  else if (x == "y")
  {
    dimx = (parametri->dimmi_dimy());
    xmin = (parametri->ymin) - dimx;
    xmax = (parametri->ymin);
    nbinx = (parametri->nbin_y + 3);
  }
  else if (x == "z")
  {
    dimx = (parametri->dimmi_dimz());
    xmin = (parametri->zmin) - dimx;
    xmax = (parametri->zmin);
    nbinx = (parametri->nbin_z + 3);
  }
  else if (x == "px")
  {
    dimx = (parametri->dimmi_dimpx());
    xmin = (parametri->pxmin) - dimx;
    xmax = (parametri->pxmin);
    nbinx = (parametri->nbin_px + 3);
  }
  else if (x == "py")
  {
    dimx = (parametri->dimmi_dimpy());
    xmin = (parametri->pymin) - dimx;
    xmax = (parametri->pymin);
    nbinx = (parametri->nbin_py + 3);
  }
  else if (x == "pz")
  {
    dimx = (parametri->dimmi_dimpz());
    xmin = (parametri->pzmin) - dimx;
    xmax = (parametri->pzmin);
    nbinx = (parametri->nbin_pz + 3);
  }
  else if (x == "gamma")
  {
    dimx = (parametri->dimmi_dimgamma());
    xmin = (parametri->gammamin) - dimx;
    xmax = (parametri->gammamin);
    nbinx = (parametri->nbin_gamma + 3);
  }
  else if (x == "theta")
  {
    dimx = (parametri->dimmi_dimtheta());
    xmin = (parametri->thetamin) - dimx;
    xmax = (parametri->thetamin);
    nbinx = (parametri->nbin_theta + 3);
  }
  else if (x == "E")
  {
    dimx = (parametri->dimmi_dimE());
    xmin = (parametri->Emin) - dimx;
    xmax = (parametri->Emin);
    nbinx = (parametri->nbin_E + 3);
  }
  else if (x == "thetaT")
  {
    dimx = (parametri->dimmi_dimthetaT());
    xmin = (parametri->thetaTmin) - dimx;
    xmax = (parametri->thetaTmin);
    nbinx = (parametri->nbin_thetaT + 3);
  }
  else if (x == "ty")
  {
    dimx = (parametri->dimmi_dimty());
    xmin = (parametri->tymin) - dimx;
    xmax = (parametri->tymin);
    nbinx = (parametri->nbin_ty + 3);
  }
  else if (x == "tz")
  {
    dimx = (parametri->dimmi_dimtz());
    xmin = (parametri->tzmin) - dimx;
    xmax = (parametri->tzmin);
    nbinx = (parametri->nbin_tz + 3);
  }
  else if (x == "w")
  {
    dimx = (parametri->dimmi_dimw());
    xmin = (parametri->wmin) - dimx;
    xmax = (parametri->wmin);
    nbinx = (parametri->nbin_w + 3);
  }
  else if (x == "ch")
  {
    dimx = (parametri->dimmi_dimch());
    xmin = (parametri->chmin) - dimx;
    xmax = (parametri->chmin);
    nbinx = (parametri->nbin_ch + 3);
  }
  else printf("variabile x non riconosciuta\n");


  for (int i = 0; i < nbinx; i++)
  {
    file_out << std::setprecision(6) << xmin << "\t" << xmax << "\t" << data_binned[i] << std::endl;
    xmin += dimx;
    xmax += dimx;
  }
  file_out.close();
}

