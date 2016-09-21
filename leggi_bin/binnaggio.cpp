
#include "binnaggio.h"

_Binnaggio::_Binnaggio(float * particelle, int npart, int ndv, Parametri * parametri, float ** data_binned, std::string binx, std::string biny)
{
  int binnare_su_x = 0, binnare_su_y = 0;
  int whichbin_x = 0, whichbin_y = 0;

  if (binx == "x") binnare_su_x = 0;
  else if (binx == "y") binnare_su_x = 1;
  else if (binx == "z") binnare_su_x = 2;
  else if (binx == "px") binnare_su_x = 3;
  else if (binx == "py") binnare_su_x = 4;
  else if (binx == "pz") binnare_su_x = 5;
  else if (binx == "gamma") binnare_su_x = 6;
  else if (binx == "theta") binnare_su_x = 7;
  else if (binx == "E") binnare_su_x = 8;
  else if (binx == "thetaT") binnare_su_x = 9;
  else if (binx == "ty") binnare_su_x = 10;
  else if (binx == "tz") binnare_su_x = 11;
  else if (binx == "w") binnare_su_x = 12;
  else if (binx == "ch") binnare_su_x = 13;
  else printf("variabile x non riconosciuta\n");

  if (biny == "x") binnare_su_y = 0;
  else if (biny == "y") binnare_su_y = 1;
  else if (biny == "z") binnare_su_y = 2;
  else if (biny == "px") binnare_su_y = 3;
  else if (biny == "py") binnare_su_y = 4;
  else if (biny == "pz") binnare_su_y = 5;
  else if (biny == "gamma") binnare_su_y = 6;
  else if (biny == "theta") binnare_su_y = 7;
  else if (biny == "E") binnare_su_y = 8;
  else if (biny == "thetaT") binnare_su_y = 9;
  else if (biny == "ty") binnare_su_y = 10;
  else if (biny == "tz") binnare_su_y = 11;
  else if (biny == "w") binnare_su_y = 12;
  else if (biny == "ch") binnare_su_y = 13;
  else printf("variabile y non riconosciuta\n");


  float x, y, z, px, py, pz, w, ch, gamma, theta, thetaT, E, ty, tz;
  float dato_da_binnare_x = 0., dato_da_binnare_y = 0.;
  for (int i = 0; i < npart; i++)
  {
    if (((ndv == 4 || ndv == 5) && parametri->file_version < 3) || (ndv == 6 && parametri->file_version >= 3))
    {
      x = *(particelle + i*ndv);
      y = *(particelle + i*ndv + 1);
      px = *(particelle + i*ndv + 2);
      py = *(particelle + i*ndv + 3);
      if (parametri->p[WEIGHT] && !parametri->overwrite_weight)
        w = *(particelle + i*ndv + 4);
      else
        w = parametri->overwrite_weight_value;
      if (parametri->file_version >= 3 && !parametri->overwrite_charge)
        ch = *(particelle + i*ndv + 5);
      else
        ch = parametri->overwrite_charge_value;
      gamma = (float)(sqrt(1. + px*px + py*py) - 1.);      //gamma-1
      theta = (float)(atan2(py, px)*180. / M_PI);          //theta sgatto
      thetaT = (float)atan(sqrt((py*py / (px*px))));       //theta turch
      E = (float)(gamma*parametri->massa_particella_MeV);  //energia
      if (px > 0.0) ty = py / px;
      else ty = FLT_MAX;
      tz = 0.0;
      if (binnare_su_x == 0) dato_da_binnare_x = x;
      else if (binnare_su_x == 1) dato_da_binnare_x = y;
      else if (binnare_su_x == 2) std::cout << "Unable to bin on z in 2D" << std::endl;
      else if (binnare_su_x == 3) dato_da_binnare_x = px;
      else if (binnare_su_x == 4) dato_da_binnare_x = py;
      else if (binnare_su_x == 5) std::cout << "Unable to bin on pz in 2D" << std::endl;
      else if (binnare_su_x == 6) dato_da_binnare_x = gamma;
      else if (binnare_su_x == 7) dato_da_binnare_x = theta;
      else if (binnare_su_x == 8) dato_da_binnare_x = E;
      else if (binnare_su_x == 9) dato_da_binnare_x = thetaT;
      else if (binnare_su_x == 10) dato_da_binnare_x = ty;
      else if (binnare_su_x == 11) std::cout << "Unable to bin on tz in 2D" << std::endl;
      else if (binnare_su_x == 12) dato_da_binnare_x = w;
      else if (binnare_su_x == 13) dato_da_binnare_x = ch;
      if (binnare_su_y == 0) dato_da_binnare_y = x;
      else if (binnare_su_y == 1) dato_da_binnare_y = y;
      else if (binnare_su_y == 2) std::cout << "Unable to bin on z in 2D" << std::endl;
      else if (binnare_su_y == 3) dato_da_binnare_y = px;
      else if (binnare_su_y == 4) dato_da_binnare_y = py;
      else if (binnare_su_y == 5) std::cout << "Unable to bin on pz in 2D" << std::endl;
      else if (binnare_su_y == 6) dato_da_binnare_y = gamma;
      else if (binnare_su_y == 7) dato_da_binnare_y = theta;
      else if (binnare_su_y == 8) dato_da_binnare_y = E;
      else if (binnare_su_y == 9) dato_da_binnare_y = thetaT;
      else if (binnare_su_y == 10) dato_da_binnare_y = ty;
      else if (binnare_su_y == 11) std::cout << "Unable to bin on tz in 2D" << std::endl;
      else if (binnare_su_y == 12) dato_da_binnare_y = w;
      else if (binnare_su_y == 13) dato_da_binnare_y = ch;
    }
    else if (((ndv == 6 || ndv == 7) && parametri->file_version < 3) || (ndv == 8 && parametri->file_version >= 3))
    {
      x = *(particelle + i*ndv);
      y = *(particelle + i*ndv + 1);
      z = *(particelle + i*ndv + 2);
      px = *(particelle + i*ndv + 3);
      py = *(particelle + i*ndv + 4);
      pz = *(particelle + i*ndv + 5);
      if (parametri->p[WEIGHT] && !parametri->overwrite_weight)
        w = *(particelle + i*ndv + 6);
      else
        w = parametri->overwrite_weight_value;
      if (parametri->file_version >= 3 && !parametri->overwrite_charge)
        ch = *(particelle + i*ndv + 7);
      else
        ch = parametri->overwrite_charge_value;
      gamma = (float)(sqrt(1. + px*px + py*py + pz*pz) - 1.);             //gamma-1
      theta = (float)(atan2(sqrt(py*py + pz*pz), px)*180. / M_PI);        //theta sgatto
      thetaT = (float)atan(sqrt((py*py / (px*px)) + (pz*pz / (px*px))));  //theta turch
      E = (float)(gamma*parametri->massa_particella_MeV);                 //energia
      if (px > 0.0)
      {
        ty = py / px;
        tz = pz / px;
      }
      else
      {
        ty = FLT_MAX;
        tz = FLT_MAX;
      }
      if (binnare_su_x < 6) dato_da_binnare_x = *(particelle + i*ndv + binnare_su_x);
      else if (binnare_su_x == 6) dato_da_binnare_x = gamma;
      else if (binnare_su_x == 7) dato_da_binnare_x = theta;
      else if (binnare_su_x == 8) dato_da_binnare_x = E;
      else if (binnare_su_x == 9) dato_da_binnare_x = thetaT;
      else if (binnare_su_x == 10) dato_da_binnare_x = ty;
      else if (binnare_su_x == 11) dato_da_binnare_x = tz;
      else if (binnare_su_x == 12) dato_da_binnare_x = w;
      else if (binnare_su_x == 13) dato_da_binnare_x = ch;
      if (binnare_su_y < 6) dato_da_binnare_y = *(particelle + i*ndv + binnare_su_y);
      else if (binnare_su_y == 6) dato_da_binnare_y = gamma;
      else if (binnare_su_y == 7) dato_da_binnare_y = theta;
      else if (binnare_su_y == 8) dato_da_binnare_y = E;
      else if (binnare_su_y == 9) dato_da_binnare_y = thetaT;
      else if (binnare_su_y == 10) dato_da_binnare_y = ty;
      else if (binnare_su_y == 11) dato_da_binnare_y = tz;
      else if (binnare_su_y == 12) dato_da_binnare_y = w;
      else if (binnare_su_y == 13) dato_da_binnare_y = ch;
    }


    if (dato_da_binnare_x < parametri->minimi[binnare_su_x])
    {
      whichbin_x = 0;
    }
    else if (dato_da_binnare_x > parametri->massimi[binnare_su_x])
    {
      whichbin_x = parametri->dimmi_nbin(binnare_su_x) + 2;
    }
    else
    {
      whichbin_x = (int)(((dato_da_binnare_x - parametri->minimi[binnare_su_x]) / parametri->dimmi_dim(binnare_su_x)) + 1.0);
    }
    if (dato_da_binnare_y < parametri->minimi[binnare_su_y])
    {
      whichbin_y = 0;
    }
    else if (dato_da_binnare_y > parametri->massimi[binnare_su_y])
    {
      whichbin_y = parametri->dimmi_nbin(binnare_su_y) + 2;
    }
    else
    {
      whichbin_y = (int)(((dato_da_binnare_y - parametri->minimi[binnare_su_y]) / parametri->dimmi_dim(binnare_su_y)) + 1.0);
    }

    data_binned[whichbin_x][whichbin_y] += w;
  }
}

_Binnaggio::_Binnaggio(float * particelle, int npart, int ndv, Parametri * parametri, float * data_binned, std::string binx)
{
  int binnare_su_x = 0;
  int whichbin_x = 0;

  if (binx == "x") binnare_su_x = 0;
  else if (binx == "y") binnare_su_x = 1;
  else if (binx == "z") binnare_su_x = 2;
  else if (binx == "px") binnare_su_x = 3;
  else if (binx == "py") binnare_su_x = 4;
  else if (binx == "pz") binnare_su_x = 5;
  else if (binx == "gamma") binnare_su_x = 6;
  else if (binx == "theta") binnare_su_x = 7;
  else if (binx == "E") binnare_su_x = 8;
  else if (binx == "thetaT") binnare_su_x = 9;
  else if (binx == "ty") binnare_su_x = 10;
  else if (binx == "tz") binnare_su_x = 11;
  else if (binx == "w") binnare_su_x = 12;
  else if (binx == "ch") binnare_su_x = 13;
  else printf("variabile x non riconosciuta\n");

  float x, y, px, py, pz, w, ch, gamma, theta, thetaT, E, ty, tz;
  float dato_da_binnare_x = 0.;
  fflush(stdout);
  for (int i = 0; i < npart; i++)
  {
    if (((ndv == 4 || ndv == 5) && parametri->file_version < 3) || (ndv == 6 && parametri->file_version >= 3))
    {
      x = *(particelle + i*ndv);
      y = *(particelle + i*ndv + 1);
      px = *(particelle + i*ndv + 2);
      py = *(particelle + i*ndv + 3);
      if (parametri->p[WEIGHT] && !parametri->overwrite_weight)
        w = *(particelle + i*ndv + 4);
      else
        w = parametri->overwrite_weight_value;
      if (parametri->file_version >= 3 && !parametri->overwrite_charge)
        ch = *(particelle + i*ndv + 5);
      else
        ch = parametri->overwrite_charge_value;
      gamma = (float)(sqrt(1. + px*px + py*py) - 1.);      //gamma-1
      theta = (float)(atan2(py, px)*180. / M_PI);          //theta sgatto
      thetaT = (float)atan(sqrt((py*py / (px*px))));       //theta turch
      E = (float)(gamma*parametri->massa_particella_MeV);  //energia
      if (px > 0.0) ty = py / px;
      else ty = 0.0;
      tz = 0.0;
      if (binnare_su_x == 0) dato_da_binnare_x = x;
      else if (binnare_su_x == 1) dato_da_binnare_x = y;
      else if (binnare_su_x == 2) std::cout << "Unable to bin on z in 2D" << std::endl;
      else if (binnare_su_x == 3) dato_da_binnare_x = px;
      else if (binnare_su_x == 4) dato_da_binnare_x = py;
      else if (binnare_su_x == 5) std::cout << "Unable to bin on pz in 2D" << std::endl;
      else if (binnare_su_x == 6) dato_da_binnare_x = gamma;
      else if (binnare_su_x == 7) dato_da_binnare_x = theta;
      else if (binnare_su_x == 8) dato_da_binnare_x = E;
      else if (binnare_su_x == 9) dato_da_binnare_x = thetaT;
      else if (binnare_su_x == 10) dato_da_binnare_x = ty;
      else if (binnare_su_x == 11) std::cout << "Unable to bin on tz in 2D" << std::endl;
      else if (binnare_su_x == 12) dato_da_binnare_x = w;
    }
    else if (((ndv == 6 || ndv == 7) && parametri->file_version < 3) || (ndv == 8 && parametri->file_version >= 3))
    {
      px = *(particelle + i*ndv + 3);
      py = *(particelle + i*ndv + 4);
      pz = *(particelle + i*ndv + 5);
      if (parametri->p[WEIGHT] && !parametri->overwrite_weight)
        w = *(particelle + i*ndv + 6);
      else
        w = parametri->overwrite_weight_value;
      if (parametri->file_version >= 3 && !parametri->overwrite_charge)
        ch = *(particelle + i*ndv + 7);
      else
        ch = parametri->overwrite_charge_value;
      gamma = (float)(sqrt(1. + px*px + py*py + pz*pz) - 1.);             //gamma
      theta = (float)(atan2(sqrt(py*py + pz*pz), px)*180. / M_PI);        //theta sgatto
      thetaT = (float)atan(sqrt((py*py / (px*px)) + (pz*pz / (px*px))));  //theta turch
      E = (float)(gamma*parametri->massa_particella_MeV);                 //energia
      if (px > 0.0)
      {
        ty = py / px;
        tz = pz / px;
      }
      else
      {
        ty = 0.0;
        tz = 0.0;
      }
      if (binnare_su_x < 6) dato_da_binnare_x = *(particelle + i*ndv + binnare_su_x);
      else if (binnare_su_x == 6) dato_da_binnare_x = gamma;
      else if (binnare_su_x == 7) dato_da_binnare_x = theta;
      else if (binnare_su_x == 8) dato_da_binnare_x = E;
      else if (binnare_su_x == 9) dato_da_binnare_x = thetaT;
      else if (binnare_su_x == 10) dato_da_binnare_x = ty;
      else if (binnare_su_x == 11) dato_da_binnare_x = tz;
      else if (binnare_su_x == 12) dato_da_binnare_x = w;
      else if (binnare_su_x == 13) dato_da_binnare_x = ch;
    }

    if (dato_da_binnare_x < parametri->minimi[binnare_su_x])
    {
      whichbin_x = 0;
    }
    else if (dato_da_binnare_x > parametri->massimi[binnare_su_x])
    {
      whichbin_x = parametri->dimmi_nbin(binnare_su_x) + 2;
    }
    else
    {
      whichbin_x = (int)(((dato_da_binnare_x - parametri->minimi[binnare_su_x]) / parametri->dimmi_dim(binnare_su_x)) + 1.0);
    }

    data_binned[whichbin_x] += w;
  }
}

