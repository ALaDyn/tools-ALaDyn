
#include "binnaggio.h"

_Binning::_Binning(aladyn_float * parts, Parameters * params, aladyn_float ** data_binned, std::string binx, std::string biny)
{
  int bin_on_x = 0, bin_on_y = 0;
  int whichbin_x = 0, whichbin_y = 0;
  size_t ndv = params->ndv;

  double min_x = 0.0, min_y = 0.0, max_x = 1.0, max_y = 1.0;

  if (binx == "x") bin_on_x = COLUMN_X; 
  else if (binx == "y") bin_on_x = COLUMN_Y; 
  else if (binx == "z") bin_on_x = COLUMN_Z; 
  else if (binx == "px") bin_on_x = COLUMN_PX; 
  else if (binx == "py") bin_on_x = COLUMN_PY; 
  else if (binx == "pz") bin_on_x = COLUMN_PZ; 
  else if (binx == "gamma") bin_on_x = COLUMN_GAMMA; 
  else if (binx == "theta") bin_on_x = COLUMN_THETA; 
  else if (binx == "E") bin_on_x = COLUMN_E; 
  else if (binx == "thetaT") bin_on_x = COLUMN_THETAT; 
  else if (binx == "ty") bin_on_x = COLUMN_TY; 
  else if (binx == "tz") bin_on_x = COLUMN_TZ; 
  else if (binx == "w") bin_on_x = COLUMN_W; 
  else if (binx == "ch") bin_on_x = COLUMN_CH; 
  else std::cerr << "Unrecognized x variable" << std::endl;

  if (biny == "x") bin_on_y = COLUMN_X;
  else if (biny == "y") bin_on_y = COLUMN_Y;
  else if (biny == "z") bin_on_y = COLUMN_Z;
  else if (biny == "px") bin_on_y = COLUMN_PX;
  else if (biny == "py") bin_on_y = COLUMN_PY;
  else if (biny == "pz") bin_on_y = COLUMN_PZ;
  else if (biny == "gamma") bin_on_y = COLUMN_GAMMA;
  else if (biny == "theta") bin_on_y = COLUMN_THETA;
  else if (biny == "E") bin_on_y = COLUMN_E;
  else if (biny == "thetaT") bin_on_y = COLUMN_THETAT;
  else if (biny == "ty") bin_on_y = COLUMN_TY;
  else if (biny == "tz") bin_on_y = COLUMN_TZ;
  else if (biny == "w") bin_on_y = COLUMN_W;
  else if (biny == "ch") bin_on_y = COLUMN_CH;
  else std::cerr << "Unrecognized y variable" << std::endl;

  min_x = params->get_min(bin_on_x);
  max_x = params->get_max(bin_on_x);
  min_y = params->get_min(bin_on_y);
  max_y = params->get_max(bin_on_y);

  aladyn_float x, y, z, px, py, pz, w, ch, gamma, theta, thetaT, E, ty, tz;
  aladyn_float going_to_bin_on_x = 0., going_to_bin_on_y = 0.;
  for (int i = 0; i < params->nparts; i++)
  {
    if (params->sim_is_2d)
    {
      x = *(parts + i*ndv);
      y = *(parts + i*ndv + 1);
      px = *(parts + i*ndv + 2);
      py = *(parts + i*ndv + 3);
      if (params->file_has_weight && !params->overwrite_weight)
        w = *(parts + i*ndv + 4);
      else
        w = params->overwrite_weight_value;
      if (params->file_has_charge && !params->overwrite_charge)
        ch = *(parts + i*ndv + 5);
      else
        ch = params->overwrite_charge_value;
      gamma = (aladyn_float)(sqrt(1. + px*px + py*py) - 1.);      //gamma-1
      theta = (aladyn_float)(atan2(py, px)*180. / M_PI);          //theta sgatto
      thetaT = (aladyn_float)atan(sqrt((py*py / (px*px))));       //theta turch
      E = (aladyn_float)(gamma*params->mass_MeV);
      if (px > 0.0) ty = py / px;
      else ty = FLT_MAX;
      tz = 0.0;
      if (bin_on_x == COLUMN_X) going_to_bin_on_x = x;
      else if (bin_on_x == COLUMN_Y) going_to_bin_on_x = y;
      else if (bin_on_x == COLUMN_Z) std::cout << "Unable to bin on z in 2D" << std::endl;
      else if (bin_on_x == COLUMN_PX) going_to_bin_on_x = px;
      else if (bin_on_x == COLUMN_PY) going_to_bin_on_x = py;
      else if (bin_on_x == COLUMN_PZ) std::cout << "Unable to bin on pz in 2D" << std::endl;
      else if (bin_on_x == COLUMN_GAMMA) going_to_bin_on_x = gamma;
      else if (bin_on_x == COLUMN_THETA) going_to_bin_on_x = theta;
      else if (bin_on_x == COLUMN_E) going_to_bin_on_x = E;
      else if (bin_on_x == COLUMN_THETAT) going_to_bin_on_x = thetaT;
      else if (bin_on_x == COLUMN_TY) going_to_bin_on_x = ty;
      else if (bin_on_x == COLUMN_TZ) std::cout << "Unable to bin on tz in 2D" << std::endl;
      else if (bin_on_x == COLUMN_W) going_to_bin_on_x = w;
      else if (bin_on_x == COLUMN_CH) going_to_bin_on_x = ch;
      if (bin_on_y == COLUMN_X) going_to_bin_on_y = x;
      else if (bin_on_y == COLUMN_Y) going_to_bin_on_y = y;
      else if (bin_on_y == COLUMN_Z) std::cout << "Unable to bin on z in 2D" << std::endl;
      else if (bin_on_y == COLUMN_PX) going_to_bin_on_y = px;
      else if (bin_on_y == COLUMN_PY) going_to_bin_on_y = py;
      else if (bin_on_y == COLUMN_PZ) std::cout << "Unable to bin on pz in 2D" << std::endl;
      else if (bin_on_y == COLUMN_GAMMA) going_to_bin_on_y = gamma;
      else if (bin_on_y == COLUMN_THETA) going_to_bin_on_y = theta;
      else if (bin_on_y == COLUMN_E) going_to_bin_on_y = E;
      else if (bin_on_y == COLUMN_THETAT) going_to_bin_on_y = thetaT;
      else if (bin_on_y == COLUMN_TY) going_to_bin_on_y = ty;
      else if (bin_on_y == COLUMN_TZ) std::cout << "Unable to bin on tz in 2D" << std::endl;
      else if (bin_on_y == COLUMN_W) going_to_bin_on_y = w;
      else if (bin_on_y == COLUMN_CH) going_to_bin_on_y = ch;
    }
    else
    {
      x = *(parts + i*ndv);
      y = *(parts + i*ndv + 1);
      z = *(parts + i*ndv + 2);
      px = *(parts + i*ndv + 3);
      py = *(parts + i*ndv + 4);
      pz = *(parts + i*ndv + 5);
      if (params->file_has_weight && !params->overwrite_weight)
        w = *(parts + i*ndv + 6);
      else
        w = params->overwrite_weight_value;
      if (params->file_has_charge && !params->overwrite_charge)
        ch = *(parts + i*ndv + 7);
      else
        ch = params->overwrite_charge_value;
      gamma = (aladyn_float)(sqrt(1. + px*px + py*py + pz*pz) - 1.);             //gamma-1
      theta = (aladyn_float)(atan2(sqrt(py*py + pz*pz), px)*180. / M_PI);        //theta sgatto
      thetaT = (aladyn_float)atan(sqrt((py*py / (px*px)) + (pz*pz / (px*px))));  //theta turch
      E = (aladyn_float)(gamma*params->mass_MeV);
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
      if (bin_on_x < COLUMN_GAMMA) going_to_bin_on_x = *(parts + i*ndv + bin_on_x);
      else if (bin_on_x == COLUMN_GAMMA) going_to_bin_on_x = gamma;
      else if (bin_on_x == COLUMN_THETA) going_to_bin_on_x = theta;
      else if (bin_on_x == COLUMN_E) going_to_bin_on_x = E;
      else if (bin_on_x == COLUMN_THETAT) going_to_bin_on_x = thetaT;
      else if (bin_on_x == COLUMN_TY) going_to_bin_on_x = ty;
      else if (bin_on_x == COLUMN_TZ) going_to_bin_on_x = tz;
      else if (bin_on_x == COLUMN_W) going_to_bin_on_x = w;
      else if (bin_on_x == COLUMN_CH) going_to_bin_on_x = ch;
      if (bin_on_y < COLUMN_GAMMA) going_to_bin_on_y = *(parts + i*ndv + bin_on_y);
      else if (bin_on_y == COLUMN_GAMMA) going_to_bin_on_y = gamma;
      else if (bin_on_y == COLUMN_THETA) going_to_bin_on_y = theta;
      else if (bin_on_y == COLUMN_E) going_to_bin_on_y = E;
      else if (bin_on_y == COLUMN_THETAT) going_to_bin_on_y = thetaT;
      else if (bin_on_y == COLUMN_TY) going_to_bin_on_y = ty;
      else if (bin_on_y == COLUMN_TZ) going_to_bin_on_y = tz;
      else if (bin_on_y == COLUMN_W) going_to_bin_on_y = w;
      else if (bin_on_y == COLUMN_CH) going_to_bin_on_y = ch;
    }


    if (going_to_bin_on_x < min_x)
    {
      whichbin_x = 0;
    }
    else if (going_to_bin_on_x > max_x)
    {
      whichbin_x = params->get_nbin(bin_on_x) + 2;
    }
    else
    {
      whichbin_x = (int)(((going_to_bin_on_x - min_x) / params->get_size_(bin_on_x)) + 1.0);
    }
    if (going_to_bin_on_y < min_y)
    {
      whichbin_y = 0;
    }
    else if (going_to_bin_on_y > max_x)
    {
      whichbin_y = params->get_nbin(bin_on_y) + 2;
    }
    else
    {
      whichbin_y = (int)(((going_to_bin_on_y - min_y) / params->get_size_(bin_on_y)) + 1.0);
    }

    data_binned[whichbin_x][whichbin_y] += w;
  }
}

_Binning::_Binning(aladyn_float * parts, Parameters * params, aladyn_float * data_binned, std::string binx)
{
  int bin_on_x = 0;
  int whichbin_x = 0;
  size_t ndv = params->ndv;
  double min_x = 0.0, max_x = 1.0;

  if (binx == "x") bin_on_x = COLUMN_X;
  else if (binx == "y") bin_on_x = COLUMN_Y;
  else if (binx == "z") bin_on_x = COLUMN_Z;
  else if (binx == "px") bin_on_x = COLUMN_PX;
  else if (binx == "py") bin_on_x = COLUMN_PY;
  else if (binx == "pz") bin_on_x = COLUMN_PZ;
  else if (binx == "gamma") bin_on_x = COLUMN_GAMMA;
  else if (binx == "theta") bin_on_x = COLUMN_THETA;
  else if (binx == "E") bin_on_x = COLUMN_E;
  else if (binx == "thetaT") bin_on_x = COLUMN_THETAT;
  else if (binx == "ty") bin_on_x = COLUMN_TY;
  else if (binx == "tz") bin_on_x = COLUMN_TZ;
  else if (binx == "w") bin_on_x = COLUMN_W;
  else if (binx == "ch") bin_on_x = COLUMN_CH;
  else std::cerr << "Unrecognized x variable" << std::endl;

  min_x = params->get_min(bin_on_x);
  max_x = params->get_max(bin_on_x);


  aladyn_float x, y, px, py, pz, w, ch, gamma, theta, thetaT, E, ty, tz;
  aladyn_float going_to_bin_on_x = 0.;
  for (int i = 0; i < params->nparts; i++)
  {
    if (params->sim_is_2d)
    {
      x = *(parts + i*ndv);
      y = *(parts + i*ndv + 1);
      px = *(parts + i*ndv + 2);
      py = *(parts + i*ndv + 3);
      if (params->file_has_weight && !params->overwrite_weight)
        w = *(parts + i*ndv + 4);
      else
        w = params->overwrite_weight_value;
      if (params->file_has_charge && !params->overwrite_charge)
        ch = *(parts + i*ndv + 5);
      else
        ch = params->overwrite_charge_value;
      gamma = (aladyn_float)(sqrt(1. + px*px + py*py) - 1.);      //gamma-1
      theta = (aladyn_float)(atan2(py, px)*180. / M_PI);          //theta sgatto
      thetaT = (aladyn_float)atan(sqrt((py*py / (px*px))));       //theta turch
      E = (aladyn_float)(gamma*params->mass_MeV);
      if (px > 0.0) ty = py / px;
      else ty = 0.0;
      tz = 0.0;
      if (bin_on_x == COLUMN_X) going_to_bin_on_x = x;
      else if (bin_on_x == COLUMN_Y) going_to_bin_on_x = y;
      else if (bin_on_x == COLUMN_Z) std::cout << "Unable to bin on z in 2D" << std::endl;
      else if (bin_on_x == COLUMN_PX) going_to_bin_on_x = px;
      else if (bin_on_x == COLUMN_PY) going_to_bin_on_x = py;
      else if (bin_on_x == COLUMN_PZ) std::cout << "Unable to bin on pz in 2D" << std::endl;
      else if (bin_on_x == COLUMN_GAMMA) going_to_bin_on_x = gamma;
      else if (bin_on_x == COLUMN_THETA) going_to_bin_on_x = theta;
      else if (bin_on_x == COLUMN_E) going_to_bin_on_x = E;
      else if (bin_on_x == COLUMN_THETAT) going_to_bin_on_x = thetaT;
      else if (bin_on_x == COLUMN_TY) going_to_bin_on_x = ty;
      else if (bin_on_x == COLUMN_TZ) std::cout << "Unable to bin on tz in 2D" << std::endl;
      else if (bin_on_x == COLUMN_W) going_to_bin_on_x = w;
      else if (bin_on_x == COLUMN_CH) going_to_bin_on_x = ch;
    }
    else
    {
      px = *(parts + i*ndv + 3);
      py = *(parts + i*ndv + 4);
      pz = *(parts + i*ndv + 5);
      if (params->file_has_weight && !params->overwrite_weight)
        w = *(parts + i*ndv + 6);
      else
        w = params->overwrite_weight_value;
      if (params->file_has_charge && !params->overwrite_charge)
        ch = *(parts + i*ndv + 7);
      else
        ch = params->overwrite_charge_value;
      gamma = (aladyn_float)(sqrt(1. + px*px + py*py + pz*pz) - 1.);             //gamma
      theta = (aladyn_float)(atan2(sqrt(py*py + pz*pz), px)*180. / M_PI);        //theta sgatto
      thetaT = (aladyn_float)atan(sqrt((py*py / (px*px)) + (pz*pz / (px*px))));  //theta turch
      E = (aladyn_float)(gamma*params->mass_MeV);                 //energia
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
      if (bin_on_x < COLUMN_GAMMA) going_to_bin_on_x = *(parts + i*ndv + bin_on_x);
      else if (bin_on_x == COLUMN_GAMMA) going_to_bin_on_x = gamma;
      else if (bin_on_x == COLUMN_THETA) going_to_bin_on_x = theta;
      else if (bin_on_x == COLUMN_E) going_to_bin_on_x = E;
      else if (bin_on_x == COLUMN_THETAT) going_to_bin_on_x = thetaT;
      else if (bin_on_x == COLUMN_TY) going_to_bin_on_x = ty;
      else if (bin_on_x == COLUMN_TZ) going_to_bin_on_x = tz;
      else if (bin_on_x == COLUMN_W) going_to_bin_on_x = w;
      else if (bin_on_x == COLUMN_CH) going_to_bin_on_x = ch;
    }

    if (going_to_bin_on_x < min_x)
    {
      whichbin_x = 0;
    }
    else if (going_to_bin_on_x > max_x)
    {
      whichbin_x = params->get_nbin(bin_on_x) + 2;
    }
    else
    {
      whichbin_x = (int)(((going_to_bin_on_x - min_x) / params->get_size_(bin_on_x)) + 1.0);
    }

    data_binned[whichbin_x] += w;
  }
}

