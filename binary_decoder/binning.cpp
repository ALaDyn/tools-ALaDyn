
#include "binning.h"

_Binning::_Binning(aladyn_float * parts, size_t npart, Parameters * params, densityplot * plot)
{
  size_t bin_on_x = plot->get_x_column_to_bin(), bin_on_y = plot->get_y_column_to_bin();
  size_t whichbin_x = 0, whichbin_y = 0;
  size_t ndv = params->ndv;

  aladyn_float min_x = 0.0, min_y = 0.0, max_x = 1.0, max_y = 1.0;

  min_x = plot->min_value_x;
  max_x = plot->max_value_x;
  min_y = plot->min_value_y;
  max_y = plot->max_value_y;

  aladyn_float x, y, z, px, py, pz, w, ch, gamma, theta, thetaT, E, ty, tz;
  aladyn_float going_to_bin_on_x = 0., going_to_bin_on_y = 0.;
  for (int i = 0; i < npart; i++)
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
      whichbin_x = plot->nbin_x + 2;
    }
    else
    {
      whichbin_x = (int)(((going_to_bin_on_x - min_x) / plot->get_bin_size_x()) + 1.0);
    }
    if (going_to_bin_on_y < min_y)
    {
      whichbin_y = 0;
    }
    else if (going_to_bin_on_y > max_x)
    {
      whichbin_y = plot->nbin_y + 2;
    }
    else
    {
      whichbin_y = (int)(((going_to_bin_on_y - min_y) / plot->get_bin_size_y()) + 1.0);
    }

    plot->data[whichbin_x][whichbin_y] += w;
  }
}


_Binning::_Binning(aladyn_float * parts, size_t npart, Parameters * params, histo * plot)
{
  size_t bin_on_x = plot->get_column_to_bin();
  size_t whichbin_x = 0;
  size_t ndv = params->ndv;
  aladyn_float min_x = 0.0, max_x = 1.0;
  min_x = plot->min_value;
  max_x = plot->max_value;

  aladyn_float x, y, px, py, pz, w, ch, gamma, theta, thetaT, E, ty, tz;
  aladyn_float going_to_bin_on_x = 0.;
  for (int i = 0; i < npart; i++)
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
      whichbin_x = plot->nbin + 2;
    }
    else
    {
      whichbin_x = (int)(((going_to_bin_on_x - min_x) / plot->get_bin_size()) + 1.0);
    }

    plot->data[whichbin_x] += w;
  }
}

