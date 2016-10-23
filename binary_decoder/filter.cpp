
#include "filter.h"


#ifndef _WIN32
#include <strings.h>
#else
int    strcasecmp(const char* s1, const char* s2)
{
  for (;;)
  {
    int c1 = tolower(*((unsigned char*)s1++));
    int c2 = tolower(*((unsigned char*)s2++));

    if ((c1 != c2) || (c1 == '\0'))
    {
      return(c1 - c2);
    }
  }
}
#endif

namespace cost
{
  unsigned int xmin = __0X00;
  unsigned int ymin = __0X01;
  unsigned int zmin = __0X02;
  unsigned int pxmin = __0X03;
  unsigned int pymin = __0X04;
  unsigned int pzmin = __0X05;
  unsigned int xmax = __0X06;
  unsigned int ymax = __0X07;
  unsigned int zmax = __0X08;
  unsigned int pxmax = __0X09;
  unsigned int pymax = __0X10;
  unsigned int pzmax = __0X11;
  unsigned int emin = __0X12;
  unsigned int emax = __0X13;
  unsigned int thetamin = __0X14;
  unsigned int thetamax = __0X15;
  unsigned int thetaTmin = __0X16;
  unsigned int thetaTmax = __0X17;
  unsigned int tymin = __0X18;
  unsigned int tymax = __0X19;
  unsigned int tzmin = __0X20;
  unsigned int tzmax = __0X21;
  unsigned int wmin = __0X22;
  unsigned int wmax = __0X23;
  unsigned int chmin = __0X24;
  unsigned int chmax = __0X25;

  unsigned int all_filter_index[] =
  {
    xmin, ymin, zmin,
    pxmin, pymin, pzmin,
    xmax, ymax, zmax,
    pxmax, pymax, pzmax,
    emin, emax,
    thetamin, thetamax,
    thetaTmin, thetaTmax,
    tymin, tymax,
    tzmin, tzmax,
    wmin, wmax,
    chmin, chmax
  };
}


_Filter::_Filter(Parameters * params, aladyn_float *data, unsigned int n_data[], aladyn_float *val, unsigned int filter_mask)
{
  aladyn_float * pntt_loc, p[] = { 0, 0, 0 }, E = 0., theta = 0., thetaT = 0., ty = 0., tz = 0., w = 0., ch = 0.;
  unsigned int current_index = 0, tests[32];
  bool flag;
  unsigned char tot_test = 0;

  flag_filters = 0;
  if (!filter_mask) filter_mask = internal_mask;
  if (!filter_mask)
  {
    return;
  }
  for (unsigned char c = 0; c < 32; ++c)
  {
    unsigned int r = filter_mask & cost[c];
    if (!r) continue;
    tests[tot_test++] = cost[c];
  }
  for (unsigned int i = 0; i < n_data[0]; ++i)
  {
    pntt_loc = data + i*n_data[1];
    flag = true;

    if (params->sim_is_2d)
    {
      p[0] = pntt_loc[2], p[1] = pntt_loc[3], p[2] = 0.0;
      if (params->file_has_weight && !params->overwrite_weight)
        w = pntt_loc[4];
      else
        w = params->overwrite_weight_value;
      if (params->file_has_charge && !params->overwrite_charge)
        ch = pntt_loc[5];
      else
        ch = params->overwrite_charge_value;
    }
    else
    {
      p[0] = pntt_loc[3], p[1] = pntt_loc[4], p[2] = pntt_loc[5];
      if (params->file_has_weight && !params->overwrite_weight)
        w = pntt_loc[6];
      else
        w = params->overwrite_weight_value;
      if (params->file_has_charge && !params->overwrite_charge)
        ch = pntt_loc[7];
      else
        ch = params->overwrite_charge_value;
    }



    for (unsigned char c = 0; c < tot_test; ++c)
    {
      if (!flag) break;

      if (tests[c] == __0X12 || tests[c] == __0X13)
        E = (aladyn_float)(params->mass_MeV * (sqrtf(1.0f + p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) - 1.f));

      else if (tests[c] == __0X14 || tests[c] == __0X15)
        theta = (aladyn_float)(atan2(sqrt(p[1] * p[1] + p[2] * p[2]), p[0])*180. / M_PI);

      else if (tests[c] == __0X16 || tests[c] == __0X17)
        thetaT = (aladyn_float)atan(sqrt((p[1] * p[1] / (p[0] * p[0])) + (p[2] * p[2] / (p[0] * p[0]))));

      else if (tests[c] == __0X18 || tests[c] == __0X19)
        ty = p[1] / p[0];

      else if (tests[c] == __0X20 || tests[c] == __0X21)
        tz = p[2] / p[0];

      switch (tests[c])
      {
      case __0X00: // cost::xmin
        name = xmin;
        flag_filters.minus_xmin = pntt_loc[(int)name] >= val[(int)name];
        flag = flag && flag_filters.minus_xmin;
        break;
      case __0X01: // cost::ymin
        name = ymin;
        flag_filters.minus_ymin = pntt_loc[(int)name] >= val[(int)name];
        flag = flag && flag_filters.minus_ymin;
        break;
      case __0X02: // cost::zmin
        name = zmin;
        flag_filters.minus_zmin = pntt_loc[(int)name] >= val[(int)name];
        flag = flag && flag_filters.minus_zmin;
        break;
      case __0X03: // cost::pxmin
        name = pxmin;
        flag_filters.minus_pxmin = pntt_loc[(int)name] >= val[(int)name];
        flag = flag && flag_filters.minus_pxmin;
        break;
      case __0X04: // cost::pymin
        name = pymin;
        flag_filters.minus_pymin = pntt_loc[(int)name] >= val[(int)name];
        flag = flag && flag_filters.minus_pymin;
        break;
      case __0X05: // cost::pzmin
        name = pzmin;
        flag_filters.minus_pzmin = pntt_loc[(int)name] >= val[(int)name];
        flag = flag && flag_filters.minus_pzmin;
        break;
      case __0X06: // cost::xmax
        name = xmax;
        flag_filters.plus_xmax = pntt_loc[(int)name - 6] <= val[(int)name];
        flag = flag && flag_filters.plus_xmax;
        break;
      case __0X07: // cost::ymax
        name = ymax;
        flag_filters.plus_ymax = pntt_loc[(int)name - 6] <= val[(int)name];
        flag = flag && flag_filters.plus_ymax;
        break;
      case __0X08: // cost::zmax
        name = zmax;
        flag_filters.plus_zmax = pntt_loc[(int)name - 6] <= val[(int)name];
        flag = flag && flag_filters.plus_zmax;
        break;
      case __0X09: // cost::pxmax
        name = pxmax;
        flag_filters.plus_pxmax = pntt_loc[(int)name - 6] <= val[(int)name];
        flag = flag && flag_filters.plus_pxmax;
        break;
      case __0X10: // cost::pymax
        name = pymax;
        flag_filters.plus_pymax = pntt_loc[(int)name - 6] <= val[(int)name];
        flag = flag && flag_filters.plus_pymax;
        break;
      case __0X11: // cost::pzmax
        name = pzmax;
        flag_filters.plus_pzmax = pntt_loc[(int)name - 6] <= val[(int)name];
        flag = flag && flag_filters.plus_pzmax;
        break;
      case __0X12: // cost::emin
        name = emin;
        flag_filters.minus_Emin = E >= val[12];
        flag = flag && flag_filters.minus_Emin;
        break;
      case __0X13:  // cost::emax
        name = emax;
        flag_filters.plus_Emax = E <= val[13];
        flag = flag && flag_filters.plus_Emax;
        break;
      case __0X14: // cost::thetamin
        name = thetamin;
        flag_filters.minus_thetamin = theta >= val[14];
        flag = flag && flag_filters.minus_thetamin;
        break;
      case __0X15: // cost::thetamax
        name = thetamax;
        flag_filters.plus_thetamax = theta <= val[15];
        flag = flag && flag_filters.plus_thetamax;
        break;
      case __0X16: // cost::thetaTmin
        name = thetaTmin;
        flag_filters.minus_thetaTmin = thetaT >= val[16];
        flag = flag && flag_filters.minus_thetaTmin;
        break;
      case __0X17: // cost::thetaTmax
        name = thetaTmax;
        flag_filters.plus_thetaTmax = thetaT <= val[17];
        flag = flag && flag_filters.plus_thetaTmax;
        break;
      case __0X18: // cost::tymin
        name = tymin;
        flag_filters.minus_tymin = ty >= val[18];
        flag = flag && flag_filters.minus_tymin;
        break;
      case __0X19: // cost::tymax
        name = tymax;
        flag_filters.plus_tymax = ty <= val[19];
        flag = flag && flag_filters.plus_tymax;
        break;
      case __0X20: // cost::tzmin
        name = tzmin;
        flag_filters.minus_tzmin = tz >= val[20];
        flag = flag && flag_filters.minus_tzmin;
        break;
      case __0X21: // cost::tzmax
        name = tzmax;
        flag_filters.plus_tzmax = tz <= val[21];
        flag = flag && flag_filters.plus_tzmax;
        break;
      case __0X22: // cost::wmin
        name = wmin;
        flag_filters.minus_wmin = w >= val[22];
        flag = flag && flag_filters.minus_wmin;
        break;
      case __0X23: // cost::wmax
        name = wmax;
        flag_filters.plus_wmax = w <= val[23];
        flag = flag && flag_filters.plus_wmax;
        break;
      case __0X24: // cost::chmin
        name = chmin;
        flag_filters.minus_chmin = ch >= val[24];
        flag = flag && flag_filters.minus_chmin;
        break;
      case __0X25: // cost::chmax
        name = chmax;
        flag_filters.plus_chmax = ch <= val[25];
        flag = flag && flag_filters.plus_chmax;
        break;
      }

    }
    if (!flag) continue;
    for (unsigned int k = 0; k < n_data[1]; ++k) data[n_data[1] * current_index + k] = pntt_loc[k];
    current_index++;
  }
  n_data[0] = current_index;
  internal_mask = 0;
}

const char * _Filter::descr[] =
{
  "+xmin",
  "+ymin",
  "+zmin",
  "+xmax",
  "+ymax",
  "+zmax",
  "+pxmin",
  "+pymin",
  "+pzmin",
  "+pxmax",
  "+pymax",
  "+pzmax",
  "+Emin",
  "+Emax",
  "+thetamin",
  "+thetamax",
  "+thetaTmin",
  "+thetaTmax",
  "+tymin",
  "+tymax",
  "+tzmin",
  "+tzmax",
  "+wmin",
  "+wmax",
  "+chmin",
  "+chmax"
  //add here others if needed
};

const unsigned int _Filter::cost[] =
{
  __0X00, __0X01, __0X02, __0X03, __0X04, __0X05, __0X06, __0X07,
  __0X08, __0X09, __0X10, __0X11, __0X12, __0X13, __0X14, __0X15,
  __0X16, __0X17, __0X18, __0X19, __0X20, __0X21, __0X22, __0X23,
  __0X24, __0X25, __0X26, __0X27, __0X28, __0X29, __0X30, __0X31
};

aladyn_float * _Filter::build_filter(Parameters * params)
{
  char ** my_args;
  aladyn_float * my_vals;
  int index[ENABLED_FILTERS], counter = 0;
  for (int i = 1; i < params->argc; ++i)
  {
    if (params->argv[i][0] == '+')
      index[counter++] = i;
  }
  index[counter] = -1;
  if (!counter) return (aladyn_float *)nullptr;
  my_args = new char *[ENABLED_FILTERS + 1], my_args[ENABLED_FILTERS] = 0;
  my_vals = new aladyn_float[ENABLED_FILTERS + 1];
  for (int i = 0; i < counter; ++i)
    my_args[i] = const_cast<char*>(params->argv[index[i]].c_str()),
    my_vals[i] = (aladyn_float)atof(params->argv[index[i] + 1].c_str());
  my_args[counter] = 0;
  return build_filter(
    my_args[0], my_vals[0],
    my_args[1], my_vals[1],
    my_args[2], my_vals[2],
    my_args[3], my_vals[3],
    my_args[4], my_vals[4],
    my_args[5], my_vals[5],
    my_args[6], my_vals[6],
    my_args[7], my_vals[7],
    my_args[8], my_vals[8],
    my_args[9], my_vals[9],
    my_args[10], my_vals[10],
    my_args[11], my_vals[11],
    my_args[12], my_vals[12],
    my_args[13], my_vals[13],
    my_args[14], my_vals[14],
    my_args[15], my_vals[15],
    my_args[16], my_vals[16],
    my_args[17], my_vals[17],
    my_args[18], my_vals[18],
    my_args[19], my_vals[19],
    my_args[20], my_vals[20],
    my_args[21], my_vals[21],
    my_args[22], my_vals[22],
    my_args[23], my_vals[23],
    my_args[24], my_vals[24],
    my_args[25], my_vals[25],
    my_args[ENABLED_FILTERS]);
}


aladyn_float * _Filter::build_filter(const char *p, ...)
{
  va_list app;
  char * buff = new char[256], *z = buff;
  aladyn_float val, *all_vals = new aladyn_float[ENABLED_FILTERS];
  va_start(app, p);
  strcpy(buff, p);
  do
  {
    val = (aladyn_float)va_arg(app, double);
    find_filter(buff, val, all_vals);
    if (!all_vals) return (aladyn_float*)nullptr;
    z = va_arg(app, char *);
    if (!z) break;
    strcpy(buff, z);
  } while (1);
  va_end(app);
  return all_vals;
}


void _Filter::find_filter(char *b, aladyn_float v, aladyn_float *& V)
{
  int i;
  for (i = 0; i < ENABLED_FILTERS; ++i) if (!strcasecmp(b, descr[i])) break;
  if (i >= ENABLED_FILTERS)
  {
    V = (aladyn_float*)nullptr; // non è detto che sia la cosa migliore da farsi
    return;
  }
  V[i] = v;
  internal_mask |= cost::all_filter_index[i];
}

unsigned int _Filter::internal_mask = 0;

