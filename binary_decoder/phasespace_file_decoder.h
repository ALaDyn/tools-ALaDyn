#pragma once

#include "binary_decoder.h"
#include "binning.h"
#include "filter.h"
#include "output.h"
#include "swap_tools.h"


class Part_v1 {
  double x, y, z;
  double px, py, pz;
};
class Part_v2 {
  double x, y, z;
  double px, py, pz;
  double weight;
};
class Part_v3 {
  double x, y, z;
  double px, py, pz;
  float weight;
  int charge;
};
union Part_vX {
  Part_v1 part_v1;
  Part_v2 part_v2;
  Part_v3 part_v3;
};
union double_as_two_float {
  double d;
  float f[2]; //f[0] = peso, f[1] = carica
};




int read_phase_space_file(Parameters * );


