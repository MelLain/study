#pragma once

#include <memory>

#include "common.h"
#include "matrix.h"
#include "grid.h"

namespace DTS {

std::shared_ptr<DM> FivePointsLaplass(const DM& src, const Grid& grid);
double ProductByPointAndSum(const DM& src_1, const DM& src_2, const Grid& grid);

}  // namespace DTS
