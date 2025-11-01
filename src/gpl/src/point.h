// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2018-2025, The OpenROAD Authors

#pragma once

#include <sstream>
#include <string>

namespace gpl {

class FloatPoint
{
 public:
  float x = 0;
  float y = 0;
  FloatPoint() = default;
  FloatPoint(float x, float y) : x(x), y(y) {}
  std::string to_string() const
  {
    std::ostringstream oss;
    oss << "(" << x << ", " << y << ")";
    return oss.str();
  }
};

}  // namespace gpl
