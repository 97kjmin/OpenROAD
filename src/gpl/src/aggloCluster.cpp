#include "aggloCluster.h"

// Standard Library
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

// OpenMP
#include <omp.h>

// Boost Geometry
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>

// OR-Tools
#include "ortools/graph/graph.h"
#include "ortools/graph/linear_assignment.h"

// OpenDB
#include "odb/db.h"

// OpenSTA
#include "db_sta/dbNetwork.hh"
#include "db_sta/dbSta.hh"
#include "sta/ExceptionPath.hh"
#include "sta/FuncExpr.hh"
#include "sta/Graph.hh"
#include "sta/Liberty.hh"
#include "sta/MinMax.hh" 
#include "sta/Network.hh"
#include "sta/Path.hh"
#include "sta/PathAnalysisPt.hh"
#include "sta/PathEnd.hh"
#include "sta/PathExpanded.hh"
#include "sta/Search.hh"
#include "sta/Sequential.hh"
#include "sta/Sta.hh"
#include "sta/Units.hh"

// OpenROAD Modules
#include "rsz/Resizer.hh"
#include "utl/Logger.h"

// GPL Internal
#include "placerBase.h" 

namespace gpl {

//==============================================================================
// MasterMask Implementation
//==============================================================================

std::string 
MasterMask::to_string() const
{
  std::ostringstream oss;
  oss << "[F:" << func_id_ << "]"
      << " [Clk:" << (has_positive_clock_edge_ ? "+" : "-") << "]"
      << " [Clr:" << (has_clear_ ? "Y" : "N") << "]"
      << " [Pre:" << (has_preset_ ? "Y" : "N") << "]"
      << " [Q:" << (has_q_pin_ ? "Y" : "N") << "]"
      << " [QN:" << (has_qn_pin_ ? "Y" : "N") << "]"
      << " [Scan:" << (has_scan_ ? "Y" : "N") << "]";
  return oss.str();
}

//==============================================================================
// InstMask Implementation
//==============================================================================

std::string 
InstMask::to_string() const
{
  std::ostringstream oss;
  oss << "[Clk:" << (clock_net_ ? clock_net_->getName() : "N/A") << "]"
      << " [Clr:" << (clear_net_ ? clear_net_->getName() : "N/A") << "]"
      << " [Pre:" << (preset_net_ ? preset_net_->getName() : "N/A") << "]"
      << " [SE:" << (scan_enable_net_ ? scan_enable_net_->getName() : "N/A") << "]"
      << " [SI:" << (scan_in_net_ ? scan_in_net_->getName() : "N/A") << "]";
  return oss.str();
}

//==============================================================================
// FlopUnit Implementation
//==============================================================================

std::string 
FlopUnit::to_string() const
{
  std::ostringstream oss;
  oss << "Inst: " << inst_->getName() << "\n"
      << "  - Master Mask: " << master_mask_.to_string() << "\n"
      << "  - Inst Mask:   " << inst_mask_.to_string() << "\n"
      << "  - Orig Pt:     " << orig_pt_.to_string();
  return oss.str();
}

//==============================================================================
// VirtualBin Implementation
//==============================================================================

int64_t
VirtualBin::getOverflowArea() const
{
  const float result = static_cast<float>(inst_placed_area_)
                     + static_cast<float>(macro_placed_area_) * target_density_
                     + static_cast<float>(non_place_area_) * target_density_
                     - static_cast<float>(bin_area_) * target_density_;
  return static_cast<int64_t>(std::max(0.0f, result));
}

//==============================================================================
// VirtualBinGrid Implementation
//==============================================================================

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

VirtualBinGrid::VirtualBinGrid(int lx, int ly, int ux, int uy,
                               int bin_cnt_x, int bin_cnt_y,
                               double bin_size_x, double bin_size_y,
                               float target_density, float target_overflow)
    : lx_(lx),
      ly_(ly),
      ux_(ux),
      uy_(uy),
      bin_cnt_x_(bin_cnt_x),
      bin_cnt_y_(bin_cnt_y),
      bin_size_x_(bin_size_x),
      bin_size_y_(bin_size_y),
      target_density_(target_density),
      target_overflow_(target_overflow)
{
  bins_.reserve(bin_cnt_x_ * bin_cnt_y_);
  for (int y = 0; y < bin_cnt_y_; ++y) {
    for (int x = 0; x < bin_cnt_x_; ++x) {
      const int bin_lx = lx_ + std::lround(x * bin_size_x_);
      const int bin_ly = ly_ + std::lround(y * bin_size_y_);
      const int bin_ux = lx_ + std::lround((x + 1) * bin_size_x_);
      const int bin_uy = ly_ + std::lround((y + 1) * bin_size_y_);
      bins_.emplace_back(bin_lx, bin_ly, bin_ux, bin_uy, target_density_);
    }
  }
}

//------------------------------------------------------------------------------
// Overflow Check
//------------------------------------------------------------------------------

bool 
VirtualBinGrid::checkOverflow()
{
  int64_t total_overflow_area = 0;
  int64_t total_inst_placed_area = 0;
  int64_t total_macro_placed_area = 0;
  int64_t total_non_place_area = 0;

  for (const auto& bin : bins_) {
    total_overflow_area += bin.getOverflowArea();
    total_inst_placed_area += bin.getInstPlacedArea();
    total_macro_placed_area += bin.getMacroPlacedArea();
    total_non_place_area += bin.getNonPlaceArea();
  }

  const float total_occupied_area = static_cast<float>(total_inst_placed_area)
                                   + static_cast<float>(total_macro_placed_area) * target_density_
                                   + static_cast<float>(total_non_place_area) * target_density_;

  const float overflow_ratio = static_cast<float>(total_overflow_area) / total_occupied_area;

  return overflow_ratio > target_overflow_;
}

//------------------------------------------------------------------------------
// Bin Index Calculation
//------------------------------------------------------------------------------

std::pair<int, int> 
VirtualBinGrid::getMinMaxIdxX(const odb::Rect& box) const
{
  const int lower_idx = static_cast<int>(std::floor((box.xMin() - lx_) / bin_size_x_));
  const int upper_idx = static_cast<int>(std::ceil((box.xMax() - lx_) / bin_size_x_));
  
  return {std::max(0, lower_idx), std::min(bin_cnt_x_, upper_idx)};
}

std::pair<int, int> 
VirtualBinGrid::getMinMaxIdxY(const odb::Rect& box) const
{
  const int lower_idx = static_cast<int>(std::floor((box.yMin() - ly_) / bin_size_y_));
  const int upper_idx = static_cast<int>(std::ceil((box.yMax() - ly_) / bin_size_y_));
  
  return {std::max(0, lower_idx), std::min(bin_cnt_y_, upper_idx)};
}

//------------------------------------------------------------------------------
// Instance Placement Update
//------------------------------------------------------------------------------

void
VirtualBinGrid::accumulateInstPlacement(odb::dbMaster* master,
                                        const Point& center_pos,
                                        std::vector<int64_t>& area_deltas,
                                        bool is_add) const
{
  if (!master || bins_.empty()) {
    return;
  }

  if (area_deltas.size() != bins_.size()) {
    return;
  }

  const int w = master->getWidth();
  const int h = master->getHeight();
  const int lx = center_pos.get<0>() - w / 2;
  const int ly = center_pos.get<1>() - h / 2;
  const odb::Rect bbox(lx, ly, lx + w, ly + h);

  const auto [min_x_idx, max_x_idx] = getMinMaxIdxX(bbox);
  const auto [min_y_idx, max_y_idx] = getMinMaxIdxY(bbox);

  const int sign = is_add ? 1 : -1;

  for (int y = min_y_idx; y < max_y_idx; ++y) {
    for (int x = min_x_idx; x < max_x_idx; ++x) {
      const int bin_idx = y * bin_cnt_x_ + x;
      if (bin_idx < 0 || bin_idx >= static_cast<int>(area_deltas.size())) {
        continue;
      }

      const VirtualBin& bin = bins_[bin_idx];

      const int overlap_lx = std::max(bbox.xMin(), bin.lx());
      const int overlap_ly = std::max(bbox.yMin(), bin.ly());
      const int overlap_ux = std::min(bbox.xMax(), bin.ux());
      const int overlap_uy = std::min(bbox.yMax(), bin.uy());

      if (overlap_ux > overlap_lx && overlap_uy > overlap_ly) {
        const int64_t overlap_area = static_cast<int64_t>(overlap_ux - overlap_lx)
                                    * static_cast<int64_t>(overlap_uy - overlap_ly);
        area_deltas[bin_idx] += sign * overlap_area;
      }
    }
  }
}

void 
VirtualBinGrid::updateInstPlacement(odb::dbMaster* master, 
                                    const Point& center_pos, 
                                    bool is_add)
{
  const int w = master->getWidth();
  const int h = master->getHeight();
  const int lx = center_pos.get<0>() - w / 2;
  const int ly = center_pos.get<1>() - h / 2;
  const odb::Rect bbox(lx, ly, lx + w, ly + h);

  const auto [min_x_idx, max_x_idx] = getMinMaxIdxX(bbox);
  const auto [min_y_idx, max_y_idx] = getMinMaxIdxY(bbox);

  for (int y = min_y_idx; y < max_y_idx; ++y) {
    for (int x = min_x_idx; x < max_x_idx; ++x) {
      const int bin_idx = y * bin_cnt_x_ + x;
      VirtualBin& bin = bins_[bin_idx];

      const int overlap_lx = std::max(bbox.xMin(), bin.lx());
      const int overlap_ly = std::max(bbox.yMin(), bin.ly());
      const int overlap_ux = std::min(bbox.xMax(), bin.ux());
      const int overlap_uy = std::min(bbox.yMax(), bin.uy());

      if (overlap_ux > overlap_lx && overlap_uy > overlap_ly) {
        const int64_t overlap_area = static_cast<int64_t>(overlap_ux - overlap_lx)
                                    * static_cast<int64_t>(overlap_uy - overlap_ly);

        if (is_add) {
          bin.addInstPlacedArea(overlap_area);
        } else {
          bin.subInstPlacedArea(overlap_area);
        }
      }
    }
  }
}

void 
VirtualBinGrid::addInst(odb::dbMaster* master, const Point& center_pos)
{
  updateInstPlacement(master, center_pos, true);
}

void 
VirtualBinGrid::removeInst(odb::dbMaster* master, const Point& center_pos)
{
  updateInstPlacement(master, center_pos, false);
}

bool
VirtualBinGrid::wouldOverflow(odb::dbMaster* master1, const Point& pos1,
                              odb::dbMaster* master2, const Point& pos2,
                              odb::dbMaster* new_master, const Point& new_pos) const
{
  if (bins_.empty()) {
    return false;
  }

  std::vector<int64_t> area_deltas(bins_.size(), 0);
  accumulateInstPlacement(master1, pos1, area_deltas, false);
  accumulateInstPlacement(master2, pos2, area_deltas, false);
  accumulateInstPlacement(new_master, new_pos, area_deltas, true);

  int64_t total_overflow_area = 0;
  int64_t total_inst_area = 0;
  int64_t total_macro_area = 0;
  int64_t total_non_place_area = 0;

  for (size_t idx = 0; idx < bins_.size(); ++idx) {
    const VirtualBin& bin = bins_[idx];
    const int64_t inst_area = std::max<int64_t>(0, bin.getInstPlacedArea() + area_deltas[idx]);

    total_inst_area += inst_area;
    total_macro_area += bin.getMacroPlacedArea();
    total_non_place_area += bin.getNonPlaceArea();

    const float overflow = static_cast<float>(inst_area)
                           + static_cast<float>(bin.getMacroPlacedArea()) * bin.getTargetDensity()
                           + static_cast<float>(bin.getNonPlaceArea()) * bin.getTargetDensity()
                           - static_cast<float>(bin.getBinArea()) * bin.getTargetDensity();
    total_overflow_area += static_cast<int64_t>(std::max(0.0f, overflow));
  }

  const float total_occupied_area = static_cast<float>(total_inst_area)
                                   + static_cast<float>(total_macro_area) * target_density_
                                   + static_cast<float>(total_non_place_area) * target_density_;

  if (total_occupied_area <= 0.0f) {
    return false;
  }

  const float overflow_ratio = static_cast<float>(total_overflow_area) / total_occupied_area;
  return overflow_ratio > target_overflow_;
}

void
VirtualBinGrid::applyMerge(odb::dbMaster* master1, const Point& pos1,
                           odb::dbMaster* master2, const Point& pos2,
                           odb::dbMaster* new_master, const Point& new_pos)
{
  if (master1) {
    removeInst(master1, pos1);
  }

  if (master2) {
    removeInst(master2, pos2);
  }

  if (new_master) {
    addInst(new_master, new_pos);
  }
}

//==============================================================================
// AggloCluster Implementation
//==============================================================================

// Constructor
AggloCluster::AggloCluster(odb::dbDatabase* db,
                           sta::dbSta* sta,
                           utl::Logger* log,
                           rsz::Resizer* resizer,
                           float target_density,
                           float target_overflow,
                           float region_scale_factor,
                           int num_paths_per_endpoint,
                           int threads,
                           int num_samples,
                           bool verbose)
    : db_(db),
      block_(db->getChip()->getBlock()),
      sta_(sta),
      network_(sta_->getDbNetwork()),
      corner_(sta_->cmdCorner()),
      resizer_(resizer),
      log_(log),
      threads_(threads),
      num_paths_per_endpoint_(num_paths_per_endpoint),
      num_samples_(num_samples),
      target_density_(target_density),
      target_overflow_(target_overflow),
      region_scale_factor_(region_scale_factor),
      verbose_(verbose)
{
}

// Destructor
AggloCluster::~AggloCluster() = default;

//------------------------------------------------------------------------------
// Main Entry Point
//------------------------------------------------------------------------------

void 
AggloCluster::doAggloCluster()
{
  using Clock = std::chrono::steady_clock;

  const auto overall_start = Clock::now();
  std::vector<std::pair<std::string, double>> phase_timings;
  phase_timings.reserve(11);

  const int64_t initial_hpwl = computeTotalHpwl();

  auto measure_phase = [&](const std::string& label, auto&& phase_fn) {
    const auto phase_start = Clock::now();
    phase_fn();
    const auto phase_end = Clock::now();
    const double seconds = std::chrono::duration<double>(phase_end - phase_start).count();
    phase_timings.emplace_back(label, seconds);
  };

  measure_phase("Phase 1: initVirtualBinGrid", [&] { initVirtualBinGrid(); });
  measure_phase("Phase 2: readCompatibleMasters", [&] { readCompatibleMasters(); });
  measure_phase("Phase 3: readFlopUnits", [&] { readFlopUnits(); });
  measure_phase("Phase 4: createCompatibleGroups", [&] { createCompatibleGroups(); });
  measure_phase("Phase 5: readTimingPaths", [&] { readTimingPaths(); });
  measure_phase("Phase 6: analyzeTimingPaths", [&] { analyzeTimingPaths(); });
  measure_phase("Phase 7: calcFeasibleRegions", [&] { calcFeasibleRegions(); });
  measure_phase("Phase 8: createFlopClusters", [&] { createFlopClusters(); });
  measure_phase("Phase 9: createCompatibilityGraph", [&] { createCompatibilityGraph(); });
  measure_phase("Phase 10: runAgglomerativeClustering", [&] { runAgglomerativeClustering(); });
  measure_phase("Phase 11: implementClusters", [&] { implementClusters(); });

  const auto overall_end = Clock::now();
  const double total_seconds = std::chrono::duration<double>(overall_end - overall_start).count();

  const int64_t final_hpwl = computeTotalHpwl();
  const int64_t delta_hpwl = final_hpwl - initial_hpwl;

  if (verbose_) {
    const auto old_flags = std::cout.flags();
    const auto old_precision = std::cout.precision();

    std::cout << "\n[doAggloCluster] Phase Timing Summary" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    for (const auto& [label, seconds] : phase_timings) {
      std::cout << "  - " << label << ": " << seconds << " s" << std::endl;
    }
    std::cout << "  - Total runtime: " << total_seconds << " s" << std::endl;

    std::cout << std::setprecision(2);
    const double initial_hpwl_um = block_->dbuToMicrons(initial_hpwl);
    const double final_hpwl_um = block_->dbuToMicrons(final_hpwl);
    const double delta_hpwl_um = block_->dbuToMicrons(delta_hpwl);

    std::cout << "\n[doAggloCluster] HPWL Summary" << std::endl;
    std::cout << "  - Initial HPWL: " << initial_hpwl_um << " um" << std::endl;
    std::cout << "  - Final HPWL:   " << final_hpwl_um << " um" << std::endl;
    std::cout << "  - Delta HPWL:   " << delta_hpwl_um << " um" << std::endl;

    std::cout.flags(old_flags);
    std::cout.precision(old_precision);
  }
}

//==============================================================================
// Phase 1: initVirtualBinGrid()
//==============================================================================

void 
AggloCluster::initVirtualBinGrid()
{
  if (verbose_) {
    std::cout << "\n[initVirtualBinGrid] Starting virtual bin grid initialization..." << std::endl;
  }

  // Step 1: Get core area information
  const odb::Rect& core_area = block_->getCoreArea();
  const int64_t total_core_area = static_cast<int64_t>(core_area.dx()) 
                                  * static_cast<int64_t>(core_area.dy());

  if (verbose_) {
    std::cout << "[Step 1] Core Area Calculation" << std::endl;
    std::cout << "  - Core bounds: (" << core_area.xMin() << ", " << core_area.yMin() 
              << ") -> (" << core_area.xMax() << ", " << core_area.yMax() << ")" << std::endl;
    std::cout << "  - Core dimensions: " << core_area.dx() << " x " << core_area.dy() << std::endl;
  }

  // Step 2: Calculate average placeable instance area
  int64_t placeable_area = 0;
  int placeable_count = 0;
  
  for (odb::dbInst* inst : block_->getInsts()) {
    if (!inst->getPlacementStatus().isFixed()) {
      placeable_area += inst->getBBox()->getBox().area();
      placeable_count++;
    }
  }

  const int64_t avg_inst_area = placeable_area / placeable_count;

  if (verbose_) {
    std::cout << "[Step 2] Placeable Instance Analysis" << std::endl;
    std::cout << "  - Placeable instance count: " << placeable_count << std::endl;
    std::cout << "  - Total placeable area: " << placeable_area << " DBU^2" << std::endl;
    std::cout << "  - Average instance area: " << avg_inst_area << " DBU^2" << std::endl;
  }

  // Step 3: Determine ideal bin count based on target density
  const int64_t ideal_bin_area = (target_density_ > 0)
      ? std::round(static_cast<float>(avg_inst_area) / target_density_)
      : 0;

  int ideal_bin_cnt = (ideal_bin_area > 0) 
      ? total_core_area / ideal_bin_area 
      : 4;
  ideal_bin_cnt = std::max(ideal_bin_cnt, 4);

  if (verbose_) {
    std::cout << "[Step 3] Ideal Bin Count Determination" << std::endl;
    std::cout << "  - Target density: " << target_density_ << std::endl;
    std::cout << "  - Ideal bin area: " << ideal_bin_area << " DBU^2" << std::endl;
    std::cout << "  - Ideal bin count: " << ideal_bin_cnt << std::endl;
  }

  // Step 4: Calculate bin grid dimensions (power-of-2 based for efficiency)
  const int width = core_area.dx();
  const int height = core_area.dy();
  const int aspect_ratio = roundDownToPowerOfTwo(std::max(width, height) 
                                                  / std::min(width, height));

  int base_bin_cnt = 2;
  for (; base_bin_cnt <= 1024; base_bin_cnt *= 2) {
    const int total_bins = base_bin_cnt * (base_bin_cnt * aspect_ratio);
    if ((base_bin_cnt == 2 || total_bins <= ideal_bin_cnt) 
        && 4 * total_bins > ideal_bin_cnt) {
      break;
    }
  }

  const int bin_cnt_x = (width > height) ? base_bin_cnt * aspect_ratio : base_bin_cnt;
  const int bin_cnt_y = (width > height) ? base_bin_cnt : base_bin_cnt * aspect_ratio;

  const double bin_size_x = static_cast<double>(width) / bin_cnt_x;
  const double bin_size_y = static_cast<double>(height) / bin_cnt_y;

  if (verbose_) {
    std::cout << "[Step 4] Bin Grid Dimensions Calculation" << std::endl;
    std::cout << "  - Core aspect ratio (rounded to power-of-2): " << aspect_ratio << std::endl;
    std::cout << "  - Base bin count: " << base_bin_cnt << std::endl;
    std::cout << "  - Final bin grid: " << bin_cnt_x << " x " << bin_cnt_y 
              << " = " << (bin_cnt_x * bin_cnt_y) << " bins" << std::endl;
    std::cout << "  - Bin size: " << bin_size_x << " x " << bin_size_y << " DBU" << std::endl;
  }

  // Step 5: Create virtual bin grid
  virtual_bin_grid_ = VirtualBinGrid(core_area.xMin(),
                                     core_area.yMin(),
                                     core_area.xMax(),
                                     core_area.yMax(),
                                     bin_cnt_x,
                                     bin_cnt_y,
                                     bin_size_x,
                                     bin_size_y,
                                     target_density_,
                                     target_overflow_);

  if (verbose_) {
    std::cout << "[Step 5] Virtual Bin Grid Created" << std::endl;
    std::cout << "  - Target overflow threshold: " << target_overflow_ << std::endl;
    std::cout << "  - Grid initialized with " << virtual_bin_grid_.getBins().size() << " bins" << std::endl;
  }

  // Step 6: Populate bins with existing instance areas
  populateBinGrid();

  // Step 7: Check current overflow status
  if (verbose_) {
    std::cout << "\n[Step 7] Current Overflow Analysis" << std::endl;
    
    int64_t total_overflow_area = 0;
    int64_t total_inst_placed_area = 0;
    int64_t total_macro_placed_area = 0;
    int64_t total_non_place_area = 0;
    int overflow_bin_count = 0;
    
    const auto& bins = virtual_bin_grid_.getBins();
    for (const auto& bin : bins) {
      const int64_t overflow = bin.getOverflowArea();
      if (overflow > 0) {
        overflow_bin_count++;
      }
      total_overflow_area += overflow;
      total_inst_placed_area += bin.getInstPlacedArea();
      total_macro_placed_area += bin.getMacroPlacedArea();
      total_non_place_area += bin.getNonPlaceArea();
    }
    
    const float total_occupied_area = static_cast<float>(total_inst_placed_area)
                                     + static_cast<float>(total_macro_placed_area) * target_density_
                                     + static_cast<float>(total_non_place_area) * target_density_;
    
    const float overflow_ratio = (total_occupied_area > 0) 
        ? static_cast<float>(total_overflow_area) / total_occupied_area 
        : 0.0f;
    
    std::cout << "  - Total overflow area: " << total_overflow_area << " DBU^2" << std::endl;
    std::cout << "  - Total occupied area: " << static_cast<int64_t>(total_occupied_area) << " DBU^2" << std::endl;
    std::cout << "  - Current overflow ratio: " << std::fixed << std::setprecision(4) 
              << overflow_ratio << " (" << (overflow_ratio * 100) << "%)" << std::endl;
    std::cout << "  - Target overflow threshold: " << target_overflow_ 
              << " (" << (target_overflow_ * 100) << "%)" << std::endl;
    std::cout << "  - Bins with overflow: " << overflow_bin_count << " / " << bins.size()
              << " (" << std::fixed << std::setprecision(1) 
              << (100.0 * overflow_bin_count / bins.size()) << "%)" << std::endl;
    std::cout << "  - Overflow status: " << (overflow_ratio > target_overflow_ ? "EXCEEDED" : "OK") << std::endl;
  }

  if (verbose_) {
    std::cout << "[initVirtualBinGrid] Completed successfully.\n" << std::endl;
  }
}

void
AggloCluster::populateBinGrid()
{
  auto& bins = virtual_bin_grid_.getBins();
  const int bin_cnt_x = virtual_bin_grid_.getBinCntX();
  
  int64_t total_inst_area = 0;
  int64_t total_macro_area = 0;
  int64_t total_fixed_area = 0;
  int processed_inst_count = 0;

  for (odb::dbInst* inst : block_->getInsts()) {
    const odb::Rect inst_bbox = inst->getBBox()->getBox();
    const auto [min_x, max_x] = virtual_bin_grid_.getMinMaxIdxX(inst_bbox);
    const auto [min_y, max_y] = virtual_bin_grid_.getMinMaxIdxY(inst_bbox);

    const bool show_debug = verbose_ && (processed_inst_count < num_samples_);

    if (show_debug) {
      std::cout << "\n[Sample Instance #" << processed_inst_count << "] " << inst->getName() << std::endl;
      std::cout << "  BBox: (" << inst_bbox.xMin() << ", " << inst_bbox.yMin() 
                << ") -> (" << inst_bbox.xMax() << ", " << inst_bbox.yMax() << ")" << std::endl;
      std::cout << "  Area: " << inst_bbox.area() << " DBU^2" << std::endl;
      std::cout << "  Status: " << (inst->getPlacementStatus().isFixed() ? "FIXED" : "PLACEABLE") << std::endl;
      std::cout << "  Type: " << inst->getMaster()->getType().getString() << std::endl;
      std::cout << "  Bin range: X[" << min_x << ", " << max_x << "), Y[" << min_y << ", " << max_y << ")" << std::endl;
    }

    int64_t inst_total_overlap = 0;
    int bin_overlap_count = 0;

    for (int y = min_y; y < max_y; ++y) {
      for (int x = min_x; x < max_x; ++x) {
        VirtualBin& bin = bins[y * bin_cnt_x + x];
        
        const odb::Rect bin_rect(bin.lx(), bin.ly(), bin.ux(), bin.uy());
        const odb::Rect overlap = inst_bbox.intersect(bin_rect);
        
        if (overlap.area() == 0) {
          continue;
        }

        const int64_t overlap_area = static_cast<int64_t>(overlap.area());
        inst_total_overlap += overlap_area;
        bin_overlap_count++;
        
        if (show_debug) {
          std::cout << "    Bin[" << x << "," << y << "]: (" 
                    << bin.lx() << ", " << bin.ly() << ") -> (" 
                    << bin.ux() << ", " << bin.uy() << "), overlap = " << overlap_area << " DBU^2" << std::endl;
        }
        
        if (inst->getPlacementStatus().isFixed()) {
          bin.addNonPlaceArea(overlap_area);
          total_fixed_area += overlap_area;
        } else if (inst->getMaster()->getType() == odb::dbMasterType::BLOCK) {
          bin.addMacroPlacedArea(overlap_area);
          total_macro_area += overlap_area;
        } else {
          bin.addInstPlacedArea(overlap_area);
          total_inst_area += overlap_area;
        }
      }
    }

    if (show_debug) {
      std::cout << "  Overlapping bins: " << bin_overlap_count << std::endl;
      std::cout << "  Total overlap: " << inst_total_overlap << " DBU^2" << std::endl;
      
      std::string classification;
      if (inst->getPlacementStatus().isFixed()) {
        classification = "FIXED/BLOCKAGE";
      } else if (inst->getMaster()->getType() == odb::dbMasterType::BLOCK) {
        classification = "MACRO";
      } else {
        classification = "STANDARD CELL";
      }
      std::cout << "  Classification: " << classification << std::endl;
    }

    processed_inst_count++;
  }

  if (verbose_) {
    const odb::Rect& core_area = block_->getCoreArea();
    const int64_t total_core_area = static_cast<int64_t>(core_area.dx()) 
                                    * static_cast<int64_t>(core_area.dy());
    
    std::cout << "\n[Step 6] Bin Population Completed" << std::endl;
    std::cout << "  - Total instances processed: " << processed_inst_count << std::endl;
    std::cout << "  - Standard cell area: " << total_inst_area << " DBU^2" << std::endl;
    std::cout << "  - Macro area: " << total_macro_area << " DBU^2" << std::endl;
    std::cout << "  - Fixed/blockage area: " << total_fixed_area << " DBU^2" << std::endl;
    std::cout << "  - Total occupied area: " << (total_inst_area + total_macro_area + total_fixed_area) << " DBU^2" << std::endl;
    std::cout << "  - Core utilization: " << std::fixed << std::setprecision(2) 
              << (100.0 * (total_inst_area + total_macro_area + total_fixed_area) / total_core_area) << "%" << std::endl;
  }
}

//==============================================================================
// Phase 2: readCompatibleMasters()
//==============================================================================

void 
AggloCluster::readCompatibleMasters()
{
  if (verbose_) {
    std::cout << "\n[readCompatibleMasters] Starting compatible master analysis..." << std::endl;
  }

  compatible_masters_.clear();
  func_str_to_func_id_.clear();
  representative_masters_.clear();

  // Step 1: Scan all library masters and group by (MasterMask, bit_width)
  if (verbose_) {
    std::cout << "[Step 1] Scanning library masters..." << std::endl;
  }

  std::map<MasterMask, std::map<int, std::vector<odb::dbMaster*>>> master_groups;
  std::map<MasterMask, std::map<int, std::set<std::string>>> master_name_tracker;
  const char* temp_inst_name = "_temp_master_check";

  int total_masters = 0;
  int valid_flop_masters = 0;
  int duplicate_masters = 0;

  for (odb::dbLib* lib : db_->getLibs()) {
    for (odb::dbMaster* master : lib->getMasters()) {
      total_masters++;

      // Create temporary instance to check if master is a valid flop
      odb::dbInst* temp_inst = odb::dbInst::create(block_, master, temp_inst_name);
      if (!temp_inst) {
        continue;
      }

      if (isValidFlop(temp_inst)) {
        const MasterMask mask = createMasterMask(temp_inst);
        const int bit_width = getNumDPins(temp_inst);
        const std::string master_name = master->getName();
        
        // Check for duplicate master names
        if (master_name_tracker[mask][bit_width].count(master_name) > 0) {
          duplicate_masters++;
        } else {
          // Add only if not duplicate
          valid_flop_masters++;
          master_groups[mask][bit_width].push_back(master);
          master_name_tracker[mask][bit_width].insert(master_name);
        }
      }

      odb::dbInst::destroy(temp_inst);
    }
  }

  if (verbose_) {
    std::cout << "  - Total masters scanned: " << total_masters << std::endl;
    std::cout << "  - Valid flop masters found: " << valid_flop_masters << std::endl;
    std::cout << "  - Duplicate masters skipped: " << duplicate_masters << std::endl;
    std::cout << "  - Unique master masks: " << master_groups.size() << std::endl;
  }

  // Step 2: Filter groups that support multi-bit clustering (>= 2 bit widths)
  if (verbose_) {
    std::cout << "\n[Step 2] Filtering multi-bit compatible groups..." << std::endl;
  }

  for (const auto& [mask, bit_map] : master_groups) {
    if (bit_map.size() >= 2) {
      compatible_masters_[mask] = bit_map;
    }
  }

  if (verbose_) {
    std::cout << "  - Compatible master groups: " << compatible_masters_.size() << std::endl;
    
    // Show bit width distribution
    std::map<size_t, int> bit_width_dist;
    for (const auto& [mask, bit_map] : compatible_masters_) {
      bit_width_dist[bit_map.size()]++;
    }
    
    std::cout << "  - Bit width support distribution:" << std::endl;
    for (const auto& [num_widths, count] : bit_width_dist) {
      std::cout << "    Groups supporting " << num_widths << " bit widths: " << count << std::endl;
    }
  }

  // Step 3: Select representative master (largest area) for each bit width
  if (verbose_) {
    std::cout << "\n[Step 3] Selecting representative masters..." << std::endl;
  }

  for (const auto& [mask, bit_map] : compatible_masters_) {
    std::map<int, odb::dbMaster*> largest_masters;

    for (const auto& [bit_width, masters] : bit_map) {
      if (masters.empty()) {
        continue;
      }

      // Find master with largest area (width * height)
      auto largest_it = std::max_element(
          masters.begin(), 
          masters.end(),
          [](odb::dbMaster* a, odb::dbMaster* b) {
            return a->getWidth() * a->getHeight() < b->getWidth() * b->getHeight();
          });

      largest_masters[bit_width] = *largest_it;
    }

    representative_masters_[mask] = largest_masters;
  }

  if (verbose_) {
    std::cout << "  - Total representative groups: " << representative_masters_.size() << std::endl;

    // Sample output for first few groups
    int group_count = 0;

    for (const auto& [mask, bit_map] : compatible_masters_) {
      if (group_count >= num_samples_) {
        break;
      }

      std::cout << "\n[Sample Group #" << group_count << "]" << std::endl;
      std::cout << "  Master Mask: " << mask.to_string() << std::endl;
      std::cout << "  Supported bit widths: ";
      for (const auto& [bit_width, masters] : bit_map) {
        std::cout << bit_width << "(" << masters.size() << " masters) ";
      }
      std::cout << std::endl;

      // Show all masters and highlight the representative
      const auto& repr_map = representative_masters_[mask];
      std::cout << "  Masters by bit width:" << std::endl;
      
      for (const auto& [bit_width, masters] : bit_map) {
        std::cout << "    " << bit_width << "-bit (" << masters.size() << " total):" << std::endl;
        
        // Get representative master for this bit width
        odb::dbMaster* repr_master = repr_map.at(bit_width);
        const int repr_area = repr_master->getWidth() * repr_master->getHeight();
        
        // Show all masters with their areas
        for (odb::dbMaster* master : masters) {
          const int area = master->getWidth() * master->getHeight();
          const bool is_repr = (master == repr_master);
          
          std::cout << "      " << (is_repr ? "[REPR] " : "       ") 
                    << master->getName() 
                    << " (Area: " << area << " DBU^2";
          
          if (!is_repr && area < repr_area) {
            const int percent_smaller = (repr_area - area) * 100 / repr_area;
            std::cout << ", " << percent_smaller << "% smaller";
          }
          
          std::cout << ")" << std::endl;
        }
      }

      group_count++;
    }

    if (compatible_masters_.size() > static_cast<size_t>(num_samples_)) {
      std::cout << "\n  ... and " << (compatible_masters_.size() - num_samples_) 
                << " more groups" << std::endl;
    }

    std::cout << "[readCompatibleMasters] Completed successfully.\n" << std::endl;
  }
}

//==============================================================================
// Phase 3: readFlopUnits()
//==============================================================================

void 
AggloCluster::readFlopUnits()
{
  if (verbose_) {
    std::cout << "\n[readFlopUnits] Starting flop unit extraction..." << std::endl;
  }

  flop_units_.clear();
  inst_to_flop_id_.clear();

  // Step 1: Extract valid flop instances and create flop units
  int flop_id = 0;
  int total_instances = 0;
  
  for (odb::dbInst* inst : block_->getInsts()) {
    total_instances++;
    
    if (!isValidFlop(inst)) {
      continue;
    }

    const MasterMask master_mask = createMasterMask(inst);
    const InstMask inst_mask = createInstMask(inst);
    const odb::Rect bbox = inst->getBBox()->getBox();

    flop_units_.emplace_back(flop_id, inst, bbox, master_mask, inst_mask);
    inst_to_flop_id_[inst] = flop_id++;
  }

  if (verbose_) {
    // Step 2: Report extraction statistics
    std::cout << "[Step 1] Flop Unit Extraction Completed" << std::endl;
    std::cout << "  - Total instances checked: " << total_instances << std::endl;
    std::cout << "  - Valid flop instances found: " << flop_units_.size() << std::endl;

    if (flop_units_.empty()) {
      std::cout << "[readFlopUnits] No valid flop units found.\n" << std::endl;
      return;
    }

    // Step 4: Show sample flop units with detailed information
    std::cout << "\n[Step 2] Sample Flop Units" << std::endl;
    const size_t sample_count = std::min(num_samples_, static_cast<int>(flop_units_.size()));
    
    for (size_t i = 0; i < sample_count; ++i) {
      const FlopUnit& flop = flop_units_[i];
      std::cout << "  Flop #" << i << ": " << flop.inst_->getName() << std::endl;
      std::cout << "    Master: " << flop.inst_->getMaster()->getName() << std::endl;
      std::cout << "    Master Mask: " << flop.master_mask_.to_string() << std::endl;
      std::cout << "    Inst Mask: " << flop.inst_mask_.to_string() << std::endl;
      std::cout << "    Location: (" << flop.orig_pt_.x << ", " 
                << flop.orig_pt_.y << ")" << std::endl;
    }
    
    if (flop_units_.size() > sample_count) {
      std::cout << "  ... and " << (flop_units_.size() - sample_count) 
                << " more flop units" << std::endl;
    }

    std::cout << "[readFlopUnits] Completed successfully.\n" << std::endl;
  }
}

//==============================================================================
// Phase 4: createCompatibleGroups()
//==============================================================================

void 
AggloCluster::createCompatibleGroups()
{
  if (verbose_) {
    std::cout << "\n[createCompatibleGroups] Starting compatible group creation..." << std::endl;
  }

  // Step 1: Group flops by identical mask pairs
  compatible_groups_.clear();

  for (size_t i = 0; i < flop_units_.size(); ++i) {
    const FlopUnit& unit = flop_units_[i];
    compatible_groups_[{unit.master_mask_, unit.inst_mask_}].push_back(i);
  }

  if (verbose_) {
    std::cout << "[Step 1] Group Creation Completed" << std::endl;
    std::cout << "  - Total compatible groups created: " << compatible_groups_.size() << std::endl;

    // Step 2: Analyze group size distribution
    std::cout << "\n[Step 2] Group Size Distribution" << std::endl;
    
    std::map<size_t, int> size_dist;
    for (const auto& [masks, flop_ids] : compatible_groups_) {
      size_dist[flop_ids.size()]++;
    }
    
    for (const auto& [size, count] : size_dist) {
      std::cout << "  - Groups with " << size << " flop(s): " << count << std::endl;
    }

    // Step 3: Show sample groups (sorted by size for better visibility)
    if (!compatible_groups_.empty()) {
      std::cout << "\n[Step 3] Sample Compatible Groups" << std::endl;

      // Sort groups by size (largest first)
      std::vector<std::pair<std::pair<MasterMask, InstMask>, std::vector<int>>> sorted_groups;
      sorted_groups.reserve(compatible_groups_.size());
      for (const auto& group : compatible_groups_) {
        sorted_groups.push_back(group);
      }
      
      std::sort(sorted_groups.begin(), sorted_groups.end(),
                [](const auto& a, const auto& b) {
                  return a.second.size() > b.second.size();
                });

      // Display top groups
      const size_t num_groups_to_show = std::min(static_cast<size_t>(num_samples_), 
                                                   sorted_groups.size());
      
      for (size_t group_idx = 0; group_idx < num_groups_to_show; ++group_idx) {
        const auto& [masks, flop_ids] = sorted_groups[group_idx];
        const auto& [master_mask, inst_mask] = masks;
        
        std::cout << "\n  Group #" << group_idx << ": " << flop_ids.size() << " flops" << std::endl;
        std::cout << "    Master Mask: " << master_mask.to_string() << std::endl;
        std::cout << "    Inst Mask:   " << inst_mask.to_string() << std::endl;
        
        // Show sample flop instances (first 3)
        std::cout << "    Sample flops (first 3): ";
        const size_t num_flops_to_show = std::min(size_t(3), flop_ids.size());
        
        for (size_t i = 0; i < num_flops_to_show; ++i) {
          if (i > 0) std::cout << ", ";
          std::cout << flop_units_[flop_ids[i]].inst_->getName();
        }
        
        if (flop_ids.size() > 3) {
          std::cout << " ... (+" << (flop_ids.size() - 3) << " more)";
        }
        std::cout << std::endl;
      }

      if (sorted_groups.size() > num_samples_) {
        std::cout << "\n  ... and " << (sorted_groups.size() - num_samples_) 
                  << " more groups" << std::endl;
      }
    }

    std::cout << "[createCompatibleGroups] Completed successfully.\n" << std::endl;
  }
}

//==============================================================================
// Phase 5: readTimingPaths()
//==============================================================================

void 
AggloCluster::readTimingPaths()
{
  if (verbose_) {
    std::cout << "\n[readTimingPaths] Starting timing path extraction..." << std::endl;
  }

  timing_paths_.clear();

  // Step 1: Initialize STA and find timing paths
  sta_->ensureGraph();
  sta_->ensureLevelized();
  sta_->searchPreamble();

  sta::PathEndSeq path_ends = sta_->search()->findPathEnds(
      /* from */ nullptr,
      /* thrus */ nullptr,
      /* to */ nullptr,
      /* unconstrained */ false,
      /* corner */ sta_->cmdCorner(),
      /* min_max */ sta::MinMaxAll::max(),
      /* group_path_count */ INT_MAX,
      /* endpoint_path_count */ num_paths_per_endpoint_,
      /* unique_pins */ false,
      /* slack_min */ -sta::INF,
      /* slack_max */ sta::INF,
      /* sort_by_slack */ true,
      /* groups */ nullptr,
      /* setup */ true,
      /* hold */ false,
      /* recovery */ false,
      /* removal */ false,
      /* clk_gating_setup */ false,
      /* clk_gating_hold */ false);

  if (verbose_) {
    std::cout << "[Step 1] STA Path Search Completed" << std::endl;
    std::cout << "  - Total path ends found: " << path_ends.size() << std::endl;
    std::cout << "  - Paths per endpoint: " << num_paths_per_endpoint_ << std::endl;
  }

  // Step 2: Filter and extract flop-to-flop timing paths
  int total_paths_analyzed = 0;
  int valid_flop_paths = 0;
  int skipped_short_paths = 0;
  int skipped_non_flop_paths = 0;

  for (sta::PathEnd* path_end : path_ends) {
    total_paths_analyzed++;
    
    const sta::Slack slack = path_end->slack(sta_);
    sta::Path* path = path_end->path();
    sta::PathExpanded expanded(path, sta_);

    // Skip paths with insufficient vertices
    if (expanded.size() < 2) {
      skipped_short_paths++;
      continue;
    }

    // Extract start flop instance (launch flop)
    const sta::Path* start_path = expanded.path(expanded.startIndex());
    sta::Vertex* sta_start_vertex = start_path ? start_path->vertex(sta_) : nullptr;
    sta::Pin* sta_start_pin = sta_start_vertex ? sta_start_vertex->pin() : nullptr;
    sta::Instance* sta_start_inst = sta_start_pin ? network_->instance(sta_start_pin) : nullptr;
    odb::dbInst* db_start_inst = sta_start_inst ? network_->staToDb(sta_start_inst) : nullptr;

    // Extract end flop instance (capture flop)
    sta::Vertex* sta_end_vertex = path_end->vertex(sta_);
    sta::Pin* sta_end_pin = sta_end_vertex ? sta_end_vertex->pin() : nullptr;
    sta::Instance* sta_end_inst = sta_end_pin ? network_->instance(sta_end_pin) : nullptr;
    odb::dbInst* db_end_inst = sta_end_inst ? network_->staToDb(sta_end_inst) : nullptr;

    // Validate both endpoints are registered flops
    const bool start_is_valid = db_start_inst 
        && inst_to_flop_id_.find(db_start_inst) != inst_to_flop_id_.end();
    const bool end_is_valid = db_end_inst 
        && inst_to_flop_id_.find(db_end_inst) != inst_to_flop_id_.end();

    if (!start_is_valid || !end_is_valid) {
      skipped_non_flop_paths++;
      continue;
    }

    const int start_flop_id = inst_to_flop_id_[db_start_inst];
    const int end_flop_id = inst_to_flop_id_[db_end_inst];

    // Extract pin information (Q/QN for start, D for end)
    odb::dbITerm* start_pin = sta_start_pin ? network_->flatPin(sta_start_pin) : nullptr;
    odb::dbITerm* end_pin = sta_end_pin ? network_->flatPin(sta_end_pin) : nullptr;

    timing_paths_.emplace_back(path, slack, start_flop_id, end_flop_id, start_pin, end_pin);
    valid_flop_paths++;
  }

  if (verbose_) {
    std::cout << "[Step 2] Path Filtering Completed" << std::endl;
    std::cout << "  - Total paths analyzed: " << total_paths_analyzed << std::endl;
    std::cout << "  - Valid flop-to-flop paths: " << valid_flop_paths << std::endl;
    std::cout << "  - Skipped (short paths): " << skipped_short_paths << std::endl;
    std::cout << "  - Skipped (non-flop endpoints): " << skipped_non_flop_paths << std::endl;

    // Step 3: Analyze slack distribution
    if (valid_flop_paths > 0) {
      // Collect and sort slack values
      std::vector<sta::Slack> slacks;
      slacks.reserve(timing_paths_.size());
      for (const auto& tp : timing_paths_) {
        slacks.push_back(tp.slack_);
      }
      std::sort(slacks.begin(), slacks.end());

      // Calculate statistics
      const sta::Slack min_slack = slacks.front();
      const sta::Slack max_slack = slacks.back();
      const sta::Slack median_slack = slacks[slacks.size() / 2];
      const sta::Slack avg_slack = std::accumulate(slacks.begin(), slacks.end(), 0.0) 
                                    / slacks.size();
      const int negative_slack_count = std::count_if(slacks.begin(), slacks.end(), 
                                                      [](sta::Slack s) { return s < 0; });

      // Get time units for display
      sta::Unit* time_unit = sta_->units()->timeUnit();
      const char* time_suffix = time_unit->scaledSuffix();

      std::cout << "\n[Step 3] Slack Distribution Analysis" << std::endl;
      std::cout << "  - Min slack: " << time_unit->asString(min_slack, 3) << time_suffix << std::endl;
      std::cout << "  - Max slack: " << time_unit->asString(max_slack, 3) << time_suffix << std::endl;
      std::cout << "  - Median slack: " << time_unit->asString(median_slack, 3) << time_suffix << std::endl;
      std::cout << "  - Average slack: " << time_unit->asString(avg_slack, 3) << time_suffix << std::endl;
      std::cout << "  - Negative slack paths: " << negative_slack_count 
                << " (" << std::fixed << std::setprecision(1) 
                << (100.0 * negative_slack_count / valid_flop_paths) << "%)" << std::endl;

      // Display sample critical paths
      std::cout << "\n[Sample Critical Paths (worst slack)]" << std::endl;
      const int num_paths_to_show = std::min(num_samples_, static_cast<int>(timing_paths_.size()));
      
      for (int i = 0; i < num_paths_to_show; ++i) {
        const auto& tp = timing_paths_[i];
        const FlopUnit& start_flop = flop_units_[tp.start_flop_idx_];
        const FlopUnit& end_flop = flop_units_[tp.end_flop_idx_];
        
        // Show path endpoints
        std::cout << "  Path #" << i << ": " << start_flop.inst_->getName();
        if (tp.start_pin_) {
          std::cout << "/" << tp.start_pin_->getMTerm()->getName();
        }
        std::cout << " -> " << end_flop.inst_->getName();
        if (tp.end_pin_) {
          std::cout << "/" << tp.end_pin_->getMTerm()->getName();
        }
        std::cout << std::endl;
        
        // Show slack value
        std::cout << "    Slack: " << time_unit->asString(slacks[i], 3) << time_suffix << std::endl;
        
        // Show detailed pin information
        if (tp.start_pin_ || tp.end_pin_) {
          std::cout << "    Pin details:" << std::endl;
          
          if (tp.start_pin_) {
            std::cout << "      Start: " << tp.start_pin_->getInst()->getName() 
                      << "/" << tp.start_pin_->getMTerm()->getName();
            if (tp.start_pin_->getNet()) {
              std::cout << " (net: " << tp.start_pin_->getNet()->getName() << ")";
            }
            std::cout << std::endl;
          }
          
          if (tp.end_pin_) {
            std::cout << "      End:   " << tp.end_pin_->getInst()->getName() 
                      << "/" << tp.end_pin_->getMTerm()->getName();
            if (tp.end_pin_->getNet()) {
              std::cout << " (net: " << tp.end_pin_->getNet()->getName() << ")";
            }
            std::cout << std::endl;
          }
        }
      }
      
      if (timing_paths_.size() > static_cast<size_t>(num_samples_)) {
        std::cout << "  ... and " << (timing_paths_.size() - num_samples_) 
                  << " more paths" << std::endl;
      }
    }

    std::cout << "[readTimingPaths] Completed successfully.\n" << std::endl;
  }
}

//==============================================================================
// Phase 6: analyzeTimingPaths()
//==============================================================================

void 
AggloCluster::analyzeTimingPaths()
{
  if (verbose_) {
    std::cout << "\n[analyzeTimingPaths] Starting timing path analysis..." << std::endl;
  }

  // Step 0: Clear existing budgets
  for (FlopUnit& flop : flop_units_) {
    flop.pin_budgets_.clear();
  }

  // Step 1: Distribute slack budgets to flop pins
  int start_pin_budget_count = 0;
  int end_pin_budget_count = 0;

  for (size_t i = 0; i < timing_paths_.size(); ++i) {
    const auto& path = timing_paths_[i];

    // Calculate slack budget (split equally, including negative slack)
    const sta::Slack budget = path.slack_ / 2.0f;

    // Start Flop: Store budget for its output pin (Q or QN)
    if (path.start_flop_idx_ >= 0 && path.start_pin_) {
      FlopUnit& start_flop = flop_units_[path.start_flop_idx_];
      start_flop.pin_budgets_[path.start_pin_].emplace_back(i, budget);
      start_pin_budget_count++;
    }

    // End Flop: Store budget for its input pin (D)
    if (path.end_flop_idx_ >= 0 && path.end_pin_) {
      FlopUnit& end_flop = flop_units_[path.end_flop_idx_];
      end_flop.pin_budgets_[path.end_pin_].emplace_back(i, budget);
      end_pin_budget_count++;
    }
  }

  if (verbose_) {
    std::cout << "[Step 1] Slack Budget Distribution Completed" << std::endl;
    std::cout << "  - Total timing paths processed: " << timing_paths_.size() << std::endl;
    std::cout << "  - Start pin budgets assigned: " << start_pin_budget_count << std::endl;
    std::cout << "  - End pin budgets assigned: " << end_pin_budget_count << std::endl;
  }

  // Step 2: Sort pin budgets in ascending order (most critical first)
  int total_pins_sorted = 0;
  
  for (FlopUnit& flop : flop_units_) {
    for (auto& [pin, budget_list] : flop.pin_budgets_) {
      std::sort(budget_list.begin(),
                budget_list.end(),
                [](const auto& a, const auto& b) { 
                  return a.second < b.second;  // Sort by slack budget (ascending)
                });
      total_pins_sorted++;
    }
  }

  if (verbose_) {
    std::cout << "\n[Step 2] Pin Budget Sorting Completed" << std::endl;
    std::cout << "  - Total pins with budgets sorted: " << total_pins_sorted << std::endl;
  }

  //----------------------------------------------------------------------------
  // Step 3: Analyze constraint distribution
  //----------------------------------------------------------------------------
  if (verbose_) {
    int flops_with_constraints = 0;
    int flops_without_constraints = 0;
    int total_pin_constraints = 0;
    int d_pin_constraints = 0;
    int q_pin_constraints = 0;
    int qn_pin_constraints = 0;
    
    for (const auto& flop : flop_units_) {
      if (!flop.pin_budgets_.empty()) {
        flops_with_constraints++;
        for (const auto& [pin, budget_list] : flop.pin_budgets_) {
          total_pin_constraints++;
          if (isDPin(pin)) {
            d_pin_constraints++;
          } else if (isQPin(pin)) {
            q_pin_constraints++;
          } else if (isQNPin(pin)) {
            qn_pin_constraints++;
          }
        }
      } else {
        flops_without_constraints++;
      }
    }
    
    std::cout << "\n[Step 3] Timing Constraint Statistics" << std::endl;
    std::cout << "  - Flops with timing constraints: " << flops_with_constraints 
              << " / " << flop_units_.size() 
              << " (" << std::fixed << std::setprecision(1) 
              << (100.0 * flops_with_constraints / flop_units_.size()) << "%)" << std::endl;
    std::cout << "  - Flops without constraints: " << flops_without_constraints << std::endl;
    std::cout << "  - Total pin-level constraints: " << total_pin_constraints << std::endl;
    std::cout << "  - Constraint breakdown:" << std::endl;
    std::cout << "    D pins (input):  " << d_pin_constraints << std::endl;
    std::cout << "    Q pins (output): " << q_pin_constraints << std::endl;
    std::cout << "    QN pins (output): " << qn_pin_constraints << std::endl;

    //--------------------------------------------------------------------------
    // Step 4: Sample flop pin budget details
    //--------------------------------------------------------------------------
    std::cout << "\n[Step 4] Sample Flop Pin Budgets" << std::endl;
    
    sta::Unit* time_unit = sta_->units()->timeUnit();
    const char* time_suffix = time_unit->scaledSuffix();
    
    int sample_count = 0;
    for (size_t i = 0; i < flop_units_.size(); ++i) {
      const FlopUnit& flop = flop_units_[i];
      
      if (sample_count >= num_samples_) {
        break;
      }
      
      if (flop.pin_budgets_.empty()) {
        continue;
      }
      
      std::cout << "\n  Flop #" << sample_count << ": " << flop.inst_->getName() << std::endl;
      std::cout << "    Master: " << flop.inst_->getMaster()->getName() << std::endl;
      std::cout << "    Total pins with constraints: " << flop.pin_budgets_.size() << std::endl;
      
      for (const auto& [pin, budget_list] : flop.pin_budgets_) {
        std::string pin_type;
        if (isDPin(pin)) {
          pin_type = "D (input)";
        } else if (isQPin(pin)) {
          pin_type = "Q (output)";
        } else if (isQNPin(pin)) {
          pin_type = "QN (output)";
        } else {
          pin_type = "Unknown";
        }
        
        const sta::Slack min_budget = budget_list.front().second;  // Already sorted ascending
        const sta::Slack max_budget = budget_list.back().second;
        
        std::cout << "      Pin: " << pin->getMTerm()->getName() 
                  << " [" << pin_type << "]" << std::endl;
        std::cout << "        Associated paths: " << budget_list.size() << std::endl;
        std::cout << "        Min slack budget (most critical): " 
                  << time_unit->asString(min_budget, 3) << time_suffix << std::endl;
        std::cout << "        Max slack budget: " 
                  << time_unit->asString(max_budget, 3) << time_suffix << std::endl;
      }
      
      sample_count++;
    }
    
    if (flops_with_constraints > num_samples_) {
      std::cout << "\n  ... and " << (flops_with_constraints - num_samples_) 
                << " more flops with constraints" << std::endl;
    }
    
    std::cout << "[analyzeTimingPaths] Completed successfully.\n" << std::endl;
  }
}

//==============================================================================
// Phase 7: calcFeasibleRegions()
//==============================================================================

void 
AggloCluster::calcFeasibleRegions()
{
  if (verbose_) {
    std::cout << "\n[calcFeasibleRegions] Starting feasible region calculation..." << std::endl;
  }

  const int total_flops = flop_units_.size();

  if (threads_ > 1) {
    // Multi-threaded mode: Parallel processing for performance
    if (verbose_) {
      std::cout << "[Step 1] Calculating feasible regions (Parallel)" << std::endl;
      std::cout << "  - Mode: Parallel" << std::endl;
      std::cout << "  - Threads: " << threads_ << std::endl;
      std::cout << "  - Total flops: " << total_flops << std::endl;
    }

    #pragma omp parallel for schedule(static) num_threads(threads_)
    for (int i = 0; i < total_flops; ++i) {
      calcFeasibleRegion(flop_units_[i], false);
    }

    if (verbose_) {
      std::cout << "[calcFeasibleRegions] Completed successfully." << std::endl;
      std::cout << "  - Total flops processed: " << total_flops << std::endl;
      
      // Output feasible regions for first num_samples_ flops
      std::cout << "\n[Verification] Feasible regions for first " << std::min(num_samples_, total_flops) << " flops:" << std::endl;
      for (int i = 0; i < std::min(num_samples_, total_flops); ++i) {
        const FlopUnit& flop = flop_units_[i];
        const Box& fr = flop.feasible_region_;
        
        std::cout << "  Flop #" << i << " (" << flop.inst_->getName() << "):" << std::endl;
        
        if (boost::geometry::is_empty(fr)) {
          std::cout << "    Feasible region: EMPTY" << std::endl;
        } else {
          const Point min_corner = fr.min_corner();
          const Point max_corner = fr.max_corner();
          std::cout << "    Feasible region (UV): [(" 
                    << min_corner.get<0>() << ", " << min_corner.get<1>() << ") -> ("
                    << max_corner.get<0>() << ", " << max_corner.get<1>() << ")]" << std::endl;
          
          // Calculate area
          const int width = max_corner.get<0>() - min_corner.get<0>();
          const int height = max_corner.get<1>() - min_corner.get<1>();
          std::cout << "    Size: " << width << " x " << height << " = " << (width * height) << " DBU²" << std::endl;
        }
      }
    }
  } else {
    // Single-threaded mode: Sequential processing with optional debug output
    if (verbose_) {
      std::cout << "[Step 1] Calculating feasible regions (Sequential)" << std::endl;
      std::cout << "  - Mode: Sequential" << std::endl;
      std::cout << "  - Total flops: " << total_flops << std::endl;
      std::cout << "  - Sample size: " << num_samples_ << std::endl;
    }

    for (int i = 0; i < total_flops; ++i) {
      FlopUnit& flop = flop_units_[i];
      const bool enable_verbose = verbose_ && (i < num_samples_);

      if (enable_verbose) {
        std::cout << "\n[Sample Flop #" << i << "] " << flop.inst_->getName() << std::endl;
        std::cout << "  Master: " << flop.inst_->getMaster()->getName() << std::endl;
        std::cout << "  Location: (" << flop.orig_pt_.x << ", " << flop.orig_pt_.y << ")" << std::endl;
      }

      calcFeasibleRegion(flop, enable_verbose);
    }
  }
}

//==============================================================================
// Phase 8: createFlopClusters()
//==============================================================================

void 
AggloCluster::createFlopClusters()
{
  if (verbose_) {
    std::cout << "\n[createFlopClusters] Starting flop cluster initialization..." << std::endl;
  }

  // Step 0: Initialize data structures
  flop_clusters_.clear();
  flop_clusters_.reserve(flop_units_.size());
  feasible_regions_.clear();
  flop_cluster_is_valid_.assign(flop_units_.size(), true);
  flop_cluster_no_further_merge_.assign(flop_units_.size(), false);

  if (verbose_) {
    std::cout << "[Step 1] Data Structure Initialization Completed" << std::endl;
    std::cout << "  - Cluster capacity reserved: " << flop_units_.size() << std::endl;
    std::cout << "  - Validation flags initialized: " << flop_units_.size() << std::endl;
  }

  // Step 1: Create initial single-flop clusters
  for (size_t i = 0; i < flop_units_.size(); ++i) {
    flop_clusters_.emplace_back(i, flop_units_[i]);
    flop_units_[i].cluster_idx_ = i;

    const Box& feasible_region = flop_units_[i].feasible_region_;
    feasible_regions_.insert(std::make_pair(feasible_region, i));
  }

  if (verbose_) {
    std::cout << "\n[Step 2] Initial Cluster Creation Completed" << std::endl;
    std::cout << "  - Total flop clusters created: " << flop_clusters_.size() << std::endl;
    std::cout << "  - Feasible regions registered: " << feasible_regions_.size() << std::endl;

    // Step 2: Calculate overlap statistics
    std::vector<int> overlap_counts;
    overlap_counts.reserve(flop_clusters_.size());
    
    for (size_t idx = 0; idx < flop_clusters_.size(); ++idx) {
      const std::vector<FlopClusterEntry> intersecting_clusters = getIntersectedCluster(flop_clusters_[idx]);
      
      int count = 0;
      for (const auto& [other_box, other_idx] : intersecting_clusters) {
        if (static_cast<size_t>(other_idx) != idx) {
          count++;
        }
      }
      overlap_counts.push_back(count);
    }

    if (!overlap_counts.empty()) {
      std::sort(overlap_counts.begin(), overlap_counts.end());
      
      int min_overlap = overlap_counts.front();
      int max_overlap = overlap_counts.back();
      int median_overlap = overlap_counts[overlap_counts.size() / 2];
      double avg_overlap = std::accumulate(overlap_counts.begin(), overlap_counts.end(), 0.0) / overlap_counts.size();
      
      std::cout << "\n[Step 3] Feasible Region Overlap Analysis" << std::endl;
      std::cout << "  - Min overlaps per cluster: " << min_overlap << std::endl;
      std::cout << "  - Max overlaps per cluster: " << max_overlap << std::endl;
      std::cout << "  - Median overlaps per cluster: " << median_overlap << std::endl;
      std::cout << "  - Average overlaps per cluster: " << std::fixed << std::setprecision(1) << avg_overlap << std::endl;
    }

    // Step 3: Display sample clusters
    if (!flop_clusters_.empty()) {
      std::cout << "\n[Step 4] Sample Cluster Details" << std::endl;
      const size_t num_samples = std::min(num_samples_, static_cast<int>(flop_clusters_.size()));
      
      for (size_t i = 0; i < num_samples; ++i) {
        const FlopCluster& cluster = flop_clusters_[i];
        const FlopUnit& flop = flop_units_[i];
        
        std::cout << "\n  Cluster #" << i << ":" << std::endl;
        std::cout << "    Flop instance: " << flop.inst_->getName() << std::endl;
        std::cout << "    Master: " << flop.inst_->getMaster()->getName() << std::endl;
        std::cout << "    Original location: (" << flop.orig_pt_.x << ", " << flop.orig_pt_.y << ")" << std::endl;
        
        // Feasible region details (45-degree rotated coordinate system)
        const Box& fr = flop.feasible_region_;
        const bool is_empty = boost::geometry::is_empty(fr);
        
        if (is_empty) {
          std::cout << "    Feasible region: EMPTY (no valid placement)" << std::endl;
        } else {
          std::cout << "    Feasible region (45° rotated): (" 
                    << fr.min_corner().get<0>() << ", " << fr.min_corner().get<1>() << ") -> ("
                    << fr.max_corner().get<0>() << ", " << fr.max_corner().get<1>() << ")" << std::endl;
          
          int fr_width = fr.max_corner().get<0>() - fr.min_corner().get<0>();
          int fr_height = fr.max_corner().get<1>() - fr.min_corner().get<1>();
          std::cout << "    FR dimensions: " << fr_width << " x " << fr_height << " (rotated coords)" << std::endl;
        }
        
        // Overlapping clusters analysis
        const int overlap_count = overlap_counts[i];
        std::cout << "    Overlapping clusters: " << overlap_count << std::endl;
        
        if (overlap_count > 0 && overlap_count <= 3) {
          const std::vector<FlopClusterEntry> intersecting = getIntersectedCluster(cluster);
          std::cout << "      Examples: ";
          int shown = 0;
          for (const auto& [other_box, other_idx] : intersecting) {
            if (static_cast<int>(other_idx) == static_cast<int>(i)) continue;
            if (shown > 0) std::cout << ", ";
            std::cout << "Cluster#" << other_idx;
            shown++;
            if (shown >= 3) break;
          }
          std::cout << std::endl;
        } else if (overlap_count > 3) {
          const std::vector<FlopClusterEntry> intersecting = getIntersectedCluster(cluster);
          std::cout << "      Examples: ";
          int shown = 0;
          for (const auto& [other_box, other_idx] : intersecting) {
            if (static_cast<int>(other_idx) == static_cast<int>(i)) continue;
            if (shown > 0) std::cout << ", ";
            std::cout << "Cluster#" << other_idx;
            shown++;
            if (shown >= 3) break;
          }
          std::cout << " ... (+" << (overlap_count - 3) << " more)" << std::endl;
        }
        
        // Cluster state
        std::cout << "    State: " 
                  << (flop_cluster_is_valid_[i] ? "Valid" : "Invalid") 
                  << ", "
                  << (flop_cluster_no_further_merge_[i] ? "Final" : "Mergeable") 
                  << std::endl;
      }
      
      if (flop_clusters_.size() > num_samples) {
        std::cout << "\n  ... and " << (flop_clusters_.size() - num_samples) 
                  << " more clusters" << std::endl;
      }
    }

    std::cout << "[createFlopClusters] Completed successfully.\n" << std::endl;
  }
}

//==============================================================================
// Phase 9: createCompatibilityGraph()
//==============================================================================

void
AggloCluster::createCompatibilityGraph()
{
  if (verbose_) {
    std::cout << "\n[createCompatibilityGraph] Starting compatibility graph construction..." << std::endl;
  }

  // Step 0: Initialize data structures
  edge_pq_.clear();
  adj_list_.clear();

  if (verbose_) {
    std::cout << "[Step 1] Data Structure Initialization" << std::endl;
    std::cout << "  - Edge priority queue cleared" << std::endl;
    std::cout << "  - Adjacency list cleared" << std::endl;
  }

  // Step 1: Build edges for all clusters
  const int total_clusters = flop_clusters_.size();

  if (threads_ > 1) {
    // Multi-threaded mode: Parallel processing for performance
    if (verbose_) {
      std::cout << "\n[Step 2] Building Edges (Parallel)" << std::endl;
      std::cout << "  - Mode: Parallel" << std::endl;
      std::cout << "  - Threads: " << threads_ << std::endl;
      std::cout << "  - Total clusters: " << total_clusters << std::endl;
    }

    // Collect edges from all clusters in parallel
    std::vector<std::vector<Edge>> all_edges(flop_clusters_.size());

    #pragma omp parallel for schedule(static) num_threads(threads_)
    for (size_t i = 0; i < flop_clusters_.size(); ++i) {
      all_edges[i] = updateEdges(flop_clusters_[i], false);
    }

    // Merge all edges into global data structures (sequential merge)
    for (const auto& edges : all_edges) {
      for (const Edge& e : edges) {
        edge_pq_.insert(e);
        adj_list_[e.n1].insert(e);
        adj_list_[e.n2].insert(e);
      }
    }

    if (verbose_) {
      std::cout << "\n[Step 3] Edge Construction Completed" << std::endl;
      std::cout << "  - Total edges created: " << edge_pq_.size() << std::endl;
      std::cout << "  - Adjacency list entries: " << adj_list_.size() << std::endl;

      // Analyze edge weight distribution
      if (!edge_pq_.empty()) {
        std::vector<double> edge_weights;
        edge_weights.reserve(edge_pq_.size());
        for (const auto& edge : edge_pq_) {
          edge_weights.push_back(edge.weight);
        }
        std::sort(edge_weights.begin(), edge_weights.end());

        const double min_weight = edge_weights.front();
        const double max_weight = edge_weights.back();
        const double median_weight = edge_weights[edge_weights.size() / 2];
        const double avg_weight = std::accumulate(edge_weights.begin(), edge_weights.end(), 0.0) / edge_weights.size();
        const int negative_count = std::count_if(edge_weights.begin(), edge_weights.end(),
                                                  [](double w) { return w < 0; });
        const float negative_pct = 100.0f * negative_count / edge_weights.size();

        std::cout << "\n[Step 4] Edge Weight Distribution" << std::endl;
        std::cout << "  - Min weight (best): " << std::fixed << std::setprecision(2) << min_weight << std::endl;
        std::cout << "  - Max weight: " << max_weight << std::endl;
        std::cout << "  - Median weight: " << median_weight << std::endl;
        std::cout << "  - Average weight: " << avg_weight << std::endl;
        std::cout << "  - HPWL improving edges: " << negative_count
                  << " (" << std::setprecision(1) << negative_pct << "%)" << std::endl;
      }

      std::cout << "[createCompatibilityGraph] Completed successfully.\n" << std::endl;
    }
  } else {
    // Single-threaded mode: Sequential processing with optional debug output
    if (verbose_) {
      std::cout << "\n[Step 2] Building Edges (Sequential)" << std::endl;
      std::cout << "  - Mode: Sequential" << std::endl;
      std::cout << "  - Total clusters: " << total_clusters << std::endl;
      std::cout << "  - Sample size: " << num_samples_ << std::endl;
    }

    for (size_t i = 0; i < flop_clusters_.size(); ++i) {
      const bool enable_verbose = verbose_ && (i < static_cast<size_t>(num_samples_));
      std::vector<Edge> edges = updateEdges(flop_clusters_[i], enable_verbose);
      
      // Insert edges into global data structures
      for (const Edge& e : edges) {
        edge_pq_.insert(e);
        adj_list_[e.n1].insert(e);
        adj_list_[e.n2].insert(e);
      }
    }

    if (verbose_) {
      std::cout << "\n[Step 3] Edge Construction Completed" << std::endl;
      std::cout << "  - Total edges created: " << edge_pq_.size() << std::endl;
      std::cout << "  - Adjacency list entries: " << adj_list_.size() << std::endl;
      
      if (total_clusters > num_samples_) {
        std::cout << "  - Detailed output shown for first " << num_samples_ << " clusters" << std::endl;
      }

      std::cout << "[createCompatibilityGraph] Completed successfully.\n" << std::endl;
    }
  }
}

//==============================================================================
// Phase 10: runAgglomerativeClustering()
//==============================================================================

void 
AggloCluster::runAgglomerativeClustering()
{
  if (verbose_) {
    std::cout << "\n[runAgglomerativeClustering] Starting agglomerative clustering..." << std::endl;
    std::cout << "  - Initial edge queue size: " << edge_pq_.size() << std::endl;
    std::cout << "  - Initial cluster count: " << flop_clusters_.size() << std::endl;
  }

  // Step 0: Initialize counters
  int iteration = 0;
  int total_merges = 0;
  int intermediate_merges = 0;
  int final_merges = 0;
  int skipped_invalid = 0;
  int skipped_density = 0;
  std::map<int, int> merge_size_dist;

  // Step 1: Main clustering loop
  while (!edge_pq_.empty()) {
    const Edge best_edge = *edge_pq_.begin();
    edge_pq_.erase(edge_pq_.begin());

    const int n1_idx = best_edge.n1;
    const int n2_idx = best_edge.n2;

    if (!flop_cluster_is_valid_[n1_idx] || !flop_cluster_is_valid_[n2_idx]) {
      skipped_invalid++;
      continue;
    }

    const FlopCluster& c1 = flop_clusters_[n1_idx];
    const FlopCluster& c2 = flop_clusters_[n2_idx];
    const bool show_debug = verbose_ && (total_merges < num_samples_);

    // Step 1.1: Display merge information
    if (show_debug) {
      auto print_feasible_region = [](const Box& fr) {
        if (boost::geometry::is_empty(fr)) {
          std::cout << "EMPTY";
        } else {
          std::cout << "[(" << fr.min_corner().get<0>() << "," << fr.min_corner().get<1>() 
                    << ") -> (" << fr.max_corner().get<0>() << "," << fr.max_corner().get<1>() << ")]";
        }
      };

      std::cout << "\n[Merge #" << total_merges << " | Iteration " << iteration << "]" << std::endl;
      std::cout << "  ✓ Selected edge: Cluster[" << n1_idx << "] + Cluster[" << n2_idx << "]" << std::endl;
      std::cout << "    C1: size=" << c1.flops_.size() << ", pos=(" 
                << c1.curr_pt_.x << ", " << c1.curr_pt_.y << ")" << std::endl;
      std::cout << "    C2: size=" << c2.flops_.size() << ", pos=(" 
                << c2.curr_pt_.x << ", " << c2.curr_pt_.y << ")" << std::endl;
      std::cout << "    Edge weight: " << std::fixed << std::setprecision(2) 
                << best_edge.weight << (best_edge.weight < 0 ? " ✓ improvement" : " ✗ degradation") << std::endl;
      std::cout << "    New position: (" << best_edge.pos.x << ", " << best_edge.pos.y << ")" << std::endl;
      
      std::cout << "    C1 feasible region: ";
      print_feasible_region(c1.feasible_region_);
      std::cout << std::endl;
      
      std::cout << "    C2 feasible region: ";
      print_feasible_region(c2.feasible_region_);
      std::cout << std::endl;
    }

    // Step 1.2: Verify density feasibility with current occupancy
    bool density_ok = true;
    odb::dbMaster* master1 = getClusterMaster(c1);
    odb::dbMaster* master2 = getClusterMaster(c2);
    odb::dbMaster* merged_master = nullptr;
    const int merged_bits = static_cast<int>(c1.flops_.size() + c2.flops_.size());
    const auto rep_it = representative_masters_.find(c1.master_mask_);
    if (rep_it != representative_masters_.end()) {
      const auto bit_it = rep_it->second.find(merged_bits);
      if (bit_it != rep_it->second.end()) {
        merged_master = bit_it->second;
      }
    }

    if (!master1 || !master2 || !merged_master) {
      density_ok = false;
    } else {
      const Point c1_pos(std::lround(c1.curr_pt_.x), std::lround(c1.curr_pt_.y));
      const Point c2_pos(std::lround(c2.curr_pt_.x), std::lround(c2.curr_pt_.y));
      const Point merged_pos(std::lround(best_edge.pos.x), std::lround(best_edge.pos.y));

      density_ok = !virtual_bin_grid_.wouldOverflow(master1, c1_pos,
                                                    master2, c2_pos,
                                                    merged_master, merged_pos);
    }

    if (!density_ok) {
      skipped_density++;

      if (show_debug) {
        std::cout << "  ✗ Merge skipped: density overflow" << std::endl;
      }

      removeEdges(c1, show_debug);
      removeEdges(c2, show_debug);

      for (int cluster_idx : {n1_idx, n2_idx}) {
        if (!flop_cluster_is_valid_[cluster_idx]
            || flop_cluster_no_further_merge_[cluster_idx]) {
          continue;
        }

        const std::vector<Edge> refreshed_edges = updateEdges(flop_clusters_[cluster_idx], show_debug);
        for (const Edge& e : refreshed_edges) {
          edge_pq_.insert(e);
          adj_list_[e.n1].insert(e);
          adj_list_[e.n2].insert(e);
        }
      }

      continue;
    }

    // Step 1.2: Execute merge
    const int new_cluster_idx = mergeClusters(best_edge, show_debug);
    const FlopCluster& new_cluster = flop_clusters_[new_cluster_idx];

    if (show_debug) {
      std::cout << "  ✓ Created Cluster[" << new_cluster_idx << "]: size=" << new_cluster.flops_.size() << std::endl;
      std::cout << "    Feasible region: ";
      if (boost::geometry::is_empty(new_cluster.feasible_region_)) {
        std::cout << "EMPTY ⚠ WARNING!" << std::endl;
      } else {
        const Box& fr = new_cluster.feasible_region_;
        std::cout << "[(" << fr.min_corner().get<0>() << "," << fr.min_corner().get<1>() 
                  << ") -> (" << fr.max_corner().get<0>() << "," << fr.max_corner().get<1>() << ")]" << std::endl;
      }
    }

    total_merges++;
    merge_size_dist[new_cluster.flops_.size()]++;

    // Step 1.3: Determine merge type and handle accordingly
    const bool can_merge_further = isFurtherMergeable(new_cluster, show_debug);

    if (can_merge_further) {
      intermediate_merges++;
      
      if (show_debug) {
        std::cout << "  ✓ Decision: INTERMEDIATE (can grow further)" << std::endl;
      }

      const std::vector<Edge> new_edges = updateEdges(new_cluster);
      
      for (const Edge& e : new_edges) {
        edge_pq_.insert(e);
        adj_list_[e.n1].insert(e);
        adj_list_[e.n2].insert(e);
      }
      
      if (show_debug) {
        std::cout << "    New edges created: " << new_edges.size() << std::endl;
        
        const int neighbor_count = adj_list_.count(new_cluster_idx) ? adj_list_[new_cluster_idx].size() : 0;
        if (neighbor_count > 0 && neighbor_count <= 5) {
          std::cout << "    Top neighbors (showing up to 3):" << std::endl;
          int shown = 0;
          for (const Edge& e : adj_list_[new_cluster_idx]) {
            if (shown >= 3) break;
            const int neighbor_idx = (e.n1 == new_cluster_idx) ? e.n2 : e.n1;
            std::cout << "      Cluster[" << neighbor_idx << "]: size=" 
                      << flop_clusters_[neighbor_idx].flops_.size() 
                      << ", weight=" << std::fixed << std::setprecision(2) << e.weight << std::endl;
            shown++;
          }
        }
      }
    } 
    else {
      final_merges++;
      flop_cluster_no_further_merge_[new_cluster_idx] = true;

      if (show_debug) {
        std::cout << "  ✓ Decision: FINAL (reached maximum bit-width)" << std::endl;
      }

      const std::set<int> affected_flop_units = distributeSlack(new_cluster, show_debug);
      
      // Calculate feasible regions for affected flop units
      std::set<int> affected_flop_clusters;
      for (int u_idx : affected_flop_units) {
        calcFeasibleRegion(flop_units_[u_idx]);
        
        const int c_idx = flop_units_[u_idx].cluster_idx_;
        if (flop_cluster_is_valid_[c_idx] && !flop_cluster_no_further_merge_[c_idx]) {
          affected_flop_clusters.insert(c_idx);
        }
      }

      if (show_debug) {
        std::cout << "    Affected flops: " << affected_flop_units.size() << std::endl;
        std::cout << "    Affected clusters: " << affected_flop_clusters.size() << std::endl;
      }

      // Update feasible regions for affected clusters
      for (int c_idx : affected_flop_clusters) {
        updateFeasibleRegion(flop_clusters_[c_idx], show_debug);
      }

      // Update edges for affected clusters
      std::vector<Edge> all_new_edges;
      for (int c_idx : affected_flop_clusters) {
        std::vector<Edge> edges = updateEdges(flop_clusters_[c_idx]);
        all_new_edges.insert(all_new_edges.end(), edges.begin(), edges.end());
      }
      
      for (const Edge& e : all_new_edges) {
        edge_pq_.insert(e);
        adj_list_[e.n1].insert(e);
        adj_list_[e.n2].insert(e);
      }

      if (show_debug) {
        std::cout << "    Total new edges from affected clusters: " << all_new_edges.size() << std::endl;
      }
    }

    iteration++;

    if (verbose_ && iteration % 100 == 0) {
      std::cout << "\n[Progress Update] Iteration " << iteration << std::endl;
      std::cout << "  - Total merges: " << total_merges 
                << " (intermediate: " << intermediate_merges 
                << ", final: " << final_merges << ")" << std::endl;
      std::cout << "  - Skipped (invalid): " << skipped_invalid << std::endl;
  std::cout << "  - Skipped (density): " << skipped_density << std::endl;
      std::cout << "  - Remaining edges: " << edge_pq_.size() << std::endl;
      std::cout << "  - Valid clusters: " << std::count(flop_cluster_is_valid_.begin(), 
                                                        flop_cluster_is_valid_.end(), true) << std::endl;
      
      if (!merge_size_dist.empty()) {
        std::cout << "  - Current merge distribution:" << std::endl;
        for (const auto& [size, count] : merge_size_dist) {
          std::cout << "    " << size << "-bit: " << count << " merge(s)" << std::endl;
        }
      }
    }
  }

  // Step 2: Final statistics
  if (verbose_) {
    std::cout << "\n[Step 2] Final Clustering Statistics" << std::endl;
    std::cout << "  Execution summary:" << std::endl;
    std::cout << "    Total iterations: " << iteration << std::endl;
    std::cout << "    Total merges: " << total_merges 
              << " (intermediate: " << intermediate_merges 
              << ", final: " << final_merges << ")" << std::endl;
  std::cout << "    Skipped (invalid): " << skipped_invalid << std::endl;
  std::cout << "    Skipped (density): " << skipped_density << std::endl;

    std::map<int, int> cluster_size_dist;
    std::map<int, std::vector<int>> size_to_clusters;
    int single_flop_clusters = 0;
    int merged_clusters = 0;
    int total_flops_in_mbff = 0;

    for (size_t i = 0; i < flop_clusters_.size(); ++i) {
      if (flop_cluster_is_valid_[i]) {
        const int size = flop_clusters_[i].flops_.size();
        cluster_size_dist[size]++;
        size_to_clusters[size].push_back(i);
        
        if (size == 1) {
          single_flop_clusters++;
        } else {
          merged_clusters++;
          total_flops_in_mbff += size;
        }
      }
    }

    const int total_valid = single_flop_clusters + merged_clusters;
    const int total_flops = flop_units_.size();
    const float mbff_ratio = total_flops > 0 ? (100.0f * total_flops_in_mbff / total_flops) : 0.0f;

    std::cout << "\n  Final cluster status:" << std::endl;
    std::cout << "    Total valid clusters: " << total_valid << std::endl;
    std::cout << "    Single-flop clusters: " << single_flop_clusters 
              << " (" << std::fixed << std::setprecision(1) 
              << (100.0f * single_flop_clusters / total_valid) << "%)" << std::endl;
    std::cout << "    Multi-bit clusters (MBFF): " << merged_clusters 
              << " (" << (100.0f * merged_clusters / total_valid) << "%)" << std::endl;
    std::cout << "    Flops in MBFF: " << total_flops_in_mbff << " / " << total_flops 
              << " (" << mbff_ratio << "%)" << std::endl;

    if (!cluster_size_dist.empty()) {
      std::cout << "\n  Cluster size distribution:" << std::endl;
      int max_size = 0;
      for (const auto& [size, count] : cluster_size_dist) {
        max_size = std::max(max_size, size);
        std::cout << "    " << size << "-bit: " << count << " cluster(s)";
        
        if (size > 1) {
          const float pct = 100.0f * count / merged_clusters;
          std::cout << " (" << std::fixed << std::setprecision(1) << pct << "% of MBFF)";
        }
        std::cout << std::endl;
      }
      
      std::cout << "    Maximum bit-width achieved: " << max_size << std::endl;
    }

    std::cout << "\n  Termination analysis:" << std::endl;
    if (edge_pq_.empty()) {
      std::cout << "    Reason: Edge queue exhausted (no more mergeable pairs)" << std::endl;
      
      int remaining_intermediate = 0;
      int remaining_final = 0;
      for (size_t i = 0; i < flop_clusters_.size(); ++i) {
        if (flop_cluster_is_valid_[i]) {
          if (flop_cluster_no_further_merge_[i]) {
            remaining_final++;
          } else {
            remaining_intermediate++;
          }
        }
      }
      
      std::cout << "    Remaining clusters that could grow: " << remaining_intermediate << std::endl;
      std::cout << "    Remaining final clusters: " << remaining_final << std::endl;
      
      if (remaining_intermediate > 0) {
        std::cout << "    ⚠ WARNING: " << remaining_intermediate 
                  << " cluster(s) marked as intermediate but no edges found!" << std::endl;
        std::cout << "      This suggests edge creation or compatibility issues." << std::endl;
      }
    }

    if (merged_clusters > 0) {
      std::cout << "\n  Sample multi-bit clusters:" << std::endl;
      
      std::vector<std::pair<int, int>> size_count_pairs;
      for (const auto& [size, count] : cluster_size_dist) {
        if (size > 1) {
          size_count_pairs.emplace_back(size, count);
        }
      }
      std::sort(size_count_pairs.begin(), size_count_pairs.end(), 
                [](const auto& a, const auto& b) { return a.first > b.first; });
      
      int samples_shown = 0;
      constexpr int MAX_SAMPLES = 5;
      
      for (const auto& [size, count] : size_count_pairs) {
        if (samples_shown >= MAX_SAMPLES) break;
        
        const auto& cluster_indices = size_to_clusters[size];
        for (size_t i = 0; i < std::min(size_t(2), cluster_indices.size()) && samples_shown < MAX_SAMPLES; ++i) {
          const int c_idx = cluster_indices[i];
          const FlopCluster& cluster = flop_clusters_[c_idx];
          
          std::cout << "    [" << (samples_shown + 1) << "] Cluster #" << c_idx 
                    << ": " << size << "-bit MBFF" << std::endl;
          std::cout << "        Position: (" << cluster.curr_pt_.x << ", " << cluster.curr_pt_.y << ")" << std::endl;
          
          std::cout << "        Feasible region: ";
          if (boost::geometry::is_empty(cluster.feasible_region_)) {
            std::cout << "EMPTY";
          } else {
            const Box& fr = cluster.feasible_region_;
            std::cout << "[(" << fr.min_corner().get<0>() << "," << fr.min_corner().get<1>() 
                      << ") -> (" << fr.max_corner().get<0>() << "," << fr.max_corner().get<1>() << ")]";
          }
          std::cout << std::endl;
          
          std::cout << "        Member flops: ";
          int member_count = 0;
          for (int flop_idx : cluster.flops_) {
            if (member_count > 0) std::cout << ", ";
            std::cout << flop_units_[flop_idx].inst_->getName();
            member_count++;
            if (member_count >= 3 && cluster.flops_.size() > 3) {
              std::cout << " ... (+" << (cluster.flops_.size() - 3) << " more)";
              break;
            }
          }
          std::cout << std::endl;
          
          samples_shown++;
        }
      }
    }

    std::cout << "[runAgglomerativeClustering] Completed successfully.\n" << std::endl;
  }
}

//==============================================================================
// Phase 11: implementClusters()
//==============================================================================

void 
AggloCluster::implementClusters()
{
  if (verbose_) {
    std::cout << "\n[implementClusters] Starting cluster implementation..." << std::endl;
  }

  // Identify final clusters to implement
  std::vector<int> final_cluster_indices;
  for (int i = 0; i < flop_clusters_.size(); ++i) {
    if (flop_cluster_is_valid_[i] && flop_clusters_[i].flops_.size() > 1) {
      final_cluster_indices.push_back(i);
    }
  }

  if (verbose_) {
    std::cout << "[Step 1] Identified " << final_cluster_indices.size() 
              << " clusters to implement" << std::endl;
  }

  // Implement each cluster
  int implemented = 0;
  int skipped = 0;

  for (int idx = 0; idx < final_cluster_indices.size(); ++idx) {
    const int cluster_idx = final_cluster_indices[idx];
    const bool show_debug = (idx < num_samples_) && verbose_;
    
    if (implementSingleCluster(flop_clusters_[cluster_idx], show_debug, idx)) {
      implemented++;
    } else {
      skipped++;
    }
  }

  if (verbose_) {
    if (final_cluster_indices.size() > num_samples_) {
      std::cout << "\n  ... and " << (final_cluster_indices.size() - num_samples_) 
                << " more clusters processed" << std::endl;
    }
    std::cout << "\n[Step 2] Implementation Summary" << std::endl;
    std::cout << "  Total: " << final_cluster_indices.size() << std::endl;
    std::cout << "  Implemented: " << implemented << std::endl;
    std::cout << "  Skipped: " << skipped << std::endl;
    std::cout << "[implementClusters] Completed successfully.\n" << std::endl;
  }
}

//==============================================================================
// Phase 2-3 Helper Functions: Mask Creation & Flop Validation
//==============================================================================

// -----------------------------------------------------------------------------
// [Level 1] Flop Validation
// -----------------------------------------------------------------------------

bool 
AggloCluster::isValidFlop(odb::dbInst* inst) const
{
  // Step 1: Check Liberty cell availability
  const sta::LibertyCell* lib_cell = getLibertyCell(inst);
  if (!lib_cell || lib_cell->sequentials().empty()) {
    return false;
  }

  // Step 2: Filter out invalid cell types
  if (lib_cell->isClockGate() || resizer_->dontUse(lib_cell)) {
    return false;
  }

  // Step 3: Filter out latches (only flip-flops allowed)
  const bool has_latch = std::any_of(
      lib_cell->sequentials().begin(),
      lib_cell->sequentials().end(),
      [](const sta::Sequential* seq) { return seq->isLatch(); });
  
  if (has_latch) {
    return false;
  }

  // Step 4: Validate pin configuration
  int d_pins = 0;
  int q_pins = 0;
  int qn_pins = 0;

  for (odb::dbITerm* iterm : inst->getITerms()) {
    if (isQPin(iterm)) {
      q_pins++;
    } else if (isQNPin(iterm)) {
      qn_pins++;
    } else if (isDPin(iterm)) {
      d_pins++;
    } else if (isClockPin(iterm) || isPowerPin(iterm) || isPresetPin(iterm)
               || isClearPin(iterm) || isScanInPin(iterm)
               || isScanEnablePin(iterm)) {
      // Known control/power pins - valid
      continue;
    } else {
      // Unrecognized pin type - invalid
      return false;
    }
  }

  // Step 5: Verify D/Q pin count matches (D == max(Q, QN))
  return d_pins == std::max(q_pins, qn_pins);
}

// -----------------------------------------------------------------------------
// [Level 1] Mask Creation (Top-level)
// -----------------------------------------------------------------------------

InstMask 
AggloCluster::createInstMask(odb::dbInst* inst)
{
  odb::dbNet* clock_net = nullptr;
  odb::dbNet* clear_net = nullptr;
  odb::dbNet* preset_net = nullptr;
  odb::dbNet* scan_enable_net = nullptr;
  odb::dbNet* scan_in_net = nullptr;

  // Helper lambda to assign net with conflict detection
  auto assignNet = [&](odb::dbNet*& target_net, odb::dbNet* new_net, const char* net_type) {
    if (!target_net) {
      target_net = new_net;
    } else if (target_net != new_net && verbose_) {
      std::cout << "[createInstMask] Warning: Instance " << inst->getName() 
                << " has multiple different " << net_type << " nets. "
                << "Using the first one (" << target_net->getName() << ")" << std::endl;
    }
  };

  for (odb::dbITerm* iterm : inst->getITerms()) {
    odb::dbNet* net = iterm->getNet();
    if (!net) {
      continue;
    }

    if (isClockPin(iterm)) {
      assignNet(clock_net, net, "clock");
    } else if (isClearPin(iterm)) {
      assignNet(clear_net, net, "clear");
    } else if (isPresetPin(iterm)) {
      assignNet(preset_net, net, "preset");
    } else if (isScanEnablePin(iterm)) {
      assignNet(scan_enable_net, net, "scan enable");
    } else if (isScanInPin(iterm)) {
      assignNet(scan_in_net, net, "scan in");
    }
  }

  return InstMask(clock_net, clear_net, preset_net, scan_enable_net, scan_in_net);
}

MasterMask 
AggloCluster::createMasterMask(odb::dbInst* inst)
{
  // Step 1: Get Liberty cell
  const sta::LibertyCell* lib_cell = getLibertyCell(inst);
  if (!lib_cell) {
    if (verbose_) {
      std::cout << "[createMasterMask] Warning: No Liberty cell found for instance " 
                << inst->getName() << std::endl;
    }
    return MasterMask();
  }

  // Step 2: Extract function expression from sequential element
  const sta::FuncExpr* expr = nullptr;
  if (!lib_cell->sequentials().empty()) {
    expr = lib_cell->sequentials().front()->data();
  }
  
  if (!expr) {
    if (verbose_) {
      std::cout << "[createMasterMask] Warning: No function expression for sequential in instance " 
                << inst->getName() << " (master: " << inst->getMaster()->getName() << ")" << std::endl;
    }
    return MasterMask();
  }

  // Step 3: Extract master properties
  const int func_id = getFuncId(expr, inst);
  const bool has_pos_clk = hasPositiveClockEdge(inst);
  const bool has_clear = hasClear(inst);
  const bool has_preset = hasPreset(inst);
  const bool has_scan = hasScan(inst);

  // Step 4: Check for Q/QN output pins
  bool has_q = false;
  bool has_qn = false;
  
  for (odb::dbITerm* iterm : inst->getITerms()) {
    if (!has_q && isQPin(iterm)) {
      has_q = true;
      if (has_qn) break;  // Found both, early exit
    }
    if (!has_qn && isQNPin(iterm)) {
      has_qn = true;
      if (has_q) break;   // Found both, early exit
    }
  }

  return MasterMask(func_id, has_pos_clk, has_clear, has_preset, has_q, has_qn, has_scan);
}

// -----------------------------------------------------------------------------
// [Level 2] Liberty Cell and Function Analysis
// -----------------------------------------------------------------------------

const sta::LibertyCell* 
AggloCluster::getLibertyCell(odb::dbInst* inst) const
{
  const sta::Cell* cell = network_->dbToSta(inst->getMaster());
  if (!cell) {
    if (verbose_) {
      std::cout << "[getLibertyCell] Warning: Cannot convert master " 
                << inst->getMaster()->getName() << " to STA cell" << std::endl;
    }
    return nullptr;
  }
  
  const sta::LibertyCell* lib_cell = network_->libertyCell(cell);
  if (!lib_cell) {
    if (verbose_) {
      std::cout << "[getLibertyCell] Warning: No Liberty cell for master " 
                << inst->getMaster()->getName() << std::endl;
    }
    return nullptr;
  }
  
  if (const sta::TestCell* test_cell = lib_cell->testCell()) {
    lib_cell = test_cell;
  }
  return lib_cell;
}

int 
AggloCluster::getFuncId(const sta::FuncExpr* expr, odb::dbInst* inst)
{
  std::string func_str = getFuncStr(expr, inst);

  auto it = func_str_to_func_id_.find(func_str);
  if (it != func_str_to_func_id_.end()) {
    return it->second;
  }

  int new_id = func_str_to_func_id_.size();
  func_str_to_func_id_[func_str] = new_id;
  return new_id;
}

std::string 
AggloCluster::getFuncStr(const sta::FuncExpr* expr, odb::dbInst* inst) const
{
  if (!expr) {
    return "()";
  }

  const sta::FuncExpr::Operator op = expr->op();

  // Handle port reference
  if (op == sta::FuncExpr::op_port) {
    const FlopPort p_type = getPortType(expr->port(), inst);
    return "p(" + std::to_string(static_cast<int>(p_type)) + ")";
  }

  const std::string op_str = std::to_string(op);

  // Handle unary operator (NOT)
  if (op == sta::FuncExpr::op_not) {
    return "op" + op_str + "(" + getFuncStr(expr->left(), inst) + ")";
  }

  // Handle binary operators (AND, OR, XOR, etc.)
  return "op" + op_str + "(" + getFuncStr(expr->left(), inst) + "," 
         + getFuncStr(expr->right(), inst) + ")";
}

// -----------------------------------------------------------------------------
// [Level 2] Master Property Checkers
// -----------------------------------------------------------------------------

bool 
AggloCluster::hasClear(odb::dbInst* inst) const
{
  const sta::LibertyCell* lib_cell = getLibertyCell(inst);
  if (!lib_cell) {
    return false;
  }

  return std::any_of(
      lib_cell->sequentials().begin(),
      lib_cell->sequentials().end(),
      [](const sta::Sequential* seq) { return seq->clear() != nullptr; });
}

bool 
AggloCluster::hasPreset(odb::dbInst* inst) const
{
  const sta::LibertyCell* lib_cell = getLibertyCell(inst);
  if (!lib_cell) {
    return false;
  }

  return std::any_of(
      lib_cell->sequentials().begin(),
      lib_cell->sequentials().end(),
      [](const sta::Sequential* seq) { return seq->preset() != nullptr; });
}

bool 
AggloCluster::hasPositiveClockEdge(odb::dbInst* inst) const
{
  const sta::LibertyCell* lib_cell = getLibertyCell(inst);
  if (!lib_cell) {
    return false;
  }

  // Check if any sequential has negative edge (!CLK pattern)
  const bool has_negative_edge = std::any_of(
      lib_cell->sequentials().begin(),
      lib_cell->sequentials().end(),
      [](const sta::Sequential* seq) {
        const sta::FuncExpr* clk = seq->clock();
        if (!clk) return false;
        const sta::FuncExpr* left = clk->left();
        const sta::FuncExpr* right = clk->right();
        return left && !right;  // !CLK pattern (negative edge)
      });

  return !has_negative_edge;
}

bool 
AggloCluster::hasScan(odb::dbInst* inst) const
{
  const sta::LibertyCell* lib_cell = getLibertyCell(inst);
  if (!lib_cell) {
    return false;
  }

  return getLibertyScanIn(lib_cell) && getLibertyScanEnable(lib_cell);
}

// -----------------------------------------------------------------------------
// [Level 2] Pin Counting Utilities
// -----------------------------------------------------------------------------

int 
AggloCluster::getNumDPins(odb::dbInst* inst) const
{
  return std::count_if(
      inst->getITerms().begin(),
      inst->getITerms().end(),
      [this](odb::dbITerm* iterm) { return isDPin(iterm); });
}

int 
AggloCluster::getNumQPins(odb::dbInst* inst) const
{
  return std::count_if(
      inst->getITerms().begin(),
      inst->getITerms().end(),
      [this](odb::dbITerm* iterm) { return isQPin(iterm); });
}

int 
AggloCluster::getNumQNPins(odb::dbInst* inst) const
{
  return std::count_if(
      inst->getITerms().begin(),
      inst->getITerms().end(),
      [this](odb::dbITerm* iterm) { return isQNPin(iterm); });
}

// -----------------------------------------------------------------------------
// [Level 3] Pin Type Checkers
// -----------------------------------------------------------------------------

bool 
AggloCluster::isClockPin(odb::dbITerm* iterm) const
{
  return iterm->getSigType() == odb::dbSigType::CLOCK
         || network_->isRegClkPin(network_->dbToSta(iterm));
}

bool 
AggloCluster::isDPin(odb::dbITerm* iterm) const
{
  const sta::LibertyPort* lib_port = getLibertyPort(iterm);
  if (lib_port == nullptr) {
    return false;
  }

  return iterm->getIoType() == odb::dbIoType::INPUT
         && !(isClockPin(iterm) || isPowerPin(iterm) || isPresetPin(iterm)
              || isClearPin(iterm) || isScanInPin(iterm)
              || isScanEnablePin(iterm));
}

bool 
AggloCluster::isPowerPin(odb::dbITerm* iterm) const
{
  return iterm->getSigType()
      .isSupply();  // dbSigType::isSupply returns true for POWER or GROUND
}

bool 
AggloCluster::isPresetPin(odb::dbITerm* iterm) const
{
  const sta::LibertyPort* lib_port = getLibertyPort(iterm);
  if (!lib_port) {
    return false;
  }

  const sta::LibertyCell* lib_cell = getLibertyCell(iterm->getInst());
  if (!lib_cell) {
    return false;
  }

  return std::any_of(
      lib_cell->sequentials().begin(),
      lib_cell->sequentials().end(),
      [lib_port](const sta::Sequential* seq) {
        return seq->preset() && seq->preset()->hasPort(lib_port);
      });
}

bool 
AggloCluster::isClearPin(odb::dbITerm* iterm) const
{
  const sta::LibertyPort* lib_port = getLibertyPort(iterm);
  if (!lib_port) {
    return false;
  }

  const sta::LibertyCell* lib_cell = getLibertyCell(iterm->getInst());
  if (!lib_cell) {
    return false;
  }

  return std::any_of(
      lib_cell->sequentials().begin(),
      lib_cell->sequentials().end(),
      [lib_port](const sta::Sequential* seq) {
        return seq->clear() && seq->clear()->hasPort(lib_port);
      });
}

bool 
AggloCluster::isQPin(odb::dbITerm* iterm) const
{
  // Step 1: Basic type checks
  if (iterm->getIoType() != odb::dbIoType::OUTPUT 
      || isClockPin(iterm) || isPowerPin(iterm)) {
    return false;
  }

  // Step 2: Get Liberty port and function
  const sta::LibertyPort* lib_port = getLibertyPort(iterm);
  if (!lib_port) {
    return false;
  }

  const sta::FuncExpr* func = getPortFunc(lib_port);
  if (!func) {
    return false;
  }

  const sta::LibertyPort* func_port = func->port();
  if (!func_port) {
    return false;
  }

  // Step 3: Check if this matches a sequential output
  const sta::LibertyCell* lib_cell = getLibertyCell(iterm->getInst());
  if (!lib_cell) {
    return false;
  }

  const std::string pin_func_name = func_port->name();
  return std::any_of(
      lib_cell->sequentials().begin(),
      lib_cell->sequentials().end(),
      [&pin_func_name](const sta::Sequential* seq) {
        const sta::LibertyPort* output = seq->output();
        return output && output->name() == pin_func_name;
      });
}

bool 
AggloCluster::isQNPin(odb::dbITerm* iterm) const
{
  // Step 1: Basic type checks
  if (iterm->getIoType() != odb::dbIoType::OUTPUT 
      || isClockPin(iterm) || isPowerPin(iterm)) {
    return false;
  }

  // Step 2: Get Liberty port and function
  const sta::LibertyPort* lib_port = getLibertyPort(iterm);
  if (!lib_port) {
    return false;
  }

  const sta::FuncExpr* func = getPortFunc(lib_port);
  if (!func) {
    return false;
  }

  const sta::LibertyPort* func_port = func->port();
  if (!func_port) {
    return false;
  }

  // Step 3: Check if this matches a sequential inverted output
  const sta::LibertyCell* lib_cell = getLibertyCell(iterm->getInst());
  if (!lib_cell) {
    return false;
  }

  const std::string pin_func_name = func_port->name();
  return std::any_of(
      lib_cell->sequentials().begin(),
      lib_cell->sequentials().end(),
      [&pin_func_name](const sta::Sequential* seq) {
        const sta::LibertyPort* output_inv = seq->outputInv();
        return output_inv && output_inv->name() == pin_func_name;
      });
}

bool 
AggloCluster::isScanInPin(odb::dbITerm* iterm) const
{
  const sta::LibertyCell* lib_cell = getLibertyCell(iterm->getInst());
  if (!lib_cell) {
    return false;
  }

  const sta::LibertyPort* scan_in_port = sta::getLibertyScanIn(lib_cell);
  if (!scan_in_port) {
    return false;
  }

  odb::dbMTerm* mterm = network_->staToDb(scan_in_port);
  return mterm && iterm->getInst()->getITerm(mterm) == iterm;
}

bool 
AggloCluster::isScanEnablePin(odb::dbITerm* iterm) const
{
  const sta::LibertyCell* lib_cell = getLibertyCell(iterm->getInst());
  if (!lib_cell) {
    return false;
  }

  const sta::LibertyPort* scan_enable_port = sta::getLibertyScanEnable(lib_cell);
  if (!scan_enable_port) {
    return false;
  }

  odb::dbMTerm* mterm = network_->staToDb(scan_enable_port);
  return mterm && iterm->getInst()->getITerm(mterm) == iterm;
}

// -----------------------------------------------------------------------------
// [Level 3] Port Analysis Utilities
// -----------------------------------------------------------------------------

const sta::LibertyPort* 
AggloCluster::getLibertyPort(odb::dbITerm* iterm) const
{
  const sta::Pin* pin = network_->dbToSta(iterm);
  if (pin == nullptr) {
    if (verbose_) {
      const std::string pin_name = iterm->getMTerm()->getName();
      // Skip VDD/VSS warnings
      if (pin_name != "VDD" && pin_name != "VSS") {
        std::cout << "[getLibertyPort] Warning: Cannot convert pin " 
                  << pin_name << " of instance " 
                  << iterm->getInst()->getName() << " to STA pin" << std::endl;
      }
    }
    return nullptr;
  }
  
  const sta::LibertyPort* lib_port = network_->libertyPort(pin);
  if (lib_port == nullptr) {
    if (verbose_) {
      const std::string pin_name = iterm->getMTerm()->getName();
      // Skip VDD/VSS warnings
      if (pin_name != "VDD" && pin_name != "VSS") {
        std::cout << "[getLibertyPort] Warning: No Liberty port for pin " 
                  << pin_name << " of instance " 
                  << iterm->getInst()->getName() << std::endl;
      }
    }
    return nullptr;
  }
  return lib_port;
}

const sta::FuncExpr* 
AggloCluster::getPortFunc(const sta::LibertyPort* port) const
{
  if (!port) {
    if (verbose_) {
      std::cout << "[getPortFunc] Warning: Null port provided" << std::endl;
    }
    return nullptr;
  }

  const sta::FuncExpr* function = port->function();
  if (function) {
    return function;
  }

  // Check parent bundle/bus
  sta::LibertyCellPortIterator port_iter(port->libertyCell());
  while (port_iter.hasNext()) {
    const sta::LibertyPort* next_port = port_iter.next();
    function = next_port->function();
    if (!function) {
      continue;
    }
    if (next_port->hasMembers()) {
      std::unique_ptr<sta::ConcretePortMemberIterator> mem_iter(
          next_port->memberIterator());
      while (mem_iter->hasNext()) {
        const sta::ConcretePort* mem_port = mem_iter->next();
        if (mem_port == port) {
          return function;
        }
      }
    }
  }
  
  if (verbose_) {
    std::cout << "[getPortFunc] Warning: No function found for port " 
              << port->name() << std::endl;
  }
  return nullptr;
}

FlopPort 
AggloCluster::getPortType(const sta::LibertyPort* lib_port, odb::dbInst* inst) const
{
  if (!lib_port) {
    if (verbose_) {
      std::cout << "[getPortType] Warning: Null Liberty port provided" << std::endl;
    }
    return unknown;
  }

  if (!inst) {
    if (verbose_) {
      std::cout << "[getPortType] Warning: Null instance provided" << std::endl;
    }
    return unknown;
  }

  odb::dbMTerm* mterm = network_->staToDb(lib_port);
  if (mterm != nullptr) {
    odb::dbITerm* iterm = inst->getITerm(mterm);
    if (iterm != nullptr) {
      if (isDPin(iterm)) {
        return d;
      }
      if (isQPin(iterm)) {
        return q;
      }
      if (isQNPin(iterm)) {
        return qn;
      }
      if (isClearPin(iterm)) {
        return clear;
      }
      if (isPresetPin(iterm)) {
        return preset;
      }
      if (isScanInPin(iterm)) {
        return si;
      }
      if (isScanEnablePin(iterm)) {
        return se;
      }
      if (isPowerPin(iterm)) {
        if (iterm->getSigType() == odb::dbSigType::GROUND) {
          return vss;
        }
        return vdd;
      }
    }
  }

  const sta::LibertyCell* lib_cell = getLibertyCell(inst);
  if (lib_cell == nullptr) {
    if (verbose_) {
      std::cout << "[getPortType] Warning: Could not find Liberty cell for instance " 
                << inst->getName() << std::endl;
    }
    return unknown;
  }

  for (const sta::Sequential* seq : lib_cell->sequentials()) {
    if (sta::LibertyPort::equiv(lib_port, seq->output())) {
      return func;
    }
    if (sta::LibertyPort::equiv(lib_port, seq->outputInv())) {
      return ifunc;
    }
  }

  if (verbose_) {
    std::cout << "[getPortType] Warning: Could not classify port " << lib_port->name() 
              << " in instance " << inst->getName() << std::endl;
  }
  return unknown;
}

//==============================================================================
// Phase 7 Helper Functions: Feasible Region Calculation 
//==============================================================================

// -----------------------------------------------------------------------------
// [Level 1] Top-level: Compute feasible region for each flop unit
// -----------------------------------------------------------------------------

void
AggloCluster::calcFeasibleRegion(FlopUnit& flop, bool verbose)
{
  // Save original precision settings for high-precision debug output
  const std::streamsize original_precision = std::cout.precision();
  const std::ios_base::fmtflags original_flags = std::cout.flags();
  
  if (verbose) {
    std::cout << std::fixed << std::setprecision(18);
    std::cout << "\n=== Calculating Feasible Region for Flop: " << flop.inst_->getName() << " ===" << std::endl;
  }

  // Step 1: Extract wire parasitic parameters
  const auto est = resizer_->getEstimateParasitics();
  const double unit_c = est->wireSignalCapacitance(corner_);
  const double unit_r = est->wireSignalResistance(corner_);

  if (verbose) {
    std::cout << "  [Step 1] Wire Parameters" << std::endl;
    std::cout << "    Unit resistance: " << unit_r << std::endl;
    std::cout << "    Unit capacitance: " << unit_c << std::endl;
  }

  // Step 2: Classify pins by type
  std::vector<odb::dbITerm*> d_pins, q_pins, qn_pins;
  odb::dbITerm* clk_pin = nullptr;
  
  for (odb::dbITerm* iterm : flop.inst_->getITerms()) {
    if (isDPin(iterm)) {
      d_pins.push_back(iterm);
    } else if (isQPin(iterm)) {
      q_pins.push_back(iterm);
    } else if (isQNPin(iterm)) {
      qn_pins.push_back(iterm);
    } else if (isClockPin(iterm)) {
      clk_pin = iterm;
    }
  }
  
  if (verbose) {
    std::cout << "  [Step 2] Pin Classification" << std::endl;
    std::cout << "    D pins: " << d_pins.size() << std::endl;
    for (size_t i = 0; i < d_pins.size(); ++i) {
      std::cout << "      [" << i << "] " << d_pins[i]->getMTerm()->getName() << std::endl;
    }
    std::cout << "    Q pins: " << q_pins.size() << std::endl;
    for (size_t i = 0; i < q_pins.size(); ++i) {
      std::cout << "      [" << i << "] " << q_pins[i]->getMTerm()->getName() << std::endl;
    }
    std::cout << "    QN pins: " << qn_pins.size() << std::endl;
    for (size_t i = 0; i < qn_pins.size(); ++i) {
      std::cout << "      [" << i << "] " << qn_pins[i]->getMTerm()->getName() << std::endl;
    }
    std::cout << "    Clock pin: " << (clk_pin ? clk_pin->getMTerm()->getName() : "NOT FOUND") << std::endl;
  }
  
  // Validate clock pin existence
  if (!clk_pin) {
    if (verbose) {
      std::cout << "    [WARNING] No clock pin found - skipping this flop" << std::endl;
    }
    return;
  }
  
  // Step 3: Process output pins (Q/QN - fanout constraints)
  std::vector<odb::dbITerm*> output_pins;
  output_pins.reserve(q_pins.size() + qn_pins.size());
  output_pins.insert(output_pins.end(), q_pins.begin(), q_pins.end());
  output_pins.insert(output_pins.end(), qn_pins.begin(), qn_pins.end());
  
  for (odb::dbITerm* out_pin : output_pins) {
    processFanOutPin(flop, out_pin, clk_pin->getMTerm(), est, unit_r, unit_c, verbose);
  }

  // Step 4: Process input pins (D - fanin constraints)
  for (odb::dbITerm* d_pin : d_pins) {
    processFanInPin(flop, d_pin, est, unit_r, unit_c, verbose);
  }
  
  // Step 5: Compute final feasible region by intersecting all pin constraints
  std::vector<odb::dbITerm*> all_pins;
  all_pins.reserve(d_pins.size() + q_pins.size() + qn_pins.size());
  all_pins.insert(all_pins.end(), d_pins.begin(), d_pins.end());
  all_pins.insert(all_pins.end(), q_pins.begin(), q_pins.end());
  all_pins.insert(all_pins.end(), qn_pins.begin(), qn_pins.end());
  
  computeFinalFeasibleRegion(flop, all_pins, verbose);
  
  // Restore original precision settings
  if (verbose) {
    std::cout.precision(original_precision);
    std::cout.flags(original_flags);
  }
}

// -----------------------------------------------------------------------------
// [Level 2] Main processing functions for different pin types
// -----------------------------------------------------------------------------

void
AggloCluster::processFanOutPin(FlopUnit& flop,
                               odb::dbITerm* out_pin,
                               odb::dbMTerm* clk_pin_lib,
                               est::EstimateParasitics* est,
                               double unit_r,
                               double unit_c,
                               bool verbose)
{
  if (verbose) {
    std::cout << "\n  [Step 3.1] Processing Fan-Out Pin: " << out_pin->getMTerm()->getName() << std::endl;
  }
  
  // Step 1: Check timing constraints
  const auto budget_it = flop.pin_budgets_.find(out_pin);
  if (budget_it == flop.pin_budgets_.end() || budget_it->second.empty()) {
    if (verbose) {
      std::cout << "    No timing constraints found" << std::endl;
      std::cout << "    Result: Using unconstrained region (inverse box)" << std::endl;
    }
    flop.pin_feasible_regions_[out_pin] = boost::geometry::make_inverse<Box>();
    return;
  }
  
  if (verbose) {
    std::cout << "    Total timing paths for this pin: " << budget_it->second.size() << std::endl;
  }
  
  // Step 2: Extract most critical path (already sorted in ascending order)
  const auto& [critical_path_idx, slack_budget] = budget_it->second[0];
  
  if (verbose) {
    std::cout << "    Most critical path:" << std::endl;
    std::cout << "      Path index: " << critical_path_idx << std::endl;
    std::cout << "      Slack budget: " << slack_budget << std::endl;
  }
  
  // Step 3: Validate net connection
  odb::dbNet* fanout_net = out_pin->getNet();
  if (!fanout_net) {
    if (verbose) {
      std::cout << "    Pin not connected to any net" << std::endl;
      std::cout << "    Result: Skipping this pin" << std::endl;
    }
    return;
  }
  
  if (verbose) {
    std::cout << "    Fan-out net: " << fanout_net->getName() << std::endl;
  }
  
  // Step 4: Extract cell delay characterization
  const auto cap_delay = extractCapacitanceDelayPoints(
      flop.inst_, 
      clk_pin_lib->getName(), 
      out_pin->getMTerm()->getName(), 
      0);
  
  std::vector<float> coeffs;
  coeffs.reserve(cap_delay.size() - 1);
  for (size_t i = 1; i < cap_delay.size(); ++i) {
    const float dy = cap_delay[i].second - cap_delay[i-1].second;
    const float dx = cap_delay[i].first - cap_delay[i-1].first;
    coeffs.push_back(dy / dx);
  }
  
  if (verbose) {
    std::cout << "    Cell delay analysis:" << std::endl;
    std::cout << "      Cap-delay points: " << cap_delay.size() << std::endl;
    for (size_t i = 0; i < cap_delay.size(); ++i) {
      std::cout << "      Point[" << i << "]: cap=" << cap_delay[i].first
                << ", delay=" << cap_delay[i].second << std::endl;
    }
    std::cout << "      Delay coefficients (∂delay/∂cap):" << std::endl;
    for (size_t i = 0; i < coeffs.size(); ++i) {
      std::cout << "        Segment [" << i << "-" << (i+1) << "]: " << coeffs[i] << std::endl;
    }
  }
  
  // Step 5: Build Steiner tree and compute parasitics
  const sta::Pin* driver_pin_sta = network_->dbToSta(out_pin);
  est::SteinerTree* driver_steiner_tree = est->makeSteinerTree(driver_pin_sta);
  
  const sta::Net* fanout_net_sta = network_->dbToSta(fanout_net);
  float pin_capacitance = 0.0f, wire_capacitance_total = 0.0f;
  const sta::MinMax* mm = sta::MinMax::max();
  sta_->connectedCap(fanout_net_sta, corner_, mm, pin_capacitance, wire_capacitance_total);
  const float total_net_capacitance = pin_capacitance + wire_capacitance_total;
  
  const auto top_steiner_point = driver_steiner_tree->top();
  const auto driver_point = driver_steiner_tree->drvrPt();
  const auto steiner_location = driver_steiner_tree->location(top_steiner_point);
  const auto driver_location = driver_steiner_tree->location(driver_point);
  const double l1 = dbuToMeters(driver_steiner_tree->distance(driver_point, top_steiner_point));
  
  const float wire_cap = l1 * unit_c;
  const float wire_res = l1 * unit_r;
  const float capacitance_beyond_first_steiner = total_net_capacitance - wire_cap;
  
  if (verbose) {
    std::cout << "    Steiner tree analysis:" << std::endl;
    std::cout << "      Driver location: (" << driver_location.getX() << ", " << driver_location.getY() << ")" << std::endl;
    std::cout << "      Top Steiner point: (" << steiner_location.getX() << ", " << steiner_location.getY() << ")" << std::endl;
    std::cout << "      Distance to Steiner (l1): " << l1 << std::endl;
    std::cout << "    Net capacitance:" << std::endl;
    std::cout << "      Pin capacitance: " << pin_capacitance << std::endl;
    std::cout << "      Wire capacitance: " << wire_capacitance_total << std::endl;
    std::cout << "      Total capacitance: " << total_net_capacitance << std::endl;
    std::cout << "    Wire parasitics (to Steiner):" << std::endl;
    std::cout << "      Wire capacitance: " << wire_cap << std::endl;
    std::cout << "      Wire resistance: " << wire_res << std::endl;
    std::cout << "      Remaining capacitance: " << capacitance_beyond_first_steiner << std::endl;
  }
  
  // Step 6: Solve for maximum distance using coefficient segments
  float max_dist = -1.0f;
  bool found_valid_segment = false;
  
  if (verbose) {
    std::cout << "    Solving for maximum distance:" << std::endl;
  }
  
  for (size_t seg_idx = 0; seg_idx < coeffs.size(); ++seg_idx) {
    const float coeff = coeffs[seg_idx];
    const float segment_cap_min = cap_delay[seg_idx].first;
    const float segment_cap_max = cap_delay[seg_idx + 1].first;
    const bool is_last_iteration = (seg_idx == coeffs.size() - 1);
    
    if (verbose) {
      std::cout << "      Trying segment [" << seg_idx << "-" << (seg_idx+1) << "]:" << std::endl;
      std::cout << "        Coefficient: " << coeff << std::endl;
      std::cout << "        Cap range: [" << segment_cap_min << ", " << segment_cap_max << "]" << std::endl;
    }
    
    const float candidate_dist = solveMaxDistanceFanOut(
        l1, unit_r, unit_c, slack_budget, coeff, capacitance_beyond_first_steiner, verbose);
    
    if (candidate_dist < 0) {
      if (verbose) {
        std::cout << "        ✗ Invalid solution (negative or no real root)" << std::endl;
      }
      continue;
    }
    
    const float resulting_cap = candidate_dist * unit_c + capacitance_beyond_first_steiner;
    
    if (verbose) {
      std::cout << "        Candidate distance: " << candidate_dist << std::endl;
      std::cout << "        Resulting capacitance: " << resulting_cap << std::endl;
    }
    
    if (resulting_cap >= segment_cap_min && resulting_cap <= segment_cap_max) {
      max_dist = candidate_dist;
      found_valid_segment = true;
      
      if (verbose) {
        std::cout << "        ✓ Valid! Capacitance is within segment range." << std::endl;
      }
      break;
    } else {
      if (verbose) {
        std::cout << "        ✗ Invalid. Capacitance outside segment range." << std::endl;
      }
      
      if (is_last_iteration) {
        max_dist = candidate_dist;
        found_valid_segment = true;
        if (verbose) {
          std::cout << "        ⚠ Last iteration: Using candidate_dist despite out-of-range cap: " << candidate_dist << std::endl;
        }
      }
    }
  }
  
  if (!found_valid_segment) {
    max_dist = l1;
    if (verbose) {
      std::cout << "      Warning: No valid segment found, using l1 as default distance: " << l1 << std::endl;
    }
  }
  
  max_dist = metersToDbu(max_dist);

  // Step 7: Apply bin size constraint
  const double bin_size_x = virtual_bin_grid_.getBinSizeX();
  const double bin_size_y = virtual_bin_grid_.getBinSizeY();
  const double avg_bin_size = (bin_size_x + bin_size_y) / 2.0;
  const float l1_dbu = metersToDbu(l1);
  const float lower_limit = l1_dbu;
  const float upper_limit = l1_dbu + (avg_bin_size * region_scale_factor_);
  
  if (verbose) {
    std::cout << "      Applying distance constraint:" << std::endl;
    std::cout << "        Calculated max_dist: " << max_dist << " DBU" << std::endl;
    std::cout << "        Lower bound (l1): " << lower_limit << " DBU" << std::endl;
    std::cout << "        Bin size (avg): " << avg_bin_size << " DBU" << std::endl;
    std::cout << "        Region scale factor: " << region_scale_factor_ << std::endl;
    std::cout << "        Upper bound (l1 + bin*k): " << upper_limit << " DBU" << std::endl;
  }
  
  max_dist = std::clamp(max_dist, lower_limit, upper_limit);

  // Step 8: Create feasible region box
  const int steiner_x = steiner_location.getX();
  const int steiner_y = steiner_location.getY();
  
  if (verbose) {
    std::cout << "    Maximum placement distance:" << std::endl;
    std::cout << "      Final max_dist: " << max_dist << " DBU" << std::endl;
    std::cout << "    Result: Feasible region centered at (" << steiner_x << ", " << steiner_y << ")" << std::endl;
  }
  
  flop.pin_feasible_regions_[out_pin] = createFeasibleBox(steiner_x, steiner_y, max_dist);
}

void
AggloCluster::processFanInPin(FlopUnit& flop,
                              odb::dbITerm* d_pin,
                              est::EstimateParasitics* est,
                              double unit_r,
                              double unit_c,
                              bool verbose)
{
  if (verbose) {
    std::cout << "\n  [Step 4.1] Processing Fan-In Pin: " << d_pin->getMTerm()->getName() << std::endl;
  }
  
  // Step 1: Check timing constraints
  const auto budget_it = flop.pin_budgets_.find(d_pin);
  if (budget_it == flop.pin_budgets_.end() || budget_it->second.empty()) {
    if (verbose) {
      std::cout << "    No timing constraints found" << std::endl;
      std::cout << "    Result: Using unconstrained region (inverse box)" << std::endl;
    }
    flop.pin_feasible_regions_[d_pin] = boost::geometry::make_inverse<Box>();
    return;
  }
  
  if (verbose) {
    std::cout << "    Total timing paths for this pin: " << budget_it->second.size() << std::endl;
  }
  
  // Step 2: Extract most critical path (already sorted in ascending order)
  const auto& [critical_path_idx, slack_budget] = budget_it->second[0];
  
  if (verbose) {
    std::cout << "    Most critical path:" << std::endl;
    std::cout << "      Path index: " << critical_path_idx << std::endl;
    std::cout << "      Slack budget: " << slack_budget << std::endl;
  }
  
  // Step 3: Validate net connection
  odb::dbNet* fi_net = d_pin->getNet();
  if (!fi_net) {
    if (verbose) {
      std::cout << "    Pin not connected to any net" << std::endl;
      std::cout << "    Result: Skipping this pin" << std::endl;
    }
    return;
  }
  
  if (verbose) {
    std::cout << "    Fan-in net: " << fi_net->getName() << std::endl;
  }
  
  // Step 4: Get driver pin
  odb::dbITerm* fi_net_drvr_pin = fi_net->get1stITerm();
  if (!fi_net_drvr_pin) {
    if (verbose) {
      std::cout << "    No driver pin found on net" << std::endl;
      std::cout << "    Result: Skipping this pin" << std::endl;
    }
    return;
  }
  
  if (verbose) {
    std::cout << "    Driver instance: " << fi_net_drvr_pin->getInst()->getName() << std::endl;
    std::cout << "    Driver pin: " << fi_net_drvr_pin->getMTerm()->getName() << std::endl;
  }
  
  // Step 5: Calculate D pin capacitance
  const sta::Pin* ipin_sta = network_->dbToSta(d_pin);
  const float ipin_cap = getPinCapacitance(ipin_sta);
  
  if (verbose) {
    std::cout << "    D pin input capacitance: " << ipin_cap << std::endl;
  }
  
  // Step 6: Extract cell delay characterization from driver
  odb::dbInst* fi_inst = fi_net_drvr_pin->getInst();
  const auto slews = getInstanceInputSlews(fi_inst);
  
  std::vector<std::pair<float, float>> cap_delay;
  float worst_delay = 0.0f;
  
  if (verbose) {
    std::cout << "    Driver cell delay analysis:" << std::endl;
    std::cout << "      Input pins with slews: " << slews.size() << std::endl;
  }
  
  for (const auto& [iterm, slew] : slews) {
    const auto pts = extractCapacitanceDelayPoints(
        fi_inst,
        iterm->getMTerm()->getName(),
        fi_net_drvr_pin->getMTerm()->getName(),
        slew);
    
    if (!pts.empty()) {
      const float tail_second = pts.back().second;
      if (tail_second >= worst_delay) {
        worst_delay = tail_second;
        cap_delay = pts;
      }
    }
  }
  
  if (cap_delay.empty() || cap_delay.size() < 2) {
    if (verbose) {
      std::cout << "    Insufficient cap-delay data points" << std::endl;
      std::cout << "    Result: Skipping this pin" << std::endl;
    }
    return;
  }
  
  // Step 7: Calculate delay coefficients
  std::vector<float> coeffs;
  coeffs.reserve(cap_delay.size() - 1);
  for (size_t i = 1; i < cap_delay.size(); ++i) {
    const float dy = cap_delay[i].second - cap_delay[i-1].second;
    const float dx = cap_delay[i].first - cap_delay[i-1].first;
    coeffs.push_back(dy / dx);
  }
  
  if (verbose) {
    std::cout << "      Cap-delay points: " << cap_delay.size() << std::endl;
    for (size_t i = 0; i < cap_delay.size(); ++i) {
      std::cout << "      Point[" << i << "]: cap=" << cap_delay[i].first
                << ", delay=" << cap_delay[i].second << std::endl;
    }
    std::cout << "      Delay coefficients (∂delay/∂cap):" << std::endl;
    for (size_t i = 0; i < coeffs.size(); ++i) {
      std::cout << "        Segment [" << i << "-" << (i+1) << "]: " << coeffs[i] << std::endl;
    }
  }
  
  // Step 8: Build Steiner tree and compute net capacitance
  const sta::Pin* fi_net_drvr_pin_sta = network_->dbToSta(fi_net_drvr_pin);
  est::SteinerTree* fi_tree = est->makeSteinerTree(fi_net_drvr_pin_sta);
  
  float pin_cap = 0.0f, wire_cap = 0.0f;
  const sta::MinMax* mm = sta::MinMax::max();
  const sta::Net* fi_net_drvr_sta = network_->dbToSta(fi_net_drvr_pin->getNet());
  sta_->connectedCap(fi_net_drvr_sta, corner_, mm, pin_cap, wire_cap);
  const float total_cap = pin_cap + wire_cap;

  if (verbose) {
    std::cout << "    Net capacitance:" << std::endl;
    std::cout << "      Pin capacitance: " << pin_cap << std::endl;
    std::cout << "      Wire capacitance: " << wire_cap << std::endl;
    std::cout << "      Total capacitance: " << total_cap << std::endl;
  }

  const auto top_pt = fi_tree->top();
  const auto top_loc = fi_tree->location(top_pt);
  const int top_x = top_loc.getX();
  const int top_y = top_loc.getY();
  
  if (verbose) {
    std::cout << "    Steiner tree analysis:" << std::endl;
    std::cout << "      Top Steiner point: (" << top_x << ", " << top_y << ")" << std::endl;
    std::cout << "      Branch count: " << fi_tree->branchCount() << std::endl;
    std::cout << "      Pin count: " << fi_tree->pinCount() << std::endl;
  }
  
  // Step 9: Find path in Steiner tree to target pin
  const int branch_count = fi_tree->branchCount();
  const int pin_count = fi_tree->pinCount();
  
  int target_pt = -1;
  for (int i = 0; i < pin_count; i++) {
    if (fi_tree->pin(i) == ipin_sta) {
      target_pt = i;
      break;
    }
  }
  
  if (target_pt == -1) {
    if (verbose) {
      std::cout << "    Target pin not found in Steiner tree" << std::endl;
      std::cout << "    Result: Skipping this pin" << std::endl;
    }
    return;
  }
  
  const int drvr_pt = fi_tree->drvrPt();
  std::vector<int> node_path;
  
  if (!findSteinerPathRecursive(fi_tree, drvr_pt, target_pt, node_path)) {
    if (verbose) {
      std::cout << "    Failed to find path in Steiner tree" << std::endl;
      std::cout << "    Result: Skipping this pin" << std::endl;
    }
    return;
  }
  
  if (node_path.empty() || node_path.size() < 2) {
    if (verbose) {
      std::cout << "    Node path too short (size: " << node_path.size() << ")" << std::endl;
      std::cout << "    Result: Skipping this pin" << std::endl;
    }
    return;
  }
  
  if (verbose) {
    std::cout << "    Path finding:" << std::endl;
    std::cout << "      Target pin index: " << target_pt << std::endl;
    std::cout << "      Driver point index: " << drvr_pt << std::endl;
    std::cout << "      Path nodes: " << node_path.size() << std::endl;
  }
  
  // Step 10: Calculate path segments
  const size_t total_segments = node_path.size() - 1;
  float pre_leaf_length_dbu = 0.0f;
  float final_segment_length_dbu = 0.0f;
  
  for (size_t k = 0; k < total_segments; ++k) {
    const int path_node1 = node_path[k];
    const int path_node2 = node_path[k + 1];
    
    int current_wire_length = 0;
    
    for (int i = 0; i < branch_count; i++) {
      odb::Point pt1, pt2;
      int steiner_pt1, steiner_pt2, wire_length;
      fi_tree->branch(i, pt1, steiner_pt1, pt2, steiner_pt2, wire_length);
      
      if ((steiner_pt1 == path_node1 && steiner_pt2 == path_node2) ||
          (steiner_pt1 == path_node2 && steiner_pt2 == path_node1)) {
        current_wire_length = wire_length;
        break;
      }
    }
    
    if (k < total_segments - 1) {
      pre_leaf_length_dbu += current_wire_length;
    } else {
      final_segment_length_dbu = current_wire_length;
    }
  }
  
  // Step 11: Solve for maximum distance using coefficient segments
  const float l1 = dbuToMeters(final_segment_length_dbu);
  const float on_path_R_wo_last = unit_r * dbuToMeters(pre_leaf_length_dbu);
  
  if (verbose) {
    std::cout << "    Wire path segmentation:" << std::endl;
    std::cout << "      Total segments: " << total_segments << std::endl;
    std::cout << "      Pre-leaf segments: " << (total_segments - 1) << std::endl;
    std::cout << "      Pre-leaf length: " << dbuToMeters(pre_leaf_length_dbu) << std::endl;
    std::cout << "      Final segment length: " << l1 << std::endl;
    std::cout << "      Pre-leaf resistance: " << on_path_R_wo_last << std::endl;
  }
  
  float max_dist = -1.0f;
  bool found_valid_segment = false;
  
  if (verbose) {
    std::cout << "    Solving for maximum distance:" << std::endl;
  }
  
  for (size_t seg_idx = 0; seg_idx < coeffs.size(); ++seg_idx) {
    const float coeff = coeffs[seg_idx];
    const float segment_cap_min = cap_delay[seg_idx].first;
    const float segment_cap_max = cap_delay[seg_idx + 1].first;
    const bool is_last_iteration = (seg_idx == coeffs.size() - 1);
    
    if (verbose) {
      std::cout << "      Trying segment [" << seg_idx << "-" << (seg_idx+1) 
                << "], coeff=" << coeff 
                << ", cap_range=[" << segment_cap_min << ", " << segment_cap_max << "]" << std::endl;
    }
    
    const float candidate_dist = solveMaxDistanceFanIn(
        l1, unit_r, unit_c, slack_budget, coeff, on_path_R_wo_last, ipin_cap, verbose);
    
    if (candidate_dist < 0) {
      if (verbose) {
        std::cout << "        ✗ Invalid solution (negative or no real root)" << std::endl;
      }
      continue;
    }
    
    const float resulting_cap = total_cap + (candidate_dist - l1) * unit_c;
    
    if (verbose) {
      std::cout << "      Candidate distance: " << candidate_dist << ", resulting_cap: " << resulting_cap << std::endl;
    }
    
    if (resulting_cap >= segment_cap_min && resulting_cap <= segment_cap_max) {
      max_dist = candidate_dist;
      found_valid_segment = true;
      if (verbose) {
        std::cout << "      ✓ Valid segment found! Using max_dist = " << max_dist << std::endl;
      }
      break;
    } else {
      if (verbose) {
        std::cout << "      ✗ Resulting cap out of range for this segment" << std::endl;
      }
      
      if (is_last_iteration) {
        max_dist = candidate_dist;
        found_valid_segment = true;
        if (verbose) {
          std::cout << "      ⚠ Last iteration: Using candidate_dist despite out-of-range cap: " << candidate_dist << std::endl;
        }
      }
    }
  }
  
  if (!found_valid_segment) {
    max_dist = l1;
    if (verbose) {
      std::cout << "      Warning: No valid segment found, using l1 as default distance: " << l1 << std::endl;
    }
  }
  
  max_dist = metersToDbu(max_dist);

  // Step 12: Apply bin size constraint
  const double bin_size_x = virtual_bin_grid_.getBinSizeX();
  const double bin_size_y = virtual_bin_grid_.getBinSizeY();
  const double avg_bin_size = (bin_size_x + bin_size_y) / 2.0;
  const float l1_dbu = metersToDbu(l1);
  const float lower_limit = l1_dbu;
  const float upper_limit = l1_dbu + (avg_bin_size * region_scale_factor_);
  
  if (verbose) {
    std::cout << "      Applying distance constraint:" << std::endl;
    std::cout << "        Calculated max_dist: " << max_dist << " DBU" << std::endl;
    std::cout << "        Lower bound (l1): " << lower_limit << " DBU" << std::endl;
    std::cout << "        Bin size (avg): " << avg_bin_size << " DBU" << std::endl;
    std::cout << "        Region scale factor: " << region_scale_factor_ << std::endl;
    std::cout << "        Upper bound (l1 + bin*k): " << upper_limit << " DBU" << std::endl;
  }
  
  max_dist = std::clamp(max_dist, lower_limit, upper_limit);

  // Step 13: Create feasible region box
  if (verbose) {
    std::cout << "    Maximum placement distance:" << std::endl;
    std::cout << "      Final max_dist: " << max_dist << " DBU" << std::endl;
    std::cout << "    Result: Feasible region centered at (" << top_x << ", " << top_y << ")" << std::endl;
  }
  
  flop.pin_feasible_regions_[d_pin] = createFeasibleBox(top_x, top_y, max_dist);
}

void
AggloCluster::computeFinalFeasibleRegion(FlopUnit& flop,
                                         const std::vector<odb::dbITerm*>& all_pins,
                                         bool verbose)
{
  if (verbose) {
    std::cout << "\n  [Step 5] Computing Final Feasible Region" << std::endl;
    std::cout << "    Total pins to process: " << all_pins.size() << std::endl;
  }
  
  Box final_region = boost::geometry::make_inverse<Box>();
  bool has_valid_region = false;
  
  // Intersect all pin feasible regions
  for (odb::dbITerm* pin : all_pins) {
    const auto it = flop.pin_feasible_regions_.find(pin);
    if (it == flop.pin_feasible_regions_.end()) {
      continue;
    }
    
    const Box& pin_box = it->second;
    
    // Skip inverse boxes (unconstrained regions)
    const bool is_inverse_box = pin_box.min_corner().get<0>() > pin_box.max_corner().get<0>() ||
                                 pin_box.min_corner().get<1>() > pin_box.max_corner().get<1>();
    
    if (is_inverse_box) {
      if (verbose) {
        std::cout << "    Pin " << pin->getMTerm()->getName() 
                  << " has inverse box (no constraints) - skipping" << std::endl;
      }
      continue;
    }
    
    if (verbose) {
      std::cout << "    Pin " << pin->getMTerm()->getName() 
                << " has valid region: [(" 
                << pin_box.min_corner().get<0>() << ", " << pin_box.min_corner().get<1>() 
                << ") - (" 
                << pin_box.max_corner().get<0>() << ", " << pin_box.max_corner().get<1>() 
                << ")]" << std::endl;
    }
    
    if (!has_valid_region) {
      final_region = pin_box;
      has_valid_region = true;
    } else {
      Box intersection;
      boost::geometry::intersection(final_region, pin_box, intersection);
      
      if (boost::geometry::is_empty(intersection)) {
        if (verbose) {
          std::cout << "      WARNING: No intersection - constraints conflict!" << std::endl;
          std::cout << "      Setting feasible region to empty box (no valid placement)" << std::endl;
        }
        flop.feasible_region_ = Box();
        return;
      }
      
      final_region = intersection;
      if (verbose) {
        std::cout << "      Intersection result: [(" 
                  << final_region.min_corner().get<0>() << ", " << final_region.min_corner().get<1>() 
                  << ") - (" 
                  << final_region.max_corner().get<0>() << ", " << final_region.max_corner().get<1>() 
                  << ")]" << std::endl;
      }
    }
  }
  
  // If all pins have no constraints, use the flop's current location as feasible region
  if (!has_valid_region) {
    // Get flop's original location
    const int flop_x = static_cast<int>(flop.orig_pt_.x);
    const int flop_y = static_cast<int>(flop.orig_pt_.y);
    
    // Transform to rotated coordinate system
    const Point flop_pos(flop_x, flop_y);
    const Point flop_pos_transformed = transformCoords(flop_pos);
    
    const int u = flop_pos_transformed.get<0>();
    const int v = flop_pos_transformed.get<1>();
    
    // Create a point box at the flop's location
    final_region = Box(Point(u, v), Point(u, v));
    
    if (verbose) {
      std::cout << "    All pins have no constraints - using flop's current location as feasible region" << std::endl;
      std::cout << "    Flop location (XY): (" << flop_x << ", " << flop_y << ")" << std::endl;
      std::cout << "    Transformed (UV): (" << u << ", " << v << ")" << std::endl;
    }
  }
  
  flop.feasible_region_ = final_region;
  
  if (verbose) {
    if (has_valid_region) {
      std::cout << "    Final feasible region: [(" 
                << final_region.min_corner().get<0>() << ", " 
                << final_region.min_corner().get<1>() << ") - (" 
                << final_region.max_corner().get<0>() << ", " 
                << final_region.max_corner().get<1>() << ")]" << std::endl;
    } else {
      std::cout << "    Final feasible region (flop location): [(" 
                << final_region.min_corner().get<0>() << ", " 
                << final_region.min_corner().get<1>() << ") - (" 
                << final_region.max_corner().get<0>() << ", " 
                << final_region.max_corner().get<1>() << ")]" << std::endl;
    }
  }
}

// -----------------------------------------------------------------------------
// [Level 3] Utility functions for timing and parasitic analysis
// -----------------------------------------------------------------------------

float
AggloCluster::getPinCapacitance(const sta::Pin* pin) const
{
  // Step 1: Get Liberty port for the pin
  sta::LibertyPort* port = network_->libertyPort(pin);
  if (!port) {
    return 0.0f;  // No Liberty model available
  }
  
  // Step 2: Get corner-specific port (corner 0 = default corner)
  sta::LibertyPort* corner_port = port->cornerPort(0);
  if (!corner_port) {
    return 0.0f;  // No corner-specific data
  }
  
  // Step 3: Return capacitance value
  return corner_port->capacitance();
}

std::vector<std::pair<float, float>>
AggloCluster::extractCapacitanceDelayPoints(odb::dbInst* inst,
                                            const std::string& input_pin_name,
                                            const std::string& output_pin_name,
                                            float input_slew) const
{
  std::vector<std::pair<float, float>> cap_delay_points;

  // Step 1: Get Liberty cell representation
  const sta::LibertyCell* liberty_cell = getLibertyCell(inst);
  if (!liberty_cell) {
    return cap_delay_points;  // Empty result if cell not found
  }

  // Step 2: Find timing arc and capacitance axis for input->output path
  const sta::TableAxis* capacitance_axis = nullptr;
  sta::GateTableModel* timing_model = nullptr;

  if (!findTimingArcModel(liberty_cell, input_pin_name, output_pin_name, 
                          capacitance_axis, timing_model)) {
    return cap_delay_points;  // No matching timing arc found
  }

  // Step 3: Get PVT conditions and calculate delays
  const sta::Pvt* pvt_conditions = liberty_cell->libertyLibrary()->defaultOperatingConditions();
  const auto* capacitance_values = capacitance_axis->values();
  
  cap_delay_points.reserve(capacitance_values->size());
  
  for (float output_capacitance : *capacitance_values) {
    float delay = calculateGateDelay(timing_model, pvt_conditions, 
                                     input_slew, output_capacitance);
    cap_delay_points.emplace_back(output_capacitance, delay);
  }

  return cap_delay_points;
}

std::vector<std::pair<odb::dbITerm*, float>>
AggloCluster::getInstanceInputSlews(odb::dbInst* inst) const
{
  std::vector<std::pair<odb::dbITerm*, float>> input_slews;
  
  // Step 1: Get timing graph and analysis parameters
  sta::Graph* timing_graph = sta_->graph();
  const sta::MinMax* worst_case = sta::MinMax::max();

  // Step 2: Iterate through all terminals of the instance
  for (odb::dbITerm* terminal : inst->getITerms()) {
    // Step 2a: Filter input pins only (skip outputs and power pins)
    if (terminal->getIoType() != odb::dbIoType::INPUT) {
      continue;
    }

    // Step 2b: Get slew value for this input pin
    const float slew = getTerminalSlew(terminal, timing_graph, worst_case);
    
    // Step 2c: Include only pins with valid (non-zero) slew values
    // Zero slew indicates unconnected or unanalyzed pins
    if (slew > 0.0f) {
      input_slews.emplace_back(terminal, slew);
    }
  }

  return input_slews;
}

bool
AggloCluster::findSteinerPathRecursive(est::SteinerTree* tree,
                                       int current_pt,
                                       int target_pt,
                                       std::vector<int>& path)
{
  // Step 1: Check for invalid node (base case)
  if (current_pt == -1) {
    return false;
  }

  // Step 2: Add current node to path
  path.push_back(current_pt);

  // Step 3: Check if target found (base case)
  if (current_pt == target_pt) {
    return true;
  }

  // Step 4: Search left subtree (recursive case)
  if (findSteinerPathRecursive(tree, tree->left(current_pt), target_pt, path)) {
    return true;
  }

  // Step 5: Search right subtree (recursive case)
  if (findSteinerPathRecursive(tree, tree->right(current_pt), target_pt, path)) {
    return true;
  }

  // Step 6: Backtrack - remove current node if target not found in subtrees
  path.pop_back();
  return false;
}

// -----------------------------------------------------------------------------
// [Level 3] Quadratic equation solving for max distance calculation
// -----------------------------------------------------------------------------

float
AggloCluster::solveMaxDistanceFanOut(float l1,
                               float unit_r,
                               float unit_c,
                               float slack_budget,
                               float coeff,
                               float wo_fst_stt_cap,
                               bool verbose) const
{
  // Step 1: Build quadratic equation coefficients
  // Solve: a*x^2 + b*x + c = 0 where x is new Manhattan distance (meters)
  const auto a = unit_r * unit_c / 2;
  const auto b = wo_fst_stt_cap * unit_r + coeff * unit_c / 2;
  const auto c = -std::pow(l1, 2) * unit_r * unit_c / 2 
                 - slack_budget 
                 - wo_fst_stt_cap * l1 * unit_r 
                 - coeff * l1 * unit_c / 2;
  
  // Step 2: Calculate discriminant
  const float D = std::pow(b, 2) - 4 * a * c;
  
  if (verbose) {
    std::cout << "    Solving quadratic equation for maximum distance (FanOut):" << std::endl;
    std::cout << "      Coefficients: a = " << a << ", b = " << b << ", c = " << c << std::endl;
    std::cout << "      Discriminant: " << D << std::endl;
  }
  
  // Step 3: Check for real roots
  if (D < 0) {
    if (verbose) {
      std::cout << "      No real roots - discriminant is negative" << std::endl;
      std::cout << "      Returning invalid value (-1.0)" << std::endl;
    }
    return -1.0f;
  }
  
  // Step 4: Calculate roots using quadratic formula
  const float root1 = (-b + sqrt(D)) / (2 * a);
  const float root2 = (-b - sqrt(D)) / (2 * a);
  
  if (verbose) {
    std::cout << "      Root 1: " << root1 << std::endl;
    std::cout << "      Root 2: " << root2 << std::endl;
  }
  
  // Step 5: Select larger positive root
  if (root1 > 0 || root2 > 0) {
    const float max_dist = std::max(root1, root2);
    if (verbose) {
      std::cout << "      Selected maximum distance: " << max_dist << std::endl;
    }
    return max_dist;
  }
  
  if (verbose) {
    std::cout << "      Both roots are negative - returning invalid value (-1.0)" << std::endl;
  }
  return -1.0f;
}

float
AggloCluster::solveMaxDistanceFanIn(float l1,
                                   float unit_r,
                                   float unit_c,
                                   float slack_budget,
                                   float coeff,
                                   float on_path_R_wo_last,
                                   float ipin_cap,
                                   bool verbose) const
{
  // Step 1: Build quadratic equation coefficients for Fan-In case
  // Solve: a*x^2 + b*x + c = 0 where x is new Manhattan distance (meters)
  const auto a = unit_r * unit_c / 2;
  const auto b = on_path_R_wo_last * unit_c + coeff * unit_c / 2 + unit_r * ipin_cap;
  const auto c = -pow(l1, 2) * unit_r * unit_c / 2 
                 - on_path_R_wo_last * unit_c * l1 
                 - coeff * l1 * unit_c / 2 
                 - ipin_cap * l1 * unit_r
                 - slack_budget;
  
  // Step 2: Calculate discriminant
  const float D = pow(b, 2) - 4 * a * c;
  
  if (verbose) {
    std::cout << "    Solving quadratic equation for maximum distance (FanIn):" << std::endl;
    std::cout << "      Coefficients: a = " << a << ", b = " << b << ", c = " << c << std::endl;
    std::cout << "      Discriminant: " << D << std::endl;
  }
  
  // Step 3: Check for real roots
  if (D < 0) {
    if (verbose) {
      std::cout << "      No real roots - discriminant is negative" << std::endl;
      std::cout << "      Returning invalid value (-1.0)" << std::endl;
    }
    return -1.0f;
  }
  
  // Step 4: Calculate roots using quadratic formula
  const float root1 = (-b + sqrt(D)) / (2 * a);
  const float root2 = (-b - sqrt(D)) / (2 * a);
  
  if (verbose) {
    std::cout << "      Root 1: " << root1 << std::endl;
    std::cout << "      Root 2: " << root2 << std::endl;
  }
  
  // Step 5: Select larger positive root
  if (root1 > 0 || root2 > 0) {
    const float max_dist = std::max(root1, root2);
    if (verbose) {
      std::cout << "      Selected maximum distance: " << max_dist << std::endl;
    }
    return max_dist;
  }
  
  if (verbose) {
    std::cout << "      Both roots are negative - returning invalid value (-1.0)" << std::endl;
  }
  return -1.0f;
}

Box
AggloCluster::createFeasibleBox(int steiner_x, int steiner_y, float max_dist) const
{
  // Step 1: Create 4 corner points at Manhattan distance max_dist from Steiner point
  const Point p1(steiner_x - static_cast<int>(max_dist), steiner_y);
  const Point p2(steiner_x + static_cast<int>(max_dist), steiner_y);
  const Point p3(steiner_x, steiner_y - static_cast<int>(max_dist));
  const Point p4(steiner_x, steiner_y + static_cast<int>(max_dist));
  
  // Step 2: Transform to 45-degree rotated coordinates
  const Point t1 = transformCoords(p1);
  const Point t2 = transformCoords(p2);
  const Point t3 = transformCoords(p3);
  const Point t4 = transformCoords(p4);
  
  // Step 3: Find bounding box in transformed coordinates
  const int min_x = std::min({t1.get<0>(), t2.get<0>(), t3.get<0>(), t4.get<0>()});
  const int max_x = std::max({t1.get<0>(), t2.get<0>(), t3.get<0>(), t4.get<0>()});
  const int min_y = std::min({t1.get<1>(), t2.get<1>(), t3.get<1>(), t4.get<1>()});
  const int max_y = std::max({t1.get<1>(), t2.get<1>(), t3.get<1>(), t4.get<1>()});
  
  const Point box_min(min_x, min_y);
  const Point box_max(max_x, max_y);
  
  return Box(box_min, box_max);
}

// -----------------------------------------------------------------------------
// [Level 4] Helper functions for extractCapacitanceDelayPoints
// -----------------------------------------------------------------------------

bool
AggloCluster::findTimingArcModel(const sta::LibertyCell* liberty_cell,
                                 const std::string& input_pin_name,
                                 const std::string& output_pin_name,
                                 const sta::TableAxis*& capacitance_axis,
                                 sta::GateTableModel*& timing_model) const
{
  // Step 1: Search through all timing arcs in the cell
  for (sta::TimingArcSet* arc_set : liberty_cell->timingArcSets()) {
    const sta::LibertyPort* from_port = arc_set->from();
    const sta::LibertyPort* to_port = arc_set->to();
    
    if (!from_port || !to_port) {
      continue;
    }
    
    // Step 2: Check if this arc matches the requested input->output path
    if (std::string(from_port->name()) != input_pin_name || 
        std::string(to_port->name()) != output_pin_name) {
      continue;
    }

    // Step 3: Look for gate table model in matching arc set
    for (sta::TimingArc* arc : arc_set->arcs()) {
      auto* gate_model = dynamic_cast<sta::GateTableModel*>(arc->model());
      if (!gate_model) {
        continue;
      }

      // Step 4: Extract delay model and find output capacitance axis
      if (findCapacitanceAxis(gate_model, capacitance_axis)) {
        timing_model = gate_model;
        return true;  // Success
      }
    }
  }

  return false;  // No matching timing model found
}

bool
AggloCluster::findCapacitanceAxis(sta::GateTableModel* gate_model,
                                  const sta::TableAxis*& capacitance_axis) const
{
  // Step 1: Get delay model from gate model
  auto* delay_model = gate_model->delayModel();
  
  // Step 2: Search all three axes for output capacitance
  for (const sta::TableAxis* axis : {delay_model->axis1(), 
                                     delay_model->axis2(), 
                                     delay_model->axis3()}) {
    if (axis && axis->variable() == sta::TableAxisVariable::total_output_net_capacitance) {
      capacitance_axis = axis;
      return true;
    }
  }
  
  return false;
}

float
AggloCluster::calculateGateDelay(sta::GateTableModel* timing_model,
                                 const sta::Pvt* pvt_conditions,
                                 float input_slew,
                                 float output_capacitance) const
{
  // Step 1: Set up query parameters
  const bool pocv_enabled = false;  // POCV (Parametric On-Chip Variation) disabled
  sta::Slew output_slew;  // Output slew (not used but required by API)
  sta::ArcDelay arc_delay;
  
  // Step 2: Query timing model for delay at this capacitance point
  timing_model->gateDelay(pvt_conditions, 
                         input_slew, 
                         output_capacitance, 
                         pocv_enabled, 
                         arc_delay, 
                         output_slew);

  // Step 3: Convert to float and ensure non-negative
  const float delay = static_cast<float>(arc_delay);
  return std::max(0.0f, delay);
}

// -----------------------------------------------------------------------------
// [Level 4] Helper function for getInstanceInputSlews
// -----------------------------------------------------------------------------

float
AggloCluster::getTerminalSlew(odb::dbITerm* terminal,
                              sta::Graph* timing_graph,
                              const sta::MinMax* min_max) const
{
  // Step 1: Convert physical pin to timing analysis pin
  sta::Pin* sta_pin = network_->dbToSta(terminal);
  
  // Step 2: Get load vertex (input side where signal arrives)
  sta::Vertex* load_vertex = timing_graph->pinLoadVertex(sta_pin);
  
  // Step 3: Extract slew at this vertex (worst-case transition time)
  return sta_->vertexSlew(load_vertex, min_max);
}

//==============================================================================
// Phase 9 Helper Functions: Compatibility Graph Generation
//==============================================================================

// -----------------------------------------------------------------------------
// [Level 1] Top-level: Build/update compatibility graph edges
// -----------------------------------------------------------------------------

std::vector<Edge>
AggloCluster::updateEdges(const FlopCluster& cluster, bool verbose)
{
  std::vector<Edge> new_edges;
  const int cluster_id = cluster.id_;

  // Early exit: Invalid or no-merge clusters
  if (!flop_cluster_is_valid_[cluster_id] || flop_cluster_no_further_merge_[cluster_id]) {
    if (verbose) {
      const char* reason = !flop_cluster_is_valid_[cluster_id] ? "invalid" : "no further merge";
      std::cout << "\n[updateEdges] Cluster #" << cluster_id << ": Skipped (" << reason << ")" << std::endl;
    }
    return new_edges;
  }

  // Early exit: Empty feasible region
  if (boost::geometry::is_empty(cluster.feasible_region_)) {
    if (verbose) {
      std::cout << "\n[updateEdges] Cluster #" << cluster_id << ": Skipped (empty feasible region)" << std::endl;
    }
    return new_edges;
  }

  if (verbose) {
    std::cout << "\n[updateEdges] Cluster #" << cluster_id << std::endl;
    std::cout << "  Size: " << cluster.flops_.size() << "-bit" << std::endl;
    std::cout << "  Position: (" << cluster.curr_pt_.x << ", " << cluster.curr_pt_.y << ")" << std::endl;
    
    // Show feasible region
    const Box& fr = cluster.feasible_region_;
    const Point min_uv = fr.min_corner();
    const Point max_uv = fr.max_corner();
    std::cout << "  Feasible region (UV): [(" 
              << min_uv.get<0>() << ", " << min_uv.get<1>() << ") - ("
              << max_uv.get<0>() << ", " << max_uv.get<1>() << ")]" << std::endl;
  }

  // Query R-Tree for intersecting clusters
  const std::vector<FlopClusterEntry> intersecting = getIntersectedCluster(cluster, verbose);

  if (verbose) {
    std::cout << "  Intersecting clusters: " << intersecting.size() << std::endl;
    
    // Show size distribution
    std::map<int, int> size_dist;
    for (const auto& [box, idx] : intersecting) {
      if (idx != cluster_id && flop_cluster_is_valid_[idx]) {
        size_dist[flop_clusters_[idx].flops_.size()]++;
      }
    }
    
    if (!size_dist.empty()) {
      std::cout << "    Size distribution: ";
      bool first = true;
      for (const auto& [size, count] : size_dist) {
        if (!first) std::cout << ", ";
        std::cout << size << "-bit(" << count << ")";
        first = false;
      }
      std::cout << std::endl;
    }
  }

  // Track duplicates and statistics
  std::set<std::pair<int, int>> edge_pairs;
  int skipped_duplicate = 0;
  int skipped_invalid = 0;
  int skipped_incompatible = 0;
  int skipped_no_placement = 0;
  int compatible_count = 0;

  // Evaluate each intersecting cluster
  for (const auto& [box, neighbor_id] : intersecting) {
    if (neighbor_id == cluster_id) continue;  // Skip self

    // Check for duplicate edge
    const auto edge_pair = std::minmax(cluster_id, neighbor_id);
    if (!edge_pairs.insert(edge_pair).second) {
      skipped_duplicate++;
      continue;
    }

    // Check validity
    if (!flop_cluster_is_valid_[neighbor_id] || flop_cluster_no_further_merge_[neighbor_id]) {
      skipped_invalid++;
      continue;
    }

    const FlopCluster& neighbor = flop_clusters_[neighbor_id];

    // Check compatibility
    const bool show_compat_debug = verbose && (compatible_count < 3);
    if (!checkCompatibility(cluster, neighbor, show_compat_debug)) {
      skipped_incompatible++;
      continue;
    }

    compatible_count++;

    // Calculate placement
    const bool show_place_debug = verbose && (new_edges.size() < 3);
    const PlacementCandidate pc = calcPlacementCandidate(cluster, neighbor, show_place_debug);

    if (!pc.is_valid_) {
      skipped_no_placement++;
      continue;
    }

    // Create edge
    new_edges.emplace_back(cluster_id, neighbor_id, pc.merge_gain_, pc.placement_pos_);

    if (verbose && new_edges.size() <= 3) {
      const int merged_size = cluster.flops_.size() + neighbor.flops_.size();
      std::cout << "    Edge #" << new_edges.size() << " with Cluster[" << neighbor_id << "]" << std::endl;
      std::cout << "      Weight: " << std::fixed << std::setprecision(2) << pc.merge_gain_
                << (pc.merge_gain_ < 0 ? " (improves)" : " (degrades)") << std::endl;
      std::cout << "      Sizes: " << cluster.flops_.size() << " + " << neighbor.flops_.size()
                << " → " << merged_size << " bits" << std::endl;
      std::cout << "      Position: (" << pc.placement_pos_.x << ", " << pc.placement_pos_.y << ")" << std::endl;
    }
  }

  if (verbose) {
    std::cout << "  Summary:" << std::endl;
    std::cout << "    Candidates evaluated: " << intersecting.size() << std::endl;
    std::cout << "    Skipped (duplicate): " << skipped_duplicate << std::endl;
    std::cout << "    Skipped (invalid): " << skipped_invalid << std::endl;
    std::cout << "    Compatible: " << compatible_count << std::endl;
    std::cout << "    Skipped (incompatible): " << skipped_incompatible << std::endl;
    std::cout << "    Skipped (no placement): " << skipped_no_placement << std::endl;
    std::cout << "    Edges created: " << new_edges.size() << std::endl;
    if (new_edges.size() > 3) {
      std::cout << "    (showing first 3 only)" << std::endl;
    }
  }
  
  return new_edges;
}

// -----------------------------------------------------------------------------
// [Level 2] Main logic for edge creation
// -----------------------------------------------------------------------------

std::vector<FlopClusterEntry>
AggloCluster::getIntersectedCluster(const FlopCluster& cluster, bool verbose) const
{
  std::vector<FlopClusterEntry> results;

  if (boost::geometry::is_empty(cluster.feasible_region_)) {
    if (verbose) {
      std::cout << "      [getIntersectedCluster] Empty feasible region" << std::endl;
    }
    return results;
  }

  feasible_regions_.query(bgi::intersects(cluster.feasible_region_),
                          std::back_inserter(results));

  if (verbose) {
    std::cout << "      [getIntersectedCluster] Found " << results.size() << " clusters" << std::endl;
    if (results.size() > 0 && results.size() <= 5) {
      std::cout << "        IDs: ";
      for (size_t i = 0; i < results.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << results[i].second;
      }
      std::cout << std::endl;
    }
  }

  return results;
}

bool 
AggloCluster::checkCompatibility(const FlopCluster& c1,
                                 const FlopCluster& c2,
                                 bool verbose) const
{
  if (verbose) {
    std::cout << "      [checkCompatibility] Clusters " << c1.id_ << " and " << c2.id_ << std::endl;
  }

  // Check mask compatibility
  const bool master_match = (c1.master_mask_ == c2.master_mask_);
  const bool inst_match = (c1.inst_mask_ == c2.inst_mask_);
  
  if (verbose) {
    std::cout << "        Master mask: " << (master_match ? "✓" : "✗") << std::endl;
    std::cout << "        Inst mask: " << (inst_match ? "✓" : "✗") << std::endl;
  }

  if (!master_match || !inst_match) {
    if (verbose) {
      std::cout << "        Result: INCOMPATIBLE (mask mismatch)" << std::endl;
    }
    return false;
  }

  // Check merged size support
  const size_t merged_size = c1.flops_.size() + c2.flops_.size();
  
  if (verbose) {
    std::cout << "        Merged size: " << c1.flops_.size() << " + " << c2.flops_.size() 
              << " = " << merged_size << " bits" << std::endl;
  }

  const auto it = compatible_masters_.find(c1.master_mask_);
  if (it == compatible_masters_.end()) {
    if (verbose) {
      std::cout << "        Result: INCOMPATIBLE (mask not in library)" << std::endl;
    }
    return false;
  }

  const bool has_master = (it->second.count(merged_size) > 0);
  
  if (verbose) {
    std::cout << "        Master for " << merged_size << "-bit: " << (has_master ? "✓" : "✗") << std::endl;
    if (!has_master) {
      std::cout << "        Available: ";
      for (const auto& [bits, masters] : it->second) {
        std::cout << bits << " ";
      }
      std::cout << "bits" << std::endl;
    }
    std::cout << "        Result: " << (has_master ? "COMPATIBLE" : "INCOMPATIBLE") << std::endl;
  }
  
  return has_master;
}

odb::dbMaster*
AggloCluster::getClusterMaster(const FlopCluster& cluster) const
{
  const int bits = static_cast<int>(cluster.flops_.size());
  if (bits <= 0) {
    return nullptr;
  }

  if (bits == 1) {
    const int flop_idx = *cluster.flops_.begin();
    if (flop_idx < 0 || flop_idx >= static_cast<int>(flop_units_.size())) {
      return nullptr;
    }
    return flop_units_[flop_idx].inst_->getMaster();
  }

  const auto mask_it = representative_masters_.find(cluster.master_mask_);
  if (mask_it == representative_masters_.end()) {
    return nullptr;
  }

  const auto bit_it = mask_it->second.find(bits);
  if (bit_it == mask_it->second.end()) {
    return nullptr;
  }

  return bit_it->second;
}

PlacementCandidate
AggloCluster::calcPlacementCandidate(const FlopCluster& c1, const FlopCluster& c2, bool verbose)
{
  if (verbose) {
    std::cout << "      [calcPlacementCandidate] Clusters " << c1.id_ << " and " << c2.id_ << std::endl;
  }

  // Get masters
  odb::dbMaster* master1 = getClusterMaster(c1);
  odb::dbMaster* master2 = getClusterMaster(c2);
  
  const int merged_bits = c1.flops_.size() + c2.flops_.size();
  odb::dbMaster* merged_master = nullptr;
  const auto mask_it = representative_masters_.find(c1.master_mask_);
  if (mask_it != representative_masters_.end()) {
    const auto bit_it = mask_it->second.find(merged_bits);
    if (bit_it != mask_it->second.end()) {
      merged_master = bit_it->second;
    }
  }

  if (verbose) {
    std::cout << "        Masters:" << std::endl;
    std::cout << "          C1 (" << c1.flops_.size() << "-bit): " 
              << (master1 ? master1->getName() : "NULL") << std::endl;
    std::cout << "          C2 (" << c2.flops_.size() << "-bit): " 
              << (master2 ? master2->getName() : "NULL") << std::endl;
    std::cout << "          Merged (" << merged_bits << "-bit): " 
              << (merged_master ? merged_master->getName() : "NULL") << std::endl;
  }

  if (!master1 || !master2 || !merged_master) {
    if (verbose) {
      std::cout << "        Result: INVALID (missing master)" << std::endl;
    }
    return PlacementCandidate(c1.curr_pt_, 0.0, false);
  }

  // Calculate HPWL-optimal box
  const Box hpwl_box = calcMedianBox(c1, c2);

  if (verbose) {
    std::cout << "        HPWL box: [(" << hpwl_box.min_corner().get<0>() << ", "
              << hpwl_box.min_corner().get<1>() << ") - (" << hpwl_box.max_corner().get<0>()
              << ", " << hpwl_box.max_corner().get<1>() << ")]" << std::endl;
  }

  // Get feasible region intersection
  const Box feasible_uv = getFeasibleRegionIntersection(c1, c2);
  if (boost::geometry::is_empty(feasible_uv)) {
    if (verbose) {
      std::cout << "        Result: INVALID (no feasible region overlap)" << std::endl;
    }
    return PlacementCandidate(c1.curr_pt_, 0.0, false);
  }

  if (verbose) {
    const Point min_uv = feasible_uv.min_corner();
    const Point max_uv = feasible_uv.max_corner();
    std::cout << "        Feasible (UV): [(" << min_uv.get<0>() << ", " << min_uv.get<1>()
              << ") - (" << max_uv.get<0>() << ", " << max_uv.get<1>() << ")]" << std::endl;
  }

  // Project HPWL-optimal point onto feasible region
  const Point projected = project(hpwl_box, feasible_uv);

  if (verbose) {
    std::cout << "        Projected optimal: (" << projected.get<0>() << ", " 
              << projected.get<1>() << ")" << std::endl;
  }

  // Generate candidate positions
  std::vector<Point> candidates = generateUniformSamples(feasible_uv, 16);
  candidates.push_back(projected);

  if (verbose) {
    std::cout << "        Candidates: " << candidates.size() << std::endl;
  }

  // Sort by Manhattan distance to projected point
  auto manhattan = [](const Point& a, const Point& b) {
    const int64_t dx = std::abs(static_cast<int64_t>(a.get<0>()) - static_cast<int64_t>(b.get<0>()));
    const int64_t dy = std::abs(static_cast<int64_t>(a.get<1>()) - static_cast<int64_t>(b.get<1>()));
    return dx + dy;
  };

  std::sort(candidates.begin(), candidates.end(),
            [&](const Point& a, const Point& b) {
              return manhattan(a, projected) < manhattan(b, projected);
            });

  // Evaluate candidates
  for (size_t i = 0; i < candidates.size(); ++i) {
    const Point& cand = candidates[i];
    
    if (verbose && i < 3) {
      std::cout << "          #" << i << ": (" << cand.get<0>() << ", " << cand.get<1>()
                << "), dist=" << manhattan(cand, projected) << std::endl;
    }

    // Check density constraint
    if (checkPlacementDensityConstraint(c1, master1, c2, master2, cand, merged_master)) {
      const FloatPoint pos(static_cast<float>(cand.get<0>()), static_cast<float>(cand.get<1>()));
      
      if (verbose && i < 3) {
        std::cout << "            Density: ✓" << std::endl;
      }
      
      const double gain = calcMergeHPWLGain(c1, c2, pos, verbose);

      if (verbose) {
        std::cout << "        Result: VALID (candidate #" << i << ")" << std::endl;
      }

      return PlacementCandidate(pos, gain, true);
    } else if (verbose && i < 3) {
      std::cout << "            Density: ✗" << std::endl;
    }
  }

  // All failed
  if (verbose) {
    std::cout << "        Result: INVALID (all " << candidates.size() 
              << " candidates failed density)" << std::endl;
  }
  return PlacementCandidate(c1.curr_pt_, 0.0, false);
}

// -----------------------------------------------------------------------------
// [Level 3] Helpers for placement candidate calculation
// -----------------------------------------------------------------------------

Box 
AggloCluster::calcMedianBox(const FlopCluster& c1, const FlopCluster& c2) const
{
  // Collect all external pins connected to both clusters
  std::set<std::variant<odb::dbITerm*, odb::dbBTerm*>> connected_pins = getConnectedPins(c1);
  std::set<std::variant<odb::dbITerm*, odb::dbBTerm*>> pins_from_c2 = getConnectedPins(c2);
  connected_pins.insert(pins_from_c2.begin(), pins_from_c2.end());

  // Extract pin coordinates
  std::vector<int> x_coords;
  std::vector<int> y_coords;

  for (const auto& pin_variant : connected_pins) {
    std::visit(
        [&](auto&& pin) {
          using T = std::decay_t<decltype(pin)>;

          if constexpr (std::is_same_v<T, odb::dbITerm*>) {
            int pin_x, pin_y;
            if (pin->getAvgXY(&pin_x, &pin_y)) {
              x_coords.push_back(pin_x);
              y_coords.push_back(pin_y);
            }
          } else if constexpr (std::is_same_v<T, odb::dbBTerm*>) {
            odb::Rect bbox = pin->getBBox();
            if (!bbox.isInverted()) {
              x_coords.push_back(bbox.xCenter());
              y_coords.push_back(bbox.yCenter());
            }
          }
        },
        pin_variant);
  }

  // Calculate median box (lower and upper medians)
  size_t n = x_coords.size();
  std::sort(x_coords.begin(), x_coords.end());
  std::sort(y_coords.begin(), y_coords.end());

  int min_x = x_coords[(n - 1) / 2];
  int max_x = x_coords[n / 2];
  int min_y = y_coords[(n - 1) / 2];
  int max_y = y_coords[n / 2];

  return Box(Point(min_x, min_y), Point(max_x, max_y));
}

Box
AggloCluster::getFeasibleRegionIntersection(const FlopCluster& c1,
                                            const FlopCluster& c2) const
{
  Box intersection_box;
  boost::geometry::intersection(c1.feasible_region_,
                                c2.feasible_region_,
                                intersection_box);

  return intersection_box;
}

Point 
AggloCluster::project(const Box& hpwl_box, const Box& feasible_box) const
{
  // Calculate center of HPWL box in XY coordinates
  const Point& min_xy = hpwl_box.min_corner();
  const Point& max_xy = hpwl_box.max_corner();

  int target_x = (min_xy.get<0>() + max_xy.get<0>()) / 2;
  int target_y = (min_xy.get<1>() + max_xy.get<1>()) / 2;
  Point target_xy(target_x, target_y);

  // Transform to UV coordinates
  Point target_uv = transformCoords(target_xy);

  // Clamp to feasible region bounds
  const Point& min_uv = feasible_box.min_corner();
  const Point& max_uv = feasible_box.max_corner();

  int u_final = std::clamp(target_uv.get<0>(), min_uv.get<0>(), max_uv.get<0>());
  int v_final = std::clamp(target_uv.get<1>(), min_uv.get<1>(), max_uv.get<1>());

  Point final_uv(u_final, v_final);

  // Transform back to XY coordinates
  Point final_xy = inverseTransformCoords(final_uv);

  return final_xy;
}

std::vector<gpl::Point> 
gpl::AggloCluster::generateUniformSamples(const Box& box, int p) const
{
  std::vector<Point> samples;
  if (p <= 0) {
    return samples;
  }

  int u_min = box.min_corner().get<0>();
  int u_max = box.max_corner().get<0>();
  int v_min = box.min_corner().get<1>();
  int v_max = box.max_corner().get<1>();

  // Special case: single sample at center
  if (p == 1) {
    int u = (u_min + u_max) / 2;
    int v = (v_min + v_max) / 2;
    Point center_p_uv(u, v);
    Point center_p_xy = inverseTransformCoords(center_p_uv);
    samples.push_back(center_p_xy);
    return samples;
  }

  // Multiple samples: create grid including boundaries
  int num_steps = static_cast<int>(std::ceil(std::sqrt(p)));
  
  double u_step = (u_max == u_min) ? 0.0 : static_cast<double>(u_max - u_min) / (num_steps - 1.0);
  double v_step = (v_max == v_min) ? 0.0 : static_cast<double>(v_max - v_min) / (num_steps - 1.0);

  for (int i = 0; i < num_steps; ++i) {
    for (int j = 0; j < num_steps; ++j) {
      
      // Ensure exact boundaries at last step
      int u = (i == num_steps - 1) ? u_max : (u_min + static_cast<int>(std::round(i * u_step)));
      int v = (j == num_steps - 1) ? v_max : (v_min + static_cast<int>(std::round(j * v_step)));
      
      Point sample_p_uv(u, v); 
      Point sample_p_xy = inverseTransformCoords(sample_p_uv);
      samples.push_back(sample_p_xy);
    }
  }
  
  return samples;
}

bool 
gpl::AggloCluster::checkPlacementDensityConstraint(
    const FlopCluster& c1,
    odb::dbMaster* master1, 
    const FlopCluster& c2,
    odb::dbMaster* master2, 
    const Point& new_pos,
    odb::dbMaster* new_master)     
{
  if (!master1 || !master2 || !new_master) {
    if (verbose_) {
      std::cout << "[checkPlacementDensityConstraint] Null master. Skipping check." << std::endl;
    }
    return false;
  }

  // Convert cluster positions to integer points
  Point c1_pos(std::lround(c1.curr_pt_.x), std::lround(c1.curr_pt_.y));
  Point c2_pos(std::lround(c2.curr_pt_.x), std::lround(c2.curr_pt_.y));

  const bool has_overflow = virtual_bin_grid_.wouldOverflow(master1, c1_pos,
                                                            master2, c2_pos,
                                                            new_master, new_pos);

  if (verbose_ && has_overflow) {
    std::cout << "[checkPlacementDensityConstraint] candidate ✗ violates density" << std::endl;
  }

  return !has_overflow; 
}

double 
AggloCluster::calcMergeHPWLGain(const FlopCluster& c1, 
                                const FlopCluster& c2, 
                                const FloatPoint& merge_pos,
                                bool verbose) const
{
  if (verbose) {
    std::cout << "      [calcMergeHPWLGain]" << std::endl;
    std::cout << "        C1: (" << static_cast<int>(c1.curr_pt_.x) << ", "
              << static_cast<int>(c1.curr_pt_.y) << "), " << c1.flops_.size() << " bits" << std::endl;
    std::cout << "        C2: (" << static_cast<int>(c2.curr_pt_.x) << ", "
              << static_cast<int>(c2.curr_pt_.y) << "), " << c2.flops_.size() << " bits" << std::endl;
    std::cout << "        Merge: (" << static_cast<int>(merge_pos.x) << ", "
              << static_cast<int>(merge_pos.y) << ")" << std::endl;
  }
  
  // Build instance sets for quick lookup
  std::set<odb::dbInst*> c1_insts, c2_insts;
  for (int idx : c1.flops_) c1_insts.insert(flop_units_[idx].inst_);
  for (int idx : c2.flops_) c2_insts.insert(flop_units_[idx].inst_);
  
  // Collect unique nets
  std::set<odb::dbNet*> nets;
  auto collect_nets = [&](const FlopCluster& c) {
    for (int idx : c.flops_) {
      for (odb::dbITerm* iterm : flop_units_[idx].inst_->getITerms()) {
        odb::dbNet* net = iterm->getNet();
        if (net && (net->getSigType() == odb::dbSigType::SIGNAL || 
                    net->getSigType() == odb::dbSigType::CLOCK)) {
          nets.insert(net);
        }
      }
    }
  };
  collect_nets(c1);
  collect_nets(c2);
  
  if (verbose) {
    std::cout << "        Nets: " << nets.size() << std::endl;
  }
  
  // Positions
  const int c1_x = static_cast<int>(c1.curr_pt_.x);
  const int c1_y = static_cast<int>(c1.curr_pt_.y);
  const int c2_x = static_cast<int>(c2.curr_pt_.x);
  const int c2_y = static_cast<int>(c2.curr_pt_.y);
  const int merge_x = static_cast<int>(merge_pos.x);
  const int merge_y = static_cast<int>(merge_pos.y);
  
  double orig_hpwl = 0.0;
  double new_hpwl = 0.0;
  int shown = 0;
  
  // Process each net
  for (odb::dbNet* net : nets) {
    bool has_c1 = false;
    bool has_c2 = false;
    std::vector<std::pair<int, int>> external;
    
    // Check connections
    for (odb::dbITerm* iterm : net->getITerms()) {
      odb::dbInst* inst = iterm->getInst();
      if (c1_insts.count(inst)) {
        has_c1 = true;
      } else if (c2_insts.count(inst)) {
        has_c2 = true;
      } else {
        const FloatPoint coord = getPinCoordinate(iterm);
        external.emplace_back(static_cast<int>(coord.x), static_cast<int>(coord.y));
      }
    }
    
    // Include top-level ports
    for (odb::dbBTerm* bterm : net->getBTerms()) {
      const FloatPoint coord = getPinCoordinate(bterm);
      external.emplace_back(static_cast<int>(coord.x), static_cast<int>(coord.y));
    }
    
    // Calculate original HPWL
    int min_x = std::numeric_limits<int>::max();
    int max_x = std::numeric_limits<int>::min();
    int min_y = std::numeric_limits<int>::max();
    int max_y = std::numeric_limits<int>::min();
    
    for (const auto& [x, y] : external) {
      min_x = std::min(min_x, x);
      max_x = std::max(max_x, x);
      min_y = std::min(min_y, y);
      max_y = std::max(max_y, y);
    }
    
    if (has_c1) {
      min_x = std::min(min_x, c1_x);
      max_x = std::max(max_x, c1_x);
      min_y = std::min(min_y, c1_y);
      max_y = std::max(max_y, c1_y);
    }
    
    if (has_c2) {
      min_x = std::min(min_x, c2_x);
      max_x = std::max(max_x, c2_x);
      min_y = std::min(min_y, c2_y);
      max_y = std::max(max_y, c2_y);
    }
    
    const double net_orig = (max_x >= min_x && max_y >= min_y) 
                            ? (max_x - min_x) + (max_y - min_y) 
                            : 0.0;
    orig_hpwl += net_orig;
    
    // Calculate merged HPWL
    min_x = std::numeric_limits<int>::max();
    max_x = std::numeric_limits<int>::min();
    min_y = std::numeric_limits<int>::max();
    max_y = std::numeric_limits<int>::min();
    
    for (const auto& [x, y] : external) {
      min_x = std::min(min_x, x);
      max_x = std::max(max_x, x);
      min_y = std::min(min_y, y);
      max_y = std::max(max_y, y);
    }
    
    if (has_c1 || has_c2) {
      min_x = std::min(min_x, merge_x);
      max_x = std::max(max_x, merge_x);
      min_y = std::min(min_y, merge_y);
      max_y = std::max(max_y, merge_y);
    }
    
    const double net_new = (max_x >= min_x && max_y >= min_y) 
                           ? (max_x - min_x) + (max_y - min_y) 
                           : 0.0;
    new_hpwl += net_new;
    
    if (verbose && shown < 3) {
      std::cout << "        Net '" << net->getName() << "':" << std::endl;
      std::cout << "          C1: " << (has_c1 ? "✓" : "✗")
                << ", C2: " << (has_c2 ? "✓" : "✗")
                << ", External: " << external.size() << std::endl;
      std::cout << "          Original: " << net_orig << " → Merged: " << net_new
                << " (Δ " << (net_new - net_orig) << ")" << std::endl;
      shown++;
    }
  }
  
  const double gain = new_hpwl - orig_hpwl;
  
  if (verbose) {
    std::cout << "        Total HPWL:" << std::endl;
    std::cout << "          Original: " << orig_hpwl << std::endl;
    std::cout << "          Merged: " << new_hpwl << std::endl;
    std::cout << "          Gain: " << gain << (gain < 0 ? " (improves)" : " (degrades)") << std::endl;
  }
  
  return gain;
}

// -----------------------------------------------------------------------------
// [Level 4] Lower-level helpers
// -----------------------------------------------------------------------------

std::set<std::variant<odb::dbITerm*, odb::dbBTerm*>>
AggloCluster::getConnectedPins(const FlopCluster& cluster) const
{
  std::set<std::variant<odb::dbITerm*, odb::dbBTerm*>> connected_pins;
  
  // Build set of instances in cluster for fast lookup
  std::set<odb::dbInst*> cluster_insts;
  for (int flop_idx : cluster.flops_) {
    cluster_insts.insert(flop_units_[flop_idx].inst_);
  }
  
  // Collect external pins connected to cluster
  for (int flop_idx : cluster.flops_) {
    const FlopUnit& flop = flop_units_[flop_idx];
    
    for (odb::dbITerm* iterm : flop.inst_->getITerms()) {
      odb::dbNet* net = iterm->getNet();
      
      // Skip invalid or non-signal/clock nets
      if (net == nullptr
          || (net->getSigType() != odb::dbSigType::SIGNAL
              && net->getSigType() != odb::dbSigType::CLOCK)) {
        continue;
      }
      
      // Add external instance terminals
      for (odb::dbITerm* net_iterm : net->getITerms()) {
        if (cluster_insts.find(net_iterm->getInst()) == cluster_insts.end()) {
          connected_pins.insert(net_iterm);
        }
      }
      
      // Add top-level block terminals
      for (odb::dbBTerm* net_bterm : net->getBTerms()) {
        connected_pins.insert(net_bterm);
      }
    }
  }
  
  return connected_pins;
}

FloatPoint 
AggloCluster::getPinCoordinate(const std::variant<odb::dbITerm*, odb::dbBTerm*>& pin_variant) const
{
  int x = 0, y = 0;

  if (std::holds_alternative<odb::dbITerm*>(pin_variant)) {
    odb::dbITerm* iterm = std::get<odb::dbITerm*>(pin_variant);
    
    if (!iterm->getAvgXY(&x, &y)) {
      if (verbose_) {
        std::cout << "[WARNING] Failed to get coordinates for ITerm: " 
                  << iterm->getName() << ", using (0, 0)" << std::endl;
      }
      return FloatPoint(0.0f, 0.0f);
    }
    return FloatPoint(static_cast<float>(x), static_cast<float>(y));
  } 
  
  odb::dbBTerm* bterm = std::get<odb::dbBTerm*>(pin_variant);
  odb::Rect bbox = bterm->getBBox();
  
  if (bbox.isInverted()) {
    if (verbose_) {
      std::cout << "[WARNING] Invalid BBox for BTerm: " 
                << bterm->getName() << ", using (0, 0)" << std::endl;
    }
    return FloatPoint(0.0f, 0.0f);
  }
  
  return FloatPoint(static_cast<float>(bbox.xCenter()),
                    static_cast<float>(bbox.yCenter()));
}

//==============================================================================
// Phase 10 Helper Functions: Agglomerative Clustering
//==============================================================================

//------------------------------------------------------------------------------
// [Level 1] Top-level: Execute merge operations and maintain graph
//------------------------------------------------------------------------------

int 
AggloCluster::mergeClusters(const Edge& edge, bool verbose)
{
  const FlopCluster& c1 = flop_clusters_[edge.n1];
  const FlopCluster& c2 = flop_clusters_[edge.n2];

  if (verbose) {
    std::cout << "\n      [mergeClusters] Merging Cluster[" << edge.n1 << "] + Cluster[" << edge.n2 << "]" << std::endl;
    std::cout << "        C1: size=" << c1.flops_.size() << ", pos=(" 
              << c1.curr_pt_.x << ", " << c1.curr_pt_.y << ")" << std::endl;
    std::cout << "        C2: size=" << c2.flops_.size() << ", pos=(" 
              << c2.curr_pt_.x << ", " << c2.curr_pt_.y << ")" << std::endl;
    std::cout << "        New position: (" << edge.pos.x << ", " << edge.pos.y << ")" << std::endl;
  }

  // Invalidate old clusters
  flop_cluster_is_valid_[edge.n1] = false;
  flop_cluster_is_valid_[edge.n2] = false;

  // Remove old edges and feasible regions
  removeEdges(c1, verbose);
  removeEdges(c2, verbose);
  feasible_regions_.remove(std::make_pair(c1.feasible_region_, edge.n1));
  feasible_regions_.remove(std::make_pair(c2.feasible_region_, edge.n2));

  // Create new merged cluster
  const int new_cluster_idx = flop_clusters_.size();
  flop_clusters_.emplace_back(new_cluster_idx, c1, c2, edge);
  flop_cluster_is_valid_.push_back(true);
  flop_cluster_no_further_merge_.push_back(false);
  
  FlopCluster& new_cluster = flop_clusters_.back();

  if (verbose) {
    std::cout << "        ✓ Created Cluster[" << new_cluster_idx << "]: size=" << new_cluster.flops_.size() << std::endl;
    
    if (boost::geometry::is_empty(new_cluster.feasible_region_)) {
      std::cout << "          Feasible region: EMPTY ⚠ WARNING!" << std::endl;
    } else {
      const Box& fr = new_cluster.feasible_region_;
      std::cout << "          Feasible region: (" << fr.min_corner().get<0>() << ", " 
                << fr.min_corner().get<1>() << ") -> (" << fr.max_corner().get<0>() << ", " 
                << fr.max_corner().get<1>() << ")" << std::endl;
    }
  }

  // Update flop unit assignments
  for (int flop_idx : new_cluster.flops_) {
    flop_units_[flop_idx].cluster_idx_ = new_cluster_idx;
    flop_units_[flop_idx].curr_pt_ = new_cluster.curr_pt_;
  }

  // Insert new feasible region into R-Tree
  if (!boost::geometry::is_empty(new_cluster.feasible_region_)) {
    feasible_regions_.insert(std::make_pair(new_cluster.feasible_region_, new_cluster_idx));
    
    if (verbose) {
      std::cout << "          R-Tree: Inserted successfully" << std::endl;
    }
  } else if (verbose) {
    std::cout << "          R-Tree: Skipped (empty region - won't be discoverable!)" << std::endl;
  }

  odb::dbMaster* master1 = getClusterMaster(c1);
  odb::dbMaster* master2 = getClusterMaster(c2);
  odb::dbMaster* merged_master = getClusterMaster(new_cluster);

  if (master1 && master2 && merged_master) {
    const Point c1_pos(std::lround(c1.curr_pt_.x), std::lround(c1.curr_pt_.y));
    const Point c2_pos(std::lround(c2.curr_pt_.x), std::lround(c2.curr_pt_.y));
    const Point merged_pos(std::lround(new_cluster.curr_pt_.x), std::lround(new_cluster.curr_pt_.y));

    virtual_bin_grid_.applyMerge(master1, c1_pos, master2, c2_pos, merged_master, merged_pos);

    if (verbose && virtual_bin_grid_.checkOverflow()) {
      std::cout << "          WARNING: Density overflow detected after merge (unexpected)" << std::endl;
    }
  } else if (verbose) {
    std::cout << "          WARNING: Missing master information for density update" << std::endl;
  }

  return new_cluster_idx;
}

bool
AggloCluster::isFurtherMergeable(const FlopCluster& cluster, bool verbose) const
{
  const auto it = compatible_masters_.find(cluster.master_mask_);
  if (it == compatible_masters_.end()) {
    if (verbose) {
      std::cout << "    ✗ Mask not found in compatible_masters_ (possible bug)" << std::endl;
    }
    return false;
  }

  const auto& bits_to_masters = it->second;
  if (bits_to_masters.empty()) {
    if (verbose) {
      std::cout << "    ✗ No masters available for this mask" << std::endl;
    }
    return false;
  }

  const int max_bits = bits_to_masters.rbegin()->first;
  const int current_size = static_cast<int>(cluster.flops_.size());
  const bool can_merge = current_size < max_bits;
  
  if (verbose) {
    std::cout << "    isFurtherMergeable: " << (can_merge ? "✓ YES" : "✗ NO") 
              << " (size=" << current_size << "/" << max_bits << ")" << std::endl;
  }
  
  return can_merge;
}  

std::set<int>
AggloCluster::distributeSlack(const FlopCluster& cluster, bool verbose)
{
  if (verbose) {
    std::cout << "        [distributeSlack] Cluster[" << cluster.id_ 
              << "]: " << cluster.flops_.size() << " flop(s)" << std::endl;
  }

  std::set<int> affected_flop_units;
  int total_paths = 0;
  int internal_paths = 0;
  int redistributed = 0;
  int shown_details = 0;
  constexpr int MAX_DETAILS = 5;

  for (int flop_idx : cluster.flops_) {
    const FlopUnit& flop = flop_units_[flop_idx];
    
    if (verbose && shown_details < MAX_DETAILS) {
      std::cout << "          Flop[" << flop_idx << "]: " << flop.inst_->getName() << std::endl;
      const float dx = std::abs(flop.curr_pt_.x - flop.orig_pt_.x);
      const float dy = std::abs(flop.curr_pt_.y - flop.orig_pt_.y);
      std::cout << "            Movement: " << (dx + dy) << " DBU" << std::endl;
    }

    const auto used_slacks = calcUsedSlacks(flop, verbose && shown_details < MAX_DETAILS);

    for (const auto& [pin, used_slack_list] : used_slacks) {
      for (const auto& [path_idx, used_slack] : used_slack_list) {
        total_paths++;
        const TimingPath& path = timing_paths_[path_idx];

        const int other_flop_idx = (path.start_flop_idx_ == flop_idx) 
                                     ? path.end_flop_idx_ 
                                     : path.start_flop_idx_;

        if (flop_units_[other_flop_idx].cluster_idx_ == cluster.id_) {
          internal_paths++;
          continue;
        }

        FlopUnit& other_flop = flop_units_[other_flop_idx];
        affected_flop_units.insert(other_flop_idx);

        odb::dbITerm* other_pin = (path.start_flop_idx_ == flop_idx) 
                                    ? path.end_pin_ 
                                    : path.start_pin_;
        
        if (!other_pin || other_flop.pin_budgets_.find(other_pin) == other_flop.pin_budgets_.end()) {
          continue;
        }

        sta::Slack original_budget = 0.0;
        const auto& current_budgets = flop.pin_budgets_.at(pin);
        for (const auto& [budget_path_idx, budget_slack] : current_budgets) {
          if (budget_path_idx == path_idx) {
            original_budget = budget_slack;
            break;
          }
        }

        auto& other_budgets = other_flop.pin_budgets_.at(other_pin);
        
        for (auto& [budget_path_idx, budget_slack] : other_budgets) {
          if (budget_path_idx == path_idx) {
            const sta::Slack old_slack = budget_slack;
            const sta::Slack slack_delta = original_budget - used_slack;
            budget_slack += slack_delta;
            
            if (verbose && shown_details < MAX_DETAILS) {
              std::cout << "            Path[" << path_idx << "]: Flop[" << flop_idx 
                        << "] → Flop[" << other_flop_idx << "]" << std::endl;
              std::cout << "              Unused slack: " << std::fixed << std::setprecision(12) 
                        << slack_delta << std::endl;
              std::cout << "              Other slack: " << old_slack << " → " << budget_slack 
                        << (slack_delta > 0 ? " ✓" : (slack_delta < 0 ? " ⚠" : "")) << std::endl;
              shown_details++;
            }
            
            redistributed++;
            break;
          }
        }

        std::sort(other_budgets.begin(), other_budgets.end(),
                  [](const auto& a, const auto& b) { return a.second < b.second; });
      }
    }
  }
  
  if (verbose) {
    std::cout << "          Summary: paths=" << total_paths 
              << ", internal=" << internal_paths 
              << ", redistributed=" << redistributed 
              << ", affected=" << affected_flop_units.size() << std::endl;
    if (shown_details >= MAX_DETAILS) {
      std::cout << "          (First " << MAX_DETAILS << " redistributions shown)" << std::endl;
    }
  }

  return affected_flop_units;
}

void 
AggloCluster::removeEdges(const FlopCluster& cluster, bool verbose)
{
  const int cluster_idx = cluster.id_;

  if (adj_list_.find(cluster_idx) == adj_list_.end()) {
    if (verbose) {
      std::cout << "          removeEdges(Cluster[" << cluster_idx << "]): No edges" << std::endl;
    }
    return;
  }

  const auto& edges_to_remove = adj_list_.at(cluster_idx);
  const int edge_count = edges_to_remove.size();
  
  if (verbose) {
    std::cout << "          removeEdges(Cluster[" << cluster_idx << "]): " 
              << edge_count << " edge(s)" << std::endl;
  }

  for (const Edge& edge : edges_to_remove) {
    edge_pq_.erase(edge);

    const int neighbor_idx = (edge.n1 == cluster_idx) ? edge.n2 : edge.n1;
    if (adj_list_.count(neighbor_idx)) {
      adj_list_.at(neighbor_idx).erase(edge);
    }
  }

  adj_list_.erase(cluster_idx);
}

void 
AggloCluster::updateFeasibleRegion(FlopCluster& cluster, bool verbose)
{
  if (verbose) {
    std::cout << "          [updateFeasibleRegion] Cluster[" << cluster.id_ 
              << "]: " << cluster.flops_.size() << " flop(s)" << std::endl;
  }

  const Box old_feasible_region = cluster.feasible_region_;
  const bool old_was_empty = boost::geometry::is_empty(old_feasible_region);

  if (cluster.flops_.empty()) {
    cluster.feasible_region_ = Box();
    if (verbose) {
      std::cout << "            WARNING: Cluster has no flops!" << std::endl;
    }
  } else {
    auto it = cluster.flops_.begin();
    cluster.feasible_region_ = flop_units_[*it].feasible_region_;
    
    int intersections = 0;
    for (++it; it != cluster.flops_.end(); ++it) {
      Box temp_result;
      boost::geometry::intersection(cluster.feasible_region_, 
                                    flop_units_[*it].feasible_region_, 
                                    temp_result);
      cluster.feasible_region_ = temp_result;
      intersections++;
      
      if (boost::geometry::is_empty(cluster.feasible_region_)) {
        if (verbose) {
          std::cout << "            WARNING: Feasible region became EMPTY after " 
                    << intersections << " intersection(s)" << std::endl;
        }
        break;
      }
    }
    
    if (verbose && !boost::geometry::is_empty(cluster.feasible_region_)) {
      const Box& fr = cluster.feasible_region_;
      std::cout << "            New region: (" << fr.min_corner().get<0>() << ", " 
                << fr.min_corner().get<1>() << ") -> (" << fr.max_corner().get<0>() << ", " 
                << fr.max_corner().get<1>() << ") [" << intersections << " intersections]" << std::endl;
    }
  }

  if (!old_was_empty) {
    feasible_regions_.remove(std::make_pair(old_feasible_region, cluster.id_));
  }
  
  if (!boost::geometry::is_empty(cluster.feasible_region_)) {
    feasible_regions_.insert(std::make_pair(cluster.feasible_region_, cluster.id_));
    if (verbose) {
      std::cout << "            ✓ R-Tree updated" << std::endl;
    }
  } else if (verbose) {
    std::cout << "            R-Tree: Skipped (empty region)" << std::endl;
  }
}

//==============================================================================
// Reverse Engineering: calcUsedSlacks()
//==============================================================================

//------------------------------------------------------------------------------
// [Level 1] Main Entry Point: Calculate used slacks based on flop movement
//------------------------------------------------------------------------------

std::unordered_map<odb::dbITerm*, std::vector<std::pair<int, sta::Slack>>>
AggloCluster::calcUsedSlacks(const FlopUnit& flop, bool verbose)
{
  std::unordered_map<odb::dbITerm*, std::vector<std::pair<int, sta::Slack>>> used_slacks;
  
  if (verbose) {
    std::cout << "\n[calcUsedSlacks] Flop: " << flop.inst_->getName() << std::endl;
    const float dx = std::abs(flop.curr_pt_.x - flop.orig_pt_.x);
    const float dy = std::abs(flop.curr_pt_.y - flop.orig_pt_.y);
    std::cout << "  Movement: " << (dx + dy) << " DBU" << std::endl;
  }
  
  const auto est = resizer_->getEstimateParasitics();
  const double unit_c = est->wireSignalCapacitance(corner_);
  const double unit_r = est->wireSignalResistance(corner_);
  
  std::vector<odb::dbITerm*> d_pins, q_pins, qn_pins;
  odb::dbITerm* clk_pin = nullptr;
  
  for (odb::dbITerm* iterm : flop.inst_->getITerms()) {
    if (isDPin(iterm)) {
      d_pins.push_back(iterm);
    } else if (isQPin(iterm)) {
      q_pins.push_back(iterm);
    } else if (isQNPin(iterm)) {
      qn_pins.push_back(iterm);
    } else if (isClockPin(iterm)) {
      clk_pin = iterm;
    }
  }
  
  if (!clk_pin) {
    if (verbose) {
      std::cout << "  WARNING: No clock pin - returning empty" << std::endl;
    }
    return used_slacks;
  }
  
  std::vector<odb::dbITerm*> output_pins;
  output_pins.reserve(q_pins.size() + qn_pins.size());
  output_pins.insert(output_pins.end(), q_pins.begin(), q_pins.end());
  output_pins.insert(output_pins.end(), qn_pins.begin(), qn_pins.end());
  
  if (verbose && !output_pins.empty()) {
    std::cout << "  Processing " << output_pins.size() << " output pin(s)..." << std::endl;
  }
  
  for (odb::dbITerm* out_pin : output_pins) {
    calcUsedSlacksFanOut(flop, out_pin, clk_pin->getMTerm(), 
                         est, unit_r, unit_c, used_slacks, verbose);
  }
  
  if (verbose && !d_pins.empty()) {
    std::cout << "  Processing " << d_pins.size() << " input pin(s)..." << std::endl;
  }
  
  for (odb::dbITerm* d_pin : d_pins) {
    calcUsedSlacksFanIn(flop, d_pin, est, unit_r, unit_c, used_slacks, verbose);
  }
  
  if (verbose) {
    std::cout << "[calcUsedSlacks] Completed\n" << std::endl;
  }
  
  return used_slacks;
}

//------------------------------------------------------------------------------
// [Level 2] FanOut Pin Processing (Q/QN pins)
//------------------------------------------------------------------------------

void
AggloCluster::calcUsedSlacksFanOut(
    const FlopUnit& flop,
    odb::dbITerm* out_pin,
    odb::dbMTerm* clk_pin_lib,
    est::EstimateParasitics* est,
    double unit_r,
    double unit_c,
    std::unordered_map<odb::dbITerm*, std::vector<std::pair<int, sta::Slack>>>& used_slacks,
    bool verbose)
{
  if (verbose) {
    std::cout << "\n  [FanOut: " << out_pin->getMTerm()->getName() << "]" << std::endl;
  }
  
  const auto budget_it = flop.pin_budgets_.find(out_pin);
  if (budget_it == flop.pin_budgets_.end() || budget_it->second.empty()) {
    if (verbose) {
      std::cout << "    No timing constraints - skipping" << std::endl;
    }
    return;
  }
  
  odb::dbNet* fanout_net = out_pin->getNet();
  if (!fanout_net) {
    if (verbose) {
      std::cout << "    Not connected to net - skipping" << std::endl;
    }
    return;
  }
  
  const auto cap_delay = extractCapacitanceDelayPoints(
      flop.inst_, clk_pin_lib->getName(), out_pin->getMTerm()->getName(), 0);
  
  if (cap_delay.empty() || cap_delay.size() < 2) {
    if (verbose) {
      std::cout << "    Insufficient cap-delay data - skipping" << std::endl;
    }
    return;
  }
  
  std::vector<float> coeffs;
  coeffs.reserve(cap_delay.size() - 1);
  for (size_t i = 1; i < cap_delay.size(); ++i) {
    coeffs.push_back((cap_delay[i].second - cap_delay[i-1].second) / 
                     (cap_delay[i].first - cap_delay[i-1].first));
  }
  
  const sta::Pin* driver_pin_sta = network_->dbToSta(out_pin);
  est::SteinerTree* driver_steiner_tree = est->makeSteinerTree(driver_pin_sta);
  
  const sta::Net* fanout_net_sta = network_->dbToSta(fanout_net);
  float pin_cap = 0.0f, wire_cap_total = 0.0f;
  const sta::MinMax* mm = sta::MinMax::max();
  sta_->connectedCap(fanout_net_sta, corner_, mm, pin_cap, wire_cap_total);
  const float total_net_cap = pin_cap + wire_cap_total;
  
  const auto top_steiner_point = driver_steiner_tree->top();
  const auto driver_point = driver_steiner_tree->drvrPt();
  const auto steiner_location = driver_steiner_tree->location(top_steiner_point);
  const double l1 = dbuToMeters(driver_steiner_tree->distance(driver_point, top_steiner_point));
  
  const float wire_cap = l1 * unit_c;
  const float wo_fst_stt_cap = total_net_cap - wire_cap;
  
  const int steiner_x = steiner_location.getX();
  const int steiner_y = steiner_location.getY();
  const float actual_dist = dbuToMeters(std::abs(flop.curr_pt_.x - steiner_x) + 
                                        std::abs(flop.curr_pt_.y - steiner_y));
  
  if (verbose) {
    std::cout << "    Distance: " << l1 << " → " << actual_dist << " m (Δ=" 
              << (actual_dist - l1) << ")" << std::endl;
  }
  
  std::vector<std::pair<int, sta::Slack>> pin_used_slacks;
  pin_used_slacks.reserve(budget_it->second.size());
  
  for (const auto& [path_idx, slack_budget] : budget_it->second) {
    float coeff = coeffs.empty() ? 0.0f : coeffs.back();
    const float new_cap = actual_dist * unit_c + wo_fst_stt_cap;
    
    for (size_t seg_idx = 0; seg_idx < coeffs.size(); ++seg_idx) {
      if (new_cap >= cap_delay[seg_idx].first && 
          new_cap <= cap_delay[seg_idx + 1].first) {
        coeff = coeffs[seg_idx];
        break;
      }
    }
    
    const float rc_delay = (actual_dist * actual_dist - l1 * l1) * unit_r * unit_c / 2;
    const float wire_delay = (actual_dist - l1) * wo_fst_stt_cap * unit_r;
    const float cell_delay = coeff * (actual_dist - l1) * unit_c / 2;
    const sta::Slack used_slack = rc_delay + wire_delay + cell_delay;
    
    pin_used_slacks.emplace_back(path_idx, used_slack);
    
    if (verbose) {
      std::cout << "    Path[" << path_idx << "]: budget=" << slack_budget 
                << ", used=" << used_slack 
                << ", remaining=" << (slack_budget - used_slack) << std::endl;
    }
  }
  
  used_slacks[out_pin] = std::move(pin_used_slacks);
}

//------------------------------------------------------------------------------
// [Level 2] FanIn Pin Processing (D pins)
//------------------------------------------------------------------------------

void
AggloCluster::calcUsedSlacksFanIn(
    const FlopUnit& flop,
    odb::dbITerm* d_pin,
    est::EstimateParasitics* est,
    double unit_r,
    double unit_c,
    std::unordered_map<odb::dbITerm*, std::vector<std::pair<int, sta::Slack>>>& used_slacks,
    bool verbose)
{
  if (verbose) {
    std::cout << "\n  [FanIn: " << d_pin->getMTerm()->getName() << "]" << std::endl;
  }
  
  const auto budget_it = flop.pin_budgets_.find(d_pin);
  if (budget_it == flop.pin_budgets_.end() || budget_it->second.empty()) {
    if (verbose) {
      std::cout << "    No timing constraints - skipping" << std::endl;
    }
    return;
  }
  
  odb::dbNet* fi_net = d_pin->getNet();
  if (!fi_net) {
    if (verbose) {
      std::cout << "    Not connected to net - skipping" << std::endl;
    }
    return;
  }
  
  odb::dbITerm* fi_net_drvr_pin = fi_net->get1stITerm();
  if (!fi_net_drvr_pin) {
    if (verbose) {
      std::cout << "    No driver pin - skipping" << std::endl;
    }
    return;
  }
  
  const sta::Pin* ipin_sta = network_->dbToSta(d_pin);
  const float ipin_cap = getPinCapacitance(ipin_sta);
  
  odb::dbInst* fi_inst = fi_net_drvr_pin->getInst();
  const auto slews = getInstanceInputSlews(fi_inst);
  
  std::vector<std::pair<float, float>> cap_delay;
  float worst_delay = 0.0f;
  
  for (const auto& [iterm, slew] : slews) {
    const auto pts = extractCapacitanceDelayPoints(
        fi_inst, iterm->getMTerm()->getName(), 
        fi_net_drvr_pin->getMTerm()->getName(), slew);
    
    if (!pts.empty() && pts.back().second >= worst_delay) {
      worst_delay = pts.back().second;
      cap_delay = pts;
    }
  }
  
  if (cap_delay.empty() || cap_delay.size() < 2) {
    if (verbose) {
      std::cout << "    Insufficient cap-delay data - skipping" << std::endl;
    }
    return;
  }
  
  std::vector<float> coeffs;
  coeffs.reserve(cap_delay.size() - 1);
  for (size_t i = 1; i < cap_delay.size(); ++i) {
    coeffs.push_back((cap_delay[i].second - cap_delay[i-1].second) / 
                     (cap_delay[i].first - cap_delay[i-1].first));
  }
  
  const sta::Pin* fi_net_drvr_pin_sta = network_->dbToSta(fi_net_drvr_pin);
  est::SteinerTree* fi_tree = est->makeSteinerTree(fi_net_drvr_pin_sta);
  
  float pin_cap = 0.0f, wire_cap = 0.0f;
  const sta::MinMax* mm = sta::MinMax::max();
  const sta::Net* fi_net_drvr_sta = network_->dbToSta(fi_net_drvr_pin->getNet());
  sta_->connectedCap(fi_net_drvr_sta, corner_, mm, pin_cap, wire_cap);
  const float total_cap = pin_cap + wire_cap;

  const auto top_pt = fi_tree->top();
  const auto top_loc = fi_tree->location(top_pt);
  const int top_x = top_loc.getX();
  const int top_y = top_loc.getY();
  
  const int pin_count = fi_tree->pinCount();
  int target_pt = -1;
  for (int i = 0; i < pin_count; i++) {
    if (fi_tree->pin(i) == ipin_sta) {
      target_pt = i;
      break;
    }
  }
  
  if (target_pt == -1) {
    if (verbose) {
      std::cout << "    Target pin not in Steiner tree - skipping" << std::endl;
    }
    return;
  }
  
  std::vector<int> node_path;
  if (!findSteinerPathRecursive(fi_tree, fi_tree->drvrPt(), target_pt, node_path)) {
    if (verbose) {
      std::cout << "    Failed to find path - skipping" << std::endl;
    }
    return;
  }
  
  if (node_path.size() < 2) {
    if (verbose) {
      std::cout << "    Path too short - skipping" << std::endl;
    }
    return;
  }
  
  const size_t total_segments = node_path.size() - 1;
  float pre_leaf_length_dbu = 0.0f;
  float final_segment_length_dbu = 0.0f;
  const int branch_count = fi_tree->branchCount();
  
  for (size_t k = 0; k < total_segments; ++k) {
    const int path_node1 = node_path[k];
    const int path_node2 = node_path[k + 1];
    int current_wire_length = 0;
    
    for (int i = 0; i < branch_count; i++) {
      odb::Point pt1, pt2;
      int steiner_pt1, steiner_pt2, wire_length;
      fi_tree->branch(i, pt1, steiner_pt1, pt2, steiner_pt2, wire_length);
      
      if ((steiner_pt1 == path_node1 && steiner_pt2 == path_node2) ||
          (steiner_pt1 == path_node2 && steiner_pt2 == path_node1)) {
        current_wire_length = wire_length;
        break;
      }
    }
    
    if (k < total_segments - 1) {
      pre_leaf_length_dbu += current_wire_length;
    } else {
      final_segment_length_dbu = current_wire_length;
    }
  }
  
  const float l1 = dbuToMeters(final_segment_length_dbu);
  const float on_path_R_wo_last = unit_r * dbuToMeters(pre_leaf_length_dbu);
  const float actual_dist = dbuToMeters(std::abs(flop.curr_pt_.x - top_x) + 
                                        std::abs(flop.curr_pt_.y - top_y));
  
  if (verbose) {
    std::cout << "    Distance: " << l1 << " → " << actual_dist << " m (Δ=" 
              << (actual_dist - l1) << ")" << std::endl;
  }
  
  std::vector<std::pair<int, sta::Slack>> pin_used_slacks;
  pin_used_slacks.reserve(budget_it->second.size());
  
  for (const auto& [path_idx, slack_budget] : budget_it->second) {
    float coeff = coeffs.empty() ? 0.0f : coeffs.back();
    const float new_cap = total_cap + (actual_dist - l1) * unit_c;
    
    for (size_t seg_idx = 0; seg_idx < coeffs.size(); ++seg_idx) {
      if (new_cap >= cap_delay[seg_idx].first && 
          new_cap <= cap_delay[seg_idx + 1].first) {
        coeff = coeffs[seg_idx];
        break;
      }
    }
    
    const float rc_delay = (actual_dist * actual_dist - l1 * l1) * unit_r * unit_c / 2;
    const float on_path_delay = (actual_dist - l1) * on_path_R_wo_last * unit_c;
    const float cell_delay = coeff * (actual_dist - l1) * unit_c / 2;
    const float ipin_delay = (actual_dist - l1) * unit_r * ipin_cap;
    const sta::Slack used_slack = rc_delay + on_path_delay + cell_delay + ipin_delay;
    
    pin_used_slacks.emplace_back(path_idx, used_slack);
    
    if (verbose) {
      std::cout << "    Path[" << path_idx << "]: budget=" << slack_budget 
                << ", used=" << used_slack 
                << ", remaining=" << (slack_budget - used_slack) << std::endl;
    }
  }
  
  used_slacks[d_pin] = std::move(pin_used_slacks);
}


//==============================================================================
// Phase 11: Implementation of Clusters Helper Functions
//==============================================================================

// -----------------------------------------------------------------------------
// [Level 1] Top-level: Convert clusters to MBFF instances
// -----------------------------------------------------------------------------

bool 
AggloCluster::implementSingleCluster(const FlopCluster& cluster, bool verbose, int debug_idx)
{
  if (verbose) {
    std::cout << "\n[Cluster #" << debug_idx << " | ID=" << cluster.id_ << "]" << std::endl;
    std::cout << "  Size: " << cluster.flops_.size() << " flop(s)" << std::endl;
    std::cout << "  Position: (" << cluster.curr_pt_.x << ", " << cluster.curr_pt_.y << ")" << std::endl;
    
    std::cout << "  Flop IDs: [";
    int shown = 0;
    for (int flop_id : cluster.flops_) {
      if (shown > 0) std::cout << ", ";
      std::cout << flop_id;
      if (++shown >= 5) {
        std::cout << ", ...";
        break;
      }
    }
    std::cout << "]" << std::endl;
  }

  // Step 1: Extract net bundles
  if (verbose) {
    std::cout << "\n  [Step 1] Extracting net bundles..." << std::endl;
  }
  
  std::vector<NetBundle> net_bundles = getNetBundles(cluster);
  
  if (net_bundles.empty()) {
    if (verbose) {
      std::cout << "    ✗ SKIP: No valid net bundles" << std::endl;
    }
    return false;
  }
  
  if (verbose) {
    std::cout << "    ✓ Extracted " << net_bundles.size() << " bundle(s)" << std::endl;
    
    for (int i = 0; i < net_bundles.size(); ++i) {
      const NetBundle& nb = net_bundles[i];
      std::cout << "      [" << i << "] Flop_" << nb.flop_id_ 
                << " (" << nb.flop_inst_->getName() << ")" << std::endl;
      
      if (nb.d_net_) {
        std::cout << "        D:  " << nb.d_net_->getName() 
                  << " (fanout=" << nb.d_net_->getITermCount() << ")" << std::endl;
      }
      if (nb.q_net_) {
        std::cout << "        Q:  " << nb.q_net_->getName() 
                  << " (fanout=" << nb.q_net_->getITermCount() << ")" << std::endl;
      }
      if (nb.qn_net_) {
        std::cout << "        QN: " << nb.qn_net_->getName() 
                  << " (fanout=" << nb.qn_net_->getITermCount() << ")" << std::endl;
      }
    }
  }
  
  // Step 2: Get port bundles for candidate masters
  if (verbose) {
    std::cout << "\n  [Step 2] Analyzing candidate MBFF masters..." << std::endl;
  }
  
  std::map<odb::dbMaster*, std::vector<PortBundle>> master_port_bundles 
      = getAllPortBundles(cluster, net_bundles);
  
  if (master_port_bundles.empty()) {
    if (verbose) {
      std::cout << "    ✗ SKIP: No compatible " << net_bundles.size() 
                << "-bit MBFF masters" << std::endl;
    }
    return false;
  }
  
  if (verbose) {
    std::cout << "    ✓ Found " << master_port_bundles.size() << " compatible master(s)" << std::endl;
    
    int master_idx = 0;
    for (const auto& [master, port_bundles] : master_port_bundles) {
      std::cout << "      [" << master_idx << "] " << master->getName() << std::endl;
      std::cout << "        Size: " << master->getWidth() << "×" << master->getHeight() 
                << " DBU (area=" << (master->getWidth() * master->getHeight()) << ")" << std::endl;
      std::cout << "        Port bundles: " << port_bundles.size() << std::endl;
      
      for (int i = 0; i < port_bundles.size(); ++i) {
        const PortBundle& pb = port_bundles[i];
        std::cout << "          [" << i << "] ";
        
        if (pb.d_mterm_) {
          std::cout << "D=" << pb.d_mterm_->getName() 
                    << "@(" << pb.d_local_pos_.get<0>() << "," << pb.d_local_pos_.get<1>() << ")";
        }
        if (pb.q_mterm_) {
          std::cout << ", Q=" << pb.q_mterm_->getName()
                    << "@(" << pb.q_local_pos_.get<0>() << "," << pb.q_local_pos_.get<1>() << ")";
        }
        if (pb.qn_mterm_) {
          std::cout << ", QN=" << pb.qn_mterm_->getName()
                    << "@(" << pb.qn_local_pos_.get<0>() << "," << pb.qn_local_pos_.get<1>() << ")";
        }
        std::cout << std::endl;
      }
      
      master_idx++;
    }
  }
  
  // Step 3: Find optimal assignment using Cost-Scaling Push-Relabel algorithm
  if (verbose) {
    std::cout << "\n  [Step 3] Computing optimal port assignment..." << std::endl;
  }
  
  MasterPortAssignment best_assignment = assignPorts(cluster, net_bundles, master_port_bundles);
  
  if (best_assignment.best_master_ == nullptr) {
    if (verbose) {
      std::cout << "    ✗ SKIP: Assignment algorithm failed" << std::endl;
    }
    return false;
  }
  
  if (verbose) {
    std::cout << "    ✓ Optimal assignment found" << std::endl;
    std::cout << "      Master: " << best_assignment.best_master_->getName() << std::endl;
    std::cout << "      Total HPWL: " << best_assignment.min_cost_ << " DBU" << std::endl;
    
    const auto& port_bundles = master_port_bundles.at(best_assignment.best_master_);
    const Point inst_center(std::lround(cluster.curr_pt_.x), std::lround(cluster.curr_pt_.y));
    
    double total_cost_check = 0.0;
    
    for (int i = 0; i < best_assignment.net_to_port_assignment_.size(); ++i) {
      const int port_idx = best_assignment.net_to_port_assignment_[i];
      const NetBundle& nb = net_bundles[i];
      const PortBundle& pb = port_bundles[port_idx];
      
      std::cout << "      [" << i << "] NetBundle(Flop_" << nb.flop_id_ 
                << ") → PortBundle[" << port_idx << "]" << std::endl;
      
      double d_cost = 0.0, q_cost = 0.0, qn_cost = 0.0;
      
      if (nb.d_net_ && pb.d_mterm_) {
        odb::Rect bbox = getNetBBoxWithoutPin(nb.d_net_, nb.d_iterm_);
        const Point global_pos = getGlobalMTermPos(pb.d_local_pos_, 
                                                   best_assignment.best_master_, 
                                                   inst_center);
        bbox.merge(odb::Point(global_pos.get<0>(), global_pos.get<1>()));
        d_cost = bbox.dx() + bbox.dy();
        std::cout << "        D:  " << nb.d_net_->getName() << " → " 
                  << pb.d_mterm_->getName() << " (HPWL=" << d_cost << ")" << std::endl;
      }
      
      if (nb.q_net_ && pb.q_mterm_) {
        odb::Rect bbox = getNetBBoxWithoutPin(nb.q_net_, nb.q_iterm_);
        const Point global_pos = getGlobalMTermPos(pb.q_local_pos_, 
                                                   best_assignment.best_master_, 
                                                   inst_center);
        bbox.merge(odb::Point(global_pos.get<0>(), global_pos.get<1>()));
        q_cost = bbox.dx() + bbox.dy();
        std::cout << "        Q:  " << nb.q_net_->getName() << " → " 
                  << pb.q_mterm_->getName() << " (HPWL=" << q_cost << ")" << std::endl;
      }
      
      if (nb.qn_net_ && pb.qn_mterm_) {
        odb::Rect bbox = getNetBBoxWithoutPin(nb.qn_net_, nb.qn_iterm_);
        const Point global_pos = getGlobalMTermPos(pb.qn_local_pos_, 
                                                   best_assignment.best_master_, 
                                                   inst_center);
        bbox.merge(odb::Point(global_pos.get<0>(), global_pos.get<1>()));
        qn_cost = bbox.dx() + bbox.dy();
        std::cout << "        QN: " << nb.qn_net_->getName() << " → " 
                  << pb.qn_mterm_->getName() << " (HPWL=" << qn_cost << ")" << std::endl;
      }
      
      total_cost_check += d_cost + q_cost + qn_cost;
    }
    
    std::cout << "      Cost verification: computed=" << total_cost_check 
              << ", assignment_solver=" << best_assignment.min_cost_;
    
    if (std::abs(total_cost_check - best_assignment.min_cost_) < 1.0) {
      std::cout << " ✓" << std::endl;
    } else {
      std::cout << " ⚠ mismatch!" << std::endl;
    }
  }
  
  // Step 4: Apply implementation
  if (verbose) {
    std::cout << "\n  [Step 4] Creating MBFF instance..." << std::endl;
  }
  
  applyImplementation(cluster, best_assignment, net_bundles);
  
  if (verbose) {
    std::cout << "    ✓ Created: mbff_cluster_" << cluster.id_ << std::endl;
    std::cout << "    ✓ Destroyed " << net_bundles.size() << " original flop(s)" << std::endl;
    std::cout << "  ✓ Cluster #" << debug_idx << " implementation COMPLETE" << std::endl;
  }
  
  return true;
}

// -----------------------------------------------------------------------------
// [Level 2] Helpers for master selection and port assignment
// -----------------------------------------------------------------------------

std::vector<NetBundle> 
AggloCluster::getNetBundles(const FlopCluster& cluster) const
{
  std::vector<NetBundle> net_bundles;
  net_bundles.reserve(cluster.flops_.size());

  for (int flop_id : cluster.flops_) {
    const FlopUnit& flop = flop_units_[flop_id];
    NetBundle bundle;
    bundle.flop_id_ = flop_id;
    bundle.flop_inst_ = flop.inst_;

    // Collect D, Q, QN nets and pins from this flop instance
    for (odb::dbITerm* iterm : flop.inst_->getITerms()) {
      if (isDPin(iterm)) {
        bundle.d_net_ = iterm->getNet();
        bundle.d_iterm_ = iterm;
      } else if (isQPin(iterm)) {
        bundle.q_net_ = iterm->getNet();
        bundle.q_iterm_ = iterm;
      } else if (isQNPin(iterm)) {
        bundle.qn_net_ = iterm->getNet();
        bundle.qn_iterm_ = iterm;
      }
    }
    
    // Valid bundle requires D net and at least one output (Q or QN)
    if (bundle.d_net_ && (bundle.q_net_ || bundle.qn_net_)) {
      net_bundles.push_back(bundle);
    }
  }
  
  return net_bundles;
}

std::map<odb::dbMaster*, std::vector<PortBundle>> 
AggloCluster::getAllPortBundles(const FlopCluster& cluster, 
                                 const std::vector<NetBundle>& net_bundles) const
{
  std::map<odb::dbMaster*, std::vector<PortBundle>> master_port_bundles;
  
  const int num_nets = net_bundles.size();
  
  // Check if compatible masters exist for this cluster's mask and bit-width
  const auto mask_it = compatible_masters_.find(cluster.master_mask_);
  if (mask_it == compatible_masters_.end()) {
    return master_port_bundles;
  }
  
  const auto bits_it = mask_it->second.find(num_nets);
  if (bits_it == mask_it->second.end()) {
    return master_port_bundles;
  }
  
  const auto& candidate_masters = bits_it->second;
  
  // Extract port bundles from each candidate master
  for (odb::dbMaster* master : candidate_masters) {
    std::vector<PortBundle> port_bundles = getPortBundles(master);
    
    // Skip masters with insufficient port bundles
    if (static_cast<int>(port_bundles.size()) < num_nets) {
      continue;
    }
    
    master_port_bundles[master] = std::move(port_bundles);
  }
  
  return master_port_bundles;
}

MasterPortAssignment 
AggloCluster::assignPorts(
    const FlopCluster& cluster,
    const std::vector<NetBundle>& net_bundles,
    const std::map<odb::dbMaster*, std::vector<PortBundle>>& master_port_bundles)
{
  typedef util::StaticGraph<> Graph;

  MasterPortAssignment best_result;
  const int num_nets = net_bundles.size();
  
  if (master_port_bundles.empty()) {
    return best_result;
  }
  
  const Point inst_center(std::lround(cluster.curr_pt_.x), 
                          std::lround(cluster.curr_pt_.y));

  // Evaluate each candidate master
  for (const auto& [master, port_bundles] : master_port_bundles) {
    const int num_ports = port_bundles.size();
    
    if (num_ports < num_nets) {
      continue;
    }
    
    // Build bipartite graph: nets (left) <-> ports (right)
    const int num_left_nodes = num_nets;
    const int num_nodes = num_nets + num_ports;
    const int num_arcs = num_nets * num_ports; 
    
    Graph graph(num_nodes, num_arcs);
    std::vector<int64_t> arc_costs;
    arc_costs.reserve(num_arcs);

    // Create edges from each net to each port with HPWL cost
    for (int i = 0; i < num_nets; ++i) {
      for (int j = 0; j < num_ports; ++j) {
        const int tail = i; 
        const int head = num_left_nodes + j;
        graph.AddArc(tail, head);
        
        const double cost_double = calcAssignmentCost(net_bundles[i], 
                                                      port_bundles[j], 
                                                      master,
                                                      inst_center);
        arc_costs.push_back(static_cast<int64_t>(std::round(cost_double)));
      }
    }

    // Solve assignment using Cost-Scaling Push-Relabel algorithm
    // (Goldberg & Kennedy, 1995) - O(n*m*log(nC))
    graph.Build();
    ::operations_research::LinearSumAssignment assignment_solver(graph, num_left_nodes);

    for (int arc = 0; arc < num_arcs; ++arc) {
      assignment_solver.SetArcCost(arc, arc_costs[arc]);
    }
    
    if (!assignment_solver.ComputeAssignment()) {
      continue;
    }
    
    const double total_cost = assignment_solver.GetCost();
    
    // Update best result if this master yields lower cost
    if (total_cost < best_result.min_cost_) {
      best_result.min_cost_ = total_cost;
      best_result.best_master_ = master;
      best_result.net_to_port_assignment_.resize(num_nets);
      
      for (int i = 0; i < num_left_nodes; ++i) {
        const int assigned_global_node_idx = assignment_solver.GetMate(i);
        const int j = assigned_global_node_idx - num_left_nodes;
        best_result.net_to_port_assignment_[i] = j;
      }
    }
  }
  
  return best_result;
}


void 
AggloCluster::applyImplementation(const FlopCluster& cluster,
                                  const MasterPortAssignment& result,
                                  const std::vector<NetBundle>& net_bundles)
{ 
  odb::dbMaster* best_master = result.best_master_;
  if (!best_master) {
    return;
  }
  
  const int num_bits = net_bundles.size();
  
  // Step 1: Create new MBFF instance
  const std::string new_inst_name = "mbff_cluster_" + std::to_string(cluster.id_);
  odb::dbInst* new_inst = odb::dbInst::create(block_, best_master, new_inst_name.c_str());
  
  if (!new_inst) {
    log_->error(utl::GPL, 9974, "Failed to create new instance: {}", new_inst_name);
    return;
  }

  // Step 2: Set instance placement (centered at cluster position)
  const int center_x = std::lround(cluster.curr_pt_.x);
  const int center_y = std::lround(cluster.curr_pt_.y);
  const int origin_x = center_x - best_master->getWidth() / 2;
  const int origin_y = center_y - best_master->getHeight() / 2;
  
  new_inst->setLocation(origin_x, origin_y);
  new_inst->setPlacementStatus(odb::dbPlacementStatus::PLACED);

  // Step 3: Connect common nets (CLK, CLR, PRE, SE, SI)
  // All flops in cluster have identical InstMask
  const InstMask& inst_mask = cluster.inst_mask_;
  
  for (odb::dbITerm* iterm : new_inst->getITerms()) {
    if (isClockPin(iterm) && inst_mask.clock_net_) {
      iterm->connect(inst_mask.clock_net_);
    } else if (isClearPin(iterm) && inst_mask.clear_net_) {
      iterm->connect(inst_mask.clear_net_);
    } else if (isPresetPin(iterm) && inst_mask.preset_net_) {
      iterm->connect(inst_mask.preset_net_);
    } else if (isScanEnablePin(iterm) && inst_mask.scan_enable_net_) {
      iterm->connect(inst_mask.scan_enable_net_);
    } else if (isScanInPin(iterm) && inst_mask.scan_in_net_) {
      iterm->connect(inst_mask.scan_in_net_);
    }
  }
  
  // Step 4: Connect data/scan nets based on Hungarian assignment
  const std::vector<PortBundle> port_bundles = getPortBundles(best_master);

  for (int i = 0; i < num_bits; ++i) {
    const int j = result.net_to_port_assignment_[i];
    
    const NetBundle& net_b = net_bundles[i];
    const PortBundle& port_b = port_bundles[j];

    // Connect D pin
    if (net_b.d_net_ && port_b.d_mterm_) {
      net_b.d_iterm_->disconnect();
      new_inst->getITerm(port_b.d_mterm_)->connect(net_b.d_net_);
    }
    
    // Connect Q pin
    if (net_b.q_net_ && port_b.q_mterm_) {
      net_b.q_iterm_->disconnect(); 
      new_inst->getITerm(port_b.q_mterm_)->connect(net_b.q_net_);
    }
    
    // Connect QN pin
    if (net_b.qn_net_ && port_b.qn_mterm_) {
      net_b.qn_iterm_->disconnect();
      new_inst->getITerm(port_b.qn_mterm_)->connect(net_b.qn_net_);
    }
  }

  // Step 5: Destroy original 1-bit instances
  for (const auto& net_b : net_bundles) {
    odb::dbInst::destroy(net_b.flop_inst_);
  }
}

// -----------------------------------------------------------------------------
// [Level 3] Lower-level helpers for implementation
// -----------------------------------------------------------------------------

double 
AggloCluster::calcAssignmentCost(const NetBundle& net_bundle,
                                 const PortBundle& port_bundle,
                                 odb::dbMaster* master,
                                 const Point& new_inst_center) const
{
  double total_hpwl_cost = 0.0;

  // Calculate HPWL cost for D pin
  if (net_bundle.d_net_ && port_bundle.d_mterm_) {
    odb::Rect bbox = getNetBBoxWithoutPin(net_bundle.d_net_, net_bundle.d_iterm_);
    const Point global_pos = getGlobalMTermPos(port_bundle.d_local_pos_, 
                                               master, 
                                               new_inst_center);
    bbox.merge(odb::Point(global_pos.get<0>(), global_pos.get<1>()));
    total_hpwl_cost += bbox.dx() + bbox.dy();
  }

  // Calculate HPWL cost for Q pin
  if (net_bundle.q_net_ && port_bundle.q_mterm_) {
    odb::Rect bbox = getNetBBoxWithoutPin(net_bundle.q_net_, net_bundle.q_iterm_);
    const Point global_pos = getGlobalMTermPos(port_bundle.q_local_pos_, 
                                               master, 
                                               new_inst_center);
    bbox.merge(odb::Point(global_pos.get<0>(), global_pos.get<1>()));
    total_hpwl_cost += bbox.dx() + bbox.dy();
  }

  // Calculate HPWL cost for QN pin
  if (net_bundle.qn_net_ && port_bundle.qn_mterm_) {
    odb::Rect bbox = getNetBBoxWithoutPin(net_bundle.qn_net_, net_bundle.qn_iterm_);
    const Point global_pos = getGlobalMTermPos(port_bundle.qn_local_pos_, 
                                               master, 
                                               new_inst_center);
    bbox.merge(odb::Point(global_pos.get<0>(), global_pos.get<1>()));
    total_hpwl_cost += bbox.dx() + bbox.dy();
  }
  
  return total_hpwl_cost;
}

odb::Rect 
AggloCluster::getNetBBoxWithoutPin(odb::dbNet* net, odb::dbITerm* pin_to_ignore) const
{
  odb::Rect bbox;
  bbox.mergeInit();
  
  if (!net) {
    return bbox;
  }

  for (odb::dbITerm* iterm : net->getITerms()) {
    if (iterm == pin_to_ignore) {
      continue;
    }
    FloatPoint pin_coord = getPinCoordinate(iterm);
    bbox.merge(odb::Point(static_cast<int>(pin_coord.x), static_cast<int>(pin_coord.y)));
  }

  for (odb::dbBTerm* bterm : net->getBTerms()) {
    FloatPoint pin_coord = getPinCoordinate(bterm);
    bbox.merge(odb::Point(static_cast<int>(pin_coord.x), static_cast<int>(pin_coord.y)));
  }
  
  if (bbox.isInverted()) {
      bbox.set_xlo(0); bbox.set_ylo(0); bbox.set_xhi(0); bbox.set_yhi(0);
  }

  return bbox;
}

Point 
AggloCluster::getGlobalMTermPos(const Point& local_port_pos, 
                                odb::dbMaster* master,
                                const Point& inst_center) const
{
  // Calculate instance origin from center position
  const int center_x = inst_center.get<0>();
  const int center_y = inst_center.get<1>();
  const int origin_x = center_x - master->getWidth() / 2;
  const int origin_y = center_y - master->getHeight() / 2;

  // Transform local port position to global coordinates
  return Point(origin_x + local_port_pos.get<0>(), 
               origin_y + local_port_pos.get<1>());
}

std::vector<PortBundle> 
AggloCluster::getPortBundles(odb::dbMaster* master) const
{
  std::vector<PortBundle> port_bundles;

  // Create temporary instance to analyze pin properties
  const char* temp_inst_name = "_temp_port_bundle_check";
  odb::dbInst* temp_inst = odb::dbInst::create(block_, master, temp_inst_name);
  
  if (!temp_inst) {
    return port_bundles;
  }

  // Collect D, Q, QN pins separately using isDPin/isQPin/isQNPin functions
  std::vector<odb::dbITerm*> d_pins;
  std::vector<odb::dbITerm*> q_pins;
  std::vector<odb::dbITerm*> qn_pins;

  for (odb::dbITerm* iterm : temp_inst->getITerms()) {
    if (isDPin(iterm)) {
      d_pins.push_back(iterm);
    } else if (isQPin(iterm)) {
      q_pins.push_back(iterm);
    } else if (isQNPin(iterm)) {
      qn_pins.push_back(iterm);
    }
  }

  // Assume the order is preserved: (D0, Q0, QN0), (D1, Q1, QN1), ...
  // The number of bundles is determined by the number of D pins
  const int num_bundles = d_pins.size();
  port_bundles.reserve(num_bundles);

  for (int i = 0; i < num_bundles; ++i) {
    PortBundle bundle;
    bundle.bundle_idx_ = i;
    bundle.master_ = master;

    // D pin (must exist)
    if (i < d_pins.size()) {
      odb::dbITerm* d_iterm = d_pins[i];
      bundle.d_mterm_ = d_iterm->getMTerm();
      
      // Get D pin local position
      odb::Rect bbox = bundle.d_mterm_->getBBox();
      bundle.d_local_pos_.set<0>(bbox.xCenter());
      bundle.d_local_pos_.set<1>(bbox.yCenter());
    }

    // Q pin (may not exist for all bundles)
    if (i < q_pins.size()) {
      odb::dbITerm* q_iterm = q_pins[i];
      bundle.q_mterm_ = q_iterm->getMTerm();
      
      // Get Q pin local position
      odb::Rect bbox = bundle.q_mterm_->getBBox();
      bundle.q_local_pos_.set<0>(bbox.xCenter());
      bundle.q_local_pos_.set<1>(bbox.yCenter());
    }

    // QN pin (may not exist for all bundles)
    if (i < qn_pins.size()) {
      odb::dbITerm* qn_iterm = qn_pins[i];
      bundle.qn_mterm_ = qn_iterm->getMTerm();
      
      // Get QN pin local position
      odb::Rect bbox = bundle.qn_mterm_->getBBox();
      bundle.qn_local_pos_.set<0>(bbox.xCenter());
      bundle.qn_local_pos_.set<1>(bbox.yCenter());
    }

    port_bundles.push_back(bundle);
  }

  // Destroy temporary instance
  odb::dbInst::destroy(temp_inst);

  return port_bundles;
}

//==============================================================================
// Unit Conversion Helper Functions
//==============================================================================

double 
AggloCluster::dbuToMeters(int dist) const
{
  int dbu = db_->getTech()->getDbUnitsPerMicron();
  return dist / (dbu * 1e+6);
}

int 
AggloCluster::metersToDbu(double dist) const
{
  int dbu = db_->getTech()->getDbUnitsPerMicron();
  return dist * dbu * 1e+6;
}

int64_t
AggloCluster::computeTotalHpwl() const
{
  int64_t total_hpwl = 0;

  for (odb::dbNet* net : block_->getNets()) {
    if (net == nullptr) {
      continue;
    }

    const odb::dbSigType sig_type = net->getSigType();
    if (sig_type != odb::dbSigType::SIGNAL && sig_type != odb::dbSigType::CLOCK) {
      continue;
    }

    int min_x = std::numeric_limits<int>::max();
    int max_x = std::numeric_limits<int>::min();
    int min_y = std::numeric_limits<int>::max();
    int max_y = std::numeric_limits<int>::min();
    bool has_pin = false;

    for (odb::dbITerm* iterm : net->getITerms()) {
      int x = 0;
      int y = 0;
      if (iterm->getAvgXY(&x, &y)) {
        min_x = std::min(min_x, x);
        max_x = std::max(max_x, x);
        min_y = std::min(min_y, y);
        max_y = std::max(max_y, y);
        has_pin = true;
      }
    }

    for (odb::dbBTerm* bterm : net->getBTerms()) {
      odb::Rect bbox = bterm->getBBox();
      if (!bbox.isInverted()) {
        const int x = bbox.xCenter();
        const int y = bbox.yCenter();
        min_x = std::min(min_x, x);
        max_x = std::max(max_x, x);
        min_y = std::min(min_y, y);
        max_y = std::max(max_y, y);
        has_pin = true;
      }
    }

    if (has_pin) {
      total_hpwl += static_cast<int64_t>(max_x - min_x)
                    + static_cast<int64_t>(max_y - min_y);
    }
  }

  return total_hpwl;
}

//==============================================================================
// Coordinate Transformation Helper Functions
//==============================================================================

Point 
AggloCluster::transformCoords(const Point& p) const
{
  const int x = p.get<0>();
  const int y = p.get<1>();
  return Point(x + y, y - x);
}

Point 
AggloCluster::inverseTransformCoords(const Point& p) const
{
  const int x_prime = p.get<0>();
  const int y_prime = p.get<1>();
  return Point((x_prime - y_prime) / 2, (x_prime + y_prime) / 2);
}

//==============================================================================
// Other Utility Helper Functions
//==============================================================================

unsigned int 
AggloCluster::roundDownToPowerOfTwo(unsigned int x)

{
  x |= (x >> 1);
  x |= (x >> 2);
  x |= (x >> 4);
  x |= (x >> 8);
  x |= (x >> 16);
  return x ^ (x >> 1);
}





}  // end namespace gpl
