#include "aggloCluster.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>

#include "ortools/graph/linear_assignment.h"
#include "ortools/graph/graph.h"
#include "db_sta/dbNetwork.hh"
#include "db_sta/dbSta.hh"
#include "odb/db.h"
#include "rsz/Resizer.hh"
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
#include "utl/Logger.h"
#include "placerBase.h" 

namespace gpl {

//==============================================================================
// MasterMask Implementation
//==============================================================================

std::string MasterMask::to_string() const
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

std::string InstMask::to_string() const
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

std::string FlopUnit::to_string() const
{
  std::ostringstream oss;
  oss << "Inst: " << inst_->getName() << "\n";
  oss << "  - Master Mask: " << master_mask_.to_string() << "\n";
  oss << "  - Inst Mask:   " << inst_mask_.to_string() << "\n";
  oss << "  - Orig Pt:     " << orig_pt_.to_string();
  return oss.str();
}

//==============================================================================
// VirtualBinGrid Implementation
//==============================================================================

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

  // Accumulate area statistics from all bins
  for (const auto& bin : bins_) {
    total_overflow_area += bin.getOverflowArea();
    total_inst_placed_area += bin.getInstPlacedArea();
    total_macro_placed_area += bin.getMacroPlacedArea();
    total_non_place_area += bin.getNonPlaceArea();
  }

  // Calculate total occupied area considering target density
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
  // Calculate bin indices using floating-point division
  const int lower_idx = static_cast<int>(std::floor((box.xMin() - lx_) / bin_size_x_));
  const int upper_idx = static_cast<int>(std::ceil((box.xMax() - lx_) / bin_size_x_));
  
  // Clamp to valid range [0, bin_cnt_x_) - upper_idx is exclusive
  return {std::max(0, lower_idx), std::min(bin_cnt_x_, upper_idx)};
}

std::pair<int, int> 
VirtualBinGrid::getMinMaxIdxY(const odb::Rect& box) const
{
  // Calculate bin indices using floating-point division
  const int lower_idx = static_cast<int>(std::floor((box.yMin() - ly_) / bin_size_y_));
  const int upper_idx = static_cast<int>(std::ceil((box.yMax() - ly_) / bin_size_y_));
  
  // Clamp to valid range [0, bin_cnt_y_) - upper_idx is exclusive
  return {std::max(0, lower_idx), std::min(bin_cnt_y_, upper_idx)};
}

//------------------------------------------------------------------------------
// Instance Placement Update (Helper)
//------------------------------------------------------------------------------

void 
VirtualBinGrid::updateInstPlacement(odb::dbMaster* master, 
                                    const Point& center_pos, 
                                    bool is_add)
{
  // Calculate instance bounding box
  const int w = master->getWidth();
  const int h = master->getHeight();
  const int lx = center_pos.get<0>() - w / 2;
  const int ly = center_pos.get<1>() - h / 2;
  const odb::Rect bbox(lx, ly, lx + w, ly + h);

  // Find bins that overlap with instance (upper_idx is exclusive)
  const auto [min_x_idx, max_x_idx] = getMinMaxIdxX(bbox);
  const auto [min_y_idx, max_y_idx] = getMinMaxIdxY(bbox);

  // Update each overlapping bin
  for (int y = min_y_idx; y < max_y_idx; ++y) {
    for (int x = min_x_idx; x < max_x_idx; ++x) {
      const int bin_idx = y * bin_cnt_x_ + x;
      VirtualBin& bin = bins_[bin_idx];

      // Calculate overlap area between instance and bin
      const int overlap_lx = std::max(bbox.xMin(), bin.lx());
      const int overlap_ly = std::max(bbox.yMin(), bin.ly());
      const int overlap_ux = std::min(bbox.xMax(), bin.ux());
      const int overlap_uy = std::min(bbox.yMax(), bin.uy());

      // Only update if there's actual overlap
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

//------------------------------------------------------------------------------
// Public API
//------------------------------------------------------------------------------

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
                           int num_paths_per_endpoint,
                           int threads,
                           float feasible_region_bin_multiplier,
                           bool verbose)
    : db_(db),
      block_(db->getChip()->getBlock()),
      sta_(sta),
      network_(sta_->getDbNetwork()),
      corner_(sta_->cmdCorner()),
      log_(log),
      resizer_(resizer),
      target_density_(target_density),
      target_overflow_(target_overflow),
      num_paths_per_endpoint_(num_paths_per_endpoint),
      threads_(threads),
      feasible_region_bin_multiplier_(feasible_region_bin_multiplier),
      verbose_(verbose)
{
}

// Destructor
AggloCluster::~AggloCluster() = default;

//------------------------------------------------------------------------------
// Main Entry Point
//------------------------------------------------------------------------------

void AggloCluster::doAggloCluster()
{
  initVirtualBinGrid();

  readCompatibleMasters();

  readFlopUnits();

  createCompatibleGroups();

  readTimingPaths();

  analyzeTimingPaths();

  calcFeasibleRegions();

  createFlopClusters();

  createCompatibilityGraph();

  runAgglomerativeClustering();

  exit(0); // Temporary exit for debugging

  implementClusters();
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

  //----------------------------------------------------------------------------
  // Step 1: Calculate total core area
  //----------------------------------------------------------------------------
  const odb::Rect& core_area = block_->getCoreArea();
  const int64_t total_core_area = static_cast<int64_t>(core_area.dx()) 
                                  * static_cast<int64_t>(core_area.dy());

  if (verbose_) {
    std::cout << "[Step 1] Core Area Calculation" << std::endl;
    std::cout << "  - Core bounds: (" << core_area.xMin() << ", " << core_area.yMin() 
              << ") -> (" << core_area.xMax() << ", " << core_area.yMax() << ")" << std::endl;
    std::cout << "  - Core dimensions: " << core_area.dx() << " x " << core_area.dy() << std::endl;
    std::cout << "  - Total core area: " << total_core_area << " DBU^2" << std::endl;
  }

  //----------------------------------------------------------------------------
  // Step 2: Calculate average placeable instance area
  //----------------------------------------------------------------------------
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

  //----------------------------------------------------------------------------
  // Step 3: Determine ideal bin count
  //----------------------------------------------------------------------------
  const int64_t ideal_bin_area = (target_density_ > 0)
      ? std::round(static_cast<float>(avg_inst_area) / target_density_)
      : 0;

  int ideal_bin_cnt = (ideal_bin_area > 0) 
      ? total_core_area / ideal_bin_area 
      : 4;
  ideal_bin_cnt = std::max(ideal_bin_cnt, 4);  // Minimum 2x2 grid

  if (verbose_) {
    std::cout << "[Step 3] Ideal Bin Count Determination" << std::endl;
    std::cout << "  - Target density: " << target_density_ << std::endl;
    std::cout << "  - Ideal bin area: " << ideal_bin_area << " DBU^2" << std::endl;
    std::cout << "  - Ideal bin count: " << ideal_bin_cnt << std::endl;
  }

  //----------------------------------------------------------------------------
  // Step 4: Calculate bin grid dimensions (power-of-2 based)
  //----------------------------------------------------------------------------
  const int width = core_area.dx();
  const int height = core_area.dy();
  const int aspect_ratio = roundDownToPowerOfTwo(std::max(width, height) 
                                                  / std::min(width, height));

  // Find optimal bin count (power of 2) that approximates ideal_bin_cnt
  int base_bin_cnt = 2;
  for (; base_bin_cnt <= 1024; base_bin_cnt *= 2) {
    const int total_bins = base_bin_cnt * (base_bin_cnt * aspect_ratio);
    if ((base_bin_cnt == 2 || total_bins <= ideal_bin_cnt) 
        && 4 * total_bins > ideal_bin_cnt) {
      break;
    }
  }

  // Apply aspect ratio to bin counts
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

  //----------------------------------------------------------------------------
  // Step 5: Create virtual bin grid
  //----------------------------------------------------------------------------
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

  //----------------------------------------------------------------------------
  // Step 6: Populate bins with existing instance areas
  //----------------------------------------------------------------------------
  auto& bins = virtual_bin_grid_.getBins();
  
  int64_t total_inst_area = 0;
  int64_t total_macro_area = 0;
  int64_t total_fixed_area = 0;
  int processed_inst_count = 0;

  // Sample instances for detailed debug output (first 5)
  constexpr int DEBUG_SAMPLE_SIZE = 5;

  for (odb::dbInst* inst : block_->getInsts()) {
    const odb::Rect inst_bbox = inst->getBBox()->getBox();
    const auto [min_x, max_x] = virtual_bin_grid_.getMinMaxIdxX(inst_bbox);
    const auto [min_y, max_y] = virtual_bin_grid_.getMinMaxIdxY(inst_bbox);

    const bool show_debug = verbose_ && (processed_inst_count < DEBUG_SAMPLE_SIZE);

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

    // Update all bins that overlap with this instance
    for (int y = min_y; y < max_y; ++y) {
      for (int x = min_x; x < max_x; ++x) {
        VirtualBin& bin = bins[y * bin_cnt_x + x];
        
        // Calculate overlap between instance and bin
        const odb::Rect bin_rect(bin.lx(), bin.ly(), bin.ux(), bin.uy());
        const odb::Rect overlap = inst_bbox.intersect(bin_rect);
        
        if (overlap.area() == 0) {
          continue;
        }

        // Classify and add area based on instance type
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
    std::cout << "\n[Step 6] Bin Population Completed" << std::endl;
    std::cout << "  - Total instances processed: " << processed_inst_count << std::endl;
    std::cout << "  - Standard cell area: " << total_inst_area << " DBU^2" << std::endl;
    std::cout << "  - Macro area: " << total_macro_area << " DBU^2" << std::endl;
    std::cout << "  - Fixed/blockage area: " << total_fixed_area << " DBU^2" << std::endl;
    std::cout << "  - Total occupied area: " << (total_inst_area + total_macro_area + total_fixed_area) << " DBU^2" << std::endl;
    std::cout << "  - Core utilization: " << std::fixed << std::setprecision(2) 
              << (100.0 * (total_inst_area + total_macro_area + total_fixed_area) / total_core_area) << "%" << std::endl;
    std::cout << "[initVirtualBinGrid] Completed successfully.\n" << std::endl;
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

  //----------------------------------------------------------------------------
  // Step 1: Scan all library masters and group by (MasterMask, bit_width)
  //----------------------------------------------------------------------------
  std::map<MasterMask, std::map<int, std::vector<odb::dbMaster*>>> master_groups;
  std::map<MasterMask, std::map<int, std::set<std::string>>> master_name_tracker;  // Track unique names
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
          if (verbose_) {
            std::cout << "[DEBUG] Duplicate master detected: " << master_name 
                      << " (bit_width: " << bit_width << ")" << std::endl;
          }
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
    std::cout << "[Step 1] Library Scan Completed" << std::endl;
    std::cout << "  - Total masters scanned: " << total_masters << std::endl;
    std::cout << "  - Valid flop masters found: " << valid_flop_masters << std::endl;
    std::cout << "  - Duplicate masters skipped: " << duplicate_masters << std::endl;
    std::cout << "  - Unique master masks: " << master_groups.size() << std::endl;
  }

  //----------------------------------------------------------------------------
  // Step 2: Filter groups that support multi-bit clustering (>= 2 bit widths)
  //----------------------------------------------------------------------------
  for (const auto& [mask, bit_map] : master_groups) {
    if (bit_map.size() >= 2) {
      compatible_masters_[mask] = bit_map;
    }
  }

  if (verbose_) {
    std::cout << "[Step 2] Multi-bit Compatibility Filter" << std::endl;
    std::cout << "  - Compatible master groups: " << compatible_masters_.size() << std::endl;
  }

  //----------------------------------------------------------------------------
  // Step 3: Select representative master (largest area) for each bit width
  //----------------------------------------------------------------------------
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

  //----------------------------------------------------------------------------
  // Debug output: Display compatible master groups and representatives
  //----------------------------------------------------------------------------
  if (verbose_) {
    std::cout << "[Step 3] Representative Masters Selected" << std::endl;
    std::cout << "  - Total representative groups: " << representative_masters_.size() << std::endl;

    // Sample output for first few groups
    constexpr int DEBUG_SAMPLE_SIZE = 3;
    int group_count = 0;

    for (const auto& [mask, bit_map] : compatible_masters_) {
      if (group_count >= DEBUG_SAMPLE_SIZE) {
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
        int repr_area = repr_master->getWidth() * repr_master->getHeight();
        
        // Show all masters with their areas
        for (odb::dbMaster* master : masters) {
          int area = master->getWidth() * master->getHeight();
          bool is_repr = (master == repr_master);
          
          std::cout << "      " << (is_repr ? "[REPR] " : "       ") 
                    << master->getName() 
                    << " (Area: " << area << " DBU^2";
          
          if (!is_repr && area < repr_area) {
            std::cout << ", " << ((repr_area - area) * 100 / repr_area) << "% smaller";
          }
          
          std::cout << ")" << std::endl;
        }
      }

      group_count++;
    }

    if (compatible_masters_.size() > DEBUG_SAMPLE_SIZE) {
      std::cout << "\n  ... and " << (compatible_masters_.size() - DEBUG_SAMPLE_SIZE) 
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
    std::cout << "[Step 1] Flop Unit Extraction Completed" << std::endl;
    std::cout << "  - Total instances checked: " << total_instances << std::endl;
    std::cout << "  - Valid flop instances found: " << flop_units_.size() << std::endl;

    // Sample output for first few flop units
    constexpr int DEBUG_SAMPLE_SIZE = 3;
    if (flop_units_.size() > 0) {
      std::cout << "\n[Sample Flop Units]" << std::endl;
      for (size_t i = 0; i < std::min(DEBUG_SAMPLE_SIZE, static_cast<int>(flop_units_.size())); ++i) {
        const FlopUnit& flop = flop_units_[i];
        std::cout << "  Flop #" << i << ": " << flop.inst_->getName() << std::endl;
        std::cout << "    Master: " << flop.inst_->getMaster()->getName() << std::endl;
        std::cout << "    Master Mask: " << flop.master_mask_.to_string() << std::endl;
        std::cout << "    Inst Mask: " << flop.inst_mask_.to_string() << std::endl;
        std::cout << "    Location: (" << flop.orig_pt_.x << ", " 
                  << flop.orig_pt_.y << ")" << std::endl;
      }
      if (flop_units_.size() > DEBUG_SAMPLE_SIZE) {
        std::cout << "  ... and " << (flop_units_.size() - DEBUG_SAMPLE_SIZE) 
                  << " more flop units" << std::endl;
      }
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

  compatible_groups_.clear();

  for (size_t i = 0; i < flop_units_.size(); ++i) {
    const FlopUnit& unit = flop_units_[i];
    compatible_groups_[{unit.master_mask_, unit.inst_mask_}].push_back(i);
  }

  if (verbose_) {
    std::cout << "[Step 1] Group Creation Completed" << std::endl;
    std::cout << "  - Total compatible groups created: " << compatible_groups_.size() << std::endl;

    // Sort groups by size (largest first) for better visibility
    std::vector<std::pair<std::pair<MasterMask, InstMask>, std::vector<int>>> sorted_groups;
    for (const auto& group : compatible_groups_) {
      sorted_groups.push_back(group);
    }
    std::sort(sorted_groups.begin(), sorted_groups.end(),
              [](const auto& a, const auto& b) {
                return a.second.size() > b.second.size();
              });

    // Show statistics
    std::cout << "\n[Step 2] Group Size Distribution" << std::endl;
    std::map<size_t, int> size_dist;
    for (const auto& group : sorted_groups) {
      size_dist[group.second.size()]++;
    }
    for (const auto& [size, count] : size_dist) {
      std::cout << "  - Groups with " << size << " flop(s): " << count << std::endl;
    }

    // Show details of top groups
    constexpr int DEBUG_SAMPLE_SIZE = 5;
    if (!sorted_groups.empty()) {
      std::cout << "\n[Sample Compatible Groups]" << std::endl;
      int group_idx = 0;
      for (const auto& [masks, flop_ids] : sorted_groups) {
        if (group_idx >= DEBUG_SAMPLE_SIZE) break;
        
        const auto& [master_mask, inst_mask] = masks;
        std::cout << "\n  Group #" << group_idx << ": " << flop_ids.size() << " flops" << std::endl;
        std::cout << "    Master Mask: " << master_mask.to_string() << std::endl;
        std::cout << "    Inst Mask:   " << inst_mask.to_string() << std::endl;
        
        // Show first 3 flop instances in this group
        std::cout << "    Sample flops (first 3): ";
        for (size_t i = 0; i < std::min(size_t(3), flop_ids.size()); ++i) {
          if (i > 0) std::cout << ", ";
          std::cout << flop_units_[flop_ids[i]].inst_->getName();
        }
        if (flop_ids.size() > 3) {
          std::cout << " ... (+" << (flop_ids.size() - 3) << " more)";
        }
        std::cout << std::endl;
        
        group_idx++;
      }

      if (sorted_groups.size() > DEBUG_SAMPLE_SIZE) {
        std::cout << "\n  ... and " << (sorted_groups.size() - DEBUG_SAMPLE_SIZE) 
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

  //----------------------------------------------------------------------------
  // Step 1: Initialize STA and find timing paths
  //----------------------------------------------------------------------------
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

  //----------------------------------------------------------------------------
  // Step 2: Filter and extract flop-to-flop timing paths
  //----------------------------------------------------------------------------
  int total_paths_analyzed = 0;
  int valid_flop_paths = 0;
  int skipped_short_paths = 0;
  int skipped_non_flop_paths = 0;

  for (sta::PathEnd* path_end : path_ends) {
    total_paths_analyzed++;
    
    sta::Slack slack = path_end->slack(sta_);
    sta::Path* path = path_end->path();
    sta::PathExpanded expanded(path, sta_);

    if (expanded.size() < 2) {
      skipped_short_paths++;
      continue;
    }

    // Start FF
    const sta::Path* start_path = expanded.path(expanded.startIndex());
    sta::Vertex* sta_start_vertex
        = start_path ? start_path->vertex(sta_) : nullptr;
    sta::Pin* sta_start_pin
        = sta_start_vertex ? sta_start_vertex->pin() : nullptr;
    sta::Instance* sta_start_inst
        = sta_start_pin ? network_->instance(sta_start_pin) : nullptr;
    odb::dbInst* db_start_inst
        = sta_start_inst ? network_->staToDb(sta_start_inst) : nullptr;

    // End FF
    sta::Vertex* sta_end_vertex = path_end->vertex(sta_);
    sta::Pin* sta_end_pin = sta_end_vertex ? sta_end_vertex->pin() : nullptr;
    sta::Instance* sta_end_inst
        = sta_end_pin ? network_->instance(sta_end_pin) : nullptr;
    odb::dbInst* db_end_inst
        = sta_end_inst ? network_->staToDb(sta_end_inst) : nullptr;

    if (!db_start_inst || !db_end_inst
        || inst_to_flop_id_.find(db_start_inst) == inst_to_flop_id_.end()
        || inst_to_flop_id_.find(db_end_inst) == inst_to_flop_id_.end()) {
      skipped_non_flop_paths++;
      continue;
    }

    int start_flop_id = inst_to_flop_id_[db_start_inst];
    int end_flop_id = inst_to_flop_id_[db_end_inst];

    // Extract the actual pins involved in this timing path
    // Start pin: output pin of the start flop (Q or QN)
    odb::dbITerm* start_pin = sta_start_pin ? network_->flatPin(sta_start_pin) : nullptr;
    
    // End pin: input pin of the end flop (D)
    odb::dbITerm* end_pin = sta_end_pin ? network_->flatPin(sta_end_pin) : nullptr;

    timing_paths_.emplace_back(path, slack, start_flop_id, end_flop_id, start_pin, end_pin);
    int path_idx = timing_paths_.size() - 1;

    valid_flop_paths++;
  }

  if (verbose_) {
    std::cout << "[Step 2] Path Filtering Completed" << std::endl;
    std::cout << "  - Total paths analyzed: " << total_paths_analyzed << std::endl;
    std::cout << "  - Valid flop-to-flop paths: " << valid_flop_paths << std::endl;
    std::cout << "  - Skipped (short paths): " << skipped_short_paths << std::endl;
    std::cout << "  - Skipped (non-flop endpoints): " << skipped_non_flop_paths << std::endl;

    //--------------------------------------------------------------------------
    // Step 3: Analyze slack distribution
    //--------------------------------------------------------------------------
    if (valid_flop_paths > 0) {
      std::vector<sta::Slack> slacks;
      slacks.reserve(timing_paths_.size());
      for (const auto& tp : timing_paths_) {
        slacks.push_back(tp.slack_);
      }
      std::sort(slacks.begin(), slacks.end());

      sta::Slack min_slack = slacks.front();
      sta::Slack max_slack = slacks.back();
      sta::Slack median_slack = slacks[slacks.size() / 2];
      sta::Slack avg_slack = std::accumulate(slacks.begin(), slacks.end(), 0.0) / slacks.size();

      int negative_slack_count = std::count_if(slacks.begin(), slacks.end(), 
                                                [](sta::Slack s) { return s < 0; });

      // Convert to user units for display
      sta::Unit* time_unit = sta_->units()->timeUnit();
      const char* time_suffix = time_unit->scaledSuffix();

      std::cout << "\n[Step 3] Slack Distribution Analysis" << std::endl;
      std::cout << "  - Min slack: " << time_unit->asString(min_slack, 3) << time_suffix << std::endl;
      std::cout << "  - Max slack: " << time_unit->asString(max_slack, 3) << time_suffix << std::endl;
      std::cout << "  - Median slack: " << time_unit->asString(median_slack, 3) << time_suffix << std::endl;
      std::cout << "  - Average slack: " << time_unit->asString(avg_slack, 3) << time_suffix << std::endl;
      std::cout << "  - Negative slack paths: " << negative_slack_count 
                << " (" << std::fixed << std::setprecision(1) << (100.0 * negative_slack_count / valid_flop_paths) << "%)" << std::endl;

      // Sample critical paths
      constexpr int DEBUG_SAMPLE_SIZE = 3;
      std::cout << "\n[Sample Critical Paths (worst slack)]" << std::endl;
      for (int i = 0; i < std::min(DEBUG_SAMPLE_SIZE, static_cast<int>(timing_paths_.size())); ++i) {
        const auto& tp = timing_paths_[i];
        const FlopUnit& start_flop = flop_units_[tp.start_flop_idx_];
        const FlopUnit& end_flop = flop_units_[tp.end_flop_idx_];
        
        std::cout << "  Path #" << i << ": " << start_flop.inst_->getName();
        if (tp.start_pin_) {
          std::cout << "/" << tp.start_pin_->getMTerm()->getName();
        }
        std::cout << " -> " << end_flop.inst_->getName();
        if (tp.end_pin_) {
          std::cout << "/" << tp.end_pin_->getMTerm()->getName();
        }
        std::cout << std::endl;
        std::cout << "    Slack: " << time_unit->asString(slacks[i], 3) << time_suffix << std::endl;
        
        // Display pin details
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
      
      if (timing_paths_.size() > DEBUG_SAMPLE_SIZE) {
        std::cout << "  ... and " << (timing_paths_.size() - DEBUG_SAMPLE_SIZE) 
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

  // Clear existing budgets
  for (FlopUnit& flop : flop_units_) {
    flop.pin_budgets_.clear();
  }

  //----------------------------------------------------------------------------
  // Step 1: Analyze each timing path and distribute slack budgets
  //----------------------------------------------------------------------------
  for (size_t i = 0; i < timing_paths_.size(); ++i) {
    const auto& path = timing_paths_[i];

    // Calculate slack budget (split equally, including negative slack)
    const sta::Slack budget = path.slack_ / 2.0f;

    // Start Flop: Store budget for its output pin (Q or QN)
    if (path.start_flop_idx_ >= 0 && path.start_pin_) {
      FlopUnit& start_flop = flop_units_[path.start_flop_idx_];
      start_flop.pin_budgets_[path.start_pin_].emplace_back(i, budget);
    }

    // End Flop: Store budget for its input pin (D)
    if (path.end_flop_idx_ >= 0 && path.end_pin_) {
      FlopUnit& end_flop = flop_units_[path.end_flop_idx_];
      end_flop.pin_budgets_[path.end_pin_].emplace_back(i, budget);
    }
  }

  if (verbose_) {
    std::cout << "[Step 1] Slack Budget Distribution Completed" << std::endl;
    std::cout << "  - Total timing paths processed: " << timing_paths_.size() << std::endl;
  }

  //----------------------------------------------------------------------------
  // Step 2: Sort pin budgets in ascending order of slack budget
  //----------------------------------------------------------------------------
  for (FlopUnit& flop : flop_units_) {
    for (auto& [pin, budget_list] : flop.pin_budgets_) {
      std::sort(budget_list.begin(),
                budget_list.end(),
                [](const auto& a, const auto& b) { 
                  return a.second < b.second;  // Sort by slack budget (ascending)
                });
    }
  }

  if (verbose_) {
    std::cout << "[Step 2] Pin Budget Sorting Completed" << std::endl;
  }

  //----------------------------------------------------------------------------
  // Step 3: Statistics and sample flop analysis
  //----------------------------------------------------------------------------
  if (verbose_) {
    int flops_with_constraints = 0;
    int total_pin_constraints = 0;
    int total_d_pins = 0;
    int total_q_qn_pins = 0;
    
    for (const auto& flop : flop_units_) {
      if (!flop.pin_budgets_.empty()) {
        flops_with_constraints++;
        for (const auto& [pin, budget_list] : flop.pin_budgets_) {
          total_pin_constraints++;
          if (isDPin(pin)) {
            total_d_pins++;
          } else {
            total_q_qn_pins++;
          }
        }
      }
    }
    
    std::cout << "\n[Step 3] Timing Constraint Statistics" << std::endl;
    std::cout << "  - Flops with timing constraints: " << flops_with_constraints 
              << " / " << flop_units_.size() << std::endl;
    std::cout << "  - Total pin-level constraints: " << total_pin_constraints << std::endl;
    std::cout << "  - D pin constraints: " << total_d_pins << std::endl;
    std::cout << "  - Q/QN pin constraints: " << total_q_qn_pins << std::endl;

    // Sample flop analysis - show first few flops with constraints
    constexpr int DEBUG_SAMPLE_SIZE = 3;
    std::cout << "\n[Sample Flop Pin Budgets]" << std::endl;
    
    sta::Unit* time_unit = sta_->units()->timeUnit();
    const char* time_suffix = time_unit->scaledSuffix();
    
    int sample_count = 0;
    for (const auto& flop : flop_units_) {
      if (sample_count >= DEBUG_SAMPLE_SIZE) {
        break;
      }
      
      if (flop.pin_budgets_.empty()) {
        continue;
      }
      
      std::cout << "  Flop #" << sample_count << ": " << flop.inst_->getName() << std::endl;
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
        std::cout << "        Paths: " << budget_list.size() << std::endl;
        std::cout << "        Min slack budget: " << time_unit->asString(min_budget, 3) << time_suffix << std::endl;
        std::cout << "        Max slack budget: " << time_unit->asString(max_budget, 3) << time_suffix << std::endl;
      }
      
      sample_count++;
    }
    
    if (flops_with_constraints > DEBUG_SAMPLE_SIZE) {
      std::cout << "  ... and " << (flops_with_constraints - DEBUG_SAMPLE_SIZE) 
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

  int flop_count = 0;
  constexpr int DEBUG_SAMPLE_SIZE = 530;

  for (FlopUnit& flop : flop_units_) {
    bool enable_verbose = verbose_ && (flop_count < DEBUG_SAMPLE_SIZE);

    if (enable_verbose) {
      std::cout << "\n[Sample Flop #" << flop_count << "] " << flop.inst_->getName() << std::endl;
      std::cout << "  Master: " << flop.inst_->getMaster()->getName() << std::endl;
      std::cout << "  Location: (" << flop.orig_pt_.x << ", " << flop.orig_pt_.y << ")" << std::endl;
    }

    calcFeasibleRegion(flop, enable_verbose);

    flop_count++;
  }

  if (verbose_) {
    std::cout << "\n[calcFeasibleRegions] Completed successfully." << std::endl;
    std::cout << "  - Total flops processed: " << flop_count << std::endl;
    
    if (flop_count > DEBUG_SAMPLE_SIZE) {
      std::cout << "  - Detailed output shown for first " << DEBUG_SAMPLE_SIZE << " flops" << std::endl;
    }
    std::cout << std::endl;
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

  flop_clusters_.clear();
  flop_clusters_.reserve(flop_units_.size());
  feasible_regions_.clear();
  flop_cluster_is_valid_.assign(flop_units_.size(), true);
  flop_cluster_no_further_merge_.assign(flop_units_.size(), false);

  if (verbose_) {
    std::cout << "[Step 1] Data Structure Initialization" << std::endl;
    std::cout << "  - Reserved capacity: " << flop_units_.size() << " clusters" << std::endl;
  }

  //----------------------------------------------------------------------------
  // Step 2: Create initial single-flop clusters
  //----------------------------------------------------------------------------
  for (size_t i = 0; i < flop_units_.size(); ++i) {
    flop_clusters_.emplace_back(i, flop_units_[i]);
    flop_units_[i].cluster_idx_ = i;

    const Box& feasible_region = flop_units_[i].feasible_region_;
    feasible_regions_.insert(std::make_pair(feasible_region, i));
  }

  if (verbose_) {
    std::cout << "[Step 2] Initial Cluster Creation Completed" << std::endl;
    std::cout << "  - Total flop clusters created: " << flop_clusters_.size() << std::endl;
    std::cout << "  - Feasible regions registered: " << feasible_regions_.size() << std::endl;

    //--------------------------------------------------------------------------
    // Step 3: Sample cluster details with overlap analysis
    //--------------------------------------------------------------------------
    constexpr int DEBUG_SAMPLE_SIZE = 5;
    if (!flop_clusters_.empty()) {
      std::cout << "\n[Sample Flop Clusters with Overlap Analysis]" << std::endl;
      for (size_t i = 0; i < std::min(DEBUG_SAMPLE_SIZE, static_cast<int>(flop_clusters_.size())); ++i) {
        const FlopCluster& cluster = flop_clusters_[i];
        const FlopUnit& flop = flop_units_[i];
        
        std::cout << "  Cluster #" << i << ":" << std::endl;
        std::cout << "    Flop instance: " << flop.inst_->getName() << std::endl;
        std::cout << "    Master: " << flop.inst_->getMaster()->getName() << std::endl;
        std::cout << "    Original location: (" << flop.orig_pt_.x << ", " << flop.orig_pt_.y << ")" << std::endl;
        
        // Feasible region details (45-degree rotated coordinate system)
        const Box& fr = flop.feasible_region_;
        std::cout << "    Feasible region (45° rotated): (" 
                  << fr.min_corner().get<0>() << ", " << fr.min_corner().get<1>() << ") -> ("
                  << fr.max_corner().get<0>() << ", " << fr.max_corner().get<1>() << ")" << std::endl;
        
        int fr_width = fr.max_corner().get<0>() - fr.min_corner().get<0>();
        int fr_height = fr.max_corner().get<1>() - fr.min_corner().get<1>();
        std::cout << "    FR dimensions: " << fr_width << " x " << fr_height << " (rotated coords)" << std::endl;
        
        // Count overlapping clusters using RTree query
        const std::vector<FlopClusterEntry> intersecting_clusters = getIntersectedCluster(cluster);
        
        int overlap_count = 0;
        std::vector<int> overlapping_cluster_ids;
        
        for (const auto& [other_box, other_idx] : intersecting_clusters) {
          if (static_cast<int>(other_idx) == static_cast<int>(i)) {
            continue;  // Skip self
          }
          
          overlap_count++;
          if (overlapping_cluster_ids.size() < 3) {  // Store first 3 for display
            overlapping_cluster_ids.push_back(other_idx);
          }
        }
        
        std::cout << "    Overlapping clusters: " << overlap_count << std::endl;
        if (!overlapping_cluster_ids.empty()) {
          std::cout << "      Examples: ";
          for (size_t j = 0; j < overlapping_cluster_ids.size(); ++j) {
            std::cout << "Cluster#" << overlapping_cluster_ids[j];
            if (j < overlapping_cluster_ids.size() - 1) {
              std::cout << ", ";
            }
          }
          if (overlap_count > 3) {
            std::cout << " ... (+" << (overlap_count - 3) << " more)";
          }
          std::cout << std::endl;
        }
        
        // Cluster state
        std::cout << "    Is valid: " << (flop_cluster_is_valid_[i] ? "Yes" : "No") << std::endl;
        std::cout << "    Further mergeable: " << (flop_cluster_no_further_merge_[i] ? "No" : "Yes") << std::endl;
      }
      
      if (flop_clusters_.size() > DEBUG_SAMPLE_SIZE) {
        std::cout << "  ... and " << (flop_clusters_.size() - DEBUG_SAMPLE_SIZE) 
                  << " more clusters" << std::endl;
      }
    }

    //--------------------------------------------------------------------------
    // Step 4: Feasible region overlap statistics
    //--------------------------------------------------------------------------
    if (!feasible_regions_.empty()) {
      std::cout << "\n[Feasible Region Overlap Statistics]" << std::endl;
      std::cout << "  - Total entries in feasible_regions_: " << feasible_regions_.size() << std::endl;
      
      // Calculate overlap distribution using RTree queries
      std::vector<int> overlap_counts;
      overlap_counts.reserve(feasible_regions_.size());
      
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
        
        std::cout << "  - Min overlaps per FR: " << min_overlap << std::endl;
        std::cout << "  - Max overlaps per FR: " << max_overlap << std::endl;
        std::cout << "  - Median overlaps per FR: " << median_overlap << std::endl;
        std::cout << "  - Average overlaps per FR: " << std::fixed << std::setprecision(1) << avg_overlap << std::endl;
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

  edge_pq_.clear();
  adj_list_.clear();

  if (verbose_) {
    std::cout << "[Step 1] Data Structure Initialization" << std::endl;
    std::cout << "  - Edge priority queue cleared" << std::endl;
    std::cout << "  - Adjacency list cleared" << std::endl;
  }

  //----------------------------------------------------------------------------
  // Step 2: Build edges for all clusters
  //----------------------------------------------------------------------------
  constexpr int DEBUG_SAMPLE_SIZE = 5;
  
  for (size_t i = 0; i < flop_clusters_.size(); ++i) {
    bool enable_verbose = verbose_ && (i < DEBUG_SAMPLE_SIZE);
    updateEdges(flop_clusters_[i], enable_verbose);
  }

  if (verbose_) {
    std::cout << "\n[Step 2] Edge Construction Completed" << std::endl;
    std::cout << "  - Total clusters processed: " << flop_clusters_.size() << std::endl;
    std::cout << "  - Total edges created: " << edge_pq_.size() << std::endl;
    std::cout << "  - Adjacency list entries: " << adj_list_.size() << std::endl;
    
    if (flop_clusters_.size() > DEBUG_SAMPLE_SIZE) {
      std::cout << "  - Detailed output shown for first " << DEBUG_SAMPLE_SIZE << " clusters" << std::endl;
    }

    //--------------------------------------------------------------------------
    // Step 3: Edge statistics
    //--------------------------------------------------------------------------
    if (!edge_pq_.empty()) {
      std::vector<double> edge_weights;
      edge_weights.reserve(edge_pq_.size());
      for (const auto& edge : edge_pq_) {
        edge_weights.push_back(edge.weight);
      }
      std::sort(edge_weights.begin(), edge_weights.end());

      double min_weight = edge_weights.front();
      double max_weight = edge_weights.back();
      double median_weight = edge_weights[edge_weights.size() / 2];
      double avg_weight = std::accumulate(edge_weights.begin(), edge_weights.end(), 0.0) / edge_weights.size();

      int negative_weight_count = std::count_if(edge_weights.begin(), edge_weights.end(),
                                                 [](double w) { return w < 0; });

      std::cout << "\n[Step 3] Edge Weight Distribution" << std::endl;
      std::cout << "  - Min weight (best improvement): " << std::fixed << std::setprecision(2) << min_weight << std::endl;
      std::cout << "  - Max weight: " << max_weight << std::endl;
      std::cout << "  - Median weight: " << median_weight << std::endl;
      std::cout << "  - Average weight: " << avg_weight << std::endl;
      std::cout << "  - Edges with HPWL improvement (negative): " << negative_weight_count
                << " (" << std::fixed << std::setprecision(1) << (100.0 * negative_weight_count / edge_weights.size()) << "%)" << std::endl;

      // Show top edges (best improvements)
      std::cout << "\n[Sample Top Edges (best HPWL improvement)]" << std::endl;
      constexpr int EDGE_SAMPLE_SIZE = 3;
      int edge_count = 0;
      for (const auto& edge : edge_pq_) {
        if (edge_count >= EDGE_SAMPLE_SIZE) {
          break;
        }
        
        const FlopCluster& c1 = flop_clusters_[edge.n1];
        const FlopCluster& c2 = flop_clusters_[edge.n2];
        
        std::cout << "  Edge #" << edge_count << ": Cluster[" << edge.n1 << "] <-> Cluster[" << edge.n2 << "]" << std::endl;
        std::cout << "    Weight (HPWL gain): " << std::fixed << std::setprecision(2) << edge.weight << std::endl;
        std::cout << "    Cluster[" << edge.n1 << "] size: " << c1.flops_.size() << " flop(s)" << std::endl;
        std::cout << "    Cluster[" << edge.n2 << "] size: " << c2.flops_.size() << " flop(s)" << std::endl;
        std::cout << "    Merged size would be: " << (c1.flops_.size() + c2.flops_.size()) << " flop(s)" << std::endl;
        std::cout << "    Merge position: (" << edge.pos.x << ", " << edge.pos.y << ")" << std::endl;
        
        edge_count++;
      }
      
      if (edge_pq_.size() > EDGE_SAMPLE_SIZE) {
        std::cout << "  ... and " << (edge_pq_.size() - EDGE_SAMPLE_SIZE) << " more edges" << std::endl;
      }
    }

    std::cout << "[createCompatibilityGraph] Completed successfully.\n" << std::endl;
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

  int iteration = 0;
  int total_merges = 0;
  int intermediate_merges = 0;
  int final_merges = 0;
  int skipped_invalid = 0;

  // Track cluster size distribution for debugging
  std::map<int, int> merge_size_dist;  // size -> count

  constexpr int DEBUG_SAMPLE_SIZE = std::numeric_limits<int>::max();  // Show all samples

  while (!edge_pq_.empty()) {
    // Get edge with minimum cost (most negative = best HPWL improvement)
    Edge best_edge = *edge_pq_.begin();
    edge_pq_.erase(edge_pq_.begin());

    const int n1_idx = best_edge.n1;
    const int n2_idx = best_edge.n2;

    // Skip if either cluster is already merged
    if (!flop_cluster_is_valid_[n1_idx] || !flop_cluster_is_valid_[n2_idx]) {
      skipped_invalid++;
      continue;
    }

    const FlopCluster& c1 = flop_clusters_[n1_idx];
    const FlopCluster& c2 = flop_clusters_[n2_idx];

    // Determine if we should show debug output for this merge
    const bool show_debug = verbose_ && (total_merges < DEBUG_SAMPLE_SIZE);

    if (show_debug) {
      std::cout << "\n╔═══════════════════════════════════════════════════════════════╗" << std::endl;
      std::cout << "║  MERGE #" << total_merges << " (Iteration " << iteration << ")" << std::endl;
      std::cout << "╚═══════════════════════════════════════════════════════════════╝" << std::endl;
      std::cout << "  Edge: Cluster[" << n1_idx << "] + Cluster[" << n2_idx << "]" << std::endl;
      std::cout << "    • C1: size=" << c1.flops_.size() << ", pos=(" 
                << c1.curr_pt_.x << ", " << c1.curr_pt_.y << ")" << std::endl;
      std::cout << "    • C2: size=" << c2.flops_.size() << ", pos=(" 
                << c2.curr_pt_.x << ", " << c2.curr_pt_.y << ")" << std::endl;
      std::cout << "    • Edge weight (HPWL): " << std::fixed << std::setprecision(2) 
                << best_edge.weight << (best_edge.weight < 0 ? " (improvement)" : " (degradation)") << std::endl;
      std::cout << "    • New position: (" << best_edge.pos.x << ", " << best_edge.pos.y << ")" << std::endl;
      
      // Show feasible region info
      const Box& fr1 = c1.feasible_region_;
      const Box& fr2 = c2.feasible_region_;
      std::cout << "    • C1 feasible region: ";
      if (boost::geometry::is_empty(fr1)) {
        std::cout << "EMPTY";
      } else {
        std::cout << "[(" << fr1.min_corner().get<0>() << "," << fr1.min_corner().get<1>() 
                  << ") -> (" << fr1.max_corner().get<0>() << "," << fr1.max_corner().get<1>() << ")]";
      }
      std::cout << std::endl;
      
      std::cout << "    • C2 feasible region: ";
      if (boost::geometry::is_empty(fr2)) {
        std::cout << "EMPTY";
      } else {
        std::cout << "[(" << fr2.min_corner().get<0>() << "," << fr2.min_corner().get<1>() 
                  << ") -> (" << fr2.max_corner().get<0>() << "," << fr2.max_corner().get<1>() << ")]";
      }
      std::cout << std::endl;
    }

    // Create new merged cluster
    const int new_cluster_idx = mergeClusters(best_edge, show_debug);
    const FlopCluster& new_cluster = flop_clusters_[new_cluster_idx];

    if (show_debug) {
      std::cout << "  ▸ New Cluster[" << new_cluster_idx << "]: size=" << new_cluster.flops_.size() << std::endl;
      const Box& new_fr = new_cluster.feasible_region_;
      std::cout << "    • Feasible region: ";
      if (boost::geometry::is_empty(new_fr)) {
        std::cout << "EMPTY (WARNING!)" << std::endl;
      } else {
        std::cout << "[(" << new_fr.min_corner().get<0>() << "," << new_fr.min_corner().get<1>() 
                  << ") -> (" << new_fr.max_corner().get<0>() << "," << new_fr.max_corner().get<1>() << ")]" << std::endl;
      }
    }

    total_merges++;
    merge_size_dist[new_cluster.flops_.size()]++;

    // Check if cluster can grow further
    const bool can_merge_further = isFurtherMergeable(new_cluster, show_debug);

    if (can_merge_further) {
      // Intermediate merge: can still grow larger
      intermediate_merges++;
      
      if (show_debug) {
        std::cout << "  ▸ Decision: INTERMEDIATE (can grow to larger bit-width)" << std::endl;
        std::cout << "    • Action: Update edges only" << std::endl;
      }

      updateEdges(new_cluster, show_debug);
      
      if (show_debug) {
        const int new_edge_count = adj_list_.count(new_cluster_idx) ? adj_list_[new_cluster_idx].size() : 0;
        std::cout << "    • New edges created: " << new_edge_count << std::endl;
        
        if (new_edge_count > 0 && new_edge_count <= 5) {
          std::cout << "    • Top neighbors:" << std::endl;
          int shown = 0;
          for (const Edge& e : adj_list_[new_cluster_idx]) {
            if (shown >= 3) break;
            int neighbor_idx = (e.n1 == new_cluster_idx) ? e.n2 : e.n1;
            std::cout << "      - Cluster[" << neighbor_idx << "]: size=" 
                      << flop_clusters_[neighbor_idx].flops_.size() 
                      << ", weight=" << std::fixed << std::setprecision(2) << e.weight << std::endl;
            shown++;
          }
        }
      }
    } 
    else {
      // Final merge: reached maximum bit-width
      final_merges++;
      flop_cluster_no_further_merge_[new_cluster_idx] = true;

      if (show_debug) {
        std::cout << "  ▸ Decision: FINAL (reached maximum bit-width)" << std::endl;
        std::cout << "    • Action: Distribute slack to neighbors" << std::endl;
      }

      // Distribute unused slack to connected flops
      const std::set<int> affected_flop_units = distributeSlack(new_cluster, show_debug);
      
      std::set<int> affected_flop_clusters;
      
      // Update feasible regions of affected flop units
      for (int u_idx : affected_flop_units) {
        calcFeasibleRegion(flop_units_[u_idx]);
        
        // Collect affected clusters for later updates
        const int c_idx = flop_units_[u_idx].cluster_idx_;
        if (flop_cluster_is_valid_[c_idx] && !flop_cluster_no_further_merge_[c_idx]) {
          affected_flop_clusters.insert(c_idx);
        }
      }

      if (show_debug) {
        std::cout << "    • Affected flops: " << affected_flop_units.size() << std::endl;
        std::cout << "    • Affected clusters: " << affected_flop_clusters.size() << std::endl;
      }

      // Update feasible regions of affected clusters
      for (int c_idx : affected_flop_clusters) {
        updateFeasibleRegion(flop_clusters_[c_idx], show_debug);
      }

      // Update edges of affected clusters
      int total_new_edges = 0;
      for (int c_idx : affected_flop_clusters) {
        int before_count = adj_list_.count(c_idx) ? adj_list_[c_idx].size() : 0;
        updateEdges(flop_clusters_[c_idx]);
        int after_count = adj_list_.count(c_idx) ? adj_list_[c_idx].size() : 0;
        total_new_edges += (after_count - before_count);
      }

      if (show_debug) {
        std::cout << "    • Total new edges from affected clusters: " << total_new_edges << std::endl;
      }
    }

    iteration++;

    // Periodic progress report
    if (verbose_ && iteration % 100 == 0) {
      std::cout << "\n  [Progress] Iteration " << iteration << std::endl;
      std::cout << "    - Total merges: " << total_merges 
                << " (intermediate: " << intermediate_merges 
                << ", final: " << final_merges << ")" << std::endl;
      std::cout << "    - Skipped (invalid): " << skipped_invalid << std::endl;
      std::cout << "    - Remaining edges: " << edge_pq_.size() << std::endl;
      std::cout << "    - Valid clusters: " << std::count(flop_cluster_is_valid_.begin(), 
                                                          flop_cluster_is_valid_.end(), true) << std::endl;
      
      // Show current size distribution
      std::cout << "    - Current merge size distribution:" << std::endl;
      for (const auto& [size, count] : merge_size_dist) {
        std::cout << "      * " << size << "-bit: " << count << " merge(s)" << std::endl;
      }
    }
  }

  if (verbose_) {
    std::cout << "\n╔═══════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║  FINAL CLUSTERING STATISTICS                                  ║" << std::endl;
    std::cout << "╚═══════════════════════════════════════════════════════════════╝" << std::endl;
    std::cout << "  Execution Summary:" << std::endl;
    std::cout << "    • Total iterations: " << iteration << std::endl;
    std::cout << "    • Total merges: " << total_merges 
              << " (intermediate: " << intermediate_merges 
              << ", final: " << final_merges << ")" << std::endl;
    std::cout << "    • Skipped (invalid): " << skipped_invalid << std::endl;

    // Count final valid clusters by size
    std::map<int, int> cluster_size_dist;
    std::map<int, std::vector<int>> size_to_clusters;  // For showing examples
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

    std::cout << "\n  Final Cluster Status:" << std::endl;
    std::cout << "    • Total valid clusters: " << total_valid << std::endl;
    std::cout << "    • Single-flop clusters: " << single_flop_clusters 
              << " (" << std::fixed << std::setprecision(1) 
              << (100.0f * single_flop_clusters / total_valid) << "%)" << std::endl;
    std::cout << "    • Multi-bit clusters (MBFF): " << merged_clusters 
              << " (" << (100.0f * merged_clusters / total_valid) << "%)" << std::endl;
    std::cout << "    • Flops in MBFF: " << total_flops_in_mbff << " / " << total_flops 
              << " (" << mbff_ratio << "%)" << std::endl;

    if (!cluster_size_dist.empty()) {
      std::cout << "\n  Cluster Size Distribution:" << std::endl;
      int max_size = 0;
      for (const auto& [size, count] : cluster_size_dist) {
        max_size = std::max(max_size, size);
        std::cout << "    • " << size << "-bit: " << count << " cluster(s)";
        
        // Show percentage for merged clusters
        if (size > 1) {
          float pct = 100.0f * count / merged_clusters;
          std::cout << " (" << std::fixed << std::setprecision(1) << pct << "% of MBFF)";
        }
        std::cout << std::endl;
      }
      
      std::cout << "    • Maximum bit-width achieved: " << max_size << std::endl;
    }

    // Show why clustering stopped
    std::cout << "\n  Termination Analysis:" << std::endl;
    if (edge_pq_.empty()) {
      std::cout << "    • Reason: Edge queue exhausted (no more mergeable pairs)" << std::endl;
      
      // Analyze remaining clusters
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
      std::cout << "    • Remaining clusters that could grow: " << remaining_intermediate << std::endl;
      std::cout << "    • Remaining final clusters: " << remaining_final << std::endl;
      
      if (remaining_intermediate > 0) {
        std::cout << "    • ⚠ WARNING: " << remaining_intermediate 
                  << " cluster(s) marked as intermediate but no edges found!" << std::endl;
        std::cout << "      This suggests edge creation or compatibility issues." << std::endl;
      }
    }

    // Show sample merged clusters
    if (merged_clusters > 0) {
      std::cout << "\n  Sample Multi-bit Clusters:" << std::endl;
      
      // Show largest clusters first
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
          int c_idx = cluster_indices[i];
          const FlopCluster& cluster = flop_clusters_[c_idx];
          
          std::cout << "    [" << (samples_shown + 1) << "] Cluster #" << c_idx 
                    << ": " << size << "-bit MBFF" << std::endl;
          std::cout << "        Position: (" << cluster.curr_pt_.x << ", " << cluster.curr_pt_.y << ")" << std::endl;
          
          const Box& fr = cluster.feasible_region_;
          std::cout << "        Feasible region: ";
          if (boost::geometry::is_empty(fr)) {
            std::cout << "EMPTY";
          } else {
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

    std::cout << "\n" << std::endl;
  }
}

//==============================================================================
// Phase 11: implementClusters()
//==============================================================================

void 
AggloCluster::implementClusters()
{
  // 1. 최종 클러스터 식별
  std::vector<int> final_cluster_indices;
  for (int i = 0; i < flop_clusters_.size(); ++i) {
    // 유효하고(valid), 병합된(size > 1) 클러스터
    if (flop_cluster_is_valid_[i] && flop_clusters_[i].flops_.size() > 1) {
      final_cluster_indices.push_back(i);
    }
  }

  // 2. 각 클러스터 구현
  for (int cluster_idx : final_cluster_indices) {
    implementSingleCluster(flop_clusters_[cluster_idx]);
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
  const sta::LibertyCell* lib_cell = getLibertyCell(inst);
  if (lib_cell == nullptr || lib_cell->sequentials().empty()) {
    return false;
  }

  if (lib_cell->isClockGate() || resizer_->dontUse(lib_cell)
      || std::any_of(
          lib_cell->sequentials().begin(),
          lib_cell->sequentials().end(),
          [](const sta::Sequential* seq) { return seq->isLatch(); })) {
    return false;
  }

  int d_pins = 0;
  int q_pins = 0;
  int qn_pins = 0;
  std::vector<std::string> unrecognized_pins;

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
    } else {
      // Unrecognized pin type
      unrecognized_pins.push_back(iterm->getMTerm()->getName());
    }
  }

  // Debug output for invalid flops
  if (!unrecognized_pins.empty()) {
    if (verbose_) {
      std::cout << "[isValidFlop] Invalid flop: " << inst->getName() 
                << " (Master: " << inst->getMaster()->getName() << ")" << std::endl;
      std::cout << "  Reason: Unrecognized pin(s): ";
      for (size_t i = 0; i < unrecognized_pins.size(); ++i) {
        std::cout << unrecognized_pins[i];
        if (i < unrecognized_pins.size() - 1) {
          std::cout << ", ";
        }
      }
      std::cout << std::endl;
    }
    return false;
  }

  if (d_pins != std::max(q_pins, qn_pins)) {
    if (verbose_) {
      std::cout << "[isValidFlop] Invalid flop: " << inst->getName()
                << " (Master: " << inst->getMaster()->getName() << ")" << std::endl;
      std::cout << "  Reason: Pin count mismatch - D:" << d_pins 
                << ", Q:" << q_pins << ", QN:" << qn_pins << std::endl;
    }
    return false;
  }

  return true;
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

  for (odb::dbITerm* iterm : inst->getITerms()) {
    odb::dbNet* net = iterm->getNet();
    if (net == nullptr) {
      continue;
    }

    if (isClockPin(iterm)) {
      if (clock_net == nullptr) {
        clock_net = net;
      } else if (clock_net != net) {
        if (verbose_) {
          std::cout << "[createInstMask] Instance " << inst->getName() 
                    << " has multiple different clock nets. Using the first one found (" 
                    << clock_net->getName() << ")." << std::endl;
        }
      }
    } else if (isClearPin(iterm)) {
      if (clear_net == nullptr) {
        clear_net = net;
      } else if (clear_net != net) {
        if (verbose_) {
          std::cout << "[createInstMask] Instance " << inst->getName()
                    << " has multiple different clear nets. Using the first one found ("
                    << clear_net->getName() << ")." << std::endl;
        }
      }
    } else if (isPresetPin(iterm)) {
      if (preset_net == nullptr) {
        preset_net = net;
      } else if (preset_net != net) {
        if (verbose_) {
          std::cout << "[createInstMask] Instance " << inst->getName()
                    << " has multiple different preset nets. Using the first one found ("
                    << preset_net->getName() << ")." << std::endl;
        }
      }
    } else if (isScanEnablePin(iterm)) {
      if (scan_enable_net == nullptr) {
        scan_enable_net = net;
      } else if (scan_enable_net != net) {
        if (verbose_) {
          std::cout << "[createInstMask] Instance " << inst->getName()
                    << " has multiple different scan enable nets. Using the first one found ("
                    << scan_enable_net->getName() << ")." << std::endl;
        }
      }
    } else if (isScanInPin(iterm)) {
      if (scan_in_net == nullptr) {
        scan_in_net = net;
      } else if (scan_in_net != net) {
        if (verbose_) {
          std::cout << "[createInstMask] Instance " << inst->getName()
                    << " has multiple different scan in nets. Using the first one found ("
                    << scan_in_net->getName() << ")." << std::endl;
        }
      }
    }
  }

  return InstMask(
      clock_net, clear_net, preset_net, scan_enable_net, scan_in_net);
}

MasterMask 
AggloCluster::createMasterMask(odb::dbInst* inst)
{
  const sta::LibertyCell* lib_cell = getLibertyCell(inst);
  if (lib_cell == nullptr) {
    return MasterMask();
  }

  const sta::FuncExpr* expr = nullptr;
  if (!lib_cell->sequentials().empty()) {
    expr = lib_cell->sequentials().front()->data();
  }
  if (expr == nullptr) {
    return MasterMask();
  }

  int func_id = getFuncId(expr, inst);
  bool clk = hasPositiveClockEdge(inst);
  bool clear = hasClear(inst);
  bool preset = hasPreset(inst);
  bool scan = hasScan(inst);
  bool qpin = false;
  bool qnpin = false;
  for (odb::dbITerm* iterm : inst->getITerms()) {
    if (!qpin && isQPin(iterm)) {
      qpin = true;
    }
    if (!qnpin && isQNPin(iterm)) {
      qnpin = true;
    }
    if (qpin && qnpin) {
      break;
    }
  }
  return MasterMask(func_id, clk, clear, preset, qpin, qnpin, scan);
}

// -----------------------------------------------------------------------------
// [Level 2] Liberty Cell and Function Analysis
// -----------------------------------------------------------------------------

const sta::LibertyCell* 
AggloCluster::getLibertyCell(odb::dbInst* inst) const
{
  const sta::Cell* cell = network_->dbToSta(inst->getMaster());
  const sta::LibertyCell* lib_cell = network_->libertyCell(cell);
  if (!lib_cell) {
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
  if (expr == nullptr) {
    return "()";
  }

  if (expr->op() == sta::FuncExpr::op_port) {
    FlopPort p_type = getPortType(expr->port(), inst);
    return "p(" + std::to_string(int(p_type)) + ")";
  }

  std::string op_str = std::to_string(expr->op());

  if (expr->op() == sta::FuncExpr::op_not) {
    return "op" + op_str + "(" + getFuncStr(expr->left(), inst) + ")";
  }

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
  if (lib_cell == nullptr) {
    return false;
  }
  for (sta::Sequential* seq : lib_cell->sequentials()) {
    if (seq->clear()) {
      return true;
    }
  }
  return false;
}

bool 
AggloCluster::hasPreset(odb::dbInst* inst) const
{
  const sta::LibertyCell* lib_cell = getLibertyCell(inst);
  if (lib_cell == nullptr) {
    return false;
  }
  for (sta::Sequential* seq : lib_cell->sequentials()) {
    if (seq->preset()) {
      return true;
    }
  }
  return false;
}

bool 
AggloCluster::hasPositiveClockEdge(odb::dbInst* inst) const
{
  const sta::LibertyCell* lib_cell = getLibertyCell(inst);
  if (lib_cell == nullptr) {
    return false;
  }

  for (sta::Sequential* seq : lib_cell->sequentials()) {
    const sta::FuncExpr* left = seq->clock()->left();
    const sta::FuncExpr* right = seq->clock()->right();
    if (left && !right) {
      return false;  // !CLK - negative edge
    }
  }
  return true;
}

bool 
AggloCluster::hasScan(odb::dbInst* inst) const
{
  const sta::LibertyCell* lib_cell = getLibertyCell(inst);
  if (lib_cell == nullptr) {
    return false;
  }
  return lib_cell && getLibertyScanIn(lib_cell)
         && getLibertyScanEnable(lib_cell);
}

// -----------------------------------------------------------------------------
// [Level 2] Pin Counting Utilities
// -----------------------------------------------------------------------------

int 
AggloCluster::getNumDPins(odb::dbInst* inst) const
{
  int num = 0;
  for (odb::dbITerm* iterm : inst->getITerms()) {
    num += isDPin(iterm);
  }
  return num;
}

int 
AggloCluster::getNumQPins(odb::dbInst* inst) const
{
  int num = 0;
  for (odb::dbITerm* iterm : inst->getITerms()) {
    num += isQPin(iterm);
  }
  return num;
}

int 
AggloCluster::getNumQNPins(odb::dbInst* inst) const
{
  int num = 0;
  for (odb::dbITerm* iterm : inst->getITerms()) {
    num += isQNPin(iterm);
  }
  return num;
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
  if (lib_port == nullptr) {
    return false;
  }
  const sta::LibertyCell* lib_cell = getLibertyCell(iterm->getInst());
  if (lib_cell == nullptr) {
    return false;
  }
  for (const sta::Sequential* seq : lib_cell->sequentials()) {
    if (seq->preset() && seq->preset()->hasPort(lib_port)) {
      return true;
    }
  }
  return false;
}

bool 
AggloCluster::isClearPin(odb::dbITerm* iterm) const
{
  const sta::LibertyPort* lib_port = getLibertyPort(iterm);
  if (lib_port == nullptr) {
    return false;
  }
  const sta::LibertyCell* lib_cell = getLibertyCell(iterm->getInst());
  if (lib_cell == nullptr) {
    return false;
  }
  for (const sta::Sequential* seq : lib_cell->sequentials()) {
    if (seq->clear() && seq->clear()->hasPort(lib_port)) {
      return true;
    }
  }
  return false;
}

bool 
AggloCluster::isQPin(odb::dbITerm* iterm) const
{
  if (iterm->getIoType() != odb::dbIoType::OUTPUT) {
    return false;
  }
  if (isClockPin(iterm) || isPowerPin(iterm)) {
    return false;
  }

  const sta::LibertyPort* lib_port = getLibertyPort(iterm);
  if (lib_port == nullptr) {
    return false;
  }

  const sta::FuncExpr* func = getPortFunc(lib_port);
  if (func == nullptr) {
    return false;
  }

  const sta::LibertyPort* funcPort = func->port();
  if (funcPort == nullptr) {
    return false;
  }

  const std::string pinFuncName = funcPort->name();

  const sta::LibertyCell* lib_cell = getLibertyCell(iterm->getInst());
  if (lib_cell == nullptr) {
    return false;
  }

  for (const sta::Sequential* seq : lib_cell->sequentials()) {
    const sta::LibertyPort* output = seq->output();
    if (output != nullptr && output->name() == pinFuncName) {
      return true;
    }
  }

  return false;
}

bool 
AggloCluster::isQNPin(odb::dbITerm* iterm) const
{
  if (iterm->getIoType() != odb::dbIoType::OUTPUT) {
    return false;
  }
  if (isClockPin(iterm) || isPowerPin(iterm)) {
    return false;
  }

  const sta::LibertyPort* lib_port = getLibertyPort(iterm);
  if (lib_port == nullptr) {
    return false;
  }

  const sta::FuncExpr* func = getPortFunc(lib_port);
  if (func == nullptr) {
    return false;
  }

  const sta::LibertyPort* funcPort = func->port();
  if (funcPort == nullptr) {
    return false;
  }

  const std::string pinFuncName = funcPort->name();

  const sta::LibertyCell* lib_cell = getLibertyCell(iterm->getInst());
  if (lib_cell == nullptr) {
    return false;
  }

  for (const sta::Sequential* seq : lib_cell->sequentials()) {
    const sta::LibertyPort* outputInv = seq->outputInv();
    if (outputInv != nullptr && outputInv->name() == pinFuncName) {
      return true;
    }
  }

  return false;
}

bool 
AggloCluster::isScanInPin(odb::dbITerm* iterm) const
{
  const sta::LibertyCell* lib_cell = getLibertyCell(iterm->getInst());
  if (lib_cell == nullptr) {
    return false;
  }
  if (lib_cell && sta::getLibertyScanIn(lib_cell)) {
    odb::dbMTerm* mterm = network_->staToDb(sta::getLibertyScanIn(lib_cell));
    return iterm->getInst()->getITerm(mterm) == iterm;
  }
  return false;
}

bool 
AggloCluster::isScanEnablePin(odb::dbITerm* iterm) const
{
  const sta::LibertyCell* lib_cell = getLibertyCell(iterm->getInst());
  if (lib_cell == nullptr) {
    return false;
  }
  if (lib_cell && sta::getLibertyScanEnable(lib_cell)) {
    odb::dbMTerm* mterm
        = network_->staToDb(sta::getLibertyScanEnable(lib_cell));
    return iterm->getInst()->getITerm(mterm) == iterm;
  }
  return false;
}

// -----------------------------------------------------------------------------
// [Level 3] Port Analysis Utilities
// -----------------------------------------------------------------------------

const sta::LibertyPort* 
AggloCluster::getLibertyPort(odb::dbITerm* iterm) const
{
  const sta::Pin* pin = network_->dbToSta(iterm);
  if (pin == nullptr) {
    return nullptr;
  }
  const sta::LibertyPort* lib_port = network_->libertyPort(pin);
  if (lib_port == nullptr) {
    return nullptr;
  }
  return lib_port;
}

const sta::FuncExpr* 
AggloCluster::getPortFunc(const sta::LibertyPort* port) const
{
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
  return nullptr;
}

FlopPort 
AggloCluster::getPortType(const sta::LibertyPort* lib_port, odb::dbInst* inst) const
{
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
      std::cout << "[getPortType] Could not find Liberty cell for instance: " 
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
    std::cout << "[getPortType] Could not recognize port: " << lib_port->name() << std::endl;
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
  // Save original precision settings
  std::streamsize original_precision = std::cout.precision();
  std::ios_base::fmtflags original_flags = std::cout.flags();
  
  if (verbose) {
    std::cout << std::fixed << std::setprecision(18);
  }

  const auto est = resizer_->getEstimateParasitics();
  const double unit_c = est->wireSignalCapacitance(corner_);
  const double unit_r = est->wireSignalResistance(corner_);

  if (verbose) {
    std::cout << "\n=== Calculating Feasible Region for Flop: " << flop.inst_->getName() << " ===" << std::endl;
    std::cout << "  [Step 1] Wire Parameters" << std::endl;
    std::cout << "    Unit resistance: " << unit_r << std::endl;
    std::cout << "    Unit capacitance: " << unit_c << std::endl;
  }

  // Classify pins
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
  
  if (!clk_pin) {
    if (verbose) {
      std::cout << "    [WARNING] No clock pin found - skipping this flop" << std::endl;
    }
    return;
  }
  
  // Process output pins
  std::vector<odb::dbITerm*> output_pins;
  output_pins.reserve(q_pins.size() + qn_pins.size());
  output_pins.insert(output_pins.end(), q_pins.begin(), q_pins.end());
  output_pins.insert(output_pins.end(), qn_pins.begin(), qn_pins.end());
  
  for (odb::dbITerm* out_pin : output_pins) {
    processFanOutPin(flop, out_pin, clk_pin->getMTerm(), est, unit_r, unit_c, verbose);
  }

  // Process input pins
  for (odb::dbITerm* d_pin : d_pins) {
    processFanInPin(flop, d_pin, est, unit_r, unit_c, verbose);
  }
  
  // Compute final feasible region from all pins
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
  
  // Check timing constraints
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
  
  // Get most critical path
  const auto& [critical_path_idx, slack_budget] = budget_it->second[0];
  
  if (verbose) {
    std::cout << "    Most critical path:" << std::endl;
    std::cout << "      Path index: " << critical_path_idx << std::endl;
    std::cout << "      Slack budget: " << slack_budget << std::endl;
  }
  
  // Validate net connection
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
  
  // Extract cell delay characterization
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
  
  // Build Steiner tree and compute parasitics
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
  
  // Solve for maximum distance using coefficient segments
  float max_dist = -1.0f;  // Invalid value if no valid segment found
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
      // Special handling for last iteration: use candidate_dist even if out of range
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

  // Apply bin size constraint: ensure max_dist is in [l1, l1 + bin_size * k]
  const double bin_size_x = virtual_bin_grid_.getBinSizeX();
  const double bin_size_y = virtual_bin_grid_.getBinSizeY();
  const double avg_bin_size = (bin_size_x + bin_size_y) / 2.0;
  const float l1_dbu = metersToDbu(l1);
  const float lower_limit = l1_dbu;
  const float upper_limit = l1_dbu + (avg_bin_size * feasible_region_bin_multiplier_);
  
  if (verbose) {
    std::cout << "      Applying distance constraint:" << std::endl;
    std::cout << "        Calculated max_dist: " << max_dist << " DBU" << std::endl;
    std::cout << "        Lower bound (l1): " << lower_limit << " DBU" << std::endl;
    std::cout << "        Bin size (avg): " << avg_bin_size << " DBU" << std::endl;
    std::cout << "        Bin multiplier: " << feasible_region_bin_multiplier_ << std::endl;
    std::cout << "        Upper bound (l1 + bin*k): " << upper_limit << " DBU" << std::endl;
  }
  
  max_dist = std::clamp(max_dist, lower_limit, upper_limit);

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
  
  // Check timing constraints
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
  
  // Get most critical path
  const auto& [critical_path_idx, slack_budget] = budget_it->second[0];
  
  if (verbose) {
    std::cout << "    Most critical path:" << std::endl;
    std::cout << "      Path index: " << critical_path_idx << std::endl;
    std::cout << "      Slack budget: " << slack_budget << std::endl;
  }
  
  // Validate net connection
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
  
  // Get driver pin
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
  
  // Calculate D pin capacitance
  const sta::Pin* ipin_sta = network_->dbToSta(d_pin);
  const float ipin_cap = getPinCapacitance(ipin_sta);
  
  if (verbose) {
    std::cout << "    D pin input capacitance: " << ipin_cap << std::endl;
  }
  
  // Extract cell delay characterization
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
  
  // Calculate coefficients
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
  
  // Build Steiner tree
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
  
  // Find path in Steiner tree
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
  
  // Calculate path segments
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
  
  // Solve for maximum distance using coefficient segments
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
  
  float max_dist = -1.0f;  // Invalid value if no valid segment found
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
      // Special handling for last iteration: use candidate_dist even if out of range
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

  // Apply bin size constraint: ensure max_dist is in [l1, l1 + bin_size * k]
  const double bin_size_x = virtual_bin_grid_.getBinSizeX();
  const double bin_size_y = virtual_bin_grid_.getBinSizeY();
  const double avg_bin_size = (bin_size_x + bin_size_y) / 2.0;
  const float l1_dbu = metersToDbu(l1);
  const float lower_limit = l1_dbu;
  const float upper_limit = l1_dbu + (avg_bin_size * feasible_region_bin_multiplier_);
  
  if (verbose) {
    std::cout << "      Applying distance constraint:" << std::endl;
    std::cout << "        Calculated max_dist: " << max_dist << " DBU" << std::endl;
    std::cout << "        Lower bound (l1): " << lower_limit << " DBU" << std::endl;
    std::cout << "        Bin size (avg): " << avg_bin_size << " DBU" << std::endl;
    std::cout << "        Bin multiplier: " << feasible_region_bin_multiplier_ << std::endl;
    std::cout << "        Upper bound (l1 + bin*k): " << upper_limit << " DBU" << std::endl;
  }
  
  max_dist = std::clamp(max_dist, lower_limit, upper_limit);

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
    
    // Skip inverse boxes (no constraints)
    if (pin_box.min_corner().get<0>() > pin_box.max_corner().get<0>() ||
        pin_box.min_corner().get<1>() > pin_box.max_corner().get<1>()) {
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
  
  flop.feasible_region_ = final_region;
  
  if (verbose) {
    std::cout << "    Final feasible region: [(" 
              << final_region.min_corner().get<0>() << ", " 
              << final_region.min_corner().get<1>() << ") - (" 
              << final_region.max_corner().get<0>() << ", " 
              << final_region.max_corner().get<1>() << ")]" << std::endl;
  }
}

// -----------------------------------------------------------------------------
// [Level 3] Utility functions for timing and parasitic analysis
// -----------------------------------------------------------------------------

float
AggloCluster::getPinCapacitance(const sta::Pin* pin) const
{
  // Get the Liberty port for the pin
  sta::LibertyPort* port = network_->libertyPort(pin);
  if (!port) {
    return 0.0f;  // No Liberty model available
  }
  
  // Get corner-specific port (corner 0 = default corner)
  sta::LibertyPort* corner_port = port->cornerPort(0);
  if (!corner_port) {
    return 0.0f;  // No corner-specific data
  }
  
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
  
  // Get timing graph and analysis parameters
  sta::Graph* timing_graph = sta_->graph();
  const sta::MinMax* worst_case = sta::MinMax::max();

  // Iterate through all terminals of the instance
  for (odb::dbITerm* terminal : inst->getITerms()) {
    // Only process input pins (skip outputs and power pins)
    if (terminal->getIoType() != odb::dbIoType::INPUT) {
      continue;
    }

    // Get slew value for this input pin
    float slew = getTerminalSlew(terminal, timing_graph, worst_case);
    
    // Only include pins with valid (non-zero) slew values
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
  // Base case: Invalid node
  if (current_pt == -1) {
    return false;
  }

  // Add current node to path
  path.push_back(current_pt);

  // Base case: Found target
  if (current_pt == target_pt) {
    return true;
  }

  // Recursive case: Search left subtree
  if (findSteinerPathRecursive(tree, tree->left(current_pt), target_pt, path)) {
    return true;
  }

  // Recursive case: Search right subtree
  if (findSteinerPathRecursive(tree, tree->right(current_pt), target_pt, path)) {
    return true;
  }

  // Backtrack: Remove current node from path if target not found in subtrees
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
  // Solve quadratic equation: a*x^2 + b*x + c = 0
  // where x is the new Manhattan distance from Steiner point to flop (in meters)
  auto a = unit_r * unit_c;
  auto b = wo_fst_stt_cap * unit_r + coeff * unit_c;
  auto c = -std::pow(l1, 2) * unit_r * unit_c 
           - slack_budget 
           - wo_fst_stt_cap * l1 * unit_r 
           - coeff * l1 * unit_c;
  
  float D = std::pow(b, 2) - 4 * a * c;
  
  if (verbose) {
    std::cout << "    Solving quadratic equation for maximum distance (FanOut):" << std::endl;
    std::cout << "      Coefficients: a = " << a << ", b = " << b << ", c = " << c << std::endl;
    std::cout << "      Discriminant: " << D << std::endl;
  }
  
  if (D < 0) {
    if (verbose) {
      std::cout << "      No real roots - discriminant is negative" << std::endl;
      std::cout << "      Returning invalid value (-1.0)" << std::endl;
    }
    return -1.0f;
  }
  
  // Calculate roots using quadratic formula (result in meters)
  float root1 = (-b + sqrt(D)) / (2 * a);
  float root2 = (-b - sqrt(D)) / (2 * a);
  
  if (verbose) {
    std::cout << "      Root 1: " << root1 << std::endl;
    std::cout << "      Root 2: " << root2 << std::endl;
  }
  
  // Use larger positive root
  if (root1 > 0 || root2 > 0) {
    float max_dist = std::max(root1, root2);
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
  // Solve quadratic equation for Fan-In case
  // where x is the new Manhattan distance from Steiner point to flop (in meters)
  auto a = unit_r * unit_c;
  auto b = on_path_R_wo_last * unit_c + coeff * unit_c + unit_r * ipin_cap;
  auto c = -pow(l1, 2) * unit_r * unit_c 
           - on_path_R_wo_last * unit_c * l1 
           - coeff * l1 * unit_c 
           - ipin_cap * l1 * unit_r
           - slack_budget;
  
  float D = pow(b, 2) - 4 * a * c;
  
  if (verbose) {
    std::cout << "    Solving quadratic equation for maximum distance (FanIn):" << std::endl;
    std::cout << "      Coefficients: a = " << a << ", b = " << b << ", c = " << c << std::endl;
    std::cout << "      Discriminant: " << D << std::endl;
  }
  
  if (D < 0) {
    if (verbose) {
      std::cout << "      No real roots - discriminant is negative" << std::endl;
      std::cout << "      Returning invalid value (-1.0)" << std::endl;
    }
    return -1.0f;
  }
  
  // Calculate roots using quadratic formula (result in meters)
  float root1 = (-b + sqrt(D)) / (2 * a);
  float root2 = (-b - sqrt(D)) / (2 * a);
  
  if (verbose) {
    std::cout << "      Root 1: " << root1 << std::endl;
    std::cout << "      Root 2: " << root2 << std::endl;
  }
  
  // Use larger positive root
  if (root1 > 0 || root2 > 0) {
    float max_dist = std::max(root1, root2);
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
  // Create 4 corner points at max_dist from Steiner point (Manhattan distance)
  Point p1(steiner_x - static_cast<int>(max_dist), steiner_y);
  Point p2(steiner_x + static_cast<int>(max_dist), steiner_y);
  Point p3(steiner_x, steiner_y - static_cast<int>(max_dist));
  Point p4(steiner_x, steiner_y + static_cast<int>(max_dist));
  
  // Transform to 45-degree rotated coordinates
  Point t1 = transformCoords(p1);
  Point t2 = transformCoords(p2);
  Point t3 = transformCoords(p3);
  Point t4 = transformCoords(p4);
  
  // Find bounding box in transformed coordinates
  int min_x = std::min({t1.get<0>(), t2.get<0>(), t3.get<0>(), t4.get<0>()});
  int max_x = std::max({t1.get<0>(), t2.get<0>(), t3.get<0>(), t4.get<0>()});
  int min_y = std::min({t1.get<1>(), t2.get<1>(), t3.get<1>(), t4.get<1>()});
  int max_y = std::max({t1.get<1>(), t2.get<1>(), t3.get<1>(), t4.get<1>()});
  
  Point box_min(min_x, min_y);
  Point box_max(max_x, max_y);
  
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
  // Search through all timing arcs in the cell
  for (sta::TimingArcSet* arc_set : liberty_cell->timingArcSets()) {
    const sta::LibertyPort* from_port = arc_set->from();
    const sta::LibertyPort* to_port = arc_set->to();
    
    if (!from_port || !to_port) {
      continue;
    }
    
    // Check if this arc matches the requested input->output path
    if (std::string(from_port->name()) != input_pin_name || 
        std::string(to_port->name()) != output_pin_name) {
      continue;
    }

    // Found matching arc set, now look for the gate table model
    for (sta::TimingArc* arc : arc_set->arcs()) {
      auto* gate_model = dynamic_cast<sta::GateTableModel*>(arc->model());
      if (!gate_model) {
        continue;
      }

      // Extract the delay model and find the output capacitance axis
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
  auto* delay_model = gate_model->delayModel();
  
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
  const bool pocv_enabled = false;  // POCV (Parametric On-Chip Variation) disabled
  sta::Slew output_slew;  // Output slew (not used but required by API)
  sta::ArcDelay arc_delay;
  
  // Query the timing model for delay at this capacitance point
  timing_model->gateDelay(pvt_conditions, 
                         input_slew, 
                         output_capacitance, 
                         pocv_enabled, 
                         arc_delay, 
                         output_slew);

  // Convert to float and ensure non-negative
  float delay = static_cast<float>(arc_delay);
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
  // Convert physical pin to timing analysis pin
  sta::Pin* sta_pin = network_->dbToSta(terminal);
  
  // Get the load vertex (input side of the pin where signal arrives)
  sta::Vertex* load_vertex = timing_graph->pinLoadVertex(sta_pin);
  
  // Extract slew at this vertex (worst-case transition time)
  return sta_->vertexSlew(load_vertex, min_max);
}

//==============================================================================
// Phase 9 Helper Functions: Compatibility Graph Generation
//==============================================================================

// -----------------------------------------------------------------------------
// [Level 1] Top-level: Build/update compatibility graph edges
// -----------------------------------------------------------------------------

void 
AggloCluster::updateEdges(const FlopCluster& cluster, bool verbose)
{
  const int i = cluster.id_;

  if (verbose) {
    std::cout << "\n[updateEdges] Processing Cluster #" << i << std::endl;
    std::cout << "  Cluster size: " << cluster.flops_.size() << " flop(s)" << std::endl;
    std::cout << "  Current position: (" << cluster.curr_pt_.x << ", " << cluster.curr_pt_.y << ")" << std::endl;
  }

  if (!flop_cluster_is_valid_[i] || flop_cluster_no_further_merge_[i]) {
    if (verbose) {
      std::cout << "  Status: Skipped (";
      if (!flop_cluster_is_valid_[i]) {
        std::cout << "invalid";
      } else {
        std::cout << "no further merge";
      }
      std::cout << ")" << std::endl;
    }
    return;
  }

  if (boost::geometry::is_empty(cluster.feasible_region_)) {
    if (verbose) {
      std::cout << "  Status: Skipped (empty feasible region)" << std::endl;
      std::cout << "  ⚠ This cluster cannot create edges - feasible region is empty!" << std::endl;
    }
    return;
  }

  if (verbose) {
    const Box& fr = cluster.feasible_region_;
    Point min_uv = fr.min_corner();
    Point max_uv = fr.max_corner();
    int u_min = min_uv.get<0>();
    int v_min = min_uv.get<1>();
    int u_max = max_uv.get<0>();
    int v_max = max_uv.get<1>();
    
    // Four corners of the UV box
    Point corner1_xy = inverseTransformCoords(Point(u_min, v_min));
    Point corner2_xy = inverseTransformCoords(Point(u_max, v_min));
    Point corner3_xy = inverseTransformCoords(Point(u_max, v_max));
    Point corner4_xy = inverseTransformCoords(Point(u_min, v_max));
    
    std::cout << "  Feasible region (UV rotated): [(" 
              << u_min << ", " << v_min << ") - ("
              << u_max << ", " << v_max << ")]" << std::endl;
    std::cout << "  Feasible region (XY original, diamond shape):" << std::endl;
    std::cout << "    Corner1: (" << corner1_xy.get<0>() << ", " << corner1_xy.get<1>() << ")" << std::endl;
    std::cout << "    Corner2: (" << corner2_xy.get<0>() << ", " << corner2_xy.get<1>() << ")" << std::endl;
    std::cout << "    Corner3: (" << corner3_xy.get<0>() << ", " << corner3_xy.get<1>() << ")" << std::endl;
    std::cout << "    Corner4: (" << corner4_xy.get<0>() << ", " << corner4_xy.get<1>() << ")" << std::endl;
    
    // Verify R-Tree query before using it
    std::cout << "  ▸ Verifying R-Tree query for this cluster..." << std::endl;
  }

  std::vector<FlopClusterEntry> intersecting_entries = getIntersectedCluster(cluster, verbose);

  if (verbose) {
    std::cout << "  ▸ R-Tree query result: " << intersecting_entries.size() << " intersecting cluster(s)" << std::endl;
    
    // Check if this cluster finds itself in the R-Tree
    bool found_self = false;
    for (const auto& [box, idx] : intersecting_entries) {
      if (idx == i) {
        found_self = true;
        break;
      }
    }
    
    if (found_self) {
      std::cout << "    ✓ Self-intersection detected (normal - cluster finds itself)" << std::endl;
    } else {
      std::cout << "    ⚠ WARNING: Cluster does NOT find itself in R-Tree!" << std::endl;
      std::cout << "    This may indicate R-Tree indexing issue" << std::endl;
    }
    
    // Show distribution of intersecting cluster sizes
    if (!intersecting_entries.empty()) {
      std::map<int, int> size_dist;
      for (const auto& [box, idx] : intersecting_entries) {
        if (idx != i && flop_cluster_is_valid_[idx]) {
          int size = flop_clusters_[idx].flops_.size();
          size_dist[size]++;
        }
      }
      
      if (!size_dist.empty()) {
        std::cout << "    • Intersecting cluster sizes: ";
        bool first = true;
        for (const auto& [size, count] : size_dist) {
          if (!first) std::cout << ", ";
          std::cout << size << "-bit(" << count << ")";
          first = false;
        }
        std::cout << std::endl;
      }
    }
  }

  int valid_candidates = 0;
  int skipped_lower_id = 0;
  int skipped_invalid = 0;
  int skipped_incompatible = 0;
  int skipped_no_placement = 0;
  int edges_created = 0;

  for (const auto& entry : intersecting_entries) {
    const int j = entry.second;

    // Skip self-edges
    if (i == j) {
      continue;
    }

    // Skip if edge already exists between i and j
    bool edge_exists = false;
    for (const auto& edge : adj_list_[i]) {
      if ((edge.n1 == i && edge.n2 == j) || (edge.n1 == j && edge.n2 == i)) {
        edge_exists = true;
        break;
      }
    }
    
    if (edge_exists) {
      skipped_lower_id++;  // Reuse this counter for "already has edge"
      continue;
    }

    if (!flop_cluster_is_valid_[j] || flop_cluster_no_further_merge_[j]) {
      skipped_invalid++;
      continue;
    }

    const FlopCluster& adj_cluster = flop_clusters_[j];

    if (!checkCompatibility(cluster, adj_cluster, verbose && valid_candidates < 3)) {
      skipped_incompatible++;
      continue;
    }

    valid_candidates++;

    PlacementCandidate pc = calcPlacementCandidate(cluster, adj_cluster, verbose && edges_created < 3);

    if (!pc.is_valid_) {
      skipped_no_placement++;
      continue;
    }

    Edge new_edge(i, j, pc.merge_gain_, pc.placement_pos_);

    edge_pq_.insert(new_edge);
    adj_list_[i].insert(new_edge);
    adj_list_[j].insert(new_edge);
    edges_created++;

    if (verbose && edges_created <= 3) {
      std::cout << "    Edge #" << edges_created << " created with Cluster[" << j << "]" << std::endl;
      std::cout << "      HPWL gain: " << std::fixed << std::setprecision(2) << pc.merge_gain_ 
                << (pc.merge_gain_ < 0 ? " (improvement)" : " (degradation)") << std::endl;
      std::cout << "      Merge position: (" << pc.placement_pos_.x << ", " << pc.placement_pos_.y << ")" << std::endl;
      std::cout << "      Adjacent cluster size: " << adj_cluster.flops_.size() << " flop(s)" << std::endl;
      std::cout << "      Merged size would be: " << (cluster.flops_.size() + adj_cluster.flops_.size()) << " flop(s)" << std::endl;
    }
  }

  if (verbose) {
    std::cout << "  Edge creation summary:" << std::endl;
    std::cout << "    - Intersecting clusters: " << intersecting_entries.size() << std::endl;
    std::cout << "    - Skipped (edge already exists): " << skipped_lower_id << std::endl;
    std::cout << "    - Skipped (invalid/no-merge): " << skipped_invalid << std::endl;
    std::cout << "    - Compatible candidates: " << valid_candidates << std::endl;
    std::cout << "    - Skipped (incompatible masks): " << skipped_incompatible << std::endl;
    std::cout << "    - Skipped (no valid placement): " << skipped_no_placement << std::endl;
    std::cout << "    - Edges successfully created: " << edges_created << std::endl;
    if (edges_created > 3) {
      std::cout << "      (showing first 3 edges only)" << std::endl;
    }
  }
}

// -----------------------------------------------------------------------------
// [Level 2] Main logic for edge creation
// -----------------------------------------------------------------------------

std::vector<FlopClusterEntry>
AggloCluster::getIntersectedCluster(const FlopCluster& cluster, bool verbose) const
{
  std::vector<FlopClusterEntry> intersecting_clusters;

  if (boost::geometry::is_empty(cluster.feasible_region_)) {
    if (verbose) {
      std::cout << "      [getIntersectedCluster] Empty feasible region - returning 0 clusters" << std::endl;
    }
    return intersecting_clusters;
  }

  feasible_regions_.query(bgi::intersects(cluster.feasible_region_),
                          std::back_inserter(intersecting_clusters));

  if (verbose) {
    std::cout << "      [getIntersectedCluster] RTree query completed" << std::endl;
    std::cout << "        Intersecting clusters found: " << intersecting_clusters.size() << std::endl;
    
    if (intersecting_clusters.size() > 0 && intersecting_clusters.size() <= 5) {
      std::cout << "        Cluster IDs: ";
      for (size_t i = 0; i < intersecting_clusters.size(); ++i) {
        std::cout << intersecting_clusters[i].second;
        if (i < intersecting_clusters.size() - 1) {
          std::cout << ", ";
        }
      }
      std::cout << std::endl;
    }
  }

  return intersecting_clusters;
}

bool 
AggloCluster::checkCompatibility(const FlopCluster& c1,
                                 const FlopCluster& c2,
                                 bool verbose) const
{
  if (verbose) {
    std::cout << "      [checkCompatibility] Checking clusters " << c1.id_ << " and " << c2.id_ << std::endl;
  }

  // Check mask compatibility
  bool masks_match = (c1.master_mask_ == c2.master_mask_
                      && c1.inst_mask_ == c2.inst_mask_);
  
  if (verbose) {
    std::cout << "        Master mask match: " << (c1.master_mask_ == c2.master_mask_ ? "YES" : "NO") << std::endl;
    std::cout << "        Inst mask match: " << (c1.inst_mask_ == c2.inst_mask_ ? "YES" : "NO") << std::endl;
  }

  if (!masks_match) {
    if (verbose) {
      std::cout << "        Result: INCOMPATIBLE (mask mismatch)" << std::endl;
    }
    return false;
  }

  // Check if merged size is supported
  size_t total_flops = c1.flops_.size() + c2.flops_.size();
  
  if (verbose) {
    std::cout << "        Cluster sizes: " << c1.flops_.size() << " + " << c2.flops_.size() 
              << " = " << total_flops << " flop(s)" << std::endl;
  }

  auto it = compatible_masters_.find(c1.master_mask_);

  if (it != compatible_masters_.end()) {
    const auto& bits_to_masters = it->second;
    bool has_master = bits_to_masters.count(total_flops) > 0;
    
    if (verbose) {
      std::cout << "        Master available for " << total_flops << "-bit: " 
                << (has_master ? "YES" : "NO") << std::endl;
      if (!has_master) {
        std::cout << "        Available bit widths: ";
        for (const auto& [bits, masters] : bits_to_masters) {
          std::cout << bits << " ";
        }
        std::cout << std::endl;
      }
      std::cout << "        Result: " << (has_master ? "COMPATIBLE" : "INCOMPATIBLE (no master)") << std::endl;
    }
    
    return has_master;
  }

  if (verbose) {
    std::cout << "        Result: INCOMPATIBLE (mask not found in compatible_masters)" << std::endl;
  }

  return false;
}

PlacementCandidate
AggloCluster::calcPlacementCandidate(const FlopCluster& c1, const FlopCluster& c2, bool verbose)
{
  if (verbose) {
    std::cout << "      [calcPlacementCandidate] Computing placement for clusters " 
              << c1.id_ << " and " << c2.id_ << std::endl;
  }

  // Get masters for both clusters and the merged result
  const MasterMask& mask = c1.master_mask_;
  
  auto get_cluster_master = [&, this](const FlopCluster& c) -> odb::dbMaster* {
    int bits = c.flops_.size();
    if (bits == 1) {
      int flop_id = *c.flops_.begin();
      return flop_units_[flop_id].inst_->getMaster();
    } else {
      auto mask_it = representative_masters_.find(c.master_mask_);
      if (mask_it == representative_masters_.end()) {
        return nullptr;
      }
      auto bit_it = mask_it->second.find(bits);
      if (bit_it == mask_it->second.end()) {
        return nullptr;
      }
      return bit_it->second;
    }
  };
  
  odb::dbMaster* master1 = get_cluster_master(c1);
  odb::dbMaster* master2 = get_cluster_master(c2);
  
  // Find the master for merged cluster
  int new_bits = c1.flops_.size() + c2.flops_.size();
  odb::dbMaster* new_master = nullptr;
  auto mask_it = representative_masters_.find(mask);
  if (mask_it != representative_masters_.end()) {
    auto bit_it = mask_it->second.find(new_bits);
    if (bit_it != mask_it->second.end()) {
      new_master = bit_it->second;
    }
  }

  if (verbose) {
    std::cout << "        Master lookup:" << std::endl;
    std::cout << "          Cluster[" << c1.id_ << "] (" << c1.flops_.size() << "-bit): " 
              << (master1 ? master1->getName() : "NULL") << std::endl;
    std::cout << "          Cluster[" << c2.id_ << "] (" << c2.flops_.size() << "-bit): " 
              << (master2 ? master2->getName() : "NULL") << std::endl;
    std::cout << "          Merged (" << new_bits << "-bit): " 
              << (new_master ? new_master->getName() : "NULL") << std::endl;
  }

  // Return invalid if required masters are unavailable
  if (master1 == nullptr || master2 == nullptr || new_master == nullptr) {
    if (verbose) {
      std::cout << "        Result: INVALID (missing master)" << std::endl;
    }
    return PlacementCandidate(c1.curr_pt_, 0.0, false);
  }

  // Calculate HPWL-optimal bounding box in XY coordinates
  const Box hpwl_box_xy = calcMedianBox(c1, c2);

  if (verbose) {
    std::cout << "        HPWL-optimal box (XY): [(" 
              << hpwl_box_xy.min_corner().get<0>() << ", "
              << hpwl_box_xy.min_corner().get<1>() << ") - ("
              << hpwl_box_xy.max_corner().get<0>() << ", "
              << hpwl_box_xy.max_corner().get<1>() << ")]" << std::endl;
  }

  // Get timing-feasible region intersection in UV coordinates
  const Box fr_box_uv = getFeasibleRegionIntersection(c1, c2);
  if (boost::geometry::is_empty(fr_box_uv)) {
    if (verbose) {
      std::cout << "        Result: INVALID (empty feasible region intersection)" << std::endl;
    }
    return PlacementCandidate(c1.curr_pt_, 0.0, false);
  }

  if (verbose) {
    Point fr_min_uv = fr_box_uv.min_corner();
    Point fr_max_uv = fr_box_uv.max_corner();
    int u_min = fr_min_uv.get<0>();
    int v_min = fr_min_uv.get<1>();
    int u_max = fr_max_uv.get<0>();
    int v_max = fr_max_uv.get<1>();
    
    // Four corners of the UV box
    Point corner1_xy = inverseTransformCoords(Point(u_min, v_min));
    Point corner2_xy = inverseTransformCoords(Point(u_max, v_min));
    Point corner3_xy = inverseTransformCoords(Point(u_max, v_max));
    Point corner4_xy = inverseTransformCoords(Point(u_min, v_max));
    
    std::cout << "        Feasible region (UV rotated): [(" 
              << u_min << ", " << v_min << ") - ("
              << u_max << ", " << v_max << ")]" << std::endl;
    std::cout << "        Feasible region (XY original, diamond shape):" << std::endl;
    std::cout << "          Corner1: (" << corner1_xy.get<0>() << ", " << corner1_xy.get<1>() << ")" << std::endl;
    std::cout << "          Corner2: (" << corner2_xy.get<0>() << ", " << corner2_xy.get<1>() << ")" << std::endl;
    std::cout << "          Corner3: (" << corner3_xy.get<0>() << ", " << corner3_xy.get<1>() << ")" << std::endl;
    std::cout << "          Corner4: (" << corner4_xy.get<0>() << ", " << corner4_xy.get<1>() << ")" << std::endl;
  }

  // Project HPWL-optimal point onto feasible region
  Point P_proj_xy = project(hpwl_box_xy, fr_box_uv);

  if (verbose) {
    std::cout << "        Projected optimal point (XY): (" 
              << P_proj_xy.get<0>() << ", " << P_proj_xy.get<1>() << ")" << std::endl;
  }

  // Generate uniform samples within feasible region
  std::vector<Point> xy_candidates = generateUniformSamples(fr_box_uv, 16);
  xy_candidates.push_back(P_proj_xy);

  if (verbose) {
    std::cout << "        Generated " << xy_candidates.size() << " candidate positions" << std::endl;
  }

  // Sort candidates by Manhattan distance to projected point
  auto manhattan_dist = [](const Point& a, const Point& b) {
    int64_t dx = std::abs(static_cast<int64_t>(a.get<0>()) - static_cast<int64_t>(b.get<0>()));
    int64_t dy = std::abs(static_cast<int64_t>(a.get<1>()) - static_cast<int64_t>(b.get<1>()));
    return dx + dy;
  };

  std::sort(xy_candidates.begin(),
            xy_candidates.end(),
            [&](const Point& a, const Point& b) {
              return manhattan_dist(a, P_proj_xy) < manhattan_dist(b, P_proj_xy);
            });

  if (verbose) {
    std::cout << "        Evaluating candidates (sorted by Manhattan distance):" << std::endl;
  }

  // Evaluate candidates in order of proximity to optimal point
  int candidate_idx = 0;
  for (const auto& xy_cand_int : xy_candidates) { 
    
    if (verbose && candidate_idx < 3) {
      std::cout << "          Candidate #" << candidate_idx << ": (" 
                << xy_cand_int.get<0>() << ", " << xy_cand_int.get<1>() << ")" << std::endl;
      int64_t dist = manhattan_dist(xy_cand_int, P_proj_xy);
      std::cout << "            Manhattan distance to optimal: " << dist << std::endl;
    }

    // Check if placement satisfies density constraint
    if (checkPlacementDensityConstraint(c1, master1, c2, master2, xy_cand_int, new_master)) {
      
      // Convert to float coordinates
      FloatPoint xy_cand_float(static_cast<float>(xy_cand_int.get<0>()),
                               static_cast<float>(xy_cand_int.get<1>()));

      if (verbose) {
        std::cout << "            Density check: PASSED" << std::endl;
      }
      
      // Calculate HPWL gain: negative value indicates improvement
      double gain = calcMergeHPWLGain(c1, c2, xy_cand_float, verbose);

      if (verbose) {
        std::cout << "        Result: VALID (selected candidate #" << candidate_idx << ")" << std::endl;
      }

      // Return first valid candidate
      return PlacementCandidate(xy_cand_float, gain, true);
    } else {
      if (verbose && candidate_idx < 3) {
        std::cout << "            Density check: FAILED (overflow)" << std::endl;
      }
    }

    candidate_idx++;
  }

  // All candidates failed density check
  if (verbose) {
    std::cout << "        Result: INVALID (all " << xy_candidates.size() 
              << " candidates failed density check)" << std::endl;
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

  // Remove old clusters and add merged cluster
  virtual_bin_grid_.removeInst(master1, c1_pos);
  virtual_bin_grid_.removeInst(master2, c2_pos);
  virtual_bin_grid_.addInst(new_master, new_pos);

  bool has_overflow = virtual_bin_grid_.checkOverflow();

  // Rollback if overflow detected
  if (has_overflow) {
    virtual_bin_grid_.removeInst(new_master, new_pos);
    virtual_bin_grid_.addInst(master2, c2_pos);
    virtual_bin_grid_.addInst(master1, c1_pos);
    return false;
  }
  
  return true; 
}

double 
AggloCluster::calcMergeHPWLGain(const FlopCluster& c1, 
                                const FlopCluster& c2, 
                                const FloatPoint& merge_position,
                                bool verbose) const
{
  if (verbose) {
    std::cout << "      [calcMergeHPWLGain] Calculating HPWL gain for merge" << std::endl;
    std::cout << "        c1 current position: (" 
              << static_cast<int>(c1.curr_pt_.x) << ", "
              << static_cast<int>(c1.curr_pt_.y) << "), flops: " << c1.flops_.size() << std::endl;
    std::cout << "        c2 current position: (" 
              << static_cast<int>(c2.curr_pt_.x) << ", "
              << static_cast<int>(c2.curr_pt_.y) << "), flops: " << c2.flops_.size() << std::endl;
    std::cout << "        Merge position: (" 
              << static_cast<int>(merge_position.x) << ", "
              << static_cast<int>(merge_position.y) << ")" << std::endl;
  }
  
  // Collect all instances in both clusters
  std::set<odb::dbInst*> c1_insts, c2_insts;
  for (int flop_idx : c1.flops_) {
    c1_insts.insert(flop_units_[flop_idx].inst_);
  }
  for (int flop_idx : c2.flops_) {
    c2_insts.insert(flop_units_[flop_idx].inst_);
  }
  
  // Collect all unique nets connected to either cluster
  std::set<odb::dbNet*> all_nets;
  for (int flop_idx : c1.flops_) {
    for (odb::dbITerm* iterm : flop_units_[flop_idx].inst_->getITerms()) {
      odb::dbNet* net = iterm->getNet();
      if (net && (net->getSigType() == odb::dbSigType::SIGNAL || 
                  net->getSigType() == odb::dbSigType::CLOCK)) {
        all_nets.insert(net);
      }
    }
  }
  for (int flop_idx : c2.flops_) {
    for (odb::dbITerm* iterm : flop_units_[flop_idx].inst_->getITerms()) {
      odb::dbNet* net = iterm->getNet();
      if (net && (net->getSigType() == odb::dbSigType::SIGNAL || 
                  net->getSigType() == odb::dbSigType::CLOCK)) {
        all_nets.insert(net);
      }
    }
  }
  
  if (verbose) {
    std::cout << "        Total unique nets: " << all_nets.size() << std::endl;
  }
  
  double original_hpwl = 0.0;
  double merged_hpwl = 0.0;
  int net_count = 0;
  
  int merge_x = static_cast<int>(merge_position.x);
  int merge_y = static_cast<int>(merge_position.y);
  int c1_x = static_cast<int>(c1.curr_pt_.x);
  int c1_y = static_cast<int>(c1.curr_pt_.y);
  int c2_x = static_cast<int>(c2.curr_pt_.x);
  int c2_y = static_cast<int>(c2.curr_pt_.y);
  
  // Process each net
  for (odb::dbNet* net : all_nets) {
    // Check connections to each cluster
    bool has_c1_connection = false;
    bool has_c2_connection = false;
    
    // Collect external pin coordinates (not in c1 or c2)
    std::vector<std::pair<int, int>> external_pins;
    
    for (odb::dbITerm* net_iterm : net->getITerms()) {
      odb::dbInst* inst = net_iterm->getInst();
      
      if (c1_insts.find(inst) != c1_insts.end()) {
        has_c1_connection = true;
      } else if (c2_insts.find(inst) != c2_insts.end()) {
        has_c2_connection = true;
      } else {
        FloatPoint pin_coord = getPinCoordinate(net_iterm);
        external_pins.emplace_back(static_cast<int>(pin_coord.x), 
                                   static_cast<int>(pin_coord.y));
      }
    }
    
    // Include top-level ports (BTerms)
    for (odb::dbBTerm* net_bterm : net->getBTerms()) {
      FloatPoint pin_coord = getPinCoordinate(net_bterm);
      external_pins.emplace_back(static_cast<int>(pin_coord.x), 
                                 static_cast<int>(pin_coord.y));
    }
    
    // Calculate original HPWL (before merge)
    double net_original_hpwl = 0.0;
    int min_x = std::numeric_limits<int>::max();
    int max_x = std::numeric_limits<int>::min();
    int min_y = std::numeric_limits<int>::max();
    int max_y = std::numeric_limits<int>::min();
    
    for (const auto& [x, y] : external_pins) {
      min_x = std::min(min_x, x);
      max_x = std::max(max_x, x);
      min_y = std::min(min_y, y);
      max_y = std::max(max_y, y);
    }
    
    if (has_c1_connection) {
      min_x = std::min(min_x, c1_x);
      max_x = std::max(max_x, c1_x);
      min_y = std::min(min_y, c1_y);
      max_y = std::max(max_y, c1_y);
    }
    
    if (has_c2_connection) {
      min_x = std::min(min_x, c2_x);
      max_x = std::max(max_x, c2_x);
      min_y = std::min(min_y, c2_y);
      max_y = std::max(max_y, c2_y);
    }
    
    if (max_x >= min_x && max_y >= min_y) {
      net_original_hpwl = (max_x - min_x) + (max_y - min_y);
      original_hpwl += net_original_hpwl;
    }
    
    // Calculate merged HPWL (after merge)
    double net_merged_hpwl = 0.0;
    min_x = std::numeric_limits<int>::max();
    max_x = std::numeric_limits<int>::min();
    min_y = std::numeric_limits<int>::max();
    max_y = std::numeric_limits<int>::min();
    
    for (const auto& [x, y] : external_pins) {
      min_x = std::min(min_x, x);
      max_x = std::max(max_x, x);
      min_y = std::min(min_y, y);
      max_y = std::max(max_y, y);
    }
    
    // Both clusters now at merge position
    if (has_c1_connection || has_c2_connection) {
      min_x = std::min(min_x, merge_x);
      max_x = std::max(max_x, merge_x);
      min_y = std::min(min_y, merge_y);
      max_y = std::max(max_y, merge_y);
    }
    
    if (max_x >= min_x && max_y >= min_y) {
      net_merged_hpwl = (max_x - min_x) + (max_y - min_y);
      merged_hpwl += net_merged_hpwl;
    }
    
    net_count++;
    if (verbose && net_count <= 3) {
      std::cout << "        Net #" << net_count << " '" << net->getName() << "':" << std::endl;
      std::cout << "          c1 connected: " << (has_c1_connection ? "YES" : "NO") << std::endl;
      std::cout << "          c2 connected: " << (has_c2_connection ? "YES" : "NO") << std::endl;
      std::cout << "          External pins: " << external_pins.size() << std::endl;
      std::cout << "          Original HPWL: " << net_original_hpwl << " DBU" << std::endl;
      std::cout << "          Merged HPWL: " << net_merged_hpwl << " DBU" << std::endl;
      std::cout << "          Net gain: " << (net_merged_hpwl - net_original_hpwl) << " DBU" << std::endl;
    }
  }
  
  double gain = merged_hpwl - original_hpwl;
  
  if (verbose) {
    std::cout << "        Nets processed: " << net_count << std::endl;
    std::cout << "        Original HPWL (c1, c2 separate): " << original_hpwl << " DBU" << std::endl;
    std::cout << "        Merged HPWL (at merge position): " << merged_hpwl << " DBU" << std::endl;
    std::cout << "        HPWL gain: " << gain << (gain < 0 ? " (improvement)" : " (degradation)") << std::endl;
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
      std::cerr << "[FATAL] Failed to get coordinates for ITerm: " 
                << iterm->getName() << std::endl;
      std::exit(1);
    }
    return FloatPoint(static_cast<float>(x), static_cast<float>(y));
  } 
  
  odb::dbBTerm* bterm = std::get<odb::dbBTerm*>(pin_variant);
  odb::Rect bbox = bterm->getBBox();
  
  if (bbox.isInverted()) {
    std::cerr << "[FATAL] Invalid BBox for BTerm: " 
              << bterm->getName() << std::endl;
    std::exit(1);
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
    std::cout << "        C1 size: " << c1.flops_.size() << " flop(s), position: (" 
              << c1.curr_pt_.x << ", " << c1.curr_pt_.y << ")" << std::endl;
    std::cout << "        C2 size: " << c2.flops_.size() << " flop(s), position: (" 
              << c2.curr_pt_.x << ", " << c2.curr_pt_.y << ")" << std::endl;
    std::cout << "        New position: (" << edge.pos.x << ", " << edge.pos.y << ")" << std::endl;
  }

  // Step 1: Mark old clusters as invalid
  flop_cluster_is_valid_[edge.n1] = false;
  flop_cluster_is_valid_[edge.n2] = false;

  if (verbose) {
    std::cout << "        [Step 1] Marked old clusters as invalid" << std::endl;
  }

  // Step 2: Remove edges associated with old clusters
  removeEdges(c1, verbose);
  removeEdges(c2, verbose);

  if (verbose) {
    std::cout << "        [Step 2] Removed edges from old clusters" << std::endl;
  }

  // Step 3: Remove feasible regions of old clusters from R-Tree
  feasible_regions_.remove(std::make_pair(c1.feasible_region_, edge.n1));
  feasible_regions_.remove(std::make_pair(c2.feasible_region_, edge.n2));

  if (verbose) {
    std::cout << "        [Step 3] Removed feasible regions from R-Tree" << std::endl;
  }

  // Step 4: Create new merged cluster
  const int new_cluster_idx = flop_clusters_.size();
  flop_clusters_.emplace_back(new_cluster_idx, c1, c2, edge);
  flop_cluster_is_valid_.push_back(true);
  flop_cluster_no_further_merge_.push_back(false);
  
  FlopCluster& new_cluster = flop_clusters_.back();

  if (verbose) {
    std::cout << "        [Step 4] Created new cluster #" << new_cluster_idx << std::endl;
    std::cout << "          New cluster size: " << new_cluster.flops_.size() << " flop(s)" << std::endl;
    
    const Box& fr = new_cluster.feasible_region_;
    const bool fr_empty = boost::geometry::is_empty(fr);
    if (fr_empty) {
      std::cout << "          WARNING: Feasible region is EMPTY!" << std::endl;
    } else {
      std::cout << "          Feasible region: (" << fr.min_corner().get<0>() << ", " 
                << fr.min_corner().get<1>() << ") -> (" << fr.max_corner().get<0>() << ", " 
                << fr.max_corner().get<1>() << ")" << std::endl;
    }
  }

  // Step 5: Update all flop units that belong to this new cluster
  for (int flop_idx : new_cluster.flops_) {
    flop_units_[flop_idx].cluster_idx_ = new_cluster_idx;
    flop_units_[flop_idx].curr_pt_ = new_cluster.curr_pt_;
  }

  if (verbose) {
    std::cout << "        [Step 5] Updated " << new_cluster.flops_.size() << " flop unit(s)" << std::endl;
  }

  // Step 6: Insert new cluster's feasible region into R-Tree
  if (!boost::geometry::is_empty(new_cluster.feasible_region_)) {
    feasible_regions_.insert(std::make_pair(new_cluster.feasible_region_, new_cluster_idx));
    
    if (verbose) {
      std::cout << "        [Step 6] Inserted feasible region into R-Tree" << std::endl;
      
      // Verify insertion by querying the R-Tree
      std::vector<FlopClusterEntry> query_result;
      feasible_regions_.query(
          boost::geometry::index::intersects(new_cluster.feasible_region_),
          std::back_inserter(query_result)
      );
      
      // Check if our newly inserted cluster appears in the query
      bool found = false;
      int self_count = 0;
      for (const auto& [box, idx] : query_result) {
        if (idx == new_cluster_idx) {
          found = true;
          self_count++;
        }
      }
      
      if (found) {
        std::cout << "          ✓ Verification: Found in R-Tree query" << std::endl;
        if (self_count > 1) {
          std::cout << "          ⚠ WARNING: Found " << self_count 
                    << " times (should be 1) - duplicate insertion!" << std::endl;
        }
      } else {
        std::cout << "          ✗ ERROR: NOT found in R-Tree after insertion!" << std::endl;
        std::cout << "          This is a critical bug - cluster won't be discoverable" << std::endl;
      }
      
      // Show overlapping clusters
      int overlap_count = query_result.size() - self_count;
      std::cout << "          • Overlapping clusters: " << overlap_count << std::endl;
      
      if (overlap_count > 0 && overlap_count <= 5) {
        std::cout << "          • Overlapping cluster IDs: ";
        bool first = true;
        for (const auto& [box, idx] : query_result) {
          if (idx != new_cluster_idx) {
            if (!first) std::cout << ", ";
            std::cout << idx;
            first = false;
          }
        }
        std::cout << std::endl;
      } else if (overlap_count > 5) {
        std::cout << "          • (showing first 5 overlapping clusters)" << std::endl;
        std::cout << "          • Overlapping cluster IDs: ";
        int shown = 0;
        for (const auto& [box, idx] : query_result) {
          if (idx != new_cluster_idx && shown < 5) {
            if (shown > 0) std::cout << ", ";
            std::cout << idx;
            shown++;
          }
        }
        std::cout << " ... (+" << (overlap_count - 5) << " more)" << std::endl;
      }
    }
  } else if (verbose) {
    std::cout << "        [Step 6] Skipped R-Tree insertion (empty feasible region)" << std::endl;
    std::cout << "          ⚠ WARNING: Empty feasible region means this cluster" << std::endl;
    std::cout << "          won't be discoverable for future edge creation!" << std::endl;
  }

  return new_cluster_idx;
}

bool
AggloCluster::isFurtherMergeable(const FlopCluster& cluster, bool verbose) const
{
  // Check if compatible masters exist for this cluster's mask
  const auto it = compatible_masters_.find(cluster.master_mask_);
  if (it == compatible_masters_.end()) {
    if (verbose) {
      std::cout << "    ✗ Mask not found in compatible_masters_" << std::endl;
      std::cout << "      (This should not happen - possible bug)" << std::endl;
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

  // Get maximum bit-width available for this mask
  const int max_bits = bits_to_masters.rbegin()->first;
  const int current_size = static_cast<int>(cluster.flops_.size());

  // Can merge further if current size is less than max available bits
  const bool can_merge = current_size < max_bits;
  
  if (verbose) {
    std::cout << "    ▸ isFurtherMergeable: " << (can_merge ? "✓ YES" : "✗ NO") << std::endl;
    std::cout << "      Current size: " << current_size << " bit(s)" << std::endl;
    std::cout << "      Maximum available: " << max_bits << " bit(s)" << std::endl;
    
    if (!can_merge && current_size == max_bits) {
      std::cout << "      Reason: Already at maximum bit-width" << std::endl;
    } else if (!can_merge) {
      std::cout << "      Reason: Current size (" << current_size 
                << ") >= max (" << max_bits << ")" << std::endl;
    }
    
    // Show available bit-widths
    std::cout << "      Available bit-widths: ";
    bool first = true;
    for (const auto& [bits, masters] : bits_to_masters) {
      if (!first) std::cout << ", ";
      std::cout << bits;
      first = false;
    }
    std::cout << std::endl;
  }
  
  return can_merge;
}  

std::set<int>
AggloCluster::distributeSlack(const FlopCluster& cluster, bool verbose)
{
  if (verbose) {
    std::cout << "        [distributeSlack] Processing cluster #" << cluster.id_ 
              << " with " << cluster.flops_.size() << " flop(s)" << std::endl;
  }

  std::set<int> affected_flop_units;
  int total_paths_processed = 0;
  int internal_paths_skipped = 0;
  int slack_redistributed_count = 0;
  int detailed_outputs = 0;
  constexpr int MAX_DETAILED_OUTPUTS = 5;  // Show first 5 redistributions in detail

  // Process each flop in the cluster
  for (int flop_idx : cluster.flops_) {
    const FlopUnit& flop = flop_units_[flop_idx];
    
    if (verbose && detailed_outputs < MAX_DETAILED_OUTPUTS) {
      std::cout << "          ▸ Processing flop #" << flop_idx 
                << " (inst: " << flop.inst_->getName() << ")" << std::endl;
      std::cout << "            Original position: (" << flop.orig_pt_.x << ", " << flop.orig_pt_.y << ")" << std::endl;
      std::cout << "            Current position: (" << flop.curr_pt_.x << ", " << flop.curr_pt_.y << ")" << std::endl;
      
      const float dx = std::abs(flop.curr_pt_.x - flop.orig_pt_.x);
      const float dy = std::abs(flop.curr_pt_.y - flop.orig_pt_.y);
      const float manhattan_distance = dx + dy;
      std::cout << "            Manhattan distance: " << manhattan_distance << " DBU" << std::endl;
    }

    // Calculate used slacks for each pin of this flop
    const auto used_slacks = calcUsedSlacks(flop, verbose && detailed_outputs < MAX_DETAILED_OUTPUTS);

    // Distribute slack to connected flops outside this cluster
    for (const auto& [pin, used_slack_list] : used_slacks) {
      for (const auto& [path_idx, used_slack] : used_slack_list) {
        total_paths_processed++;
        const TimingPath& path = timing_paths_[path_idx];

        // Identify the flop on the other end of this timing path
        const int other_flop_idx = (path.start_flop_idx_ == flop_idx) 
                                     ? path.end_flop_idx_ 
                                     : path.start_flop_idx_;

        // Skip if the other flop is also in the same cluster (internal path)
        if (flop_units_[other_flop_idx].cluster_idx_ == cluster.id_) {
          internal_paths_skipped++;
          continue;
        }

        FlopUnit& other_flop = flop_units_[other_flop_idx];
        affected_flop_units.insert(other_flop_idx);

        // Get the corresponding pin on the other flop
        odb::dbITerm* other_pin = (path.start_flop_idx_ == flop_idx) 
                                    ? path.end_pin_ 
                                    : path.start_pin_;
        
        if (!other_pin || other_flop.pin_budgets_.find(other_pin) == other_flop.pin_budgets_.end()) {
          continue;
        }

        // Find the original slack budget for this path on current flop
        sta::Slack original_budget = 0.0;
        const auto& current_budgets = flop.pin_budgets_.at(pin);
        for (const auto& [budget_path_idx, budget_slack] : current_budgets) {
          if (budget_path_idx == path_idx) {
            original_budget = budget_slack;
            break;
          }
        }

        // Update the slack budget on the other flop
        auto& other_budgets = other_flop.pin_budgets_.at(other_pin);
        sta::Slack old_other_budget = 0.0;
        sta::Slack new_other_budget = 0.0;
        
        for (auto& [budget_path_idx, budget_slack] : other_budgets) {
          if (budget_path_idx == path_idx) {
            old_other_budget = budget_slack;
            
            // Redistribute the unused slack to the other flop
            const sta::Slack slack_delta = original_budget - used_slack;
            budget_slack += slack_delta;
            new_other_budget = budget_slack;
            
            if (verbose && detailed_outputs < MAX_DETAILED_OUTPUTS) {
              std::cout << "            ━━ Slack Redistribution ━━" << std::endl;
              std::cout << "               Path #" << path_idx << ": Flop[" << flop_idx << "] → Flop[" << other_flop_idx << "]" << std::endl;
              std::cout << "               Pin: " << pin->getMTerm()->getName() << " → " << other_pin->getMTerm()->getName() << std::endl;
              std::cout << "               Original slack budget (current flop): " << std::fixed << std::setprecision(12) << original_budget << std::endl;
              std::cout << "               Used slack (after movement): " << std::fixed << std::setprecision(12) << used_slack << std::endl;
              std::cout << "               Unused slack (to redistribute): " << std::fixed << std::setprecision(12) << slack_delta << std::endl;
              std::cout << "               Other flop slack: " << std::fixed << std::setprecision(12) << old_other_budget << " → " << new_other_budget << std::endl;
              
              if (slack_delta > 0) {
                std::cout << "               ✓ Positive redistribution (slack increased)" << std::endl;
              } else if (slack_delta < 0) {
                std::cout << "               ⚠ Negative redistribution (slack decreased)" << std::endl;
              } else {
                std::cout << "               ○ No change (exact usage)" << std::endl;
              }
              detailed_outputs++;
            }
            
            slack_redistributed_count++;
            break;
          }
        }

        // Re-sort budgets for the updated pin (ascending order of slack)
        std::sort(other_budgets.begin(), other_budgets.end(),
                  [](const auto& a, const auto& b) { 
                    return a.second < b.second; 
                  });
      }
    }
  }
  
  if (verbose) {
    std::cout << "          ━━━━━━━━━━━━━━━━━━━━━━━━━━━" << std::endl;
    std::cout << "          Total paths processed: " << total_paths_processed << std::endl;
    std::cout << "          Internal paths skipped: " << internal_paths_skipped << std::endl;
    std::cout << "          Slack redistributed: " << slack_redistributed_count << " time(s)" << std::endl;
    std::cout << "          Affected flop units: " << affected_flop_units.size() << std::endl;
    if (detailed_outputs >= MAX_DETAILED_OUTPUTS) {
      std::cout << "          (Showing first " << MAX_DETAILED_OUTPUTS << " detailed redistributions)" << std::endl;
    }
  }

  return affected_flop_units;
}

void 
AggloCluster::removeEdges(const FlopCluster& cluster, bool verbose)
{
  const int cluster_idx = cluster.id_;

  // Check if this cluster has any edges in adjacency list
  if (adj_list_.find(cluster_idx) == adj_list_.end()) {
    if (verbose) {
      std::cout << "          [removeEdges] Cluster #" << cluster_idx 
                << ": No edges to remove" << std::endl;
    }
    return;
  }

  const auto& edges_to_remove = adj_list_.at(cluster_idx);
  const int edge_count = edges_to_remove.size();
  
  if (verbose) {
    std::cout << "          [removeEdges] Cluster #" << cluster_idx 
              << ": Removing " << edge_count << " edge(s)" << std::endl;
  }

  // Remove each edge from priority queue and neighbor's adjacency list
  for (const Edge& edge : edges_to_remove) {
    // Remove from global edge priority queue
    edge_pq_.erase(edge);

    // Remove from neighbor's adjacency list
    const int neighbor_idx = (edge.n1 == cluster_idx) ? edge.n2 : edge.n1;
    if (adj_list_.count(neighbor_idx)) {
      adj_list_.at(neighbor_idx).erase(edge);
    }
  }

  // Remove this cluster's entry from adjacency list
  adj_list_.erase(cluster_idx);
}


void 
AggloCluster::updateFeasibleRegion(FlopCluster& cluster, bool verbose)
{
  if (verbose) {
    std::cout << "          [updateFeasibleRegion] Cluster #" << cluster.id_ 
              << " with " << cluster.flops_.size() << " flop(s)" << std::endl;
  }

  const Box old_feasible_region = cluster.feasible_region_;
  const bool old_was_empty = boost::geometry::is_empty(old_feasible_region);

  // Calculate new feasible region as intersection of all member flops
  if (cluster.flops_.empty()) {
    cluster.feasible_region_ = Box();
    
    if (verbose) {
      std::cout << "            WARNING: Cluster has no flops!" << std::endl;
    }
  } else {
    // Initialize with first flop's region
    auto it = cluster.flops_.begin();
    cluster.feasible_region_ = flop_units_[*it].feasible_region_;
    
    if (verbose) {
      std::cout << "            Starting with flop #" << *it << "'s feasible region" << std::endl;
    }

    int intersection_count = 0;
    
    // Intersect with remaining flops' regions
    for (++it; it != cluster.flops_.end(); ++it) {
      const int flop_idx = *it;
      const Box& flop_feasible_region = flop_units_[flop_idx].feasible_region_;
      
      Box temp_result;
      boost::geometry::intersection(cluster.feasible_region_, 
                                    flop_feasible_region, 
                                    temp_result);
      
      intersection_count++;
      
      cluster.feasible_region_ = temp_result;
      
      // Early exit if intersection becomes empty
      if (boost::geometry::is_empty(cluster.feasible_region_)) {
        if (verbose) {
          std::cout << "            WARNING: Feasible region became EMPTY after " 
                    << intersection_count << " intersection(s)" << std::endl;
        }
        break;
      }
    }
    
    if (verbose && !boost::geometry::is_empty(cluster.feasible_region_)) {
      const Box& fr = cluster.feasible_region_;
      std::cout << "            New feasible region: (" << fr.min_corner().get<0>() << ", " 
                << fr.min_corner().get<1>() << ") -> (" << fr.max_corner().get<0>() << ", " 
                << fr.max_corner().get<1>() << ")" << std::endl;
      std::cout << "            Performed " << intersection_count << " intersection(s)" << std::endl;
    }
  }

  // Update R-Tree: remove old entry, insert new entry
  if (!old_was_empty) {
    feasible_regions_.remove(std::make_pair(old_feasible_region, cluster.id_));
    
    if (verbose) {
      std::cout << "            Removed old feasible region from R-Tree" << std::endl;
    }
  }
  
  if (!boost::geometry::is_empty(cluster.feasible_region_)) {
    feasible_regions_.insert(std::make_pair(cluster.feasible_region_, cluster.id_));
    
    if (verbose) {
      std::cout << "            Inserted new feasible region into R-Tree" << std::endl;
    }
  } else if (verbose) {
    std::cout << "            Skipped R-Tree insertion (empty feasible region)" << std::endl;
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
    std::cout << "\n[calcUsedSlacks] Analyzing flop: " << flop.inst_->getName() << std::endl;
    std::cout << "  Original position: (" << flop.orig_pt_.x << ", " << flop.orig_pt_.y << ")" << std::endl;
    std::cout << "  Current position:  (" << flop.curr_pt_.x << ", " << flop.curr_pt_.y << ")" << std::endl;
    
    const float dx = flop.curr_pt_.x - flop.orig_pt_.x;
    const float dy = flop.curr_pt_.y - flop.orig_pt_.y;
    const float manhattan_movement = std::abs(dx) + std::abs(dy);
    
    std::cout << "  Movement (Manhattan): " << manhattan_movement << " DBU" << std::endl;
  }
  
  // Get wire parameters
  const auto est = resizer_->getEstimateParasitics();
  const double unit_c = est->wireSignalCapacitance(corner_);
  const double unit_r = est->wireSignalResistance(corner_);
  
  // Classify pins
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
      std::cout << "  [WARNING] No clock pin found - returning empty result" << std::endl;
    }
    return used_slacks;
  }
  
  // Process output pins (Q/QN)
  std::vector<odb::dbITerm*> output_pins;
  output_pins.reserve(q_pins.size() + qn_pins.size());
  output_pins.insert(output_pins.end(), q_pins.begin(), q_pins.end());
  output_pins.insert(output_pins.end(), qn_pins.begin(), qn_pins.end());
  
  if (verbose && !output_pins.empty()) {
    std::cout << "\n  Processing " << output_pins.size() << " output pin(s)..." << std::endl;
  }
  
  for (odb::dbITerm* out_pin : output_pins) {
    calcUsedSlacksFanOut(flop, out_pin, clk_pin->getMTerm(), 
                         est, unit_r, unit_c, used_slacks, verbose);
  }
  
  // Process input pins (D)
  if (verbose && !d_pins.empty()) {
    std::cout << "\n  Processing " << d_pins.size() << " input pin(s)..." << std::endl;
  }
  
  for (odb::dbITerm* d_pin : d_pins) {
    calcUsedSlacksFanIn(flop, d_pin, est, unit_r, unit_c, used_slacks, verbose);
  }
  
  if (verbose) {
    std::cout << "[calcUsedSlacks] Completed.\n" << std::endl;
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
    std::cout << "\n  [FanOut Pin: " << out_pin->getMTerm()->getName() << "]" << std::endl;
  }
  
  // Check if this pin has timing constraints
  const auto budget_it = flop.pin_budgets_.find(out_pin);
  if (budget_it == flop.pin_budgets_.end() || budget_it->second.empty()) {
    if (verbose) {
      std::cout << "    No timing constraints - skipping" << std::endl;
    }
    return;
  }
  
  if (verbose) {
    std::cout << "    Total timing paths: " << budget_it->second.size() << std::endl;
  }
  
  // Validate net connection
  odb::dbNet* fanout_net = out_pin->getNet();
  if (!fanout_net) {
    if (verbose) {
      std::cout << "    Pin not connected to net - skipping" << std::endl;
    }
    return;
  }
  
  // Extract cell delay characterization
  const auto cap_delay = extractCapacitanceDelayPoints(
      flop.inst_, 
      clk_pin_lib->getName(), 
      out_pin->getMTerm()->getName(), 
      0);
  
  if (cap_delay.empty() || cap_delay.size() < 2) {
    if (verbose) {
      std::cout << "    Insufficient cap-delay data - skipping" << std::endl;
    }
    return;
  }
  
  std::vector<float> coeffs;
  coeffs.reserve(cap_delay.size() - 1);
  for (size_t i = 1; i < cap_delay.size(); ++i) {
    const float dy = cap_delay[i].second - cap_delay[i-1].second;
    const float dx = cap_delay[i].first - cap_delay[i-1].first;
    coeffs.push_back(dy / dx);
  }
  
  // Build Steiner tree and compute parasitics
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
  const double l1 = dbuToMeters(driver_steiner_tree->distance(driver_point, top_steiner_point));
  
  const float wire_cap = l1 * unit_c;
  const float wo_fst_stt_cap = total_net_capacitance - wire_cap;
  
  // Calculate actual MANHATTAN distance to Steiner point after movement
  const int steiner_x = steiner_location.getX();
  const int steiner_y = steiner_location.getY();
  const float dx_new = std::abs(flop.curr_pt_.x - steiner_x);
  const float dy_new = std::abs(flop.curr_pt_.y - steiner_y);
  const float actual_dist_dbu = dx_new + dy_new;
  const float actual_dist = dbuToMeters(actual_dist_dbu);
  
  if (verbose) {
    std::cout << "    Steiner point: (" << steiner_x << ", " << steiner_y << ")" << std::endl;
    std::cout << "    Original distance (l1): " << l1 << " m" << std::endl;
    std::cout << "    Actual distance: " << actual_dist << " m" << std::endl;
    std::cout << "    Distance change: " << (actual_dist - l1) << " m" << std::endl;
  }
  
  // Calculate used slack for each timing path
  std::vector<std::pair<int, sta::Slack>> pin_used_slacks;
  pin_used_slacks.reserve(budget_it->second.size());
  
  for (const auto& [path_idx, slack_budget] : budget_it->second) {
    // Find appropriate coefficient segment based on new capacitance
    float coeff = coeffs.empty() ? 0.0f : coeffs.back();
    const float new_cap = actual_dist * unit_c + wo_fst_stt_cap;
    
    for (size_t seg_idx = 0; seg_idx < coeffs.size(); ++seg_idx) {
      if (new_cap >= cap_delay[seg_idx].first && 
          new_cap <= cap_delay[seg_idx + 1].first) {
        coeff = coeffs[seg_idx];
        break;
      }
    }
    
    // Calculate delay increase components
    // 1. RC delay (quadratic term)
    const float rc_delay_increase = 
        (actual_dist * actual_dist - l1 * l1) * unit_r * unit_c;
    
    // 2. Wire-load interaction delay (linear term)
    const float wire_delay_increase = 
        (actual_dist - l1) * wo_fst_stt_cap * unit_r;
    
    // 3. Cell delay increase due to capacitance change
    const float cell_delay_increase = 
        coeff * (actual_dist - l1) * unit_c;
    
    // Total delay increase = used slack
    const float total_delay_increase = 
        rc_delay_increase + wire_delay_increase + cell_delay_increase;
    
    const sta::Slack used_slack = total_delay_increase;
    
    pin_used_slacks.emplace_back(path_idx, used_slack);
    
    if (verbose) {
      std::cout << "    Path " << path_idx << ":" << std::endl;
      std::cout << "      Original slack budget: " << slack_budget << std::endl;
      std::cout << "      RC delay increase: " << rc_delay_increase << std::endl;
      std::cout << "      Wire delay increase: " << wire_delay_increase << std::endl;
      std::cout << "      Cell delay increase: " << cell_delay_increase << std::endl;
      std::cout << "      Total used slack: " << used_slack << std::endl;
      std::cout << "      Remaining slack: " << (slack_budget - used_slack) << std::endl;
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
    std::cout << "\n  [FanIn Pin: " << d_pin->getMTerm()->getName() << "]" << std::endl;
  }
  
  // Check timing constraints
  const auto budget_it = flop.pin_budgets_.find(d_pin);
  if (budget_it == flop.pin_budgets_.end() || budget_it->second.empty()) {
    if (verbose) {
      std::cout << "    No timing constraints - skipping" << std::endl;
    }
    return;
  }
  
  if (verbose) {
    std::cout << "    Total timing paths: " << budget_it->second.size() << std::endl;
  }
  
  // Validate net connection
  odb::dbNet* fi_net = d_pin->getNet();
  if (!fi_net) {
    if (verbose) {
      std::cout << "    Pin not connected to net - skipping" << std::endl;
    }
    return;
  }
  
  // Get driver pin
  odb::dbITerm* fi_net_drvr_pin = fi_net->get1stITerm();
  if (!fi_net_drvr_pin) {
    if (verbose) {
      std::cout << "    No driver pin found - skipping" << std::endl;
    }
    return;
  }
  
  // Calculate D pin capacitance
  const sta::Pin* ipin_sta = network_->dbToSta(d_pin);
  const float ipin_cap = getPinCapacitance(ipin_sta);
  
  // Extract cell delay characterization
  odb::dbInst* fi_inst = fi_net_drvr_pin->getInst();
  const auto slews = getInstanceInputSlews(fi_inst);
  
  std::vector<std::pair<float, float>> cap_delay;
  float worst_delay = 0.0f;
  
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
      std::cout << "    Insufficient cap-delay data - skipping" << std::endl;
    }
    return;
  }
  
  // Calculate coefficients
  std::vector<float> coeffs;
  coeffs.reserve(cap_delay.size() - 1);
  for (size_t i = 1; i < cap_delay.size(); ++i) {
    const float dy = cap_delay[i].second - cap_delay[i-1].second;
    const float dx = cap_delay[i].first - cap_delay[i-1].first;
    coeffs.push_back(dy / dx);
  }
  
  // Build Steiner tree
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
  
  // Find path in Steiner tree
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
      std::cout << "    Target pin not found in Steiner tree - skipping" << std::endl;
    }
    return;
  }
  
  const int drvr_pt = fi_tree->drvrPt();
  std::vector<int> node_path;
  
  if (!findSteinerPathRecursive(fi_tree, drvr_pt, target_pt, node_path)) {
    if (verbose) {
      std::cout << "    Failed to find path in Steiner tree - skipping" << std::endl;
    }
    return;
  }
  
  if (node_path.empty() || node_path.size() < 2) {
    if (verbose) {
      std::cout << "    Path too short - skipping" << std::endl;
    }
    return;
  }
  
  // Calculate path segments
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
  
  const float l1 = dbuToMeters(final_segment_length_dbu);
  const float on_path_R_wo_last = unit_r * dbuToMeters(pre_leaf_length_dbu);
  
  // Calculate actual MANHATTAN distance to top Steiner point after movement
  const float dx_new = std::abs(flop.curr_pt_.x - top_x);
  const float dy_new = std::abs(flop.curr_pt_.y - top_y);
  const float actual_dist_dbu = dx_new + dy_new;
  const float actual_dist = dbuToMeters(actual_dist_dbu);
  
  if (verbose) {
    std::cout << "    Top Steiner point: (" << top_x << ", " << top_y << ")" << std::endl;
    std::cout << "    Original final segment (l1): " << l1 << " m" << std::endl;
    std::cout << "    Actual distance: " << actual_dist << " m" << std::endl;
    std::cout << "    Distance change: " << (actual_dist - l1) << " m" << std::endl;
    std::cout << "    Pre-leaf resistance: " << on_path_R_wo_last << std::endl;
  }
  
  // Calculate used slack for each timing path
  std::vector<std::pair<int, sta::Slack>> pin_used_slacks;
  pin_used_slacks.reserve(budget_it->second.size());
  
  for (const auto& [path_idx, slack_budget] : budget_it->second) {
    // Find appropriate coefficient
    float coeff = coeffs.empty() ? 0.0f : coeffs.back();
    const float new_cap = total_cap + (actual_dist - l1) * unit_c;
    
    for (size_t seg_idx = 0; seg_idx < coeffs.size(); ++seg_idx) {
      if (new_cap >= cap_delay[seg_idx].first && 
          new_cap <= cap_delay[seg_idx + 1].first) {
        coeff = coeffs[seg_idx];
        break;
      }
    }
    
    // Calculate delay increase components
    // 1. RC delay (quadratic term)
    const float rc_delay_increase = 
        (actual_dist * actual_dist - l1 * l1) * unit_r * unit_c;
    
    // 2. On-path resistance interaction with new wire capacitance
    const float on_path_delay_increase = 
        (actual_dist - l1) * on_path_R_wo_last * unit_c;
    
    // 3. Cell delay increase due to capacitance change
    const float cell_delay_increase = 
        coeff * (actual_dist - l1) * unit_c;
    
    // 4. Input pin capacitance interaction with new wire resistance
    const float ipin_delay_increase = 
        (actual_dist - l1) * unit_r * ipin_cap;
    
    // Total delay increase = used slack
    const float total_delay_increase = 
        rc_delay_increase + on_path_delay_increase + 
        cell_delay_increase + ipin_delay_increase;
    
    const sta::Slack used_slack = total_delay_increase;
    
    pin_used_slacks.emplace_back(path_idx, used_slack);
    
    if (verbose) {
      std::cout << "    Path " << path_idx << ":" << std::endl;
      std::cout << "      Original slack budget: " << slack_budget << std::endl;
      std::cout << "      RC delay increase: " << rc_delay_increase << std::endl;
      std::cout << "      On-path delay increase: " << on_path_delay_increase << std::endl;
      std::cout << "      Cell delay increase: " << cell_delay_increase << std::endl;
      std::cout << "      Input pin delay increase: " << ipin_delay_increase << std::endl;
      std::cout << "      Total used slack: " << used_slack << std::endl;
      std::cout << "      Remaining slack: " << (slack_budget - used_slack) << std::endl;
    }
  }
  
  used_slacks[d_pin] = std::move(pin_used_slacks);
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

// Implement

void 
AggloCluster::implementSingleCluster(const FlopCluster& cluster)
{
  // 1. 클러스터에 속한 1-bit FF들의 (D, Q, QN) 넷 번들 정보를 가져옵니다.
  std::vector<NetBundle> net_bundles = getNetBundles(cluster);

  if (net_bundles.empty()) {
    if (verbose_) {
      std::cout << "[implementSingleCluster] Cluster " << cluster.id_ 
                << " (" << cluster.flops_.size() << " flops) has no valid net bundles. Skipping." << std::endl;
    }
    return;
  }
  
  // 2. 헝가리안 알고리즘으로 HPWL이 최소가 되는 
  //    '최적 마스터'와 '넷-핀 매핑'을 찾습니다.
  MasterPortAssignment best_assignment = findBestMasterAndAssignment(cluster, net_bundles);

  if (best_assignment.best_master_ == nullptr) {
    if (verbose_) {
      std::cout << "[implementSingleCluster] Cluster " << cluster.id_
                << " (" << cluster.flops_.size() << " flops) found no valid master or assignment. Skipping." << std::endl;
    }
    return;
  }

  if (verbose_) {
    std::cout << "[implementSingleCluster] Cluster " << cluster.id_ 
              << ": Best master found: " << best_assignment.best_master_->getName()
              << ". Assignment cost: " << best_assignment.min_cost_ << std::endl;
  }

  // 3. 찾은 최적의 결과를 바탕으로 실제 DB를 수정합니다.
  //    (기존 Inst 삭제, 새 Inst 생성, 넷 연결)
  applyImplementation(cluster, best_assignment, net_bundles);
}

std::vector<NetBundle> 
AggloCluster::getNetBundles(const FlopCluster& cluster) const
{
  throw std::logic_error("getNetBundles not implemented");
  // std::vector<NetBundle> net_bundles;
  // net_bundles.reserve(cluster.flops_.size());

  // for (int flop_id : cluster.flops_) {
  //   const FlopUnit& flop = flop_units_[flop_id];
  //   NetBundle bundle;
  //   bundle.flop_id_ = flop_id;
  //   bundle.flop_inst_ = flop.inst_;

  //   for (odb::dbITerm* iterm : flop.inst_->getITerms()) {
  //     if (isDPin(iterm)) {
  //       bundle.d_net_ = iterm->getNet();
  //       bundle.d_iterm_ = iterm;
  //     } else if (isQPin(iterm)) {
  //       bundle.q_net_ = iterm->getNet();
  //       bundle.q_iterm_ = iterm;
  //     } else if (isQNPin(iterm)) {
  //       bundle.qn_net_ = iterm->getNet();
  //       bundle.qn_iterm_ = iterm;
  //     }
  //   }
    
  //   // D와 Q(또는 QN) 넷이 모두 존재해야 유효한 번들로 간주
  //   if (bundle.d_net_ && (bundle.q_net_ || bundle.qn_net_)) {
  //       net_bundles.push_back(bundle);
  //   } else {
  //       log_->warn(utl::GPL, 9970, "Flop {} ({}) is missing D/Q/QN nets. Excluding from bundle.",
  //                  flop.inst_->getName(), flop_id);
  //   }
  // }
  // return net_bundles;
}

MasterPortAssignment 
AggloCluster::findBestMasterAndAssignment(
    const FlopCluster& cluster,
    const std::vector<NetBundle>& net_bundles)
{
  throw std::logic_error("findBestMasterAndAssignment not implemented");
  // typedef util::StaticGraph<> Graph;

  // MasterPortAssignment best_result;
  // const int num_nets = net_bundles.size();
  
  // if (compatible_masters_.find(cluster.master_mask_) == compatible_masters_.end() ||
  //     compatible_masters_.at(cluster.master_mask_).find(num_nets) == compatible_masters_.at(cluster.master_mask_).end()) {
  //     log_->warn(utl::GPL, 9971, "No compatible {}-bit masters found for cluster {}.", num_nets, cluster.id_);
  //     return best_result;
  // }
  
  // const auto& candidate_masters = compatible_masters_.at(cluster.master_mask_).at(num_nets);
  
  // Point inst_center(std::lround(cluster.curr_pt_.x), std::lround(cluster.curr_pt_.y));

  // for (odb::dbMaster* master : candidate_masters) {
    
  //   std::vector<PortBundle> port_bundles = getPortBundles(master);
  //   const int num_ports = port_bundles.size();

  //   if (num_ports < num_nets) {
  //     log_->warn(utl::GPL, 9972, "Master {} has only {} port bundles, but cluster needs {}. Skipping.",
  //                master->getName(), num_ports, num_nets);
  //     continue;
  //   }
    
  //   const int num_left_nodes = num_nets;
  //   const int num_nodes = num_nets + num_ports;
  //   const int num_arcs = num_nets * num_ports; 
    
  //   Graph graph(num_nodes, num_arcs);
  //   std::vector<int64_t> arc_costs;
  //   arc_costs.reserve(num_arcs);

  //   for (int i = 0; i < num_nets; ++i) {
  //     for (int j = 0; j < num_ports; ++j) {
  //       int tail = i; 
  //       int head = num_left_nodes + j;
  //       graph.AddArc(tail, head);
        
  //       // [호출부 수정] 'master' 포인터를 calcAssignmentCost로 전달
  //       double cost_double = calcAssignmentCost(net_bundles[i], 
  //                                               port_bundles[j], 
  //                                               master,           // <-- [수정] master 전달
  //                                               inst_center);
        
  //       arc_costs.push_back(static_cast<int64_t>(std::round(cost_double)));
  //     }
  //   }

  //   graph.Build();
  //   ::operations_research::LinearSumAssignment assignment_solver(graph, num_left_nodes);

  //   for (int arc = 0; arc < num_arcs; ++arc) {
  //     assignment_solver.SetArcCost(arc, arc_costs[arc]);
  //   }
    
  //   if (assignment_solver.ComputeAssignment()) {
  //     double total_cost = assignment_solver.GetCost();
      
  //     if (total_cost < best_result.min_cost_) {
  //       best_result.min_cost_ = total_cost;
  //       best_result.best_master_ = master;
  //       best_result.net_to_port_assignment_.resize(num_nets);
        
  //       for (int i = 0; i < num_left_nodes; ++i) {
  //         int assigned_global_node_idx = assignment_solver.GetMate(i);
  //         int j = assigned_global_node_idx - num_left_nodes;
  //         best_result.net_to_port_assignment_[i] = j;
  //       }
  //     }
  //   } else {
  //       log_->warn(utl::GPL, 109, "Hungarian assignment failed for master {}.", master->getName());
  //   }
  // }
  
  // return best_result;
}

/**
 * @brief (헝가리안 비용함수) 
 * 특정 넷 번들을 특정 포트 번들에 할당했을 때의 HPWL 비용을 계산합니다.
 */
double 
AggloCluster::calcAssignmentCost(const NetBundle& net_bundle,
                                 const PortBundle& port_bundle,
                                 odb::dbMaster* master,         // <-- [수정] master 받기
                                 const Point& new_inst_center) const
{
  throw std::logic_error("calcAssignmentCost not implemented");
  // double total_hpwl_cost = 0.0;

  // // 1. 'D' 핀 비용 계산
  // if (net_bundle.d_net_ && port_bundle.d_mterm_) {
  //   odb::Rect bbox_d = getNetBBoxWithoutPin(net_bundle.d_net_, net_bundle.d_iterm_);
  //   Point global_d_pos = getGlobalMTermPos(port_bundle.d_local_pos_, master, new_inst_center);
  //   bbox_d.merge(global_d_pos);
  //   total_hpwl_cost += bbox_d.dx() + bbox_d.dy();
  // }

  // // 2. 'Q' 핀 비용 계산
  // if (net_bundle.q_net_ && port_bundle.q_mterm_) {
  //   odb::Rect bbox_q = getNetBBoxWithoutPin(net_bundle.q_net_, net_bundle.q_iterm_);
  //   Point global_q_pos = getGlobalMTermPos(port_bundle.q_local_pos_, master, new_inst_center);
  //   bbox_q.merge(global_q_pos);
  //   total_hpwl_cost += bbox_q.dx() + bbox_q.dy();
  // }

  // // 3. 'QN' 핀 비용 계산
  // if (net_bundle.qn_net_ && port_bundle.qn_mterm_) {
  //   odb::Rect bbox_qn = getNetBBoxWithoutPin(net_bundle.qn_net_, net_bundle.qn_iterm_);
  //   Point global_qn_pos = getGlobalMTermPos(port_bundle.qn_local_pos_, master, new_inst_center);
  //   bbox_qn.merge(global_qn_pos);
  //   total_hpwl_cost += bbox_qn.dx() + bbox_qn.dy();
  // }
  
  // return total_hpwl_cost;
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
    // [수정] 핀 위치가 아닌 인스턴스의 BBox 중심으로 HPWL 계산
    // (또는 getAvgXY 사용)
    int x, y;
    if (iterm->getAvgXY(&x, &y)) {
       bbox.merge(odb::Point(x,y));
    }
  }

  for (odb::dbBTerm* bterm : net->getBTerms()) {
    odb::Rect bterm_bbox = bterm->getBBox();
    if (!bterm_bbox.isInverted()) {
        bbox.merge(bterm_bbox);
    }
  }
  
  if (bbox.isInverted()) {
      bbox.set_xlo(0); bbox.set_ylo(0); bbox.set_xhi(0); bbox.set_yhi(0);
  }

  return bbox;
}

Point 
AggloCluster::getGlobalMTermPos(const Point& local_port_pos, 
                                odb::dbMaster* master,         // <-- [수정] master 받기
                                const Point& inst_center) const
{
  // [로직 수정]
  // local_port_pos는 마스터의 원점(0,0) 기준 핀의 로컬 좌표입니다.
  // inst_center는 새 인스턴스의 중심 좌표입니다.

  throw std::logic_error("getGlobalMTermPos not implemented");

  // // 1. 인스턴스의 원점(Origin)을 계산합니다.
  // int center_x = inst_center.get<0>();
  // int center_y = inst_center.get<1>();
  // int origin_x = center_x - master->getWidth() / 2;
  // int origin_y = center_y - master->getHeight() / 2;

  // // 2. 핀의 글로벌 좌표 = 인스턴스 원점 + 핀의 로컬 좌표
  // return Point(origin_x + local_port_pos.get<0>(), 
  //              origin_y + local_port_pos.get<1>());
}

// /**
//  * @brief 마스터의 로컬 MTerm 좌표를 글로벌 좌표(DBU)로 변환합니다.
//  */
// Point 
// AggloCluster::getGlobalMTermPos(const Point& local_port_pos, 
//                                 const Point& inst_center) const
// {
//   throw std::logic_error("getGlobalMTermPos not implemented");
//   // [주의] inst_center는 인스턴스의 *중심*입니다. 
//   //        local_port_pos는 마스터의 *원점(0,0)* 기준입니다.
//   //        정확한 계산을 위해서는 마스터의 (width/2, height/2) 오프셋이 필요합니다.
//   //        하지만, 여기서는 모든 MTerm에 동일한 오프셋이 적용되므로
//   //        inst_center를 임시 원점으로 사용해도 비용 계산(상대 비교)은 유효합니다.
  
//   // // [간단한 구현]
//   // // (inst_center가 (0,0)이라 가정할 때의 mterm 위치) + (실제 inst_center 위치)
//   // return Point(inst_center.get<0>() + local_port_pos.get<0>(),
//   //              inst_center.get<1>() + local_port_pos.get<1>());

//   // [더 정확한 구현 (master 포인터 필요)]
//   // int center_x = inst_center.get<0>();
//   // int center_y = inst_center.get<1>();
//   // int origin_x = center_x - master->getWidth() / 2;
//   // int origin_y = center_y - master->getHeight() / 2;
//   // return Point(origin_x + local_port_pos.get<0>(), 
//   //              origin_y + local_port_pos.get<1>());
// }

void 
AggloCluster::applyImplementation(const FlopCluster& cluster,
                                  const MasterPortAssignment& result,
                                  const std::vector<NetBundle>& net_bundles)
{

  throw std::logic_error("applyImplementation not implemented");
  // --- 이 함수는 실제 DB를 수정하므로 매우 주의해야 합니다 ---
  
  // odb::dbMaster* best_master = result.best_master_;
  // if (!best_master) {
  //     return;
  // }
  
  // int num_bits = net_bundles.size();
  
  // // 1. 새 MBFF 인스턴스 생성
  // std::string new_inst_name = "mbff_cluster_" + std::to_string(cluster.id_);
  // odb::dbInst* new_inst = odb::dbInst::create(block_, best_master, new_inst_name.c_str());
  
  // if (!new_inst) {
  //     log_->error(utl::GPL, 9974, "Failed to create new instance: {}", new_inst_name);
  //     return;
  // }

  // // 2. 새 인스턴스 위치 설정 (중심점 기준)
  // int center_x = std::lround(cluster.curr_pt_.x);
  // int center_y = std::lround(cluster.curr_pt_.y);
  // int origin_x = center_x - best_master->getWidth() / 2;
  // int origin_y = center_y - best_master->getHeight() / 2;
  // new_inst->setLocation(origin_x, origin_y);
  // new_inst->setPlacementStatus(odb::dbPlacementStatus::PLACED);

  // // 3. 공통 넷 연결 (CLK, CLR, PRE, SE 등)
  // // (InstMask는 클러스터 내 모든 FF가 동일하다고 가정)
  // const InstMask& inst_mask = cluster.inst_mask_;
  
  // if (inst_mask.clock_net_) {
  //   odb::dbMTerm* clk_mterm = best_master->findMTerm("CK"); // "CK" 또는 "CLK" 등
  //   if (clk_mterm) new_inst->getITerm(clk_mterm)->connect(inst_mask.clock_net_);
  // }
  // if (inst_mask.clear_net_) {
  //   odb::dbMTerm* clr_mterm = best_master->findMTerm("CLR"); // "CLR" 또는 "R" 등
  //   if (clr_mterm) new_inst->getITerm(clr_mterm)->connect(inst_mask.clear_net_);
  // }
  // if (inst_mask.preset_net_) {
  //   odb::dbMTerm* pre_mterm = best_master->findMTerm("PRE"); // "PRE" 또는 "S" 등
  //   if (pre_mterm) new_inst->getITerm(pre_mterm)->connect(inst_mask.preset_net_);
  // }
  // if (inst_mask.scan_enable_net_) {
  //   odb::dbMTerm* se_mterm = best_master->findMTerm("SE"); // "SE" 등
  //   if (se_mterm) new_inst->getITerm(se_mterm)->connect(inst_mask.scan_enable_net_);
  // } 
  // if (inst_mask.scan_in_net_) {
  //   odb::dbMTerm* si_mterm = best_master->findMTerm("SI"); // "SE" 등
  //   if (si_mterm) new_inst->getITerm(si_mterm)->connect(inst_mask.scan_in_net_);    
  // }
  
  // // 4. 데이터/스캔 넷 연결 (헝가리안 할당 결과 기반)
  // std::vector<PortBundle> port_bundles = getPortBundles(best_master);

  // for (int i = 0; i < num_bits; ++i) { // i = 넷 번들 인덱스
  //   int j = result.net_to_port_assignment_[i]; // j = 포트 번들 인덱스
    
  //   const NetBundle& net_b = net_bundles[i];
  //   const PortBundle& port_b = port_bundles[j]; // [j]가 포트 번들 인덱스

  //   // D 핀 연결
  //   if (net_b.d_net_ && port_b.d_mterm_) {
  //     net_b.d_iterm_->disconnect(); // 기존 핀 연결 해제
  //     new_inst->getITerm(port_b.d_mterm_)->connect(net_b.d_net_);
  //   }
  //   // Q 핀 연결
  //   if (net_b.q_net_ && port_b.q_mterm_) {
  //     net_b.q_iterm_->disconnect(); 
  //     new_inst->getITerm(port_b.q_mterm_)->connect(net_b.q_net_);
  //   }
  //   // QN 핀 연결
  //   if (net_b.qn_net_ && port_b.qn_mterm_) {
  //     net_b.qn_iterm_->disconnect();
  //     new_inst->getITerm(port_b.qn_mterm_)->connect(net_b.qn_net_);
  //   }
  //   // (필요시 SI/SO 핀 로직 추가)
  // }

  // // 5. 기존 1-bit 인스턴스 삭제
  // for (const auto& net_b : net_bundles) {
  //   odb::dbInst::destroy(net_b.inst_);
  // }
}

std::vector<PortBundle> 
AggloCluster::getPortBundles(odb::dbMaster* master) const
{
  throw std::logic_error("getPortBundles not implemented");
  // // 정규식: (D, Q, QN, SI, SO)로 시작하고 (그룹 1)
  // //         (\d+) : 1개 이상의 숫자로 끝남 (그룹 2)
  // std::regex pin_regex(R"(^(D|Q|QN)(\d+)$)");
  // std::smatch match;

  // // Key: 핀 이름에서 추출한 인덱스 (e.g., "QN0" -> 0, "D1" -> 1)
  // std::map<int, PortBundle> port_map; 

  // for (odb::dbMTerm* mterm : master->getMTerms()) {
  //   std::string mterm_name = mterm->getName();
    
  //   if (std::regex_match(mterm_name, match, pin_regex) && match.size() == 3) {
      
  //     std::string type = match[1].str(); // "D", "Q", "QN" 등
      
  //     // [수정] 핀 이름의 숫자를 정수로 변환하여 인덱스로 바로 사용
  //     int index = std::stoi(match[2].str()); // e.g., "0", "1", "2" ...

  //     PortBundle& bundle = port_map[index]; // 맵에 접근 (없으면 생성)
  //     bundle.bundle_idx_ = index;

  //     // 핀의 로컬 좌표 계산 (MTerm의 BBox 중심 사용)
  //     Point local_pos(0, 0);
  //     odb::Rect bbox;
  //     if (mterm->getBBox(bbox)) { 
  //         local_pos.set<0>(bbox.xCenter());
  //         local_pos.set<1>(bbox.yCenter());
  //     } 

  //     if (type == "D") {
  //       bundle.d_mterm_ = mterm;
  //       bundle.d_local_pos_ = local_pos;
  //     } else if (type == "Q") {
  //       bundle.q_mterm_ = mterm;
  //       bundle.q_local_pos_ = local_pos;
  //     } else if (type == "QN") {
  //       bundle.qn_mterm_ = mterm;
  //       bundle.qn_local_pos_ = local_pos;
  //     }
  //   }
  // }

  // std::vector<PortBundle> port_bundles;
  // port_bundles.reserve(port_map.size());
  // for (auto const& [index, bundle] : port_map) {
  //   port_bundles.push_back(bundle);
  // }
  // return port_bundles;
}





// float AggloCluster::calcHPWL(const std::vector<int>& flop_indices,
//                                   const std::vector<odb::Point>& cluster_centers) const
// {
//   std::set<odb::dbNet*> nets;
//   std::set<odb::dbInst*> cluster_insts;
//   for (int flop_idx : flop_indices) {
//     cluster_insts.insert(flop_units_[flop_idx].inst_);
//   }

//   auto pins = ConnectedPins(flop_indices);
//   for (const auto& pin_variant : pins) {
//     std::visit(
//         [&](auto&& arg) {
//           if (arg->getNet()) {
//             nets.insert(arg->getNet());
//           }
//         },
//         pin_variant);
//   }

//   float total_hpwl = 0;
//   for (odb::dbNet* net : nets) {
//     odb::Rect net_bbox;
//     net_bbox.mergeInit();

//     for (odb::dbITerm* iterm : net->getITerms()) {
//       // Exclude ITerms belonging to the flops in the cluster
//       if (cluster_insts.find(iterm->getInst()) == cluster_insts.end()) {
//         // Use instance's BBox for HPWL calculation, not pin's BBox
//         net_bbox.merge(iterm->getInst()->getBBox()->getBox());
//       }
//     }
//     for (odb::dbBTerm* bterm : net->getBTerms()) {
//       for (auto bpin : bterm->getBPins()) {
//         net_bbox.merge(bpin->getBBox());
//       }
//     }

//     for (const auto& center : cluster_centers) {
//       net_bbox.merge(center);
//     }

//     total_hpwl += net_bbox.dx() + net_bbox.dy();
//   }

//   return total_hpwl;
// }

// std::set<std::variant<odb::dbITerm*, odb::dbBTerm*>>
// AggloCluster::getConnectedPins(const std::vector<int>& flop_indices) const
// {
//   std::set<std::variant<odb::dbITerm*, odb::dbBTerm*>> connected_pins;
//   std::set<odb::dbInst*> cluster_insts;
//   for (int flop_idx : flop_indices) {
//     cluster_insts.insert(flop_units_[flop_idx].inst_);
//   }

//   for (int flop_idx : flop_indices) {
//     const FlopUnit& flop = flop_units_[flop_idx];
//     for (odb::dbITerm* iterm : flop.inst_->getITerms()) {
//       odb::dbNet* net = iterm->getNet();
//       if (net == nullptr
//           || (net->getSigType() != odb::dbSigType::SIGNAL
//               && net->getSigType() != odb::dbSigType::CLOCK)) {
//         continue;
//       }
//       for (odb::dbITerm* net_iterm : net->getITerms()) {
//         if (cluster_insts.find(net_iterm->getInst()) == cluster_insts.end()) {
//           connected_pins.insert(net_iterm);
//         }
//       }
//       for (odb::dbBTerm* net_bterm : net->getBTerms()) {
//         connected_pins.insert(net_bterm);
//       }
//     }
//   }
//   return connected_pins;
// }



// PlacementResult AggloCluster::calcPlacementLoc(int cluster_idx1, int cluster_idx2)
// {
//   PlacementResult result;
  
//   const Box& box1 = flop_clusters_[cluster_idx1].feasible_region_;
//   const Box& box2 = flop_clusters_[cluster_idx2].feasible_region_;

//   Box intersection_box;
//   boost::geometry::intersection(box1, box2, intersection_box);

//   if (boost::geometry::is_empty(intersection_box)) {
//     return result;
//   }

//   std::vector<int> merged_flops;
//   merged_flops.insert(merged_flops.end(),
//                       flop_clusters_[cluster_idx1].flops_.begin(),
//                       flop_clusters_[cluster_idx1].flops_.end());
//   merged_flops.insert(merged_flops.end(),
//                       flop_clusters_[cluster_idx2].flops_.begin(),
//                       flop_clusters_[cluster_idx2].flops_.end());
//   auto connected_pins = getConnectedPins(merged_flops);

//   // Calculate original HPWL
//   const auto& c1 = flop_clusters_[cluster_idx1];
//   const auto& c2 = flop_clusters_[cluster_idx2];
//   std::vector<odb::Point> original_centers;
//   original_centers.emplace_back(c1.curr_pt_.x, c1.curr_pt_.y);
//   original_centers.emplace_back(c2.curr_pt_.x, c2.curr_pt_.y);
//   float original_hpwl = calcHPWL(merged_flops, original_centers);

//   // Remove cluster's own pins to get external pins for median calculation
//   // This part of your logic seems to be for finding the best placement location
//   // by considering external connections.

//   std::vector<odb::Point> connected_pin_coords;
//   for (const auto& pin_variant : connected_pins) {
//     std::visit(
//         [&](auto&& arg) {
//           using T = std::decay_t<decltype(arg)>;
//           odb::Rect bbox;
//           if constexpr (std::is_same_v<T, odb::dbITerm*>) {
//             bbox = arg->getInst()->getBBox()->getBox();
//           } else if constexpr (std::is_same_v<T, odb::dbBTerm*>) {
//             for (auto bpin : arg->getBPins()) {
//               bbox.merge(bpin->getBBox());
//             }
//           }
//           connected_pin_coords.emplace_back(bbox.xCenter(), bbox.yCenter());
//         },
//         pin_variant);
//   }

//   // 3. Calculate the HPWL-minimizing median point
//   FloatPoint median_point;
//   if (connected_pin_coords.empty()) {
//     median_point = {0.0, 0.0};
//   } else {
//     std::vector<float> x_coords, y_coords;
//     x_coords.reserve(connected_pin_coords.size());
//     y_coords.reserve(connected_pin_coords.size());

//     for (const auto& pt : connected_pin_coords) {
//       x_coords.push_back(pt.getX());
//       y_coords.push_back(pt.getY());
//     }

//     size_t mid_idx = x_coords.size() / 2;
//     std::nth_element(x_coords.begin(), x_coords.begin() + mid_idx, x_coords.end());
//     std::nth_element(y_coords.begin(), y_coords.begin() + mid_idx, y_coords.end());
//     median_point = {x_coords[mid_idx], y_coords[mid_idx]};
//   }

//   // 4. Find the closest point in the intersection_box to the median_point
//   // Transform median_point to the rotated coordinate system
//   Point median_pt_transformed = transformCoords(Point(median_point.x, median_point.y));
//   Point closest_pt_bg;

//   // Find the closest point in the rotated space
//   if (boost::geometry::within(median_pt_transformed, intersection_box)) {
//     closest_pt_bg = median_pt_transformed;
//   } else {
//     float closest_x = std::max((float) intersection_box.min_corner().get<0>(),
//                    std::min((float) median_pt_transformed.get<0>(), (float) intersection_box.max_corner().get<0>()));
//     float closest_y = std::max((float) intersection_box.min_corner().get<1>(),
//                    std::min((float) median_pt_transformed.get<1>(), (float) intersection_box.max_corner().get<1>()));
//     closest_pt_bg.set<0>(closest_x);
//     closest_pt_bg.set<1>(closest_y);
//   }

//   // Transform the result back to the original coordinate system
//   closest_pt_bg = inverseTransformCoords(closest_pt_bg);
//   result.placement_loc = {static_cast<float>(closest_pt_bg.get<0>()),
//                           static_cast<float>(closest_pt_bg.get<1>())};
//   // Calculate merged HPWL
//   float merged_hpwl = calcHPWL(merged_flops,
//                                {odb::Point(result.placement_loc.x, result.placement_loc.y)}); 

//   result.gain = static_cast<double>(original_hpwl - merged_hpwl);
//   result.valid = true;
//   return result;
// }





}  // end namespace gpl
