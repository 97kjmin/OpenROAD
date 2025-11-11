#include "aggloCluster.h"

#include <algorithm>
#include <cmath>
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
  
  exit(0);

  createFlopClusters();

  createCompatibilityGraph();

  runAgglomerativeClustering();

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
  constexpr int DEBUG_SAMPLE_SIZE = 10;

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
  flop_clusters_.clear();
  flop_clusters_.reserve(flop_units_.size());
  feasible_regions_.clear();
  flop_cluster_is_valid_.assign(flop_units_.size(), true);
  flop_cluster_no_further_merge_.assign(flop_units_.size(), false);

  for (size_t i = 0; i < flop_units_.size(); ++i) {
    flop_clusters_.emplace_back(i, flop_units_[i]);
    flop_units_[i].cluster_idx_ = i;

    const Box& feasible_region = flop_units_[i].feasible_region_;
    feasible_regions_.insert(std::make_pair(feasible_region, i));
  }
}

//==============================================================================
// Phase 9: createCompatibilityGraph()
//==============================================================================

void
AggloCluster::createCompatibilityGraph()
{
  edge_pq_.clear();
  adj_list_.clear();

  for (size_t i = 0; i < flop_clusters_.size(); ++i) {
    updateEdges(flop_clusters_[i]);
  }
}

//==============================================================================
// Phase 10: runAgglomerativeClustering()
//==============================================================================

void 
AggloCluster::runAgglomerativeClustering()
{
  while (!edge_pq_.empty()) {
    Edge best_edge = *edge_pq_.rbegin();
    edge_pq_.erase(std::prev(edge_pq_.end()));

    int n1_idx = best_edge.n1;
    int n2_idx = best_edge.n2;

    if (!flop_cluster_is_valid_[n1_idx] || !flop_cluster_is_valid_[n2_idx]) {
      continue;
    }

    // mergeClusters create new cluster and mark n1 and n2 as invalid
    int n_new_idx = mergeClusters(best_edge);
    const FlopCluster& n_new = flop_clusters_[n_new_idx];

    if (isFurtherMergeable(n_new)) {
      // Slack distribution is not performed, so just update edges
      updateEdges(n_new); 
    } 
    else {
      flop_cluster_no_further_merge_[n_new_idx] = true;

      // Perform slack distribution if final cluster is formed      
      std::set<int> affected_flop_units = distributeSlack(n_new); 
      
      std::set<int> affected_flop_clusters;
      
      for (int u_idx : affected_flop_units) {
        // Update feasible region of affected 'FlopUnit's
        calcFeasibleRegion(flop_units_[u_idx]); 
        
        // Get affected 'Cluster's for later FR and edge updates
        int c_idx = flop_units_[u_idx].cluster_idx_;
        if (flop_cluster_is_valid_[c_idx] && !flop_cluster_no_further_merge_[c_idx]) {
          affected_flop_clusters.insert(c_idx);
        }
      }

      // Update feasible regions of affected 'Cluster's
      for (int c_idx : affected_flop_clusters) {
        updateFeasibleRegion(flop_clusters_[c_idx]);
      }

      // Update edges of affected 'Cluster's
      for (int c_idx : affected_flop_clusters) {
        updateEdges(flop_clusters_[c_idx]);
      }
    }
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
AggloCluster::calcFeasibleRegion(FlopUnit& flop, 
                                 bool verbose,
                                 double unit_r, 
                                 double unit_c)
{
  //----------------------------------------------------------------------------
  // Step 1: Get parasitic estimator
  //----------------------------------------------------------------------------
  auto est_ = resizer_->getEstimateParasitics();
  
  if (verbose) {
    std::cout << "  [Step 1] Wire Parameters" << std::endl;
    std::cout << "    Unit resistance: " << unit_r << " Ohm/DBU" << std::endl;
    std::cout << "    Unit capacitance: " << unit_c << " F/DBU" << std::endl;
  }
  
  //----------------------------------------------------------------------------
  // Step 2: Classify pins
  //----------------------------------------------------------------------------
  PinClassification pins = classifyFlopPins(flop.inst_);
  
  if (verbose) {
    std::cout << "  [Step 2] Pin Classification" << std::endl;
    std::cout << "    D pins: " << pins.d_pins.size() << std::endl;
    for (size_t i = 0; i < pins.d_pins.size(); ++i) {
      std::cout << "      [" << i << "] " << pins.d_pins[i]->getMTerm()->getName() << std::endl;
    }
    std::cout << "    Q pins: " << pins.q_pins.size() << std::endl;
    for (size_t i = 0; i < pins.q_pins.size(); ++i) {
      std::cout << "      [" << i << "] " << pins.q_pins[i]->getMTerm()->getName() << std::endl;
    }
    std::cout << "    QN pins: " << pins.qn_pins.size() << std::endl;
    for (size_t i = 0; i < pins.qn_pins.size(); ++i) {
      std::cout << "      [" << i << "] " << pins.qn_pins[i]->getMTerm()->getName() << std::endl;
    }
    std::cout << "    Clock pin: " << (pins.clk_pin ? pins.clk_pin->getMTerm()->getName() : "NOT FOUND") << std::endl;
  }
  
  if (!pins.clk_pin) {
    if (verbose) {
      std::cout << "    [WARNING] No clock pin found - skipping this flop" << std::endl;
    }
    return;
  }
  
  odb::dbMTerm* clk_pin_lib = pins.clk_pin->getMTerm();
  
  //----------------------------------------------------------------------------
  // Step 3: Process Fan-Out pins (Q and QN)
  //----------------------------------------------------------------------------
  std::vector<odb::dbITerm*> output_pins;
  output_pins.insert(output_pins.end(), pins.q_pins.begin(), pins.q_pins.end());
  output_pins.insert(output_pins.end(), pins.qn_pins.begin(), pins.qn_pins.end());
  
  for (odb::dbITerm* out_pin : output_pins) {
    processFanOutPin(flop, out_pin, clk_pin_lib, est_, unit_r, unit_c, verbose);
  }

  //----------------------------------------------------------------------------
  // Step 4: Process Fan-In pins (D pins)
  //----------------------------------------------------------------------------
  for (odb::dbITerm* d_pin : pins.d_pins) {
    processFanInPin(flop, d_pin, est_, unit_r, unit_c, verbose);
  }
  
  //----------------------------------------------------------------------------
  // Step 5: Compute final feasible region
  //----------------------------------------------------------------------------
  std::vector<odb::dbITerm*> all_pins;
  all_pins.insert(all_pins.end(), pins.d_pins.begin(), pins.d_pins.end());
  all_pins.insert(all_pins.end(), pins.q_pins.begin(), pins.q_pins.end());
  all_pins.insert(all_pins.end(), pins.qn_pins.begin(), pins.qn_pins.end());
  
  computeFinalFeasibleRegion(flop, all_pins, verbose);
}


// -----------------------------------------------------------------------------
// [Level 2] Main processing functions for different pin types
// -----------------------------------------------------------------------------

AggloCluster::PinClassification
AggloCluster::classifyFlopPins(odb::dbInst* inst) const
{
  PinClassification result;
  result.clk_pin = nullptr;
  
  for (odb::dbITerm* iterm : inst->getITerms()) {
    if (isDPin(iterm)) {
      result.d_pins.push_back(iterm);
    } else if (isQPin(iterm)) {
      result.q_pins.push_back(iterm);
    } else if (isQNPin(iterm)) {
      result.qn_pins.push_back(iterm);
    } else if (isClockPin(iterm)) {
      result.clk_pin = iterm;
    }
  }
  return result;
}

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
  
  // Step 1: Check if pin has timing constraints
  auto budget_it = flop.pin_budgets_.find(out_pin);
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
  
  // Step 2: Get most critical path (first element after sorting)
  const auto& [critical_path_idx, slack_budget] = budget_it->second[0];
  
  if (verbose) {
    std::cout << "    Most critical path:" << std::endl;
    std::cout << "      Path index: " << critical_path_idx << std::endl;
    
    sta::Unit* time_unit = sta_->units()->timeUnit();
    const char* time_suffix = time_unit->scaledSuffix();
    std::cout << "      Slack budget: " << time_unit->asString(slack_budget, 3) << time_suffix << std::endl;
  }
  
  // Get the fan-out net
  odb::dbNet* fo_net = out_pin->getNet();
  if (!fo_net) {
    if (verbose) {
      std::cout << "    Pin not connected to any net" << std::endl;
      std::cout << "    Result: Skipping this pin" << std::endl;
    }
    return;
  }
  
  if (verbose) {
    std::cout << "    Fan-out net: " << fo_net->getName() << std::endl;
  }
  
  // Step 3: Extract cell delay coefficient
  odb::dbMTerm* opin_lib = out_pin->getMTerm();
  auto cap_delay = extractCapacitanceDelayPoints(flop.inst_, 
                                                   clk_pin_lib->getName(), 
                                                   opin_lib->getName(), 
                                                   0);
  float coeff = (cap_delay[1].second - cap_delay[0].second) / 
                (cap_delay[1].first - cap_delay[0].first);
  
  if (verbose) {
    sta::Unit* time_unit = sta_->units()->timeUnit();
    sta::Unit* cap_unit = sta_->units()->capacitanceUnit();
    const char* time_suffix = time_unit->scaledSuffix();
    const char* cap_suffix = cap_unit->scaledSuffix();
    
    std::cout << "    Cell delay analysis:" << std::endl;
    std::cout << "      Cap-delay points: " << cap_delay.size() << std::endl;
    std::cout << "      Point[0]: cap=" << cap_unit->asString(cap_delay[0].first) << cap_suffix
              << ", delay=" << time_unit->asString(cap_delay[0].second) << time_suffix << std::endl;
    std::cout << "      Point[1]: cap=" << cap_unit->asString(cap_delay[1].first) << cap_suffix
              << ", delay=" << time_unit->asString(cap_delay[1].second) << time_suffix << std::endl;
    std::cout << "      Delay coefficient (∂delay/∂cap): " << coeff << std::endl;
  }
  
  // Step 4: Build Steiner tree and calculate net delay
  sta::Pin* drvr_pin = network_->dbToSta(out_pin);
  est::SteinerTree* drvr_tree = est->makeSteinerTree(drvr_pin);
  
  // Net capacitance computation
  const sta::Net* drvr_dbnet = network_->dbToSta(fo_net);
  auto p_cap = 0.0f, w_cap = 0.0f;
  const sta::MinMax* mm = sta::MinMax::max();
  sta_->connectedCap(drvr_dbnet, corner_, mm, p_cap, w_cap);
  auto c_total = p_cap + w_cap;
  
  // Get distance from driver to first Steiner point
  auto top_pt = drvr_tree->top();
  auto drvr_pt = drvr_tree->drvrPt();
  auto top_loc = drvr_tree->location(top_pt);
  auto drvr_loc = drvr_tree->location(drvr_pt);
  float l1 = drvr_tree->distance(drvr_pt, top_pt) / pow(10, 6);
  
  // Wire resistance and capacitance up to first Steiner point
  auto wire_cap = l1 * unit_c;
  auto wire_res = l1 * unit_r;
  auto wo_fst_stt_cap = c_total - wire_cap;
  
  if (verbose) {
    sta::Unit* cap_unit = sta_->units()->capacitanceUnit();
    sta::Unit* res_unit = sta_->units()->resistanceUnit();
    const char* cap_suffix = cap_unit->scaledSuffix();
    const char* res_suffix = res_unit->scaledSuffix();
    
    std::cout << "    Steiner tree analysis:" << std::endl;
    std::cout << "      Driver location: (" << drvr_loc.getX() << ", " << drvr_loc.getY() << ")" << std::endl;
    std::cout << "      Top Steiner point: (" << top_loc.getX() << ", " << top_loc.getY() << ")" << std::endl;
    std::cout << "      Distance to Steiner (l1): " << l1 << " mm" << std::endl;
    std::cout << "    Net capacitance:" << std::endl;
    std::cout << "      Pin capacitance: " << cap_unit->asString(p_cap) << cap_suffix << std::endl;
    std::cout << "      Wire capacitance: " << cap_unit->asString(w_cap) << cap_suffix << std::endl;
    std::cout << "      Total capacitance: " << cap_unit->asString(c_total) << cap_suffix << std::endl;
    std::cout << "    Wire parasitics (to Steiner):" << std::endl;
    std::cout << "      Wire capacitance: " << cap_unit->asString(wire_cap) << cap_suffix << std::endl;
    std::cout << "      Wire resistance: " << res_unit->asString(wire_res) << res_suffix << std::endl;
    std::cout << "      Remaining capacitance: " << cap_unit->asString(wo_fst_stt_cap) << cap_suffix << std::endl;
  }
  
  // Step 5: Solve for max distance
  float max_dist = solveMaxDistanceFanOut(l1, unit_r, unit_c, slack_budget, coeff, wo_fst_stt_cap, verbose);
  
  // Step 6: Create feasible region box
  int steiner_x = top_loc.getX();
  int steiner_y = top_loc.getY();
  
  if (verbose) {
    std::cout << "    Maximum placement distance:" << std::endl;
    std::cout << "      Max Manhattan distance from Steiner: " << max_dist << " DBU" << std::endl;
    std::cout << "      Max distance in mm: " << (max_dist / pow(10, 6)) << " mm" << std::endl;
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
  
  // Step 1: Check if pin has timing constraints
  auto budget_it = flop.pin_budgets_.find(d_pin);
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
  
  // Step 2: Get most critical path
  const auto& [critical_path_idx, slack_budget] = budget_it->second[0];
  
  if (verbose) {
    std::cout << "    Most critical path:" << std::endl;
    std::cout << "      Path index: " << critical_path_idx << std::endl;
    
    sta::Unit* time_unit = sta_->units()->timeUnit();
    const char* time_suffix = time_unit->scaledSuffix();
    std::cout << "      Slack budget: " << time_unit->asString(slack_budget, 3) << time_suffix << std::endl;
  }
  
  // Get the fan-in net
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
  
  // Get driver pin of fan-in net
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
  
  // Step 3: Calculate D pin capacitance
  const sta::Pin* ipin_sta = network_->dbToSta(d_pin);
  float ipin_cap = getPinCapacitance(ipin_sta);
  
  if (verbose) {
    sta::Unit* cap_unit = sta_->units()->capacitanceUnit();
    const char* cap_suffix = cap_unit->scaledSuffix();
    std::cout << "    D pin input capacitance: " << cap_unit->asString(ipin_cap) << cap_suffix << std::endl;
  }
  
  // Step 4: Extract cell delay coefficient
  odb::dbInst* fi_inst = fi_net_drvr_pin->getInst();
  std::vector<std::pair<odb::dbITerm*, float>> slews = getInstanceInputSlews(fi_inst);
  
  std::vector<std::pair<float, float>> cap_delay;
  float worst_delay = 0.0;
  
  if (verbose) {
    std::cout << "    Driver cell delay analysis:" << std::endl;
    std::cout << "      Input pins with slews: " << slews.size() << std::endl;
  }
  
  for (auto& [iterm, slew] : slews) {
    auto pts = extractCapacitanceDelayPoints(
        fi_inst,
        iterm->getMTerm()->getName(),
        fi_net_drvr_pin->getMTerm()->getName(),
        slew);
    
    if (!pts.empty()) {
      float tail_second = pts.back().second;
      if (tail_second >= worst_delay) {
        worst_delay = tail_second;
        cap_delay = std::move(pts);
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
  
  float coeff = (cap_delay[1].second - cap_delay[0].second) / 
                (cap_delay[1].first - cap_delay[0].first);
  
  if (verbose) {
    sta::Unit* time_unit = sta_->units()->timeUnit();
    sta::Unit* cap_unit = sta_->units()->capacitanceUnit();
    const char* time_suffix = time_unit->scaledSuffix();
    const char* cap_suffix = cap_unit->scaledSuffix();
    
    std::cout << "      Cap-delay points: " << cap_delay.size() << std::endl;
    std::cout << "      Point[0]: cap=" << cap_unit->asString(cap_delay[0].first) << cap_suffix
              << ", delay=" << time_unit->asString(cap_delay[0].second) << time_suffix << std::endl;
    std::cout << "      Point[1]: cap=" << cap_unit->asString(cap_delay[1].first) << cap_suffix
              << ", delay=" << time_unit->asString(cap_delay[1].second) << time_suffix << std::endl;
    std::cout << "      Delay coefficient (∂delay/∂cap): " << coeff << std::endl;
  }
  
  // Step 5: Build Steiner tree for fan-in net
  const sta::Pin* fi_net_drvr_pin_sta = network_->dbToSta(fi_net_drvr_pin);
  est::SteinerTree* fi_tree = est->makeSteinerTree(fi_net_drvr_pin_sta);
  
  auto top_pt = fi_tree->top();
  auto top_loc = fi_tree->location(top_pt);
  const int top_x = top_loc.getX();
  const int top_y = top_loc.getY();
  
  if (verbose) {
    std::cout << "    Steiner tree analysis:" << std::endl;
    std::cout << "      Top Steiner point: (" << top_x << ", " << top_y << ")" << std::endl;
    std::cout << "      Branch count: " << fi_tree->branchCount() << std::endl;
    std::cout << "      Pin count: " << fi_tree->pinCount() << std::endl;
  }
  
  // Step 6: Find path from driver to target D pin in Steiner tree
  int branch_count = fi_tree->branchCount();
  int pin_count = fi_tree->pinCount();
  
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
  
  int drvr_pt = fi_tree->drvrPt();
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
  
  // Calculate pre-leaf and final segment lengths
  const size_t total_segments = node_path.size() - 1;
  std::vector<int> pre_leaf_branch_indices;
  int final_branch_index = -1;
  float pre_leaf_length_dbu = 0.0;
  float final_segment_length_dbu = 0.0;
  
  for (size_t k = 0; k < total_segments; ++k) {
    int path_node1 = node_path[k];
    int path_node2 = node_path[k + 1];
    
    int current_wire_length = 0;
    int current_branch_index = -1;
    
    for (int i = 0; i < branch_count; i++) {
      odb::Point pt1, pt2;
      int steiner_pt1, steiner_pt2, wire_length;
      fi_tree->branch(i, pt1, steiner_pt1, pt2, steiner_pt2, wire_length);
      
      if ((steiner_pt1 == path_node1 && steiner_pt2 == path_node2) ||
          (steiner_pt1 == path_node2 && steiner_pt2 == path_node1)) {
        current_branch_index = i;
        current_wire_length = wire_length;
        break;
      }
    }
    
    if (k < total_segments - 1) {
      pre_leaf_length_dbu += current_wire_length;
      pre_leaf_branch_indices.push_back(current_branch_index);
    } else {
      final_segment_length_dbu = current_wire_length;
      final_branch_index = current_branch_index;
    }
  }
  
  // Step 7: Solve for max distance
  float l1 = final_segment_length_dbu / pow(10, 6);
  float on_path_R_wo_last = unit_r * pre_leaf_length_dbu / pow(10, 6);
  
  if (verbose) {
    sta::Unit* res_unit = sta_->units()->resistanceUnit();
    const char* res_suffix = res_unit->scaledSuffix();
    
    std::cout << "    Wire path segmentation:" << std::endl;
    std::cout << "      Total segments: " << total_segments << std::endl;
    std::cout << "      Pre-leaf segments: " << (total_segments - 1) << std::endl;
    std::cout << "      Pre-leaf length: " << pre_leaf_length_dbu << " DBU (" << (pre_leaf_length_dbu / pow(10, 6)) << " mm)" << std::endl;
    std::cout << "      Final segment length (l1): " << final_segment_length_dbu << " DBU (" << l1 << " mm)" << std::endl;
    std::cout << "      Pre-leaf resistance: " << res_unit->asString(on_path_R_wo_last) << res_suffix << std::endl;
  }
  
  float max_dist = solveMaxDistanceFanIn(l1, unit_r, unit_c, slack_budget, 
                                         coeff, on_path_R_wo_last, ipin_cap, verbose);
  
  // Step 8: Create feasible region box
  if (verbose) {
    std::cout << "    Maximum placement distance:" << std::endl;
    std::cout << "      Max Manhattan distance from Steiner: " << max_dist << " DBU" << std::endl;
    std::cout << "      Max distance in mm: " << (max_dist / pow(10, 6)) << " mm" << std::endl;
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
  
  // Start with inverse box (empty intersection)
  Box final_region = boost::geometry::make_inverse<Box>();
  bool has_valid_region = false;
  
  // Intersect all pin feasible regions
  for (odb::dbITerm* pin : all_pins) {
    auto it = flop.pin_feasible_regions_.find(pin);
    if (it == flop.pin_feasible_regions_.end()) {
      continue;  // Pin not found
    }
    
    const Box& pin_box = it->second;
    
    // Skip inverse boxes (no constraints) - they don't affect the intersection
    if (pin_box.min_corner().get<0>() > pin_box.max_corner().get<0>() ||
        pin_box.min_corner().get<1>() > pin_box.max_corner().get<1>()) {
      if (verbose) {
        std::cout << "    Pin " << pin->getMTerm()->getName() 
                  << " has inverse box (no constraints) - skipping" << std::endl;
      }
      continue;
    }
    
    // Valid box - compute intersection with final region
    if (verbose) {
      std::cout << "    Pin " << pin->getMTerm()->getName() 
                << " has valid region: [(" 
                << pin_box.min_corner().get<0>() << ", " << pin_box.min_corner().get<1>() 
                << ") - (" 
                << pin_box.max_corner().get<0>() << ", " << pin_box.max_corner().get<1>() 
                << ")]" << std::endl;
    }
    
    if (!has_valid_region) {
      // First valid box - initialize final_region
      final_region = pin_box;
      has_valid_region = true;
    } else {
      // Compute intersection (all constraints must be satisfied)
      Box intersection;
      boost::geometry::intersection(final_region, pin_box, intersection);
      
      // Check if intersection is empty
      if (boost::geometry::is_empty(intersection)) {
        if (verbose) {
          std::cout << "      WARNING: No intersection - constraints conflict!" << std::endl;
          std::cout << "      Setting feasible region to inverse box (no valid placement)" << std::endl;
        }
        // No valid placement exists - set to inverse box
        flop.feasible_region_ = boost::geometry::make_inverse<Box>();
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
  
  // If all pins had inverse boxes, use core area
  if (!has_valid_region) {
    if (verbose) {
      std::cout << "    No valid pin constraints - using core area" << std::endl;
    }
    
    const odb::Rect& core_area = block_->getCoreArea();
    Point core_min_orig(core_area.xMin(), core_area.yMin());
    Point core_max_orig(core_area.xMax(), core_area.yMax());
    Point core_min_tf = transformCoords(core_min_orig);
    Point core_max_tf = transformCoords(core_max_orig);
    
    int min_x = std::min(core_min_tf.get<0>(), core_max_tf.get<0>());
    int max_x = std::max(core_min_tf.get<0>(), core_max_tf.get<0>());
    int min_y = std::min(core_min_tf.get<1>(), core_max_tf.get<1>());
    int max_y = std::max(core_min_tf.get<1>(), core_max_tf.get<1>());
    
    Point box_min(min_x, min_y);
    Point box_max(max_x, max_y);
    final_region = Box(box_min, box_max);
    
    if (verbose) {
      std::cout << "    Core area transformed: [(" << min_x << ", " << min_y 
                << ") - (" << max_x << ", " << max_y << ")]" << std::endl;
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
  // where x is the new Manhattan distance from Steiner point to flop
  auto a = unit_r * unit_c;
  auto b = wo_fst_stt_cap * unit_r + coeff * unit_c;
  auto c = -std::pow(l1, 2) * unit_r * unit_c 
           - slack_budget 
           - wo_fst_stt_cap * l1 * unit_r 
           - coeff * l1 * unit_c;
  
  float D = std::pow(b, 2) - 4 * a * c;
  float max_dist = l1 * pow(10, 6);  // Default to original distance
  
  if (verbose) {
    std::cout << "    Solving quadratic equation for maximum distance (FanOut):" << std::endl;
    std::cout << "      Coefficients: a = " << a << ", b = " << b << ", c = " << c << std::endl;
    std::cout << "      Discriminant: " << D << std::endl;
  }
  
  if (D < 0) {
    if (verbose) {
      std::cout << "      No real roots - discriminant is negative" << std::endl;
      std::cout << "      Using original distance (l1): " << (l1 * pow(10, 6)) << " DBU (" 
                << l1 << " mm)" << std::endl;
    }
  } else {
    // Calculate roots using quadratic formula
    float root1 = (-b + sqrt(D)) / (2 * a) * pow(10, 6);
    float root2 = (-b - sqrt(D)) / (2 * a) * pow(10, 6);
    
    if (verbose) {
      std::cout << "      Root 1: " << root1 << " DBU (" 
                << (root1 / pow(10, 6)) << " mm)" << std::endl;
      std::cout << "      Root 2: " << root2 << " DBU (" 
                << (root2 / pow(10, 6)) << " mm)" << std::endl;
    }
    
    // Use larger positive root, or original distance if both negative
    if (root1 > 0 || root2 > 0) {
      max_dist = std::max(root1, root2);
      if (verbose) {
        std::cout << "      Selected maximum distance: " << max_dist << " DBU (" 
                  << (max_dist / pow(10, 6)) << " mm)" << std::endl;
      }
    } else {
      if (verbose) {
        std::cout << "      Both roots are negative - using original distance" << std::endl;
      }
    }
  }
  
  return max_dist;
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
  auto a = unit_r * unit_c;
  auto b = on_path_R_wo_last * unit_c + coeff * unit_c + on_path_R_wo_last * ipin_cap;
  auto c = -pow(l1, 2) * unit_r * unit_c 
           - on_path_R_wo_last * unit_c * l1 
           - coeff * l1 * unit_c 
           - ipin_cap * l1 * unit_r
           - slack_budget;
  
  float D = pow(b, 2) - 4 * a * c;
  float max_dist = l1 * pow(10, 6);  // Default to original distance
  
  if (verbose) {
    std::cout << "    Solving quadratic equation for maximum distance (FanIn):" << std::endl;
    std::cout << "      Coefficients: a = " << a << ", b = " << b << ", c = " << c << std::endl;
    std::cout << "      Discriminant: " << D << std::endl;
  }
  
  if (D < 0) {
    if (verbose) {
      std::cout << "      No real roots - discriminant is negative" << std::endl;
      std::cout << "      Using original distance (l1): " << (l1 * pow(10, 6)) << " DBU (" 
                << l1 << " mm)" << std::endl;
    }
  } else {
    // Calculate roots using quadratic formula
    float root1 = (-b + sqrt(D)) / (2 * a) * pow(10, 6);
    float root2 = (-b - sqrt(D)) / (2 * a) * pow(10, 6);
    
    if (verbose) {
      std::cout << "      Root 1: " << root1 << " DBU (" 
                << (root1 / pow(10, 6)) << " mm)" << std::endl;
      std::cout << "      Root 2: " << root2 << " DBU (" 
                << (root2 / pow(10, 6)) << " mm)" << std::endl;
    }
    
    // Use larger positive root, or original distance if both negative
    if (root1 > 0 || root2 > 0) {
      max_dist = std::max(root1, root2);
      if (verbose) {
        std::cout << "      Selected maximum distance: " << max_dist << " DBU (" 
                  << (max_dist / pow(10, 6)) << " mm)" << std::endl;
      }
    } else {
      if (verbose) {
        std::cout << "      Both roots are negative - using original distance" << std::endl;
      }
    }
  }
  
  return max_dist;
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

// -----------------------------------------------------------------------------
// [Level 3] Quadratic equation solving for max distance calculation
// -----------------------------------------------------------------------------

std::unordered_map<odb::dbITerm*, std::pair<int, sta::Slack>> 
AggloCluster::calcUsedSlacks(const FlopUnit& flop_unit) const 
{
  throw std::logic_error("calcUsedSlacks not implemented");
}

//==============================================================================
// FlopCluster Helper Functions
//==============================================================================

Box 
AggloCluster::calcMedianBox(const FlopCluster& c1,
                            const FlopCluster& c2) const

{
  std::set<std::variant<odb::dbITerm*, odb::dbBTerm*>> connected_pins
      = getConnectedPins(c1);

  std::set<std::variant<odb::dbITerm*, odb::dbBTerm*>> pins_from_c2
      = getConnectedPins(c2);

  connected_pins.insert(pins_from_c2.begin(), pins_from_c2.end());

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
          }
          else if constexpr (std::is_same_v<T, odb::dbBTerm*>) {
            odb::Rect bbox = pin->getBBox();
            if (!bbox.isInverted()) {
              x_coords.push_back(bbox.xCenter());
              y_coords.push_back(bbox.yCenter());
            }
          }
        },
        pin_variant);
  }

  size_t n = x_coords.size();
  std::sort(x_coords.begin(), x_coords.end());
  std::sort(y_coords.begin(), y_coords.end());

  int min_x = x_coords[(n - 1) / 2];
  int max_x = x_coords[n / 2];
  int min_y = y_coords[(n - 1) / 2];
  int max_y = y_coords[n / 2];

  return Box(Point(min_x, min_y), Point(max_x, max_y));
}

PlacementCandidate
AggloCluster::calcPlacementCandidate(const FlopCluster& c1, const FlopCluster& c2)
{
  // --- 0. 병합될 마스터 정보 가져오기 ---
  const MasterMask& mask = c1.master_mask_; // (c1, c2 마스크는 동일)

  // 클러스터의 현재 마스터를 가져오는 헬퍼 람다
  auto get_cluster_master = [&, this](const FlopCluster& c) -> odb::dbMaster* {
    int bits = c.flops_.size();
    if (bits == 1) {
      // 1-bit (Single cluster) -> FlopUnit에서 직접 가져옴
      int flop_id = *c.flops_.begin();
      return flop_units_[flop_id].inst_->getMaster();
    } else {
      // 2-bit 이상 (Merged cluster) -> 대표 마스터 맵에서 가져옴
      auto mask_it = representative_masters_.find(c.master_mask_);
      if (mask_it == representative_masters_.end()) return nullptr;
      auto bit_it = mask_it->second.find(bits);
      if (bit_it == mask_it->second.end()) return nullptr;
      return bit_it->second;
    }
  };
  
  odb::dbMaster* master1 = get_cluster_master(c1);
  odb::dbMaster* master2 = get_cluster_master(c2);
  
  // 병합될 새 마스터는 항상 representative_masters_ 맵에서 찾아야 함
  int new_bits = c1.flops_.size() + c2.flops_.size();
  odb::dbMaster* new_master = nullptr;
  auto mask_it = representative_masters_.find(mask);
  if (mask_it != representative_masters_.end()) {
    auto bit_it = mask_it->second.find(new_bits);
    if (bit_it != mask_it->second.end()) {
      new_master = bit_it->second;
    }
  }

  // 병합에 필요한 마스터 정보가 없는 경우 (e.g., 3-bit, 5-bit 마스터가 없음)
  if (master1 == nullptr || master2 == nullptr || new_master == nullptr) {
    return PlacementCandidate(c1.curr_pt_, 0.0, false);
  }

  // --- 1. HPWL 최적 박스 (XY 공간) ---
  const Box hpwl_box_xy = calcMedianBox(c1, c2);

  // --- 2. 타이밍 Feasible Region (UV 공간) ---
  const Box fr_box_uv = getFeasibleRegionIntersection(c1, c2);
  if (boost::geometry::is_empty(fr_box_uv)) {
    // 공통된 Feasible Region이 없으면 병합 불가
    return PlacementCandidate(c1.curr_pt_, 0.0, false);
  }

  // --- 3. HPWL 최적점 프로젝션 (project 함수가 XY 좌표 반환) ---
  Point P_proj_xy = project(hpwl_box_xy, fr_box_uv);

  // --- 4. Feasible Region 샘플링 (generateUniformSamples가 XY 좌표 반환) ---
  std::vector<Point> xy_candidates = generateUniformSamples(fr_box_uv, 16);
  xy_candidates.push_back(P_proj_xy); // 프로젝션 포인트도 후보에 추가

  // --- 5. 후보군들을 P_proj_xy (XY 좌표)에 가까운 순서로 정렬 ---
  // (정렬 기준: XY 공간에서의 유클리드 거리 제곱)
  auto dist_sq_xy = [](const Point& a, const Point& b) {
    // [수정] .get<>() 멤버 함수 사용
    int64_t dx = a.get<0>() - b.get<0>();
    int64_t dy = a.get<1>() - b.get<1>();
    return dx * dx + dy * dy;
  };

  std::sort(xy_candidates.begin(),
            xy_candidates.end(),
            [&](const Point& a, const Point& b) {
              // P_proj_xy를 기준으로 정렬
              return dist_sq_xy(a, P_proj_xy) < dist_sq_xy(b, P_proj_xy);
            });

  // --- 6. 정렬된 후보군(XY)을 순회하며 검사 ---
  for (const auto& xy_cand_int : xy_candidates) { 
    
    // 밀도 제약 검사
    if (checkPlacementDensityConstraint(c1, master1, c2, master2, xy_cand_int, new_master)) {
      
      // 밀도 통과! -> FloatPoint로 변환
      FloatPoint xy_cand_float(static_cast<float>(xy_cand_int.get<0>()),
                               static_cast<float>(xy_cand_int.get<1>()));

      // HPWL diff 계산
      double hpwl_diff = calcHPWLDiff(c1, c2, xy_cand_float);
      double gain = -hpwl_diff; // HPWL 감소량이 gain

      // 배치 가능한 첫 번째 후보를 즉시 반환
      return PlacementCandidate(xy_cand_float, gain, true);
    }
  }

  // --- 7. 모든 후보가 밀도 제약에 실패한 경우 ---
  return PlacementCandidate(c1.curr_pt_, 0.0, false);
}

double 
AggloCluster::calcHPWLDiff(const FlopCluster& c1,
                            const FlopCluster& c2,
                            const FloatPoint& new_pos) const
  {
  double old_hpwl = 0.0;
  double new_hpwl = 0.0;

  auto pins1 = getConnectedPins(c1);
  for (const auto& pin_variant : pins1) {
    FloatPoint pin_loc = getPinCoordinate(pin_variant);
    old_hpwl += calcManhattanDistance(pin_loc, c1.curr_pt_); 
  }

  auto pins2 = getConnectedPins(c2);
  for (const auto& pin_variant : pins2) {
    FloatPoint pin_loc = getPinCoordinate(pin_variant);
    old_hpwl += calcManhattanDistance(pin_loc, c2.curr_pt_); 
  }

  std::set<std::variant<odb::dbITerm*, odb::dbBTerm*>> all_pins = pins1;
  all_pins.insert(pins2.begin(), pins2.end()); 

  for (const auto& pin_variant : all_pins) {
    FloatPoint pin_loc = getPinCoordinate(pin_variant);
    new_hpwl += calcManhattanDistance(pin_loc, new_pos); 
  }

  return new_hpwl - old_hpwl; 
}

bool 
AggloCluster::checkCompatibility(const FlopCluster& c1,
                                      const FlopCluster& c2) const
{
  if (!(c1.master_mask_ == c2.master_mask_
        && c1.inst_mask_ == c2.inst_mask_)) {
    return false;
  }

  size_t total_flops = c1.flops_.size() + c2.flops_.size();
  auto it = compatible_masters_.find(c1.master_mask_);

  if (it != compatible_masters_.end()) {
    const auto& bits_to_masters = it->second;
    if (bits_to_masters.count(total_flops) > 0) {
      return true;
    }
  }

  return false;
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
      std::cout << "[checkPlacementDensityConstraint] Null master passed. Skipping density check." << std::endl;
    }
    return false;
  }

  Point c1_pos(std::lround(c1.curr_pt_.x), std::lround(c1.curr_pt_.y));
  Point c2_pos(std::lround(c2.curr_pt_.x), std::lround(c2.curr_pt_.y));

  virtual_bin_grid_.removeInst(master1, c1_pos);
  virtual_bin_grid_.removeInst(master2, c2_pos);
  virtual_bin_grid_.addInst(new_master, new_pos);

  bool has_overflow = virtual_bin_grid_.checkOverflow();

  if (has_overflow) {
    virtual_bin_grid_.removeInst(new_master, new_pos);
    virtual_bin_grid_.addInst(master2, c2_pos);
    virtual_bin_grid_.addInst(master1, c1_pos);

    return false;
  }
  return true; 
}

std::set<int>
AggloCluster::distributeSlack(const FlopCluster& cluster)
{
  throw std::logic_error("distributeSlack not implemented");
  // std::set<int> affected_flop_units;

  // for (int flop_idx : cluster.flops_) {
  //   const FlopUnit& flop = flop_units_[flop_idx];

  //   // Calculate used slacks for each pin of the flop 
  //   std::unordered_map<odb::dbITerm*, std::pair<int, sta::Slack>> used_slacks = calcUsedSlacks(flop);

  //   for (const auto& [pin, used_slack_pair] : used_slacks) {
  //     const auto& [path_idx, used_slack] = used_slack_pair;
  //     const TimingPath& path = timing_paths_[path_idx];

  //     // Find the flop on the other side of the timing path
  //     int other_flop_idx = (path.start_flop_idx_ == flop_idx)
  //                              ? path.end_flop_idx_
  //                              : path.start_flop_idx_;

  //     // If the other flop is also in the same cluster, no need to distribute slack
  //     if (flop_units_[other_flop_idx].cluster_idx_ == cluster.id_) {
  //       continue;
  //     }

  //     FlopUnit& other_flop = flop_units_[other_flop_idx];
  //     affected_flop_units.insert(other_flop_idx);

  //     // Find the connected pin on the other flop for this path
  //     odb::dbITerm* other_pin = getConnectedPinOnPath(path.path_, other_flop);
  //     if (other_pin == nullptr
  //         || other_flop.pin_budgets_.find(other_pin)
  //                == other_flop.pin_budgets_.end()) {
  //       continue;
  //     }

  //     auto& other_budgets = other_flop.pin_budgets_.at(other_pin);
  //     for (auto& budget_pair : other_budgets) {
  //       if (budget_pair.first == path_idx) {
  //         // Find the original budget allocated to this path on the current flop
  //         sta::Slack original_budget = 0.0; 
  //         const auto& budgets = flop.pin_budgets_.at(pin);
  //         for (const auto& orig_budget_pair : budgets) {
  //           if (orig_budget_pair.first == path_idx) {
  //             original_budget = orig_budget_pair.second;
  //             break;
  //           }
  //         }

  //         // Update the budget on the other flop by adding the slack delta
  //         sta::Slack slack_delta = original_budget - used_slack;
  //         budget_pair.second += slack_delta;
  //         break;
  //       }
  //     }

  //     // Sort the budgets for the updated pin
  //     std::sort(other_budgets.begin(),
  //               other_budgets.end(),
  //               [](const auto& a, const auto& b) { return a.second < b.second; });
  //   }
  // }
  // return affected_flop_units;
}

std::set<std::variant<odb::dbITerm*, odb::dbBTerm*>>
AggloCluster::getConnectedPins(const FlopCluster& cluster) const
{
  std::set<std::variant<odb::dbITerm*, odb::dbBTerm*>> connected_pins;
  
  std::set<odb::dbInst*> cluster_insts;
  for (int flop_idx : cluster.flops_) {
    cluster_insts.insert(flop_units_[flop_idx].inst_);
  }
  
  for (int flop_idx : cluster.flops_) {
    const FlopUnit& flop = flop_units_[flop_idx];
    for (odb::dbITerm* iterm : flop.inst_->getITerms()) {
      odb::dbNet* net = iterm->getNet();
      if (net == nullptr
          || (net->getSigType() != odb::dbSigType::SIGNAL
              && net->getSigType() != odb::dbSigType::CLOCK)) {
        continue;
      }
      for (odb::dbITerm* net_iterm : net->getITerms()) {
        if (cluster_insts.find(net_iterm->getInst()) == cluster_insts.end()) {
          connected_pins.insert(net_iterm);
        }
      }
      for (odb::dbBTerm* net_bterm : net->getBTerms()) {
        connected_pins.insert(net_bterm);
      }
    }
  }
  return connected_pins;
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

std::vector<FlopClusterEntry>
AggloCluster::getIntersectedCluster(const FlopCluster& cluster) const
{
  std::vector<FlopClusterEntry> intersecting_clusters;

  if (boost::geometry::is_empty(cluster.feasible_region_)) {
    return intersecting_clusters;
  }

  feasible_regions_.query(bgi::intersects(cluster.feasible_region_),
                          std::back_inserter(intersecting_clusters));

  return intersecting_clusters;
}

bool
AggloCluster::isFurtherMergeable(const FlopCluster& cluster) const
{
  auto it = compatible_masters_.find(cluster.master_mask_);
  if (it == compatible_masters_.end()) {
    return false;
  }

  const auto& bits_to_masters = it->second;
  if (bits_to_masters.empty()) {
    return false;
  }

  int max_bits = bits_to_masters.rbegin()->first;

  return cluster.flops_.size() < max_bits;
}  

int 
AggloCluster::mergeClusters(const Edge& edge)
{
  const FlopCluster& c1 = flop_clusters_[edge.n1];
  const FlopCluster& c2 = flop_clusters_[edge.n2];

  // Marking old clusters as invalid
  flop_cluster_is_valid_[edge.n1] = false;
  flop_cluster_is_valid_[edge.n2] = false;

  // Remove edges associated with old clusters
  removeEdges(c1);
  removeEdges(c2);

  // Remove feasible regions of old clusters from R-Tree
  feasible_regions_.remove(std::make_pair(c1.feasible_region_, edge.n1));
  feasible_regions_.remove(std::make_pair(c2.feasible_region_, edge.n2));

  // Create new cluster by merging c1 and c2
  int n_new_idx = flop_clusters_.size();
  flop_clusters_.emplace_back(n_new_idx, c1, c2, edge);
  flop_cluster_is_valid_.push_back(true);
  flop_cluster_no_further_merge_.push_back(false);
  FlopCluster& n_new = flop_clusters_.back();

  // Update FlopUnit's info
  for (int flop_idx : n_new.flops_) {
    flop_units_[flop_idx].cluster_idx_ = n_new_idx;
    flop_units_[flop_idx].curr_pt_ = n_new.curr_pt_;
  }

  // Insert feasible region of new cluster into R-Tree
  if (!boost::geometry::is_empty(n_new.feasible_region_)) {
    feasible_regions_.insert(std::make_pair(n_new.feasible_region_, n_new_idx));
  }

  return n_new_idx;
}

void 
AggloCluster::removeEdges(const FlopCluster& cluster)
{
  const int cluster_idx = cluster.id_;

  if (adj_list_.find(cluster_idx) == adj_list_.end()) {
    return;
  }

  auto& edges_to_remove = adj_list_.at(cluster_idx);
  
  for (const Edge& edge : edges_to_remove) {
    edge_pq_.erase(edge);

    int neighbor_idx = (edge.n1 == cluster_idx) ? edge.n2 : edge.n1;

    if (adj_list_.count(neighbor_idx)) {
      adj_list_.at(neighbor_idx).erase(edge);
    }
  }

  adj_list_.erase(cluster_idx);
}

void 
AggloCluster::updateEdges(const FlopCluster& cluster)
{
  const int i = cluster.id_;

  if (!flop_cluster_is_valid_[i] || flop_cluster_no_further_merge_[i]) {
    return;
  }

  if (boost::geometry::is_empty(cluster.feasible_region_)) {
    return;
  }

  std::vector<FlopClusterEntry> intersecting_entries = getIntersectedCluster(cluster);

  for (const auto& entry : intersecting_entries) {
    const int j = entry.second;

    if (i >= j) {
      continue;
    }

    if (!flop_cluster_is_valid_[j] || flop_cluster_no_further_merge_[j]) {
      continue;
    }

    const FlopCluster& adj_cluster = flop_clusters_[j];

    if (!checkCompatibility(cluster, adj_cluster)) {
      continue;
    }

    PlacementCandidate pc = calcPlacementCandidate(cluster, adj_cluster);

    if (!pc.is_valid_) {
      continue;
    }

    Edge new_edge(i, j, pc.merge_gain_, pc.placement_pos_);

    edge_pq_.insert(new_edge);
    adj_list_[i].insert(new_edge);
    adj_list_[j].insert(new_edge);
  }
}

void 
AggloCluster::updateFeasibleRegion(FlopCluster& cluster)
{
  const Box old_fr = cluster.feasible_region_;

  if (cluster.flops_.empty()) {
    cluster.feasible_region_ = Box();
  } else {
    auto it = cluster.flops_.begin();
    cluster.feasible_region_ = flop_units_[*it].feasible_region_;
    for (++it; it != cluster.flops_.end(); ++it) {
      const int flop_idx = *it;
      const Box& flop_fr = flop_units_[flop_idx].feasible_region_;
      boost::geometry::intersection(
          cluster.feasible_region_, flop_fr, cluster.feasible_region_);
      if (boost::geometry::is_empty(cluster.feasible_region_)) {
        break;
      }
    }
  }

  if (!boost::geometry::is_empty(old_fr)) {
    feasible_regions_.remove(std::make_pair(old_fr, cluster.id_));
  }
  if (!boost::geometry::is_empty(cluster.feasible_region_)) {
    feasible_regions_.insert(std::make_pair(cluster.feasible_region_, cluster.id_));
  }
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

Point 
AggloCluster::project(const Box& hpwl_box, const Box& feasible_box) const
{
    const Point& min_xy = hpwl_box.min_corner();
    const Point& max_xy = hpwl_box.max_corner();

    int target_x = (min_xy.get<0>() + max_xy.get<0>()) / 2;
    int target_y = (min_xy.get<1>() + max_xy.get<1>()) / 2;
    Point target_xy(target_x, target_y);

    Point target_uv = transformCoords(target_xy);

    const Point& min_uv = feasible_box.min_corner();
    const Point& max_uv = feasible_box.max_corner();

    int u_final = std::clamp(target_uv.get<0>(),  // u_target
                             min_uv.get<0>(),     // u_min
                             max_uv.get<0>());    // u_max

    int v_final = std::clamp(target_uv.get<1>(),  // v_target
                             min_uv.get<1>(),     // v_min
                             max_uv.get<1>());    // v_max

    Point final_uv(u_final, v_final);

    Point final_xy = inverseTransformCoords(final_uv);

    return final_xy;
}

//==============================================================================
// Other Utility Helper Functions
//==============================================================================

FloatPoint 
AggloCluster::getPinCoordinate(const std::variant<odb::dbITerm*, odb::dbBTerm*>& pin_variant) const
{
  int x = 0, y = 0;

  if (std::holds_alternative<odb::dbITerm*>(pin_variant)) {
    odb::dbITerm* iterm = std::get<odb::dbITerm*>(pin_variant);
    
    if (iterm->getAvgXY(&x, &y)) {
      return FloatPoint(static_cast<float>(x), static_cast<float>(y));
    } else {
      // Critical error - cannot proceed without pin coordinates
      std::cerr << "[FATAL] Failed to get average XY for ITerm: " << iterm->getName() << std::endl;
      std::exit(1);
    }
  } else {
    odb::dbBTerm* bterm = std::get<odb::dbBTerm*>(pin_variant);
    
    odb::Rect bbox = bterm->getBBox();
    if (!bbox.isInverted()) {
      return FloatPoint(static_cast<float>(bbox.xCenter()),
                        static_cast<float>(bbox.yCenter()));
    } else {
      // Critical error - cannot proceed without pin coordinates
      std::cerr << "[FATAL] Failed to get valid BBox for BTerm (is inverted): " 
                << bterm->getName() << std::endl;
      std::exit(1);
    }
  }
}

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

double // 이거 수정해야하네 
AggloCluster::calcManhattanDistance(const FloatPoint& p1, const FloatPoint& p2) const
{
  double dx = std::abs(p1.x - p2.x);
  double dy = std::abs(p1.y - p2.y);
  return dx + dy;
}

std::vector<gpl::Point> 
gpl::AggloCluster::generateUniformSamples(
    const Box& box, // 입력 Box는 UV 공간의 Feasible Region
    int p) const
{
  std::vector<Point> samples;
  if (p <= 0) {
    return samples;
  }

  int u_min = box.min_corner().get<0>();
  int u_max = box.max_corner().get<0>();
  int v_min = box.min_corner().get<1>();
  int v_max = box.max_corner().get<1>();

  // [수정] p=1인 경우 (Box의 중심점) 예외 처리
  if (p == 1) {
    int u = (u_min + u_max) / 2;
    int v = (v_min + v_max) / 2;
    Point center_p_uv(u, v);
    Point center_p_xy = inverseTransformCoords(center_p_uv);
    samples.push_back(center_p_xy);
    return samples;
  }

  // --- p > 1 인 경우 (경계 포함 그리드) ---
  
  // p=2~4일 때 2, p=5~9일 때 3, ... (p > 1이므로 num_steps >= 2 보장)
  int num_steps = static_cast<int>(std::ceil(std::sqrt(p)));
  
  // 경계를 포함하기 위해 (num_steps - 1)로 나눔
  double u_step = (u_max == u_min) ? 0.0 : static_cast<double>(u_max - u_min) / (num_steps - 1.0);
  double v_step = (v_max == v_min) ? 0.0 : static_cast<double>(v_max - v_min) / (num_steps - 1.0);

  // 루프를 0부터 num_steps-1까지
  for (int i = 0; i < num_steps; ++i) {
    for (int j = 0; j < num_steps; ++j) {
      
      // std::round를 사용하고, 마지막 스텝은 max값으로 명시적 보장
      int u = (i == num_steps - 1) ? u_max : (u_min + static_cast<int>(std::round(i * u_step)));
      int v = (j == num_steps - 1) ? v_max : (v_min + static_cast<int>(std::round(j * v_step)));
      
      Point sample_p_uv(u, v); 

      Point sample_p_xy = inverseTransformCoords(sample_p_uv);
      samples.push_back(sample_p_xy);
    }
  }
  
  return samples; // XY 좌표 리스트 반환
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
