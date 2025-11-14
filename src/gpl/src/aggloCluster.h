#pragma once

//==============================================================================
// File Overview
//==============================================================================
// 
// This header defines the Agglomerative Clustering engine for converting
// single-bit flip-flops into multi-bit flip-flops (MBFFs).
//
// Key Responsibilities:
// 1. Data structures for representing flops, clusters, and timing information
// 2. Main clustering algorithm orchestration
// 3. Helper functions for mask creation, feasible region calculation, etc.
//
// Main Components:
// 1. Mask info: MasterMask, InstMask (flop characteristics for grouping)
// 2. Timing info: TimingPath (timing path and slack between flops)
// 3. Spatial info: VirtualBin, VirtualBinGrid (density constraint checking)
// 4. Graph info: Edge, EdgeHash (compatibility graph edges)
// 5. Helper structures: NetBundle, PortBundle, MasterPortAssignment, PlacementCandidate
// 6. Core structures: FlopUnit, FlopCluster (individual and merged flops)
// 7. Main engine: AggloCluster (orchestrates the entire clustering flow)
//==============================================================================

// Standard library
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <variant>
#include <vector>

// Boost geometry
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

// OpenROAD
#include "odb/db.h"
#include "sta/Delay.hh"
#include "point.h"

//==============================================================================
// Namespace Aliases
//==============================================================================

namespace bgi = boost::geometry::index;
namespace bg = boost::geometry::model;

//==============================================================================
// Forward Declarations
//==============================================================================

namespace utl {
class Logger;
}

namespace rsz {
class Resizer;
}

namespace est {
class EstimateParasitics;
class SteinerTree;
}

namespace sta {
class dbNetwork;
class dbSta;
class Corner;
class FuncExpr;
class LibertyCell;
class LibertyPort;
class Path;
class Pin;
class TableAxis;
class GateTableModel;
class Pvt;
class Graph;
class MinMax;
}

//==============================================================================
// GPL Namespace
//==============================================================================

namespace gpl {

//==============================================================================
// Type Aliases
//==============================================================================

using Point = bg::point<int, 2, boost::geometry::cs::cartesian>;
using Box = bg::box<Point>;
using FlopClusterEntry = std::pair<Box, int>;
using FlopClusterRTree = bgi::rtree<FlopClusterEntry, bgi::rstar<16, 4, 4>>;

//==============================================================================
// Enumerations
//==============================================================================

enum FlopPort
{
  d, q, qn, clear, preset, si, se, vdd, vss, func, ifunc, unknown
};

//==============================================================================
// 1. Mask Structures
//==============================================================================

/**
 * @brief Master cell characteristics for flop grouping
 * 
 * Encapsulates pin characteristics of the master cell (function ID, clock edge,
 * CLR, PRE, Q, QN, Scan). Flops with identical masks can be grouped together.
 */
struct MasterMask
{
  // Member variables
  int func_id_{-1};
  bool has_positive_clock_edge_{false};
  bool has_clear_{false};
  bool has_preset_{false};
  bool has_q_pin_{false};
  bool has_qn_pin_{false};
  bool has_scan_{false};

  // Constructors
  MasterMask() = default;
  
  MasterMask(int func_id,
             bool has_positive_clock_edge,
             bool has_clear,
             bool has_preset,
             bool has_q_pin,
             bool has_qn_pin,
             bool has_scan)
      : func_id_(func_id),
        has_positive_clock_edge_(has_positive_clock_edge),
        has_clear_(has_clear),
        has_preset_(has_preset),
        has_q_pin_(has_q_pin),
        has_qn_pin_(has_qn_pin),
        has_scan_(has_scan)
  {
  }

  // Comparison operators
  auto _tie() const
  {
    return std::tie(func_id_,
                    has_positive_clock_edge_,
                    has_clear_,
                    has_preset_,
                    has_q_pin_,
                    has_qn_pin_,
                    has_scan_);
  }

  bool operator<(const MasterMask& rhs) const { return _tie() < rhs._tie(); }
  bool operator==(const MasterMask& rhs) const { return _tie() == rhs._tie(); }

  // Utilities
  std::string to_string() const;
};

/**
 * @brief Instance net connections for flop grouping
 * 
 * Stores actual net connections for a flip-flop instance
 * (CLK, CLR, PRE, SE, SI nets). Flops with identical connections
 * can potentially be merged.
 */
struct InstMask
{
  // Member variables
  odb::dbNet* clock_net_{nullptr};
  odb::dbNet* clear_net_{nullptr};
  odb::dbNet* preset_net_{nullptr};
  odb::dbNet* scan_enable_net_{nullptr};
  odb::dbNet* scan_in_net_{nullptr};

  // Constructors
  InstMask() = default;
  
  InstMask(odb::dbNet* clock_net,
           odb::dbNet* clear_net,
           odb::dbNet* preset_net,
           odb::dbNet* scan_enable_net,
           odb::dbNet* scan_in_net)
      : clock_net_(clock_net),
        clear_net_(clear_net),
        preset_net_(preset_net),
        scan_enable_net_(scan_enable_net),
        scan_in_net_(scan_in_net)
  {
  }

  // Comparison operators
  auto _tie() const
  {
    return std::tie(clock_net_,
                    clear_net_,
                    preset_net_,
                    scan_enable_net_,
                    scan_in_net_);
  }

  bool operator<(const InstMask& rhs) const { return _tie() < rhs._tie(); }
  bool operator==(const InstMask& rhs) const { return _tie() == rhs._tie(); }

  // Utilities
  std::string to_string() const;
};

//==============================================================================
// 2. Timing Information
//==============================================================================

/**
 * @brief Timing path between two flip-flops
 * 
 * Encapsulates timing path and slack information for timing-driven clustering.
 */
struct TimingPath
{
  // Member variables
  sta::Path* path_{nullptr};
  sta::Slack slack_{sta::INF};
  int start_flop_idx_{-1};
  int end_flop_idx_{-1};
  odb::dbITerm* start_pin_{nullptr};  // Output pin of start flop (Q or QN)
  odb::dbITerm* end_pin_{nullptr};    // Input pin of end flop (D)

  // Constructor
  TimingPath(sta::Path* path, 
             sta::Slack slack, 
             int start_idx, 
             int end_idx,
             odb::dbITerm* start_pin = nullptr,
             odb::dbITerm* end_pin = nullptr)
      : path_(path),
        slack_(slack),
        start_flop_idx_(start_idx),
        end_flop_idx_(end_idx),
        start_pin_(start_pin),
        end_pin_(end_pin)
  {
  }

  // Comparison operator
  auto _tie() const
  {
    return std::tie(slack_, start_flop_idx_, end_flop_idx_, start_pin_, end_pin_, path_);
  }

  bool operator<(const TimingPath& other) const { return _tie() < other._tie(); }
};

//==============================================================================
// 3. Spatial Information (Density Constraint Grid)
//==============================================================================

/**
 * @brief Virtual bin for density constraint checking
 * 
 * Represents a single bin in the virtual grid for tracking placement density.
 */
class VirtualBin
{
public:
  // Constructors
  VirtualBin() = default;
  
  VirtualBin(int lx, int ly, int ux, int uy, float target_density)
      : lx_(lx), ly_(ly), ux_(ux), uy_(uy), target_density_(target_density)
  {
    bin_area_ = static_cast<int64_t>(ux_ - lx_) * static_cast<int64_t>(uy_ - ly_);
  }

  // Boundary accessors
  int lx() const { return lx_; }
  int ly() const { return ly_; }
  int ux() const { return ux_; }
  int uy() const { return uy_; }

  // Area query methods
  int64_t getInstPlacedArea() const { return inst_placed_area_; }
  int64_t getMacroPlacedArea() const { return macro_placed_area_; }
  int64_t getNonPlaceArea() const { return non_place_area_; }
  int64_t getOverflowArea() const;
  int64_t getBinArea() const { return bin_area_; }
  float getTargetDensity() const { return target_density_; }

  // Area update methods
  void addInstPlacedArea(int64_t area) { inst_placed_area_ += area; }
  void addMacroPlacedArea(int64_t area) { macro_placed_area_ += area; }
  void addNonPlaceArea(int64_t area) { non_place_area_ += area; }
  void subInstPlacedArea(int64_t area) { inst_placed_area_ -= area; }
  void subMacroPlacedArea(int64_t area) { macro_placed_area_ -= area; }
  void subNonPlaceArea(int64_t area) { non_place_area_ -= area; }

private:
  // Bin boundaries
  int lx_{0}, ly_{0}, ux_{0}, uy_{0};
  
  // Areas
  int64_t bin_area_{0};
  int64_t inst_placed_area_{0};
  int64_t macro_placed_area_{0};
  int64_t non_place_area_{0};
  
  // Configuration
  float target_density_{0.0f};
};

/**
 * @brief Virtual bin grid for density management
 * 
 * Manages a uniform grid of bins covering the chip area for density checking.
 */
class VirtualBinGrid
{
public:
  // Constructors
  VirtualBinGrid() = default;
  
  VirtualBinGrid(int lx, int ly, int ux, int uy,
                 int bin_cnt_x, int bin_cnt_y,
                 double bin_size_x, double bin_size_y,
                 float target_density, float target_overflow);

  // Query methods
  bool checkOverflow();
  bool wouldOverflow(odb::dbMaster* master1, const Point& pos1,
                     odb::dbMaster* master2, const Point& pos2,
                     odb::dbMaster* new_master, const Point& new_pos) const;
  std::vector<VirtualBin>& getBins() { return bins_; }
  std::pair<int, int> getMinMaxIdxX(const odb::Rect& box) const;
  std::pair<int, int> getMinMaxIdxY(const odb::Rect& box) const;
  int getBinCntX() const { return bin_cnt_x_; }
  int getBinCntY() const { return bin_cnt_y_; }
  double getBinSizeX() const { return bin_size_x_; }
  double getBinSizeY() const { return bin_size_y_; }

  // Update methods
  void addInst(odb::dbMaster* master, const Point& center_pos);
  void removeInst(odb::dbMaster* master, const Point& center_pos);
  void applyMerge(odb::dbMaster* master1, const Point& pos1,
                  odb::dbMaster* master2, const Point& pos2,
                  odb::dbMaster* new_master, const Point& new_pos);

private:
  void accumulateInstPlacement(odb::dbMaster* master, const Point& center_pos,
                               std::vector<int64_t>& area_deltas, bool is_add) const;
  void updateInstPlacement(odb::dbMaster* master, const Point& center_pos, bool is_add);

  // Grid boundaries
  int lx_{0}, ly_{0}, ux_{0}, uy_{0};
  
  // Grid dimensions
  int bin_cnt_x_{0}, bin_cnt_y_{0};
  double bin_size_x_{0.0}, bin_size_y_{0.0};
  
  // Configuration
  float target_density_{0.0f};
  float target_overflow_{0.0f};
  
  // Bins
  std::vector<VirtualBin> bins_;
};

//==============================================================================
// 4. Graph Information (Compatibility Graph)
//==============================================================================

/**
 * @brief Edge between mergeable clusters
 * 
 * Represents a potential merge between two compatible clusters.
 */
struct Edge
{
  // Member variables
  int n1, n2;         // Cluster indices (n1 < n2)
  double weight;      // Merge cost
  FloatPoint pos;     // Proposed placement after merge

  // Constructor
  Edge(int node1, int node2, double w, FloatPoint p)
      : weight(w), pos(p)
  {
    n1 = std::min(node1, node2);
    n2 = std::max(node1, node2);
  }

  // Comparison operators
  auto _tie() const { return std::tie(weight, n1, n2, pos.x, pos.y); }
  bool operator<(const Edge& other) const { return _tie() < other._tie(); }
  bool operator==(const Edge& other) const { return _tie() == other._tie(); }
};

/**
 * @brief Hash functor for Edge
 */
struct EdgeHash
{
  std::size_t operator()(const Edge& e) const
  {
    return std::hash<int>{}(e.n1) 
         ^ (std::hash<int>{}(e.n2) << 1)
         ^ (std::hash<double>{}(e.weight) << 2) 
         ^ (FloatPoint::Hash{}(e.pos) << 3);
  }
};

//==============================================================================
// 5. Helper Structures (MBFF Implementation Support)
//==============================================================================

/**
 * @brief Net bundle for a single 1-bit flop
 * 
 * Stores D, Q, QN net connections and their corresponding pins
 * for a single flop instance in a cluster.
 */
struct NetBundle
{
  int flop_id_{-1};                    // Flop unit ID
  odb::dbInst* flop_inst_{nullptr};    // Flop instance
  odb::dbNet* d_net_{nullptr};         // D net
  odb::dbNet* q_net_{nullptr};         // Q net
  odb::dbNet* qn_net_{nullptr};        // QN net
  odb::dbITerm* d_iterm_{nullptr};     // D pin
  odb::dbITerm* q_iterm_{nullptr};     // Q pin
  odb::dbITerm* qn_iterm_{nullptr};    // QN pin
};

/**
 * @brief Port (D, Q, QN) information and local coordinates of MBFF master
 * 
 * Stores each port bundle (D0/Q0/QN0, D1/Q1/QN1, etc.) of the MBFF master cell
 * and their local coordinates within the master.
 */
struct PortBundle
{
  int bundle_idx_{-1};              // Port bundle index (0, 1, 2, ...)
  odb::dbMaster* master_{nullptr};  // Associated master
  odb::dbMTerm* d_mterm_{nullptr};  // D port master terminal
  odb::dbMTerm* q_mterm_{nullptr};  // Q port master terminal
  odb::dbMTerm* qn_mterm_{nullptr}; // QN port master terminal
  Point d_local_pos_;               // D port local coordinate
  Point q_local_pos_;               // Q port local coordinate
  Point qn_local_pos_;              // QN port local coordinate
};

/**
 * @brief Optimal master and net-pin assignment 
 * 
 * Stores the result of optimally assigning multiple 1-bit flops in a cluster
 * to MBFF master ports. Includes optimal master selection and cost.
 */
struct MasterPortAssignment
{
  odb::dbMaster* best_master_{nullptr};        // Selected MBFF master
  std::vector<int> net_to_port_assignment_;    // net_to_port_assignment_[i] = j
                                               // means i-th NetBundle assigned to j-th PortBundle
  double min_cost_{std::numeric_limits<double>::max()};  // Total HPWL cost
};

/**
 * @brief Possible placement location and gain for cluster merge
 * 
 * Stores the new placement location after merging two clusters and the resulting
 * gain (cost reduction).
 */
struct PlacementCandidate
{
  FloatPoint placement_pos_;    // New placement position after merge
  double merge_gain_;           // Gain (cost reduction) from this merge
  bool is_valid_;               // Whether this placement candidate is feasible

  PlacementCandidate(FloatPoint pos, double gain, bool valid)
      : placement_pos_(pos), merge_gain_(gain), is_valid_(valid)
  {
  }
};

//==============================================================================
// 6. Core Structures (Flop Unit and Cluster)
//==============================================================================

/**
 * @brief Single 1-bit flip-flop unit
 * 
 * Encapsulates all information of an individual flop instance.
 * Includes master/instance masks, placement info, timing constraints, feasible regions, etc.
 */
struct FlopUnit
{
  // Identity
  int id_{-1};                            // Flop unit ID
  odb::dbInst* inst_{nullptr};            // Flop instance pointer

  // Cell Characteristics
  MasterMask master_mask_;                // Master cell characteristics
  InstMask inst_mask_;                    // Net connections (CLK, CLR, PRE, SE, SI)

  // Placement (DBU)
  FloatPoint orig_pt_;                    // Original placement position
  FloatPoint curr_pt_;                    // Current placement position

  // Timing Constraints
  std::unordered_map<odb::dbITerm*, std::vector<std::pair<int, sta::Slack>>> pin_budgets_;
                                          // Key: pin of this flop (Q/QN for fan-out, D for fan-in)
                                          // Value: vector of (timing_path_idx, slack_budget)
                                          //        sorted in ascending order of slack_budget

  // Feasible Region (45-degree transformed coordinates)
  std::unordered_map<odb::dbITerm*, Box> pin_feasible_regions_;
                                          // Key: pin (ITerm)
                                          // Value: feasible region for this pin
  Box feasible_region_;                   // Intersection of all pin regions

  // Cluster Assignment
  int cluster_idx_{-1};                   // Index in flop_clusters_

  FlopUnit(int id,
           odb::dbInst* inst,
           const odb::Rect& bbox,
           const MasterMask& master_mask,
           const InstMask& inst_mask)
      : id_(id),
        inst_(inst),
        master_mask_(master_mask),
        inst_mask_(inst_mask),
        cluster_idx_(-1)
  {
    orig_pt_.x = curr_pt_.x = static_cast<float>(bbox.xCenter());
    orig_pt_.y = curr_pt_.y = static_cast<float>(bbox.yCenter());
  }

  std::string to_string() const;
};

/**
 * @brief Merged cluster of multiple flip-flops
 * 
 * A cluster of 1-bit flops suitable for merging to implement as a single MBFF.
 * All flops in a cluster have identical masks, and feasible region is the
 * intersection of all flop regions.
 */
struct FlopCluster
{
  // Identity
  int id_{-1};                            // Cluster ID

  // Members
  std::set<int> flops_;                   // FlopUnit IDs in this cluster

  // Cell Characteristics (all flops have same masks)
  MasterMask master_mask_;                // Master cell characteristics
  InstMask inst_mask_;                    // Net connections

  // Placement (DBU)
  FloatPoint curr_pt_;                    // Current placement position

  // Timing Constraints
  Box feasible_region_;                   // Intersection of all member flops' regions

  // Constructor from a single FlopUnit
  FlopCluster(int id, const FlopUnit& flop_unit)
      : id_(id),
        master_mask_(flop_unit.master_mask_),
        inst_mask_(flop_unit.inst_mask_),
        curr_pt_(flop_unit.curr_pt_),
        feasible_region_(flop_unit.feasible_region_)
  {
    flops_.insert(flop_unit.id_);
  }

  // Constructor from merging two FlopClusters
  FlopCluster(int id,
              const FlopCluster& c1,
              const FlopCluster& c2,
              const Edge& edge)
      : id_(id),
        master_mask_(c1.master_mask_),    // c1 and c2 have same masks
        inst_mask_(c1.inst_mask_),        // c1 and c2 have same masks
        curr_pt_(edge.pos)
  {
    // Compute feasible region as intersection of c1 and c2 feasible regions
    boost::geometry::intersection(c1.feasible_region_,
                                  c2.feasible_region_,
                                  feasible_region_);

    // Combine flops from both clusters
    flops_.insert(c1.flops_.begin(), c1.flops_.end());
    flops_.insert(c2.flops_.begin(), c2.flops_.end());
  }
};

//==============================================================================
// 7. Main Class: Agglomerative Clustering Engine
//==============================================================================

/**
 * @brief Hierarchical flip-flop clustering engine for MBFF synthesis
 * 
 * Converts single-bit flip-flops into multi-bit flip-flops (MBFFs) through
 * agglomerative clustering, optimizing for timing, density, and power.
 *
 * Algorithm Overview:
 * 1. Initialize density-aware virtual bin grid
 * 2. Read compatible MBFF masters from Liberty
 * 3. Extract flip-flop instances and create flop units
 * 4. Group compatible flops by masks
 * 5. Extract timing paths from STA
 * 6. Analyze timing and compute slack budgets
 * 7. Calculate feasible placement regions
 * 8. Initialize single-flop clusters
 * 9. Build compatibility graph with merge candidates
 * 10. Greedily merge clusters by cost
 * 11. Implement final clusters as MBFF instances
 */
class AggloCluster
{
public:
  // Constructor & Destructor
  AggloCluster(odb::dbDatabase* db,
               sta::dbSta* sta,
               utl::Logger* log,
               rsz::Resizer* resizer,
               float target_density,
               float target_overflow,
               float region_scale_factor,
               int num_paths_per_endpoint,
               int threads,
               int num_samples,
               bool verbose);
  ~AggloCluster();

  // Main entry point
  void doAggloCluster();

private:
  //==========================================================================
  // External References
  //==========================================================================
  
  odb::dbDatabase* db_;
  odb::dbBlock* block_;
  sta::dbSta* sta_;
  sta::dbNetwork* network_;
  sta::Corner* corner_;
  rsz::Resizer* resizer_;
  utl::Logger* log_;

  //==========================================================================
  // Configuration
  //==========================================================================
  
  // Parallelization
  int threads_;
  
  // Timing analysis
  int num_paths_per_endpoint_;
  
  // Density constraints
  float target_density_;
  float target_overflow_;
  
  // Feasible region constraint
  float region_scale_factor_;
  
  // Debug & output
  int num_samples_;
  bool verbose_;

  //==========================================================================
  // Data Structures
  //==========================================================================
  
  // Flop & cluster data
  std::vector<FlopUnit> flop_units_;
  std::vector<FlopCluster> flop_clusters_;
  std::vector<bool> flop_cluster_is_valid_;
  std::vector<bool> flop_cluster_no_further_merge_;
  
  // Timing information
  std::vector<TimingPath> timing_paths_;
  
  // Master & compatibility mappings
  std::map<MasterMask, std::map<int, std::vector<odb::dbMaster*>>> compatible_masters_;
  std::map<MasterMask, std::map<int, odb::dbMaster*>> representative_masters_;
  std::map<std::pair<MasterMask, InstMask>, std::vector<int>> compatible_groups_;
  
  // Compatibility graph
  std::set<Edge> edge_pq_;
  std::unordered_map<int, std::unordered_set<Edge, EdgeHash>> adj_list_;
  FlopClusterRTree feasible_regions_;
  
  // Helper mappings
  std::unordered_map<std::string, int> func_str_to_func_id_;
  std::unordered_map<odb::dbInst*, int> inst_to_flop_id_;
  
  // Spatial grid
  VirtualBinGrid virtual_bin_grid_;

  //==========================================================================
  // Main Algorithm Pipeline
  //==========================================================================
  
  void initVirtualBinGrid();
  void readCompatibleMasters();
  void readFlopUnits();
  void readTimingPaths();
  void createCompatibleGroups();
  void analyzeTimingPaths();
  void calcFeasibleRegions();
  void createFlopClusters();
  void createCompatibilityGraph();
  void runAgglomerativeClustering();
  void implementClusters();

  //==========================================================================
  // Helper Functions
  //==========================================================================

  // --- Phase 1: Virtual Bin Grid Initialization ---
  
  void populateBinGrid();

  // --- Phase 2-3: Master & Flop Analysis ---
  
  bool isValidFlop(odb::dbInst* inst) const;
  InstMask createInstMask(odb::dbInst* inst);
  MasterMask createMasterMask(odb::dbInst* inst);
  
  const sta::LibertyCell* getLibertyCell(odb::dbInst* inst) const;
  int getFuncId(const sta::FuncExpr* expr, odb::dbInst* inst);
  std::string getFuncStr(const sta::FuncExpr* expr, odb::dbInst* inst) const;
  
  bool hasClear(odb::dbInst* inst) const;
  bool hasPreset(odb::dbInst* inst) const;
  bool hasPositiveClockEdge(odb::dbInst* inst) const;
  bool hasScan(odb::dbInst* inst) const;
  
  int getNumDPins(odb::dbInst* inst) const;
  int getNumQPins(odb::dbInst* inst) const;
  int getNumQNPins(odb::dbInst* inst) const;
  
  bool isClockPin(odb::dbITerm* iterm) const;
  bool isDPin(odb::dbITerm* iterm) const;
  bool isQPin(odb::dbITerm* iterm) const;
  bool isQNPin(odb::dbITerm* iterm) const;
  bool isClearPin(odb::dbITerm* iterm) const;
  bool isPresetPin(odb::dbITerm* iterm) const;
  bool isPowerPin(odb::dbITerm* iterm) const;
  bool isScanEnablePin(odb::dbITerm* iterm) const;
  bool isScanInPin(odb::dbITerm* iterm) const;
  
  const sta::LibertyPort* getLibertyPort(odb::dbITerm* iterm) const;
  const sta::FuncExpr* getPortFunc(const sta::LibertyPort* port) const;
  FlopPort getPortType(const sta::LibertyPort* lib_port, odb::dbInst* inst) const;

  // --- Phase 7: Feasible Region Calculation ---
  
  void calcFeasibleRegion(FlopUnit& flop, bool verbose = false);
  void processFanOutPin(FlopUnit& flop, odb::dbITerm* out_pin, odb::dbMTerm* clk_pin_lib,
                        est::EstimateParasitics* est, double unit_r, double unit_c, bool verbose = false);
  void processFanInPin(FlopUnit& flop, odb::dbITerm* in_pin, 
                       est::EstimateParasitics* est, double unit_r, double unit_c, bool verbose = false);
  void computeFinalFeasibleRegion(FlopUnit& flop, const std::vector<odb::dbITerm*>& all_pins, bool verbose = false);
  
  float getPinCapacitance(const sta::Pin* pin) const;
  std::vector<std::pair<float, float>> extractCapacitanceDelayPoints(odb::dbInst* inst,
      const std::string& input_pin_name, const std::string& output_pin_name, float input_slew) const;
  std::vector<std::pair<odb::dbITerm*, float>> getInstanceInputSlews(odb::dbInst* inst) const;
  bool findSteinerPathRecursive(est::SteinerTree* tree, int current_pt, int target_pt, std::vector<int>& path);
  
  float solveMaxDistanceFanOut(float l1, float unit_r, float unit_c, float slack_budget,
                               float coeff, float wo_fst_stt_cap, bool verbose = false) const;
  float solveMaxDistanceFanIn(float l1, float unit_r, float unit_c, float slack_budget,
                              float coeff, float on_path_R_wo_last, float ipin_cap, bool verbose = false) const;
  Box createFeasibleBox(int steiner_x, int steiner_y, float max_dist) const;
  
  bool findTimingArcModel(const sta::LibertyCell* liberty_cell, const std::string& input_pin_name,
      const std::string& output_pin_name, const sta::TableAxis*& capacitance_axis, 
      sta::GateTableModel*& timing_model) const;
  bool findCapacitanceAxis(sta::GateTableModel* gate_model, const sta::TableAxis*& capacitance_axis) const;
  float calculateGateDelay(sta::GateTableModel* timing_model, const sta::Pvt* pvt_conditions,
                           float input_slew, float output_capacitance) const;
  float getTerminalSlew(odb::dbITerm* terminal, sta::Graph* timing_graph, const sta::MinMax* min_max) const;

  // --- Phase 9: Compatibility Graph Construction ---
  
  std::vector<Edge> updateEdges(const FlopCluster& cluster, bool verbose = false);
  std::vector<FlopClusterEntry> getIntersectedCluster(const FlopCluster& cluster, bool verbose = false) const;
  bool checkCompatibility(const FlopCluster& c1, const FlopCluster& c2, bool verbose = false) const;
  odb::dbMaster* getClusterMaster(const FlopCluster& cluster) const;
  PlacementCandidate calcPlacementCandidate(const FlopCluster& c1, const FlopCluster& c2, bool verbose = false);
  
  Box calcMedianBox(const FlopCluster& c1, const FlopCluster& c2) const;
  Box getFeasibleRegionIntersection(const FlopCluster& c1, const FlopCluster& c2) const;
  Point project(const Box& hpwl_box, const Box& feasible_box) const;
  std::vector<Point> generateUniformSamples(const Box& box, int p) const;
  bool checkPlacementDensityConstraint(const FlopCluster& c1, odb::dbMaster* master1,
                                       const FlopCluster& c2, odb::dbMaster* master2,
                                       const Point& new_pos, odb::dbMaster* new_master);
  double calcMergeHPWLGain(const FlopCluster& c1, const FlopCluster& c2,
                           const FloatPoint& merge_position, bool verbose = false) const;
  std::set<std::variant<odb::dbITerm*, odb::dbBTerm*>> getConnectedPins(const FlopCluster& cluster) const;
  FloatPoint getPinCoordinate(const std::variant<odb::dbITerm*, odb::dbBTerm*>& pin_variant) const;

  // --- Phase 10: Agglomerative Clustering ---
  
  int mergeClusters(const Edge& merge_edge, bool verbose = false);
  bool isFurtherMergeable(const FlopCluster& cluster, bool verbose = false) const;
  void removeEdges(const FlopCluster& cluster, bool verbose = false);
  std::set<int> distributeSlack(const FlopCluster& cluster, bool verbose = false);
  void updateFeasibleRegion(FlopCluster& cluster, bool verbose = false);
  
  std::unordered_map<odb::dbITerm*, std::vector<std::pair<int, sta::Slack>>>
  calcUsedSlacks(const FlopUnit& flop, bool verbose = false);
  void calcUsedSlacksFanOut(const FlopUnit& flop, odb::dbITerm* out_pin, odb::dbMTerm* clk_pin_lib,
      est::EstimateParasitics* est, double unit_r, double unit_c,
      std::unordered_map<odb::dbITerm*, std::vector<std::pair<int, sta::Slack>>>& used_slacks, bool verbose = false);
  void calcUsedSlacksFanIn(const FlopUnit& flop, odb::dbITerm* d_pin,
      est::EstimateParasitics* est, double unit_r, double unit_c,
      std::unordered_map<odb::dbITerm*, std::vector<std::pair<int, sta::Slack>>>& used_slacks, bool verbose = false);

  // --- Phase 11: MBFF Implementation ---
  
  bool implementSingleCluster(const FlopCluster& cluster, bool verbose = false, int debug_idx = -1);
  std::vector<NetBundle> getNetBundles(const FlopCluster& cluster) const;
  std::map<odb::dbMaster*, std::vector<PortBundle>> getAllPortBundles(
      const FlopCluster& cluster, const std::vector<NetBundle>& net_bundles) const;
  MasterPortAssignment assignPorts(const FlopCluster& cluster,
                                   const std::vector<NetBundle>& net_bundles,
                                   const std::map<odb::dbMaster*, std::vector<PortBundle>>& master_port_bundles);
  void applyImplementation(const FlopCluster& cluster, const MasterPortAssignment& result,
                           const std::vector<NetBundle>& net_bundles);
  double calcAssignmentCost(const NetBundle& net_bundle, const PortBundle& port_bundle,
                            odb::dbMaster* master, const Point& new_inst_origin) const;
  odb::Rect getNetBBoxWithoutPin(odb::dbNet* net, odb::dbITerm* pin_to_ignore) const;
  Point getGlobalMTermPos(const Point& local_port_pos, odb::dbMaster* master, const Point& inst_center) const;
  std::vector<PortBundle> getPortBundles(odb::dbMaster* master) const;

  // --- General Utilities ---
  
  double dbuToMeters(int dist) const;
  int metersToDbu(double dist) const;
  Point transformCoords(const Point& p) const;
  Point inverseTransformCoords(const Point& p) const;
  unsigned int roundDownToPowerOfTwo(unsigned int x);
  int64_t computeTotalHpwl() const;
};

}  // namespace gpl
