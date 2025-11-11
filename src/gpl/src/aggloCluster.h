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

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <cstdint>
#include <map>
#include <ostream>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <variant>
#include <vector>
#include "odb/db.h"
#include "sta/Delay.hh"
#include "point.h"

namespace bgi = boost::geometry::index;
namespace bg = boost::geometry::model;

namespace utl {
class Logger;
}  // namespace utl

namespace rsz {
class Resizer;
}  // namespace rsz

namespace est {
class EstimateParasitics;
class SteinerTree;
}  // namespace est

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
}  // namespace sta

namespace gpl {

using Point = bg::point<int, 2, boost::geometry::cs::cartesian>;
using Box = bg::box<Point>;
using FlopClusterEntry = std::pair<Box, int>;  
using FlopClusterRTree = bgi::rtree<FlopClusterEntry, bgi::rstar<16, 4, 4>>;

enum FlopPort
{
  d, q, qn, clear, preset, si, se, vdd, vss, func, ifunc, unknown
};

//==============================================================================
// 1. Mask Structures
//==============================================================================

/**
 * @brief Represents the characteristics of an MBFF master cell
 * 
 * Stores various pin characteristics of the master cell (Function ID, clock edge,
 * CLR, PRE, Q, QN, Scan). Used for grouping flops with the same characteristics.
 */
struct MasterMask
{
  int func_id_{-1};
  bool has_positive_clock_edge_{false};
  bool has_clear_{false};
  bool has_preset_{false};
  bool has_q_pin_{false};
  bool has_qn_pin_{false};
  bool has_scan_{false};

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

  bool operator<(const MasterMask& rhs) const
  {
    return _tie() < rhs._tie();
  }

  bool operator==(const MasterMask& rhs) const
  {
    return _tie() == rhs._tie();
  }

  std::string to_string() const;
};

/**
 * @brief Represents the net connections of a flip-flop instance
 * 
 * Stores which nets the individual flip-flop instance is actually connected to
 * (CLK, CLR, PRE, SE, SI nets).
 */
struct InstMask
{
  odb::dbNet* clock_net_{nullptr};
  odb::dbNet* clear_net_{nullptr};
  odb::dbNet* preset_net_{nullptr};
  odb::dbNet* scan_enable_net_{nullptr};
  odb::dbNet* scan_in_net_{nullptr};

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

  auto _tie() const
  {
    return std::tie(clock_net_,
                    clear_net_,
                    preset_net_,
                    scan_enable_net_,
                    scan_in_net_);
  }

  bool operator<(const InstMask& rhs) const
  {
    return _tie() < rhs._tie();
  }

  bool operator==(const InstMask& rhs) const
  {
    return _tie() == rhs._tie();
  }

  std::string to_string() const;
};

//==============================================================================
// 2. Timing Information
//==============================================================================

/**
 * @brief Timing path information between two flip-flops
 * 
 * Stores the timing path and slack information between the source and destination flops.
 * Used for timing constraint analysis.
 */
struct TimingPath
{
  sta::Path* path_{nullptr};
  sta::Slack slack_{sta::INF};
  int start_flop_idx_{-1};
  int end_flop_idx_{-1};
  odb::dbITerm* start_pin_{nullptr};  // Output pin of start flop (Q or QN)
  odb::dbITerm* end_pin_{nullptr};    // Input pin of end flop (D)

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

  auto _tie() const
  {
    return std::tie(slack_, start_flop_idx_, end_flop_idx_, start_pin_, end_pin_, path_);
  }

  bool operator<(const TimingPath& other) const
  {
    return _tie() < other._tie();
  }
};

//==============================================================================
// 3. Spatial Information (Density Constraint Grid)
//==============================================================================

/**
 * @brief Single virtual bin for density constraint checking
 * 
 * Stores the boundary and placed area information of each virtual bin.
 * Used for checking density constraints.
 */
class VirtualBin
{
public:
  VirtualBin() = default;

  VirtualBin(int lx, int ly, int ux, int uy, float target_density)
      : lx_(lx),
        ly_(ly),
        ux_(ux),
        uy_(uy),
        target_density_(target_density)
  {
    bin_area_ = static_cast<int64_t>(ux_ - lx_) * static_cast<int64_t>(uy_ - ly_);
  }

  // Boundary accessors
  int lx() const { return lx_; }
  int ly() const { return ly_; }
  int ux() const { return ux_; }
  int uy() const { return uy_; }

  // Area update methods
  void addInstPlacedArea(int64_t area) { inst_placed_area_ += area; }
  void addMacroPlacedArea(int64_t area) { macro_placed_area_ += area; }
  void addNonPlaceArea(int64_t area) { non_place_area_ += area; }
  void subInstPlacedArea(int64_t area) { inst_placed_area_ -= area; }
  void subMacroPlacedArea(int64_t area) { macro_placed_area_ -= area; }
  void subNonPlaceArea(int64_t area) { non_place_area_ -= area; }

  // Area getters 
  int64_t getInstPlacedArea() const { return inst_placed_area_; }
  int64_t getMacroPlacedArea() const { return macro_placed_area_; }
  int64_t getNonPlaceArea() const { return non_place_area_; }

  int64_t getOverflowArea() const
  {
    // Calculate overflow as: inst_area + macro_area*density + non_place*density - bin_area*density
    const float result = static_cast<float>(inst_placed_area_)
                       + static_cast<float>(macro_placed_area_) * target_density_
                       + static_cast<float>(non_place_area_) * target_density_
                       - static_cast<float>(bin_area_) * target_density_;
    return static_cast<int64_t>(std::max(0.0f, result));
  }

private:
  int lx_, ly_, ux_, uy_;
  float target_density_;
  int64_t bin_area_ = 0;
  int64_t inst_placed_area_ = 0;
  int64_t macro_placed_area_ = 0;
  int64_t non_place_area_ = 0;
};

/**
 * @brief Virtual bin grid: manages all virtual bins in the chip area
 * 
 * Divides the chip area into a uniform grid and manages density information
 * of each bin. Used for checking density constraints and validating placement.
 */
class VirtualBinGrid
{
public:
  VirtualBinGrid() = default;

  VirtualBinGrid(int lx,
                 int ly,
                 int ux,
                 int uy,
                 int bin_cnt_x,
                 int bin_cnt_y,
                 double bin_size_x,
                 double bin_size_y,
                 float target_density,
                 float target_overflow)
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
        int bin_lx = lx_ + std::lround(x * bin_size_x_);
        int bin_ly = ly_ + std::lround(y * bin_size_y_);
        int bin_ux = lx_ + std::lround((x + 1) * bin_size_x_);
        int bin_uy = ly_ + std::lround((y + 1) * bin_size_y_);
        bins_.emplace_back(bin_lx, bin_ly, bin_ux, bin_uy, target_density_);
      }
    }
  }

  // Accessors
  bool checkOverflow();
  std::vector<VirtualBin>& getBins() { return bins_; }
  std::pair<int, int> getMinMaxIdxX(const odb::Rect& box) const;
  std::pair<int, int> getMinMaxIdxY(const odb::Rect& box) const;

  // Modifiers
  void addInst(odb::dbMaster* master, const Point& center_pos);
  void removeInst(odb::dbMaster* master, const Point& center_pos);

private:
  int lx_, ly_, ux_, uy_;
  int bin_cnt_x_, bin_cnt_y_;
  double bin_size_x_, bin_size_y_;
  float target_density_;
  float target_overflow_;
  std::vector<VirtualBin> bins_;

  void updateInstPlacement(odb::dbMaster* master, const Point& center_pos, bool is_add);
};

//==============================================================================
// 4. Graph Information (Compatibility Graph)
//==============================================================================

/**
 * @brief Edge between two mergeable clusters
 * 
 * Represents an edge between two mergeable flip-flop clusters.
 * Contains merge cost (weight) and proposed placement location after merge.
 */
struct Edge
{
  int n1, n2;           // Two FlopCluster indices
  double weight;        // Edge weight (cost to merge)
  FloatPoint pos;       // Proposed placement location after merging

  Edge(int node1, int node2, double w, FloatPoint p)
      : weight(w), pos(p)
  {
    n1 = std::min(node1, node2);
    n2 = std::max(node1, node2);
  }

  auto _tie() const { return std::tie(weight, n1, n2, pos.x, pos.y); }

  bool operator<(const Edge& other) const { return _tie() < other._tie(); }
  bool operator==(const Edge& other) const { return _tie() == other._tie(); }
};

/**
 * @brief Hash function for Edge
 * 
 * Hash function to use edges as keys in unordered_set/unordered_map.
 */
struct EdgeHash
{
  std::size_t operator()(const Edge& e) const
  {
    return std::hash<int>{}(e.n1) ^ (std::hash<int>{}(e.n2) << 1)
           ^ (std::hash<double>{}(e.weight) << 2) ^ (FloatPoint::Hash{}(e.pos) << 3);
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
        cluster_idx_(-1),
        master_mask_(master_mask),
        inst_mask_(inst_mask)
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
        curr_pt_(flop_unit.curr_pt_),
        feasible_region_(flop_unit.feasible_region_),
        master_mask_(flop_unit.master_mask_),
        inst_mask_(flop_unit.inst_mask_)
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
 * @brief Merge-based flip-flop clustering engine
 * 
 * Hierarchically merges compatible flops to convert them into multi-bit flip-flops (MBFFs).
 * Performs optimal clustering considering timing constraints, density constraints, and energy efficiency.
 *
 * Algorithm Flow:
 * 1. Initialize virtual bin grid for density checking
 * 2. Read compatible masters from liberty library
 * 3. Read flip-flop instances and create flop units
 * 4. Create compatible groups based on masks
 * 5. Read timing paths from STA
 * 6. Analyze timing paths and compute slack budgets
 * 7. Calculate feasible regions for each flop unit
 * 8. Initialize flop clusters (one per flop)
 * 9. Create compatibility graph (edges between mergeable clusters)
 * 10. Run agglomerative clustering (greedily merge clusters)
 * 11. Implement final clusters by creating MBFF instances
 */
class AggloCluster
{
  public:
    AggloCluster(odb::dbDatabase* db,
                 sta::dbSta* sta,
                 utl::Logger* log,
                 rsz::Resizer* resizer,
                 float target_density,
                 float target_overflow,
                 int num_paths_per_endpoint,
                 int threads,
                 bool verbose = false);
    ~AggloCluster();

    // Main clustering flow
    void doAggloCluster();

  private:
    //========================================================================
    // External References (OpenROAD)
    //========================================================================
    odb::dbDatabase* db_;
    odb::dbBlock* block_;
    sta::dbSta* sta_;
    sta::dbNetwork* network_;
    sta::Corner* corner_;
    rsz::Resizer* resizer_;
    utl::Logger* log_;

    //========================================================================
    // Configuration Parameters
    //========================================================================
    int threads_;
    int num_paths_per_endpoint_; // Number of timing paths to consider per flop endpoint
    float target_density_;
    float target_overflow_;
    bool verbose_;

    //========================================================================
    // Master & Compatibility Information
    //========================================================================
    // Key: MasterMask, Value: map of (bit-width -> vector of compatible dbMasters)
    std::map<MasterMask, std::map<int, std::vector<odb::dbMaster*>>> compatible_masters_;
    // Key: MasterMask, Value: map of (bit-width -> largest compatible dbMaster)
    std::map<MasterMask, std::map<int, odb::dbMaster*>> representative_masters_;
    // Key: (MasterMask, InstMask), Value: indices of compatible FlopUnits
    std::map<std::pair<MasterMask, InstMask>, std::vector<int>> compatible_groups_;

    //========================================================================
    // Flop & Cluster Data
    //========================================================================
    std::vector<FlopUnit> flop_units_;
    std::vector<FlopCluster> flop_clusters_;
    std::vector<bool> flop_cluster_is_valid_;
    std::vector<bool> flop_cluster_no_further_merge_;

    //========================================================================
    // Timing Information
    //========================================================================
    std::vector<TimingPath> timing_paths_;

    //========================================================================
    // Compatibility Graph (for Agglomerative Clustering)
    //========================================================================
    std::set<Edge> edge_pq_; // Priority queue of mergeable edges (sorted by weight)
    std::unordered_map<int, std::unordered_set<Edge, EdgeHash>> adj_list_; // Adjacency list: Key=FlopCluster index, Value=edges connected to it
    FlopClusterRTree feasible_regions_; // R-Tree for fast spatial queries on feasible regions (45-degree transformed)

    //========================================================================
    // Helper Mappings & Spatial Grid
    //========================================================================
    std::unordered_map<std::string, int> func_str_to_func_id_;
    std::unordered_map<odb::dbInst*, int> inst_to_flop_id_;
    VirtualBinGrid virtual_bin_grid_;

    //========================================================================
    // Main Algorithm Steps (in execution order)
    //========================================================================
    
    // ---- Initialization & Data Collection ----
    void initVirtualBinGrid();           // Initialize virtual bin grid for density checking
    void readCompatibleMasters();        // Read compatible MBFF masters from Liberty
    void readFlopUnits();                // Read flip-flop instances and create FlopUnits
    void readTimingPaths();              // Extract timing paths from STA
    
    // ---- Grouping & Analysis ----
    void createCompatibleGroups();       // Group flops by (MasterMask, InstMask)
    void analyzeTimingPaths();           // Compute slack budgets for each pin
    void calcFeasibleRegions();          // Calculate feasible regions from slack budgets
    
    // ---- Graph Construction ----
    void createFlopClusters();           // Initialize clusters (one per flop)
    void createCompatibilityGraph();     // Build graph with mergeable cluster edges
    
    // ---- Clustering & Implementation ----
    void runAgglomerativeClustering();   // Greedily merge clusters by edge weights
    void implementClusters();            // Convert final clusters into MBFF instances

    //========================================================================
    // Helper Functions (grouped by algorithm phase, organized by hierarchy)
    //========================================================================

    // ╔═══════════════════════════════════════════════════════════════════╗
    // ║ Phase 2-3: readCompatibleMasters() & readFlopUnits()              ║
    // ╚═══════════════════════════════════════════════════════════════════╝
    
    // [Level 1] Flop validation
    bool isValidFlop(odb::dbInst* inst) const;
    
    // [Level 1] Mask creation (top-level)
    InstMask createInstMask(odb::dbInst* inst);
    MasterMask createMasterMask(odb::dbInst* inst);
    
    // [Level 2] Liberty cell and function analysis
    const sta::LibertyCell* getLibertyCell(odb::dbInst* inst) const;
    int getFuncId(const sta::FuncExpr* expr, odb::dbInst* inst);
    std::string getFuncStr(const sta::FuncExpr* expr, odb::dbInst* inst) const;
    
    // [Level 2] Master property checkers
    bool hasClear(odb::dbInst* inst) const;
    bool hasPreset(odb::dbInst* inst) const;
    bool hasPositiveClockEdge(odb::dbInst* inst) const;
    bool hasScan(odb::dbInst* inst) const;
    
    // [Level 2] Pin counting utilities
    int getNumDPins(odb::dbInst* inst) const;
    int getNumQPins(odb::dbInst* inst) const;
    int getNumQNPins(odb::dbInst* inst) const;
    
    // [Level 3] Pin type checkers (used by mask creation and validation)
    bool isClockPin(odb::dbITerm* iterm) const;
    bool isDPin(odb::dbITerm* iterm) const;
    bool isQPin(odb::dbITerm* iterm) const;
    bool isQNPin(odb::dbITerm* iterm) const;
    bool isClearPin(odb::dbITerm* iterm) const;
    bool isPresetPin(odb::dbITerm* iterm) const;
    bool isPowerPin(odb::dbITerm* iterm) const;
    bool isScanEnablePin(odb::dbITerm* iterm) const;
    bool isScanInPin(odb::dbITerm* iterm) const;
    
    // [Level 3] Port analysis utilities
    const sta::LibertyPort* getLibertyPort(odb::dbITerm* iterm) const;
    const sta::FuncExpr* getPortFunc(const sta::LibertyPort* port) const;
    FlopPort getPortType(const sta::LibertyPort* lib_port, odb::dbInst* inst) const;

    // ╔═══════════════════════════════════════════════════════════════════╗
    // ║ Phase 7: calcFeasibleRegions()                                    ║
    // ╚═══════════════════════════════════════════════════════════════════╝
    
    // [Level 1] Top-level: Compute feasible region for each flop unit
    void calcFeasibleRegion(FlopUnit& flop, bool verbose = false, double unit_r=7.5e5, double unit_c=16e-12);
    
    // [Level 2] Main processing functions for different pin types
    struct PinClassification {
      std::vector<odb::dbITerm*> d_pins;
      std::vector<odb::dbITerm*> q_pins;
      std::vector<odb::dbITerm*> qn_pins;
      odb::dbITerm* clk_pin;
    };
    
    PinClassification classifyFlopPins(odb::dbInst* inst) const;
    void processFanOutPin(FlopUnit& flop, odb::dbITerm* out_pin, odb::dbMTerm* clk_pin_lib, est::EstimateParasitics* est, double unit_r, double unit_c, bool verbose = false);
    void processFanInPin(FlopUnit& flop, odb::dbITerm* in_pin, est::EstimateParasitics* est, double unit_r, double unit_c, bool verbose = false);
    void computeFinalFeasibleRegion(FlopUnit& flop, const std::vector<odb::dbITerm*>& all_pins, bool verbose = false);
    
    // [Level 3] Utility functions for timing and parasitic analysis
    float getPinCapacitance(const sta::Pin* pin) const;
    std::vector<std::pair<float, float>> extractCapacitanceDelayPoints(odb::dbInst* inst, const std::string& input_pin_name, const std::string& output_pin_name, float input_slew) const;
    std::vector<std::pair<odb::dbITerm*, float>> getInstanceInputSlews(odb::dbInst* inst) const;
    bool findSteinerPathRecursive(est::SteinerTree* tree, int current_pt, int target_pt, std::vector<int>& path);
    
    // [Level 3] Quadratic equation solving for max distance calculation
    float solveMaxDistanceFanOut(float l1, float unit_r, float unit_c, float slack_budget, float coeff, float wo_fst_stt_cap, bool verbose = false) const;
    float solveMaxDistanceFanIn(float l1, float unit_r, float unit_c, float slack_budget, float coeff, float on_path_R_wo_last, float ipin_cap, bool verbose = false) const;
    Box createFeasibleBox(int steiner_x, int steiner_y, float max_dist) const;
    
    // [Level 4] Helper functions for extractCapacitanceDelayPoints
    bool findTimingArcModel(const sta::LibertyCell* liberty_cell, const std::string& input_pin_name, const std::string& output_pin_name, const sta::TableAxis*& capacitance_axis, sta::GateTableModel*& timing_model) const;
    bool findCapacitanceAxis(sta::GateTableModel* gate_model, const sta::TableAxis*& capacitance_axis) const;
    float calculateGateDelay(sta::GateTableModel* timing_model, const sta::Pvt* pvt_conditions, float input_slew, float output_capacitance) const;
    
    // [Level 4] Helper function for getInstanceInputSlews
    float getTerminalSlew(odb::dbITerm* terminal, sta::Graph* timing_graph, const sta::MinMax* min_max) const;

    // ╔═══════════════════════════════════════════════════════════════════╗
    // ║ Phase 8: createFlopClusters()                                    ║
    // ╚═══════════════════════════════════════════════════════════════════╝
    // (No helper functions needed, uses initialized FlopUnits)

    // ╔═══════════════════════════════════════════════════════════════════╗
    // ║ Phase 9: createCompatibilityGraph()                              ║
    // ╚═══════════════════════════════════════════════════════════════════╝
    // [Level 1] Top-level: Check if two clusters can be merged
    bool checkCompatibility(const FlopCluster& c1, const FlopCluster& c2) const;
    PlacementCandidate calcPlacementCandidate(const FlopCluster& c1, const FlopCluster& c2);
    double calcHPWLDiff(const FlopCluster& c1, const FlopCluster& c2, const FloatPoint& new_pos) const;
    
    //   [Level 2] Helpers for placement analysis
    Box calcMedianBox(const FlopCluster& c1, const FlopCluster& c2) const;
    Box getFeasibleRegionIntersection(const FlopCluster& c1, const FlopCluster& c2) const;
    bool checkPlacementDensityConstraint(const FlopCluster& c1, odb::dbMaster* master1, const FlopCluster& c2, odb::dbMaster* master2, const Point& new_pos, odb::dbMaster* new_master);

    // ╔═══════════════════════════════════════════════════════════════════╗
    // ║ Phase 10: runAgglomerativeClustering()                           ║
    // ╚═══════════════════════════════════════════════════════════════════╝
    // [Level 1] Top-level: Execute merge operations and maintain graph
    int mergeClusters(const Edge& merge_edge);
    bool isFurtherMergeable(const FlopCluster& cluster) const;
    
    //   [Level 2] Helpers for edge/region management
    void removeEdges(const FlopCluster& cluster);
    void updateEdges(const FlopCluster& cluster);
    void updateFeasibleRegion(FlopCluster& cluster);
    std::set<std::variant<odb::dbITerm*, odb::dbBTerm*>> getConnectedPins(const FlopCluster& cluster) const;
    std::vector<FlopClusterEntry> getIntersectedCluster(const FlopCluster& cluster) const;
    std::set<int> distributeSlack(const FlopCluster& cluster);
    // [Level 1] Top-level: Calculate slack budgets for each pin
    std::unordered_map<odb::dbITerm*, std::pair<int, sta::Slack>> calcUsedSlacks(const FlopUnit& flop_unit) const;

    // ╔═══════════════════════════════════════════════════════════════════╗
    // ║ Phase 11: implementClusters()                                    ║
    // ╚═══════════════════════════════════════════════════════════════════╝
    // [Level 1] Top-level: Convert clusters to MBFF instances
    void implementSingleCluster(const FlopCluster& cluster);
    
    //   [Level 2] Helpers for master selection and port assignment
    MasterPortAssignment findBestMasterAndAssignment(const FlopCluster& cluster, const std::vector<NetBundle>& net_bundles);
    std::vector<NetBundle> getNetBundles(const FlopCluster& cluster) const;
    std::vector<PortBundle> getPortBundles(odb::dbMaster* master) const;
    double calcAssignmentCost(const NetBundle& net_bundle, const PortBundle& port_bundle, odb::dbMaster* master, const Point& new_inst_origin) const;
    
    //     [Level 3] Lower-level helpers for implementation
    void applyImplementation(const FlopCluster& cluster, const MasterPortAssignment& result, const std::vector<NetBundle>& net_bundles);
    odb::Rect getNetBBoxWithoutPin(odb::dbNet* net, odb::dbITerm* pin_to_ignore) const;
    Point getGlobalMTermPos(const Point& local_port_pos, odb::dbMaster* master, const Point& inst_center) const;
    FloatPoint getPinCoordinate(const std::variant<odb::dbITerm*, odb::dbBTerm*>& pin_variant) const;

    // ╔═══════════════════════════════════════════════════════════════════╗
    // ║ General Utilities (used across multiple phases)                   ║
    // ╚═══════════════════════════════════════════════════════════════════╝

    // Coordinate transformation
    Point transformCoords(const Point& p) const;
    Point inverseTransformCoords(const Point& p) const;
    Point project(const Box& hpwl_box, const Box& feasible_box) const;
    
    // Math utilities
    unsigned int roundDownToPowerOfTwo(unsigned int x);
    double calcManhattanDistance(const FloatPoint& p1, const FloatPoint& p2) const;
    std::vector<Point> generateUniformSamples(const Box& box, int p) const;

};

}  // namespace gpl
