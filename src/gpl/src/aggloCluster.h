#pragma once

#include <map>
#include <ostream>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <variant>
#include <vector>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
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

namespace sta {
class dbNetwork;
class dbSta;
class Corner;
class FuncExpr;
class LibertyCell;
class LibertyPort;
class Path;
}  // namespace sta

namespace gpl {

using Point = bg::point<int, 2, boost::geometry::cs::cartesian>;
using Box = bg::box<Point>;
using FlopEntry = std::pair<Box, int>; // Feasible region & index of FlopCluster
using FlopRTree = bgi::rtree<FlopEntry, bgi::rstar<16, 4, 4>>;

enum FlopPort
{
  d, q, qn, clear, preset, si, se, vdd, vss, func, ifunc, unknown
};

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

  bool operator<(const MasterMask& rhs) const
  {
    return std::tie(func_id_,
                    has_positive_clock_edge_,
                    has_clear_,
                    has_preset_,
                    has_q_pin_,
                    has_qn_pin_,
                    has_scan_)
           < std::tie(rhs.func_id_,
                      rhs.has_positive_clock_edge_,
                      rhs.has_clear_,
                      rhs.has_preset_,
                      rhs.has_q_pin_,
                      rhs.has_qn_pin_,
                      rhs.has_scan_);
  }

  bool operator==(const MasterMask& rhs) const
  {
    return std::tie(func_id_,
                    has_positive_clock_edge_,
                    has_clear_,
                    has_preset_,
                    has_q_pin_,
                    has_qn_pin_,
                    has_scan_)
           == std::tie(rhs.func_id_,
                       rhs.has_positive_clock_edge_,
                       rhs.has_clear_,
                       rhs.has_preset_,
                       rhs.has_q_pin_,
                       rhs.has_qn_pin_,
                       rhs.has_scan_);
  }

  std::string to_string() const
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
};

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

  bool operator<(const InstMask& rhs) const
  {
    return std::tie(clock_net_,
                    clear_net_,
                    preset_net_,
                    scan_enable_net_,
                    scan_in_net_)
           < std::tie(rhs.clock_net_,
                      rhs.clear_net_,
                      rhs.preset_net_,
                      rhs.scan_enable_net_,
                      rhs.scan_in_net_);
  }

  bool operator==(const InstMask& rhs) const
  {
    return std::tie(clock_net_,
                    clear_net_,
                    preset_net_,
                    scan_enable_net_,
                    scan_in_net_)
           == std::tie(rhs.clock_net_,
                       rhs.clear_net_,
                       rhs.preset_net_,
                       rhs.scan_enable_net_,
                       rhs.scan_in_net_);
  }

  std::string to_string() const
  {
    std::ostringstream oss;
    oss << "[Clk:" << (clock_net_ ? clock_net_->getName() : "N/A") << "]"
        << " [Clr:" << (clear_net_ ? clear_net_->getName() : "N/A") << "]"
        << " [Pre:" << (preset_net_ ? preset_net_->getName() : "N/A") << "]"
        << " [SE:" << (scan_enable_net_ ? scan_enable_net_->getName() : "N/A")
        << "]"
        << " [SI:" << (scan_in_net_ ? scan_in_net_->getName() : "N/A") << "]";
    return oss.str();
  }
};

struct TimingPath
{
  sta::Path* path_{nullptr};
  sta::Slack slack_{sta::INF};
  int start_flop_idx_{-1};
  int end_flop_idx_{-1};

  TimingPath(sta::Path* path, sta::Slack slack, int start_idx, int end_idx)
      : path_(path),
        slack_(slack),
        start_flop_idx_(start_idx),
        end_flop_idx_(end_idx)
  {
  }

  bool operator<(const TimingPath& other) const
  {
    return std::tie(slack_, start_flop_idx_, end_flop_idx_, path_)
           < std::tie(other.slack_,
                      other.start_flop_idx_,
                      other.end_flop_idx_,
                      other.path_);
  }
};

struct PlacementCandidate
{
  FloatPoint pos_;
  double gain_;
  bool valid_;

  PlacementCandidate(FloatPoint pos, double gain, bool valid)
      : pos_(pos), gain_(gain), valid_(valid)
  {
  }
};

struct Edge {
  // Two FlopCluste indices 
  int n1, n2;

  // Edge weight (e.g., cost to merge) 
  double weight;

  // Proposed placement location after merging
  FloatPoint pos; 

  Edge(int node1, int node2, double w, FloatPoint p) 
      : weight(w), pos(p) {
      n1 = std::min(node1, node2);
      n2 = std::max(node1, node2);
  }

  bool operator<(const Edge& other) const
  {
      return std::tie(weight, n1, n2, pos.x, pos.y)
              < std::tie(other.weight, other.n1, other.n2, other.pos.x, other.pos.y);
  }

  bool operator==(const Edge& other) const
  {
      return std::tie(n1, n2, weight, pos.x, pos.y)
              == std::tie(other.n1, other.n2, other.weight, other.pos.x, other.pos.y);
  }
};

struct EdgeHash {
  std::size_t operator()(const Edge& e) const {
      std::size_t h1 = std::hash<int>{}(e.n1);
      std::size_t h2 = std::hash<int>{}(e.n2);
      std::size_t h3 = std::hash<double>{}(e.weight);
      std::size_t h4 = FloatPoint::Hash{}(e.pos); 

      return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3);
  }
};

struct FlopUnit
{
  // Index of this flop unit
  int id_{-1};

  // Pointer to the dbInst representing this flop
  odb::dbInst* inst_{nullptr};
  
  // Masks
  MasterMask master_mask_;
  InstMask inst_mask_;

  // Placement information (DBU)
  FloatPoint orig_pt_;
  FloatPoint curr_pt_;

  // Timing Information
  std::vector<int> timing_paths_;  // Indices to timing_paths_ in AggloCluster
  std::unordered_map<odb::dbITerm*, std::vector<std::pair<int, sta::Slack>>> pin_budgets_; // Key: ITerm, Value: vector of (path index, slack budget)
  
  // Feasible region information (Transformed to 45 degree coordinate system) 
  std::unordered_map<odb::dbITerm*, Box> pin_feasible_regions_; // Key: ITerm, Value: feasible region box of the pin
  Box feasible_region_; // Intersection of all pin feasible regions

  // Cluster information
  int cluster_idx_{-1};

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

  std::string to_string() const
  {
    std::ostringstream oss;
    oss << "Inst: " << inst_->getName() << "\n";
    oss << "  - Master Mask: " << master_mask_.to_string() << "\n";
    oss << "  - Inst Mask:   " << inst_mask_.to_string() << "\n";
    oss << "  - Orig Pt:     " << orig_pt_.to_string();
    return oss.str();
  }
};

struct FlopCluster
{
  // Index of this cluster
  int id_{-1};

  // Indices of flops in this cluster
  std::set<int> flops_;

  // Masks
  MasterMask master_mask_;
  InstMask inst_mask_;

  // Placement Information (DBU) 
  FloatPoint curr_pt_;

  // Feasible region (Intersection of all flops' feasible regions)
  Box feasible_region_;

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
        master_mask_(c1.master_mask_), // c1 and c2 should have the same masks
        inst_mask_(c1.inst_mask_), // c1 and c2 should have the same masks
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

class AggloCluster
{
  public:
    AggloCluster(odb::dbDatabase* db,
                  sta::dbSta* sta,
                  utl::Logger* log,
                  rsz::Resizer* resizer,
                  int num_paths_per_endpoint, 
                  int threads);
    ~AggloCluster();

    void doAggloCluster(bool debug);

  private:
    /* OpenROAD */
    odb::dbDatabase* db_;
    odb::dbBlock* block_;
    sta::dbSta* sta_;
    sta::dbNetwork* network_;
    sta::Corner* corner_;
    rsz::Resizer* resizer_;
    utl::Logger* log_;

    /* Config */
    int threads_;
    int num_paths_per_endpoint_; // Number of timing paths to consider per flop endpoint

    /* Key Data */
    // Key: MasterMask, Value: Corresponding compatible dbMasters by bit-width
    std::map<MasterMask, std::map<int, std::vector<odb::dbMaster*>>> compatible_masters_;
    // Key: (MasterMask, InstMask), Value: Indices of compatible FlopUnits
    std::map<std::pair<MasterMask, InstMask>, std::vector<int>> compatible_groups_;
    // Timing paths 
    std::vector<TimingPath> timing_paths_;
    // Flop units
    std::vector<FlopUnit> flop_units_;
    // Flop clusters
    std::vector<FlopCluster> flop_clusters_;
    std::vector<bool> flop_cluster_is_valid_; 
    std::vector<bool> flop_cluster_is_final_;
    // Compatibility graph
    std::set<Edge> edge_pq_; // Min-heap priority queue of Edges
    std::unordered_map<int, std::unordered_set<Edge, EdgeHash>> adj_list_; // Key: FlopCluster index, Value: set of Edges connected to the FlopCluster
    // Feasible regions R-Tree (FlopCluster feasible regions) - 45 degree transformed
    FlopRTree feasible_regions_;
    // Helper maps
    std::unordered_map<std::string, int> func_str_to_func_id_; // Key: function string, Value: function ID 
    std::unordered_map<odb::dbInst*, int> inst_to_flop_id_; // Key: dbInst pointer, Value: index of FlopUnit in flop_units_

    /* Main Functions (Sorted by Steps) */ 
    // (1) Update compatible_masters_ by reading liberty information
    void readCompatibleMasters(); 
    // (2) Intialize flop_units_ 
    void readFlopUnits();
    // (3) Create compatible_groups_ from flop_units_
    void createCompatibleGroups();
    // (4) Read timing paths and associate them with flop_units_
    void readTimingPaths();
    // (5) Analyze timing paths to compute pin budgets for each flop unit
    void analyzeTimingPaths();
    // (6) Calculate feasible regions for each flop unit
    void calcFeasibleRegions();
    // (7) Initialize flop_clusters_ from flop_units_
    void createFlopClusters();
    // (8) Create compatibility graph (edge_pq_ and adj_list_)
    void createCompatibilityGraph();
    // (9) Run the agglomerative clustering process
    void runAgglomerativeClustering();

    /* Helper Functions (Sorted by Usage) */

    // (1) Pin/Port related functions
    const sta::LibertyPort* getLibertyPort(odb::dbITerm* iterm) const;
    const sta::FuncExpr* getPortFunc(const sta::LibertyPort* port) const;
    FlopPort getPortType(const sta::LibertyPort* lib_port, odb::dbInst* inst) const;
    bool isClearPin(odb::dbITerm* iterm) const;
    bool isClockPin(odb::dbITerm* iterm) const;
    bool isDPin(odb::dbITerm* iterm) const;
    bool isPowerPin(odb::dbITerm* iterm) const;
    bool isPresetPin(odb::dbITerm* iterm) const;
    bool isQPin(odb::dbITerm* iterm) const;
    bool isQNPin(odb::dbITerm* iterm) const;
    bool isScanEnablePin(odb::dbITerm* iterm) const;
    bool isScanInPin(odb::dbITerm* iterm) const;

    // (2) Inst/Master related functions
    InstMask createInstMask(odb::dbInst* inst);
    MasterMask createMasterMask(odb::dbInst* inst);
    int getFuncId(const sta::FuncExpr* expr, odb::dbInst* inst);
    std::string getFuncStr(const sta::FuncExpr* expr, odb::dbInst* inst) const;
    const sta::LibertyCell* getLibertyCell(odb::dbInst* inst) const;
    int getNumDPins(odb::dbInst* inst) const;
    int getNumQNPins(odb::dbInst* inst) const;
    int getNumQPins(odb::dbInst* inst) const;
    bool hasClear(odb::dbInst* inst) const;
    bool hasPositiveClockEdge(odb::dbInst* inst) const;
    bool hasPreset(odb::dbInst* inst) const;
    bool hasScan(odb::dbInst* inst) const;
    bool isValidFlop(odb::dbInst* inst) const;

    // (3) FlopUnit related functions
    void calcFeasibleRegion(FlopUnit& flop); // Calculate feasible region for a given FlopUnit
    std::unordered_map<odb::dbITerm*, std::pair<int, sta::Slack>> calcUsedSlacks(const FlopUnit& flop_unit) const;
    odb::dbITerm* getConnectedPinOnPath(const sta::Path* path, const FlopUnit& flop) const; // For a given FlopUnit, find the connected pin on the given timing path

    // (4) FlopCluster related functions
    PlacementCandidate calcPlacementCandidate(const FlopCluster& c1, const FlopCluster& c2) const; // Calculate placement location and gain for merging two FlopClusters
    bool checkCompatibility(const FlopCluster& c1, const FlopCluster& c2) const; // Check if two FlopClusters are compatible for merging (i.e., same masks + available bits)
    std::set<int> distributeSlack(const FlopCluster& cluster); // Distribute slack budget and returns set of affected FlopUnit indices
    std::vector<FlopEntry> getIntersectedCluster(const FlopCluster& cluster) const; // Get FlopClusters whose feasible regions intersect with the given FlopCluster
    bool isFurtherMergeable(const FlopCluster& cluster) const; // Check if a FlopCluster can be further merged with any other cluster
    int mergeClusters(const Edge& merge_edge);
    void removeEdges(const FlopCluster& cluster); // Remove all edges connected to a given FlopCluster
    void updateEdges(const FlopCluster& cluster); // Update edges of a given FlopCluster 
    void updateFeasibleRegion(FlopCluster& cluster); // Calculate feasible region for a given FlopCluster and update feasible_regions_ R-Tree

    // (5) Coordinate transformation functions
    Point transformCoords(const Point& p) const;
    Point inverseTransformCoords(const Point& p) const;  


    
    // float calcHPWL(
    //   const std::vector<int>& flop_indices, 
    //   const std::vector<odb::Point>& cluster_centers) const;
    
    // std::set<std::variant<odb::dbITerm*, odb::dbBTerm*>> getConnectedPins(
    //   const std::vector<int>& flop_indices) const;

    // // helper functions
    // Point transformCoords(const Point& p) const;
    // Point inverseTransformCoords(const Point& p) const;  



    // Report functions
    void reportCompatibleMasters() const;
    void reportCompatibleGroups() const;
    void reportFlopUnits() const;
    void reportFlopPaths() const;
    void reportFeasibleRegions() const;
    void reportTimingPaths() const;
    void reportCompatibilityGraph() const;
};

}  // namespace gpl
