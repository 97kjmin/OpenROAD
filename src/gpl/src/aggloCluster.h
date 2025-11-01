#pragma once

#include <map>
#include <unordered_map>
#include <vector>
#include <string>
#include <ostream>

#include "odb/db.h"
#include "point.h"
#include "sta/Delay.hh"

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

enum FlopPort
{
  d,
  q,
  qn,
  clear,
  preset,
  si,
  se,
  vdd,
  vss,
  func,
  ifunc,
  unknown
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

  MasterMask(int func,
             bool clk,
             bool clear,
             bool preset,
             bool qpin,
             bool qnpin,
             bool scan)
      : func_id_(func),
        has_positive_clock_edge_(clk),
        has_clear_(clear),
        has_preset_(preset),
        has_q_pin_(qpin),
        has_qn_pin_(qnpin),
        has_scan_(scan)
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

  InstMask(odb::dbNet* clk,
           odb::dbNet* clear,
           odb::dbNet* preset,
           odb::dbNet* scan_enable,
           odb::dbNet* scan_in)
      : clock_net_(clk),
        clear_net_(clear),
        preset_net_(preset),
        scan_enable_net_(scan_enable),
        scan_in_net_(scan_in)
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
  int start_flop_idx_{-1};
  int end_flop_idx_{-1};
  sta::Slack slack_{ sta::INF };
  sta::Path* path_{nullptr};

  TimingPath(int start_idx, int end_idx, sta::Slack slack, sta::Path* path)
      : start_flop_idx_(start_idx),
        end_flop_idx_(end_idx),
        slack_(slack),
        path_(path)
  {
  }

  bool operator<(const TimingPath& other) const
  {
    return slack_ < other.slack_;
  }
};

struct FlopUnit
{
  odb::dbInst* inst_{nullptr};
  FloatPoint orig_pt_;
  FloatPoint curr_pt_;
  int cluster_idx_{-1};
  MasterMask master_mask_;
  InstMask inst_mask_;

  FlopUnit(odb::dbInst* inst,
           const odb::Rect& bbox,
           const MasterMask& master_mask,
           const InstMask& inst_mask)
      : inst_(inst),
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
  // OpenROAD vars
  odb::dbDatabase* db_;
  odb::dbBlock* block_;
  sta::dbSta* sta_;
  sta::dbNetwork* network_;
  sta::Corner* corner_;
  rsz::Resizer* resizer_;
  utl::Logger* log_;

  // Config vars
  int threads_;
  int num_paths_per_endpoint_;
  float multiplier_;

  // Internal vars
  std::map<MasterMask, std::map<int, std::vector<odb::dbMaster*>>>
      compatible_masters_;
  std::map<std::pair<MasterMask, InstMask>, std::vector<int>> compatible_groups_;
  std::unordered_map<std::string, int> func_str_to_func_id_;
  std::unordered_map<odb::dbInst*, int> inst_to_flop_id_;
  std::vector<FlopUnit> flop_units_;
  std::vector<TimingPath> timing_paths_;
  std::unordered_map<int, std::vector<int>> paths_by_start_flop_id_;
  std::unordered_map<int, std::vector<int>> paths_by_end_flop_id_;

  // Main functions
  void createCompatibleGroups();
  void readCompatibleMasters();
  void readFlopUnits();
  void readTimingPaths();
  void reportTimingPaths() const;

  // Pin/Port related functions
  const sta::LibertyPort* getLibertyPort(odb::dbITerm* iterm) const;
  const sta::FuncExpr* getPortFunc(const sta::LibertyPort* port) const;
  FlopPort getPortType(const sta::LibertyPort* lib_port,
                       odb::dbInst* inst) const;
  bool isClearPin(odb::dbITerm* iterm) const;
  bool isClockPin(odb::dbITerm* iterm) const;
  bool isDPin(odb::dbITerm* iterm) const;
  bool isPowerPin(odb::dbITerm* iterm) const;
  bool isPresetPin(odb::dbITerm* iterm) const;
  bool isScanInPin(odb::dbITerm* iterm) const;
  bool isScanEnablePin(odb::dbITerm* iterm) const;
  bool isQPin(odb::dbITerm* iterm) const;
  bool isQNPin(odb::dbITerm* iterm) const;

  // Inst/Master related functions
  MasterMask createMasterMask(odb::dbInst* inst);
  InstMask createInstMask(odb::dbInst* inst);
  int getFuncId(const sta::FuncExpr* expr, odb::dbInst* inst);
  std::string getFuncStr(const sta::FuncExpr* expr, odb::dbInst* inst) const;
  const sta::LibertyCell* getLibertyCell(
      odb::dbInst* inst) const;  
  int getNumDPins(odb::dbInst* inst) const;
  int getNumQPins(odb::dbInst* inst) const;
  int getNumQNPins(odb::dbInst* inst) const;
  bool hasClear(odb::dbInst* inst) const;
  bool hasPositiveClockEdge(odb::dbInst* inst) const;
  bool hasPreset(odb::dbInst* inst) const;
  bool hasScan(odb::dbInst* inst) const;
  bool isValidFlop(odb::dbInst* inst) const;

  // Report functions 
  void reportCompatibleMasters() const;
  void reportCompatibleGroups() const;
  void reportFlopUnits() const;
};

}  // namespace gpl
