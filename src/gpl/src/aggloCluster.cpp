#include "aggloCluster.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>

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

namespace gpl {

// Constructor
AggloCluster::AggloCluster(odb::dbDatabase* db,
                           sta::dbSta* sta,
                           utl::Logger* log,
                           rsz::Resizer* resizer,
                           int num_paths_per_endpoint,
                           int threads)
    : db_(db),
      block_(db->getChip()->getBlock()),
      sta_(sta),
      network_(sta_->getDbNetwork()),
      corner_(sta_->cmdCorner()),
      log_(log),
      resizer_(resizer),
      num_paths_per_endpoint_(num_paths_per_endpoint),
      threads_(threads)
{
}

// Destructor
AggloCluster::~AggloCluster() = default;

// Core Agglomerative Clustering Function
void AggloCluster::doAggloCluster(bool debug)
{
  readCompatibleMasters();
  if (debug) {
    reportCompatibleMasters();
  }

  readFlopUnits();
  if (debug) {
    reportFlopUnits();
  }

  createCompatibleGroups();
  if (debug) {
    reportCompatibleGroups();
  }

  readTimingPaths();
  if (debug) {
    reportTimingPaths();
  }

  analyzeTimingPaths();
  if (debug) {
    reportFlopPaths();
  }

  calcFeasibleRegions();
  if (debug) {
    reportFeasibleRegions();
  }

  createFlopClusters();

  createCompatibilityGraph();
  if (debug) {
    reportCompatibilityGraph();
  }
}

/* Main Functions */
void 
AggloCluster::readCompatibleMasters()
{
  compatible_masters_.clear();
  func_str_to_func_id_.clear();

  std::map<MasterMask, std::map<int, std::vector<odb::dbMaster*>>>
      temp_master_map;

  const char* temp_inst_name = "_temp_master_check";

  for (odb::dbLib* lib : db_->getLibs()) {
    for (odb::dbMaster* master : lib->getMasters()) {
      odb::dbInst* temp_inst
          = odb::dbInst::create(block_, master, temp_inst_name);
      if (temp_inst == nullptr) {
        continue;
      }

      if (isValidFlop(temp_inst)) {
        MasterMask master_mask = createMasterMask(temp_inst);
        int bits = getNumDPins(temp_inst);
        temp_master_map[master_mask][bits].push_back(master);
      }

      odb::dbInst::destroy(temp_inst);
    }
  }

  // Update compatible_masters_ with groups that have multi-bit masters
  std::copy_if(temp_master_map.begin(),
               temp_master_map.end(),
               std::inserter(compatible_masters_, compatible_masters_.end()),
               [](const auto& pair) { return pair.second.size() >= 2; });
}

void 
AggloCluster::readFlopUnits()
{
  flop_units_.clear();
  inst_to_flop_id_.clear();

  int flop_id = 0;
  for (odb::dbInst* inst : block_->getInsts()) {
    if (!isValidFlop(inst)) {
      continue;
    }

    flop_units_.emplace_back(flop_id,
                             inst,
                             inst->getBBox()->getBox(),
                             createMasterMask(inst),
                             createInstMask(inst));
    inst_to_flop_id_[inst] = flop_id++;
  }
}

void 
AggloCluster::createCompatibleGroups()
{
  compatible_groups_.clear();

  for (size_t i = 0; i < flop_units_.size(); ++i) {
    const FlopUnit& unit = flop_units_[i];
    compatible_groups_[{unit.master_mask_, unit.inst_mask_}].push_back(i);
  }
}

void 
AggloCluster::readTimingPaths()
{
  timing_paths_.clear();
  for (auto& flop_unit : flop_units_) {
    flop_unit.timing_paths_.clear();
  }

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

  for (sta::PathEnd* path_end : path_ends) {
    sta::Slack slack = path_end->slack(sta_);
    sta::Path* path = path_end->path();
    sta::PathExpanded expanded(path, sta_);

    if (expanded.size() < 2) {
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
      continue;
    }

    int start_flop_id = inst_to_flop_id_[db_start_inst];
    int end_flop_id = inst_to_flop_id_[db_end_inst];

    timing_paths_.emplace_back(path, slack, start_flop_id, end_flop_id);
    int path_idx = timing_paths_.size() - 1;

    flop_units_[start_flop_id].timing_paths_.push_back(path_idx);
    flop_units_[end_flop_id].timing_paths_.push_back(path_idx);
  }
}

void 
AggloCluster::analyzeTimingPaths()
{
  for (FlopUnit& flop : flop_units_) {
    flop.pin_budgets_.clear();
  }

  for (size_t i = 0; i < timing_paths_.size(); ++i) {
    const auto& path = timing_paths_[i];

    const sta::Slack budget = (path.slack_ > 0) ? path.slack_ / 2.0f : 0.0f; // Assign budget 0 if slack is non-positive

    // Start Flop
    FlopUnit& start_flop = flop_units_[path.start_flop_idx_];
    if (odb::dbITerm* pin = getConnectedPinOnPath(path.path_, start_flop)) {
      start_flop.pin_budgets_[pin].emplace_back(i, budget);
    }

    // End Flop
    FlopUnit& end_flop = flop_units_[path.end_flop_idx_];
    if (odb::dbITerm* pin = getConnectedPinOnPath(path.path_, end_flop)) {
      end_flop.pin_budgets_[pin].emplace_back(i, budget);
    }
  } 

  // Sort pin budgets for each flop unit (ascending order of budget)
  for (FlopUnit& flop : flop_units_) {
    for (auto& pair : flop.pin_budgets_) {
      std::sort(pair.second.begin(),
                pair.second.end(),
                [](const auto& a, const auto& b) { return a.second < b.second; });
    }
  }
}

void 
AggloCluster::calcFeasibleRegions()
{
    for (FlopUnit& flop : flop_units_) {
        calcFeasibleRegion(flop);
    }
}

void 
AggloCluster::createFlopClusters()
{
  flop_clusters_.clear();
  flop_clusters_.reserve(flop_units_.size());
  feasible_regions_.clear();
  flop_cluster_is_valid_.assign(flop_units_.size(), true);
  flop_cluster_is_final_.assign(flop_units_.size(), false);

  for (size_t i = 0; i < flop_units_.size(); ++i) {
    flop_clusters_.emplace_back(i, flop_units_[i]);
    flop_units_[i].cluster_idx_ = i;

    const Box& feasible_region = flop_units_[i].feasible_region_;
    feasible_regions_.insert(std::make_pair(feasible_region, i));
  }
}

void
AggloCluster::createCompatibilityGraph()
{
  edge_pq_.clear();
  adj_list_.clear();

  for (size_t i = 0; i < flop_clusters_.size(); ++i) {
    updateEdges(flop_clusters_[i]);
  }
}

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
      flop_cluster_is_final_[n_new_idx] = true;

      // Perform slack distribution if final cluster is formed      
      std::set<int> affected_flop_units = distributeSlack(n_new); 
      
      std::set<int> affected_flop_clusters;
      
      for (int u_idx : affected_flop_units) {
        // Update feasible region of affected 'FlopUnit's
        calcFeasibleRegion(flop_units_[u_idx]); 
        
        // Get affected 'Cluster's for later FR and edge updates
        int c_idx = flop_units_[u_idx].cluster_idx_;
        if (flop_cluster_is_valid_[c_idx] && !flop_cluster_is_final_[c_idx]) {
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

/* Helper Functions */
// (1) Pin/Port related functions
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
AggloCluster::getPortType(const sta::LibertyPort* lib_port,odb::dbInst* inst) const
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
    log_->warn(utl::GPL,
               9993,
               "Could not find Liberty cell for instance {}",
               inst->getName());
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

  log_->warn(utl::GPL, 9992, "Could not recognize port {}", lib_port->name());
  return unknown;
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

// (2) Inst/Master related functions

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
        log_->warn(utl::GPL,
                   9998,
                   "Instance {} has multiple different clock nets. Using the "
                   "first one found ({}).",
                   inst->getName(),
                   clock_net->getName());
      }
    } else if (isClearPin(iterm)) {
      if (clear_net == nullptr) {
        clear_net = net;
      } else if (clear_net != net) {
        log_->warn(utl::GPL,
                   9997,
                   "Instance {} has multiple different clear nets. Using the "
                   "first one found ({}).",
                   inst->getName(),
                   clear_net->getName());
      }
    } else if (isPresetPin(iterm)) {
      if (preset_net == nullptr) {
        preset_net = net;
      } else if (preset_net != net) {
        log_->warn(utl::GPL,
                   9996,
                   "Instance {} has multiple different preset nets. Using the "
                   "first one found ({}).",
                   inst->getName(),
                   preset_net->getName());
      }
    } else if (isScanEnablePin(iterm)) {
      if (scan_enable_net == nullptr) {
        scan_enable_net = net;
      } else if (scan_enable_net != net) {
        log_->warn(utl::GPL,
                   9995,
                   "Instance {} has multiple different scan enable nets. Using "
                   "the first one found ({}).",
                   inst->getName(),
                   scan_enable_net->getName());
      }
    } else if (isScanInPin(iterm)) {
      if (scan_in_net == nullptr) {
        scan_in_net = net;
      } else if (scan_in_net != net) {
        log_->warn(utl::GPL,
                   9994,
                   "Instance {} has multiple different scan in nets. Using the "
                   "first one found ({}).",
                   inst->getName(),
                   scan_in_net->getName());
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
AggloCluster::hasScan(odb::dbInst* inst) const
{
  const sta::LibertyCell* lib_cell = getLibertyCell(inst);
  if (lib_cell == nullptr) {
    return false;
  }
  return lib_cell && getLibertyScanIn(lib_cell)
         && getLibertyScanEnable(lib_cell);
}

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
    } else {
      log_->warn(utl::GPL,
                 9991,
                 "Instance {} / Master {} has an unrecognized pin {} - "
                 "considering invalid.",
                 inst->getName(),
                 inst->getMaster()->getName(),
                 iterm->getMTerm()->getName());
      return false;
    }
  }

  if (d_pins != std::max(q_pins, qn_pins)) {
    return false;
  }

  return true;
}

// (3) FlopUnit related functions

void
AggloCluster::calcFeasibleRegion(FlopUnit& flop)
{
  throw std::logic_error("calcFeasibleRegion not implemented");
}

std::unordered_map<odb::dbITerm*, std::pair<int, sta::Slack>> 
AggloCluster::calcUsedSlacks(const FlopUnit& flop_unit) const 
{
  throw std::logic_error("calcUsedSlacks not implemented");
}

odb::dbITerm*
AggloCluster::getConnectedPinOnPath(const sta::Path* path,
                                  const FlopUnit& flop) const
{
  sta::PathExpanded expanded(path, sta_);

  if (expanded.size() < 2) {
    return nullptr;
  }

  const sta::Path* start_path = expanded.path(expanded.startIndex());
  sta::Pin* start_pin = (start_path) ? start_path->vertex(sta_)->pin() : nullptr;
  odb::dbInst* start_inst
      = (start_pin) ? network_->staToDb(network_->instance(start_pin)) : nullptr;

  if (start_inst != nullptr && start_inst == flop.inst_) {
    const sta::Path* next_path = expanded.path(expanded.startIndex() + 1);
    sta::Pin* next_pin = (next_path) ? next_path->vertex(sta_)->pin() : nullptr;
    
    return (next_pin) ? network_->flatPin(next_pin) : nullptr;
  }

  const sta::Path* end_path = expanded.path(expanded.size() - 1);
  sta::Pin* end_pin = (end_path) ? end_path->vertex(sta_)->pin() : nullptr;
  odb::dbInst* end_inst
      = (end_pin) ? network_->staToDb(network_->instance(end_pin)) : nullptr;

  if (end_inst != nullptr && end_inst == flop.inst_) {
    const sta::Path* prev_path = expanded.path(expanded.size() - 2);
    sta::Pin* prev_pin = (prev_path) ? prev_path->vertex(sta_)->pin() : nullptr;

    return (prev_pin) ? network_->flatPin(prev_pin) : nullptr;
  }

  return nullptr;
}

// (4) FlopCluster related functions

PlacementCandidate
AggloCluster::calcPlacementCandidate(const FlopCluster& c1,
                                     const FlopCluster& c2) const
{
  throw std::logic_error("calcMergeLocation not implemented");
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

std::set<int>
AggloCluster::distributeSlack(const FlopCluster& cluster)
{
  std::set<int> affected_flop_units;

  for (int flop_idx : cluster.flops_) {
    const FlopUnit& flop = flop_units_[flop_idx];

    // Calculate used slacks for each pin of the flop 
    std::unordered_map<odb::dbITerm*, std::pair<int, sta::Slack>> used_slacks = calcUsedSlacks(flop);

    for (const auto& [pin, used_slack_pair] : used_slacks) {
      const auto& [path_idx, used_slack] = used_slack_pair;
      const TimingPath& path = timing_paths_[path_idx];

      // Find the flop on the other side of the timing path
      int other_flop_idx = (path.start_flop_idx_ == flop_idx)
                               ? path.end_flop_idx_
                               : path.start_flop_idx_;

      // If the other flop is also in the same cluster, no need to distribute slack
      if (flop_units_[other_flop_idx].cluster_idx_ == cluster.id_) {
        continue;
      }

      FlopUnit& other_flop = flop_units_[other_flop_idx];
      affected_flop_units.insert(other_flop_idx);

      // Find the connected pin on the other flop for this path
      odb::dbITerm* other_pin = getConnectedPinOnPath(path.path_, other_flop);
      if (other_pin == nullptr
          || other_flop.pin_budgets_.find(other_pin)
                 == other_flop.pin_budgets_.end()) {
        continue;
      }

      auto& other_budgets = other_flop.pin_budgets_.at(other_pin);
      for (auto& budget_pair : other_budgets) {
        if (budget_pair.first == path_idx) {
          // Find the original budget allocated to this path on the current flop
          sta::Slack original_budget = 0.0; 
          const auto& budgets = flop.pin_budgets_.at(pin);
          for (const auto& orig_budget_pair : budgets) {
            if (orig_budget_pair.first == path_idx) {
              original_budget = orig_budget_pair.second;
              break;
            }
          }

          // Update the budget on the other flop by adding the slack delta
          sta::Slack slack_delta = original_budget - used_slack;
          budget_pair.second += slack_delta;
          break;
        }
      }

      // Sort the budgets for the updated pin
      std::sort(other_budgets.begin(),
                other_budgets.end(),
                [](const auto& a, const auto& b) { return a.second < b.second; });
    }
  }
  return affected_flop_units;
}

std::vector<FlopEntry>
AggloCluster::getIntersectedCluster(const FlopCluster& cluster) const
{
  std::vector<FlopEntry> intersecting_clusters;

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
  flop_cluster_is_final_.push_back(false);
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

  if (!flop_cluster_is_valid_[i] || flop_cluster_is_final_[i]) {
    return;
  }

  if (boost::geometry::is_empty(cluster.feasible_region_)) {
    return;
  }

  std::vector<FlopEntry> intersecting_entries = getIntersectedCluster(cluster);

  for (const auto& entry : intersecting_entries) {
    const int j = entry.second;

    if (i >= j) {
      continue;
    }

    if (!flop_cluster_is_valid_[j] || flop_cluster_is_final_[j]) {
      continue;
    }

    const FlopCluster& adj_cluster = flop_clusters_[j];

    if (!checkCompatibility(cluster, adj_cluster)) {
      continue;
    }

    PlacementCandidate pc = calcPlacementCandidate(cluster, adj_cluster);

    if (!pc.valid_) {
      continue;
    }

    Edge new_edge(i, j, pc.gain_, pc.pos_);

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

// (5) Coordinate transformation functions

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



// float AggloCluster::calcHPWL(const std::vector<int>& flop_indices,
//                                   const std::vector<odb::Point>& cluster_centers) const
// {
//   std::set<odb::dbNet*> nets;
//   std::set<odb::dbInst*> cluster_insts;
//   for (int flop_idx : flop_indices) {
//     cluster_insts.insert(flop_units_[flop_idx].inst_);
//   }

//   auto pins = getConnectedPins(flop_indices);
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





void 
AggloCluster::reportCompatibleMasters() const
{
  std::cout << "--- Compatible Masters Report ---" << std::endl;
  int group_idx = 1;
  for (const auto& master_mask_pair : compatible_masters_) {
    const MasterMask& master_mask = master_mask_pair.first;
    const auto& bits_to_masters = master_mask_pair.second;

    std::cout << "Group " << group_idx++ << ": " << master_mask.to_string()
              << std::endl;

    for (const auto& bits_pair : bits_to_masters) {
      std::cout << "  " << bits_pair.first << "-bit Masters:" << std::endl;
      for (odb::dbMaster* master : bits_pair.second) {
        std::cout << "    - " << master->getName() << std::endl;
      }
    }
    std::cout << std::endl;
  }
  std::cout << "--- End of Report ---" << std::endl;
}

void 
AggloCluster::reportFlopUnits() const
{
  std::cout << "--- Flop Units Report ---" << std::endl;
  std::cout << "Total Flop Units: " << flop_units_.size() << std::endl;

  if (!flop_units_.empty()) {
    std::cout << "\nSample Flop Units (up to 5):" << std::endl;
    for (size_t i = 0; i < std::min(flop_units_.size(), (size_t) 5); ++i) {
      std::cout << flop_units_[i].to_string() << std::endl;
    }
  }
  std::cout << "\n--- End of Report ---" << std::endl;
}

void 
AggloCluster::reportCompatibleGroups() const
{
  std::cout << "--- Compatible Groups Report ---" << std::endl;
  std::cout << "Total Compatible Groups: " << compatible_groups_.size()
            << std::endl;

  int group_idx = 1;
  for (const auto& group_pair : compatible_groups_) {
    const auto& master_mask = group_pair.first.first;
    const auto& inst_mask = group_pair.first.second;
    const auto& flops = group_pair.second;

    std::cout << "\nGroup " << group_idx++ << " (" << flops.size()
              << " flops):" << std::endl;
    std::cout << "  - Master Mask: " << master_mask.to_string() << std::endl;
    std::cout << "  - Inst Mask:   " << inst_mask.to_string() << std::endl;
  }
  std::cout << "\n--- End of Report ---" << std::endl;
}

void
AggloCluster::reportTimingPaths() const
{
  log_->report("--- Timing Paths Report ---");
  log_->report("Total paths checked: {}", timing_paths_.size());

  if (!timing_paths_.empty()) {
    log_->report("\nSample Timing Paths (up to 5 worst slack):");
    int count = 0;
    for (const auto& path : timing_paths_) {
      if (count >= 5) {
        break;
      }
      const FlopUnit& start_flop = flop_units_[path.start_flop_idx_];
      const FlopUnit& end_flop = flop_units_[path.end_flop_idx_];
      log_->report("  - Path {}: Start FF: {}, End FF: {}, Slack: {:.4f}",
                   count + 1,
                   start_flop.inst_->getName(),
                   end_flop.inst_->getName(),
                   sta_->units()->timeUnit()->staToUser(path.slack_));
      count++;
    }
  }
  log_->report("\n--- End of Report ---");
}

void
AggloCluster::reportFlopPaths() const
{
  log_->report("--- Flop Timing Paths Report ---");
  log_->report("Total Flop Units: {}", flop_units_.size());

  if (!flop_units_.empty()) {
    log_->report("\nSample Flop Units (up to 5):");
    for (size_t i = 0; i < std::min(flop_units_.size(), (size_t) 5); ++i) {
      const FlopUnit& flop = flop_units_[i];
      log_->report("  - Flop: {}", flop.inst_->getName());
      log_->report("    - Total timing paths: {}", flop.timing_paths_.size());
      if (!flop.pin_budgets_.empty()) {
        log_->report("    - Pin Budgets:");
        for (const auto& pair : flop.pin_budgets_) {
          odb::dbITerm* pin = pair.first;
          const auto& budgets = pair.second;
          if (!budgets.empty()) {
            log_->report("      - Pin: {}, Num Paths: {}, Smallest Budget: {:.4f}ns",
                         pin->getMTerm()->getName(),
                         budgets.size(),
                         sta_->units()->timeUnit()->staToUser(budgets[0].second));
          }
        }
      }
    }
  }
}

void
AggloCluster::reportFeasibleRegions() const
{
  log_->report("--- Feasible Regions Report ---");
  log_->report("Total Flop Units: {}", flop_units_.size());

  if (!flop_units_.empty()) {
    log_->report("\nSample Flop Units (up to 5):");
    for (size_t i = 0; i < std::min(flop_units_.size(), (size_t) 5); ++i) {
      const FlopUnit& flop = flop_units_[i];
      log_->report("  - Flop: {}", flop.inst_->getName());
      if (boost::geometry::is_empty(flop.feasible_region_)) {
        log_->report("    - Feasible Region: Empty");
      } else {
        log_->report("    - Feasible Region: min({} {}), max({} {})",
                     flop.feasible_region_.min_corner().get<0>(),
                     flop.feasible_region_.min_corner().get<1>(),
                     flop.feasible_region_.max_corner().get<0>(),
                     flop.feasible_region_.max_corner().get<1>());
      }
    }
  }
}

void
AggloCluster::reportCompatibilityGraph() const
{
  log_->report("--- Compatibility Graph Report ---");
  log_->report("Total edges in priority queue: {}", edge_pq_.size());

  if (!edge_pq_.empty()) {
    log_->report("\nSample Edges (up to 10 with highest gain):");
    int count = 0;
    for (auto it = edge_pq_.rbegin(); it != edge_pq_.rend() && count < 10; ++it, ++count) {
      const Edge& edge = *it;
      log_->report("  - Edge {}:", count + 1);
      log_->report("    - Cluster 1 ({} flops):", flop_clusters_[edge.n1].flops_.size());
      for (int flop_idx : flop_clusters_[edge.n1].flops_) {
        log_->report("      - {}", flop_units_[flop_idx].inst_->getName());
      }
      log_->report("    - Cluster 2 ({} flops):", flop_clusters_[edge.n2].flops_.size());
      for (int flop_idx : flop_clusters_[edge.n2].flops_) {
        log_->report("      - {}", flop_units_[flop_idx].inst_->getName());
      }
      log_->report("    - Gain (HPWL reduction): {:.2f}", edge.weight);
      log_->report("    - Suggested Placement: {}", edge.pos.to_string());
    }
  }
  log_->report("\n--- End of Report ---");
}





}  // end namespace gpl
