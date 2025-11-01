#include "aggloCluster.h"

#include <algorithm>
#include <string>
#include <iostream>
#include <vector>

#include "sta/Sta.hh"
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
#include "sta/Units.hh"
#include "utl/Logger.h"

namespace gpl {

void AggloCluster::createCompatibleGroups()
{
  compatible_groups_.clear();

  for (size_t i = 0; i < flop_units_.size(); ++i) {
    const FlopUnit& unit = flop_units_[i];
    compatible_groups_[{unit.master_mask_, unit.inst_mask_}].push_back(i);
  }
}

InstMask AggloCluster::createInstMask(odb::dbInst* inst)
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

MasterMask AggloCluster::createMasterMask(odb::dbInst* inst)
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

int AggloCluster::getFuncId(const sta::FuncExpr* expr, odb::dbInst* inst)
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

std::string AggloCluster::getFuncStr(const sta::FuncExpr* expr,
                                     odb::dbInst* inst) const
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

const sta::LibertyCell* AggloCluster::getLibertyCell(odb::dbInst* inst) const
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

const sta::LibertyPort* AggloCluster::getLibertyPort(odb::dbITerm* iterm) const
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

const sta::FuncExpr* AggloCluster::getPortFunc(
    const sta::LibertyPort* port) const
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

int AggloCluster::getNumDPins(odb::dbInst* inst) const
{
  int num = 0;
  for (odb::dbITerm* iterm : inst->getITerms()) {
    num += isDPin(iterm);
  }
  return num;
}

int AggloCluster::getNumQPins(odb::dbInst* inst) const
{
  int num = 0;
  for (odb::dbITerm* iterm : inst->getITerms()) {
    num += isQPin(iterm);
  }
  return num;
}

int AggloCluster::getNumQNPins(odb::dbInst* inst) const
{
  int num = 0;
  for (odb::dbITerm* iterm : inst->getITerms()) {
    num += isQNPin(iterm);
  }
  return num;
}

FlopPort AggloCluster::getPortType(const sta::LibertyPort* lib_port,
                                   odb::dbInst* inst) const
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

bool AggloCluster::hasClear(odb::dbInst* inst) const
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

bool AggloCluster::hasPositiveClockEdge(odb::dbInst* inst) const
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

bool AggloCluster::hasPreset(odb::dbInst* inst) const
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

bool AggloCluster::hasScan(odb::dbInst* inst) const
{
  const sta::LibertyCell* lib_cell = getLibertyCell(inst);
  if (lib_cell == nullptr) {
    return false;
  }
  return lib_cell && getLibertyScanIn(lib_cell)
         && getLibertyScanEnable(lib_cell);
}

bool AggloCluster::isClearPin(odb::dbITerm* iterm) const
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

bool AggloCluster::isClockPin(odb::dbITerm* iterm) const
{
  return iterm->getSigType() == odb::dbSigType::CLOCK
         || network_->isRegClkPin(network_->dbToSta(iterm));
}

bool AggloCluster::isDPin(odb::dbITerm* iterm) const
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

bool AggloCluster::isPowerPin(odb::dbITerm* iterm) const
{
  return iterm->getSigType()
      .isSupply();  // dbSigType::isSupply returns true for POWER or GROUND
}

bool AggloCluster::isPresetPin(odb::dbITerm* iterm) const
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

bool AggloCluster::isScanInPin(odb::dbITerm* iterm) const
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

bool AggloCluster::isScanEnablePin(odb::dbITerm* iterm) const
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

bool AggloCluster::isQPin(odb::dbITerm* iterm) const
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

bool AggloCluster::isQNPin(odb::dbITerm* iterm) const
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

bool AggloCluster::isValidFlop(odb::dbInst* inst) const
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
      log_->warn(
          utl::GPL,
          9991,
          "Instance {} / Master {} has an unrecognized pin {} - considering invalid.",
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

void AggloCluster::readCompatibleMasters()
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

  std::copy_if(temp_master_map.begin(),
               temp_master_map.end(),
               std::inserter(compatible_masters_, compatible_masters_.end()),
               [](const auto& pair) { return pair.second.size() >= 2; });
}

void AggloCluster::readFlopUnits()
{
  flop_units_.clear();
  inst_to_flop_id_.clear();

  int flop_id = 0;
  for (odb::dbInst* inst : block_->getInsts()) {
    if (!isValidFlop(inst)) {
      continue;
    }

    flop_units_.emplace_back(inst,
                             inst->getBBox()->getBox(),
                             createMasterMask(inst),
                             createInstMask(inst));
    inst_to_flop_id_[inst] = flop_id++;
  }
}


void AggloCluster::readTimingPaths()
{
  timing_paths_.clear();
  paths_by_start_flop_id_.clear();
  paths_by_end_flop_id_.clear();

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
    sta::Vertex* sta_start_vertex = start_path ? start_path->vertex(sta_) : nullptr;
    sta::Pin* sta_start_pin = sta_start_vertex ? sta_start_vertex->pin() : nullptr;
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

    timing_paths_.emplace_back(start_flop_id, end_flop_id, slack, path);
    int path_idx = timing_paths_.size() - 1;

    paths_by_start_flop_id_[start_flop_id].push_back(path_idx);
    paths_by_end_flop_id_[end_flop_id].push_back(path_idx);
  }
}

void AggloCluster::reportCompatibleMasters() const
{
  std::cout << "--- Compatible Masters Report ---" << std::endl;
  int group_idx = 1;
  for (const auto& master_mask_pair : compatible_masters_) {
    const MasterMask& master_mask = master_mask_pair.first;
    const auto& bits_to_masters = master_mask_pair.second;

    std::cout << "Group " << group_idx++ << ": " << master_mask.to_string() << std::endl;

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

void AggloCluster::reportFlopUnits() const
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

void AggloCluster::reportCompatibleGroups() const
{
  std::cout << "--- Compatible Groups Report ---" << std::endl;
  std::cout << "Total Compatible Groups: " << compatible_groups_.size()
            << std::endl;

  int group_idx = 1;
  for (const auto& group_pair : compatible_groups_) {
    const auto& master_mask = group_pair.first.first;
    const auto& inst_mask = group_pair.first.second;
    const auto& flops = group_pair.second;

    std::cout << "\nGroup " << group_idx++ << " (" << flops.size() << " flops):" << std::endl;
    std::cout << "  - Master Mask: " << master_mask.to_string() << std::endl;
    std::cout << "  - Inst Mask:   " << inst_mask.to_string() << std::endl;
  }
  std::cout << "\n--- End of Report ---" << std::endl;
}

void AggloCluster::reportTimingPaths() const
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

}

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
      threads_(threads),
      multiplier_(block_->getDbUnitsPerMicron())
{
}

AggloCluster::~AggloCluster() = default;

}  // end namespace gpl
