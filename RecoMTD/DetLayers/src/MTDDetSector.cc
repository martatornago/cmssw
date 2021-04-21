#define EDM_ML_DEBUG

#include "RecoMTD/DetLayers/interface/MTDDetSector.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/DetLayers/interface/MeasurementEstimator.h"
#include "MTDDiskSectorBuilderFromDet.h"
#include "DataFormats/GeometrySurface/interface/RectangularPlaneBounds.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

MTDDetSector::MTDDetSector(vector<const GeomDet*>::const_iterator first, vector<const GeomDet*>::const_iterator last)
    : GeometricSearchDet(false), theDets(first, last) {
  init();
}

MTDDetSector::MTDDetSector(const vector<const GeomDet*>& vdets) : GeometricSearchDet(false), theDets(vdets) { init(); }

void MTDDetSector::init() {
  // Add here the sector build based on a collection of GeomDets, mimic what done in ForwardDetRingOneZ
  // using the code from tracker BladeShapeBuilderFromDet
  // simple initial version, no sorting for the time being
  setDisk(MTDDiskSectorBuilderFromDet()(theDets));
}

const vector<const GeometricSearchDet*>& MTDDetSector::components() const {
  // FIXME dummy impl.
  edm::LogError("MTDDetLayers") << "temporary dummy implementation of MTDDetSector::components()!!";
  static const vector<const GeometricSearchDet*> result;
  return result;
}

pair<bool, TrajectoryStateOnSurface> MTDDetSector::compatible(const TrajectoryStateOnSurface& ts,
                                                              const Propagator& prop,
                                                              const MeasurementEstimator& est) const {
  TrajectoryStateOnSurface ms = prop.propagate(ts, specificSurface());

#ifdef EDM_ML_DEBUG
  LogTrace("MTDDetLayers") << "MTDDetSector::compatible, sector: \n"
                           << (*this) << "\n  TS at Z,R,phi: " << std::fixed << std::setw(14) << ts.globalPosition().z()
                           << " , " << std::setw(14) << ts.globalPosition().perp() << " , " << std::setw(14)
                           << ts.globalPosition().phi();
  if (ms.isValid()) {
    LogTrace("MTDDetLayers") << " DEST at Z,R,phi: " << std::fixed << std::setw(14) << ms.globalPosition().z() << " , "
                             << std::setw(14) << ms.globalPosition().perp() << " , " << std::setw(14)
                             << ms.globalPosition().phi() << " local Z: " << std::setw(14) << ms.localPosition().z();
  } else {
    LogTrace("MTDDetLayers") << " DEST: not valid";
  }
#endif

  return make_pair(ms.isValid() and est.estimate(ms, specificSurface()) != 0, ms);
}

vector<GeometricSearchDet::DetWithState> MTDDetSector::compatibleDets(const TrajectoryStateOnSurface& startingState,
                                                                      const Propagator& prop,
                                                                      const MeasurementEstimator& est) const {
  LogTrace("MTDDetLayers") << "MTDDetSector::compatibleDets, sector: \n"
                           << (*this) << "\n  TS at Z,R,phi: " << std::fixed << std::setw(14)
                           << startingState.globalPosition().z() << " , " << std::setw(14)
                           << startingState.globalPosition().perp() << " , " << std::setw(14)
                           << startingState.globalPosition().phi();

  vector<DetWithState> result;

  // Propagate and check that the result is within bounds
  pair<bool, TrajectoryStateOnSurface> compat = compatible(startingState, prop, est);
  if (!compat.first) {
    LogTrace("MTDDetLayers") << "    MTDDetSector::compatibleDets: not compatible"
                             << "    (should not have been selected!)";
    return result;
  }

  TrajectoryStateOnSurface& tsos = compat.second;
  GlobalPoint startPos = tsos.globalPosition();

  LogTrace("MTDDetLayers") << "Starting position: " << startPos << " starting p/pT: " << tsos.globalMomentum().mag()
                           << " / " << tsos.globalMomentum().perp();

  // determine distance of det center from extrapolation on the surface, sort dets accordingly

  size_t idetMin = basicComponents().size();
  double dist2Min = std::numeric_limits<double>::max();
  std::vector<std::pair<double, size_t> > tmpDets;
  tmpDets.reserve(basicComponents().size());

  for (size_t idet = 0; idet < basicComponents().size(); idet++) {
    double dist2 = (startPos - theDets[idet]->position()).mag2();
    tmpDets.emplace_back(dist2, idet);
    if (dist2 < dist2Min) {
      dist2Min = dist2;
      idetMin = idet;
    }
  }

  //look for the compatibledets considering each line of the sector

  if (add(idetMin, result, tsos, prop, est)) {
    compatibleDetsLine(idetMin, result, tsos, prop, est, startPos);

    for (int iside = -1; iside <= 1; iside += 2) {
      size_t shift(1);
      bool isCompatible(true);
      size_t idetNew(theDets.size());

      while (isCompatible) {
        idetNew = ShiftedModuleIndex(theDets[idetMin]->geographicalId().rawId(), 0, iside * shift);
        if (idetNew == theDets.size()) {
          break;
        }
        isCompatible = add(idetNew, result, tsos, prop, est);
        if (isCompatible) {
          compatibleDetsLine(idetNew, result, tsos, prop, est, startPos);
          shift++;
        }
      }
    }
  }

#ifdef EDM_ML_DEBUG
  if (result.empty()) {
    LogTrace("MTDDetLayers") << "MTDDetSector::compatibleDets, closest not compatible!";
  } else {
    LogTrace("MTDDetLayers") << "MTDDetSector::compatibleDets, found " << result.size() << " compatible dets";
  }
#endif

  return result;
}

void MTDDetSector::compatibleDetsV(const TrajectoryStateOnSurface&,
                                   const Propagator&,
                                   const MeasurementEstimator&,
                                   std::vector<DetWithState>&) const {
  edm::LogError("MTDDetLayers") << "At the moment not a real implementation";
}

vector<DetGroup> MTDDetSector::groupedCompatibleDets(const TrajectoryStateOnSurface& startingState,
                                                     const Propagator& prop,
                                                     const MeasurementEstimator& est) const {
  // FIXME should be implemented to allow returning  overlapping chambers
  // as separate groups!
  edm::LogInfo("MTDDetLayers") << "dummy implementation of MTDDetSector::groupedCompatibleDets()";
  vector<DetGroup> result;
  return result;
}

bool MTDDetSector::add(size_t idet,
                       vector<DetWithState>& result,
                       const TrajectoryStateOnSurface& tsos,
                       const Propagator& prop,
                       const MeasurementEstimator& est) const {
  pair<bool, TrajectoryStateOnSurface> compat = theCompatibilityChecker.isCompatible(theDets[idet], tsos, prop, est);

  if (compat.first) {
    result.push_back(DetWithState(theDets[idet], compat.second));
    LogTrace("MTDDetLayers") << "MTDDetSector::compatibleDets found compatible det idetMin " << idet
                             << " detId = " << theDets[idet]->geographicalId().rawId() << " at "
                             << theDets[idet]->position()
                             << " dist = " << std::sqrt((tsos.globalPosition() - theDets[idet]->position()).mag2());
  }

  return compat.first;
}

std::ostream& operator<<(std::ostream& os, const MTDDetSector& id) {
  os << " MTDDetSector at " << std::fixed << id.specificSurface().position() << std::endl
     << " L/W/T   : " << std::setw(14) << id.specificSurface().bounds().length() << " / " << std::setw(14)
     << id.specificSurface().bounds().width() << " / " << std::setw(14) << id.specificSurface().bounds().thickness()
     << std::endl
     << " rmin    : " << std::setw(14) << id.specificSurface().innerRadius() << std::endl
     << " rmax    : " << std::setw(14) << id.specificSurface().outerRadius() << std::endl
     << " phi ref : " << std::setw(14) << id.specificSurface().position().phi() << std::endl
     << " phi w/2 : " << std::setw(14) << id.specificSurface().phiHalfExtension() << std::endl;
  return os;
}

std::pair<uint32_t, uint32_t> MTDDetSector::ModuleIndex(uint32_t DetId) const {
  ETLDetId start_mod(DetId);

  uint32_t c = start_mod.module();
  uint32_t sensor_type = start_mod.modType();  //1=left; 2=right;
  uint32_t disc_side = start_mod.discSide();   //0=front; 1=back;
  uint32_t column, row, global_row;
  row = halfRowRange - 1;

  std::array<uint32_t, halfRowRange> start_copy;

  if (disc_side == 0) {
    if (sensor_type == 2) {
      start_copy = start_copy_FR;
    } else {
      start_copy = start_copy_FL;
    }
  } else {
    if (sensor_type == 2) {
      start_copy = start_copy_BR;
    } else {
      start_copy = start_copy_BL;
    }
  }

  for (size_t i = 0; i < halfRowRange - 1; i++) {
    if (start_copy[i] <= c && start_copy[i + 1] > c) {
      row = i;
      break;
    }
  }
  column = c - start_copy[row];

  global_row = row * 2;
  if ((disc_side == 0 && sensor_type == 1) || (disc_side == 1 && sensor_type == 2)) {
    global_row++;
  }

  return std::make_pair(column, global_row);
}

size_t MTDDetSector::ShiftedModuleIndex(uint32_t DetId, int horizontalShift, int verticalShift) const {
  ETLDetId start_mod(DetId);

  uint32_t module = start_mod.module();
  uint32_t modtyp = start_mod.modType();
  uint32_t discside = start_mod.discSide();
  std::pair<uint32_t, uint32_t> pair = ModuleIndex(DetId);
  int row_init, row, global_row_init, global_row, column;
  double offset_init;
  column = pair.first + horizontalShift;
  global_row_init = pair.second;
  global_row = global_row_init + verticalShift;

  if (global_row < 0 || global_row >= static_cast<int>(halfRowRange) * 2) {
    return theDets.size();
  }

  std::array<uint32_t, halfRowRange> start_copy;
  std::array<double, halfRowRange> offset;
  size_t lastModule;

  if (discside == 0) {
    if (global_row_init % 2 == 0) {
      row_init = global_row_init / 2;  //front right
      offset_init = offset_FR[row_init];
    } else {
      row_init = (global_row_init - 1) / 2;
      offset_init = offset_FL[row_init];
    }
    if (global_row % 2 == 0) {
      row = global_row / 2;  //front right
      modtyp = 2;
      start_copy = start_copy_FR;
      offset = offset_FR;
      lastModule = lastModule_frontRight;
    } else {
      row = (global_row - 1) / 2;
      modtyp = 1;
      start_copy = start_copy_FL;
      offset = offset_FL;
      lastModule = lastModule_frontLeft;
      
    }
  } else {
    if (global_row_init % 2 == 0) {
      row_init = global_row_init / 2;  //back left
      offset_init = offset_BL[row_init];
    } else {
      row_init = (global_row_init - 1) / 2;
      offset_init = offset_BR[row_init];
    }
    if (global_row % 2 == 0) {
      row = global_row / 2;  //back left
      modtyp = 1;
      start_copy = start_copy_BL;
      offset = offset_BL;
      lastModule = lastModule_backLeft;
    } else {
      row = (global_row - 1) / 2;
      modtyp = 2;
      start_copy = start_copy_BR;
      offset = offset_BR;
      lastModule = lastModule_backRight;
    }
  }

  if (row >= static_cast<int>(halfRowRange) || row < 0) {
    return theDets.size();
  }

  module = start_copy[row] + column;
  if (offset_init != offset[row]) {
    module += (offset_init - offset[row]) / (sensor_module_x + deltaX);
  }
    if (row == halfRowRange - 1) {
        if (module > lastModule || module < start_copy[row]) {
          return theDets.size();
        }
    } else if (module >= start_copy[row + 1] || module < start_copy[row]) {
    return theDets.size();
  }

  size_t id;
  id = (modtyp == 1) ? module : module + lastModule;

  if ((discside == 0 && module > lastModule_frontLeft + lastModule_frontRight) ||
      (discside == 1 && module > lastModule_backLeft + lastModule_backRight)) {
    return theDets.size();
  }

  return id - 1;
}

void MTDDetSector::compatibleDetsLine(const size_t idetMin,
                                      vector<DetWithState>& result,
                                      const TrajectoryStateOnSurface& tsos,
                                      const Propagator& prop,
                                      const MeasurementEstimator& est,
                                      GlobalPoint startPos) const {
  for (int iside = -1; iside <= 1; iside += 2) {
    size_t shift(1);
    bool isCompatible(true);
    size_t idetNew(theDets.size());

    while (isCompatible) {
      idetNew = ShiftedModuleIndex(theDets[idetMin]->geographicalId().rawId(), iside * shift, 0);
      if (idetNew == theDets.size()) {
        break;
      }
      isCompatible = add(idetNew, result, tsos, prop, est);
      shift++;
    }
  }

  return;
}
