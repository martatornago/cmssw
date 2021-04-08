#ifndef RecoMTD_DetLayers_MTDDetSector_H
#define RecoMTD_DetLayers_MTDDetSector_H

#include "TrackingTools/DetLayers/interface/GeometricSearchDet.h"
#include "DataFormats/GeometrySurface/interface/BoundDiskSector.h"

#include <ostream>

class GeomDet;

class MTDDetSector : public GeometricSearchDet {
public:
  using GeometricSearchDet::GeometricSearchDet;

  /// Construct from iterators on GeomDet*
  MTDDetSector(std::vector<const GeomDet*>::const_iterator first, std::vector<const GeomDet*>::const_iterator last);

  /// Construct from a vector of GeomDet*
  MTDDetSector(const std::vector<const GeomDet*>& dets);

  ~MTDDetSector() override{};

  // GeometricSearchDet structure

  const std::vector<const GeomDet*>& basicComponents() const override { return theDets; }

  const BoundSurface& surface() const final { return *theDiskS; }

  const std::vector<const GeometricSearchDet*>& components() const override;

  std::pair<bool, TrajectoryStateOnSurface> compatible(const TrajectoryStateOnSurface& ts,
                                                       const Propagator& prop,
                                                       const MeasurementEstimator& est) const override;

  std::vector<DetWithState> compatibleDets(const TrajectoryStateOnSurface& startingState,
                                           const Propagator& prop,
                                           const MeasurementEstimator& est) const override;

  void compatibleDetsV(const TrajectoryStateOnSurface& startingState,
                       const Propagator& prop,
                       const MeasurementEstimator& est,
                       std::vector<DetWithState>& result) const override;

  std::vector<DetGroup> groupedCompatibleDets(const TrajectoryStateOnSurface& startingState,
                                              const Propagator& prop,
                                              const MeasurementEstimator& est) const override;

  // GeometricSearchDet extension

  const BoundDiskSector& specificSurface() const { return *theDiskS; }
    
  std::pair<uint32_t, uint32_t> ModuleIndex(uint32_t DetId) const;
   
  uint32_t ShiftedModuleIndex(uint32_t, int, int) const;

  size_t compatibleDetsLine(size_t idetMin,
                          std::vector<DetWithState>& result,
                          const TrajectoryStateOnSurface& tsos,
                          const Propagator& prop,
                          const MeasurementEstimator& est,
                          GlobalPoint startPos) const;

protected:
  void setDisk(BoundDiskSector* diskS) { theDiskS = diskS; }

  bool add(size_t idet,
           std::vector<DetWithState>& result,
           const TrajectoryStateOnSurface& tsos,
           const Propagator& prop,
           const MeasurementEstimator& est) const;

private:
  ReferenceCountingPointer<BoundDiskSector> theDiskS;
  std::vector<const GeomDet*> theDets;

  // Window of detid ordered modules around that closest to the track extrapolation on the sector surface
  // needed to limit the size of the vector of distances to sort
  // value 50 based on the possible mismatch of module number between adjacent
  // modules, due to left-right type imparity

  void init();
    
  const double x_offset = 22.55;  //all constants in mm
  const double sensor_module_x = 43.1;
  const double deltaX = 0.5;
  const size_t lastModule_backLeft = 512;
  const size_t lastModule_backRight = 1025;
  const size_t lastModule_frontLeft = 516;
  const size_t lastModule_frontRight = 1032;

    
  std::vector<double> offset_FR = {x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset+5*(sensor_module_x+deltaX), x_offset+6*(sensor_module_x+deltaX), x_offset+7*(sensor_module_x+deltaX), x_offset+8*(sensor_module_x+deltaX), x_offset+8*(sensor_module_x+deltaX), x_offset+7*(sensor_module_x+deltaX), x_offset+6*(sensor_module_x+deltaX), x_offset+2*(sensor_module_x+deltaX), x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset,0};
  std::vector<double> offset_FL = {0, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset+2*(sensor_module_x+deltaX), x_offset+6*(sensor_module_x+deltaX), x_offset+7*(sensor_module_x+deltaX), x_offset+8*(sensor_module_x+deltaX), x_offset+8*(sensor_module_x+deltaX), x_offset+7*(sensor_module_x+deltaX), x_offset+6*(sensor_module_x+deltaX), x_offset+5*(sensor_module_x+deltaX), x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, 0};
  std::vector<double> offset_BR = {0, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset+3*(sensor_module_x+deltaX), x_offset+6*(sensor_module_x+deltaX), x_offset+7*(sensor_module_x+deltaX), x_offset+7*(sensor_module_x+deltaX), x_offset+8*(sensor_module_x+deltaX), x_offset+7*(sensor_module_x+deltaX), x_offset+6*(sensor_module_x+deltaX), x_offset+4*(sensor_module_x+deltaX), x_offset,x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset,0};
  std::vector<double> offset_BL = {0, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset+4*(sensor_module_x+deltaX), x_offset+6*(sensor_module_x+deltaX), x_offset+7*(sensor_module_x+deltaX), x_offset+8*(sensor_module_x+deltaX), x_offset+7*(sensor_module_x+deltaX), x_offset+7*(sensor_module_x+deltaX), x_offset+6*(sensor_module_x+deltaX), x_offset+3*(sensor_module_x+deltaX), x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset, x_offset,0};
    
  std::vector<uint32_t> start_copy_FR = {1, 7, 18, 33, 50, 69, 90, 112, 136, 161, 186, 207, 227, 247, 266, 285, 305, 325, 349, 374, 398, 421, 443, 463, 481, 497, 510};
  std::vector<uint32_t> start_copy_FL = {1, 8, 21, 37, 55, 75, 97, 120, 144, 169, 193, 213, 233, 252, 271, 291, 311, 332, 357, 382, 406, 428, 449, 468, 485, 500, 511};
  std::vector<uint32_t> start_copy_BR = {1, 10, 23, 39, 57, 77, 99, 122, 146, 171, 194, 214, 234, 254, 273, 293, 313, 335, 360, 384, 407, 430, 451, 470, 487, 501, 511};
  std::vector<uint32_t> start_copy_BL = {1, 4, 14, 28, 45, 64, 85, 107, 130, 154, 179, 201, 221, 241, 260, 280, 300, 320, 343, 368, 392, 415, 437, 457, 475, 491, 504};
    
};

std::ostream& operator<<(std::ostream&, const MTDDetSector&);

#endif
