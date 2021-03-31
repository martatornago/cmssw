//#define EDM_ML_DEBUG

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


std::pair<uint32_t, uint32_t> MTDDetSector::ModuleIndex(uint32_t DetId) const {

 ETLDetId start_mod(DetId);

 uint32_t c = start_mod.module();
 uint32_t sensor_type = start_mod.modType();  //1=left; 2=right;
 uint32_t disc_side = start_mod.discSide(); //0=front; 1=back;
 uint32_t column, row, global_row;
 row=0;
    
 if(disc_side == 0){
  if( sensor_type == 2 ){
   for( int i=0; i<27; i++ ){
       if( start_copy_FR[i]<=c && start_copy_FR[i+1]>c ){
           row = i;
       }
   }
   column = c - start_copy_FR[row];
  }else{
   for( int i=0; i<27; i++ ){
       if( start_copy_FL[i]<=c && start_copy_FL[i+1]>c ){
           row = i;
       }
   }
   column = c - start_copy_FL[row];
  }
 }else{
  if( sensor_type == 2 ){
   for( int i=0; i<27; i++ ){
       if( start_copy_BR[i]<=c && start_copy_BR[i+1]>c ){
           row = i;
       }
   }
   column = c - start_copy_BR[row];
  }else{
   for( int i=0; i<27; i++ ){
       if( start_copy_BL[i]<=c && start_copy_BL[i+1]>c ){
           row = i;
       }
   }
   column = c - start_copy_BL[row];
  }
 }
    
     
 if(disc_side == 0){
     if(sensor_type == 1){
         global_row = row*2+1;
     }else{
         global_row = row*2;
     }
 }else{
     if(sensor_type == 1){
         global_row = row*2;
     }else{
         global_row = row*2+1;
     }
 }

    
    
    return std::make_pair (column, global_row);
    

}

uint32_t MTDDetSector::ShiftedModuleIndex(uint32_t DetId, int horizontalShift, int verticalShift) const {
  
    ETLDetId start_mod(DetId);
    
    uint32_t module = start_mod.module();
    uint32_t modtyp = start_mod.modType();
    uint32_t discside = start_mod.discSide();
    
    std::pair<uint32_t, uint32_t> pair = ModuleIndex(DetId);
    
    uint32_t row_init, row, global_row_init, global_row, column;

    if(horizontalShift<0 && abs(horizontalShift)>pair.first){
        return 0;
    }else{
        column = pair.first + horizontalShift;
    }

    
    double offset_init;
    
    global_row_init = pair.second;
    global_row = global_row_init + verticalShift;
    
    if(discside == 0){
        if(global_row_init%2==0){
            row_init = global_row_init/2; //front right
            offset_init = offset_FR[row_init];
        }else{
            row_init = (global_row_init-1)/2;
            offset_init = offset_FL[row_init+1];
        }
        if(global_row%2==0){
            row=global_row/2; //front right
            modtyp = 2;
        }else{
            row=(global_row-1)/2;
            modtyp = 1;
        }
    }else{
        if(global_row_init%2==0){
            row_init = global_row_init/2; //back left
            offset_init = offset_BL[row_init+1];
        }else{
            row_init = (global_row_init-1)/2;
            offset_init = offset_BR[row_init+1];
        }
        if(global_row%2==0){
            row=global_row/2; //back left
            modtyp = 1;
        }else{
            row=(global_row-1)/2;
            modtyp = 2;
        }
    }
    
    if(discside==0){
        if(modtyp==1){
            if(row>=27) return 0;
            else{
                module = start_copy_FL[row] + column;
                if(offset_init!=offset_FL[row+1]){
                    module += (offset_init - offset_FL[row+1])/(sensor_module_x + deltaX);
                }
                if(module>=start_copy_FL[row+1]) return 0;
            }
        }else{
            if(row>=27) return 0;
            else{
                module = start_copy_FR[row] + column;
                if(offset_init!=offset_FR[row]){
                    module += (offset_init - offset_FR[row])/(sensor_module_x + deltaX);
                }
                if(module>=start_copy_FR[row+1]) return 0;
            }
        }
    }else{
        if(modtyp==1){
            if(row>=27) return 0;
            else{
                module = start_copy_BL[row] + column;
                if(offset_init!=offset_BL[row+1]){
                    module += (offset_init - offset_BL[row+1])/(sensor_module_x + deltaX);
                }
                if(module>=start_copy_BL[row+1]) return 0;
            }
        }else{
            if(row>=27) return 0;
            else{
                module = start_copy_BR[row] + column;
                if(offset_init!=offset_BR[row+1]){
                    module += (offset_init - offset_BR[row+1])/(sensor_module_x + deltaX);
                }
                if(module>=start_copy_BR[row+1]) return 0;
            }
        }
    }

    uint32_t id;
    
    if(modtyp==1){
        id = module;
    }else{
        if(discside==0){
            id = module+517;
        }else{
            id = module+514;
            }
        }
    
    return id;
    
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

  LogTrace("MTDDetLayers") << "Starting position: " << startPos;

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
    LogTrace("MTDDetLayers") << "MTDDetSector::compatibleDets " << std::fixed << std::setw(8) << idet << " "
                             << theDets[idet]->geographicalId().rawId() << " dist = " << std::setw(10)
                             << std::sqrt(dist2) << " Min idet/dist = " << std::setw(8) << idetMin << " "
                             << std::setw(10) << std::sqrt(dist2Min) << " " << theDets[idet]->position();
  }

  // loop on an interval od ordered detIds around the minimum
  // set a range of GeomDets around the minimum compatible with the geometry of ETL

    std::vector<size_t> detsMatrix;

//look for the compatibledets inside a matrix with matrixSizexmatrixSize dimension
    for(int vertShift = -matrixSize; vertShift<int(matrixSize+1); vertShift++){
        for(int horShift = -matrixSize; horShift<int(matrixSize+1); horShift++){
            if(vertShift!=0 || horShift!=0){
                detsMatrix.push_back(ShiftedModuleIndex(theDets[idetMin]->geographicalId().rawId(), horShift, vertShift));
            }
        }
    }
    
    for (size_t det=0; det<detsMatrix.size(); det++) {
      if (add(detsMatrix[det], result, tsos, prop, est)) {
        LogTrace("MTDDetLayers") << "MTDDetSector::compatibleDets found compatible det " << detsMatrix[det]
                                 << " detId = " << theDets[detsMatrix[det]]->geographicalId().rawId() << " at "
                                 << theDets[detsMatrix[det]]->position() << " dist = " << std::sqrt(tmpDets[detsMatrix[det]].first);
      }
    }
    
//look for the compatibledets considering each line of the sector
    
//    size_t upShift=1;
//    size_t downShift=-1;
//    bool isCompatibleUp=1;
//    bool isCompatibleDown=1;
//
//    size_t idetMinNew = 0;
//    size_t idetMinUp = 0;
//    size_t idetMinDown = 0;
//
//    compatibleDetsLine(idetMin, detsMatrix, result, tsos, prop, est, idetMinNew);
//
//    while(isCompatibleUp){
//        idetMinUp = ShiftedModuleIndex(idetMinNew,0,upShift);
//        if(add(idetMinUp, result, tsos, prop, est)){
//            detsMatrix.push_back(idetMinUp);
//            compatibleDetsLine(idetMinUp, detsMatrix, result, tsos, prop, est, idetMinNew);
//            upShift ++;
//        }else{
//            isCompatibleUp = 0;
//        }
//      }
//    while(isCompatibleDown){
//        idetMinDown = ShiftedModuleIndex(idetMinNew,0,downShift);
//        if(add(idetMinDown, result, tsos, prop, est)){
//            detsMatrix.push_back(idetMinDown);
//            compatibleDetsLine(idetMinDown, detsMatrix, result, tsos, prop, est, idetMinNew);
//            downShift +=-1;
//        }else{
//            isCompatibleDown = 0;
//        }
//      }
    
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

//void MTDDetSector::compatibleDetsLine(size_t idetMin,
//                                      std::vector<size_t> detsMatrix,
//                                      vector<DetWithState>& result,
//                                      const TrajectoryStateOnSurface& tsos,
//                                      const Propagator& prop,
//                                      const MeasurementEstimator& est,
//                                      size_t idetMinNew) const {
//    size_t negShift=-1;
//    size_t posShift=1;
//    size_t maxNegidet=0;
//    size_t maxPosidet=0;
//    bool isCompatiblePos=1;
//    bool isCompatibleNeg=1;
//
//    while(isCompatiblePos){
//        if (add(idetMin+posShift, result, tsos, prop, est)){
//            detsMatrix.push_back(idetMin+posShift);
//            maxPosidet = idetMin+posShift;
//            posShift ++;
//        }else{
//            isCompatiblePos = 0;
//        }
//    }
//    
//    while(isCompatibleNeg){
//        if (add(idetMin+negShift, result, tsos, prop, est)){
//            detsMatrix.push_back(idetMin+negShift);
//            maxNegidet = idetMin+negShift;
//            negShift +=-1;
//        }else{
//            isCompatibleNeg = 0;
//        }
//    }
//    
//    idetMinNew = size_t((maxPosidet-maxNegidet)/2);
//
//    
//}
