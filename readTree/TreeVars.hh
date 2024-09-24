#ifndef _TREEVARS_HH_
#define _TREEVARS_HH_
#include <vector>
int EventID, NHits;
std::vector<int> *PDGID = 0, *TrackID = 0, *ParentTrackID = 0, *DetID = 0, *ProcessID = 0;
std::vector<double> *HitInPositionX = 0, *HitInPositionY = 0, *HitInPositionZ = 0;
std::vector<double> *HitInMomentumX = 0, *HitInMomentumY = 0, *HitInMomentumZ = 0;
std::vector<double> *VertexPositionX = 0, *VertexPositionY = 0, *VertexPositionZ = 0;
std::vector<double> *HitOutPositionX = 0, *HitOutPositionY = 0, *HitOutPositionZ = 0;
std::vector<double> *HitGlobalTime = 0, *HitEdep = 0, *TrackLength = 0;
int TOF_NHits;
std::vector<int> *TOF_PDGID = 0, *TOF_TrackID = 0, *TOF_ParentTrackID = 0, *TOF_DetID = 0, *TOF_ProcessID = 0;
std::vector<double> *TOF_HitInPositionX = 0, *TOF_HitInPositionY = 0, *TOF_HitInPositionZ = 0;
std::vector<double> *TOF_HitInMomentumX = 0, *TOF_HitInMomentumY = 0, *TOF_HitInMomentumZ = 0;
std::vector<double> *TOF_VertexPositionX = 0, *TOF_VertexPositionY = 0, *TOF_VertexPositionZ = 0;
std::vector<double> *TOF_HitOutPositionX = 0, *TOF_HitOutPositionY = 0, *TOF_HitOutPositionZ = 0;
std::vector<double> *TOF_HitGlobalTime = 0, *TOF_HitEdep = 0, *TOF_TrackLength = 0;
double TOF_TotalEdep, TOF_TotalTrackLength;
std::vector<double> *TOF_ModuleEdep = 0, *TOF_ModuleTrackLength = 0;
std::vector<int> *TOF_DID = 0;
#endif