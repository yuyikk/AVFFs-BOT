#include "TFile.h"
#include "TString.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TSystem.h"
// #include "/opt/anaconda3/envs/hep/macros/SetDefaultStyle.cxx"

#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "TreeVars.hh"
#include "Histograms.hh"
// #include "../include/Constants.hh"
using namespace std;
const int DetIdVrt1 = 10001;
const int DetIdVrt2 = 10002;
const double scaleFactor = 10.0;
constexpr double kFieldSizeX = 34;
constexpr double kFieldSizeY = 142;
constexpr double kFieldSizeZ = 100;
constexpr double timescale = 16; // 16ns for 10e7 electrons with electron beam 100µA
double EdepCut[10] = {0, 0.1, 0.5, 1, 2, 5, 10, 15, 20, 30};
int pid_list[8] = {11, 22, -11, 2112, 2212, 211, -211, 111};
double pid_mass[8] = {0.511, 0, 0.511, 939.565, 938.272, 139.57, 139.57, 134.98};
TString pid_name[8] = {"e^{-}", "#gamma", "e^{+}", "neutron", "proton", "#pi^{+}", "#pi^{-}", "#pi^{0}"};
TVector3 targetPosition(0, 0, -850);
TVector3 NeutronArm_ZAxis(1, 0, 1);
void Binning(const int &nbins, double *binning, const double &min, const double &max);
int getIdx(int x, int *vec, int size = 8);
TTree *MakeTree(TFile *);
void MakeChain(TChain *);
void InitHist();
void ActivateBranch(TTree *);
void ActivateBranch(TChain *);
void Rotation(double *vec);
double getEk(const double &mom, const int &pid_idx);
int getEkIdx(const double &Ek);
void VectorToArray(std::vector<double> &vec, double *p);
int main(int argc, char **argv)
{
    // Binning(100, 1e-6, 2.5);
    // Open the ROOT file
    // TFile *file = TFile::Open(argv[1]);
    // if (!file || file->IsZombie())
    // {
    //     std::cerr << "Error opening file!" << std::endl;
    //     return -1;
    // }
    TFile *ofile = new TFile(argv[2], "recreate");
    InitHist();
    // Get the tree
    // TTree *tree = MakeTree(file);
    TChain *tree = new TChain("T", "BOT_Chain"); //(TTree *)file->Get("T");
    tree->Add(argv[1]);
    MakeChain(tree);
    if (!tree)
    {
        cerr << "Making tree failed." << endl;
        return -1;
    }
    // Loop over entries

    TVector3 hitMom;
    Long64_t nEntries = tree->GetEntries();
    Long64_t percentage = nEntries / 50;
    cout << "Total events: " << nEntries << endl;
    TString per = "%";
    for (Long64_t i = 0; i < nEntries; i++)
    {
        if (i % percentage == 0)
        {
            cout << Form("%.1f%s", i * 100. / nEntries, per.Data()) << endl;
        }
        tree->GetEntry(i);
        int n = HitInPositionX->size();
        // Loop over hits on virtual detectors
        for (int ihit = 0; ihit < n; ++ihit)
        {
            hitMom.SetXYZ(HitInMomentumX->at(ihit), HitInMomentumY->at(ihit), HitInMomentumZ->at(ihit));
            bool isBackward = NeutronArm_ZAxis.Dot(hitMom) <= 0;
            if (isBackward)
            {
                continue;
            }
            int id = DetID->at(ihit);
            int pid = PDGID->at(ihit);
            int trkId = TrackID->at(ihit);
            int pid_idx = getIdx(pid, pid_list);
            // The z-axis is choosen to be the axis of neutron arm
            double xz[2] = {HitInPositionX->at(ihit) / 10, HitInPositionZ->at(ihit) / 10 + 850};
            double Vxz[2] = {VertexPositionX->at(ihit) / 10, VertexPositionZ->at(ihit) / 10 + 850};
            Rotation(xz);
            Rotation(Vxz);
            double hitX = xz[0];
            double hitY = HitInPositionY->at(ihit) / 10;
            double hitZ = xz[1] - 850;
            double hitVx = Vxz[0];
            double hitVy = VertexPositionY->at(ihit) / 10;
            double hitVz = Vxz[1] - 850;
            double Ek = getEk(hitMom.Mag(), pid_idx);
            int Ek_idx_max = getEkIdx(Ek);
            bool inField = (fabs(hitX) < 10) && (fabs(hitY) < 42);
            bool goodPidIdx = (pid_idx < 8 && pid_idx >= 0);
            bool goodEkIdx = (Ek_idx_max >= 0 && Ek_idx_max < 10);
            for (int idx = 0; idx <= Ek_idx_max; ++idx)
            {
                if (id == DetIdVrt1)
                {
                    if (!inField)
                    {
                        continue;
                    }
                    if (goodPidIdx && goodEkIdx)
                    {
                        h_hitXY_Vrt1[pid_idx][idx]->Fill(scaleFactor * hitX, scaleFactor * hitY);
                        h_hitZ_Vrt1[pid_idx][idx]->Fill(hitZ);
                        h_hitVx_Vrt1[pid_idx][idx]->Fill(hitVx);
                        h_hitVy_Vrt1[pid_idx][idx]->Fill(hitVy);
                        h_hitVz_Vrt1[pid_idx][idx]->Fill(hitVz);
                        h_hitVxVy_Vrt1[pid_idx][idx]->Fill(hitVx, hitVy);
                        h_hitVxVz_Vrt1[pid_idx][idx]->Fill(hitVx, hitVz);
                        h_hitVyVz_Vrt1[pid_idx][idx]->Fill(hitVy, hitVz);
                    }
                    // if ((hitZ - 850) > -680)
                    // {
                    //     cout << "Mom: " << hitMom.X() << "\t" << hitMom.Y() << "\t" << hitMom.Z() << endl;
                    //     cout << "Pos: " << hitVx << "\t" << hitVy << "\t" << hitVz << endl;
                    //     cout << NeutronArm_ZAxis.Dot(hitMom) << "\t" << isBackward << endl;
                    // }
                }
                else if (id == DetIdVrt2)
                {
                    if (goodPidIdx && goodEkIdx)
                    {
                        h_hitXY_Vrt2[pid_idx][idx]->Fill(hitX, hitY);
                        h_hitZ_Vrt2[pid_idx][idx]->Fill(hitZ);
                        h_hitVx_Vrt2[pid_idx][idx]->Fill(hitVx);
                        h_hitVy_Vrt2[pid_idx][idx]->Fill(hitVy);
                        h_hitVz_Vrt2[pid_idx][idx]->Fill(hitVz);
                        h_hitVxVy_Vrt2[pid_idx][idx]->Fill(hitVx, hitVy);
                        h_hitVxVz_Vrt2[pid_idx][idx]->Fill(hitVx, hitVz);
                        h_hitVyVz_Vrt2[pid_idx][idx]->Fill(hitVy, hitVz);
                        vector<int>::iterator itr = find(TrackID->begin(), TrackID->end(), trkId);
                        bool onVirtual1 = false;
                        double hitPos[2] = {-9999, -9999};
                        double hitPosY = -9999;
                        for (int k = 0; k < TrackID->size(); ++k)
                        {
                            if (TrackID->at(k) == trkId && DetID->at(k) == DetIdVrt1)
                            {
                                onVirtual1 = true;
                                hitPos[0] = HitInPositionX->at(k) / 10;
                                hitPos[1] = HitInPositionZ->at(k) / 10 + 850;
                                hitPosY = HitInPositionY->at(k) / 10;
                                break;
                            }
                        }
                        if (onVirtual1)
                        {
                            Rotation(hitPos);
                            h_hitXY_Vrt12_1[pid_idx][idx]->Fill(scaleFactor * hitPos[0], scaleFactor * hitPosY);
                            h_hitXY_Vrt12_2[pid_idx][idx]->Fill(hitX, hitY);
                            h_hitVx_Vrt12[pid_idx][idx]->Fill(hitVx);
                            h_hitVy_Vrt12[pid_idx][idx]->Fill(hitVy);
                            h_hitVz_Vrt12[pid_idx][idx]->Fill(hitVz);
                            h_hitVxVy_Vrt12[pid_idx][idx]->Fill(hitVx, hitVy);
                            h_hitVxVz_Vrt12[pid_idx][idx]->Fill(hitVx, hitVz);
                            h_hitVyVz_Vrt12[pid_idx][idx]->Fill(hitVy, hitVz);
                        }
                    }
                }
            }

            if (id == DetIdVrt1)
            {
                if (!inField)
                {
                    continue;
                }
                if (goodPidIdx)
                {
                    h_hitMom_Vrt1[pid_idx]->Fill(hitMom.Mag());
                    h_hitEk_Vrt1[pid_idx]->Fill(Ek);
                }
            }
            else if (id == DetIdVrt2)
            {
                if (goodPidIdx)
                {
                    h_hitMom_Vrt2[pid_idx]->Fill(hitMom.Mag());
                    h_hitEk_Vrt2[pid_idx]->Fill(Ek);
                }
            }
        }
    } // End of event loop
    ofile->Write();
    ofile->Close();
    // file->Close();
    return 0;
}
TTree *MakeTree(TFile *file)
{
    if (!file || file->IsZombie())
    {
        std::cerr << "Error opening file!" << std::endl;
        return 0;
    }
    TTree *tree = (TTree *)file->Get("T");
    if (!tree)
    {
        std::cerr << "Error getting tree!" << std::endl;
        return 0;
    }

    // Define variables to hold the branch data
    // Set branch addresses
    tree->SetBranchAddress("EventID", &EventID);
    tree->SetBranchAddress("NHits", &NHits);
    tree->SetBranchAddress("PDGID", &PDGID);
    tree->SetBranchAddress("TrackID", &TrackID);
    tree->SetBranchAddress("ParentTrackID", &ParentTrackID);
    tree->SetBranchAddress("DetID", &DetID);
    tree->SetBranchAddress("ProcessID", &ProcessID);
    tree->SetBranchAddress("HitInPosition.X", &HitInPositionX);
    tree->SetBranchAddress("HitInPosition.Y", &HitInPositionY);
    tree->SetBranchAddress("HitInPosition.Z", &HitInPositionZ);
    tree->SetBranchAddress("HitInMomentum.X", &HitInMomentumX);
    tree->SetBranchAddress("HitInMomentum.Y", &HitInMomentumY);
    tree->SetBranchAddress("HitInMomentum.Z", &HitInMomentumZ);
    tree->SetBranchAddress("VertexPosition.X", &VertexPositionX);
    tree->SetBranchAddress("VertexPosition.Y", &VertexPositionY);
    tree->SetBranchAddress("VertexPosition.Z", &VertexPositionZ);
    tree->SetBranchAddress("HitOutPosition.X", &HitOutPositionX);
    tree->SetBranchAddress("HitOutPosition.Y", &HitOutPositionY);
    tree->SetBranchAddress("HitOutPosition.Z", &HitOutPositionZ);
    tree->SetBranchAddress("HitGlobalTime", &HitGlobalTime);
    tree->SetBranchAddress("HitEdep", &HitEdep);
    tree->SetBranchAddress("TrackLength", &TrackLength);
    tree->SetBranchAddress("TOF.NHits", &TOF_NHits);
    tree->SetBranchAddress("TOF.PDGID", &TOF_PDGID);
    tree->SetBranchAddress("TOF.TrackID", &TOF_TrackID);
    tree->SetBranchAddress("TOF.ParentTrackID", &TOF_ParentTrackID);
    tree->SetBranchAddress("TOF.DetID", &TOF_DetID);
    tree->SetBranchAddress("TOF.ProcessID", &TOF_ProcessID);
    tree->SetBranchAddress("TOF.HitInPosition.X", &TOF_HitInPositionX);
    tree->SetBranchAddress("TOF.HitInPosition.Y", &TOF_HitInPositionY);
    tree->SetBranchAddress("TOF.HitInPosition.Z", &TOF_HitInPositionZ);
    tree->SetBranchAddress("TOF.HitInMomentum.X", &TOF_HitInMomentumX);
    tree->SetBranchAddress("TOF.HitInMomentum.Y", &TOF_HitInMomentumY);
    tree->SetBranchAddress("TOF.HitInMomentum.Z", &TOF_HitInMomentumZ);
    tree->SetBranchAddress("TOF.VertexPosition.X", &TOF_VertexPositionX);
    tree->SetBranchAddress("TOF.VertexPosition.Y", &TOF_VertexPositionY);
    tree->SetBranchAddress("TOF.VertexPosition.Z", &TOF_VertexPositionZ);
    tree->SetBranchAddress("TOF.HitOutPosition.X", &TOF_HitOutPositionX);
    tree->SetBranchAddress("TOF.HitOutPosition.Y", &TOF_HitOutPositionY);
    tree->SetBranchAddress("TOF.HitOutPosition.Z", &TOF_HitOutPositionZ);
    tree->SetBranchAddress("TOF.HitGlobalTime", &TOF_HitGlobalTime);
    tree->SetBranchAddress("TOF.HitEdep", &TOF_HitEdep);
    tree->SetBranchAddress("TOF.TrackLength", &TOF_TrackLength);
    tree->SetBranchAddress("TOF.TotalEdep", &TOF_TotalEdep);
    tree->SetBranchAddress("TOF.TotalTrackLength", &TOF_TotalTrackLength);
    tree->SetBranchAddress("TOF.ModuleEdep", &TOF_ModuleEdep);
    tree->SetBranchAddress("TOF.ModuleTrackLength", &TOF_ModuleTrackLength);
    // tree->SetBranchAddress("TOF.DID", &TOF_DID);
    ActivateBranch(tree);
    return tree;
}

void MakeChain(TChain *tree)
{
    if (!tree)
    {
        std::cerr << "Error getting tree!" << std::endl;
        return;
    }

    // Define variables to hold the branch data
    // Set branch addresses
    tree->SetBranchAddress("EventID", &EventID);
    tree->SetBranchAddress("NHits", &NHits);
    tree->SetBranchAddress("PDGID", &PDGID);
    tree->SetBranchAddress("TrackID", &TrackID);
    tree->SetBranchAddress("ParentTrackID", &ParentTrackID);
    tree->SetBranchAddress("DetID", &DetID);
    tree->SetBranchAddress("ProcessID", &ProcessID);
    tree->SetBranchAddress("HitInPosition.X", &HitInPositionX);
    tree->SetBranchAddress("HitInPosition.Y", &HitInPositionY);
    tree->SetBranchAddress("HitInPosition.Z", &HitInPositionZ);
    tree->SetBranchAddress("HitInMomentum.X", &HitInMomentumX);
    tree->SetBranchAddress("HitInMomentum.Y", &HitInMomentumY);
    tree->SetBranchAddress("HitInMomentum.Z", &HitInMomentumZ);
    tree->SetBranchAddress("VertexPosition.X", &VertexPositionX);
    tree->SetBranchAddress("VertexPosition.Y", &VertexPositionY);
    tree->SetBranchAddress("VertexPosition.Z", &VertexPositionZ);
    tree->SetBranchAddress("HitOutPosition.X", &HitOutPositionX);
    tree->SetBranchAddress("HitOutPosition.Y", &HitOutPositionY);
    tree->SetBranchAddress("HitOutPosition.Z", &HitOutPositionZ);
    tree->SetBranchAddress("HitGlobalTime", &HitGlobalTime);
    tree->SetBranchAddress("HitEdep", &HitEdep);
    tree->SetBranchAddress("TrackLength", &TrackLength);
    tree->SetBranchAddress("TOF.NHits", &TOF_NHits);
    tree->SetBranchAddress("TOF.PDGID", &TOF_PDGID);
    tree->SetBranchAddress("TOF.TrackID", &TOF_TrackID);
    tree->SetBranchAddress("TOF.ParentTrackID", &TOF_ParentTrackID);
    tree->SetBranchAddress("TOF.DetID", &TOF_DetID);
    tree->SetBranchAddress("TOF.ProcessID", &TOF_ProcessID);
    tree->SetBranchAddress("TOF.HitInPosition.X", &TOF_HitInPositionX);
    tree->SetBranchAddress("TOF.HitInPosition.Y", &TOF_HitInPositionY);
    tree->SetBranchAddress("TOF.HitInPosition.Z", &TOF_HitInPositionZ);
    tree->SetBranchAddress("TOF.HitInMomentum.X", &TOF_HitInMomentumX);
    tree->SetBranchAddress("TOF.HitInMomentum.Y", &TOF_HitInMomentumY);
    tree->SetBranchAddress("TOF.HitInMomentum.Z", &TOF_HitInMomentumZ);
    tree->SetBranchAddress("TOF.VertexPosition.X", &TOF_VertexPositionX);
    tree->SetBranchAddress("TOF.VertexPosition.Y", &TOF_VertexPositionY);
    tree->SetBranchAddress("TOF.VertexPosition.Z", &TOF_VertexPositionZ);
    tree->SetBranchAddress("TOF.HitOutPosition.X", &TOF_HitOutPositionX);
    tree->SetBranchAddress("TOF.HitOutPosition.Y", &TOF_HitOutPositionY);
    tree->SetBranchAddress("TOF.HitOutPosition.Z", &TOF_HitOutPositionZ);
    tree->SetBranchAddress("TOF.HitGlobalTime", &TOF_HitGlobalTime);
    tree->SetBranchAddress("TOF.HitEdep", &TOF_HitEdep);
    tree->SetBranchAddress("TOF.TrackLength", &TOF_TrackLength);
    tree->SetBranchAddress("TOF.TotalEdep", &TOF_TotalEdep);
    tree->SetBranchAddress("TOF.TotalTrackLength", &TOF_TotalTrackLength);
    tree->SetBranchAddress("TOF.ModuleEdep", &TOF_ModuleEdep);
    tree->SetBranchAddress("TOF.ModuleTrackLength", &TOF_ModuleTrackLength);
    // tree->SetBranchAddress("TOF.DID", &TOF_DID);
    ActivateBranch(tree);
    // return tree;
}
void Rotation(double *vec)
{
    double x = vec[0];
    double y = vec[1];
    double cos = std::sqrt(2) / 2;
    double sin = std::sqrt(2) / 2;
    vec[0] = cos * x - sin * y;
    vec[1] = sin * x + cos * y;
}
int getIdx(int x, int *vec, int size)
{
    int rt = 8;
    for (int i = 0; i < size; ++i)
    {
        if (x == vec[i])
        {
            rt = i;
            break;
        }
    }

    return rt;
}
void ActivateBranch(TTree *tree)
{
    if (!tree)
    {
        return;
    }
    tree->SetBranchStatus("*", false);
    tree->SetBranchStatus("NHits", true);
    tree->SetBranchStatus("PDGID", true);
    tree->SetBranchStatus("TrackID", true);
    // tree->SetBranchStatus("ParentTrackID", &ParentTrackID);
    tree->SetBranchStatus("DetID", true);
    // tree->SetBranchStatus("ProcessID", &ProcessID);
    tree->SetBranchStatus("HitInPosition.X", true);
    tree->SetBranchStatus("HitInPosition.Y", true);
    tree->SetBranchStatus("HitInPosition.Z", true);
    tree->SetBranchStatus("HitInMomentum.X", true);
    tree->SetBranchStatus("HitInMomentum.Y", true);
    tree->SetBranchStatus("HitInMomentum.Z", true);
    tree->SetBranchStatus("VertexPosition.X", true);
    tree->SetBranchStatus("VertexPosition.Y", true);
    tree->SetBranchStatus("VertexPosition.Z", true);
    // tree->SetBranchStatus("HitOutPosition.X", true);
    // tree->SetBranchStatus("HitOutPosition.Y", true);
    // tree->SetBranchStatus("HitOutPosition.Z", true);
    tree->SetBranchStatus("HitGlobalTime", true);
    // tree->SetBranchStatus("HitEdep", true);
    // tree->SetBranchStatus("TrackLength", &TrackLength);
    // tree->SetBranchStatus("TOF.NHits", true);
    // tree->SetBranchStatus("TOF.PDGID", true);
    // tree->SetBranchStatus("TOF.TrackID", &TOF_TrackID);
    // tree->SetBranchStatus("TOF.ParentTrackID", &TOF_ParentTrackID);
    // tree->SetBranchStatus("TOF.DetID", true);
    // tree->SetBranchStatus("TOF.ProcessID", &TOF_ProcessID);
    // tree->SetBranchStatus("TOF.HitInPosition.X", true);
    // tree->SetBranchStatus("TOF.HitInPosition.Y", true);
    // tree->SetBranchStatus("TOF.HitInPosition.Z", true);
    // tree->SetBranchStatus("TOF.HitInMomentum.X", true);
    // tree->SetBranchStatus("TOF.HitInMomentum.Y", true);
    // tree->SetBranchStatus("TOF.HitInMomentum.Z", true);
    // tree->SetBranchStatus("TOF.VertexPosition.X", true);
    // tree->SetBranchStatus("TOF.VertexPosition.Y", true);
    // tree->SetBranchStatus("TOF.VertexPosition.Z", true);
    // tree->SetBranchStatus("TOF.HitOutPosition.X", &TOF_HitOutPositionX);
    // tree->SetBranchStatus("TOF.HitOutPosition.Y", &TOF_HitOutPositionY);
    // tree->SetBranchStatus("TOF.HitOutPosition.Z", &TOF_HitOutPositionZ);
    // tree->SetBranchStatus("TOF.HitGlobalTime", true);
    // tree->SetBranchStatus("TOF.HitEdep", true);
    // tree->SetBranchStatus("TOF.TrackLength", &TOF_TrackLength);
    // tree->SetBranchStatus("TOF.TotalEdep", true);
    // tree->SetBranchStatus("TOF.TotalTrackLength", &TOF_TotalTrackLength);
    // tree->SetBranchStatus("TOF.ModuleEdep", true);
    // tree->SetBranchStatus("TOF.ModuleTrackLength", &TOF_ModuleTrackLength);
    // tree->SetBranchStatus("TOF.DID", &TOF_DID);
}

void ActivateBranch(TChain *tree)
{
    if (!tree)
    {
        return;
    }
    tree->SetBranchStatus("*", false);
    tree->SetBranchStatus("NHits", true);
    tree->SetBranchStatus("PDGID", true);
    tree->SetBranchStatus("TrackID", &TrackID);
    // tree->SetBranchStatus("ParentTrackID", &ParentTrackID);
    tree->SetBranchStatus("DetID", true);
    // tree->SetBranchStatus("ProcessID", &ProcessID);
    tree->SetBranchStatus("HitInPosition.X", true);
    tree->SetBranchStatus("HitInPosition.Y", true);
    tree->SetBranchStatus("HitInPosition.Z", true);
    tree->SetBranchStatus("HitInMomentum.X", true);
    tree->SetBranchStatus("HitInMomentum.Y", true);
    tree->SetBranchStatus("HitInMomentum.Z", true);
    tree->SetBranchStatus("VertexPosition.X", true);
    tree->SetBranchStatus("VertexPosition.Y", true);
    tree->SetBranchStatus("VertexPosition.Z", true);
    // tree->SetBranchStatus("HitOutPosition.X", true);
    // tree->SetBranchStatus("HitOutPosition.Y", true);
    // tree->SetBranchStatus("HitOutPosition.Z", true);
    tree->SetBranchStatus("HitGlobalTime", true);
    // tree->SetBranchStatus("HitEdep", true);
    // tree->SetBranchStatus("TrackLength", &TrackLength);
    // tree->SetBranchStatus("TOF.NHits", true);
    // tree->SetBranchStatus("TOF.PDGID", true);
    // tree->SetBranchStatus("TOF.TrackID", &TOF_TrackID);
    // tree->SetBranchStatus("TOF.ParentTrackID", &TOF_ParentTrackID);
    // tree->SetBranchStatus("TOF.DetID", true);
    // tree->SetBranchStatus("TOF.ProcessID", &TOF_ProcessID);
    // tree->SetBranchStatus("TOF.HitInPosition.X", true);
    // tree->SetBranchStatus("TOF.HitInPosition.Y", true);
    // tree->SetBranchStatus("TOF.HitInPosition.Z", true);
    // tree->SetBranchStatus("TOF.HitInMomentum.X", true);
    // tree->SetBranchStatus("TOF.HitInMomentum.Y", true);
    // tree->SetBranchStatus("TOF.HitInMomentum.Z", true);
    // tree->SetBranchStatus("TOF.VertexPosition.X", true);
    // tree->SetBranchStatus("TOF.VertexPosition.Y", true);
    // tree->SetBranchStatus("TOF.VertexPosition.Z", true);
    // tree->SetBranchStatus("TOF.HitOutPosition.X", &TOF_HitOutPositionX);
    // tree->SetBranchStatus("TOF.HitOutPosition.Y", &TOF_HitOutPositionY);
    // tree->SetBranchStatus("TOF.HitOutPosition.Z", &TOF_HitOutPositionZ);
    // tree->SetBranchStatus("TOF.HitGlobalTime", true);
    // tree->SetBranchStatus("TOF.HitEdep", true);
    // tree->SetBranchStatus("TOF.TrackLength", &TOF_TrackLength);
    // tree->SetBranchStatus("TOF.TotalEdep", true);
    // tree->SetBranchStatus("TOF.TotalTrackLength", &TOF_TotalTrackLength);
    // tree->SetBranchStatus("TOF.ModuleEdep", true);
    // tree->SetBranchStatus("TOF.ModuleTrackLength", &TOF_ModuleTrackLength);
    // tree->SetBranchStatus("TOF.DID", &TOF_DID);
}
void InitHist()
{
    TString hn = "", ht = "";

    // for (int i = 0; i < 10; ++i)
    // {
    //     hn = Form("h_nHits_Vrt1_%d", i + 1);
    //     ht = Form("Number of hits on virtual 1 with E_{k} > %.1f MeV;nHits", EdepCut[i]);
    //     h_nHits_Vrt1[i] = new TH1D(hn, ht, 9, 0, 9);
    //     hn = Form("h_nHits_Vrt2_%d", i + 1);
    //     ht = Form("Number of hits on virtual 2 with E_{k} > %.1f MeV;nHits", EdepCut[i]);
    //     h_nHits_Vrt2[i] = new TH1D(hn, ht, 9, 0, 9);
    // }
    for (int i = 0; i < 8; ++i)
    {
        for (int j = 0; j < 10; ++j)
        {
            hn = Form("h_hitXY_Vrt1_%d_%d", i, j);
            ht = Form("virtual 1: hit position x vs y for %s with E_{k} > %.1f (MeV);x (cm);y (cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitXY_Vrt1[i][j] = new TH2D(hn, ht, 180, -180, 180, 750, -750, 750);

            hn = Form("h_hitXY_Vrt2_%d_%d", i, j);
            ht = Form("virtual 2: hit position x vs y for %s with E_{k} > %.1f (MeV);x (cm);y (cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitXY_Vrt2[i][j] = new TH2D(hn, ht, 180, -180, 180, 750, -750, 750);

            hn = Form("h_hitZ_Vrt1_%d_%d", i, j);
            ht = Form("virtual 1: hit position z for %s with E_{k} > %.1f (MeV);z (cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitZ_Vrt1[i][j] = new TH1D(hn, ht, 200, -1000, -600);

            hn = Form("h_hitZ_Vrt2_%d_%d", i, j);
            ht = Form("virtual 2: hit position z for %s with E_{k} > %.1f (MeV);z (cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitZ_Vrt2[i][j] = new TH1D(hn, ht, 900, -1000, 800);

            hn = Form("h_hitVx_Vrt1_%d_%d", i, j);
            ht = Form("virtual 1: hit vertex x for %s with E_{k} > %.1f MeV;x(cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitVx_Vrt1[i][j] = new TH1D(hn, ht, 200, -200, 200);

            hn = Form("h_hitVy_Vrt1_%d_%d", i, j);
            ht = Form("virtual 1: hit vertex y for %s with E_{k} > %.1f MeV;y(cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitVy_Vrt1[i][j] = new TH1D(hn, ht, 150, -150, 150);

            hn = Form("h_hitVz_Vrt1_%d_%d", i, j);
            ht = Form("virtual 1: hit vertex z for %s with E_{k} > %.1f MeV;z(cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitVz_Vrt1[i][j] = new TH1D(hn, ht, 200, -1000, -600);

            hn = Form("h_hitVxVy_Vrt1_%d_%d", i, j);
            ht = Form("virtual 1: hit vertex x vs y for %s with E_{k} > %.1f (MeV);x (cm);y (cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitVxVy_Vrt1[i][j] = new TH2D(hn, ht, 200, -200, 200, 150, -150, 150);

            hn = Form("h_hitVxVz_Vrt1_%d_%d", i, j);
            ht = Form("virtual 1: hit vertex x vs z for %s with E_{k} > %.1f (MeV);x (cm);z (cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitVxVz_Vrt1[i][j] = new TH2D(hn, ht, 200, -200, 200, 200, -1000, -600);

            hn = Form("h_hitVyVz_Vrt1_%d_%d", i, j);
            ht = Form("virtual 1: hit vertex y vs z for %s with E_{k} > %.1f (MeV);y (cm);z (cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitVyVz_Vrt1[i][j] = new TH2D(hn, ht, 150, -150, 150, 200, -1000, -600);

            hn = Form("h_hitVx_Vrt2_%d_%d", i, j);
            ht = Form("virtual 2: hit vertex x for %s with E_{k} > %.1f MeV;x(cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitVx_Vrt2[i][j] = new TH1D(hn, ht, 300, -300, 300);

            hn = Form("h_hitVy_Vrt2_%d_%d", i, j);
            ht = Form("virtual 2: hit vertex y for %s with E_{k} > %.1f MeV;y(cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitVy_Vrt2[i][j] = new TH1D(hn, ht, 500, -500, 500);

            hn = Form("h_hitVz_Vrt2_%d_%d", i, j);
            ht = Form("virtual 2: hit vertex z for %s with E_{k} > %.1f MeV;z(cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitVz_Vrt2[i][j] = new TH1D(hn, ht, 900, -1000, 800);

            hn = Form("h_hitVxVy_Vrt2_%d_%d", i, j);
            ht = Form("virtual 2: hit vertex x vs y for %s with E_{k} > %.1f (MeV);x (cm);y (cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitVxVy_Vrt2[i][j] = new TH2D(hn, ht, 300, -300, 300, 500, -500, 500);

            hn = Form("h_hitVxVz_Vrt2_%d_%d", i, j);
            ht = Form("virtual 2: hit vertex x vs z for %s with E_{k} > %.1f (MeV);x (cm);z (cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitVxVz_Vrt2[i][j] = new TH2D(hn, ht, 300, -300, 300, 200, -1000, -600);

            hn = Form("h_hitVyVz_Vrt2_%d_%d", i, j);
            ht = Form("virtual 2: hit vertex y vs z for %s with E_{k} > %.1f (MeV);y (cm);z (cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitVyVz_Vrt2[i][j] = new TH2D(hn, ht, 500, -500, 500, 900, -1000, 800);


            hn = Form("h_hitXY_Vrt12_1_%d_%d", i, j);
            ht = Form("virtual 12: hit position x vs y on virtual 1 for %s with E_{k} > %.1f (MeV);x (cm);y (cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitXY_Vrt12_1[i][j] = new TH2D(hn, ht, 180, -180, 180, 750, -750, 750);

            hn = Form("h_hitXY_Vrt12_2_%d_%d", i, j);
            ht = Form("virtual 12: hit position x vs y on virtual 2 for %s with E_{k} > %.1f (MeV);x (cm);y (cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitXY_Vrt12_2[i][j] = new TH2D(hn, ht, 180, -180, 180, 750, -750, 750);

            hn = Form("h_hitVx_Vrt12_%d_%d", i, j);
            ht = Form("virtual 1&2: hit vertex x for %s with E_{k} > %.1f MeV;x(cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitVx_Vrt12[i][j] = new TH1D(hn, ht, 200, -200, 200);



            hn = Form("h_hitVy_Vrt12_%d_%d", i, j);
            ht = Form("virtual 1&2: hit vertex y for %s with E_{k} > %.1f MeV;y(cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitVy_Vrt12[i][j] = new TH1D(hn, ht, 150, -150, 150);

            hn = Form("h_hitVz_Vrt12_%d_%d", i, j);
            ht = Form("virtual 1&2: hit vertex z for %s with E_{k} > %.1f MeV;z(cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitVz_Vrt12[i][j] = new TH1D(hn, ht, 200, -1000, -600);

            hn = Form("h_hitVxVy_Vrt12_%d_%d", i, j);
            ht = Form("virtual 1&2: hit vertex x vs y for %s with E_{k} > %.1f (MeV);x (cm);y (cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitVxVy_Vrt12[i][j] = new TH2D(hn, ht, 200, -200, 200, 150, -150, 150);

            hn = Form("h_hitVxVz_Vrt12_%d_%d", i, j);
            ht = Form("virtual 1&2: hit vertex x vs z for %s with E_{k} > %.1f (MeV);x (cm);z (cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitVxVz_Vrt12[i][j] = new TH2D(hn, ht, 200, -200, 200, 200, -1000, -600);

            hn = Form("h_hitVyVz_Vrt12_%d_%d", i, j);
            ht = Form("virtual 1&2: hit vertex y vs z for %s with E_{k} > %.1f (MeV);y (cm);z (cm)", pid_name[i].Data(), EdepCut[j]);
            h_hitVyVz_Vrt12[i][j] = new TH2D(hn, ht, 150, -150, 150, 200, -1000, -600);
        }
    }
    // vector<double> vec;
    // Binning(vec, 1e-5, 100);
    // cout << "size: " << vec.size() << endl;
    // int nbins = vec.size() - 1;
    // double binning[vec.size()];
    // VectorToArray(vec, binning);
    int nbins = 300;
    double binning[nbins + 1];
    Binning(nbins, binning, 1e-5, 10000);
    for (int i = 0; i < 8; ++i)
    {
        hn = Form("h_hitMom_Vrt1_%d", i + 1);
        ht = Form("virtual 1: hit momuntum for %s;p (MeV)", pid_name[i].Data());
        h_hitMom_Vrt1[i] = new TH1D(hn, ht, nbins, binning);
        // h_hitMom_Vrt1[i] = new TH1D(hn, ht, 70, -4, 3);

        hn = Form("h_hitEk_Vrt1_%d", i + 1);
        ht = Form("virtual 1: hit E_{k} for %s;E_{k} (MeV)", pid_name[i].Data());
        h_hitEk_Vrt1[i] = new TH1D(hn, ht, nbins, binning);
        // h_hitEk_Vrt1[i] = new TH1D(hn, ht, 70, -4, 3);

        hn = Form("h_hitMom_Vrt2_%d", i + 1);
        ht = Form("virtual 2: hit momuntum for %s;p (MeV)", pid_name[i].Data());
        h_hitMom_Vrt2[i] = new TH1D(hn, ht, nbins, binning);
        // h_hitMom_Vrt2[i] = new TH1D(hn, ht, 70, -4, 3);

        hn = Form("h_hitEk_Vrt2_%d", i + 1);
        ht = Form("virtual 2: hit E_{k} for %s;E_{k} (MeV)", pid_name[i].Data());
        h_hitEk_Vrt2[i] = new TH1D(hn, ht, nbins, binning);
        // h_hitEk_Vrt2[i] = new TH1D(hn, ht, 70, -4, 3);

        // hn = Form("h_hitE_Vrt1_%d", i + 1);
        // ht = Form("virtual 1: hit energy for %s;E (MeV)", pid_name[i].Data());
        // h_hitE_Vrt1[i] = new TH1D(hn, ht, 160, 0, 80);
    }
}
void Binning(const int &nbins, double *binning, const double &min, const double &max)
{
    double logMin = std::log10(min);
    double logMax = std::log10(max);
    double logBinWidth = (logMax - logMin) / nbins;

    // Fill bin edges
    for (int i = 0; i <= nbins; ++i)
    {
        binning[i] = std::pow(10, logMin + i * logBinWidth);
    }
}
double getEk(const double &mom, const int &pid_idx)
{
    double mass = pid_mass[pid_idx];
    double E = sqrt(mom * mom + mass * mass);
    return E - mass;
}
int getEkIdx(const double &Ek)
{
    int rt = -1;
    for (int i = 9; i >= 0; --i)
    {
        if (Ek > EdepCut[i])
        {
            rt = i;
            break;
        }
    }
    return rt;
}
void VectorToArray(std::vector<double> &vec, double *p)
{
    for (int i = 0; i < vec.size(); ++i)
    {
        p[i] = vec.at(i);
    }
}