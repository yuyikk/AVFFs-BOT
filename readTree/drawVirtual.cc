#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TKey.h>
#include "TH2D.h"
#include <TString.h>
#include <iostream>
#include <map>
#include <vector>
#include "Histograms.hh"

double timescale = timescale = 1. * atof(getenv("MULT")) * pow(10, atoi(getenv("POW"))) * 1.6e-9; // unit us // in the unit of us
const double defaultLegXY[4] = {0.15, 0.3, 0.3, 0.85};
const char *pid_name[8] = {"e^{-}", "#gamma", "e^{+}", "neutron", "proton", "#pi^{+}", "#pi^{-}", "#pi^{0}"};
const char *EkCut[10] = {"E_{k} (MeV) > 0",
                        "E_{k} (MeV) > 0.1",
                        "E_{k} (MeV) > 0.5",
                        "E_{k} (MeV) > 1",
                        "E_{k} (MeV) > 2",
                        "E_{k} (MeV) > 5",
                        "E_{k} (MeV) > 10",
                        "E_{k} (MeV) > 15",
                        "E_{k} (MeV) > 20",
                        "E_{k} (MeV) > 30"};
const double EkCutVal[10] = {0, 0.1, 0.5, 1, 2, 5, 10, 15, 20, 30};
int pid_color[8] = {kRed - 4, kOrange - 2, kBlue, kCyan, kMagenta, kGreen, kBlack, kGreen + 2};
const double areaVrt1 = 20 * 84;
const double areaVrt2 = 200 * 840;
TFile *file = 0;
void InitHist(TFile *f);
// void drawEk(TH1D *h);
void drawEk();
void drawVertex();
using namespace std;
void drawVirtual()
{
    cout << Form("Number of events %1f", 1. * atof(getenv("MULT")) * pow(10, atoi(getenv("POW")))) << "\t" << "timescale: " << timescale << endl;
    gStyle->SetPadLeftMargin(0.16);
    gStyle->SetPadBottomMargin(0.17);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetOptStat(false);
    gStyle->SetLabelSize(0.05, "XY");
    // Open the ROOT file containing the histograms
    // file = TFile::Open("./0907-0908.BOT.Vritual.Vacuum.root"); // Replace with your actual file name
    file = TFile::Open(Form("%s", getenv("INPUT"))); // Replace with your actual file name
    if (!file || file->IsZombie())
    {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    InitHist(file);
    drawEk();
    drawVertex();
    file->Close();
}
void drawEk()
{

    // Draw Ek distributions
    TCanvas *c = new TCanvas("c", "c", 1200, 900);
    c->SetLogx();
    c->SetLogy();
    TLegend *legend = new TLegend({0.2, 0.3, 0.35, 0.85});
    legend->SetTextSize(0.04);
    legend->SetFillColorAlpha(0, 0);
    int rebin[8] = {2, 2, 4, 6, 6, 6, 6, 6};
    for (int i = 0; i < 8; ++i)
    {
        h_hitEk_Vrt1[i]->Rebin(rebin[i]);
        for (int ib = 0; ib < h_hitEk_Vrt1[i]->GetNbinsX(); ++ib)
        {
            double scale = 1. / (timescale * h_hitEk_Vrt1[i]->GetBinWidth(ib + 1) * areaVrt1);
            h_hitEk_Vrt1[i]->SetBinContent(ib + 1, scale * h_hitEk_Vrt1[i]->GetBinContent(ib + 1));
        }
        h_hitEk_Vrt1[i]->SetMaximum(1e7);
        h_hitEk_Vrt1[i]->SetMinimum(1e-10);
        h_hitEk_Vrt1[i]->SetLineColor(pid_color[i]);
        h_hitEk_Vrt1[i]->GetYaxis()->SetTitle("particle rate (MHz/cm^{2}/MeV)");
        h_hitEk_Vrt1[i]->SetTitleOffset(1.25);
        h_hitEk_Vrt1[i]->SetStats(false);
        h_hitEk_Vrt1[i]->Draw("SAME HIST");
        legend->AddEntry(h_hitEk_Vrt1[i], pid_name[i], "l");
    }
    legend->Draw("same");
    c->Draw();
    c->Update();
    c->Print(Form("%s/Vrt1_Ek.pdf", getenv("OUTDIR")));
    legend->Clear();
    TCanvas *c2 = new TCanvas("c2", "c2", 1200, 900);
    c2->SetLogx();
    c2->SetLogy();
    for (int i = 0; i < 8; ++i)
    {
        h_hitEk_Vrt2[i]->Rebin(rebin[i]);
        for (int ib = 0; ib < h_hitEk_Vrt2[i]->GetNbinsX(); ++ib)
        {
            double scale = 1. / (timescale * h_hitEk_Vrt2[i]->GetBinWidth(ib + 1) * areaVrt2);
            h_hitEk_Vrt2[i]->SetBinContent(ib + 1, scale * h_hitEk_Vrt2[i]->GetBinContent(ib + 1));
        }
        h_hitEk_Vrt2[i]->SetMaximum(1e4);
        h_hitEk_Vrt2[i]->SetMinimum(1e-10);
        h_hitEk_Vrt2[i]->SetLineColor(pid_color[i]);
        h_hitEk_Vrt2[i]->GetYaxis()->SetTitle("particle rate (MHz/cm^{2}/MeV)");
        h_hitEk_Vrt2[i]->SetTitleOffset(1.25);
        h_hitEk_Vrt2[i]->SetStats(false);
        h_hitEk_Vrt2[i]->Draw("SAME HIST");
        legend->AddEntry(h_hitEk_Vrt2[i], pid_name[i], "l");
    }
    legend->Draw("same");

    c2->Update();
    c2->Print(Form("%s/Vrt2_Ek.pdf", getenv("OUTDIR")));
}
void drawVertex()
{

    
    // for (int i = 0; i < 8; ++i) // loop for each particle
    // {
    //     TCanvas *c = new TCanvas("c", "c", 1200, 900);
    //     c->cd();
    //     c->SetLogy();
    //     TLegend *legend = new TLegend({0.15, 0.5, 0.3, 0.9});
    //     legend->SetTextSize(0.04);
    //     int rebin[8] = {2, 2, 2, 2, 2, 2, 2, 2};
    //     for (int j = 0; j < 8; ++j) // loop for each Ek cut
    //     {
    //         h_hitVz_Vrt1[i][j]->SetLineColor(pid_color[j]);
    //         h_hitVz_Vrt1[i][0]->SetMaximum(1e7);
    //         h_hitVz_Vrt1[i][j]->SetMinimum(3e-5);
    //         h_hitVz_Vrt1[i][j]->GetYaxis()->SetTitle("particle rate (MHz)");
    //         legend->AddEntry(h_hitVz_Vrt1[i][j], EkCut[j], "l");
    //         h_hitVz_Vrt1[i][j]->Draw("SAME HIST");
    //     }
    //     legend->Draw("SAME");
    //     gPad->Update();
    //     c->Print(Form("%s/h_hitVz_Vrt1_particle_%s.pdf", getenv("OUTDIR"), pid_name[i]));
        
    // }

    // for (int i = 0; i < 8; ++i) // loop for each particle
    // {
    //     TCanvas *c = new TCanvas("c", "c", 1200, 900);
    //     c->cd();
    //     c->SetLogy();
    //     TLegend *legend = new TLegend({0.5, 0.5, 0.65, 0.9});
    //     legend->SetTextSize(0.04);
    //     int rebin[8] = {20, 20, 20, 20, 20, 20, 20, 20};
    //     for (int j = 0; j < 8; ++j) // loop for each Ek cut
    //     {
    //         h_hitVz_Vrt2[i][j]->SetLineColor(pid_color[j]);
    //         h_hitVz_Vrt2[i][j]->Rebin(rebin[j]);
    //         h_hitVz_Vrt2[i][j]->Scale(1 / timescale);
    //         h_hitVz_Vrt2[i][j]->SetMinimum(3e-5);
    //         h_hitVz_Vrt2[i][j]->SetMaximum(1e7);
    //         // int bin = h_hitVz_Vrt2[i][j]->FindBin(655);
    //         // h_hitVz_Vrt2[i][j]->SetBinContent(bin, 0);
    //         // cout << "Bin Content: " << h_hitVz_Vrt2[i][j]->GetBinContent(bin) << endl;
    //         h_hitVz_Vrt2[i][j]->GetYaxis()->SetTitle("particle rate (MHz)");
    //         legend->AddEntry(h_hitVz_Vrt2[i][j], EkCut[j], "l");
    //         h_hitVz_Vrt2[i][j]->Draw("SAME HIST");
    //     }
    //     legend->Draw("SAME");
    //     gPad->Update();
    //     c->Print(Form("%s/h_hitVz_Vrt2_particle_%s.pdf", getenv("OUTDIR"),pid_name[i]));
    // }

    for (int j = 0; j < 10; ++j) // loop for each Ekcut
    {
        TCanvas *c = new TCanvas("c", "c", 1200, 900);
        c->cd();
        c->SetLogy();
        TLegend *legend = new TLegend({0.2, 0.5, 0.35, 0.9});
        legend->SetTextSize(0.04);
        legend->SetFillColorAlpha(0, 0);
        int rebin[8] = {2, 2, 2, 2, 2, 2, 2, 2};
        for (int i = 0; i < 8; ++i) // loop for each particle
        {
            h_hitVz_Vrt1[i][j]->SetLineColor(pid_color[i]);
            // h_hitVz_Vrt1[i][j]->Rebin(rebin[i]);
            h_hitVz_Vrt1[i][j]->SetMaximum(1e8);
            h_hitVz_Vrt1[i][j]->SetMinimum(1e-4);
            h_hitVz_Vrt1[i][j]->GetYaxis()->SetTitle("particle rate (MHz)");
            legend->AddEntry(h_hitVz_Vrt1[i][j], pid_name[i], "l");
            h_hitVz_Vrt1[i][j]->Draw("SAME HIST");
        }
        legend->Draw("SAME");
        gPad->Update();
        c->Print(Form("%s/h_hitVz_Vrt1_Ekcut%.1f.pdf", getenv("OUTDIR"), EkCutVal[j]));
    }

    for (int j = 0; j < 10; ++j) // loop for each Ekcut
    {
        TCanvas *c = new TCanvas("c", "c", 1200, 900);
        c->cd();
        c->SetLogy();
        TLegend *legend = new TLegend({0.5, 0.5, 0.65, 0.9});
        legend->SetTextSize(0.04);
        legend->SetFillColorAlpha(0, 0);
        int rebin[8] = {2, 2, 2, 2, 2, 2, 2, 2};
        for (int i = 0; i < 8; ++i) // loop for each particle
        {
            h_hitVz_Vrt2[i][j]->SetLineColor(pid_color[i]);
            // h_hitVz_Vrt2[i][j]->Rebin(rebin[i]);
            h_hitVz_Vrt2[i][j]->SetMinimum(1e-4);
            h_hitVz_Vrt2[i][j]->SetMaximum(1e8);
            h_hitVz_Vrt2[i][j]->GetYaxis()->SetTitle("particle rate (MHz)");
            legend->AddEntry(h_hitVz_Vrt2[i][j], pid_name[i], "l");
            h_hitVz_Vrt2[i][j]->Draw("SAME HIST");
        }
        legend->Draw("SAME");
        gPad->Update();
        c->Print(Form("%s/h_hitVz_Vrt2_Ekcut%.1f.pdf", getenv("OUTDIR"), EkCutVal[j]));
    }

}
void InitHist(TFile *f)
{
    TString hn = "";
    for (int i = 0; i < 8; ++i)
    {
        hn = Form("h_hitMom_Vrt1_%d", i + 1);
        h_hitMom_Vrt1[i] = (TH1D *)f->Get(hn.Data());

        hn = Form("h_hitEk_Vrt1_%d", i + 1);
        h_hitEk_Vrt1[i] = (TH1D *)f->Get(hn.Data());

        hn = Form("h_hitMom_Vrt2_%d", i + 1);
        h_hitMom_Vrt2[i] = (TH1D *)f->Get(hn.Data());

        hn = Form("h_hitEk_Vrt2_%d", i + 1);
        h_hitEk_Vrt2[i] = (TH1D *)f->Get(hn.Data());
    }
    for (int i = 0; i < 8; ++i)
    {
        for (int j = 0; j < 10; ++j)
        {
            hn = Form("h_hitXY_Vrt1_%d_%d", i, j);
            h_hitXY_Vrt1[i][j] = (TH2D *)f->Get(hn.Data());

            hn = Form("h_hitXY_Vrt2_%d_%d", i, j);
            h_hitXY_Vrt2[i][j] = (TH2D *)f->Get(hn.Data());

            hn = Form("h_hitZ_Vrt1_%d_%d", i, j);
            h_hitZ_Vrt1[i][j] = (TH1D *)f->Get(hn.Data());

            hn = Form("h_hitZ_Vrt2_%d_%d", i, j);
            h_hitZ_Vrt2[i][j] = (TH1D *)f->Get(hn.Data());

            hn = Form("h_hitVx_Vrt1_%d_%d", i, j);
            h_hitVx_Vrt1[i][j] = (TH1D *)f->Get(hn.Data());

            hn = Form("h_hitVy_Vrt1_%d_%d", i, j);
            h_hitVy_Vrt1[i][j] = (TH1D *)f->Get(hn.Data());

            hn = Form("h_hitVz_Vrt1_%d_%d", i, j);
            h_hitVz_Vrt1[i][j] = (TH1D *)f->Get(hn.Data());
            h_hitVz_Vrt1[i][j]->Rebin(2);
            h_hitVz_Vrt1[i][j]->Scale(1 / timescale);

            hn = Form("h_hitVxVy_Vrt1_%d_%d", i, j);
            h_hitVxVy_Vrt1[i][j] = (TH2D *)f->Get(hn.Data());

            hn = Form("h_hitVxVz_Vrt1_%d_%d", i, j);
            h_hitVxVz_Vrt1[i][j] = (TH2D *)f->Get(hn.Data());

            hn = Form("h_hitVyVz_Vrt1_%d_%d", i, j);
            h_hitVyVz_Vrt1[i][j] = (TH2D *)f->Get(hn.Data());

            hn = Form("h_hitVx_Vrt2_%d_%d", i, j);
            h_hitVx_Vrt2[i][j] = (TH1D *)f->Get(hn.Data());

            hn = Form("h_hitVy_Vrt2_%d_%d", i, j);
            h_hitVy_Vrt2[i][j] = (TH1D *)f->Get(hn.Data());

            hn = Form("h_hitVz_Vrt2_%d_%d", i, j);
            h_hitVz_Vrt2[i][j] = (TH1D *)f->Get(hn.Data());
            h_hitVz_Vrt2[i][j]->Rebin(9);
            h_hitVz_Vrt2[i][j]->Scale(1 / timescale);

            hn = Form("h_hitVxVy_Vrt2_%d_%d", i, j);
            h_hitVxVy_Vrt2[i][j] = (TH2D *)f->Get(hn.Data());

            hn = Form("h_hitVxVz_Vrt2_%d_%d", i, j);
            h_hitVxVz_Vrt2[i][j] = (TH2D *)f->Get(hn.Data());

            hn = Form("h_hitVyVz_Vrt2_%d_%d", i, j);
            h_hitVyVz_Vrt2[i][j] = (TH2D *)f->Get(hn.Data());
        }
    }
}