#ifndef _HISTOGRAMS_HH_
#define _HISTOGRAMS_HH_
class TH1D;
class TH2D;

TH1D *h_hitVx_Vrt1[8][10];
TH1D *h_hitVy_Vrt1[8][10];
TH1D *h_hitVz_Vrt1[8][10];
TH1D *h_hitVx_Vrt2[8][10];
TH1D *h_hitVy_Vrt2[8][10];
TH1D *h_hitVz_Vrt2[8][10];

TH2D *h_hitVxVy_Vrt1[8][10];
TH2D *h_hitVxVz_Vrt1[8][10];
TH2D *h_hitVyVz_Vrt1[8][10];

TH2D *h_hitVxVy_Vrt2[8][10];
TH2D *h_hitVxVz_Vrt2[8][10];
TH2D *h_hitVyVz_Vrt2[8][10];

TH1D *h_hitVx_Vrt12[8][10];
TH1D *h_hitVy_Vrt12[8][10];
TH1D *h_hitVz_Vrt12[8][10];

TH2D *h_hitVxVy_Vrt12[8][10];
TH2D *h_hitVxVz_Vrt12[8][10];
TH2D *h_hitVyVz_Vrt12[8][10];

TH2D *h_hitXY_Vrt12_1[8][10];
TH2D *h_hitXY_Vrt12_2[8][10];

TH2D *h_hitXY_Vrt1[8][10];
TH2D *h_hitXY_Vrt2[8][10];
TH1D *h_hitZ_Vrt1[8][10];
TH1D *h_hitZ_Vrt2[8][10];
TH1D *h_hitMom_Vrt1[8];
TH1D *h_hitMom_Vrt2[8];
TH1D *h_hitEk_Vrt1[8];
TH1D *h_hitEk_Vrt2[8];
// TH1D *h_hitE_Vrt1[8];
// TH1D *h_hitE_Vrt2[8];

#endif