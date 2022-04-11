#include <iostream>
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TGraphPolar.h"
#include <math.h>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "THStack.h"
#include "TSpline.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TROOT.h"

int pt_bins = 51;
int phi_bins = 51;



using namespace std;


struct pt_phi_shift{
  TSpline3 *spline;
  TH1D *hist;
};

void crop_hist(TH1D *hist){
  int nbins = hist->GetNbinsX();
  int bin_low = -999;
  int bin_high = -999;
  for (size_t i = 1; i <= nbins; i++) {
    int cont = hist->GetBinContent(i);
    if(bin_low == -999 && cont > 0) bin_low = i;
    if(cont > 0) bin_high = i;
  }
  //std::cout << "bin_low = " << bin_low << ", bin_high = " << bin_high <<'\n';
  double binsize = hist->GetBinWidth(1);
  double width = binsize*(bin_high - bin_low);
  double lowedge = hist->GetBinCenter(bin_low) - 0.5*binsize;
  hist->GetXaxis()->SetRangeUser(lowedge - 0.2*width, lowedge + 1.2*width);
}

string name_add(const char* part1, int mass){
  string fullPath = "plots/";
  fullPath += part1;
  fullPath += "_";
  fullPath += to_string(mass);
  fullPath += ".png";
  return fullPath;
}

string name_add_HME(const char* part1, int mass){
  string fullPath = "HME/";
  fullPath += part1;
  fullPath += "_";
  fullPath += to_string(mass);
  fullPath += ".png";
  return fullPath;
}

void Calculations(Double_t GenMET_pt, Double_t GenMET_phi, TH2D *Two_Roots, TH1D *H_reco1, TH1D *H_reco2, TLorentzVector p_com, TLorentzVector h_bb) {
  TLorentzVector nu, nu1, h2_1, h2_2, HH_1, HH_2;
  nu.SetPx(GenMET_pt*cos(GenMET_phi));
  nu.SetPy(GenMET_pt*sin(GenMET_phi));

  nu1.SetPx(GenMET_pt*cos(GenMET_phi));
  nu1.SetPy(GenMET_pt*sin(GenMET_phi));

  Double_t Mom_calc = (125*125 + (nu.Px() + p_com.Px())*(nu.Px() + p_com.Px()) + (nu.Py()+p_com.Py())*(nu.Py()+p_com.Py()) + p_com.Pz()*p_com.Pz() - p_com.E()*p_com.E() - GenMET_pt*GenMET_pt)/GenMET_pt;
  Double_t ac, bc, cc, D, sol1, sol2;

  ac = p_com.E() - p_com.Pz();
  bc = -Mom_calc;
  cc = p_com.E() + p_com.Pz();
  D = bc*bc - 4*ac*cc;


  if (D > 0) {
    sol1 = (-bc + sqrt(D))/(2.*ac);
    sol2 = (-bc - sqrt(D))/(2.*ac);

    Double_t Sh1, Ch1, Sh2, Ch2;

    Sh1 = 0.5*(sol1 - 1/sol1);
    Ch1 = 0.5*(sol1 + 1/sol1);
    nu.SetPz(GenMET_pt*Sh1);
    nu.SetE(GenMET_pt*Ch1);

    Sh2 = 0.5*(sol2 - 1/sol2);
    Ch2 = 0.5*(sol2 + 1/sol2);

    nu1.SetPz(GenMET_pt*Sh2);
    nu1.SetE(GenMET_pt*Ch2);

    h2_1 = p_com + nu;
    h2_2 = p_com + nu1;

    HH_1 = h2_1 + h_bb;
    HH_2 = h2_2 + h_bb;

    //std::cout << "H2_1 mass " << h2_1.M() << '\n';
    //std::cout << "H2_2 mass " << h2_2.M() << '\n' << '\n';

    Two_Roots->Fill(HH_1.M(), HH_2.M());
    H_reco1->Fill(HH_1.M());
    H_reco2->Fill(HH_2.M());

  }
}

double Determinant(Double_t GenMET_pt, Double_t GenMET_phi, TLorentzVector p_com, TLorentzVector h_bb) {
  TLorentzVector nu;
  nu.SetPx(GenMET_pt*cos(GenMET_phi));
  nu.SetPy(GenMET_pt*sin(GenMET_phi));

  Double_t Mom_calc = (125*125 + (nu.Px() + p_com.Px())*(nu.Px() + p_com.Px()) + (nu.Py()+p_com.Py())*(nu.Py()+p_com.Py()) + p_com.Pz()*p_com.Pz() - p_com.E()*p_com.E() - GenMET_pt*GenMET_pt)/GenMET_pt;
  Double_t ac, bc, cc, D;

  ac = p_com.E() - p_com.Pz();
  bc = -Mom_calc;
  cc = p_com.E() + p_com.Pz();
  D = bc*bc - 4*ac*cc;

  return D;
}

void CalculationsClose(Double_t GenMET_pt, Double_t GenMET_phi, TH1D *H_reco, TLorentzVector p_com, TLorentzVector h_bb, Double_t closest) {
  TLorentzVector nu, nu1, h2_1, h2_2, HH_1, HH_2;
  nu.SetPx(GenMET_pt*cos(GenMET_phi));
  nu.SetPy(GenMET_pt*sin(GenMET_phi));

  nu1.SetPx(GenMET_pt*cos(GenMET_phi));
  nu1.SetPy(GenMET_pt*sin(GenMET_phi));

  Double_t Mom_calc = (125*125 + (nu.Px() + p_com.Px())*(nu.Px() + p_com.Px()) + (nu.Py()+p_com.Py())*(nu.Py()+p_com.Py()) + p_com.Pz()*p_com.Pz() - p_com.E()*p_com.E() - GenMET_pt*GenMET_pt)/GenMET_pt;
  Double_t ac, bc, cc, D, sol1, sol2;

  ac = p_com.E() - p_com.Pz();
  bc = -Mom_calc;
  cc = p_com.E() + p_com.Pz();
  D = bc*bc - 4*ac*cc;

  if (D >= 0) {
    sol1 = (-bc + sqrt(D))/(2.*ac);
    sol2 = (-bc - sqrt(D))/(2.*ac);

    Double_t Sh1, Ch1, Sh2, Ch2;

    Sh1 = 0.5*(sol1 - 1/sol1);
    Ch1 = 0.5*(sol1 + 1/sol1);
    nu.SetPz(GenMET_pt*Sh1);
    nu.SetE(GenMET_pt*Ch1);

    Sh2 = 0.5*(sol2 - 1/sol2);
    Ch2 = 0.5*(sol2 + 1/sol2);

    nu1.SetPz(GenMET_pt*Sh2);
    nu1.SetE(GenMET_pt*Ch2);

    h2_1 = p_com + nu;
    h2_2 = p_com + nu1;

    HH_1 = h2_1 + h_bb;
    HH_2 = h2_2 + h_bb;

    if (abs(HH_1.M() - closest) <  abs(HH_2.M() - closest)){
        H_reco->Fill(HH_1.M());
    }
    else{
        H_reco->Fill(HH_2.M());
    }

  }
}

TH1D* Getcdf_inv(TH1D *histo){
  vector<double> cdf_v;
  double nbins = histo->GetNbinsX();
  cdf_v.push_back(0.);
  double nevents = 0;
  for (size_t i = 1; i <= nbins; i++) {
    nevents += histo->GetBinContent(i);
  }
  for (int i = 1; i < nbins + 1; i++) {
    cdf_v.push_back(histo->GetBinContent(i)/nevents + cdf_v[i-1]);
  }
  TH1D *cdf = new TH1D("cdf1","cdf_gen", nbins, cdf_v.data());
  for (size_t i = 0; i < cdf_v.size() - 1; i++) {
    cdf->SetBinContent(i+1, histo->GetBinCenter(i+1));
  }
  return cdf;
}

TGraph* Getcdf_inv_graph(TH1D *histo){
  vector<double> cdf_v, x_vector;
  double nbins = histo->GetNbinsX();
  Int_t nevents = histo->GetEntries() - histo->GetBinContent(0) - histo->GetBinContent(nbins + 1);
  cdf_v.push_back(0.);
  x_vector.push_back(histo->GetBinLowEdge(1));
  for (int i = 1; i < nbins + 1; i++) {
    cdf_v.push_back(histo->GetBinContent(i)/nevents + cdf_v[i-1]);
    x_vector.push_back(histo->GetBinCenter(i));
  }
  cdf_v.push_back(1.);
  x_vector.push_back(histo->GetBinLowEdge(nbins) + histo->GetBinWidth(nbins));
  TGraph *cdf = new TGraph(cdf_v.size(), cdf_v.data(), x_vector.data());
  return cdf;
}

TSpline3* Getcdf_inv_spline(TH1D *histo){
  vector<double> cdf_v, x_vector, fitx, fity;
  double nbins = histo->GetNbinsX();
  Int_t nevents = 0;
  for (int i = 1; i <= nbins; i++) {
    nevents += histo->GetBinContent(i);
  }
  double cumm_sum = 0;
  double db = histo->GetBinWidth(1)/2.;
  for (int i = 1; i < nbins + 1; i++) {
    cumm_sum += histo->GetBinContent(i)/nevents;
    if(cumm_sum >= 0.1 && cumm_sum <= 0.9){
      cdf_v.push_back(cumm_sum);
      x_vector.push_back(histo->GetBinCenter(i) + db);
    }
  }
  double x[1000];
  double y[1000];

  for (size_t i = 0; i < cdf_v.size(); i++) {
    x[i] = cdf_v[i];
    y[i] = x_vector[i];
    //std::cout << "x = " << x[i] << ",y = " << y[i] << '\n';
  }
  std::cout << "Finish output" << '\n';
  TSpline3 *cdf = new TSpline3("cgf", x, y, cdf_v.size());
  return cdf;
}

TH1D* Getcdf(TH1D *histo){
  vector<double> cdf_v;
  double nbins = histo->GetNbinsX();
  cdf_v.push_back(0.);
  Int_t nevents = histo->GetEntries() - histo->GetBinContent(0) - histo->GetBinContent(nbins + 1);
  for (int i = 1; i < nbins + 1; i++) {
    cdf_v.push_back(histo->GetBinContent(i)/nevents + cdf_v[i-1]);
  }
  TH1D *cdf = (TH1D*) histo->Clone();
  for (size_t i = 0; i < cdf_v.size() - 1; i++) {
    cdf->SetBinContent(i+1, cdf_v[i]);
  }
  return cdf;
}

void Get_pdf(int mass, TH2D *Pt_Phi_diff_pdf) {

  for (size_t i = 1; i <= 4; i++) {
    string fullPath = "GluGluToRadionToHHTo2B2WToLNu2J_M-";
    fullPath += to_string(mass);
    fullPath += "_narrow_13TeV-madgraph/NanoAODproduction_2016_cfg_NANO_";
    fullPath += to_string(i);
    fullPath += "_Skim.root";
    TFile *myFile = TFile::Open(fullPath.data());
    TTree *t = (TTree*)myFile->Get("Events");

    Float_t GenMET_phi;
    Float_t GenMET_pt;
    Float_t GenNuFromWFromHiggs_phi[2];
    Float_t GenNuFromWFromHiggs_pt[2];

    t->SetBranchAddress("GenMET_phi", &GenMET_phi);
    t->SetBranchAddress("GenMET_pt",  &GenMET_pt);
    t->SetBranchAddress("GenNuFromWFromHiggs_phi",    &GenNuFromWFromHiggs_phi);
    t->SetBranchAddress("GenNuFromWFromHiggs_pt",     &GenNuFromWFromHiggs_pt);
    Long64_t nevent = t->GetEntries();
    std::cout << "nevent " << nevent <<'\n';
    for (size_t i = 0; i < nevent; i++) {
      t->GetEntry(i);
      Pt_Phi_diff_pdf->Fill(GenNuFromWFromHiggs_pt[0] - GenMET_pt, GenNuFromWFromHiggs_phi[0] - GenMET_phi);
      //std::cout << "GenNuFromWFromHiggs_pt[0] " << GenNuFromWFromHiggs_pt[0] << '\n';
    }
    myFile->Close();
  }

  auto c = new TCanvas("c","c");
  c->SetCanvasSize(2000, 2000);
  c->SetWindowSize(1500, 900);

  Pt_Phi_diff_pdf->GetXaxis()->SetTitle("GenNu_pt - GenMET_pt, GeV");
  Pt_Phi_diff_pdf->GetYaxis()->SetTitle("GenNu_phi - GenMET_phi, rad");
  Pt_Phi_diff_pdf->Draw("colz");
  c->SetLogz();
  c->SaveAs(name_add_HME("CDF/Pt_Phi_diff_pdf", mass).data());


  delete c;
}

void Get_cdf(int mass, vector<pt_phi_shift> &cdf) {

  double pt_low = -mass/4.;
  double pt_high = mass/6.5;
  double phi_low = -TMath::Pi()/2.;
  double phi_high = TMath::Pi()/2.;

  TH2D *Pt_Phi_diff_pdf = new TH2D("Pt_Phi_diff_pdf2", "Pt Phi diff pdf", pt_bins, pt_low, pt_high, phi_bins, phi_low, phi_high);
  Pt_Phi_diff_pdf->GetXaxis()->SetTitle("GenNu_pt - GenMET_pt, GeV");
  Pt_Phi_diff_pdf->GetYaxis()->SetTitle("GenNu_phi - GenMET_phi, rad");

  Get_pdf(mass, Pt_Phi_diff_pdf);


  auto c = new TCanvas("c","c");
  c->SetCanvasSize(2000, 2000);
  c->SetWindowSize(1500, 900);

  TH1D *X_proj =  Pt_Phi_diff_pdf->ProjectionX();
  X_proj->GetXaxis()->SetTitle("GenNu_pt - GenMET_pt, GeV");
  X_proj->Draw();
  c->SaveAs(name_add_HME("CDF/Pt_diff_pdf", mass).data());


  TH1D *Y_proj =  Pt_Phi_diff_pdf->ProjectionY();
  Y_proj->GetXaxis()->SetTitle("GenNu_phi - GenMET_phi, rad");
  Y_proj->Draw();
  c->SaveAs(name_add_HME("CDF/Phi_diff_pdf", mass).data());

  auto Pt_cdf = Getcdf_inv(X_proj);
  Pt_cdf->SetStats(0);
  Pt_cdf->Draw();


  TSpline3 *Pt_cdf_s = Getcdf_inv_spline(X_proj);
  Pt_cdf_s->SetLineColor(kRed);
  Pt_cdf_s->Draw("same");
  c->SaveAs(name_add_HME("CDF/Pt_cdf", mass).data());
  std::cout << "wombat" << '\n';
  pt_phi_shift cdf_single;

  cdf_single.spline = Pt_cdf_s;
  cdf_single.hist = Pt_cdf;
  cdf.push_back(cdf_single);

  for (int i = 0; i < pt_bins; i++) {
    double pt_bin_width = X_proj->GetBinWidth(1);
    double pt =  (i + 0.5)*pt_bin_width + pt_low;
    string hist_name = "Phi_pdf_";
    hist_name += to_string(pt);
    TH1D* phi_hist = new TH1D(hist_name.data(), hist_name.data(), phi_bins, phi_low, phi_high);
    for (size_t j = 0; j < phi_bins; j++) {
      double phi_bin_width = Y_proj->GetBinWidth(1);
      double phi = (j + 0.5)*phi_bin_width + phi_low;
      int bin2d = Pt_Phi_diff_pdf->FindBin(pt, phi);
      phi_hist->SetBinContent(j+1, Pt_Phi_diff_pdf->GetBinContent(bin2d));
    }
    TH1D* phi_cdf = Getcdf_inv(phi_hist);
    phi_cdf->SetTitle(hist_name.data());
    phi_cdf->Draw();

    TSpline3 *phi_cdf_s = Getcdf_inv_spline(phi_hist);
    phi_cdf_s->SetLineColor(kRed);
    phi_cdf_s->Draw("same");
    c->SaveAs(name_add_HME("CDF/Pt_cdf", mass*100 + i).data());

    cdf_single.spline = phi_cdf_s;
    cdf_single.hist = phi_cdf;
    cdf.push_back(cdf_single);

  }


  delete c;
}

void Get_pt_phi_shift(double &pt, double &phi, int mass, TRandom3 &r, vector<pt_phi_shift> &cdf){
  double pt_low = -mass/4.;
  double pt_high = mass/6.5;
  double phi_low = -TMath::Pi()/2.;
  double phi_high = TMath::Pi()/2.;
  double r_pt = r.Rndm();
  double r_phi = r.Rndm();
  double pt_bin_width = (pt_high - pt_low)/pt_bins;

  if (abs(r_pt - 0.5) < 0.4) {
    pt = cdf[0].spline->Eval(r_pt);
  }
  else{
    pt = cdf[0].hist->Interpolate(r_pt);
  }

  int Phi_hist = (pt - pt_low)/pt_bin_width + 1;
  if (Phi_hist < 1) Phi_hist == 1;
  if (Phi_hist > phi_bins ) Phi_hist == phi_bins ;
  if (abs(Phi_hist - 30) <= 10) {
    if (abs(r_phi - 0.5) <= 0.4) {
      phi = cdf[Phi_hist].spline->Eval(r_phi);
    }
    else{
      phi = cdf[Phi_hist].hist->Interpolate(r_phi);
    }
  }
  else{
    phi = cdf[Phi_hist].hist->Interpolate(r_phi);
  }

}

double CalculationsHME_blind(TLorentzVector p_com, TLorentzVector h_bb, int mass, int num) {

  TRandom3 r;
  auto root_blind = new TH1D("TwHME", "Both Roots HME", 301, 0.9*mass, 1.2*mass);

  for (size_t j = 0; j < 20000; j++) {
    double angle = r.Uniform(-1, 1);
    double eta_gen = -log(tan(acos(angle)/2));
    double phi_gen = r.Uniform(TMath::Pi(), TMath::Pi());
    TLorentzVector nu, hh;
    double pt = (125*125 - p_com.M()*p_com.M())/(2*(p_com.E()*cosh(eta_gen) - p_com.Pz()*sinh(eta_gen) - p_com.Px()*cos(phi_gen) - p_com.Py()*sin(phi_gen)));

    nu.SetPtEtaPhiM(pt, eta_gen, phi_gen, 0.);
    hh = nu + p_com + h_bb;

    root_blind->Fill(hh.M());
  }


  Int_t MaxBin = root_blind->GetMaximumBin();
  Int_t x,y,z;
  root_blind->GetBinXYZ(MaxBin, x, y, z);


  Double_t closest = ((TAxis*)root_blind->GetXaxis())->GetBinCenter(x); //*0.5 + ((TAxis*)Two_Roots_true->GetYaxis())->GetBinCenter(y)*0.5;

  delete root_blind;
  return closest;

}

double CalculationsHME_eta(TRandom3 r, Double_t GenMET_pt, Double_t GenMET_phi, TLorentzVector p_com, TLorentzVector h_bb, int mass, int num) {
  auto hme_eta_single = new TH1D("TwHME_eta", "HME uniform eta", 820,230, 1050);

  for (size_t j = 0; j < 10000; j++) {
    //double angle = r.Uniform(-1, 1);
    double eta_gen = r.Uniform(-6, 6);  //-log(tan(acos(angle)/2));
    double phi = GenMET_phi;
    double pt = GenMET_pt; //(125*125 - p_com.M()*p_com.M())/(2*(p_com.E()*cosh(eta_gen) - p_com.Pz()*sinh(eta_gen) - p_com.Px()*cos(phi) - p_com.Py()*sin(phi)));
    TLorentzVector nu, hh;
    nu.SetPtEtaPhiM(pt, eta_gen, phi, 0.);
    hh = nu + p_com + h_bb;
    if(abs((nu + p_com).M() - 125) < 1.){
        hme_eta_single->Fill(hh.M());
    }
  }


  Int_t MaxBin = hme_eta_single->GetMaximumBin();
  Int_t x,y,z;
  hme_eta_single->GetBinXYZ(MaxBin, x, y, z);
  Double_t closest = ((TAxis*)hme_eta_single->GetXaxis())->GetBinCenter(x);
  if((abs(closest - mass) >= sqrt(mass)) && num%100 == 50){
    std::cout << abs(closest - mass) << '\n';
    auto c = new TCanvas("c","c");
    c->SetCanvasSize(1500, 1500);
    c->SetWindowSize(1300, 900);

    crop_hist(hme_eta_single);
    hme_eta_single->Draw();
    c->SaveAs(name_add_HME("No_spread/HME_eta_single", num).data());
    delete c;
  }



   //*0.5 + ((TAxis*)Two_Roots_true->GetYaxis())->GetBinCenter(y)*0.5;

  delete hme_eta_single;
  return closest;

}

double CalculationsHME_eta_spread(TRandom3 r, Double_t GenMET_pt, Double_t GenMET_phi, TLorentzVector p_com, TLorentzVector h_bb, int mass, int num, vector<pt_phi_shift> &cdf) {
  auto hme_eta_single = new TH1D("TwHME_eta", "HME uniform eta", 820, 230, 1050);
  int min_bin = 1;
  //int max_bin = Pt_cdf->GetNbinsX();
  //std::cout << "nbins " << max_bin << '\n';
  for (size_t j = 0; j < 10000; j++) {
    double eta_gen = r.Uniform(-6, 6);  //-log(tan(acos(angle)/2));
    Double_t pt, phi;
    Get_pt_phi_shift(pt, phi, mass, r, cdf);
    pt += GenMET_pt;
    phi += GenMET_phi;

    TLorentzVector nu, hh;
    nu.SetPtEtaPhiM(pt, eta_gen, phi, 0.);
    hh = nu + p_com + h_bb;
    if(abs((nu + p_com).M() - 125) < 1. && abs(h_bb.M() - 125) < 1.){
        hme_eta_single->Fill(hh.M());
    }
  }


  Int_t MaxBin = hme_eta_single->GetMaximumBin();
  Int_t x,y,z;
  hme_eta_single->GetBinXYZ(MaxBin, x, y, z);
  Double_t closest = ((TAxis*)hme_eta_single->GetXaxis())->GetBinCenter(x); //*0.5 + ((TAxis*)Two_Roots_true->GetYaxis())->GetBinCenter(y)*0.5;


  if(abs(closest - mass) >= sqrt(mass) && r.Rndm() < 0.1){
    auto c = new TCanvas("c","c");
    c->SetCanvasSize(1500, 1500);
    c->SetWindowSize(1300, 900);
    crop_hist(hme_eta_single);
    hme_eta_single->Draw();
    c->SaveAs(name_add_HME("Spread/HME_eta_spread", num).data());
    delete c;
  }



  delete hme_eta_single;
  return closest;

}

void CalculationsHME(TRandom3 r, Double_t GenMET_pt, Double_t GenMET_phi, TLorentzVector p_com, TLorentzVector h_bb, double mass, int num, vector<double> &firstroot, vector<double> &secondroot, vector<double> &noselection, vector<pt_phi_shift> &cdf) {
  THStack hs("hs","Both Roots");
  auto *First_Root = new TH1D("Fr_HME", "First Roots HME", 820, 230, 1050);
  auto *Second_Root = new TH1D("Sr_HME", "Second Roots HME", 820, 230, 1050);
  auto *Two_Roots_true = new TH1D("Sr_HME_b", "Second Roots HME_v", 820, 230, 1050);
  auto *Two_Roots_2D = new TH2D("Sr_HME_2D", "Second Roots 2D", 820, 230, 1050, 820, 230, 1050);

  double spt  = sqrt(GenMET_pt);
  double sphi =  sqrt(2*TMath::Pi());

  for (size_t j = 0; j < 10000; j++) {
    TLorentzVector nu, nu1, h2_1, h2_2, HH_1, HH_2;
    Double_t pt, phi;
    Get_pt_phi_shift(pt, phi, mass, r, cdf);
    pt += GenMET_pt;
    phi += GenMET_phi;

    nu.SetPx(pt*cos(phi));
    nu.SetPy(pt*sin(phi));

    nu1.SetPx(pt*cos(phi));
    nu1.SetPy(pt*sin(phi));

    Double_t Mom_calc = (125*125 + (nu.Px() + p_com.Px())*(nu.Px() + p_com.Px()) + (nu.Py()+p_com.Py())*(nu.Py()+p_com.Py()) + p_com.Pz()*p_com.Pz() - p_com.E()*p_com.E() - pt*pt)/pt;
    Double_t ac, bc, cc, D, sol1, sol2;

    ac = p_com.E() - p_com.Pz();
    bc = -Mom_calc;
    cc = p_com.E() + p_com.Pz();
    D = bc*bc - 4*ac*cc;

    if (D >= 0) {
      sol1 = (-bc + sqrt(D))/(2.*ac);
      sol2 = (-bc - sqrt(D))/(2.*ac);

      Double_t Sh1, Ch1, Sh2, Ch2;

      Sh1 = 0.5*(sol1 - 1/sol1);
      Ch1 = 0.5*(sol1 + 1/sol1);
      nu.SetPz(pt*Sh1);
      nu.SetE(pt*Ch1);

      Sh2 = 0.5*(sol2 - 1/sol2);
      Ch2 = 0.5*(sol2 + 1/sol2);

      nu1.SetPz(pt*Sh2);
      nu1.SetE(pt*Ch2);

      h2_1 = p_com + nu;
      h2_2 = p_com + nu1;

      HH_1 = h2_1 + h_bb;
      HH_2 = h2_2 + h_bb;

      double m1 = HH_1.M();
      double m2 = HH_2.M();
      if(abs(h2_1.M() - 125) < 1. && abs(h2_2.M() - 125) < 1. && abs(h_bb.M() - 125) < 1.){
        First_Root->Fill(m1);
        Second_Root->Fill(m2);
        Two_Roots_true->Fill(m1);
        Two_Roots_true->Fill(m2);
        Two_Roots_2D->Fill(m1, m2);
      }


    }
  }

  Int_t MaxBin = First_Root->GetMaximumBin();
  Int_t x,y,z;
  First_Root->GetBinXYZ(MaxBin, x, y, z);
  double f_root = ((TAxis*)First_Root->GetXaxis())->GetBinCenter(x);
  firstroot.push_back(f_root);

  MaxBin = Second_Root->GetMaximumBin();
  Second_Root->GetBinXYZ(MaxBin, x, y, z);
  double s_root = ((TAxis*)Second_Root->GetXaxis())->GetBinCenter(x);
  secondroot.push_back(s_root);

  crop_hist(First_Root);
  crop_hist(Second_Root);

  hs.Add(First_Root);
  hs.Add(Second_Root);

  MaxBin = Two_Roots_true->GetMaximumBin();
  Two_Roots_true->GetBinXYZ(MaxBin, x, y, z);
  double ns_root = ((TAxis*)Two_Roots_true->GetXaxis())->GetBinCenter(x);

  noselection.push_back(ns_root);

  if (num % 100 == 99) {
    auto c = new TCanvas("c","c");
    c->SetCanvasSize(1500, 1500);
    c->SetWindowSize(1300, 900);

    First_Root->SetLineColor(kRed);
    Second_Root->SetLineColor(kBlue);

    hs.Draw("nostack");
    c->SaveAs(name_add_HME("HME2D/Both_roots", num).data());

    Two_Roots_2D->Draw();
    c->SaveAs(name_add_HME("HME_both/Both_roots", num).data());
    delete c;
  }




  delete First_Root;
  delete Second_Root;
  delete Two_Roots_true;
  delete Two_Roots_2D;
}

void Plot_hist(Double_t mass, const char* name, vector<pt_phi_shift> &cdf, TRandom3 &r) {


    TFile *myFile = TFile::Open(name);
    TTree *t = (TTree*)myFile->Get("Events");

    Float_t GenMET_phi;
    Float_t GenMET_pt;

    Float_t GenLepFromWFromHiggs_phi[2];
    Float_t GenLepFromWFromHiggs_eta[2];
    Float_t GenLepFromWFromHiggs_mass[2];
    Float_t GenLepFromWFromHiggs_pt[2];
    Int_t   GenLepFromWFromHiggs_charge[2];
    Int_t   GenLepFromWFromHiggs_pdgId;

    Float_t GenNuFromWFromHiggs_phi[2];
    Float_t GenNuFromWFromHiggs_eta[2];
    Float_t GenNuFromWFromHiggs_mass[2];
    Float_t GenNuFromWFromHiggs_pt[2];
    Int_t   GenNuFromWFromHiggs_charge[2];

    Float_t GenBQuarkFromHiggs_phi[2];
    Float_t GenBQuarkFromHiggs_eta[2];
    Float_t GenBQuarkFromHiggs_mass[2];
    Float_t GenBQuarkFromHiggs_pt[2];
    Int_t   GenBQuarkFromHiggs_charge[2];

    Float_t GenQuarkFromWFromHiggs_phi[2];
    Float_t GenQuarkFromWFromHiggs_eta[2];
    Float_t GenQuarkFromWFromHiggs_mass[2];
    Float_t GenQuarkFromWFromHiggs_pt[2];
    Int_t   GenQuarkFromWFromHiggs_charge[2];

    Float_t GenQuarkBarFromWFromHiggs_phi[2];
    Float_t GenQuarkBarFromWFromHiggs_eta[2];
    Float_t GenQuarkBarFromWFromHiggs_mass[2];
    Float_t GenQuarkBarFromWFromHiggs_pt[2];
    Int_t   GenQuarkBarFromWFromHiggs_charge[2];

    t->SetBranchAddress("GenMET_phi", &GenMET_phi);
    t->SetBranchAddress("GenMET_pt",  &GenMET_pt);

    t->SetBranchAddress("GenLepFromWFromHiggs_phi",    &GenLepFromWFromHiggs_phi);
    t->SetBranchAddress("GenLepFromWFromHiggs_eta",    &GenLepFromWFromHiggs_eta);
    t->SetBranchAddress("GenLepFromWFromHiggs_mass",   &GenLepFromWFromHiggs_mass);
    t->SetBranchAddress("GenLepFromWFromHiggs_pt",     &GenLepFromWFromHiggs_pt);
    t->SetBranchAddress("GenLepFromWFromHiggs_charge", &GenLepFromWFromHiggs_charge);
    t->SetBranchAddress("GenLepFromWFromHiggs_pdgId",  &GenLepFromWFromHiggs_pdgId);

    t->SetBranchAddress("GenNuFromWFromHiggs_phi",    &GenNuFromWFromHiggs_phi);
    t->SetBranchAddress("GenNuFromWFromHiggs_eta",    &GenNuFromWFromHiggs_eta);
    t->SetBranchAddress("GenNuFromWFromHiggs_mass",   &GenNuFromWFromHiggs_mass);
    t->SetBranchAddress("GenNuFromWFromHiggs_pt",     &GenNuFromWFromHiggs_pt);
    t->SetBranchAddress("GenNuFromWFromHiggs_charge", &GenNuFromWFromHiggs_charge);

    t->SetBranchAddress("GenBQuarkFromHiggs_phi",    &GenBQuarkFromHiggs_phi);
    t->SetBranchAddress("GenBQuarkFromHiggs_eta",    &GenBQuarkFromHiggs_eta);
    t->SetBranchAddress("GenBQuarkFromHiggs_mass",   &GenBQuarkFromHiggs_mass);
    t->SetBranchAddress("GenBQuarkFromHiggs_pt",     &GenBQuarkFromHiggs_pt);
    t->SetBranchAddress("GenBQuarkFromHiggs_charge", &GenBQuarkFromHiggs_charge);

    t->SetBranchAddress("GenQuarkFromWFromHiggs_phi",    &GenQuarkFromWFromHiggs_phi);
    t->SetBranchAddress("GenQuarkFromWFromHiggs_eta",    &GenQuarkFromWFromHiggs_eta);
    t->SetBranchAddress("GenQuarkFromWFromHiggs_mass",   &GenQuarkFromWFromHiggs_mass);
    t->SetBranchAddress("GenQuarkFromWFromHiggs_pt",     &GenQuarkFromWFromHiggs_pt);
    t->SetBranchAddress("GenQuarkFromWFromHiggs_charge", &GenQuarkFromWFromHiggs_charge);

    t->SetBranchAddress("GenQuarkBarFromWFromHiggs_phi",    &GenQuarkBarFromWFromHiggs_phi);
    t->SetBranchAddress("GenQuarkBarFromWFromHiggs_eta",    &GenQuarkBarFromWFromHiggs_eta);
    t->SetBranchAddress("GenQuarkBarFromWFromHiggs_mass",   &GenQuarkBarFromWFromHiggs_mass);
    t->SetBranchAddress("GenQuarkBarFromWFromHiggs_pt",     &GenQuarkBarFromWFromHiggs_pt);
    t->SetBranchAddress("GenQuarkBarFromWFromHiggs_charge", &GenQuarkBarFromWFromHiggs_charge);

    Long64_t nevent = t->GetEntries() - 1;
    std::cout << "nevent " << nevent <<'\n';

    auto c = new TCanvas("c","c");
    c->SetCanvasSize(1500, 1500);
    c->SetWindowSize(1300, 900);


    TH1D *H_reco1 = new TH1D("H1_reco", "HH- mass reco", 820, 230, 1050);
    TH1D *H_reco2 = new TH1D("H2_reco", "HH+ mass reco", 820, 230, 1050);
    TH1D *H_reco1_true = new TH1D("H1_reco_true", "HH- mass reco true", 820, 230, 1050);
    TH1D *H_reco2_true = new TH1D("H2_reco_true", "HH+ mass reco true", 820, 230, 1050);
    TH1D *h_lnqq_true = new TH1D("h_lnqq_true", "Higgs mass from true neutrino", 820, 230, 1050);
    TH1D *h_bb_hist = new TH1D("h_bb_hist", "Higgs mass from bb", 101, 124, 126);
    TH1D *HH_true_hist = new TH1D("h_bb_hist", "HH mass from true neutrino", 101, 0.995*mass, 1.005*mass);
    TH1D *PT_Diff_Miss =  new TH1D("PT_Diff_Miss", "Missing Pt - True pt", 101, -mass/6.5, mass/6.5);
    TH1D *Phi_Diff_Miss =  new TH1D("Phi_Diff_Miss", "Missing Phi - True pt", 101,  -TMath::Pi()/3.,  TMath::Pi()/3.);
    TH2D *Pt_Phi_correlation = new TH2D("Pt_Phi_correlation", "Pt Phi diff correlation", 101, -mass/5., mass/5., 101, -TMath::Pi()/2.,  TMath::Pi()/2.);
    TH2D *Pt_Phi_correlation_gen = new TH2D("Pt_Phi_correlation", "Pt Phi diff correlation", 101, -mass/5., mass/5., 101, -TMath::Pi()/2.,  TMath::Pi()/2.);

    TH1D *ML_root =  new TH1D("Most likley root", "Most likley root", 101, 0.995*mass, 1.005*mass);
    TH2D *Two_Roots = new TH2D("Tw", "Both Roots", 820, 230, 1050, 820, 230, 1050);
    TH2D *Two_Roots_true = new TH2D("Tw", "Both Roots", 820, 230, 1050, 820, 230, 1050);


    TH1D *H_hme = new TH1D("Hh_HME", "Hh HME", 820, 230, 1050);
    TH1D *H_hme_blind = new TH1D("HH_HME_blind", "Hh HME, random eta, phi", 820, 230, 1050);
    TH1D *H_hme_eta = new TH1D("HH_HME_eta", "Hh HME, random eta", 820, 230, 1050);
    TH1D *H_hme_eta_spread = new TH1D("HH_HME_eta_spread", "Hh HME, random eta, with spread", 820, 230, 1050);
    TH1D *H_hme_noselect = new TH1D("HH_nos", "Hh HME, no selection", 820, 230, 1050);


    for (size_t i = 0; i < nevent; i++) {
      TLorentzVector b1, b2, h_bb, l, nu_true, q, q_bar, w_qq, w_nl_true, p_com, h2_true, HH_true, pt_miss_differance;;
      t->GetEntry(i);
      b1.SetPtEtaPhiM(GenBQuarkFromHiggs_pt[0],GenBQuarkFromHiggs_eta[0],GenBQuarkFromHiggs_phi[0],GenBQuarkFromHiggs_mass[0]);
      b2.SetPtEtaPhiM(GenBQuarkFromHiggs_pt[1],GenBQuarkFromHiggs_eta[1],GenBQuarkFromHiggs_phi[1],GenBQuarkFromHiggs_mass[1]);
      nu_true.SetPtEtaPhiM(GenNuFromWFromHiggs_pt[0],GenNuFromWFromHiggs_eta[0],GenNuFromWFromHiggs_phi[0],GenNuFromWFromHiggs_mass[0]);
      l.SetPtEtaPhiM(GenLepFromWFromHiggs_pt[0],GenLepFromWFromHiggs_eta[0],GenLepFromWFromHiggs_phi[0],GenLepFromWFromHiggs_mass[0]);
      q.SetPtEtaPhiM(GenQuarkFromWFromHiggs_pt[0],GenQuarkFromWFromHiggs_eta[0],GenQuarkFromWFromHiggs_phi[0],GenQuarkFromWFromHiggs_mass[0]);
      q_bar.SetPtEtaPhiM(GenQuarkBarFromWFromHiggs_pt[0],GenQuarkBarFromWFromHiggs_eta[0],GenQuarkBarFromWFromHiggs_phi[0],GenQuarkBarFromWFromHiggs_mass[0]);

      std::vector<float> b1_v, b2_v;


      h_bb = b1 + b2;
      //std::cout << "Higgs Mass = " << () << '\n';
      //std::cout << "Higgs Mass = " << h_bb.M() << '\n';               //h_bb.SetPxPyPzE(b1.Px() + b2.Px(), b1.Py() + b2.Py(), b1.Pz() + b2.Pz(),b1.E() + b2.E());
      w_qq = q + q_bar;             //.SetPxPyPzE(q.Px() + q_bar.Px(), q.Py() + q_bar.Py(), q.Pz() + q_bar.Pz(),q.E() + q_bar.E());
      w_nl_true = l + nu_true;      //.SetPxPyPzE(l.Px() + nu_true.Px(), l.Py() + nu_true.Py(), l.Pz() + nu_true.Pz(), l.E() + nu_true.E());
      p_com = w_qq + l;             // .SetPxPyPzE(w_qq.Px() + l.Px(), w_qq.Py() + l.Py(), w_qq.Pz() + l.Pz(),w_qq.E() + l.E());
      h2_true = w_qq + w_nl_true;   //.SetPxPyPzE(w_qq.Px() + w_nl_true.Px(), w_qq.Py() + w_nl_true.Py(), w_qq.Pz() + w_nl_true.Pz(),w_qq.E() + w_nl_true.E());
      HH_true = h_bb + h2_true;

      pt_miss_differance = b1 + b2 + l + q + q_bar;

      //std::cout << "Diff Pt = " << pt_miss_differance.Pt() << " Phi = " << pt_miss_differance.Phi() << '\n';
      //std::cout << "Nu   Pt = " << nu_true.Pt() << " Phi = " << nu_true.Phi() << '\n';
      //std::cout << "Miss Pt = " << GenMET_pt << " Phi = " << GenMET_phi << '\n';


      //PT_miss->Fill(GenMET_pt);
      //PT_Diff->Fill(pt_miss_differance.Pt());
      PT_Diff_Miss->Fill(nu_true.Pt() - GenMET_pt);
      Phi_Diff_Miss->Fill(nu_true.Phi() - GenMET_phi);


      h_lnqq_true->Fill(h2_true.M());
      h_bb_hist->Fill(h_bb.M());
      HH_true_hist->Fill(HH_true.M());

      Pt_Phi_correlation->Fill(nu_true.Pt() - GenMET_pt, nu_true.Phi() - GenMET_phi);

      TH2D *Pt_Phi_Det = new TH2D("Pt_Phi_Det", "Pt_Phi_Det", 101, -mass/5., mass/5., 101, -TMath::Pi()/2.,  TMath::Pi()/2.);
      if (i%100 == 99) {
        for (double dpt = -mass/5.; dpt <= mass/5.; dpt += mass/2.5/100) {
          for (double dphi = -TMath::Pi()/2.; dphi <= TMath::Pi()/2.; dphi += TMath::Pi()/100.) {
            double d = Determinant(nu_true.Pt() + dpt, nu_true.Phi() + dphi, p_com, h_bb);
            if (d > 0) {
              Pt_Phi_Det->SetBinContent(Pt_Phi_Det->FindBin(dpt, dphi), d);
            }
            //std::cout << "dpt = " << dpt << ", dphi = " << dphi << ", D = " << d <<'\n';
          }
        }
        c->SetLogz();
        Pt_Phi_Det->Draw("colz");
        c->SaveAs(name_add_HME("Det/Det_scan", i + 100000*mass).data());
      }



      delete Pt_Phi_Det;
      Calculations(GenMET_pt, GenMET_phi, Two_Roots, H_reco1, H_reco2, p_com, h_bb);
      Calculations(nu_true.Pt(), nu_true.Phi(), Two_Roots_true, H_reco1_true, H_reco2_true, p_com, h_bb);

      //std::cout << "Higgs mass = " << h_bb.M() << '\n';
    }



    TGraph2D *plot3D = new TGraph2D();

    Int_t MaxBin = Two_Roots_true->GetMaximumBin();
    Int_t x,y,z;
    Two_Roots_true->GetBinXYZ(MaxBin, x, y, z);

    std::cout << "X = "<< ((TAxis*)Two_Roots_true->GetXaxis())->GetBinCenter(x) << ", Y = " <<  ((TAxis*)Two_Roots_true->GetYaxis())->GetBinCenter(y)  <<'\n';
    Double_t closest = ((TAxis*)Two_Roots_true->GetXaxis())->GetBinCenter(x)*0.5 + ((TAxis*)Two_Roots_true->GetYaxis())->GetBinCenter(y)*0.5;
    vector<double> firstroot, secondroot, bothroots, noselection;
    for (size_t i = 0; i < nevent; i++) {
      TLorentzVector b1, b2, h_bb, l, nu_true, q, q_bar, w_qq, w_nl_true, p_com, h2_true, HH_true, pt_miss_differance;;
      t->GetEntry(i);
      b1.SetPtEtaPhiM(GenBQuarkFromHiggs_pt[0],GenBQuarkFromHiggs_eta[0],GenBQuarkFromHiggs_phi[0],GenBQuarkFromHiggs_mass[0]);
      b2.SetPtEtaPhiM(GenBQuarkFromHiggs_pt[1],GenBQuarkFromHiggs_eta[1],GenBQuarkFromHiggs_phi[1],GenBQuarkFromHiggs_mass[1]);
      nu_true.SetPtEtaPhiM(GenNuFromWFromHiggs_pt[0],GenNuFromWFromHiggs_eta[0],GenNuFromWFromHiggs_phi[0],GenNuFromWFromHiggs_mass[0]);
      l.SetPtEtaPhiM(GenLepFromWFromHiggs_pt[0],GenLepFromWFromHiggs_eta[0],GenLepFromWFromHiggs_phi[0],GenLepFromWFromHiggs_mass[0]);
      q.SetPtEtaPhiM(GenQuarkFromWFromHiggs_pt[0],GenQuarkFromWFromHiggs_eta[0],GenQuarkFromWFromHiggs_phi[0],GenQuarkFromWFromHiggs_mass[0]);
      q_bar.SetPtEtaPhiM(GenQuarkBarFromWFromHiggs_pt[0],GenQuarkBarFromWFromHiggs_eta[0],GenQuarkBarFromWFromHiggs_phi[0],GenQuarkBarFromWFromHiggs_mass[0]);

      std::vector<float> b1_v, b2_v;

      h_bb = b1 + b2;
      //std::cout << "Higgs Mass = " << () << '\n';
      //std::cout << "Higgs Mass = " << h_bb.M() << '\n';               //h_bb.SetPxPyPzE(b1.Px() + b2.Px(), b1.Py() + b2.Py(), b1.Pz() + b2.Pz(),b1.E() + b2.E());
      w_qq = q + q_bar;             //.SetPxPyPzE(q.Px() + q_bar.Px(), q.Py() + q_bar.Py(), q.Pz() + q_bar.Pz(),q.E() + q_bar.E());
      w_nl_true = l + nu_true;      //.SetPxPyPzE(l.Px() + nu_true.Px(), l.Py() + nu_true.Py(), l.Pz() + nu_true.Pz(), l.E() + nu_true.E());
      p_com = w_qq + l;             // .SetPxPyPzE(w_qq.Px() + l.Px(), w_qq.Py() + l.Py(), w_qq.Pz() + l.Pz(),w_qq.E() + l.E());
      h2_true = w_qq + w_nl_true;   //.SetPxPyPzE(w_qq.Px() + w_nl_true.Px(), w_qq.Py() + w_nl_true.Py(), w_qq.Pz() + w_nl_true.Pz(),w_qq.E() + w_nl_true.E());
      HH_true = h_bb + h2_true;

      pt_miss_differance = b1 + b2 + l + q + q_bar;

      CalculationsHME(r, GenMET_pt, GenMET_phi, p_com, h_bb, mass, i, firstroot, secondroot, noselection, cdf);
      H_hme_eta->Fill(CalculationsHME_eta(r, GenMET_pt, GenMET_phi, p_com, h_bb, mass, i));
      H_hme_eta_spread->Fill(CalculationsHME_eta_spread(r, GenMET_pt, GenMET_phi,  p_com, h_bb, mass, i, cdf));
      CalculationsClose(nu_true.Pt(), nu_true.Phi(), ML_root, p_com, h_bb, closest);
      if (i % 10 == 9) {
        plot3D->SetPoint(2*i, nu_true.Pt() - GenMET_pt, nu_true.Phi() - GenMET_phi, firstroot[i-1]);
        plot3D->SetPoint(2*i+1,nu_true.Pt() - GenMET_pt, nu_true.Phi() - GenMET_phi, secondroot[i-1]);
      }



      std::cout << "Even processsing = " << i << '\n';
    }

    bothroots = firstroot;
    bothroots.insert(std::end(bothroots), std::begin(secondroot), std::end(secondroot));

    sort(bothroots.begin(), bothroots.end());
    double closest_HME = bothroots[bothroots.size()/2];
    for (size_t i = 0; i < firstroot.size(); i++) {
      H_hme_noselect->Fill(noselection[i]);
      if (abs(closest_HME - firstroot[i]) < abs(closest_HME - secondroot[i])) {
        H_hme->Fill(firstroot[i]);
      }
      else{
        H_hme->Fill(secondroot[i]);
      }
    }

    auto c1 = new TCanvas("c1","c1");
    c1->SetCanvasSize(1500, 1500);
    c1->SetWindowSize(1300, 900);

    crop_hist(ML_root);
    ML_root->Draw();
    c1->SaveAs(name_add("Most_likley_root", mass).data());

    crop_hist(H_reco1);
    H_reco1->Draw();
    c1->SaveAs(name_add("H_reco1", mass).data());

    crop_hist(H_reco2);
    H_reco2->Draw();
    c1->SaveAs(name_add("H_reco2", mass).data());

    crop_hist(HH_true_hist);
    HH_true_hist->Draw();
    c1->SaveAs(name_add("HH_true_hist", mass).data());

    crop_hist(H_hme);
    H_hme->Draw();
    c1->SaveAs(name_add_HME("HH_hme", mass).data());

    crop_hist(H_hme_noselect);
    H_hme_noselect->Draw();
    c1->SaveAs(name_add_HME("H_hme_noselect", mass).data());

    crop_hist(H_hme_blind);
    H_hme_blind->Draw();
    c1->SaveAs(name_add_HME("HH_hme_blind", mass).data());

    crop_hist(H_hme_eta);
    H_hme_eta->Draw();
    c1->SaveAs(name_add_HME("HH_hme_eta", mass).data());


    crop_hist(H_hme_eta_spread);
    H_hme_eta_spread->Draw();
    c1->SaveAs(name_add_HME("HH_HME_eta_spread", mass).data());


    Two_Roots_true->Draw("COL");
    c1->SaveAs(name_add("Two_Roots_true", mass).data());

    TFile *WriteFile = new TFile("raweventdisplay.root", "RECREATE");

    plot3D->SetMarkerStyle(20);
    plot3D->Draw("pcol");
    c1->SaveAs(name_add_HME("plot3D", mass).data());
    WriteFile->WriteObject(plot3D, "3Dplot");


}


void Plot_hist_TT(int mass, const char* name, vector<pt_phi_shift> &cdf, TRandom3 &r) {


    TFile *myFile = TFile::Open(name);
    TTree *t = (TTree*)myFile->Get("Events");

    Float_t GenMET_phi;
    Float_t GenMET_pt;

    Float_t GenLepFromWFromTop_phi[2];
    Float_t GenLepFromWFromTop_eta[2];
    Float_t GenLepFromWFromTop_mass[2];
    Float_t GenLepFromWFromTop_pt[2];
    Int_t   GenLepFromWFromTop_charge[2];
    Int_t   GenLepFromWFromTop_pdgId;

    Float_t GenTauFromTop_phi[2];
    Float_t GenTauFromTop_eta[2];
    Float_t GenTauFromTop_mass[2];
    Float_t GenTauFromTop_pt[2];
    Int_t   GenTauFromTop_charge[2];
    Int_t   GenTauFromTop_pdgId;

    Float_t GenNuTauFromTop_phi[2];
    Float_t GenNuTauFromTop_eta[2];
    Float_t GenNuTauFromTop_mass[2];
    Float_t GenNuTauFromTop_pt[2];
    Int_t   GenNuTauFromTop_charge[2];
    Int_t   GenNuTauFromTop_pdgId;

    Float_t GenNuFromWFromTop_phi[2];
    Float_t GenNuFromWFromTop_eta[2];
    Float_t GenNuFromWFromTop_mass[2];
    Float_t GenNuFromWFromTop_pt[2];
    Int_t   GenNuFromWFromTop_charge[2];

    Float_t GenBQuarkFromTop_phi[2];
    Float_t GenBQuarkFromTop_eta[2];
    Float_t GenBQuarkFromTop_mass[2];
    Float_t GenBQuarkFromTop_pt[2];
    Int_t   GenBQuarkFromTop_charge[2];

    Float_t GenQuarkFromWFromTop_phi[2];
    Float_t GenQuarkFromWFromTop_eta[2];
    Float_t GenQuarkFromWFromTop_mass[2];
    Float_t GenQuarkFromWFromTop_pt[2];
    Int_t   GenQuarkFromWFromTop_charge[2];



    t->SetBranchAddress("GenMET_phi", &GenMET_phi);
    t->SetBranchAddress("GenMET_pt",  &GenMET_pt);

    t->SetBranchAddress("GenTauFromTop_phi",    &GenTauFromTop_phi);
    t->SetBranchAddress("GenTauFromTop_eta",    &GenTauFromTop_eta);
    t->SetBranchAddress("GenTauFromTop_mass",   &GenTauFromTop_mass);
    t->SetBranchAddress("GenTauFromTop_pt",     &GenTauFromTop_pt);
    t->SetBranchAddress("GenTauFromTop_charge", &GenTauFromTop_charge);
    t->SetBranchAddress("GenTauFromTop_pdgId",  &GenTauFromTop_pdgId);

    t->SetBranchAddress("GenNuTauFromTop_phi",    &GenNuTauFromTop_phi);
    t->SetBranchAddress("GenNuTauFromTop_eta",    &GenNuTauFromTop_eta);
    t->SetBranchAddress("GenNuTauFromTop_mass",   &GenNuTauFromTop_mass);
    t->SetBranchAddress("GenNuTauFromTop_pt",     &GenNuTauFromTop_pt);
    t->SetBranchAddress("GenNuTauFromTop_charge", &GenNuTauFromTop_charge);
    t->SetBranchAddress("GenNuTauFromTop_pdgId",  &GenNuTauFromTop_pdgId);

    t->SetBranchAddress("GenLepFromWFromTop_phi",    &GenLepFromWFromTop_phi);
    t->SetBranchAddress("GenLepFromWFromTop_eta",    &GenLepFromWFromTop_eta);
    t->SetBranchAddress("GenLepFromWFromTop_mass",   &GenLepFromWFromTop_mass);
    t->SetBranchAddress("GenLepFromWFromTop_pt",     &GenLepFromWFromTop_pt);
    t->SetBranchAddress("GenLepFromWFromTop_charge", &GenLepFromWFromTop_charge);
    t->SetBranchAddress("GenLepFromWFromTop_pdgId",  &GenLepFromWFromTop_pdgId);

    t->SetBranchAddress("GenNuFromWFromTop_phi",    &GenNuFromWFromTop_phi);
    t->SetBranchAddress("GenNuFromWFromTop_eta",    &GenNuFromWFromTop_eta);
    t->SetBranchAddress("GenNuFromWFromTop_mass",   &GenNuFromWFromTop_mass);
    t->SetBranchAddress("GenNuFromWFromTop_pt",     &GenNuFromWFromTop_pt);
    t->SetBranchAddress("GenNuFromWFromTop_charge", &GenNuFromWFromTop_charge);

    t->SetBranchAddress("GenBQuarkFromTop_phi",    &GenBQuarkFromTop_phi);
    t->SetBranchAddress("GenBQuarkFromTop_eta",    &GenBQuarkFromTop_eta);
    t->SetBranchAddress("GenBQuarkFromTop_mass",   &GenBQuarkFromTop_mass);
    t->SetBranchAddress("GenBQuarkFromTop_pt",     &GenBQuarkFromTop_pt);
    t->SetBranchAddress("GenBQuarkFromTop_charge", &GenBQuarkFromTop_charge);

    t->SetBranchAddress("GenQuarkFromWFromTop_phi",    &GenQuarkFromWFromTop_phi);
    t->SetBranchAddress("GenQuarkFromWFromTop_eta",    &GenQuarkFromWFromTop_eta);
    t->SetBranchAddress("GenQuarkFromWFromTop_mass",   &GenQuarkFromWFromTop_mass);
    t->SetBranchAddress("GenQuarkFromWFromTop_pt",     &GenQuarkFromWFromTop_pt);
    t->SetBranchAddress("GenQuarkFromWFromTop_charge", &GenQuarkFromWFromTop_charge);



    Long64_t nevent = t->GetEntries() - 1;
    std::cout << "nevent " << nevent <<'\n';

    auto c = new TCanvas("c","c");
    c->SetCanvasSize(1500, 1500);
    c->SetWindowSize(1300, 900);


    TH1D *H_reco1 = new TH1D("H1_reco", "HH- mass reco", 820, 230, 1050);
    TH1D *H_reco2 = new TH1D("H2_reco", "HH+ mass reco", 820, 230, 1050);
    TH1D *H_reco1_true = new TH1D("H1_reco_true", "HH- mass reco true", 820, 230, 1050);
    TH1D *H_reco2_true = new TH1D("H2_reco_true", "HH+ mass reco true", 820, 230, 1050);
    TH1D *h_lnqq_true = new TH1D("h_lnqq_true", "Higgs mass from true neutrino", 820, 230, 1050);
    TH1D *h_bb_hist = new TH1D("h_bb_hist", "Higgs mass from bb", 101, 124, 126);
    TH1D *HH_true_hist = new TH1D("h_bb_hist", "HH mass from true neutrino", 820, 230, 1050);


    TH1D *ML_root =  new TH1D("Most likley root", "Most likley root", 820, 230, 1050);
    TH2D *Two_Roots = new TH2D("Tw", "Both Roots", 820, 230, 1050, 820, 230, 1050);
    TH2D *Two_Roots_true = new TH2D("Tw", "Both Roots", 820, 230, 1050, 820, 230, 1050);


    TH1D *H_hme = new TH1D("Hh_HME", "Hh HME", 820, 230, 1050);
    TH1D *H_hme_blind = new TH1D("HH_HME_blind", "Hh HME, random eta, phi", 820, 230, 1050);
    TH1D *H_hme_eta = new TH1D("HH_HME_eta", "Hh HME, random eta", 820, 230, 1050);
    TH1D *H_hme_eta_spread = new TH1D("HH_HME_eta_spread", "Hh HME, random eta, with spread", 820, 230, 1050);
    TH1D *H_hme_noselect = new TH1D("HH_nos", "Hh HME, no selection", 820, 230, 1050);

    double GenNuTop_pt = 0;
    for (size_t i = 0; i < nevent; i++) {
      TLorentzVector b1, b2, h_bb, l, nu_true, q, q_bar, w_qq, w_nl_true, p_com, h2_true, HH_true, pt_miss_differance;;
      t->GetEntry(i);
      b1.SetPtEtaPhiM(GenBQuarkFromTop_pt[0],GenBQuarkFromTop_eta[0],GenBQuarkFromTop_phi[0],GenBQuarkFromTop_mass[0]);
      b2.SetPtEtaPhiM(GenBQuarkFromTop_pt[1],GenBQuarkFromTop_eta[1],GenBQuarkFromTop_phi[1],GenBQuarkFromTop_mass[1]);
      if (abs(GenNuTop_pt - GenNuFromWFromTop_pt[0]) <= 0.001) {
        nu_true.SetPtEtaPhiM(GenNuTauFromTop_pt[0],GenNuTauFromTop_eta[0],GenNuTauFromTop_phi[0],GenNuTauFromTop_mass[0]);
        l.SetPtEtaPhiM(GenTauFromTop_pt[0],GenTauFromTop_eta[0],GenTauFromTop_phi[0],GenTauFromTop_mass[0]);
        //std::cout << "Chose tau" << '\n';
        //std::cout << "Mass = " << GenTauFromTop_mass[0] <<'\n';
      }
      else{
        nu_true.SetPtEtaPhiM(GenNuFromWFromTop_pt[0],GenNuFromWFromTop_eta[0],GenNuFromWFromTop_phi[0],GenNuFromWFromTop_mass[0]);
        l.SetPtEtaPhiM(GenLepFromWFromTop_pt[0],GenLepFromWFromTop_eta[0],GenLepFromWFromTop_phi[0],GenLepFromWFromTop_mass[0]);
      }
      q.SetPtEtaPhiM(GenQuarkFromWFromTop_pt[0],GenQuarkFromWFromTop_eta[0],GenQuarkFromWFromTop_phi[0],GenQuarkFromWFromTop_mass[0]);
      q_bar.SetPtEtaPhiM(GenQuarkFromWFromTop_pt[1],GenQuarkFromWFromTop_eta[1],GenQuarkFromWFromTop_phi[1],GenQuarkFromWFromTop_mass[1]);
      GenNuTop_pt = GenLepFromWFromTop_pt[0];

      std::vector<float> b1_v, b2_v;
      std::cout << "neutrino mass = " << nu_true.M() <<'\n';

      h_bb = b1 + b2;
      //std::cout << "Higgs Mass = " << () << '\n';
      //std::cout << "Higgs Mass = " << h_bb.M() << '\n';               //h_bb.SetPxPyPzE(b1.Px() + b2.Px(), b1.Py() + b2.Py(), b1.Pz() + b2.Pz(),b1.E() + b2.E());
      w_qq = q + q_bar;             //.SetPxPyPzE(q.Px() + q_bar.Px(), q.Py() + q_bar.Py(), q.Pz() + q_bar.Pz(),q.E() + q_bar.E());
      w_nl_true = l + nu_true;      //.SetPxPyPzE(l.Px() + nu_true.Px(), l.Py() + nu_true.Py(), l.Pz() + nu_true.Pz(), l.E() + nu_true.E());
      p_com = w_qq + l;             // .SetPxPyPzE(w_qq.Px() + l.Px(), w_qq.Py() + l.Py(), w_qq.Pz() + l.Pz(),w_qq.E() + l.E());
      h2_true = w_qq + w_nl_true;   //.SetPxPyPzE(w_qq.Px() + w_nl_true.Px(), w_qq.Py() + w_nl_true.Py(), w_qq.Pz() + w_nl_true.Pz(),w_qq.E() + w_nl_true.E());
      HH_true = h_bb + h2_true;

      pt_miss_differance = b1 + b2 + l + q + q_bar;

      //std::cout << "Diff Pt = " << pt_miss_differance.Pt() << " Phi = " << pt_miss_differance.Phi() << '\n';
      //std::cout << "Nu   Pt = " << nu_true.Pt() << " Phi = " << nu_true.Phi() << '\n';
      //std::cout << "Miss Pt = " << GenMET_pt << " Phi = " << GenMET_phi << '\n';


      //PT_miss->Fill(GenMET_pt);
      //PT_Diff->Fill(pt_miss_differance.Pt());


      h_lnqq_true->Fill(h2_true.M());
      h_bb_hist->Fill(h_bb.M());
      HH_true_hist->Fill(HH_true.M());

      Calculations(GenMET_pt, GenMET_phi, Two_Roots, H_reco1, H_reco2, p_com, h_bb);
      Calculations(nu_true.Pt(), nu_true.Phi(), Two_Roots_true, H_reco1_true, H_reco2_true, p_com, h_bb);

      //std::cout << "Higgs mass = " << h_bb.M() << '\n';
    }



    TGraph2D *plot3D = new TGraph2D();

    Int_t MaxBin = Two_Roots_true->GetMaximumBin();
    Int_t x,y,z;
    Two_Roots_true->GetBinXYZ(MaxBin, x, y, z);

    std::cout << "X = "<< ((TAxis*)Two_Roots_true->GetXaxis())->GetBinCenter(x) << ", Y = " <<  ((TAxis*)Two_Roots_true->GetYaxis())->GetBinCenter(y)  <<'\n';
    Double_t closest = ((TAxis*)Two_Roots_true->GetXaxis())->GetBinCenter(x)*0.5 + ((TAxis*)Two_Roots_true->GetYaxis())->GetBinCenter(y)*0.5;
    vector<double> firstroot, secondroot, bothroots, noselection;
    for (size_t i = 0; i < nevent; i++) {
      TLorentzVector b1, b2, h_bb, l, nu_true, q, q_bar, w_qq, w_nl_true, p_com, h2_true, HH_true, pt_miss_differance;;
      t->GetEntry(i);
      b1.SetPtEtaPhiM(GenBQuarkFromTop_pt[0],GenBQuarkFromTop_eta[0],GenBQuarkFromTop_phi[0],GenBQuarkFromTop_mass[0]);
      b2.SetPtEtaPhiM(GenBQuarkFromTop_pt[1],GenBQuarkFromTop_eta[1],GenBQuarkFromTop_phi[1],GenBQuarkFromTop_mass[1]);
      nu_true.SetPtEtaPhiM(GenNuFromWFromTop_pt[0],GenNuFromWFromTop_eta[0],GenNuFromWFromTop_phi[0],GenNuFromWFromTop_mass[0]);
      l.SetPtEtaPhiM(GenLepFromWFromTop_pt[0],GenLepFromWFromTop_eta[0],GenLepFromWFromTop_phi[0],GenLepFromWFromTop_mass[0]);
      q.SetPtEtaPhiM(GenQuarkFromWFromTop_pt[0],GenQuarkFromWFromTop_eta[0],GenQuarkFromWFromTop_phi[0],GenQuarkFromWFromTop_mass[0]);
      q_bar.SetPtEtaPhiM(GenQuarkFromWFromTop_pt[1],GenQuarkFromWFromTop_eta[1],GenQuarkFromWFromTop_phi[1],GenQuarkFromWFromTop_mass[1]);

      std::vector<float> b1_v, b2_v;

      h_bb = b1 + b2;
      //std::cout << "Higgs Mass = " << () << '\n';
      //std::cout << "Higgs Mass = " << h_bb.M() << '\n';               //h_bb.SetPxPyPzE(b1.Px() + b2.Px(), b1.Py() + b2.Py(), b1.Pz() + b2.Pz(),b1.E() + b2.E());
      w_qq = q + q_bar;             //.SetPxPyPzE(q.Px() + q_bar.Px(), q.Py() + q_bar.Py(), q.Pz() + q_bar.Pz(),q.E() + q_bar.E());
      w_nl_true = l + nu_true;      //.SetPxPyPzE(l.Px() + nu_true.Px(), l.Py() + nu_true.Py(), l.Pz() + nu_true.Pz(), l.E() + nu_true.E());
      p_com = w_qq + l;             // .SetPxPyPzE(w_qq.Px() + l.Px(), w_qq.Py() + l.Py(), w_qq.Pz() + l.Pz(),w_qq.E() + l.E());
      h2_true = w_qq + w_nl_true;   //.SetPxPyPzE(w_qq.Px() + w_nl_true.Px(), w_qq.Py() + w_nl_true.Py(), w_qq.Pz() + w_nl_true.Pz(),w_qq.E() + w_nl_true.E());
      HH_true = h_bb + h2_true;

      pt_miss_differance = b1 + b2 + l + q + q_bar;

      CalculationsHME(r, GenMET_pt, GenMET_phi, p_com, h_bb, mass, i, firstroot, secondroot, noselection, cdf);
      H_hme_eta->Fill(CalculationsHME_eta(r, GenMET_pt, GenMET_phi, p_com, h_bb, mass, i));
      H_hme_eta_spread->Fill(CalculationsHME_eta_spread(r, GenMET_pt, GenMET_phi,  p_com, h_bb, mass, i, cdf));
      CalculationsClose(nu_true.Pt(), nu_true.Phi(), ML_root, p_com, h_bb, closest);
      if (i % 10 == 9) {
        plot3D->SetPoint(2*i, nu_true.Pt() - GenMET_pt, nu_true.Phi() - GenMET_phi, firstroot[i-1]);
        plot3D->SetPoint(2*i+1,nu_true.Pt() - GenMET_pt, nu_true.Phi() - GenMET_phi, secondroot[i-1]);
      }



      std::cout << "Even processsing = " << i << '\n';
    }

    bothroots = firstroot;
    bothroots.insert(std::end(bothroots), std::begin(secondroot), std::end(secondroot));

    sort(bothroots.begin(), bothroots.end());
    double closest_HME = bothroots[bothroots.size()/2];
    for (size_t i = 0; i < firstroot.size(); i++) {
      H_hme_noselect->Fill(noselection[i]);
      if (abs(closest_HME - firstroot[i]) < abs(closest_HME - secondroot[i])) {
        H_hme->Fill(firstroot[i]);
      }
      else{
        H_hme->Fill(secondroot[i]);
      }
    }

    auto c1 = new TCanvas("c1","c1");
    c1->SetCanvasSize(1500, 1500);
    c1->SetWindowSize(1300, 900);

    crop_hist(ML_root);
    ML_root->Draw();
    c1->SaveAs(name_add("Most_likley_root_TT", mass).data());

    crop_hist(H_reco1);
    H_reco1->Draw();
    c1->SaveAs(name_add("H_reco1_TT", mass).data());

    crop_hist(H_reco2);
    H_reco2->Draw();
    c1->SaveAs(name_add("H_reco2_TT", mass).data());

    crop_hist(HH_true_hist);
    HH_true_hist->Draw();
    c1->SaveAs(name_add("HH_true_hist_TT", mass).data());

    crop_hist(H_hme);
    H_hme->Draw();
    c1->SaveAs(name_add_HME("HH_hme_TT", mass).data());

    crop_hist(H_hme_noselect);
    H_hme_noselect->Draw();
    c1->SaveAs(name_add_HME("H_hme_noselect_TT", mass).data());

    crop_hist(H_hme_blind);
    H_hme_blind->Draw();
    c1->SaveAs(name_add_HME("HH_hme_blind_TT", mass).data());

    crop_hist(H_hme_eta);
    H_hme_eta->Draw();
    c1->SaveAs(name_add_HME("HH_hme_eta_TT", mass).data());


    crop_hist(H_hme_eta_spread);
    H_hme_eta_spread->Draw();
    c1->SaveAs(name_add_HME("HH_HME_eta_spread_TT", mass).data());


    Two_Roots_true->Draw("COL");
    c1->SaveAs(name_add("Two_Roots_true-TT", mass).data());
    /*
    TFile *WriteFile = new TFile("raweventdisplay.root", "RECREATE");

    plot3D->SetMarkerStyle(20);
    plot3D->Draw("pcol");
    c1->SaveAs(name_add_HME("plot3D_TT", mass).data());
    WriteFile->WriteObject(plot3D, "3Dplot");
    */

}

int main(){
  // Create a histogram for the values we read.

  // Open the file containing the tree.
  // TFile *myFile = TFile::Open("250_lnqq_2016_Skim.root");
  // TFile *myFile = TFile::Open("NanoAODproduction_2016_cfg_NANO_Skim.root");
  //TFile *myFile = TFile::Open("GluGluToRadionToHHTo2B2WToLNu2J_M-250_narrow_13TeV-madgraph/NanoAODproduction_2016_cfg_NANO_3_Skim.root");

  vector<const char*> name;
  name.push_back("GluGluToRadionToHHTo2B2WToLNu2J_M-250_narrow_13TeV-madgraph/NanoAODproduction_2016_cfg_NANO_3_Skim.root");
  name.push_back("GluGluToRadionToHHTo2B2WToLNu2J_M-260_narrow_13TeV-madgraph/NanoAODproduction_2016_cfg_NANO_3_Skim.root");
  name.push_back("GluGluToRadionToHHTo2B2WToLNu2J_M-270_narrow_13TeV-madgraph/NanoAODproduction_2016_cfg_NANO_3_Skim.root");
  name.push_back("GluGluToRadionToHHTo2B2WToLNu2J_M-300_narrow_13TeV-madgraph/NanoAODproduction_2016_cfg_NANO_3_Skim.root");
  name.push_back("GluGluToRadionToHHTo2B2WToLNu2J_M-350_narrow_13TeV-madgraph/NanoAODproduction_2016_cfg_NANO_3_Skim.root");
  name.push_back("GluGluToRadionToHHTo2B2WToLNu2J_M-400_narrow_13TeV-madgraph/NanoAODproduction_2016_cfg_NANO_3_Skim.root");
  name.push_back("GluGluToRadionToHHTo2B2WToLNu2J_M-450_narrow_13TeV-madgraph/NanoAODproduction_2016_cfg_NANO_3_Skim.root");
  name.push_back("GluGluToRadionToHHTo2B2WToLNu2J_M-500_narrow_13TeV-madgraph/NanoAODproduction_2016_cfg_NANO_3_Skim.root");
  name.push_back("GluGluToRadionToHHTo2B2WToLNu2J_M-550_narrow_13TeV-madgraph/NanoAODproduction_2016_cfg_NANO_3_Skim.root");
  name.push_back("GluGluToRadionToHHTo2B2WToLNu2J_M-600_narrow_13TeV-madgraph/NanoAODproduction_2016_cfg_NANO_3_Skim.root");
  name.push_back("GluGluToRadionToHHTo2B2WToLNu2J_M-700_narrow_13TeV-madgraph/NanoAODproduction_2016_cfg_NANO_3_Skim.root");
  name.push_back("GluGluToRadionToHHTo2B2WToLNu2J_M-800_narrow_13TeV-madgraph/NanoAODproduction_2016_cfg_NANO_3_Skim.root");
  name.push_back("GluGluToRadionToHHTo2B2WToLNu2J_M-900_narrow_13TeV-madgraph/NanoAODproduction_2016_cfg_NANO_3_Skim.root");
  name.push_back("GluGluToRadionToHHTo2B2WToLNu2J_M-1000_narrow_13TeV-madgraph/NanoAODproduction_2016_cfg_NANO_3_Skim.root");
  vector<int> mass;
  mass.push_back(250);
  mass.push_back(260);
  mass.push_back(270);
  mass.push_back(300);
  mass.push_back(350);
  mass.push_back(400);
  mass.push_back(450);
  mass.push_back(500);
  mass.push_back(550);
  mass.push_back(600);
  mass.push_back(700);
  mass.push_back(800);
  mass.push_back(900);
  mass.push_back(1000);



  /*
  for (size_t i = 0; i < mass.size(); i++) {
    Plot_hist(mass[i], name[i]);
  }
  */
  vector<pt_phi_shift> cdf;
  Get_cdf(mass[3], cdf);
  TRandom3 r;
  double pt_low = -mass[3]/4.;
  double pt_high = mass[3]/6.5;
  double phi_low = -TMath::Pi()/2.;
  double phi_high = TMath::Pi()/2.;

  auto c = new TCanvas("c","c");
  c->SetCanvasSize(2000, 2000);
  c->SetWindowSize(1500, 900);
  c->SetLogz();

  TH2D *Pt_Phi_diff_pdf_gen = new TH2D("Pt_Phi_diff_pdf_gen", "Pt Phi diff gen pdf", pt_bins, pt_low, pt_high, phi_bins, phi_low, phi_high);
  double pt, phi;
  for (size_t i = 0; i < 300000; i++) {
    Get_pt_phi_shift(pt, phi, mass[3], r, cdf);
    Pt_Phi_diff_pdf_gen->Fill(pt, phi);
    //std::cout << "Pt = " << pt << ", Phi = " << phi <<'\n';
  }

  Pt_Phi_diff_pdf_gen->GetXaxis()->SetTitle("GenNu_pt - GenMET_pt, GeV");
  Pt_Phi_diff_pdf_gen->GetYaxis()->SetTitle("GenNu_phi - GenMET_phi, rad");
  Pt_Phi_diff_pdf_gen->Draw("colz");

  c->SaveAs(name_add_HME("CDF/Pt_Phi_diff_gen_pdf", mass[3]).data());

  Plot_hist(mass[3], name[3], cdf, r);

  const char* TT_name = "TT_bar_Skim.root";
  Plot_hist_TT(mass[3], TT_name, cdf, r);

  delete c;
  return 0;
}
