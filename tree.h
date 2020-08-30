#include <TTree.h>
#include <TFile.h>
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
#include <fstream>

using namespace std;

typedef struct dssd
{
  Double_t e;
  Int_t id;
  ULong64_t ts;
} dssd;

typedef struct Recoil
{
  Double_t e;
  Int_t x, y;
  ULong64_t ts;
  Double_t me;
  ULong64_t mts;

  vector<Double_t> ge;
  vector<ULong64_t> gts;
  vector<Int_t> gid;
} recoil;

typedef struct Decay
{
  Double_t e, xe, ye, boxe;
  Int_t x, y;
  ULong64_t ts;
  Double_t dt;

  vector<Double_t> ge; //xa能量
  vector<ULong64_t> gts;
  vector<Int_t> gid;
} decay;

class tree
{
public:
  vector<dssd> *br_x, *br_y, *br_box, *br_mw, *br_xa, *br_xaa, *br_gs, *br_de;
  Double_t xesum, yesum, desum, mesum, pde; //desum -dssd正背面较低的能量; pde -ppac阳极信号

  tree()
  {
    ipt = NULL;
    opt = NULL;
  }
  tree(TTree *ipt_)
  {
    ipt = ipt_;
    Init();
  }
  void Init();
  void Loop(TTree *opt_);
  void BranchOpt();
  void GetEvent();

  TTree *ipt;
  TTree *opt;

  map<ULong64_t, Recoil> maprec; //recoil
  map<ULong64_t, Decay> mapdec;  //decay
  vector<Recoil> rvec;
  vector<Decay> dvec;
};
