#include "tree.h"
using namespace std;

void tree::Init()
{
  if (ipt == NULL)
    return;

  br_gs = 0;
  br_mw = 0;
  br_x = 0;
  br_y = 0;
  br_box = 0;
  br_xa = 0;
  br_de = 0;

  ipt->SetBranchAddress("gs", &br_gs);
  ipt->SetBranchAddress("mw", &br_mw);
  ipt->SetBranchAddress("x", &br_x);
  ipt->SetBranchAddress("y", &br_y);
  ipt->SetBranchAddress("de", &br_de);
  ipt->SetBranchAddress("box", &br_box);
  ipt->SetBranchAddress("xa", &br_xa);
  ipt->SetBranchAddress("xesum", &xesum);
  ipt->SetBranchAddress("yesum", &yesum);
  ipt->SetBranchAddress("desum", &desum);
  ipt->SetBranchAddress("mesum", &mesum);
}

void tree::Loop(TTree *opt_)
{
  if (opt_ == NULL)
    return;
  opt = opt_;
  BranchOpt();

  maprec.clear();
  mapdec.clear();

  Long64_t nentries = ipt->GetEntriesFast();
  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    ipt->GetEntry(jentry);
    if ((*br_x).size() == 0 || (*br_y).size() == 0)
      continue;
    //front-back correlation
    if (abs(xesum - yesum) > 500)
      continue;

    pde = 0;
    if ((*br_de).size() > 0)
      pde = (*br_de)[0].e;

    if ((mesum > 10 || pde > 4000) && desum > 10000 && desum < 50000)
    {
      Recoil re;
      re.e = desum;
      re.x = (*br_x)[0].id;
      re.y = (*br_y)[0].id;
      re.ts = (*br_x)[0].ts;
      re.mts = (*br_mw)[0].ts;
      if ((*br_mw).size() > 1)
      {
        for (int i = 1; i < (*br_mw).size(); i++)
          re.mts = TMath::Min(re.mts, (*br_mw)[i].ts);
      }
      re.me = mesum;
      for (int i = 0; i < (*br_gs).size(); i++)
      {
        re.ge.push_back((*br_gs)[i].e);
        re.gid.push_back((*br_gs)[i].id);
        re.gts.push_back((*br_gs)[i].ts);
      }
      maprec.insert(make_pair(re.ts, re));
    }  

    if (mesum < 10 && pde < 1000 && desum < 10000)
    { //decay
      Decay de;
      de.e = desum;
      de.xe = (*br_x)[0].e;
      de.ye = (*br_y)[0].e;
      de.x = (*br_x)[0].id;
      de.y = (*br_y)[0].id;
      de.ts = (*br_x)[0].ts;
      de.boxe = 0;
      if ((*br_box).size() > 0)
        de.boxe = (*br_box)[0].e;
      for (int i = 0; i < (*br_xa).size(); i++)
      {
        de.ge.push_back((*br_xa)[i].e);
        de.gid.push_back((*br_xa)[i].id);
        de.gts.push_back((*br_xa)[i].ts);
      }
      mapdec.insert(make_pair(de.ts, de));
    }
    if (jentry % 10000 == 0)
      cout << "Process " << jentry << " / " << nentries << endl;
  }
  cout << "recoil" << maprec.size() << "\t"
       << "decay" << mapdec.size() << endl;

  ULong64_t us = 100;
  ULong64_t ms = 1000 * us;
  ULong64_t sec = 1000 * ms;
  ULong64_t twindow = 200 * sec; //200s
  int i = 0;
  cout << "recoil-decay " << endl;
  for (auto ia = maprec.begin(); ia != maprec.end(); ia++)
  {
    rvec.clear();
    dvec.clear();
    //recoil-decay correlation
    auto ib = mapdec.lower_bound(ia->first - 20 * sec);
    for (auto ic = ib; ic != mapdec.end(); ic++)
    {
      if (ic->first >= ia->first + twindow)
        break;
      Int_t delx = abs(ic->second.x - ia->second.x);
      Int_t dely = abs(ic->second.y - ia->second.y);
      if (delx > 0 || dely > 0)
        continue;
      ic->second.dt = Long64_t(ic->first - ia->first) / (1.0e5);
      dvec.push_back(ic->second);
    }

    if (dvec.size() > 0)
    {
      rvec.push_back(ia->second);
      opt->Fill();
    }
    if (i % 5000 == 0)
      cout << "Process " << i << " / " << maprec.size() << endl;
    i++;
  }

  cout << endl;
}

void tree::BranchOpt()
{

  opt->Branch("recoil", &rvec);
  opt->Branch("decay", &dvec);
}
