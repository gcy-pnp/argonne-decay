
#include "tree.h"
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TString.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
  if(argc < 2)
    {
      std::cout<< "eg. ./xxx  [run number]" <<std::endl;
      exit(1);
    }
  
  int run = atoi(argv[1]);
  cout<<"input:"<<Form("dssd%d.root",run)<<endl;
  auto ipf = new TFile(Form("/data/d2/zhli_agafa/agafa/dssd/dssd%d.root",run),"READ");
  auto ipt = (TTree*)ipf->Get("opt");

  auto opf = new TFile(Form("decay%d.root",run),"RECREATE");
  auto opt = new TTree("opt","opt");
  auto it = new tree(ipt);
  it->Loop(opt);
  opt->Write();
  opf->Close();
  
  return 0;
}

