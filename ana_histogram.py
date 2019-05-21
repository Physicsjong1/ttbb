from ROOT import *
import ROOT
gROOT.SetBatch()

c1 = TCanvas('c','c',800,1000)
c1.SetMargin(.15, .1, .1, .1)

f1 = TFile('ttbbo.root','read')
f2 = TFile('ttho.root','read')

names = [x.GetName() for x in f1.GetListOfKeys()]
names.remove('tree')

for x in names:
	  
  gStyle.SetOptStat(kFALSE)

  hist_1 = f1.Get(x)
  hist_2 = f2.Get(x)

  maxi = max(hist_1.GetMaximum(), hist_2.GetMaximum())
  hist_1.SetMaximum(maxi*1.2)
  c1.cd()
  hist_1.SetLineColor(kRed+2)
  hist_1.SetLineWidth(2)
  hist_2.SetLineColor(kBlue+2)
  hist_2.SetLineWidth(2)
  
  if x == 'nbjet' or x == 'gennbjet' or x == 'matchednbjet':
    hist_1.GetXaxis().SetTitle('Number of b-jets')
    hist_1.GetYaxis().SetTitle('Normalized Entries')
    hist_1.GetXaxis().SetTitleSize(0.04)
    hist_1.GetYaxis().SetTitleSize(0.04) 
  elif x == 'mbb' or x == 'genmbb' or x == 'matchedmbb':
    hist_1.GetXaxis().SetTitle('Invariant Mass [GeV]')
    hist_1.GetYaxis().SetTitle('Normalized Entries')
    hist_1.GetXaxis().SetTitleSize(0.04)
    hist_1.GetYaxis().SetTitleSize(0.04) 
  
  elif x == 'dRbb' or x == 'gendRbb' or x == 'matcheddRbb':
    hist_1.GetXaxis().SetTitle('dR=#sqrt{(#Delta#phi)^{2}+(#Delta#eta)^{2}}')
    hist_1.GetYaxis().SetTitle('Normalized Entries')
    hist_1.GetXaxis().SetTitleSize(0.04)
    hist_1.GetYaxis().SetTitleSize(0.04) 
    hist_1.SetAxisRange(0.,4.,"x")

  hist_1.DrawNormalized('hist')
  hist_2.DrawNormalized('same')
  
  leg = TLegend(0.15, 0.8, 0.3, 0.9)
  leg.AddEntry(hist_1, 'ttbb','l')
  leg.AddEntry(hist_2, 'tthbb','l')
  leg.Draw()

  c1.Print('figures/'+x+'.pdf')
