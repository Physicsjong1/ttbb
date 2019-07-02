from ROOT import *
import ROOT
gROOT.SetBatch()

c1 = TCanvas('c','c',800,1000)
c1.SetMargin(.2, .1, .1, .1)

f1 = TFile('analysis/ttbbdi.root','read')
f2 = TFile('analysis/tthbbdi.root','read')

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


  if x.find('nbjet') != -1:
    hist_1.GetXaxis().SetTitle('Number of b-jets')
    hist_1.GetXaxis().SetTitleSize(0.05)

  elif x.find('mbb') != -1 or x.find('pt') != -1 or x == bjet1e or x == bjet2e or x == bjet3e:
    hist_1.GetXaxis().SetTitle('GeV')
    hist_1.GetXaxis().SetTitleSize(0.05)
  
  elif x.find('dRbb') != -1:
    hist_1.GetXaxis().SetTitle('dR=#sqrt{(#Delta#phi)^{2}+(#Delta#eta)^{2}}')
    hist_1.GetXaxis().SetTitleSize(0.05)

  elif x.find('eta') != -1:
    hist_1.GetXaxis().SetTitle('#eta')
    hist_1.GetXaxis().SetTitleSize(0.05)

  elif x.find('phi') != -1:
    hist_1.GetXaxis().SetTitle('#phi')
    hist_1.GetXaxis().SetTitleSize(0.05)
  
  hist_1.GetYaxis().SetTitle('Normalized Entries') 
  hist_1.GetYaxis().SetTitleSize(0.05)
  
  if x.find('nbjet') != -1 or x.find('njet') != -1 or x.find('nelectron') != -1 or x.find('nmuon') != -1:
    hist_1.Draw()
    hist_2.Draw('same')
    hist_1.GetYaxis().SetTitle('Entries') 
    hist_1.GetYaxis().SetTitleSize(0.05)

  else:
    hist_1.DrawNormalized('hist')
    hist_2.DrawNormalized('same') 

  leg = TLegend(0.23, 0.79, 0.38, 0.89)
  leg.AddEntry(hist_1, 'ttbb','l')
  leg.AddEntry(hist_2, 'tthbb','l')
  leg.SetTextSize(0.05)
  leg.SetBorderSize(0)
  leg.Draw()

  c1.Print('analysis/hist/'+x+'.pdf')
