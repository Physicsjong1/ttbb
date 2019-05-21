from ROOT import *
import ROOT

c1 = TCanvas('c1','c1',800,1000)
c2 = TCanvas('c2','c2',800,1000)
c3 = TCanvas('c3','c3',800,1000)

gStyle.SetOptStat(kFALSE)

f1 = TFile('out_prefixttbb.root','read')
f2 = TFile('out_prefixtthbb.root','read')

t1 = f1.Get("tree")
t2 = f2.Get("tree")

names = [x.GetName() for x in t1.GetListOfLeaves()]

#print(names)

for x in names:
  print x
  c1.cd()
  
  t1.SetLineColor(kRed+1)
  t1.SetLineWidth(2)
  t2.SetLineColor(kBlue+1)
  t2.SetLineWidth(2)
  
  n=t1.GetEntries()

  #Make temporary histogram1 at c1
  t1.Draw(x,'')
  h1 = gPad.GetPrimitive('htemp')  
  m = h1.GetMaximum()
  h1.SetMaximum(1.25*m)
  h1.GetXaxis().SetTitleSize(0.04)
  h1.GetYaxis().SetTitleSize(0.04)
  h1.GetYaxis().SetTitle('Normalized Entries') 
  h1.DrawNormalized()

  #Move to c2 to generate another temp histogram2
  c2.cd()
  t2.Draw(x,'')
  h2 = gPad.GetPrimitive('htemp')
  
  #Move again to c1 to draw histogram2 on same canvas
  c1.cd()
  h2.DrawNormalized('same')

  leg = TLegend(0.75, 0.8, 0.9, 0.9)
  leg.AddEntry(t1, 'tthbb','l')
  leg.AddEntry(t2, 'tthbb','l')
  leg.Draw()
  
  c1.Print('figures/'+x+'.pdf')
