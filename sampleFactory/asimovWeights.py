from ConfigParser import ConfigParser
from pprint import pprint
from glob import glob
from ROOT import *

# this is a very haphazard script for making asimov datasets

config = ConfigParser()
config.read("hatsXsects.ini")
crossSections = dict([sample, float(xsec)] for sample, xsec in config.items('Xsects'))
pprint(crossSections)

maxXsect = max(crossSections.values())
relativeXsects = {}
for sample in crossSections.keys():
  relativeXsects[sample] = crossSections[sample]/maxXsect
pprint(relativeXsects)

files = []
for sampleName in glob("sigs/sig*.root"):
  files.append(TFile(sampleName))

nEvents = {}
for key in relativeXsects.keys():
  nTotalAnalyzed = 0
  nTotalAccepted = 0
  for sampleFile in files:
    print "key: %s, file: %s" % (key.lower() , sampleFile.GetName().lower())
    if key.lower() in sampleFile.GetName().lower():
      print "found sample file %s corresponding to key %s" % (sampleFile.GetName(), key)
      nTotalAnalyzed += sampleFile.Get("hCounter").GetBinContent(1)
      nTotalAccepted += sampleFile.Get("tree").GetEntriesFast()
  nEvents[key] = {"analyzed": nTotalAnalyzed, "accepted": nTotalAccepted}
  print nEvents[key]
      
print nEvents     

weights = {}
for key in relativeXsects.keys():
  weights[key] = relativeXsects[key]/(nEvents[key]["analyzed"]/10e7)
maxWeight = max(weights.values())
for key in weights.keys():
  weights[key] = weights[key]/maxWeight
print weights
  
  


gSystem.CompileMacro("culler.C", "gOck")
gSystem.Load("culler_C")
sigFiles = []
for key in relativeXsects.keys():
 for sigName in glob("sigs/*.root"):
   if key.lower() in sigName.lower():
     print "working on file %s" % sigName
     sigFiles.append(TFile(sigName))
     nEvents = int(sigFiles[-1].Get("tree").GetEntriesFast()*weights[key])
     print " --> will cull %i events from it" % (nEvents)
     tree = sigFiles[-1].Get("tree")
     instance = culler(tree)
     instance.Loop("culledSigs/culled_%s" % sigName.replace("sigs/",""), nEvents)
     del instance
     del tree

