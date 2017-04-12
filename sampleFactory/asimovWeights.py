from ConfigParser import ConfigParser
from pprint import pprint
from glob import glob
from ROOT import *
from sys import argv

# this is a very haphazard script for making asimov datasets
sigOrData = argv[1]

config = ConfigParser()
config.read("hatsXsects.ini")
crossSections = dict([sample, float(xsec)] for sample, xsec in config.items('Xsects_%s' % sigOrData))
#pprint(crossSections)

maxXsect = max(crossSections.values())
relativeXsects = {}
for sample in crossSections.keys():
  relativeXsects[sample] = crossSections[sample]/maxXsect
print "-----\nrelative cross sections:"
pprint(relativeXsects)
print "-----\n"

files = []
sampleMatchString = "/mnt/hadoop/users/hakala/hatsSamples_preliminary/%s*.root" % sigOrData
for sampleName in glob(sampleMatchString):
  files.append(TFile(sampleName))
#print files
files[0].ls()

nEvents = {}
for key in relativeXsects.keys():
  nTotalAnalyzed = 0
  nTotalAccepted = 0
  for sampleFile in files:
    if key.lower() in sampleFile.GetName().lower():
      nTotalAnalyzed += sampleFile.Get("hCounter").GetBinContent(1)
      nTotalAccepted += sampleFile.Get("tree").GetEntriesFast()
  nEvents[key] = {"analyzed": nTotalAnalyzed, "accepted": nTotalAccepted}
print "-----\nnEvents map:"
pprint(nEvents)
print "-----\n"

cullFractions={}
for key in nEvents.keys():
  #print "calculating cullFraction for ", key
  #print "relativeXsect is", relativeXsects[key]
  #print "nEvents[key]['analyzed'] is", nEvents[key]["analyzed"]
  cullFractions[key] = relativeXsects[key]/(nEvents[key]["analyzed"]/10e8)
maxFraction = float(max(cullFractions.values()))
pprint (cullFractions)
for key in cullFractions.keys():
  cullFractions[key] = cullFractions[key]/maxFraction
print "-----\ncullFractions:"
pprint(cullFractions)
print "-----\n"
  
   

gSystem.CompileMacro("culler.C", "gOck")
gSystem.Load("culler_C")
sampleFiles = []
sampleEventsMap = {}
for key in relativeXsects.keys():
  sampleEventsMap[key]=0
for key in sampleEventsMap.keys():
  for sampleName in glob(sampleMatchString):
    if key.lower() in sampleName.lower():
      sampleFiles.append(TFile(sampleName))
      cullHowMany = int(sampleFiles[-1].Get("tree").GetEntries()*cullFractions[key])
      sampleEventsMap[key]+= cullHowMany
print "-----\nsampleEventsMap:"
pprint(sampleEventsMap)
print "-----\n:"


for inFile in sampleFiles:
  print inFile
  inFileName = inFile.GetName()
  for key in cullFractions.keys():
    #print "key.lower() is", key.lower()
    #print "name.lower() is", inFileName.lower()
    #print "key.lower() in inFileName.lower() is:", (key.lower() in inFileName.lower())
    if key.lower() in inFileName.lower():
      tree = inFile.Get("tree")
      print "inFile %s with number of events %i corresponds to key %s" % (inFile.GetName(), tree.GetEntries(),  key)
      cullEvents=int(cullFractions[key]*tree.GetEntries())
      print "will cull %i events from it" % cullEvents
      instance = culler(tree)
      instance.Loop("culled_%s/culled_%s" % (sigOrData, inFile.GetName().replace("/mnt/hadoop/users/hakala/hatsSamples_preliminary/","")), cullEvents)
      del instance

