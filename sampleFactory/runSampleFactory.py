from os import listdir
from os.path import isfile, join

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-l", "--load", dest = "load",
                  action = "store_true", default = False,
                  help = "do not recompile the macro, instead load it a compiled library" )
parser.add_option("-i", "--inDir", dest = "inDir",
                  help = "the input directory of root files" )
parser.add_option("-o", "--outSuffix", dest = "outSuffix",
                  help = "the file name suffix for of the output root files" )
(options, args) = parser.parse_args()

from ROOT import *
chain = TChain("ntuplizer/tree")
inFiles = []
for inFile in listdir(options.inDir):
  if isfile(join(options.inDir, inFile)):
    inFiles.append(join(options.inDir, inFile))
for sample in inFiles:
  chain.Add(sample)

if not options.load:
  gSystem.CompileMacro("sampleFactory.C", "gOck")
gSystem.Load("sampleFactory_C")
instance = sampleFactory(chain)
instance.Loop("output/mc_%s.root" % options.outSuffix, "output/data_%s.root" % options.outSuffix, "output/sig_%s.root" % options.outSuffix)

newHcounter = TH1I("hCounter", "Events counter", 5,0,5)
for inFileName in inFiles:
  inFile = TFile(inFileName)
  newHcounter.SetBinContent(1, newHcounter.GetBinContent(1) + inFile.Get("ntuplizer/hCounter").GetBinContent(1))
sigFile = TFile("output/sig_%s.root" % options.outSuffix, "UPDATE")
newHcounter.Write()
sigFile.Close()
dataFile = TFile("output/data_%s.root" % options.outSuffix, "UPDATE")
newHcounter.Write()
dataFile.Close()
