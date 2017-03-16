from optparse import OptionParser
parser = OptionParser()
parser.add_option("-l", "--load", dest = "load",
                  action = "store_true", default = False,
                  help = "do not recompile the macro, instead load it a compiled library" )
(options, args) = parser.parse_args()

from ROOT import *
if not options.load:
  gSystem.CompileMacro("sampleFactory.C", "gOck")
gSystem.Load("sampleFactory_C")
instance = sampleFactory()
instance.Loop("mc_test.root", "data_test.root", "sig_test.root")
