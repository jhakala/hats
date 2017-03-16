
def pyMakeClass(options, arguments):
  from ROOT import TFile, TTree

  inFile = TFile(options.inputSample)
  fileKeyNames = []
  for key in inFile.GetListOfKeys():
    fileKeyNames.append(key.GetName())
    if options.debug:
      pprint(fileKeyNames)
  if "/" in options.treeName:
    if options.debug:
      print "looks like the tree name specified is within some TDirectory structure"
    treePath = options.treeName.split("/")
    treeName = treePath.pop()
    for tdirectory in treePath:
      for dirKeys in inFile.Get(tdirectory).GetListOfKeys():
        fileKeyNames.append(dirKeys.GetName())
      inFile.cd(tdirectory)
    if options.debug:
      print "found these keys in the TDirectory structure:"
      pprint(fileKeyNames)
      
  else:
    treeName = options.treeName
  
  if not treeName in fileKeyNames:
    print "Error: could not find a tree named %s in %s" % (options.treeName, options.inputSample)
    inFile.ls()
    exit(1)
  else:
    inTree = inFile.Get(options.treeName)
    inTree.MakeClass(options.className)
  

from pprint import pprint

if __name__ == "__main__":
  from os import path
  from optparse import OptionParser

  defaultSample = "sampleInput/flatTuple_999.root"
  defaultTreeName = "ntuplizer/tree"

  parser = OptionParser()
  parser.add_option("-i", "--inputSample", dest="inputSample",
                    default = defaultSample,
                    help = "the input root file to make a class for"  )
  parser.add_option("-t", "--treeName",   dest="treeName",
                    default = defaultTreeName,
                    help = "the name of the class that gets made"     )
  parser.add_option("-c", "--className",   dest="className",
                    default = "sampleFactory",
                    help = "the name of the class that gets made"     )
  parser.add_option("-d", "--debug",  dest="debug",
                    action="store_true", default=False,
                    help = "turn on debugging info"                   )
  (opts, args) = parser.parse_args()
  if opts.debug:
    print "options:" 
    pprint(opts)
    print "args:"
    pprint(args)
  
  
  if not path.exists(opts.inputSample):
    print "error: could not find input sample %s" % opts.inputSample 
    if opts.inputSample == defaultSample:
      print "you can specify a different sample from the default one using the -i option."
    exit(1)
  else:
    pyMakeClass(opts, args)
 
  
