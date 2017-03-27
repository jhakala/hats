## This is an example of what you might do when you first look at some root files 
## that you've never seen before. First thing is to take a look at the basic
## structure of the file and the branches therein. Let's say the first thing
## we want to do is to look at our signal file.
##
## These will be like commands you run from an interactive pyroot session.

from ROOT import * 

firstLookFile = TFile("hatsSamples/sigFile1.root")
firstLookFile.ls()

## The output looks something like this. So now we can get the tree and look at it
####  TFile**    hatsSamples/sigFile1.root
####   TFile*    hatsSamples/sigFile1.root
####    KEY: TTree  tree;1  tree

hatsTree = firstLookFile.Get("tree")
for branch in hatsTree.GetListOfBranches():
  print branch.GetName()

## This shows us a bunch of the branch names like e.g.:
####  jetAK4_N
####  jetAK4_pt
####  jetAK4_eta
####  jetAK4_mass
####  jetAK4_phi
####  jetAK4_e
## now we might draw  the AK4 jet pT distribution

hatsTree.Draw("jetAK4_pt")

## Interesting! It looks like the AK4 jet pT spectrum has a peak around 400 in this signal sample
## 
## It is possible to loop over trees and compute more complicated quantities inside pyroot...
## This is covered in other exercises. 
##
## Here we will use a C++ class for number crunching, since it will save time.
## But we'll mostly just C++ classes as a plugin to fill the gaps where python slows us down 
## but leverage python for the myriad of tasks where it WILL save us lots of time.
##
## Let's call it "hatsTrees"

hatsTree.MakeClass("hatsTrees")

## The macro won't do anything for us now, so we need to edit it to compute interesting numbers
## In this repo there is an unmodified version of what MakeClass gives you
## You can see what I've added to the c++ class by running (in a shell)
##  vimdiff firstStepMacro.C hatsTrees.C
##   vimdiff firstStepMacro.h hatsTrees.h
