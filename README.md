# Internal web app for Kunkels 

## Goal 

My goal is to build an app that we can deploy that automates the following tasks:

+ the creation of a set of oligos that make a set of mutants. The app will have to be given the coding nucleotide sequence, and a wild type PDB structure. 
+ the ordering of those mutants via Kunkel from Transcriptic via the Kunkel API
+ it will then create a sequencing order (since we can't use Transcriptic sequencing since it is too expensive) 
+ it will then ship the plate (TS API)
+ [how can it possibly submit an order to Genscript?]
+ submit a order to GenScript
+ retrieve sequencing results 
+ algorithmically compare those results to the wild type sequence 
+ collect all the mutations in a given clone into a list 
+ compare that list against the list that was ordered 
  + [submit the mutants that were not found through the whole thing again]
+ generate a map that gives the well in which to find each mutant 
+ print this map out and give to a human 
