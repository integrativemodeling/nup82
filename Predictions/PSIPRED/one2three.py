
f=open("oneletter.Nsp1.txt","r")
for l in f:
 for aa in l:
  if "A"==aa : print "ALA"
  elif "C"==aa : print "CYS"
  elif "D"==aa : print "ASP"
  elif "E"==aa : print "GLU"
  elif "F"==aa : print "PHE"
  elif "G"==aa : print "GLY"
  elif "H"==aa : print "HIS"
  elif "I"==aa : print "ILE"
  elif "K"==aa : print "LYS"
  elif "L"==aa : print "LEU"

  elif "N"==aa : print "ASN"
  elif "M"==aa : print "MET"
  elif "P"==aa : print "PRO"
  elif "Q"==aa : print "GLN"
  elif "R"==aa : print "ARG"
  elif "S"==aa : print "SER"
  elif "T"==aa : print "THR"
  elif "V"==aa : print "VAL"
  elif "W"==aa : print "TRP"
  elif "Y"==aa : print "TYR"
  
