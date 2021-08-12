cat /scratch/aob2x/moments/*.meta | awk '{
  split($2, sp, ".")
  print "nEff."sp[3]"."$0
}' | sed 's/|/./g' > /project/berglandlab/moments/moments.genomalicious.delim
