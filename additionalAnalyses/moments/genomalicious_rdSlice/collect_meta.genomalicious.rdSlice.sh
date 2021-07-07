cat /scratch/aob2x/moments_rdSlice/*.meta | awk '{
  split($2, sp, ".")
  print "counts.rdSlice."sp[3]"."$0
}' | sed 's/|/./g' > /project/berglandlab/moments/moments.genomalicious.rdSlice.delim
