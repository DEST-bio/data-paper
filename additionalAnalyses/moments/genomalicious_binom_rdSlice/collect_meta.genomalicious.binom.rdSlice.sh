
cat /scratch/aob2x/moments_binom_rdSlice/*.meta | awk '{
  split($2, sp, ".")
  print "binom_rdSlice."sp[3]"."$0
}' | sed 's/|/./g' > /project/berglandlab/moments/moments.genomalicious.binom.rdSlice.delim
