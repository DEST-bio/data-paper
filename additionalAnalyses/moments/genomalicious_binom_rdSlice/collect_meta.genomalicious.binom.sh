
cat /scratch/aob2x/moments_binom/*.meta | awk '{
  split($2, sp, ".")
  print "binom."sp[3]"."$0
}' | sed 's/|/./g' > /project/berglandlab/moments/moments.genomalicious.binom.delim
