cd /project/berglandlab/moments/SNAPE/
ls * | sed 's/$/\t83960116/g' | awk '{
  n=split($0, sp, "/")
  x=split(sp[n], spp, ".")
  print spp[1]"."spp[2]".PoolSNP\t/project/berglandlab/moments/PoolSNP/"spp[1]"."spp[2]".unfolded.PoolSNP.fs\t"$2
  print spp[1]"."spp[2]".SNAPE\t/project/berglandlab/moments/SNAPE/"spp[1]"."spp[2]".unfolded.SNAPE.fs\t"$2
}' > /project/berglandlab/moments/moments_jobs.delim


cat /scratch/aob2x/pairs.csv | sed '1d' | sed 's/"//g' | awk -F',' '{
  print $5"."$4"."$3"\t/project/berglandlab/moments/"$3"/"$5"."$4".unfolded."$3".fs\t83960116"
}' > /project/berglandlab/moments/moments_jobs.delim









ls -d /project/berglandlab/moments/PoolSNP/* | sed 's/$/\t83960116/g' | awk '{
  n=split($0, sp, "/")
  x=split(sp[n], spp, ".")
  print spp[1]"."spp[2]"."spp[4]"\t"$0
}' > /project/berglandlab/moments/moments_jobs.delim

ls -d /project/berglandlab/moments/SNAPE/* | sed 's/$/\t83960116/g' | awk '{
  n=split($0, sp, "/")
  x=split(sp[n], spp, ".")
  print spp[1]"."spp[2]"."spp[4]"\t"$0
}' >> /project/berglandlab/moments/moments_jobs.delim





#### small test
ls -d /project/berglandlab/moments/PoolSNP/* | sed 's/$/\t83960116/g' | head | awk '{
  n=split($0, sp, "/")
  x=split(sp[n], spp, ".")
  print spp[1]"."spp[2]"."spp[4]"\t"$0
}' > /project/berglandlab/moments/moments_jobs.head.delim

ls -d /project/berglandlab/moments/SNAPE/* | sed 's/$/\t83960116/g' | head | awk '{
  n=split($0, sp, "/")
  x=split(sp[n], spp, ".")
  print spp[1]"."spp[2]"."spp[4]"\t"$0
}' >> /project/berglandlab/moments/moments_jobs.head.delim
