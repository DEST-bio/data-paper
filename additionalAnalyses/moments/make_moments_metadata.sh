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
