RPG = newPop(founderPop)

dat <- pullQtlHaplo(RPG)
info <- hapsim::haplodata(dat)
info$freqs <- rep(start_allele_freq,n_chr) #frequencies for 0 allele

hap=info

haplos <- hapsim::haplosim(n_founders*2,hap)
haplos$freqs

for (each in 1:n_chr){
  assign(paste0("chr",each),matrix(as.integer(haplos$data[,each])))
}

haplotypes <- list(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10)

genMapRep = rep(list(seq(0,1,length.out=n_qtl)),n_chr)



