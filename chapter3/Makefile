#!/bin/bash
#need to generate 100 folders of GS, GS_trim, JsonD beforehand.

default:all
all:

.PHONY: para_space gill100 gsam100 rm_gap gimp100 plot_sim plot_im plot_sim_err plot_im_err

#command
RSCRIPT= Rscript --vanilla
RM     = rm -i

#path
SCR  = ../Script/chapter3

len        := 999
omega_z    := 1
sim_size   := 500
sample_size:= 1000
ncore      := 4
target     := NULL

#################################################100 simulations
#Pick a parameter space (100)
para_space:$(SCR)/para_space.R
	$(RSCRIPT) $<

#simulations
nameL  := $(shell test -f index.txt && cat index.txt)
gill100:  $(patsubst %, Gs/%, $(nameL))
gsam100:  $(patsubst %, Results/Gse/%.est.json, $(nameL))
gimp100:  $(patsubst %, Results/Pise/%.est.json,$(nameL))
Gs:          gill100
Results/Gse: gsam100
Results/Pise:gimp100

#100 gillespie simulations
Gs/%: $(SCR)/haha_gillespie100.R trueP.100.txt
	$(RSCRIPT) $< $@ $(word 2,$^) $(len) $(omega_z) $(sim_size)

#100 gillespie simu estimates
Results/Gse/%.est.json: $(SCR)/gillespie_sim_est2.R Gs/%
	$(RSCRIPT) $< $(word 2,$^) $@

#generate original seqs
rm_gap:runbash100.sh
	bash $< 100 500

#100/10k switch
Results/Pise/%.est.json: $(SCR)/phase_im_sampling_1e4.R Gs_trim/% JsonD/%/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $@ $(sample_size) 64


#plot the true vs sim estimates
plot_sim: $(SCR)/plot_trueVsSim.R
	$(RSCRIPT) $< $(target) 1

plot_sim_err: $(SCR)/plot_err.R
	$(RSCRIPT) $< $(target) 1

plot_im: $(SCR)/plot_trueVsSim.R
	$(RSCRIPT) $< $(target) 0

plot_im_err: $(SCR)/plot_err.R
	$(RSCRIPT) $< $(target) 0



help:
	@echo "make para_space:	  generate 100 parameter sets as our true values."	
	@echo "make gill100:	    run 100 gillespie simulations to generate 100x500 alignments."	
	@echo "make gsam100:	    run EM on generated alignments to obtain para estimates"
	@echo "make rm_gap:	      remove the gaps from simulated alignments"
	@echo "make gimp100:  	  run importance sampling algo on non-gap alignments to obtain para estimates"
	@echo "make plot_sim:  	  plot the true vs sim estimates" 
	@echo "make plot_im:  	  plot the true vs im estimates"
	@echo "make plot_sim_err: plot the error percentage of parameters from simulations"
	@echo "make plot_im_err:  plot the error percentage of parameters from sampling results"
	
	
clean:
	@rm -i Results/*
	@rm -i Figure/*
