TESTSET=minlplib
SOLVER = WORHP
OPTFILE = 0
TIME = 30

TESTINSTANCES = $(shell cat $(TESTSET).test)

.PHONY: test
test : $(TESTSET)_$(SOLVER)_$(OPTFILE).trc

trc/%_$(SOLVER)_$(OPTFILE).trc : opt/%_$(SOLVER)_$(OPTFILE)/examiner2.opt
	@mkdir -p trc lst log
	../gams/gams instances/$*.gms SOLVER=EXAMINER2 --MINLP=RMINLP --MIQCP=RMIQCP LO=3 O=lst/$*_$(SOLVER)_$(OPTFILE).lst RESLIM=$(TIME) OPTFILE=1 OPTDIR=$(<D) #> log/$*_$(SOLVER)_$(OPTFILE).log 2>&1

opt/%_$(SOLVER)_$(OPTFILE)/examiner2.opt :
	@mkdir -p $(@D)
	@echo "subsolver $(SOLVER)" > $@
	@echo "subsolveropt $(OPTFILE)" >> $@
	@echo "trace $(subst opt/,trc/,$(@D)).trc" >> $@
	@echo "tracestyle 1" >> $@
	@echo "examineInitPoint 0" >> $@
	@echo "examineSolvPoint 1" >> $@

$(TESTSET)_$(SOLVER)_$(OPTFILE).trc : $(addprefix trc/, $(addsuffix _$(SOLVER)_$(OPTFILE).trc, $(TESTINSTANCES)))
	cat $^ > $@
