TESTSET=minlplib
SOLVER = WORHP
OPTFILE = 0
TIME = 30

TESTINSTANCES = $(shell cat $(TESTSET).test)

.PHONY: test
test : $(TESTSET)_$(SOLVER)_$(OPTFILE).trc

trc/%_$(SOLVER)_$(OPTFILE).trc :
	@mkdir -p trc lst log
	../gams/gams instances/$*.gms SOLVER=$(SOLVER) --MINLP=RMINLP --MIQCP=RMIQCP TRACEOPT=3 TRACE=$@ LO=3 O=lst/$*_$(SOLVER)_$(OPTFILE).lst RESLIM=$(TIME) OPTFILE=$(OPTFILE) > log/$*_$(SOLVER)_$(OPTFILE).log 2>&1

$(TESTSET)_$(SOLVER)_$(OPTFILE).trc : $(addprefix trc/, $(addsuffix _$(SOLVER)_$(OPTFILE).trc, $(TESTINSTANCES)))
	cat $^ > $@
