#########################################################
#    _ _       __                __         		#
#   (_) |     / /___ _____  ____/ /__  _____		#
#  / /| | /| / / __ `/ __ \/ __  / _ \/ ___/		#
# / / | |/ |/ / /_/ / / / / /_/ /  __/ /    		#
#/_/  |__/|__/\__,_/_/ /_/\__,_/\___/_/     		#
# Dynamics of Interestellar Wanderers			#
# Jorge I. Zuluaga et al. [)] 2017			#
# http://github.com/seap-udea/iWander.git		#
#########################################################
# Makefile
#########################################################
include compiler.in
BRANCH=$(shell bash .getbranch)
PROGRAMS=wanderer encounters probability reconstruct

#########################################################
# MAIN RULES
#########################################################
test:
	@make .test.exe 
	@echo "It works!"

all:
	@echo "Compiling main programs..."
	for program in $(PROGRAMS);do make $$program.exe;done

analysis:all
	./wanderer.exe
	./encounters.exe
	./probability.exe
	$(PYTHON) bin/progenitors.py

edit:
	emacs -nw makefile *.cpp *.conf bin/*.py bin/*.sh *.sh

#########################################################
# COMPILATION
#########################################################
%.exe:%.opp
	$(CPP) $^ $(LFLAGS) -o $@

%.opp:%.cpp %.conf iwander.conf
	$(CPP) -c $(CFLAGS) $< -o $@

#########################################################
# INSTALL
#########################################################
pack:
	@echo "Packing data..."
	@bash .store/pack.sh pack

unpack:
	@echo "Unpacking data..."
	@bash .store/pack.sh unpack

#########################################################
# CLEAN
#########################################################
cleancrap:
	@echo "Cleaning crap..."
	@find . -name "*~" -exec rm -rf {} \;
	@find . -name "#*#" -exec rm -rf {} \;

cleanexe:
	@echo "Cleaning executable..."
	@rm -rf *.pyc *.out *.exe *.o *opp *.log .[a-zA-Z0-9]*.exe
	@rm -rf __pycache__

cleandata:
	@echo "Cleaning data..."
	@touch scratch/foo
	@-rm -r scratch/* 
	@touch log/foo
	@-rm -r log/* 

clean:cleancrap cleanexe
	@echo "Cleaning..."
	@rm -rf *.png *.dat

cleanall:clean cleandata

#########################################################
# GIT RELATED
#########################################################
merge:	
	@echo "Merging branches..."
	@-make master
	@-git merge dev

commit:
	@echo "Commiting..."
	@-git commit -am "Commit"
	@-git push origin $(BRANCH)

pull:
	@echo "Pulling latest version..."
	@-git reset --hard HEAD
	@-git pull origin $(BRANCH)

master:
	@-git checkout master

dev:
	@-git checkout dev

branch:
	@-echo $(BRANCH)

