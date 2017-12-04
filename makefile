include compiler.in
BRANCH=$(shell bash .getbranch)
PROGRAMS=wanderer encounters probability

test:
	@make .test.exe 
	@echo "It works!"

branch:
	@echo $(BRANCH)

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
	@rm *.csv

clean:cleancrap cleanexe
	@echo "Cleaning..."
	@rm -rf *.png *.dat

cleanall:clean cleandata

all:
	@echo "Compiling main programs..."
	for program in $(PROGRAMS);do make $$program.exe;done

analysis:all
	./wanderer.exe
	./encounters.exe
	./probability.exe
	python3 progenitors.py

%.exe:%.opp
	$(CPP) $^ $(LFLAGS) -o $@

%.opp:%.cpp %.conf iwander.conf
	$(CPP) -c $(CFLAGS) $< -o $@

dev:
	@-git checkout dev

master:
	@-git checkout master

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

pack:
	@echo "Packing data..."
	@bash .store/pack.sh pack

unpack:
	@echo "Unpacking data..."
	@bash .store/pack.sh unpack
