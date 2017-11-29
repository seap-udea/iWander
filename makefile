include compiler.in
BRANCH=$(shell bash .getbranch)

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
	@rm -rf *.pyc *.out *.{exe,o,opp} *.log .[a-zA-Z0-9]*.{exe,o,opp}

clean:cleancrap cleanexe
	@echo "Cleaning..."
	@rm -rf *.png *.dat

%.exe:%.opp
	$(CPP) $^ $(LFLAGS) -o $@

%.opp:%.cpp %.conf
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
