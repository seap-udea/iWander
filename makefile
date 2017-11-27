include compiler.in

test:
	@make .test.exe 
	@echo "It works!"

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

commit:
	@echo "Commiting..."
	@-git commit -am "Commit"
	@-git push origin master

pull:
	@-git reset --hard HEAD
	@-git pull

pack:
	@echo "Packing data..."
	@bash .store/pack.sh pack

unpack:
	@echo "Unpacking data..."
	@bash .store/pack.sh unpack
