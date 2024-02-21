GEFLASH.exe: EqOfSt.o  mixRuleSubroutines.o  Input.o main.o  RachfordRice.o NewtonRaphson.o FuelVap.o
	g++ -g EqOfSt.o  mixRuleSubroutines.o  Input.o main.o   RachfordRice.o  NewtonRaphson.o FuelVap.o  -o GEFLASH.exe

NewtonRaphson.o: NewtonRaphson.C
	g++ -c -g NewtonRaphson.C

EqOfSt.o: EqOfSt.C 
	g++ -c -g EqOfSt.C  

mixRuleSubroutines.o: mixRuleSubroutines.C 
	g++ -c -g mixRuleSubroutines.C  

Input.o: Input.C 
	g++ -c -g Input.C

RachfordRice.o: RachfordRice.C
	g++ -c -g RachfordRice.C

FuelVap.o: FuelVap.C
	g++ -c -g FuelVap.C

main.o: main.C  
	g++ -c -g main.C


clean:
	rm -f PHVap.exe *.o  

 
