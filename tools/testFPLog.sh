
tests=1000


for e in 8  11 ; do 
  for f in 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 ; do

target=Virtex4
name=FPLog_$e'_'$f'_'$target
cd /home/fdedinec/Boulot/Recherche/Hard/FloPoCo.bibine/trunk
./flopoco -target=$target -name=$name FPLog $e $f TestBench $tests  
 mv flopoco.vhdl ghdl
cd ghdl
echo
echo Testing $e $f for $target
   ghdl -a --ieee=synopsys -fexplicit flopoco.vhdl
   ghdl -e --ieee=synopsys -fexplicit TestBench_$name
   ghdl -r --ieee=synopsys TestBench_$name --vcd=TestBench_$name.vcd

target=StratixII
name=FPLog_$e'_'$f'_'$target
cd /home/fdedinec/Boulot/Recherche/Hard/FloPoCo.bibine/trunk
./flopoco -target=$target -name=$name FPLog $e $f TestBench $tests  
 mv flopoco.vhdl ghdl
cd ghdl
echo
echo Testing $e $f for $target
   ghdl -a --ieee=synopsys -fexplicit flopoco.vhdl
   ghdl -e --ieee=synopsys -fexplicit TestBench_$name
   ghdl -r --ieee=synopsys TestBench_$name --vcd=TestBench_$name.vcd


  done ;
done

