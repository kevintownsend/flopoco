#testing Shifters
flopoco -target=Virtex4 -frequency=400 LeftShifter 40 5
flopoco -target=Virtex4 -frequency=400 LeftShifter 40 45
flopoco -target=Virtex4 -frequency=400 RightShifter 40 5
flopoco -target=Virtex4 -frequency=400 RightShifter 5 45


#testing LeadingZeroCounting
flopoco -target=Virtex4 -frequency=400 LZOC 55
flopoco -target=Virtex4 -frequency=400 LZOCShifter 79 45
flopoco -target=Virtex5 -frequency=500 LZCShifter 34 45
flopoco -target=Virtex5 -frequency=500 LOCShifter 12 8
flopoco -target=Virtex5 -frequency=500 LZOCShifterSticky 79 45
flopoco -target=Virtex4 -frequency=400 LZCShifterSticky 24 18
flopoco -target=Virtex4 -frequency=400 LOCShifterSticky 22 56


#Testing IntAdder
flopoco -target=Virtex4 -frequency=400 IntAdder 1
flopoco -target=Virtex4 -frequency=400 IntAdder 16
flopoco -target=Virtex4 -frequency=400 IntAdder 128
flopoco -target=Virtex5 -frequency=500 IntAdder 128


#flexible IntAdder
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 0 0 -1 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 1 0 -1 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 2 0 -1 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 3 0 -1 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 0 1 -1 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 1 1 -1 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 2 1 -1 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 3 1 -1 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 0 0 0 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 1 0 0 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 2 0 0 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 3 0 0 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 0 1 0 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 1 1 0 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 2 1 0 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 3 1 0 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 0 0 1 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 1 0 1 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 2 0 1 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 3 0 1 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 0 1 1 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 1 1 1 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 2 1 1 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 3 1 1 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 0 0 2 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 1 0 2 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 2 0 2 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 3 0 2 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 0 1 2 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 1 1 2 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 2 1 2 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 3 1 2 0
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 0 0 -1 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 1 0 -1 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 2 0 -1 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 3 0 -1 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 0 1 -1 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 1 1 -1 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 2 1 -1 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 3 1 -1 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 0 0 0 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 1 0 0 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 2 0 0 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 3 0 0 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 0 1 0 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 1 1 0 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 2 1 0 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 3 1 0 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 0 0 1 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 1 0 1 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 2 0 1 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 3 0 1 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 0 1 1 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 1 1 1 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 2 1 1 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 3 1 1 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 0 0 2 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 1 0 2 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 2 0 2 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 3 0 2 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 0 1 2 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 1 1 2 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 2 1 2 1
flopoco -target=Virtex4 -frequency=400 MyIntAdder 128 3 1 2 1

#intDualSub

flopoco -frequency=500 IntDualSub 26 0
flopoco -frequency=500 IntDualSub 26 1
flopoco -frequency=500 IntDualSub 216 0
flopoco -frequency=500 IntDualSub 216 1
flopoco -frequency=500 IntDualSub 1 0
flopoco -frequency=500 IntDualSub 1 1

#IntNAdder

flopoco -frequency=500 IntNAdder 1 1
flopoco -frequency=500 IntNAdder 1 2
flopoco -frequency=500 IntNAdder 10 1
flopoco -frequency=500 IntNAdder 100 1
flopoco -frequency=500 IntNAdder 100 2
flopoco -frequency=500 IntNAdder 100 10





