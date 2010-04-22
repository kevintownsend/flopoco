import random
import math
import os
import subprocess
import string
import commands
import pipes
import time
#-------------------------------------------------------------------------------

#function declarations
#prints the operator parameters and their types nicely
def printParameters(list):
	i=0
	print "      Name|      Type"
	print "---------------------"
	while (i)<len(list):
		print '%10s|%10s' % (list[i][0], list[i][1]) 
		i=i+1
	return

#prints the operator parameters and their types nicely
def printResults(list, lfile):
	i=0
	print '%93s|%12s|%11s|%4s' % ("TEST", "GENERATE_HDL", "COMPILE HDL", "PASS") 
	lfile.write( '%93s|%12s|%11s|%4s\n' % ("TEST", "GENERATE_HDL", "COMPILE HDL", "PASS")) 
	print "----------------------------------------------------------------------------------------------------"
	lfile.write("----------------------------------------------------------------------------------------------------\n")
	while (i)<len(list):
		print '%93s|%12s|%11s|%4s' % (list[i][0], list[i][1], list[i][2], list[i][3]) 
		lfile.write('%93s|%12s|%11s|%4s' % (list[i][0], list[i][1], list[i][2], list[i][3]))
		lfile.write("\n")
		i=i+1
	return

#maps a list of random numbers to a list of operator parameters
#return a list containing tuples of the form (name, random_value)
def mapRandoms(paramList):
	i=0
	randomMapList=[]
	while (i)<len(paramList):
		if paramList[i][2] != paramList[i][3]:
			randomMapList.append(random.randrange(paramList[i][2],paramList[i][3],1))
		else:
			randomMapList.append(paramList[i][3])
		i=i+1
	return randomMapList

#nicely prints the random mappings on the parameters without their paranthesis
def prinRandomMappings(paramList, RandomList):
	i=0
	randomMapListPrint=""
	while (i)<len(paramList):
		randomMapListPrint=randomMapListPrint + " " + paramList[i][0] + " " +`RandomList[i]`
		i=i+1
	return randomMapListPrint

#writes the random numbers in a string ready appended to the operator description
def cmdLine(paramList):
	cmd_list=""
	i=0
	while (i)<len(paramList):
		cmd_list = cmd_list + " " + `paramList[i]`
		i=i+1
	return cmd_list
#-------------------------------------------------------------------------------



print "----------------------------------------------------------------------------------------------------"
print "--------------------------- FloPoCo PreRelease Tests -----------------------------------------------"
executable_name = "./flopoco"
operators = [
#             ["LeftShifter",       [ ["wIn"   ,"in",   4  , 120] ,["MaxShift"  ,"in",   2, 64 ]  ]], 
#             ["RightShifter",      [ ["wIn"   ,"in",   4  , 120] ,["MaxShift"  ,"in",   2, 64 ]  ]], 
#             ["LZOC",              [ ["wIn"   ,"in",   4  , 120]                                 ]], 
#             ["LZOCShifter",       [ ["wIn"   ,"in",   4  , 120], ["wOut"      ,"in",   4, 120]  ]], 
#             ["LZCShifter",        [ ["wIn"   ,"in",   4  , 120], ["wOut"      ,"in",   4, 120]  ]], 
#             ["LOCShifter",        [ ["wIn"   ,"in",   4  , 120], ["wOut"      ,"in",   4, 120]  ]],    
#             ["LZOCShifterSticky", [ ["wIn"   ,"in",   4  , 120], ["wOut"      ,"in",   4, 120]  ]], 
#             ["LZCShifterSticky",  [ ["wIn"   ,"in",   4  , 120], ["wOut"      ,"in",   4, 120]  ]],  
#             ["LOCShifterSticky",  [ ["wIn"   ,"in",   4  , 120], ["wOut"      ,"in",   4, 120]  ]],  
#             ["IntAdder",          [ ["wIn"   ,"in",   1  ,  64]                                 ]],   
#             ["IntNAdder",         [ ["wIn"   ,"in",   1  , 120], ["N"         ,"in",   1,  30]  ]],
#             ["IntCompressorTree", [ ["wIn"   ,"in",   1  , 120], ["N"         ,"in",   1,  30]  ]],      
#             ["IntMultiplier",     [ ["wInX"  ,"in",   1  ,  64], ["wInY"      ,"in",   1,  64]  ]],  
#             ["IntKaratsuba",      [ ["wIn"   ,"in",  18  ,  51]                                 ]], 
#             ["IntSquarer",        [ ["wIn"   ,"in",   2  ,  68]                                 ]],    
#             ["Fix2FP",            [ ["LSB"   ,"in",-100  , 100], ["MSB"       ,"in",-100, 100], ["Signed","in", 0, 1], ["wE","in",4, 15], ["wF","in",8, 64] ]], 
#			 ["FPAdder",           [ ["wE"    ,"in",   4  ,  15], ["wF"        ,"in",   8,  64]  ]],
#			 ["FPMultiplier",      [ ["wE"    ,"in",   4  ,  15], ["wF"        ,"in",   8,  64]  ]],
#			 ["FPDiv",             [ ["wE"    ,"in",   4  ,  15], ["wF"        ,"in",   8,  64]  ]],
#			 ["FPExp",             [ ["wE"    ,"in",   4  ,  15], ["wF"        ,"in",   8,  64]  ]], 
#             ["FPLog",             [ ["wE"    ,"in",   6  ,  11], ["wF"        ,"in",   8,  52]  ]], 
#             ["FPSqrt",            [ ["wE"    ,"in",   8  ,   8], ["wF"        ,"in",  23,  23], ["useDSP","in", 1, 1], ["correctlyRounded","in", 0, 1] ]],
#             ["FPSqrt",            [ ["wE"    ,"in",   4  ,  15], ["wF"        ,"in",   8,  64], ["useDSP","in", 0, 0], ["correctlyRounded","in", 1, 1] ]], 
#             ["IntConstMult",      [ ["w"     ,"in",   2  , 120], ["c"         ,"in",   0,32000] ]], 
#             ["FPConstMult",       [ ["wE_in" ,"in",   4  ,  15], ["wF_in"     ,"in",   4,  64], ["wE_out","out",4, 15], ["wF_out","out",4, 64], ["cst_sgn","in",0, 1], ["cst_exp","in",-128, 128], ["cst_int_sig","in",0, 32000] ]],
#             ["LongAcc",           [ ["wE"    ,"in",   4  ,  15], ["wF"        ,"in",   8,  64], ["MaxMSB","in", -100, 100], ["LSB","in", -100, 100], ["MSB","in", -100, 100]]],
#             ["OutputIEEE",        [ ["wE_in" ,"in"   ,4  ,  15], ["wF_in"     ,"in",   4,  64], ["wE_out","out",4, 15], ["wF_out","out",4, 64]]],
#             ["InputIEEE",         [ ["wE_in" ,"in"   ,4  ,  15], ["wF_in"     ,"in",   4,  64], ["wE_out","out",4, 15], ["wF_out","out",4, 64]] ]
            ] #TODO Add the rest of operators
res = []
input_combinations_per_operator = 1;
test_cases_per_combination = 100;

#REMOVE TEMPORARY MODELSIM FILES
os.system("rm -f wlf*")

#CREATE A LOG FILE 
logfile = open( "release_test.log","w")

#logfile.write("Parsing operators from internal list: \n")
##parse operator list and test each operator
#i=0
#pass_all = True
#while i<len(operators):
#	j=input_combinations_per_operator
#	while j>0:
#		commands.getoutput("rm vsim*")
#		commands.getoutput("killall vsimk")
#		commands.getoutput("killall -9 vish")
#		print "--------------------------- ", operators[i][0]   ," ---------------------------"
#		print "Parameters: "
#		printParameters( operators[i][1] )
#		uut = operators[i][0]
#		run_cmd = executable_name + " " + uut + " " + cmdLine(mapRandoms(operators[i][1])) + " TestBench " + `test_cases_per_combination`
#		print run_cmd
#		print logfile.write(run_cmd + "\n")
#		modelsim_food = commands.getoutput(run_cmd)
#		print modelsim_food
#		logfile.write(modelsim_food+"\n")

#		did_generate_vhdl = True
#		status = string.find(modelsim_food, "Pipeline depth = 42")
#		if status < 0:
#			did_generate_vhdl = False

#		modelsim_food = modelsim_food[string.find(modelsim_food, "vdel") : string.find(modelsim_food, "To run the simulation using gHDL,")-1 ]

#		finished = False
#		pass_test = True
#		did_compile = True

#		if did_generate_vhdl:
#		#start modelsim
#			p = subprocess.Popen("vsim -c", shell=True, bufsize=1, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
#			(child_stdin, child_stdout, child_stderr) = (p.stdin, p.stdout, p.stderr)
#			child_stdin.write("vdel -all -lib work\n")
#			child_stdin.write('vlib work\n')
#			child_stdin.write('vcom flopoco.vhdl \n')
#			child_stdin.write( modelsim_food[string.find(modelsim_food,"vsim"):string.find(modelsim_food,"add")+1]+"\n" )
#			child_stdin.write('add wave -r * \n')
#			child_stdin.write(modelsim_food[string.find(modelsim_food,"run"):]+"\n")
#			child_stdin.write('exit \n')
#			#-------------------------------------------------------------------		
#			while ((not finished) and (did_compile)):
#				st = child_stdout.readline()
#				print st[:len(st)-1]
#				logfile.write(st[:len(st)-1]+"\n")

#				status = string.find(st, "Error:")
#				if status > 0:
#					did_compile = False
#					pass_test = False

#				status = string.find(st, "Incorrect")
#				if status > 0:
#					pass_test = False

#				status = string.find(st, "Stopped")
#				if status > 0:
#					finished = True
#					did_compile = False

#				status = string.find(st, "Error loading design")
#				if status > 0:
#					finished = True

#			#-------------------------------------------------------------------		
#			if did_compile:
#				child_stdin.write("exit \n")

#			child_stdout.close()
#			child_stdin.close()
#			child_stderr.close()
#			#p.wait()
#			if pass_test:
#				print "Test: ", run_cmd , " has SUCCEDED "
#			else:
#				print "Test: ", run_cmd , " has FAILED "
#		
#		pass_test = pass_test and did_generate_vhdl
#		res.append( [run_cmd,`pass_test`])
#		pass_all = pass_all and pass_test
#		j = j -1
#	i = i+1


pass_all = True
#PARSE EXTERNAL TEST FILE
print("Parsing operators from external list: \n")
logfile.write("Parsing operators from external list: \n")

fd = open ("flopoco_test.cmd")
if fd < 0:
	print ("Unable to open file flopoco_test.cmd")
	logfile.write("Unable to open file flopoco_test.cmd")
for line in fd:
	if ((line[0]!='#') and (len(line)>1)):
		commands.getoutput("rm -f vsim*")
		commands.getoutput("killall -9 vsimk")
		run_cmd = line[:len(line)-1] + " TestBench " + `test_cases_per_combination`
		print run_cmd
		logfile.write(run_cmd+"\n")
		modelsim_food = commands.getoutput(run_cmd)

		did_generate_vhdl = True
		status = string.find(modelsim_food, "Output file: flopoco.vhdl")
		if status < 0:
			did_generate_vhdl = False
			print("Did not generate VHDL");

		modelsim_food = modelsim_food[string.find(modelsim_food, "vdel") : string.find(modelsim_food, "To run the simulation using gHDL,")-1 ]

		print("The modelsim food is:\n"+modelsim_food);

		finished = False
		pass_test = True
		did_compile = True

		if did_generate_vhdl:
			print("It did generate VHDL");
			#start modelsim
			p = subprocess.Popen("vsim -c", shell=True, bufsize=1, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
			(child_stdin, child_stdout, child_stderr) = (p.stdin, p.stdout, p.stderr)
			child_stdin.write("vdel -all -lib work\n")
			child_stdin.write('vlib work\n')
			child_stdin.write('vcom flopoco.vhdl \n')
			child_stdin.write('vcom flopoco.vhdl \n')
			child_stdin.write( modelsim_food[string.find(modelsim_food,"vsim"):string.find(modelsim_food,"add")-1]+"\n" )
			child_stdin.write('add wave -r * \n')
			child_stdin.write(modelsim_food[string.find(modelsim_food,"run"):]+"\n")
			child_stdin.write('exit \n')
	
			while ((not finished) and (did_compile)):
				st = child_stdout.readline()
				#print st[0:len(st)-2]
				logfile.write(st[0:len(st)-2]+"\n")
				status = string.find(st, "Error:")
				if status > 0:
					pass_test = False
#					did_compile = False
		
				status = string.find(st, "Incorrect")
				if status > 0:
					pass_test = False
		
				status = string.find(st, "Stopped")
				if status > 0:
					finished = True
#					did_compile = False
		
				status = string.find(st, "Error loading design")
				if status > 0:
					did_compile = False
					finished = True
	
			if did_compile:
				child_stdin.write('exit \n')

			child_stdout.close()
			child_stdin.close()
			child_stderr.close()
			#p.wait()
#		did_compile = not did_compile
		pass_test = pass_test and did_generate_vhdl
		res.append( [run_cmd, `did_generate_vhdl`, `did_compile`, `pass_test`])
		pass_all = pass_all and pass_test
	
		
		if pass_test:
			print "Test: ", run_cmd , " has SUCCEDED "
		else:
			print "Test: ", run_cmd , " has FAILED "
fd.close()


printResults( res , logfile)

print "FINAL PASS STATUS: ", `pass_all`
logfile.write("FINAL PASS STATUS:" +  `pass_all`)

logfile.close()
print "Logfile was written to: release_test.log"

