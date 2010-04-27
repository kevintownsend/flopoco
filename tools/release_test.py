import random
import math
import os
import subprocess
import string
import commands
import pipes
import time



input_combinations_per_operator = 1
test_cases_per_combination = 100
useModelSim=True # if True, use modelsim; if False, use ghdl


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


res = []
#REMOVE TEMPORARY MODELSIM FILES
os.system("rm -f wlf*")

#CREATE A LOG FILE 
logfile = open( "release_test.log","w")

pass_all = True
#PARSE EXTERNAL TEST FILE
print("Parsing operators from external list: \n")
logfile.write("Parsing operators from external list: \n")

for filename in os.listdir("tests/"):
	if (filename[len(filename)-1]!='~'):
		print "processing " + filename + " ... "
		
		fd = open ("tests/"+filename)
		if fd < 0:
			print "Unable to open file ".filename
			logfile.write("Unable to open file " + filename)

		skipfile = False
		for line in fd:
			if (line[0]=='!'):
				skipfile = True
				print "Skipping the rest of the file "+filename+"; ! encountered"
				
			if ((line[0]!='#') and (len(line)>1) and (skipfile == False)):
				run_cmd = line[:len(line)-1] + " TestBench " + `test_cases_per_combination`
				print run_cmd
				logfile.write(run_cmd+"\n")
				modelsim_food = commands.getoutput(run_cmd)
				did_generate_vhdl = True
				status = string.find(modelsim_food, "Output file: flopoco.vhdl")

				if status < 0:
					did_generate_vhdl = False
					print("Did not generate VHDL");

				if useModelSim:
					commands.getoutput("rm -f vsim*")
					commands.getoutput("killall -9 vsimk")
					modelsim_food = modelsim_food[string.find(modelsim_food, "vdel") : string.find(modelsim_food, "To run the simulation using gHDL,")-1 ]
				else:
					ghdl_food = modelsim_food[string.find(modelsim_food, "   ghdl") : string.find(modelsim_food, "Final")-1 ]
	
				finished = False
				pass_test = True
				did_compile = True

				if did_generate_vhdl:
					print("It did generate VHDL");

					if useModelSim:
					#start modelsim
						p = subprocess.Popen("vsim -c", shell=True, bufsize=1, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
						(child_stdin, child_stdout, child_stderr) = (p.stdin, p.stdout, p.stderr)
						child_stdin.write("vdel -all -lib work\n")
						child_stdin.write('vlib work\n')
						child_stdin.write('vcom flopoco.vhdl \n')
						child_stdin.write('vcom flopoco.vhdl \n')
						child_stdin.write(modelsim_food[string.find(modelsim_food,"vsim"):string.find(modelsim_food,"add")-1]+"\n" )
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
								#did_compile = False
						
								status = string.find(st, "Incorrect")
								if status > 0:
									pass_test = False
		
							status = string.find(st, "Stopped")
							if status > 0:
								finished = True
								#did_compile = False
		
								status = string.find(st, "Error loading design")
								if status > 0:
									did_compile = False
									finished = True

						if did_compile:
							child_stdin.write('exit \n')
							child_stdout.close()
							child_stdin.close()
							child_stderr.close()

					else: # use ghdl

						cmd=ghdl_food[string.find(ghdl_food,"ghdl -a"):string.find(ghdl_food,".vhdl")+5]
						logfile.write(cmd+"\n")
						print(cmd)
		 				status=commands.getoutput(cmd)
		 				if(status):
							print "ghdl -a error:"
							print status
							logfile.write(status+"\n")
							did_compile=False
							pass_test=False

						cmd=ghdl_food[string.find(ghdl_food,"ghdl -e"):string.find(ghdl_food,"   ghdl -r")-1]
						logfile.write(cmd+"\n")
						print(cmd)
		 				status=commands.getoutput(cmd)
		 				if(status):
							print "ghdl -e error:"
							print status
							logfile.write(status+"\n")
							pass_test=False
							did_compile=False

						cmd=ghdl_food[string.find(ghdl_food,"ghdl -r"):string.find(ghdl_food,".vcd")+4]
						logfile.write(cmd+"\n")
						print(cmd)
		 				status=commands.getoutput(cmd)
		 				if string.find(status, "Incorrect") !=-1:
							print "ghdl -r error:"
							print status
							logfile.write(status+"\n")
							pass_test=False
						commands.getoutput("rm *.vcd e~testbench* testbench*")

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

