#! /opt/apps/python/2.7.1/bin/python

# This scripts checks over all of the directories to see which input files have been completed (by 
# looking to see if the last input file has been generated and is not empty) to build a list of 
# input files which need to be run.  It then divides that list into 50 sublists, and runs the 
# sublist specified by <job n>.


import os
import sys
from numpy import array, ceil
import threading
import Queue

directories = array([ '0-E30D_K31E+N27C','15-E30D+Y31C','20-K31E+N29C','37-Monomer+R20C','42-Monomer+K32C','7-WT+G28C','10-WT+K32C','16-E30D+K32C','21-K31E+Y31C','38-Monomer+N27C','43-Monomer+S33C','8-WT+N29C','11-WT+N54C','17-E30D+N54C','22-K31E+K32C','39-Monomer+G28C','44-Monomer+N54C','9-WT+Y31C','12-E30D+N27C','18-K31E+N27C','23-K31E+N54C','3-E30D_K31E+Y31C','4-E30D_K31E+K32C','13-E30D+G28C','19-K31E+G28C','2-E30D_K31E+N29C','40-Monomer+N29C','5-E30D_K31E+N54C','14-E30D+N29C','1-E30D_K31E+G28C','36-Monomer+I18C','41-Monomer+Y31C','6-WT+N27C' ])

#directories = array([ '0-E30D_K31E+N27C' ])

mg10_affix = array([ 'ncd_auto_2_10_193.txt', 'rcom_auto_78_10_193.txt', 'rcom_auto_2_10_193.txt', 'ncd_auto_78_10_193.txt' ])
mg19_affix = array([ 'ncd_auto_2_19_193.txt', 'rcom_auto_78_19_193.txt', 'rcom_auto_2_19_193.txt', 'ncd_auto_78_19_193.txt' ])
affix = array([ '2_60_257.txt', '2_10_193.txt', '2_19_193.txt', '2_10_193_1.txt', '2_19_193_1.txt', '2_10_193_0.txt', '2_19_193_0.txt', '2_30_129.txt' ])

R1S=0  
R1E=359   
R2S=0  
R2E=359   
FRAMES=80
cwd=os.getcwd()
exitFlag = 0

input_files = []

def get_tasks( threadName, q ) : 
	while not exitFlag : 
		queueLock.acquire()
		if not workQueue.empty() : 
			data = q.get()
			queueLock.release()
			print "%s processing %s" %(threadName,data)
			d=data[0]
			r1=data[1]
			r2=data[2]
			curdir="%s/%s/%i-%i"%(cwd,d,r1,r2)
			for frame in range(FRAMES+1) : 
				in1="%s/%s_%i-%i-frame%i_2.in" %(curdir,d,r1,r2,frame) 
		    		in2="%s/%s_%i-%i-frame%i_2_mg10.in" %(curdir,d,r1,r2,frame) 
		     		in3="%s/%s_%i-%i-frame%i_2_mg19.in" %(curdir,d,r1,r2,frame)
				if os.path.isfile(in1) and os.path.getsize(in1) != 0 :
					for i in affix : 
						out1="%s/%s_%i-%i-frame%i_%s"%(curdir,d,r1,r2,frame,i)
						if not os.path.isfile(out1) or os.path.getsize(out1) == 0 : 
							input_files.append(in1)
#							if i != affix[0] : 
#								print i,affix[0],out1
#								sys.exit()
							break
				if os.path.isfile(in2) and os.path.getsize(in2) != 0 :
					for i in mg10_affix : 
						out2="%s/%s_%i-%i-frame%i_%s"%(curdir,d,r1,r2,frame,i)
						#try : print out2, os.path.isfile(out2),os.path.getsize(out2)
						#except : print "\n"
						if not os.path.isfile(out2) or os.path.getsize(out2) == 0 : 
							input_files.append(in2)
#							if i != mg10_affix[0] : 
#								print i,mg10_affix[0],out2
#								sys.exit()
							break
				if os.path.isfile(in3) and os.path.getsize(in3) != 0 :
					for i in mg19_affix : 
						out3="%s/%s_%i-%i-frame%i_%s"%(curdir,d,r1,r2,frame,i)
						if not os.path.isfile(out3) or os.path.getsize(out3) == 0 : 
							input_files.append(in3)
#							if i != mg19_affix[0] : 
#								print i,mg19_affix[0],out3
#								sys.exit()
							break
		else : 
			queueLock.release()
#	input_files=array(input_files)
#	return input_files

def sort_tasks( i, tasks ) : 
	n=len(tasks)
	length=ceil(n/50.)
	sorted_tasks=tasks[length*i:length*(i+1)]
	return sorted_tasks

def do_apbs( input ) : 
	os.chdir(os.path.dirname(input))	
	cmd = "/home1/01635/awr397/apbs-1.3/bin/apbs %s" %input
	os.system(cmd)


TASKS = []
for d in directories : 
	for r1 in range(R1S,R1E,30) :
		for r2 in range(R2S,R2E,30) : 
			TASKS.append((d,r1,r2))

class myThread( threading.Thread ) : 
	def __init__(self, threadID, name, q) : 
		threading.Thread.__init__(self)
		self.threadID = threadID
		self.name = name
		self.q = q
	def run(self) : 
		print "Starting " + self.name
		get_tasks( self.name, self.q ) 
		print "Exiting " + self.name

threadList = []
for i in range(48) : 
	threadList.append("Thread-%i"%i)
queueLock = threading.Lock()
workQueue = Queue.Queue(0)
threads = []
threadID = 1

for tName in threadList : 
	thread = myThread(threadID,tName,workQueue)
	thread.start()
	threads.append(thread)
	threadID += 1

queueLock.acquire()
for task in TASKS : 
	workQueue.put(task)
queueLock.release()

while not workQueue.empty() :
	pass
exitFlag = 1
for t in threads : 
	t.join()
print "Exiting Main Thread"

input_files=array(input_files)
#filelist = open('inputlist.txt','w')
filelist = open('t.txt','w')
for input in input_files : filelist.write("%s\n"%input)
filelist.close()

for i in range(50) : 
	cmd="qsub ./run_apbs.py %i" %i
	os.system(cmd)

