#ABAQUS PDE - use the command print to print variables in ABAQUS
#
#Abaqus Python Library: 
from odbAccess import*
import os
import sys
import numpy as np
#Open Abaqus file: 
strWorkPath=os.getcwd()
odbPath = os.path.join(strWorkPath,r'E1.odb' )
odbFile = openOdb(path=odbPath)
#
#1ST TEXT FILE:
#
#What it does: chosen_node_stresses, chosen_node_disp1, chosen_node_disp2 are input parameters
# 		to decide which nodes we want to analyze to extract the stresses, and displacements. 
#
#Define input nodes:
chosen_node_stresses = 36007 #int(raw_input("choose a node for stresses:"))
chosen_node_disp1 = 1006	 #int(raw_input("choose the 1st node - disp:"))
chosen_node_disp2 = 1018     #int(raw_input("choose the 2nd node - disp:"))
inc_time=0 
if (chosen_node_disp1>chosen_node_disp2):
	print 'ERROR: node_disp1 label must be lesser than node_disp2 label'
	sys.exit()
nodes_disp = odbFile.rootAssembly.instances['PART-1-1'].nodeSets['ALL'] 
outputFile1 = open(os.path.join(strWorkPath, r'%s%s%s' % ('Node_data_', os.path.basename(os.path.normpath(strWorkPath)),'.dat')),'w')
outputFile1.write('Node for stresses: %d Nodes chosen to calculate the sheath radius: %d %d \n' %(chosen_node_stresses,chosen_node_disp1,chosen_node_disp2) )
outputFile1.write('time Mises Press Dispxx%d Dispxx%d \n' %(chosen_node_disp1,chosen_node_disp2))
for iStep in odbFile.steps.values(): 
	for iFrame in iStep.frames:	
		time = iFrame.frameValue+inc_time 
		nodes_strain = iFrame.fieldOutputs['E'].getSubset(position=ELEMENT_NODAL, region=nodes_disp) 
		nodes_stress = iFrame.fieldOutputs['S'].getSubset(position=ELEMENT_NODAL, region=nodes_disp) 
		nodes_displacement = iFrame.fieldOutputs['U']
		outputFile1.write('%10.4E ' % (time))
		for ivalue in nodes_stress.values: 
			if(chosen_node_stresses==ivalue.nodeLabel): 
				sig1=ivalue.data[0]
				sig2=ivalue.data[1]
				Mises=ivalue.mises
				Press=ivalue.press
				outputFile1.write('%10.4E %10.4E ' % (Mises,Press))
				break		
		for n,ivalue in enumerate(nodes_displacement.values): 
			if(chosen_node_disp1==ivalue.nodeLabel):
				initial_coordinate=nodes_disp.nodes[n].coordinates[0]
				Uxx=ivalue.data[0]
				outputFile1.write('%10.6E ' %(Uxx+initial_coordinate))
			if(chosen_node_disp2==ivalue.nodeLabel):
				initial_coordinate=nodes_disp.nodes[n].coordinates[0]
				Uxx=ivalue.data[0]
				outputFile1.write('%10.6E \n' %(Uxx+initial_coordinate))
				break
	inc_time=iFrame.frameValue+inc_time
print '1st txt file COMPLETED'
outputFile1.close()  
#
#2ST TEXT FILE:
#
#Important: Disregard the first line of code. Each point has two stresses assigned -- one for each element. 
#	Although the first node appears with two values, only the second one counts. 
#
#What it does: lFrame decides which steps we wanna plot. last_step defines the step. 
#	For each frame, it writes the stresses or displacement (or whatever is chosen)
#	along a path chosen by nodes_disp. In this example the 'BOTTOM_NODES' was the node
# 	set chosen. The 1st row gives the node number, the 2nd the coordinates, the 3rd on 
#	gives the stresses etc for each time frame.
#
# Header:
# Nodel Label; Node Coordinate; Pressure, Mises etc. for each frame.  
#
#Define variables - nodes of interest, step of interest, number of frames in the last step, Frames, Frames of interest:
nodes_disp = odbFile.rootAssembly.instances['PART-1-1'].nodeSets['BOTTOM_NODES'] 
Step=odbFile.steps.keys()[-1]
Number_frames= len(odbFile.steps[Step].frames)
Frame=odbFile.steps[Step].frames 
lFrame=[Frame[1],Frame[int(Number_frames/2.5)],Frame[int(Number_frames/2.0)],Frame[-1]]
#Iniciate Outputs
outputFile2 = open(os.path.join(strWorkPath, r'%s%s%s' % ('Path_data_', os.path.basename(os.path.normpath(strWorkPath)),'.dat')),'w')
#Header
outputFile2.write('Node_Number Node_Coordinate ')  
for i,iFrame in enumerate(lFrame):
	outputFile2.write('t=%0.1f ' %(lFrame[i].frameValue))
outputFile2.write(' \n')
#Iniciate Lists
for n,iNode in enumerate(nodes_disp.nodes):
	label=iNode.label
	initial_coordinate=iNode.coordinates[0]
	for j in xrange(1,-1,-1):
		outputFile2.write('%4d %10.4E ' % (label,initial_coordinate))
		for iFrame in lFrame:
			Stress=iFrame.fieldOutputs['S'].getSubset(position=ELEMENT_NODAL, region=nodes_disp).values[2*n-j]	
			outputFile2.write('%10.4E %10.4E %10.4E %10.4E %10.4E %10.4E ' % (Stress.mises,Stress.press, Stress.data[0],Stress.data[1],Stress.data[2],Stress.data[3]))
		outputFile2.write('\n')
print '2nd txt file COMPLETED'
outputFile2.close()