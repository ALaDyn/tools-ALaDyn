#!/usr/bin/python
######################################################################
# Name:         ALaDyn_generates_vts.py
# Author:       A. Marocchino
# Date:			2014-02-18
# Purpose:      generates vts file for 3D plots
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, base64
import numpy as np
###>>>
###>>>
home_path = os.path.expanduser('~')
print home_path
sys.path.append(os.path.join(home_path,'Codes/ALaDyn_Code/tools-ALaDyn/ALaDyn_Pythons'))
###>>>
### --- ###
from read_ALaDyn_bin import *
from ALaDyn_plot_utilities_1 import *
from ALaDyn_write_vts import *
### --- ###



### --- ### shell inputs
if(len(sys.argv)<2):
	print 'Input [1]: frame_begin'
	print 'Input [2]: frame_end'
	print 'Input [3]: yz-shaving'

if sys.argv[1] == -1:
	frame_begin		  = 0
	frame_end         = last_output(os.getcwd())
	cell_cut		  = 0
else:
	frame_begin 		= int(		sys.argv[1])	
	frame_end			= int(		sys.argv[2])
	cell_cut		 	= int(		sys.argv[3])
### --- ###


#--- *** ---#
if __name__ == '__main__':
	
	#-path
	path = os.getcwd()
	
	#-folder output structure
	generate_folder_vts(path)


#	write_vts(path,2)

	frame =0; s='%2.2i'%frame
 	rhobunch, X,Y,Z	= read_ALaDyn_bin(path,'Bdenout'+s+'.bin','grid')

	print '--- Generating VTS ---'
	for i in range(frame_begin, min(frame_end,last_output(os.getcwd())) + 1 ):
		print '--- frame >',i
		write_vts(path,i,X,Y,Z)
		
# 
# # 	frame = 0
# # 	s='%2.2i'%frame 				#conversion to 2-character-long-string
# # 	file_name = 'Bdenout'+s+'.bin'
# # 	matrix,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
# # # 	file_name = 'Edenout'+s+'.bin'
# # # 	matrix2,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
# # 	#- cut & sign
# # 	matrix = np.abs( matrix )
# # # 	matrix2 = np.abs( matrix2 )
# # 	size = matrix.shape
# 	
# 	NI=2
# 	NJ=2
# 	NK=2
# 	size=[2,2,2]
# 
# 	
# 
# 	#- writing vts header
# 	f = open('test.vts','w+')
# 	f.write('<?xml version="1.0"?>' + '\n')
# 	f.write('<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">' + '\n')
# 	f.write('<StructuredGrid WholeExtent=" %d %d %d %d %d %d "> \n' % (0,size[1],0,size[0],0,size[2]) )
# 	f.write('<Piece Extent=" %d %d %d %d %d %d "> \n' % (0,size[1],0,size[0],0,size[2]) )
# 
# 
# 
# 	f.write('<Points> \n')
#  	f.write('<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="binary"> \n')
#  	buffer=[]
#  	for i in range(1,NI+2):
#  		for j in range(1,NJ+2):
#  			for k in range(1,NK+2):
#  				buffer.append(1./float(NJ)*(float(j)-1.))
#  				buffer.append(1./float(NI)*(float(i)-1.))
#  				buffer.append(1./float(NK)*(float(k)-1.))
# 
# 	s = base64.b64encode(np.array(buffer,dtype=np.float32))
#  	f.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
#  	f.write(s)
# 	f.write('</DataArray> \n')
# 	f.write('</Points> \n')
# 
# 
# 
# 
# 	f.write('<PointData>\n')
# 	f.write('<DataArray type="Float32" Name="rhoEdge" format="binary"> \n')
# 	buffer=[]
# 	n=0
# 	for i in range(0,NI+1):
# 		for j in range(0,NJ+1):
# 			for k in range(0,NK+1):
# 				buffer.append(n*1.)
# 				n+=1
# 
# 	s = base64.b64encode(np.array(buffer,dtype=np.float32))
#  	f.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
#  	f.write(s)
# # 	for item in buffer:
# # 		f.write(str(item)+'\n')	
# 	f.write('</DataArray> \n')
# 	f.write('</PointData> \n')
# 	
# 	
# 	
# 
# 	
# 	
# 	
# 	
# 	f.write('<CellData> \n')
# 	f.write('<DataArray type="Float32" Name="rho" format="ascii"> \n')
# 
# 	buffer=[]
# 	n=0
# 	for i in range(0,NI):
# 		for j in range(0,NJ):
# 			for k in range(0,NK):
# 				buffer.append(str(n*1.))
# 				n+=1
# 	
# 	for item in buffer:
# 		f.write(str(item)+'\n')
# 	
# 	f.write('</DataArray> \n')
# 
# 	f.write('</CellData> \n')
# 
# 	f.write('</Piece> \n')
# 	f.write('</StructuredGrid> \n')
# 	f.write('</VTKFile>')
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # 	#-writing MESH
# # 	f.write('<Points> \n')
# # 	f.write('<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii"> \n')
# # 
# # 	buffer=[]
# # 	n=1
# # 	for i in range(0,NI+1):
# # 		for j in range(0,NJ+1):
# # 			for k in range(0,NK+1):
# # 				buffer.append(1./NJ*(j-1))
# # 				buffer.append(1./NI*(i-1))
# # 				buffer.append(1./NK*(k-1))
# # 				n=n+3
# # 	n=n-1
# # # 	for k in z:
# # # 		for j in y:
# # # 			for i in x:
# # # 				f.write(str(i)+' '+str(j)+' '+str(k)+'\n')
# # 	for item in buffer:
# # 		f.write(str(item)+'\n')
# # 
# # 	f.write('</DataArray> \n')
# # 	f.write('</Points> \n')
# # 	f.write('<PointData></PointData> \n')
# # 
# # 	f.write('</Piece> \n')
# # 	f.write('</StructuredGrid> \n')
# # 	f.write('</VTKFile> \n')
# # 
# # 	
# # # 
# # # 	#-writing a scalar quantity: DENSITY
# # # 	f.write('<CellData>')
# # # 	f.write('<DataArray type="Float32" Name="rho" format="ascii">')
# # # 	f.write('</DataArray>')
# # # 
# # # 
# # # 	f.write('</DataArray>')
# # # 	
# # # 
# # # n=NI*NJ*NK
# # # call writeLine(lun,'<DataArray type="Float32" Name="Te" format="ascii">')
# # # call random_number(buffer(1:n))
# # # 
# # # i=1
# # # do while(i.le.n)
# # # 	write(str,'(10g)') buffer(i:min(i+9,n))
# # # 	call writeline(lun,str)
# # # 	i=i+10	
# # # enddo
# # # call writeLine(lun,'</DataArray>')
# # # 
# # # call writeLine(lun,'</CellData>')
# # # 
# # # !
# # # 
# # # !    <CellData Scalars="Te">
# # # !      <DataArray type="Float32" Name="Te" format="ascii" RangeMin="0.085676759481" RangeMax="0.79933160543">
# # # !        1 4 4 1
# # # !      </DataArray>
# # # !      <DataArray type="Float32" Name="Ti" format="ascii" RangeMin="0.085676759481" RangeMax="0.79933160543">
# # # !        4 1 1 4
# # # !      </DataArray>
# # # !    </CellData>
# # # !
# # # !write(str,'(a,i)') 'POINT_DATA',NNODES
# # # !call writeLine(lun,str)
# # # !call writeLine(lun,'SCALARS pippo float')
# # # !call writeLine(lun,'LOOKUP_TABLE default')
# # # !
# # # !do i=1,NI+1
# # # !do j=1,NJ+1
# # # !write(str,'(3g)') (j-1)*NI+(i-1)
# # # !call writeLine(lun,str)
# # # !enddo
# # # !enddo
# # # 
# # # call writeLine(lun,'</Piece>')
# # # call writeLine(lun,'</StructuredGrid>')
# # # call writeLine(lun,'</VTKFile>')
# # # 
# # # 
# # # 
# # # 
# # # 
# # # 
# # # 
# # # 	f.close()
# # # 
# # # 
# # # 
# # # # 	for i in range(frame_begin, min(frame_end,last_output(os.getcwd())) + 1 ):
# # # # 		print '-------------------'
# # # # 		
# # # # 		
# # # # 		if output_exists(path,'rho',i) == True:
# # # # 			print 'rho --- frame >>> ',i
# # # # 			plot_density_sections(path,i,rho_min,rho_max,isolines,celltocut,magnification_fig,savedata)
# # # # 
# # # # 		if output_exists(path,'E',i) == True:
# # # # 			print 'E --- frame >>> ',i
# # # # 			plot_Efield_sections(path,i,magnification_fig,savedata)
# # # # 
# # # # 			
# # # # 		if output_exists(path,'B',i) == True:
# # # # 			print 'B --- frame >>> ',i		
# # # # 			plot_Bfield_sections(path,i,magnification_fig,savedata)
# # # # 		
# # # # 		
# # # # 			
# # # # 		if output_exists(path,'Moving_window_axes',i) == True:
# # # # 			print 'Moving window axes --- frame >>> ',i
# # # # 			if savedata == 'True':
# # # #  				print 'Moving Window Coordinates data --- frame >>> ',i		
# # # #  				save_moving_window_coordinates(path,i)
# # # # 
# # # # 			
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
# # # # 
