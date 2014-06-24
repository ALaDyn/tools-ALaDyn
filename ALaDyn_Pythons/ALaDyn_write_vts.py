#!/usr/bin/python
######################################################################
# Name:         ALaDyn_write_vts.py
# Author:       A. Marocchino
# Date:			2014-02-18
# Purpose:      write the VTS file
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, base64
import numpy as np
###>>>
###>>>
home_path = os.path.expanduser('~')
sys.path.append(os.path.join(home_path,'Codes/ALaDyn_Code/tools-ALaDyn/ALaDyn_Pythons'))
###>>>
### --- ###
from read_ALaDyn_bin import *
### --- ###


def matrix2vector( M ):
	s = M.shape
	v = []
	if len(s) == 3:
		for i in range(0,s[0]):
			for j in range(0,s[1]):
				for k in range(0,s[2]):
					v.append(M[i,j,k])
	if len(s) == 2:
		for i in range(0,s[0]):
			for j in range(0,s[1]):
				v.append(M[i,j])	
	return v

def matrix2vectorField( M1, M2, M3):
	s = M1.shape
	v = []
	if len(s) == 3:
		for i in range(0,s[0]):
			for j in range(0,s[1]):
				for k in range(0,s[2]):
					v.append(M1[i,j,k])
					v.append(M2[i,j,k])
					v.append(M3[i,j,k])
	if len(s) == 2:
		for i in range(0,s[0]):
			for j in range(0,s[1]):
				v.append(M1[i,j])
				v.append(M2[i,j])
				v.append(M3[i,j])
	return v


#--- *** ---#
def write_vts(path,frame,X,Y,Z,cell_cut):
# 	for i in range(frame_begin, min(frame_end,last_output(os.getcwd())) + 1 ):
# 		write_VTS(path,i)

 	sf='%2.2i'%frame 				#conversion to 2-character-long-string
 	file_name 		= 'Bdenout'+sf+'.bin'
 	rhobunch, x,y,z	= read_ALaDyn_bin(path,file_name,'grid')
 	rhobunch		= np.abs( rhobunch )
 	file_name 		= 'Edenout'+sf+'.bin'
 	rhobck, x,y,z	= read_ALaDyn_bin(path,file_name,'grid')
 	rhobck			= np.abs( rhobck )
 	
 	#- matrix shaving
 	if cell_cut > 0:
	 	rhobunch = rhobunch[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
 		rhobck   = rhobck[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
	 	Y		 = Y[cell_cut:-cell_cut]
 		Z        = Z[cell_cut:-cell_cut]
 	#-
 	size 			= rhobunch.shape
 	#- 

 	
	#- writing vts header
	f = open(os.path.join(path,'VTS_files','ALaDyn_output_'+sf+'.vts'),'w+')
	f.write('<?xml version="1.0"?>' + '\n')
	f.write('<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">' + '\n')
	f.write('<StructuredGrid WholeExtent=" %d %d %d %d %d %d "> \n' % (0,size[1]-1,0,size[2]-1,0,size[0]-1) )
	f.write('<Piece Extent=" %d %d %d %d %d %d "> \n' % (0,size[1]-1,0,size[2]-1,0,size[0]-1) )
	#--- ---#


	#- generating vector-MESH
 	mesh=[]
 	for i in range(0,size[0]):
 		for j in range(0,size[1]):
 			for k in range(0,size[2]):
				mesh.append( Y[j] )
				mesh.append( Z[k] )
				mesh.append( X[i] )
	#- Writing MESH -#
	f.write('<Points> \n')
 	f.write('<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="binary"> \n')
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f.write(  base64.b64encode(np.array(len(s),dtype=np.int32))  )
 	f.write(  s  )
	f.write('</DataArray> \n')
	f.write('</Points> \n')
	#- -#

	#- Point Data Begin-#
 	f.write('<PointData>\n')
#	f.write('<CellData> \n')
	#- -#

	#- Writing BUNCH Density -#
	f.write('<DataArray type="Float32" Name="rho_bunch" format="binary"> \n')
	mesh = matrix2vector( rhobunch )
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	f.write(s)
	f.write('</DataArray> \n')

	#- Writing Tot-Density -#
	f.write('<DataArray type="Float32" Name="rho" format="binary"> \n')
	mesh = matrix2vector( rhobunch+rhobck )
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	f.write(s)
	f.write('</DataArray> \n')
		#-Deallocate memory
	rhobunch = [0.]; rhobck = [0.]

	#- Writing E-field -#
			#- background -#
 	file_name = 'Exfout'+sf+'.bin'
	Ex,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Eyfout'+sf+'.bin'
	Ey,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Ezfout'+sf+'.bin'
	Ez,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
			#- bunch(es) -#
	file_name = 'Exbout'+sf+'.bin'
	Exb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Eybout'+sf+'.bin'
	Eyb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Ezbout'+sf+'.bin'
	Ezb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')

 	#- matrix shaving
 	if cell_cut > 0:
	 	Ex  = Ex[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
	 	Ey  = Ey[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
	 	Ez  = Ez[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
	 	Exb = Exb[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
	 	Eyb = Eyb[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
	 	Ezb = Ezb[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
 	#- 

	
	#- E-field bunch
	f.write('<DataArray type="Float32" Name="E_bunch" NumberOfComponents="3" format="binary"> \n')
	mesh = matrix2vectorField(Exb,Eyb,Ezb)
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	f.write(s)
	f.write('</DataArray> \n')
	#- E-field total
	f.write('<DataArray type="Float32" Name="E" NumberOfComponents="3" format="binary"> \n')
	mesh = matrix2vectorField(Ex+Exb,Ey+Eyb,Ez+Ezb)
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	f.write(s)
	f.write('</DataArray> \n')
		#-Deallocate memory
	Ex=[0.];Ey=[0.];Ez=[0.];Exb=[0.];Eyb=[0.];Ezb=[0.];


	#- Writing B-field -#
			#- background -#
 	file_name = 'Bxfout'+sf+'.bin'
	Bx,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Byfout'+sf+'.bin'
	By,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Bzfout'+sf+'.bin'
	Bz,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
# 			#- bunch(es) -#
# 	file_name = 'Bxbout'+sf+'.bin'
# 	Bxb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
# 	file_name = 'Bybout'+sf+'.bin'
# 	Byb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
# 	file_name = 'Bzbout'+sf+'.bin'
# 	Bzb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')

 	#- matrix shaving
 	if cell_cut > 0:
	 	Bx  = Bx[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
	 	By  = By[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
	 	Bz  = Bz[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
 	#- 

	
	#- B-field bunch
# 	f.write('<DataArray type="Float32" Name="B_bunch" NumberOfComponents="3" format="binary"> \n')
# 	mesh = matrix2vectorField(Bxb,Byb,Bzb)
# 	s = base64.b64encode(np.array(mesh,dtype=np.float32))
#  	f.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
#  	f.write(s)
# 	f.write('</DataArray> \n')
	#- B-field total
	f.write('<DataArray type="Float32" Name="B" NumberOfComponents="3" format="binary"> \n')
	mesh = matrix2vectorField(Bx,By,Bz)
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	f.write(s)
	f.write('</DataArray> \n')
		#-Deallocate memory
	Bx=[0.];By=[0.];Bz=[0.];#Bxb=[0.];Byb=[0.];Bzb=[0.];



	#- Point Data End-#
 	f.write('</PointData> \n')
#	f.write('</CellData> \n')
	#- -#







		
# 	frame = 0
# 	s='%2.2i'%frame 				#conversion to 2-character-long-string
# 	file_name = 'Bdenout'+s+'.bin'
# 	matrix,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
# # 	file_name = 'Edenout'+s+'.bin'
# # 	matrix2,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
# 	#- cut & sign
# 	matrix = np.abs( matrix )
# # 	matrix2 = np.abs( matrix2 )
# 	size = matrix.shape
	

	

	#- writing vts header
# 	f = open('test.vts','w+')
# 	f.write('<?xml version="1.0"?>' + '\n')
# 	f.write('<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">' + '\n')
# 	f.write('<StructuredGrid WholeExtent=" %d %d %d %d %d %d "> \n' % (0,size[1],0,size[0],0,size[2]) )
# 	f.write('<Piece Extent=" %d %d %d %d %d %d "> \n' % (0,size[1],0,size[0],0,size[2]) )



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

# 	f.write('</CellData> \n')

	f.write('</Piece> \n')
	f.write('</StructuredGrid> \n')
	f.write('</VTKFile>')






#--- *** ---# DENSITY VTS: bunch+capillary+total
def write_density_vts(path,frame,X,Y,Z,cell_cut):

 	sf='%2.2i'%frame 				#conversion to 2-character-long-string
 	file_name 		= 'Bdenout'+sf+'.bin'
 	rhobunch, x,y,z	= read_ALaDyn_bin(path,file_name,'grid')
 	rhobunch		= np.abs( rhobunch )
 	file_name 		= 'Edenout'+sf+'.bin'
 	rhobck, x,y,z	= read_ALaDyn_bin(path,file_name,'grid')
 	rhobck			= np.abs( rhobck )
 	
 	#- matrix shaving
 	if cell_cut > 0:
	 	rhobunch = rhobunch[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
 		rhobck   = rhobck[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
	 	Y		 = Y[cell_cut:-cell_cut]
 		Z        = Z[cell_cut:-cell_cut]
 	#-
 	size 			= rhobunch.shape
 	#- 

 	
	#- writing vts header
		#total
	f_total = open(os.path.join(path,'VTS_files','ALaDyn_rho_'+sf+'.vts'),'w+')
	f_total.write('<?xml version="1.0"?>' + '\n')
	f_total.write('<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">' + '\n')
	f_total.write('<StructuredGrid WholeExtent=" %d %d %d %d %d %d "> \n' % (0,size[1]-1,0,size[2]-1,0,size[0]-1) )
	f_total.write('<Piece Extent=" %d %d %d %d %d %d "> \n' % (0,size[1]-1,0,size[2]-1,0,size[0]-1) )
		#bunch(es)
	f_bunch = open(os.path.join(path,'VTS_files','ALaDyn_rho_bunch_'+sf+'.vts'),'w+')
	f_bunch.write('<?xml version="1.0"?>' + '\n')
	f_bunch.write('<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">' + '\n')
	f_bunch.write('<StructuredGrid WholeExtent=" %d %d %d %d %d %d "> \n' % (0,size[1]-1,0,size[2]-1,0,size[0]-1) )
	f_bunch.write('<Piece Extent=" %d %d %d %d %d %d "> \n' % (0,size[1]-1,0,size[2]-1,0,size[0]-1) )
		#background
	f_bck = open(os.path.join(path,'VTS_files','ALaDyn_rho_bck_'+sf+'.vts'),'w+')
	f_bck.write('<?xml version="1.0"?>' + '\n')
	f_bck.write('<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">' + '\n')
	f_bck.write('<StructuredGrid WholeExtent=" %d %d %d %d %d %d "> \n' % (0,size[1]-1,0,size[2]-1,0,size[0]-1) )
	f_bck.write('<Piece Extent=" %d %d %d %d %d %d "> \n' % (0,size[1]-1,0,size[2]-1,0,size[0]-1) )
	#--- ---#


	#- generating vector-MESH
 	mesh=[]
 	for i in range(0,size[0]):
 		for j in range(0,size[1]):
 			for k in range(0,size[2]):
				mesh.append( Y[j] )
				mesh.append( Z[k] )
				mesh.append( X[i] )
	
	#- Writing MESH -#
		#it's the same for the whole 3 cases
		#total
	f_total.write('<Points> \n')
 	f_total.write('<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="binary"> \n')
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f_total.write(  base64.b64encode(np.array(len(s),dtype=np.int32))  )
 	f_total.write(  s  )
	f_total.write('</DataArray> \n')
	f_total.write('</Points> \n')
		#bunch(es)
	f_bunch.write('<Points> \n')
 	f_bunch.write('<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="binary"> \n')
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f_bunch.write(  base64.b64encode(np.array(len(s),dtype=np.int32))  )
 	f_bunch.write(  s  )
	f_bunch.write('</DataArray> \n')
	f_bunch.write('</Points> \n')
		#background
	f_bck.write('<Points> \n')
 	f_bck.write('<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="binary"> \n')
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f_bck.write(  base64.b64encode(np.array(len(s),dtype=np.int32))  )
 	f_bck.write(  s  )
	f_bck.write('</DataArray> \n')
	f_bck.write('</Points> \n')
	#- -#

	#- Point Data Begin-#
 	f_total.write('<PointData>\n')
 	f_bunch.write('<PointData>\n')
 	f_bck.write('<PointData>\n')
	#- -#

	#- Writing DensitIES -#
		#total
	f_total.write('<DataArray type="Float32" Name="rho" format="binary"> \n')
	mesh = matrix2vector( rhobunch+rhobck )
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f_total.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	f_total.write(s)
	f_total.write('</DataArray> \n')
		#bunch(es)
	f_bunch.write('<DataArray type="Float32" Name="rho_bunch" format="binary"> \n')
	mesh = matrix2vector( rhobunch )
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f_bunch.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	f_bunch.write(s)
	f_bunch.write('</DataArray> \n')
		#background
	f_bck.write('<DataArray type="Float32" Name="rho_bck" format="binary"> \n')
	mesh = matrix2vector( rhobck )
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f_bck.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	f_bck.write(s)
	f_bck.write('</DataArray> \n')

		#-Deallocate memory
	rhobunch = [0.]; rhobck = [0.]


	#- Point Data End-#
 	f_total.write('</PointData> \n')
 	f_bunch.write('</PointData> \n')
 	f_bck.write('</PointData> \n')
	#- -#

	# --- #
	f_total.write('</Piece> \n')
	f_bunch.write('</Piece> \n')
	f_bck.write('</Piece> \n')
	f_total.write('</StructuredGrid> \n')
	f_bunch.write('</StructuredGrid> \n')
	f_bck.write('</StructuredGrid> \n')
	f_total.write('</VTKFile>')
	f_bunch.write('</VTKFile>')
	f_bck.write('</VTKFile>')











#--- *** ---#
def write_E_vts(path,frame,X,Y,Z,cell_cut):
# 	for i in range(frame_begin, min(frame_end,last_output(os.getcwd())) + 1 ):
# 		write_VTS(path,i)

 	sf='%2.2i'%frame 				#conversion to 2-character-long-string
 	file_name = 'Exfout'+sf+'.bin'
	Ex,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Eyfout'+sf+'.bin'
	Ey,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Ezfout'+sf+'.bin'
	Ez,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
			#- bunch(es) -#
	file_name = 'Exbout'+sf+'.bin'
	Exb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Eybout'+sf+'.bin'
	Eyb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Ezbout'+sf+'.bin'
	Ezb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')

 	#- matrix shaving
 	if cell_cut > 0:
	 	Ex  = Ex[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
	 	Ey  = Ey[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
	 	Ez  = Ez[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
	 	Exb = Exb[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
	 	Eyb = Eyb[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
	 	Ezb = Ezb[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
	 	Y		 = Y[cell_cut:-cell_cut]
 		Z        = Z[cell_cut:-cell_cut]
 	#-
 	size = Ex.shape
 	#- 

 	
	#- writing vts header
	f = open(os.path.join(path,'VTS_files','ALaDyn_E_'+sf+'.vts'),'w+')
	f.write('<?xml version="1.0"?>' + '\n')
	f.write('<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">' + '\n')
	f.write('<StructuredGrid WholeExtent=" %d %d %d %d %d %d "> \n' % (0,size[1]-1,0,size[2]-1,0,size[0]-1) )
	f.write('<Piece Extent=" %d %d %d %d %d %d "> \n' % (0,size[1]-1,0,size[2]-1,0,size[0]-1) )
	#--- ---#

	#- writing Ex-vts header
	fEx = open(os.path.join(path,'VTS_files','ALaDyn_Ex_'+sf+'.vts'),'w+')
	fEx.write('<?xml version="1.0"?>' + '\n')
	fEx.write('<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">' + '\n')
	fEx.write('<StructuredGrid WholeExtent=" %d %d %d %d %d %d "> \n' % (0,size[1]-1,0,size[2]-1,0,size[0]-1) )
	fEx.write('<Piece Extent=" %d %d %d %d %d %d "> \n' % (0,size[1]-1,0,size[2]-1,0,size[0]-1) )
	#- writing Ey-vts header
	fEy = open(os.path.join(path,'VTS_files','ALaDyn_Ey_'+sf+'.vts'),'w+')
	fEy.write('<?xml version="1.0"?>' + '\n')
	fEy.write('<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">' + '\n')
	fEy.write('<StructuredGrid WholeExtent=" %d %d %d %d %d %d "> \n' % (0,size[1]-1,0,size[2]-1,0,size[0]-1) )
	fEy.write('<Piece Extent=" %d %d %d %d %d %d "> \n' % (0,size[1]-1,0,size[2]-1,0,size[0]-1) )
	#- writing Ex-vts header
	fEz = open(os.path.join(path,'VTS_files','ALaDyn_Ez_'+sf+'.vts'),'w+')
	fEz.write('<?xml version="1.0"?>' + '\n')
	fEz.write('<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">' + '\n')
	fEz.write('<StructuredGrid WholeExtent=" %d %d %d %d %d %d "> \n' % (0,size[1]-1,0,size[2]-1,0,size[0]-1) )
	fEz.write('<Piece Extent=" %d %d %d %d %d %d "> \n' % (0,size[1]-1,0,size[2]-1,0,size[0]-1) )
	#--- ---#




	#- generating vector-MESH
 	mesh=[]
 	for i in range(0,size[0]):
 		for j in range(0,size[1]):
 			for k in range(0,size[2]):
				mesh.append( Y[j] )
				mesh.append( Z[k] )
				mesh.append( X[i] )
	#- Writing MESH -#
	f.write('<Points> \n')
 	f.write('<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="binary"> \n')
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f.write(  base64.b64encode(np.array(len(s),dtype=np.int32))  )
 	f.write(  s  )
	f.write('</DataArray> \n')
	f.write('</Points> \n')
	#- -#

	#- Writing ExMESH -#
	fEx.write('<Points> \n')
 	fEx.write('<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="binary"> \n')
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	fEx.write(  base64.b64encode(np.array(len(s),dtype=np.int32))  )
 	fEx.write(  s  )
	fEx.write('</DataArray> \n')
	fEx.write('</Points> \n')
	#- Writing ExMESH -#
	fEy.write('<Points> \n')
 	fEy.write('<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="binary"> \n')
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	fEy.write(  base64.b64encode(np.array(len(s),dtype=np.int32))  )
 	fEy.write(  s  )
	fEy.write('</DataArray> \n')
	fEy.write('</Points> \n')
	#- Writing ExMESH -#
	fEz.write('<Points> \n')
 	fEz.write('<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="binary"> \n')
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	fEz.write(  base64.b64encode(np.array(len(s),dtype=np.int32))  )
 	fEz.write(  s  )
	fEz.write('</DataArray> \n')
	fEz.write('</Points> \n')
	#- -#


	#- Point Data Begin-#
 	f.write('<PointData>\n')
 	fEx.write('<PointData>\n')
 	fEy.write('<PointData>\n')
 	fEz.write('<PointData>\n')
	#- -#

	#- Writing E-field -#
	#- E-field total
	f.write('<DataArray type="Float32" Name="E" NumberOfComponents="3" format="binary"> \n')
	mesh = matrix2vectorField(Ex+Exb,Ey+Eyb,Ez+Ezb)
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	f.write(s)
	f.write('</DataArray> \n')

	#- Writing Ex-field -#
	fEx.write('<DataArray type="Float32" Name="Ex" NumberOfComponents="1" format="binary"> \n')
	mesh = matrix2vector(Ex+Exb)
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	fEx.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	fEx.write(s)
	fEx.write('</DataArray> \n')
	#- Writing Ex-field -#
	fEy.write('<DataArray type="Float32" Name="Ey" NumberOfComponents="1" format="binary"> \n')
	mesh = matrix2vector(Ey+Eyb)
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	fEy.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	fEy.write(s)
	fEy.write('</DataArray> \n')
	#- Writing Ex-field -#
	fEz.write('<DataArray type="Float32" Name="Ez" NumberOfComponents="1" format="binary"> \n')
	mesh = matrix2vector(Ez+Ezb)
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	fEz.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	fEz.write(s)
	fEz.write('</DataArray> \n')

		#-Deallocate memory
	Ex=[0.];Ey=[0.];Ez=[0.];Exb=[0.];Eyb=[0.];Ezb=[0.];


	#- Point Data End-#
 	f.write('</PointData> \n')
 	fEx.write('</PointData> \n')
 	fEy.write('</PointData> \n')
 	fEz.write('</PointData> \n')
	#- -#
	f.write('</Piece> \n')
	f.write('</StructuredGrid> \n')
	f.write('</VTKFile>')

	fEx.write('</Piece> \n')
	fEx.write('</StructuredGrid> \n')
	fEx.write('</VTKFile>')

	fEy.write('</Piece> \n')
	fEy.write('</StructuredGrid> \n')
	fEy.write('</VTKFile>')

	fEz.write('</Piece> \n')
	fEz.write('</StructuredGrid> \n')
	fEz.write('</VTKFile>')




#--- *** ---#
def write_B_vts(path,frame,X,Y,Z,cell_cut):
# 	for i in range(frame_begin, min(frame_end,last_output(os.getcwd())) + 1 ):
# 		write_VTS(path,i)

 	sf='%2.2i'%frame 				#conversion to 2-character-long-string
 	file_name = 'Bxfout'+sf+'.bin'
	Bx,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Byfout'+sf+'.bin'
	By,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Bzfout'+sf+'.bin'
	Bz,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
 	size 			= Bx.shape


 	#- matrix shaving
 	if cell_cut > 0:
	 	Bx  = Bx[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
	 	By  = By[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
	 	Bz  = Bz[:,cell_cut:-cell_cut,cell_cut:-cell_cut]
	 	Y		 = Y[cell_cut:-cell_cut]
 		Z        = Z[cell_cut:-cell_cut]
 	#-
 	size = Bx.shape
 	#- 

 	
	#- writing vts header
	f = open(os.path.join(path,'VTS_files','ALaDyn_B_'+sf+'.vts'),'w+')
	f.write('<?xml version="1.0"?>' + '\n')
	f.write('<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">' + '\n')
	f.write('<StructuredGrid WholeExtent=" %d %d %d %d %d %d "> \n' % (0,size[1]-1,0,size[2]-1,0,size[0]-1) )
	f.write('<Piece Extent=" %d %d %d %d %d %d "> \n' % (0,size[1]-1,0,size[2]-1,0,size[0]-1) )
	#--- ---#


	#- generating vector-MESH
 	mesh=[]
 	for i in range(0,size[0]):
 		for j in range(0,size[1]):
 			for k in range(0,size[2]):
				mesh.append( Y[j] )
				mesh.append( Z[k] )
				mesh.append( X[i] )
	#- Writing MESH -#
	f.write('<Points> \n')
 	f.write('<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="binary"> \n')
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f.write(  base64.b64encode(np.array(len(s),dtype=np.int32))  )
 	f.write(  s  )
	f.write('</DataArray> \n')
	f.write('</Points> \n')
	#- -#

	#- Point Data Begin-#
 	f.write('<PointData>\n')
	#- -#

	#- Writing B-field -#
	f.write('<DataArray type="Float32" Name="B" NumberOfComponents="3" format="binary"> \n')
	mesh = matrix2vectorField(Bx,By,Bz)
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	f.write(s)
	f.write('</DataArray> \n')
		#-Deallocate memory
	Bx=[0.];By=[0.];Bz=[0.];#Bxb=[0.];Byb=[0.];Bzb=[0.];

	#- Point Data End-#
 	f.write('</PointData> \n')
	#- -#
	f.write('</Piece> \n')
	f.write('</StructuredGrid> \n')
	f.write('</VTKFile>')










#--- *** ---#
def write_vts_section_longitudinal(path,frame,X,Y,Z,cell_cut):

 	sf='%2.2i'%frame 				#conversion to 2-character-long-string
 	file_name 		= 'Bdenout'+sf+'.bin'
 	rhobunch, x,y,z	= read_ALaDyn_bin(path,file_name,'grid')
 	rhobunch		= np.abs( rhobunch )
 	file_name 		= 'Edenout'+sf+'.bin'
 	rhobck, x,y,z	= read_ALaDyn_bin(path,file_name,'grid')
 	rhobck			= np.abs( rhobck )

	#-original size
	sizeO = rhobunch.shape
 	
 	#- matrix shaving
 	if cell_cut > 0:
	 	rhobunch = rhobunch[:, sizeO[1]/2, cell_cut:-cell_cut]
 		rhobck   = rhobck[:, sizeO[1]/2, cell_cut:-cell_cut]
	 	Y		 = Y[ sizeO[1]/2 ]
 		Z        = Z[cell_cut:-cell_cut]
	else:
	 	rhobunch = rhobunch[:, sizeO[1]/2,:]
 		rhobck   = rhobck[:, sizeO[1]/2,:]
	 	Y		 = Y[ sizeO[1]/2 ]
 		Z        = Z[:]
 	#-
 	size 		 = rhobunch.shape
 	#- 

 	
	#- writing vts header
	f = open(os.path.join(path,'VTS_files','ALaDyn_section_longitudinal_output_'+sf+'.vts'),'w+')
	f.write('<?xml version="1.0"?>' + '\n')
	f.write('<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">' + '\n')
	f.write('<StructuredGrid WholeExtent=" %d %d %d %d %d %d "> \n' % (0,0,0,size[1]-1,0,size[0]-1) )
	f.write('<Piece Extent=" %d %d %d %d %d %d "> \n' % (0,0,0,size[1]-1,0,size[0]-1) )
	#--- ---#


	#- generating vector-MESH
 	mesh=[]
 	for i in range(0,size[0]):
		for k in range(0,size[1]):
			mesh.append( Y    )
			mesh.append( Z[k] )
			mesh.append( X[i] )
	#- Writing MESH -#
	f.write('<Points> \n')
 	f.write('<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="binary"> \n')
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f.write(  base64.b64encode(np.array(len(s),dtype=np.int32))  )
 	f.write(  s  )
	f.write('</DataArray> \n')
	f.write('</Points> \n')
	#- -#

	#- Point Data Begin-#
 	f.write('<PointData>\n')
#	f.write('<CellData> \n')
	#- -#

	#- Writing BUNCH Density -#
	f.write('<DataArray type="Float32" Name="rho_bunch" format="binary"> \n')
	mesh = matrix2vector( rhobunch )
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	f.write(s)
	f.write('</DataArray> \n')

	#- Writing Tot-Density -#
	f.write('<DataArray type="Float32" Name="rho" format="binary"> \n')
	mesh = matrix2vector( rhobunch+rhobck )
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	f.write(s)
	f.write('</DataArray> \n')
		#-Deallocate memory
	rhobunch = [0.]; rhobck = [0.]

	#- Writing E-field -#
			#- background -#
 	file_name = 'Exfout'+sf+'.bin'
	Ex,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Eyfout'+sf+'.bin'
	Ey,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Ezfout'+sf+'.bin'
	Ez,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
			#- bunch(es) -#
	file_name = 'Exbout'+sf+'.bin'
	Exb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Eybout'+sf+'.bin'
	Eyb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Ezbout'+sf+'.bin'
	Ezb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')

 	#- matrix shaving
 	if cell_cut > 0:
	 	Ex  = Ex[:,sizeO[1]/2,cell_cut:-cell_cut]
	 	Ey  = Ey[:,sizeO[1]/2,cell_cut:-cell_cut]
	 	Ez  = Ez[:,sizeO[1]/2,cell_cut:-cell_cut]
	 	Exb = Exb[:,sizeO[1]/2,cell_cut:-cell_cut]
	 	Eyb = Eyb[:,sizeO[1]/2,cell_cut:-cell_cut]
	 	Ezb = Ezb[:,sizeO[1]/2,cell_cut:-cell_cut]
 	else:
	 	Ex  = Ex[:,sizeO[1]/2,:]
	 	Ey  = Ey[:,sizeO[1]/2,:]
	 	Ez  = Ez[:,sizeO[1]/2,:]
	 	Exb = Exb[:,sizeO[1]/2,:]
	 	Eyb = Eyb[:,sizeO[1]/2,:]
	 	Ezb = Ezb[:,sizeO[1]/2,:]
 	#- 

	
	#- E-field bunch
	f.write('<DataArray type="Float32" Name="E_bunch" NumberOfComponents="3" format="binary"> \n')
	mesh = matrix2vectorField(Exb,Eyb,Ezb)
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	f.write(s)
	f.write('</DataArray> \n')
	#- E-field total
	f.write('<DataArray type="Float32" Name="E" NumberOfComponents="3" format="binary"> \n')
	mesh = matrix2vectorField(Ex+Exb,Ey+Eyb,Ez+Ezb)
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	f.write(s)
	f.write('</DataArray> \n')
	#- Ex-total
	f.write('<DataArray type="Float32" Name="Ex" NumberOfComponents="1" format="binary"> \n')
	mesh = matrix2vector(Ex+Exb)
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	f.write(s)
	f.write('</DataArray> \n')
	#- Ey-total
	f.write('<DataArray type="Float32" Name="Ey" NumberOfComponents="1" format="binary"> \n')
	mesh = matrix2vector(Ey+Eyb)
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	f.write(s)
	f.write('</DataArray> \n')
	#- Ez-total
	f.write('<DataArray type="Float32" Name="Ez" NumberOfComponents="1" format="binary"> \n')
	mesh = matrix2vector(Ez+Ezb)
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	f.write(s)
	f.write('</DataArray> \n')
		#-Deallocate memory
	Ex=[0.];Ey=[0.];Ez=[0.];Exb=[0.];Eyb=[0.];Ezb=[0.];


	#- Writing B-field -#
			#- background -#
 	file_name = 'Bxfout'+sf+'.bin'
	Bx,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Byfout'+sf+'.bin'
	By,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Bzfout'+sf+'.bin'
	Bz,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
# 			#- bunch(es) -#
# 	file_name = 'Bxbout'+sf+'.bin'
# 	Bxb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
# 	file_name = 'Bybout'+sf+'.bin'
# 	Byb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
# 	file_name = 'Bzbout'+sf+'.bin'
# 	Bzb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')

 	#- matrix shaving
 	if cell_cut > 0:
	 	Bx  = Bx[:,sizeO[1]/2,cell_cut:-cell_cut]
	 	By  = By[:,sizeO[1]/2,cell_cut:-cell_cut]
	 	Bz  = Bz[:,sizeO[1]/2,cell_cut:-cell_cut]
 	else:
	 	Bx  = Bx[:,sizeO[1]/2,:]
	 	By  = By[:,sizeO[1]/2,:]
	 	Bz  = Bz[:,sizeO[1]/2,:]
 	#- 

	
	#- B-field bunch
# 	f.write('<DataArray type="Float32" Name="B_bunch" NumberOfComponents="3" format="binary"> \n')
# 	mesh = matrix2vectorField(Bxb,Byb,Bzb)
# 	s = base64.b64encode(np.array(mesh,dtype=np.float32))
#  	f.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
#  	f.write(s)
# 	f.write('</DataArray> \n')
	#- B-field total
	f.write('<DataArray type="Float32" Name="B" NumberOfComponents="3" format="binary"> \n')
	mesh = matrix2vectorField(Bx,By,Bz)
	s = base64.b64encode(np.array(mesh,dtype=np.float32))
 	f.write(base64.b64encode(np.array(len(s),dtype=np.int32)))
 	f.write(s)
	f.write('</DataArray> \n')
		#-Deallocate memory
	Bx=[0.];By=[0.];Bz=[0.];#Bxb=[0.];Byb=[0.];Bzb=[0.];



	#- Point Data End-#
 	f.write('</PointData> \n')
#	f.write('</CellData> \n')
	#- -#


	f.write('</Piece> \n')
	f.write('</StructuredGrid> \n')
	f.write('</VTKFile>')
