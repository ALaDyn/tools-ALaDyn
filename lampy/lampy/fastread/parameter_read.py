import os
import struct


def export_parameters(dir_path, file_name):
    path = os.path.join(os.path.join(dir_path, file_name))
    file_read = open(path+'.bin', 'rb')

    file_read.seek(4)
    Nparam = struct.unpack('=i', file_read.read(4))[0]
    file_read.seek(16)
    integerdata_temp = struct.unpack('='+Nparam*'i', file_read.read(Nparam*4))
    file_read.seek(104)
    realdata_temp = struct.unpack('='+Nparam*'f', file_read.read(Nparam*4))

    file_read.close()
    integerdata = dict()
    integerdata['proc_y'] = integerdata_temp[0]
    integerdata['proc_z'] = integerdata_temp[1]
    integerdata['proc_x'] = integerdata_temp[2]
    integerdata['nx_tot'] = integerdata_temp[3]
    integerdata['ny_tot'] = integerdata_temp[4]
    integerdata['ny_per_proc'] = integerdata_temp[5]
    integerdata['nz_tot'] = integerdata_temp[6]
    integerdata['nz_per_proc'] = integerdata_temp[7]
    integerdata['boundary_x'] = integerdata_temp[8]
    integerdata['boundary_y'] = integerdata_temp[9]
    integerdata['boundary_z'] = integerdata_temp[10]
    integerdata['env_non_env'] = integerdata_temp[11]
    integerdata['unknown1'] = integerdata_temp[12]
    integerdata['unknown2'] = integerdata_temp[13]
    integerdata['ndim'] = integerdata_temp[14]
    integerdata['unknown3'] = integerdata_temp[15]
    integerdata['unknown4'] = integerdata_temp[16]
    integerdata['unknown5'] = integerdata_temp[17]
    integerdata['unknown6'] = integerdata_temp[18]
    integerdata['unknown7'] = integerdata_temp[19]

    realdata = dict()
    realdata['time'] = realdata_temp[0]
    realdata['x_min'] = realdata_temp[1]
    realdata['x_max'] = realdata_temp[2]
    realdata['y_min'] = realdata_temp[3]
    realdata['y_max'] = realdata_temp[4]
    realdata['z_min'] = realdata_temp[5]
    realdata['z_max'] = realdata_temp[6]
    realdata['pulse_duration'] = realdata_temp[7]
    realdata['waist'] = realdata_temp[8]
    realdata['n_over_nc'] = realdata_temp[9]
    realdata['a0'] = realdata_temp[10]
    realdata['lambda_0'] = realdata_temp[11]
    realdata['E0'] = realdata_temp[12]
    realdata['unknown1'] = realdata_temp[13]
    realdata['np_per_cell'] = realdata_temp[14]
    realdata['unknown2'] = realdata_temp[15]
    realdata['unknown3'] = realdata_temp[16]
    realdata['unknown4'] = realdata_temp[17]
    realdata['unknown5'] = realdata_temp[18]
    realdata['unknown6'] = realdata_temp[19]

    return (integerdata, realdata)


def _read_file_timestep(file_path):

    file_read = open(file_path, 'rb')

    file_read.seek(4)
    Nparam = struct.unpack('=i', file_read.read(4))[0]
    file_read.seek(16)
    _ = struct.unpack('='+Nparam*'i', file_read.read(Nparam*4))
    file_read.seek(104)
    realdata_temp = struct.unpack('='+Nparam*'f', file_read.read(Nparam*4))
    file_read.close()
    realdata = dict()
    realdata['time'] = realdata_temp[0]

    return realdata['time']


def _output_directories(dir_path):

    listdir = [o for o in os.listdir(dir_path)
               if os.path.isdir(os.path.join(dir_path, o)) and len(o) == 4]

    def intval(val):
        return int(val)

    listdir.sort(key=intval)
    return listdir


def _read_box_limits(file_path):

    file_read = open(file_path, 'rb')

    file_read.seek(4)
    Nparam = struct.unpack('=i', file_read.read(4))[0]
    file_read.seek(16)
    _ = struct.unpack('='+Nparam*'i', file_read.read(Nparam*4))
    file_read.seek(104)
    realdata_temp = struct.unpack('='+Nparam*'f', file_read.read(Nparam*4))
    file_read.close()
    realdata = dict()
    realdata['x_min'] = realdata_temp[1]
    realdata['x_max'] = realdata_temp[2]
    realdata['y_min'] = realdata_temp[3]
    realdata['y_max'] = realdata_temp[4]
    realdata['z_min'] = realdata_temp[5]
    realdata['z_max'] = realdata_temp[6]

    return realdata
