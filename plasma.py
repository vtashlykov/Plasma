import ctypes
from typing import List, Union, Optional

def list_to_c_double_pointer(arr: Optional[List[Union[int, float]]] = None):
    if arr is not None:
        return (ctypes.c_double * len(arr))(*arr)
    return (ctypes.c_double * len([]))(*[])

class PlasmaStruct(ctypes.Structure):
    _fields_ = [('lambda', ctypes.c_double),
                ('Ne', ctypes.c_double),
                ('Te', ctypes.c_double),
                ('Ti', ctypes.c_double),
                ('nu_i', ctypes.c_double),
                ('nu_e', ctypes.c_double),
                ('Con', ctypes.c_double*103)
    ]


# Prepare variables.
FLENGTH=5000
con = [0] * 103
con[7] = 100
c_con = list_to_c_double_pointer(con)
arr1=[0]*2*FLENGTH
arr2=[0]*2*FLENGTH
c_x_list = (ctypes.c_double * len(arr1))(*arr1)
c_y_list = (ctypes.c_double * len(arr2))(*arr2)

# struct = PlasmaStruct(
#     ctypes.c_double(2),
#     ctypes.c_double(1),
#     ctypes.c_double(2000),
#     ctypes.c_double(1000),
#     ctypes.c_double(1000),
#     ctypes.c_double(100),
#     c_con
# )

lib = ctypes.CDLL('./build_lib/libPlasma')
lib.Set_pars.argtypes = [ctypes.c_char_p, ctypes.POINTER(PlasmaStruct)]
lib.Set_pars.restype = None
lib.Spectrum.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(PlasmaStruct)]
lib.Spectrum.restype = None

filename=b'Pars.config'
struct = PlasmaStruct()
lib.Set_pars(filename, struct)
lib.Spectrum(c_x_list, c_y_list, ctypes.byref(struct))
file=open("Spectrum.dat", 'w')
for i in range(2*FLENGTH):
    file.write(str(c_x_list[i])+'\t'+str(c_y_list[i])+'\n')
print([{'x': k, 'y': v} for k, v in zip(list(c_x_list), list(c_y_list))] )