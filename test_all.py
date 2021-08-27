from simulation import functions as fc
import numpy as np

def create_mesh_test():
    
    """
    check the create_mesh function by checking the length of the name list and the mesh created.
    """
    print('testing function create_mesh')
    m_list = fc.create_mesh(0.025, 0.1)
    assert len(m_list)==55
    assert m_list[1]=='2 hr=0.025_dfc0.05_000 degree.msh'
    

def stress_determination_test():
    """
    check the stress determination function by comparing the example with the reference number obtained in the history data.
    check if the failure parameter is greater than one when the mesh breaks.
    check if the failure parameter residual is smaller than the set value.
    """
    print('testing stress determination')
    mesh = '4 hr=0.025_dfc0.05_170 degree.msh'
    stress, failure = fc.find_py(mesh)
    assert np.isclose(stress,44641113.28125)
    assert failure > 1
    assert (failure -1)<1e-4
    
create_mesh_test()
stress_determination_test()
print('Test Complete!')


