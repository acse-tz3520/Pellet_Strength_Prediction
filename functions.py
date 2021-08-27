import gmsh
import matplotlib.pyplot as plt
import pandas as pd
import pickle
from firedrake import *
import sklearn


def create_mesh(hole_r, pellet_radius):
    """
    Create the mesh set with the give hole radius and pellet radius,
    output the meshes and the name list of the meshes created.

    Parameters:
    -----------
    hole_r : the hole radius(unit:m)
    pellet_radius: the pellet radius(unit:m)
    """

    # Create a range of meshes with different orientations
    angle_range = [10*i for i in range(18)]
    # Create name list for the meshes generated
    mesh_list = []

    # make sure the mesh size small enough to simulate the infinitesimal elements
    mesh_size = hole_r/5
    # radius for the pellet raduys
    outer_radius = pellet_radius
    # assign hole radius
    hole_radius = hole_r
    # initial value for the spacing between pellet centre and hole centre to be twice the hole radius
    d_from_c = hole_r*2

    # Mesh for single hole at the centre
    # Mesh for 1 hole at the centre
    if hole_radius < outer_radius:

        gmsh.initialize()
        # add a new model
        gmsh.model.add("t1")

        # assume contact angle = 1 degree,the angle for the sector will be (pi/180)*1 in radian
        # horizontal projection for the radius of the pellet
        R_h = outer_radius*sin(pi/(180*2))
        # vertical projection for the radius of the pellet
        R_v = outer_radius*cos(pi/(180*2))

        # bottom right vertex
        gmsh.model.geo.addPoint(R_h, -R_v, 0, mesh_size, 1)
        # bottom left vertex
        gmsh.model.geo.addPoint(-R_h, -R_v, 0, mesh_size, 2)
        # top left vertex
        gmsh.model.geo.addPoint(-R_h, R_v, 0, mesh_size, 3)
        # top right vertex
        gmsh.model.geo.addPoint(R_h, R_v, 0, mesh_size, 4)
        # centre point
        gmsh.model.geo.addPoint(0, 0, 0, 0.5, 5)

        # bottom arc of the pellet
        gmsh.model.geo.addCircleArc(1, 5, 2, 1)
        # left arc of the pellet
        gmsh.model.geo.addCircleArc(2, 5, 3, 2)
        # top arc of the pellet
        gmsh.model.geo.addCircleArc(3, 5, 4, 3)
        # right arc of the pellet
        gmsh.model.geo.addCircleArc(4, 5, 1, 4)
        gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)

        # centre
        gmsh.model.geo.addPoint(
            0, 0, 0, mesh_size, 6)

        # four vertices in the circle
        # top right vertex
        gmsh.model.geo.addPoint(+hole_radius/sqrt(2),
                                +hole_radius/sqrt(2), 0, mesh_size, 7)
        # bottom left vertex
        gmsh.model.geo.addPoint(-hole_radius/sqrt(2),
                                -hole_radius/sqrt(2), 0, mesh_size, 8)
        # top left vertex
        gmsh.model.geo.addPoint(-hole_radius/sqrt(2),
                                +hole_radius/sqrt(2),  0, mesh_size, 9)
        # bottom right vertex
        gmsh.model.geo.addPoint(+hole_radius/sqrt(2),
                                -hole_radius/sqrt(2),  0, mesh_size, 10)

        # top arc of the hole
        gmsh.model.geo.addCircleArc(9, 6, 7, 6)
        # right arc of the hole
        gmsh.model.geo.addCircleArc(7, 6, 10, 7)
        # bottom arc of the hole
        gmsh.model.geo.addCircleArc(10, 6, 8, 8)
        # left arc of the hole
        gmsh.model.geo.addCircleArc(8, 6, 9, 9)
        # complete the hole
        gmsh.model.geo.addCurveLoop([9, 6, 7, 8], 2)

        # construct the plane surface
        gmsh.model.geo.addPlaneSurface([1, 2], 1)
        gmsh.model.geo.synchronize()

        # Labelling the top and bottom arcs of the pellet
        # bottom labelled 5
        gmsh.model.addPhysicalGroup(1, [1], 5)
        # top labelled 5
        gmsh.model.addPhysicalGroup(1, [3], 6)
        # adding physical group
        ps = gmsh.model.addPhysicalGroup(2, [1])
        gmsh.model.setPhysicalName(2, ps, "PunchedDom")

        # generate a 2D mesh
        gmsh.model.mesh.generate(2)
        # name and save
        # make the angle three digits

        an = '000'
        # convert the hole radius to string
        hr = str(float(hole_radius))
        # convert the distance from centre to stirng
        d_f_c = '0'
        # assign the name for the msh file saved
        name = '01 hr='+hr+'_dfc'+d_f_c+'_'+an+' degree.msh'
        # adding the generated mesh to the list
        mesh_list.append(name)
        # save
        gmsh.write(name)

        gmsh.finalize()

    # Starting meshes generation with number of holes >=2
    number_of_holes = 2

    while True:
        # checking if the holes are completely inside the pellet
        while 1.5*hole_radius+d_from_c <= outer_radius:
            # calculating the spacing between two adjacent hole' centres
            residual_angle = (180-(360/number_of_holes))/2
            cc_distance = d_from_c*cos(residual_angle*pi/180)*2
            if 2.5*hole_radius <= cc_distance:

                # Creating mesh set for different orientation
                for angle_step in angle_range:
                    gmsh.initialize()
                    # add a new model
                    gmsh.model.add("t"+str(angle_step))

                    # assume contact angle = 1 degree,the angle for the sector will be (pi/180)*1 in radian
                    # horizontal projection for the radius of the pellet
                    R_h = outer_radius*sin(pi/(180*2))
                    # vertical projection for the radius of the pellet
                    R_v = outer_radius*cos(pi/(180*2))

                    # bottom right vertex
                    gmsh.model.geo.addPoint(R_h, -R_v, 0, mesh_size, 1)
                    # bottom left vertex
                    gmsh.model.geo.addPoint(-R_h, -R_v, 0, mesh_size, 2)
                    # top left vertex
                    gmsh.model.geo.addPoint(-R_h, R_v, 0, mesh_size, 3)
                    # top right vertex
                    gmsh.model.geo.addPoint(R_h, R_v, 0, mesh_size, 4)
                    # centre point
                    gmsh.model.geo.addPoint(0, 0, 0, 0.5, 5)

                    # bottom arc of the pellet
                    gmsh.model.geo.addCircleArc(1, 5, 2, 1)
                    # left arc of the pellet
                    gmsh.model.geo.addCircleArc(2, 5, 3, 2)
                    # top arc of the pellet
                    gmsh.model.geo.addCircleArc(3, 5, 4, 3)
                    # right arc of the pellet
                    gmsh.model.geo.addCircleArc(4, 5, 1, 4)
                    gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)

                    # 5 points per circle(4 verticies + one centre),therefore consum 5 point labels for
                    # each holes
                    j = 5
                    # angle between two adjacent hole' centres in radian
                    angle = 2*pi/number_of_holes
                    # number of elements
                    n_element = [1]

                    for i in range(number_of_holes):

                        # vertical projection for the radius of the hole
                        r_h = d_from_c*cos(angle*i+angle_step*pi/180)
                        # vertical projection for the radius of the hole
                        r_v = d_from_c*sin(angle*i+angle_step*pi/180)

                        # centre
                        gmsh.model.geo.addPoint(
                            r_h, r_v, 0, mesh_size, 6+i*j)

                        # four vertices in the circle
                        # top right vertex
                        gmsh.model.geo.addPoint(r_h+hole_radius/sqrt(2),
                                                r_v+hole_radius/sqrt(2), 0, mesh_size, 7+i*j)
                        # bottom left vertex
                        gmsh.model.geo.addPoint(r_h-hole_radius/sqrt(2),
                                                r_v-hole_radius/sqrt(2), 0, mesh_size, 8+i*j)
                        # top left vertex
                        gmsh.model.geo.addPoint(r_h-hole_radius/sqrt(2),
                                                r_v+hole_radius/sqrt(2),  0, mesh_size, 9+i*j)
                        # bottom right vertex
                        gmsh.model.geo.addPoint(r_h+hole_radius/sqrt(2),
                                                r_v-hole_radius/sqrt(2),  0, mesh_size, 10+i*j)

                        # top arc of the hole
                        gmsh.model.geo.addCircleArc(9+i*j, 6+i*j, 7+i*j, 6+i*j)
                        # right arc of the hole
                        gmsh.model.geo.addCircleArc(
                            7+i*j, 6+i*j, 10+i*j, 7+i*j)
                        # bottom arc of the hole
                        gmsh.model.geo.addCircleArc(
                            10+i*j, 6+i*j, 8+i*j, 8+i*j)
                        # left arc of the hole
                        gmsh.model.geo.addCircleArc(8+i*j, 6+i*j, 9+i*j, 9+i*j)
                        # complete the hole
                        gmsh.model.geo.addCurveLoop(
                            [9+i*j, 6+i*j, 7+i*j, 8+i*j], 2+i)
                        n_element.append(2+i)

                    # construct the plane surface
                    gmsh.model.geo.addPlaneSurface(n_element, 1)
                    gmsh.model.geo.synchronize()

                    # Labelling the top and bottom arcs of the pellet
                    # bottom labelled 5
                    gmsh.model.addPhysicalGroup(1, [1], 5)
                    # top labelled 5
                    gmsh.model.addPhysicalGroup(1, [3], 6)
                    # adding physical group
                    ps = gmsh.model.addPhysicalGroup(2, [1])
                    gmsh.model.setPhysicalName(2, ps, "PunchedDom")

                    # generate a 2D mesh
                    gmsh.model.mesh.generate(2)
                    # name and save
                    # make the angle three digits
                    if angle_step >= 100:
                        an = str(angle_step)
                    elif angle_step >= 10:
                        an = '0'+str(angle_step)
                    else:
                        an = '000'
                    # convert the hole radius to string
                    hr = str(float(hole_radius))
                    # convert the distance from centre to stirng
                    d_f_c = str(d_from_c)
                    # assign the name for the msh file saved
                    name = str(number_of_holes)+' hr='+hr + \
                        '_dfc'+d_f_c+'_'+an+' degree.msh'
                    # adding the generated mesh to the list
                    mesh_list.append(name)
                    # save
                    gmsh.write(name)
                    gmsh.finalize()

            # move the hole away from the pellet centre by one unit of the hole radius
            d_from_c += hole_radius

        # adding one hole to the pellet
        number_of_holes += 1
        # reset the distance from centre to twice the hole radius
        d_from_c = hole_radius*2
        # the spacing between two adjacent holes' centres
        residual_angle = (180-(360/number_of_holes))/2
        while 2.5*hole_radius >= (d_from_c*cos(residual_angle*pi/180)*2):
            # the d_from_c is not sufficient to maintain all holes without merge, move the hole
            # away by one unit of the hole radius
            d_from_c += hole_radius

        # checking if the holes are still completely in the pellet
        if hole_radius+d_from_c >= outer_radius:
            break
    return mesh_list


def compression(mesh_file, py):
    """
    Simulate the compression test on the mesh imported, with the lower plate be stationary
    and the upper plate compressing downwards. 

    Parameters:
    -----------
    mesh_file : the mesh file for the pellet
    py: stress loaded on the top of the plate
    """

    # value set
    young_modulus = 40.5e9
    poission_ratio = 0.17
    rho = 2310
    # loading the mesh
    mesh_g = Mesh(mesh_file)
    # construct the vector function space
    V = VectorFunctionSpace(mesh_g, "CG", 1)
    # Dirichlet boundary condition at the bottom, both horizontal and vertical component 0
    bc = DirichletBC(V, Constant([0, 0]), 5)

    # Enter constant value
    ym = Constant(young_modulus)  # young modulus
    ps = Constant(poission_ratio)  # possion ratio
    rho = Constant(rho)  # density
    g = Constant(1)
    f = as_vector([0, -rho*g])  # the test function(body force per unit volume)

    mu = ym/(2*(Constant(1)+ps))  # shear modulus
    lambda_ = ym*ps/((Constant(1)+ps)*(Constant(1)-2*ps))  # strain
    Id = Identity(mesh_g.geometric_dimension())  # 2x2 Identity tensor

    # only setting py as the compression is vertically downwards
    px = 0
    py = py
    # Neumman boundary condition
    tra_x = as_vector([Constant(px), 0])
    tra_y = as_vector([0, Constant(-py)])

    # the symmetric strain rate tensor
    def epsilon(u):
        return 0.5*(grad(u) + grad(u).T)

    # the stress tensor
    def sigma(u):
        return lambda_*div(u)*Id + 2*mu*epsilon(u)
    # solving the equations
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(sigma(u), epsilon(v))*dx
    L = dot(f, v)*dx+dot(tra_x, v)*ds(6)+dot(tra_y, v)*ds(6)

    uh = Function(V)
    solve(a == L, uh,
          bcs=bc,
          solver_parameters={"ksp_type": "cg",
                             "ksp_max_it": 200,
                             "pc_type": "gamg",
                             "mat_type": "aij",
                             "ksp_monitor": None})

    # output the stress
    P0_tensor = TensorFunctionSpace(mesh_g, 'DG', 0)
    stress = Function(P0_tensor)
    stress.interpolate(sigma(uh))

    return stress


def extract_stress_components(mesh_file, stress, plotting=False):
    """
    Extact the maxmium principle stress from the stress obtained from the compression function.

    Parameters:
    -----------
    mesh_file : the mesh file for the pellet
    stress: stress components obtained from function "compression"
    plotting: Whether to enable the plotting for the maximum principal stress contour
    """

    # load the mesh file
    mesh = Mesh(mesh_file)
    # Defining a discontinuous scalar field constant on each element
    P0_scalar = FunctionSpace(mesh, 'DG', 0)
    # Calculating sigma_xx
    stress00 = Function(P0_scalar)
    stress00.assign(stress.sub(0))
    # Calculating sigma_yy
    stress11 = Function(P0_scalar)
    stress11.assign(stress.sub(3))
    # Calculating sigma_xy
    stress01 = Function(P0_scalar)
    stress01.assign(stress.sub(1))
    sigmap1 = Function(P0_scalar)
    sigmap1.assign(((stress00 + stress11) / 2) +
                   (((stress00 - stress11) / 2)**2 + stress01**2)**.5)

    # plotting function
    if plotting == True:
        name = mesh_file[:-3]+'png'
        fig, axes = plt.subplots()
        tc = tricontourf(sigmap1, levels=50, axes=axes)
        cb = plt.colorbar(tc, ax=axes)
        plt.savefig(name)
        # clear the output to avoid lack of memory
        plt.close(fig)
        plt.close("all")
    return sigmap1


def compute_failure_modes(meshfile, sigmap1):
    """
    compute the mode 1 failure parameter.

    Parameters:
    -----------
    mesh_file : the mesh file for the pellet
    sigmap1: the max principal stress  obtained from function "extract_stress_components"
    """

    # define the tensile stress value
    ft = 5.07e6
    # load the mesh file
    mesh = Mesh(meshfile)
    # Defining a discontinuous scalar field constant on each element
    P0_scalar = FunctionSpace(mesh, 'DG', 0)
    # Compute the mode 1 failure parameters. If it is < 1, the mesh will not break.
    mode1fail = Function(P0_scalar)
    mode1fail.assign(((sigmap1 >= 0) * sigmap1) / ft)

    return mode1fail


def checking_load_failure(mode1fail):
    """
    check the mode 1 failure parameter on the mesh elements to check if the pellet breaks.

    Parameters:
    -----------
    mode1fail: the mode 1 failure parameters obtained from function "compute_failure_modes"
    """

    # varaible set to recording the maximum mode 1 failure parameter
    max_mode1fail = 0
    mode1fail = mode1fail
    mode1fail_dat = mode1fail.dat.data
    # counter for the number of element with mode 1 failure parameter > 1
    ct = 0
    # check every point's mode 1 failure parameter
    for i in range(len(mode1fail_dat)):
        if (mode1fail_dat[i] >= 1):
            ct += 1
            # update the maximum mode 1 failure number
            if mode1fail_dat[i] > max_mode1fail:
                max_mode1fail = mode1fail_dat[i]
    # no element with mode 1 failure parameter > 1
    if ct == 0:
        return False
    else:
        return max_mode1fail


def find_py(mesh_file, initial_value=10, tol=1e-4, plotting=False):
    """
    Find the critical stress for the mode 1 failure for the pellet in the given mesh file.

    Parameters:
    -----------
    mesh_file : the mesh file for the pellet
    initial_value: the initial value of stress for the breaking test
    tol: tolerance for the mode 1 failure parameter residual 
    plotting: enable the plotting of the contour at the critical stress
    """

    # setting upper and lower bound for the binary search
    lower_bound = 0
    upper_bound = initial_value
    # load the mesh
    meshfile = Mesh(mesh_file)

    # compute the stress, sigma_max and the mode 1 failure value
    st = compression(mesh_file, upper_bound)
    sig = extract_stress_components(mesh_file, st)
    mode1fail = compute_failure_modes(mesh_file, sig)

    # searching the order of the magnitude of the py
    while not checking_load_failure(mode1fail):
        upper_bound *= 10
        st = compression(mesh_file, upper_bound)
        sig = extract_stress_components(mesh_file, st)
        mode1fail = compute_failure_modes(mesh_file, sig)

    # Now the upper limit has been obtained

    while True:
        # compare the mid point between lower bound and upper bound
        test_value = lower_bound+0.5*(upper_bound-lower_bound)

        st = compression(mesh_file, test_value)
        sig = extract_stress_components(mesh_file, st)
        mode1fail = compute_failure_modes(mesh_file, sig)

        # if test value smaller than the true value, test value becomes the lower bound
        if not checking_load_failure(mode1fail):
            lower_bound = test_value
        # if test value greater than the true value, test value becomes the uuper bound
        else:
            upper_bound = test_value

        # check if the accuracy met(compare to the tolerance) and if the pellet breaks
        if (checking_load_failure(mode1fail)-1) <= tol and checking_load_failure(mode1fail) > 1:
            break

    # Plot the coutour if required
    if plotting == True:
        st = compression(mesh_file, upper_bound)
        # run the function extract_stress_components with plotting enabled
        extract_stress_components(mesh_file, st, plotting=True)

    return upper_bound, checking_load_failure(mode1fail)


def calculation(outer_r, hole_r, mesh_file_list, plotting=False):
    """
    Create the two tables, one with the geometric features and the angle rotated, the critical loading and the mode 1 failure parameter.
    The other one is the aggregate version for the former one, where only the the geometric features and the maxmimum/minimum critical loading would be recorded.

    Parameters:
    -----------
    outer_r: the pellet radius
    hole_r: the hole radius
    mesh_file_list: name list of the mesh files
    plotting: enable the plotting of the contour at the critical stress
    """

    # Create the lists to record the parameters
    hole_number = []
    hole_radius = []
    dis_from_centre = []
    angle_list = []
    outer_r = outer_r
    # list for the stress
    py_list = []
    # list for the mode 1 failure parameter
    fail_list = []
    # list for the force
    force_list = []
    # calculating the small top arc length contacted with the loading plate
    length_arc = outer_r*pi/180

    ct = 0
    # find the maximum loading stress for every mesh
    for i in mesh_file_list:

        # hole_number_end, end point for the hole number in the name
        hne = i.find(' ')

        # hole_radius_start & hole_radius_end, start and end points for the hole number in the name
        hrs = i.find('=')+1
        hre = i.find('_')
        # distance_from_centre_start & distance_from_centre_end, start and end points
        # for the distance_from_centre in the name
        dfcs = i.find('c')+1
        dfce = i.find('_', hre+1)

        # append the parameters
        hole_number.append(int(i[:hne]))
        hole_radius.append(float(i[hrs:hre]))
        dis_from_centre.append(float(i[dfcs:dfce]))
        angle_list.append(int(i[-14:-11]))
        pma, failma = find_py(i, plotting=plotting)
        py_list.append(pma)
        fail_list.append(failma)

    # calculating the force
    force_list = [length_arc*py_list[i] for i in range(len(py_list))]
    outer_radius_list = [outer_r]*len(py_list)

    table = pd.DataFrame({'pellet_radius': outer_radius_list, 'hole_number': hole_number, 'hole_radius': hole_radius, 'distance_from_centre': dis_from_centre,
                          'angle': angle_list, 'force': force_list, 'max_failure_parameter': fail_list})
    table.to_csv('raw result.csv')

    # below will aggregate all the combinations with same geometric features thus find the max and min loading

    # list for all possible hole numbers
    hn = table['hole_number'].value_counts().index.tolist()
    # list for all possible hole radii
    hr = table['hole_radius'].value_counts().index.tolist()
    # list for all possible distances from centre
    dfc = table['distance_from_centre'].value_counts().index.tolist()
    # list for all possible pellet radius
    pr = table['pellet_radius'].value_counts().index.tolist()

    # making list to append and to make the final table
    hn_final = []
    hr_final = []
    dfc_final = []
    max_f = []
    min_f = []
    max_f_fail = []
    min_f_fail = []
    pr_final = []
    # for each combination
    for pra in pr:
        for num in hn:
            for ra in hr:
                for dis in dfc:
                    # filter the table with the current option of parameters
                    temp = table[table['hole_number'] == num]
                    temp = temp[temp['hole_radius'] == ra]
                    temp = temp[temp['distance_from_centre'] == dis]
                    temp = temp[temp['pellet_radius'] == pra]
                    # check if there is any value in the temporary table
                    if len(temp) != 0:
                        hn_final.append(num)
                        hr_final.append(ra)
                        dfc_final.append(dis)
                        pr_final.append(pra)
                        # Find the max loading and the corresponding failure parameter
                        max_f.append(temp['force'].max())
                        max_f_fail.append(
                            temp['max_failure_parameter'][temp["force"].idxmax()])
                        # Find the min loading and the corresponding failure parameter
                        min_f.append(temp['force'].min())
                        min_f_fail.append(
                            temp['max_failure_parameter'][temp["force"].idxmin()])

    # make new table showing only the max and min loading (no information about the orientation/angle )
    data = pd.DataFrame({'pellet_radius': pr_final, 'hole_number': hn_final, 'hole_radius': hr_final, 'distance_from_centre': dfc_final,
                         'maximum_load': max_f, 'failure parameter at maximum_load': max_f_fail,
                         'minimum_load': min_f, 'failure parameter at minimum_load': min_f_fail})
    # sort the values
    data = data.sort_values(['pellet_radius', 'hole_number', 'hole_radius',
                            'distance_from_centre'], ascending=[True, True, True, True])
    data.reset_index(drop=True, inplace=True)
    data.to_csv('aggregate result.csv')


def prediction(outer_r, surface_area):
    """
    create possible combinations for the given pellet radius and the surface area, predict the strength for every combination 
    and display geometric features for the strongest shape.

    Parameters:
    -----------
    outer_r: the pellet radius
    surface_area : the sum of the perimeter (pellet+all the holes)
    """

    # The sum of the perimeters(pellet + holes' sum)
    residual = surface_area/(2*pi)-outer_r
    # if residual smller than 0, not applicable
    if residual < 0:
        print('Inputs are invalid for the pellet!')
    else:
        # create list for geometric features
        hole_number_list = []
        ratio_pellet_hole_list = []
        hole_radius_list = []
        dfc_list = []

        # one hole
        hole_radius = residual
        if hole_radius < outer_r:
            hole_number_list.append(1)
            ratio_pellet_hole_list.append(outer_r/residual)
            hole_radius_list.append(residual)
            dfc_list.append(0)

        # find maximum number of holes
        # setting constraint: hole radius must be greater than pellet radius/10
        n = 2
        while residual/n > outer_r/10:
            n += 1

        # looping to find the all possible combinations
        for number_of_holes in range(2, n):

            # calculating the spacing betwwen adjacent holes sufficient
            hole_radius = residual/number_of_holes
            residual_angle = (180-(360/number_of_holes))/2
            d_from_c = hole_radius*2
            cc_distance = d_from_c*cos(residual_angle*pi/180)*2
            # check the spacing betwwen adjacent holes sufficient
            while 1.5*hole_radius+d_from_c <= outer_r:

                if 2.5*hole_radius <= cc_distance:
                    hole_number_list.append(number_of_holes)
                    ratio_pellet_hole_list.append(outer_r/hole_radius)
                    hole_radius_list.append(hole_radius)
                    dfc_list.append(d_from_c)
                d_from_c += hole_radius
                cc_distance = d_from_c*cos(residual_angle*pi/180)*2

        surface_area_list = [surface_area]*len(hole_radius_list)

        # make dataframe for the input of the model
        xx = pd.DataFrame({'surface_area': surface_area, 'hole_number': hole_number_list,
                           'ratio pellet_hole': ratio_pellet_hole_list, 'hole_radius': hole_radius_list,
                           'distance_from_centre': dfc_list})
        # using max and min model to make prediction
        max_model = pickle.load(open('max_model.pkl', 'rb'))
        min_model = pickle.load(open('min_model.pkl', 'rb'))

        # Find the maximum value and the index to locate the position of the geometric features in the list
        max_predict = max(max_model.predict(xx))
        max_index = list(max_model.predict(xx)).index(max_predict)

        # printing results
        print('Geometric faetures for the maximum strength:')
        print('hole number = ', hole_number_list[max_index])
        print('hole radius = ', round(
            hole_radius_list[max_index], 2), '(', hole_radius_list[max_index], ')')
        print('distance from centre = ', round(
            dfc_list[max_index], 2), '(', dfc_list[max_index], ')')
        print('Maximum Predicted Strength : ', round(
            max_model.predict(xx)[max_index]))
        print('Minimum Predicted Strength : ', round(
            min_model.predict(xx)[max_index]))