from functions import *
import pickle
print('========================================================')
print('please select the functionality')
print('1. Mesh Generation')
print('2. Mesh Generation and Maximum Loading Calculation')
print('3. Maximum and Minimum Loading Prediction')
print('4. Quit')
x = input('Please enter the option: ')
print('========================================================')

if x == '1':
    print('Mesh Generation chosen!')
    out_r = float(input('Please enter the pellet radius(unit:m):'))
    hol_r = float(input('Please enter the hole radius(unit:m):'))
    m_list = create_mesh(hol_r, out_r)

    print('Success! Please find the folder for the output meshes')


elif x == '2':
    print('Mesh Generation and Maximum Loading Calculation chosen!')
    out_r = float(input('Please enter the pellet radius(unit:m):'))
    hol_r = float(input('Please enter the hole radius(unit:m):'))
    m_list = create_mesh(hol_r, out_r)
    calculation(out_r, hol_r, m_list)
    print('Success! Please find the folder for the output tables')

elif x == '3':
    print('Maximum and Minimum Loading Prediction chosen!')
    surface_area = float(input('Please enter the surface area(unit:m):'))
    out_r = float(input('Please enter the pellet radius(unit:m):'))
    prediction(out_r,surface_area)

elif x == '4':
    pass
else:
    print('Invalid Input!')