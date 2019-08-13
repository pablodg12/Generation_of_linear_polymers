from sympy import *
import numpy as np
import move as mv
import copy
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import distance_matrix
from numpy import inf
import sys
import random
import math

def calculate_end_to_end(data,L,w=0,get_chain=False):
    L = L[1] - L[0]
    #L1 = [-2.3708598967242978e+00,6.9762359896724035e+01]
    tmp_chain = []
    for i,v in enumerate(data):
        
        x = v[0]
        y = v[1]
        z = v[2]
        if i > 0 :
            a = int(np.abs(tmp_chain[-1][0]-x)/L)
            b = int(np.abs(tmp_chain[-1][1]-y)/L)
            c = int(np.abs(tmp_chain[-1][2]-z)/L)
            for m in range(0,a+1):
                if tmp_chain[-1][0]-x > L/2:
                    x = x + L
                elif tmp_chain[-1][0]-x < -L/2:
                    x = x - L
            for m in range(0,b+1):
                if tmp_chain[-1][1]-y > L/2:
                    y = y + L
                elif tmp_chain[-1][1]-y < -L/2:
                    y = y - L
            for m in range(0,c+1):
                if tmp_chain[-1][2]-z > L/2:
                    z = z + L
                elif tmp_chain[-1][2]-z < -L/2:
                    z = z - L
            tmp_chain.append([x,y,z])
        else:
            tmp_chain.append([x,y,z])
    if w==0:
        ee = np.linalg.norm(np.array(tmp_chain[0])-np.array(tmp_chain[-1]))**2
    else:
        ee = np.linalg.norm(np.array(tmp_chain[0])-np.array(tmp_chain[w]))**2
    if get_chain == False:
        return(ee)
    else:
        return(tmp_chain)
        


def algoritmo_test(chain,d):
    vm = calculate_end_to_end(chain,d,0,True)
    theta = []
    angles = []
    radius = []
    for k in range(0,len(vm)-3):
        g = vm.copy()
        point_a = np.array([g[k][0],g[k][1],g[k][2]])
        point_b = np.array([g[k+1][0],g[k+1][1],g[k+1][2]])
        point_c = np.array([g[k+2][0],g[k+2][1],g[k+2][2]])
        point_d = np.array([g[k+3][0],g[k+3][1],g[k+3][2]])
        
        point_d = mv.traslation(point_d,-1*point_c)
        point_b = mv.traslation(point_b,-1*point_c)
        point_a = mv.traslation(point_a,-1*point_c)
        point_c = mv.traslation(point_c,-1*point_c)

        rot_matrix_x = mv.rotation_matrix_x(point_b,0,center=True)
        point_a = np.dot(rot_matrix_x,point_a)
        point_b = np.dot(rot_matrix_x,point_b)
        point_c = np.dot(rot_matrix_x,point_c)
        point_d = np.dot(rot_matrix_x,point_d)

        rot_matrix_y = mv.rotation_matrix_y(point_b,0,center=True)
        point_a = np.dot(rot_matrix_y,point_a)
        point_b = np.dot(rot_matrix_y,point_b)
        point_c = np.dot(rot_matrix_y,point_c)
        point_d = np.dot(rot_matrix_y,point_d)


        rot_matrix_z = mv.rotation_matrix_z(point_a,0,center=True)
        point_a = np.dot(rot_matrix_z,point_a)
        point_b = np.dot(rot_matrix_z,point_b)
        point_c = np.dot(rot_matrix_z,point_c)
        point_d = np.dot(rot_matrix_z,point_d)
        
        a,b,c = mv.cart2sph(point_d[0], point_d[1], point_d[2])
        radius.append(a)
        theta.append(b)
        angles.append(c)
        #angles.append(cart2sph(point_d[0], point_d[1], point_d[2])[2])
    
    return radius,theta,angles



def Fm(x):
    
    phi0 = x[0]
    phi1 = x[1]
    phi2 = x[2]
    phi4 = x[3]
    phi5 = x[4]
    phi6 = x[5]
    
    _,_,_,_,_,_,_,p6a,p6b,p6c =  mv.get_spheric_2(cadena[1:5])
    _,_,_,_,_,_,_,p0a,p0b,p0c =  mv.get_spheric_2([cadena[k] for k in reversed(range(8,12))])
    _,_,_,X0,Y0,Z0,T0 = mv.get_spheric([cadena[k] for k in reversed(range(8,12))])
    _,_,_,X6,Y6,Z6,T6 = mv.get_spheric(cadena[1:5])
    
    t = cadena[2:7]
    
    w_t = np.array([1.54*np.sin(1.911)*np.cos(phi6),1.54*np.sin(1.911)*np.sin(phi6),1.54*np.cos(1.911)])
    t[2] = np.ndarray.tolist(np.dot(X6.T,np.dot(Y6.T,np.dot(Z6.T,w_t)))+T6)
    
    _,_,_,_,_,_,_,p5a,p5b,p5c = mv.get_spheric_2(t[0:4])
    _,_,_,X5,Y5,Z5,T5 = mv.get_spheric(t[0:4])
    
    w_t = np.array([1.54*np.sin(1.911)*np.cos(phi5),1.54*np.sin(1.911)*np.sin(phi5),1.54*np.cos(1.911)])
    t[3] = np.ndarray.tolist(np.dot(X5.T,np.dot(Y5.T,np.dot(Z5.T,w_t)))+T5)    
    _,_,_,_,_,_,_,p4a,p4b,p4c = mv.get_spheric_2(t[1:])

    t = cadena[6:11]
    
    w_t = np.array([1.54*np.sin(1.911)*np.cos(phi0),1.54*np.sin(1.911)*np.sin(phi0),1.54*np.cos(1.911)])
    t[2] = np.ndarray.tolist(np.dot(X0.T,np.dot(Y0.T,np.dot(Z0.T,w_t)))+T0)
    
    _,_,_,_,_,_,_,p1a,p1b,p1c = mv.get_spheric_2([t[k] for k in reversed(range(1,len(t)))])
    _,_,_,X1,Y1,Z1,T1 = mv.get_spheric([t[k] for k in reversed(range(1,len(t)))])
    
    w_t = np.array([1.54*np.sin(1.911)*np.cos(phi1),1.54*np.sin(1.911)*np.sin(phi1),1.54*np.cos(1.911)])
    t[1] = np.ndarray.tolist(np.dot(X1.T,np.dot(Y1.T,np.dot(Z1.T,w_t)))+T1)
    _,_,_,_,_,_,_,p2a,p2b,p2c = mv.get_spheric_2([t[k] for k in reversed(range(0,len(t)-1))])
    

    
    phi0 = phi0
    phi1 = phi1
    phi2 = phi2
    phi4 = phi4
    phi5 = phi5
    phi6 = phi6
    
    theta1_p0 = p0a
    theta2_p0 = p0b
    theta3_p0 = p0c
    
    theta1_p1 = p1a
    theta2_p1 = p1b
    theta3_p1 = p1c
    
    theta1_p2 = p2a
    theta2_p2 = p2b
    theta3_p2 = p2c
    
    theta1_p4 = p4a
    theta2_p4 = p4b
    theta3_p4 = p4c
    
    theta1_p5 = p5a
    theta2_p5 = p5b
    theta3_p5 = p5c
    
    theta1_p6 = p6a
    theta2_p6 = p6b
    theta3_p6 = p6c
    
    p0x=cadena[9][0]
    p0y=cadena[9][1]
    p0z=cadena[9][2]
    
    p6x=cadena[3][0]
    p6y=cadena[3][1]
    p6z=cadena[3][2]
    
    a,b = F(phi0,phi1,phi2,phi4,phi5,phi6,theta1_p0,theta2_p0,theta3_p0,theta1_p6,theta2_p6,theta3_p6,
               theta1_p1,theta2_p1,theta3_p1,theta1_p2,theta2_p2,theta3_p2,
               theta1_p5,theta2_p5,theta3_p5,theta1_p4,theta2_p4,theta3_p4,p0x,p0y,p0z,p6x,p6y,p6z)
    return(np.linalg.norm(np.array([a[0][0],a[1][0],a[2][0],b])))


def convert_angles(angle):
    if angle >= 0 and angle <= np.pi/2:
        return(np.pi/2-angle)
    if angle > np.pi/2 and angle <= np.pi:
        return((angle-np.pi/2))
    if angle > np.pi and angle <= 3*np.pi/2:
        return(angle-np.pi/2)
    if angle > 3*np.pi/2 and angle <= 2*np.pi:
        return((2*np.pi-angle)+np.pi/2)

def potential(all_data,L,pr=False,save=-1):
    data = copy.deepcopy(all_data)
    energyA = 0
    energyB = 0
    R0 = 0.112
    E = 4.01
    r_cut = 10.5
    #L = 26.3829

    phi = []

    for t,k in enumerate(data):
        k = np.array(k)
        _,_,c = algoritmo_test(k,[0,L])
        phi = phi + c
    for k in phi:
        k = k - np.pi/2
        if k < 0:
            a = 2*np.pi + k
        else:
            a = k
            
        energyA += 1.736*np.cos(a)**0 - 4.490*np.cos(a)**1 + 0.776*np.cos(a)**2 + 6.990*np.cos(a)**3

    data = np.array(data)

    if save != -1:
        np.save("dihedral_data_"+str(save)+".npy",phi)
        
    for a in range(0,data.shape[0]):
        for b in range(0,data.shape[1]): 
            for c in range(0,data.shape[0]):
                for d in range(0,data.shape[1]):
                    if (a == c) and (d == b-3 or d == b-2 or d == b-1 or d == b or d == b+1 or d == b+2 or d == b+3):
                        continue
                    else:
                        p1 = copy.deepcopy(data[a,b])
                        p2 = copy.deepcopy(data[c,d])
                        if (p1[0]-p2[0]) > L/2:
                            p2[0] = p2[0] + L 
                        elif (p1[0]-p2[0]) < -L/2:
                            p2[0] = p2[0] - L
                        if (p1[1]-p2[1]) > L/2:
                            p2[1] = p2[1] + L 
                        elif (p1[1]-p2[1]) < -L/2:
                            p2[1] = p2[1] - L 
                        if (p1[2]-p2[2]) > L/2:
                            p2[2] = p2[2] + L 
                        elif (p1[2]-p2[2]) < -L/2:
                            p2[2] = p2[2] - L                       
                        r = np.linalg.norm(p1-p2)
                        if r < 10.5:
                            energyB += 0.448*((4.01/r)**12 - (4.01/r)**6)
    energynofake = energyA
    energyA = energyA*10**(int(math.log10(abs(energyB))-1))
    if pr == True:
        print("Real Dihedral Energy: ",energynofake)
        print("Dihedral Energy Fake: ",energyA)
        print("LJ Energy: ",energyB)
        print("Total Energy ",energyA+energyB)
        #print("...................")

    if energyB > 0:
        return(energyA+energyB)
    else:
        return(energyA)

if __name__ == '__main__':
    #np.random.seed(0)
    #random.seed(0)
    phi, theta = symbols('phi theta') #Angulos genericos
    phi0, phi1, phi2, phi4, phi5, phi6 = symbols('phi_0 phi_1 phi_2 phi_4 phi_5 phi_6') #Angulos
    theta1_p0,theta2_p0,theta3_p0 = symbols('theta1_p0 theta2_p0 theta3_p0')
    theta1_p1,theta2_p1,theta3_p1 = symbols('theta1_p1 theta2_p1 theta3_p1')
    theta1_p2,theta2_p2,theta3_p2 = symbols('theta1_p2 theta2_p2 theta3_p2')
    theta1_p4,theta2_p4,theta3_p4 = symbols('theta1_p4 theta2_p4 theta3_p4')
    theta1_p5,theta2_p5,theta3_p5 = symbols('theta1_p5 theta2_p5 theta3_p5')
    theta1_p6,theta2_p6,theta3_p6 = symbols('theta1_p6 theta2_p6 theta3_p6')
    p0x,p0y,p0z = symbols('p0x p0y p0z')
    p6x,p6y,p6z = symbols('p6x p6y p6z')

    w  = Matrix ([1.54*sin(1.911)*cos(phi),1.54*sin(1.911)*sin(phi),1.54*np.cos(1.911) ]) #param del cono
    Ry = Matrix([[cos(theta),0,sin(theta)],[0,1,0],[-1*sin(theta),0,cos(theta)]]) #matriz de rotación xz
    Rx = Matrix([[1,0,0],[0,cos(theta),-sin(theta)],[0,sin(theta),cos(theta)]]) #matriz de rotación xz
    Rz = Matrix([[cos(theta),-sin(theta),0],[sin(theta),cos(theta),0],[0,0,1]])

    p1 = Rx.T.subs(theta,theta1_p0)*(Ry.T.subs(theta,theta2_p0)*(Rz.T.subs(theta,theta3_p0)*w.subs(phi,phi0))) + Matrix([p0x,p0y,p0z])
    p5 = Rx.T.subs(theta,theta1_p6)*(Ry.T.subs(theta,theta2_p6)*(Rz.T.subs(theta,theta3_p6)*w.subs(phi,phi6))) +  Matrix([p6x,p6y,p6z])

    p2 = Rx.T.subs(theta,theta1_p1)*(Ry.T.subs(theta,theta2_p1)*(Rz.T.subs(theta,theta3_p1)*w.subs(phi,phi1))) + p1
    p4 = Rx.T.subs(theta,theta1_p5)*(Ry.T.subs(theta,theta2_p5)*(Rz.T.subs(theta,theta3_p5)*w.subs(phi,phi5))) + p5

    p3_mas = Rx.T.subs(theta,theta1_p4)*(Ry.T.subs(theta,theta2_p4)*(Rz.T.subs(theta,theta3_p4)*w.subs(phi,phi4))) + p4
    p3_menos = Rx.T.subs(theta,theta1_p2)*(Ry.T.subs(theta,theta2_p2)*(Rz.T.subs(theta,theta3_p2)*w.subs(phi,phi2))) + p2

    restriction = cos(1.911) - (p3_mas-p4).dot(p3_menos-p2)/((p3_mas-p4).norm() * (p3_menos-p2).norm())

    F = [simplify(trigsimp(p3_mas - p3_menos)),simplify(trigsimp(restriction))]

    F =  lambdify((phi0,phi1,phi2,phi4,phi5,phi6,theta1_p0,theta2_p0,theta3_p0,theta1_p6,theta2_p6,theta3_p6,
               theta1_p1,theta2_p1,theta3_p1,theta1_p2,theta2_p2,theta3_p2,
               theta1_p5,theta2_p5,theta3_p5,theta1_p4,theta2_p4,theta3_p4,p0x,p0y,p0z,p6x,p6y,p6z),F,'numpy')

    #datas = np.load("data1.npy")
    #all_data = [datas[0:300],datas[300:]]
    datas = np.load("data100_300_300K.npy")

    total_data = copy.deepcopy(datas)

    steps = 10000000000000
    energia_p = []
    end_to_end = []

    n_atoms = 300 + 6
    tol = 1e-5

    L = 26.3829
    ee_target = 7586.304102021168

    gt_fast = 1
    lo = 1.54
    Tho = 1.911

    for j in range(0,steps):

        print("step ",j)

        if gt_fast == 1:
            if j % 50 != 0:
                a = calculate_end_to_end(total_data[0],[0,L])
                b = calculate_end_to_end(total_data[1],[0,L])
                print("End-to-End Promedio: ",(a+b)/2)
                print("..............................")
            else:
                a = calculate_end_to_end(total_data[0],[0,L])
                b = calculate_end_to_end(total_data[1],[0,L])
                end_to_end.append((a+b)/2)
                np.save("end_to_end_data.npy",end_to_end)
                print("End-to-End Promedio: ",(a+b)/2)
                print("..............................")

        if gt_fast == 1 and ee_target*1.05>(a+b)/2 and ee_target*0.95<(a+b)/2:
            print("Active minimization")
            all_data = copy.deepcopy(total_data)
            gt_fast = 0
            n_atoms = 300

        if gt_fast == 1:
            for k_i in range(0,2):
                for i in range(0,n_atoms-12):

                    all_data = np.zeros((total_data.shape[0],total_data.shape[1]+6,3))
                    all_data[k_i][3:total_data.shape[1]+3] = copy.deepcopy(total_data[k_i])

                    Th = np.random.uniform(0,1)*np.pi
                    Ph = np.random.uniform(0,2*np.pi,1)

                    all_data[k_i][2][0] = all_data[k_i][3][0] + lo*np.sin(Th)*np.cos(Ph) 
                    all_data[k_i][2][1] = all_data[k_i][3][1] + lo*np.sin(Th)*np.sin(Ph)  
                    all_data[k_i][2][2] = all_data[k_i][3][2] + lo*np.cos(Th)   

                    temporal_0 = mv.traslation(all_data[k_i][3],-all_data[k_i][2])
                    rot_matrix_x = mv.rotation_matrix_x(temporal_0,0,center=True)
                    temporal_1 = np.dot(rot_matrix_x,temporal_0)
                    rot_matrix_y = mv.rotation_matrix_y(temporal_1,0,center=True)

                    r_phi = np.random.uniform(0,2*np.pi,1)
                    all_data[k_i][1][0] = lo*np.sin(Tho)*np.cos(r_phi) 
                    all_data[k_i][1][1] = lo*np.sin(Tho)*np.sin(r_phi) 
                    all_data[k_i][1][2] = lo*np.cos(Tho) 
                    all_data[k_i][1] = mv.traslation(np.dot(np.linalg.inv(rot_matrix_x),np.dot(np.linalg.inv(rot_matrix_y),all_data[k_i][1])),all_data[k_i][2])

                    point_a = all_data[k_i][3]
                    point_b = all_data[k_i][2]
                    point_c = all_data[k_i][1]

                    temporal_0 = mv.traslation(point_b,-point_c)
                    temporala_0 = mv.traslation(point_a,-point_c)
                    rot_matrix_x = mv.rotation_matrix_x(temporal_0,0,center=True)
                    temporal_1 = np.dot(rot_matrix_x,temporal_0)
                    temporala_1 = np.dot(rot_matrix_x,temporala_0)
                    rot_matrix_y = mv.rotation_matrix_y(temporal_1,0,center=True)
                    temporal_2 = np.dot(rot_matrix_y,temporal_1)
                    temporala_2 = np.dot(rot_matrix_y,temporala_1)
                    rot_matrix_z = mv.rotation_matrix_z(temporala_2,0,center=True)

                    r_phi = np.random.uniform(0,2*np.pi,1)[0]
                    all_data[k_i][0][0] = lo*np.sin(Tho)*np.cos(r_phi) 
                    all_data[k_i][0][1] = lo*np.sin(Tho)*np.sin(r_phi) 
                    all_data[k_i][0][2] = lo*np.cos(Tho)
                    all_data[k_i][0] = mv.traslation(np.dot(np.linalg.inv(rot_matrix_x),np.dot(np.linalg.inv(rot_matrix_y),np.dot(np.linalg.inv(rot_matrix_z),all_data[k_i][0]))),all_data[k_i][1])

                    for asd in range(0,3):
                        for dsa in range(0,3):
                            if all_data[k_i][asd][dsa] > L:
                                all_data[k_i][asd][dsa] = all_data[k_i][asd][dsa] - L
                            elif all_data[k_i][asd][dsa] < 0:
                                all_data[k_i][asd][dsa] = all_data[k_i][asd][dsa] + L

                    #########################

                    Th = np.random.uniform(0,1)*np.pi
                    Ph = np.random.uniform(0,2*np.pi,1)

                    all_data[k_i][total_data.shape[1]+3][0] = all_data[k_i][total_data.shape[1]-1+3][0] + lo*np.sin(Th)*np.cos(Ph) 
                    all_data[k_i][total_data.shape[1]+3][1] = all_data[k_i][total_data.shape[1]-1+3][1] + lo*np.sin(Th)*np.sin(Ph)  
                    all_data[k_i][total_data.shape[1]+3][2] = all_data[k_i][total_data.shape[1]-1+3][2] + lo*np.cos(Th)   

                    temporal_0 = mv.traslation(all_data[k_i][total_data.shape[1]-1+3],-all_data[k_i][total_data.shape[1]+3])
                    rot_matrix_x = mv.rotation_matrix_x(temporal_0,0,center=True)
                    temporal_1 = np.dot(rot_matrix_x,temporal_0)
                    rot_matrix_y = mv.rotation_matrix_y(temporal_1,0,center=True)

                    r_phi = np.random.uniform(0,2*np.pi,1)
                    all_data[k_i][total_data.shape[1]+1+3][0] = lo*np.sin(Tho)*np.cos(r_phi) 
                    all_data[k_i][total_data.shape[1]+1+3][1] = lo*np.sin(Tho)*np.sin(r_phi) 
                    all_data[k_i][total_data.shape[1]+1+3][2] = lo*np.cos(Tho) 
                    all_data[k_i][total_data.shape[1]+1+3] = mv.traslation(np.dot(np.linalg.inv(rot_matrix_x),np.dot(np.linalg.inv(rot_matrix_y),all_data[k_i][total_data.shape[1]+1+3])),all_data[k_i][total_data.shape[1]+3])

                    point_a = all_data[k_i][total_data.shape[1]-1+3]
                    point_b = all_data[k_i][total_data.shape[1]+3]
                    point_c = all_data[k_i][total_data.shape[1]+1+3]

                    temporal_0 = mv.traslation(point_b,-point_c)
                    temporala_0 = mv.traslation(point_a,-point_c)
                    rot_matrix_x = mv.rotation_matrix_x(temporal_0,0,center=True)
                    temporal_1 = np.dot(rot_matrix_x,temporal_0)
                    temporala_1 = np.dot(rot_matrix_x,temporala_0)
                    rot_matrix_y = mv.rotation_matrix_y(temporal_1,0,center=True)
                    temporal_2 = np.dot(rot_matrix_y,temporal_1)
                    temporala_2 = np.dot(rot_matrix_y,temporala_1)
                    rot_matrix_z = mv.rotation_matrix_z(temporala_2,0,center=True)

                    r_phi = np.random.uniform(0,2*np.pi,1)[0]
                    all_data[k_i][total_data.shape[1]+2+3][0] = lo*np.sin(Tho)*np.cos(r_phi) 
                    all_data[k_i][total_data.shape[1]+2+3][1] = lo*np.sin(Tho)*np.sin(r_phi) 
                    all_data[k_i][total_data.shape[1]+2+3][2] = lo*np.cos(Tho)
                    all_data[k_i][total_data.shape[1]+2+3] = mv.traslation(np.dot(np.linalg.inv(rot_matrix_x),
                        np.dot(np.linalg.inv(rot_matrix_y),np.dot(np.linalg.inv(rot_matrix_z),all_data[k_i][total_data.shape[1]+2+3]))),all_data[k_i][total_data.shape[1]+1+3])

                    for asd in range(total_data.shape[1]-1,total_data.shape[1]+2+3+1):
                        for dsa in range(0,3):
                            if all_data[k_i][asd][dsa] > L:
                                all_data[k_i][asd][dsa] = all_data[k_i][asd][dsa] - L
                            elif all_data[k_i][asd][dsa] < 0:
                                all_data[k_i][asd][dsa] = all_data[k_i][asd][dsa] + L
                    #print("Energía Potencial..: ",potential_energy_1)
                    p = 1
                    v = []
                    #i = np.random.randint(0,n_atoms-12)
                    data = copy.deepcopy(all_data[k_i])

                    tmp_chain = []
                    cadena = copy.deepcopy(data[i:i+13])

                    for ii,vv in enumerate(cadena):
                        x = vv[0]
                        y = vv[1]
                        z = vv[2]
                        if ii > 0 :
                            a = int(np.abs(tmp_chain[-1][0]-x)/L)
                            b = int(np.abs(tmp_chain[-1][1]-y)/L)
                            c = int(np.abs(tmp_chain[-1][2]-z)/L)
                            for m in range(0,a+1):
                                if tmp_chain[-1][0]-x > L/2:
                                    x = x + L
                                if tmp_chain[-1][0]-x < -L/2:
                                    x = x - L
                            for m in range(0,b+1):
                                if tmp_chain[-1][1]-y > L/2:
                                    y = y + L
                                if tmp_chain[-1][1]-y < -L/2:
                                    y = y - L
                            for m in range(0,c+1):
                                if tmp_chain[-1][2]-z > L/2:
                                    z = z + L
                                if tmp_chain[-1][2]-z < -L/2:
                                    z = z - L
                            tmp_chain.append([x,y,z])
                        else:
                            tmp_chain.append([x,y,z])
                    cadena = copy.deepcopy(tmp_chain)
                    counter = 0
                    while(p>tol):
                        counter+=1
                        if counter == 10:
                            break
                        x0 = np.random.uniform(0,2*np.pi,6)
                        try:
                            gg = minimize(Fm,x0, tol=1e-6)
                            p = gg.fun
                            for h in gg.x:
                                if h < 0:
                                    v.append(2*np.pi+h)
                                elif h > 2*np.pi:
                                    v.append(h-2*np.pi)
                                else: 
                                    v.append(h)
                        except RuntimeWarning:
                            continue
                    if counter == 10:
                        continue
                    _,_,_,X0,Y0,Z0,T0 = mv.get_spheric([cadena[k] for k in reversed(range(8,12))])
                    _,_,_,X6,Y6,Z6,T6 = mv.get_spheric(cadena[1:5])

                    w_t = np.array([1.54*np.sin(1.911)*np.cos(v[5]),1.54*np.sin(1.911)*np.sin(v[5]),1.54*np.cos(1.911)])
                    cadena[4] = np.ndarray.tolist(np.dot(X6.T,np.dot(Y6.T,np.dot(Z6.T,w_t)))+T6)

                    w_t = np.array([1.54*np.sin(1.911)*np.cos(v[0]),1.54*np.sin(1.911)*np.sin(v[0]),1.54*np.cos(1.911)])
                    cadena[8] = np.ndarray.tolist(np.dot(X0.T,np.dot(Y0.T,np.dot(Z0.T,w_t)))+T0)

                    _,_,_,X1,Y1,Z1,T1 = mv.get_spheric([cadena[k] for k in reversed(range(7,11))])
                    _,_,_,X5,Y5,Z5,T5 = mv.get_spheric(cadena[2:6])

                    w_t = np.array([1.54*np.sin(1.911)*np.cos(v[1]),1.54*np.sin(1.911)*np.sin(v[1]),1.54*np.cos(1.911)])
                    cadena[7] = np.ndarray.tolist(np.dot(X1.T,np.dot(Y1.T,np.dot(Z1.T,w_t)))+T1)

                    w_t = np.array([1.54*np.sin(1.911)*np.cos(v[4]),1.54*np.sin(1.911)*np.sin(v[4]),1.54*np.cos(1.911)])
                    cadena[5] = np.ndarray.tolist(np.dot(X5.T,np.dot(Y5.T,np.dot(Z5.T,w_t)))+T5)

                    _,_,_,X4,Y4,Z4,T4 = mv.get_spheric(cadena[3:7])
                    w_t = np.array([1.54*np.sin(1.911)*np.cos(v[3]),1.54*np.sin(1.911)*np.sin(v[3]),1.54*np.cos(1.911)])
                    cadena[6] = np.ndarray.tolist(np.dot(X4.T,np.dot(Y4.T,np.dot(Z4.T,w_t)))+T4)

                    for mnq in range(0,13):
                        for uj in range(0,3):
                            if cadena[mnq][uj] > L:
                                cadena[mnq][uj] = cadena[mnq][uj]- L
                            elif cadena[mnq][uj] < 0:
                                cadena[mnq][uj]= cadena[mnq][uj] + L


                    all_data_copy = copy.deepcopy(all_data)

                    data_2 = copy.deepcopy(data)

                    data_2[i:i+13] = copy.deepcopy(cadena)

                    all_data_copy[k_i] = copy.deepcopy(data_2)

                    if gt_fast == 1:
                        eeee1 = calculate_end_to_end(all_data[k_i][3:total_data.shape[1]+3],[0,L],i+5-4)
                        eeee2 = calculate_end_to_end(all_data_copy[k_i][3:total_data.shape[1]+3],[0,L],i+5-4)
                        if (eeee2-(i+2)*ee_target/(n_atoms-6))**2 < (eeee1-(i+2)*ee_target/(n_atoms-6))**2:
                            data = copy.deepcopy(data_2)
                            del all_data
                            total_data[k_i] = copy.deepcopy(all_data_copy[k_i][3:total_data.shape[1]+3])
                            #all_data = all_data_copy
                        del data,cadena,all_data_copy,data_2
                    else:
                        potential_energy_2 = potential(all_data_copy,L)
                        if potential_energy_2 <= potential_energy_1:
                            data = copy.deepcopy(data_2)
                            del all_data
                            all_data = copy.deepcopy(all_data_copy)
                            potential_energy_1 = potential_energy_2

                        del data,cadena,all_data_copy,data_2,potential_energy_2


                    da = copy.deepcopy(total_data[k_i]) 
                    chain = copy.deepcopy(da[len(da)-4:])

                    ee1 = calculate_end_to_end(da,[0,L])

                    tmp_chain = []
                    for ii,vv in enumerate(chain):
                        x = vv[0]
                        y = vv[1]
                        z = vv[2]
                        if ii > 0 :
                            a = int(np.abs(tmp_chain[-1][0]-x)/L)
                            b = int(np.abs(tmp_chain[-1][1]-y)/L)
                            c = int(np.abs(tmp_chain[-1][2]-z)/L)
                            for m in range(0,a+1):
                                if tmp_chain[-1][0]-x > L/2:
                                    x = x + L
                                if tmp_chain[-1][0]-x < -L/2:
                                    x = x - L
                            for m in range(0,b+1):
                                if tmp_chain[-1][1]-y > L/2:
                                    y = y + L
                                if tmp_chain[-1][1]-y < -L/2:
                                    y = y - L
                            for m in range(0,c+1):
                                if tmp_chain[-1][2]-z > L/2:
                                    z = z + L
                                if tmp_chain[-1][2]-z < -L/2:
                                    z = z - L
                            tmp_chain.append([x,y,z])
                        else:
                            tmp_chain.append([x,y,z])

                    for _ in range(0,1):
                        chain = np.array(tmp_chain)
                        _,_,_,X0,Y0,Z0,T0 = mv.get_spheric(chain)
                        tmp = np.random.uniform(0,2*np.pi,1)
                        w_t = np.array([1.54*np.sin(1.911)*np.cos(tmp),1.54*np.sin(1.911)*np.sin(tmp),1.54*np.cos(1.911)])
                        
                        chain[3] = np.ndarray.tolist(np.dot(X0.T,np.dot(Y0.T,np.dot(Z0.T,w_t)))+T0)
                        for tk in range(0,4):
                            for uj in range(0,3):
                                if chain[tk][uj] > L:
                                    chain[tk][uj]  = chain[tk][uj] - L
                                elif chain[tk][uj]  < 0:
                                    chain[tk][uj] = chain[tk][uj]  + L   

                        da[len(da)-4:] = copy.deepcopy(chain)

                        ee2 = calculate_end_to_end(da,[0,L])

                        all_data_copy = copy.deepcopy(total_data)
                        all_data_copy[k_i] = copy.deepcopy(da)
                        
                        #potential_energy_2 = potential(all_data_copy)

                        
                        #if potential_energy_2 <= potential_energy_1 and (ee2 - 7586)**2 < (ee1 - 7586)**2:
                        if  (ee2 - ee_target )**2 <= (ee1 - ee_target )**2:
                            #del total_data
                            total_data[k_i] = copy.deepcopy(all_data_copy[k_i])
                            break
                            #potential_energy_1 = potential_energy_2
                            #print(ee2)
                            
                    del chain,all_data_copy,da
                        
                    data = copy.deepcopy(total_data[k_i]) 
                    cadena = copy.deepcopy(data[0:4])

                    ee1 = calculate_end_to_end(data,[0,L])

                    tmp_chain = []
                    for ii,vv in enumerate(cadena):
                        x = vv[0]
                        y = vv[1]
                        z = vv[2]
                        if ii > 0 :
                            a = int(np.abs(tmp_chain[-1][0]-x)/L)
                            b = int(np.abs(tmp_chain[-1][1]-y)/L)
                            c = int(np.abs(tmp_chain[-1][2]-z)/L)
                            for m in range(0,a+1):
                                if tmp_chain[-1][0]-x > L/2:
                                    x = x + L
                                if tmp_chain[-1][0]-x < -L/2:
                                    x = x - L
                            for m in range(0,b+1):
                                if tmp_chain[-1][1]-y > L/2:
                                    y = y + L
                                if tmp_chain[-1][1]-y < -L/2:
                                    y = y - L
                            for m in range(0,c+1):
                                if tmp_chain[-1][2]-z > L/2:
                                    z = z + L
                                if tmp_chain[-1][2]-z < -L/2:
                                    z = z - L
                            tmp_chain.append([x,y,z])
                        else:
                            tmp_chain.append([x,y,z])
                    
                    for _ in range(0,1):
                        cadena = copy.deepcopy(tmp_chain)
                        _,_,_,X0,Y0,Z0,T0 = mv.get_spheric([cadena[k] for k in reversed(range(0,4))])
                        tmp = np.random.uniform(0,2*np.pi,1)
                        w_t = np.array([1.54*np.sin(1.911)*np.cos(tmp),1.54*np.sin(1.911)*np.sin(tmp),1.54*np.cos(1.911)])
                        
                        cadena[0] = np.ndarray.tolist(np.dot(X0.T,np.dot(Y0.T,np.dot(Z0.T,w_t)))+T0)
                        for tk in range(0,4):
                            for uj in range(0,3):
                                if cadena[tk][uj] > L:
                                    cadena[tk][uj]  = cadena[tk][uj] - L
                                elif cadena[tk][uj]  < 0:
                                    cadena[tk][uj] = cadena[tk][uj]  + L 

                        data_2 = copy.deepcopy(data)
                        data_2[0:4] = copy.deepcopy(cadena)

                        ee2 = calculate_end_to_end(data_2,[0,L])

                        all_data_copy = copy.deepcopy(total_data)
                        all_data_copy[k_i] = copy.deepcopy(data_2)
                        #potential_energy_2 = potential(all_data_copy)
                        
                        #if potential_energy_2 <= potential_energy_1 and (ee2 - 7586)**2 < (ee1 - 7586)**2:
                        if  (ee2 - ee_target )**2 <= (ee1 - ee_target)**2:
                            #del all_data
                            total_data[k_i] = copy.deepcopy(all_data_copy[k_i])
                            #potential_energy_1 = potential_energy_2
                            #print(ee2)
                            break
        
        if gt_fast == 0:

            k_i = np.random.randint(0,2)

            if j % 50 != 0:
                potential_energy_1 = potential(all_data,L,True)
                a = calculate_end_to_end(all_data[0],[0,L])
                b = calculate_end_to_end(all_data[1],[0,L])
                print("End-to-End Promedio: ",(a+b)/2)
                print("..............................")
            else:
                potential_energy_1 = potential(all_data,L,True,j)
                a = calculate_end_to_end(all_data[0],[0,L])
                b = calculate_end_to_end(all_data[1],[0,L])
                end_to_end.append((a+b)/2)
                np.save("end_to_end_data.npy",end_to_end)
                print("End-to-End Promedio: ",(a+b)/2)
                print("..............................")

            energia_p.append(potential_energy_1)
            p = 1
            v = []
            i = np.random.randint(0,n_atoms-12)
            data = copy.deepcopy(all_data[k_i])

            tmp_chain = []
            cadena = copy.deepcopy(data[i:i+13])

            for ii,vv in enumerate(cadena):
                x = vv[0]
                y = vv[1]
                z = vv[2]
                if ii > 0 :
                    a = int(np.abs(tmp_chain[-1][0]-x)/L)
                    b = int(np.abs(tmp_chain[-1][1]-y)/L)
                    c = int(np.abs(tmp_chain[-1][2]-z)/L)
                    for m in range(0,a+1):
                        if tmp_chain[-1][0]-x > L/2:
                            x = x + L
                        if tmp_chain[-1][0]-x < -L/2:
                            x = x - L
                    for m in range(0,b+1):
                        if tmp_chain[-1][1]-y > L/2:
                            y = y + L
                        if tmp_chain[-1][1]-y < -L/2:
                            y = y - L
                    for m in range(0,c+1):
                        if tmp_chain[-1][2]-z > L/2:
                            z = z + L
                        if tmp_chain[-1][2]-z < -L/2:
                            z = z - L
                    tmp_chain.append([x,y,z])
                else:
                    tmp_chain.append([x,y,z])
            cadena = copy.deepcopy(tmp_chain)
            counter = 0
            while(p>tol):
                counter+=1
                if counter == 10:
                    break
                x0 = np.random.uniform(0,2*np.pi,6)
                try:
                    gg = minimize(Fm,x0, tol=1e-6)
                    p = gg.fun
                    for h in gg.x:
                        if h < 0:
                            v.append(2*np.pi+h)
                        elif h > 2*np.pi:
                            v.append(h-2*np.pi)
                        else: 
                            v.append(h)
                except RuntimeWarning:
                    continue
            if counter == 10:
                continue
            _,_,_,X0,Y0,Z0,T0 = mv.get_spheric([cadena[k] for k in reversed(range(8,12))])
            _,_,_,X6,Y6,Z6,T6 = mv.get_spheric(cadena[1:5])

            w_t = np.array([1.54*np.sin(1.911)*np.cos(v[5]),1.54*np.sin(1.911)*np.sin(v[5]),1.54*np.cos(1.911)])
            cadena[4] = np.ndarray.tolist(np.dot(X6.T,np.dot(Y6.T,np.dot(Z6.T,w_t)))+T6)

            w_t = np.array([1.54*np.sin(1.911)*np.cos(v[0]),1.54*np.sin(1.911)*np.sin(v[0]),1.54*np.cos(1.911)])
            cadena[8] = np.ndarray.tolist(np.dot(X0.T,np.dot(Y0.T,np.dot(Z0.T,w_t)))+T0)

            _,_,_,X1,Y1,Z1,T1 = mv.get_spheric([cadena[k] for k in reversed(range(7,11))])
            _,_,_,X5,Y5,Z5,T5 = mv.get_spheric(cadena[2:6])

            w_t = np.array([1.54*np.sin(1.911)*np.cos(v[1]),1.54*np.sin(1.911)*np.sin(v[1]),1.54*np.cos(1.911)])
            cadena[7] = np.ndarray.tolist(np.dot(X1.T,np.dot(Y1.T,np.dot(Z1.T,w_t)))+T1)

            w_t = np.array([1.54*np.sin(1.911)*np.cos(v[4]),1.54*np.sin(1.911)*np.sin(v[4]),1.54*np.cos(1.911)])
            cadena[5] = np.ndarray.tolist(np.dot(X5.T,np.dot(Y5.T,np.dot(Z5.T,w_t)))+T5)

            _,_,_,X4,Y4,Z4,T4 = mv.get_spheric(cadena[3:7])
            w_t = np.array([1.54*np.sin(1.911)*np.cos(v[3]),1.54*np.sin(1.911)*np.sin(v[3]),1.54*np.cos(1.911)])
            cadena[6] = np.ndarray.tolist(np.dot(X4.T,np.dot(Y4.T,np.dot(Z4.T,w_t)))+T4)

            for mnq in range(0,13):
                for uj in range(0,3):
                    if cadena[mnq][uj] > L:
                        cadena[mnq][uj] = cadena[mnq][uj]- L
                    elif cadena[mnq][uj] < 0:
                        cadena[mnq][uj]= cadena[mnq][uj] + L


            all_data_copy = copy.deepcopy(all_data)

            data_2 = copy.deepcopy(data)

            data_2[i:i+13] = copy.deepcopy(cadena)

            all_data_copy[k_i] = copy.deepcopy(data_2)

            potential_energy_2 = potential(all_data_copy,L)

            if potential_energy_2 <= potential_energy_1:
                data = copy.deepcopy(data_2)
                del all_data
                all_data = copy.deepcopy(all_data_copy)
                potential_energy_1 = potential_energy_2

            #del data,cadena,all_data_copy,data_2,potential_energy_2

            np.save("energy.npy",energia_p)
            np.save("data_final.npy",all_data)
