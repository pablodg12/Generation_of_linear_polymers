import numpy as np

#def car2sph(x,y,z):
#    azimuth = np.arctan2(y,x)
#    elevation = np.arctan2(z,np.sqrt(x**2 + y**2))
#    r = np.sqrt(x**2 + y**2 + z**2)
#    return(r,elevation,azimuth)
def cart2sph(x, y, z):
    x_0 = np.array([x,y,z])
    r = np.sqrt(x**2+y**2+z**2)
   
    if z == 0:
        theta = np.pi/2
    if z > 0:
        theta = np.arctan(np.sqrt(x**2 + y**2)/z)
    if z < 0:
        theta = np.pi + np.arctan(np.sqrt(x**2 + y**2)/z)
    if x > 0 and y > 0:
        phi = np.arctan(y/x)
    if x > 0 and y < 0:
        phi = 2*np.pi + np.arctan(y/x)
    if x == 0:
        phi = (np.pi/2)*(np.sign(y))
    if x < 0:
        phi = np.pi + np.arctan(y/x) 
    return(r,theta,phi)

#def cart2sph(x,y,z):
#    radius = np.sqrt(x * x + y * y + z * z)
#    theta = np.arctan2(y, x)
#    phi = np.arccos(z / radius)
#    
#    return(radius,theta,phi)

#def cart2sph(x,y,z):
#    azimuth = np.arctan2(y,x)
#    elevation = np.arctan2(z,np.sqrt(x**2 + y**2))
#    r = np.sqrt(x**2 + y**2 + z**2)
#    return(r,elevation,azimuth)

#Traslation Matrix
def traslation(x_0,t):
    x_0 = np.append(x_0,[1])
    traslation_matrix = np.array([[1,0,0,t[0]],
                                  [0,1,0,t[1]],
                                  [0,0,1,t[2]],
                                  [0,0,0,1]])
    return(np.dot(traslation_matrix,x_0)[:-1])

def rotation_matrix_y(x_0,theta,center=False,g=False):
    if center == True:
        projection = x_0 - np.dot(x_0,np.array([0,1,0]))*np.array([0,1,0])
        rule = np.dot(np.array([0,0,1]),projection)/(np.linalg.norm(projection))
        
        s = np.dot(np.cross(np.array([0,1,0]),x_0),np.array([0,0,1]))/(np.linalg.norm(x_0))
        if s < 0:
            theta = 2*np.pi -np.arccos(rule)
        else:
            theta = np.arccos(rule)
        
    #rotation_matrix = np.array([[1,0,0],
    #                           [0,np.cos(theta),-np.sin(theta)],
    #                              [0,np.sin(theta),np.cos(theta)]])
            
    rotation_matrix = np.array([[np.cos(theta),0,np.sin(theta)],
                              [0,1,0],
                              [-np.sin(theta),0,np.cos(theta)]])
    if g == True:
        return(rotation_matrix,theta)
    return(rotation_matrix)

def rotation_matrix_x(x_0,theta,center=False,g=False):
    if center == True:
        projection = x_0 - np.dot(x_0,np.array([1,0,0]))*np.array([1,0,0])
        rule = np.dot(np.array([0,0,1]),projection)/(np.linalg.norm(projection))
        
        s = np.dot(np.cross(np.array([1,0,0]),x_0),np.array([0,0,1]))/(np.linalg.norm(x_0))
        if s < 0:
            theta = 2*np.pi -np.arccos(rule)
        else:
                theta = np.arccos(rule)
        
    rotation_matrix = np.array([[1,0,0],
                               [0,np.cos(theta),-np.sin(theta)],
                               [0,np.sin(theta),np.cos(theta)]])
    if g == True:
        return(rotation_matrix,theta)
    return(rotation_matrix)
    
def rotation_matrix_z(x_0,theta,center=False,g=False):
    if center == True:
        projection = x_0 - np.dot(x_0,np.array([0,0,1]))*np.array([0,0,1])
        rule = np.dot(np.array([0,1,0]),projection)/(np.linalg.norm(projection))
        
        s = np.dot(np.cross(np.array([0,0,1]),x_0),np.array([0,1,0]))/(np.linalg.norm(x_0))
        if s < 0:
            theta = 2*np.pi -np.arccos(rule)
        else:
            theta = np.arccos(rule)
    rotation_matrix = np.array([[np.cos(theta),-np.sin(theta),0],
                                   [np.sin(theta),np.cos(theta),0],
                                   [0,0,1]])
    if g == True:
        return(rotation_matrix,theta)
    return(rotation_matrix)


def f_correction(point_a,point_b,d):
    dif = np.abs((point_a-point_b)) > np.sqrt((d[1]-d[0])**2)/2
    pg0 = point_b[0]
    pg1 = point_b[1]
    pg2 = point_b[2]
    if np.sum(dif) > 0:
        if dif[0] == True:
            min_distance_p3 = np.sqrt((point_a[0]-d)**2)
            min_distance_p4 = np.sqrt((point_b[0]-d)**2)
            #min_distance_p3 = np.abs(point_a[0]-d)
            #min_distance_p4 = np.abs(point_b[0]-d)
            if np.argmin(min_distance_p3) == 0:
                pg0 = point_a[0] - np.min(min_distance_p3) - np.min(min_distance_p4)
                #pg0 = point_a[0] - point_a[0] - (d[1]-point_b[0])
            else:
                pg0 = point_a[0] + np.min(min_distance_p3) + np.min(min_distance_p4)
                
        if dif[1] == True:
            min_distance_p3 = np.sqrt((point_a[1]-d)**2)
            min_distance_p4 = np.sqrt((point_b[1]-d)**2)
            #min_distance_p3 = np.abs(point_a[1]-d)
            #min_distance_p4 = np.abs(point_b[1]-d)
            if np.argmin(min_distance_p3) == 0:
                pg1 = point_a[1] - np.min(min_distance_p3) - np.min(min_distance_p4)
            else:
                pg1 = point_a[1] + np.min(min_distance_p3) + np.min(min_distance_p4)
        if dif[2] == True:
            min_distance_p3 = np.sqrt((point_a[2]-d)**2)
            min_distance_p4 = np.sqrt((point_b[2]-d)**2)
            #min_distance_p3 = np.abs(point_a[2]-d)
            #min_distance_p4 = np.abs(point_b[2]-d)
            if np.argmin(min_distance_p3) == 0:
                pg2 = point_a[2] - np.min(min_distance_p3) - np.min(min_distance_p4)
            else:
                pg2 = point_a[2] + np.min(min_distance_p3) + np.min(min_distance_p4)
    return(pg0,pg1,pg2)

def b2_correction(point_a,point_b,point_c,d):
    dif = np.abs((point_a-point_b)) > np.sqrt((d[1]-d[0])**2)/2
    dif_2 = point_c-point_b
    pg0=point_b[0]
    pg1=point_b[1]
    pg2=point_b[2]
    if np.sum(dif) > 0:
        if dif[0] == True:
            min_distance_p3 = np.sqrt((point_a[0]-d)**2)
            min_distance_p4 = np.sqrt((point_b[0]-d)**2)
            #min_distance_p3 = np.abs(point_a[0]-d)
            #min_distance_p4 = np.abs(point_b[0]-d)
            if np.argmin(min_distance_p3) == 0:
                pg0 = point_a[0] - np.min(min_distance_p3) - np.min(min_distance_p4)
            else:
                pg0 = point_a[0]+ np.min(min_distance_p3) + np.min(min_distance_p4)
        if dif[1] == True:
            min_distance_p3 = np.sqrt((point_a[1]-d)**2)
            min_distance_p4 = np.sqrt((point_b[1]-d)**2)
            #min_distance_p3 = np.abs(point_a[1]-d)
            #min_distance_p4 = np.abs(point_b[1]-d)
            if np.argmin(min_distance_p3) == 0:
                pg1 = point_a[1] - np.min(min_distance_p3) - np.min(min_distance_p4)
            else:
                pg1 = point_a[1]+ np.min(min_distance_p3) + np.min(min_distance_p4)
        if dif[2] == True:
            min_distance_p3 = np.sqrt((point_a[2]-d)**2)
            min_distance_p4 = np.sqrt((point_b[2]-d)**2)
            #min_distance_p3 = np.abs(point_a[2]-d)
            #min_distance_p4 = np.abs(point_b[2]-d)
            if np.argmin(min_distance_p3) == 0:
                pg2 = point_a[2] - np.min(min_distance_p3) - np.min(min_distance_p4)
            else:
                pg2 = point_a[2]+ np.min(min_distance_p3) + np.min(min_distance_p4)
        return(pg0,pg1,pg2,pg0+dif_2[0],pg1+dif_2[1],pg2+dif_2[2])
    return(pg0,pg1,pg2,point_c[0],point_c[1],point_c[2])

def get_spheric(cadena):
    mini_chain = []
    mini_chain.append(traslation(np.array(cadena[0]),-1*np.array(cadena[2])))
    mini_chain.append(traslation(np.array(cadena[1]),-1*np.array(cadena[2])))
    mini_chain.append(traslation(np.array(cadena[2]),-1*np.array(cadena[2])))
    mini_chain.append(traslation(np.array(cadena[3]),-1*np.array(cadena[2])))
    
    rot_matrix_x = rotation_matrix_x(mini_chain[1],0,center=True)
    
    mini_chain[0] = np.dot(rot_matrix_x,mini_chain[0])
    mini_chain[1] = np.dot(rot_matrix_x,mini_chain[1])
    mini_chain[2] = np.dot(rot_matrix_x,mini_chain[2])
    mini_chain[3] = np.dot(rot_matrix_x,mini_chain[3])
    
    rot_matrix_y = rotation_matrix_y(mini_chain[1],0,center=True)
    
    mini_chain[0] = np.dot(rot_matrix_y,mini_chain[0])
    mini_chain[1] = np.dot(rot_matrix_y,mini_chain[1])
    mini_chain[2] = np.dot(rot_matrix_y,mini_chain[2])
    mini_chain[3] = np.dot(rot_matrix_y,mini_chain[3])
    
    rot_matrix_z = rotation_matrix_z(mini_chain[0],0,center=True)
    
    mini_chain[0] = np.dot(rot_matrix_z,mini_chain[0])
    mini_chain[1] = np.dot(rot_matrix_z,mini_chain[1])
    mini_chain[2] = np.dot(rot_matrix_z,mini_chain[2])
    mini_chain[3] = np.dot(rot_matrix_z,mini_chain[3])
    
    r,theta,phi = cart2sph(mini_chain[3][0], mini_chain[3][1],  mini_chain[3][2])

    return(r,theta,phi,rot_matrix_x,rot_matrix_y,rot_matrix_z,cadena[2])

def get_spheric_2(cadena):
    mini_chain = []
    mini_chain.append(traslation(np.array(cadena[0]),-1*np.array(cadena[2])))
    mini_chain.append(traslation(np.array(cadena[1]),-1*np.array(cadena[2])))
    mini_chain.append(traslation(np.array(cadena[2]),-1*np.array(cadena[2])))
    mini_chain.append(traslation(np.array(cadena[3]),-1*np.array(cadena[2])))
    
    rot_matrix_x,a = rotation_matrix_x(mini_chain[1],0,center=True,g=True)
    
    mini_chain[0] = np.dot(rot_matrix_x,mini_chain[0])
    mini_chain[1] = np.dot(rot_matrix_x,mini_chain[1])
    mini_chain[2] = np.dot(rot_matrix_x,mini_chain[2])
    mini_chain[3] = np.dot(rot_matrix_x,mini_chain[3])
    
    rot_matrix_y,b = rotation_matrix_y(mini_chain[1],0,center=True,g=True)
    
    mini_chain[0] = np.dot(rot_matrix_y,mini_chain[0])
    mini_chain[1] = np.dot(rot_matrix_y,mini_chain[1])
    mini_chain[2] = np.dot(rot_matrix_y,mini_chain[2])
    mini_chain[3] = np.dot(rot_matrix_y,mini_chain[3])
    
    rot_matrix_z,c = rotation_matrix_z(mini_chain[0],0,center=True,g=True)
    
    mini_chain[0] = np.dot(rot_matrix_z,mini_chain[0])
    mini_chain[1] = np.dot(rot_matrix_z,mini_chain[1])
    mini_chain[2] = np.dot(rot_matrix_z,mini_chain[2])
    mini_chain[3] = np.dot(rot_matrix_z,mini_chain[3])
    
    r,theta,phi = cart2sph(mini_chain[3][0], mini_chain[3][1],  mini_chain[3][2])

    return(r,theta,phi,rot_matrix_x,rot_matrix_y,rot_matrix_z,cadena[2],a,b,c)

