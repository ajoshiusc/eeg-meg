
# coding: utf-8

# In[13]:

from dfsio import writedfs
from surfproc import (mean_curvature, view_patch, smooth_surf_function,
                      get_sphere, add_normals, face_v_conn, face_areas)
from numpy.linalg import norm as norm2
import scipy as sp
brain_name = 'C:\\Users\\ajoshi\\Downloads\\C113AV_scan1\\\
c7452AD_AV_ave.brain.dfs'
skull_name = 'C:\\Users\\ajoshi\\Downloads\\C113AV_scan1\\\
c7452AD_AV_ave.outer_skull.dfs'
brain = get_sphere([0, 0, 0], radius=.06000, res=150)
skull = get_sphere([0, 0, 0], radius=.12000, res=50)

# brain = readdfs(brain_name)
# brain=smooth_patch(brain,100)
# brain = reducepatch(brain,.9,0)
# brain=smooth_patch(brain,10)

# skull = readdfs(skull_name)
# skull = reducepatch(skull,.9,0)
brain = add_normals(brain)
# brain.normals=brain.vertices/5.0
skull = add_normals(skull)
# view_patch(brain)
# view_patch(skull)
writedfs('brain.dfs', brain)
writedfs('skull.dfs', skull)
print "# vertices in Brain: " + str(brain.vertices.shape[0])
print "# vertices in Skull: " + str(skull.vertices.shape[0])
m = mean_curvature(brain)
# m = smooth_surf_function(brain, m, 10, 100)
br = view_patch(brain, attrib=m, opacity=1, show=0)
view_patch(skull, opacity=.1, fig=br, show=1)
Tri = face_v_conn(brain)

Q = 1.0+sp.zeros((brain.faces.shape[0], 1))
Q = (1.0/3.0)*Tri*Q
Q = Q[:, 0]
br_face_area = face_areas(brain)
Q = 1.0+sp.zeros(brain.vertices.shape[0])
# Q=0.0*Q
# Q[1800]=1
# Q=smooth_surf_function(brain,Q,10,10)
area_v = (1.0/3.0)*Tri*br_face_area
v = sp.zeros(skull.vertices.shape[0])  # eq 8 from Sarvas
v_aaj = sp.zeros(skull.vertices.shape[0])  # joshi
view_patch(brain, attrib=Q, opacity=1)
for i in range(skull.vertices.shape[0]):
    r0 = brain.vertices
    r_r0 = skull.vertices[i, ]-r0
    norm_r_r0 = norm2(r_r0, ord=2, axis=1)
    den = norm_r_r0 ** 3.0
    num = r_r0
    v1 = num/sp.expand_dims(den, axis=1)
    Qvec = sp.expand_dims(Q, axis=1)*brain.normals
    v1_aaj = (2.0*m)/norm_r_r0
    v1 = sp.sum(v1 * Qvec, axis=1)
#    v1_aaj=v1_aaj*sp.sign(v1)
    v[i] = sp.sum(v1*area_v)
    v_aaj[i] = sp.sum(Q*v1_aaj*area_v)
    if sp.mod(i, 100) == 0:
        print 'i='+str(i)


# In[15]:

view_patch(skull, attrib=v,show=1)
#view_patch(skull, attrib=v_aaj,show=1)
