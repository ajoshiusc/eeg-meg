'''||AUM||'''
import numpy as np
from scipy.sparse.linalg import cg  # ,spilu,LinearOperator
from scipy.sparse import eye, spdiags, csc_matrix, vstack  # , diags
from vtk import (vtkSphereSource, vtkPolyData, vtkDecimatePro, vtkPoints, vtkCleanPolyData, vtkCellArray,
                 vtkPolyDataConnectivityFilter, vtkSmoothPolyDataFilter, vtkPolyDataNormals, vtkCurvatures)
from dfsio import readdfs, writedfs
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray, vtk_to_numpy
import scipy as sp
from mayavi import mlab
import trimesh as tm

__author__ = "Anand A. Joshi"
__copyright__ = "University of Southern California, Los Angeles"
__email__ = "ajoshi@sipi.usc.edu"


def face_areas(s):
    r12 = s.vertices[s.faces[:, 0],:]
    r13 = s.vertices[s.faces[:, 2],:] - r12
    r12 = s.vertices[s.faces[:, 1],:] - r12
    area1 = sp.sqrt(sp.sum(sp.cross(r12,r13)** 2, axis=1)) / 2.0
    return area1


def smooth_surf_function(s, f0, a1=None, a2=None, aniso=None, normalize=None):
    if a1 is None:
        a1 = 3.1

    if a2 is None:
        a2 = 3.1

    if aniso is None:
        aniso = np.ones((len(s.vertices), 1))

    if normalize is None:
        normalize = 0

    S, Dx, Dy = get_stiffness_matrix_tri_wt(s, aniso)
    M = vstack((a1*S, a2*vstack((Dx, Dy))))
    A = vstack((eye(len(s.vertices)), M))
    b = np.concatenate((f0, np.zeros(M.shape[0])))
    AtA = A.T*A

    M = spdiags(1.0/AtA.diagonal(), [0], AtA.shape[0], AtA.shape[1])
#    M=np.diag(np.diag(A.transpose()*A))
    f, nits = cg(AtA, A.T*b, tol=1e-100, maxiter=600, M=M)
#    print('f',f)
    if normalize > 0:
        f = f*np.linalg.norm(f0)/np.linalg.norm(f)

    return f


def get_stiffness_matrix_tri_wt(surf1, W):
    W = np.squeeze(W)
    X = surf1.vertices[:, 0]
    Y = surf1.vertices[:, 1]
    Z = surf1.vertices[:, 2]
    NumTri = surf1.faces.shape[0]
#    NumVertx = X.shape[0]
    vertx_1 = surf1.faces[:, 0]
    vertx_2 = surf1.faces[:, 1]
    vertx_3 = surf1.faces[:, 2]
    V1 = np.column_stack((X[vertx_1], Y[vertx_1], Z[vertx_1]))
    V2 = np.column_stack((X[vertx_2], Y[vertx_2], Z[vertx_2]))
    V3 = np.column_stack((X[vertx_3], Y[vertx_3], Z[vertx_3]))
    x1 = np.zeros((NumTri))
    y1 = np.zeros((NumTri))
    v2_v1temp = V2-V1
    x2 = np.linalg.norm(v2_v1temp, axis=1)
    y2 = np.zeros((NumTri))
    x3 = np.einsum('ij,ij->i', (V3-V1),
                   (v2_v1temp/np.column_stack((x2, x2, x2))))
    mynorm = np.cross((v2_v1temp), V3-V1, axis=1)
    yunit = np.cross(mynorm, v2_v1temp, axis=1)
    y3 = np.einsum('ij,ij->i', yunit, (V3-V1))/np.linalg.norm(yunit, axis=1)
    sqrt_DT = (np.abs((x1*y2 - y1*x2)+(x2*y3 - y2*x3)+(x3*y1 - y3*x1)))
    Ar = 0.5*(np.abs((x1*y2 - y1*x2)+(x2*y3 - y2*x3)+(x3*y1 - y3*x1)))
#    Ar=1.0+0.0*Ar
    y1 = y1/sqrt_DT
    y2 = y2/sqrt_DT
    y3 = y3/sqrt_DT
    x1 = x1/sqrt_DT
    x2 = x2/sqrt_DT
    x3 = x3/sqrt_DT
    tmp_A = np.concatenate((y2-y3, y3-y1, y1-y2), axis=0)
    tmp_B = np.concatenate((x3-x2, x1-x3, x2-x1), axis=0)

    rowno = np.arange(0, NumTri)
    rowno_all = np.concatenate((rowno, rowno, rowno))
    vertx_all = np.concatenate((vertx_1, vertx_2, vertx_3))
    Dx = csc_matrix((tmp_A, (rowno_all, vertx_all)), (NumTri, len(X)))
    Dy = csc_matrix((tmp_B, (rowno_all, vertx_all)), (NumTri, len(X)))

    TC = face_v_conn(surf1)
    Wt = (1.0/3.0)*(TC.T*W)
    Wt = spdiags(Wt*Ar, (0), NumTri, NumTri)
    S = Dx.T*Wt*Dx + Dy.T*Wt*Dy

    return S, Dx, Dy


def face_v_conn(FV):

    nFaces = FV.faces.shape[0]
#    nVertices = FV.vertices.shape[0]

    rows = np.concatenate((FV.faces[:, 0], FV.faces[:, 1], FV.faces[:, 2]))
    cols = np.concatenate((np.arange(0, nFaces), np.arange(0, nFaces),
                           np.arange(0, nFaces)))
    data = np.ones((len(rows)))
    TriConn = csc_matrix((data, (rows, cols)))

    return TriConn


def mkVtkIdList(it):
    vil = vtk.vtkIdList()
    for i in it:
        vil.InsertNextId(int(i))
    return vil


def createPolyData(verts, faces):

    poly = vtkPolyData()

    points = vtkPoints()
    for i in range(0,verts.shape[0]):
        points.InsertPoint(i, verts[i,])
#    points.SetData(numpy_to_vtk(verts))

    poly.SetPoints(points)

    tri = vtkCellArray()
    for i in range(0, faces.shape[0]):
        tri.InsertNextCell(3)
        tri.InsertCellPoint(faces[i, 0])
        tri.InsertCellPoint(faces[i, 1])
        tri.InsertCellPoint(faces[i, 2])

    poly.SetPolys(tri)

    return poly


def reducepatch(surf, ratio=0.90, VERBOSITY=0):
    ratio = 1.0 - ratio
    pol = createPolyData(surf.vertices, surf.faces)

    decimate = vtkDecimatePro()
    decimate.SetInputConnection(pol.GetProducerPort())
    decimate.SetTargetReduction(ratio)
    decimate.Update()

    decimatedPoly = vtkPolyData()
    decimatedPoly.ShallowCopy(decimate.GetOutput())


    connectivity2 = vtkPolyDataConnectivityFilter()
    connectivity2.SetInput(decimatedPoly)
    connectivity2.SetExtractionModeToLargestRegion()
    connectivity2.Update()
    decimatedPoly = connectivity2.GetOutput()

    cleaner = vtkCleanPolyData()
    cleaner.SetInput(decimatedPoly)
    cleaner.Update()
    decimatedPoly = cleaner.GetOutput()

    if VERBOSITY >0 :
        print("\n # extracted regions: " +
            str(connectivity2.GetNumberOfExtractedRegions()))
        print("Before decimation: " + str(pol.GetNumberOfPoints()) + " points;" +
              str(pol.GetNumberOfPolys()) + " polygons.")
        print("After decimation: " + str(decimatedPoly.GetNumberOfPoints()) + " points;"+
              str(decimatedPoly.GetNumberOfPolys()) + " polygons .")

    pts = decimatedPoly.GetPoints()
    surf.vertices = vtk_to_numpy(pts.GetData())
    faces1 = decimatedPoly.GetPolys()
    f1 = faces1.GetData()
    f2 = vtk_to_numpy(f1)
    f2 = f2.reshape(len(f2)/4,4)
    surf.faces = f2[:, 1:]
     
    return surf

def add_normals(s1):
    normals = vtkPolyDataNormals()
    normals.SetInput(createPolyData(s1.vertices, s1.faces))
    normals.ComputePointNormalsOn()
    normals.ComputeCellNormalsOn()
    normals.SplittingOff()
    normals.Update()
    normals.AutoOrientNormalsOn()
    normals.ConsistencyOn()
    # normals.Update()
    n1 = normals.GetOutput()
    p1 = n1.GetPointData()
    s1.normals = vtk_to_numpy(p1.GetNormals())
    s1.vertices = vtk_to_numpy(n1.GetPoints().GetData())
    faces = vtk_to_numpy(n1.GetPolys().GetData())
    faces = faces.reshape(len(faces)/4,4)
    s1.faces = faces[:, 1:]

    return s1
def mean_curvature(s):
    curve1 = vtkCurvatures()
    s=createPolyData(s.vertices, s.faces)
    curve1.SetInput(s)
    curve1.SetCurvatureTypeToMean()
    curve1.Update()
    m=curve1.GetOutput()
    m=vtk_to_numpy(m.GetPointData().GetScalars())
    return m

def get_sphere(center= [0,0,0], radius= 5.0, res=100):
    source = vtkSphereSource()
    source.SetCenter(center[0], center[1], center[2])
    source.SetRadius(radius)
    source.SetThetaResolution(res)
    source.SetPhiResolution(res)
    source.Update()
    surf1= source.GetOutput()
    pts = surf1.GetPoints()
    vert1 = vtk_to_numpy(pts.GetData())
    faces1 = surf1.GetPolys()
    f1 = faces1.GetData()
    f2 = vtk_to_numpy(f1)
    f2 = f2.reshape(len(f2) / 4, 4)
    class surf:
        pass
    surf.faces = f2[:, 1:]
    surf.vertices = vert1
    return surf

#

def view_patch(r,attrib=[]):
    mlab.figure()
    if len(attrib)>0:
        mlab.triangular_mesh(r.vertices[:, 0], r.vertices[:, 1], r.vertices[:, 2],
                             r.faces, representation='surface', opacity=1, scalars=attrib)
    else:
        mlab.triangular_mesh(r.vertices[:, 0], r.vertices[:, 1], r.vertices[:, 2],
                             r.faces, representation='surface', opacity=1)

    mlab.gcf().scene.parallel_projection = True
    mlab.view(azimuth=0, elevation=90)
    mlab.colorbar(orientation='horizontal')

    mlab.draw()
    mlab.show()


def readdfsVTK(fname):
    s=readdfs(fname)
    poly = createPolyData(s.vertices, s.faces)
    return poly


def smooth_patch(surf, iter1=15, relax1=0.1):
    smoothFilter = vtkSmoothPolyDataFilter()
    smoothFilter.SetInput(createPolyData(surf.vertices, surf.faces))
    smoothFilter.SetNumberOfIterations(iter1)
    smoothFilter.SetRelaxationFactor(relax1)
    smoothFilter.Update()
    surf1= smoothFilter.GetOutput()
    pts = surf1.GetPoints()
    vert1 = vtk_to_numpy(pts.GetData())
    faces1 = surf1.GetPolys()
    f1 = faces1.GetData()
    f2 = vtk_to_numpy(f1)
    f2 = f2.reshape(len(f2) / 4, 4)
    surf.faces = f2[:, 1:]
    surf.vertices = vert1
    return surf
