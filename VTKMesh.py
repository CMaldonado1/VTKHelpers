import vtk
import numpy as np
import os
import meshio # tested with 2.3.0
from pyevtk.hl import pointsToVTK

'''
This module is aimed to simplify the implementation of common tasks on VTK meshes,
that result overly convoluted if the usual VTK Python wrapper for C++ is used,
and render the code difficult to follow.
'''

class VTKObject():

    def __init__(self, filename=None):

        self.filename = filename
        self.n_points = None
        self.n_cells = None
        self.points = None
        self.edges = []
        self.triangles = []
        self.neighbors_dict = {}

        if filename is not None:
            self.reader = vtk.vtkPolyDataReader()
            self.reader.SetFileName(self.filename)
            self.reader.Update()
            self.getMesh()
            self.extractPointsAndEdges()
            self.generateNeighborsDict()
            self.extractTriangles()
            self.extractPartitionIDs()


    def getMesh(self):

        '''
        :return: numpy array where each element is a triple of (x, y, z)
        coordinates and a set of indices representing the links to that point
        '''

        output = self.reader.GetOutput()
        self.n_points = output.GetNumberOfPoints()
        self.n_cells = output.GetNumberOfCells()
        self.points = np.array([output.GetPoint(i) for i in range(self.n_points)])
        # return output


    def extractPointsAndEdges(self):

        output = self.reader.GetOutput()
        for i in range(self.n_cells):
            pts_cell = output.GetCell(i)
            for j in range(pts_cell.GetNumberOfEdges()):
                # print("%s %s" % (pts_cell.GetEdge(j).GetPointId(0), pts_cell.GetEdge(j).GetPointId(1)))
                self.edges.append(([int(pts_cell.GetEdge(j).GetPointId(i)) for i in (0, 1)]))


    def generateNeighborsDict(self):
        for nb_pair in self.edges:
            self.neighbors_dict[nb_pair[0]] = self.neighbors_dict.get(nb_pair[0], []) + [nb_pair[1]]


    def extractTriangles(self):
        output = self.reader.GetOutput()
        # self.triangles = tuple([[int(output.GetCell(j).GetPointId(i)] for i in (0, 1, 2)) for j in range(self.n_cells)])
        self.triangles = [[int(output.GetCell(j).GetPointId(i)) for i in (0, 1, 2)] for j in range(self.n_cells)]


    def extractPartitionIDs(self):

        '''
        Extracting Subpart ID
        '''

        output = self.reader.GetOutput()

        pp = output.GetPointData().GetArray(0)
        self.subpartID = [int(pp.GetComponent(i, 0)) for i in range(self.n_points)]


    def extractSubpart(self, ids):

        subvtk = VTKObject()
        subvtk.points = np.array([self.points[i] for i, x in enumerate(self.points) if self.subpartID[i] in ids])
        subvtk.point_ids = { i for i, id in enumerate(self.subpartID) if id in ids }
        subvtk.n_points = len(subvtk.points)

        #
        subvtk.triangles = [tuple(triang) for triang in self.triangles if all([kk in subvtk.point_ids for kk in triang])]

        # : self.neighbors_dict[i] for i, _ in enumerate(self.subpartID) if self.subpartID[i] in ids}

        return subvtk


    # point cloud to vtk
    def save_utk(self, filename):
        # for frame in range(data.shape[0]):
        # x, y, z = data[frame][:, 0], data[frame][:, 1], data[frame][:, 2]
        x, y, z = (self.points[:,i] for i in (0, 1, 2))

        # It doesn't work sending directly the coordinates
        index = np.random.choice(x.shape[0], size=x.shape[0], replace=False)
        x, y, z = x[index], y[index], z[index]

        img_dir = os.path.dirname(filename)
        if len(img_dir) != 0 and not os.path.exists(img_dir):
             os.mkdir(img_dir)

        # saving the utk files
        pointsToVTK(filename, x, y, z, data=None)


    # mesh to vtk
    def SaveMeshToVTK(self, filename):
        meshio.write_points_cells(
            filename,
            self.points,
            {'triangle': np.array(self.triangles)}
        )
        #mesh = meshio.Mesh(self.points, )
        #meshio.write(filename, mesh)


    # def save(self, filename):
    #
    #     #for frame in range(self.points.shape[0]):
    #         # x, y, z = data[frame][:, 0], data[frame][:, 1], data[frame][:, 2]
    #         # x, y, z = (self.points[frame][:, i] for i in (0, 1, 2))
    #
    #
    #         # It doesn't work sending directly the coordinates
    #         # index = np.random.choice(x.shape[0], size=x.shape[0], replace=False)
    #         # x, y, z = x[index], y[index], z[index]
    #
    #     x = self.points[:, 0]
    #     y = self.points[:, 1]
    #     z = self.points[:, 2]
    #
    #     # for i in (0, 1, 2))
    #     # print(x)
    #
    #     odir = os.path.dirname(filename)
    #     if not len(odir) == 0 and not os.path.exists(odir):
    #         os.mkdir(odir)
    #
    #     # saving the utk files
    #     pointsToVTK(filename, x, y, z, data=None)