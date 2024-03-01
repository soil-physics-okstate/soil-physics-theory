# assimilation.py
# ---------------------------------------------------------
# This module is part of the CRITERIA3D_LAB distribution
# https://github.com/ARPA-SIMC/CRITERIA3D_LAB
# ---------------------------------------------------------

from dataStructures import *
import rectangularMesh
import criteria3D
import soil
import numpy as np
from scipy.interpolate import interpn


def buildDataStructuresForInterpolation(initialState):
    # Select just the vertices that compose a hyper-plane
    initialState = initialState.sort_values(by=['x', 'y', 'z'])
    # Understand which interpolation to apply (1D, 2D, 3D)
    x = initialState["x"].unique()
    y = initialState["y"].unique()
    z = initialState["z"].unique()
    # Instantiate and populate the data structure where all the known coordinates are stored
    points = []
    if len(x) > 1:
        points += [x]
    if len(y) > 1:
        points += [y]
    if len(z) > 1:
        points += [z]
    # Instantiate the data structure where all the values of the known points are stored
    if len(points) == 1:
        values = np.zeros(len(points[0]))
    if len(points) == 2:
        values = np.zeros((len(points[0]), len(points[1])))
    if len(points) == 3:
        values = np.zeros((len(points[0]), len(points[1]), len(points[2])))
    # Populate the data structure where all the values of the known points are stored

    for i in range(x.shape[-1]):
        for j in range(y.shape[-1]):
            for k in range(z.shape[-1]):
                psi = initialState[(initialState["x"] == x[i]) & (initialState["y"] == y[j])
                                   & (initialState["z"] == z[k])]["value"]
                # from kPa to meters
                observedPsi = float(psi / 9.81)
                curve = C3DParameters.waterRetentionCurve
                horizonIndex = soil.getHorizonIndex(z[k])
                theta = soil.thetaFromPsi(curve, observedPsi, soil.horizons[horizonIndex])
                index = rectangularMesh.getCellIndex(x[i], y[j], z[k])
                if index != NODATA:
                    currentTheta = soil.getVolumetricWaterContent(index)
                    # check validity range of sensors
                    if (abs(theta) > 0.3) and (abs(currentTheta) > 0.3):
                        value = 0
                    else:
                        value = theta - currentTheta
                else:
                    value = 0

                if len(points) == 1:
                    # Even though we do not know which is the coordinate that has not one unique value,
                    # we can sum all of the coordinates because the index of the ones that has just one unique value would be 0
                    values[i + j + k] = value
                if len(points) == 2:
                    # We understand which is the pair of coordinates to consider, based on the one to not consider (lenght = 1)
                    if len(x) == 1:
                        values[j, k] = value
                    if len(y) == 1:
                        values[i, k] = value
                    if len(z) == 1:
                        values[i, j] = value
                if len(points) == 3:
                    values[i, j, k] = value
    return points, values, [x, y, z]


def getVertices(domain):
    mins = [min(axis) for axis in domain]
    maxs = [max(axis) for axis in domain]
    origin = (mins[0], mins[1], mins[2])
    maxFirstCoordinate = (maxs[0], mins[1], mins[2])
    maxSecondCoordinate = (mins[0], maxs[1], mins[2])
    maxThirdCoordinate = (mins[0], mins[1], maxs[2])
    return origin, maxFirstCoordinate, maxSecondCoordinate, maxThirdCoordinate


def interpolate(initialState):
    points, values, domain = buildDataStructuresForInterpolation(initialState)
    origin, maxFirstCoordinate, maxSecondCoordinate, maxThirdCoordinate = getVertices(domain)

    originIndex = rectangularMesh.getCellIndex(origin[0], origin[1], origin[2])
    maxFirstCoordinateIndex = rectangularMesh.getCellIndex(maxFirstCoordinate[0], maxFirstCoordinate[1],
                                                           maxFirstCoordinate[2])
    maxSecondCoordinateIndex = rectangularMesh.getCellIndex(maxSecondCoordinate[0], maxSecondCoordinate[1],
                                                            maxSecondCoordinate[2])
    maxThirdCoordinateIndex = rectangularMesh.getCellIndex(maxThirdCoordinate[0], maxThirdCoordinate[1],
                                                           maxThirdCoordinate[2])

    interpolated_points = {}
    for x_inc in range(0, maxFirstCoordinateIndex - originIndex + 1, 1):
        for y_inc in range(0, maxSecondCoordinateIndex - originIndex + C3DStructure.nrRectanglesInXAxis,
                           C3DStructure.nrRectanglesInXAxis):
            for z_inc in range(0, maxThirdCoordinateIndex - originIndex + C3DStructure.nrRectangles,
                               C3DStructure.nrRectangles):

                index = originIndex + x_inc + y_inc + z_inc
                coords = [round(num, 2) for num in rectangularMesh.getXYDepth(index)]
                x, y, z = coords[0], coords[1], coords[2]

                if len(points) == 1:
                    point = np.array([x + y + z])
                if len(points) == 2:
                    if len(domain[0]) == 1:
                        point = np.array([y, z])
                    if len(domain[1]) == 1:
                        point = np.array([x, z])
                    if len(domain[2]) == 1:
                        point = np.array([x, y])
                if len(points) == 3:
                    point = np.array([x, y, z])

                value = interpn(points, values, point, bounds_error=False)
                interpolated_points[index] = value[0]
    return interpolated_points


def getCloserIndex(i, indices):
    min_index = float('inf')
    min_distance = float('inf')
    # x_i, y_i, z_i = rectangularMesh.getXYDepth(i)
    for index in indices:
        # x_index, y_index, z_index = rectangularMesh.getXYDepth(index)
        distance = criteria3D.getCellDistance(i, index)
        if distance < min_distance:
            min_index = index
            min_distance = distance
    return min_index


def assimilate(initialState):
    interpolated_points = interpolate(initialState)

    if not ((initialState['x'] < 0).any()):
        initialState[['x', 'y']] = -initialState[['x', 'y']]
        new_interpolated_points = interpolate(initialState)
        interpolated_points = {**interpolated_points, **new_interpolated_points}

    indices = interpolated_points.keys()
    for i in range(C3DStructure.nrCells):
        x, y, depth = rectangularMesh.getXYDepth(i)
        if depth != 0:
            if i in indices:
                value = interpolated_points[i]
            else:
                min_index = getCloserIndex(i, indices)
                value = interpolated_points[min_index]
            # assign residual of volumetric water content
            currentTheta = soil.getVolumetricWaterContent(i)
            theta = currentTheta + value
            curve = C3DParameters.waterRetentionCurve
            horizonIndex = soil.getHorizonIndex(depth)
            psi = soil.psiFromTheta(curve, theta, soil.horizons[horizonIndex])
            criteria3D.setMatricPotential(i, psi)
