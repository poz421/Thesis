#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <tuple>
#include <unordered_map>
#include "base.h"
#include "distanceUtils.h"
#include <cmath>

class TrajGrid {

    /*
    A class representing a simple grid to index trajectory points.
    
    Atrributes:
        data { (int,int): {int: [int]} }: dictionary mapping a grid cell
                to a dictionary containing traj points lying in the cell.
        startLat (float): smallest y-value of indexed points.
        startLon (float): smallest x-value of indexed points.
        delta (float): size of grid cell.
        numX (int): number of cells in the x-direction.
        numY (int): number of cells in the y-direction.
    */
public:
    std::map<std::pair<int, int>, std::unordered_map<int, std::vector<int>>> data;
    double startLat;
    double startLon;
    double delta;
    int numX;
    int numY;
    //Simple initialization of class. 
    TrajGrid(std::map<std::pair<int, int>, std::unordered_map<int, std::vector<int>>> data, double startLat, double startLon, double delta, int numX, int numY);
};

/*
Create an index of trajectory points.

    Args:
        trajs ({int : traj}): dictionary of traj objects.
        simpTajs ([simpleTraj]): list of simplified trajectories whose points
                are to be indexed.
        delta (float): size of grid cell.

    Returns:
        trajGrid object containing the points feom the simplified trajectories.
*/


TrajGrid gridData(std::unordered_map<int, Traj> trajs, std::vector<SimpleTraj> simpTrajs, double delta);
       
