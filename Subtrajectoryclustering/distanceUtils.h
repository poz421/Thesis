#pragma once
#include <vector>
#include "base.h"
#include <algorithm>

/*

Class representing a 'simplified' (sub)trajectory.

    Instead of storing the actual points, stores the indices of points left
    over from the parent (sub)trajectory after simplification.

    Attributes:
        trajID: ID of parent (sub)trajectory.
        indices: indices of points left over from the parent trajectory. Indices
                 are based on the list of points sorted by timestamp.
*/
class SimpleTraj {
public:

    int trajID;
    std::vector<int> indices;

    // Default constructor
    SimpleTraj();

    // Parameterized constructor
    SimpleTraj(int trajID, const std::vector<int>& indices);


    
};

/*
Simplify a (sub)trajectory in the forward direction.

    Traverse the points in ascending order of timestamps, only retaining a point if
    it is more than a certain distance away from the last retained point. The first
    and last points are always retained.

    Args:
        uTraj (traj): trajectory to be simplified (points are sorted by timestamp).
        tau (float): distance for simplification.
        start (int): starting index of trajectory point to be simplified.
        end (int): ending index of trajectory point to be simplified.

    Returns:
        A forward-simplified version of the subtrajectory lying between start and end.

*/
SimpleTraj fSimplify(const Traj& uTraj, double tau, int start, int end);

/*
Simplify a (sub)trajectory in the reverse direction.

    Traverse the points in descending order of timestamps, only retaining a point if
    it is more than a certain distance away from the last retained point. The first
    and last points are always retained.

    Args:
        uTraj (traj): trajectory to be simplified (points are sorted by timestamp).
        tau (float): distance for simplification.
        start (int): starting index of trajectory point to be simplified.
        end (int): ending index of trajectory point to be simplified.

    Returns:
        A backward-simplified version of the subtrajectory lying between start and end.

*/

SimpleTraj bSimplify(const Traj& uTraj, double tau, int start, int end);


/*

Return canonical intervals of the input interval.

    Args:
        i (int): left endpoint of input interval (inclusive).
        j (int): right endpoint of input interval (inclusive).

    Returns:
        List of canonical intervals. Each interval is of the form (start,end)
        where start and end are integers.
*/

std::vector<std::pair<int, int>> canonise(int i, int j);