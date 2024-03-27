#include "base.h"
#include "grid.h"
#include "frechet.h"
#include "distanceUtils.h"
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>

/*

Check if a point lies inside an axis-aligned square.

    Args:
        center (pt): center of the square.
        point (pt): point to be checked.
        r (float): half side lengtg of square.

    Returns:
        True if point lies inside the square, False otherwise.

*/
bool inSquare(const Pt& center, const Pt& point, double r);


/*

Find points in grid that are 'close' to the input point.

    Compute all points in the (at most 9) neighboring cells of the grid cell
    containing the input point (inluding the input point's cell).

    Args:
        grid (trajGrid): grid containing trajectory points.
        point (pt): input point.

    Returns:
        Dictionay of the form {int : [int]}, mapping trajID to the indices of
        its points that lie in the neighboring cells of input point's cell.
*/
//tested in c++ but not in python
std::unordered_map<int, std::vector<int>> findNeighboringPts(const TrajGrid& grid, const Pt& point);

/*
Compute all subtrajectories, with both endpoints coming from a given set of
        points, that have small Frechet distance to a pathlet.

        Args:
            trajs ({int : traj}): dictionary mapping IDs to trajectories.
            simpTraj (simpleTraj): simplified trajectory whose subtrajectories are
                    to be checked for proximity to input pathlet.
            pth (simpleTraj): input simplified pathlet.
            r (float): guess on the Frechet distance.
            possStarts ([int]): indices of possible starting positions for subtrajectories
                    of simpTraj.
            possEnds ([int]): indices of possible ending positions for subtrajectories
                    of simpTraj.
            upth (pathlet): original unsimplified pathlet.
            pathDic ({(pathlet, subTraj) : float}): dictionary mapping pathlet-subtrajectory
                    pair to the approximate Frechet distance b/w them.

        Returns:
            Nothing is returned. However, pathDic is updated to store the newly computed
            distances.

*/

void computeDistances(
    std::unordered_map<int, Traj>& trajs,
    const SimpleTraj& simpTraj,
    const SimpleTraj& pth,
    double r,
    std::vector<int>& possStarts,
    std::vector<int>& possEnds,
    const Pathlet& upth,
    std::unordered_map<std::pair<Pathlet, SubTraj>, double>& pathDic
);


/*

Compute approximate Frechet distances b/w subtrajectories and pathlets.

    These distances are used as inputs to the greedy algorithm. We avoid computing
    all pairwise distances by overlaying a grid over the points, and discarding
    pathlet-subtrajectiry pairs whose endpoints do not lie in neighboring cells.

    Args:
        trajs ({int: traj}): dictionary of trajectories.
        rmin, rmax (float): lower and upper bounds on the Frechet distances.

    Returns:
        A dictionary of the form {(pathlet, subTraj) : float} storing the distance
        b/w pathlet and subtrajectory.

*/
std::unordered_map<std::pair<Pathlet, SubTraj>, double> process(
    std::unordered_map<int, Traj>& trajs,
    double rmin,
    double rmax,
    int simplify);