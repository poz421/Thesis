#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "base.h"
#include "heapdict.cpp"


/*

Compute pre-requisite data structures for the greedy algorithm.

        Args:
            trajs ({int : traj}): dict mapping ID to traj objects.
            distPairs ({(pathlet, int) : [(subTraj, float)]): dict mapping a
                    pathlet-trajID pair to a list of subtraj-float pairs, where the
                    subtraj belong to the resp. traj, and the float values are computed
                    Frechet distances.

        Returns:
            A 4-tuple (strajCov, ptStraj, strajPth, trajCov), where
            strajCov ({subtraj : int}) : dict storing coverage (#points) in all subtrajs.<
                    in distPairs.
            ptStraj ({pt : {subtraj}}) : dict storing for each point, the set of subtrajs.
                    in distPairs containing it.
            strajPth ({subtraj: [pathlet]) : dict storing for each pathlet (subtrajectory?) in distPairs,
                    the list of pathlets associated with it.
            trajCov ({int : int}) : dict storing the #points in each trajectory.

*/

std::tuple<std::unordered_map<SubTraj, int, SubTraj::Hash >,
    std::unordered_map<Pt, std::unordered_set<SubTraj, SubTraj::Hash>, Pt::Hash>,
    std::unordered_map<SubTraj, std::vector<Pathlet>, SubTraj::Hash>,
    std::unordered_map<int, int>>

    preprocessGreedy(
        const std::unordered_map<int, Traj>& trajs,
        const std::unordered_map<std::pair<Pathlet, int>, std::vector<std::pair<SubTraj, double>>> distPairs);



/*

Process a point picked in an interation of the greedy algorithm.

        Args:
            p (pt): Point to be processed.
            ptStraj ({pt : {subtraj}}) : dict storing for each point, the set of subtrajs
                    in distPairs containing it.
            strajCov ({subtraj : int}) : dict storing coverage (#points) in all subtrajs
                    in distPairs.
            strajPth ({subtraj: [pathlet]) : dict storing for each subtraj in distPairs,
                    the list of pathlets associated with it.
            trajs ({int : traj}) : dict mapping IDs to trajectories.
            trajCov ({int : int}) : dict storing the #points in each trajectory.
            distPairs ({(pathlet, int) : [(subTraj, float)]): dict mapping a
                    pathlet-trajID pair to a list of subtraj-float pairs, where the
                    subtraj belong to the resp. traj, and the float values are computed
                    Frechet distances.
            numUnprocessedPts (int) : no. of points left to be processed.
            queue : priority queue

        Returns:
            Set of the form {(pathlet, int}) containing pathlet-trajID pairs that the point
            can be assigned to.
*/

int processPoint(
    const Pt& p,
    std::unordered_map<Pt, std::unordered_set<SubTraj, SubTraj::Hash>, Pt::Hash>& ptStraj,
    std::unordered_map<SubTraj, int, SubTraj::Hash>& strajCov,
    std::unordered_map<SubTraj, std::vector<Pathlet>, SubTraj::Hash>& strajPth,
    const std::unordered_map<int, Traj>& trajs,
    std::unordered_map<int, int>& trajCov,
    const std::unordered_map<std::pair<Pathlet, int>, std::vector<std::pair<SubTraj, double>>>& distPairs,
    int& numUnprocessedPts,
    HeapDict<int, double>& queue,
    std::unordered_map<Pt, bool, Pt::Hash>& processed
);

/*

Process the points of a subtrajectory in an iteration of the greedy algorithm.

        Args:
            straj (subtraj): subtrajectory whose points are to be processed.
            strajCov ({subtraj : int}) : dict storing coverage (#points) in all subtrajs
                    in distPairs.
            trajs ({int : traj}) : dict mapping IDs to trajectories.
            trajCov ({int : int}) : dict storing the #points in each trajectory.
            ptStraj ({pt : {subtraj}}) : dict storing for each point, the set of subtrajs
                    in distPairs containing it.
            strajPth ({subtraj: [pathlet]) : dict storing for each subtraj in distPairs,
                    the list of pathlets associated with it.
            distPairs ({(pathlet, int) : [(subTraj, float)]): dict mapping a
                    pathlet-trajID pair to a list of subtraj-float pairs, where the
                    subtraj belong to the resp. traj, and the float values are computed
                    Frechet distances.
            numUnprocessedPts (int) : no. of points left to be processed.
            queue : priority queue

        Returns:
            Set of the form {(pathlet, int}) containing pathlet-trajID pairs that the point
            can be assigned to.

*/

std::unordered_set<std::pair<Pathlet, int>> processSubtraj(
    SubTraj straj,
    std::unordered_map<SubTraj, int, SubTraj::Hash>& strajCov,
    std::unordered_map<int, Traj>& trajs,
    std::unordered_map<int, int>& trajCov,
    std::unordered_map<Pt, std::unordered_set<SubTraj, SubTraj::Hash>, Pt::Hash>& ptStraj,
    std::unordered_map<SubTraj, std::vector<Pathlet>, SubTraj::Hash>& strajPth,
    std::unordered_map<std::pair<Pathlet, int>, std::vector<std::pair<SubTraj, double>>>& distPairs,
    int& numUnprocessedPts,
    HeapDict<int, double>& queue,
    std::unordered_map<Pt, bool, Pt::Hash>& processed
);


/*

Process the unprocessed points of a trajectory.

        This function is called when an interation of the greedy algorithm decides to leave the unprocessed
        points of the trajectory unassigned.

        Args:
            trID (int): ID of the traj whose points are to be processed.
            ptStraj ({pt : {subtraj}}) : dict storing for each point, the set of subtrajs
                    in distPairs containing it.
            strajCov ({subtraj : int}) : dict storing coverage (#points) in all subtrajs
                    in distPairs.
            strajPth ({subtraj: [pathlet]) : dict storing for each subtraj in distPairs,
                    the list of pathlets associated with it.
            trajs ({int : traj}) : dict mapping IDs to trajectories.
            trajCov ({int : int}) : dict storing the #points in each trajectory.
            distPairs ({(pathlet, int) : [(subTraj, float)]): Dictionary mapping a
                    pathlet-trajID pair to a list of subtraj-float pairs, where the
                    subtraj belong to the resp. traj, and the float values are computed
                    Frechet distances.
            numUnprocessedPts (int) : no. of points left to be processed.
            queue : priority queue

        Returns:
            Set of the form {(pathlet, int}) containing pathlet-trajID pairs that the point
            can be assigned to.

*/


std::unordered_set<std::pair<Pathlet, int>>
processTraj(int trID,
    std::unordered_map<Pt, std::unordered_set<SubTraj, SubTraj::Hash>, Pt::Hash>& ptStraj,
    std::unordered_map<SubTraj, int, SubTraj::Hash>& strajCov,
    std::unordered_map<SubTraj, std::vector<Pathlet>, SubTraj::Hash>& strajPth,
    const std::unordered_map<int, Traj>& trajs,
    std::unordered_map<int, int>& trajCov,
    const std::unordered_map<std::pair<Pathlet, int>, std::vector<std::pair<SubTraj, double>>>& distPairs,
    int& numUnprocessedPts,
    HeapDict<int, double>& queue,
    std::unordered_map<int, int>& unassignedPts,
    std::unordered_map<Pt, bool, Pt::Hash>& processed);


/*
    Compute coverage-cost ratio for a pathlet.

        Args:
            trajStrajDist {int : (subtraj, float)} : dict containing a subtrajectory for
                    a trajectory ID alongwith the distance from the pathlet.
            c1, c3 (float): parameters of the greedy algorithm.
            strajCov ({subtraj : int}) : dict storing coverage (#points) of subtrajs.

        Returns:
            The coverage-cost ratio of the pathlet.


*/

double computeCovCostRatio(const std::unordered_map<int,
    std::pair<SubTraj, double>>&trajStrajDist,
    float c1,
    float c3,
    const std::unordered_map<SubTraj, int, SubTraj::Hash>& strajCov);


/*

Help in computing the subtraj of a trajectory with the maximum coverage-cost ratio.

        For a guess 'r' on the coverage-cost ratio, computes the subtraj s that maximizes the quantity
        (coverage - c3*r*distance), and if it is > 0, returns the (s,distance), else returns (None, None).

        Args:
            strajDists ([(subtraj, float)]): list of subtraj-distance pairs from a single trajectory.
            strajCov ({subtraj : int}): dict storing coverage (#points) of subtrajs.
            r (float): guess on coverage-cost ratio.
            c3 : parameter of greedy algorithm.

        Returns:
            (subtraj, float) or (None,None).

*/

std::pair<SubTraj, double> optStrajAdvancedHelper(const std::vector<std::pair<SubTraj, double>>& strajDists,
    std::unordered_map<SubTraj, int, SubTraj::Hash>& strajCov,
    double r,
    float c3);

/*


Compute subtrajectories for a pathlet with optimal coverage-cost ratio.

        Args:
            pth (pathlet) : pathlet of concern.
            distPairs ({(pathlet, int) : [(subTraj, float)]): dict mapping a
                    pathlet-trajID pair to a list of subtraj-float pairs, where the
                    subtraj belong to the resp. traj, and the float values are computed
                    Frechet distances.
            pthOptStrajs ({pathlet : {int : (subtraj, float)}}): dict storing for a pathlet,
                    the (at most one) optimal subtraj. from  each traj. along with the distance.
                    This is updated by the function.
            strajCov ({subtraj : int}) : dict storing coverage (#points) of subtrajs.
            c1, c3 (float) : parameters of the greedy algorithm.
            m (int) : total no. of points.
            affectedTrajs [int] : list of trajIDs with newly covered points.

        Returns:
            Nothin, but updates pthOptStrajs.

*/

void computeOptStrajsAdvanced(
    Pathlet pth,
    const std::unordered_map<std::pair<Pathlet, int>, std::vector<std::pair<SubTraj, double>>>& distPairs,
    std::unordered_map<Pathlet, std::unordered_map<int, std::pair<SubTraj, double>>, Pathlet::Hash>& pthOptStrajs,
    std::unordered_map<SubTraj, int, SubTraj::Hash >& strajCov,
    float c1,
    float c3,
    int m,
    const std::vector<int>& affectedTrajs,
    int maxIndex);


/*


Run the greedy algorithm for pathlet cover.

    At each step, the algorithm either chooses to leave a point unassigned, or picks
    a pathlet and a set of subtrajectories assigned to the pathlet (at most one from
    each subtrajectory) depending on whichever has the highest coverage-cost ratio.
    The points that are covered by the sets picked up in each greedy step are said to
    be "processed".

    Args:
        trajs ({int : traj}): dict mapping ID to traj objects.
        distPairs ({(pathlet, int) : [(subTraj, float)]): dict mapping a
                    pathlet-trajID pair to a list of subtraj-float pairs, where the
                    subtraj belong to the resp. traj, and the float values are computed
                    Frechet distances.
        strajCov ({subtraj : int}) : dict storing coverage (#points) in all subtrajs
                    in distPairs.
        ptStraj ({pt : {subtraj}}) : dict storing for each point, the set of subtrajs
                    in distPairs containing it.
        strajPth ({subtraj: [pathlet]) : dict storing for each subtraj in distPairs,
                the list of pathlets associated with it.
        trajCov ({int : int}) : dict storing the #points in each trajectory.
        c1,c2,c3 (float): parameters of the greedy algorithm.

    Returns:
        Pathlet assignments and unassigned points as determined by the greedy algorithm,
        alongwith other relevant info about the pathlets picked.

*/

std::tuple<std::unordered_map<Pathlet, std::vector<SubTraj>, Pathlet::Hash>,
    std::vector<std::tuple<Pathlet, float, int, float>>,
    std::unordered_map<int, int>>
    runGreedy(
        std::unordered_map<int, Traj>& trajs,
        std::unordered_map<std::pair<Pathlet, int>, std::vector<std::pair<SubTraj, double>>>& distPairs,
        std::unordered_map<SubTraj, int, SubTraj::Hash >& strajCov,
        std::unordered_map < Pt, std::unordered_set<SubTraj, SubTraj::Hash>, Pt::Hash >& ptStraj,
        std::unordered_map<SubTraj, std::vector<Pathlet>, SubTraj::Hash>& strajPth,
        std::unordered_map<int, int >& trajCov,
        float c1, float c2, float c3);