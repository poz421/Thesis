#pragma once
#include "base.h"
#include <vector>
#include <set>
#include <cmath>

/*
Compute squared Euclidean distance b/w points.

    Args:
        pt1, pt2 (pt): Input points.

    Returns:
        Squared Euclidean distance b/w pt1 and pt2.

*/

double sqDist(const Pt& pt1, const Pt& pt2);



/*
Compute the distance b/w a point and a line segment.

    The distance is defined as the shortest over all distances b/w the point
    and all the points of the segment.

    Args:
        p (pt): Input point.
        seg (pt,pt): Endpoints of input segment.

    Returns:
        The distance b/w p and seg.
*/

double distPtSegment(const Pt& p, const Pt& q1, const Pt& q2);


/*

Decide if the discrete Frechet distance b/w trajectories is at most delta.

    Uses the classic dynamic programming algorithm to find the optimal correspondence.

    Args:
        trajA, trajB ([pt]): input trajectories, represented as list of pt objects.
        delta (float): guess on the discrete Frechet distance.

    Returns:
        True if the discrete Frechet distance b/w trajA and trajB is at most delta,
        False otherwise

*/
bool frechetDec(const std::vector<Pt>& trajA, const std::vector<Pt>& trajB, double delta);


/*
Decide if the semi-contiuous Frechet distance b/w trajectories is at most delta.

    The semi-continuous Frechet distance is defined similarly to its discrete cousin,
    except that the correspondences are defined b/w points on one trajectory and segments
    on the other. The cost of a correspondence pair is calculated as the distance b/w the
    point and segment. The decision procedure is based on dynamic programming.

    Args:
        trajA, trajB ([pt]): input trajectories, represented as list of pt objects.
        delta (float): guess on the discrete Frechet distance.

    Returns:
        True if the semi-continuous Frechet distance b/w trajA and trajB is at most delta,
        False otherwise.

*/

bool semiContFrechetDec(std::vector<Pt>& trajA, std::vector<Pt>& trajB, double delta);