#include "distanceUtils.h"
#include "frechet.h"

//default constructor
SimpleTraj::SimpleTraj() : trajID(-1) {}

// Parameterized constructor
SimpleTraj::SimpleTraj(int trajID, const std::vector<int>& indices) : trajID(trajID), indices(indices) {}



SimpleTraj fSimplify(const Traj& uTraj, double tau, int start, int end) {
    SimpleTraj sTraj;
    sTraj.trajID = uTraj.pts[0].trajID;
    sTraj.indices.push_back(start);
    int curPoint = start;

    for (int i = start + 1; i < end + 1; i++) {
        if (i == end) {
            sTraj.indices.push_back(i);
            continue;
        }

        if (sqDist(uTraj.pts[curPoint], uTraj.pts[i]) >= tau * tau) {
            sTraj.indices.push_back(i);
            curPoint = i;
        }
    }

    return sTraj;
}


SimpleTraj bSimplify(const Traj& uTraj, double tau, int start, int end) {
    SimpleTraj sTraj;
    sTraj.trajID = uTraj.pts[0].trajID;
    sTraj.indices.push_back(end);
    int curPoint = end;

    for (int i = end - 1; i > start-1; --i) {
        if (i == start) {
            sTraj.indices.push_back(i);
            continue;
        }

        if (sqDist(uTraj.pts[curPoint], uTraj.pts[i]) >= tau * tau) {
            sTraj.indices.push_back(i);
            curPoint = i;
        }
    }

    std::reverse(sTraj.indices.begin(), sTraj.indices.end());
    return sTraj;
}

std::vector<std::pair<int, int>> canonise(int i, int j) {
    // Hard-coded lower bound on the length of canonical pathlets
    if (j - i <= 50) {
        return {};
    }

    std::vector<std::pair<int, int>> retval;
    retval.push_back(std::make_pair(i, j));
    int midpoint = (i + j) / 2;
    auto leftIntervals = canonise(i, midpoint);
    auto rightIntervals = canonise(midpoint + 1, j);

    retval.insert(retval.end(), leftIntervals.begin(), leftIntervals.end());
    retval.insert(retval.end(), rightIntervals.begin(), rightIntervals.end());

    return retval;
}




