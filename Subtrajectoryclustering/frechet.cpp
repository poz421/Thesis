#include "frechet.h"

double sqDist(const Pt& pt1, const Pt& pt2) {
    double latDiff = pt1.lat - pt2.lat;
    double lonDiff = pt1.lon - pt2.lon;
    return (latDiff * latDiff) + (lonDiff * lonDiff);
}

double distPtSegment(const Pt& p, const Pt& q1, const Pt& q2) {
    double x = p.lon;
    double y = p.lat;
    double x1 = q1.lon;
    double y1 = q1.lat;
    double x2 = q2.lon;
    double y2 = q2.lat;

    if (x1 == x2 && y1 == y2) { // Degenerate segment.
        return std::sqrt(sqDist(p, q1));
    }

    if (x1 == x2) { // Vertical segment.
        if ((y1 <= y && y <= y2) || (y2 <= y && y <= y1)) {
            return std::abs(x - x1);
        }
        else {
            return std::min(std::sqrt(sqDist(p, q1)), std::sqrt(sqDist(p, q2)));
        }
    }
    else if (y1 == y2) { // Horizontal segment.
        if ((x1 <= x && x <= x2) || (x2 <= x && x <= x1)) {
            return std::abs(y - y1);
        }
        else {
            return std::min(std::sqrt(sqDist(p, q1)), std::sqrt(sqDist(p, q2)));
        }
    }
    else {
        //Translate so that (x1,y1) is at the origin.
        x = x - x1;
        y = y - y1;
        x2 = x2 - x1;
        y2 = y2 - y1;
        double m = y2 / x2;
        double c = y + x / m;

        //Projection of(x, y) on line passing through origin and (x2, y2).
        double x3 = c / (m + 1 / m);
        double y3 = m * c / (m + 1 / m);

        if (x2 * x3 + y2 * y3 >= 0) { // (x3, y3) between origin and (x2, y2).
            return std::sqrt((x3 - x) * (x3 - x) + (y3 - y) * (y3 - y));
        }
        else {
            return std::min(std::sqrt(sqDist(p, q1)), std::sqrt(sqDist(p, q2)));
        }
    }
}



bool frechetDec(const std::vector<Pt>& trajA, const std::vector<Pt>& trajB, double delta) {
    std::vector<std::pair<int, int>> ptQueue;
    std::set<std::pair<int, int>> visited;
    ptQueue.push_back(std::make_pair(0, 0));

    while (!ptQueue.empty()) {
        std::pair<int, int> current = ptQueue[0];
        ptQueue.erase(ptQueue.begin());

        if (visited.find(current) != visited.end()) {
            continue;
        }

        if (current == std::make_pair(int(trajA.size() - 1), int( trajB.size() - 1))) {
            return true;
        }

        visited.insert(current);
        int i = current.first;
        int j = current.second;
        if (i == 0 && j == 0)
            if (sqDist(trajA[0], trajB[0]) > delta * delta)
                return false;

        // Bounds check, add points that are within delta (using squared distance)
        if (i + 1 < trajA.size() && sqDist(trajA[i + 1], trajB[j]) <= delta * delta) {
            ptQueue.push_back(std::make_pair(i + 1, j));
        }

        if (j + 1 < trajB.size() && sqDist(trajA[i], trajB[j + 1]) <= delta * delta) {
            ptQueue.push_back(std::make_pair(i, j + 1));
        }

        if (i + 1 < trajA.size() && j + 1 < trajB.size() && sqDist(trajA[i + 1], trajB[j + 1]) <= delta * delta) {
            ptQueue.push_back(std::make_pair(i + 1, j + 1));
        }
    }

    return false;
}



// Function to compute the semi-continuous Frechet distance
bool semiContFrechetDec(std::vector<Pt>& trajA, std::vector<Pt>& trajB, double delta) {
    trajA.push_back(trajA.back());//traj[trajA.size()-1]
    trajB.push_back(trajB.back());

    std::vector<std::pair<int, int>> ptQueue;
    ptQueue.push_back(std::make_pair(0, 0));
    std::set<std::pair<int, int>> visited;

    while (!ptQueue.empty()) {
        std::pair<int, int> current = ptQueue[0];
        ptQueue.erase(ptQueue.begin());

        if (visited.find(current) != visited.end()) {
            continue;
        }

        if (current.first == (trajA.size() - 1) && current.second == (trajB.size() - 1)) {
            return true;
        }

        visited.insert(current);

        int i = current.first;
        int j = current.second;

        if (i + 1 < trajA.size()) {
            if (j != 0) {
                Pt seg1 = trajA[i];
                Pt seg2 = trajA[i + 1];
                Pt seg3 = trajB[j - 1];
                Pt seg4 = trajB[j];
                if (std::max(distPtSegment(seg2, seg3, seg4), distPtSegment(seg4, seg1, seg2)) <= delta) {
                    ptQueue.push_back({ i + 1, j });
                }
            }
            else {
                if (sqDist(trajA[i + 1], trajB[j]) <= delta * delta) {
                    ptQueue.push_back({ i + 1, j });
                }
            }
        }

        if (j + 1 < trajB.size()) {
            if (i != 0) {
                Pt seg1 = trajA[i - 1];
                Pt seg2 = trajA[i];
                Pt seg3 = trajB[j];
                Pt seg4 = trajB[j + 1];
                if (std::max(distPtSegment(seg2, seg3, seg4), distPtSegment(seg4, seg1, seg2)) <= delta) {
                    ptQueue.push_back({ i, j + 1 });
                }
            }
            else {
                if (sqDist(trajA[i], trajB[j + 1]) <= delta * delta) {
                    ptQueue.push_back({ i, j + 1 });
                }
            }
        }

        if ((i + 1) < trajA.size() && (j + 1) < trajB.size()) {
            Pt seg1 = trajA[i];
            Pt seg2 = trajA[i + 1];
            Pt seg3 = trajB[j];
            Pt seg4 = trajB[j + 1];
            if (std::max(distPtSegment(seg2, seg3, seg4), distPtSegment(seg4, seg1, seg2)) <= delta) {
                ptQueue.push_back({ i + 1, j + 1 });
            }
        }
    }
    return false;
}
