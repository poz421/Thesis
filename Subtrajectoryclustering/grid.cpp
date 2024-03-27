#include "grid.h"

TrajGrid::TrajGrid(std::map<std::pair<int, int>, std::unordered_map<int, std::vector<int>>> data, double startLat, double startLon, double delta, int numX, int numY)
    : data(data), startLat(startLat), startLon(startLon), delta(delta), numX(numX), numY(numY) {}




TrajGrid gridData(std::unordered_map<int, Traj> trajs, std::vector<SimpleTraj> simpTrajs, double delta) {
    //Compute the x and y extent of the points.
    
    double xMin = std::numeric_limits<double>::max();
    double yMin = std::numeric_limits<double>::max();
    double xMax = std::numeric_limits<double>::min();
    double yMax = std::numeric_limits<double>::min();

    for (const SimpleTraj& simpTraj : simpTrajs) {
        int trID = simpTraj.trajID;
        for (int index : simpTraj.indices) {
            const Pt& p = trajs.at(trID).pts.at(index);
            xMin = std::min(xMin, p.lon);
            xMax = std::max(xMax, p.lon);
            yMin = std::min(yMin, p.lat);
            yMax = std::max(yMax, p.lat);
        }
    }

    /*
     Initialize grid which is a dictionary mapping from the cell (xIndex, yIndex)
     to the contents of the cell.
     grid is of the form {(xIndex, yIndex) : {trajID : [indices of points in simpTraj]}}.
    
    */
    int xCells = static_cast<int>(std::ceil((xMax - xMin) / delta));
    int yCells = static_cast<int>(std::ceil((yMax - yMin) / delta));

    std::map<std::pair<int, int>, std::unordered_map<int, std::vector<int>>> grid;

    for (const SimpleTraj& simpTraj : simpTrajs) {
        int trID = simpTraj.trajID;
        for (int index : simpTraj.indices) {
            const Pt& p = trajs.at(trID).pts.at(index);
            int xCell = int(std::floor((p.lon - xMin) / delta));
            int yCell = int(std::floor((p.lat - yMin) / delta));

            if (grid.find(std::make_pair(xCell, yCell)) == grid.end()) {
                grid[std::make_pair(xCell, yCell)] = {};
            }

            if (grid[std::make_pair(xCell, yCell)].find(trID) == grid[std::make_pair(xCell, yCell)].end()) {
                grid[std::make_pair(xCell, yCell)][trID] = {};
            }

            grid[std::make_pair(xCell, yCell)][trID].push_back(index);
        }
    }

    return TrajGrid(grid, yMin, xMin, delta, xCells, yCells);
}