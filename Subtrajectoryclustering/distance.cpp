#include "distance.h"
#include <cmath>
#include <chrono>
#include "unordered_set"

// Function to convert degrees to radians
double deg2rad(double deg) {
    return (deg * 3.14 / 180.0);
}

// Calculate Haversine distance between two points
double haversineDistance(double lat1, double lon1, double lat2, double lon2) {
    // Radius of the Earth in kilometers
    const double R = 6371.0;

    // Convert latitude and longitude from degrees to radians
    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    lat2 = deg2rad(lat2);
    lon2 = deg2rad(lon2);

    // Haversine formula
    double dlon = lon2 - lon1;
    double dlat = lat2 - lat1;
    double a = std::sin(dlat / 2) * std::sin(dlat / 2) +
        std::cos(lat1) * std::cos(lat2) *
        std::sin(dlon / 2) * std::sin(dlon / 2);
    double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));
    double distance = R * c;

    return distance;
}


bool inSquare(const Pt& center, const Pt& point, double r) {
    return std::abs(center.lat - point.lat) <= r && std::abs(center.lon - point.lon) <= r;
}

std::vector<int> selectEquidistantPoints(int lowerBounds,int upperBounds,int d)
{
    std::vector<int> ret;
    int num = int((upperBounds - lowerBounds) / d);
    if (num == 0)
        num = 1;
    for (int i = lowerBounds; i < upperBounds; i += num)
        ret.push_back(i);
    return ret;
}

std::vector<bool> removeMarkedPathlets(std::vector<std::vector<Pt>>& pathlets, double omega)
{
    std::vector<bool> ret(pathlets.size(), true);
    for(int i=0;i<pathlets.size();i++)
        for (int j = i + 1; j < pathlets.size(); j++)
        { 
            std::vector<Pt> first = pathlets[i];
            std::vector<Pt> second = pathlets[j];
            bool found = false;

            for (int k = 0; k < first.size(); k++)
            {
                for (int l = 0; l < second.size(); l++)
                    if (haversineDistance(first[k].lat,first[k].lon,second[l].lat,second[l].lon) < omega)
                    {
                        
                        ret[j] = false;
                        found = true;
                        break;
                    }
                if (found)
                    break;
            }
        }
    return ret;
}

std::unordered_map<int, std::vector<int>> findNeighboringPts(const TrajGrid& grid, const Pt& point) {
    // Find grid cell containing input point.

    
    int xCell = static_cast<int>(std::floor((point.lon - grid.startLon) / grid.delta));
    int yCell = static_cast<int>(std::floor((point.lat - grid.startLat) / grid.delta));

    std::unordered_map<int, std::vector<int>> possibleTrajs;

    for (int i = -1; i <2 ; i++) {
        for (int j = -1; j < 2; j++) {
            if (0 <= xCell + i && xCell + i < grid.numX && 0 <= yCell + j && yCell + j < grid.numY) {  // Boundary check
                auto gridCellIter = grid.data.find(std::make_pair(xCell + i, yCell + j));
                //If there are no pts in the cell it won't be in data dictionary.
                if (gridCellIter != grid.data.end()) {
                    for (const auto& kv : grid.data.at(std::make_pair(xCell + i, yCell + j))) {
                        int trajID = kv.first;
                        const std::vector<int>& indices = kv.second;
                       
                        // If this is the first time you're adding points from this trajectory.
                        if (possibleTrajs.find(trajID) == possibleTrajs.end()) {
                            possibleTrajs[trajID] = std::vector<int>();
                        }

                        possibleTrajs[trajID].insert(possibleTrajs[trajID].end(), indices.begin(), indices.end());
                    }
                }
            }
        }
    }

    return possibleTrajs;
}

void computeDistances(
    std::unordered_map<int, Traj>& trajs,
    const SimpleTraj& simpTraj,
    const SimpleTraj& pth,
    double r,
    std::vector<int>& possStarts,
    std::vector<int>& possEnds,
    const Pathlet& upth,
    std::unordered_map<std::pair<Pathlet, SubTraj>, double>& pathDic
)

{
    if (possStarts.empty() || possEnds.empty()) {
        return;
    }



    std::sort(possStarts.begin(), possStarts.end());
    std::sort(possEnds.begin(), possEnds.end());


    std::vector<Pt> pathletPoints;

    for (int index : pth.indices) {
        pathletPoints.push_back(trajs[pth.trajID].pts[index]);
    }



    //std::vector<SubTraj> subTrajMatches;
    //std::vector<int> startSkips;

    const std::vector<int> indices = simpTraj.indices;
    ;
    int curBound = possStarts[0];

    for (int i = 0; i < possStarts.size(); i++) {
        for (int j = int(possEnds.size() - 1); j > -1; j--) {
            if (possEnds[j] < possStarts[i]) {
                break;
            }
            if (possEnds[j] <= curBound) {
                break;
            }
            
            SubTraj subTrajectory(simpTraj.trajID, std::make_pair(possStarts[i], possEnds[j]));
             auto pair = std::make_pair(upth, subTrajectory);

            if (pathDic.find(pair) != pathDic.end()) {
                 curBound = possEnds[j];
                 break;
             }

            std::vector<Pt> subTrajPoints;
            for (int k = 0; k < indices.size(); ++k) {
                if (possStarts[i] <= indices[k] && indices[k] <= possEnds[j]) {
                    subTrajPoints.push_back(trajs[simpTraj.trajID].pts[indices[k]]);
                }
            }

            if (frechetDec(subTrajPoints, pathletPoints,  3*r)) {
                pathDic[std::make_pair(upth, subTrajectory)] =  r;
                curBound = possEnds[j];
                break;
            }

        }
    }
}



std::unordered_map<std::pair<Pathlet, SubTraj>, double> process(
    std::unordered_map<int, Traj>& trajs,
    double rmin,
    double rmax,
    int simplify) {
    auto start_time = std::chrono::high_resolution_clock::now();

    double r = rmin;
    std::unordered_map<std::pair<Pathlet, SubTraj>, double> pathDic;
    std::vector<Pathlet> pathlets;


    // Create "unsimplified" canonical pathlets.
    for (const auto& tr : trajs) {
        if (tr.second.pts.size() >= 32)
            for (const auto& b : canonise(0, int(tr.second.pts.size() - 1))) {
                pathlets.push_back(Pathlet(tr.first, std::make_pair(b.first, b.second)));
            }
        else
            pathlets.push_back(Pathlet(tr.first,std::make_pair(0,tr.second.pts.size()-1)));
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //std::cout << "Computing pathlet: " << duration.count() << " milliseconds" << std::endl;

    start_time = std::chrono::high_resolution_clock::now();
    //std::cout << "size before sparsify" << std::endl;
    //std::cout << pathlets.size() << std::endl;
    
   
    
    
    if (simplify == 1) {
        int d = 2;




        std::vector<std::vector<Pt>> simplePoints;
        for (auto p : pathlets)
        {
            auto simpleIndex = selectEquidistantPoints(p.bounds.first, p.bounds.second, d);
            std::vector<Pt> points;

            for (int index : simpleIndex)
            {
                Traj t = trajs[p.trajID];
                points.push_back(t.pts[index]);

            }
            simplePoints.push_back(points);
        }


        
        double omega = 0.0005;
        std::vector<bool> patToMantain = removeMarkedPathlets(simplePoints, omega);

        std::vector<Pathlet> pathletNew;

        for (int i = 0; i < patToMantain.size(); i++)
        {
            if (patToMantain[i])
            {
                pathletNew.push_back(pathlets[i]);
            }
        }

        pathlets = pathletNew;

        //std::cout << "size after sparsify" << std::endl;


        end_time = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        //std::cout << "samplying pathlets execution time: " << duration.count() << " milliseconds" << std::endl;

    }
    std::cout <<"NpathIn;" << pathlets.size() << std::endl;

    auto start1=std::chrono::high_resolution_clock::now();
    std::chrono::milliseconds distanceTime(0), candidateTime(0), simplyTime(0);
    while (r <= rmax) {
        start_time = std::chrono::high_resolution_clock::now();
        std::unordered_map<int, SimpleTraj> fSimpTrajs, bSimpTrajs;

        // (r/2)-simplify trajectories, store in a map.
        for (const auto& tr : trajs) {
           
            SimpleTraj fTraj = fSimplify(tr.second, r, 0, int(tr.second.pts.size()) - 1);
            SimpleTraj bTraj = bSimplify(tr.second, r, 0, int(tr.second.pts.size()) - 1);
            fSimpTrajs[tr.first] = fTraj;
            bSimpTrajs[tr.first] = bTraj;
        }

        std::vector<SimpleTraj> fSimpTrajsVector(fSimpTrajs.size());
        std::transform(fSimpTrajs.begin(), fSimpTrajs.end(), fSimpTrajsVector.begin(),
            [](const auto& pair) { return pair.second; });

        std::vector<SimpleTraj> bSimpTrajsVector(bSimpTrajs.size());
        std::transform(bSimpTrajs.begin(), bSimpTrajs.end(), bSimpTrajsVector.begin(),
            [](const auto& pair) { return pair.second; });

       


        // Create grid of resolution r for fSimpTrajs and bSimpTrajs.
        TrajGrid fGrid = gridData(trajs, fSimpTrajsVector, r);
        TrajGrid bGrid = gridData(trajs, bSimpTrajsVector, r);


        end_time = std::chrono::high_resolution_clock::now();
        simplyTime+= std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        start_time= std::chrono::high_resolution_clock::now();
        // Find candidate subtrajectories for each pathlet.
        for (const auto& pth : pathlets) {
            // Simplify pathlet.
           auto start_time2 = std::chrono::high_resolution_clock::now();
            
            SimpleTraj forwardPathlet = fSimplify(trajs.at(pth.trajID), r, pth.bounds.first, pth.bounds.second);
            SimpleTraj backwardPathlet = bSimplify(trajs.at(pth.trajID), r, pth.bounds.first, pth.bounds.second);

            // Forward simplified subtrajectories with endpoints close to those of pth.
            std::unordered_map<int, std::vector<int>> startFNeighbors = findNeighboringPts(fGrid, trajs.at(pth.trajID).pts[pth.bounds.first]);
            std::unordered_map<int, std::vector<int>> endFNeighbors = findNeighboringPts(fGrid, trajs.at(pth.trajID).pts[pth.bounds.second]);
            std::vector<int> possibleFSubTrajs;

            
            for (const auto& entry : startFNeighbors) {
                if (endFNeighbors.find(entry.first) != endFNeighbors.end()) {
                    possibleFSubTrajs.push_back(entry.first);
                }
            }

            
           
            
            
            // Backward simplified subtrajectories with endpoints close to those of pth.
            std::unordered_map<int, std::vector<int>> startBNeighbors = findNeighboringPts(bGrid, trajs.at(pth.trajID).pts[pth.bounds.first]);
            std::unordered_map<int, std::vector<int>> endBNeighbors = findNeighboringPts(bGrid, trajs.at(pth.trajID).pts[pth.bounds.second]);
            std::vector<int> possibleBSubTrajs;

            for (const auto& entry : startBNeighbors) {
                if (endBNeighbors.find(entry.first) != endBNeighbors.end()) {
                    possibleBSubTrajs.push_back(entry.first);
                }
            }
          


           
            
            // Record distance between pathlet and itself to be 0.
            SubTraj subpath = SubTraj(pth.trajID, pth.bounds);
            if (pathDic.find(std::make_pair(pth,subpath)) == pathDic.end()) {
                pathDic[std::make_pair(pth, subpath)] = 0.0;
            }
            
            auto end_time2 = std::chrono::high_resolution_clock::now();
            candidateTime+= std::chrono::duration_cast<std::chrono::milliseconds>(end_time2 - start_time2);

            //start_time = std::chrono::high_resolution_clock::now();
            // Compute distance between forwardPathlet and possibleFSubTrajs.
            for (int trajID : possibleFSubTrajs) {
                SimpleTraj simpTraj = fSimpTrajs.at(trajID);
                computeDistances(trajs, simpTraj, forwardPathlet, r, startFNeighbors[trajID], endFNeighbors[trajID], pth, pathDic);
            }
            
           

            // Compute distance between backwardPathlet and possibleBSubTrajs.
            for (int trajID : possibleBSubTrajs) {
                SimpleTraj simpTraj = bSimpTrajs.at(trajID);
                computeDistances(trajs, simpTraj, backwardPathlet, r, startBNeighbors[trajID], endBNeighbors[trajID], pth, pathDic);
            }
            //end_time = std::chrono::high_resolution_clock::now();
            //distanceTime += std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        }
        end_time = std::chrono::high_resolution_clock::now();
        distanceTime += std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        r *=2 ;
    }
    
    auto end1 = std::chrono::high_resolution_clock::now();
    auto prova= std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1);

    //std::cout << "total " << prova.count() << " milliseconds" << std::endl;
    //std::cout << "simplyfy trajectory execution time: " << simplyTime.count() << " milliseconds" << std::endl;
   // std::cout << "Find candidate pathlet for distance execution time: " << candidateTime.count() << " milliseconds" << std::endl;
    //std::cout << "calculate distances +find candidate execution time: " << distanceTime.count() << " milliseconds" << std::endl;

    return pathDic;
}


