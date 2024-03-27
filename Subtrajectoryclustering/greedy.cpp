#include "greedy.h"
#include <iomanip>
#include <fstream>
#include <chrono>
int iter = 0;
int iterGreedy = 0, iterProcc = 0;
std::vector<std::pair<Pathlet, int>> retValProva;
std::vector<std::pair<Pathlet, int>> retValProctraj;

std::tuple<std::unordered_map<SubTraj, int, SubTraj::Hash >,
    std::unordered_map<Pt, std::unordered_set<SubTraj, SubTraj::Hash>, Pt::Hash>,
    std::unordered_map<SubTraj, std::vector<Pathlet>, SubTraj::Hash>,
    std::unordered_map<int, int>>

    preprocessGreedy(

        const std::unordered_map<int, Traj>& trajs,

        const std::unordered_map<std::pair<Pathlet, int>, std::vector<std::pair<SubTraj, double>>> distPairs)
{
    std::unordered_map<SubTraj, int, SubTraj::Hash> strajCov;
    std::unordered_map<Pt, std::unordered_set<SubTraj, SubTraj::Hash>, Pt::Hash> ptStraj;
    std::unordered_map<SubTraj, std::vector<Pathlet>, SubTraj::Hash> strajPth;
    std::unordered_map<int, int> trajCov;

    auto start_time = std::chrono::high_resolution_clock::now();
    std::chrono::milliseconds execution(0), execution1(0);

    for (const auto& entry : distPairs) {
        const Pathlet& pth = entry.first.first;
        int trID = entry.first.second;
        const std::vector<std::pair<SubTraj, double>>& value = entry.second;

        for (const auto& subtrajDistPair : value) {
            const SubTraj& straj = subtrajDistPair.first;
            double dist = subtrajDistPair.second;

            strajCov[straj] = straj.bounds.second - straj.bounds.first + 1;

            if (strajPth.find(straj) != strajPth.end()) {
                strajPth[straj].push_back(pth);
            }
            else {
                strajPth[straj] = { pth };
            }
            auto start_time1 = std::chrono::high_resolution_clock::now();
            for (int j = straj.bounds.first; j < straj.bounds.second + 1; j++) {
                const Pt& p = trajs.at(trID).pts.at(j);
                if (ptStraj.find(p) != ptStraj.end()) {
                    ptStraj[p].insert(straj);
                }
                else {
                    ptStraj[p] = { straj };
                }
            }
            auto end_time1 = std::chrono::high_resolution_clock::now();
            execution1 += std::chrono::duration_cast<std::chrono::milliseconds>(end_time1 - start_time1);

        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    execution = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time) - execution1;
    //std::cout << "Straj execution time: " << execution.count() << " milliseconds" << std::endl;
    //std::cout << "Computing ptStraj execution time: " << execution1.count() << " milliseconds" << std::endl;

    start_time = std::chrono::high_resolution_clock::now();
    for (const auto& entry : trajs) {
        int trID = entry.first;
        const Traj& tra = entry.second;
        trajCov[trID] = int(tra.pts.size());

        for (const Pt& p : tra.pts) {
            if (ptStraj.find(p) == ptStraj.end()) {
                ptStraj[p] = {};
            }
        }
    }

    end_time = std::chrono::high_resolution_clock::now();
    execution = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    //std::cout << "Computing Traj execution time: " << execution.count() << " milliseconds" << std::endl;

    return std::make_tuple(strajCov, ptStraj, strajPth, trajCov);

}


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
)
{

    iterProcc++;

    int index = 0;
    for (const auto& straj : ptStraj[p]) {
        int trID = straj.trajID;
        strajCov[straj] = strajCov[straj] - 1;

        for (const auto& pth : strajPth[straj]) {
            retValProva[index++] = std::make_pair(pth, trID);
            //retVal.insert({ pth, trID });
        }
    }

    ptStraj[p].clear();// this is also marking that p is processed. 
    processed[p] = true;
    numUnprocessedPts--;

    // Check if we also have singleton sets.

    if (queue.size() > 0) {//DOUBT
        trajCov[p.trajID]--;

        // Change priority of trajID to zero if its coverage is 0.
        if (trajCov[p.trajID] == 0) {
            queue.insert(p.trajID, 0);
        }
    }






    return index;
}


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
) {
    int trID = straj.trajID;
    Traj& tr = trajs[trID];
    std::unordered_set<std::pair<Pathlet, int>> retVal;

    for (int i = straj.bounds.first; i < straj.bounds.second + 1; i++) {
        Pt& p = tr.pts[i];

        // Check if point p is unprocessed.
        if (processed[p] == false) {
            int maxIndex = processPoint(p, ptStraj, strajCov, strajPth, trajs, trajCov, distPairs, numUnprocessedPts, queue, processed);
            for (int i = 0; i < maxIndex; i++)
                retVal.insert(retValProva[i]);
            //retVal.insert(result.begin(), result.end());
        }
    }

    if (strajCov[straj] != 0) {
        std::cout << "Error!! Coverage should have been 0, instead of " << strajCov[straj] << std::endl;
    }

    return retVal;
}



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
    std::unordered_map<Pt, bool, Pt::Hash>& processed)
{
    const std::vector<Pt>& points = trajs.at(trID).pts;
    std::unordered_set<std::pair<Pathlet, int>> retVal;
    int count = 0;
    std::unordered_set<std::pair<Pathlet, int>> result;
    int index = 0;
    for (const Pt& p : points) {
        // Check if p is unprocessed.

        if (processed[p] == false) {
            // Add p to the list of unassigned points.
            if (unassignedPts.find(trID) != unassignedPts.end()) {
                unassignedPts[trID]++;
            }
            else {
                unassignedPts[trID] = 1;
            }

            int maxIndex = processPoint(p, ptStraj, strajCov, strajPth, trajs, trajCov, distPairs, numUnprocessedPts, queue, processed);
            for (int i = 0; i < maxIndex; i++)
                retVal.insert(retValProva[i]);
            // retVal.insert(result.begin(), result.end());
            
        }
    }
    
    if (trajCov[trID] != 0) {
        std::cout << "Error!! Coverage should have been zero."<<trajCov[trID]<<" "<<trID << std::endl;
    }

    return retVal;
}


double computeCovCostRatio(const std::unordered_map<int, std::pair<SubTraj, double>>& trajStrajDist, float c1, float c3, const std::unordered_map<SubTraj, int, SubTraj::Hash>& strajCov) {
    float curCov = 0.0f;
    double curCost = c1;

    for (const auto& kv : trajStrajDist) {
        const SubTraj& straj = kv.second.first;
        double dist = kv.second.second;

        //DOUBT don't know how to put this condition
        if (straj.trajID == -1) {
            continue;
        }

        curCov += strajCov.at(straj);
        curCost += 1.0 * c3 * dist;
    }

    if (curCov == 0) {
        return 0;
    }
    else {
        return (1.0 * curCov) / (1.0 * curCost);
    }
}

std::pair<SubTraj, double> optStrajAdvancedHelper(const std::vector<std::pair<SubTraj, double>>& strajDists,
    std::unordered_map<SubTraj, int, SubTraj::Hash>& strajCov,
    double r,
    float c3) {
    double temp = 0;
    SubTraj stra;
    double dista = -1;


    for (int i = 0; i < strajDists.size(); i++) {
        SubTraj straj = strajDists[i].first;
        double dist = strajDists[i].second;
        int cov = strajCov.at(straj);

        if (temp < cov - c3 * r * dist) {
            temp = cov - c3 * r * dist;
            stra = straj;
            dista = dist;
        }
    }


    return std::make_pair(stra, dista);


}


void computeOptStrajsAdvanced(
    Pathlet pth,
    const std::unordered_map<std::pair<Pathlet, int>, std::vector<std::pair<SubTraj, double>>>& distPairs,
    std::unordered_map<Pathlet, std::unordered_map<int, std::pair<SubTraj, double>>, Pathlet::Hash>& pthOptStrajs,
    std::unordered_map<SubTraj, int, SubTraj::Hash >& strajCov,
    float c1,
    float c3,
    int m,
    const std::vector<int>& affectedTrajs,
    int maxIndex)
{
   
    std::unordered_map<int, std::pair<SubTraj, double>> optStrajs = pthOptStrajs.at(pth);
    std::unordered_map<int, std::pair<SubTraj, double>> ret, temp;
    if (c3 == 0) {
        // Just try to maximize coverage when c3 = 0.
        for (int trID : affectedTrajs) {
            std::vector<std::pair<SubTraj, double>> strajDists = distPairs.at(std::make_pair(pth, trID));
            // Helper's result does not depend on r when c3 = 0.
            std::pair<SubTraj, float> result = optStrajAdvancedHelper(strajDists, strajCov, 1, c3);
            ret[trID] = result;
        }
    }
    else {
        // r is a guess on the coverage-cost ratio.
        double summation = 0;  // Quantity to check when to stop iterating through values of r.
        double rmin = 1.0 / (c1 + c3);
        double r = rmin;
        double rmax = m * rmin;
        while (r <= rmax) {
            for (int i = 0; i < affectedTrajs.size();i++) {
                if (affectedTrajs[i] == iter) {
                    int trID = i;
                    if (distPairs.find(std::make_pair(pth, trID)) != distPairs.end()) {
                        std::vector<std::pair<SubTraj, double>> strajDists = distPairs.at(std::make_pair(pth, trID));
                        std::pair<SubTraj, float> result = optStrajAdvancedHelper(strajDists, strajCov, r, c3);
                        temp[trID] = result;
                        //DOUBT of this condition
                        if (result.first.trajID != -1) {
                            summation += (strajCov.at(result.first) - c3 * r * result.second);  // summation += cov - c3*r*d.
                        }
                    }
                }
            }

            if (summation < c1 * r) {
                break;
            }
            else {
                ret = temp;  // Replace ret with current subtrajs.
                temp.clear();
                summation = 0;
            }
            r *= 2;
        }
    }

   
    // Update pthOptStrajs.
    pthOptStrajs[pth] = ret;
}

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
        float c1, float c2, float c3) {



    auto start_time = std::chrono::high_resolution_clock::now();


    std::unordered_map<Pathlet, double, Pathlet::Hash> pthOptCovCost;
    std::unordered_map<Pathlet, std::unordered_map<int, std::pair<SubTraj, double>>, Pathlet::Hash> pthOptStrajs;

    // Build skeleton of pthOptStrajs.
    for (const auto& kv : distPairs) {
        Pathlet pth = kv.first.first;
        int trID = kv.first.second;
        if (pthOptStrajs.find(pth) == pthOptStrajs.end()) {
            pthOptStrajs[pth] = {};
        }
        pthOptStrajs[pth][trID] = std::make_pair(SubTraj(), -1); //default subtraj==nullSubtraj
    }
    
    std::vector<int> affectedTrajs;
    affectedTrajs.resize(trajs.size()+1);
    for (int i = 0; i < affectedTrajs.size(); i++)
        affectedTrajs[i] = -1;

    // Compute pthOptStrajs.
    for (const auto& kv : distPairs) {
        Pathlet pth = kv.first.first;
        for (const auto& subTrajDist : pthOptStrajs[pth]) {
            affectedTrajs[subTrajDist.first] = iter;
        }
         computeOptStrajsAdvanced(pth, distPairs, pthOptStrajs, strajCov, c1, c3, int(ptStraj.size()), affectedTrajs, iter);
         iter++;
    }

   

    int count = 0;

    //std::unordered_map<Pathlet, float,Pathlet::Hash> distances;
    Pathlet pth;
    // Compute pthOptCovCost from pthOptStrajs.
    for (const auto& kv : pthOptStrajs) {
        pth = kv.first;
        pthOptCovCost[pth] = computeCovCostRatio(kv.second, c1, c3, strajCov);
        if (pthOptCovCost[pth] == 0) {
            //std::cout << pth << std::endl;
            std::cout << "Error" << std::endl;
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

   // std::cout << "initialize pathlet coverage cost ratio execution time: " << duration.count() << " milliseconds" << std::endl;


    auto start_time2 = std::chrono::high_resolution_clock::now();
    // Initialize a max priority queue of pathlets ordered by coverage cost ratios.
    HeapDict<Pathlet, double> queue1;
    for (const auto& kv : pthOptCovCost) {
        queue1.insert(kv.first, -1 * kv.second);
    }



    // Initialize a priority queue of trajs, with coverage to cost ratio of a singleton set, i.e., |T|/c2.
    HeapDict<int, double> queue2;
    for (const auto& kv : trajCov) {
        int trID = kv.first;
        int cov = kv.second;
        queue2.insert(trID, -(1.0 * cov) / (1.0 * c2));
    }

    std::unordered_map<Pathlet, std::vector<SubTraj>, Pathlet::Hash> pthAssignments;
    std::vector<std::tuple<Pathlet, float, int, float>> pthStats;
    std::unordered_map<int, int> unassignedPts;

    std::unordered_map<Pt, bool, Pt::Hash> processed;
    std::unordered_map<Pt, bool, Pt::Hash> coverPoint;

    for (auto p : ptStraj)
    {
        processed.insert({ p.first, false });
        coverPoint.insert({ p.first, false });
    }
    int numUnprocessedPts = int(ptStraj.size());
    


    auto end_time2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time2 - start_time2);

    //std::cout << "creating priority queue execution time: " << duration.count() << " milliseconds" << std::endl;



    start_time = std::chrono::high_resolution_clock::now();
    std::chrono::milliseconds execution1(0);
    std::chrono::microseconds execution2(0);

    std::pair<Pathlet, double> x1;
    std::pair<int, double> x2;
    std::pair<Pathlet, double> x;
    std::pair<int, double> y;
    std::unordered_set<std::pair<Pathlet, int>> affectedPths;
    std::unordered_set<std::pair<Pathlet, int>> begin;


    int max = 0, countRet = 0;

    for (const auto p : ptStraj) {
        for (const auto& straj : ptStraj[p.first]) {
            for (const auto& pth : strajPth[straj]) {
                countRet++;
            }
        }
        if (countRet > max)
            max = countRet;
        countRet = 0;
    }




    //std::cout << "countRetval:" << trajCov.size() << std::endl;
    retValProva.reserve(max);



    int countProcessTraj = 0;
    std::unordered_map<Pathlet,int> affectedPathlets;
    affectedPathlets.reserve(pthOptCovCost.size());
    for (auto p : pthOptCovCost)
        affectedPathlets[p.first] = -1;


    for (int i = 0; i < affectedTrajs.size(); i++)
        affectedTrajs[i] = -1;
    iter = 0;
    double dist = 0;
    int coveragePoint = 0;
    int numPoint = numUnprocessedPts;
    std::unordered_map<Pathlet, int> pathFreq;

    while (numUnprocessedPts > 0) {
        iterGreedy++;
       // std::cout << "num of points is " << numUnprocessedPts << std::endl;
        
        if (queue1.size() > 0) {
            x1 = queue1.peekItem();
            x2 = queue2.peekItem();
        }
        else {
            x2 = queue2.peekItem();
            x1 = std::make_pair(Pathlet(), x2.second + 1);
        }

        affectedPths.clear();
        begin.clear();
        // Set of the form {(pth, trID)} whose coverage changes after new points are
        // processed in the current iteration.


        // If the most beneficial pathlet is more beneficial than leaving a point unassigned.
        if (x1.second <= x2.second) {
            x = queue1.popItem();
            pth = x.first;
            std::unordered_set<int> freqId;

            for (auto sub : pthOptStrajs[pth])
                
            {
                
                   

                if (sub.second.first.trajID != -1 && sub.second.second!=-1)
                {
                    if (pth.trajID == 719 && pth.bounds == std::make_pair(932, 1863))
                        std::cout << sub.second.first.trajID << std::endl;
                   
                    freqId.insert(sub.second.first.trajID);

                   // std::cout << sub.second.first << std::endl;
                    //std::cout << sub.second.second << std::endl;
                    dist += sub.second.second;


                    int id = sub.second.first.trajID;
                    int lowBd = sub.second.first.bounds.first;
                    int upBd = sub.second.first.bounds.second;
                    for (int i = lowBd; i < upBd; i++)
                    {
                        Pt point = trajs[id].pts.at(i);
                        coverPoint[point] = true;
                        coveragePoint++;
                        
                    }
                        
                }     
                   
            }
            //std::cout << "size " << freqId.size() << std::endl;

            if (pathFreq.find(pth) == pathFreq.end())
                pathFreq[pth] = freqId.size();
           // else
             //S   pathFreq[pth] += freqId.size();
            double fracThickness = 0;
            bool condition = pthAssignments.find(pth) != pthAssignments.end();
            for (const auto& kv : pthOptStrajs[pth]) {
                int trID = kv.first;
                std::pair<SubTraj, float> pair = kv.second;
                SubTraj straj = pair.first;
                if (straj.trajID == -1)
                    continue;
                if (condition) {
                    pthAssignments[pth].push_back(straj);
                }
                else {
                    pthAssignments[pth] = { straj };
                }
                fracThickness += 1.0f * strajCov[straj] / trajs[straj.trajID].pts.size();

                begin = processSubtraj(straj, strajCov, trajs, trajCov, ptStraj, strajPth, distPairs, numUnprocessedPts, queue2, processed);

                affectedPths.insert(begin.begin(), begin.end());


            }

            pthStats.push_back(std::make_tuple(pth, pthOptCovCost[pth], count, fracThickness));

        }
        else {
            countProcessTraj++;
            y = queue2.popItem();
            int trID = y.first;
            begin = processTraj(trID, ptStraj, strajCov, strajPth, trajs, trajCov, distPairs, numUnprocessedPts, queue2, unassignedPts, processed);

            affectedPths.insert(begin.begin(), begin.end());
        }

        auto start_time1 = std::chrono::high_resolution_clock::now();

       


        // Update coverage-cost ratio of affected pathlets.
        //int maxIndex = 0;
        for (const auto& pthTraID : affectedPths) {
            Pathlet path = std::get<0>(pthTraID);
            affectedPathlets[path] = iter;
            affectedTrajs[std::get<1>(pthTraID)] = iter;
            //maxIndex++;
        }

        for (auto p: affectedPathlets) {

            if (p.second == iter)
            {
                Pathlet path = p.first;
                computeOptStrajsAdvanced(path, distPairs, pthOptStrajs, strajCov, c1, c3, int(ptStraj.size()), affectedTrajs, iter);
                pthOptCovCost[path] = computeCovCostRatio(pthOptStrajs[path], c1, c3, strajCov);
                queue1.insert(path, -1.0 * pthOptCovCost[path]);
            }
        }

        auto end_time1 = std::chrono::high_resolution_clock::now();
        execution1 += std::chrono::duration_cast<std::chrono::milliseconds>(end_time1 - start_time1);

        count++;
        iter++;
    }

    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

   // std::cout << "running greedy algorithm execution time: " << duration.count() - execution1.count() << " milliseconds" << std::endl;
    //std::cout << "update queue of coverage cost ratio execution time: " << execution1.count() << " milliseconds" << std::endl;
   
    //std::cout << "coveragePoint: " << coveragePoint << std::endl;

    
    coveragePoint = 0;
    for (auto p : coverPoint)
        if (p.second)
            coveragePoint++;

    double percent = (double)coveragePoint / (double)numPoint;
    std::cout<<"CovPoint%;" << percent << std::endl;


    std::unordered_set<Pathlet> patRet;
    for (auto val : pthStats)
    {
       // std::cout << std::get<0>(val) << std::endl;
        patRet.insert(std::get<0>(val));
    }
   // std::cout << patRet.size() << std::endl;
    std::cout <<"NpathOut;" << pathFreq.size() << std::endl;
    int pathPointCount = 0;
    for (auto elem : pathFreq)
    {
       
            std::cout << elem.first.trajID<<";"<<elem.first.bounds.first<<";"<<elem.first.bounds.second<<";"<< elem.second << std::endl;
            pathPointCount += elem.first.bounds.second - elem.first.bounds.first;
  
    }
    //std::cout << "PointPt: " << pathPointCount << std::endl;
   
    //std::cout << "PointPt: " << pathPointCount << std::endl;

    double ratio = 0;
    for (auto elem : unassignedPts)
        ratio += elem.second /trajs[elem.first].pts.size();
       
    double objective = c1*pathFreq.size() + c2* ratio + c3 * dist;

    std::cout << "Obj;"<<objective << std::endl;


    return std::make_tuple(pthAssignments, pthStats, unassignedPts);
}


