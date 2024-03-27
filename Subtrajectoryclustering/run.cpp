#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>
#include <chrono>
#include "base.h"
#include "ioUtils.h"
#include "distance.h"
#include "greedy.h"
#include "heapdict.h"
#include "frechet.h"
#include <random>

int main(int argc, char* argv[]) {


    //phase1: loading trajectories:
    std::string fileName = argv[1];
    int simplify = std::stoi(argv[2]);

    auto start_time = std::chrono::high_resolution_clock::now();
    auto start_time_total = std::chrono::steady_clock::now();

	
    std::unordered_map<int, Traj> tr = readTrajsFromTxtFile(fileName);

 

    auto end_time = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //std::cout << "Loading trajectories execution time: " << duration.count() << " milliseconds\n" << std::endl;
    std::cout << "Phase;1" << std::endl;

    //phase 2: computing frechet distance:
    start_time = std::chrono::high_resolution_clock::now();
	std::unordered_map<std::pair<Pathlet, SubTraj>, double> proc = process(tr,0.0005,0.002, simplify);

    std::unordered_map<std::pair<Pathlet, int>, std::vector<std::pair<SubTraj, double>>> distPairs2;

    // Iterate through distPairs1 and populate distPairs2
    for (const auto& entry : proc) {
        Pathlet pth = entry.first.first;
        int trID = entry.first.second.trajID;
        double dist = entry.second;
        SubTraj straj = entry.first.second;

        std::pair<Pathlet, int> key(pth, trID);
        if (distPairs2.find(key) != distPairs2.end()) {
            distPairs2[key].emplace_back(straj, dist);
        }
        else {
            distPairs2[key] = { {straj, dist} };
        }
    }
    
   



    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //std::cout << "Computing Frechet distances execution time: " << duration.count() << " milliseconds\n" << std::endl;
    std::cout << "Phase;2" << std::endl;
    
    //phase 3: compute prerequisite data structures

    start_time = std::chrono::high_resolution_clock::now();
    HeapDict<Pathlet, int> hp;

    std::unordered_map<SubTraj, int, SubTraj::Hash > strajCov;
    std::unordered_map < Pt, std::unordered_set<SubTraj, SubTraj::Hash>,Pt::Hash >  ptStraj;
    std::unordered_map<SubTraj, std::vector<Pathlet>, SubTraj::Hash> strajPth;
    std::unordered_map<int, int > trajCov;
    
    auto tup = preprocessGreedy(tr, distPairs2);

    
    strajCov = std::get<0>(tup);
    ptStraj = std::get<1>(tup);

    strajPth = std::get<2>(tup);
    trajCov = std::get<3>(tup);
    float c1, c2 , c3 ;
    c1 = std::stof(argv[3]);
    c2 = std::stof(argv[4]);
    c3 = std::stof(argv[5]);

    std::cout << c1 << std::endl;

    
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //std::cout << "Computing prerequisite data structures execution time: " << duration.count() << " milliseconds\n" << std::endl;
    std::cout << "Phase;3" << std::endl;


    //phase 4: running greedy algorithm

    start_time = std::chrono::high_resolution_clock::now();
    std::cout << "Phase;4" << std::endl;
    auto retval = runGreedy(tr, distPairs2, strajCov, ptStraj, strajPth, trajCov, c1, c2, c3);

    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);



   // std::cout << "Running greedy algorithm execution time:  " << duration.count() << " milliseconds" << std::endl;
    auto end_time_total = std::chrono::steady_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::minutes>(end_time_total - start_time_total);


    std::cout << "RunningTime;" << duration2.count()<< std::endl;


    return 0;
}




