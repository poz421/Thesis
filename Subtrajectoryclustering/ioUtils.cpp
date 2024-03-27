#include "ioUtils.h"


std::unordered_map<int, Traj> readTrajsFromTxtFile(const std::string& fName) {
    std::unordered_map<int, Traj> trajs;

    std::ifstream file(fName);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open the file." << std::endl;
        return trajs;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.find("#") != std::string::npos) {
            continue; // Skip lines containing #
        }

        std::replace(line.begin(), line.end(), '\n', ';');
        std::istringstream linestream(line);
        std::string token;
        std::stringstream ss(line);
        int trajID;
        double timestamp, lat, lon;
        
        
        if (std::getline(ss, token, ';')) {
            trajID = std::stoi(token);
        }
        
        if (std::getline(ss, token, ';')) {
            // Ignore objID
        }

        if (std::getline(ss, token, ';')) {
            timestamp = std::stod(token);
        }
        if (std::getline(ss, token, ';')) {
            lat = std::stod(token);
       
        }
        if (std::getline(ss, token, ';')) {
            lon = std::stod(token);
        }

        if (trajs.find(trajID) == trajs.end()) {
            trajs[trajID] = Traj();
        }

       
        trajs[trajID].addPt(lat, lon, trajID, timestamp);
    }

    file.close();

    for (auto& trajPair : trajs) {
        trajPair.second.sortPts();
    }

    return trajs;
}