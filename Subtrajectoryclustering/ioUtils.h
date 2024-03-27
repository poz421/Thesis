#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "base.h"
/*
Read trajectory points from a text file.

    Each line of the file contains one point in the following format
    (angled braces provided for clarity).
    <trajID>;<objID>;<timestamp>;<lat>;<lon>
    objID is ignored.

    Args:
        fName (str): name of the file.

    Returns:
        A dictionary of trajectories of the form {trID : traj}.
        Points within a trajectory are sorted by timestamp.

*/
std::unordered_map<int, Traj> readTrajsFromTxtFile(const std::string& fName);