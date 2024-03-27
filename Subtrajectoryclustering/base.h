#pragma once
#include <functional>
#include <iostream>
#include <vector>
#include <tuple>
#ifndef PT_H
#define PT_H
class Pt
{
    /* A class representing a point in a trajectory.

        Attributes:
        lat(float) : latitude.
        lon(float) : longitude.
        trajID(int) : ID of the point's trajectory.
        t(float) : timestamp associated with the point.
        */

public:
    double lat;
    double lon;
    int trajID;
    double t;
    //Default constructor with dummy initialization values. 
    Pt();
    Pt(double lat, double lon, int trajID, double t);
    //Computes a hash so that pt objects can be used as dict keys.
    struct Hash {
        size_t operator()(const Pt& point) const;
    };
    //Define == operator for pt objects.
    bool operator==(const Pt& other) const;
    //Return string to be output while printing a pt object. 
    friend std::ostream& operator<<(std::ostream& os, const Pt& point);
    bool operator<(const Pt& other) const;
};
#endif

class Traj
{
    /*
     A class representing a trajectory. 
    
    Attributes:
        pts: list of points in the trajectory.
        trajID: int that represents the ID of the trajectory
    */
public:
    std::vector<Pt> pts;
    int trajID;
    //Initialize trajectory with empty list of points.
    Traj();
    /*
     Add a pt to the trajectory.
        
        Args:
            lat (float): latitude of point.
            lon (float): longitude of point.
            trajID (int): trajID of the point (all points of a traj will have same ID).
    */
    void addPt(double lat, double lon, int trajID, double t);
    //Sort points of trajectory in ascending order of timestamp.
    void sortPts();
};

class Pathlet {
    /*
    A class representing a pathlet (essentially a subtrajectory).
    
        Similar to traj, however avoids storing the points explicitly by referring to the
        trajectory to which the points belong, along with the starting and ending indices
        of the points in the list of points of the trajectory, sorted by timestamp.
        
        Attributes:
            trajID (int): ID of the trajectory to which points of the pathlet belong.
            bounds (int,int): start and end indices of points in the timestamp-sorted list of
                              trajectory points.
    */
public:
    int trajID;
    std::pair<int, int> bounds;
    Pathlet();
    /*
    Initialize pathlet with points from a trajectory.
        
            Args:
                trajID (int): ID of trajectory.
                bounds (int,int): start and end indices.
    */
    Pathlet(int trajID, std::pair<int, int> bounds);
    // Define == operator for pathlet objects. 
    bool operator==(const Pathlet& other) const;
    //Define a hash function so that pathlets can be used as keys in a dict. 
    struct Hash {
        std::size_t operator()(const Pathlet& path) const;
    };
    //Return string to be output while printing a pathlet object.
    friend std::ostream& operator<<(std::ostream& os, const Pathlet& pathlet);

};


class SubTraj {
    /*
    A class representing a subtrajectory.
    
    Exactly identical to a pathlet class. Defined as a class of its own for conceptual reasons.
    
    */
public:
    int trajID;
    std::pair<int, int> bounds;
    SubTraj();
    SubTraj(int trajID, std::pair<int, int> bounds);
    bool operator==(const SubTraj& other) const;
    
    struct Hash{
        std::size_t operator()(const SubTraj& st) const;
    };
   

    friend std::ostream& operator<<(std::ostream& os, const SubTraj& subTraj);
};

namespace std {
    template<>
    struct hash<Pt>
    {
        size_t operator()(const Pt& point) const {
            std::size_t h1 = std::hash<double>()(point.lat);
            std::size_t h2 = std::hash<double>()(point.lon);
            std::size_t h3 = std::hash<int>()(point.trajID);
            std::size_t h4 = std::hash<double>()(point.t);
            return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3);
        }
    };
}
namespace std {
    template <>
    struct hash<Pathlet> {
        std::size_t operator()(const Pathlet& path) const {
            std::size_t hash = std::hash<int>()(path.trajID);
            std::size_t h1 = std::hash<int>()(path.bounds.first);
            std::size_t h2 = std::hash<int>()(path.bounds.second);
            return hash ^ ( h1 << 1 ) ^ ( h2 <<2 );
        }
    };
}

namespace std {
    template <>
    struct hash<SubTraj>
    {
        std::size_t operator()(const SubTraj& st) const {
            std::size_t h1 = std::hash<int>()(st.trajID);
            std::size_t h2 = std::hash<int>()(st.bounds.first);
            std::size_t h3 = std::hash<int>()(st.bounds.second);
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }

    };
}


namespace std {
    template <>
    struct hash<std::pair<Pathlet, SubTraj>> {
        std::size_t operator()(const std::pair<Pathlet, SubTraj>& pair) const {
            std::size_t h1 = std::hash<Pathlet>()(pair.first);
            std::size_t h2 = std::hash<SubTraj>()(pair.second);
            return h1 ^ (h2 << 1);
        }
    };
}

namespace std {
    template <>
    struct hash<std::pair<Pathlet, int>> {
        std::size_t operator()(const std::pair<Pathlet, int>& pair) const {
            std::size_t h1 = std::hash<Pathlet>()(pair.first);
            std::size_t h2 = std::hash<int>()(pair.second);
            return h1 ^ (h2 << 1);
        }
    };
}

