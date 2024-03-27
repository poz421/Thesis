#include "base.h"
#include <algorithm>



    Pt::Pt() : lat(0.0), lon(0.0), trajID(-1), t(-1.0) {}
    Pt::Pt(double lat, double lon, int trajID, double t)
        : lat(lat), lon(lon), trajID(trajID), t(t) {}

    // Define a custom hash function for Pt objects
    //so that pt objects can be used as dict keys
    size_t Pt::Hash::operator()(const Pt& point) const {
            return std::hash<double>()(point.lat) ^
                std::hash<double>()(point.lon) ^
                std::hash<int>()(point.trajID) ^
                std::hash<double>()(point.t);
    }

    // Define a custom equality operator for Pt objects
    bool Pt::operator==(const Pt& other) const {
        return (lat == other.lat && lon == other.lon && trajID == other.trajID && t == other.t);
    }

    // Define a custom string representation for Pt objects
    std::ostream& operator<<(std::ostream& os, const Pt& point) {
        os << "Point TrajID " << point.trajID << " ; lat-long (" << point.lat << "," << point.lon << "); time " << point.t;
        return os;
    }
    bool Pt::operator<(const Pt& other) const {
        return t < other.t;
    }

    Traj::Traj()
    {

    }
    void Traj::addPt(double lat, double lon, int trajID, double t)
    {
        this->trajID = trajID;
        Pt p(lat, lon, trajID, t);
        pts.push_back(p);
    }

    void Traj::sortPts() {
        std::sort(pts.begin(), pts.end());
    }



    Pathlet::Pathlet()
    {
        trajID = 0;
    }


    Pathlet::Pathlet(int trajID, std::pair<int, int> bounds)
        : trajID(trajID), bounds(bounds) {}

    bool Pathlet::operator==(const Pathlet& other) const {
        return (trajID == other.trajID && bounds == other.bounds);
    }

   
    std::size_t Pathlet::Hash::operator()(const Pathlet& path) const {
         std::size_t hash = std::hash<int>()(path.trajID);
         hash ^= (std::hash<int>()(path.bounds.first)<<1) ;
         hash ^= (std::hash<int>()(path.bounds.second)<<2);
         return hash;
    }
    

    std::ostream& operator<<(std::ostream& os, const Pathlet& pathlet) {
        os << "Pathlet TrajID " << pathlet.trajID << " ; bounds (" << pathlet.bounds.first << ", " << pathlet.bounds.second << ")";
        return os;
    }


    SubTraj::SubTraj()
    {
        trajID = -1;
    }


    SubTraj::SubTraj(int trajID, std::pair<int, int> bounds)
        : trajID(trajID), bounds(bounds) {}

    bool SubTraj::operator==(const SubTraj& other) const {
        return (trajID == other.trajID && bounds == other.bounds);
    }

   
    std::size_t SubTraj::Hash::operator()(const SubTraj& st) const {
         std::size_t hash = std::hash<int>()(st.trajID);
         hash ^= (std::hash<int>()(st.bounds.first)<<1) ;
         hash ^= (std::hash<int>()(st.bounds.second)<<2);
         return hash;
        }
   

    std::ostream& operator<<(std::ostream& os, const SubTraj& subTraj) {
        os << "Subtraj TrajID " << subTraj.trajID << " ; bounds (" << subTraj.bounds.first << ", " << subTraj.bounds.second << ")";
        return os;
    }





