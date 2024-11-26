#include <iostream>
#include <cmath>
#include <algorithm>

struct Vec3 {
    double x, y, z;
};

// Helper function to find intersection range
bool intersectSlab(double p0, double d, double min, double max, double& tmin, double& tmax) {
    if (std::abs(d) < 1e-8) {
        // Line is parallel to the slab
        return p0 >= min && p0 <= max;
    }
    double t1 = (min - p0) / d;
    double t2 = (max - p0) / d;
    if (t1 > t2) std::swap(t1, t2);
    tmin = std::max(tmin, t1);
    tmax = std::min(tmax, t2);
    return tmin <= tmax;
}

bool doesLineIntersectBox(const Vec3& p0, const Vec3& d, const Vec3& boxMin, const Vec3& boxMax) {
    double tmin = -INFINITY, tmax = INFINITY;

    // Check x-axis slab
    if (!intersectSlab(p0.x, d.x, boxMin.x, boxMax.x, tmin, tmax)) return false;

    // Check y-axis slab
    if (!intersectSlab(p0.y, d.y, boxMin.y, boxMax.y, tmin, tmax)) return false;

    // Check z-axis slab
    if (!intersectSlab(p0.z, d.z, boxMin.z, boxMax.z, tmin, tmax)) return false;

    return true;
}

int main() {
    Vec3 p0 = {2, 2, 1}; // Line point
    Vec3 d = {-1,-1,-1};  // Line direction
    Vec3 boxMin = {-1, -1, -1}; // Box min corner
    Vec3 boxMax = {1, 1, 1};    // Box max corner

    if (doesLineIntersectBox(p0, d, boxMin, boxMax)) {
        std::cout << "The line intersects the box." << std::endl;
    } else {
        std::cout << "The line does not intersect the box." << std::endl;
    }

    return 0;
}
