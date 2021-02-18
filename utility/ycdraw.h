#include <vector>


namespace YCdraw {
	void DrawSphere(int nla, int nlo);
	void DrawTrajectory(std::deque<double> trajectory);
	void DrawSphereAt(int nla, int nlo, double rad, double x, double y, double z);
	void DrawTrajectory_points(std::deque<double> trajectory);
}

