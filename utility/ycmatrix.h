#include <vector>


namespace YCutility {

	std::vector<double> matrix_m_matrix(const std::vector<double>& matrix_a, const std::vector<double>& matrix_b);

	std::vector<double> matrix_mutiple(const std::vector<double>& affine_m, const std::vector<double>& coord);

	std::vector<double> a_rotate_degree(const std::vector<double>& coordinate, int degree, double length);

	std::vector<double> b_rotate_degree(const std::vector<double>& coordinate, int degree_a, int degree_b, double length_a, double length_b);

	void Copy_Mat4(double m1[16], const std::vector<double> m0);

	void MatMat4(double m01[16],const double m0[16],const std::vector<double> m1);

	template <typename REAL> std::vector<double> Quat(const REAL* q);

	template <typename REAL> void QuatQuat(REAL r[],const REAL p[],const REAL q[]);

	template <typename REAL> void Copy_Quat(REAL r[],const REAL p[]);
}

