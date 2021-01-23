#include "ycmatrix.h"
//using namespace YCutility;

std::vector<double> YCutility::matrix_m_matrix(const std::vector<double>& matrix_a, const std::vector<double>& matrix_b) {
	std::vector<double> result;
	for (int i = 0; i < 4; i += 1) {
		for (int j = 0; j < 4; j += 1) {
			double tempt = 0.f;
			for (int k = 0; k < 4; k += 1) {
				tempt = tempt + matrix_a[i * 4 + k] * matrix_b[k * 4 + j];
			}
			result.push_back(tempt);
		}
	}
	return result;
}

std::vector<double> YCutility::matrix_mutiple(const std::vector<double>& affine_m, const std::vector<double>& coord) {
	std::vector<double> result;
	for (int i = 0; i < 4; i += 1) {
		double tempt_sum = 0;
		for (int j = 0; j < 4; j += 1) {
			tempt_sum = tempt_sum + coord[j] * affine_m[i * 4 + j];
		}
		result.push_back(tempt_sum);
	}
	return result;
}


std::vector<double> YCutility::a_rotate_degree(const std::vector<double>& coordinate, int degree, double length) {
	std::vector<double> transformed_coordinate;
	//affine matrix Z rotation 
	std::vector<double> affine_tran_a = { 1,0,0, length / 2,
							   0,1,0,0,
							   0,0,1,0,
							   0,0,0,1 };

	std::vector<double> affineZ_a = { cos(degree * 3.14159 / 180),-sin(degree * 3.14159 / 180),0,0,
							   sin(degree * 3.14159 / 180), cos(degree * 3.14159 / 180),0, 0 ,
							  0, 0, 0, 0 ,
							   0, 0, 0, 1 };

	std::vector<double> affineX = { 1,0,0,0,
						   0, cos(20 * 3.14159 / 180),-sin(20 * 3.14159 / 180), 0 ,
						  0, sin(20 * 3.14159 / 180), cos(20 * 3.14159 / 180), 0 ,
						   0, 0, 0, 1 };

	std::vector<double> affineY = { cos(45 * 3.14159 / 180),0,sin(45 * 3.14159 / 180),0,
								 0, 1, 0, 0 ,
								 -sin(45 * 3.14159 / 180), 0,cos(45 * 3.14159 / 180), 0 ,
								 0, 0, 0, 1 };

	std::vector<double> tempt_affine;
	std::vector<double> tempt_affine2;

	tempt_affine = matrix_m_matrix(affineZ_a, affine_tran_a);
	//tempt_affine2 = matrix_m_matrix(tempt_affine, affineY);

	transformed_coordinate = matrix_mutiple(tempt_affine, coordinate);
	return transformed_coordinate;
}

std::vector<double> YCutility::b_rotate_degree(const std::vector<double>& coordinate, int degree_a, int degree_b, double length_a, double length_b) {
	std::vector<double> transformed_coordinate;
	//affine matrix Z rotation 
	std::vector<double> affine_tran_a = { 1,0,0, length_a,
							   0,1,0,0,
							   0,0,1,0,
							   0,0,0,1 };

	std::vector<double> affine_tran_b = { 1,0,0, length_b / 2,
						   0,1,0,0,
						   0,0,1,0,
						   0,0,0,1 };

	std::vector<double> affineX = { 1,0,0,0,
							   0, cos(20 * 3.14159 / 180),-sin(20 * 3.14159 / 180), 0 ,
							  0, sin(20 * 3.14159 / 180), cos(20 * 3.14159 / 180), 0 ,
							   0, 0, 0, 1 };

	std::vector<double> affineY = { cos(45 * 3.14159 / 180),0,sin(45 * 3.14159 / 180),0,
							 0, 1, 0, 0 ,
							 -sin(45 * 3.14159 / 180), 0,cos(45 * 3.14159 / 180), 0 ,
							 0, 0, 0, 1 };

	std::vector<double> affineZ_a = { cos(degree_a * 3.14159 / 180),-sin(degree_a * 3.14159 / 180),0,0,
							   sin(degree_a * 3.14159 / 180), cos(degree_a * 3.14159 / 180),0, 0 ,
							  0, 0, 0, 0 ,
							   0, 0, 0, 1 };

	std::vector<double> affineZ_b = { cos(degree_b * 3.14159 / 180),-sin(degree_b * 3.14159 / 180),0,0,
						   sin(degree_b * 3.14159 / 180), cos(degree_b * 3.14159 / 180),0, 0 ,
						  0, 0, 0, 0 ,
						   0, 0, 0, 1 };

	std::vector<double> tempt_affine;
	std::vector<double> tempt_affine2;
	std::vector<double> tempt_affine3;

	//std::vector<double> tempt_affine4;

	tempt_affine = matrix_m_matrix(affineZ_a, affine_tran_a);  // R2(theta1) * T (l1x)
	tempt_affine2 = matrix_m_matrix(tempt_affine, affineZ_b); // R2(theta1) * T (l1x) * R(theta2)
	tempt_affine3 = matrix_m_matrix(tempt_affine2, affine_tran_b);  // R2(theta1) * T (l1x) * R(theta2) * T(l2x/2)

	//tempt_affine4 = matrix_m_matrix(tempt_affine3, affineY);

	transformed_coordinate = matrix_mutiple(tempt_affine3, coordinate);
	return transformed_coordinate;
}

void YCutility::Copy_Mat4(double m1[16], const std::vector<double> m0)
{
	for (int i = 0; i < 16; ++i) { m1[i] = m0[i]; }
}

void YCutility::MatMat4(
	double m01[16],
	const double m0[16],
	const std::vector<double> m1)
{
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			m01[i * 4 + j] = m0[i * 4 + 0] * m1[0 * 4 + j] + m0[i * 4 + 1] * m1[1 * 4 + j] + m0[i * 4 + 2] * m1[2 * 4 + j] + m0[i * 4 + 3] * m1[3 * 4 + j];
		}
	}
}

template <typename REAL>
std::vector<double> YCutility::Quat(const REAL* q)
{
	std::vector<double> output;
	REAL x2 = q[1] * q[1] * 2.0;
	REAL y2 = q[2] * q[2] * 2.0;
	REAL z2 = q[3] * q[3] * 2.0;
	REAL xy = q[1] * q[2] * 2.0;
	REAL yz = q[2] * q[3] * 2.0;
	REAL zx = q[3] * q[1] * 2.0;
	REAL xw = q[1] * q[0] * 2.0;
	REAL yw = q[2] * q[0] * 2.0;
	REAL zw = q[3] * q[0] * 2.0;

	output.push_back(1.0 - y2 - z2);
	output.push_back(xy - zw);
	output.push_back(zx + yw);
	output.push_back(.0);

	output.push_back(xy + zw);
	output.push_back(1.0 - z2 - x2);
	output.push_back(yz - xw);
	output.push_back(.0);

	output.push_back(zx - yw);
	output.push_back(yz + xw);
	output.push_back(1.0 - x2 - y2);
	output.push_back(.0);

	output.push_back(.0);
	output.push_back(.0);
	output.push_back(.0);
	output.push_back(1.0);


	//m.SetZero();
	//m.mat[0 * 4 + 0] = 1.0 - y2 - z2; m.mat[0 * 4 + 1] = xy - zw;         m.mat[0 * 4 + 2] = zx + yw;
	//m.mat[1 * 4 + 0] = xy + zw;       m.mat[1 * 4 + 1] = 1.0 - z2 - x2;   m.mat[1 * 4 + 2] = yz - xw;
	//m.mat[2 * 4 + 0] = zx - yw;       m.mat[2 * 4 + 1] = yz + xw;         m.mat[2 * 4 + 2] = 1.0 - x2 - y2;
	//m.mat[3 * 4 + 3] = 1.0;
	return output;
}


template <typename REAL>
void YCutility::QuatQuat(
	REAL r[],
	const REAL p[],
	const REAL q[])
{
	r[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3];
	r[1] = p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2];
	r[2] = p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1];
	r[3] = p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0];
}

template <typename REAL>
void YCutility::Copy_Quat(
	REAL r[],
	const REAL p[])
{
	r[0] = p[0];
	r[1] = p[1];
	r[2] = p[2];
	r[3] = p[3];
}


