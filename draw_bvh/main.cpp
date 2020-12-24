#include <GLFW/glfw3.h>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <stdlib.h>
//#include "ycmatrix.h"


class RIGbone
{
public:
	RIGbone() {
		for (int i = 0; i < 16; ++i) { invBindMat[i] = 0.0; }
		invBindMat[0] = 1.0;
		invBindMat[5] = 1.0;
		invBindMat[10] = 1.0;
		invBindMat[15] = 1.0;
		//
		scale = 1;
		quatRelativeRot[0] = 1;
		quatRelativeRot[1] = 0;
		quatRelativeRot[2] = 0;
		quatRelativeRot[3] = 0;
		transRelative[0] = 0;
		transRelative[1] = 0;
		transRelative[2] = 0;
		ibone_parent = -1;
	}
	std::vector<double> Pos() const {
		std::vector<double> result = {affmat3Global[3], affmat3Global[7], affmat3Global[11]};
		return result;
	}
	void DeformSkin(double pos2[3],
		const double pos0[3]) const;
	void SetRotationBryant(double rx, double ry, double rz);
	void SetTranslation(double tx, double ty, double tz);
	/*
	int PickHandler(const delfem2::CVec3d& org,
					const delfem2::CVec3d& dir,
					double rad_handlr,
					double tol) const;
	 */
	void AffineJoint(const double a[16]) const;


public:
	std::string bname; 

	/**
	 * @details Inverse of Affine matrix to send the skin to the bone reference config. The joint position of this bone will be mapped to the origin
	 */

	//use for GLTF
	double invBindMat[16];

	int ibone_parent; // initialized and stay constant

	double transRelative[3]; // position of the joint position from parent joint

	double scale; 

	/**
	 * @brief rotation of the joint from parent joint (quaternion w,x,y,z).
	 * @details this value will be changed  when pose is edited
	 */
	double quatRelativeRot[4];

	/**
	 * @brief affine matrix  to send bone from the origin to the deformed pose
	 * @details this value will be set when the pose is edited using the function
	 * "void UpdateBoneRotTrans(std::vector<CRigBone>& aBone)"
	 */
	double affmat3Global[16];
};

class CChannel_bvh
{
public:
	CChannel_bvh(int ib, int ia, bool br) {
		this->ibone = ib;
		this->iaxis = ia;
		this->isrot = br;
	}
public:
	int ibone;
	int iaxis;
	bool isrot;
};

//what is this for ?
void CalcInvMat(double* a, const int n, int& info)
{
	double tmp1;

	info = 0;
	int i, j, k;
	for (i = 0; i < n; i++) {
		if (fabs(a[i * n + i]) < 1.0e-30) {
			info = 1;
			return;
		}
		if (a[i * n + i] < 0.0) {
			info--;
		}
		tmp1 = 1.0 / a[i * n + i];
		a[i * n + i] = 1.0;
		for (k = 0; k < n; k++) {
			a[i * n + k] *= tmp1;
		}
		for (j = 0; j < n; j++) {
			if (j != i) {
				tmp1 = a[j * n + i];
				a[j * n + i] = 0.0;
				for (k = 0; k < n; k++) {
					a[j * n + k] -= tmp1 * a[i * n + k];
				}
			}
		}
	}
}

std::vector<double> matrix_m_matrix(const std::vector<double>& matrix_a, const std::vector<double>& matrix_b) {
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
//[1]Utilities functions including copy matrix, copy matrix x matrix to a 4x4 matrix and 4x4 matrix return a 4x4 matrix---------------

void Copy_Mat4(double m1[16], const std::vector<double> m0)
{
	for (int i = 0; i < 16; ++i) { m1[i] = m0[i]; }
}

void MatMat4(
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
std::vector<double> Quat(const REAL* q)
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

void UpdateBoneRotTrans
(std::vector<RIGbone>& aBone)
{
	for (std::size_t ibone = 0; ibone < aBone.size(); ++ibone) {

		std::vector<double> tempt_tranform4d;
		for (int i = 0; i < 16; i++) {
			if ((i == 0) || (i == 5) || (i == 10) || (i == 15)) {
				tempt_tranform4d.push_back(1.0f);
			}
			else if (i == 3) {
				tempt_tranform4d.push_back(aBone[ibone].transRelative[0]);
			}
			else if (i == 7) {
				tempt_tranform4d.push_back(aBone[ibone].transRelative[1]);
			}
			else if (i == 11) {
				tempt_tranform4d.push_back(aBone[ibone].transRelative[2]);
			}
			else {
				tempt_tranform4d.push_back(0.0f);
			}
		}
		double X_degree = aBone[ibone].quatRelativeRot[1];
		double Y_degree = aBone[ibone].quatRelativeRot[2];
		double Z_degree = aBone[ibone].quatRelativeRot[3];

		std::vector<double> tempt_Xrotate_4d = { 1,0,0,0,
					   0, cos(X_degree * 3.14159 / 180),-sin(X_degree * 3.14159 / 180), 0 ,
					  0, sin(X_degree * 3.14159 / 180), cos(X_degree * 3.14159 / 180), 0 ,
					   0, 0, 0, 1 };

		std::vector<double> tempt_Yrotate_4d = { cos(Y_degree * 3.14159 / 180),0,sin(Y_degree * 3.14159 / 180),0,
							 0, 1, 0, 0 ,
							 -sin(Y_degree * 3.14159 / 180), 0,cos(Y_degree * 3.14159 / 180), 0 ,
							 0, 0, 0, 1 };

		std::vector<double> tempt_Zrotate_4d = { cos(Z_degree * 3.14159 / 180),-sin(Z_degree * 3.14159 / 180),0,0,
						   sin(Z_degree * 3.14159 / 180), cos(Z_degree * 3.14159 / 180),0, 0 ,
						  0, 0, 0, 0 ,
						   0, 0, 0, 1 };

		std::vector<double> temp_result;
		std::vector<double> temp_result2;
		std::vector<double> temp_result3;



		std::vector<double> quat_matrix = Quat(aBone[ibone].quatRelativeRot);

		temp_result = matrix_m_matrix(tempt_tranform4d, tempt_Xrotate_4d);
		//temp_result2 = matrix_m_matrix(temp_result, tempt_Yrotate_4d);
		//temp_result3 = matrix_m_matrix(temp_result2, tempt_Zrotate_4d);
		temp_result3 = matrix_m_matrix(temp_result, quat_matrix);

		// quertion way, and use GLM matrix 4D
		//CMat4d m01 = CMat4d::Translate(aBone[ibone].transRelative);
		//m01 = m01 * CMat4d::Quat(aBone[ibone].quatRelativeRot);
		//m01 = m01 * CMat4d::Scale(aBone[ibone].scale);


		const int ibone_p = aBone[ibone].ibone_parent;
		if (ibone_p < 0 || ibone_p >= (int)aBone.size()) { // root bone
			Copy_Mat4(aBone[ibone].affmat3Global, temp_result3);
			continue;
		}
		MatMat4(aBone[ibone].affmat3Global,
			aBone[ibone_p].affmat3Global, temp_result3);
	}
}

template <typename REAL>
void QuatQuat(
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
void Copy_Quat(
	REAL r[],
	const REAL p[])
{
	r[0] = p[0];
	r[1] = p[1];
	r[2] = p[2];
	r[3] = p[3];
}

void SetPose_BioVisionHierarchy
(std::vector<RIGbone>& aBone,
	const std::vector<CChannel_bvh>& aChannelRotTransBone,
	const double* aVal)
{
	for (auto& bone : aBone) {
		bone.quatRelativeRot[0] = 1.0;
		bone.quatRelativeRot[1] = 0.0;
		bone.quatRelativeRot[2] = 0.0;
		bone.quatRelativeRot[3] = 0.0;
	}
	const int nch = aChannelRotTransBone.size();
	for (int ich = 0; ich < nch; ++ich) {
		const int ibone = aChannelRotTransBone[ich].ibone;
		const int iaxis = aChannelRotTransBone[ich].iaxis;
		const bool isrot = aChannelRotTransBone[ich].isrot;
		const double val = aVal[ich];
		assert(ibone < (int)aBone.size());
		assert(iaxis >= 0 && iaxis < 3);
		if (!isrot) {
			aBone[ibone].transRelative[iaxis] = val;


			//old method, degree transformation ==========
			//this is for transformation 
			//if((iaxis == 0) || (iaxis == 1))
				// move the character to the center a little bit 
			//	aBone[ibone].transRelative[iaxis] = val;
			//else
			//	aBone[ibone].transRelative[iaxis] = val;
		}
		else {
			// this is for rotation  ==========
			//my way store degree in quatRelativeRot instead of quenteron 
			//aBone[ibone].quatRelativeRot[1 + iaxis] = val;

			// the way of quertino method 
			const double ar = val * 3.1415926 / 180.0;
			double v0[3] = { 0,0,0 };
			v0[iaxis] = 1.0;
			double dq[4] = { cos(ar * 0.5), v0[0] * sin(ar * 0.5), v0[1] * sin(ar * 0.5), v0[2] * sin(ar * 0.5) };
			double qtmp[4]; QuatQuat(qtmp,
				aBone[ibone].quatRelativeRot, dq);
			Copy_Quat(aBone[ibone].quatRelativeRot, qtmp);

		}
	}
	UpdateBoneRotTrans(aBone);
}



//[1]ABOVE ARE UTILITY FUNCTIONS -------------------------

//[2]load BVH file main program
bool loadBVH(std::vector<RIGbone>& aBone,
	std::vector<CChannel_bvh>& aChannelRotTransBone,
	int& nframe,
	std::vector<double>& aValueRotTransBone,
	const std::string& path_bvh)
{
	std::ifstream fin;
	fin.open(path_bvh.c_str());
	if (!fin.is_open()) {
		std::cout << "cannot open file" << std::endl;
		return false;
	}
	aBone.clear();
	aChannelRotTransBone.clear();
	//
	std::string line;
	std::vector<int> stackIndBone;
	while (std::getline(fin, line))
	{
		if (line[line.size() - 1] == '\n') line.erase(line.size() - 1); // remove the newline code
		if (line[line.size() - 1] == '\r') line.erase(line.size() - 1); // remove the newline code

		std::vector<std::string> aToken;
		std::istringstream iss(line);
		for (std::string line; iss >> line; )
			aToken.push_back(line);

		//for (int i = 0; i < aToken.size(); i += 1) {
		//	std::cout << aToken[i];
		//}
		//std::cout << std::endl;

		if (aToken[0] == "HIERARCHY") {
			assert(aBone.empty());
		}
		else if (aToken[0] == "ROOT") {
			assert(aBone.size() == 0);
			RIGbone br;
			assert(aToken.size() == 2);
			br.bname = aToken[1];
			aBone.push_back(br);
		}
		else if (aToken[0] == "{") {
			stackIndBone.push_back(aBone.size() - 1);
			if (stackIndBone.size() > 1) {
				int ibp = stackIndBone[stackIndBone.size() - 2];
				int ib = aBone.size() - 1;
				aBone[ib].ibone_parent = ibp;
			}
		}
		else if (aToken[0] == "}") {
			stackIndBone.resize(stackIndBone.size() - 1);
		}
		else if (aToken[0] == "OFFSET") {
			assert(aToken.size() == 4);
			int ib = aBone.size() - 1;
			double org_x = atof(aToken[1].c_str());
			double org_y = atof(aToken[2].c_str());
			double org_z = atof(aToken[3].c_str());
			aBone[ib].invBindMat[3] = -org_x;
			aBone[ib].invBindMat[7] = -org_y;
			aBone[ib].invBindMat[11] = -org_z;
			if (stackIndBone.size() > 1) {
				const int ibp = stackIndBone[stackIndBone.size() - 2];
				assert(ibp < (int)aBone.size());
				aBone[ib].invBindMat[3] += aBone[ibp].invBindMat[3];
				aBone[ib].invBindMat[7] += aBone[ibp].invBindMat[7];
				aBone[ib].invBindMat[11] += aBone[ibp].invBindMat[11];
			}
		}
		else if (aToken[0] == "CHANNELS") {
			assert(aToken.size() >= 2);

			int nch = stoi(aToken[1]);
			assert((int)aToken.size() == nch + 2);
			assert(!aBone.empty());
			const std::size_t ib = aBone.size() - 1;
			for (int ich = 0; ich < nch; ++ich) {
				const std::string& type_ch = aToken[ich + 2];
				if (type_ch == "Xposition") { aChannelRotTransBone.emplace_back(ib, 0, false); }     // emplace_back just faster
				else if (type_ch == "Yposition") { aChannelRotTransBone.emplace_back(ib, 1, false); }
				else if (type_ch == "Zposition") { aChannelRotTransBone.emplace_back(ib, 2, false); }
				// profs' method make the order back to natural order XYZ
				else if (type_ch == "Xrotation") { aChannelRotTransBone.emplace_back(ib, 0, true); }
				else if (type_ch == "Yrotation") { aChannelRotTransBone.emplace_back(ib, 1, true); }
				else if (type_ch == "Zrotation") { aChannelRotTransBone.emplace_back(ib, 2, true); }
				else {
					std::cout << "ERROR-->undefiend type" << std::endl;
				}
			}
		}
		else if (aToken[0] == "JOINT") {
			RIGbone br;
			assert(aToken.size() == 2);
			br.bname = aToken[1];
			aBone.push_back(br);
		}
		else if (aToken[0] == "End") {
			assert(aToken[1] == "Site");
			RIGbone br;
			assert(aToken.size() == 2);
			br.bname = aToken[1];
			aBone.push_back(br);
		}
		else if (aToken[0] == "MOTION") {
			break;
		}
	}
	nframe = 0;
	//parse frame information 
	{
		std::string tempt_str0;
		{
			std::getline(fin, line);
			std::stringstream ss(line);
			ss >> tempt_str0 >> nframe;
			//      std::cout << "frame: " << nframe << std::endl;
		}
		std::getline(fin, line);
		//    std::cout << "frametime: " << line << std::endl;
	}
	const int nchannel = aChannelRotTransBone.size();
	aValueRotTransBone.resize(nframe* nchannel);
	for (int iframe = 0; iframe < nframe; ++iframe) {
		std::getline(fin, line);
		//line = rig_v3q::MyReplace(line, '\t', ' ');   Prof.Umetani replace tab to space
		if (line[line.size() - 1] == '\n') line.erase(line.size() - 1); // remove the newline code
		if (line[line.size() - 1] == '\r') line.erase(line.size() - 1); // remove the newline code
		
		//my split 
		std::stringstream ss(line);
		std::istream_iterator<std::string> begin(ss);
		std::istream_iterator<std::string> end;
		std::vector<std::string> aToken(begin, end);

		//print out the result
		//std::copy(aToken.begin(), aToken.end(), std::ostream_iterator<std::string>(std::cout, "\n"));

		// prof.umetani's split
		//std::vector<std::string> aToken = rig_v3q::MySplit(line, ' ');
		//std::cout << aToken.size() << " " << aChannelRotTransBone.size() << std::endl;


		assert(aToken.size() == aChannelRotTransBone.size());
		for (int ich = 0; ich < nchannel; ++ich) {
			//aValueRotTransBone[iframe * nchannel + ich] = rig_v3q::myStod(aToken[ich]);
			aValueRotTransBone[iframe * nchannel + ich] = atof(aToken[ich].c_str());
		}
	}
	// --------------- bones initilization ?
	for (std::size_t ibone = 0; ibone < aBone.size(); ++ibone) {
		RIGbone& bone = aBone[ibone];
		bone.scale = 1.0;
		bone.quatRelativeRot[0] = 1.0;
		bone.quatRelativeRot[1] = 0.0;
		bone.quatRelativeRot[2] = 0.0;
		bone.quatRelativeRot[3] = 0.0;
		bone.transRelative[0] = 0.0;
		bone.transRelative[1] = 0.0;
		bone.transRelative[2] = 0.0;
		// WHY ? why transrelative takes the minus substraction 
		if (bone.ibone_parent != -1) {
			const RIGbone& bone_p = aBone[bone.ibone_parent];
			bone.transRelative[0] = (-bone.invBindMat[3]) - (-bone_p.invBindMat[3]);
			bone.transRelative[1] = (-bone.invBindMat[7]) - (-bone_p.invBindMat[7]);
			bone.transRelative[2] = (-bone.invBindMat[11]) - (-bone_p.invBindMat[11]);
		}
	}
	for (auto& bone : aBone) {
		for (int i = 0; i < 16; ++i) {
			bone.affmat3Global[i] = bone.invBindMat[i];
		}
		int info; 
		//what is this for ?
		CalcInvMat(bone.affmat3Global, 4, info);
	}
	return true;
}
//[2]ABOVE is BVH parser

//[3]professor's draw function 
void DrawSphere
(int nla, int nlo)
{
	if (nla <= 1 || nlo <= 2) { return; }
	const double pi = 3.1415926535;
	double dla = 2.0 * pi / nla;
	double dlo = pi / nlo;
	::glBegin(GL_QUADS);
	for (int ila = 0; ila < nla; ila++) {
		for (int ilo = 0; ilo < nlo; ilo++) {
			double rla0 = (ila + 0) * dla;
			double rla1 = (ila + 1) * dla;
			double rlo0 = (ilo + 0) * dlo;
			double rlo1 = (ilo + 1) * dlo;
			double p0[3] = { cos(rla0) * cos(rlo0), cos(rla0) * sin(rlo0), sin(rla0) };
			double p1[3] = { cos(rla0) * cos(rlo1), cos(rla0) * sin(rlo1), sin(rla0) };
			double p2[3] = { cos(rla1) * cos(rlo1), cos(rla1) * sin(rlo1), sin(rla1) };
			double p3[3] = { cos(rla1) * cos(rlo0), cos(rla1) * sin(rlo0), sin(rla1) };
			::glVertex3dv(p0);
			::glVertex3dv(p1);
			::glVertex3dv(p2);
			::glVertex3dv(p3);
		}
	}
	::glEnd();
}

void DrawSphereAt
(int nla, int nlo, double rad, double x, double y, double z)
{
	glTranslated(+x, +y, +z);
	glScaled(rad, rad, rad);
	DrawSphere(nla, nlo);
	glScaled(1.0 / rad, 1.0 / rad, 1.0 / rad);
	glTranslated(-x, -y, -z);
}

void Draw_RigBone(
	int ibone,
	bool is_selected,
	int ielem_selected,
	const std::vector<RIGbone>& aBone,
	double rad_bone_sphere,
	double rad_rot_hndlr)
{
	{ // draw sphere joint point
		if (is_selected) {glColor3d(1, 0, 0); }
		else { glColor3d(1, 1, 1); }
		const std::vector<double> pos = aBone[ibone].Pos();
		// undertand how it works 
		DrawSphereAt(32, 32, rad_bone_sphere, pos[0], pos[1], pos[2]);
	}
	//if (is_selected) {
		// understand how rotate works and WHY  this one should rotate the sphere let me just temporarily ignore this 
		//opengl::DrawHandlerRotation_Mat4(aBone[ibone].affmat3Global, rad_rot_hndlr, ielem_selected);
		//int ibone_parent = aBone[ibone].ibone_parent;
		//if (ibone_parent >= 0 && ibone_parent < (int)aBone.size()) {
		//	const std::vector<double> pp(aBone[ibone_parent].Pos());
		//}
		//else {
		//}
	//}
}

void DrawBone(
	const std::vector<RIGbone>& aBone,
	int ibone_selected,
	int ielem_selected,
	double rad_bone_sphere,
	double rad_rot_hndlr)
{
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glPointSize(3);
	for (unsigned int iskel = 0; iskel < aBone.size(); ++iskel) {
		const bool is_selected = (int)iskel == ibone_selected;
		Draw_RigBone(iskel,
			is_selected, ielem_selected, aBone,
			rad_bone_sphere, rad_rot_hndlr);
	}
	// draw edges whilte
	for (unsigned int ibone = 0; ibone < aBone.size(); ++ibone) {
		const RIGbone& bone = aBone[ibone];
		const int ibone_p = aBone[ibone].ibone_parent;
		if (ibone_p < 0 || ibone_p >= (int)aBone.size()) { continue; }
		const RIGbone& bone_p = aBone[ibone_p];
		bool is_selected_p = (ibone_p == ibone_selected);
		if (is_selected_p) { glColor3d(1.0, 1.0, 1.0); }
		else { glColor3d(1.0, 0.0, 0.0); }
		glBegin(GL_LINES);
		// understand what is GLVertex
		std::vector<double> tempt_cur = bone.Pos();
		std::vector<double> tempt_par = bone_p.Pos();
		glVertex3f(tempt_cur[0], tempt_cur[1], tempt_cur[2]);
		glVertex3f(tempt_par[0], tempt_par[1], tempt_par[2]);
		//opengl::myGlVertex(bone.Pos());
		//opengl::myGlVertex(bone_p.Pos());
		glEnd();
	}
}

//[3]ABOVE are draw functions


void geneCUBE_length(double x, double y, double z, std::vector<double>& out_vert, std::vector<unsigned int>& out_tri){
	out_vert = {
	-x / 2.0f,-y / 2.0f,-z / 2.0f,  //0
	-x / 2.0f,-y / 2.0f, z / 2.0f,  //1
	-x / 2.0f, y / 2.0f, z / 2.0f,  //2
	 x / 2.0f, y / 2.0f,-z / 2.0f,  //3 
	-x / 2.0f, y / 2.0f,-z / 2.0f,  //4
	 x / 2.0f,-y / 2.0f, z / 2.0f,  //5
	 x / 2.0f,-y / 2.0f,-z / 2.0f,  //6
	 x / 2.0f, y / 2.0f, z / 2.0f   //7
	};

	out_tri = {
		0,1,2,
		3,0,4,
		5,0,6,
		3,6,0,
		0,2,4,
		5,1,0,
		2,1,5,
		7,6,3,
		6,7,5,
		7,3,4,
		7,4,2,
		7,2,5
	};
}

void geneCUBE(std::vector<double>& out_vertices, std::vector<unsigned int>& out_triangles) {
	out_vertices = {
		-2.0f,-1.0f,-1.0f,  //0
		-2.0f,-1.0f, 1.0f,  //1
		-2.0f, 1.0f, 1.0f,  //2
		 2.0f, 1.0f,-1.0f,  //3 
		-2.0f, 1.0f,-1.0f,  //4
		 2.0f,-1.0f, 1.0f,  //5
		 2.0f,-1.0f,-1.0f,  //6
		 2.0f, 1.0f, 1.0f   //7
	};

	out_triangles = {
		0,1,2,
		3,0,4,
		5,0,6,
		3,6,0,
		0,2,4,
		5,1,0,
		2,1,5,
		7,6,3,
		6,7,5,
		7,3,4,
		7,4,2,
		7,2,5
	};

}

bool loadOBJ(
	const char* path,

	std::vector<double>& out_vertices,
	std::vector<unsigned int>& out_triangles
) {
	printf("Loading OBJ file %s...\n", path);


	//std::vector<double> temp_vertices;
	std::vector<unsigned int> vertexIndices;

	FILE* file = fopen(path, "r");
	if (file == NULL) {
		printf("cannot open the file !\n");
		getchar();
		return false;
	}

	while (1) {

		char lineHeader[128];
		// read the first word 
		int res = fscanf(file, "%s", lineHeader);
		if (res == EOF)
			break; // EOF = End Of File Quit 
		// else : parse lineHeader

		if (strcmp(lineHeader, "v") == 0) {
			double x, y, z;
			fscanf(file, "%lf %lf %lf\n", &x, &y, &z);
			std::cout << x << " " << y << " " << z << " " << std::endl;
			out_vertices.push_back(x);
			out_vertices.push_back(y);
			out_vertices.push_back(z);
		}
		else if (strcmp(lineHeader, "f") == 0) {
			unsigned int vertexIndex[3];
			int matches = fscanf(file, "%d %d %d/\n", &vertexIndex[0], &vertexIndex[1], &vertexIndex[2]);
			if (matches != 3) {
				printf("File can't be read");
				fclose(file);
				return false;
			}
			//std::cout << vertexIndex[0] << " " << vertexIndex[1] << " " << vertexIndex[2] << " " << std::endl;
			vertexIndices.push_back(vertexIndex[0] - 1 );
			vertexIndices.push_back(vertexIndex[1] - 1);
			vertexIndices.push_back(vertexIndex[2] - 1);
		}
		else {
			// comment
			char stupidBuffer[1000];
			fgets(stupidBuffer, 1000, file);
		}

	}
	//std::cout << temp_vertices[0] << " " << temp_vertices[1] << " " << temp_vertices[2] << std::endl;
	std::cout << "\n print vertices  index is \n";
	std::cout << vertexIndices.size() << std::endl;

	out_triangles = vertexIndices;
	fclose(file);
	return true;
}


std::vector<double> matrix_mutiple(const std::vector<double>& affine_m, const std::vector<double>& coord) {
	std::vector<double> result;
	for (int i = 0; i < 4; i += 1) {
		double tempt_sum = 0;
		for (int j = 0; j < 4; j += 1) {
			tempt_sum = tempt_sum + coord[j] * affine_m[i * 4 +j];
		}
		result.push_back(tempt_sum);
	}
	return result;
}


std::vector<double> a_rotate_degree(const std::vector<double>& coordinate, int degree, double length) {
	std::vector<double> transformed_coordinate;
	//affine matrix Z rotation 
	std::vector<double> affine_tran_a = { 1,0,0, length/2,
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

std::vector<double> b_rotate_degree(const std::vector<double>& coordinate, int degree_a,int degree_b, double length_a, double length_b) {
	std::vector<double> transformed_coordinate;
	//affine matrix Z rotation 
	std::vector<double> affine_tran_a = { 1,0,0, length_a,
							   0,1,0,0,
							   0,0,1,0,
							   0,0,0,1 };

	std::vector<double> affine_tran_b = { 1,0,0, length_b/2,
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

static void error_callback(int error, const char* description)
{
  fputs(description, stderr);
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    glfwSetWindowShouldClose(window, GL_TRUE);
}

int main(void)
{

  GLFWwindow* window;
  glfwSetErrorCallback(error_callback);
  if (!glfwInit())
    exit(EXIT_FAILURE);
  window = glfwCreateWindow(1080, 880, "Assignment_5_read_BVH", NULL, NULL);
  if (!window)
  {
    glfwTerminate();
    exit(EXIT_FAILURE);
  }

  glfwMakeContextCurrent(window);
  glfwSetKeyCallback(window, key_callback);

  glEnable(GL_DEPTH_TEST);
  // Accept fragment if it closer to the camera than the former one
  glDepthFunc(GL_LESS);


  std::vector<RIGbone> aBone;
  std::vector<CChannel_bvh> aChannelRotTransBone;
  int nframe = 0;
  std::vector<double> aValRotTransBone;

  std::string path_bvh = "jump.bvh";

  loadBVH(aBone, aChannelRotTransBone, nframe, aValRotTransBone,path_bvh);

  std::cout << "nBone:" << aBone.size() << "   aCh:" << aChannelRotTransBone.size() << std::endl;
  for (unsigned int ib = 0; ib < aBone.size(); ++ib) {
	  std::cout << ib << " " << aBone[ib].bname << std::endl;
  }


  double scale_ratio = 390.f;
  while (!glfwWindowShouldClose(window))
  {
	  {
		  static int iframe = 0;
		  const int nch = aChannelRotTransBone.size();
		  SetPose_BioVisionHierarchy(aBone, aChannelRotTransBone,  aValRotTransBone.data() + iframe * nch);

		  iframe = (iframe + 1) % nframe;  // repeat playing this character animation 
	  }
    float ratio;
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    ratio = width / (float) height;

    glViewport(0, 0, width, height);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity(); 
	// last two parameters change the clipping planes 
    glOrtho(-ratio* scale_ratio, ratio * scale_ratio, -1.f * scale_ratio, 1.f * scale_ratio, 10.f * scale_ratio, -10.f * scale_ratio);
	//gluPerspective(120,ratio,0.001f,1.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    //glRotatef((float) glfwGetTime() * 50.f, 0.f, 1.f, 0.f);

	DrawBone(aBone,-1, -1,0.1, 1.0);

	//glBegin(GL_TRIANGLES);
	//glColor3f(0.7f, 1.f, 0.f);
	//glVertex3f(-0.6f, -0.4f, 0.f);
	//glColor3f(0.7f, 1.f, 0.f);
	//glVertex3f(0.6f, -0.4f, 0.f);
	//glColor3f(0.1f, 1.f, 1.f);
	//glVertex3f(0.0f, 0.8f, 0.f);
	//glEnd();

    glfwSwapBuffers(window);
    glfwPollEvents();
  }
  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}