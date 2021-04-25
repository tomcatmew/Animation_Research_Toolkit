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
#include <deque>
// my utilities include folder files
#include "ycmatrix.cpp"
#include "ycdraw.cpp"

using namespace YCutility;
using namespace YCdraw;

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

//[1]BELOW Utilities functions updatebone and setpose ----------------[1]
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

		std::vector<double> temp_result;


		std::vector<double> quat_matrix = Quat(aBone[ibone].quatRelativeRot);

		temp_result = matrix_m_matrix(tempt_tranform4d, quat_matrix);

		// quertion way, and use GLM matrix 4D - Prof.Umetani's way
		//CMat4d m01 = CMat4d::Translate(aBone[ibone].transRelative);
		//m01 = m01 * CMat4d::Quat(aBone[ibone].quatRelativeRot);
		//m01 = m01 * CMat4d::Scale(aBone[ibone].scale);

		const int ibone_p = aBone[ibone].ibone_parent;
		if (ibone_p < 0 || ibone_p >= (int)aBone.size()) { // root bone
			Copy_Mat4(aBone[ibone].affmat3Global, temp_result);
			continue;
		}
		MatMat4(aBone[ibone].affmat3Global,
			aBone[ibone_p].affmat3Global, temp_result);
	   
	}
}

void UpdateBoneRotTrans_XYZ
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
												 0, cos(X_degree),-sin(X_degree), 0 ,
												 0, sin(X_degree), cos(X_degree), 0,
												 0, 0, 0, 1 };

		std::vector<double> tempt_Yrotate_4d = { cos(Y_degree), 0,sin(Y_degree), 0,
												 0, 1, 0, 0 ,
												 -sin(Y_degree), 0,cos(Y_degree), 0,
												 0, 0, 0, 1 };

		std::vector<double> tempt_Zrotate_4d = { cos(Z_degree),-sin(Z_degree),0, 0,
												 sin(Z_degree), cos(Z_degree),0, 0,
												 0, 0, 1, 0 ,
												 0, 0, 0, 1 };

		std::vector<double> temp_result;
		std::vector<double> temp_result2;
		std::vector<double> temp_result3;

		std::vector<double> quat_matrix = Quat(aBone[ibone].quatRelativeRot);

		temp_result = matrix_m_matrix(tempt_tranform4d, tempt_Xrotate_4d);
		temp_result2 = matrix_m_matrix(temp_result, tempt_Yrotate_4d);
		temp_result3 = matrix_m_matrix(temp_result2, tempt_Zrotate_4d);

		const int ibone_p = aBone[ibone].ibone_parent;
		if (ibone_p < 0 || ibone_p >= (int)aBone.size()) { // root bone
			Copy_Mat4(aBone[ibone].affmat3Global, temp_result3);
			continue;
		}
		MatMat4(aBone[ibone].affmat3Global,
			aBone[ibone_p].affmat3Global, temp_result3);
	}
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
		}
		else {
			// this is for rotation XYZ ========== (my way store degree in quatRelativeRot instead of quenteron )
			//double ar = val * 3.1415926 / 180.0;
			//aBone[ibone].quatRelativeRot[1 + iaxis] = ar;
			// this is for rotation XYZ ==========

			// the way of quaternion method   =======
			
			//calculate z y X rotation's cooresponding [w,x,y,z] quaternion one by one and use quat*quat to mutiple them together. 
			//finally we will get the [w,x,y,z] quaternion representation and store it into the memory. 
			//when we draw the bones, we use [x,x,y,z] and transform to 4x4 rotation matrix then simply mutiple to our object.
			const double ar = val * 3.1415926 / 180.0;
			double v0[3] = { 0,0,0 };
			v0[iaxis] = 1.0;
			double dq[4] = { cos(ar * 0.5), v0[0] * sin(ar * 0.5), v0[1] * sin(ar * 0.5), v0[2] * sin(ar * 0.5) };
			double qtmp[4]; QuatQuat(qtmp,
				aBone[ibone].quatRelativeRot, dq);
			Copy_Quat(aBone[ibone].quatRelativeRot, qtmp);
			
			// the way of quaternion method   =======
		}
	}
	//UpdateBoneRotTrans_XYZ(aBone);
	UpdateBoneRotTrans(aBone);
}

void SetPose_BioVisionHierarchy_Quat
(std::vector<RIGbone>& aBone,
	const std::vector<CChannel_bvh>& aChannelRotTransBone,
	const double* aVal)
{
	const int nch = aChannelRotTransBone.size();
	for (int ich = 0; ich < nch; ++ich) {
		const int ibone = aChannelRotTransBone[ich].ibone;
		const int iaxis = aChannelRotTransBone[ich].iaxis;
		const bool isrot = aChannelRotTransBone[ich].isrot;
		const double val = aVal[ich];
		std::cout << val << " ";
		if (ich == nch - 1)
		{
			std::cout << std::endl;
		}

		if (ich == 0 || ich == 1 || ich == 2) {
			aBone[ibone].transRelative[iaxis] = val;
		}
		else {

			int check_point = ich % 3;
			if (check_point == 0)
			{
				int mul = ich / 3;
				int start = mul * 4 - 1;
				double dq[4] = { aVal[start],aVal[start + 1] ,aVal[start + 2] ,aVal[start + 3] };
				Copy_Quat(aBone[ibone].quatRelativeRot, dq);
			}
			// the way of quaternion method   =======
		}
	}
	//UpdateBoneRotTrans_XYZ(aBone);
	UpdateBoneRotTrans(aBone);
}

//[1]ABOVE ARE UTILITY FUNCTIONS -------------------------------------[1]


//[2]load BVH file main program --------------------------------------[2]

bool loadQ(std::vector<double>& aValueQ, int& nframe, const std::string& path_q)
{
	int nchannel = 127;
	std::ifstream fin;
	fin.open(path_q.c_str());
	if (!fin.is_open()) {
		std::cout << "cannot open file" << std::endl;
		return false;
	}
	//
	aValueQ.resize(nframe * nchannel);
	std::string line;
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

		for (int ich = 0; ich < nchannel; ++ich) {
			//aValueRotTransBone[iframe * nchannel + ich] = rig_v3q::myStod(aToken[ich]);
			if (ich == 1)
			{
				aValueQ[iframe * nchannel + ich] = atof(aToken[ich].c_str()) - 180;
			}
			else {
				aValueQ[iframe * nchannel + ich] = atof(aToken[ich].c_str());
			}
		}
	}
}

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
			if (ich == 1)
			{
				aValueRotTransBone[iframe * nchannel + ich] = atof(aToken[ich].c_str()) - 180;
			}
			else {
				aValueRotTransBone[iframe * nchannel + ich] = atof(aToken[ich].c_str());
			}
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
			bone.transRelative[0] = (-bone.invBindMat[3]) - (-bone_p.invBindMat[3]) ;
			bone.transRelative[1] = (-bone.invBindMat[7]) - (-bone_p.invBindMat[7]) ;
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
//[2]ABOVE is BVH parser ---------------------------------------------[2]



//[3]professor's draw bones function ---------------------------------[3]
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
		glLineWidth(4);
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
//[3]ABOVE are draw functions ---------------------------------------[3]


//[4]BELOW are generating cubes and loadOBJ functions ---------------[4]
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
	//std::cout << "\n print vertices  index is \n";
	//std::cout << vertexIndices.size() << std::endl;

	out_triangles = vertexIndices;
	fclose(file);
	return true;
}
//[4]ABOVE are generating cubes and loadOBJ functions ---------------[4]

double cam_y = 170.f;
double camd_x = 0.f;

void move_cam_up() {
	cam_y = cam_y - 2.f;
}

void move_cam_down() {
	cam_y = cam_y + 2.f;
}

void rotate_cam_up() {
	camd_x = camd_x + 2.f;
}

void rotate_cam_down() {
	camd_x = camd_x - 2.f;
}

//[5]BELOW is main() call -------------------------------------------[5]
static void error_callback(int error, const char* description)
{
  fputs(description, stderr);
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    glfwSetWindowShouldClose(window, GL_TRUE);
  if (key == GLFW_KEY_W && action == GLFW_PRESS)
	  move_cam_up();
  if (key == GLFW_KEY_S && action == GLFW_PRESS)
	  move_cam_down();
  if (key == GLFW_KEY_I && action == GLFW_PRESS)
	  rotate_cam_up();
  if (key == GLFW_KEY_K && action == GLFW_PRESS)
	  rotate_cam_down();
}

int main(void)
{
	std::deque<double> pos_history;

  
  // contains 10 frames + data
  // [x1,y1,z1,x2,y2,z2] 
  // using std::deque 
  // each keyframe add new things , remove old things,
  
  //5 
  // current in 10th frame, 5 - 10 frames in history.
  // first frame, size of 0,  2nd frame - add into pos_history, 3rd frame - add into pos_history, ....... keep doing until reach a limit = 100 ,  
  // then remove first, add from back.


  GLFWwindow* window;
  glfwSetErrorCallback(error_callback);
  if (!glfwInit())
    exit(EXIT_FAILURE);
  window = glfwCreateWindow(1080, 880, "5_read_BVH", NULL, NULL);
  if (!window)
  {
    glfwTerminate();
    exit(EXIT_FAILURE);
  }

  glfwMakeContextCurrent(window);
  glfwSetKeyCallback(window, key_callback);

  // Accept fragment if it closer to the camera than the former one
  glDepthFunc(GL_LESS);


  std::vector<RIGbone> aBone;
  std::vector<CChannel_bvh> aChannelRotTransBone;
  int nframe = 0;
  std::vector<double> aValRotTransBone;

  std::string path_bvh = "model/LocomotionFlat01_000.bvh";
  //std::string path_bvh = "10_01.bvh";

  loadBVH(aBone, aChannelRotTransBone, nframe, aValRotTransBone,path_bvh);


  std::vector<double> aValueQ;
  int qframe = 8169;
  std::string path_q = "model/output_direct_radian_q.bvh";

  loadQ(aValueQ, qframe, path_q);

  std::cout << "nBone:" << aBone.size() << "   aCh:" << aChannelRotTransBone.size() << std::endl;
  std::cout << "qSize: " << aValueQ.size() << std::endl;
  int lefts;
  int rights;
  for (unsigned int ib = 0; ib < aBone.size(); ++ib) {
	  if (aBone[ib].bname == "LeftShoulder")
		  lefts = ib;
	  else if (aBone[ib].bname == "RightShoulder")
		  rights = ib;
	  //std::cout << ib << " " << aBone[ib].bname << std::endl;
  }
  //std::cout << rights << " " << lefts << std::endl;


  //initialization
  SetPose_BioVisionHierarchy(aBone, aChannelRotTransBone, aValRotTransBone.data());
  std::vector<double> start_pos;
  start_pos = aBone[0].Pos();
  for (int i = 0; i < 300; i++) {
	  if ((i + 1) % 3 == 2) {
		  pos_history.push_back(-180.0f);
	  }
	  else if ((i + 1) % 3 == 1) {
		  pos_history.push_back(start_pos[0]);
	  }
	  else{
		  pos_history.push_back(start_pos[2]);
	  }
  }


  //uncomment for PFNN data
  double scale_ratio = 30.f;

  //jump.bvh
  //double scale_ratio = 250.f;   

  //CMU bvh
  //double scale_ratio = 40.f; 
  while (!glfwWindowShouldClose(window))
  {
	  //manipulate the deque hip_pos_history

	  std::vector<double> trajectory;
	  std::vector<double> left_spos;
	  std::vector<double> right_spos;
	  double cross_p[3];
	  double upward[3] = {0, 10.0, 0};
	  double horizontal_s[3];
	  {
		  static int iframe = 0;
		  //const int nch = aChannelRotTransBone.size();
		  const int nch = 127;
		  std::cout << "===frame: " <<iframe << std::endl;
		  //SetPose_BioVisionHierarchy(aBone, aChannelRotTransBone,  aValRotTransBone.data() + iframe * nch);
		  SetPose_BioVisionHierarchy_Quat(aBone, aChannelRotTransBone, aValueQ.data() + iframe * nch);
		  trajectory = aBone[0].Pos();

		  //update position history
		  for (int i = 0; i < 3; i++) {
			  pos_history.pop_back();
		  }
		  pos_history.push_front(trajectory[2]);
		  pos_history.push_front(-180.f);
		  pos_history.push_front(trajectory[0]);
		  std::cout << pos_history.size() << std::endl;
		  left_spos = aBone[lefts + 1].Pos();

		  right_spos = aBone[rights + 1].Pos();

		  horizontal_s[0] = right_spos[0] - left_spos[0];
		  horizontal_s[1] = right_spos[1] - left_spos[1];
		  horizontal_s[2] = right_spos[2] - left_spos[2];

		  //crossproduct with upward direction
		  cross_p[0] = upward[1] * horizontal_s[2] - upward[2] * horizontal_s[1];
		  cross_p[1] = upward[2] * horizontal_s[0] - upward[0] * horizontal_s[2];
	      cross_p[2] = upward[0] * horizontal_s[1] - upward[1] * horizontal_s[0];

		  double arrow[3] = {cross_p[0] - trajectory[0],cross_p[1] - trajectory[1], cross_p[2] - trajectory[2] };
		  // we need to unit it based on character's position, not origin.
		  double unit = sqrt(arrow[0] * arrow[0] + arrow[1] * arrow[1] + arrow[2] * arrow[2]);
		  int arrow_scale = 15;
		  arrow[0] = arrow[0] / unit * arrow_scale;
		  arrow[1] = arrow[1] / unit * arrow_scale;
		  arrow[2] = arrow[2] / unit * arrow_scale;
	
		  cross_p[0] = trajectory[0] + arrow[0];
		  cross_p[1] = trajectory[1] + arrow[1];
		  cross_p[2] = trajectory[2] + arrow[2];

		  // repeat playing this character animation 
		  //iframe = (iframe + 1) % nframe;  
		  iframe = (iframe + 1) % qframe;
		  std::cout << " Frame : " << iframe << std::endl;
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
	//mat4x4_perspective(p, 1.57, width / (float)height, 1, 10);
    glOrtho(-ratio* scale_ratio, ratio * scale_ratio, -1.f * scale_ratio, 1.f * scale_ratio, 100.f * scale_ratio, -100.f * scale_ratio);
	//gluPerspective(120,ratio,0.001f,1.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	// uncomment for PFNN data
	glTranslated(0.f, cam_y, 0.f);

    // uncomment for 02_01.bvh data
	//glTranslated(0.f, 170.f, 0.f);

    //glRotatef((float) glfwGetTime() * 50.f, 0.f, 1.f, 0.f);

	glRotatef(camd_x, 1.0f, 0.f, 0.f);
	glEnable(GL_DEPTH_TEST);
	DrawBone(aBone,-1, -1,0.1, 1.0);

	// projection the hip point to ground
	glColor3f(0.7f, 0.3f, 0.2f);
	glEnable(GL_DEPTH_TEST);
	DrawSphereAt(32, 32, 0.2,trajectory[0], -180.f, trajectory[2]);

	glEnable(GL_DEPTH_TEST);
	glBegin(GL_LINE_STRIP);
	glLineWidth(3.0f);
	glColor3f(0.7f, 0.3f, 0.2f);
	glVertex3f(trajectory[0], -180.f, trajectory[2]);
	glVertex3f(cross_p[0] , -180.f, cross_p[2] );
	glEnd();

	DrawTrajectory(pos_history);
	DrawTrajectory_points(pos_history);
    glfwSwapBuffers(window);
    glfwPollEvents();
  }
  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}