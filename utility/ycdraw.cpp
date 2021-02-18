#include "ycdraw.h"
#include <GLFW/glfw3.h>
#include <deque>

void YCdraw::DrawSphere
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

void YCdraw::DrawSphereAt
(int nla, int nlo, double rad, double x, double y, double z)
{
	glTranslated(+x, +y, +z);
	glScaled(rad, rad, rad);
	DrawSphere(nla, nlo);
	glScaled(1.0 / rad, 1.0 / rad, 1.0 / rad);
	glTranslated(-x, -y, -z);
}

void YCdraw::DrawTrajectory
(std::deque<double> trajectory)
{
	::glBegin(GL_LINE_STRIP);
	::glLineWidth(1.5f);
	::glColor3f(0.0f, 0.5f, 0.9f);
	for (int i = 0; i < trajectory.size()/3; i++)
	{

		::glVertex3f(trajectory[i * 3], trajectory[i * 3 + 1], trajectory[i * 3 + 2]);
		//DrawSphereAt(32, 32, 0.2, trajectory[i * 3], trajectory[i * 3 + 1], trajectory[i * 3 + 2]);
	}
	::glEnd();
}

void YCdraw::DrawTrajectory_points
(std::deque<double> trajectory)
{
	::glColor3f(1.0f, 1.0f, 1.0f);
	for (int i = 0; i < trajectory.size() / 3; i = i + 9)
	{
		DrawSphereAt(32, 32, 0.2, trajectory[i * 3], trajectory[i * 3 + 1], trajectory[i * 3 + 2]);
	}
	::glEnd();
}