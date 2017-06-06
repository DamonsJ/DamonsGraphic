#include <iostream>
#include "DamonsVector.h"
#include "DamonsMatrix.h"
#include "DamonsQuaternion.h"
#include "DamonsDistance.h"

using namespace DMath;
using namespace DGraphic;

int main() {

	// 	DVector<float, 3> v1(1.0f,2.0f,3.0f);
	// 	DVector<float, 3> v2(11.0f, 12.0f, 13.0f);
	// 	DVector<float, 3> v = DVector<float,3>::RandomInRange(v1,v2);
	// 	DVector<float, 3> v3 = DMath::RoundUpToPowerOf2(v1);
	// 	std::string str = v3.ToString();
	// 	DVector<float, 3> *v4 = new DVector<float, 3>(5.2f,1.3f,3.2f);
	// 	std::string str4 = v4->ToString();
	// 	delete v4;
	// 
	// 	DMatrix<float, 4> m1(1.0f,12.0f,3.0f,  4.0f,
	// 						 4.0f, 5.0f, 16.0f, 7.0f,
	// 						 8.0f, 9.0f, 10.0f, 11.0f,
	// 						11.0f, 12.0f, 13.0f,14.0f);
	// 
	// 	DMatrix<float, 4> m2 = m1.Transpose();
	// 	DMatrix<float, 4> m3 = DMatrix<float, 4>::Identity();
	// 
	// 	bool isC = m1.InverseWithDeterminantCheck(&m3);
	// 
	// 	std::string strm1 = m1.ToString();
	// 	std::string strm2 = m2.ToString();
	// 	std::string strm3 = m3.ToString();
	// 
	// 	std::cout << strm1 << std::endl<<strm2 << std::endl << strm3;
	// 
	// 	DQuaternion<float> q1(0.50f, 0.76f, 0.38f,0.19f);
	// 	DQuaternion<float> q2(1.40f, 6.3f,8.5f, 5.9f);
	// 	q1.Normalize();
	// 	q2.Normalize();
	// 
	// 	DQuaternion<float> q3 = DQuaternion<float>::Slerp(q1, q2, 1.0f);
	// 
	// 	std::string strq1 = q1.ToString();
	// 	std::string strq2 = q2.ToString();
	// 	std::string strq3 = q3.ToString();
	// 
	// 	std::cout << strq1 << std::endl << strq2 << std::endl << strq3;
	// 	//////////////////////////////////////////////////////////////////////////

	typedef DPoint<float > DPoint;
	DPoint pt1(1.2, 3.1f, 4);
	DPoint pt2(pt1);
	float fv = pt1(0);
	float sv = pt2[2];
	std::string ptstr1 = pt1.ToString();

	DPoint *pt3 = new DPoint(1.2f, 3.1f, 4.5f);
	delete pt3;

	DPoint pt4(pt1 + pt2);
	pt4.Normalize();
	DPoint pt5 = pt4.Normalized();
	DPoint pt6 = -pt4;
	DPoint pt7 = 1.0f + pt4;;
	DPoint pt8 = pt4 + 1.0f;
	float d = pt7.DotProduct(pt8);
	DPoint pt11 = pt7.CrossProduct(pt8);
	DPoint pt10 = pt7*pt8;
	DPoint pt9 = pt7 / pt8;
	pt9 += 1.0f;

	DSegment<float > dl1(pt1, pt9);
	DSegment<float > dl2(pt1, pt9);

	std::wstring wst = dl1.ToWString();
	std::string str = dl1.ToString();

	std::cout << str;
	//DPoint dld = dl1.Direction();
	DDirection<float> dld = dl1.Direction();

	bool iss = (dl1 == dl2);

	DDirection<float> dir(dl2);
	dir = -dir;

	DPoint pt13(0.0f, 0.0f, 10.0f);
	DPoint pt14(20.0f, 0.0f, 10.0f);
	DPoint pt15(10.00001f, 10.0000001f, 0.0000001f);
	DSegment<float > dl3(pt13, pt14);
	bool isOn = dl3.Has_on(pt15);

	DBox<float > box1(pt13, pt14);
	DBox<float > box2(box1);
	bool isit = box1.IsIntersect(box2);
	bool isiw = box1 == (box2);

	DPoint pt16(30.00001f, 0.0000001f, 0.0000000f);
	DTriangle<float > tri(pt13, pt14, pt15);
	bool iso = tri.Has_on(pt16);
	float area = tri.Area();
	DBox<float > box = tri.Box();
	DDirection<float > dri = tri.Normal();

	DDirection<float > dri2(1, 0, 0);
	DPlane<float > plane1(pt13, dri2);
	bool iss1 = plane1 == plane1;
	DPoint pt17 = plane1.Projection(pt15);
	DPlane<float > plane2 = plane1.Opposite();
	bool ison = plane1.Has_on(pt16);

	DPoint pt18(0.0f, 0.0f, 0.0f);
	DPoint pt19(20.0f, 20.0f, 20.0f);
	DSegment<float > dl4(pt18, pt19);
	float ds = DDistance::PointToSegment(pt16, dl4);

	DRay<float > dl5(pt18, pt19);
	float dt = DDistance::PointToRay(pt16, dl5);
	float dt1 = DDistance::PointToTriangle(pt16,tri);

	system("pause");
	return 0;
}