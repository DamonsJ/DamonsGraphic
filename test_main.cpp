#include <iostream>
#include "DamonsVector.h"
#include "DamonsMatrix.h"
#include "DamonsQuaternion.h"
#include "DamonsDistance.h"
#include "DamonsIntersect.h"
#include "DamonsPolygon.h"

using namespace DMath;
using namespace DGraphic;

int main() {

	////DVector<float, 3> v1(1.0f,2.0f,3.0f);
	//DVector<float, 3> v2(11.0f, 12.0f, 13.0f);
	//DVector<float, 3> v = DVector<float,3>::RandomInRange(v1,v2);
	//DVector<float, 3> v3 = DMath::RoundUpToPowerOf2(v1);
	//std::string str = v3.ToString();
	//DVector<float, 3> *v4 = new DVector<float, 3>(5.2f,1.3f,3.2f);
	//std::string str4 = v4->ToString();
	//delete v4;

	DMatrix<float, 4> m1(1.0f, 12.0f, 3.0f, 4.0f,
		4.0f, 5.0f, 16.0f, 7.0f,
		8.0f, 9.0f, 10.0f, 11.0f,
		11.0f, 12.0f, 13.0f, 14.0f);

	DMatrix<float, 4> m2 = m1.Transpose();
	DMatrix<float, 4> m3 = DMatrix<float, 4>::Identity();
	DMatrix<float, 4> m4 = m1 - m3;

	bool isC = m1.InverseWithDeterminantCheck(&m3);

	std::string strm1 = m1.ToString();
	std::string strm2 = m2.ToString();
	std::string strm3 = m3.ToString();

	std::cout << strm1 << std::endl << strm2 << std::endl << strm3;

	DQuaternion<float> q1(0.50f, 0.76f, 0.38f, 0.19f);
	DQuaternion<float> q2(1.40f, 6.3f, 8.5f, 5.9f);
	q1.Normalize();
	q2.Normalize();

	DQuaternion<float> q3 = DQuaternion<float>::Slerp(q1, q2, 1.0f);

	std::string strq1 = q1.ToString();
	std::string strq2 = q2.ToString();
	std::string strq3 = q3.ToString();

	std::cout << strq1 << std::endl << strq2 << std::endl << strq3;

	DPolygon<double > poly;
	poly.AddPoint(2.0, 5.0);
	poly.AddPoint(8.0, 7.0);
	poly.AddPoint(4.0, 10.0);
	poly.AddPoint(3.0, 7.4999);
	bool issimple = poly.IsSimple();
	bool iscounter = poly.IsCounterClockWise();
	bool isconvex = poly.IsConvex();
	int isin = poly.IsPointInPolygon(3,7.499);
 	//////////////////////////////////////////////////////////////////////////
	/*
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

	DLine<float > dl6(pt18, pt19);
	float dt2 = DDistance::PointToLine(pt16,dl6);

	DPoint pt20(00.0f, 20.0f, 20.0f);
	DPoint pt21(0.0f, 0.0f, 20.0f);
	DLine<float > dl7(pt18, pt16);
	DLine<float > dl8(pt21, pt20);
	float dt3 = DDistance::LineToLine(dl7, dl8);

	DPoint pt22(0.0f, 0.0f, 0.0f);
	DPoint pt23(20.0f, 20.0f, 20.0f);
	DPoint pt24(20.0f, 0.0f, 0.0f);
	DPoint pt25(0.0f, 20.0f, 0.0f);
	DLine<float > dl9(pt22, pt23);
	DLine<float > d20(pt24, pt25);
	float dt4 = DDistance::LineToLine(dl9, d20);

	DSegment<float > s1(pt22, pt23);
	DSegment<float > s2(pt24, pt25);
	float dt5 = DDistance::SegmentToSegment(s1, s2);

	DRay<float > r1(pt22, pt23);
	DRay<float > r2(pt24, pt25);
	float dt6 = DDistance::RayToRay(r1, r2);

	float dt7 = DDistance::LineToRay(dl9, r2);
	float dt8 = DDistance::LineToSegment(dl9, s2);
	float dt9 = DDistance::RayToSegment(r1, s2);

	DPoint res;
	DDirection<float > dri3(1, 0, 0);
	DPoint pt26(-20.0f, 20.0f, 0.0f);
	DPlane<float > plane3(pt26, dri3);

	bool isinter = DIntersection::LineWithPlane(d20,plane3,res);
	bool isinte1 = DIntersection::RayWithPlane(r2, plane3, res);
	bool isinte2 = DIntersection::SegmentWithPlane(s2, plane3, res);

	DPoint pt27(-20.0f, 0.0f, 0.0f);
	DPoint pt28(-20.0f, 60.0f, 0.0f);
	DPoint pt29(-20.0f, 60.0f, 20.0f);
	DTriangle<float> tri1(pt27,pt28,pt29);

	bool isinter3 = DIntersection::LineWithTriangle(d20, tri1, res);
	bool isinter4 = DIntersection::RayWithTriangle(r2, tri1, res);
	bool isinter5 = DIntersection::SegmentWithTriangle(s2, tri1, res);

	DDirection<float > dri4(0, 1, 0);
	DPlane<float > plane4(pt22, dri4);

	DObject obj;
	DIntersection::TrianglePlaneIntersectType tp = DIntersection::PlaneWithTriangle(plane4, tri1, obj);

	DPlane<float > plane5(pt26, dri4);
	DIntersection::TrianglePlaneIntersectType tp1 = DIntersection::PlaneWithTriangle(plane5, tri1, obj);

	DBox<float> box_r(pt22, pt23);
	DPoint ray_pt1(0.0f,20.0001f,0.0f);
	DDirection<float > ray_dir1(1.0f, 1.0f, 1.0f);
	DPoint  ray_pt2(0.0f, -0.001f, 0.0f);
	DDirection<float > ray_dir2(1.0f, 0.0f, 0.0f);
	DRay<float> ray_b1(ray_pt1, ray_dir1);
	DRay<float> ray_b2(ray_pt2, ray_dir2);

	bool overlap1 = DIntersection::RayWithAABB(ray_b1, box_r);
	bool overlap2 = DIntersection::RayWithAABB(ray_b2, box_r);


	DPoint pt30(0.0f, 100.0f, 0.0f);
	DPoint pt31(0.0f, 0.0f, 0.0f);
	DPoint pt32(0.0f, 100.0f, 80.0f);
	DTriangle<float> trit(pt30, pt31, pt32);

	DPoint pt33(-1.5f, 96.5f, 77.5f);
	DPoint pt34(15.5f, 96.5f, 77.5f);
	DLine<float > dlt(pt33, pt34);
	bool isintert = DIntersection::LineWithTriangle(dlt, trit, res);
	*/
	system("pause");
	return 0;
}