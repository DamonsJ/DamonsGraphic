#ifndef _DAMONS_INTERSECT_H_
#define _DAMONS_INTERSECT_H_

#include "DamonsPoint.h"
#include "DamonsSegment.h"
#include "DamonsDirection.h"
#include "DamonsBox.h"
#include "DamonsTriangle.h"
#include "DamonsPlane.h"
#include "DamonsLine.h"

namespace DGraphic {

	/// class  DIntersection
	/// @breif DIntersection compute 3 dimension entity intersection between another
	///        it is an operation class that contains no data

	class DIntersection {

		public:
			///@breif define a uninitialize constructor
			DIntersection() {

			}
			///@breif define a uninitialize deallocator
			~DIntersection() {

			}

		public:
			/// @brief Calculate the intersection between line and plane
			///
			/// @param p DPlane of type T
			/// @param l DLine  of type T.
			/// @param result is intersect point as return
			/// @return a bool value whether  plane  and line is intersected
			template<class T>
			static bool LineWithPlane(DLine<T> &l,DPlane<T> &p, DPoint<T > &result) {
				DDirection<T> ld = l.Direction();
				DDirection<T> pn = p.GetDirection();

				T det = ld.DotProduct(pn);
				/// line is parallel with plane
				if (std::abs(det) < DEplision) {
				///check whether the line is lie on plane
					if (p.Has_on(l.GetSourcePoint())) {
						result = l.GetSourcePoint();
						return true;
					}
					else {
						return false;
					}
				}

				T t = -(pn.DotProduct(l.GetSourcePoint())- pn.DotProduct(p.GetOrigin())) / det;

				result = l.GetPointOnLine(t);

				return true;
			}

			/// @brief Calculate the intersection between ray and plane
			///
			/// @param p DPlane of type T
			/// @param l DRay  of type T.
			/// @param result is intersect point as return
			/// @return a bool value whether  plane  and ray is intersected
			template<class T>
			static bool RayWithPlane(DRay<T> &l, DPlane<T> &p, DPoint<T > &result) {
				DDirection<T> ld = l.Direction();
				DDirection<T> pn = p.GetDirection();

				T det = ld.DotProduct(pn);
				/// line is parallel with plane
				if (std::abs(det) < DEplision) {
					///check whether the line is lie on plane
					if (p.Has_on(l.GetSourcePoint())) {
						result = l.GetSourcePoint();
						return true;
					}
					else {
						return false;
					}
				}

				T t = -(pn.DotProduct(l.GetSourcePoint()) - pn.DotProduct(p.GetOrigin())) / det;
				if (t > 0) {
					result = l.GetSourcePoint()+DPoint<T>(t*ld.x(), t*ld.y(), t*ld.z());
					return true;
				}
				else
					return false;

			}

			/// @brief Calculate the intersection between segment and plane
			///
			/// @param p DPlane of type T
			/// @param l DSegment  of type T.
			/// @param result is intersect point as return
			/// @return a bool value whether  plane  and segment is intersected
			template<class T>
			static bool SegmentWithPlane(DSegment<T> &l, DPlane<T> &p, DPoint<T > &result) {
				DDirection<T> ld = l.Direction();
				DDirection<T> pn = p.GetDirection();

				T det = ld.DotProduct(pn);
				/// line is parallel with plane
				if (std::abs(det) < DEplision) {
					///check whether the line is lie on plane
					if (p.Has_on(l.GetStartPoint())) {
						result = l.GetStartPoint();
						return true;
					}
					else {
						return false;
					}
				}

				T t = -(pn.DotProduct(l.GetStartPoint()) - pn.DotProduct(p.GetOrigin())) / det;
				if (t > 0 && t < l.Length()) {
					result = l.GetStartPoint() + DPoint<T>(t*ld.x(), t*ld.y(), t*ld.z());
					return true;
				}
				else
					return false;
			}

			/// @brief Calculate the intersection between line and triangle
			///
			/// @param tri DTriangle of type T
			/// @param l DLine  of type T.
			/// @param result is intersect point as return
			/// @return a bool value whether triangle and line is intersected
			template<class T>
			static bool LineWithTriangle(DLine<T> &l, DTriangle<T> &tri, DPoint<T > &result) {
				DDirection<T> tn = tri.Normal();
				DDirection<T> ld = l.Direction();

				DPoint<T> e1 = tri[1] - tri[0];
				DPoint<T> e2 = tri[2] - tri[0];

				DPoint<T> p = ld.CrossProduct(e2);
				T det = e1.DotProduct(p);

				if (det < DEplision && det > -DEplision)
					return false;

				T OneOverDet = 1.0 / det;
	
				DPoint<T> s = l.GetSourcePoint() - tri[0];
				T u = OneOverDet*(s.DotProduct(p));
				if (u < 0.0 || u > 1.0)
					return false;

				DPoint<T> q = s.CrossProduct(e1);

				T v = OneOverDet*(ld.DotProduct(q));
				if (v < 0.0 || v > 1.0)
					return false;

				T t = (e2.DotProduct(q))* OneOverDet;
				
				result = l.GetPointOnLine(t);

				return true;
			}

			/// @brief Calculate the intersection between ray and triangle
			///
			/// @param tri DTriangle of type T
			/// @param l DRay  of type T.
			/// @param result is intersect point as return
			/// @return a bool value whether triangle and ray is intersected
			template<class T>
			static bool RayWithTriangle(DRay<T> &l, DTriangle<T> &tri, DPoint<T > &result) {
				DDirection<T> tn = tri.Normal();
				DDirection<T> ld = l.Direction();

				DPoint<T> e1 = tri[1] - tri[0];
				DPoint<T> e2 = tri[2] - tri[0];

				DPoint<T> p = ld.CrossProduct(e2);
				T det = e1.DotProduct(p);

				if (det < DEplision && det > -DEplision)
					return false;

				T OneOverDet = 1.0 / det;

				DPoint<T> s = l.GetSourcePoint() - tri[0];
				T u = OneOverDet*(s.DotProduct(p));
				if (u < 0.0 || u > 1.0)
					return false;

				DPoint<T> q = s.CrossProduct(e1);

				T v = OneOverDet*(ld.DotProduct(q));
				if (v < 0.0 || v > 1.0)
					return false;

				T t = (e2.DotProduct(q))* OneOverDet;

				if (t > 0) {
					result = l.GetSourcePoint() + DPoint<T>(t*ld.x(), t*ld.y(), t*ld.z());
					return true;
				}
				else
					return false;

			}

			/// @brief Calculate the intersection between segment and triangle
			///
			/// @param tri DTriangle of type T
			/// @param l DSegment  of type T.
			/// @param result is intersect point as return
			/// @return a bool value whether triangle and segment is intersected
			template<class T>
			static bool SegmentWithTriangle(DSegment<T> &l, DTriangle<T> &tri, DPoint<T > &result) {
				DDirection<T> tn = tri.Normal();
				DDirection<T> ld = l.Direction();

				DPoint<T> e1 = tri[1] - tri[0];
				DPoint<T> e2 = tri[2] - tri[0];

				DPoint<T> p = ld.CrossProduct(e2);
				T det = e1.DotProduct(p);

				if (det < DEplision && det > -DEplision)
					return false;

				T OneOverDet = 1.0 / det;

				DPoint<T> s = l.GetStartPoint() - tri[0];
				T u = OneOverDet*(s.DotProduct(p));
				if (u < 0.0 || u > 1.0)
					return false;

				DPoint<T> q = s.CrossProduct(e1);

				T v = OneOverDet*(ld.DotProduct(q));
				if (v < 0.0 || v > 1.0)
					return false;

				T t = (e2.DotProduct(q))* OneOverDet;

				if (t > 0 && t < l.Length()) {
					result = l.GetStartPoint() + DPoint<T>(t*ld.x(), t*ld.y(), t*ld.z());
					return true;
				}
				else
					return false;

			}
	};
};

#endif