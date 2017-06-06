#ifndef _DAMONS_DISTANCE_H_
#define _DAMONS_DISTANCE_H_

#include "DamonsPoint.h"
#include "DamonsSegment.h"
#include "DamonsDirection.h"
#include "DamonsBox.h"
#include "DamonsTriangle.h"
#include "DamonsPlane.h"

namespace DGraphic {

	/// class  DDistance
	/// @breif DDistance compute 3 dimension entity distance between another
	///        it is an operation class that contains no data

	class DDistance
	{
	public:
		///@breif define a uninitialize constructor
		DDistance() {

		}
		///@breif define a uninitialize deallocator
		~DDistance() {

		}

	public:
		/// @brief Calculate the distance between point and line
		///
		/// @param p DPoint of type T
		/// @param l DSegment  of type T.
		/// @param The distance between point  and line

		template <class T>
		static inline T PointToSegment(const DPoint<T> &p, const DSegment<T> &l) {
			/// assume the line is Ps + (Pe-Ps)*t,where t >= 0 and t <= 1
			/// if the nearest point of p lie between the start and end point of line
			/// then we can calculate t,and the nearest point ,so we can calculate distance
			/// otherwise the distance is distance between p and the two end points of line
			DPoint<T> vec = p - l.GetStartPoint();
			DPoint<T> v = DPoint<T>(l.to_vector());

			T t = vec.DotProduct(v) / l.SquaredLength();
			if (t < static_cast<T>(DEplision)) {
				return vec.Length();
			}
			else if ((t + DEplision) > static_cast<T>(1)) {
				return (p - l.GetEndPoint()).Length();
			}
			else {
				DPoint<T> r = (p - (l.GetStartPoint() + t*v));
				return r.Length();
			}

		}

		/// @brief Calculate the distance between point and ray
		///
		/// @param p DPoint of type T
		/// @param l DRay  of type T.
		/// @return The distance between point and ray

		template <class T>
		static inline T PointToRay(const DPoint<T> &p, const DRay<T> &r) {
			///same method with PointToLine but different t value
			DPoint<T> vec = p - r.GetSourcePoint();
			const DDirection<T> dir = r.Direction();
			DPoint<T> v = DPoint<T>(dir.x(), dir.y(), dir.z());

			T t = vec.DotProduct(v);
			if (t < static_cast<T>(DEplision)) {
				return vec.Length();
			}
			else {
				DPoint<T> res = (p - (r.GetSourcePoint() + t*v));
				return res.Length();
			}
		}

		/// @brief Calculate the distance between point and plane
		///
		/// @param p DPoint of type T
		/// @param pl DPlane  of type T.
		/// @return The distance between point and plane

		template <class T>
		static inline T PointToPlane(const DPoint<T> &p, const DPlane<T> &pl) {

			DDirection<T> dir = pl.GetDirection();
			DPoint<T> norm = DPoint<T>(dir.x(), dir.y(), dir.z());
			DPoint<T> orgpt = pl.GetOrigin();

			T s1 = norm.DotProduct(p);
			T s2 = norm.DotProduct(orgpt);

			return (s1 - s2);
		}

		/// @brief Calculate the distance between point and triangle
		///
		/// @param p DPoint of type T
		/// @param tri DTriangle  of type T.
		/// @return The distance between point and triangle

		template <class T>
		static inline T PointToTriangle(const DPoint<T> &p, const DTriangle<T> &tri) {

			///1. calculate the nearest point lie on triangle
			DDirection<T> dir = tri.Normal();
			DPoint<T> norm = DPoint<T>(dir.x(), dir.y(), dir.z());
			DPoint<T> orgpt = tri[0];

			T s1 = norm.DotProduct(p);
			T s2 = norm.DotProduct(orgpt);
			DPoint<T> npt = (p - (s1-s2)*norm);
			
			///2. calculate whether the nearest point is lie in triangle
			/// if lie in triangle,the distance is p and the nearest point
			/// otherwise calculate the nearest distance between p and three lines of triangle
			if (tri.Has_on(npt)) {
				return p.Distance(npt);
			}
			else {
				DSegment<T> l1 = DSegment<T>(tri[0], tri[1]);
				DSegment<T> l2 = DSegment<T>(tri[1], tri[2]);
				DSegment<T> l3 = DSegment<T>(tri[2], tri[0]);

				T d1 = PointToSegment(p, l1);
				T d2 = PointToSegment(p, l2);
				T d3 = PointToSegment(p, l3);

				return std::min(d1, std::min(d2,d3));
			}
		}

	protected:
	};
};
#endif
