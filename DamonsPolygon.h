/*
* learning from google's opensource project :mathfu
*/

#ifndef _DAMONSPOLYGON_HEADER_
#define _DAMONSPOLYGON_HEADER_

#include "DamonsObject.h"
#include "DamonsVector.h"
#include <vector>

using namespace DMath;
namespace DGraphic {
	///class DPolygon
	///@breif polygon with plane coordinates with Type T
	/// 
	/// polygon stores arbitrary number of 2d points of type <b>T</b> and provides a set
	/// functions to perform operations on polygon.
	///
	/// @tparam T type of polygon point.
	/// @note note that polygon is default closed, if input polygon has different
	/// first and last point,then the first point will be pushed .
	template< class T = double>
	class DPolygon : public DObject
	{
	public:
		/// @brief Create an uninitialized polygon.
		inline DPolygon() {}

		/// @brief create a polygon with points vector
		///
		/// @param points points of polygon

		inline DPolygon(std::vector< DVector<T, 2> > &points) {
			std::copy(points.begin(), points.end(), std::back_inserter(m_polygonPoints));
		}
		/// @brief deallocate memory
		~DPolygon() {
			std::vector<PolygonPoint >().swap(m_polygonPoints);
			m_polygonPoints.clear();
		}
	public:
		/// @brief add a point to polygon
		///
		/// @return void
		void AddPoint(const PolygonPoint &_point) { m_polygonPoints.push_back(_point); }
		void AddPoint(const T &x, const T &y) { PolygonPoint _pp(x, y); m_polygonPoints.push_back(_pp); }
	public:
		bool IsSimple() const { ; }
		bool IsConvex() const { ; }
		bool IsCounterClockWise() const { ; }
		bool IsPointInPolygon() const { ; }
	public:
		//多边形的核
		//多边形的三角剖分
	protected:
		typedef DVector<T, 2> PolygonPoint;
		std::vector<PolygonPoint > m_polygonPoints;
	};

};

#endif// 2018/09/17