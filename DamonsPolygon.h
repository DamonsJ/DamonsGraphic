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
		void AddPoint(const DVector<T, 2> &_point) { m_polygonPoints.push_back(_point); }
		void AddPoint(const T &x, const T &y) { PolygonPoint _pp(x, y); m_polygonPoints.push_back(_pp); }
	public:
		/// @brief get point size of polygon
		/// @note only count point size,if first point same with last point
		/// only one counted
		/// @return unsigned int : point size 
		unsigned int GetPolygonPointSize() const{
			return IsClosed()? m_polygonPoints.size() - 1:m_polygonPoints.size();
		}
		/// @brief close polygon 
		/// @return void
		void ClosePolygon() {
			if(!IsClosed()){
				auto p = m_polygonPoints.front();
				m_polygonPoints.push_back(p);
			}
		}
	public:
		/// @brief tell whether polygon is simple
		/// @note a polygon is simple when there are two edges intersected
		/// which means the intersected point is not end point.
		/// @return  true if simple otherwise false 
		bool IsSimple() { 
			ClosePolygon();
			unsigned int sz = GetPolygonPointSize();
			for(unsigned int i = 0;i < sz - 1;++i){
				if(IsSegmentIntersect(m_polygonPoints[i],m_polygonPoints[i+1],m_polygonPoints[i+1],m_polygonPoints[i+2]))
					return true;	
			}
			return false;
		}
		bool IsConvex() const { ; }
		bool IsCounterClockWise() const { ; }
		bool IsPointInPolygon() const { ; }
	public:
		//����εĺ�
		//����ε������ʷ�

	protected:
		/// @brief tell whether polygon is close
		/// @note a polygon is close when the first and last point is same
		/// @return  true if close otherwise false 
	    bool IsClosed() const{
			auto &p1 = m_polygonPoints.front();
			auto &p2 = m_polygonPoints.back();
			auto pd = p1 - p2;
			return pd.LengthSquared() < DEplision*DEplision;
		}
		/// @brief tell whether two segment is intersect
		/// @return  true if intersect otherwise false 		
		bool IsSegmentIntersect(DVector<T, 2> &V0,DVector<T, 2> &V1,DVector<T, 2> &U0,DVector<T, 2> &U1) {
				 							
				const T Ax = V1[0] - V0[0];				
				const T Ay = V1[1] - V0[1];
				const T Bx = U0[0] - U1[0];								
				const T By = U0[1] - U1[1];							
				const T Cx = V0[0] - U0[0];							
				const T Cy = V0[1] - U0[1];							
				const T f  = Ay*Bx - Ax*By;								
				const T d  = By*Cx - Bx*Cy;									
				if((f>0.0 && d>=0.0 && d<=f) || (f<0.0 && d<=0.0 && d>=f))	
				{													
					const T e = Ax*Cy - Ay*Cx;					
					if(f>0.0)										
					{												
						if(e>=0.0 && e<=f) return true;			
					}												
					else											
					{												
						if(e<=0.0 && e>=f) return true;			
					}												
				}

				return false;
		}
	protected:
		typedef DVector<T, 2> PolygonPoint;
		std::vector<PolygonPoint > m_polygonPoints;
	};

};

#endif// 2018/09/17