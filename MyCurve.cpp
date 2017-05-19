// MyCurve.cpp: implementation of the MyCurve class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Curve.h"
#include "matrix.h"
#include "MyCurve.h"


#include <fstream>
/*
 *	Template for using matrix.h.
 */
#ifndef _NO_NAMESPACE 
using namespace std;
using namespace math;
#define STD std
#else
#define STD
#endif

#define DIVISIONS 20
#define GL_PI  3.1415926535f

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

MyCurve::MyCurve()
{
	totalPoints = 0;
	picked = NULL;
	style = 0;
	showCtrl = false;
}

MyCurve::~MyCurve()
{

}

void MyCurve::AddPoint(float x, float y)
{
	totalPoints++;
	Point tmp(x, y);
	interpPoints.push_back(tmp);
	float n;
	if (totalPoints >= 2){
		//If there are more than 1 interpolation point, set up the 2 end points to help determine the curve.
		//They lie on the tangent of the first and last interpolation points.
		tmp = interpPoints[0] - interpPoints[1];
		n = tmp.Norm();
		endPoints[0] = interpPoints[0] + tmp / n * 50;
		tmp = interpPoints[totalPoints-1] - interpPoints[totalPoints-2];
		n = tmp.Norm();
		endPoints[1] = interpPoints[totalPoints-1] + tmp / n * 50;
	}
	
}

void MyCurve::DrawCircle(Point p)
{
	int vertices = 20;
	double r = 5.0;
	double tmpX,tmpY;
	double Angle, Angle0;
	
	Angle = (2*GL_PI)/vertices;
	Angle0 = 0.0;
	
	tmpX=p.x;
	tmpY=p.y;
	
	glBegin(GL_POLYGON);
	
	for (int i = 0; i < vertices; i++) {
		glVertex2d(tmpX + r * cos(i*Angle+Angle0), tmpY + r * sin(i*Angle+Angle0));
	}
	glEnd();
}

void MyCurve::PickPoint(float x, float y)
{
	float radius = 5.0;
	picked = NULL;
	Point tmp = Point(x, y);
	if (dist(tmp, endPoints[0]) < radius){
		picked = endPoints;
		return;
	}
	if (dist(tmp, endPoints[1]) < radius){
		picked = endPoints+1;
		return;
	}
	for (int i = 0; i < interpPoints.size(); i++){
		if (dist(tmp, interpPoints[i]) < radius){
			picked = &(interpPoints[i]);
			return;
		}
	}
}

void MyCurve::MovePicked(float x, float y)
{
	if (picked != NULL){
		(*picked).x = x;
		(*picked).y = y;
	}
	
}

void MyCurve::ClearAll()
{
	totalPoints = 0;
	interpPoints.clear();
	ctrlPoints.clear();
	curve.clear();
}

void MyCurve::DrawCurve()
{
	int i;
	if (totalPoints <= 0)
		return;
	if (totalPoints == 1){//if there is only 1 interpolation point, draw it.
		glColor3f(0,0,1);
		DrawCircle(interpPoints[0]);
	}else{//if there are more than 1 point, draw the curve.
		if (style != BSPLINE && style != HERMITE){
			glColor3f(0.0,1.0,0.0);
			//Draw the two end points.
			DrawCircle(endPoints[0]);
			DrawCircle(endPoints[1]);
			//Connect end points with interpolation points with straight lines.
			glColor3f(0.0,1.0,0.0);
			glBegin(GL_LINES);
			glVertex2d(endPoints[0].x, endPoints[0].y);
			glVertex2d(interpPoints[0].x, interpPoints[0].y);		
			glVertex2d(interpPoints[totalPoints - 1].x, interpPoints[totalPoints - 1].y);
			glVertex2d(endPoints[1].x, endPoints[1].y);
			glEnd();
		}
		//Calculate control points
		ControlPoints();
		//Interpolate the curve
		Interpolate();
		//Draw the curve
		glColor3f(1.0, 0.0, 0.0);
		for (i = 0; i < int(curve.size() - 1); i++){
			glBegin(GL_LINES);
			glVertex2d(curve[i].x, curve[i].y);
			glVertex2d(curve[i + 1].x, curve[i + 1].y);
			glEnd();
		}
		
		//Draw interpolation points.
		glColor3f(0.0,0.0,1.0);
		for (i = 0; i < totalPoints; i++){ 
			DrawCircle(interpPoints[i]);
		}
		//Draw control points and lines that connect them
		if (showCtrl && ctrlPoints.size() > 0){
			glColor3f(1.0, 1.0, 0.0);
			//B-spline, draw its control points
			if (style == BSPLINE){
				for (i = 0; i < totalPoints + 1; i++){
					DrawCircle(ctrlPoints[i]);
					glBegin(GL_LINES);
					glVertex2d(ctrlPoints[i].x, ctrlPoints[i].y);
					glVertex2d(ctrlPoints[i + 1].x, ctrlPoints[i + 1].y);
					glEnd();
				}
				DrawCircle(ctrlPoints[totalPoints + 1]);
			}else
				if (style == HERMITE){
					for (i = 0; i < totalPoints; i++)
					{
						DrawCircle(interpPoints[i] + ctrlPoints[i]);
						glBegin(GL_LINES);
						glVertex2d(interpPoints[i].x, interpPoints[i].y);
						glVertex2d(interpPoints[i].x + ctrlPoints[i].x, interpPoints[i].y + ctrlPoints[i].y);
						glEnd();
					}
				}
				else{
				//Bezier curve, draw its control points
					for (i = 0; i < totalPoints - 1; i++){
						DrawCircle(ctrlPoints[i * 2]);
						DrawCircle(ctrlPoints[i * 2 + 1]);
						glBegin(GL_LINES);
						glVertex2d(interpPoints[i].x, interpPoints[i].y);
						glVertex2d(ctrlPoints[i * 2].x, ctrlPoints[i * 2].y);
						
						glVertex2d(ctrlPoints[i * 2].x, ctrlPoints[i * 2].y);
						glVertex2d(ctrlPoints[i * 2 + 1].x, ctrlPoints[i * 2 + 1].y);
						
						glVertex2d(ctrlPoints[i * 2 + 1].x, ctrlPoints[i * 2 + 1].y);
						glVertex2d(interpPoints[i + 1].x, interpPoints[i + 1].y);
						glEnd();
					}
			}			
		}				
	}
}

void MyCurve::Interpolate()
{
	//Clear the old curve points
	curve.clear();
	//Depending on the selected style, interpolate the curve.
	switch(style){
		case BERSTEIN:	InterpBerstein(); 
			break;
		case CASTELJAU:	InterpCasteljau(); 
			break;
		case MATRIX:	InterpMatrix(); 
			break;
		case HERMITE:   InterpHermite(); 
			break;
	}
}

//////////////////////////////////////////////////////////////////////////
// Calculate the control points
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables 
// interpPoints	- type: vector<Point>
//				  discription: stores all the interpolation points
// endPoints	- type: Point[2]
//				  discription: stores the two end points tangent to the first and last interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// ctrlPoints	- type: vector<Point>
//				  discription: stores all the control points for curve interpolation and display.
//                             For Bezier curve, between very pair of consecutive interpolation points,
//                             there should be two control points. These four points determins the curve interpolation.
//                             For B-Spline, there should be totalPoints + 2 control points calculated from Ac = p.
// Hint: If you want to implement B-Spline, you need to write functions to create the A matrix as in the handouts.
//       Then you solve a linear system Ac = p, where p is the interpolation points vector and c are the control points.
//       We have provided you with a datastructure to store and solve the linear system.
//       Below is an example code, read the understand it.
//
//	matrix<float> A(3,3);
//  matrix<float> c(3,1);
//  matrix<float> p(3,1);
//  A(0,0) = 1.0; A(0,1) = 0.0; A(0,2) = 0.0;
//  A(1,0) = 0.0; A(1,1) = 1.0; A(1,2) = 0.0;
//  A(2,0) = 0.0; A(2,1) = 0.0; A(2,2) = 1.0;
//  p(0,0) = 1.0; p(1,0) = 2.0; p(2,0) = 3.0;
//  c = A.Solve(p);
//
//  The result in c is c(0,0) = 1.0; c(1,0) = 2.0; c(3,0) = 3.0, which satisfies Ac = p.



void MyCurve::ControlPoints(){
	// Prepare data
	if (totalPoints < 2)
		return;
	Point tan1, tan2, tmp;
	ctrlPoints.clear();
	int degree = 3;
	int dim;
	if (style == BSPLINE)
		dim = totalPoints + 2;
	else
		dim = totalPoints;
	matrix<double> A(dim, dim);
	matrix<double> C(dim, 2);
	matrix<double> P(dim, 2);
			
	switch(style)
	{
		case BERSTEIN:
		case CASTELJAU:
		case MATRIX:
			for (int i = 0; i < totalPoints-1; i++){
				if (i == 0)
					tan1 = interpPoints[0] + (interpPoints[0]-endPoints[0]);
				else
					tan1 = interpPoints[i] + (interpPoints[i+1]-interpPoints[i-1])/6;
				if (i == totalPoints - 2)
					tan2 = interpPoints[i+1] + (interpPoints[i+1]-endPoints[1]);
				else
					tan2 = interpPoints[i+1] + (interpPoints[i]-interpPoints[i+2])/6;
					ctrlPoints.push_back(tan1);
					ctrlPoints.push_back(tan2);
			}
			break;
		
		case HERMITE:
			//Additional constraints: second derivative are zero at end points
			//first row	
			A(0,0) = 2; A(0,1) = 1;
			for (int i = 2; i < dim; i++){
				A(0, i) = 0;
			}

			//middle rows
			for (int i = 1; i < dim - 1; i++){
				for (int j = 0; j < dim; j++)
					A(i,j) = 0;
				A(i, i-1) = 1; A(i, i) = 4; A(i, i+1) = 1;
			}

			//last row
			for (int i = 0; i < dim - 2; i++)
				A(dim-1, i) = 0;
			A(dim-1, dim-2) = 1; A(dim-1,dim-1) = 2;


			P(0,0) = 3 * (interpPoints[1].x - interpPoints[0].x);
			P(0,1) = 3 * (interpPoints[1].y - interpPoints[0].y);
			for (int i = 1; i < dim - 1; i++){
				P(i, 0) = 3 * (interpPoints[i+1].x - interpPoints[i-1].x);
				P(i, 1) = 3 * (interpPoints[i+1].y - interpPoints[i-1].y);
			}
			P(dim-1,0) =  3 * (interpPoints[dim-1].x - interpPoints[dim-2].x);
			P(dim-1,1) =  3 * (interpPoints[dim-1].y - interpPoints[dim-2].y);


			C = A.Solve(P);
			for (int i = 0; i < totalPoints; i++){
				ctrlPoints.push_back(Point(C(i,0), C(i,1)));
			}

			// Based on showCtrl, determine the boundary slope either by endPoints or automatically
			if (showCtrl){
				ctrlPoints[0] = endPoints[0] - interpPoints[0];
				ctrlPoints[totalPoints - 1] = endPoints[1] - interpPoints[totalPoints - 1];
			}		
			else{
				endPoints[0] = interpPoints[0] + ctrlPoints[0];
				endPoints[1] = interpPoints[totalPoints - 1] + ctrlPoints[totalPoints - 1];
			}
	
			break;

		}
		
}	


//////////////////////////////////////////////////////////////////////////
// Cubic Berstein Bezier Spline
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables 
// interpPoints	- type: vector<Point> -> user points
//				  discription: stores all the interpolation points
// ctrlPoints	- type: vector<Point> -> creating the control points
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// endPoints	- type: Point[2] ->    
//				  discription: stores the two end points tangent to the first and last interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curve		- type: vector<Point>
//				  discription: stores all the points that form the curve, including all interpolation points 



/* interpoation point for Berstein*/
void MyCurve::InterpBerstein(){

	int coeff[4] = {1, 3, 3, 1};

	curve.clear();
	//since you won't calc for last interpolation point/given point, totalPoints-1
	for (int i = 0; i < totalPoints-1; i++){
		vector<Point> p;
		p.push_back(interpPoints[i]);
		p.push_back(ctrlPoints[2*i]); // first ctrl point between interpPoint i and i+1
		p.push_back(ctrlPoints[2*i+1]); // second ctrl point between interpPoint i and i+1
		p.push_back(interpPoints[i+1]);

		// now calc berstein based these 4 points

		for(float k=0; k < DIVISIONS; k++){
			float t = (float)k/(float)DIVISIONS;// ei value ta tumi set koiro
			Point curve_point;
			for(int j=0; j<=3; j++){
				float x = coeff[j] * pow((1-t), (3-j)) * pow(t, j) * p[j].x;
				float y = coeff[j] * pow((1-t), (3-j)) * pow(t, j) * p[j].y;
				curve_point.x += x; // slide er calc, tumi syntax thik kore nio. 
				curve_point.y += y;
		
			}
			
			curve.push_back(curve_point);
		}

	}
	
}

//////////////////////////////////////////////////////////////////////////
// Cubic de Casteljau Bezier Spline
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables 
// interpPoints	- type: vector<Point>
//				  discription: stores all the interpolation points
// ctrlPoints	- type: vector<Point>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// endPoints	- type: Point[2]
//				  discription: stores the two end points tangent to the first and last interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curve		- type: vector<Point>
//				  discription: stores all the points that form the curve, including all interpolation points 
void MyCurve::InterpCasteljau(){
	
	curve.clear();
	//since you won't calc for last interpolation point/given point, totalPoints-1
	for (int i = 0; i < totalPoints-1; i++){
		vector<Point> p;
		p.push_back(interpPoints[i]);
		p.push_back(ctrlPoints[i*2]); // first ctrl point between interpPoint i and i+1
		p.push_back(ctrlPoints[i*2+1]); // second ctrl point between interpPoint i and i+1
		p.push_back(interpPoints[i+1]);
		
		for(int j = 0; j < DIVISIONS; j++){
			vector<Point> q;
			float t = (float)j/(float)DIVISIONS;
			q.push_back(p[0]);
			q.push_back(p[1]);
			q.push_back(p[2]);
			q.push_back(p[3]);

			//now calc point for each value of t, using castaljau algorithm

			//algorithm taken from:
			//http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/de-casteljau.html
			int n = 4; // 4 points for each segment
			//vector<Point> q;
			for(int k = 1; k <= n; k++){
				for(int ii = 0; ii < n - k; ii++){
					q[ii].x = (1 - t)*q[ii].x + t*q[ii + 1].x; 
					q[ii].y = (1 - t)*q[ii].y + t*q[ii + 1].y; 
				}
			}
			curve.push_back(q[0]);
		}	
	}
	
}

//////////////////////////////////////////////////////////////////////////
// Cubic Matrix Form Bezier Spline
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables 
// interpPoints	- type: vector<Point>
//				  discription: stores all the interpolation points
// ctrlPoints	- type: vector<Point>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// endPoints	- type: Point[2]
//				  discription: stores the two end points tangent to the first and last interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curve		- type: vector<Point>
//				  discription: stores all the points that form the curve, including all interpolation points 
void MyCurve::InterpMatrix(){
	
	curve.clear();
	
	//since you won't calc for last interpolation point/given point, totalPoints-1
	
	for (int i = 0; i < totalPoints-1; i++){
		
		// now calc matrix based 

		for(int k=0; k < DIVISIONS; k++){
			Point p0, p1, p2, p3;
			float t = (float)k/(float)DIVISIONS;// ei value ta tumi set koiro
			float t3 = t*t*t;
			float t2 = t*t;

			float b0 = -t*t*t + 3*t*t -3*t + 1;
			float b1 = 3*t*t*t - 6*t*t + 3*t;
			float b2 = -3*t*t*t + 3*t*t;
			float b3 = t*t*t;

			p0.x = b0 * interpPoints[i].x, p0.y = b0 * interpPoints[i].y;
			p1.x = b2 * ctrlPoints[2*i].x, p1.y = b2 * ctrlPoints[2*i].y;
			p2.x = b3 * ctrlPoints[2*i+1].x, p2.y = b3 * ctrlPoints[2*i+1].y;
			p3.x = b1 * interpPoints[i+1].x, p3.y = b1 * interpPoints[i+1].y;

			Point P1; 
			P1.x = p0.x + p1.x + p2.x + p3.x;
			P1.y = p0.y + p1.y + p2.y + p3.y;
			curve.push_back(P1);
				
		}
		
	}

}
	

//////////////////////////////////////////////////////////////////////////
// Bonus Points: Hermite Spline curve
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables 
// interpPoints	- type: vector<Point>
//				  discription: stores all the interpolation points
// ctrlPoints	- type: vector<Point>
//				  discription: stores the control points that helps to determine the curve.
//                             There should be totalPoints control points.
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curve		- type: vector<Point>
//				  discription: stores all the points that form the curve, including all interpolation points
void MyCurve::InterpHermite(){
	
	curve.clear();	
	for (int i =0; i< totalPoints-1; i++){
		
		for(int k=0; k < DIVISIONS; k++){
			Point P0, P1, T1, T2;
			Point p;
			float t = (float)k/(float)DIVISIONS;// ei value ta tumi set koiro
		
			float b0 =  2*t*t*t - 3*t*t + 1;
			float b1 = -2*t*t*t + 3*t*t;
			float b2 = t*t*t - 2*t*t + t;
			float b3 = t*t*t - t*t;

			P0.x = b0 * interpPoints[i].x, P0.y = b0 * interpPoints[i].y;
			P1.x = b1 * interpPoints[i+1].x, P1.y = b1 * interpPoints[i+1].y;
			T1.x = b2 * ctrlPoints[i].x, T1.y = b2 * ctrlPoints[i].y;
			T2.x = b3 * ctrlPoints[i+1].x, T2.y = b3 * ctrlPoints[i+1].y;

			p.x = P0.x + P1.x + T1.x + T2.x;
			p.y = P0.y + P1.y + T1.y + T2.y;

			curve.push_back(p);
	
		
		}

	}
}

