// PFgrid.h
// Class definitions used in sparse grid scheme for phase field
// Questions/commentS to jgruber@andrew.cmu.edu
//
// PFgrid objects are declared with the following syntax (note
// that "double" here could also be "float", etc.):
//
//		PFgrid1D gridname1(xmax);
//		PFgrid2D gridname2(xmax,ymax);
//		PFgrid3D gridname3(xmax,ymax,zmax);
//
//
// Order parameter values can be used as lvalues with the syntax
//
//		grid[x][y][z][nu] = somevalue;
//
//
// Order parameter values can be used as rvalues with the syntax
//
//		somevalue = grid[x][y][z][nu];
//
//
// Most calculations are most efficient when performed using only 
// the nonzero order parameters at a point.  The function nonzero()
// can be used to determine the correct end condition for the loop,
// and the functions value() and index() can be used to set or get
// information, e.g.
//
//		int stop = grid[x][y][z].nonzero();
//		for (int i=0; i<stop; i++) {
//			double value = grid[x][y][z].value(i);
//			int index = grid[x][y][z].index(i);
//			// use "value" and "index" in some computation
//		}
//
//
// Note that the argument to either value() or index() is the 
// position within the PFdata vector
//
// An even more efficient use of this idea is to declare a reference
// to a PFdata, PFgrid1D, or PFgrid2D class and use this in the loop,
// thus avoiding unnecessary subscript operator calls.
//

#ifndef PFGRID
#define PFGRID
#include<vector>
#include<stdlib.h>

#define EPSILON (5.0E-3)



class PFdata{
	struct PFitem{double value; int index;};
public:
	PFdata() {}
	double& operator[](int index);
	const double operator[](int index) const;
	const int nonzero() const {return v.size();}
	double& value(int i) {return v[i].value;}
	const double value(int i) const {return v[i].value;}
	int& index(int i) {return v[i].index;}
	const int index(int i) const {return v[i].index;}
	void clear() {v.clear();}
private:
	std::vector<PFitem> v;
};

inline double& PFdata::operator[](int index)
{
	std::vector<PFitem>::iterator it;
	std::vector<PFitem>::iterator end = v.end();
	for (it=v.begin(); it!=end; ++it)
		if ((*it).index==index) return (*it).value;

	if (v.size()>0) {
		PFitem& back = v.back();
		if (back.value<EPSILON) {
			back.index = index;
			return back.value;
		}
	}

	PFitem item;
	item.value = 0.0;
	item.index = index;
	v.push_back(item);
	return v.back().value;	
}

inline const double PFdata::operator[](int index) const
{
	std::vector<PFitem>::const_iterator it;
	std::vector<PFitem>::const_iterator end = v.end();
	for (it=v.begin(); it!=end; ++it)
		if ((*it).index==index) return (*it).value;
	return 0.0;
}


class PFgrid1D{
public:
	PFgrid1D(int x) {v.assign(x,PFdata());}
	PFdata& operator[](int i) {return v[i];}
	const PFdata& operator[](int i) const {return v[i];}
	void swap(PFgrid1D& grid) {v.swap(grid.v);}
private:
	std::vector<PFdata > v;
};

class PFgrid2D{
public:
	PFgrid2D(int y, int x) {v.assign(y,PFgrid1D(x));}
	PFgrid1D& operator[](int i) {return v[i];}
	const PFgrid1D& operator[](int i) const {return v[i];}
	void swap(PFgrid2D& grid) {v.swap(grid.v);}
private:
	std::vector<PFgrid1D > v;
};


class PFgrid3D{
public:
	PFgrid3D(int x, int y, int z) {v.assign(x,PFgrid2D(y,z));}
	PFgrid2D& operator[](int i) {return v[i];}
	const PFgrid2D operator[](int i) const {return v[i];}
	void swap(PFgrid3D& grid) {v.swap(grid.v);}
private:
	std::vector<PFgrid2D > v;
};

#endif
