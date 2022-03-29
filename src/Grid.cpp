#include "Grid.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
////Initializer
template<int d> Grid<d>::Grid(const VectorDi& _cell_counts,const real _dx,const VectorD& _domain_min)
:cell_counts(_cell_counts),node_counts(_cell_counts+VectorDi::Ones()),dx(_dx),domain_min(_domain_min){domain_max=domain_min+cell_counts.template cast<real>()*dx;}

template<int d> Grid<d>& Grid<d>::operator=(const Grid& copy)
{cell_counts=copy.cell_counts;node_counts=copy.node_counts;dx=copy.dx;domain_min=copy.domain_min;domain_max=copy.domain_max;return *this;}

template<int d> void Grid<d>::Initialize(const VectorDi& _cell_counts,const real _dx,const VectorD& _domain_min)
{cell_counts=_cell_counts;node_counts=cell_counts+VectorDi::Ones();dx=_dx;domain_min=_domain_min;domain_max=domain_min+cell_counts.template cast<real>()*dx;}

////Coordinate operations
template<int d> Vector<int,d> Grid<d>::Coord(const int index,const VectorDi& counts) {/*not impl*/return Vector<int,d>::Zero();}
template<> Vector<int,2> Grid<2>::Coord(const int index,const Vector<int,2>& counts) {return Vector2i(index/counts[1],index%counts[1]);}
template<> Vector<int,3> Grid<3>::Coord(const int index,const Vector<int,3>& counts) {int n_yz=counts[1]*counts[2];int mod_yz=index%n_yz;return Vector3i(index/n_yz,mod_yz/counts[2],mod_yz%counts[2]);}

template<int d> Vector<int,d> Grid<d>::Node_Coord(const int index) const {return Coord(index,node_counts);}

template<int d> Vector<int,d> Grid<d>::Cell_Coord(const int index) const {return Coord(index,cell_counts);}

template<int d> Vector<int,d> Grid<d>::Cell_Coord(const VectorD& pos) const {VectorD coord_with_frac=(pos-domain_min)/dx;return coord_with_frac.template cast<int>();}

////////////////////////////////////////////////////////////////////////////////////////////////////
////Position operations
template<int d> Vector<real,d> Grid<d>::Node(const VectorDi& node) const 
{return domain_min+node.template cast<real>()*dx;}

template<int d> Vector<real,d> Grid<d>::Center(const VectorDi& cell) const 
{return domain_min+(cell.template cast<real>()+(real).5*VectorD::Ones())*dx;}

template<int d> Vector<real,d> Grid<d>::Node(const int node_index) const
{return Node(Node_Coord(node_index));}

template<int d> Vector<real,d> Grid<d>::Center(const int cell_index) const
{return Center(Cell_Coord(cell_index));}

////////////////////////////////////////////////////////////////////////////////////////////////////
////Index operations
template<int d> int Grid<d>::Index(const Vector<int,d>& coord,const Vector<int,d>& counts){/*not impl*/return Vector<int,d>::Zero();}
template<> int Grid<1>::Index(const Vector<int,1>& coord,const Vector<int,1>& counts){return coord[0];}
template<> int Grid<2>::Index(const Vector<int,2>& coord,const Vector<int,2>& counts){return coord[0]*counts[1]+coord[1];}
template<> int Grid<3>::Index(const Vector<int,3>& coord,const Vector<int,3>& counts){return coord[0]*counts[1]*counts[2]+coord[1]*counts[2]+coord[2];}

template<int d> int Grid<d>::Cell_Index(const VectorDi& cell) const {return Index(cell,cell_counts);}

template<int d> int Grid<d>::Node_Index(const VectorDi& node) const {return Index(node,node_counts);}

////////////////////////////////////////////////////////////////////////////////////////////////////
////Validate index and coord
template<class ArrayT> bool All_Less(const ArrayT& a0,const ArrayT& a1){for(auto i=0;i<a0.size();i++){if(a0[i]>=a1[i])return false;}return true;}
template<class ArrayT> bool All_Greater_Equal(const ArrayT& a0,const ArrayT& a1){for(auto i=0;i<a0.size();i++){if(a0[i]<a1[i])return false;}return true;}

template<int d> bool Grid<d>::Valid(const VectorDi& coord,const VectorDi& counts)
{static const VectorDi zero=VectorDi::Zero();return All_Greater_Equal(coord,zero)&&All_Less(coord,counts);}

template<int d> bool Grid<d>::Valid_Node(const VectorDi& node) const {return Valid(node,node_counts);}

template<int d> bool Grid<d>::Valid_Cell(const VectorDi& cell) const {return Valid(cell,cell_counts);}

template class Grid<2>;
template class Grid<3>;