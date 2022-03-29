#ifndef __Grid_h__
#define __Grid_h__
#include "Common.h"

template<int d> class Grid
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;
public:
	VectorDi cell_counts;
	VectorDi node_counts;
	real dx;
	VectorD domain_min;
	VectorD domain_max;

	////Initializer
	Grid(const VectorDi& _cell_counts=VectorDi::Zero(),const real _dx=(real)0,const VectorD& _domain_min=VectorD::Zero());
	Grid& operator=(const Grid& copy);
	Grid(const Grid& copy){*this=copy;}
    void Initialize(const VectorDi& _cell_counts=VectorDi::Zero(),const real _dx=(real)0,const VectorD& _domain_min=VectorD::Zero());

	////Coordinate operations
    static VectorDi Coord(const int index,const VectorDi& counts);                      ////general index->coord
    VectorDi Node_Coord(const int index) const;			                                ////node index->node coord
	VectorDi Cell_Coord(const int index) const;                                         ////cell index->cell coord
	VectorDi Cell_Coord(const VectorD& pos) const;                                      ////pos->cell coord

    ////Position operations
    VectorD Node(const VectorDi& node) const;                                           ////node coord->node pos
    VectorD Center(const VectorDi& cell) const;                                         ////cell coord->cell center pos
	VectorD Node(const int node_index) const;
	VectorD Center(const int cell_index) const;

    ////Index operations
    static int Index(const Vector<int,d>& coord,const Vector<int,d>& counts);           ////general coord->index
    int Cell_Index(const VectorDi& cell) const;                                         ////cell coord->cell index
    int Node_Index(const VectorDi& node) const;                                         ////node coord->node index

    ////Validate index and coord
	static bool Valid(const VectorDi& coord,const VectorDi& counts);                    ////check if coord is within [zero,counts)
	bool Valid_Node(const VectorDi& node) const;                                        ////check if node is within [zero,node_counts)
    bool Valid_Cell(const VectorDi& cell) const;                                        ////check if node is within [zero,cell_counts)
};

#endif
