#ifndef __Particles_h__
#define __Particles_h__
#include "Common.h"

template<int d> class Particles
{using VectorD=Vector<real,d>;
using MatrixD=Matrix<real,d>;
public:
	////attributes
	ArrayPtr<VectorD> x;		////position
	ArrayPtr<VectorD> v;		////velocity
	ArrayPtr<VectorD> f;		////force
	ArrayPtr<real> m;			////mass
	ArrayPtr<real> den;			////density
	ArrayPtr<int> idx;			////index, for rigid body
	ArrayPtr<real> vol;			////particle volume
	ArrayPtr<MatrixD> f_E;		////Elestic deformation gradient
	ArrayPtr<MatrixD> f_P;		////Plastic deformation gradient
	ArrayPtr<VectorD> delta_f;  ////Delta_f
	ArrayPtr<MatrixD> Ap;		////Plastic deformation gradient

	//////////////////////////////////////////////////////////////////////////
	////common functions
	Particles()
	{
		if(x==nullptr)x.reset(new Array<VectorD>());	
		if(v==nullptr)v.reset(new Array<VectorD>());	
		if(f==nullptr)f.reset(new Array<VectorD>());	
		if(m==nullptr)m.reset(new Array<real>());		
		if(den==nullptr)den.reset(new Array<real>());	
		if(idx==nullptr)idx.reset(new Array<int>());
		if(vol==nullptr)vol.reset(new Array<real>());
		if(f_E==nullptr)f_E.reset(new Array<MatrixD>());
		if(f_P==nullptr)f_P.reset(new Array<MatrixD>());
		if(Ap==nullptr)Ap.reset(new Array<MatrixD>());
		if(delta_f==nullptr)delta_f.reset(new Array<VectorD>());
	}

	void Resize(const int size)
	{
		x->resize((size_type)size,VectorD::Zero());
		v->resize((size_type)size,VectorD::Zero());
		f->resize((size_type)size,VectorD::Zero());
		m->resize((size_type)size,(real)0);
		den->resize((size_type)size,(real)0);
		idx->resize((size_type)size,0);
		vol->resize((size_type)size,0);
		f_E->resize((size_type)size,MatrixD::Identity());
		f_P->resize((size_type)size,MatrixD::Identity());
		Ap->resize((size_type)size,MatrixD::Zero());
		delta_f->resize((size_type)size,VectorD::Zero());
	}

	int Add_Element()
	{
		x->push_back(VectorD::Zero());
		v->push_back(VectorD::Zero());
		f->push_back(VectorD::Zero());
		m->push_back((real)0);
		den->push_back((real)0);
		idx->push_back(0);
		vol->push_back(0);
		f_E->push_back(MatrixD::Identity());
		f_P->push_back(MatrixD::Identity());
		Ap->push_back(MatrixD::Zero());
		delta_f->push_back(VectorD::Zero());
		return (int)x->size()-1;
	}

	int Size() const {return (int)(*x).size();}

	//////////////////////////////////////////////////////////////////////////
	////functions for separate attributes
	////functions for x
	VectorD& X(const int i)
	{return (*x)[i];}
	
	const VectorD& X(const int i) const 
	{return (*x)[i];}

	Array<VectorD>* X()
	{return x.get();}

	const Array<VectorD>* X() const 
	{return x.get();}
	
	ArrayPtr<VectorD> XPtr()
	{return x;}
	
	const ArrayPtr<VectorD> XPtr() const
	{return x;}
	
	Array<VectorD>& XRef()
	{return *x;}

	const Array<VectorD>& XRef() const 
	{return *x;}

	////functions for x
	VectorD& Delta_F(const int i)
	{return (*delta_f)[i];}
	
	const VectorD& Delta_F(const int i) const 
	{return (*delta_f)[i];}

	Array<VectorD>* Delta_F()
	{return delta_f.get();}

	const Array<VectorD>* Delta_F() const 
	{return delta_f.get();}
	
	ArrayPtr<VectorD> Delta_FPtr()
	{return delta_f;}
	
	const ArrayPtr<VectorD> Delta_FPtr() const
	{return delta_f;}
	
	Array<VectorD>& Delta_FRef()
	{return *delta_f;}

	const Array<VectorD>& Delta_FRef() const 
	{return *delta_f;}
	
	//////////////////////////////////////////////////////////////////////////
	////functions for v
	VectorD& V(const int i)
	{return (*v)[i];}
	
	const VectorD& V(const int i) const 
	{return (*v)[i];}

	Array<VectorD>* V()
	{return v.get();}

	const Array<VectorD>* V() const 
	{return v.get();}
	
	ArrayPtr<VectorD> VPtr()
	{return v;}
	
	const ArrayPtr<VectorD> VPtr() const
	{return v;}
	
	Array<VectorD>& VRef()
	{return *v;}

	const Array<VectorD>& VRef() const 
	{return *v;}

	//////////////////////////////////////////////////////////////////////////
	////functions for f
	VectorD& F(const int i)
	{return (*f)[i];}
	
	const VectorD& F(const int i) const 
	{return (*f)[i];}

	Array<VectorD>* F()
	{return f.get();}

	const Array<VectorD>* F() const 
	{return f.get();}
	
	ArrayPtr<VectorD> FPtr()
	{return f;}
	
	const ArrayPtr<VectorD> FPtr() const
	{return f;}
	
	Array<VectorD>& FRef()
	{return *f;}

	const Array<VectorD>& FRef() const 
	{return *f;}

	//////////////////////////////////////////////////////////////////////////
	////functions for m
	real& M(const int i)
	{return (*m)[i];}
	
	const real& M(const int i) const 
	{return (*m)[i];}

	Array<real>* M()
	{return m.get();}

	const Array<real>* M() const 
	{return m.get();}
	
	ArrayPtr<real> MPtr()
	{return m;}
	
	const ArrayPtr<real> MPtr() const
	{return m;}
	
	Array<real>& MRef()
	{return *m;}

	const Array<real>& MRef() const 
	{return *m;}
	
	//////////////////////////////////////////////////////////////////////////
	////functions for f_E
	MatrixD& F_E(const int i)
	{return (*f_E)[i];}
	
	const MatrixD& F_E(const int i) const 
	{return (*f_E)[i];}

	Array<MatrixD>* F_E()
	{return f_E.get();}

	const Array<MatrixD>* F_E() const 
	{return f_E.get();}
	
	ArrayPtr<MatrixD> F_EPtr()
	{return f_E;}
	
	const ArrayPtr<MatrixD> F_EPtr() const
	{return f_E;}
	
	Array<MatrixD>& F_ERef()
	{return *f_E;}

	const Array<MatrixD>& F_ERef() const 
	{return *f_E;}

	//////////////////////////////////////////////////////////////////////////
	////functions for Ap
	MatrixD& A(const int i)
	{return (*Ap)[i];}
	
	const MatrixD& A(const int i) const 
	{return (*Ap)[i];}

	Array<MatrixD>* A()
	{return Ap.get();}

	const Array<MatrixD>* A() const 
	{return Ap.get();}
	
	ArrayPtr<MatrixD> APtr()
	{return Ap;}
	
	const ArrayPtr<MatrixD> APtr() const
	{return Ap;}
	
	Array<MatrixD>& ARef()
	{return *Ap;}

	const Array<MatrixD>& ARef() const 
	{return *Ap;}

	//////////////////////////////////////////////////////////////////////////
	////functions for f_P
	MatrixD& F_P(const int i)
	{return (*f_P)[i];}
	
	const MatrixD& F_P(const int i) const 
	{return (*f_P)[i];}

	Array<MatrixD>* F_P()
	{return f_P.get();}

	const Array<MatrixD>* F_P() const 
	{return f_P.get();}
	
	ArrayPtr<MatrixD> F_PPtr()
	{return f_P;}
	
	const ArrayPtr<MatrixD> F_PPtr() const
	{return f_P;}
	
	Array<MatrixD>& F_PRef()
	{return *f_P;}

	const Array<MatrixD>& F_PRef() const 
	{return *f_P;}

	//////////////////////////////////////////////////////////////////////////
	////functions for vol
	real& Vol(const int i)
	{return (*vol)[i];}
	
	const real& Vol(const int i) const 
	{return (*vol)[i];}

	Array<real>* Vol()
	{return vol.get();}

	const Array<real>* Vol() const 
	{return vol.get();}
	
	ArrayPtr<real> VolPtr()
	{return vol;}
	
	const ArrayPtr<real> VolPtr() const
	{return vol;}
	
	Array<real>& VolRef()
	{return *vol;}

	const Array<real>& VolRef() const 
	{return *vol;}

	//////////////////////////////////////////////////////////////////////////
	////functions for den
	real& D(const int i)
	{return (*den)[i];}
	
	const real& D(const int i) const 
	{return (*den)[i];}

	Array<real>* D()
	{return den.get();}

	const Array<real>* D() const 
	{return den.get();}
	
	ArrayPtr<real> DPtr()
	{return den;}
	
	const ArrayPtr<real> DPtr() const
	{return den;}
	
	Array<real>& DRef()
	{return *den;}

	const Array<real>& DRef() const 
	{return *den;}

	//////////////////////////////////////////////////////////////////////////
	////functions for idx
	int& I(const int i)
	{return (*idx)[i];}
	
	const int& I(const int i) const 
	{return (*idx)[i];}

	Array<int>* I()
	{return idx.get();}

	const Array<int>* I() const 
	{return idx.get();}
	
	ArrayPtr<int> IPtr()
	{return idx;}
	
	const ArrayPtr<int> IPtr() const
	{return idx;}
	
	Array<int>& IRef()
	{return *idx;}

	const Array<int>& IRef() const 
	{return *idx;}
};
#endif
