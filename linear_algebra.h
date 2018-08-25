//本模板可直接复制到.cpp文件中，.h指头文件。
// "linear_algebra.h" -*- C++ -*-
//*for linear_algebra*
//*has defined matrix*

#ifndef _LINEAR_ALGEBRA_H_
#define _LINEAR_ALGEBRA_H_ 1

//C++ includes
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <string>
#include <map>
#if __cplusplus>=201103L
#include <initializer_list> 
#endif

//C includes
#include <cmath>
#include <cstdio>

//some defines
#ifndef DEFINES_
#define DEFINES_ 1 
#define show_var(x) cout<<#x<<":"<<(x)
#define reg register

const int inf=0x3F3F3F3F;
const long long infl=0x3F3F3F3F3F3F3F3F;
const double e=pow(1+1.0/inf,inf),pi=std::acos(-1),eps=1e-6;

template<class mtype>mtype derivative(mtype func(mtype x),mtype xx){
	return (func(xx+eps)-func(xx))/eps;
}
#endif
//start
template<class mtype>//linear algebra
class matrix{
	public:
typedef size_t size_type;
typedef matrix<mtype> matrix_type;
typedef std::vector<mtype> vec_type;
#define rcm_type const matrix_type&
		explicit matrix(size_type h=2,size_type w=2,mtype init=zero
		):height(h),width(w){
			vec.resize(height);
			for(reg size_type i=0;i<height;++i)
				for(reg size_type j=0;j<width;++j)
					vec[i].push_back(init);
		}
		matrix(std::vector<vec_type> newvec
		):height(newvec.size()){
			vec=newvec;
			width=vec[0].size();
		}
#if __cplusplus>=201103L
		matrix(std::initializer_list<std::initializer_list<mtype> > ini)
		:height(ini.size()){
			vec.resize(height);
typedef typename std::initializer_list<std::initializer_list<mtype> >::iterator iite;
typedef typename std::initializer_list<mtype> mite;
		for(iite i=ini.begin();i<ini.end();++i)
			for(mite j=*i ->begin();j<*i ->end();++j)
				vec[i].push_back(*j);
		}
#endif
		matrix(rcm_type rhs){*this=rhs;}
		~matrix(){}
		matrix(mtype num,size_type s=2){*this=matrix_cast(num,s);}
		
		//operators
		vec_type& operator[](size_type index){return vec[index];}
		const vec_type& operator[](size_type index) const 
		{return vec[index];}
		matrix_type& operator=(rcm_type rhs){
			height=rhs.height;
			width=rhs.width;
			vec=rhs.vec;
			return *this;
		}
		matrix_type operator+(rcm_type rhs){matrix_type ans;
			if(height!=rhs.height||width!=rhs.width)
				throw "cannot_do_the_calc";
			ans.vec.resize(height);for(reg size_type i=0;i<height;++i)
				for(reg size_type j=0;j<width;++j)
					ans[i].push_back(vec[i][j]+rhs.vec[i][j]);
			return ans;
		}
		matrix_type operator-(rcm_type rhs){matrix_type ans;
			if(height!=rhs.height||width!=rhs.width)
				throw "cannot_do_the_calc";
			ans.vec.resize(height);for(reg size_type i=0;i<height;++i)
				for(reg size_type j=0;j<width;++j)
					ans[i].push_back(vec[i][j]-rhs.vec[i][j]);
			return ans;
		}
		matrix_type operator*(rcm_type rhs){matrix_type ans;
			if(width!=rhs.height)
				throw "cannot_do_the_calc";
			ans.vec.resize(height);
			for(reg size_type i=0;i<height;++i)ans.vec[i].resize(rhs.width);
			for(reg size_type i=0;i<height;++i)
				for(reg size_type j=0;j<rhs.width;++j)
					for(reg size_type k=0;k<width;k++)
						ans[i][j]+=vec[i][k]*rhs.vec[k][j];
			return ans;
		}
		matrix_type operator/(rcm_type rhs){
			if(rhs.det()==zero||width!=rhs.inv().height)
				throw "sorry, please check your matrixes!";matrix_type mat;
			return (*this)*(rhs.inv());
		}
		friend inline matrix_type operator*(mtype n,rcm_type mat){
			matrix_type ans(mat);
			for(reg size_type i=0;i<ans.height;++i)
				for(reg size_type j=0;j<ans.width;++j)
					ans[i][j]*=n;
			return ans;
		}
		inline matrix_type operator/(mtype n){
			matrix_type ans(*this);
			for(reg size_type i=0;i<ans.height;++i)
				for(reg size_type j=0;j<ans.width;++j)
					ans[i][j]/=n;
			return ans;
		}
		friend inline bool operator==(rcm_type lhs,rcm_type rhs){
			return lhs.vec==rhs.vec;
		}
		friend inline bool operator!=(rcm_type lhs,rcm_type rhs){
			return !(lhs==rhs);
		}
		friend inline std::ostream& operator<<(std::ostream& os,
		rcm_type mat){
			for(reg size_type i=0;i<mat.height;++i){
				for(reg size_type j=0;j<mat.width;++j)
					os<<mat[i][j]<<" ";puts("");}
			return os;
		}
		friend inline std::istream& operator>>(std::istream& is,
		matrix_type& mat){
			for(reg size_type i=0;i<mat.height;++i)
				for(reg size_type j=0;j<mat.width;++j)
					is>>mat[i][j];
			if(!is)mat=matrix_type(mat.height,mat.width);
			return is;
		}
		
		//some functions
		inline std::pair<size_type,size_type> size(){
			return std::make_pair(height,width);
		}
		inline matrix_type T(){matrix_type ans=*this;
			for(reg size_type i=0;i<height;++i)
				for(reg size_type j=0;j<width;++j)ans[j][i]=vec[i][j];
			return ans;
		}
		matrix_type inv(){return ortmat()/det();}//wrong answer(from ortmat())
		mtype det(){mtype ans=zero;
			if(height!=width)throw "no_dets";
			if(height==1&&width==1)return vec[0][0];
			for(reg size_type j=0;j<width;j++)
				ans+=vec[0][j]*alcom(0,j);
			return ans;
		}
		inline mtype alcom(size_type i,size_type j){
			return (i+j&1?-1:1)*com(i,j);
		}
		inline mtype com(size_type i,size_type j){return commat(i,j).det();}
		inline matrix_type commat(size_type i,size_type j){
			matrix_type ans(*this);
			if(i>=height||j>=width)throw "over_range";
			ans.erase(i,j);
			return ans;
		}
		inline matrix_type ortmat(){matrix_type ans;
			ans.vec.resize(height);
			for(reg size_type i=0;i<height;++i)
				for(reg size_type j=0;j<width;++j)
					ans[i].push_back(alcom(i,j));
			ans=ans.T();
			return ans;
		}
		inline void resize(size_type h,size_t w){height=h,width=w;
			vec.resize(height);
			for(reg size_type i=0;i<height;++i)vec[i].resize(width);
		}
		inline void erase(size_type i,size_type j){
			eraseh(i),erasew(j);
		}
		inline void eraseh(size_type i){vec.erase(vec.begin()+i);height--;}
		inline void erasew(size_type j){
			for(reg size_type ii=0;ii<height;++ii)
				vec[ii].erase(vec[ii].begin()+j);width--;
		}
		inline void each(mtype func(mtype x)){
			for(reg size_type i=0;i<height;++i)
				for(reg size_type j=0;j<width;++j)
					vec[i][j]=func(vec[i][j]);
		}
		inline void show(){
			for(reg size_type i=0;i<height;++i){
				for(reg size_type j=0;j<width;++j)
					std::cout<<vec[i][j]<<" ";puts("");}
		}
		inline matrix_type matrix_cast(mtype num,size_type s=2){
			matrix_type ans(s,s);for(reg size_type i=0;i<s;++i)ans[i][i]=num;
			return ans;
		}
		inline matrix_type pow(int n){
			if(n<0)return pow(-n).inv();
			--n;
			matrix_type ans(*this),tmp(*this);
			for(;n;)(n&1)&&(ans=ans*tmp,0),tmp=tmp*tmp,n>>=1;
			return ans;
		}
	protected:
		size_type height,width;
		static const mtype xx=(mtype)(1);
		static const mtype zero=xx-xx,one=xx/xx;
		std::vector<vec_type> vec;
};
#endif
