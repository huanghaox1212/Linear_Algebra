# Linear_Algebra
【模板】矩阵类模板
___
___
___
# Errors（2个未解决的错误）
大部分功能能够正常使用，以下函数暂时出现错误，正在努力改正。  
有dalao能够提供修改意见请私信我，万分感谢。  
【答案错误】伴随矩阵(ortmat函数)  
错误代码：
```cpp
inline matrix_type ortmat(){matrix_type ans;
	ans.vec.resize(height);
	for(reg size_type i=0;i<height;++i)
		for(reg size_type j=0;j<width;++j)
			ans[i].push_back(alcom(i,j));
	ans=ans.T();
	return ans;
}
```
其中$alcom(i,j)$是求$M_{ij}$的代数余子式。  
$T$函数是求矩阵的转置矩阵。

【运行错误】对角矩阵  
错误代码：
```cpp
virtual inline matrix_type matrix_cast(mtype num,size_type s=2){
	matrix_type ans(s,s);for(reg size_type i=0;i<s;++i)ans[i][i]=num;
}
```

# 正文
## 暂时支持的功能以及函数
- 矩阵加减乘（除法因为伴随矩阵错误暂不能使用），使用重载运算符+-\*/和数乘运算、数除运算。如9\*mat,mat/9,mat\*mat等。
- 快速幂（成员函数pow）
- 转置矩阵（成员函数T）
- 元素的余子式（成员函数commat）及其行列式（成员函数com）和代数余子式（成员函数alcom）
- 矩阵的行列式（成员函数det）
- 输入输出，具体见使用说明
- 删除行列和调整大小，具体见使用说明
- ==运算和!=运算
- 普通函数的求导运算（全局函数derivative）

代码见linear_algebra.h。
## 使用说明：
### 定义、初始化：
```cpp
int main(){
	int a[2][2]={{0,1},{1,0}};
	vector<vector<int> > vec;
	vector<int> vecint1(a[0],a[0]+2),vecint2(a[1],a[1]+2);//ok
	vec.push_back(vecint1);vec.push_back(vecint2);
	matrix<int> mat(vec),mat1(3,3,0),mat2(mat);
	return 0;
}
```
### 输入输出：
```cpp
cin>>mat;
cout<<mat;
mat.show();
```
**注：输出时结尾会换一行，有需要请自行修改。**

其它用法：
```cpp
int main(){
	if(mat==mat2)cout<<"equal!\n";
    else cout<<"not equal!\n";
	cout<<mat.alcom(0,0)<<endl;//求0行0列的元素的代数余子式
	show_var(mat.det());cout<<endl;//输出矩阵的行列式
    mat.resize(10,10);//行和列都调整为10
    mat.pow(10).show();//输出矩阵的10次方
    mat.eraseh(0);//删除第0行
    mat.erasew(1);//删除第1列
    mat.erase(1,2);//删除第1行、第2行
    //......
    //其他函数参见“暂时支持的功能以及函数”
	return 0;
}
```
