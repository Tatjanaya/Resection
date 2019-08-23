// Resection.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <iomanip>

using namespace std;

void TransposeMatrix(double *m1, double *m2, int row, int col); //矩阵转置
void AddMatrix(double *m1, double *m2, int row, int col); //矩阵相加
void InverseMatrix(double *a, int n); //矩阵求逆
void MultiplyMatrix(double *m1, double *m2, double *result, int m, int n, int l); //矩阵相乘，其中result指向输出矩阵

struct Point
{
	int num;
	double x;
	double y;
};

int main()
{

	vector<double> vec; //像点坐标保存向量
	vector<double> vec1; //物方空间点坐标保存向量
	ifstream infile;
	infile.open("912.txt");
	if (!infile.is_open()) cout << "can't open this file!" << endl;

	double temp;
	while (infile >> temp) {
		vec.push_back(temp);
	}

	infile.close(); //将像点坐标保存在vec中
	ifstream infile1;
	infile1.open("wu91.txt");
	if (!infile1.is_open()) cout << "can't open this file" << endl;

	double temp1;
	while (infile1 >> temp1) {
		vec1.push_back(temp1); //将物方空间点坐标保存在vec1中
	}
	infile1.close();
	/*for (int i = 0; i < vec1.size(); i++) {
	cout << vec1[i] << " ";
	cout << endl;
	}*/

	//设置内外方位元素初始值
	double Xs = 600; double Ys = 900; double Zs = -150; //外方位线元素
	double phi = (25 * 3.14159265 / 180); double omega = 0; double kappa = 0; //外方位角元素
	double x0 = 0; double y0 = 0; double f = 28;  //内方位元素
	double k1 = 0.0, k2 = 0.0, p1 = 0.0, p2 = 0.0; //畸变改正数
	double x, y, X, Y, Z;
	double deltax, deltay;

	int num = vec1.size() / 4; //参与计算的点的个数
	double A[2 * 13] = { 0.0 }, AT[13 * 2] = { 0.0 }, ATA[13 * 13] = { 0.0 },
		NATA[13 * 13] = { 0.0 }, L[2 * 1] = { 0.0 }, ATL[13 * 1] = { 0.0 }, XX[13 * 1] = { 0.1 };

	double *LL = new double[2 * num];
	for (int i = 0; i < 2 * num; i++) {
		LL[i] = 0.0;
	}
	double Q = 0.0; //单位权中误差
	double a1, a2, a3, b1, b2, b3, c1, c2, c3;
	double Qii[13] = { 0.0 };
	double mXs, mYs, mZs, mphi, momega, mkappa;
	double mff, mxx, myy;
	double mk1, mk2, mp1, mp2;
	double d = 0; //迭代次数
				  //cout << vec.size() << endl;

	do {
		//定义旋转矩阵
		a1 = cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
		a2 = -cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);
		a3 = -sin(phi) * cos(omega);
		b1 = cos(omega) * sin(kappa);
		b2 = cos(omega) * cos(kappa);
		b3 = -sin(omega);
		c1 = sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
		c2 = -sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
		c3 = cos(phi) * cos(omega);

		double ATA1[13 * 13] = { 0.0 }, ATL1[13 * 1] = { 0.0 };

		for (int i = 0; i < num; i++) {
			//读数据
			x = vec[3 * i + 1] * 0.0051966;
			y = vec[3 * i + 2] * 0.0051966;
			Z = -1 * vec1[4 * i + 1];
			X = vec1[4 * i + 2];
			Y = vec1[4 * i + 3];

			double rr = (x - x0)*(x - x0) + (y - y0)*(y - y0);

			//共线方程
			double Xp = a1 * (X - Xs) + b1 * (Y - Ys) + c1 * (Z - Zs);
			double Yp = a2 * (X - Xs) + b2 * (Y - Ys) + c2 * (Z - Zs);
			double Zp = a3 * (X - Xs) + b3 * (Y - Ys) + c3 * (Z - Zs);

			//外方位元素改正数系数
			A[0] = 1.0 * (a1 * f + a3 * (x - x0)) / Zp;
			A[1] = 1.0 * (b1 * f + b3 * (x - x0)) / Zp;
			A[2] = 1.0 * (c1 * f + c3 * (x - x0)) / Zp;
			A[3] = (y - y0) * sin(omega) - ((x - x0) * ((x - x0) * cos(kappa) - (y - y0) * sin(kappa)) / f + f * cos(kappa)) * cos(omega);
			A[4] = -f * sin(kappa) - (x - x0) * ((x - x0) * sin(kappa) + (y - y0) * cos(kappa)) / f;
			A[5] = (y - y0);
			A[6] = (x - x0) / f;
			A[7] = 1;
			A[8] = 0;
			A[9] = -(x - x0) * rr;
			A[10] = -(x - x0) * rr * rr;
			A[11] = -rr - 2 * (x - x0) * (x - x0);
			A[12] = -2 * (x - x0) * (y - y0);

			A[13] = 1.0 * (a2 * f + a3 * (y - y0)) / Zp;
			A[14] = 1.0 * (b2 * f + b3 * (y - y0)) / Zp;
			A[15] = 1.0 * (c2 * f + c3 * (y - y0)) / Zp;
			A[16] = -(x - x0) * sin(omega) - ((y - y0) * ((x - x0) * cos(kappa) - (y - y0) * sin(kappa)) / f - f * sin(kappa)) * cos(omega);
			A[17] = -f * cos(kappa) - (y - y0) * ((x - x0) * sin(kappa) + (y - y0) * cos(kappa)) / f;
			A[18] = -(x - x0);
			A[19] = (y - y0) / f;
			A[20] = 0;
			A[21] = 1;
			A[22] = -(y - y0) * rr;
			A[23] = -(y - y0) * rr * rr;
			A[24] = -2 * (x - x0) * (y - y0);
			A[25] = -rr - 2 * (y - y0) * (y - y0);

			deltax = (x - x0) * (k1 * rr + k2 * rr * rr) + p1 * (rr + 2 * (x - x0) * (x - x0)) + 2 * p2 * (x - x0) * (y - y0);
			deltay = (y - y0) * (k1 * rr + k2 * rr * rr) + p2 * (rr + 2 * (y - y0) * (y - y0)) + 2 * p1 * (x - x0) * (y - y0);

			L[0] = x + f * Xp / Zp - x0 + deltax;
			L[1] = y + f * Yp / Zp - y0 + deltay;

			TransposeMatrix(A, AT, 2, 13);
			MultiplyMatrix(AT, A, ATA, 13, 2, 13);
			MultiplyMatrix(AT, L, ATL, 13, 2, 1);
			AddMatrix(ATA, ATA1, 13, 13);
			AddMatrix(ATL, ATL1, 13, 1);

			LL[2 * i] = L[0];
			LL[2 * i + 1] = L[1];
		}
		InverseMatrix(ATA1, 13);
		MultiplyMatrix(ATA1, ATL1, XX, 13, 13, 1);

		for (int i = 0; i < 13; i++) {
			Qii[i] = ATA1[13 * i + i];
		}

		Xs += XX[0];
		Ys += XX[1];
		Zs += XX[2];
		phi += XX[3];
		omega += XX[4];
		kappa += XX[5];
		f += XX[6];
		x0 += XX[7];
		y0 += XX[8];
		k1 += XX[9];
		k2 += XX[10];
		p1 += XX[11];
		p2 += XX[12];

		d = d + 1;

	} while (fabs(XX[0]) > 0.000001 || fabs(XX[1]) > 0.000001 || fabs(XX[2]) > 0.000001
		|| fabs(XX[3]) > 0.0000001 || fabs(XX[4]) > 0.0000001 || fabs(XX[5]) > 0.0000001
		|| fabs(XX[6]) > 0.000001 || fabs(XX[7]) > 0.000001 || fabs(XX[8]) > 0.000001
		|| fabs(XX[9]) > 0.000001 || fabs(XX[10]) > 0.000001 || fabs(XX[11]) > 0.000001 || fabs(XX[12]) > 0.000001);

	for (int i = 0; i < 2 * num; i++) {
		Q += LL[i] * LL[i];
	}
	Q = sqrt(Q / (2 * num - 6));

	mXs = sqrt(Qii[0])*Q;
	mYs = sqrt(Qii[1])*Q;
	mZs = sqrt(Qii[2])*Q;
	mphi = sqrt(Qii[3])*Q;
	momega = sqrt(Qii[4])*Q;
	mkappa = sqrt(Qii[5])*Q;
	mff = sqrt(Qii[6])*Q;
	mxx = sqrt(Qii[7])*Q;
	myy = sqrt(Qii[8])*Q;
	mk1 = sqrt(Qii[9])*Q;
	mk2 = sqrt(Qii[10])*Q;
	mp1 = sqrt(Qii[11])*Q;
	mp2 = sqrt(Qii[12])*Q;

	cout << "迭代次数为：" << d << endl;
	cout << "计算结果如下：" << endl;	//输出结果
	cout << "外方位元素：" << endl;
	cout << "      Xs   =  " << setprecision(10) << Xs << "        ±" << fixed << mXs << " mm" << endl;
	cout << "      Ys   =  " << Ys << "     ±" << mYs << " mm" << endl;
	cout << "      Zs   =  " << Zs << "   ±" << mZs << " mm" << endl;
	cout << "     phi   =  " << phi << "       ±" << mphi << " rad" << endl;
	cout << "   omega   =  " << omega << "       ±" << momega << " rad" << endl;
	cout << "   kappa   =  " << kappa << "       ±" << mkappa << " rad" << endl;
	cout << "内方位元素：" << endl;
	cout << "       f   =  " << f << "      ±" << mff << " mm" << endl;
	cout << "      x0   =  " << x0 << "     ±" << mxx << " pixel" << endl;
	cout << "      y0   =  " << y0 << "     ±" << myy << " pixel" << endl;
	cout << "畸变改正系数：" << endl;
	cout << "      k1   =  " << k1 << "       ±" << mk1 << endl;
	cout << "      k2   =  " << k2 << "       ±" << mk2 << endl;
	cout << "      p1   =  " << p1 << "       ±" << mp1 << endl;
	cout << "      p2   =  " << p2 << "       ±" << mp2 << endl;
	cout << "计算出的单位权中误差为：" << setprecision(10) << Q << "mm" << endl;
	cout << "旋转矩阵" << endl;
	cout << "  " << a1 << "  " << a2 << "  " << a3 << "  " << endl;
	cout << "  " << b1 << "  " << b2 << "  " << b3 << "  " << endl;
	cout << "  " << c1 << "  " << c2 << "  " << c3 << "  " << endl;
	cout << "各点的坐标残差" << endl;
	for (int i = 0; i < num; i++) {
		cout << int(vec[3 * i]) << "  " << LL[2 * i] << "  " << LL[2 * i + 1] << endl;
	}
	cout << endl;

	ofstream fileout("result.txt");
	if (!fileout)
		cout << "wrong" << endl;
	fileout << "迭代次数为：" << d << endl;
	fileout << "计算结果如下：" << endl;	//输出结果
	fileout << "外方位元素：" << endl;
	fileout << "      Xs   =  " << setprecision(10) << Xs << "        ±" << fixed << mXs << " mm" << endl;
	fileout << "      Ys   =  " << Ys << "     ±" << mYs << " mm" << endl;
	fileout << "      Zs   =  " << Zs << "   ±" << mZs << " mm" << endl;
	fileout << "     phi   =  " << phi << "       ±" << mphi << " rad" << endl;
	fileout << "   omega   =  " << omega << "       ±" << momega << " rad" << endl;
	fileout << "   kappa   =  " << kappa << "       ±" << mkappa << " rad" << endl;
	fileout << "内方位元素：" << endl;
	fileout << "       f   =  " << f << "      ±" << mff << " mm" << endl;
	fileout << "      x0   =  " << x0 << "     ±" << mxx << " pixel" << endl;
	fileout << "      y0   =  " << y0 << "     ±" << myy << " pixel" << endl;
	fileout << "畸变改正系数：" << endl;
	fileout << "      k1   =  " << k1 << "       ±" << mk1 << endl;
	fileout << "      k2   =  " << k2 << "       ±" << mk2 << endl;
	fileout << "      p1   =  " << p1 << "       ±" << mp1 << endl;
	fileout << "      p2   =  " << p2 << "       ±" << mp2 << endl;
	fileout << "计算出的单位权中误差为：" << setprecision(10) << Q << "mm" << endl;
	fileout << "旋转矩阵" << endl;
	fileout << "  " << a1 << "  " << a2 << "  " << a3 << "  " << endl;
	fileout << "  " << b1 << "  " << b2 << "  " << b3 << "  " << endl;
	fileout << "  " << c1 << "  " << c2 << "  " << c3 << "  " << endl;
	fileout << "各点的坐标残差" << endl;
	for (int i = 0; i < num; i++) {
		fileout << int(vec[3 * i]) << "  " << LL[2 * i] << "  " << LL[2 * i + 1] << endl;
	}
	fileout << endl;
	system("pause"); 
	return 0;
}

void TransposeMatrix(double *m1, double *m2, int row, int col) {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			m2[j * row + i] = m1[i * col + j];
		}
	}
}

void AddMatrix(double *m1, double *m2, int row, int col) {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			m2[i * col + j] = m2[i * col + j] + m1[i * col + j];
		}
	}
}

void InverseMatrix(double *a, int n) {
	int i, j, k;
	for (k = 0; k<n; k++){
		for (i = 0; i<n; i++){
			if (i != k)
				*(a + i*n + k) = -*(a + i*n + k) / (*(a + k*n + k));
		}
		*(a + k*n + k) = 1 / (*(a + k*n + k));
		for (i = 0; i<n; i++){
			if (i != k){
				for (j = 0; j<n; j++){
					if (j != k)
						*(a + i*n + j) += *(a + k*n + j)* *(a + i*n + k);
				}
			}
		}
		for (j = 0; j<n; j++){
			if (j != k)
				*(a + k*n + j) *= *(a + k*n + k);
		}
	}
}

void MultiplyMatrix(double *m1, double *m2, double *result, int m, int n, int l) {
	for (int i = 0; i<m; i++) {
		for (int j = 0; j<l; j++) {
			result[i*l + j] = 0.0;							//输出矩阵初始化
			for (int k = 0; k<n; k++)
				result[i*l + j] += m1[i*n + k] * m2[j + k*l];		//输出矩阵赋值
		}
	}
}

