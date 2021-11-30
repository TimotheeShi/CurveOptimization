
#include <utility>
#include <vector>
#include <math.h>

/** 
 * ��ͷ�ļ����ڶ���������۰��߼��䵼���ĺ���
 */

/**
 * �ú���Ϊ���۰��ߵĲ�������
 * @params:
 *		phi: �������̵���λ�ǣ�
 *		R: ��ݷֲ�Բֱ����
 *		Z: ���߳��ֳ�����
 *		r: ���ֱ����
 *		e: ƫ�ľࡣ
 * @returns: 
 *		����һ���ɰ����ϵĵ�����꣨��������(x, y)����
 */
std::pair<double, double> cycloid_curve(const double& phi, const double& R, const double& Z, const double& r, const double& e) {
	double x = 0.0, y = 0.0;
	return std::pair<double, double>(x, y);
}

/**
 * �ú���Ϊ�����˰��ߵ�ģ�ֱ��R,r,e�󵼵ĵ���ֵ
 * @params:
 *		phi: �������̵���λ�ǣ�
 *		R: ��ݷֲ�Բֱ����
 *		Z: ���߳��ֳ�����
 *		r: ���ֱ����
 *		e: ƫ�ľࡣ
 * @returns:
 *		����һ����ά��������������ƫ������ֵ��[d��/dR, d��/dr, d��/de]��(��Ϊģ)
 */
std::vector<double> cycloid_drevative(const double& phi, const double& R, const double& Z, const double& r, const double& e) {
	double dR = 0.0, dr = 0.0, de = 0.0;
	std::pair<double, double> point = cycloid_curve(phi, R, Z, r, e);
	double row = sqrt(point.first * point.first + point.second * point.second);
	return std::vector<double> {dR, dr, de};
}