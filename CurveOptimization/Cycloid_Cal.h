
#include <utility>
#include <vector>
#include <math.h>

/** 
 * 此头文件用于定义计算理论摆线及其导数的函数
 */

/**
 * 该函数为理论摆线的参数方程
 * @params:
 *		phi: 参数方程的相位角；
 *		R: 针齿分布圆直径；
 *		Z: 摆线齿轮齿数；
 *		r: 针齿直径；
 *		e: 偏心距。
 * @returns: 
 *		返回一个由摆线上的点的坐标（有序数对(x, y)）。
 */
std::pair<double, double> cycloid_curve(const double& phi, const double& R, const double& Z, const double& r, const double& e) {
	double x = 0.0, y = 0.0;
	return std::pair<double, double>(x, y);
}

/**
 * 该函数为计算了摆线的模分别对R,r,e求导的导数值
 * @params:
 *		phi: 参数方程的相位角；
 *		R: 针齿分布圆直径；
 *		Z: 摆线齿轮齿数；
 *		r: 针齿直径；
 *		e: 偏心距。
 * @returns:
 *		返回一个三维的向量保存三个偏导数的值：[dρ/dR, dρ/dr, dρ/de]。(ρ为模)
 */
std::vector<double> cycloid_drevative(const double& phi, const double& R, const double& Z, const double& r, const double& e) {
	double dR = 0.0, dr = 0.0, de = 0.0;
	std::pair<double, double> point = cycloid_curve(phi, R, Z, r, e);
	double row = sqrt(point.first * point.first + point.second * point.second);
	return std::vector<double> {dR, dr, de};
}